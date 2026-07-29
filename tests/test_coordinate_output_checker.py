from __future__ import annotations

import json
import runpy
from pathlib import Path
import subprocess
import sys

import pytest

from virosync.output_contract import coordinate_contract_metadata


CHECKER_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "ci"
    / "check_coordinate_outputs.py"
)
CHECKER = runpy.run_path(str(CHECKER_PATH))
check_coordinate_output_roots = CHECKER["check_coordinate_output_roots"]
main = CHECKER["main"]


def test_checker_script_resolves_source_without_pythonpath(tmp_path: Path) -> None:
    result = subprocess.run(
        [sys.executable, str(CHECKER_PATH), str(tmp_path / "missing")],
        cwd=tmp_path,
        check=False,
        capture_output=True,
        text=True,
        env={"PATH": ""},
    )
    assert result.returncode == 1
    assert "no completion or summary metadata" in result.stderr
    assert "ModuleNotFoundError" not in result.stderr


def _write_metadata(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload) + "\n")


def test_checker_accepts_current_nested_completion_and_summary_metadata(
    tmp_path: Path,
) -> None:
    contract = coordinate_contract_metadata()
    for result_name in ("result_a", "result_b"):
        result_dir = tmp_path / result_name
        _write_metadata(
            result_dir / "virosync_run_complete.json",
            {"status": "success", **contract},
        )
        _write_metadata(
            result_dir / "phase3_synthesis" / "virosync_summary.json",
            {"statistics": {}, **contract},
        )
        predictions = (
            result_dir / "phase3_synthesis" / "virosync_predictions.bed"
        )
        predictions.write_text("")

    assert check_coordinate_output_roots([tmp_path]) == []


def test_checker_accepts_current_completion_only_for_early_exit(
    tmp_path: Path,
) -> None:
    _write_metadata(
        tmp_path / "virosync_run_complete.json",
        {"status": "success", **coordinate_contract_metadata()},
    )

    assert check_coordinate_output_roots([tmp_path]) == []


def test_checker_accepts_authoritative_schema3_terminal_zero(tmp_path: Path) -> None:
    _write_metadata(
        tmp_path / "virosync_run_state.json",
        {
            "schema_version": 3,
            "status": "success",
            "identities": coordinate_contract_metadata(),
            "result": {"terminal_phase": 1, "canonical_rows": 0},
        },
    )

    assert check_coordinate_output_roots([tmp_path]) == []


def test_checker_rejects_non_success_schema3_state(tmp_path: Path) -> None:
    _write_metadata(
        tmp_path / "virosync_run_state.json",
        {
            "schema_version": 3,
            "status": "failed",
            "identities": coordinate_contract_metadata(),
        },
    )

    errors = check_coordinate_output_roots([tmp_path])

    assert any("not schema-v3 success" in error for error in errors)


@pytest.mark.parametrize(
    ("filename", "payload"),
    [
        ("virosync_run_complete.json", {"status": "success"}),
        ("virosync_summary.json", {"statistics": {}}),
    ],
)
def test_checker_rejects_legacy_metadata_and_requires_clean_regeneration(
    tmp_path: Path,
    filename: str,
    payload: dict,
) -> None:
    _write_metadata(tmp_path / filename, payload)

    errors = check_coordinate_output_roots([tmp_path])

    assert any("coordinate_schema_version" in error for error in errors)
    assert any("--clean-run" in error for error in errors)


@pytest.mark.parametrize(
    ("field", "stale_value"),
    [
        ("coordinate_schema_version", 1),
        ("output_schema_version", 1),
        ("coordinate_convention", "1-based, closed [start, end]"),
    ],
)
def test_checker_rejects_each_stale_contract_field(
    tmp_path: Path,
    field: str,
    stale_value: int | str,
) -> None:
    payload = coordinate_contract_metadata()
    payload[field] = stale_value
    _write_metadata(tmp_path / "virosync_run_complete.json", payload)

    errors = check_coordinate_output_roots([tmp_path])

    assert any(field in error for error in errors)


def test_checker_rejects_corrupt_or_absent_metadata(tmp_path: Path) -> None:
    empty_root = tmp_path / "empty"
    empty_root.mkdir()
    corrupt_root = tmp_path / "corrupt"
    corrupt_root.mkdir()
    (corrupt_root / "virosync_summary.json").write_text("{not-json")

    errors = check_coordinate_output_roots([empty_root, corrupt_root])

    assert any(
        "no completion or summary metadata" in error for error in errors
    )
    assert any("invalid JSON" in error for error in errors)


@pytest.mark.parametrize(
    "old_relative",
    [
        Path("old/phase3_synthesis/virosync_predictions.bed"),
        Path("old/virosync_predictions_detailed.tsv"),
        Path("old/phase1/hhg/hhg_seeds.bed"),
        Path("old/phase1/region_assembly/marker_seed_regions.bed"),
    ],
)
def test_checker_rejects_old_output_subtree_beside_current_result(
    tmp_path: Path,
    old_relative: Path,
) -> None:
    current = tmp_path / "current"
    _write_metadata(
        current / "virosync_run_complete.json",
        {"status": "success", **coordinate_contract_metadata()},
    )
    _write_metadata(
        current / "phase3_synthesis" / "virosync_summary.json",
        {"statistics": {}, **coordinate_contract_metadata()},
    )
    old_output = tmp_path / old_relative
    old_output.parent.mkdir(parents=True)
    old_output.write_text("ctg\t0\t10\n")

    errors = check_coordinate_output_roots([tmp_path])

    assert any(
        str(old_output) in error and "completion metadata" in error
        for error in errors
    )


def test_checker_requires_summary_for_normal_coordinate_outputs(
    tmp_path: Path,
) -> None:
    _write_metadata(
        tmp_path / "virosync_run_complete.json",
        {"status": "success", **coordinate_contract_metadata()},
    )
    predictions = tmp_path / "phase3_synthesis" / "virosync_predictions.gff3"
    predictions.parent.mkdir(parents=True)
    predictions.write_text("##gff-version 3\n")

    errors = check_coordinate_output_roots([tmp_path])

    assert any("virosync_summary.json" in error for error in errors)


def test_checker_cli_reports_failure(tmp_path: Path, capsys) -> None:
    _write_metadata(
        tmp_path / "virosync_run_complete.json",
        {"status": "success"},
    )

    exit_code = main([str(tmp_path)])

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "--clean-run" in captured.err
