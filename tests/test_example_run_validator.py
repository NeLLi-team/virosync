from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts.ci import validate_example_run


def _write_batch(output_root: Path, *, accepted: int = 2, elapsed: str = "9.5") -> None:
    fields = [
        "genome_id",
        "status",
        "benchmark_eligible",
        "legacy_resume",
        "predictions",
        "accepted",
        "high_tier",
        "medium_tier",
        "low_tier",
        "ncldv",
        "mirus",
        "ppv",
        "cress",
        "phage",
        "viral_unknown",
        "unknown",
        "total_bp",
        "genes",
        "hallmarks",
        "elapsed_sec",
        "error",
    ]
    row = {
        "genome_id": "example",
        "status": "success",
        "benchmark_eligible": "true",
        "legacy_resume": "false",
        "predictions": "3",
        "accepted": str(accepted),
        "high_tier": "1",
        "medium_tier": "0",
        "low_tier": "1",
        "ncldv": "1",
        "mirus": "0",
        "ppv": "1",
        "cress": "0",
        "phage": "0",
        "viral_unknown": "0",
        "unknown": "0",
        "total_bp": "42",
        "genes": "7",
        "hallmarks": "2",
        "elapsed_sec": elapsed,
        "error": "",
    }
    with (output_root / "batch_summary.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


@pytest.fixture
def valid_example(monkeypatch, tmp_path: Path):
    output_root = tmp_path / "out"
    run_dir = output_root / "example"
    run_dir.mkdir(parents=True)
    _write_batch(output_root)
    state_bytes = b'{"schema_version":3}\n'
    (run_dir / validate_example_run.RUN_STATE_FILENAME).write_bytes(state_bytes)
    result = {
        "accepted_bp": 42,
        "benchmark_eligible": True,
        "canonical_rows": 2,
        "detailed_rows": 3,
        "terminal_phase": None,
        "tier_counts": {"HIGH": 1, "MEDIUM": 0, "LOW": 1},
        "class_counts": {
            "NCLDV": 1,
            "MIRUS": 0,
            "PPV": 1,
            "CRESS": 0,
            "PHAGE": 0,
            "VIRAL_UNKNOWN": 0,
            "UNKNOWN": 0,
        },
    }
    artifact = SimpleNamespace(
        relative_path="virosync_predictions.tsv",
        size=12,
        sha256="b" * 64,
        schema="canonical-predictions-v4",
        row_count=2,
    )
    state = SimpleNamespace(
        schema_version=3,
        status="success",
        run_fingerprint="a" * 64,
        attempt=1,
        result=result,
        artifacts=(artifact,),
    )
    monkeypatch.setattr(validate_example_run, "load_run_state", lambda _path: state)
    monkeypatch.setattr(
        validate_example_run,
        "plan_resume",
        lambda *_args, **_kwargs: SimpleNamespace(completed=True, reason=None),
    )
    return output_root, state_bytes


def test_snapshot_is_deterministic_path_free_and_count_bound(valid_example) -> None:
    output_root, state_bytes = valid_example
    snapshot = validate_example_run.build_snapshot(output_root)
    rendered = json.dumps(snapshot, sort_keys=True)
    assert str(output_root) not in rendered
    run = snapshot["runs"]["example"]
    assert run["run_state_sha256"] == hashlib.sha256(state_bytes).hexdigest()
    assert run["result"]["canonical_rows"] == 2
    assert run["artifacts"][0]["sha256"] == "b" * 64


def test_snapshot_schema_error_names_current_version(tmp_path: Path) -> None:
    snapshot_path = tmp_path / "snapshot.json"
    snapshot_path.write_text('{"schema_version": 1}\n', encoding="utf-8")

    with pytest.raises(
        validate_example_run.ExampleValidationError,
        match=f"is not schema v{validate_example_run.SNAPSHOT_SCHEMA_VERSION}",
    ):
        validate_example_run._load_snapshot(snapshot_path)


def test_validator_rejects_batch_state_count_disagreement(valid_example) -> None:
    output_root, _state_bytes = valid_example
    _write_batch(output_root, accepted=1)
    with pytest.raises(
        validate_example_run.ExampleValidationError,
        match="batch accepted=1 differs from state 2",
    ):
        validate_example_run.build_snapshot(output_root)


def test_resume_comparison_requires_zero_elapsed_and_exact_snapshot(
    valid_example,
    tmp_path: Path,
) -> None:
    output_root, _state_bytes = valid_example
    snapshot_path = tmp_path / "snapshot.json"
    assert (
        validate_example_run.main(
            [str(output_root), "--write-snapshot", str(snapshot_path)]
        )
        == 0
    )

    _write_batch(output_root, elapsed="0")
    assert (
        validate_example_run.main(
            [
                str(output_root),
                "--compare-snapshot",
                str(snapshot_path),
                "--require-resume",
            ]
        )
        == 0
    )

    payload = json.loads(snapshot_path.read_text())
    payload["runs"]["example"]["artifacts"][0]["sha256"] = "c" * 64
    snapshot_path.write_text(json.dumps(payload), encoding="utf-8")
    assert (
        validate_example_run.main(
            [
                str(output_root),
                "--compare-snapshot",
                str(snapshot_path),
                "--require-resume",
            ]
        )
        == 1
    )


def test_validator_enforces_expected_example_totals(valid_example) -> None:
    output_root, _state_bytes = valid_example
    assert (
        validate_example_run.main(
            [
                str(output_root),
                "--expect-predictions",
                "3",
                "--expect-accepted",
                "2",
            ]
        )
        == 0
    )
    assert (
        validate_example_run.main(
            [str(output_root), "--expect-accepted", "3"]
        )
        == 1
    )
