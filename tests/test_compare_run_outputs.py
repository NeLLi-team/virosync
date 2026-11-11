"""Compare the scientific outputs of two genome output folders."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts import compare_run_outputs as compare

PREDICTIONS = "phase3_synthesis/virosync_predictions.tsv"
DETAILED = "virosync_predictions_detailed.tsv"
SUMMARY = "phase3_synthesis/virosync_summary.json"
MODEL = "phase1/marker_validation/host_signature_model.json"


def _write(root: Path, relative: str, content: str) -> None:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _folder(root: Path, *, generated_at: str = "2026-11-10T10:00:00") -> Path:
    _write(root, PREDICTIONS, "eve_id\tstart\nEVE_1\t10\n")
    _write(root, DETAILED, "eve_id\tstart\nEVE_1\t10\nEVE_2\t50\n")
    _write(root, SUMMARY, json.dumps({"generated_at": generated_at, "statistics": {"canonical_predictions": 1}}))
    _write(root, MODEL, json.dumps({"token_weights": {"a": 1.5, "b": 2.5}, "token_bits": {"a": [1.5, 0.5]}}))
    _write(root, "run.log", f"# Generated: {generated_at}\n")
    return root


def test_identical_folders_have_no_difference(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")

    assert compare.compare_runs(first, second) == []


def test_main_prints_and_exits_zero_for_identical_folders(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")

    assert compare.main([str(first), str(second)]) == 0

    assert "identical" in capsys.readouterr().out


def test_a_table_difference_names_the_file_and_the_first_differing_line(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    _write(second, DETAILED, "eve_id\tstart\nEVE_1\t10\nEVE_2\t51\n")

    differences = compare.compare_runs(first, second)

    assert len(differences) == 1
    assert DETAILED in differences[0]
    assert "line 3" in differences[0]
    assert "EVE_2\t50" in differences[0]
    assert "EVE_2\t51" in differences[0]


def test_main_exits_one_and_prints_each_difference(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    _write(second, DETAILED, "changed\n")
    _write(second, PREDICTIONS, "changed\n")

    assert compare.main([str(first), str(second)]) == 1

    output = capsys.readouterr().out
    assert DETAILED in output
    assert PREDICTIONS in output


def test_the_time_stamp_of_a_json_file_is_ignored_and_any_other_key_is_not(tmp_path: Path) -> None:
    first = _folder(tmp_path / "a", generated_at="2026-11-10T10:00:00")
    second = _folder(tmp_path / "b", generated_at="2026-11-11T12:30:00")
    assert compare.compare_runs(first, second) == []

    _write(second, SUMMARY, json.dumps({"generated_at": "x", "statistics": {"canonical_predictions": 2}}))
    differences = compare.compare_runs(first, second)

    assert len(differences) == 1
    assert "statistics.canonical_predictions" in differences[0]


def test_json_key_order_is_ignored_but_list_order_is_not(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    _write(second, MODEL, json.dumps({"token_bits": {"a": [1.5, 0.5]}, "token_weights": {"b": 2.5, "a": 1.5}}))
    assert compare.compare_runs(first, second) == []

    _write(second, MODEL, json.dumps({"token_weights": {"a": 1.5, "b": 2.5}, "token_bits": {"a": [0.5, 1.5]}}))
    differences = compare.compare_runs(first, second)

    assert len(differences) == 2
    assert "token_bits.a[0]" in differences[0]
    assert "token_bits.a[1]" in differences[1]


def test_a_file_in_one_folder_only_is_a_difference(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    (second / DETAILED).unlink()

    differences = compare.compare_runs(first, second)

    assert len(differences) == 1
    assert "only in" in differences[0]
    assert DETAILED in differences[0]


def test_the_same_table_in_another_location_is_a_difference(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    (second / PREDICTIONS).rename(second / "virosync_predictions.tsv")

    differences = compare.compare_runs(first, second)

    assert len(differences) == 1
    assert "virosync_predictions.tsv" in differences[0]
    assert "location" in differences[0]


def test_files_outside_the_scientific_set_are_not_compared(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    _write(second, "run.log", "# Generated: another day\n")
    _write(second, "phase1/reuse_provenance.json", "{}")

    assert compare.compare_runs(first, second) == []


def test_a_missing_folder_exits_two(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    first = _folder(tmp_path / "a")

    assert compare.main([str(first), str(tmp_path / "missing")]) == 2

    assert "missing" in capsys.readouterr().err
