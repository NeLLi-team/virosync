from __future__ import annotations

import csv
from pathlib import Path
import pytest
from click.testing import CliRunner

from virosync.config import ApplicationConfig, FeatureResolution
from virosync.orchestration import cli as orchestration_cli
from virosync.orchestration import python_runner
from virosync.orchestration.cli import orchestrate
from virosync.validation.tsv_invariants import (
    InvariantIssue,
    InvariantReport,
    TSVInvariantError,
)


def _invoke_batch(
    tmp_path: Path,
    monkeypatch,
    outcomes: dict[str, object],
):
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    for genome_id in outcomes:
        (input_dir / f"{genome_id}.fna").write_text(">scaffold\nACGT\n")

    output_root = tmp_path / "results"

    def _fake_single(*, genome_path, output_dir, genome_id, config):
        output_dir.mkdir(parents=True, exist_ok=True)
        outcome = outcomes[genome_id]
        if isinstance(outcome, BaseException):
            raise outcome
        if outcome:
            (output_dir / "successful-output.txt").write_text("preserved\n")
            return {
                "genome_id": genome_id,
                "success": True,
                "benchmark_eligible": outcome != "ineligible",
                "legacy_resume": outcome == "legacy",
                "accepted": 0,
                "predictions": 0,
                "elapsed_sec": 0.0,
            }
        return {
            "genome_id": genome_id,
            "success": False,
            "error": f"injected failure for {genome_id}",
            "accepted": 0,
            "predictions": 0,
            "elapsed_sec": 0.0,
        }

    monkeypatch.setattr(python_runner, "_single_genome_callable", lambda: _fake_single)
    application = ApplicationConfig.from_dict({"schema_version": 1})
    monkeypatch.setattr(
        orchestration_cli,
        "_load_config",
        lambda path: application,
    )
    monkeypatch.setattr(
        orchestration_cli,
        "_resolve_pipeline_resources",
        lambda config, orchestration, path: config,
    )
    monkeypatch.setattr(
        orchestration_cli,
        "_resolve_optional_features",
        lambda config: (
            config,
            {
                name: FeatureResolution(False, False, False)
                for name in ("boltz", "tmvec", "interproscan")
            },
        ),
    )
    monkeypatch.setattr(
        orchestration_cli,
        "_validate_runtime_config",
        lambda config: None,
    )

    result = CliRunner().invoke(
        orchestrate,
        [
            "run",
            "-i",
            str(input_dir),
            "-o",
            str(output_root),
            "--max-concurrent-genomes",
            "2",
        ],
    )
    return result, output_root


def _summary_rows(output_root: Path) -> list[dict[str, str]]:
    with (output_root / "batch_summary.tsv").open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_partial_failure_writes_complete_summaries_and_exits_one(
    tmp_path: Path,
    monkeypatch,
) -> None:
    result, output_root = _invoke_batch(
        tmp_path,
        monkeypatch,
        {"good": True, "bad": False},
    )

    assert result.exit_code == 1, result.output
    assert "Batch Processing Failed" in result.output
    assert "Batch Processing Complete" not in result.output
    assert str(output_root / "batch_summary.tsv") in result.output
    assert str(output_root / "batch_report.md") in result.output
    success_output = output_root / "good" / "successful-output.txt"
    assert success_output.read_text() == "preserved\n"
    rows = _summary_rows(output_root)
    assert {row["genome_id"]: row["status"] for row in rows} == {
        "bad": "failed",
        "good": "success",
    }
    report = (output_root / "batch_report.md").read_text()
    assert "2 (1 successful, 1 failed)" in report
    assert "injected failure for bad" in report


def test_all_failed_writes_one_row_per_input_and_exits_one(
    tmp_path: Path,
    monkeypatch,
) -> None:
    outcomes = {"failed_a": False, "failed_b": False}
    result, output_root = _invoke_batch(tmp_path, monkeypatch, outcomes)

    assert result.exit_code == 1, result.output
    rows = _summary_rows(output_root)
    assert [row["genome_id"] for row in rows] == sorted(outcomes)
    assert [row["status"] for row in rows] == ["failed", "failed"]
    assert (output_root / "batch_report.md").exists()
    assert "2/2 genomes failed" in result.output


def test_benchmark_ineligible_success_is_reported_and_cli_fails_closed(
    tmp_path: Path,
    monkeypatch,
) -> None:
    result, output_root = _invoke_batch(
        tmp_path,
        monkeypatch,
        {"fallback": "ineligible"},
    )

    assert result.exit_code == 1, result.output
    assert "Batch Processing Completed with Warnings" in result.output
    assert "1/1 successful genomes are benchmark-ineligible" in result.output
    rows = _summary_rows(output_root)
    assert len(rows) == 1
    assert rows[0]["genome_id"] == "fallback"
    assert rows[0]["status"] == "success_with_warnings"
    assert rows[0]["benchmark_eligible"] == "false"
    assert rows[0]["legacy_resume"] == "false"
    report = (output_root / "batch_report.md").read_text()
    assert "1 success with warnings" in report
    assert "| fallback | success_with_warnings | no | no |" in report


def test_invariant_error_propagates_through_batch_to_cli_exit_one(
    tmp_path: Path,
    monkeypatch,
) -> None:
    error = TSVInvariantError(
        InvariantReport(
            rows_checked=1,
            issues=[InvariantIssue("EVE_1", "broken", "error", "injected")],
        ),
        tmp_path / "results" / "fatal" / "virosync_tsv_invariant_report.tsv",
    )

    result, output_root = _invoke_batch(
        tmp_path,
        monkeypatch,
        {"fatal": error},
    )

    assert result.exit_code == 1, result.output
    rows = _summary_rows(output_root)
    assert len(rows) == 1
    assert rows[0]["status"] == "failed"
    assert "Detailed TSV invariant check failed" in rows[0]["error"]
    assert "1/1 genomes failed" in result.output


@pytest.mark.parametrize(
    "outcomes",
    [
        {"zero_call": True},
        {"success_a": True, "success_b": True},
    ],
)
def test_all_success_including_zero_call_exits_zero(
    tmp_path: Path,
    monkeypatch,
    outcomes: dict[str, bool],
) -> None:
    result, output_root = _invoke_batch(tmp_path, monkeypatch, outcomes)

    assert result.exit_code == 0, result.output
    assert "Batch Processing Complete" in result.output
    assert "Batch Processing Failed" not in result.output
    assert all(row["status"] == "success" for row in _summary_rows(output_root))
    assert all(
        row["benchmark_eligible"] == "true"
        and row["legacy_resume"] == "false"
        for row in _summary_rows(output_root)
    )
    assert (output_root / "batch_report.md").exists()
