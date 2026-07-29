from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from virosync.config import MaskingConfig
from virosync.orchestration._flows.single_genome.manifest import (
    _write_completion_manifest,
    _write_empty_run_log,
)
from virosync.orchestration._flows.single_genome.orchestrator import (
    _revalidate_completed_run,
)
from virosync.orchestration._flows.single_genome import orchestrator
from virosync.orchestration._flows.single_genome.resume import (
    _completed_run_artifacts,
)
from virosync.orchestration._flows.single_genome.run_state import load_run_state
from virosync.orchestration._flows.single_genome.phase1_state import (
    write_phase1_state,
)
from virosync.orchestration._flows.single_genome.phase2_resume_state import (
    write_phase2_resume_state,
)
from virosync.orchestration._flows.single_genome.phase_state import (
    write_phase2_state,
)
from virosync.validation import tsv_invariants
from virosync.validation.tsv_invariants import (
    InvariantIssue,
    InvariantReport,
    TSVInvariantError,
    enforce_tsv_invariants,
)
from virosync.pipeline.phase0.masking import mask_genome_pipeline
from virosync.pipeline.host_signatures import HostSignatureModel
from virosync.pipeline.phase1.seed_merger import MergedSeed
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary


def _report_summary(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return next(reader)


def _seed_success_markers(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "run.log").write_text(
        "# ViroSync Run Log: demo\n\n## Results Summary\nGEVEs detected: 0\n"
    )
    _seed_masking_status(output_dir)
    _write_completion_manifest(output_dir, genome_id="demo", status="success")


def _seed_masking_status(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    input_fasta = output_dir / "input.fna"
    input_fasta.write_text(">demo\nACGT\n")
    mask_genome_pipeline(
        input_fasta,
        output_dir / "phase0" / "masking",
        config=MaskingConfig(),
    )


def test_fatal_invariant_writes_report_before_raising(tmp_path: Path) -> None:
    detailed = tmp_path / "virosync_predictions_detailed.tsv"
    detailed.write_text(
        "eve_id\ttotal_proteins\tncldv_top10_proteins\n"
        "EVE_1\t1\t2\n"
    )
    report_path = tmp_path / "virosync_tsv_invariant_report.tsv"

    with pytest.raises(TSVInvariantError, match="ncldv_top_count_out_of_range"):
        enforce_tsv_invariants(detailed, report_path)

    summary = _report_summary(report_path)
    assert summary["status"] == "FAIL"
    assert int(summary["error_count"]) >= 1
    assert summary["warning_count"] == "0"


def test_missing_detailed_tsv_writes_typed_failure_report(tmp_path: Path) -> None:
    report_path = tmp_path / "virosync_tsv_invariant_report.tsv"

    with pytest.raises(TSVInvariantError, match="detailed_tsv_missing"):
        enforce_tsv_invariants(tmp_path / "missing.tsv", report_path)

    assert _report_summary(report_path) == {
        "status": "FAIL",
        "rows_checked": "0",
        "issue_count": "1",
        "error_count": "1",
        "warning_count": "0",
    }
    assert "detailed_tsv_missing" in report_path.read_text()


def test_warning_only_report_is_nonfatal_and_distinct(
    tmp_path: Path,
    monkeypatch,
) -> None:
    detailed = tmp_path / "virosync_predictions_detailed.tsv"
    detailed.write_text("eve_id\n")
    warning_report = InvariantReport(
        rows_checked=0,
        issues=[
            InvariantIssue(
                eve_id=".",
                check="diagnostic_warning",
                severity="warning",
                message="warning-only diagnostic",
            )
        ],
    )
    monkeypatch.setattr(
        tsv_invariants,
        "run_tsv_invariant_checks",
        lambda **kwargs: warning_report,
    )
    report_path = tmp_path / "virosync_tsv_invariant_report.tsv"

    result = enforce_tsv_invariants(detailed, report_path)

    assert result.passed is True
    assert result.error_count == 0
    assert result.warning_count == 1
    assert _report_summary(report_path)["status"] == "PASS_WITH_WARNINGS"


def test_unknown_severity_fails_closed() -> None:
    report = InvariantReport(
        rows_checked=0,
        issues=[InvariantIssue(".", "future_check", "notice", "unknown severity")],
    )

    assert report.passed is False
    assert report.error_count == 1
    assert report.warning_count == 0


@pytest.mark.parametrize(
    "reason",
    [
        "no HMM hits",
        "no validated markers",
        "no seeds",
        "no refined boundaries",
    ],
)
def test_all_zero_result_reasons_validate_before_success(
    tmp_path: Path,
    reason: str,
) -> None:
    output_dir = tmp_path / reason.replace(" ", "_")
    output_dir.mkdir()
    (output_dir / "virosync_predictions.tsv").write_text("eve_id\n")
    (output_dir / "virosync_predictions_detailed.tsv").write_text("eve_id\n")
    _seed_masking_status(output_dir)

    _write_empty_run_log(
        output_dir=output_dir,
        genome_id="demo",
        reason=reason,
        elapsed_sec=0.0,
        fingerprint="zero-fp",
    )

    report_path = output_dir / "virosync_tsv_invariant_report.tsv"
    assert _report_summary(report_path) == {
        "status": "PASS",
        "rows_checked": "0",
        "issue_count": "0",
        "error_count": "0",
        "warning_count": "0",
    }
    manifest = json.loads((output_dir / "virosync_run_complete.json").read_text())
    assert manifest["status"] == "success"
    assert manifest["reason"] == reason
    assert (
        _completed_run_artifacts(
            output_dir,
            expected_fingerprint="zero-fp",
        )
        is None
    )
    assert (
        _completed_run_artifacts(
            output_dir,
            expected_fingerprint="zero-fp",
            allow_legacy_schema=True,
        )
        is not None
    )


def test_zero_result_missing_detailed_removes_stale_success_markers(
    tmp_path: Path,
) -> None:
    _seed_success_markers(tmp_path)

    with pytest.raises(TSVInvariantError, match="detailed_tsv_missing"):
        _write_empty_run_log(
            output_dir=tmp_path,
            genome_id="demo",
            reason="no HMM hits",
            elapsed_sec=0.0,
        )

    assert (tmp_path / "virosync_tsv_invariant_report.tsv").exists()
    assert not (tmp_path / "run.log").exists()
    assert not (tmp_path / "virosync_run_complete.json").exists()


def test_cached_success_is_revalidated_and_invalidated(tmp_path: Path) -> None:
    _seed_success_markers(tmp_path)
    detailed = tmp_path / "virosync_predictions_detailed.tsv"
    detailed.write_text(
        "eve_id\ttotal_proteins\tncldv_top10_proteins\n"
        "EVE_1\t1\t2\n"
    )
    artifacts = {"predictions_detailed": detailed}

    with pytest.raises(TSVInvariantError):
        _revalidate_completed_run(tmp_path, artifacts)

    assert (tmp_path / "virosync_tsv_invariant_report.tsv").exists()
    assert not (tmp_path / "run.log").exists()
    assert not (tmp_path / "virosync_run_complete.json").exists()


class _NullResourceMonitor:
    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        return False


def test_fresh_normal_path_enforces_invariants_before_success_markers(
    tmp_path: Path,
    monkeypatch,
) -> None:
    """Exercise the normal finalization block, not its low-level helper alone."""
    genome = tmp_path / "demo.fna"
    genome.write_text(">scaffold\nACGT\n")
    output_dir = tmp_path / "results" / "demo"
    monkeypatch.setattr(
        orchestrator,
        "_completed_run_artifacts",
        lambda *args, **kwargs: None,
    )

    def _fake_phase0(**kwargs):
        phase0_dir = output_dir / "phase0"
        masking_result = mask_genome_pipeline(
            genome,
            phase0_dir / "masking",
            config=MaskingConfig(),
        )
        proteome = phase0_dir / "proteome.fasta"
        genes = phase0_dir / "genes.gff"
        proteome.write_text("")
        genes.write_text("")
        return {
            "masked_path": genome,
            "repeat_regions": [],
            "proteome_path": proteome,
            "n_genes": 0,
            "elapsed": 0.0,
            "masking_result": masking_result,
        }

    monkeypatch.setattr(orchestrator, "_run_phase0_subflow", _fake_phase0)

    def _fake_phase1(**kwargs):
        host_signature_model = HostSignatureModel()
        merged_seeds = [
            MergedSeed(
                scaffold="scaffold",
                start=0,
                end=4,
                seed_id="seed-0-scaffold-0",
                sources=["test"],
            )
        ]
        write_phase1_state(
            output_dir / "phase1" / "resume_state.json",
            validated_markers=[],
            merged_seeds=merged_seeds,
            host_signature_model=host_signature_model,
            host_signatures=set(),
            host_deviation_summary=None,
        )
        return {
            "merged_seeds": merged_seeds,
            "validated_markers": [],
            "host_signature_model": host_signature_model,
            "host_signatures": {},
            "background": None,
            "gene_data": {},
            "host_deviation_summary": {},
            "elapsed": 0.0,
        }

    monkeypatch.setattr(orchestrator, "_run_phase1_subflow", _fake_phase1)

    def _fake_phase2(**kwargs):
        boundary = RefinedBoundary(
            scaffold="scaffold",
            start=0,
            end=4,
            seed_id="seed-0-scaffold-0",
            seed_sources=["test"],
            confidence=0.9,
        )
        phase2_dir = output_dir / "phase2"
        boundaries_bed = phase2_dir / "refined_boundaries.bed"
        write_phase2_state(phase2_dir / "refined_state.json", [boundary])
        write_phase2_resume_state(
            phase2_dir / "resume_state.json",
            refined_boundaries=[boundary],
            boundary_taxonomy_map={},
            boundary_control_stats=None,
            boundary_diamond_query=None,
        )
        boundaries_bed.write_text(
            "scaffold\t0\t4\tEVE_scaffold_0-4\t900\t.\n"
        )
        return {
            "refined_boundaries": [boundary],
            "boundary_taxonomy_map": {},
            "boundary_control_stats": None,
            "boundary_diamond_query": None,
            "proteome_index": {},
            "goto_phase3": True,
            "boundaries_bed": boundaries_bed,
            "elapsed": 0.0,
        }

    monkeypatch.setattr(orchestrator, "_run_phase2_subflow", _fake_phase2)
    monkeypatch.setattr(
        orchestrator,
        "_run_phase3_subflow",
        lambda **kwargs: {
            "verification_results": [],
            "accepted_results": [],
            "tier_counts": {"HIGH": 0, "MEDIUM": 0, "LOW": 0},
            "candidate_tier_counts": {"HIGH": 0, "MEDIUM": 0, "LOW": 0},
            "classification_stats": {},
            "accepted": 0,
            "accepted_bp": 0,
            "total_genes": 0,
            "total_hallmarks": 0,
            "quality_gate_dropped": 0,
            "elapsed": 0.0,
            "precomputed_tmvec": None,
        },
    )
    monkeypatch.setattr(orchestrator, "_generate_required_reports", lambda **kwargs: {})

    def _fake_call_task(task_fn, *args, **kwargs):
        if task_fn is orchestrator.generate_outputs_task:
            detailed = (
                Path(kwargs["output_dir"])
                / "virosync_predictions_detailed.tsv"
            )
            detailed.parent.mkdir(parents=True, exist_ok=True)
            detailed.write_text(
                "eve_id\ttotal_proteins\tncldv_top10_proteins\n"
                "EVE_1\t1\t2\n"
            )
            return {}
        return None

    monkeypatch.setattr(orchestrator, "call_task", _fake_call_task)
    from virosync.orchestration import resource_monitor
    from virosync.utils import gpu, provenance

    monkeypatch.setattr(resource_monitor, "ResourceMonitor", _NullResourceMonitor)
    monkeypatch.setattr(gpu, "release_gpu_memory", lambda: None)
    monkeypatch.setattr(provenance, "write_provenance", lambda *args, **kwargs: None)

    with pytest.raises(TSVInvariantError, match="ncldv_top_count_out_of_range"):
        orchestrator.single_genome_flow(
            genome_path=genome,
            output_dir=output_dir,
            genome_id="demo",
            device="cpu",
            resume=True,
        )

    assert (output_dir / "virosync_tsv_invariant_report.tsv").exists()
    assert not (output_dir / "run.log").exists()
    assert not (output_dir / "virosync_run_complete.json").exists()
    assert not (output_dir / "phase3.complete.json").exists()
    assert load_run_state(output_dir).status == "failed"
