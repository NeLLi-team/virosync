from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import pytest

from virosync.ablation import (
    ABLATION_CONTRACT_SHA256,
    AblationCounters,
    AblationEvents,
    AblationID,
)

from virosync.config import MaskingConfig
from virosync.orchestration import python_runner
from virosync.orchestration._flows.single_genome import (
    _completed_run_artifacts,
    _require_phase2b_gene_taxonomy_db,
    _summarize_predictions_tsv,
    _write_completion_manifest,
)
from virosync.orchestration._flows.single_genome.manifest import (
    _summarize_prediction_outputs,
)
from virosync.orchestration._flows.single_genome.loaders import (
    _build_merged_seeds_from_regions,
)
from virosync.orchestration._flows.single_genome.orchestrator import (
    _masking_request_identity,
)
from virosync.orchestration._flows.single_genome.phase1_state import (
    PHASE1_STATE_SCHEMA,
    write_phase1_state,
)
from virosync.orchestration._flows.single_genome.phase2_resume_state import (
    PHASE2_RESUME_STATE_SCHEMA,
    write_phase2_resume_state,
)
from virosync.orchestration._flows.single_genome.phase_state import (
    PHASE2_STATE_SCHEMA,
    write_phase2_state,
)
from virosync.orchestration._flows.single_genome.run_state import (
    build_artifact_identity,
    build_input_identity,
    canonical_sha256,
    compute_run_fingerprint,
    marker_sha256,
    publish_phase_completion,
    publish_run_started,
    publish_run_success,
)
from virosync.output_contract import (
    COORDINATE_CONVENTION,
    COORDINATE_SCHEMA_VERSION,
    EFFECTIVE_EVE_CLASS_COUNT_KEYS,
    OUTPUT_SCHEMA_VERSION,
    effective_eve_class_count_total,
)
from virosync.pipeline.host_signatures import HostSignatureModel
from virosync.pipeline.phase0.masking import mask_genome_pipeline
from virosync.pipeline.phase1.seed_merger import MergedSeed
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary


_ZERO_PREDICTION_HEADER = (
    "eve_id\tscaffold\tstart\tend\tlength\tconfidence_tier\t"
    "final_confidence\teffective_eve_class\n"
)


def _seed_masking_status(output_dir: Path) -> None:
    input_fasta = output_dir / "input.fna"
    input_fasta.parent.mkdir(parents=True, exist_ok=True)
    input_fasta.write_text(">demo\nACGT\n")
    mask_genome_pipeline(
        input_fasta,
        output_dir / "phase0" / "masking",
        config=MaskingConfig(),
    )


def _closed_schema3_identity(
    output_dir: Path,
    *,
    identity_seed: str = "demo",
) -> dict[str, object]:
    """Build the complete deterministic identity used by resume fixtures."""
    output_dir.mkdir(parents=True, exist_ok=True)
    input_fasta = output_dir / "input.fna"
    input_fasta.write_text(">demo\nACGT\n")
    lock_sha256 = canonical_sha256({"fixture": identity_seed, "kind": "lock"})
    runtime_sha256 = canonical_sha256(
        {"fixture": identity_seed, "kind": "runtime"}
    )
    environment_payload = {
        "lock_sha256": lock_sha256,
        "runtime_sha256": runtime_sha256,
        "requested_device": "cpu",
        "effective_device": "cpu",
    }
    return {
        "genome_id": "demo",
        "input_path": str(input_fasta.resolve()),
        "output_dir": str(output_dir.resolve()),
        "input": asdict(build_input_identity(input_fasta)),
        "config": {
            "sha256": canonical_sha256(
                {"fixture": identity_seed, "kind": "config"}
            ),
            "ablation_id": "A0",
            "ablation_contract_sha256": ABLATION_CONTRACT_SHA256,
        },
        "code": {
            "version": "test",
            "source_sha256": canonical_sha256(
                {"fixture": identity_seed, "kind": "source"}
            ),
        },
        "environment": {
            **environment_payload,
            "sha256": canonical_sha256(environment_payload),
        },
        "coordinate_schema_version": COORDINATE_SCHEMA_VERSION,
        "coordinate_convention": COORDINATE_CONVENTION,
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "summary_schema_version": 3,
        "requested_masking": _masking_request_identity(MaskingConfig()),
        "resources": [],
    }


def _start_schema3_run(output_dir: Path, *, identity_seed: str = "demo") -> str:
    """Publish a deterministic schema-v3 running state before fixture outputs."""
    identities = _closed_schema3_identity(
        output_dir,
        identity_seed=identity_seed,
    )
    run_fingerprint = compute_run_fingerprint(identities)
    state = publish_run_started(
        output_dir,
        run_fingerprint=run_fingerprint,
        identities=identities,
    )
    return state.run_fingerprint


def _publish_schema3_success(output_dir: Path, run_fingerprint: str) -> None:
    """Publish a complete four-phase, zero-result schema-v3 success fixture."""
    ablation_content = AblationEvents(
        ablation_id=AblationID.A0,
        counters=AblationCounters(),
    ).to_bytes()
    requested_masking = _masking_request_identity(MaskingConfig())
    input_fasta = output_dir / "input.fna"
    masking_result = mask_genome_pipeline(
        input_fasta,
        output_dir / "phase0" / "masking",
        config=MaskingConfig(),
    )
    masking_status = masking_result.status_path
    status_payload = json.loads(masking_status.read_text())
    proteome = output_dir / "phase0" / "proteome.fasta"
    genes = output_dir / "phase0" / "genes.gff"
    proteome.write_text(">p1\nM\n")
    genes.write_text("demo\ttest\tCDS\t1\t4\t.\t+\t0\tID=p1\n")
    phase0_ablation = output_dir / "phase0" / "ablation_events.json"
    phase0_ablation.write_bytes(ablation_content)
    phase0_artifacts = tuple(
        build_artifact_identity(path, root=output_dir, schema=schema)
        for path, schema in (
            (phase0_ablation, "virosync.ablation_events/v1"),
            (masking_status, "masking-status-v1"),
            (proteome, "protein-fasta-v1"),
            (genes, "gene-gff-v1"),
        )
    )
    status_identity = next(
        artifact
        for artifact in phase0_artifacts
        if artifact.relative_path == "phase0/masking/masking_status.json"
    )
    publish_phase_completion(
        output_dir,
        phase=0,
        run_fingerprint=run_fingerprint,
        dependency_sha256=run_fingerprint,
        artifacts=phase0_artifacts,
        outcome="complete",
        requested_masking=requested_masking,
        actual_masking={**status_payload, "status_sha256": status_identity.sha256},
    )

    phase1_state = output_dir / "phase1" / "resume_state.json"
    write_phase1_state(
        phase1_state,
        validated_markers=[],
        merged_seeds=[
            MergedSeed(
                scaffold="demo",
                start=0,
                end=4,
                seed_id="seed-0-demo-0",
                sources=["test"],
            )
        ],
        host_signature_model=HostSignatureModel(),
        host_signatures=set(),
        host_deviation_summary=None,
    )
    phase1_ablation = output_dir / "phase1" / "ablation_events.json"
    phase1_ablation.write_bytes(ablation_content)
    phase1_artifacts = (
        build_artifact_identity(
            phase1_ablation,
            root=output_dir,
            schema="virosync.ablation_events/v1",
        ),
        build_artifact_identity(
            phase1_state,
            root=output_dir,
            schema=PHASE1_STATE_SCHEMA,
        ),
    )
    publish_phase_completion(
        output_dir,
        phase=1,
        run_fingerprint=run_fingerprint,
        dependency_sha256=marker_sha256(output_dir, 0),
        artifacts=phase1_artifacts,
        outcome="complete",
    )

    boundary = RefinedBoundary(
        scaffold="demo",
        start=0,
        end=4,
        seed_id="seed-0-demo-0",
        seed_sources=["test"],
        confidence=0.9,
    )
    refined_state = output_dir / "phase2" / "refined_state.json"
    phase2_resume_state = output_dir / "phase2" / "resume_state.json"
    refined_bed = output_dir / "phase2" / "refined_boundaries.bed"
    write_phase2_state(refined_state, [boundary])
    write_phase2_resume_state(
        phase2_resume_state,
        refined_boundaries=[boundary],
        boundary_taxonomy_map={},
        boundary_control_stats=None,
        boundary_diamond_query=None,
    )
    refined_bed.write_text("demo\t0\t4\tEVE_demo_0-4\t900\t.\n")
    phase2_ablation = output_dir / "phase2" / "ablation_events.json"
    phase2_ablation.write_bytes(ablation_content)
    phase2_artifacts = tuple(
        build_artifact_identity(path, root=output_dir, schema=schema)
        for path, schema in (
            (phase2_ablation, "virosync.ablation_events/v1"),
            (refined_state, PHASE2_STATE_SCHEMA),
            (phase2_resume_state, PHASE2_RESUME_STATE_SCHEMA),
            (refined_bed, "refined-boundaries-bed-v1"),
        )
    )
    publish_phase_completion(
        output_dir,
        phase=2,
        run_fingerprint=run_fingerprint,
        dependency_sha256=marker_sha256(output_dir, 1),
        artifacts=phase2_artifacts,
        outcome="complete",
    )

    phase3_dir = output_dir / "phase3_synthesis"
    phase3_dir.mkdir(parents=True, exist_ok=True)
    root_canonical = output_dir / "virosync_predictions.tsv"
    phase3_canonical = phase3_dir / "virosync_predictions.tsv"
    canonical = root_canonical if root_canonical.is_file() else phase3_canonical
    canonical.write_text(_ZERO_PREDICTION_HEADER)
    detailed = output_dir / "virosync_predictions_detailed.tsv"
    detailed.write_text(_ZERO_PREDICTION_HEADER)
    export_dir = canonical.parent
    bed = export_dir / "virosync_predictions.bed"
    gff = export_dir / "virosync_predictions.gff3"
    summary = export_dir / "virosync_summary.json"
    bed.write_text("")
    gff.write_text("##gff-version 3\n")
    summary.write_text(
        json.dumps(
            {
                "virosync_version": "test",
                "coordinate_schema_version": COORDINATE_SCHEMA_VERSION,
                "coordinate_convention": COORDINATE_CONVENTION,
                "output_schema_version": OUTPUT_SCHEMA_VERSION,
                "statistics": {
                    "canonical_predictions": 0,
                    "total_candidates": 0,
                    "total_accepted_length_bp": 0,
                    "high_confidence": 0,
                    "medium_confidence": 0,
                    "low_confidence": 0,
                    "promoted_low_confidence": 0,
                },
                "per_scaffold": {},
            },
            sort_keys=True,
        )
        + "\n"
    )
    run_log = output_dir / "run.log"
    if not run_log.is_file():
        run_log.write_text(
            "# ViroSync Run Log: demo\n\n## Results Summary\nCanonical GEVEs: 0\n"
        )
    invariant = output_dir / "virosync_tsv_invariant_report.tsv"
    invariant.write_text(
        "status\trows_checked\tissue_count\terror_count\twarning_count\n"
        "PASS\t0\t0\t0\t0\n\n"
        "eve_id\tcheck\tseverity\tmessage\n"
    )
    notebook = output_dir / "notebooks" / "jupyter" / "eve_analysis.ipynb"
    notebook.parent.mkdir(parents=True, exist_ok=True)
    notebook.write_text(
        '{"cells":[],"metadata":{},"nbformat":4,"nbformat_minor":5}\n'
    )
    completion = _write_completion_manifest(
        output_dir,
        genome_id="demo",
        status="success",
        fingerprint=run_fingerprint,
        output_files={
            "predictions": str(canonical),
            "predictions_detailed": str(detailed),
        },
    )
    phase3_ablation = output_dir / "phase3" / "ablation_events.json"
    phase3_ablation.parent.mkdir(parents=True, exist_ok=True)
    phase3_ablation.write_bytes(ablation_content)
    ablation_events = output_dir / "ablation_events.json"
    ablation_events.write_bytes(ablation_content)
    final_artifacts = tuple(
        build_artifact_identity(path, root=output_dir, schema=schema)
        for path, schema in (
            (ablation_events, "virosync.ablation_events/v1"),
            (canonical, "canonical-predictions-v3"),
            (detailed, "detailed-predictions-v3"),
            (bed, "canonical-predictions-bed-v1"),
            (gff, "canonical-predictions-gff3-v1"),
            (summary, "virosync-summary-v3"),
            (invariant, "tsv-invariant-report-v1"),
            (run_log, "run-log-v1"),
            (completion, "completion-manifest-v2"),
            (notebook, "eve-analysis-notebook-v1"),
        )
    )
    publish_phase_completion(
        output_dir,
        phase=3,
        run_fingerprint=run_fingerprint,
        dependency_sha256=marker_sha256(output_dir, 2),
        artifacts=(
            build_artifact_identity(
                phase3_ablation,
                root=output_dir,
                schema="virosync.ablation_events/v1",
            ),
            *final_artifacts,
        ),
        outcome="complete",
    )
    publish_run_success(
        output_dir,
        run_fingerprint=run_fingerprint,
        artifacts=final_artifacts,
        result={
            "terminal_phase": None,
            "canonical_rows": 0,
            "detailed_rows": 0,
            "accepted_bp": 0,
            "class_counts": {
                eve_class: 0 for eve_class in EFFECTIVE_EVE_CLASS_COUNT_KEYS
            },
            "tier_counts": {"HIGH": 0, "MEDIUM": 0, "LOW": 0},
            "benchmark_eligible": True,
        },
    )


def test_completed_run_artifacts_require_root_outputs(tmp_path: Path) -> None:
    run_fingerprint = _start_schema3_run(tmp_path)
    phase3_dir = tmp_path / "phase3_synthesis"
    phase3_dir.mkdir()
    (phase3_dir / "virosync_predictions.tsv").write_text("eve_id\n")

    assert _completed_run_artifacts(tmp_path) is None

    (tmp_path / "virosync_predictions_detailed.tsv").write_text("eve_id\n")
    assert _completed_run_artifacts(tmp_path) is None

    (tmp_path / "run.log").write_text(
        "# ViroSync Run Log: demo\n\n## Results Summary\nCanonical GEVEs: 0\n"
    )
    assert _completed_run_artifacts(tmp_path) is None

    _publish_schema3_success(tmp_path, run_fingerprint)
    artifacts = _completed_run_artifacts(tmp_path)

    assert artifacts is not None
    assert artifacts["phase3_predictions"] == phase3_dir / "virosync_predictions.tsv"
    assert artifacts["predictions_detailed"] == tmp_path / "virosync_predictions_detailed.tsv"
    assert artifacts["run_log"] == tmp_path / "run.log"
    assert artifacts["completion_manifest"] == tmp_path / "virosync_run_complete.json"
    assert artifacts["run_state"] == tmp_path / "virosync_run_state.json"


def test_completed_run_artifacts_reject_mutated_run_log(tmp_path: Path) -> None:
    run_fingerprint = _start_schema3_run(tmp_path)
    phase3_dir = tmp_path / "phase3_synthesis"
    phase3_dir.mkdir()
    (phase3_dir / "virosync_predictions.tsv").write_text("eve_id\n")
    (tmp_path / "virosync_predictions_detailed.tsv").write_text("eve_id\n")
    (tmp_path / "run.log").write_text(
        "# ViroSync Run Log: demo\n\n## Results Summary\nCanonical GEVEs: 0\n"
    )
    _publish_schema3_success(tmp_path, run_fingerprint)
    (tmp_path / "run.log").write_text("complete\n")

    assert _completed_run_artifacts(tmp_path) is None


def test_completed_run_artifacts_reject_missing_completion_manifest(
    tmp_path: Path,
) -> None:
    phase3_dir = tmp_path / "phase3_synthesis"
    phase3_dir.mkdir()
    (phase3_dir / "virosync_predictions.tsv").write_text("eve_id\n")
    (tmp_path / "virosync_predictions_detailed.tsv").write_text("eve_id\n")
    (tmp_path / "run.log").write_text(
        "# ViroSync Run Log: demo\n\n## Results Summary\nCanonical GEVEs: 0\n"
    )

    assert _completed_run_artifacts(tmp_path) is None


def test_completed_run_artifacts_reject_corrupt_completion_manifest(
    tmp_path: Path,
) -> None:
    phase3_dir = tmp_path / "phase3_synthesis"
    phase3_dir.mkdir()
    (phase3_dir / "virosync_predictions.tsv").write_text("eve_id\n")
    (tmp_path / "virosync_predictions_detailed.tsv").write_text("eve_id\n")
    (tmp_path / "run.log").write_text(
        "# ViroSync Run Log: demo\n\n## Results Summary\nCanonical GEVEs: 0\n"
    )
    (tmp_path / "virosync_run_complete.json").write_text("{not-json")

    assert _completed_run_artifacts(tmp_path) is None


def test_completed_run_artifacts_accept_zero_result_root_outputs(
    tmp_path: Path,
) -> None:
    run_fingerprint = _start_schema3_run(tmp_path)
    (tmp_path / "virosync_predictions.tsv").write_text("eve_id\n")
    (tmp_path / "virosync_predictions_detailed.tsv").write_text("eve_id\n")
    (tmp_path / "run.log").write_text(
        "# ViroSync Run Log: demo\n\n## Results Summary\nGEVEs detected: 0\n"
    )
    _publish_schema3_success(tmp_path, run_fingerprint)

    artifacts = _completed_run_artifacts(tmp_path)

    assert artifacts is not None
    assert artifacts["phase3_predictions"] == tmp_path / "virosync_predictions.tsv"
    assert artifacts["completion_manifest"] == tmp_path / "virosync_run_complete.json"


def test_write_completion_manifest_serializes_nested_paths(tmp_path: Path) -> None:
    _seed_masking_status(tmp_path)
    manifest = _write_completion_manifest(
        output_dir=tmp_path,
        genome_id="demo",
        status="success",
        output_files={
            "plain_path": tmp_path / "out.tsv",
            "nested": {"path": tmp_path / "nested.tsv"},
        },
    )

    text = manifest.read_text()

    assert str(tmp_path / "out.tsv") in text
    assert str(tmp_path / "nested.tsv") in text


def test_completion_manifest_records_coordinate_contract(tmp_path: Path) -> None:
    _seed_masking_status(tmp_path)
    manifest = _write_completion_manifest(
        output_dir=tmp_path,
        genome_id="demo",
        status="success",
    )

    payload = json.loads(manifest.read_text())

    assert payload["coordinate_schema_version"] == COORDINATE_SCHEMA_VERSION
    assert payload["output_schema_version"] == OUTPUT_SCHEMA_VERSION
    assert payload["coordinate_convention"] == COORDINATE_CONVENTION


def test_require_phase2b_gene_taxonomy_db_allows_no_seed_runs() -> None:
    assert _require_phase2b_gene_taxonomy_db(None, has_seeds=False) is None


def test_require_phase2b_gene_taxonomy_db_rejects_marker_db_fallback() -> None:
    with pytest.raises(ValueError, match="gene_taxonomy_faa_db"):
        _require_phase2b_gene_taxonomy_db(None, has_seeds=True)


def test_require_phase2b_gene_taxonomy_db_rejects_missing_path(tmp_path: Path) -> None:
    missing = tmp_path / "missing.dmnd"

    with pytest.raises(FileNotFoundError, match="gene taxonomy database not found"):
        _require_phase2b_gene_taxonomy_db(missing, has_seeds=True)


def test_require_phase2b_gene_taxonomy_db_accepts_existing_path(tmp_path: Path) -> None:
    db = tmp_path / "combined_proteome.dmnd"
    db.write_text("diamond")

    assert _require_phase2b_gene_taxonomy_db(db, has_seeds=True) == db


def test_phase1_resume_rebuild_preserves_zero_based_seed_ids() -> None:
    seeds = _build_merged_seeds_from_regions(
        [
            {"region_id": "r1", "scaffold": "a", "start": 10, "end": 20},
            {"region_id": "r2", "scaffold": "b", "start": 30, "end": 40},
        ],
        [],
    )

    assert [seed.seed_id for seed in seeds] == ["seed_0_a_10", "seed_1_b_30"]


def test_summarize_predictions_tsv_counts_canonical_rows_for_rollups(tmp_path: Path) -> None:
    predictions_tsv = tmp_path / "virosync_predictions.tsv"
    predictions_tsv.write_text(
        "\t".join(
            [
                "eve_id",
                "confidence_tier",
                "length",
                "gene_taxonomy_total",
                "hallmark_total",
                "classification",
            ]
        )
        + "\n"
        + "\n".join(
            [
                "eve1\tHIGH\t1000\t10\t2\tNCLDV",
                "eve2\tMEDIUM\t500\t5\t1\tPLV",
                "eve3\tLOW\t250\t3\t0\tVP",
            ]
        )
        + "\n"
    )

    stats = _summarize_predictions_tsv(predictions_tsv)

    assert stats["predictions"] == 3
    assert stats["accepted"] == 3
    assert stats["high_tier"] == 1
    assert stats["medium_tier"] == 1
    assert stats["low_tier"] == 1
    assert stats["accepted_bp"] == 1750
    assert stats["total_genes"] == 18
    assert stats["total_hallmarks"] == 3
    assert stats["ncldv_count"] == 1
    # Legacy "PLV" and "VP" rows both roll up under the unified
    # Preplasmiviricota class, so neither legacy counter fires.
    assert stats["ppv_count"] == 2
    assert stats["plv_count"] == 0
    assert stats["vp_count"] == 0


def test_summarize_predictions_tsv_counts_detailed_candidates_without_accepting(
    tmp_path: Path,
) -> None:
    detailed_tsv = tmp_path / "virosync_predictions_detailed.tsv"
    detailed_tsv.write_text(
        "eve_id\tconfidence_tier\tlength\tgene_taxonomy_total\thallmark_total\tclassification\n"
        "eve1\tHIGH\t1000\t10\t2\tNCLDV\n"
        "eve2\tLOW\t250\t3\t0\tVP\n"
    )

    stats = _summarize_predictions_tsv(detailed_tsv, canonical=False)

    assert stats["predictions"] == 2
    assert stats["accepted"] == 0
    assert stats["high_tier"] == 1
    assert stats["low_tier"] == 1


def test_summarize_predictions_prefers_persisted_effective_class(tmp_path: Path) -> None:
    predictions_tsv = tmp_path / "virosync_predictions.tsv"
    predictions_tsv.write_text(
        "eve_id\tconfidence_tier\tlength\tclassification\t"
        "region_classification\teffective_eve_class\n"
        "ppv\tHIGH\t3001\tNCLDV\tNCLDV\tPPV\n"
        "future\tMEDIUM\t3001\tNCLDV\tNCLDV\tfuture-lineage\n"
    )

    stats = _summarize_predictions_tsv(predictions_tsv)

    assert stats["accepted"] == 2
    assert stats["ppv_count"] == 1
    assert stats["unknown_count"] == 1
    assert stats["ncldv_count"] == 0
    assert effective_eve_class_count_total(stats) == stats["accepted"]


def test_summarize_legacy_predictions_uses_tier_aware_raw_fallback(
    tmp_path: Path,
) -> None:
    predictions_tsv = tmp_path / "legacy_predictions.tsv"
    predictions_tsv.write_text(
        "eve_id\tconfidence_tier\tlength\tclassification\tregion_classification\n"
        "high\tHIGH\t6001\tNCLDV\tPPV\n"
        "low\tLOW\t6001\tNCLDV\tPPV\n"
        "blank\tMEDIUM\t6001\t\t\n"
    )

    stats = _summarize_predictions_tsv(predictions_tsv)

    assert stats["accepted"] == 3
    assert stats["ppv_count"] == 1
    assert stats["ncldv_count"] == 1
    assert stats["unknown_count"] == 1
    assert effective_eve_class_count_total(stats) == stats["accepted"]
    assert set(EFFECTIVE_EVE_CLASS_COUNT_KEYS.values()).issubset(stats)


def test_fresh_and_resume_batch_summaries_use_identical_persisted_metrics(
    tmp_path: Path,
) -> None:
    canonical = tmp_path / "virosync_predictions.tsv"
    detailed = tmp_path / "virosync_predictions_detailed.tsv"
    canonical.write_text(
        "eve_id\tconfidence_tier\tlength\tgene_taxonomy_total\thallmark_total\t"
        "effective_eve_class\n"
        "ppv\tLOW\t3001\t11\t5\tPPV\n"
    )
    detailed.write_text(
        "eve_id\tconfidence_tier\tlength\tgene_taxonomy_total\thallmark_total\t"
        "effective_eve_class\n"
        "ppv\tLOW\t3001\t11\t5\tPPV\n"
    )
    in_memory_phase3 = {"total_genes": 97, "total_hallmarks": 2}

    fresh = _summarize_prediction_outputs(
        canonical,
        detailed,
        expected_accepted=1,
        expected_candidates=1,
    )
    resumed = _summarize_prediction_outputs(canonical, detailed)

    assert fresh == resumed
    assert fresh["total_genes"] == 11 != in_memory_phase3["total_genes"]
    assert fresh["total_hallmarks"] == 5 != in_memory_phase3["total_hallmarks"]
    assert fresh["accepted_bp"] == 3001
    assert fresh["ppv_count"] == 1

    fresh_result = {
        "genome_id": "demo",
        "success": True,
        "elapsed_sec": 0.0,
        **fresh,
    }
    resumed_result = {
        "genome_id": "demo",
        "success": True,
        "elapsed_sec": 0.0,
        **resumed,
    }
    fresh_dir = tmp_path / "fresh"
    resume_dir = tmp_path / "resume"
    fresh_dir.mkdir()
    resume_dir.mkdir()

    fresh_batch = python_runner._write_batch_summary(fresh_dir, [fresh_result])
    resumed_batch = python_runner._write_batch_summary(resume_dir, [resumed_result])

    assert fresh_batch.read_bytes() == resumed_batch.read_bytes()


def test_fresh_persisted_summary_fails_closed_on_in_memory_count_drift(
    tmp_path: Path,
) -> None:
    canonical = tmp_path / "virosync_predictions.tsv"
    detailed = tmp_path / "virosync_predictions_detailed.tsv"
    canonical.write_text(
        "eve_id\tconfidence_tier\tlength\teffective_eve_class\n"
        "ppv\tHIGH\t3001\tPPV\n"
    )
    detailed.write_text(
        "eve_id\tconfidence_tier\tlength\teffective_eve_class\n"
        "ppv\tHIGH\t3001\tPPV\n"
    )

    with pytest.raises(ValueError, match="persisted accepted count disagrees"):
        _summarize_prediction_outputs(
            canonical,
            detailed,
            expected_accepted=2,
            expected_candidates=1,
        )
    with pytest.raises(ValueError, match="persisted candidate count disagrees"):
        _summarize_prediction_outputs(
            canonical,
            detailed,
            expected_accepted=1,
            expected_candidates=2,
        )
