from __future__ import annotations

import logging
from pathlib import Path

import pytest

from virosync.ablation import AblationID, InterventionCounts
from virosync.orchestration._flows.single_genome import phase2
from virosync.orchestration._flows.single_genome.phase2_resume_state import (
    load_phase2_resume_state,
)
from virosync.orchestration._flows.single_genome.phase_state import (
    load_phase2_state,
)
from virosync.pipeline.phase1.seed_merger import MergedSeed


def _run_a3(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict:
    masked = tmp_path / "masked.fna"
    masked.write_text(">scaffold\n" + "ACGT" * 100 + "\n")
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">p1\nM\n")
    seed = MergedSeed(
        scaffold="scaffold",
        start=20,
        end=220,
        seed_id="seed-1",
        sources=["hhg", "marker_validation"],
        confidence="high",
    )

    def forbidden(*_args, **_kwargs):
        raise AssertionError("A3 entered Phase-2 refinement biology")

    monkeypatch.setattr(phase2, "_require_phase2b_gene_taxonomy_db", forbidden)
    monkeypatch.setattr(phase2, "call_task", forbidden)
    proteome_index = {"scaffold": ["p1"]}
    monkeypatch.setattr(
        phase2,
        "build_proteome_index",
        lambda _path: proteome_index,
    )

    return phase2._run_phase2_subflow(
        masked_path=masked,
        proteome_path=proteome,
        merged_seeds=[seed],
        validated_markers=[],
        host_signature_model=None,
        output_dir=tmp_path,
        genome_id="genome",
        resume=False,
        refined_bed=tmp_path / "phase2" / "refined_boundaries.bed",
        gene_taxonomy_faa_db=None,
        marker_db=None,
        taxonomy_labels_file=None,
        host_prefixes=["EUK__"],
        host_label="EUK",
        high_pident_host_threshold=70.0,
        boundary_host_trim_enabled=True,
        boundary_host_trim_window_bp=1000,
        boundary_host_trim_step_bp=500,
        boundary_host_trim_max_host_fraction=0.5,
        boundary_host_trim_min_viral_fraction=0.1,
        boundary_host_trim_score_threshold=0.5,
        boundary_host_trim_buffer_kb=1,
        boundary_host_trim_min_overlap_score=0.2,
        boundary_host_signature_min_token_len=3,
        taxonomy_weight_mode="rank",
        boundary_taxonomy_ml_enabled=False,
        boundary_taxonomy_ml_model="logreg",
        boundary_taxonomy_ml_threshold=0.5,
        boundary_taxonomy_ml_neighbor_window=1,
        boundary_diamond_flank_genes=10,
        boundary_diamond_control_sample_size=10,
        boundary_diamond_control_min_distance=5,
        boundary_diamond_top_k=5,
        boundary_diamond_chunk_size=100,
        boundary_diamond_random_seed=42,
        threads=1,
        gene_taxonomy_threads=1,
        extended_output=False,
        search_backend="diamond",
        genome_start_time=0.0,
        logger=logging.getLogger(__name__),
        ablation_id=AblationID.A3,
    )


def test_a3_forwards_exact_phase1_intervals_without_phase2_biology(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    result = _run_a3(tmp_path, monkeypatch)

    assert result["phase_outcome"] == "passthrough"
    assert result["ablation_counts"] == InterventionCounts(1, 1, 0)
    assert result["boundary_taxonomy_map"] == {}
    assert result["boundary_control_stats"] is None
    assert result["boundary_diamond_query"] is None
    assert result["proteome_index"] == {"scaffold": ["p1"]}
    boundary = result["refined_boundaries"][0]
    assert (boundary.scaffold, boundary.start, boundary.end) == (
        "scaffold",
        20,
        220,
    )
    assert (boundary.original_start, boundary.original_end) == (20, 220)

    report_state = load_phase2_state(tmp_path / "phase2" / "refined_state.json")
    resume_state = load_phase2_resume_state(
        tmp_path / "phase2" / "resume_state.json"
    )
    assert report_state == result["refined_boundaries"]
    assert resume_state.refined_boundaries == result["refined_boundaries"]
    assert (tmp_path / "phase2" / "refined_boundaries.bed").read_text() == (
        "scaffold\t20\t220\tEVE_scaffold_20-220\t0\t.\n"
    )
