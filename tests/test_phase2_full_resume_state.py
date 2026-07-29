from __future__ import annotations

import copy
import json
import logging
from pathlib import Path

import numpy as np
import pytest

from virosync.orchestration._flows.single_genome import phase2
from virosync.orchestration._flows.single_genome.phase2_resume_state import (
    PHASE2_RESUME_STATE_ARTIFACT_TYPE,
    PHASE2_RESUME_STATE_FILENAME,
    PHASE2_RESUME_STATE_SCHEMA_VERSION,
    Phase2ResumeState,
    Phase2ResumeStateError,
    load_phase2_resume_state,
    phase2_resume_state_from_document,
    phase2_resume_state_to_document,
    write_phase2_resume_state,
)
from virosync.orchestration._flows.single_genome.phase_state import (
    PHASE2_STATE_FILENAME,
    phase2_state_to_document,
    write_phase2_state,
)
from virosync.pipeline.phase2.boundary_diamond import (
    ControlStats,
    GeneTaxonomy,
    GenomeDiamondQuery,
    SeedGeneMapping,
)
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.taxonomy_utils import TaxonomyFingerprint


def _boundary() -> RefinedBoundary:
    return RefinedBoundary(
        scaffold="scaffold/alpha",
        start=101,
        end=999,
        seed_id="seed-b",
        original_start=90,
        original_end=1010,
        candidate_start=80,
        candidate_end=1020,
        seed_sources=["hhg", "novelty"],
        seed_confidence="high",
        seed_hhg_score=0.9123456789012345,
        confidence=0.8765432109876543,
        posterior_probability=0.9345678901234567,
        state_sequence=[4, 5],
        state_posteriors=np.array(
            [[0.01, 0.02, 0.03, 0.04, 0.40, 0.50]],
            dtype=np.float64,
        ),
        hallmark_genes=["MCP", "A32"],
    )


def _taxonomy(porf_id: str, pident: float) -> GeneTaxonomy:
    return GeneTaxonomy(
        porf_id=porf_id,
        scaffold="scaffold/alpha",
        start=100,
        end=400,
        top1_target=f"EUK__target|{porf_id}",
        top1_prefix="EUK__",
        top1_pident=pident,
        top1_evalue=1.2345678901234568e-42,
        top10_prefixes=["EUK__", "NCLDV__"],
        top10_targets=[f"EUK__target|{porf_id}", "NCLDV__target|mcp"],
        top10_bits=[333.1234567890123, 222.9876543210987],
        top10_pidents=[pident, 42.42424242424242],
        top10_evalues=[1.2345678901234568e-42, 9.876543210987655e-17],
        taxonomy_fingerprint=TaxonomyFingerprint(
            weighted_tokens={
                "Viridiplantae": 0.12345678901234566,
                "Nucleocytoviricota": 0.8765432109876543,
            },
            raw_tokens={"Viridiplantae": 7, "Nucleocytoviricota": 3},
        ),
        has_ncldv_mirus=True,
        has_vp_plv=False,
        has_viral=True,
        has_hit=True,
    )


def _seed_mapping(seed_id: str, offset: int) -> SeedGeneMapping:
    return SeedGeneMapping(
        seed_id=seed_id,
        scaffold="scaffold/alpha",
        seed_start=100 + offset,
        seed_end=900 + offset,
        eve_porf_ids=[f"{seed_id}-eve-2", f"{seed_id}-eve-1"],
        upstream_porf_ids=[f"{seed_id}-up-1", f"{seed_id}-up-2"],
        downstream_porf_ids=[f"{seed_id}-down-1", f"{seed_id}-down-2"],
        flank_start_idx=11 + offset,
        flank_end_idx=29 + offset,
        flank_start_bp=50 + offset,
        flank_end_bp=1050 + offset,
        flank_genes_config=17,
    )


def _state() -> Phase2ResumeState:
    seed_b = _seed_mapping("seed-b", 0)
    seed_a = _seed_mapping("seed-a", 1)
    return Phase2ResumeState(
        refined_boundaries=[_boundary()],
        boundary_taxonomy_map={
            "porf-over": _taxonomy("porf-over", 70.001),
            "porf-under": _taxonomy("porf-under", 69.999),
        },
        boundary_control_stats=ControlStats(
            n_genes=23,
            n_no_hits=4,
            no_hit_frequency=0.17391304347826086,
            host_frequency=0.7391304347826086,
            mean_pident=69.99999999999999,
            dominant_organism="Viridiplantae species alpha",
            host_prefix="EUK__",
        ),
        boundary_diamond_query=GenomeDiamondQuery(
            eve_porf_ids={
                "seed-b": list(seed_b.eve_porf_ids),
                "seed-a": list(seed_a.eve_porf_ids),
            },
            boundary_porf_ids={
                "seed-b": [
                    *seed_b.upstream_porf_ids,
                    *seed_b.downstream_porf_ids,
                ],
                "seed-a": [
                    *seed_a.upstream_porf_ids,
                    *seed_a.downstream_porf_ids,
                ],
            },
            control_porf_ids=["control-9", "control-1", "control-5"],
            all_porf_ids=["seed-b-eve-2", "control-9", "seed-a-eve-1"],
            seed_gene_mappings={"seed-b": seed_b, "seed-a": seed_a},
        ),
    )


def _assert_phase3_inputs_equal(
    expected: Phase2ResumeState,
    observed: Phase2ResumeState,
) -> None:
    assert phase2_state_to_document(
        observed.refined_boundaries
    ) == phase2_state_to_document(expected.refined_boundaries)
    assert observed.boundary_taxonomy_map == expected.boundary_taxonomy_map
    assert observed.boundary_control_stats == expected.boundary_control_stats
    assert observed.boundary_diamond_query == expected.boundary_diamond_query
    assert list(observed.boundary_taxonomy_map) == list(
        expected.boundary_taxonomy_map
    )
    assert list(observed.boundary_diamond_query.eve_porf_ids) == [
        "seed-b",
        "seed-a",
    ]
    assert list(observed.boundary_diamond_query.seed_gene_mappings) == [
        "seed-b",
        "seed-a",
    ]


def test_phase2_resume_state_round_trip_preserves_exact_phase3_inputs(
    tmp_path: Path,
) -> None:
    original = _state()
    state_path = tmp_path / "phase2" / PHASE2_RESUME_STATE_FILENAME

    write_phase2_resume_state(
        state_path,
        refined_boundaries=original.refined_boundaries,
        boundary_taxonomy_map=original.boundary_taxonomy_map,
        boundary_control_stats=original.boundary_control_stats,
        boundary_diamond_query=original.boundary_diamond_query,
    )
    loaded = load_phase2_resume_state(state_path)

    _assert_phase3_inputs_equal(original, loaded)
    assert loaded.boundary_taxonomy_map["porf-under"].top1_pident == 69.999
    assert loaded.boundary_taxonomy_map["porf-over"].top1_pident == 70.001
    expected_fingerprint = original.boundary_taxonomy_map[
        "porf-under"
    ].taxonomy_fingerprint
    assert (
        loaded.boundary_taxonomy_map["porf-under"].taxonomy_fingerprint
        == expected_fingerprint
    )
    payload = json.loads(state_path.read_text())
    assert payload["artifact_type"] == PHASE2_RESUME_STATE_ARTIFACT_TYPE
    assert payload["schema_version"] == PHASE2_RESUME_STATE_SCHEMA_VERSION


def test_phase2_resume_state_round_trip_preserves_optional_none_values(
) -> None:
    state = Phase2ResumeState(
        refined_boundaries=[],
        boundary_taxonomy_map={},
        boundary_control_stats=None,
        boundary_diamond_query=None,
    )

    loaded = phase2_resume_state_from_document(
        phase2_resume_state_to_document(state)
    )

    assert loaded == state


def test_phase2_resume_state_rejects_schema_drift_and_duplicate_mapping_keys(
) -> None:
    document = phase2_resume_state_to_document(_state())

    unknown_schema = copy.deepcopy(document)
    unknown_schema["schema_version"] = 2
    with pytest.raises(
        Phase2ResumeStateError,
        match="unsupported.*schema_version",
    ):
        phase2_resume_state_from_document(unknown_schema)

    extra_seed_field = copy.deepcopy(document)
    seed_document = extra_seed_field["boundary_diamond_query"][
        "seed_gene_mappings"
    ][0]["value"]
    seed_document["runtime_class"] = "arbitrary.Type"
    with pytest.raises(Phase2ResumeStateError, match="extra=.*runtime_class"):
        phase2_resume_state_from_document(extra_seed_field)

    duplicate_taxonomy = copy.deepcopy(document)
    duplicate_taxonomy["boundary_taxonomy_map"].append(
        copy.deepcopy(duplicate_taxonomy["boundary_taxonomy_map"][0])
    )
    with pytest.raises(
        Phase2ResumeStateError,
        match="duplicate key 'porf-over'",
    ):
        phase2_resume_state_from_document(duplicate_taxonomy)


def test_phase2_resume_state_rejects_nonfinite_nested_float() -> None:
    document = phase2_resume_state_to_document(_state())
    document["boundary_taxonomy_map"][0]["value"]["top1_pident"] = float("nan")

    with pytest.raises(
        Phase2ResumeStateError,
        match="top1_pident must be finite",
    ):
        phase2_resume_state_from_document(document)


def test_phase2_resume_state_failed_validation_preserves_existing_file(
    tmp_path: Path,
) -> None:
    state_path = tmp_path / PHASE2_RESUME_STATE_FILENAME
    state_path.write_text("previous valid state\n")
    state = _state()
    state.boundary_taxonomy_map["porf-under"].top1_pident = float("inf")

    with pytest.raises(
        Phase2ResumeStateError,
        match="top1_pident must be finite",
    ):
        write_phase2_resume_state(
            state_path,
            refined_boundaries=state.refined_boundaries,
            boundary_taxonomy_map=state.boundary_taxonomy_map,
            boundary_control_stats=state.boundary_control_stats,
            boundary_diamond_query=state.boundary_diamond_query,
        )

    assert state_path.read_text() == "previous valid state\n"


def test_authenticated_phase2_resume_supplies_exact_checkpoint_inputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original = _state()
    phase2_dir = tmp_path / "phase2"
    phase2_dir.mkdir()
    write_phase2_state(
        phase2_dir / PHASE2_STATE_FILENAME,
        original.refined_boundaries,
    )
    write_phase2_resume_state(
        phase2_dir / PHASE2_RESUME_STATE_FILENAME,
        refined_boundaries=original.refined_boundaries,
        boundary_taxonomy_map=original.boundary_taxonomy_map,
        boundary_control_stats=original.boundary_control_stats,
        boundary_diamond_query=original.boundary_diamond_query,
    )
    proteome_index = {
        "scaffold/alpha": ["derived-from-authenticated-proteome"]
    }
    monkeypatch.setattr(
        phase2,
        "build_proteome_index",
        lambda _path: proteome_index,
    )

    result = phase2._run_phase2_subflow(
        masked_path=tmp_path / "masked.fna",
        proteome_path=tmp_path / "proteome.faa",
        merged_seeds=[object()],
        validated_markers=[],
        host_signature_model=None,
        output_dir=tmp_path,
        genome_id="genome",
        resume=True,
        refined_bed=phase2_dir / "refined_boundaries.bed",
        gene_taxonomy_faa_db=None,
        marker_db=None,
        taxonomy_labels_file=None,
        host_prefixes=[],
        host_label="EUK",
        high_pident_host_threshold=70.0,
        boundary_host_trim_enabled=False,
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
        boundary_taxonomy_ml_model="logistic",
        boundary_taxonomy_ml_threshold=0.5,
        boundary_taxonomy_ml_neighbor_window=1,
        boundary_diamond_flank_genes=17,
        boundary_diamond_control_sample_size=23,
        boundary_diamond_control_min_distance=11,
        boundary_diamond_top_k=10,
        boundary_diamond_chunk_size=100,
        boundary_diamond_random_seed=42,
        threads=1,
        gene_taxonomy_threads=1,
        extended_output=False,
        search_backend="diamond",
        genome_start_time=0.0,
        logger=logging.getLogger(__name__),
        resume_authorized=True,
    )

    resumed = Phase2ResumeState(
        refined_boundaries=result["refined_boundaries"],
        boundary_taxonomy_map=result["boundary_taxonomy_map"],
        boundary_control_stats=result["boundary_control_stats"],
        boundary_diamond_query=result["boundary_diamond_query"],
    )
    _assert_phase3_inputs_equal(original, resumed)
    assert result["proteome_index"] is proteome_index
    assert result["goto_phase3"] is True
