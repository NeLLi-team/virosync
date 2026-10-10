from __future__ import annotations

import copy
import json
import logging
import random
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from virosync.config import PipelineConfig
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
    load_phase2_state,
    phase2_state_to_document,
    write_phase2_state,
)
from virosync.pipeline.phase2.boundary_diamond import (
    ControlStats,
    GeneTaxonomy,
    GenomeDiamondQuery,
    SeedGeneMapping,
)
from virosync.pipeline.phase2.boundary_refiner import (
    RefinedBoundary,
    assign_boundary_candidate_ids,
)
from virosync.pipeline.phase2.repeat_evidence import assess_repeat_evidence
from virosync.pipeline.phase2.terminal_repeats import find_target_site_duplication
from virosync.pipeline.taxonomy_utils import TaxonomyFingerprint


def _boundary() -> RefinedBoundary:
    return RefinedBoundary(
        scaffold="scaffold/alpha",
        start=101,
        end=999,
        candidate_id="EVE_scaffold/alpha_101-999",
        seed_id="seed-b",
        original_start=90,
        original_end=1010,
        candidate_start=80,
        candidate_end=1020,
        pre_tir_start=120,
        pre_tir_end=1000,
        tir_present=True,
        tir_status="detected",
        tir_left_start=101,
        tir_left_end=141,
        tir_right_start=959,
        tir_right_end=999,
        tir_identity=0.95,
        tir_boundary_override=True,
        tir_alignment_capped=False,
        tir_alignment_length=40,
        tir_candidate_count=1,
        tir_scan_start=50,
        tir_scan_end=1050,
        tsd_sequence="GAGGCT",
        repeat_evidence=assess_repeat_evidence(
            None,
            scaffold="scaffold/alpha",
            input_id="seed-b",
            start=101,
            end=999,
            parent_start=120,
            parent_end=1000,
            coverage_intervals=(),
            extension_bp=100,
        ),
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
    assert phase2_state_to_document(observed.refined_boundaries) == phase2_state_to_document(
        expected.refined_boundaries
    )
    assert observed.boundary_taxonomy_map == expected.boundary_taxonomy_map
    assert observed.boundary_control_stats == expected.boundary_control_stats
    assert observed.boundary_diamond_query == expected.boundary_diamond_query
    assert list(observed.boundary_taxonomy_map) == list(expected.boundary_taxonomy_map)
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
    expected_fingerprint = original.boundary_taxonomy_map["porf-under"].taxonomy_fingerprint
    assert loaded.boundary_taxonomy_map["porf-under"].taxonomy_fingerprint == expected_fingerprint
    payload = json.loads(state_path.read_text())
    assert payload["artifact_type"] == PHASE2_RESUME_STATE_ARTIFACT_TYPE
    assert payload["schema_version"] == PHASE2_RESUME_STATE_SCHEMA_VERSION


def test_phase2_checkpoints_and_bed_share_disambiguated_candidate_ids(
    tmp_path: Path,
) -> None:
    boundaries = assign_boundary_candidate_ids(
        [
            RefinedBoundary(
                scaffold="scaffold",
                start=10,
                end=30,
                seed_id="ordinary",
            ),
            RefinedBoundary(
                scaffold="scaffold",
                start=10,
                end=30,
                seed_id="rescue",
                seed_sources=["frameshift_rescue"],
            ),
        ]
    )

    bed_path = phase2._write_phase2_checkpoints(
        output_dir=tmp_path,
        refined_boundaries=boundaries,
        boundary_taxonomy_map={},
        boundary_control_stats=None,
        boundary_diamond_query=None,
    )

    bed_ids = [line.split("\t")[3] for line in bed_path.read_text().splitlines()]
    report_ids = [boundary.candidate_id for boundary in load_phase2_state(tmp_path / "phase2" / PHASE2_STATE_FILENAME)]
    resume_ids = [
        boundary.candidate_id
        for boundary in load_phase2_resume_state(tmp_path / "phase2" / PHASE2_RESUME_STATE_FILENAME).refined_boundaries
    ]
    assert bed_ids == report_ids == resume_ids
    assert len(set(bed_ids)) == 2


def test_phase2_resume_state_round_trip_preserves_optional_none_values() -> None:
    state = Phase2ResumeState(
        refined_boundaries=[],
        boundary_taxonomy_map={},
        boundary_control_stats=None,
        boundary_diamond_query=None,
    )

    loaded = phase2_resume_state_from_document(phase2_resume_state_to_document(state))

    assert loaded == state


def test_phase2_resume_state_rejects_pre_repeat_schema() -> None:
    document = phase2_resume_state_to_document(_state())
    document["schema_version"] = 3

    with pytest.raises(Phase2ResumeStateError, match="unsupported.*schema_version"):
        phase2_resume_state_from_document(document)


def test_repeat_annotation_records_missing_raw_sequence() -> None:
    boundary = RefinedBoundary(scaffold="missing", start=0, end=100, seed_id="source")

    phase2._annotate_boundary_repeats(
        [boundary],
        raw_genome_path=None,
        boundary_diamond_query=None,
        extension_bp=100,
    )

    assert boundary.repeat_evidence.input_id == "source"
    assert boundary.repeat_evidence.left.status == "not_assessed"
    assert boundary.repeat_evidence.right.status == "not_assessed"
    assert boundary.repeat_evidence.candidates == ()
    assert boundary.repeat_evidence.tsd.status == "not_assessed_no_anchor"
    assert (boundary.start, boundary.end, boundary.candidate_id) == (0, 100, "")


def test_repeat_annotation_preserves_retained_tir_and_uses_real_taxonomy_windows(tmp_path: Path) -> None:
    rng = random.Random(135)
    sequence = "".join(rng.choices("ACGT", k=600))
    reverse_arm = sequence[100:160].translate(str.maketrans("ACGT", "TGCA"))[::-1]
    sequence = sequence[:91] + "AAAGAGGCT" + sequence[100:440] + reverse_arm + "GAGGCT" + sequence[506:]
    raw_genome = tmp_path / "raw.fna"
    raw_genome.write_text(f">raw\n{sequence.lower()}\n", encoding="utf-8")
    boundary = RefinedBoundary(
        scaffold="raw",
        start=100,
        end=500,
        seed_id="unsplit-parent",
        pre_tir_start=110,
        pre_tir_end=490,
        tir_present=True,
        tir_status="detected",
        tir_left_start=100,
        tir_left_end=160,
        tir_right_start=440,
        tir_right_end=500,
        tir_identity=1.0,
        tir_alignment_length=60,
        tir_candidate_count=1,
        tir_boundary_override=True,
        tsd_sequence="GAGGCT",
    )
    boundary = assign_boundary_candidate_ids([boundary])[0]
    unchanged = replace(boundary)
    query = GenomeDiamondQuery(
        seed_gene_mappings={
            "left": SeedGeneMapping("left", "raw", 100, 160, flank_start_bp=80, flank_end_bp=190),
            "right": SeedGeneMapping("right", "raw", 440, 500, flank_start_bp=420, flank_end_bp=550),
        }
    )

    phase2._annotate_boundary_repeats(
        [boundary],
        raw_genome_path=raw_genome,
        boundary_diamond_query=query,
        extension_bp=100,
    )

    evidence = boundary.repeat_evidence
    assert evidence.input_id == "unsplit-parent"
    assert (evidence.parent_start, evidence.parent_end) == (110, 490)
    assert (evidence.assessed_start, evidence.assessed_end) == (100, 500)
    assert (evidence.left.window_start, evidence.left.window_end) == (80, 190)
    assert (evidence.right.window_start, evidence.right.window_end) == (420, 550)
    assert evidence.left.status == evidence.right.status == "incomplete"
    assert "taxonomy_clipped" in evidence.left.reasons
    assert "taxonomy_clipped" in evidence.right.reasons
    assert evidence.tsd.status == "assessed_match"
    assert (evidence.tsd.anchor_start, evidence.tsd.anchor_end) == (100, 500)
    assert evidence.tsd.anchor_source == "retained_tir"
    assert evidence.tsd.sequence == find_target_site_duplication(sequence, phase2._retained_terminal_repeat(boundary))
    assert evidence.tsd.sequence == unchanged.tsd_sequence
    assert replace(boundary, repeat_evidence=None) == unchanged


def test_phase2_resume_state_rejects_schema_drift_and_duplicate_mapping_keys() -> None:
    document = phase2_resume_state_to_document(_state())

    unknown_schema = copy.deepcopy(document)
    unknown_schema["schema_version"] = PHASE2_RESUME_STATE_SCHEMA_VERSION + 1
    with pytest.raises(
        Phase2ResumeStateError,
        match="unsupported.*schema_version",
    ):
        phase2_resume_state_from_document(unknown_schema)

    extra_seed_field = copy.deepcopy(document)
    seed_document = extra_seed_field["boundary_diamond_query"]["seed_gene_mappings"][0]["value"]
    seed_document["runtime_class"] = "arbitrary.Type"
    with pytest.raises(Phase2ResumeStateError, match="extra=.*runtime_class"):
        phase2_resume_state_from_document(extra_seed_field)

    duplicate_taxonomy = copy.deepcopy(document)
    duplicate_taxonomy["boundary_taxonomy_map"].append(copy.deepcopy(duplicate_taxonomy["boundary_taxonomy_map"][0]))
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
    proteome_index = {"scaffold/alpha": ["derived-from-authenticated-proteome"]}
    monkeypatch.setattr(
        phase2,
        "build_proteome_index",
        lambda _path: proteome_index,
    )
    config = PipelineConfig().with_overrides(
        resume=True,
        boundary_host_trim_enabled=False,
        boundary_host_trim_window_bp=1000,
        boundary_host_trim_step_bp=500,
        boundary_host_trim_max_host_fraction=0.5,
        boundary_host_trim_min_viral_fraction=0.1,
        boundary_host_trim_score_threshold=0.5,
        boundary_host_trim_buffer_kb=1,
        boundary_host_trim_min_overlap_score=0.2,
        boundary_taxonomy_ml_model="logreg",
        boundary_taxonomy_ml_neighbor_window=1,
        threads=1,
        gene_taxonomy_threads=1,
        extended_output=False,
    )
    config.host.prefixes = []
    config.phase2.diamond_flank_genes = 17
    config.phase2.diamond_control_sample_size = 23
    config.phase2.diamond_control_min_distance = 11
    config.phase2.diamond_chunk_size = 100

    result = phase2._run_phase2_subflow(
        masked_path=tmp_path / "masked.fna",
        proteome_path=tmp_path / "proteome.faa",
        merged_seeds=[object()],
        validated_markers=[],
        host_signature_model=None,
        output_dir=tmp_path,
        genome_id="genome",
        config=config,
        genome_start_time=0.0,
        logger=logging.getLogger(__name__),
        resume_authorized=True,
    )

    resumed = Phase2ResumeState(
        refined_boundaries=result.refined_boundaries,
        boundary_taxonomy_map=result.boundary_taxonomy_map,
        boundary_control_stats=result.boundary_control_stats,
        boundary_diamond_query=result.boundary_diamond_query,
    )
    _assert_phase3_inputs_equal(original, resumed)
    assert result.proteome_index is proteome_index
    assert result.goto_phase3 is True
