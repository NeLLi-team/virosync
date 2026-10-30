"""Classification follows unique proteins and the final retained marker span."""

from types import SimpleNamespace

import pytest

from virosync.pipeline.phase1.hhg_seeding import Anchor
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.phase2.seed_extension import merge_reachable_boundaries


def anchor(protein: str, model: str, score: float) -> Anchor:
    return Anchor(protein, "s", 1200, 1500, "+", model, score, 1e-20, validation_status="validated")


def marker(protein: str, model: str, start: int, end: int, score: float = 100.0) -> SimpleNamespace:
    return SimpleNamespace(
        query_porf=protein,
        hmm_target=model,
        scaffold="s",
        start=start,
        end=end,
        hmm_score=score,
        hmm_evalue=1e-20,
        validation_status="validated",
        top10_prefixes="",
        top10_pidents="",
        top10_targets="",
    )


def test_seed_counts_one_protein_across_domain_queries():
    seed = MergedSeed(
        "s", 1000, 2000, anchors=[anchor("gene|aa1-90", "GVOGm0001", 100), anchor("gene|aa100-200", "vp_polb", 50)]
    )
    seed.compute_classification()
    assert seed.predicted_family == "NCLDV"
    assert seed.region_classification_vp_plv_markers == 0


def test_seed_keeps_higher_scoring_non_mcp_classification():
    seed = MergedSeed(
        "s", 1000, 2000, anchors=[anchor("gene|aa1-90", "GVOGm0001", 100), anchor("gene|aa100-200", "vp_mcp", 50)]
    )
    seed.compute_classification()
    assert seed.predicted_family == "NCLDV"
    assert seed.region_classification_vp_plv_markers == 0


def test_final_merged_class_and_counts_use_retained_unique_markers():
    from virosync.pipeline.phase2.seed_extension import refresh_boundary_classification

    left = RefinedBoundary(
        "s",
        1000,
        3000,
        predicted_family="PPV",
        region_classification_vp_plv_markers=2,
        confirmed_seeds=(ConfirmedSeed("left", "s", 1000, 2000),),
    )
    right = RefinedBoundary(
        "s",
        2500,
        4500,
        predicted_family="NCLDV",
        region_classification_ncldv_markers=2,
        confirmed_seeds=(ConfirmedSeed("right", "s", 3500, 4500),),
    )
    boundaries = merge_reachable_boundaries([left, right])
    markers = [
        marker("ppv", "vp_polb", 1200, 1400),
        marker("large", "GVOGm0001", 2700, 2800),
        marker("large", "GVOGm0001", 2700, 2800),
        marker("large2", "GVOGm0002", 3600, 3800),
        marker("outside", "mirus_mcp", 4500, 4600),
    ]
    refresh_boundary_classification(boundaries, markers)
    assert len(boundaries) == 1
    assert boundaries[0].predicted_family == "NCLDV"
    assert boundaries[0].region_classification_ncldv_markers == 2
    assert boundaries[0].region_classification_vp_plv_markers == 1
    assert boundaries[0].region_classification_mirus_markers == 0


def test_final_classification_keeps_compact_cress_policy():
    from virosync.pipeline.phase2.seed_extension import refresh_boundary_classification

    boundary = RefinedBoundary("s", 1000, 2000, predicted_family="CRESS")
    refresh_boundary_classification([boundary], [marker("rep", "cress_rep", 1100, 1300)])
    assert boundary.predicted_family == "CRESS"


def test_merge_keeps_tightest_source_barrier_for_later_proposals():
    first = RefinedBoundary(
        "s",
        1000,
        3000,
        right_stop=4000,
        right_stop_reason="host_stop",
        left_stop=0,
        confirmed_seeds=(ConfirmedSeed("a", "s", 1000, 2000),),
    )
    second = RefinedBoundary(
        "s",
        2500,
        3500,
        right_stop=5000,
        right_stop_reason="support_gap",
        left_stop=500,
        confirmed_seeds=(ConfirmedSeed("b", "s", 2800, 3200),),
    )
    later = RefinedBoundary("s", 3400, 4500, confirmed_seeds=(ConfirmedSeed("c", "s", 4200, 4500),))
    merged = merge_reachable_boundaries([first, second])
    assert merged[0].right_stop == 4000
    assert merged[0].right_stop_reason == "host_stop"
    assert merged[0].left_stop == 500
    assert len(merge_reachable_boundaries([*merged, later])) == 2


def test_tir_child_keeps_compact_cress_classification():
    from virosync.pipeline.phase2.terminal_repeats import _child_boundary

    parent = RefinedBoundary("s", 1000, 3000, predicted_family="CRESS")
    child = _child_boundary(parent, start=1200, end=2800, groups=[], scan_start=0, scan_end=4000)
    assert child.predicted_family == "CRESS"


def test_final_classification_reuses_peve_profile_family():
    from virosync.pipeline.phase2.seed_extension import refresh_boundary_classification

    boundary = RefinedBoundary("s", 1000, 2000, predicted_family="NCLDV")
    hit = marker("peve", "EVE_CAPSID_NOVEL_0008", 1100, 1200)
    hit.top10_prefixes = "EUK-pEVE__"
    hit.top10_pidents = "80"
    refresh_boundary_classification([boundary], [hit])
    assert boundary.predicted_family == "PPV"
    assert boundary.region_classification_vp_plv_markers == 1
    assert boundary.region_classification_ncldv_markers == 0


@pytest.mark.parametrize(
    "ablation,mode,expected",
    [
        ("A0", "relaxed", "UNKNOWN"),
        ("A2", "default", "UNKNOWN"),
        ("A2", "relaxed", "NCLDV"),
        ("A3", "relaxed", "UNKNOWN"),
    ],
)
def test_final_phase2_classification_uses_phase3_marker_role_policy(tmp_path, monkeypatch, ablation, mode, expected):
    import logging

    from virosync.ablation import AblationID
    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome import phase2
    from virosync.orchestration._flows.single_genome.phase3 import _build_boundary_hallmarks
    from virosync.pipeline.phase2.boundary_diamond import GenomeDiamondQuery

    hit = marker("bypassed", "GVOGm0001", 1100, 1200, 80)
    hit.validation_status = "supported"
    seed = MergedSeed(
        "s",
        1000,
        2000,
        seed_id="seed",
        predicted_family="NCLDV",
        confirmed_seeds=(ConfirmedSeed("seed", "s", 1000, 2000),),
    )
    config = PipelineConfig().with_overrides(
        resume=False, ablation_id=AblationID(ablation), assembly_mode=mode, threads=1
    )
    monkeypatch.setattr(phase2, "build_proteome_index", lambda _: {"s": []})
    monkeypatch.setattr(phase2, "_require_phase2b_gene_taxonomy_db", lambda *args, **kwargs: tmp_path / "db")
    monkeypatch.setattr(
        phase2, "_run_boundary_diamond", lambda **kwargs: phase2._BoundaryEvidence({}, None, GenomeDiamondQuery())
    )
    monkeypatch.setattr(phase2, "_recalculate_boundary_composition", lambda *args, **kwargs: None)
    monkeypatch.setattr(phase2, "_annotate_boundary_repeats", lambda *args, **kwargs: None)
    monkeypatch.setattr(phase2, "_write_phase2_checkpoints", lambda **kwargs: tmp_path / "boundaries.bed")
    arguments = dict(
        masked_path=tmp_path / "genome",
        proteome_path=tmp_path / "proteome",
        merged_seeds=[seed],
        validated_markers=[hit],
        host_signature_model=None,
        output_dir=tmp_path,
        genome_id="g",
        genome_start_time=0,
        logger=logging.getLogger(__name__),
    )
    result = phase2._run_phase2_subflow(config=config, **arguments)
    if ablation == "A3":
        reference = phase2._run_phase2_subflow(config=config.with_overrides(ablation_id=AblationID.A0), **arguments)
        selected = result.refined_boundaries[0]
        ordinary = reference.refined_boundaries[0]
        assert (
            selected.predicted_family,
            selected.region_classification_ncldv_markers,
            selected.region_classification_vp_plv_markers,
            selected.region_classification_mirus_markers,
            selected.seed_has_mcp,
        ) == (
            ordinary.predicted_family,
            ordinary.region_classification_ncldv_markers,
            ordinary.region_classification_vp_plv_markers,
            ordinary.region_classification_mirus_markers,
            ordinary.seed_has_mcp,
        )
        assert selected.confirmed_seeds == ordinary.confirmed_seeds
        assert (result.refined_boundaries[0].start, result.refined_boundaries[0].end) == (seed.start, seed.end)
    assert result.refined_boundaries[0].predicted_family == expected
    assert bool(_build_boundary_hallmarks([hit], config)) == (expected == "NCLDV")
    assert hit.validation_status == "supported"
