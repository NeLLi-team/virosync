from __future__ import annotations

import pytest

from virosync.ablation import AblationID, InterventionCounts
from virosync.pipeline.host_signatures import HostSignatureModel
from virosync.pipeline.phase1.seed_merger import MergedSeed
from virosync.pipeline.phase2.boundary_diamond import (
    ControlStats,
    GeneTaxonomy,
    SeedGeneMapping,
)
from virosync.pipeline.phase2.boundary_refiner import (
    RefinedBoundary,
    trim_boundary_by_host_taxonomy,
)
from virosync.pipeline.phase2.host_signature_trim import (
    HostTrimParams,
    trim_seed_by_host_signature,
)
from virosync.pipeline.phase2.taxonomy_seed_refiner import (
    A4HostAwareTaxonomyMLError,
    evaluate_taxonomy_seed_refinement,
    validate_taxonomy_refinement_mode,
)


def _host_model() -> HostSignatureModel:
    return HostSignatureModel(
        token_weights={"euk": 1.0},
        token_counts={"euk": 1},
        max_weight=1.0,
        min_token_length=3,
        host_prefixes=["EUK__"],
        weight_mode="rank",
    )


def _phase2a_record(porf_id: str, start: int, end: int) -> dict[str, object]:
    return {
        "porf_id": porf_id,
        "porf_start": start,
        "porf_end": end,
        "top10_targets": "",
        "top10_bitscores": "",
        "top10_pidents": "",
        "top10_evalues": "",
        "has_ncldv_mirus": True,
        "has_vp_plv": False,
        "has_viral": True,
    }


def _taxonomy(
    porf_id: str,
    start: int,
    end: int,
    prefix: str,
    *,
    viral: bool = False,
) -> GeneTaxonomy:
    return GeneTaxonomy(
        porf_id=porf_id,
        scaffold="ctg",
        start=start,
        end=end,
        top1_target=f"{prefix}target",
        top1_prefix=prefix,
        top1_pident=90.0,
        top10_prefixes=[prefix],
        top10_pidents=[50.0],
        has_ncldv_mirus=viral,
        has_viral=viral,
        has_hit=True,
    )


def test_phase2a_a4_retains_input_and_records_normal_trim_counterfactual() -> None:
    seed = MergedSeed(scaffold="ctg", start=0, end=20)
    records = [
        _phase2a_record("left", 0, 10),
        _phase2a_record("right", 10, 20),
    ]
    params = HostTrimParams(window_bp=10, step_bp=10, buffer_kb=0)

    normal, normal_summary = trim_seed_by_host_signature(
        seed,
        records,
        _host_model(),
        params=params,
    )
    selected, a4_summary = trim_seed_by_host_signature(
        seed,
        records,
        _host_model(),
        params=params,
        ablation_id=AblationID.A4,
    )

    assert (normal.start, normal.end) == (10, 20)
    assert (selected.start, selected.end) == (0, 20)
    assert normal_summary["host_coordinate_change_opportunities"] == 0
    expected = {
        "reason": "a4_host_coordinate_change_bypass",
        "trimmed_start": 0,
        "trimmed_end": 20,
        "counterfactual_trimmed_start": 10,
        "counterfactual_trimmed_end": 20,
        "host_coordinate_change_opportunities": 1,
        "host_coordinate_change_interventions": 1,
        "host_coordinate_change_changed": 1,
    }
    assert {key: a4_summary[key] for key in expected} == expected


def test_a4_taxonomy_refinement_ignores_host_barrier_but_keeps_viral_expansion() -> None:
    seed = MergedSeed(
        scaffold="ctg",
        start=100,
        end=200,
        seed_id="seed-1",
    )
    mapping = SeedGeneMapping(
        seed_id=seed.seed_id,
        scaffold="ctg",
        seed_start=seed.start,
        seed_end=seed.end,
        upstream_porf_ids=["host", "viral"],
        flank_start_bp=0,
        flank_end_bp=300,
    )
    taxonomy_map = {
        "host": _taxonomy("host", 80, 100, "EUK__"),
        "viral": _taxonomy("viral", 40, 60, "NCLDV__", viral=True),
    }

    result = evaluate_taxonomy_seed_refinement(
        [seed],
        taxonomy_map,
        {seed.seed_id: mapping},
        expansion_kb=0,
        ablation_id=AblationID.A4,
    )

    assert (result.counterfactual_seeds[0].start, result.counterfactual_seeds[0].end) == (
        100,
        200,
    )
    assert (result.selected_seeds[0].start, result.selected_seeds[0].end) == (40, 200)
    assert result.intervention_counts == InterventionCounts(
        opportunities=1,
        interventions=1,
        changed=1,
    )


def test_a4_taxonomy_refinement_suppresses_host_edge_contraction() -> None:
    seed = MergedSeed(
        scaffold="ctg",
        start=100,
        end=200,
        seed_id="seed-edge",
    )
    mapping = SeedGeneMapping(
        seed_id=seed.seed_id,
        scaffold="ctg",
        seed_start=seed.start,
        seed_end=seed.end,
        eve_porf_ids=["host-edge", "viral-core"],
        flank_start_bp=0,
        flank_end_bp=300,
    )
    taxonomy_map = {
        "host-edge": _taxonomy("host-edge", 100, 120, "EUK__"),
        "viral-core": _taxonomy(
            "viral-core",
            140,
            160,
            "NCLDV__",
            viral=True,
        ),
    }

    result = evaluate_taxonomy_seed_refinement(
        [seed],
        taxonomy_map,
        {seed.seed_id: mapping},
        expansion_kb=0,
        ablation_id=AblationID.A4,
    )

    assert (result.counterfactual_seeds[0].start, result.counterfactual_seeds[0].end) == (
        120,
        200,
    )
    assert (result.selected_seeds[0].start, result.selected_seeds[0].end) == (
        100,
        200,
    )
    assert result.intervention_counts == InterventionCounts(1, 1, 1)


def test_a4_fails_closed_for_host_feature_dependent_taxonomy_ml() -> None:
    with pytest.raises(
        A4HostAwareTaxonomyMLError,
        match="host-derived features",
    ):
        validate_taxonomy_refinement_mode(
            ablation_id=AblationID.A4,
            taxonomy_ml_enabled=True,
        )

    validate_taxonomy_refinement_mode(
        ablation_id=AblationID.A0,
        taxonomy_ml_enabled=True,
    )


def test_phase2f_a4_retains_boundary_and_records_normal_trim_counterfactual() -> None:
    boundary = RefinedBoundary(
        scaffold="ctg",
        start=0,
        end=200,
        seed_id="seed-1",
    )
    mapping = SeedGeneMapping(
        seed_id=boundary.seed_id,
        scaffold="ctg",
        seed_start=100,
        seed_end=150,
        upstream_porf_ids=["host"],
        eve_porf_ids=["viral"],
        flank_start_bp=0,
        flank_end_bp=200,
    )
    taxonomy_map = {
        "host": _taxonomy("host", 0, 50, "EUK__"),
        "viral": _taxonomy("viral", 100, 150, "NCLDV__", viral=True),
    }

    normal_start, normal_end, normal_stats = trim_boundary_by_host_taxonomy(
        boundary=boundary,
        seed_mapping=mapping,
        taxonomy_map=taxonomy_map,
        control_stats=ControlStats(),
        host_prefix="EUK__",
    )
    selected_start, selected_end, a4_stats = trim_boundary_by_host_taxonomy(
        boundary=boundary,
        seed_mapping=mapping,
        taxonomy_map=taxonomy_map,
        control_stats=ControlStats(),
        host_prefix="EUK__",
        ablation_id=AblationID.A4,
    )

    assert (normal_start, normal_end) == (50, 200)
    assert normal_stats["trimmed"] is True
    assert normal_stats["host_coordinate_change_opportunities"] == 0
    assert (selected_start, selected_end) == (0, 200)
    assert a4_stats["trimmed"] is False
    assert a4_stats["counterfactual_trimmed"] is True
    assert (a4_stats["counterfactual_start"], a4_stats["counterfactual_end"]) == (
        50,
        200,
    )
    assert a4_stats["reason"] == "a4_host_coordinate_change_bypass"
    assert (
        a4_stats["host_coordinate_change_opportunities"],
        a4_stats["host_coordinate_change_interventions"],
        a4_stats["host_coordinate_change_changed"],
    ) == (1, 1, 1)


@pytest.mark.parametrize(
    "call",
    [
        lambda: trim_seed_by_host_signature(
            MergedSeed(scaffold="ctg", start=0, end=1),
            [],
            _host_model(),
            ablation_id="A4",  # type: ignore[arg-type]
        ),
        lambda: evaluate_taxonomy_seed_refinement(
            [],
            {},
            {},
            ablation_id="A4",  # type: ignore[arg-type]
        ),
        lambda: trim_boundary_by_host_taxonomy(
            boundary=RefinedBoundary(scaffold="ctg", start=0, end=1),
            seed_mapping=None,
            taxonomy_map={},
            control_stats=ControlStats(),
            host_prefix="EUK__",
            ablation_id="A4",  # type: ignore[arg-type]
        ),
    ],
)
def test_a4_host_coordinate_surfaces_require_closed_ablation_id(call) -> None:
    with pytest.raises(TypeError, match="AblationID"):
        call()
