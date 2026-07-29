"""Focused contract tests for the A5 composition-evidence ablation."""

from __future__ import annotations

import json
from collections.abc import Callable

import pytest

from virosync.ablation import AblationID
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    EvidenceSynthesizerConfig,
    VerificationResult,
    assign_confidence_tier,
    calculate_eve_confidence,
    evaluate_composition_ablation_effect,
)


def _candidate(*, kfd: float = 0.3, gc_deviation: float = 0.15) -> VerificationResult:
    result = VerificationResult(
        eve_id="EVE_scaffold_0-1000",
        scaffold="scaffold",
        start=0,
        end=1000,
    )
    result.kfd = kfd
    result.gc_deviation = gc_deviation
    result.gene_count = 10
    result.gene_taxonomy_viral_top10 = 2
    result.gene_taxonomy_viral_interior = 2
    result.gene_taxonomy_dominant_fraction = 0.5
    return result


def _score(result: VerificationResult, ablation_id: AblationID) -> float:
    return calculate_eve_confidence(
        result,
        crf_confidence=0.95,
        use_crf_score=True,
        ablation_id=ablation_id,
    )


def test_a5_preserves_raw_fields_but_removes_weight_and_crf_bonus() -> None:
    reference = _candidate()
    reference_score = _score(reference, AblationID.A0)

    selected = _candidate()
    selected_score = _score(selected, AblationID.A5)

    assert (selected.kfd, selected.gc_deviation, selected.cub_deviation) == (
        0.3,
        0.15,
        0.0,
    )
    assert selected.score_components["scores"]["composition"] == pytest.approx(0.7)
    assert reference.score_components["weights"]["composition"] == pytest.approx(0.18)
    assert selected.score_components["weights"]["composition"] == 0.0
    assert reference.score_components["composition_bonus"] == pytest.approx(0.03)
    assert selected.score_components["composition_bonus"] == 0.0
    assert selected.score_components["composition_evidence_active"] is False

    effect = selected.composition_ablation_effect
    assert (effect.opportunities, effect.interventions, effect.changed) == (1, 1, 1)
    assert effect.reference_confidence == pytest.approx(reference_score)
    assert effect.selected_confidence == pytest.approx(selected_score)
    assert effect.composition_score == pytest.approx(0.7)
    assert selected.to_dict()["ablation_id"] == "A5"
    assert selected.to_dict()["composition_ablation_effect"] == effect.to_dict()
    json.dumps(selected.to_dict(), allow_nan=False, sort_keys=True)


def test_reference_preserves_historical_bonus_order_at_high_threshold() -> None:
    result = _candidate()
    result.jelly_roll_confidence_bonus = 0.07
    result.ncldv_completeness_ratio = 0.5
    result.structural_score = 0.543089430894308

    score = calculate_eve_confidence(
        result,
        crf_confidence=0.95,
        use_crf_score=False,
        ablation_id=AblationID.A0,
    )

    assert result.score_components["bonus_total"] == 0.38000000000000006
    assert score == 0.7
    result.final_confidence = score
    assert assign_confidence_tier(result, high=0.7, low=0.2) == "HIGH"


@pytest.mark.parametrize(
    "ablation_id",
    [
        AblationID.A0,
        AblationID.A1,
        AblationID.A2,
        AblationID.A3,
        AblationID.A4,
        AblationID.A6,
    ],
)
def test_only_a5_disables_composition(ablation_id: AblationID) -> None:
    result = _candidate()

    _score(result, ablation_id)

    assert result.ablation_id is ablation_id
    assert result.score_components["weights"]["composition"] == pytest.approx(0.18)
    assert result.score_components["composition_bonus"] == pytest.approx(0.03)
    assert result.score_components["composition_evidence_active"] is True
    assert result.composition_ablation_effect.to_dict() == {
        "opportunities": 0,
        "interventions": 0,
        "changed": 0,
        "composition_score": 0.0,
        "reference_confidence": None,
        "selected_confidence": None,
        "reference_tier": "",
        "selected_tier": "",
    }


def test_a5_removes_composition_driven_top_cap() -> None:
    def strong_candidate() -> VerificationResult:
        result = _candidate()
        result.has_mcp = True
        result.hallmark_count = 5
        result.hallmark_diversity = 5
        result.marker_category_hits = ["mcp"]
        result.marker_complement_score = 1.0
        result.region_classification_ncldv_markers = 3
        result.gene_taxonomy_viral_top10 = 10
        result.gene_taxonomy_viral_interior = 10
        result.gene_taxonomy_dominant_fraction = 1.0
        result.ncldv_completeness_ratio = 1.0
        return result

    reference = strong_candidate()
    reference_score = _score(reference, AblationID.A0)
    selected = strong_candidate()
    selected_score = _score(selected, AblationID.A5)

    assert reference.score_components["composition_cap_active"] is True
    assert selected.score_components["composition_cap_active"] is False
    assert reference.score_components["cap"] == pytest.approx(1.0)
    assert selected.score_components["cap"] == pytest.approx(0.99)
    assert reference_score == pytest.approx(1.0)
    assert selected_score == pytest.approx(0.99)
    assert selected.composition_ablation_effect.changed == 1


def test_a5_zero_composition_has_no_opportunity_or_intervention() -> None:
    reference = _candidate(kfd=0.0, gc_deviation=0.0)
    reference_score = _score(reference, AblationID.A0)
    selected = _candidate(kfd=0.0, gc_deviation=0.0)
    selected_score = _score(selected, AblationID.A5)

    effect = selected.composition_ablation_effect
    assert reference_score == pytest.approx(selected_score)
    assert selected.score_components["scores"]["composition"] == 0.0
    assert selected.score_components["weights"]["composition"] == 0.0
    assert selected.score_components["composition_bonus"] == 0.0
    assert selected.score_components["composition_cap_active"] is False
    assert (effect.opportunities, effect.interventions, effect.changed) == (0, 0, 0)
    assert effect.composition_score == 0.0


def test_a5_records_a_reference_to_selected_tier_flip() -> None:
    result = _candidate()
    result.structural_score = 0.5
    result.clustering_bonus = 0.15

    _score(result, AblationID.A5)

    effect = result.composition_ablation_effect
    assert effect.reference_tier == "HIGH"
    assert effect.selected_tier == "MEDIUM"
    assert effect.changed == 1


def test_counterfactual_effect_is_pure_and_deterministic() -> None:
    result = _candidate()
    selected_score = _score(result, AblationID.A5)
    effect = result.composition_ablation_effect
    before = result.to_dict()

    recalculated = evaluate_composition_ablation_effect(
        ablation_id=AblationID.A5,
        composition_score=effect.composition_score,
        reference_confidence=effect.reference_confidence,
        selected_confidence=selected_score,
        high_threshold=0.7,
        low_threshold=0.2,
    )

    assert recalculated == effect
    assert result.to_dict() == before
    assert evaluate_composition_ablation_effect(
        ablation_id=AblationID.A0,
        composition_score=1.0,
        reference_confidence=1.0,
        selected_confidence=0.0,
        high_threshold=0.7,
        low_threshold=0.2,
    ).to_dict() == {
        "opportunities": 0,
        "interventions": 0,
        "changed": 0,
        "composition_score": 0.0,
        "reference_confidence": None,
        "selected_confidence": None,
        "reference_tier": "",
        "selected_tier": "",
    }


def test_changed_counter_uses_final_post_score_tier_policy() -> None:
    result = _candidate()
    result.gene_count = 2
    result.gene_taxonomy_viral_top10 = 2
    result.gene_taxonomy_viral_interior = 2
    result.gene_taxonomy_dominant_fraction = 1.0
    result.interproscan_score = 1.0
    result.region_classification_ncldv_markers = 3
    result.seed_sources = ["compositional"]
    synthesizer = EvidenceSynthesizer(
        config=EvidenceSynthesizerConfig(
            ablation_id=AblationID.A5,
            use_crf_in_final_score=True,
        )
    )

    synthesizer._calculate_final_decision(
        result,
        RefinedBoundary(
            scaffold="scaffold",
            start=0,
            end=1000,
            confidence=0.95,
        ),
        has_hhg_evidence=False,
    )

    effect = result.composition_ablation_effect
    assert result.final_confidence == pytest.approx(0.199)
    assert result.confidence_tier == "LOW"
    assert effect.reference_confidence == pytest.approx(0.199)
    assert effect.selected_confidence == pytest.approx(0.199)
    assert effect.reference_tier == effect.selected_tier == "LOW"
    assert (effect.opportunities, effect.interventions, effect.changed) == (1, 1, 0)


def _configure_override_path(result: VerificationResult, path: str) -> None:
    if path == "marker_taxonomy":
        result.gene_taxonomy_has_ncldv_mirus = True
        result.hallmark_count = 2
    elif path == "mcp":
        result.has_mcp = True
        result.hallmark_count = 1
    elif path != "final":
        raise AssertionError(f"unknown test path: {path}")


@pytest.mark.parametrize("path", ["marker_taxonomy", "mcp", "final"])
def test_a5_reaches_all_three_scoring_paths(
    monkeypatch: pytest.MonkeyPatch,
    path: str,
) -> None:
    result = _candidate()
    _configure_override_path(result, path)
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig(ablation_id=AblationID.A5))
    boundary = RefinedBoundary(
        scaffold="scaffold",
        start=0,
        end=1000,
        confidence=0.95,
    )

    monkeypatch.setattr(synthesizer, "_initialize_result", lambda _boundary: result)
    no_op: Callable[..., None] = lambda *args, **kwargs: None
    for method_name in (
        "_detect_contig_edge",
        "_process_hallmark_hits",
        "_process_gene_taxonomy",
        "_process_interproscan",
        "_apply_jelly_roll_summary",
        "_run_tiebreakers",
    ):
        monkeypatch.setattr(synthesizer, method_name, no_op)

    observed = synthesizer.verify_eve(boundary, window_features=[])

    assert observed is result
    assert observed.ablation_id is AblationID.A5
    assert observed.score_components["weights"]["composition"] == 0.0
    assert observed.score_components["composition_bonus"] == 0.0
    assert observed.composition_ablation_effect.interventions == 1


def test_config_requires_the_single_ablation_enum() -> None:
    with pytest.raises(TypeError, match="ablation_id must be an AblationID"):
        EvidenceSynthesizerConfig(ablation_id="A5")  # type: ignore[arg-type]
