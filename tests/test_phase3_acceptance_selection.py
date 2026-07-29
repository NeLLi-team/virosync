from __future__ import annotations

import pytest

from virosync.ablation import AblationID, InterventionCounts
from virosync.pipeline.phase3 import acceptance_selection as selection_module
from virosync.pipeline.phase3.acceptance_selection import select_phase3_acceptance
from virosync.pipeline.phase3.evidence_synthesizer import VerificationResult


def _candidate(
    eve_id: str,
    *,
    length: int,
    confidence_tier: str,
    eve_class: str,
) -> VerificationResult:
    return VerificationResult(
        eve_id=eve_id,
        scaffold="scaffold",
        start=100,
        end=100 + length,
        confidence_tier=confidence_tier,
        region_classification=eve_class,
    )


def _scored_candidates() -> list[VerificationResult]:
    return [
        _candidate("kept", length=6_001, confidence_tier="HIGH", eve_class="NCLDV"),
        _candidate("short", length=4_000, confidence_tier="HIGH", eve_class="NCLDV"),
        _candidate("low_unknown", length=20_000, confidence_tier="LOW", eve_class="UNKNOWN"),
    ]


@pytest.mark.parametrize("ablation_id", list(AblationID)[:-1])
def test_a0_through_a5_use_normal_gate_selection(ablation_id: AblationID) -> None:
    results = _scored_candidates()

    selection = select_phase3_acceptance(results, ablation_id)

    assert selection.ablation_id is ablation_id
    assert selection.all_results == tuple(results)
    assert selection.detailed_results == tuple(results)
    assert selection.canonical_results == (results[0],)
    assert [decision.kept for decision in selection.normal_gate_decisions] == [True, False, False]
    assert selection.candidate_count == 3
    assert selection.intervention_count == 0
    assert selection.changed_count == 0
    assert selection.intervention_counts == InterventionCounts()


def test_a6_keeps_every_scored_candidate_and_retains_counterfactual_decisions() -> None:
    results = _scored_candidates()

    selection = select_phase3_acceptance(results, AblationID.A6)

    assert selection.detailed_results == tuple(results)
    assert selection.canonical_results == tuple(results)
    assert selection.promoted_low_results == ()
    assert [decision.kept for decision in selection.normal_gate_decisions] == [True, False, False]
    assert [decision.reason for decision in selection.normal_gate_decisions] == [
        "ncldv_mirus_high_medium_pass",
        "ncldv_mirus_high_medium_gate",
        "low_unsupported_family",
    ]
    assert [candidate.intervention_applied for candidate in selection.candidates] == [False, True, True]
    assert selection.candidate_count == 3
    assert selection.intervention_count == 2
    assert selection.changed_count == 2
    assert selection.intervention_counts == InterventionCounts(
        opportunities=3,
        interventions=2,
        changed=2,
    )


def test_a6_distinguishes_normal_low_promotion_from_bypass_retention() -> None:
    promoted = _candidate(
        "promoted",
        length=6_001,
        confidence_tier="LOW",
        eve_class="NCLDV",
    )
    promoted.hallmark_count = 2
    promoted.hallmark_genes = ["MCP", "POL"]
    promoted.likely_family = "NCLDV"
    bypassed = _candidate(
        "bypassed",
        length=6_001,
        confidence_tier="LOW",
        eve_class="UNKNOWN",
    )

    selection = select_phase3_acceptance(
        [promoted, bypassed],
        AblationID.A6,
    )

    assert selection.canonical_results == (promoted, bypassed)
    assert selection.promoted_low_results == (promoted,)
    assert [
        candidate.normal_gate_decision.kept
        for candidate in selection.candidates
    ] == [True, False]
    assert selection.intervention_counts == InterventionCounts(
        opportunities=2,
        interventions=1,
        changed=1,
    )


def test_normal_gate_is_evaluated_once_per_candidate(monkeypatch: pytest.MonkeyPatch) -> None:
    results = _scored_candidates()
    real_evaluator = selection_module.evaluate_v2_quality_gate
    evaluated: list[VerificationResult] = []

    def record_evaluation(result: VerificationResult):
        evaluated.append(result)
        return real_evaluator(result)

    monkeypatch.setattr(selection_module, "evaluate_v2_quality_gate", record_evaluation)

    selection = select_phase3_acceptance(results, AblationID.A6)

    assert evaluated == results
    assert selection.canonical_results == tuple(results)
    assert selection.normal_gate_decisions[1].kept is False


def test_empty_a6_selection_has_exact_zero_counts() -> None:
    selection = select_phase3_acceptance([], AblationID.A6)

    assert selection.detailed_results == ()
    assert selection.canonical_results == ()
    assert selection.intervention_counts == InterventionCounts()


def test_selection_requires_closed_ablation_id() -> None:
    with pytest.raises(TypeError, match="AblationID"):
        select_phase3_acceptance(_scored_candidates(), "A6")  # type: ignore[arg-type]
