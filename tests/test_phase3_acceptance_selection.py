from __future__ import annotations

from dataclasses import replace
from itertools import permutations

import pytest

from virosync.ablation import AblationID, InterventionCounts
from virosync.pipeline.hallmark_evidence import HallmarkEvidence, assess_hallmarks, find_spatial_witnesses
from virosync.pipeline.phase3 import acceptance_selection as selection_module
from virosync.pipeline.phase3.acceptance_selection import select_phase3_acceptance
from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    EvidenceSynthesizerConfig,
    VerificationResult,
)


def _add_scored_core(result: VerificationResult) -> None:
    """Attach a validated, retained catalogue MCP and its native assessment."""
    result.hallmark_evidence = [
        HallmarkEvidence(
            "core|aa1-20",
            "core",
            result.scaffold,
            result.start + 1,
            result.start + 61,
            "GVOGm0003",
            1e-8,
            70.0,
            60.0,
            1,
            20,
            "validated",
            30.0,
            "retained_hmm",
        )
    ]
    result.hallmark_assessment = assess_hallmarks(result)
    assert result.hallmark_assessment.summary.qualifying_evidence_ids


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
    assert [result.normal_quality_gate_reason for result in results] == [
        decision.reason for decision in selection.normal_gate_decisions
    ]
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
    assert [result.to_dict()["normal_quality_gate_reason"] for result in results] == [
        decision.reason for decision in selection.normal_gate_decisions
    ]
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
    _add_scored_core(promoted)
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
    assert [candidate.normal_gate_decision.kept for candidate in selection.candidates] == [True, False]
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


def _rescue_candidate(
    eve_id: str,
    *,
    start: int,
    end: int,
    confidence_tier: str,
    final_confidence: float,
    rescued: bool,
) -> VerificationResult:
    result = VerificationResult(
        eve_id=eve_id,
        scaffold="scaffold",
        start=start,
        end=end,
        confidence_tier=confidence_tier,
        final_confidence=final_confidence,
        region_classification="NCLDV",
        likely_family="NCLDV",
    )
    if rescued:
        result.seed_sources = ["frameshift_rescue", "hhg", "marker_validation"]
        result.frameshift_rescue_marker_ids = [f"{eve_id}_VSR0123456789abcdef"]
        result.hallmark_count = 2
        result.hallmark_genes = ["VS000087", "POLB"]
        _add_scored_core(result)
    return result


def test_overlapping_ordinary_candidate_survives_lower_tier_rescue() -> None:
    ordinary = _rescue_candidate(
        "ordinary",
        start=100,
        end=7_000,
        confidence_tier="MEDIUM",
        final_confidence=0.30,
        rescued=False,
    )
    rescue = _rescue_candidate(
        "rescue",
        start=50,
        end=7_500,
        confidence_tier="LOW",
        final_confidence=0.40,
        rescued=True,
    )

    selection = select_phase3_acceptance([rescue, ordinary], AblationID.A0)

    assert selection.canonical_results == (ordinary,)
    assert rescue.canonical_selection_outcome == "overlap_suppressed_by:ordinary"
    assert ordinary.canonical_selection_outcome == "overlap_selected"
    assert rescue.normal_quality_gate_reason == "ncldv_mirus_low_promoted"
    assert ordinary.normal_quality_gate_reason == "ncldv_mirus_high_medium_pass"


def test_higher_tier_rescue_wins_direct_overlaps_deterministically() -> None:
    ordinary_a = _rescue_candidate(
        "ordinary_a",
        start=100,
        end=7_000,
        confidence_tier="MEDIUM",
        final_confidence=0.35,
        rescued=False,
    )
    rescue = _rescue_candidate(
        "rescue",
        start=200,
        end=7_200,
        confidence_tier="HIGH",
        final_confidence=0.25,
        rescued=True,
    )
    ordinary_b = _rescue_candidate(
        "ordinary_b",
        start=300,
        end=7_300,
        confidence_tier="MEDIUM",
        final_confidence=0.45,
        rescued=False,
    )

    first = select_phase3_acceptance(
        [ordinary_b, rescue, ordinary_a],
        AblationID.A0,
    )
    second = select_phase3_acceptance(
        [ordinary_a, ordinary_b, rescue],
        AblationID.A0,
    )

    assert first.canonical_results == (rescue,)
    assert second.canonical_results == (rescue,)
    assert ordinary_a.canonical_selection_outcome == "overlap_suppressed_by:rescue"
    assert ordinary_b.canonical_selection_outcome == "overlap_suppressed_by:rescue"


def test_weak_rescue_bridge_does_not_remove_nonoverlapping_ordinary_eves() -> None:
    ordinary_a = _rescue_candidate(
        "ordinary_a",
        start=0,
        end=7_000,
        confidence_tier="MEDIUM",
        final_confidence=0.50,
        rescued=False,
    )
    ordinary_b = _rescue_candidate(
        "ordinary_b",
        start=9_000,
        end=16_000,
        confidence_tier="MEDIUM",
        final_confidence=0.20,
        rescued=False,
    )
    rescue = _rescue_candidate(
        "rescue",
        start=5_000,
        end=11_000,
        confidence_tier="MEDIUM",
        final_confidence=0.40,
        rescued=True,
    )

    selection = select_phase3_acceptance(
        [ordinary_b, rescue, ordinary_a],
        AblationID.A0,
    )

    assert selection.canonical_results == (ordinary_b, ordinary_a)
    assert rescue.canonical_selection_outcome == "overlap_suppressed_by:ordinary_a"


def test_rescue_loser_does_not_bridge_nonoverlapping_rescue_alternative() -> None:
    ordinary = _rescue_candidate(
        "ordinary",
        start=0,
        end=7_000,
        confidence_tier="MEDIUM",
        final_confidence=0.50,
        rescued=False,
    )
    rescue_loser = _rescue_candidate(
        "rescue_loser",
        start=6_000,
        end=13_000,
        confidence_tier="MEDIUM",
        final_confidence=0.40,
        rescued=True,
    )
    rescue_alternative = _rescue_candidate(
        "rescue_alternative",
        start=12_000,
        end=19_000,
        confidence_tier="HIGH",
        final_confidence=0.80,
        rescued=True,
    )

    selection = select_phase3_acceptance(
        [rescue_loser, ordinary, rescue_alternative],
        AblationID.A0,
    )

    assert selection.canonical_results == (ordinary, rescue_alternative)
    assert rescue_loser.canonical_selection_outcome == "overlap_suppressed_by:ordinary"


def test_rescue_candidate_that_lost_its_marker_is_detailed_only() -> None:
    rescue = _rescue_candidate(
        "rescue",
        start=100,
        end=7_000,
        confidence_tier="HIGH",
        final_confidence=0.80,
        rescued=True,
    )
    rescue.frameshift_rescue_marker_ids = []
    # Coherence categories can raise the aggregate count after marker
    # processing; they are not proof that an ordinary marker protein remains.
    rescue.hallmark_count = 2
    rescue.hallmark_genes = []

    selection = select_phase3_acceptance([rescue], AblationID.A0)

    assert selection.detailed_results == (rescue,)
    assert selection.canonical_results == ()
    assert rescue.canonical_selection_outcome == "rescue_marker_excluded"


def test_rescue_seed_that_lost_rescue_marker_keeps_ordinary_marker_call() -> None:
    mixed = _rescue_candidate(
        "mixed",
        start=100,
        end=7_000,
        confidence_tier="MEDIUM",
        final_confidence=0.30,
        rescued=True,
    )
    mixed.frameshift_rescue_marker_ids = []
    mixed.hallmark_count = 1
    mixed.hallmark_genes = ["POLB"]

    selection = select_phase3_acceptance([mixed], AblationID.A0)

    assert selection.canonical_results == (mixed,)
    assert mixed.canonical_selection_outcome == "kept"


def test_a6_retains_candidates_despite_shared_output_rules() -> None:
    ordinary = _rescue_candidate(
        "ordinary",
        start=100,
        end=7_000,
        confidence_tier="MEDIUM",
        final_confidence=0.30,
        rescued=False,
    )
    rescue = _rescue_candidate(
        "rescue",
        start=50,
        end=7_500,
        confidence_tier="LOW",
        final_confidence=0.40,
        rescued=True,
    )

    selection = select_phase3_acceptance([rescue, ordinary], AblationID.A6)

    assert selection.canonical_results == (rescue, ordinary)
    assert selection.changed_count == 0
    assert selection.intervention_counts == InterventionCounts(
        opportunities=2,
        interventions=0,
        changed=0,
    )

    rescue.frameshift_rescue_marker_ids = []
    rescue.hallmark_count = 0
    rescue.hallmark_genes = []
    selection = select_phase3_acceptance([rescue], AblationID.A6)

    assert selection.canonical_results == (rescue,)


def test_rescued_protein_id_reaches_verification_result() -> None:
    result = VerificationResult(
        eve_id="rescue",
        scaffold="scaffold",
        start=100,
        end=7_000,
    )
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())

    synthesizer._process_hallmark_hits(
        result,
        [
            {
                "porf_id": "scaffold_VSR0123456789abcdef|aa1-90",
                "hallmark_gene": "VS000087",
                "hmm_score": 90.0,
                "validation_status": "validated",
            }
        ],
    )

    assert result.frameshift_rescue_marker_ids == ["scaffold_VSR0123456789abcdef"]


def _add_weak_capsid(result: VerificationResult) -> None:
    result.hallmark_evidence = [
        HallmarkEvidence(
            "capsid|aa1-20",
            "capsid",
            "scaffold",
            6_000,
            6_300,
            "VP_Penton_1",
            1.0,
            10.0,
            10.0,
            1,
            20,
            "unvalidated",
            0.0,
            "retained_hmm",
        )
    ]
    result.hallmark_witnesses = find_spatial_witnesses(result.hallmark_evidence)
    result.hallmark_assessment = assess_hallmarks(result)
    result.score_components["hallmark_floor"] = 0.0


@pytest.mark.parametrize("order", list(permutations(range(3))))
def test_overlap_retains_unique_capsid_once_in_rank_order(order: tuple[int, ...]) -> None:
    ordinary = _rescue_candidate(
        "ordinary",
        start=7_000,
        end=15_000,
        confidence_tier="MEDIUM",
        final_confidence=0.2,
        rescued=False,
    )
    rescue = _rescue_candidate(
        "rescue",
        start=1_000,
        end=9_000,
        confidence_tier="MEDIUM",
        final_confidence=0.2,
        rescued=True,
    )
    alternative = _rescue_candidate(
        "alternative",
        start=2_000,
        end=9_000,
        confidence_tier="MEDIUM",
        final_confidence=0.2,
        rescued=True,
    )
    _add_weak_capsid(rescue)
    _add_weak_capsid(alternative)
    candidates = [ordinary, rescue, alternative]

    selection = select_phase3_acceptance([candidates[index] for index in order], AblationID.A0)

    assert {result.eve_id for result in selection.canonical_results} == {"ordinary", "rescue"}
    assert rescue.canonical_selection_outcome == "overlap_retained_unique_witness"
    assert alternative.canonical_selection_outcome == "overlap_suppressed_by:ordinary"
    assert selection.candidate_count == 3
    assert (rescue.start, rescue.end, rescue.confidence_tier, rescue.final_confidence) == (1_000, 9_000, "MEDIUM", 0.2)
    assert rescue.normal_quality_gate_reason == "hallmark_spatial_rescue_pass"


@pytest.mark.parametrize("winner_has_evidence", [False, True])
def test_overlap_coverage_requires_authenticated_witness_members(winner_has_evidence: bool) -> None:
    ordinary = _rescue_candidate(
        "ordinary",
        start=0,
        end=15_000,
        confidence_tier="MEDIUM",
        final_confidence=0.2,
        rescued=False,
    )
    rescue = _rescue_candidate(
        "rescue",
        start=1_000,
        end=9_000,
        confidence_tier="MEDIUM",
        final_confidence=0.2,
        rescued=True,
    )
    _add_weak_capsid(ordinary)
    _add_weak_capsid(rescue)
    if not winner_has_evidence:
        ordinary.hallmark_evidence = []

    selection = select_phase3_acceptance([ordinary, rescue], AblationID.A0)

    assert selection.candidates[1].canonical_kept is (not winner_has_evidence)
    assert selection.candidates[0].canonical_kept


def test_overlap_retains_ordinary_loser_with_shared_and_distinct_capsids() -> None:
    ordinary = _rescue_candidate(
        "ordinary",
        start=0,
        end=15_000,
        confidence_tier="MEDIUM",
        final_confidence=0.2,
        rescued=False,
    )
    rescue = _rescue_candidate(
        "rescue",
        start=1_000,
        end=9_000,
        confidence_tier="HIGH",
        final_confidence=0.3,
        rescued=True,
    )
    _add_weak_capsid(ordinary)
    _add_weak_capsid(rescue)
    ordinary.hallmark_evidence.append(
        replace(ordinary.hallmark_evidence[0], query_porf="distinct|aa1-20", gene="distinct", start=10_000, end=10_300)
    )
    ordinary.hallmark_witnesses = find_spatial_witnesses(ordinary.hallmark_evidence)
    ordinary.hallmark_assessment = assess_hallmarks(ordinary)

    selection = select_phase3_acceptance([ordinary, rescue], AblationID.A0)

    assert selection.canonical_results == (ordinary, rescue)
    assert ordinary.canonical_selection_outcome == "overlap_retained_unique_witness"
    assert rescue.canonical_selection_outcome == "overlap_selected"


@pytest.mark.parametrize("mutation", ["missing", "stale", "partial", "forged"])
def test_overlap_does_not_restore_unauthenticated_witness(mutation: str) -> None:
    ordinary = _rescue_candidate(
        "ordinary",
        start=7_000,
        end=15_000,
        confidence_tier="MEDIUM",
        final_confidence=0.2,
        rescued=False,
    )
    rescue = _rescue_candidate(
        "rescue",
        start=1_000,
        end=9_000,
        confidence_tier="MEDIUM",
        final_confidence=0.2,
        rescued=True,
    )
    _add_weak_capsid(rescue)
    if mutation == "missing":
        rescue.hallmark_evidence = []
    elif mutation == "stale":
        rescue.hallmark_assessment = replace(rescue.hallmark_assessment, end=10_000)
    elif mutation == "partial":
        rescue.hallmark_evidence = [replace(rescue.hallmark_evidence[0], end=10_000)]
    else:
        rescue.hallmark_assessment = replace(
            rescue.hallmark_assessment,
            witnesses=(replace(rescue.hallmark_assessment.witnesses[0], rule="forged"),),
        )

    selection = select_phase3_acceptance([rescue, ordinary], AblationID.A0)

    assert selection.canonical_results == (ordinary,)
    assert selection.detailed_results == (rescue, ordinary)
    assert selection.candidates[0].normal_gate_decision.kept is True
    assert rescue.canonical_selection_outcome == "overlap_suppressed_by:ordinary"


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("coverage", ["outside", "different_identity", "one_of_two", "complete", "exact_edge"])
def test_overlap_preserves_assigned_confirmed_core_identity_and_span(coverage: str, reverse: bool) -> None:
    """Only a fully retained identical core permits an overlap loser to disappear."""
    from virosync.pipeline.phase1.seed_merger import ConfirmedSeed

    ordinary = _rescue_candidate(
        "ordinary", start=100, end=15_000, confidence_tier="HIGH", final_confidence=0.8, rescued=False
    )
    rescue = _rescue_candidate(
        "rescue", start=7_000, end=20_000, confidence_tier="HIGH", final_confidence=0.9, rescued=True
    )
    core = ConfirmedSeed("original", "scaffold", 8_000, 8_500)
    ordinary.confirmed_seeds = (core,)
    if coverage == "outside":
        ordinary.confirmed_seeds = (replace(core, start=500, end=1_000),)
    elif coverage == "different_identity":
        rescue.confirmed_seeds = (replace(core, seed_id="other"),)
    elif coverage == "one_of_two":
        ordinary.confirmed_seeds = (core, replace(core, seed_id="second", start=500, end=1_000))
        rescue.confirmed_seeds = (core,)
    elif coverage == "exact_edge":
        core = replace(core, start=rescue.start, end=ordinary.end)
        ordinary.confirmed_seeds = rescue.confirmed_seeds = (core,)
    else:
        rescue.confirmed_seeds = (core,)
    original = (
        ordinary.start,
        ordinary.end,
        ordinary.confidence_tier,
        ordinary.final_confidence,
        ordinary.confirmed_seeds,
    )
    candidates = [ordinary, rescue]
    if reverse:
        candidates.reverse()
    selection = select_phase3_acceptance(candidates, AblationID.A0)
    assert all(decision.kept for decision in selection.normal_gate_decisions)
    assert rescue in selection.canonical_results
    assert (ordinary in selection.canonical_results) is (coverage not in {"complete", "exact_edge"})
    if coverage not in {"complete", "exact_edge"}:
        assert ordinary.canonical_selection_outcome == "overlap_retained_confirmed_seed"
    assert (
        ordinary.start,
        ordinary.end,
        ordinary.confidence_tier,
        ordinary.final_confidence,
        ordinary.confirmed_seeds,
    ) == original
