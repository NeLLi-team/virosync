"""Pure selection of detailed and canonical Phase-3 result surfaces."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from virosync.ablation import AblationID, InterventionCounts

from .evidence_synthesizer import VerificationResult
from .output_generator import QualityGateDecision, evaluate_v2_quality_gate


@dataclass(frozen=True, slots=True)
class CandidateAcceptance:
    """One scored candidate and its normal and selected-arm acceptance states."""

    result: VerificationResult
    normal_gate_decision: QualityGateDecision
    canonical_kept: bool
    selection_outcome: str

    @property
    def intervention_applied(self) -> bool:
        """Return whether the selected arm bypassed a normal rejection."""

        return self.canonical_kept and not self.normal_gate_decision.kept

    @property
    def acceptance_changed(self) -> bool:
        """Return whether canonical inclusion differs from the normal gate."""

        return self.canonical_kept != self.normal_gate_decision.kept


@dataclass(frozen=True, slots=True)
class AcceptanceSelection:
    """Immutable Phase-3 acceptance result with aligned gate evidence."""

    ablation_id: AblationID
    candidates: tuple[CandidateAcceptance, ...]

    @property
    def all_results(self) -> tuple[VerificationResult, ...]:
        """Return every scored candidate in input order."""

        return tuple(candidate.result for candidate in self.candidates)

    @property
    def detailed_results(self) -> tuple[VerificationResult, ...]:
        """Return the all-candidate surface used by detailed output."""

        return self.all_results

    @property
    def canonical_results(self) -> tuple[VerificationResult, ...]:
        """Return candidates selected for canonical output in input order."""

        return tuple(candidate.result for candidate in self.candidates if candidate.canonical_kept)

    @property
    def normal_gate_decisions(self) -> tuple[QualityGateDecision, ...]:
        """Return the stored normal decisions, including A6 counterfactuals."""

        return tuple(candidate.normal_gate_decision for candidate in self.candidates)

    @property
    def promoted_low_results(self) -> tuple[VerificationResult, ...]:
        """Return canonical LOW candidates promoted by the normal gate."""

        return tuple(
            candidate.result
            for candidate in self.candidates
            if candidate.canonical_kept
            and candidate.normal_gate_decision.promoted_low
        )

    @property
    def candidate_count(self) -> int:
        """Return the exact number of scored candidates evaluated."""

        return len(self.candidates)

    @property
    def intervention_count(self) -> int:
        """Return the number of normal rejections bypassed by this selection."""

        return sum(candidate.intervention_applied for candidate in self.candidates)

    @property
    def changed_count(self) -> int:
        """Return the number of canonical inclusion outcomes changed."""

        return sum(candidate.acceptance_changed for candidate in self.candidates)

    @property
    def intervention_counts(self) -> InterventionCounts:
        """Return contract-ready A6 counts, or zeros for non-A6 arms."""

        if self.ablation_id is not AblationID.A6:
            return InterventionCounts()
        return InterventionCounts(
            opportunities=self.candidate_count,
            interventions=self.intervention_count,
            changed=self.changed_count,
        )


def _has_non_rescue_marker(result: VerificationResult) -> bool:
    """Return whether the refined result retains an ordinary marker protein."""

    # hallmark_count can later be raised by coherence evidence categories.
    # hallmark_genes is populated only from marker-bearing proteins, and this
    # helper is called after every rescued marker has left the boundary.
    return bool(result.hallmark_genes)


def _is_active_rescue_branch(result: VerificationResult) -> bool:
    """Return whether rescue provenance still has a marker in the boundary."""

    return (
        "frameshift_rescue" in (result.seed_sources or [])
        and bool(result.frameshift_rescue_marker_ids)
    )


def select_phase3_acceptance(
    results: Sequence[VerificationResult],
    ablation_id: AblationID,
) -> AcceptanceSelection:
    """Evaluate the normal gate once per candidate and select canonical output.

    A0-A5 start from the normal v2 quality-gate result, then enforce
    rescue-marker containment and arbitrate directly overlapping
    ordinary/rescue branches. A6 retains every scored candidate while
    preserving the normal decisions as counterfactual evidence.
    This function performs no I/O and never invokes the output generator's
    legacy fallback filter.
    """

    if not isinstance(ablation_id, AblationID):
        raise TypeError("ablation_id must be an AblationID")

    normal_decisions: list[QualityGateDecision] = []
    canonical_kept: list[bool] = []
    selection_outcomes: list[str] = []
    bypass_gate = ablation_id is AblationID.A6
    for result in results:
        normal_decision = evaluate_v2_quality_gate(result)
        normal_decisions.append(normal_decision)
        kept = normal_decision.kept or bypass_gate
        outcome = "kept" if kept else "normal_gate_rejected"
        if (
            not bypass_gate
            and "frameshift_rescue" in (result.seed_sources or [])
            and not result.frameshift_rescue_marker_ids
            and not _has_non_rescue_marker(result)
        ):
            kept = False
            outcome = "rescue_marker_excluded"
        canonical_kept.append(kept)
        selection_outcomes.append(outcome)

    kept_indices = [index for index, kept in enumerate(canonical_kept) if kept]
    if bypass_gate:
        ordinary: list[int] = []
        rescued: list[int] = []
    else:
        ordinary = [
            index
            for index in kept_indices
            if not _is_active_rescue_branch(results[index])
        ]
        rescued = [index for index in kept_indices if index not in ordinary]
    tier_rank = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}

    def best_index(indices: list[int]) -> int:
        return min(
            indices,
            key=lambda index: (
                -tier_rank.get(results[index].confidence_tier, 0),
                -float(results[index].final_confidence),
                results[index].scaffold,
                results[index].start,
                results[index].end,
                results[index].eve_id,
            ),
        )

    def overlaps(left: int, right: int) -> bool:
        return (
            results[left].scaffold == results[right].scaffold
            and results[left].start < results[right].end
            and results[left].end > results[right].start
        )

    def rescue_wins(rescue_index: int, ordinary_index: int) -> bool:
        rescue_result = results[rescue_index]
        ordinary_result = results[ordinary_index]
        return (
            tier_rank.get(rescue_result.confidence_tier, 0)
            > tier_rank.get(ordinary_result.confidence_tier, 0)
            or (
                rescue_result.confidence_tier
                == ordinary_result.confidence_tier
                and rescue_result.final_confidence
                > ordinary_result.final_confidence
            )
        )

    # A rescue branch must beat every ordinary branch it overlaps. One weak,
    # oversized rescue interval therefore cannot erase a separate ordinary EVE
    # through a transitive overlap chain.
    viable_rescues: list[int] = []
    suppressed_by: dict[int, int] = {}
    for rescue_index in rescued:
        overlapping_ordinary = [
            ordinary_index
            for ordinary_index in ordinary
            if overlaps(rescue_index, ordinary_index)
        ]
        blockers = [
            ordinary_index
            for ordinary_index in overlapping_ordinary
            if not rescue_wins(rescue_index, ordinary_index)
        ]
        if blockers:
            suppressed_by[rescue_index] = best_index(blockers)
        else:
            viable_rescues.append(rescue_index)

    for ordinary_index in ordinary:
        winners = [
            rescue_index
            for rescue_index in viable_rescues
            if overlaps(rescue_index, ordinary_index)
        ]
        if winners:
            suppressed_by[ordinary_index] = best_index(winners)

    suppressors = set(suppressed_by.values())
    for loser, winner in suppressed_by.items():
        canonical_kept[loser] = False
        selection_outcomes[loser] = (
            f"overlap_suppressed_by:{results[winner].eve_id}"
        )
    for winner in suppressors:
        if canonical_kept[winner]:
            selection_outcomes[winner] = "overlap_selected"

    candidates: list[CandidateAcceptance] = []
    for index, result in enumerate(results):
        result.canonical_selection_outcome = selection_outcomes[index]
        candidates.append(
            CandidateAcceptance(
                result=result,
                normal_gate_decision=normal_decisions[index],
                canonical_kept=canonical_kept[index],
                selection_outcome=selection_outcomes[index],
            )
        )
    return AcceptanceSelection(ablation_id=ablation_id, candidates=tuple(candidates))
