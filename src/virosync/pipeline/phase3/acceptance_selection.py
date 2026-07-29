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


def select_phase3_acceptance(
    results: Sequence[VerificationResult],
    ablation_id: AblationID,
) -> AcceptanceSelection:
    """Evaluate the normal gate once per candidate and select canonical output.

    A0-A5 retain the normal v2 quality-gate result. A6 retains every scored
    candidate while preserving the same normal decisions as counterfactual
    evidence. This function performs no I/O and never invokes the output
    generator's legacy fallback filter.
    """

    if not isinstance(ablation_id, AblationID):
        raise TypeError("ablation_id must be an AblationID")

    candidates: list[CandidateAcceptance] = []
    bypass_gate = ablation_id is AblationID.A6
    for result in results:
        normal_decision = evaluate_v2_quality_gate(result)
        candidates.append(
            CandidateAcceptance(
                result=result,
                normal_gate_decision=normal_decision,
                canonical_kept=normal_decision.kept or bypass_gate,
            )
        )
    return AcceptanceSelection(ablation_id=ablation_id, candidates=tuple(candidates))
