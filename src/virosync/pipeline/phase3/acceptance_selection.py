"""Select Phase-3 result surfaces and store their acceptance decisions."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace

from virosync.ablation import AblationID, InterventionCounts
from virosync.output_contract import normalize_effective_eve_class
from virosync.pipeline.hallmark_evidence import HallmarkWitness, contained_witnesses, has_spatial_rescue

from .evidence_synthesizer import VerificationResult, tsd_status
from .integration_genes import INTEGRASE_MECHANISMS
from .output_generator import (
    QualityGateDecision,
    evaluate_v2_quality_gate,
    interior_taxonomy_records,
    top10_family_gene_counts,
)

TIER_RESCUE_CLASSES = ("NCLDV", "MIRUS")
TIER_RESCUE_MIN_LENGTH_BP = 5_000
TIER_RESCUE_MIN_VIRAL_GENES = 2
TIER_RESCUE_MAX_HOST_SIGNATURE_FRACTION = 0.5
TIER_RESCUE_MIN_GC_DEVIATION = 0.10


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
        """Return canonical LOW candidates promoted by the normal gate and still LOW."""
        return tuple(
            candidate.result
            for candidate in self.candidates
            if candidate.canonical_kept
            and candidate.normal_gate_decision.promoted_low
            and candidate.result.confidence_tier == "LOW"
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
    return "frameshift_rescue" in (result.seed_sources or []) and bool(result.frameshift_rescue_marker_ids)


def _authenticated_witnesses(result: VerificationResult) -> set[HallmarkWitness]:
    """Return complete spatial witnesses authenticated at the scored coordinates."""
    if not has_spatial_rescue(result):
        return set()
    return set(
        contained_witnesses(
            result.hallmark_evidence,
            result.hallmark_assessment.witnesses,
            result.scaffold,
            result.start,
            result.end,
            result.scoring_locus_ids,
        )
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
            and not has_spatial_rescue(result)
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
        ordinary = [index for index in kept_indices if not _is_active_rescue_branch(results[index])]
        rescued = [index for index in kept_indices if index not in ordinary]
    tier_rank = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}

    def preference(index: int) -> tuple[int, float, str, int, int, str]:
        return (
            -tier_rank.get(results[index].confidence_tier, 0),
            -float(results[index].final_confidence),
            results[index].scaffold,
            results[index].start,
            results[index].end,
            results[index].eve_id,
        )

    def best_index(indices: list[int]) -> int:
        return min(indices, key=preference)

    def overlaps(left: int, right: int) -> bool:
        return (
            results[left].scaffold == results[right].scaffold
            and results[left].start < results[right].end
            and results[left].end > results[right].start
        )

    def rescue_wins(rescue_index: int, ordinary_index: int) -> bool:
        rescue_result = results[rescue_index]
        ordinary_result = results[ordinary_index]
        return tier_rank.get(rescue_result.confidence_tier, 0) > tier_rank.get(ordinary_result.confidence_tier, 0) or (
            rescue_result.confidence_tier == ordinary_result.confidence_tier
            and rescue_result.final_confidence > ordinary_result.final_confidence
        )

    # A rescue branch must beat every ordinary branch it overlaps. One weak,
    # oversized rescue interval therefore cannot erase a separate ordinary EVE
    # through a transitive overlap chain.
    viable_rescues: list[int] = []
    suppressed_by: dict[int, int] = {}
    for rescue_index in rescued:
        overlapping_ordinary = [ordinary_index for ordinary_index in ordinary if overlaps(rescue_index, ordinary_index)]
        blockers = [
            ordinary_index for ordinary_index in overlapping_ordinary if not rescue_wins(rescue_index, ordinary_index)
        ]
        if blockers:
            suppressed_by[rescue_index] = best_index(blockers)
        else:
            viable_rescues.append(rescue_index)

    for ordinary_index in ordinary:
        winners = [rescue_index for rescue_index in viable_rescues if overlaps(rescue_index, ordinary_index)]
        if winners:
            suppressed_by[ordinary_index] = best_index(winners)

    suppressors = set(suppressed_by.values())
    for loser, winner in suppressed_by.items():
        canonical_kept[loser] = False
        selection_outcomes[loser] = f"overlap_suppressed_by:{results[winner].eve_id}"
    for winner in suppressors:
        if canonical_kept[winner]:
            selection_outcomes[winner] = "overlap_selected"

    if suppressed_by:
        covered_witnesses = {
            witness
            for index in kept_indices
            if canonical_kept[index]
            for witness in _authenticated_witnesses(results[index])
        }
        for loser in sorted(suppressed_by, key=preference):
            witnesses = _authenticated_witnesses(results[loser])
            if witnesses - covered_witnesses:
                canonical_kept[loser] = True
                selection_outcomes[loser] = "overlap_retained_unique_witness"
                covered_witnesses.update(witnesses)

        covered_cores = {
            core
            for index in kept_indices
            if canonical_kept[index]
            for core in results[index].confirmed_seeds
            if core.scaffold == results[index].scaffold
            and results[index].start <= core.start
            and core.end <= results[index].end
        }
        for loser in sorted(suppressed_by, key=preference):
            if canonical_kept[loser]:
                continue
            result = results[loser]
            if set(result.confirmed_seeds) - covered_cores:
                canonical_kept[loser] = True
                selection_outcomes[loser] = "overlap_retained_confirmed_seed"
                covered_cores.update(
                    core
                    for core in result.confirmed_seeds
                    if core.scaffold == result.scaffold and result.start <= core.start and core.end <= result.end
                )

    candidates: list[CandidateAcceptance] = []
    for index, result in enumerate(results):
        result.canonical_selection_outcome = selection_outcomes[index]
        result.normal_quality_gate_reason = normal_decisions[index].reason
        candidates.append(
            CandidateAcceptance(
                result=result,
                normal_gate_decision=normal_decisions[index],
                canonical_kept=canonical_kept[index],
                selection_outcome=selection_outcomes[index],
            )
        )
    return AcceptanceSelection(ablation_id=ablation_id, candidates=tuple(candidates))


def apply_tier_rescue(selection: AcceptanceSelection) -> AcceptanceSelection:
    """Rescue NCLDV and Mirus regions without a core hallmark and lift marked LOW calls.

    A gate-rejected region that meets the rescue rule and overlaps no kept call
    becomes LOW, or MEDIUM with an integration mark. Of overlapping rescues, the
    region with more top-10 viral genes, then the longer one, then the smaller
    ``eve_id`` wins. A kept LOW call that meets the rule and has a mark becomes
    MEDIUM. Call this once, after ``select_phase3_acceptance``: it rewrites
    ``confidence_tier`` and ``canonical_selection_outcome`` of the changed
    results. The A6 arm is returned unchanged.
    """
    if selection.ablation_id is AblationID.A6:
        return selection
    candidates = [_promote_marked_low(candidate) for candidate in selection.candidates]
    taken: dict[str, list[tuple[int, int]]] = {}
    for candidate in candidates:
        if candidate.canonical_kept:
            taken.setdefault(candidate.result.scaffold, []).append((candidate.result.start, candidate.result.end))
    rescuable = [
        index
        for index, candidate in enumerate(candidates)
        if candidate.selection_outcome == "normal_gate_rejected" and _meets_rescue_rule(candidate.result)
    ]
    rescuable.sort(key=lambda index: _rescue_rank(candidates[index].result))
    for index in rescuable:
        result = candidates[index].result
        intervals = taken.setdefault(result.scaffold, [])
        if any(result.start < end and start < result.end for start, end in intervals):
            continue
        intervals.append((result.start, result.end))
        result.confidence_tier = "MEDIUM" if integration_marks(result) else "LOW"
        outcome = f"tier_rescue_{result.confidence_tier.lower()}"
        candidates[index] = replace(candidates[index], canonical_kept=True, selection_outcome=outcome)
    for candidate in candidates:
        candidate.result.canonical_selection_outcome = candidate.selection_outcome
    return AcceptanceSelection(ablation_id=selection.ablation_id, candidates=tuple(candidates))


def integration_marks(result: VerificationResult) -> tuple[str, ...]:
    """Return the integration marks that lift a rescue-rule region to MEDIUM.

    Ambiguous TIRs and internal direct repeats do not count: both are common in
    repeat-rich host sequence.
    """
    marks = []
    if result.tir_status == "detected":
        marks.append("terminal_inverted_repeat")
    if tsd_status(result) == "assessed_match":
        marks.append("target_site_duplication")
    if any(hit["mechanism"] in INTEGRASE_MECHANISMS for hit in result.integration_gene_hits):
        marks.append("integrase_gene")
    if result.gc_deviation >= TIER_RESCUE_MIN_GC_DEVIATION:
        marks.append("gc_deviation")
    return tuple(marks)


def viral_top10_genes(result: VerificationResult) -> int:
    """Return the larger of the NCLDV and MIRUS top-10 gene counts; the two are never added."""
    counts = top10_family_gene_counts(interior_taxonomy_records(result))
    return max(counts["NCLDV"], counts["MIRUS"])


def _meets_rescue_rule(result: VerificationResult) -> bool:
    return (
        normalize_effective_eve_class(result.taxonomy_class) in TIER_RESCUE_CLASSES
        and result.length >= TIER_RESCUE_MIN_LENGTH_BP
        and result.host_signature_fraction <= TIER_RESCUE_MAX_HOST_SIGNATURE_FRACTION
        and viral_top10_genes(result) >= TIER_RESCUE_MIN_VIRAL_GENES
    )


def _rescue_rank(result: VerificationResult) -> tuple[int, int, str]:
    return (-viral_top10_genes(result), -result.length, result.eve_id)


def _promote_marked_low(candidate: CandidateAcceptance) -> CandidateAcceptance:
    result = candidate.result
    if not (
        candidate.canonical_kept
        and result.confidence_tier == "LOW"
        and _meets_rescue_rule(result)
        and integration_marks(result)
    ):
        return candidate
    result.confidence_tier = "MEDIUM"
    return replace(candidate, selection_outcome="tier_promoted_medium")
