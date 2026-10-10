"""Collect bounded direct-repeat and anchored short-flank evidence without changing boundaries."""

from __future__ import annotations

import hashlib
import math
from collections import defaultdict
from dataclasses import asdict, dataclass, fields, replace
from types import UnionType
from typing import Any, get_args, get_origin, get_type_hints

from virosync.pipeline.phase2.terminal_repeats import (
    _KMER_BP,
    _LOW_COMPLEXITY_ENTROPY,
    MAX_TIR_ARM_BP,
    MIN_TIR_ARM_BP,
    MIN_TIR_IDENTITY,
    TerminalRepeat,
    _best_group_segment,
    _is_low_complexity,
    _kmer_index,
    _seed_search_intervals,
    _valid_window_groups,
    find_target_site_duplication_at,
)

__all__ = [
    "DirectRepeatCandidate",
    "EndAssessment",
    "RepeatEvidence",
    "TsdAssessment",
    "assess_repeat_evidence",
    "evidence_from_dict",
    "evidence_to_dict",
    "for_changed_endpoints",
    "SEARCH_PARAMETER_FIELDS",
    "FILTER_COUNT_FIELDS",
]

MIN_ARM_BP = MIN_TIR_ARM_BP
MIN_IDENTITY = MIN_TIR_IDENTITY
SEED_BP = _KMER_BP
MAX_SEED_CHAIN_BP = MAX_TIR_ARM_BP - SEED_BP
MIN_SEED_ENTROPY_BITS = _LOW_COMPLEXITY_ENTROPY
MAX_WINDOW_RADIUS_BP = 10_000
MAX_SEED_PAIRINGS = 65_536
MAX_DIAGONALS = 8192
MAX_CANDIDATES = 128
METHOD = f"forward-{SEED_BP}mer-ungapped-score-v1"
_DNA_BASES = frozenset("ACGT")
_NULL_OFFSETS = tuple(range(-20, 0)) + tuple(range(1, 21))
_NULL_GEOMETRY = "translate-both-termini;offsets=-20..-1,+1..+20"
_END_REASONS = frozenset(
    {
        "missing_sequence",
        "missing_coverage",
        "contig_clipped",
        "taxonomy_clipped",
        "window_limit",
        "insufficient_window",
        "ambiguous_bases",
        "low_complexity",
        "partner_unavailable",
        "insufficient_paired_window",
        "seed_pairing_limit",
        "diagonal_limit",
        "candidate_limit",
        "endpoint_changed",
        "midpoint_partition",
    }
)
_TSD_REASONS = frozenset({"contig_clipped", "ambiguous_bases"})
SEARCH_PARAMETER_FIELDS = (
    "method",
    "min_arm_bp",
    "min_identity",
    "seed_bp",
    "max_seed_chain_bp",
    "min_seed_entropy_bits",
    "max_window_radius_bp",
    "max_seed_pairings",
    "max_diagonals",
    "max_candidates",
)
FILTER_COUNT_FIELDS = (
    "filtered_low_complexity_left",
    "filtered_low_complexity_right",
    "filtered_overlapping",
)


@dataclass(frozen=True, slots=True)
class EndAssessment:
    """Coverage of one requested endpoint window, in host coordinates."""

    requested_start: int
    requested_end: int
    window_start: int | None
    window_end: int | None
    status: str
    reasons: tuple[str, ...] = ()
    ambiguous_bp: int = 0


@dataclass(frozen=True, slots=True)
class DirectRepeatCandidate:
    """A same-orientation pair with unresolved host versus viral origin.

    The outer interval includes both arms; the inner interval excludes both.
    Neither geometry has boundary authority or claims a precise insertion site.
    """

    candidate_id: str
    left_start: int
    left_end: int
    right_start: int
    right_end: int
    identity: float
    alignment_length: int
    orientation: str = "direct"
    method: str = METHOD
    interpretation: str = "unresolved"


@dataclass(frozen=True, slots=True)
class TsdAssessment:
    """Legacy exact 3–9 bp comparison and correlated local shifted controls."""

    status: str
    anchor_source: str = "none"
    anchor_start: int | None = None
    anchor_end: int | None = None
    sequence: str = ""
    left_reasons: tuple[str, ...] = ()
    right_reasons: tuple[str, ...] = ()
    null_geometry: str = _NULL_GEOMETRY
    null_assessed: int = 0
    null_matches: int = 0


@dataclass(frozen=True, slots=True)
class RepeatEvidence:
    """Immutable annotation tied to the input EVE and the interval actually assessed."""

    evidence_id: str
    scaffold: str
    input_id: str
    assessed_start: int
    assessed_end: int
    parent_start: int
    parent_end: int
    left: EndAssessment
    right: EndAssessment
    candidates: tuple[DirectRepeatCandidate, ...]
    display_candidate_id: str
    tsd: TsdAssessment
    filtered_low_complexity_left: int = 0
    filtered_low_complexity_right: int = 0
    filtered_overlapping: int = 0
    method: str = METHOD
    min_arm_bp: int = MIN_ARM_BP
    min_identity: float = MIN_IDENTITY
    seed_bp: int = SEED_BP
    max_seed_chain_bp: int = MAX_SEED_CHAIN_BP
    min_seed_entropy_bits: float = MIN_SEED_ENTROPY_BITS
    max_window_radius_bp: int = MAX_WINDOW_RADIUS_BP
    max_seed_pairings: int = MAX_SEED_PAIRINGS
    max_diagonals: int = MAX_DIAGONALS
    max_candidates: int = MAX_CANDIDATES


_EvidenceRecord = RepeatEvidence | EndAssessment | DirectRepeatCandidate | TsdAssessment
_RECORD_HINTS: dict[type[_EvidenceRecord], dict[str, object]] = {
    RepeatEvidence: get_type_hints(RepeatEvidence),
    EndAssessment: get_type_hints(EndAssessment),
    DirectRepeatCandidate: get_type_hints(DirectRepeatCandidate),
    TsdAssessment: get_type_hints(TsdAssessment),
}


def assess_repeat_evidence(
    sequence: str | None,
    *,
    scaffold: str,
    input_id: str,
    start: int,
    end: int,
    parent_start: int,
    parent_end: int,
    coverage_intervals: tuple[tuple[int, int], ...],
    extension_bp: int,
    tir: TerminalRepeat | None = None,
) -> RepeatEvidence:
    """Assess raw host windows around existing cuts; never propose an accepted boundary.

    Only the continuous taxonomy-covered component containing each cut is searched.
    Inward windows stop at the input interval midpoint, so opposite arms bracket
    that midpoint. A midpoint_partition reason records this deliberate method
    constraint and is compatible with completed assessment. Other clipping and
    search limits mark incomplete assessment. The optional TIR supplies the
    existing precise short-flank anchors.
    """
    if start < 0 or start >= end or parent_start < 0 or parent_start >= parent_end or extension_bp < 0:
        raise ValueError("repeat assessment requires nonnegative, nonempty intervals and extension")
    if sequence is not None and end > len(sequence):
        raise ValueError("repeat assessment interval exceeds host sequence")
    evidence_id = _stable_id("re", scaffold, input_id, parent_start, parent_end, start, end)
    coverage = _merge_coverage(coverage_intervals)
    midpoint = (start + end) // 2
    left = _assess_end(
        sequence, start, extension_bp, coverage, (start - extension_bp, min(start + extension_bp, midpoint))
    )
    right = _assess_end(sequence, end, extension_bp, coverage, (max(end - extension_bp, midpoint), end + extension_bp))
    candidates: tuple[DirectRepeatCandidate, ...] = ()
    counts = (0, 0, 0)
    if left.window_start is None and right.window_start is not None:
        right = replace(right, status="incomplete", reasons=tuple(sorted((*right.reasons, "partner_unavailable"))))
    if right.window_start is None and left.window_start is not None:
        left = replace(left, status="incomplete", reasons=tuple(sorted((*left.reasons, "partner_unavailable"))))
    if sequence is not None and left.window_start is not None and right.window_start is not None:
        candidates, limits, counts = _discover_pairs(sequence, left, right, evidence_id)
        if limits:
            left = replace(left, status="incomplete", reasons=tuple(sorted(set(left.reasons) | limits)))
            right = replace(right, status="incomplete", reasons=tuple(sorted(set(right.reasons) | limits)))
    return RepeatEvidence(
        evidence_id=evidence_id,
        scaffold=scaffold,
        input_id=input_id,
        assessed_start=start,
        assessed_end=end,
        parent_start=parent_start,
        parent_end=parent_end,
        left=left,
        right=right,
        candidates=candidates,
        display_candidate_id=candidates[0].candidate_id if candidates else "",
        tsd=_assess_tsd(sequence, tir),
        filtered_low_complexity_left=counts[0],
        filtered_low_complexity_right=counts[1],
        filtered_overlapping=counts[2],
    )


def for_changed_endpoints(
    evidence: RepeatEvidence | None,
    *,
    start: int,
    end: int,
) -> RepeatEvidence | None:
    """Retain parent pairs while invalidating assessments at changed endpoints."""
    if evidence is None or (start, end) == (evidence.assessed_start, evidence.assessed_end):
        return evidence
    left, right = evidence.left, evidence.right
    if start != evidence.assessed_start:
        left = replace(left, status="not_assessed", reasons=tuple(sorted(set(left.reasons) | {"endpoint_changed"})))
    if end != evidence.assessed_end:
        right = replace(right, status="not_assessed", reasons=tuple(sorted(set(right.reasons) | {"endpoint_changed"})))
    return replace(
        evidence,
        left=left,
        right=right,
        display_candidate_id="",
        tsd=TsdAssessment(status="not_assessed_endpoint_changed"),
    )


def _stable_id(prefix: str, *parts: str | int) -> str:
    """Hash length-delimited identity components without traversal-order dependence."""
    encoded = "".join(f"{len(str(part))}:{part}" for part in parts)
    return f"{prefix}-{hashlib.sha256(encoded.encode()).hexdigest()[:20]}"


def _merge_coverage(intervals: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    """Merge adjacent coverage while preserving genuine taxonomy gaps."""
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if start < 0 or start >= end:
            raise ValueError("taxonomy coverage requires nonnegative, nonempty intervals")
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return tuple(merged)


def _assess_end(
    sequence: str | None,
    endpoint: int,
    extension_bp: int,
    coverage: tuple[tuple[int, int], ...],
    partition: tuple[int, int],
) -> EndAssessment:
    """Intersect the requested window with sequence and contiguous taxonomy coverage."""
    requested_start, requested_end = endpoint - extension_bp, endpoint + extension_bp
    if sequence is None:
        return EndAssessment(requested_start, requested_end, None, None, "not_assessed", ("missing_sequence",))
    component = next(((start, end) for start, end in coverage if start <= endpoint <= end), None)
    if component is None:
        return EndAssessment(requested_start, requested_end, None, None, "not_assessed", ("missing_coverage",))
    radius = min(extension_bp, MAX_WINDOW_RADIUS_BP)
    window_start = max(0, endpoint - radius, component[0], partition[0])
    window_end = min(len(sequence), endpoint + radius, component[1], partition[1])
    reasons: list[str] = []
    if partition != (requested_start, requested_end):
        reasons.append("midpoint_partition")
    if requested_start < 0 or requested_end > len(sequence):
        reasons.append("contig_clipped")
    if component[0] > max(0, partition[0]) or component[1] < min(len(sequence), partition[1]):
        reasons.append("taxonomy_clipped")
    if radius < extension_bp:
        reasons.append("window_limit")
    if window_end - window_start < MIN_ARM_BP:
        reasons.append("insufficient_window")
    arm_window = sequence[window_start:window_end].upper()
    ambiguous_bp = sum(base not in _DNA_BASES for base in arm_window)
    if ambiguous_bp:
        reasons.append("ambiguous_bases")
    if arm_window and not ambiguous_bp and _is_low_complexity(arm_window):
        reasons.append("low_complexity")
    return EndAssessment(
        requested_start,
        requested_end,
        window_start,
        window_end,
        "incomplete" if set(reasons) - {"midpoint_partition"} else "completed",
        tuple(sorted(reasons)),
        ambiguous_bp,
    )


def _discover_pairs(
    sequence: str,
    left: EndAssessment,
    right: EndAssessment,
    evidence_id: str,
) -> tuple[tuple[DirectRepeatCandidate, ...], set[str], tuple[int, int, int]]:
    """Search distinct same-orientation arms and retain deterministic alternatives."""
    assert left.window_start is not None and right.window_start is not None
    left_sequence = sequence[left.window_start : left.window_end].upper()
    right_sequence = sequence[right.window_start : right.window_end].upper()
    if min(len(left_sequence), len(right_sequence)) < MIN_ARM_BP:
        return (), {"insufficient_paired_window"}, (0, 0, 0)
    diagonals, limits = _seed_diagonals(left_sequence, right_sequence)
    candidates: dict[tuple[int, int, int, int], DirectRepeatCandidate] = {}
    counts = [0, 0, 0]
    for diagonal, seeds in diagonals.items():
        for left_start, right_start, length, identity in _diagonal_segments(
            left_sequence, right_sequence, diagonal, seeds
        ):
            left_arm = left_sequence[left_start : left_start + length]
            right_arm = right_sequence[right_start : right_start + length]
            left_low, right_low = _is_low_complexity(left_arm), _is_low_complexity(right_arm)
            counts[0] += left_low
            counts[1] += right_low
            if left_low or right_low:
                continue
            absolute_left = left.window_start + left_start
            absolute_right = right.window_start + right_start
            if absolute_left + length >= absolute_right:
                counts[2] += 1
                continue
            coordinates = (absolute_left, absolute_left + length, absolute_right, absolute_right + length)
            candidates[coordinates] = DirectRepeatCandidate(
                _stable_id("dr", evidence_id, *coordinates),
                *coordinates,
                identity,
                length,
            )
    ordered = sorted(
        candidates.values(),
        key=lambda pair: (-pair.alignment_length, -pair.identity, pair.left_start, pair.right_start),
    )
    if len(ordered) > MAX_CANDIDATES:
        limits.add("candidate_limit")
    return tuple(ordered[:MAX_CANDIDATES]), limits, (counts[0], counts[1], counts[2])


def _seed_diagonals(left: str, right: str) -> tuple[dict[int, tuple[int, ...]], set[str]]:
    """Index bounded windows and spend the pairing budget on informative rare seeds first."""
    left_index = _kmer_index(left)
    if not left_index:
        return {}, set()
    right_index = _kmer_index(right)
    diagonals: dict[int, set[int]] = defaultdict(set)
    limits: set[str] = set()
    pairings = 0
    shared = sorted(
        left_index.keys() & right_index.keys(),
        key=lambda kmer: (len(left_index[kmer]) * len(right_index[kmer]), kmer),
    )
    for kmer in shared:
        left_positions, right_positions = left_index[kmer], right_index[kmer]
        count = len(left_positions) * len(right_positions)
        if pairings + count > MAX_SEED_PAIRINGS:
            limits.add("seed_pairing_limit")
            break
        pairings += count
        for left_position in left_positions:
            for right_position in right_positions:
                diagonal = right_position - left_position
                diagonals[diagonal].add(left_position)
    supported = {}
    for diagonal, seeds in diagonals.items():
        ordered_seeds = tuple(sorted(seeds))
        alignment_start = max(0, -diagonal)
        overlap = min(len(left) - alignment_start, len(right) - alignment_start - diagonal)
        if _seed_search_intervals(
            ordered_seeds,
            left_alignment_start=alignment_start,
            overlap_length=overlap,
        ):
            supported[diagonal] = ordered_seeds
    # Single seeds and chains too short to yield a scored arm never consume the
    # diagonal budget. Under caps the rare-seed lexical tiebreak is not strand invariant.
    ranked = sorted(supported.items(), key=lambda item: (-len(item[1]), abs(item[0]), item[0]))
    if len(ranked) > MAX_DIAGONALS:
        limits.add("diagonal_limit")
    return dict(ranked[:MAX_DIAGONALS]), limits


def _diagonal_segments(
    left: str,
    right: str,
    diagonal: int,
    seeds: tuple[int, ...],
) -> list[tuple[int, int, int, float]]:
    """Apply the existing local ungapped scorer to forward-oriented arms."""
    left_start = max(0, -diagonal)
    right_start = left_start + diagonal
    overlap = min(len(left) - left_start, len(right) - right_start)
    segments: list[tuple[int, int, int, float]] = []
    if overlap < MIN_ARM_BP:
        return segments
    for interval_start, interval_end in _seed_search_intervals(
        seeds,
        left_alignment_start=left_start,
        overlap_length=overlap,
    ):
        paired = list(
            zip(
                left[left_start + interval_start : left_start + interval_end],
                right[right_start + interval_start : right_start + interval_end],
                strict=True,
            )
        )
        for valid_start, valid_end in _unambiguous_runs(paired):
            matches = [first == second for first, second in paired[valid_start:valid_end]]
            for group_start, group_end in _valid_window_groups(matches):
                segment = _best_group_segment(matches, group_start, group_end)
                if segment is not None:
                    offset, length, identity = segment
                    position = interval_start + valid_start + offset
                    segments.append((left_start + position, right_start + position, length, identity))
    return segments


def _unambiguous_runs(paired: list[tuple[str, str]]) -> list[tuple[int, int]]:
    """Prevent ambiguous sequence from being bridged as an ordinary mismatch."""
    runs: list[tuple[int, int]] = []
    start = 0
    for position, (left, right) in enumerate(paired):
        if left not in _DNA_BASES or right not in _DNA_BASES:
            if position - start >= MIN_ARM_BP:
                runs.append((start, position))
            start = position + 1
    if len(paired) - start >= MIN_ARM_BP:
        runs.append((start, len(paired)))
    return runs


def _tsd_end_reasons(sequence: str, start: int, end: int) -> tuple[str, ...]:
    """Record full legacy 3–9 bp search coverage, even if a shorter match is found."""
    reasons: list[str] = []
    if start < 0 or end > len(sequence):
        reasons.append("contig_clipped")
    if any(base not in _DNA_BASES for base in sequence[max(0, start) : min(len(sequence), end)].upper()):
        reasons.append("ambiguous_bases")
    return tuple(reasons)


def _assess_tsd(sequence: str | None, tir: TerminalRepeat | None) -> TsdAssessment:
    """Keep the legacy anchored calculation and use shifted sites only as diagnostics."""
    if tir is None:
        return TsdAssessment("not_assessed_no_anchor")
    if sequence is None:
        return TsdAssessment("not_assessed_missing_sequence", "retained_tir", tir.left_start, tir.right_end)
    left_reasons = _tsd_end_reasons(sequence, tir.left_start - 9, tir.left_start)
    right_reasons = _tsd_end_reasons(sequence, tir.right_end, tir.right_end + 9)
    duplication = find_target_site_duplication_at(sequence, tir.left_start, tir.right_end)
    assessed = matches = 0
    for offset in _NULL_OFFSETS:
        left, right = tir.left_start + offset, tir.right_end + offset
        if _tsd_end_reasons(sequence, left - 9, left):
            continue
        if _tsd_end_reasons(sequence, right, right + 9):
            continue
        assessed += 1
        matches += bool(find_target_site_duplication_at(sequence, left, right))
    status = "incomplete" if left_reasons or right_reasons else "assessed_match" if duplication else "assessed_no_match"
    return TsdAssessment(
        status,
        "retained_tir",
        tir.left_start,
        tir.right_end,
        duplication,
        left_reasons,
        right_reasons,
        null_assessed=assessed,
        null_matches=matches,
    )


def evidence_to_dict(evidence: RepeatEvidence) -> dict[str, object]:
    """Serialize the evidence without losing candidates or assessment provenance."""
    return asdict(evidence)


def evidence_from_dict(payload: object) -> RepeatEvidence:
    """Decode an exact, validated evidence record at the resume-file boundary."""
    values = _record_fields(payload, RepeatEvidence)
    values["left"] = _decode_end(values["left"])
    values["right"] = _decode_end(values["right"])
    candidates = values["candidates"]
    if not isinstance(candidates, (list, tuple)):
        raise ValueError("repeat candidates must be an array")
    values["candidates"] = tuple(
        DirectRepeatCandidate(**_record_fields(pair, DirectRepeatCandidate)) for pair in candidates
    )
    tsd_values = _record_fields(values["tsd"], TsdAssessment)
    for name in ("left_reasons", "right_reasons"):
        tsd_values[name] = _string_tuple(tsd_values[name])
    values["tsd"] = TsdAssessment(**tsd_values)
    evidence = RepeatEvidence(**values)
    _validate_evidence(evidence)
    return evidence


def _record_fields(payload: object, record_type: type[_EvidenceRecord]) -> dict[str, Any]:
    """Check dynamic JSON fields before construction and full runtime validation.

    Values remain dynamically typed only at this external schema boundary;
    _validate_evidence checks every constructed field before decoding returns.
    """
    if not isinstance(payload, dict) or set(payload) != {field.name for field in fields(record_type)}:
        raise ValueError(f"invalid {record_type.__name__} fields")
    return dict(payload)


def _string_tuple(value: object) -> tuple[str, ...]:
    """Normalize JSON string arrays while rejecting accidental string iteration."""
    if not isinstance(value, (list, tuple)) or any(not isinstance(item, str) for item in value):
        raise ValueError("repeat assessment reasons must be a string array")
    return tuple(value)


def _decode_end(payload: object) -> EndAssessment:
    values = _record_fields(payload, EndAssessment)
    values["reasons"] = _string_tuple(values["reasons"])
    return EndAssessment(**values)


def _validate_evidence(evidence: RepeatEvidence) -> None:
    """Validate schema types and relationships once at the external file boundary."""
    records: tuple[_EvidenceRecord, ...] = (evidence, evidence.left, evidence.right, evidence.tsd, *evidence.candidates)
    for record in records:
        _validate_record_types(record)
    if not (0 <= evidence.assessed_start < evidence.assessed_end and 0 <= evidence.parent_start < evidence.parent_end):
        raise ValueError("invalid repeat assessment interval")
    expected = _stable_id(
        "re",
        evidence.scaffold,
        evidence.input_id,
        evidence.parent_start,
        evidence.parent_end,
        evidence.assessed_start,
        evidence.assessed_end,
    )
    if evidence.evidence_id != expected:
        raise ValueError("repeat evidence ID does not match its source interval")
    midpoint = (evidence.assessed_start + evidence.assessed_end) // 2
    _validate_end(
        evidence.left,
        evidence.assessed_start,
        (evidence.left.requested_start, min(evidence.left.requested_end, midpoint)),
    )
    _validate_end(
        evidence.right,
        evidence.assessed_end,
        (max(evidence.right.requested_start, midpoint), evidence.right.requested_end),
    )
    if evidence.left.window_end is not None and evidence.left.window_end > midpoint:
        raise ValueError("left repeat window crosses the assessed interval midpoint")
    if evidence.right.window_start is not None and evidence.right.window_start < midpoint:
        raise ValueError("right repeat window crosses the assessed interval midpoint")
    candidate_ids: set[str] = set()
    for pair in evidence.candidates:
        _validate_pair(pair, evidence)
        if pair.candidate_id in candidate_ids:
            raise ValueError("duplicate direct-repeat candidate ID")
        candidate_ids.add(pair.candidate_id)
    if len(evidence.candidates) > evidence.max_candidates:
        raise ValueError("direct-repeat candidates exceed the recorded retention limit")
    if evidence.display_candidate_id and evidence.display_candidate_id not in candidate_ids:
        raise ValueError("repeat display candidate is absent from evidence")
    _validate_tsd(evidence.tsd)
    if "endpoint_changed" in {*evidence.left.reasons, *evidence.right.reasons}:
        if evidence.display_candidate_id or evidence.tsd.status != "not_assessed_endpoint_changed":
            raise ValueError("changed endpoints retain stale repeat display or short-flank assessment")
    expected_parameters = (
        METHOD,
        MIN_ARM_BP,
        MIN_IDENTITY,
        SEED_BP,
        MAX_SEED_CHAIN_BP,
        MIN_SEED_ENTROPY_BITS,
        MAX_WINDOW_RADIUS_BP,
        MAX_SEED_PAIRINGS,
        MAX_DIAGONALS,
        MAX_CANDIDATES,
    )
    if tuple(getattr(evidence, name) for name in SEARCH_PARAMETER_FIELDS) != expected_parameters:
        raise ValueError("unsupported repeat discovery method or parameters")
    if min(getattr(evidence, name) for name in FILTER_COUNT_FIELDS) < 0:
        raise ValueError("negative repeat filter counts")


def _validate_record_types(record: _EvidenceRecord) -> None:
    """Reject coercible but invalid JSON scalar types, including booleans and NaN."""
    annotations = _RECORD_HINTS[type(record)]
    for field in fields(record):
        value, annotation = getattr(record, field.name), annotations[field.name]
        if get_origin(annotation) is tuple:
            continue
        allowed = get_args(annotation) if get_origin(annotation) is UnionType else (annotation,)
        if type(value) not in allowed and not (float in allowed and type(value) is int):
            raise ValueError(f"invalid type for repeat evidence {field.name}")
        if isinstance(value, float) and not math.isfinite(value):
            raise ValueError(f"non-finite repeat evidence {field.name}")


def _validate_end(end: EndAssessment, endpoint: int, method_window: tuple[int, int]) -> None:
    """Check requested geometry, available sequence, and the assessment vocabulary."""
    if end.status not in {"completed", "incomplete", "not_assessed"} or end.ambiguous_bp < 0:
        raise ValueError("invalid repeat endpoint assessment status")
    if not set(end.reasons) <= _END_REASONS:
        raise ValueError("unknown repeat endpoint assessment reason")
    if end.ambiguous_bp > 0 and "ambiguous_bases" not in end.reasons:
        raise ValueError("ambiguous repeat window lacks its assessment reason")
    if end.requested_start > end.requested_end or end.requested_start + end.requested_end != 2 * endpoint:
        raise ValueError("requested repeat window is not centered on its assessed endpoint")
    if (end.window_start is None) != (end.window_end is None):
        raise ValueError("incomplete repeat window coordinates")
    if end.window_start is not None and end.window_end is not None:
        if not (max(0, end.requested_start) <= end.window_start <= end.window_end <= end.requested_end):
            raise ValueError("invalid repeat window coordinates")
        if end.ambiguous_bp > end.window_end - end.window_start:
            raise ValueError("ambiguous base count exceeds searched repeat window")
        partitioned = method_window != (end.requested_start, end.requested_end)
        if ("midpoint_partition" in end.reasons) != partitioned:
            raise ValueError("repeat midpoint partition reason disagrees with requested geometry")
    limits = set(end.reasons) - {"midpoint_partition"}
    if end.status == "completed" and (limits or end.window_start is None):
        raise ValueError("completed repeat search cannot have assessment limits")
    if end.status == "completed" and (end.window_start, end.window_end) != method_window:
        raise ValueError("completed repeat search does not cover its whole method window")
    if end.status == "incomplete" and (not limits or end.window_start is None):
        raise ValueError("incomplete repeat search requires a searched window and a limitation")
    if end.status == "not_assessed":
        if not set(end.reasons) & {"missing_sequence", "missing_coverage", "endpoint_changed"}:
            raise ValueError("unassessed repeat endpoint requires a missing input or changed endpoint")
        if end.window_start is not None and "endpoint_changed" not in end.reasons:
            raise ValueError("unassessed repeat endpoint has an unexplained searched window")


def _validate_pair(pair: DirectRepeatCandidate, evidence: RepeatEvidence) -> None:
    """Require supported pair geometry inside both actual search windows."""
    coordinates = (pair.left_start, pair.left_end, pair.right_start, pair.right_end)
    if not (0 <= pair.left_start < pair.left_end < pair.right_start < pair.right_end):
        raise ValueError("invalid direct-repeat arm coordinates")
    if (
        pair.left_end - pair.left_start != pair.alignment_length
        or pair.right_end - pair.right_start != pair.alignment_length
    ):
        raise ValueError("inconsistent direct-repeat arm lengths")
    if pair.alignment_length < MIN_ARM_BP or not MIN_IDENTITY <= pair.identity <= 1:
        raise ValueError("direct-repeat candidate fails recorded thresholds")
    if pair.candidate_id != _stable_id("dr", evidence.evidence_id, *coordinates):
        raise ValueError("invalid direct-repeat candidate ID")
    if (pair.orientation, pair.method, pair.interpretation) != ("direct", METHOD, "unresolved"):
        raise ValueError("unsupported direct-repeat candidate interpretation or method")
    for end, start, stop in (
        (evidence.left, pair.left_start, pair.left_end),
        (evidence.right, pair.right_start, pair.right_end),
    ):
        if end.window_start is None or end.window_end is None or not end.window_start <= start < stop <= end.window_end:
            raise ValueError("direct-repeat arm is outside its searched window")


def _validate_tsd(tsd: TsdAssessment) -> None:
    """Keep missing-anchor, missing-input, partial and complete short-flank states distinct."""
    unassessed = {"not_assessed_no_anchor", "not_assessed_missing_sequence", "not_assessed_endpoint_changed"}
    if tsd.status not in unassessed | {"incomplete", "assessed_match", "assessed_no_match"}:
        raise ValueError("invalid anchored short-flank assessment status")
    if tsd.null_geometry != _NULL_GEOMETRY:
        raise ValueError("unknown shifted short-flank diagnostic geometry")
    if not 0 <= tsd.null_matches <= tsd.null_assessed <= len(_NULL_OFFSETS):
        raise ValueError("invalid shifted short-flank diagnostic counts")
    reasons = set(tsd.left_reasons) | set(tsd.right_reasons)
    if not reasons <= _TSD_REASONS:
        raise ValueError("unknown short-flank assessment reason")
    if tsd.status in unassessed and (tsd.sequence or reasons or tsd.null_assessed or tsd.null_matches):
        raise ValueError("unassessed short-flank record cannot carry matching results")
    if tsd.status in {"not_assessed_no_anchor", "not_assessed_endpoint_changed"}:
        if tsd.anchor_source != "none" or tsd.anchor_start is not None or tsd.anchor_end is not None:
            raise ValueError("unanchored short-flank record cannot carry anchor coordinates")
        return
    if tsd.anchor_source != "retained_tir" or not (
        tsd.anchor_start is not None and tsd.anchor_end is not None and 0 <= tsd.anchor_start < tsd.anchor_end
    ):
        raise ValueError("invalid short-flank anchor coordinates")
    if tsd.status == "incomplete" and not reasons:
        raise ValueError("incomplete short-flank assessment requires a limitation")
    if tsd.status.startswith("assessed_") and reasons:
        raise ValueError("complete short-flank assessment has an unexplained limitation")
    if tsd.sequence and not (
        3 <= len(tsd.sequence) <= 9 and set(tsd.sequence) <= _DNA_BASES and len(set(tsd.sequence)) > 1
    ):
        raise ValueError("invalid legacy short-flank match")
    if (tsd.status == "assessed_match" and not tsd.sequence) or (tsd.status == "assessed_no_match" and tsd.sequence):
        raise ValueError("short-flank status contradicts its matched sequence")
