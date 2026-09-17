"""Terminal inverted-repeat boundary refinement.

The scanner uses ungapped reverse-complement matches. It scores full supported
segments of at least 50 bp, requires 90% identity, and rejects low-complexity
arms. For alignments longer than 500 bp, the full segment defines the outer
endpoints while the outermost 500 bp are reported as evidence and must also
meet the identity threshold. ``tir_alignment_length`` and
``tir_alignment_capped`` preserve that distinction.
"""

from __future__ import annotations

import math
from bisect import bisect_left
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from pathlib import Path

from Bio import SeqIO

from virosync.pipeline.phase1.viral_markers import get_region_classification_summary
from virosync.pipeline.phase2.boundary_diamond import (
    GeneTaxonomy,
    GenomeDiamondQuery,
    SeedGeneMapping,
    has_identity_qualified_viral_hit,
    missing_boundary_taxonomy_ids,
    pORF,
)
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.phase3.mcp_detection import is_mcp_gene

MIN_TIR_ARM_BP = 50
MAX_TIR_ARM_BP = 500
MIN_TIR_IDENTITY = 0.90

_DNA_BASES = frozenset("ACGT")
_KMER_BP = 7
_MAX_SEED_PAIRINGS = 4096
_LOW_COMPLEXITY_ENTROPY = 1.2
_MAX_LOW_COMPLEXITY_PERIOD = 12
_TANDEM_IDENTITY = 0.8
_MATCH_SCORE = 2
_MISMATCH_SCORE = -4
_MIN_TSD_BP = 3
_MAX_TSD_BP = 9
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True, slots=True)
class TerminalRepeat:
    """One ungapped terminal inverted-repeat candidate."""

    left_start: int
    left_end: int
    right_start: int
    right_end: int
    identity: float
    alignment_length: int
    alignment_capped: bool

    @property
    def arm_length(self) -> int:
        """Return the shared arm length in base pairs."""
        return self.left_end - self.left_start


@dataclass(frozen=True, slots=True)
class _MarkerGroup:
    """Validated hits belonging to one predicted protein."""

    protein_id: str
    start: int
    end: int
    marker_names: frozenset[str]


@dataclass(frozen=True, slots=True)
class _ScanSeedIndex:
    """Informative k-mer positions shared by marker-anchor scans."""

    forward: dict[str, tuple[int, ...]]
    reverse: dict[str, tuple[int, ...]]


def _reverse_complement(sequence: str) -> str:
    """Return the reverse complement of an uppercase DNA sequence."""
    return sequence.translate(_COMPLEMENT)[::-1]


def _sequence_entropy(sequence: str) -> float:
    """Return base-composition Shannon entropy in bits."""
    counts = Counter(sequence)
    length = len(sequence)
    return -sum((count / length) * math.log2(count / length) for count in counts.values())


def _is_low_complexity(sequence: str) -> bool:
    """Return whether an arm is unsuitable as terminal-repeat evidence."""
    if not sequence or not set(sequence) <= _DNA_BASES:
        return True
    most_common_fraction = max(Counter(sequence).values()) / len(sequence)
    if most_common_fraction > 0.8 or _sequence_entropy(sequence) < _LOW_COMPLEXITY_ENTROPY:
        return True
    for period in range(1, min(_MAX_LOW_COMPLEXITY_PERIOD, len(sequence) // 2) + 1):
        periodic_matches = sum(sequence[index] == sequence[index % period] for index in range(period, len(sequence)))
        if periodic_matches / (len(sequence) - period) >= _TANDEM_IDENTITY:
            return True
    return False


def _informative_kmer(kmer: str) -> bool:
    """Return whether a seed is valid and compositionally informative."""
    return not _is_low_complexity(kmer)


def _kmer_index(sequence: str) -> dict[str, tuple[int, ...]]:
    """Index informative k-mer starts in one sequence."""
    positions: dict[str, list[int]] = defaultdict(list)
    for index in range(len(sequence) - _KMER_BP + 1):
        kmer = sequence[index : index + _KMER_BP]
        if _informative_kmer(kmer):
            positions[kmer].append(index)
    return {kmer: tuple(starts) for kmer, starts in positions.items()}


def _build_scan_seed_index(sequence: str) -> _ScanSeedIndex:
    """Build reusable forward and reverse-complement seed indexes."""
    return _ScanSeedIndex(
        forward=_kmer_index(sequence),
        reverse=_kmer_index(_reverse_complement(sequence)),
    )


def _index_prefix(
    index: dict[str, tuple[int, ...]],
    limit: int,
) -> dict[str, tuple[int, ...]]:
    """Restrict an index to k-mers starting before a coordinate."""
    restricted: dict[str, tuple[int, ...]] = {}
    for kmer, starts in index.items():
        cutoff = bisect_left(starts, limit)
        if cutoff:
            restricted[kmer] = starts[:cutoff]
    return restricted


def _matching_diagonals(
    left: str,
    reverse_right: str,
    *,
    scan_seed_index: _ScanSeedIndex | None = None,
) -> dict[int, tuple[int, ...]]:
    """Return alignment diagonals and the exact seeds that support them."""
    if scan_seed_index is None:
        left_index = _kmer_index(left)
        right_index = _kmer_index(reverse_right)
    else:
        left_limit = max(0, len(left) - _KMER_BP + 1)
        right_limit = max(0, len(reverse_right) - _KMER_BP + 1)
        left_index = _index_prefix(scan_seed_index.forward, left_limit)
        right_index = _index_prefix(scan_seed_index.reverse, right_limit)

    diagonal_seeds: dict[int, set[int]] = defaultdict(set)
    for kmer in left_index.keys() & right_index.keys():
        left_positions = left_index[kmer]
        right_positions = right_index[kmer]
        if len(left_positions) * len(right_positions) > _MAX_SEED_PAIRINGS:
            continue
        for left_position in left_positions:
            for right_position in right_positions:
                diagonal_seeds[right_position - left_position].add(left_position)

    supported: dict[int, tuple[int, ...]] = {}
    for diagonal, positions in diagonal_seeds.items():
        ordered = sorted(positions)
        if any(
            later - earlier <= MAX_TIR_ARM_BP - _KMER_BP for earlier, later in zip(ordered, ordered[1:], strict=False)
        ):
            supported[diagonal] = tuple(ordered)
    return supported


def _seed_search_intervals(
    seed_positions: tuple[int, ...],
    *,
    left_alignment_start: int,
    overlap_length: int,
) -> list[tuple[int, int]]:
    """Return merged diagonal intervals that can contain a seeded arm."""
    intervals: list[tuple[int, int]] = []
    alignment_end = left_alignment_start + overlap_length
    for earlier, later in zip(seed_positions, seed_positions[1:], strict=False):
        if later - earlier > MAX_TIR_ARM_BP - _KMER_BP:
            continue
        interval_start = max(
            left_alignment_start,
            later + _KMER_BP - MAX_TIR_ARM_BP,
        )
        interval_end = min(alignment_end, earlier + MAX_TIR_ARM_BP)
        if interval_end - interval_start < MIN_TIR_ARM_BP:
            continue
        relative_interval = (
            interval_start - left_alignment_start,
            interval_end - left_alignment_start,
        )
        if intervals and relative_interval[0] <= intervals[-1][1]:
            intervals[-1] = (
                intervals[-1][0],
                max(intervals[-1][1], relative_interval[1]),
            )
        else:
            intervals.append(relative_interval)
    return intervals


def _passes_identity(mismatches: int, length: int) -> bool:
    """Return whether integer match counts meet the identity threshold."""
    matches = length - mismatches
    return matches * 100 >= round(MIN_TIR_IDENTITY * 100) * length


def _valid_window_groups(matches: list[bool]) -> list[tuple[int, int]]:
    """Group nearby starts of qualifying minimum-length alignments."""
    mismatch_prefix = [0]
    for is_match in matches:
        mismatch_prefix.append(mismatch_prefix[-1] + (not is_match))

    valid_starts = [
        start
        for start in range(len(matches) - MIN_TIR_ARM_BP + 1)
        if _passes_identity(
            mismatch_prefix[start + MIN_TIR_ARM_BP] - mismatch_prefix[start],
            MIN_TIR_ARM_BP,
        )
    ]
    if not valid_starts:
        return []

    groups: list[tuple[int, int]] = []
    group_start = valid_starts[0]
    previous = valid_starts[0]
    for current in valid_starts[1:]:
        if current - previous > MIN_TIR_ARM_BP:
            groups.append((group_start, previous))
            group_start = current
        previous = current
    groups.append((group_start, previous))
    return groups


def _best_group_segment(
    matches: list[bool],
    group_start: int,
    group_end: int,
) -> tuple[int, int, float] | None:
    """Return the score-maximal local segment containing a valid window."""
    mismatch_prefix = [0]
    for is_match in matches:
        mismatch_prefix.append(mismatch_prefix[-1] + (not is_match))

    first_start = 0
    last_start = min(group_end, len(matches) - MIN_TIR_ARM_BP)
    best: tuple[int, int, float] | None = None
    best_key: tuple[int, float, int, int] | None = None
    for start in range(first_start, last_start + 1):
        max_length = len(matches) - start
        for length in range(MIN_TIR_ARM_BP, max_length + 1):
            contained_window_start = max(group_start, start)
            contained_window_end = min(group_end, start + length - MIN_TIR_ARM_BP)
            if contained_window_start > contained_window_end:
                continue
            if not matches[start] or not matches[start + length - 1]:
                continue
            mismatches = mismatch_prefix[start + length] - mismatch_prefix[start]
            if not _passes_identity(mismatches, length):
                continue
            match_count = length - mismatches
            identity = match_count / length
            candidate = (start, length, identity)
            score = match_count * _MATCH_SCORE + mismatches * _MISMATCH_SCORE
            candidate_key = (score, identity, length, -start)
            if best_key is None or candidate_key > best_key:
                best = candidate
                best_key = candidate_key
    return best


def _diagonal_candidates(
    left: str,
    reverse_right: str,
    diagonal: int,
    seed_positions: tuple[int, ...],
    *,
    left_offset: int,
    right_end: int,
) -> list[TerminalRepeat]:
    """Find candidate repeat pairs on one ungapped alignment diagonal."""
    left_alignment_start = max(0, -diagonal)
    right_alignment_start = left_alignment_start + diagonal
    overlap_length = min(
        len(left) - left_alignment_start,
        len(reverse_right) - right_alignment_start,
    )
    if overlap_length < MIN_TIR_ARM_BP:
        return []

    candidates: list[TerminalRepeat] = []
    for interval_start, interval_end in _seed_search_intervals(
        seed_positions,
        left_alignment_start=left_alignment_start,
        overlap_length=overlap_length,
    ):
        matches = [
            left[left_alignment_start + offset] == reverse_right[right_alignment_start + offset]
            and left[left_alignment_start + offset] in _DNA_BASES
            for offset in range(interval_start, interval_end)
        ]
        for group_start, group_end in _valid_window_groups(matches):
            segment = _best_group_segment(matches, group_start, group_end)
            if segment is None:
                continue
            segment_start, full_arm_length, identity = segment
            alignment_offset = interval_start + segment_start
            left_start = left_offset + left_alignment_start + alignment_offset
            reverse_right_start = right_alignment_start + alignment_offset
            candidate_right_end = right_end - reverse_right_start
            arm_length = min(full_arm_length, MAX_TIR_ARM_BP)
            outer_matches = matches[segment_start : segment_start + arm_length]
            outer_mismatches = arm_length - sum(outer_matches)
            if not _passes_identity(outer_mismatches, arm_length):
                continue
            identity = sum(outer_matches) / arm_length
            left_arm = left[
                left_alignment_start + alignment_offset : left_alignment_start + alignment_offset + arm_length
            ]
            if _is_low_complexity(left_arm):
                continue
            candidates.append(
                TerminalRepeat(
                    left_start=left_start,
                    left_end=left_start + arm_length,
                    right_start=candidate_right_end - arm_length,
                    right_end=candidate_right_end,
                    identity=identity,
                    alignment_length=full_arm_length,
                    alignment_capped=full_arm_length > MAX_TIR_ARM_BP,
                )
            )
    return candidates


def find_terminal_inverted_repeats(
    sequence: str,
    *,
    scan_start: int,
    scan_end: int,
    marker_start: int,
    marker_end: int,
    _scan_seed_index: _ScanSeedIndex | None = None,
) -> list[TerminalRepeat]:
    """Find non-low-complexity repeat pairs that bracket a marker span."""
    sequence = sequence.upper()
    if (
        scan_start < 0
        or scan_start > marker_start
        or marker_start >= marker_end
        or marker_end > scan_end
        or scan_end > len(sequence)
    ):
        raise ValueError("terminal-repeat scan coordinates are invalid")

    left = sequence[scan_start:marker_start]
    reverse_right = _reverse_complement(sequence[marker_end:scan_end])
    if len(left) < MIN_TIR_ARM_BP or len(reverse_right) < MIN_TIR_ARM_BP:
        return []

    candidates: dict[tuple[int, int], TerminalRepeat] = {}
    for diagonal, seed_positions in _matching_diagonals(
        left,
        reverse_right,
        scan_seed_index=_scan_seed_index,
    ).items():
        for candidate in _diagonal_candidates(
            left,
            reverse_right,
            diagonal,
            seed_positions,
            left_offset=scan_start,
            right_end=scan_end,
        ):
            _retain_best_outer_pair(candidates, candidate)
    return sorted(
        candidates.values(),
        key=lambda candidate: (
            -candidate.arm_length,
            -candidate.identity,
            candidate.left_start,
            candidate.right_end,
        ),
    )


def _shared_pair_parent(
    parent: RefinedBoundary,
    children: list[RefinedBoundary],
    colliding_intervals: set[tuple[str, int, int]],
) -> RefinedBoundary:
    """Return a parent carrying the best conflicting TIR as evidence."""
    candidate = next(
        (
            child
            for child in children
            if child.tir_boundary_override and (child.scaffold, child.start, child.end) in colliding_intervals
        ),
        next(child for child in children if child.tir_boundary_override),
    )
    return replace(
        parent,
        pre_tir_start=parent.start,
        pre_tir_end=parent.end,
        tir_present=True,
        tir_status="shared_pair",
        tir_left_start=candidate.tir_left_start,
        tir_left_end=candidate.tir_left_end,
        tir_right_start=candidate.tir_right_start,
        tir_right_end=candidate.tir_right_end,
        tir_identity=candidate.tir_identity,
        tir_boundary_override=False,
        tir_alignment_capped=candidate.tir_alignment_capped,
        tir_alignment_length=candidate.tir_alignment_length,
        tir_candidate_count=1,
        tir_scan_start=candidate.tir_scan_start,
        tir_scan_end=candidate.tir_scan_end,
        tsd_sequence=candidate.tsd_sequence,
    )


def _resolve_partition_conflicts(
    parents: list[RefinedBoundary],
    plans: list[list[RefinedBoundary]],
) -> list[RefinedBoundary]:
    """Revert split plans that would create duplicate final coordinates."""
    resolved = [list(plan) for plan in plans]
    while True:
        owners_by_interval: dict[tuple[str, int, int], set[int]] = defaultdict(set)
        for parent_index, children in enumerate(resolved):
            for child in children:
                owners_by_interval[(child.scaffold, child.start, child.end)].add(parent_index)
        colliding_by_parent: dict[int, set[tuple[str, int, int]]] = defaultdict(set)
        for interval, owners in owners_by_interval.items():
            if len(owners) < 2:
                continue
            for parent_index in owners:
                if any(child.tir_boundary_override for child in resolved[parent_index]):
                    colliding_by_parent[parent_index].add(interval)
        if not colliding_by_parent:
            return sorted(
                (child for children in resolved for child in children),
                key=lambda boundary: (boundary.scaffold, boundary.start, boundary.end),
            )
        for parent_index, colliding_intervals in colliding_by_parent.items():
            resolved[parent_index] = [
                _shared_pair_parent(
                    parents[parent_index],
                    resolved[parent_index],
                    colliding_intervals,
                )
            ]


def find_target_site_duplication(sequence: str, candidate: TerminalRepeat) -> str:
    """Return the longest exact direct repeat immediately outside a TIR pair."""
    sequence = sequence.upper()
    available = min(candidate.left_start, len(sequence) - candidate.right_end, _MAX_TSD_BP)
    for length in range(available, _MIN_TSD_BP - 1, -1):
        left = sequence[candidate.left_start - length : candidate.left_start]
        right = sequence[candidate.right_end : candidate.right_end + length]
        if left == right and set(left) <= _DNA_BASES and len(set(left)) > 1:
            return left
    return ""


def _base_protein_id(query_porf: str) -> str:
    """Return the protein ID shared by domain-suffixed marker hits."""
    return query_porf.split("|aa", 1)[0]


def _marker_groups(boundary: RefinedBoundary, markers: list[object]) -> list[_MarkerGroup]:
    """Group validated hits from the current and original parent span."""
    grouped: dict[str, list[object]] = defaultdict(list)
    parent_start = min(boundary.start, boundary.original_start)
    parent_end = max(boundary.end, boundary.original_end)
    for marker in markers:
        if marker.scaffold != boundary.scaffold:
            continue
        if marker.validation_status not in ("validated", "validated_novel"):
            continue
        midpoint = (marker.start + marker.end) // 2
        if not parent_start <= midpoint <= parent_end:
            continue
        grouped[_base_protein_id(marker.query_porf)].append(marker)
    return sorted(
        (
            _MarkerGroup(
                protein_id=protein_id,
                start=min(marker.start for marker in hits),
                end=max(marker.end for marker in hits),
                marker_names=frozenset(marker.hmm_target for marker in hits if getattr(marker, "hmm_target", "")),
            )
            for protein_id, hits in grouped.items()
        ),
        key=lambda group: (group.start, group.end, group.protein_id),
    )


def _candidate_key(candidate: TerminalRepeat) -> tuple[int, int]:
    """Return the biological boundary identity of one repeat pair."""
    return candidate.left_start, candidate.right_end


def _retain_best_outer_pair(
    candidates: dict[tuple[int, int], TerminalRepeat],
    candidate: TerminalRepeat,
) -> None:
    """Keep the longest evidence segment for one pair of outer endpoints."""
    key = _candidate_key(candidate)
    previous = candidates.get(key)
    rank = (candidate.alignment_length, candidate.arm_length, candidate.identity)
    if previous is None or rank > (
        previous.alignment_length,
        previous.arm_length,
        previous.identity,
    ):
        candidates[key] = candidate


def _discover_marker_anchored_candidates(
    sequence: str,
    *,
    scan_start: int,
    scan_end: int,
    groups: list[_MarkerGroup],
) -> list[TerminalRepeat]:
    """Return distinct repeat pairs found around individual marker proteins."""
    candidates: dict[tuple[int, int], TerminalRepeat] = {}
    scan_seed_index = _build_scan_seed_index(sequence[scan_start:scan_end])
    for group in groups:
        for candidate in find_terminal_inverted_repeats(
            sequence,
            scan_start=scan_start,
            scan_end=scan_end,
            marker_start=group.start,
            marker_end=group.end,
            _scan_seed_index=scan_seed_index,
        ):
            _retain_best_outer_pair(candidates, candidate)
    return sorted(
        candidates.values(),
        key=lambda candidate: (
            candidate.left_start,
            candidate.right_end,
            -candidate.arm_length,
            -candidate.identity,
        ),
    )


def _has_candidate_conflict(candidates: list[TerminalRepeat]) -> bool:
    """Return whether any distinct candidate intervals overlap."""
    ordered = sorted(candidates, key=lambda candidate: (candidate.left_start, candidate.right_end))
    return any(
        current.right_end > following.left_start for current, following in zip(ordered, ordered[1:], strict=False)
    )


def _relevant_seed_mappings(
    boundary: RefinedBoundary,
    query: GenomeDiamondQuery,
) -> list[SeedGeneMapping]:
    """Return Phase-2b mappings whose seed intervals feed this boundary."""
    return [
        mapping
        for mapping in query.seed_gene_mappings.values()
        if mapping.scaffold == boundary.scaffold
        and (
            mapping.seed_id == boundary.seed_id
            or (mapping.seed_start < boundary.original_end and mapping.seed_end > boundary.original_start)
        )
    ]


def _scan_bounds(
    boundary: RefinedBoundary,
    mappings: list[SeedGeneMapping],
    *,
    sequence_length: int,
    extension_bp: int,
) -> tuple[int, int] | None:
    """Return the requested flank interval clipped to Phase-2b coverage."""
    if not mappings:
        return None
    starts = [boundary.start, boundary.original_start]
    ends = [boundary.end, boundary.original_end]
    if boundary.candidate_start is not None:
        starts.append(boundary.candidate_start)
    if boundary.candidate_end is not None:
        ends.append(boundary.candidate_end)
    requested_start = max(0, min(starts) - extension_bp)
    requested_end = min(sequence_length, max(ends) + extension_bp)
    coverage_start = min(mapping.flank_start_bp for mapping in mappings)
    coverage_end = max(mapping.flank_end_bp for mapping in mappings)
    scan_start = max(requested_start, coverage_start)
    scan_end = min(requested_end, coverage_end)
    if scan_start >= scan_end:
        return None
    return scan_start, scan_end


def _candidate_metadata(
    boundary: RefinedBoundary,
    candidate: TerminalRepeat,
    *,
    status: str,
    candidate_count: int,
    scan_start: int,
    scan_end: int,
    tsd_sequence: str,
) -> RefinedBoundary:
    """Attach one candidate and its scan provenance without changing bounds."""
    return replace(
        boundary,
        pre_tir_start=boundary.start,
        pre_tir_end=boundary.end,
        tir_present=True,
        tir_status=status,
        tir_left_start=candidate.left_start,
        tir_left_end=candidate.left_end,
        tir_right_start=candidate.right_start,
        tir_right_end=candidate.right_end,
        tir_identity=candidate.identity,
        tir_alignment_capped=candidate.alignment_capped,
        tir_alignment_length=candidate.alignment_length,
        tir_candidate_count=candidate_count,
        tir_scan_start=scan_start,
        tir_scan_end=scan_end,
        tsd_sequence=tsd_sequence,
    )


def _groups_within(
    groups: list[_MarkerGroup],
    start: int,
    end: int,
) -> list[_MarkerGroup]:
    """Return marker proteins wholly contained by an interval."""
    return [group for group in groups if start <= group.start and group.end <= end]


def _child_boundary(
    parent: RefinedBoundary,
    *,
    start: int,
    end: int,
    groups: list[_MarkerGroup],
    scan_start: int,
    scan_end: int,
) -> RefinedBoundary:
    """Create a coordinate-local child without inherited aggregate evidence."""
    marker_names = set().union(*(group.marker_names for group in groups)) if groups else set()
    classification = get_region_classification_summary(marker_names, None)
    seed_sources = sorted(set(parent.seed_sources or []) & {"hhg", "marker_validation"}) if groups else []
    if "frameshift_rescue" in (parent.seed_sources or []):
        from virosync.pipeline.phase1.frameshift_screening import (
            is_rescued_protein_id,
        )

        if any(is_rescued_protein_id(group.protein_id) for group in groups):
            seed_sources.append("frameshift_rescue")

    return RefinedBoundary(
        scaffold=parent.scaffold,
        start=start,
        end=end,
        seed_id=parent.seed_id,
        original_start=start,
        original_end=end,
        pre_tir_start=parent.start,
        pre_tir_end=parent.end,
        tir_status="not_detected",
        tir_scan_start=scan_start,
        tir_scan_end=scan_end,
        seed_sources=seed_sources,
        seed_confidence="low",
        seed_has_mcp=any(is_mcp_gene(marker_name) for marker_name in marker_names),
        predicted_family=classification["classification"],
        region_classification_ncldv_markers=classification["ncldv_markers"],
        region_classification_vp_plv_markers=classification["vp_plv_markers"],
        region_classification_mirus_markers=classification["mirus_markers"],
    )


def _residual_intervals(
    parent: RefinedBoundary,
    candidates: list[TerminalRepeat],
) -> list[tuple[int, int]]:
    """Return the parent slices outside disjoint repeat-bounded children."""
    intervals: list[tuple[int, int]] = []
    cursor = parent.start
    for candidate in sorted(candidates, key=lambda item: item.left_start):
        candidate_start = max(parent.start, candidate.left_start)
        candidate_end = min(parent.end, candidate.right_end)
        if candidate_start >= candidate_end:
            continue
        if cursor < candidate_start:
            intervals.append((cursor, candidate_start))
        cursor = max(cursor, candidate_end)
    if cursor < parent.end:
        intervals.append((cursor, parent.end))
    return intervals


def _has_residual_evidence(
    *,
    scaffold: str,
    start: int,
    end: int,
    groups: list[_MarkerGroup],
    taxonomy_map: dict[str, GeneTaxonomy],
) -> bool:
    """Return whether a non-TIR slice retains local viral evidence."""
    if any(group.start < end and group.end > start for group in groups):
        return True
    return any(
        taxonomy.scaffold == scaffold
        and start <= taxonomy.start
        and taxonomy.end <= end
        and has_identity_qualified_viral_hit(
            taxonomy.top10_prefixes,
            taxonomy.top10_pidents,
        )
        for taxonomy in taxonomy_map.values()
    )


def _partition_is_complete(
    *,
    children: list[RefinedBoundary],
    groups: list[_MarkerGroup],
    taxonomy_map: dict[str, GeneTaxonomy],
    proteome_index: dict[str, list[pORF]],
) -> bool:
    """Return whether a split preserves markers and independent evidence."""
    for group in groups:
        owners = [child for child in children if child.start <= group.start and group.end <= child.end]
        if len(owners) != 1:
            return False

    return all(
        not missing_boundary_taxonomy_ids(
            scaffold=child.scaffold,
            start=child.start,
            end=child.end,
            taxonomy_map=taxonomy_map,
            proteome_index=proteome_index,
        )
        for child in children
    )


def _build_partition(
    parent: RefinedBoundary,
    candidates: list[TerminalRepeat],
    *,
    sequence: str,
    groups: list[_MarkerGroup],
    scan_start: int,
    scan_end: int,
    taxonomy_map: dict[str, GeneTaxonomy],
    proteome_index: dict[str, list[pORF]],
) -> list[RefinedBoundary] | None:
    """Build TIR children and evidence-bearing residuals for one parent."""
    children: list[RefinedBoundary] = []
    for candidate in candidates:
        child_groups = _groups_within(
            groups,
            candidate.left_start,
            candidate.right_end,
        )
        child = _child_boundary(
            parent,
            start=candidate.left_start,
            end=candidate.right_end,
            groups=child_groups,
            scan_start=scan_start,
            scan_end=scan_end,
        )
        children.append(
            replace(
                _candidate_metadata(
                    child,
                    candidate,
                    status="detected",
                    candidate_count=1,
                    scan_start=scan_start,
                    scan_end=scan_end,
                    tsd_sequence=find_target_site_duplication(sequence, candidate),
                ),
                start=candidate.left_start,
                end=candidate.right_end,
                pre_tir_start=parent.start,
                pre_tir_end=parent.end,
                tir_boundary_override=True,
            )
        )

    for start, end in _residual_intervals(parent, candidates):
        if not _has_residual_evidence(
            scaffold=parent.scaffold,
            start=start,
            end=end,
            groups=groups,
            taxonomy_map=taxonomy_map,
        ):
            continue
        children.append(
            _child_boundary(
                parent,
                start=start,
                end=end,
                groups=_groups_within(groups, start, end),
                scan_start=scan_start,
                scan_end=scan_end,
            )
        )

    children.sort(key=lambda boundary: (boundary.start, boundary.end))
    if not _partition_is_complete(
        children=children,
        groups=groups,
        taxonomy_map=taxonomy_map,
        proteome_index=proteome_index,
    ):
        return None
    return children


def _candidate_evidence_parent(
    parent: RefinedBoundary,
    candidate: TerminalRepeat,
    *,
    status: str,
    candidate_count: int,
    sequence: str,
    scan_start: int,
    scan_end: int,
) -> RefinedBoundary:
    """Attach non-authoritative candidate evidence to an unsplit parent."""
    return _candidate_metadata(
        parent,
        candidate,
        status=status,
        candidate_count=candidate_count,
        scan_start=scan_start,
        scan_end=scan_end,
        tsd_sequence=find_target_site_duplication(sequence, candidate),
    )


def _plan_boundary_partition(
    boundary: RefinedBoundary,
    *,
    sequence: str,
    validated_markers: list[object],
    boundary_diamond_query: GenomeDiamondQuery,
    boundary_taxonomy_map: dict[str, GeneTaxonomy],
    proteome_index: dict[str, list[pORF]],
    extension_bp: int,
) -> list[RefinedBoundary]:
    """Return one unchanged parent or its independent child boundaries."""
    mappings = _relevant_seed_mappings(boundary, boundary_diamond_query)
    bounds = _scan_bounds(
        boundary,
        mappings,
        sequence_length=len(sequence),
        extension_bp=extension_bp,
    )
    if bounds is None:
        return [
            replace(
                boundary,
                pre_tir_start=boundary.start,
                pre_tir_end=boundary.end,
            )
        ]
    scan_start, scan_end = bounds
    groups = _marker_groups(boundary, validated_markers)
    if not groups:
        return [
            replace(
                boundary,
                pre_tir_start=boundary.start,
                pre_tir_end=boundary.end,
                tir_status="no_marker_anchor",
                tir_scan_start=scan_start,
                tir_scan_end=scan_end,
            )
        ]
    if any(group.start < scan_start or group.end > scan_end for group in groups):
        return [
            replace(
                boundary,
                pre_tir_start=boundary.start,
                pre_tir_end=boundary.end,
                tir_status="taxonomy_incomplete",
                tir_scan_start=scan_start,
                tir_scan_end=scan_end,
            )
        ]

    candidates = _discover_marker_anchored_candidates(
        sequence,
        scan_start=scan_start,
        scan_end=scan_end,
        groups=groups,
    )
    if not candidates:
        return [
            replace(
                boundary,
                pre_tir_start=boundary.start,
                pre_tir_end=boundary.end,
                tir_status="not_detected",
                tir_scan_start=scan_start,
                tir_scan_end=scan_end,
            )
        ]

    best = min(
        candidates,
        key=lambda candidate: (
            -candidate.arm_length,
            -candidate.identity,
            candidate.left_start,
            candidate.right_end,
        ),
    )
    if _has_candidate_conflict(candidates):
        return [
            _candidate_evidence_parent(
                boundary,
                best,
                status="ambiguous",
                candidate_count=len(candidates),
                sequence=sequence,
                scan_start=scan_start,
                scan_end=scan_end,
            )
        ]

    covered_candidates = [
        candidate
        for candidate in candidates
        if not missing_boundary_taxonomy_ids(
            scaffold=boundary.scaffold,
            start=candidate.left_start,
            end=candidate.right_end,
            taxonomy_map=boundary_taxonomy_map,
            proteome_index=proteome_index,
        )
    ]
    if not covered_candidates:
        return [
            _candidate_evidence_parent(
                boundary,
                best,
                status="taxonomy_incomplete",
                candidate_count=len(candidates),
                sequence=sequence,
                scan_start=scan_start,
                scan_end=scan_end,
            )
        ]

    partition = _build_partition(
        boundary,
        covered_candidates,
        sequence=sequence,
        groups=groups,
        scan_start=scan_start,
        scan_end=scan_end,
        taxonomy_map=boundary_taxonomy_map,
        proteome_index=proteome_index,
    )
    if partition is None:
        return [
            _candidate_evidence_parent(
                boundary,
                best,
                status="ambiguous",
                candidate_count=len(candidates),
                sequence=sequence,
                scan_start=scan_start,
                scan_end=scan_end,
            )
        ]
    return partition


def refine_boundaries_with_terminal_repeats(
    boundaries: list[RefinedBoundary],
    *,
    raw_genome_path: Path,
    validated_markers: list[object],
    boundary_diamond_query: GenomeDiamondQuery,
    boundary_taxonomy_map: dict[str, GeneTaxonomy],
    proteome_index: dict[str, list[pORF]],
    extension_bp: int,
) -> list[RefinedBoundary]:
    """Partition boundaries around independent, taxonomy-covered TIR pairs."""
    scaffold_index = SeqIO.index(str(raw_genome_path), "fasta")
    try:
        plans: list[list[RefinedBoundary]] = []
        for boundary in boundaries:
            if boundary.scaffold not in scaffold_index:
                raise ValueError(f"boundary scaffold absent from raw genome: {boundary.scaffold}")
            sequence = str(scaffold_index[boundary.scaffold].seq).upper()
            plans.append(
                _plan_boundary_partition(
                    boundary,
                    sequence=sequence,
                    validated_markers=validated_markers,
                    boundary_diamond_query=boundary_diamond_query,
                    boundary_taxonomy_map=boundary_taxonomy_map,
                    proteome_index=proteome_index,
                    extension_bp=extension_bp,
                )
            )
        return _resolve_partition_conflicts(boundaries, plans)
    finally:
        scaffold_index.close()
