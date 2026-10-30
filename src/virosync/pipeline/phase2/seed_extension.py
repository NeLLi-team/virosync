"""One coordinate authority for growth from complete confirmed Phase 1 seeds."""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from dataclasses import replace
from typing import TYPE_CHECKING, TypeAlias

from virosync.ablation import AblationID
from virosync.pipeline.host_signatures import HostSignatureModel, host_lineage_share
from virosync.pipeline.phase1.marker_roles import DEFAULT_A2_SINGLE_MARKER_MIN_SCORE, decide_marker_hit_role
from virosync.pipeline.phase1.seed_merger import (
    ConfirmedSeed,
    MergedSeed,
    classify_marker_hits,
    validate_confirmed_seed_containment,
)
from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy, has_identity_qualified_viral_hit, pORF
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary

if TYPE_CHECKING:
    from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit


GeneEvidence: TypeAlias = tuple[pORF, bool, bool, bool]


def confirmed_seed_core(seed: MergedSeed) -> MergedSeed:
    """Require original source membership and return its complete coordinate hull."""
    if not seed.confirmed_seeds:
        raise ValueError("candidate has no confirmed seed membership")
    core = replace(
        seed,
        start=min(source.start for source in seed.confirmed_seeds),
        end=max(source.end for source in seed.confirmed_seeds),
    )
    validate_confirmed_seed_containment(core.scaffold, core.start, core.end, core.confirmed_seeds)
    return core


def extend_confirmed_seed(
    seed: MergedSeed,
    genes: list[pORF],
    taxonomy: dict[str, GeneTaxonomy],
    envelope: tuple[int, int],
    host_model: HostSignatureModel | None,
    host_threshold: float,
    *,
    ml_support: frozenset[str] = frozenset(),
    disable_host_stop: bool = False,
) -> RefinedBoundary:
    """Retain the core and commit only connected support inside the searched span."""
    seed = confirmed_seed_core(seed)
    available = host_model is not None and bool(host_model.token_weights)
    coverage_complete = True
    evidence: list[GeneEvidence] = []
    for gene in genes:
        if gene.start >= envelope[1] or gene.end <= envelope[0]:
            continue
        tax = taxonomy.get(gene.id)
        viral = tax is not None and has_identity_qualified_viral_hit(tax.top10_prefixes, tax.top10_pidents)
        host = bool(
            available
            and host_model is not None
            and tax is not None
            and host_lineage_share(tax, host_model) >= host_threshold
        )
        coverage_complete &= tax is not None
        evidence.append((gene, viral or gene.id in ml_support, host and not viral, tax is not None))
    left, left_reason, left_stop = _walk(seed.start, envelope[0], evidence, -1, disable_host_stop)
    right, right_reason, right_stop = _walk(seed.end, envelope[1], evidence, 1, disable_host_stop)
    return RefinedBoundary(
        scaffold=seed.scaffold,
        start=left,
        end=right,
        seed_id=seed.seed_id,
        original_start=seed.start,
        original_end=seed.end,
        confirmed_seeds=seed.confirmed_seeds,
        left_stop_reason=left_reason,
        right_stop_reason=right_reason,
        left_stop=left_stop,
        right_stop=right_stop,
        host_stop_available=available and coverage_complete,
        seed_sources=list(seed.sources),
        seed_confidence=seed.confidence,
        seed_hhg_score=seed.hhg_score,
        seed_novelty_score=seed.novelty_score,
        seed_compositional_score=seed.compositional_score,
        seed_has_mcp=seed.has_mcp,
        predicted_family=seed.predicted_family,
        region_classification_ncldv_markers=seed.region_classification_ncldv_markers,
        region_classification_vp_plv_markers=seed.region_classification_vp_plv_markers,
        region_classification_mirus_markers=seed.region_classification_mirus_markers,
    )


def _walk(
    core_edge: int, envelope_edge: int, evidence: list[GeneEvidence], direction: int, disable_host_stop: bool
) -> tuple[int, str, int | None]:
    """Classify all overlapping flank genes before accepting any long feature."""
    # Directed coordinates make left and right walks identical. Crossing genes
    # participate in the stop assessment but cannot remove any core sequence.
    edge = direction * core_edge
    cap = max(edge, direction * envelope_edge)
    ordered = sorted(
        (
            (
                min(direction * g.start, direction * g.end),
                max(direction * g.start, direction * g.end),
                support,
                host,
                covered,
            )
            for g, support, host, covered in evidence
            if max(direction * g.start, direction * g.end) > edge
        ),
        key=lambda row: (row[0], row[1]),
    )
    barrier = None
    run: list[int] = []
    if not disable_host_stop:
        for start, end, support, host, covered in ordered:
            run = [*run, start] if host else []
            if len(run) == 3:
                barrier = max(edge, run[0])
                break
    missing_edge = next((max(edge, start) for start, end, support, host, covered in ordered if not covered), None)
    limit = min(cap, barrier) if barrier is not None else cap
    if missing_edge is not None:
        limit = min(limit, missing_edge)
    accepted = edge
    pending = 0
    reason = "search_envelope"
    for start, end, support, host, covered in ordered:
        if start >= limit:
            break
        if support:
            accepted = max(accepted, min(end, limit))
            pending = 0
        else:
            pending += 1
            if pending >= 5:
                reason = "support_gap"
                break
    if barrier is not None and pending < 5:
        reason = "host_stop"
    if missing_edge is not None and (barrier is None or missing_edge < barrier) and pending < 5:
        reason = "taxonomy_unavailable"
    if accepted == edge and not ordered:
        reason = "no_exterior_evidence"
    return direction * accepted, reason, direction * barrier if barrier is not None else None


def merge_reachable_boundaries(boundaries: list[RefinedBoundary]) -> list[RefinedBoundary]:
    """Merge all compatible active overlaps, including across a nested separate seed."""
    finished: list[RefinedBoundary] = []
    active: list[RefinedBoundary] = []
    for boundary in sorted(boundaries, key=lambda item: (item.scaffold, item.start, item.end, item.seed_id)):
        finished.extend(item for item in active if item.scaffold != boundary.scaffold or item.end <= boundary.start)
        active = [item for item in active if item.scaffold == boundary.scaffold and item.end > boundary.start]
        index = 0
        while index < len(active):
            if _can_merge(active[index], boundary):
                boundary = _merge_pair(active.pop(index), boundary)
                index = 0
            else:
                index += 1
        active.append(boundary)
    return sorted(finished + active, key=lambda item: (item.scaffold, item.start, item.end, item.seed_id))


def _can_merge(previous: RefinedBoundary, boundary: RefinedBoundary) -> bool:
    union_start = min(previous.start, boundary.start)
    union_end = max(previous.end, boundary.end)
    blocked = any(
        (item.left_stop is not None and union_start < item.left_stop)
        or (item.right_stop is not None and union_end > item.right_stop)
        for item in (previous, boundary)
    )
    separate = previous.predicted_family == "CRESS" or boundary.predicted_family == "CRESS"
    separate |= ("frameshift_rescue" in previous.seed_sources) != ("frameshift_rescue" in boundary.seed_sources)
    return (
        previous.scaffold == boundary.scaffold
        and previous.start < boundary.end
        and boundary.start < previous.end
        and not blocked
        and not separate
    )


def _merge_pair(previous: RefinedBoundary, boundary: RefinedBoundary) -> RefinedBoundary:
    """Union reachable coordinates and retain every source constraint."""
    outer = boundary if boundary.end > previous.end else previous
    left_source = max((previous, boundary), key=lambda item: item.left_stop if item.left_stop is not None else -1)
    right_source = min(
        (outer, previous, boundary),
        key=lambda item: item.right_stop if item.right_stop is not None else float("inf"),
    )
    return replace(
        previous,
        start=min(previous.start, boundary.start),
        end=max(previous.end, boundary.end),
        confirmed_seeds=tuple(dict.fromkeys(previous.confirmed_seeds + boundary.confirmed_seeds)),
        seed_sources=sorted(set(previous.seed_sources + boundary.seed_sources)),
        original_start=min(previous.original_start, boundary.original_start),
        original_end=max(previous.original_end, boundary.original_end),
        seed_hhg_score=max(previous.seed_hhg_score, boundary.seed_hhg_score),
        seed_has_mcp=previous.seed_has_mcp or boundary.seed_has_mcp,
        host_stop_available=previous.host_stop_available and boundary.host_stop_available,
        left_stop=left_source.left_stop,
        left_stop_reason=left_source.left_stop_reason,
        right_stop=right_source.right_stop,
        right_stop_reason=right_source.right_stop_reason,
    )


def accept_seed_partition(parent: RefinedBoundary, proposed: list[RefinedBoundary]) -> list[RefinedBoundary] | None:
    """Accept only reachable children owning each complete core exactly once."""
    ordered = sorted(proposed, key=lambda child: (child.start, child.end))
    if any(left.end > right.start for left, right in zip(ordered, ordered[1:])):
        return None
    assigned: list[ConfirmedSeed] = []
    accepted = []
    for child in proposed:
        if child.scaffold != parent.scaffold or child.start < parent.start or child.end > parent.end:
            return None
        cores = tuple(core for core in parent.confirmed_seeds if child.start <= core.start and core.end <= child.end)
        if not cores or any(core in assigned for core in cores):
            return None
        assigned.extend(cores)
        accepted.append(
            replace(
                child,
                confirmed_seeds=cores,
                seed_id=min(core.seed_id for core in cores),
                tir_boundary_override=child.tir_boundary_override
                and (child.start, child.end) != (parent.start, parent.end),
                original_start=min(core.start for core in cores),
                original_end=max(core.end for core in cores),
                left_stop=parent.left_stop,
                right_stop=parent.right_stop,
                left_stop_reason=parent.left_stop_reason,
                right_stop_reason=parent.right_stop_reason,
                host_stop_available=parent.host_stop_available,
            )
        )
    return accepted if set(assigned) == set(parent.confirmed_seeds) else None


def refresh_boundary_classification(
    boundaries: list[RefinedBoundary],
    markers: Sequence[ValidatedMarkerHit],
    *,
    ablation_id: AblationID = AblationID.A0,
    single_marker_min_score: float = DEFAULT_A2_SINGLE_MARKER_MIN_SCORE,
) -> None:
    """Recompute family evidence using the same retained-marker roles as Phase3."""
    from virosync.pipeline.phase3.mcp_detection import is_mcp_gene

    by_scaffold: dict[str, list[ValidatedMarkerHit]] = {}
    for marker in markers:
        if decide_marker_hit_role(
            marker, ablation_id=ablation_id, single_marker_min_score=single_marker_min_score
        ).is_retained_evidence:
            by_scaffold.setdefault(marker.scaffold, []).append(marker)
    for boundary in boundaries:
        retained = [
            marker
            for marker in by_scaffold.get(boundary.scaffold, [])
            if marker.start < boundary.end and marker.end > boundary.start
        ]
        summary = classify_marker_hits(retained)
        # Compact CRESS classification has its own Phase1 identity-qualified
        # marker policy and is not assigned by the generic family classifier.
        if boundary.predicted_family != "CRESS":
            boundary.predicted_family = summary["classification"]
        boundary.region_classification_ncldv_markers = summary["ncldv_markers"]
        boundary.region_classification_vp_plv_markers = summary["vp_plv_markers"]
        boundary.region_classification_mirus_markers = summary["mirus_markers"]
        boundary.seed_has_mcp = any(is_mcp_gene(marker.hmm_target) for marker in retained)


def validate_source_membership(seeds: Sequence[MergedSeed], boundaries: Sequence[RefinedBoundary]) -> None:
    """Require exact original source membership across the Phase1/Phase2 boundary."""
    expected = Counter(core for seed in seeds for core in seed.confirmed_seeds)
    observed = Counter(core for boundary in boundaries for core in boundary.confirmed_seeds)
    if expected != observed:
        raise ValueError(
            "Phase2 source membership differs from Phase1: "
            f"missing={list((expected - observed).elements())!r}; "
            f"unexpected={list((observed - expected).elements())!r}"
        )
