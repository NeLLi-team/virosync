"""One coordinate authority for growth from protected seed-gene hulls."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, NamedTuple, TypeAlias

from virosync.ablation import AblationID
from virosync.pipeline.host_signatures import HostSignatureModel, host_lineage_share, host_stop_available
from virosync.pipeline.phase1.hhg_seeding import Anchor
from virosync.pipeline.phase1.marker_roles import DEFAULT_A2_SINGLE_MARKER_MIN_SCORE, decide_marker_hit_role
from virosync.pipeline.phase1.seed_merger import (
    ConfirmedSeed,
    MergedSeed,
    classify_marker_hits,
    validate_confirmed_seed_containment,
)
from virosync.pipeline.phase2.boundary_diamond import (
    GeneTaxonomy,
    _listify_top10_field,
    has_viral_hit_veto,
    has_viral_match,
    pORF,
)
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary

if TYPE_CHECKING:
    from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit


# Gene, positive growth support, host match, queried, independent viral/marker reset,
# viral match (validated viral marker or viral top-10 hit).
GeneEvidence: TypeAlias = tuple[pORF, bool, bool, bool, bool, bool]

# A run of three host genes is bridged when at least BRIDGE_MIN_VIRAL_GENES of the
# next BRIDGE_WINDOW_GENES walk entries have a viral match.
BRIDGE_WINDOW_GENES = 5
BRIDGE_MIN_VIRAL_GENES = 2

# A gene with hits and a host share below this value is foreign to the host.
FOREIGN_HOST_SHARE = 0.10


@dataclass(frozen=True)
class _EvidenceRules:
    """Per-seed inputs that decide what each gene adds to the walk."""

    host_model: HostSignatureModel | None
    host_stop_available: bool
    host_threshold: float
    ml_support: frozenset[str]
    viral_marker_genes: frozenset[str]


class _Entry(NamedTuple):
    """One gene of a directed walk; coordinates grow away from the seed."""

    start: int
    end: int
    support: bool
    host: bool
    covered: bool
    reset_host_run: bool
    viral_match: bool


class _HostRun(NamedTuple):
    """A run of three host genes in walk order."""

    first: _Entry
    third_index: int
    bridged: bool


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
    viral_marker_genes: frozenset[str] = frozenset(),
    disable_host_stop: bool = False,
    protected_seed_intervals: tuple[tuple[int, int], ...] = (),
    walk_limit: tuple[int, int] | None = None,
) -> RefinedBoundary:
    """Retain the core and commit only connected support inside the walk limit.

    Three host-signature matches without intervening qualifying viral support
    stop growth, unless at least two of the next five walk entries have a viral match.
    A bridge changes only where the walk may go. The run that limits the walk always
    sets a stop. A bridged run sets one only if the walk did not cross it, that is,
    if the walk accepted no supporting gene that follows the run's third host gene.
    Neutral genes do not reset the count. A viral match (a validated viral marker
    or a viral top-10 hit) supports growth at any host share. Marker-set
    membership supports it below the host-share threshold. Deep-only hits veto
    host status without extending the hull.

    The walk limit defaults to ``envelope``. Genes outside the envelope but inside
    the limit support growth only through a viral match, and a searched gene
    without hits there counts toward the support gap.

    A supporting gene is accepted when it ends before the limit that stops the
    walk. It is also accepted when it starts before a host stop or before the first
    gene without a record and ends inside that gene, even if that takes the region
    past the walk limit. The stop then moves to the end of the accepted gene, so it
    never lies inside the returned region. Any other gene that crosses the limit is
    not accepted. Only supporters before the first gene without a record can bridge
    a run of three host genes. Disabling the host stop (A4) removes the stop and its
    limit but keeps the limit and allowance of a gene without a record. It also accepts
    every gene that the walk with the stop accepts.
    """
    # Repeat assessment keeps the search context without protecting its padding.
    original_start, original_end = seed.start, seed.end
    seed = confirmed_seed_core(seed)
    available = host_stop_available(host_model)
    rules = _EvidenceRules(host_model, available, host_threshold, ml_support, viral_marker_genes)
    limit = envelope if walk_limit is None else walk_limit
    coverage_complete = True
    evidence: list[GeneEvidence] = []
    for gene in genes:
        if not _overlaps(gene, limit):
            continue
        tax = taxonomy.get(gene.id)
        inside_envelope = _overlaps(gene, envelope)
        # A gene outside the envelope never changes whether host stopping was assessable.
        coverage_complete &= tax is not None or not inside_envelope
        item = _gene_evidence(gene, tax, rules, inside_envelope=inside_envelope)
        if item is not None:
            evidence.append(item)
    left, left_reason, left_stop = _walk(seed.start, limit[0], evidence, -1, disable_host_stop)
    right, right_reason, right_stop = _walk(seed.end, limit[1], evidence, 1, disable_host_stop)
    return RefinedBoundary(
        scaffold=seed.scaffold,
        start=left,
        end=right,
        seed_id=seed.seed_id,
        original_start=original_start,
        original_end=original_end,
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


def seed_walk_limit(
    envelope: tuple[int, int], genes: Sequence[pORF], taxonomy: Mapping[str, GeneTaxonomy]
) -> tuple[int, int]:
    """Widen a seed's search envelope over consecutive genes that have a taxonomy record.

    Each edge moves outward over the genes beyond it, in scaffold order, and stops
    before the first gene without a record. A gene that crosses the edge moves it to
    its far end. Other seeds on the scaffold are searched with their own envelopes,
    so such records often exist. No search runs here.
    """
    return _outward_limit(envelope[0], genes, taxonomy, -1), _outward_limit(envelope[1], genes, taxonomy, 1)


def _outward_limit(edge: int, genes: Sequence[pORF], taxonomy: Mapping[str, GeneTaxonomy], direction: int) -> int:
    """Return the edge moved outward over consecutive genes with a record; direction is 1 or -1."""
    directed_edge = direction * edge
    outer = sorted(
        (min(direction * gene.start, direction * gene.end), max(direction * gene.start, direction * gene.end), gene.id)
        for gene in genes
        if max(direction * gene.start, direction * gene.end) > directed_edge
    )
    limit = directed_edge
    for _start, end, gene_id in outer:
        if gene_id not in taxonomy:
            break
        limit = max(limit, end)
    return direction * limit


def _overlaps(gene: pORF, interval: tuple[int, int]) -> bool:
    return gene.start < interval[1] and gene.end > interval[0]


def _gene_evidence(
    gene: pORF, tax: GeneTaxonomy | None, rules: _EvidenceRules, *, inside_envelope: bool
) -> GeneEvidence | None:
    """Return what one gene adds to the walk, or None for a searched gene to skip.

    Inside the seed's envelope a searched gene without hits or independent evidence
    is skipped. Outside it, only a viral match supports growth, and such a gene stays
    in the walk as an entry without support. Host status and the host-run reset are
    the same in both places.
    """
    marker = gene.id in rules.viral_marker_genes or (tax is not None and tax.phase1_validated_viral_marker)
    share = (
        host_lineage_share(tax, rules.host_model)
        if rules.host_stop_available and rules.host_model is not None and tax is not None
        else None
    )
    host_share = share is not None and share >= rules.host_threshold
    has_hits = tax is not None and any(str(p).strip() for p in _listify_top10_field(tax.top10_prefixes))
    reset_host_run = marker or (tax is not None and has_viral_hit_veto(tax))
    if inside_envelope and tax is not None and not has_hits and not reset_host_run and gene.id not in rules.ml_support:
        return None
    viral_match = tax is not None and has_viral_match(tax)
    support = viral_match
    if inside_envelope:
        foreign = share is not None and share < FOREIGN_HOST_SHARE and has_hits
        member = gene.id in rules.viral_marker_genes and not host_share
        support = support or foreign or member or gene.id in rules.ml_support
    host = tax is not None and is_extension_host_gene(tax, host_share_ok=host_share, marker=marker)
    return (gene, support, host, tax is not None, reset_host_run, viral_match)


def is_extension_host_gene(taxonomy: GeneTaxonomy, *, host_share_ok: bool, marker: bool) -> bool:
    """Count a host-signature match only when no qualifying viral evidence vetoes it."""
    return host_share_ok and not marker and not has_viral_hit_veto(taxonomy)


def _walk(
    core_edge: int, limit_edge: int, evidence: list[GeneEvidence], direction: int, disable_host_stop: bool
) -> tuple[int, str, int | None]:
    """Classify all overlapping flank genes before accepting any long feature.

    Returns the accepted end, the stop reason and the stop coordinate. The first unbridged
    run of three host genes limits the walk and always sets a stop. A bridged run sets a stop
    only if the walk did not cross it: it is crossed when the walk accepted a supporting entry
    that follows the run's third host gene in walk order. Walk order is directed, so the test
    holds for both sides. The first such run sets the stop, at the run's first host start or
    at the accepted end if that lies further out, so a stop never lies inside the region.
    """
    # Directed coordinates make left and right walks identical. A supporting gene
    # that crosses the limit is accepted only when it ends inside the gene that sets
    # a host stop or a missing-taxonomy edge; it never removes core sequence.
    # Disabling the host stop (A4) drops the stop but keeps that allowance, so A4
    # never accepts less than the walk with the stop.
    edge = direction * core_edge
    cap = max(edge, direction * limit_edge)
    ordered = sorted(
        (
            _Entry(
                min(direction * g.start, direction * g.end),
                max(direction * g.start, direction * g.end),
                support,
                host,
                covered,
                reset_host_run,
                viral_match,
            )
            for g, support, host, covered, reset_host_run, viral_match in evidence
            if max(direction * g.start, direction * g.end) > edge
        ),
        key=lambda entry: (entry.start, entry.end),
    )
    runs = _host_runs(ordered)
    first_host = runs[-1].first if runs and not runs[-1].bridged else None
    first_uncovered = next((entry for entry in ordered if not entry.covered), None)
    barrier = None if first_host is None else max(edge, first_host.start)
    missing_edge = None if first_uncovered is None else max(edge, first_uncovered.start)
    # A supporting entry may end at `reach`: the limit, or the end of the gene that sets it.
    # Host and missing-record ends both constrain acceptance, including nested genes.
    host_setters = [] if first_host is None else [(max(edge, first_host.start), first_host.end)]
    gap_setters = [] if first_uncovered is None else [(max(edge, first_uncovered.start), first_uncovered.end)]
    limit, reach = min([(cap, cap), *host_setters, *gap_setters])
    if first_host is not None:
        reach = min(reach, first_host.end)
    if first_uncovered is not None:
        reach = min(reach, first_uncovered.end)
    # The walk scans to `scan_limit` and accepts up to `scan_reach`. Ablation A4 drops
    # the host setter there, but it still accepts what the walk with the stop accepts.
    scan_limit, scan_reach = min([(cap, cap), *gap_setters]) if disable_host_stop else (limit, reach)
    if disable_host_stop:
        barrier = None
    accepted = edge
    last_accepted = -1
    pending = 0
    reason = "search_envelope"
    for index, entry in enumerate(ordered):
        if entry.start >= scan_limit:
            break
        if entry.support:
            if entry.end <= scan_reach or (entry.start < limit and entry.end <= reach):
                accepted = max(accepted, entry.end)
                last_accepted = index
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
    # The unbridged run always sets a stop. A bridged run sets one only if the walk did
    # not cross it. A stop must not lie inside the region, so it moves out to the accepted end.
    uncrossed = next((run for run in runs if not run.bridged or last_accepted <= run.third_index), None)
    stop = None if disable_host_stop or uncrossed is None else max(edge, uncrossed.first.start, accepted)
    return direction * accepted, reason, direction * stop if stop is not None else None


def _host_runs(ordered: list[_Entry]) -> list[_HostRun]:
    """Return the runs of three host genes in walk order, up to the first unbridged run.

    A run is bridged when at least ``BRIDGE_MIN_VIRAL_GENES`` of the next
    ``BRIDGE_WINDOW_GENES`` entries have a viral match. Foreign, marker-set and ML
    support do not count. Deep-hit vetoes reset the host counter but do not bridge it. Only entries before the first gene
    without a taxonomy record count, because the walk cannot reach later ones.
    """
    runs: list[_HostRun] = []
    run: list[_Entry] = []
    first_missing = next((index for index, entry in enumerate(ordered) if not entry.covered), len(ordered))
    for index, entry in enumerate(ordered):
        if entry.host:
            run = [*run, entry]
        elif entry.reset_host_run:
            run = []
        if len(run) < 3:
            continue
        window = ordered[index + 1 : min(index + 1 + BRIDGE_WINDOW_GENES, first_missing)]
        bridged = sum(follower.viral_match for follower in window) >= BRIDGE_MIN_VIRAL_GENES
        runs.append(_HostRun(run[0], index, bridged))
        if not bridged:
            break
        run = []
    return runs


def merge_reachable_boundaries(boundaries: list[RefinedBoundary], compatible=None) -> list[RefinedBoundary]:
    """Merge all compatible active overlaps, including across a nested separate seed."""
    finished: list[RefinedBoundary] = []
    active: list[RefinedBoundary] = []
    for boundary in sorted(boundaries, key=lambda item: (item.scaffold, item.start, item.end, item.seed_id)):
        finished.extend(item for item in active if item.scaffold != boundary.scaffold or item.end <= boundary.start)
        active = [item for item in active if item.scaffold == boundary.scaffold and item.end > boundary.start]
        index = 0
        while index < len(active):
            if _can_merge(active[index], boundary) and (compatible is None or compatible(active[index], boundary)):
                boundary = _merge_pair(active.pop(index), boundary)
                index = 0
            else:
                if compatible is not None and _can_merge(active[index], boundary):
                    from virosync.pipeline.eve_type_rules import separate_overlap

                    active[index], boundary = separate_overlap(active[index], boundary)
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


def accept_seed_partition(
    parent: RefinedBoundary, proposed: list[RefinedBoundary], *, allowed_interval: tuple[int, int]
) -> list[RefinedBoundary] | None:
    """Accept only reachable children owning each complete core exactly once."""
    if not allowed_interval[0] <= parent.start < parent.end <= allowed_interval[1]:
        raise ValueError("allowed interval must contain parent")
    ordered = sorted(proposed, key=lambda child: (child.start, child.end))
    if any(left.end > right.start for left, right in zip(ordered, ordered[1:])):
        return None
    assigned: list[ConfirmedSeed] = []
    accepted = []
    for child in proposed:
        if child.scaffold != parent.scaffold or child.start < allowed_interval[0] or child.end > allowed_interval[1]:
            return None
        if (parent.left_stop is not None and child.start < parent.left_stop) or (
            parent.right_stop is not None and child.end > parent.right_stop
        ):
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
    metadata_families: dict,
) -> None:
    """Recompute ordinary family evidence using the configured marker roles."""
    from virosync.pipeline.phase3.mcp_detection import is_mcp_gene

    by_scaffold: dict[str, list[Anchor]] = {}
    for marker in markers:
        role = decide_marker_hit_role(marker, ablation_id=ablation_id, single_marker_min_score=single_marker_min_score)
        if not (role.is_production_validated or role.is_tier1_bypassed):
            continue
        by_scaffold.setdefault(marker.scaffold, []).append(
            Anchor(
                porf_id=marker.query_porf,
                scaffold=marker.scaffold,
                start=marker.start,
                end=marker.end,
                strand=getattr(marker, "strand", ""),
                hallmark_gene=marker.hmm_target,
                score=marker.hmm_score,
                evalue=marker.hmm_evalue,
                validation_status=marker.validation_status,
                pfam_override=getattr(marker, "pfam_override", ""),
                tier1_bypassed=role.is_tier1_bypassed,
                top10_prefixes=getattr(marker, "top10_prefixes", ""),
                top10_pidents=getattr(marker, "top10_pidents", ""),
                top10_targets=getattr(marker, "top10_targets", ""),
            )
        )
    for boundary in boundaries:
        retained = [
            marker
            for marker in by_scaffold.get(boundary.scaffold, [])
            if marker.start < boundary.end and marker.end > boundary.start
        ]
        summary = classify_marker_hits(retained, metadata_families=metadata_families)
        # Compact CRESS classification has its own Phase1 identity-qualified
        # marker policy and is not assigned by the generic family classifier.
        if boundary.predicted_family != "CRESS":
            boundary.predicted_family = summary["classification"]
        boundary.region_classification_ncldv_markers = summary["ncldv_markers"]
        boundary.region_classification_vp_plv_markers = summary["vp_plv_markers"]
        boundary.region_classification_mirus_markers = summary["mirus_markers"]
        boundary.seed_has_mcp = any(is_mcp_gene(marker.hallmark_gene) for marker in retained)


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
