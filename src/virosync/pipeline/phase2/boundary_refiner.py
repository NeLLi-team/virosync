"""
Boundary Refiner for EVE Detection.

This module contains the older two-tier CRF boundary-refinement implementation:
1. Tier 1: binary screening on 1kb windows to identify a region of interest
2. Tier 2: anatomical refinement on 250bp windows for precise boundaries

The active workflow now uses gene extension, batched Diamond taxonomy, and host
taxonomy trimming. The CRF refiner remains for compatibility with older
experiments and artifacts.
"""

import logging
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from virosync.ablation import AblationID, InterventionCounts
from virosync.pipeline.host_signatures import (
    HostSignatureModel,
    score_host_signature_record,
)
from virosync.pipeline.phase1.seed_merger import MergedSeed
from virosync.pipeline.phase2.boundary_diamond import (
    MIN_VIRAL_HIT_PIDENT,
    SeedGeneMapping,
    VIRAL_PREFIXES,
    has_identity_qualified_viral_hit,
    pORF,
)
from virosync.pipeline.taxonomy_utils import calculate_fingerprint_overlap, compute_hit_weight


logger = logging.getLogger(__name__)


def extend_seeds_by_genes(
    seeds: list[MergedSeed],
    proteome_index: dict[str, list[pORF]],
    extension_genes: int = 5,
) -> list[MergedSeed]:
    """Extend each seed by ±N genes, then merge overlapping seeds on same scaffold.

    For each seed:
    - Find all pORFs on the scaffold from proteome_index
    - Identify genes overlapping the seed region
    - Extend to include ±extension_genes on each side
    - Update seed.start/end to the outermost extended gene coordinates

    Then merge overlapping seeds on the same scaffold (sorted by start).

    Args:
        seeds: List of MergedSeed from Phase 1
        proteome_index: scaffold -> sorted list of pORF objects
        extension_genes: Number of genes to extend on each side (default 5)

    Returns:
        List of (possibly merged) extended MergedSeed objects
    """
    from dataclasses import replace
    from collections import defaultdict

    if not seeds:
        return []

    # Group seeds by scaffold
    scaffold_seeds: dict[str, list[MergedSeed]] = defaultdict(list)
    for seed in seeds:
        scaffold_seeds[seed.scaffold].append(seed)

    extended_all: list[MergedSeed] = []

    for scaffold, scaffold_seed_list in scaffold_seeds.items():
        scaffold_porfs = proteome_index.get(scaffold, [])
        if not scaffold_porfs:
            # No genes on scaffold — keep original seeds unchanged
            extended_all.extend(scaffold_seed_list)
            continue

        # Extend each seed by ±N genes
        extended_scaffold: list[MergedSeed] = []
        for seed in scaffold_seed_list:
            if seed.predicted_family == "CRESS":
                extended_scaffold.append(seed)
                continue

            # Find genes overlapping the seed (midpoint within seed, or any overlap)
            overlapping_indices = []
            for i, p in enumerate(scaffold_porfs):
                if p.start < seed.end and p.end > seed.start:
                    overlapping_indices.append(i)

            if not overlapping_indices:
                # No overlapping genes — keep original seed as-is
                extended_scaffold.append(seed)
                logger.debug(
                    "Seed %s has no overlapping genes on %s; keeping original bounds",
                    seed.seed_id, scaffold,
                )
                continue

            first_idx = min(overlapping_indices)
            last_idx = max(overlapping_indices)

            # Extend by ±N genes, clamped to scaffold bounds
            ext_first = max(0, first_idx - extension_genes)
            ext_last = min(
                len(scaffold_porfs) - 1,
                last_idx + extension_genes,
            )

            new_start = scaffold_porfs[ext_first].start
            new_end = scaffold_porfs[ext_last].end

            extended_scaffold.append(replace(seed, start=new_start, end=new_end))

        # Sort by start position for merging
        extended_scaffold.sort(key=lambda s: (s.start, s.end))

        # Merge overlapping extended seeds
        merged_scaffold: list[MergedSeed] = []
        current = extended_scaffold[0]

        for next_seed in extended_scaffold[1:]:
            cress_in_pair = (
                current.predicted_family == "CRESS"
                or next_seed.predicted_family == "CRESS"
            )
            mixed_rescue_pair = (
                "frameshift_rescue" in current.sources
            ) != (
                "frameshift_rescue" in next_seed.sources
            )
            if (
                next_seed.start < current.end
                and not cress_in_pair
                and not mixed_rescue_pair
            ):
                # Overlapping — merge
                combined_sources = sorted(set(current.sources) | set(next_seed.sources))
                combined_anchors = current.anchors + [
                    a for a in next_seed.anchors if a not in current.anchors
                ]
                combined_hhg = current.hhg_anchors + [
                    a for a in next_seed.hhg_anchors if a not in current.hhg_anchors
                ]
                combined_ncldv = (
                    current.region_classification_ncldv_markers
                    + next_seed.region_classification_ncldv_markers
                )
                combined_vp_plv = (
                    current.region_classification_vp_plv_markers
                    + next_seed.region_classification_vp_plv_markers
                )
                combined_mirus = (
                    current.region_classification_mirus_markers
                    + next_seed.region_classification_mirus_markers
                )
                logger.info(
                    "Merging overlapping extended seeds on %s: [%d-%d] + [%d-%d] -> [%d-%d]",
                    scaffold, current.start, current.end,
                    next_seed.start, next_seed.end,
                    current.start, max(current.end, next_seed.end),
                )
                _conf_rank = {"high": 3, "medium": 2, "low": 1}
                best_conf = max(
                    current.confidence, next_seed.confidence,
                    key=lambda c: _conf_rank.get(c, 0),
                )
                current = replace(
                    current,
                    end=max(current.end, next_seed.end),
                    sources=combined_sources,
                    confidence=best_conf,
                    hhg_score=max(current.hhg_score, next_seed.hhg_score),
                    anchors=combined_anchors,
                    hhg_anchors=combined_hhg,
                    region_classification_ncldv_markers=combined_ncldv,
                    region_classification_vp_plv_markers=combined_vp_plv,
                    region_classification_mirus_markers=combined_mirus,
                )
            else:
                merged_scaffold.append(current)
                current = next_seed

        merged_scaffold.append(current)
        extended_all.extend(merged_scaffold)

    # Re-assign seed_ids after extension/merge (indices may have changed)
    for idx, seed in enumerate(extended_all):
        seed.seed_id = f"seed_{idx}_{seed.scaffold}_{seed.start}"

    logger.info(
        "Gene extension: %d seeds -> %d extended/merged seeds (±%d genes)",
        len(seeds), len(extended_all), extension_genes,
    )

    return extended_all


@dataclass
class RefinedBoundary:
    """
    Refined EVE boundary from Phase 2 boundary processing.

    Contains precise coordinates and confidence metrics.
    """

    scaffold: str
    start: int
    end: int

    # Original seed info
    seed_id: str = ""  # Stable ID from MergedSeed for boundary-to-seed mapping
    original_start: int = 0
    original_end: int = 0
    candidate_start: Optional[int] = None
    candidate_end: Optional[int] = None
    host_trim_reason: str = ""
    host_trim_common_euk_taxonomy: str = ""

    # Validated-marker floor (Phase-3 re-admit METADATA -- never applied to
    # start/end). When this boundary's own seed span carried >=2 validated viral
    # proteins, or one confirmed frameshift-rescued protein on a rescue-derived
    # boundary, this records their min..max coordinates as a recovery-only
    # alternative span. The Phase-3 re-admit pass (which fires only on REJECTED
    # boundaries) may synthesize and gate this floored alternative without ever
    # mutating an accepted boundary.
    marker_floor_start: Optional[int] = None
    marker_floor_end: Optional[int] = None

    # Seed evidence metadata (passed from Phase 1)
    seed_sources: list[str] = field(default_factory=list)  # ["hhg", "novelty", "compositional"]
    seed_confidence: str = "low"  # "high", "medium", "low"
    seed_hhg_score: float = 0.0
    seed_novelty_score: float = 0.0
    seed_compositional_score: float = 0.0
    seed_has_mcp: bool = False
    # Region classification from Phase 1 (NCLDV, VP, PLV, MIRUS, MIXED, UNKNOWN)
    predicted_family: str = ""
    region_classification_ncldv_markers: int = 0
    region_classification_vp_plv_markers: int = 0
    region_classification_mirus_markers: int = 0

    # Boundary confidence
    confidence: float = 0.0
    posterior_probability: float = 0.0

    # Anatomical breakdown
    core_viral_start: Optional[int] = None
    core_viral_end: Optional[int] = None
    flank_5_start: Optional[int] = None
    flank_5_end: Optional[int] = None
    flank_3_start: Optional[int] = None
    flank_3_end: Optional[int] = None

    # CRF state sequence
    state_sequence: list[int] = field(default_factory=list)
    state_posteriors: np.ndarray = field(default=None)

    # Evidence
    hallmark_genes: list[str] = field(default_factory=list)
    max_kfd: float = 0.0
    gc_deviation: float = 0.0
    cub_deviation: float = 0.0
    mean_novelty: float = 0.0

    # Window features for Phase 3 coherence analysis
    # Stored as list of WindowFeatures from Tier 2 (or Tier 1 if skip_tier2)
    window_features: list = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start

    def to_bed_line(self) -> str:
        """Format as BED line."""
        score = int(min(1000, self.confidence * 1000))
        return f"{self.scaffold}\t{self.start}\t{self.end}\tEVE_{self.scaffold}_{self.start}\t{score}\t."

    def to_gff_line(self) -> str:
        """Format as GFF3 line with attributes."""
        score = int(min(1000, self.confidence * 1000))
        attrs = [
            f"ID=EVE_{self.scaffold}_{self.start}",
            f"confidence={self.confidence:.3f}",
            f"posterior={self.posterior_probability:.3f}",
        ]
        if self.core_viral_start is not None:
            # Same 0-based half-open to 1-based inclusive shift as the
            # coordinate columns below, so the attributes match their record.
            attrs.append(f"core_start={self.core_viral_start + 1}")
            attrs.append(f"core_end={self.core_viral_end}")
        if self.hallmark_genes:
            attrs.append(f"hallmarks={','.join(self.hallmark_genes)}")

        return f"{self.scaffold}\tViroSync\tEVE\t{self.start+1}\t{self.end}\t{score}\t.\t.\t{';'.join(attrs)}"


def constrain_to_seed_bounds(
    boundary: RefinedBoundary,
    seed_mapping: SeedGeneMapping,
) -> RefinedBoundary:
    """
    Constrain refined boundary to not extend beyond ±N genes from original seed.

    This enforces the hard cap on boundary extension to ensure all genes
    within the final boundary have Diamond taxonomy data from Phase 2b.

    Args:
        boundary: Refined boundary (may extend beyond seed)
        seed_mapping: SeedGeneMapping with flanking bounds

    Returns:
        New RefinedBoundary constrained to seed flanking bounds
    """
    constrained_start = max(boundary.start, seed_mapping.flank_start_bp)
    constrained_end = min(boundary.end, seed_mapping.flank_end_bp)

    # Check if constraint was applied
    start_constrained = constrained_start > boundary.start
    end_constrained = constrained_end < boundary.end

    if start_constrained or end_constrained:
        # Log the constraint application
        constraint_msg = []
        if start_constrained:
            constraint_msg.append(
                f"start {boundary.start} -> {constrained_start} "
                f"(capped at -{seed_mapping.flank_genes_config} genes)"
            )
        if end_constrained:
            constraint_msg.append(
                f"end {boundary.end} -> {constrained_end} "
                f"(capped at +{seed_mapping.flank_genes_config} genes)"
            )

        logger.info(
            "Boundary constraint applied to %s: %s",
            seed_mapping.seed_id,
            "; ".join(constraint_msg),
        )

    # If constraint would make region invalid, use original seed bounds
    if constrained_start >= constrained_end:
        logger.warning(
            "Boundary constraint would create invalid region for %s: "
            "[%d, %d] -> [%d, %d]. Using original seed bounds [%d, %d]",
            seed_mapping.seed_id,
            boundary.start,
            boundary.end,
            constrained_start,
            constrained_end,
            seed_mapping.seed_start,
            seed_mapping.seed_end,
        )
        constrained_start = seed_mapping.seed_start
        constrained_end = seed_mapping.seed_end

    # Create new RefinedBoundary with constrained coordinates
    # Use dataclasses.replace to preserve nested objects (like WindowFeatures)
    # instead of asdict which deep-converts them to dicts
    from dataclasses import replace

    return replace(boundary, start=constrained_start, end=constrained_end)


def annotate_boundaries_with_marker_floor(
    boundaries: list,
    validated_markers: list,
) -> int:
    """Annotate each boundary with its validated-marker floor span (NO mutation).

    Phase-2 host-aware trimming (``host_signature_trim`` window scan +
    ``trim_boundary_by_host_taxonomy``) can collapse a marker-dense NCLDV/MIRUS
    seed below its validated-marker span, stripping the hallmark genes out of the
    final boundary so the v2 NCLDV/MIRUS output gate (``length > 5000 OR
    has_mcp``) then rejects a genuine EVE.

    The earlier fix MUTATED the boundary here to contain the marker span. That
    regressed recall: extending an ALREADY-ACCEPTED region's boundary pulls host
    genes in, lowers Phase-3 confidence, drops the tier MEDIUM->LOW, and the
    stricter LOW NCLDV gate then rejects a region MEDIUM had accepted on length
    alone -- so NCLDV genes were LOST on rhizophagus/tstriata. Instead this function only COMPUTES and STORES
    the floored-alternative span as metadata (``marker_floor_start`` /
    ``marker_floor_end``); ``boundary.start`` / ``boundary.end`` are left exactly
    as Phase-2 trimming produced them. The Phase-3 re-admit pass consumes this
    metadata on REJECTED boundaries only, so an accepted boundary is never
    modified and the accepted set structurally cannot regress.

    Specificity guard: ordinary markers still require two distinct proteins.
    One marker suffices only when both the seed provenance and generated protein
    ID identify a confirmed frameshift rescue. Markers are scoped to the
    boundary's own seed span
    (``original_start`` / ``original_end``) so the floor never pulls in validated
    markers from a different region on the same scaffold, and only validated
    markers (``validation_status`` in ``{validated, validated_novel}``) count --
    host genes are always "unvalidated" (no viral Diamond support), so host
    regions / spurious HMM-only hits get no floor.

    Returns the number of boundaries annotated with a floor strictly wider than
    their current span.
    """
    if not boundaries or not validated_markers:
        return 0

    n_annotated = 0
    from virosync.pipeline.phase1.frameshift_screening import (
        is_rescued_protein_id,
    )

    for boundary in boundaries:
        span_markers = [
            m for m in validated_markers
            if m.scaffold == boundary.scaffold
            and m.validation_status in ("validated", "validated_novel")
            and boundary.original_start <= (m.start + m.end) // 2 <= boundary.original_end
        ]
        # Count distinct marker-bearing proteins. One ordinary pORF can produce
        # several HMM hits, so raw hit count must not satisfy the floor.
        rescue_boundary = "frameshift_rescue" in (
            getattr(boundary, "seed_sources", []) or []
        )
        has_rescue_marker = any(
            is_rescued_protein_id(m.query_porf) for m in span_markers
        )
        if (
            len({m.query_porf for m in span_markers}) < 2
            and not (rescue_boundary and has_rescue_marker)
        ):
            continue
        marker_lo = min(m.start for m in span_markers)
        marker_hi = max(m.end for m in span_markers)
        floor_start = min(boundary.start, marker_lo)
        floor_end = max(boundary.end, marker_hi)
        # Only record a floor that would actually widen the boundary; an equal or
        # narrower span gives the Phase-3 re-admit pass nothing to recover.
        if floor_start < boundary.start or floor_end > boundary.end:
            boundary.marker_floor_start = floor_start
            boundary.marker_floor_end = floor_end
            logger.info(
                "Validated-marker floor recorded: %s %d-%d (current) -> floor %d-%d "
                "(%d validated markers, span %d-%d)",
                boundary.seed_id or boundary.scaffold,
                boundary.start, boundary.end, floor_start, floor_end,
                len(span_markers), marker_lo, marker_hi,
            )
            n_annotated += 1
    return n_annotated


def merge_adjacent_viral_boundaries(
    boundaries: list[RefinedBoundary],
    taxonomy_map: dict,
    max_gap_bp: int = 10000,
    min_viral_fraction: float = 0.3,
) -> list[RefinedBoundary]:
    """
    Merge adjacent EVE boundaries if the gap between them contains viral genes.

    After boundary constraints are applied, two EVEs might be separated by
    a small gap that contains viral genes (from the flanking regions).
    If the gap shows viral taxonomy, the EVEs should be merged.

    Args:
        boundaries: List of RefinedBoundary objects (constrained)
        taxonomy_map: Dict mapping porf_id -> GeneTaxonomy from Phase 2b
        max_gap_bp: Maximum gap between boundaries to consider merging (bp)
        min_viral_fraction: Minimum fraction of gap genes that must be viral

    Returns:
        List of merged RefinedBoundary objects
    """
    if not boundaries or len(boundaries) < 2:
        return boundaries

    from collections import defaultdict

    # Group by scaffold
    scaffold_boundaries: dict[str, list[RefinedBoundary]] = defaultdict(list)
    for b in boundaries:
        scaffold_boundaries[b.scaffold].append(b)

    merged_all = []
    total_merges = 0

    for scaffold, scaffold_list in scaffold_boundaries.items():
        # Sort by start position
        scaffold_list.sort(key=lambda b: b.start)

        merged_scaffold = []
        current = scaffold_list[0]

        for next_boundary in scaffold_list[1:]:
            gap_start = current.end
            gap_end = next_boundary.start
            cress_in_pair = (
                getattr(current, "predicted_family", "") == "CRESS"
                or getattr(next_boundary, "predicted_family", "") == "CRESS"
            )
            current_sources = getattr(current, "seed_sources", []) or []
            next_sources = getattr(next_boundary, "seed_sources", []) or []
            mixed_rescue_pair = (
                "frameshift_rescue" in current_sources
            ) != (
                "frameshift_rescue" in next_sources
            )

            # Strictly overlapping same-provenance boundaries merge. Touching
            # half-open intervals continue through the evidence-aware gap path.
            if (
                gap_end < gap_start
                and not cress_in_pair
                and not mixed_rescue_pair
            ):
                from dataclasses import replace

                combined_ncldv = (
                    getattr(current, "region_classification_ncldv_markers", 0)
                    + getattr(next_boundary, "region_classification_ncldv_markers", 0)
                )
                combined_vp_plv = (
                    getattr(current, "region_classification_vp_plv_markers", 0)
                    + getattr(next_boundary, "region_classification_vp_plv_markers", 0)
                )
                combined_mirus = (
                    getattr(current, "region_classification_mirus_markers", 0)
                    + getattr(next_boundary, "region_classification_mirus_markers", 0)
                )
                combined_sources = sorted(set(current_sources) | set(next_sources))

                logger.info(
                    "Merging overlapping EVEs on %s: [%d-%d] + [%d-%d] "
                    "(overlap=%dbp)",
                    scaffold,
                    current.start,
                    current.end,
                    next_boundary.start,
                    next_boundary.end,
                    gap_start - gap_end,
                )
                current = replace(
                    current,
                    end=max(current.end, next_boundary.end),
                    original_end=max(
                        current.original_end, next_boundary.original_end
                    ),
                    confidence=max(current.confidence, next_boundary.confidence),
                    region_classification_ncldv_markers=combined_ncldv,
                    region_classification_vp_plv_markers=combined_vp_plv,
                    region_classification_mirus_markers=combined_mirus,
                    seed_sources=combined_sources,
                )
                total_merges += 1
                continue

            # Check if gap is within merge distance
            if (
                not cress_in_pair
                and not mixed_rescue_pair
                and gap_end - gap_start <= max_gap_bp
            ):
                # Find genes in the gap and check their taxonomy
                gap_genes = []
                for porf_id, tax in taxonomy_map.items():
                    if (
                        hasattr(tax, "scaffold")
                        and tax.scaffold == scaffold
                        and hasattr(tax, "start")
                        and hasattr(tax, "end")
                        and tax.start >= gap_start
                        and tax.end <= gap_end
                    ):
                        gap_genes.append(tax)

                should_merge = False
                merge_reason = ""

                # Calculate viral fraction in gap
                if gap_genes:
                    n_viral = sum(
                        1
                        for g in gap_genes
                        if getattr(g, "has_ncldv_mirus", False)
                        or getattr(g, "has_vp_plv", False)
                    )
                    viral_fraction = n_viral / len(gap_genes)
                    if viral_fraction >= min_viral_fraction:
                        should_merge = True
                        merge_reason = (
                            f"gap=%dbp, %d/%d genes viral"
                            % (gap_end - gap_start, n_viral, len(gap_genes))
                        )

                # Flanking viral context: if both sides of the gap have
                # viral genes within the nearest N genes, merge even when
                # the gap itself contains no-hit or host genes.
                if not should_merge:
                    _flank_n = 5
                    _is_viral = lambda g: (
                        getattr(g, "has_ncldv_mirus", False)
                        or getattr(g, "has_vp_plv", False)
                    )
                    upstream_flank = sorted(
                        (
                            tax for tax in taxonomy_map.values()
                            if getattr(tax, "scaffold", None) == scaffold
                            and hasattr(tax, "end") and tax.end <= gap_start
                            and hasattr(tax, "start") and tax.start >= current.start
                        ),
                        key=lambda g: g.start,
                        reverse=True,
                    )[:_flank_n]
                    downstream_flank = sorted(
                        (
                            tax for tax in taxonomy_map.values()
                            if getattr(tax, "scaffold", None) == scaffold
                            and hasattr(tax, "start") and tax.start >= gap_end
                            and hasattr(tax, "end") and tax.end <= next_boundary.end
                        ),
                        key=lambda g: g.start,
                    )[:_flank_n]
                    if (
                        upstream_flank
                        and downstream_flank
                        and any(_is_viral(g) for g in upstream_flank)
                        and any(_is_viral(g) for g in downstream_flank)
                    ):
                        should_merge = True
                        merge_reason = (
                            "gap=%dbp, flanking viral context (%d+%d flank genes)"
                            % (gap_end - gap_start, len(upstream_flank), len(downstream_flank))
                        )

                if should_merge:
                    from dataclasses import replace

                    logger.info(
                        "Merging adjacent EVEs on %s: [%d-%d] + [%d-%d] (%s)",
                        scaffold,
                        current.start,
                        current.end,
                        next_boundary.start,
                        next_boundary.end,
                        merge_reason,
                    )
                    combined_ncldv = (
                        getattr(current, "region_classification_ncldv_markers", 0)
                        + getattr(next_boundary, "region_classification_ncldv_markers", 0)
                    )
                    combined_vp_plv = (
                        getattr(current, "region_classification_vp_plv_markers", 0)
                        + getattr(next_boundary, "region_classification_vp_plv_markers", 0)
                    )
                    combined_mirus = (
                        getattr(current, "region_classification_mirus_markers", 0)
                        + getattr(next_boundary, "region_classification_mirus_markers", 0)
                    )
                    combined_sources = sorted(set(current_sources) | set(next_sources))

                    current = replace(
                        current,
                        end=next_boundary.end,
                        original_end=next_boundary.original_end,
                        confidence=max(current.confidence, next_boundary.confidence),
                        region_classification_ncldv_markers=combined_ncldv,
                        region_classification_vp_plv_markers=combined_vp_plv,
                        region_classification_mirus_markers=combined_mirus,
                        seed_sources=combined_sources,
                    )
                    total_merges += 1
                    continue

            # No merge - finalize current and move to next
            merged_scaffold.append(current)
            current = next_boundary

        # Don't forget last boundary
        merged_scaffold.append(current)
        merged_all.extend(merged_scaffold)

    if total_merges > 0:
        logger.info(
            "Merged %d adjacent EVE pairs based on viral taxonomy in gaps",
            total_merges,
        )

    return merged_all


def should_trim_gene_as_host(
    gene_taxonomy,
    control_stats,
    host_prefix: str,
    host_baseline_fingerprint: dict = None,
    min_overlap_score: float = 0.40,
    neighbor_context: Optional[dict] = None,
    unknown_host_penalty: float = 2.0,
    unknown_viral_bonus: float = 2.0,
    host_signature_model: Optional[HostSignatureModel] = None,
    host_signature_threshold: float = 0.5,
) -> tuple[bool, float]:
    """
    Decide if gene should be trimmed as host-like using fingerprint matching.

    Args:
        gene_taxonomy: GeneTaxonomy for the gene
        control_stats: ControlStats from host region
        host_prefix: Expected host prefix (e.g., "EUK__")
        host_baseline_fingerprint: Dict of {token: weight} from control genes
        min_overlap_score: Minimum overlap (0-1) to classify as host

    Returns:
        Tuple of (should_trim, overlap_score)
    """
    # One identity-qualified viral hit is enough to protect a gene from host trimming.
    if has_identity_qualified_viral_hit(
        getattr(gene_taxonomy, "top10_prefixes", []) or [],
        getattr(gene_taxonomy, "top10_pidents", []) or [],
    ):
        return False, 0.0

    # No hit → infer from local neighborhood context.
    if not gene_taxonomy.has_hit:
        if not neighbor_context:
            return False, 0.0

        host_score = max(0.0, float(neighbor_context.get("host_score", 0.0)))
        viral_score = max(0.0, float(neighbor_context.get("viral_score", 0.0)))
        sampled_neighbors = int(neighbor_context.get("neighbors", 0))
        if sampled_neighbors <= 0:
            return False, 0.0

        adjusted_host = host_score * max(0.0, unknown_host_penalty)
        adjusted_viral = viral_score * max(0.0, unknown_viral_bonus)
        unknown_score = adjusted_host - adjusted_viral
        return unknown_score > 0.0 and adjusted_host > 0.0, unknown_score

    top10_prefixes = getattr(gene_taxonomy, "top10_prefixes", []) or []
    has_host_prefix = bool(
        gene_taxonomy.top1_prefix == host_prefix
        or any(prefix == host_prefix for prefix in top10_prefixes)
    )

    phase1_score = 0.0
    phase1_host_like = False
    if host_signature_model is not None:
        phase1_score = score_host_signature_record(gene_taxonomy, host_signature_model)
        phase1_host_like = phase1_score >= host_signature_threshold and has_host_prefix

    overlap_score = 0.0
    overlap_host_like = False
    # Use fingerprint matching if available
    if host_baseline_fingerprint and gene_taxonomy.taxonomy_fingerprint:
        overlap_host_like, overlap_score = calculate_fingerprint_overlap(
            gene_taxonomy.taxonomy_fingerprint,
            host_baseline_fingerprint,
            min_overlap_score=min_overlap_score,
        )
    elif gene_taxonomy.top1_prefix == host_prefix:
        # Fallback: prefix-only matching
        overlap_host_like, overlap_score = True, 1.0

    combined_score = max(overlap_score, phase1_score)
    return (phase1_host_like or overlap_host_like), combined_score


def _build_seed_gene_order(seed_mapping: SeedGeneMapping, taxonomy_map: dict) -> list[str]:
    """Build positional gene order for one seed using available taxonomy coordinates."""
    unique_ids = set(seed_mapping.upstream_porf_ids)
    unique_ids.update(seed_mapping.eve_porf_ids)
    unique_ids.update(seed_mapping.downstream_porf_ids)

    with_pos = []
    for porf_id in unique_ids:
        tax = taxonomy_map.get(porf_id)
        if tax:
            with_pos.append((tax.start, tax.end, porf_id))
    with_pos.sort(key=lambda item: (item[0], item[1]))
    return [item[2] for item in with_pos]


def _compute_ranked_host_viral_scores(gene_taxonomy, host_prefix: str, taxonomy_weight_mode: str) -> tuple[float, float]:
    """Compute host and viral support from ranked top-k prefixes."""
    prefixes = gene_taxonomy.top10_prefixes or []
    if not prefixes:
        host = 1.0 if gene_taxonomy.top1_prefix == host_prefix else 0.0
        viral = 1.0 if gene_taxonomy.has_viral else 0.0
        return host, viral

    bitscores = gene_taxonomy.top10_bits or []
    pidents = gene_taxonomy.top10_pidents or []
    host_score = 0.0
    viral_score = 0.0
    for rank, prefix in enumerate(prefixes):
        bits = bitscores[rank] if rank < len(bitscores) else 0.0
        weight = compute_hit_weight(rank, bits, taxonomy_weight_mode)
        if prefix == host_prefix:
            host_score += weight
        if prefix in VIRAL_PREFIXES:
            try:
                pident = float(pidents[rank])
            except (IndexError, TypeError, ValueError):
                pident = 0.0
            if pident >= MIN_VIRAL_HIT_PIDENT:
                viral_score += weight
    return host_score, viral_score


def _build_unknown_neighbor_context(
    porf_id: str,
    gene_order: list[str],
    gene_index: dict[str, int],
    taxonomy_map: dict,
    host_prefix: str,
    host_baseline_fingerprint: dict,
    min_overlap_score: float,
    taxonomy_weight_mode: str,
    neighbor_window: int,
) -> Optional[dict]:
    """Summarize host/viral neighborhood support for a no-hit gene."""
    center_idx = gene_index.get(porf_id)
    if center_idx is None or neighbor_window < 1:
        return None

    host_score = 0.0
    viral_score = 0.0
    sampled_neighbors = 0
    for delta in range(1, neighbor_window + 1):
        decay = 1.0 / float(delta)
        for idx in (center_idx - delta, center_idx + delta):
            if idx < 0 or idx >= len(gene_order):
                continue
            neighbor_tax = taxonomy_map.get(gene_order[idx])
            if not neighbor_tax or not neighbor_tax.has_hit:
                continue
            sampled_neighbors += 1
            host_rank, viral_rank = _compute_ranked_host_viral_scores(
                neighbor_tax,
                host_prefix,
                taxonomy_weight_mode,
            )
            host_score += host_rank * decay
            viral_score += viral_rank * decay
            if host_baseline_fingerprint and neighbor_tax.taxonomy_fingerprint:
                is_host_like, overlap_score = calculate_fingerprint_overlap(
                    neighbor_tax.taxonomy_fingerprint,
                    host_baseline_fingerprint,
                    min_overlap_score=min_overlap_score,
                )
                if is_host_like:
                    host_score += overlap_score * decay

    if sampled_neighbors == 0:
        return None
    return {
        "neighbors": sampled_neighbors,
        "host_score": host_score,
        "viral_score": viral_score,
    }


_DENSITY_LOOKAHEAD = 5  # check next N genes in walk direction
_DENSITY_HOST_BLOCK = 3  # N consecutive host-like genes block viral signal


def _classify_gene_label(
    tax,
    control_stats,
    host_prefix: str,
    host_baseline_fingerprint: dict,
    min_overlap_score: float,
    neighbor_context,
    unknown_host_penalty: float,
    unknown_viral_bonus: float,
    host_signature_model,
    host_signature_threshold: float,
) -> str:
    """Classify a gene as V (strong viral), H (host-like), A (ambiguous), or N (no hit).

    V: at least one viral top-10 Diamond hit with pident >= 25%.
    H: should_trim_gene_as_host() returns True.
    A: has Diamond hit but not host-like.
    N: no Diamond hit and not classified as host by neighbor context.
    """
    if has_identity_qualified_viral_hit(
        getattr(tax, "top10_prefixes", []) or [],
        getattr(tax, "top10_pidents", []) or [],
    ):
        return "V"

    trim, _ = should_trim_gene_as_host(
        tax,
        control_stats,
        host_prefix,
        host_baseline_fingerprint,
        min_overlap_score,
        neighbor_context=neighbor_context,
        unknown_host_penalty=unknown_host_penalty,
        unknown_viral_bonus=unknown_viral_bonus,
        host_signature_model=host_signature_model,
        host_signature_threshold=host_signature_threshold,
    )
    if trim:
        return "H"

    return "N" if not tax.has_hit else "A"


def _should_continue_trimming(
    labels: list[str],
    current_idx: int,
) -> bool:
    """Decide whether to continue trimming past a non-host gene using density rules.

    Called when the per-gene check says "don't trim" (label is A or N).
    Looks ahead in the walk direction (remaining labels after current_idx)
    to decide if this gene is in host territory or EVE territory.

    Rules:
    - If a V gene exists at current or ±1 position → stop (EVE territory)
    - If no V in the next 5 genes → continue trimming (host territory)
    - If V in next 5 but ≥3 consecutive H genes before it → continue (host corridor blocks)
    - If V in next 5 and <3 consecutive H before it → stop (EVE territory)
    """
    n = len(labels)

    # V at current position → definitely EVE, stop
    if labels[current_idx] == "V":
        return False

    # V at ±1 position → EVE territory, stop
    for delta in (-1, 1):
        j = current_idx + delta
        if 0 <= j < n and labels[j] == "V":
            return False

    # Look ahead: remaining genes AFTER current in walk direction
    ahead = labels[current_idx + 1: current_idx + 1 + _DENSITY_LOOKAHEAD]

    # No V in next 5 → host territory, continue trimming
    if "V" not in ahead:
        return True

    # V exists in next 5 — check if ≥3 consecutive H genes block it
    consecutive_h = 0
    for lbl in ahead:
        if lbl == "V":
            break
        if lbl == "H":
            consecutive_h += 1
        else:
            consecutive_h = 0  # A or N breaks the host corridor

    if consecutive_h >= _DENSITY_HOST_BLOCK:
        return True  # host corridor blocks viral signal

    return False  # viral is close and unblocked → EVE territory


def _select_host_taxonomy_trim(
    *,
    boundary,
    counterfactual_start: int,
    counterfactual_end: int,
    counterfactual_stats: dict,
    ablation_id: AblationID,
) -> tuple[int, int, dict]:
    """Select Phase-2f coordinates and attach candidate-level A4 evidence."""

    counterfactual_changed = (
        counterfactual_start != boundary.start
        or counterfactual_end != boundary.end
    )
    counts = (
        InterventionCounts(
            opportunities=1,
            interventions=int(counterfactual_changed),
            changed=int(counterfactual_changed),
        )
        if ablation_id is AblationID.A4
        else InterventionCounts()
    )
    if ablation_id is AblationID.A4:
        selected_start, selected_end = boundary.start, boundary.end
    else:
        selected_start, selected_end = counterfactual_start, counterfactual_end

    stats = dict(counterfactual_stats)
    stats.update(
        {
            "ablation_id": ablation_id.value,
            "counterfactual_start": counterfactual_start,
            "counterfactual_end": counterfactual_end,
            "counterfactual_trimmed": counterfactual_changed,
            "trimmed": (
                selected_start != boundary.start
                or selected_end != boundary.end
            ),
            "host_coordinate_change_opportunities": counts.opportunities,
            "host_coordinate_change_interventions": counts.interventions,
            "host_coordinate_change_changed": counts.changed,
        }
    )
    if ablation_id is AblationID.A4 and counterfactual_changed:
        stats["reason"] = "a4_host_coordinate_change_bypass"
    return selected_start, selected_end, stats


def trim_boundary_by_host_taxonomy(
    boundary,
    seed_mapping,
    taxonomy_map: dict,
    control_stats,
    host_prefix: str,
    host_baseline_fingerprint: dict = None,
    host_signature_model: Optional[HostSignatureModel] = None,
    host_signature_threshold: float = 0.5,
    min_overlap_score: float = 0.40,
    taxonomy_weight_mode: str = "rank",
    unknown_neighbor_window: int = 1,
    unknown_host_penalty: float = 2.0,
    unknown_viral_bonus: float = 2.0,
    ablation_id: AblationID = AblationID.A0,
) -> tuple:
    """
    Trim boundary inward from flanks AND inside EVE region based on host taxonomy.

    Strategy:
    - Start from boundary edges, move inward
    - For each flanking gene (closest to edge first):
        - If fingerprint matches host baseline → trim it away
        - If viral signal or low overlap → STOP trimming
    - ALSO trim genes INSIDE the EVE region (seed) that are host-like
    - Return trimmed (start, end) and trim statistics

    Uses taxonomy fingerprint matching for robust host detection:
    - Aggregates taxonomy tokens from top-10 Diamond hits
    - Compares against host baseline fingerprint from control genes
    - Fallback to prefix matching if fingerprints unavailable

    Args:
        boundary: RefinedBoundary object
        seed_mapping: SeedGeneMapping with flanking and EVE gene lists
        taxonomy_map: Dict of pORF ID → GeneTaxonomy
        control_stats: ControlStats from control region
        host_prefix: Host taxonomy prefix
        host_baseline_fingerprint: Optional dict of {token: weight} from control genes
        host_signature_model: Optional Phase 1 host-signature model
        host_signature_threshold: Host-like threshold for Phase 1 model score
        min_overlap_score: Minimum overlap (0-1) to classify as host

    Returns:
        Tuple of (trimmed_start, trimmed_end, stats_dict)
    """
    if not isinstance(ablation_id, AblationID):
        raise TypeError("ablation_id must be an AblationID")
    if not seed_mapping or not taxonomy_map:
        return _select_host_taxonomy_trim(
            boundary=boundary,
            counterfactual_start=boundary.start,
            counterfactual_end=boundary.end,
            counterfactual_stats={},
            ablation_id=ablation_id,
        )

    gene_order = _build_seed_gene_order(seed_mapping, taxonomy_map)
    gene_index = {porf_id: idx for idx, porf_id in enumerate(gene_order)}

    # Precompute density labels for all genes (H/V/A/N)
    gene_labels = {}
    for porf_id in gene_order:
        tax = taxonomy_map.get(porf_id)
        if not tax:
            gene_labels[porf_id] = "N"
            continue
        nc = _build_unknown_neighbor_context(
            porf_id=porf_id,
            gene_order=gene_order,
            gene_index=gene_index,
            taxonomy_map=taxonomy_map,
            host_prefix=host_prefix,
            host_baseline_fingerprint=host_baseline_fingerprint,
            min_overlap_score=min_overlap_score,
            taxonomy_weight_mode=taxonomy_weight_mode,
            neighbor_window=unknown_neighbor_window,
        )
        gene_labels[porf_id] = _classify_gene_label(
            tax,
            control_stats,
            host_prefix,
            host_baseline_fingerprint,
            min_overlap_score,
            nc,
            unknown_host_penalty,
            unknown_viral_bonus,
            host_signature_model,
            host_signature_threshold,
        )

    label_str = "".join(gene_labels.get(pid, "?") for pid in gene_order)
    logger.debug(
        "Density labels for %s (%d-%d): %s (upstream=%d, eve=%d, downstream=%d)",
        boundary.scaffold, boundary.start, boundary.end, label_str,
        len(seed_mapping.upstream_porf_ids),
        len(seed_mapping.eve_porf_ids),
        len(seed_mapping.downstream_porf_ids),
    )

    trimmed_start = boundary.start
    trimmed_end = boundary.end

    upstream_trimmed = 0
    upstream_stopped_by = None
    downstream_trimmed = 0
    downstream_stopped_by = None

    # Trim upstream (5' direction, left side)
    # upstream_porf_ids are ordered closest-to-EVE first (walking outward)
    upstream_walk_labels = [gene_labels.get(pid, "N") for pid in seed_mapping.upstream_porf_ids]
    for walk_i, porf_id in enumerate(seed_mapping.upstream_porf_ids):
        tax = taxonomy_map.get(porf_id)
        if not tax:
            continue

        # Skip genes that don't overlap current boundary
        # (already trimmed by previous iterations or CRF)
        if tax.end <= trimmed_start:
            continue  # Gene is already outside boundary

        label = gene_labels.get(porf_id, "N")
        if label == "H":
            # Host-like gene → trim it away
            new_start = max(trimmed_start, tax.end)
            if new_start > trimmed_start:
                upstream_trimmed += (new_start - trimmed_start)
                trimmed_start = new_start
                logger.debug(
                    "Trimming upstream gene %s (label=%s)", porf_id, label
                )
        else:
            # Non-host gene → check density rules before stopping
            if _should_continue_trimming(upstream_walk_labels, walk_i):
                # Host territory despite this gene — trim it anyway
                new_start = max(trimmed_start, tax.end)
                if new_start > trimmed_start:
                    upstream_trimmed += (new_start - trimmed_start)
                    trimmed_start = new_start
                    logger.debug(
                        "Trimming upstream gene %s (label=%s, density override)",
                        porf_id, label,
                    )
            else:
                upstream_stopped_by = f"{porf_id}:{tax.top1_prefix}"
                break

    # Trim downstream (3' direction, right side)
    # downstream_porf_ids are ordered closest-to-EVE first (walking outward)
    downstream_walk_labels = [gene_labels.get(pid, "N") for pid in seed_mapping.downstream_porf_ids]
    for walk_i, porf_id in enumerate(seed_mapping.downstream_porf_ids):
        tax = taxonomy_map.get(porf_id)
        if not tax:
            continue

        # Skip genes that don't overlap current boundary
        if tax.start >= trimmed_end:
            continue  # Gene is already outside boundary

        label = gene_labels.get(porf_id, "N")
        if label == "H":
            # Host-like gene → trim it away
            new_end = min(trimmed_end, tax.start)
            if new_end < trimmed_end:
                downstream_trimmed += (trimmed_end - new_end)
                trimmed_end = new_end
                logger.debug(
                    "Trimming downstream gene %s (label=%s)", porf_id, label
                )
        else:
            # Non-host gene → check density rules before stopping
            if _should_continue_trimming(downstream_walk_labels, walk_i):
                new_end = min(trimmed_end, tax.start)
                if new_end < trimmed_end:
                    downstream_trimmed += (trimmed_end - new_end)
                    trimmed_end = new_end
                    logger.debug(
                        "Trimming downstream gene %s (label=%s, density override)",
                        porf_id, label,
                    )
            else:
                downstream_stopped_by = f"{porf_id}:{tax.top1_prefix}"
                break

    # CRITICAL: Also trim genes INSIDE the EVE region (seed) that are host-like
    # This shrinks the EVE region from both edges by removing internal host genes
    eve_upstream_trimmed = 0
    eve_downstream_trimmed = 0
    eve_upstream_stopped_by = None
    eve_downstream_stopped_by = None

    if seed_mapping.eve_porf_ids:
        # Sort eve genes by position to process from edges inward
        eve_genes_with_pos = []
        for porf_id in seed_mapping.eve_porf_ids:
            tax = taxonomy_map.get(porf_id)
            if tax:
                eve_genes_with_pos.append((tax.start, tax.end, porf_id, tax))

        eve_genes_with_pos.sort(key=lambda x: x[0])  # Sort by start position

        # Trim from left edge of EVE region (earliest genes → inward)
        # Walk order: left-to-right (already sorted by position)
        eve_left_walk_labels = [gene_labels.get(pid, "N") for _, _, pid, _ in eve_genes_with_pos]
        for walk_i, (start, end, porf_id, tax) in enumerate(eve_genes_with_pos):
            # Skip genes already outside boundary
            if end <= trimmed_start or start >= trimmed_end:
                continue

            label = gene_labels.get(porf_id, "N")
            if label == "H":
                new_start = max(trimmed_start, end)
                if new_start > trimmed_start:
                    eve_upstream_trimmed += (new_start - trimmed_start)
                    trimmed_start = new_start
                    logger.debug(
                        "Trimming EVE upstream gene %s (label=%s)", porf_id, label
                    )
            else:
                if _should_continue_trimming(eve_left_walk_labels, walk_i):
                    new_start = max(trimmed_start, end)
                    if new_start > trimmed_start:
                        eve_upstream_trimmed += (new_start - trimmed_start)
                        trimmed_start = new_start
                        logger.debug(
                            "Trimming EVE upstream gene %s (label=%s, density override)",
                            porf_id, label,
                        )
                else:
                    eve_upstream_stopped_by = f"{porf_id}:{tax.top1_prefix}"
                    break

        # Trim from right edge of EVE region (latest genes → inward)
        # Walk order: right-to-left (reverse of position-sorted)
        eve_right_walk = list(reversed(eve_genes_with_pos))
        eve_right_walk_labels = [gene_labels.get(pid, "N") for _, _, pid, _ in eve_right_walk]
        for walk_i, (start, end, porf_id, tax) in enumerate(eve_right_walk):
            # Skip genes already outside boundary
            if end <= trimmed_start or start >= trimmed_end:
                continue

            label = gene_labels.get(porf_id, "N")
            if label == "H":
                new_end = min(trimmed_end, start)
                if new_end < trimmed_end:
                    eve_downstream_trimmed += (trimmed_end - new_end)
                    trimmed_end = new_end
                    logger.debug(
                        "Trimming EVE downstream gene %s (label=%s)", porf_id, label
                    )
            else:
                if _should_continue_trimming(eve_right_walk_labels, walk_i):
                    new_end = min(trimmed_end, start)
                    if new_end < trimmed_end:
                        eve_downstream_trimmed += (trimmed_end - new_end)
                        trimmed_end = new_end
                        logger.debug(
                            "Trimming EVE downstream gene %s (label=%s, density override)",
                            porf_id, label,
                        )
                else:
                    eve_downstream_stopped_by = f"{porf_id}:{tax.top1_prefix}"
                    break

    # Ensure boundary didn't invert
    if trimmed_start >= trimmed_end:
        # Trimming would eliminate boundary entirely - keep original
        return _select_host_taxonomy_trim(
            boundary=boundary,
            counterfactual_start=boundary.start,
            counterfactual_end=boundary.end,
            counterfactual_stats={
                "trimmed": False,
                "reason": "would_invert",
            },
            ablation_id=ablation_id,
        )

    # Total trimming from all sources
    total_upstream = upstream_trimmed + eve_upstream_trimmed
    total_downstream = downstream_trimmed + eve_downstream_trimmed

    stats = {
        "trimmed": total_upstream > 0 or total_downstream > 0,
        "upstream_bp": total_upstream,
        "downstream_bp": total_downstream,
        "upstream_trimmed": upstream_trimmed,
        "downstream_trimmed": downstream_trimmed,
        "eve_upstream_trimmed": eve_upstream_trimmed,
        "eve_downstream_trimmed": eve_downstream_trimmed,
        "upstream_stopped_by": upstream_stopped_by or eve_upstream_stopped_by,
        "downstream_stopped_by": downstream_stopped_by or eve_downstream_stopped_by,
    }

    return _select_host_taxonomy_trim(
        boundary=boundary,
        counterfactual_start=trimmed_start,
        counterfactual_end=trimmed_end,
        counterfactual_stats=stats,
        ablation_id=ablation_id,
    )
