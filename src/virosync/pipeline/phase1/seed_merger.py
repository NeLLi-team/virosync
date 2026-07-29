"""
Seed Merger for Phase 1.

Defines the MergedSeed dataclass consumed by Phase 2 and keeps a
merge_seeds() helper used by the compositional seeding sub-routine.
The active pipeline produces HHG seeds only; the legacy novelty and
compositional seed sources are retained in the dataclass schema as
zero placeholders but are no longer populated in the main flow.
"""

import logging
import math
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Optional

from .hhg_seeding import HHGSeed, Anchor

logger = logging.getLogger(__name__)


@dataclass
class MergedSeed:
    """
    Seed record consumed by Phase 2.

    The active workflow populates marker-derived HHG/validation fields. Legacy
    novelty/compositional fields remain for output and compatibility schemas.
    """

    scaffold: str
    start: int
    end: int
    seed_id: str = ""  # Stable ID assigned during merge: "seed_{idx}_{scaffold}_{start}"
    sources: list[str] = field(default_factory=list)  # Active workflow: ["hhg", "marker_validation"]
    hhg_score: float = 0.0
    novelty_score: float = 0.0
    compositional_score: float = 0.0  # KFD/CUB-based anomaly score
    mean_kfd: float = 0.0
    mean_composite: float = 0.0
    max_kfd: float = 0.0
    max_composite: float = 0.0
    gc_deviation: float = 0.0
    cub_deviation: float = 0.0
    n_windows: int = 0
    cluster_ids: list[int] = field(default_factory=list)
    anchors: list[Anchor] = field(default_factory=list)
    hhg_anchors: list[Anchor] = field(default_factory=list)  # Alias for backward compatibility
    priority: float = 0.0
    confidence: str = "low"  # "high", "medium", "low"
    score: float = 0.0  # Generic score for compatibility
    # Region classification from seed markers (NCLDV, VP, PLV, MIRUS, MIXED, UNKNOWN)
    predicted_family: str = ""
    region_classification_ncldv_markers: int = 0
    region_classification_vp_plv_markers: int = 0
    region_classification_mirus_markers: int = 0
    host_trim_original_start: Optional[int] = None
    host_trim_original_end: Optional[int] = None
    host_trimmed_start: Optional[int] = None
    host_trimmed_end: Optional[int] = None
    host_trim_reason: str = ""
    host_trim_common_euk_taxonomy: str = ""

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def num_evidence_sources(self) -> int:
        """Number of independent evidence sources."""
        return len(set(self.sources))

    @property
    def is_multi_evidence(self) -> bool:
        """True if seed has evidence from multiple sources."""
        return self.num_evidence_sources >= 2

    @property
    def is_triple_evidence(self) -> bool:
        """True if seed has evidence from all three sources."""
        return self.num_evidence_sources >= 3

    @property
    def has_mcp(self) -> bool:
        """True if seed contains a Major Capsid Protein anchor (most diagnostic marker)."""
        from virosync.pipeline.phase3.mcp_detection import is_mcp_gene
        all_anchors = self.anchors + self.hhg_anchors
        return any(is_mcp_gene(a.hallmark_gene) for a in all_anchors)

    @property
    def hallmark_genes(self) -> list[str]:
        """List of unique hallmark genes from HHG anchors."""
        all_anchors = self.anchors + self.hhg_anchors
        return list(set(a.hallmark_gene for a in all_anchors))

    def compute_classification(self, seed_marker_allowlist: list[str] | None = None) -> None:
        """
        Compute region classification based on hallmark genes.

        Sets predicted_family and marker count fields based on
        the classify_region_by_markers function.

        Deduplicates by gene (porf_id): when one gene hits multiple HMM
        models, only the best-scoring model is counted.

        Args:
            seed_marker_allowlist: List of high-purity seed markers from config
        """
        from virosync.pipeline.phase1.viral_markers import get_region_classification_summary

        # Deduplicate anchors by gene: keep best-scoring model per porf_id
        all_anchors = self.anchors + self.hhg_anchors
        best_by_gene: dict[str, str] = {}
        best_score_by_gene: dict[str, float] = {}
        for anchor in all_anchors:
            gene_id = anchor.porf_id
            if gene_id not in best_by_gene or anchor.score > best_score_by_gene[gene_id]:
                best_by_gene[gene_id] = anchor.hallmark_gene
                best_score_by_gene[gene_id] = anchor.score
        marker_names = set(best_by_gene.values())
        summary = get_region_classification_summary(marker_names, seed_marker_allowlist)

        self.predicted_family = summary["classification"]
        self.region_classification_ncldv_markers = summary["ncldv_markers"]
        self.region_classification_vp_plv_markers = summary["vp_plv_markers"]
        self.region_classification_mirus_markers = summary["mirus_markers"]

    def to_bed_line(self) -> str:
        """Format as BED line for output."""
        name = f"seed_{self.scaffold}_{self.start}_{self.end}"
        score = int(min(1000, self.priority * 100))
        return f"{self.scaffold}\t{self.start}\t{self.end}\t{name}\t{score}\t."

    def to_extended_bed_line(self) -> str:
        """Format as extended BED line with additional evidence columns."""
        name = f"seed_{self.scaffold}_{self.start}_{self.end}"
        score = int(min(1000, self.priority * 100))
        sources_str = ",".join(sorted(set(self.sources)))
        hallmarks_str = ",".join(self.hallmark_genes) if self.hallmark_genes else "."
        return (
            f"{self.scaffold}\t{self.start}\t{self.end}\t{name}\t{score}\t.\t"
            f"{sources_str}\t{self.confidence}\t{hallmarks_str}\t"
            f"{self.hhg_score:.2f}\t{self.novelty_score:.2f}\t{self.compositional_score:.2f}"
        )


def merge_seeds(
    hhg_seeds: Optional[list[HHGSeed]] = None,
    novelty_seeds: Optional[list[Any]] = None,  # legacy novelty seed list; unused in active flow
    compositional_seeds: Optional[list] = None,  # CompositionalSeed
    all_seeds_list: Optional[list[MergedSeed]] = None,  # Pre-converted seeds
    merge_distance: int = 10000,
) -> list[MergedSeed]:
    """
    Merge seed records for compatibility workflows.

    The active pipeline builds MergedSeed objects directly from marker-validated
    regions; this helper remains for legacy/offline callers.

    Args:
        hhg_seeds: Seeds from HHG seeding
        novelty_seeds: Legacy novelty seed list
        compositional_seeds: Legacy compositional seed list
        all_seeds_list: Pre-converted MergedSeed list (alternative input)
        merge_distance: Maximum gap between seeds to merge

    Returns:
        List of MergedSeed objects, sorted by priority
    """
    hhg_seeds = hhg_seeds or []
    novelty_seeds = novelty_seeds or []
    compositional_seeds = compositional_seeds or []

    logger.info(
        f"Merging seeds: {len(hhg_seeds)} HHG, "
        f"{len(novelty_seeds)} novelty, "
        f"{len(compositional_seeds)} compositional"
    )

    # Convert to common format
    all_seeds = []

    # Handle pre-converted seeds
    if all_seeds_list:
        for seed in all_seeds_list:
            all_seeds.append({
                "scaffold": seed.scaffold,
                "start": seed.start,
                "end": seed.end,
                "source": seed.sources[0] if seed.sources else "unknown",
                "hhg_score": seed.hhg_score,
                "novelty_score": seed.novelty_score,
                "compositional_score": seed.compositional_score,
                "anchors": seed.anchors + seed.hhg_anchors,
            })

    for seed in hhg_seeds:
        # Use max anchor HMM bit score as hhg_score (more representative than neighbor_score)
        # Anchor scores are typically 50-500+ for true viral hallmarks
        max_anchor_score = max(a.score for a in seed.anchors) if seed.anchors else 0.0
        all_seeds.append({
            "scaffold": seed.scaffold,
            "start": seed.start,
            "end": seed.end,
            "source": "hhg",
            "hhg_score": max_anchor_score,
            "novelty_score": 0.0,
            "compositional_score": 0.0,
            "mean_kfd": 0.0,
            "mean_composite": 0.0,
            "max_kfd": 0.0,
            "max_composite": 0.0,
            "gc_deviation": 0.0,
            "n_windows": 0,
            "cluster_ids": [],
            "anchors": seed.anchors,
        })

    for seed in novelty_seeds:
        all_seeds.append({
            "scaffold": seed.scaffold,
            "start": seed.start,
            "end": seed.end,
            "source": "novelty",
            "hhg_score": 0.0,
            "novelty_score": seed.mean_novelty,
            "compositional_score": 0.0,
            "mean_kfd": 0.0,
            "mean_composite": 0.0,
            "max_kfd": 0.0,
            "max_composite": 0.0,
            "gc_deviation": 0.0,
            "n_windows": 0,
            "cluster_ids": [],
            "anchors": [],
        })

    for seed in compositional_seeds:
        all_seeds.append({
            "scaffold": seed.scaffold,
            "start": seed.start,
            "end": seed.end,
            "source": "compositional",
            "hhg_score": 0.0,
            "novelty_score": 0.0,
            "compositional_score": seed.mean_composite if hasattr(seed, 'mean_composite') else seed.score,
            "mean_kfd": getattr(seed, "mean_kfd", 0.0),
            "mean_composite": getattr(seed, "mean_composite", 0.0),
            "max_kfd": getattr(seed, "max_kfd", 0.0),
            "max_composite": getattr(seed, "max_composite", 0.0),
            "gc_deviation": getattr(seed, "gc_deviation", 0.0),
            "cub_deviation": getattr(seed, "cub_deviation", 0.0),
            "n_windows": getattr(seed, "n_windows", 0),
            "cluster_ids": getattr(seed, "cluster_ids", []),
            "anchors": [],
        })

    if not all_seeds:
        logger.info("No seeds to merge")
        return []

    # Group by scaffold
    scaffold_seeds = defaultdict(list)
    for seed in all_seeds:
        scaffold_seeds[seed["scaffold"]].append(seed)

    # Merge overlapping seeds per scaffold
    merged = []

    for scaffold, seeds in scaffold_seeds.items():
        # Sort by start position
        seeds.sort(key=lambda s: s["start"])

        current = seeds[0].copy()
        current["sources"] = [current.pop("source")]

        for seed in seeds[1:]:
            # Check for overlap or proximity
            if seed["start"] - current["end"] <= merge_distance:
                # Merge seeds
                current["end"] = max(current["end"], seed["end"])
                current["sources"].append(seed["source"])
                current["hhg_score"] = max(current["hhg_score"], seed["hhg_score"])
                current["novelty_score"] = max(current["novelty_score"], seed["novelty_score"])
                current["compositional_score"] = max(
                    current.get("compositional_score", 0.0),
                    seed.get("compositional_score", 0.0)
                )
                current["mean_kfd"] = max(
                    current.get("mean_kfd", 0.0),
                    seed.get("mean_kfd", 0.0),
                )
                current["mean_composite"] = max(
                    current.get("mean_composite", 0.0),
                    seed.get("mean_composite", 0.0),
                )
                current["max_kfd"] = max(
                    current.get("max_kfd", 0.0),
                    seed.get("max_kfd", 0.0),
                )
                current["max_composite"] = max(
                    current.get("max_composite", 0.0),
                    seed.get("max_composite", 0.0),
                )
                current["gc_deviation"] = max(
                    current.get("gc_deviation", 0.0),
                    seed.get("gc_deviation", 0.0),
                )
                current["cub_deviation"] = max(
                    current.get("cub_deviation", 0.0),
                    seed.get("cub_deviation", 0.0),
                )
                current["n_windows"] = current.get("n_windows", 0) + seed.get("n_windows", 0)
                current["cluster_ids"] = sorted(
                    set(current.get("cluster_ids", [])) | set(seed.get("cluster_ids", []))
                )
                current["anchors"].extend(seed["anchors"])
            else:
                # Emit current and start new
                merged.append(_create_merged_seed(scaffold, current))
                current = seed.copy()
                current["sources"] = [current.pop("source")]

        # Emit final seed
        merged.append(_create_merged_seed(scaffold, current))

    # Sort by priority (descending)
    merged.sort(key=lambda s: s.priority, reverse=True)

    # Assign stable seed_id after sorting (format: "seed_{idx}_{scaffold}_{start}")
    for idx, seed in enumerate(merged):
        seed.seed_id = f"seed_{idx}_{seed.scaffold}_{seed.start}"

    # Log summary
    triple_evidence = sum(1 for s in merged if s.is_triple_evidence)
    multi_evidence = sum(1 for s in merged if s.is_multi_evidence and not s.is_triple_evidence)

    # Count by source
    hhg_only = sum(1 for s in merged if set(s.sources) == {"hhg"})
    novelty_only = sum(1 for s in merged if set(s.sources) == {"novelty"})
    compositional_only = sum(1 for s in merged if set(s.sources) == {"compositional"})
    with_mcp = sum(1 for s in merged if s.has_mcp)

    # Confidence breakdown
    high_conf = sum(1 for s in merged if s.confidence == "high")
    medium_conf = sum(1 for s in merged if s.confidence == "medium")
    low_conf = sum(1 for s in merged if s.confidence == "low")

    logger.info(f"Merged to {len(merged)} seeds:")
    logger.info(f"  Triple evidence (HHG+Novelty+Compositional): {triple_evidence}")
    logger.info(f"  Multi evidence: {multi_evidence}")
    logger.info(f"  HHG only: {hhg_only}")
    logger.info(f"  Novelty only: {novelty_only}")
    logger.info(f"  Compositional only: {compositional_only}")
    logger.info(f"  With MCP (capsid): {with_mcp}")
    logger.info(f"  Confidence: {high_conf} high, {medium_conf} medium, {low_conf} low")

    if merged:
        total_bp = sum(s.length for s in merged)
        logger.info(f"  Total coverage: {total_bp:,} bp")

    return merged


def _create_merged_seed(scaffold: str, seed_dict: dict) -> MergedSeed:
    """Create MergedSeed from intermediate dictionary."""
    sources = list(set(seed_dict["sources"]))
    anchors = seed_dict["anchors"]
    compositional_score = seed_dict.get("compositional_score", 0.0)

    # Check for MCP (Major Capsid Protein) - most diagnostic marker
    from virosync.pipeline.phase3.mcp_detection import is_mcp_gene
    has_mcp = any(is_mcp_gene(a.hallmark_gene) for a in anchors)

    # Calculate priority score
    # HHG (hallmark genes) is strongest evidence, followed by compositional, then novelty
    hhg_weight = 1.0
    novelty_weight = 0.5
    compositional_weight = 0.7  # Compositional is more reliable than novelty alone

    priority = (
        seed_dict["hhg_score"] * hhg_weight
        + seed_dict["novelty_score"] * novelty_weight
        + compositional_score * 100 * compositional_weight  # Scale up compositional (0-1 range)
    )

    marker_types = {a.hallmark_gene.lower() for a in anchors}
    marker_count = len(marker_types)
    marker_strength = 0.0
    for marker in marker_types:
        if is_mcp_gene(marker):
            # MCP is the most diagnostic viral capsid hallmark across giant-virus
            # and virophage-related families; weight all canonical MCP models
            # (GVOGm0003, VS000086/OG1352, VS000309/OG484, gamadvirusMCP, plus family-scoped
            # prefixes like plv_mcp / vp_mcp / mirus_mcp) uniformly at 3.0.
            marker_strength += 3.0
        elif marker.startswith("gvogm"):
            marker_strength += 1.5
        elif marker.startswith("og"):
            marker_strength += 1.2
        else:
            marker_strength += 1.0

    mean_kfd = seed_dict.get("mean_kfd", 0.0)
    gc_dev = seed_dict.get("gc_deviation", 0.0)
    comp_strength = max(compositional_score, mean_kfd * 3.0, gc_dev * 5.0)
    comp_strength = min(comp_strength, 1.0)

    length_kb = max(seed_dict["end"] - seed_dict["start"], 0) / 1000.0
    length_bonus = min(length_kb / 100.0, 1.0)  # cap at 100kb

    marker_bonus = 1.0 + min(marker_count, 5) * 0.1 + min(marker_strength, 6.0) * 0.1
    priority *= (1.0 + comp_strength * 0.5 + length_bonus * 0.5)
    priority *= marker_bonus

    # Bonus for multi-evidence
    num_sources = len(sources)
    if num_sources == 2:
        priority *= 1.5
    elif num_sources >= 3:
        priority *= 2.0  # Triple evidence is very strong

    # MCP bonus: seeds with MCP are much more likely to be real EVEs
    if has_mcp:
        priority *= 2.0

    # Calculate confidence level based on evidence
    has_hhg = "hhg" in sources
    has_novelty = "novelty" in sources
    has_compositional = "compositional" in sources

    # Triple evidence = high confidence
    if num_sources >= 3:
        confidence = "high"
    elif has_hhg and (has_novelty or has_compositional) and has_mcp:
        confidence = "high"
    elif has_hhg and (has_novelty or has_compositional):
        confidence = "medium"
    elif has_hhg and has_mcp:
        confidence = "medium"  # MCP is strong single evidence
    elif has_compositional and has_novelty:
        confidence = "medium"  # Two non-HHG sources together
    else:
        confidence = "low"  # Single evidence only

    confidence_score = 1.0 - math.exp(-priority / 3000.0) if priority > 0 else 0.0

    return MergedSeed(
        scaffold=scaffold,
        start=seed_dict["start"],
        end=seed_dict["end"],
        sources=sources,
        hhg_score=seed_dict["hhg_score"],
        novelty_score=seed_dict["novelty_score"],
        compositional_score=compositional_score,
        mean_kfd=seed_dict.get("mean_kfd", 0.0),
        mean_composite=seed_dict.get("mean_composite", 0.0),
        max_kfd=seed_dict.get("max_kfd", 0.0),
        max_composite=seed_dict.get("max_composite", 0.0),
        gc_deviation=seed_dict.get("gc_deviation", 0.0),
        cub_deviation=seed_dict.get("cub_deviation", 0.0),
        n_windows=seed_dict.get("n_windows", 0),
        cluster_ids=seed_dict.get("cluster_ids", []),
        anchors=anchors,
        hhg_anchors=[],  # For backward compatibility
        priority=priority,
        confidence=confidence,
        score=confidence_score,
    )


def filter_seeds_by_priority(
    seeds: list[MergedSeed],
    min_priority: float = 0.0,
    max_seeds: int = None,
) -> list[MergedSeed]:
    """
    Filter seeds by priority threshold and/or count limit.

    Args:
        seeds: List of MergedSeed objects
        min_priority: Minimum priority score to include
        max_seeds: Maximum number of seeds to return

    Returns:
        Filtered list of MergedSeed objects
    """
    filtered = [s for s in seeds if s.priority >= min_priority]

    if max_seeds is not None and len(filtered) > max_seeds:
        filtered = filtered[:max_seeds]

    return filtered


def write_seeds_bed(seeds: list[MergedSeed], output_path) -> None:
    """
    Write seeds to BED format.

    BED columns:
    1. chrom
    2. start
    3. end
    4. name (seed_scaffold_start_end)
    5. score (priority * 100, max 1000)
    6. strand (.)

    Args:
        seeds: List of MergedSeed objects
        output_path: Output BED file path
    """
    from pathlib import Path

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        for seed in seeds:
            f.write(seed.to_bed_line() + "\n")

    logger.info(f"Wrote {len(seeds)} seeds to {output_path}")
