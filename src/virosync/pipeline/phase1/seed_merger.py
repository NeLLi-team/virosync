"""
Seed Merger for Phase 1.

Defines the MergedSeed dataclass consumed by Phase 2.
The active pipeline produces HHG seeds only; the legacy novelty and
compositional seed sources are retained in the dataclass schema as
zero placeholders but are no longer populated in the main flow.
"""

import logging
from dataclasses import dataclass, field
from typing import Optional

from .hhg_seeding import Anchor

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
    def has_mcp(self) -> bool:
        """True if seed contains a Major Capsid Protein anchor (most diagnostic marker)."""
        from virosync.pipeline.phase3.mcp_detection import is_mcp_gene
        all_anchors = self.anchors + self.hhg_anchors
        return any(is_mcp_gene(a.hallmark_gene) for a in all_anchors)

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
