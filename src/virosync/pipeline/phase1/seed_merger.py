"""Seed Merger for Phase 1.

Defines the MergedSeed dataclass consumed by Phase 2.
The active pipeline produces HHG seeds only; the legacy novelty and
compositional seed sources are retained in the dataclass schema as
zero placeholders but are no longer populated in the main flow.
"""

import logging
from collections.abc import Sequence
from dataclasses import dataclass, field

from .hhg_seeding import Anchor

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class ConfirmedSeed:
    """Complete Phase 1 interval and immutable original identity."""

    seed_id: str
    scaffold: str
    start: int
    end: int


@dataclass
class MergedSeed:
    """Seed record consumed by Phase 2.

    The active workflow populates marker-derived HHG/validation fields. Legacy
    novelty/compositional fields remain for output and compatibility schemas.
    """

    scaffold: str
    start: int
    end: int
    confirmed_seeds: tuple[ConfirmedSeed, ...] = field(default=(), kw_only=True)
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
        """Compute region classification based on hallmark genes.

        Sets predicted_family and marker count fields based on
        the classify_region_by_markers function.

        Deduplicates by gene (porf_id): when one gene hits multiple HMM
        models, the highest HMM score selects one model per base protein ID.
        A selected model that is a
        structure-derived profile with only EUK-pEVE references counts as a marker
        of its training-EVE class.

        Args:
            seed_marker_allowlist: List of high-purity seed markers from config
        """
        summary = classify_marker_hits(self.anchors + self.hhg_anchors, seed_marker_allowlist)

        self.predicted_family = summary["classification"]
        self.region_classification_ncldv_markers = summary["ncldv_markers"]
        self.region_classification_vp_plv_markers = summary["vp_plv_markers"]
        self.region_classification_mirus_markers = summary["mirus_markers"]


def validate_confirmed_seed_containment(
    scaffold: str, start: int, end: int, confirmed_seeds: tuple[ConfirmedSeed, ...]
) -> None:
    """Reject missing, duplicate, relocated or truncated confirmed source intervals."""
    if not confirmed_seeds:
        raise ValueError("candidate has no confirmed seed membership")
    identities: set[str] = set()
    for seed in confirmed_seeds:
        if not seed.seed_id or seed.seed_id in identities:
            raise ValueError(f"missing or duplicate confirmed seed identity: {seed.seed_id!r}")
        identities.add(seed.seed_id)
        if seed.scaffold != scaffold or not (start <= seed.start < seed.end <= end):
            raise ValueError(f"candidate does not contain complete confirmed seed {seed.seed_id!r}")


def classify_marker_hits(hits: Sequence[object], seed_marker_allowlist: list[str] | None = None) -> dict:
    """Apply existing family rules to the highest-scoring hit per base protein."""
    from virosync.pipeline.phase1.viral_markers import base_marker_gene_id, get_region_classification_summary
    from virosync.pipeline.phase3.gene_taxonomy import peve_profile_class

    best_by_gene: dict[str, tuple[float, object]] = {}
    for hit in hits:
        protein = str(getattr(hit, "porf_id", "") or getattr(hit, "query_porf", ""))
        gene = base_marker_gene_id(protein)
        score = float(getattr(hit, "hmm_score", getattr(hit, "score", 0.0)))
        current = best_by_gene.get(gene)
        if current is None or score > current[0]:
            best_by_gene[gene] = (score, hit)
    selected = [hit for _score, hit in best_by_gene.values()]
    names_and_hits = [
        (str(getattr(hit, "hallmark_gene", "") or getattr(hit, "hmm_target", "")), hit) for hit in selected
    ]
    class_by_marker = {
        name.lower(): profile_class for name, hit in names_and_hits if (profile_class := peve_profile_class(name, hit))
    }
    return get_region_classification_summary(
        {name for name, _hit in names_and_hits}, seed_marker_allowlist, class_by_marker=class_by_marker
    )
