"""
Evidence Correlation Graph for EVE Verification.

Builds a graph representation of co-occurring evidence types
within predicted EVE regions and calculates a Coherence Score
that measures how well different evidence supports the prediction.

A high Coherence Score indicates logically consistent evidence
(e.g., hallmark genes co-occur with high KFD and novelty).
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional

import numpy as np

from virosync.pipeline.phase3.mcp_detection import is_mcp_gene

try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    nx = None

logger = logging.getLogger(__name__)


class EvidenceType(Enum):
    """Types of evidence for EVE prediction."""

    # Hallmark gene evidence
    HALLMARK_MCP = "hallmark_mcp"  # Major Capsid Protein
    HALLMARK_A32 = "hallmark_a32"  # Packaging ATPase
    HALLMARK_D5 = "hallmark_d5"  # Primase-helicase
    HALLMARK_VLTF3 = "hallmark_vltf3"  # Late transcription factor
    HALLMARK_POLB = "hallmark_polb"  # DNA Polymerase B
    HALLMARK_RNAPL = "hallmark_rnapl"  # RNA Polymerase Large
    HALLMARK_RNAPS = "hallmark_rnaps"  # RNA Polymerase Small
    HALLMARK_OTHER = "hallmark_other"  # Other hallmark genes

    # Compositional evidence
    HIGH_KFD = "high_kfd"  # K-mer frequency deviation
    HIGH_CUB = "high_cub"  # Codon usage bias
    ANOMALOUS_GC = "anomalous_gc"  # GC content anomaly

    # Homology evidence
    HIGH_NOVELTY = "high_novelty"  # ORFan/novel genes
    VIRAL_BLAST_HIT = "viral_blast_hit"  # BLAST hit to viral proteins
    NO_HOST_HIT = "no_host_hit"  # No hit to host proteins

    # Structural evidence
    VIRAL_STRUCTURE = "viral_structure"  # Structural homology to viral proteins
    CONFIDENT_FOLD = "confident_fold"  # High-confidence structure prediction

    # Legacy boundary-state evidence. The enum values keep their historical
    # names so old evidence_graph.json consumers remain compatible.
    CRF_CORE_VIRAL = "crf_core_viral"  # Boundary state classified as core viral
    CRF_VIRAL_FLANK = "crf_viral_flank"  # Boundary state classified as viral flank
    CRF_HIGH_CONFIDENCE = "crf_high_confidence"  # High boundary-state posterior

    # Genomic context
    GENE_CLUSTER = "gene_cluster"  # Multiple viral genes clustered
    SYNTENY_BREAK = "synteny_break"  # Synteny disruption with relatives


@dataclass
class WindowEvidence:
    """Evidence observed in a single genomic window."""

    scaffold: str
    start: int
    end: int
    evidence_types: set[EvidenceType] = field(default_factory=set)
    evidence_scores: dict[EvidenceType, float] = field(default_factory=dict)

    def add_evidence(self, etype: EvidenceType, score: float = 1.0) -> None:
        """Add evidence type with optional score."""
        self.evidence_types.add(etype)
        self.evidence_scores[etype] = score

    @property
    def evidence_count(self) -> int:
        return len(self.evidence_types)


@dataclass
class EvidenceProfile:
    """
    Complete evidence profile for an EVE candidate.

    Contains window-by-window evidence and aggregate metrics.
    """

    eve_id: str
    scaffold: str
    start: int
    end: int

    # Window-level evidence
    windows: list[WindowEvidence] = field(default_factory=list)

    # Aggregate evidence counts
    evidence_counts: dict[EvidenceType, int] = field(default_factory=dict)
    evidence_coverage: dict[EvidenceType, float] = field(default_factory=dict)

    # Derived metrics
    hallmark_diversity: int = 0
    has_virus_specific_marker: bool = False
    multi_evidence_windows: int = 0

    def compute_aggregates(self) -> None:
        """Compute aggregate metrics from window evidence."""
        if not self.windows:
            return

        # Count evidence occurrences. Keys are ordered by enum value so the
        # serialized dicts do not inherit set iteration order, which varies
        # with PYTHONHASHSEED.
        counts: dict[EvidenceType, int] = {}
        for window in self.windows:
            for etype in window.evidence_types:
                counts[etype] = counts.get(etype, 0) + 1
        self.evidence_counts = dict(
            sorted(counts.items(), key=lambda item: item[0].value)
        )

        # Calculate coverage (fraction of windows with each evidence)
        n_windows = len(self.windows)
        self.evidence_coverage = {
            etype: count / n_windows for etype, count in self.evidence_counts.items()
        }

        # Hallmark diversity
        hallmark_types = [
            e for e in self.evidence_counts.keys() if e.value.startswith("hallmark_")
        ]
        self.hallmark_diversity = len(hallmark_types)

        # Virus-specific markers
        virus_specific = {
            EvidenceType.HALLMARK_MCP,
            EvidenceType.HALLMARK_A32,
            EvidenceType.HALLMARK_D5,
            EvidenceType.HALLMARK_VLTF3,
        }
        self.has_virus_specific_marker = bool(
            set(self.evidence_counts.keys()) & virus_specific
        )

        # Multi-evidence windows
        self.multi_evidence_windows = sum(
            1 for w in self.windows if w.evidence_count >= 2
        )


class EvidenceCorrelationGraph:
    """
    Graph representation of evidence co-occurrence patterns.

    Nodes represent evidence types, edges represent co-occurrence
    within the same or adjacent windows. Edge weights encode
    co-occurrence strength.
    """

    def __init__(self):
        """Initialize evidence correlation graph."""
        if not HAS_NETWORKX:
            raise ImportError("networkx required. Install with: pip install networkx")

        self.graph = nx.Graph()
        self._cooccurrence_matrix = None

    def build_from_profile(
        self,
        profile: EvidenceProfile,
        adjacency_distance: int = 1,
    ) -> None:
        """
        Build graph from evidence profile.

        Args:
            profile: EvidenceProfile with window-level evidence
            adjacency_distance: How many windows apart to consider co-occurrence
        """
        self.graph.clear()

        # Add nodes for all evidence types present. Ordering by enum value keeps
        # node and edge insertion order independent of set iteration order.
        all_evidence = set()
        for window in profile.windows:
            all_evidence.update(window.evidence_types)
        evidence_list = sorted(all_evidence, key=lambda e: e.value)

        for etype in evidence_list:
            self.graph.add_node(
                etype.value,
                count=profile.evidence_counts.get(etype, 0),
                coverage=profile.evidence_coverage.get(etype, 0.0),
            )

        # Build co-occurrence matrix
        n_windows = len(profile.windows)
        n_evidence = len(evidence_list)

        if n_evidence == 0:
            return

        cooccur = np.zeros((n_evidence, n_evidence))

        for i, window in enumerate(profile.windows):
            # Evidence in this window
            window_evidence = [
                evidence_list.index(e)
                for e in window.evidence_types
                if e in evidence_list
            ]

            # Co-occurrence within window
            for e1 in window_evidence:
                for e2 in window_evidence:
                    cooccur[e1, e2] += 1

            # Co-occurrence with adjacent windows
            for d in range(1, adjacency_distance + 1):
                if i + d < n_windows:
                    adj_evidence = [
                        evidence_list.index(e)
                        for e in profile.windows[i + d].evidence_types
                        if e in evidence_list
                    ]
                    for e1 in window_evidence:
                        for e2 in adj_evidence:
                            # Decay weight with distance
                            weight = 1.0 / (d + 1)
                            cooccur[e1, e2] += weight
                            cooccur[e2, e1] += weight

        self._cooccurrence_matrix = cooccur

        # Add edges for significant co-occurrences
        for i in range(n_evidence):
            for j in range(i + 1, n_evidence):
                if cooccur[i, j] > 0:
                    self.graph.add_edge(
                        evidence_list[i].value,
                        evidence_list[j].value,
                        weight=cooccur[i, j],
                        normalized_weight=cooccur[i, j] / n_windows,
                    )

    def compute_coherence_score(self) -> float:
        """
        Compute Coherence Score measuring evidence consistency.

        The score considers:
        1. Evidence diversity (more types = higher)
        2. Co-occurrence density (more connections = higher)
        3. Expected correlations (hallmarks + compositional + structural)
        4. Penalties for contradictory evidence

        Returns:
            Coherence score between 0 and 1
        """
        if self.graph.number_of_nodes() == 0:
            return 0.0

        scores = []

        # 1. Diversity score (0-0.25)
        n_nodes = self.graph.number_of_nodes()
        diversity_score = min(0.25, n_nodes / 20)  # Cap at 5 evidence types
        scores.append(diversity_score)

        # 2. Connectivity score (0-0.25)
        if n_nodes > 1:
            max_edges = n_nodes * (n_nodes - 1) / 2
            actual_edges = self.graph.number_of_edges()
            connectivity = actual_edges / max_edges
            connectivity_score = 0.25 * connectivity
        else:
            connectivity_score = 0.0
        scores.append(connectivity_score)

        # 3. Expected correlation bonus (0-0.3)
        expected_pairs = [
            # Hallmarks should co-occur with compositional anomalies
            (EvidenceType.HALLMARK_MCP.value, EvidenceType.HIGH_KFD.value),
            (EvidenceType.HALLMARK_POLB.value, EvidenceType.HIGH_KFD.value),
            # Hallmarks should co-occur with novelty
            (EvidenceType.HALLMARK_MCP.value, EvidenceType.HIGH_NOVELTY.value),
            (EvidenceType.HALLMARK_A32.value, EvidenceType.HIGH_NOVELTY.value),
            # CRF core viral should align with hallmarks
            (EvidenceType.CRF_CORE_VIRAL.value, EvidenceType.HALLMARK_MCP.value),
            (EvidenceType.CRF_CORE_VIRAL.value, EvidenceType.HIGH_KFD.value),
            # Structural evidence should align
            (EvidenceType.VIRAL_STRUCTURE.value, EvidenceType.HALLMARK_MCP.value),
            (EvidenceType.VIRAL_STRUCTURE.value, EvidenceType.HIGH_NOVELTY.value),
        ]

        expected_bonus = 0.0
        for e1, e2 in expected_pairs:
            if self.graph.has_edge(e1, e2):
                weight = self.graph[e1][e2].get("normalized_weight", 0)
                expected_bonus += 0.0375 * min(1.0, weight)  # 0.3 / 8 pairs

        scores.append(expected_bonus)

        # 4. Virus-specific marker bonus (0-0.2)
        virus_specific = {
            EvidenceType.HALLMARK_MCP.value,
            EvidenceType.HALLMARK_A32.value,
            EvidenceType.HALLMARK_D5.value,
            EvidenceType.HALLMARK_VLTF3.value,
        }
        vs_present = virus_specific & set(self.graph.nodes())
        virus_specific_score = 0.2 * (len(vs_present) / 4)
        scores.append(virus_specific_score)

        return sum(scores)

    def get_evidence_summary(self) -> dict:
        """Get summary statistics of the evidence graph."""
        if self.graph.number_of_nodes() == 0:
            return {
                "n_evidence_types": 0,
                "n_connections": 0,
                "coherence_score": 0.0,
                "evidence_types": [],
            }

        return {
            "n_evidence_types": self.graph.number_of_nodes(),
            "n_connections": self.graph.number_of_edges(),
            "coherence_score": self.compute_coherence_score(),
            "evidence_types": list(self.graph.nodes()),
            "strongest_connections": self._get_strongest_connections(5),
        }

    def _get_strongest_connections(self, n: int = 5) -> list[tuple[str, str, float]]:
        """Get the n strongest co-occurrence connections."""
        if self.graph.number_of_edges() == 0:
            return []

        edges = [
            (u, v, d.get("normalized_weight", 0))
            for u, v, d in self.graph.edges(data=True)
        ]
        # Total key: weight ties are broken by node names so the truncation
        # below always keeps the same pairs.
        edges.sort(key=lambda x: (-x[2], x[0], x[1]))
        return edges[:n]


@dataclass
class CoherenceAnalysis:
    """
    Complete coherence analysis for an EVE candidate.

    Combines evidence profile with graph-based coherence scoring.
    """

    eve_id: str
    profile: EvidenceProfile
    graph: EvidenceCorrelationGraph
    coherence_score: float = 0.0

    # Component scores
    diversity_score: float = 0.0
    connectivity_score: float = 0.0
    expected_correlation_score: float = 0.0
    virus_specific_score: float = 0.0

    # Interpretation
    interpretation: str = ""
    confidence_level: str = ""  # high, medium, low, reject

    def compute_coherence(self) -> None:
        """Compute coherence analysis."""
        self.graph.build_from_profile(self.profile)
        self.coherence_score = self.graph.compute_coherence_score()

        # Interpret score
        if self.coherence_score >= 0.7:
            self.confidence_level = "high"
            self.interpretation = "Strong, coherent evidence from multiple sources"
        elif self.coherence_score >= 0.5:
            self.confidence_level = "medium"
            self.interpretation = "Moderate evidence with some support"
        elif self.coherence_score >= 0.3:
            self.confidence_level = "low"
            self.interpretation = "Weak evidence, may require manual review"
        else:
            self.confidence_level = "reject"
            self.interpretation = "Insufficient coherent evidence"

    def to_dict(self) -> dict:
        """Convert coherence analysis to dictionary for JSON serialization."""
        return {
            "eve_id": self.eve_id,
            "scaffold": self.profile.scaffold,
            "start": self.profile.start,
            "end": self.profile.end,
            "coherence_score": self.coherence_score,
            "confidence_level": self.confidence_level,
            "interpretation": self.interpretation,
            "component_scores": {
                "diversity": self.diversity_score,
                "connectivity": self.connectivity_score,
                "expected_correlation": self.expected_correlation_score,
                "virus_specific": self.virus_specific_score,
            },
            "evidence_profile": {
                "hallmark_diversity": self.profile.hallmark_diversity,
                "has_virus_specific_marker": self.profile.has_virus_specific_marker,
                "multi_evidence_windows": self.profile.multi_evidence_windows,
                "evidence_counts": {
                    k.value: v for k, v in self.profile.evidence_counts.items()
                },
                "evidence_coverage": {
                    k.value: v for k, v in self.profile.evidence_coverage.items()
                },
            },
            "graph_summary": self.graph.get_evidence_summary() if self.graph else {},
            "windows": [
                {
                    "start": w.start,
                    "end": w.end,
                    "evidence_types": sorted(e.value for e in w.evidence_types),
                    "evidence_scores": {
                        k.value: v for k, v in w.evidence_scores.items()
                    },
                }
                for w in self.profile.windows
            ],
        }


def build_evidence_profile(
    eve_id: str,
    scaffold: str,
    start: int,
    end: int,
    window_features: list,
    crf_states: list[int],
    crf_posteriors: Optional[np.ndarray],
    hallmark_hits: Optional[list] = None,
    novelty_scores: Optional[dict] = None,
    structural_results: Optional[list] = None,
    window_size: int = 250,
) -> EvidenceProfile:
    """
    Build evidence profile from pipeline outputs.

    Args:
        eve_id: EVE identifier
        scaffold: Scaffold name
        start: EVE start position
        end: EVE end position
        window_features: WindowFeatures from Phase 2
        crf_states: Legacy boundary-state sequence
        crf_posteriors: Legacy boundary-state posterior probabilities
        hallmark_hits: HMM hallmark hits
        novelty_scores: Novelty scores per pORF
        structural_results: Structural homology results
        window_size: Window size for evidence aggregation

    Returns:
        EvidenceProfile with complete evidence
    """
    profile = EvidenceProfile(
        eve_id=eve_id,
        scaffold=scaffold,
        start=start,
        end=end,
    )

    # Create windows
    for i, feat in enumerate(window_features):
        window = WindowEvidence(
            scaffold=scaffold,
            start=feat.start,
            end=feat.end,
        )

        # Compositional evidence
        if feat.kfd > 0.3:
            window.add_evidence(EvidenceType.HIGH_KFD, feat.kfd)
        if hasattr(feat, "cub_deviation") and feat.cub_deviation > 0.1:
            window.add_evidence(EvidenceType.HIGH_CUB, feat.cub_deviation)
        if hasattr(feat, "gc_deviation") and feat.gc_deviation > 0.1:
            window.add_evidence(EvidenceType.ANOMALOUS_GC, feat.gc_deviation)

        # Novelty evidence
        if feat.novelty_score > 0.7:
            window.add_evidence(EvidenceType.HIGH_NOVELTY, feat.novelty_score)

        # Hallmark evidence
        if feat.has_hallmark:
            window.add_evidence(EvidenceType.HALLMARK_OTHER, 1.0)

        # Legacy boundary-state evidence
        if i < len(crf_states):
            state = crf_states[i]
            if state == 5:  # CORE_VIRAL
                window.add_evidence(EvidenceType.CRF_CORE_VIRAL)
            elif state == 4:  # VIRAL_FLANK
                window.add_evidence(EvidenceType.CRF_VIRAL_FLANK)

            # High-confidence boundary-state support
            if crf_posteriors is not None and i < len(crf_posteriors):
                # Only compute CORE/FLANK posterior when 6-state Tier-2 posteriors are present.
                if crf_posteriors.shape[1] >= 6:
                    viral_prob = crf_posteriors[i, 4] + crf_posteriors[i, 5]
                    if viral_prob > 0.9:
                        window.add_evidence(EvidenceType.CRF_HIGH_CONFIDENCE, viral_prob)

        profile.windows.append(window)

    # Add hallmark hit details
    if hallmark_hits:
        for hit in hallmark_hits:
            if isinstance(hit, dict):
                hit_start = hit.get("start", 0)
                hit_end = hit.get("end", 0)
                gene = hit.get("hallmark_gene", "")
            else:
                hit_start = getattr(hit, "start", 0)
                hit_end = getattr(hit, "end", 0)
                gene = getattr(hit, "hallmark_gene", "")

            # Find overlapping windows
            for window in profile.windows:
                if hit_start < window.end and hit_end > window.start:
                    gene_lower = gene.lower()

                    if is_mcp_gene(gene_lower):
                        window.add_evidence(EvidenceType.HALLMARK_MCP)
                    elif "a32" in gene_lower:
                        window.add_evidence(EvidenceType.HALLMARK_A32)
                    elif "d5" in gene_lower:
                        window.add_evidence(EvidenceType.HALLMARK_D5)
                    elif "vltf3" in gene_lower:
                        window.add_evidence(EvidenceType.HALLMARK_VLTF3)
                    elif "polb" in gene_lower:
                        window.add_evidence(EvidenceType.HALLMARK_POLB)
                    elif "rnapl" in gene_lower:
                        window.add_evidence(EvidenceType.HALLMARK_RNAPL)
                    elif "rnaps" in gene_lower:
                        window.add_evidence(EvidenceType.HALLMARK_RNAPS)

    # Add structural evidence
    if structural_results:
        for result in structural_results:
            if result.supports_viral_origin:
                # Find window containing this pORF
                for window in profile.windows:
                    window.add_evidence(
                        EvidenceType.VIRAL_STRUCTURE,
                        result.structural_evidence_score,
                    )
                    if result.prediction and result.prediction.is_confident:
                        window.add_evidence(EvidenceType.CONFIDENT_FOLD)

    # Compute aggregates
    profile.compute_aggregates()

    return profile


def analyze_eve_coherence(
    eve_id: str,
    scaffold: str,
    start: int,
    end: int,
    window_features: list,
    crf_states: list[int],
    crf_posteriors: Optional[np.ndarray] = None,
    hallmark_hits: Optional[list] = None,
    novelty_scores: Optional[dict] = None,
    structural_results: Optional[list] = None,
) -> CoherenceAnalysis:
    """
    Perform complete coherence analysis for an EVE candidate.

    Args:
        eve_id: EVE identifier
        scaffold: Scaffold name
        start: Start position
        end: End position
        window_features: Window-level evidence features
        crf_states: Legacy boundary-state sequence
        crf_posteriors: Legacy boundary-state posteriors
        hallmark_hits: Hallmark gene hits
        novelty_scores: Novelty scores
        structural_results: Structural homology results

    Returns:
        CoherenceAnalysis with complete results
    """
    # Build evidence profile
    profile = build_evidence_profile(
        eve_id=eve_id,
        scaffold=scaffold,
        start=start,
        end=end,
        window_features=window_features,
        crf_states=crf_states,
        crf_posteriors=crf_posteriors,
        hallmark_hits=hallmark_hits,
        novelty_scores=novelty_scores,
        structural_results=structural_results,
    )

    # Build graph and compute coherence
    graph = EvidenceCorrelationGraph()
    analysis = CoherenceAnalysis(
        eve_id=eve_id,
        profile=profile,
        graph=graph,
    )
    analysis.compute_coherence()

    logger.debug(
        f"EVE {eve_id}: coherence={analysis.coherence_score:.3f}, "
        f"level={analysis.confidence_level}"
    )

    return analysis


def write_evidence_graph_json(
    analyses: list[CoherenceAnalysis],
    output_path: "Path",
    genome_id: str = "",
) -> None:
    """
    Write evidence graph analyses to JSON file.

    Args:
        analyses: List of CoherenceAnalysis objects from Phase 3
        output_path: Path to write JSON file
        genome_id: Optional genome identifier for metadata
    """
    import json

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    output_data = {
        "metadata": {
            "genome_id": genome_id,
            "n_eves": len(analyses),
            "version": "1.0",
        },
        "eves": [analysis.to_dict() for analysis in analyses],
        "summary": {
            "total_coherence_scores": [a.coherence_score for a in analyses],
            "mean_coherence": (
                sum(a.coherence_score for a in analyses) / len(analyses)
                if analyses else 0.0
            ),
            "confidence_levels": {
                "high": sum(1 for a in analyses if a.confidence_level == "high"),
                "medium": sum(1 for a in analyses if a.confidence_level == "medium"),
                "low": sum(1 for a in analyses if a.confidence_level == "low"),
                "reject": sum(1 for a in analyses if a.confidence_level == "reject"),
            },
        },
    }

    with open(output_path, "w") as f:
        json.dump(output_data, f, indent=2)

    logger.info(f"Wrote evidence graph JSON: {output_path}")
