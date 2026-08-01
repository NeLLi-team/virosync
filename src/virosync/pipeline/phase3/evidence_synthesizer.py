"""
Evidence Synthesizer for EVE Verification.

Implements the Gated Escalation Logic that routes predictions through
validation pathways based on confidence levels:

- High confidence (P > 0.95): Accept directly
- Low confidence (P < 0.60): Reject as noise
- Ambiguous (0.60 ≤ P ≤ 0.95): Send to tie-breaker modules

Tie-breaker modules include:
- Structural homology (Boltz + FoldSeek, optional)
- Evidence Correlation Graph (Coherence Score)
- Phylogenetic validation (GVClass + Diamond BLASTp)
"""

from __future__ import annotations

import csv
import logging
from collections import Counter
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING, Optional

from virosync.output_contract import canonical_family
from virosync.ablation import AblationID
from virosync.config import get_config
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.phase1.viral_markers import (
    CRESS_MARKER_MODELS,
    base_marker_gene_id,
    is_cress_specific_top1_marker,
    is_identity_qualified_cress_marker,
)
from virosync.pipeline.host_signatures import (
    HostSignatureModel,
    host_signature_density_evalue_weighted,
)
from virosync.utils.path_safety import require_strict_child, safe_filename_component

from .evidence_graph import (
    CoherenceAnalysis,
    analyze_eve_coherence,
)
from .phylogenetic_validation import (
    PhylogeneticValidator,
    PhylogeneticValidationResult,
)

if TYPE_CHECKING:
    from .structural_homology import BoltzFoldSeekAnalyzer, StructuralHomologyResult
    from .tmvec_database import TMVecDatabaseSearch

logger = logging.getLogger(__name__)


MARKER_CATEGORY_KEYWORDS = {
    "capsid": ["capsid", "mcp", "major capsid", "minor capsid", "coat"],
    "portal": ["portal"],
    "terminase": ["terminase"],
    "packaging_atpase": ["packaging atpase", "a32", "packaging", "atpase"],
    "penton": ["penton"],
    "triplex": ["triplex"],
    "protease": ["protease", "peptidase"],
    "polymerase": ["polymerase", "polb", "dna polymerase"],
    "helicase": ["helicase", "primase", "d5"],
    "transcription": ["vltf", "transcription factor", "rna polymerase", "rnap"],
}

MARKER_CATEGORY_GROUPS = {
    "structural": {"capsid", "penton", "triplex", "portal"},
    "packaging": {"terminase", "packaging_atpase"},
    "replication": {"polymerase", "helicase"},
    "transcription": {"transcription"},
    "protease": {"protease"},
}

from virosync.pipeline.phase3.mcp_detection import (
    is_mcp_gene,
)

VP_CORE_MARKER_PREFIXES = {
    "mcp": ("vp_mcp",),
    "penton": ("vp_penton",),
    "atpase": ("vp_atpase",),
    "protease": ("vp_pro",),
}

PLV_CORE_MARKER_PREFIXES = {
    "mcp": ("plv_mcp",),
    "pc": ("plv_pc",),
}

def _cress_gene_support(hallmark_hits: list) -> tuple[set[str], set[str]]:
    """Return identity-qualified and CRESS-specific-top-hit gene IDs."""

    qualified: set[str] = set()
    specific_top1: set[str] = set()
    for hit in hallmark_hits:
        if not is_identity_qualified_cress_marker(hit):
            continue
        if isinstance(hit, dict):
            porf_id = str(hit.get("porf_id") or hit.get("query_name") or "")
        else:
            porf_id = str(
                getattr(hit, "porf_id", None)
                or getattr(hit, "query_porf", "")
                or ""
            )
        if porf_id:
            gene_id = base_marker_gene_id(porf_id)
            qualified.add(gene_id)
            if is_cress_specific_top1_marker(hit):
                specific_top1.add(gene_id)
    return qualified, specific_top1


def _marker_totals_from_annotation(annotation_index: dict[str, dict]) -> dict[str, int]:
    """
    Count total unique marker categories for completeness calculation.

    Groups similar HMM variants (e.g., PLV_MCP_1-10 count as one "mcp" category)
    to match the per-protein deduplication logic in compute_marker_completeness.
    """
    ncldv_total = 0
    mirus_total = 0
    plv_categories: set[str] = set()

    for name, info in annotation_index.items():
        name_lower = name.lower()
        family = (info.get("family") or "").upper()

        if name_lower.startswith("gvogm") and family == "NCLDV":
            ncldv_total += 1

        if name_lower.startswith("mirus_") and family == "MIRUS":
            mirus_total += 1

        if name_lower.startswith("plv_"):
            # Group PLV markers by functional category
            category = None
            for label, prefixes in PLV_CORE_MARKER_PREFIXES.items():
                if name_lower.startswith(prefixes):
                    category = f"plv_{label}"
                    break
            # If not grouped, treat as unique marker
            if not category:
                category = name_lower
            plv_categories.add(category)

    return {
        "vp_total": 4,  # mcp, penton, atpase, protease (from VP_CORE_MARKER_PREFIXES)
        "ncldv_total": 9,
        "mirus_total": mirus_total,
        "plv_total": len(plv_categories),
    }


def _normalize_text(text: str) -> str:
    return text.lower().strip()


def _categorize_text(text: str) -> set[str]:
    text = _normalize_text(text)
    categories = set()
    for category, keywords in MARKER_CATEGORY_KEYWORDS.items():
        if any(k in text for k in keywords):
            categories.add(category)
    return categories


def _infer_family_from_model(model_name: str, source: str, description: str) -> str:
    name = _normalize_text(model_name)
    source = _normalize_text(source)
    desc = _normalize_text(description)
    if name.startswith("vp_") or "virophage" in desc or "virophage" in source:
        return "PPV"
    if name.startswith("plv_") or "polinton-like" in desc or "plv" in source:
        return "PPV"
    if name.startswith("mirus_") or "mirus" in source or "mirus" in desc:
        return "MIRUS"
    if name.startswith("gvogm") or name.startswith("gamadvirus") or name in {"og1352", "og484"}:
        return "NCLDV"
    if "ncldv" in desc or "ncldv" in source:
        return "NCLDV"
    return "UNKNOWN"


def load_marker_annotation_index(annotation_path: Optional[Path]) -> dict[str, dict]:
    if not annotation_path or not Path(annotation_path).exists():
        return {}
    index: dict[str, dict] = {}
    with open(annotation_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            model_name = (row.get("model_name") or "").strip()
            if not model_name:
                continue
            source = row.get("source") or ""
            description = row.get("description") or ""
            majority_annotation = row.get("majority_annotation") or ""
            text = " ".join([model_name, description, majority_annotation])
            categories = _categorize_text(text)
            family = _infer_family_from_model(model_name, source, description)
            # capscan MCP markers carry their Bellas group in the description
            # ("capscan <group> Major Capsid Protein"); surface it for likely_group.
            capscan_group = ""
            if "_caps_" in model_name.lower() and description.startswith("capscan "):
                capscan_group = description[len("capscan "):].split(" Major Capsid Protein")[0].strip()
            index[model_name.lower()] = {
                "categories": categories,
                "family": family,
                "description": description,
                "majority_annotation": majority_annotation,
                "capscan_group": capscan_group,
            }
    return index


def summarize_marker_hits(
    hallmark_genes: list[str],
    annotation_index: Optional[dict[str, dict]] = None,
) -> dict:
    annotation_index = annotation_index or {}
    categories: set[str] = set()
    families: Counter[str] = Counter()
    has_mcp = False

    for gene in hallmark_genes:
        gene_key = gene.lower()
        annotation = annotation_index.get(gene_key)
        if annotation:
            categories.update(annotation.get("categories", set()))
            family = annotation.get("family", "UNKNOWN")
            if family != "UNKNOWN":
                families[family] += 1
            if (
                "capsid" in annotation.get("categories", set())
                and gene.upper() not in CRESS_MARKER_MODELS
            ):
                has_mcp = True
            continue

        # Fallback: infer from gene name only
        categories.update(_categorize_text(gene))
        if gene_key.startswith(("vp_", "plv_")):
            # vp_ and plv_ markers are both Preplasmiviricota; the subgroup is
            # kept separately in ppv_subtype.
            families["PPV"] += 1
        elif gene_key.startswith("mirus_"):
            families["MIRUS"] += 1
        elif gene_key.startswith(("gvogm", "gamadvirus")) or gene_key in {"og1352", "og484"}:
            families["NCLDV"] += 1
        if is_mcp_gene(gene_key) or "capsid" in gene_key:
            has_mcp = True

    group_hits = 0
    for group, group_categories in MARKER_CATEGORY_GROUPS.items():
        if categories & group_categories:
            group_hits += 1

    complement_score = min(group_hits / 4.0, 1.0) if group_hits > 0 else 0.0
    dominant_family = "UNKNOWN"
    dominant_fraction = 0.0
    if families:
        dominant_family, dominant_count = families.most_common(1)[0]
        dominant_fraction = dominant_count / max(1, sum(families.values()))

    return {
        "categories": sorted(categories),
        "families": sorted(families.keys()),
        "family_counts": dict(families),
        "complement_score": complement_score,
        "dominant_family": dominant_family,
        "dominant_fraction": dominant_fraction,
        "has_mcp": has_mcp,
    }


def compute_marker_completeness(
    hallmark_genes: list[str],
    annotation_index: Optional[dict[str, dict]] = None,
    hallmark_hits: Optional[list[dict]] = None,
) -> dict[str, object]:
    """
    Compute marker completeness with per-protein deduplication.

    When a single protein hits multiple similar HMM models (e.g., PLV_MCP_1,
    PLV_MCP_2, PLV_MCP_9), only count the highest-scoring hit to avoid
    artificially inflating completeness scores.
    """
    annotation_index = annotation_index or {}
    totals = _marker_totals_from_annotation(annotation_index)

    # Build mapping from (protein_id, marker_category) -> best hit
    # This deduplicates hits from the same protein to similar HMM variants
    protein_markers: dict[tuple[str, str], dict] = {}

    if hallmark_hits:
        for i, hit in enumerate(hallmark_hits):
            # Extract protein ID and HMM model name
            if isinstance(hit, dict):
                porf_id = hit.get("porf_id") or hit.get("query_name") or f"unknown_{i}"
                gene = hit.get("hallmark_gene", "")
                score = hit.get("score", 0.0)
            else:
                porf_id = getattr(hit, "porf_id", None) or getattr(hit, "query_name", f"unknown_{i}")
                gene = getattr(hit, "hallmark_gene", "")
                score = getattr(hit, "score", 0.0)

            if not gene:
                continue

            key = gene.lower()

            # Determine marker category (group similar HMM variants)
            category = None
            family = None

            # VP markers: group by function
            for label, prefixes in VP_CORE_MARKER_PREFIXES.items():
                if key.startswith(prefixes):
                    category = f"vp_{label}"
                    family = "VP"
                    break

            # PLV markers: group by function
            if not category:
                for label, prefixes in PLV_CORE_MARKER_PREFIXES.items():
                    if key.startswith(prefixes):
                        category = f"plv_{label}"
                        family = "PLV"
                        break

            # NCLDV markers: each is unique
            if not category and key.startswith("gvogm"):
                category = key
                family = "NCLDV"

            # Mirus markers: each is unique
            if not category and key.startswith("mirus_"):
                category = key
                family = "MIRUS"

            # Other PLV markers without specific grouping
            if not category and key.startswith("plv_"):
                category = key
                family = "PLV"

            if category:
                marker_key = (porf_id, category)
                # Keep only highest-scoring hit per protein-category pair
                if marker_key not in protein_markers or score > protein_markers[marker_key]["score"]:
                    protein_markers[marker_key] = {
                        "gene": gene,
                        "category": category,
                        "family": family,
                        "score": score,
                    }

    # If no hallmark_hits provided, fall back to old behavior (no deduplication)
    if not protein_markers:
        vp_hits: set[str] = set()
        plv_hits: set[str] = set()
        ncldv_hits: set[str] = set()
        mirus_hits: set[str] = set()

        for gene in hallmark_genes:
            key = gene.lower()
            for label, prefixes in VP_CORE_MARKER_PREFIXES.items():
                if key.startswith(prefixes):
                    vp_hits.add(label)
            for label, prefixes in PLV_CORE_MARKER_PREFIXES.items():
                if key.startswith(prefixes):
                    plv_hits.add(label)
            if key.startswith("gvogm"):
                ncldv_hits.add(key)
            if key.startswith("mirus_"):
                mirus_hits.add(key)
            if key.startswith("plv_") and not any(key.startswith(p) for label, prefixes in PLV_CORE_MARKER_PREFIXES.items() for p in prefixes):
                plv_hits.add(key)
    else:
        # Count unique categories from deduplicated hits
        vp_hits: set[str] = set()
        plv_hits: set[str] = set()
        ncldv_hits: set[str] = set()
        mirus_hits: set[str] = set()

        for (porf_id, category), info in protein_markers.items():
            family = info["family"]
            if family == "VP":
                # Extract VP marker label (e.g., "vp_mcp" -> "mcp")
                vp_hits.add(category.replace("vp_", ""))
            elif family == "PLV":
                # For PLV, count grouped markers (e.g., "plv_mcp") and ungrouped (e.g., "plv_pc_054")
                if category.startswith("plv_"):
                    # For grouped markers like "plv_mcp", use the base name
                    for label, prefixes in PLV_CORE_MARKER_PREFIXES.items():
                        if category == f"plv_{label}":
                            plv_hits.add(category)
                            break
                    else:
                        # Ungrouped PLV markers: use full name
                        plv_hits.add(category)
                else:
                    plv_hits.add(category)
            elif family == "NCLDV":
                ncldv_hits.add(category)
            elif family == "MIRUS":
                mirus_hits.add(category)

    vp_total = totals["vp_total"]
    ncldv_total = totals["ncldv_total"]
    mirus_total = totals["mirus_total"]
    plv_total = totals["plv_total"]
    ppv_total = vp_total + plv_total
    ppv_hits = len(vp_hits) + len(plv_hits)

    vp_ratio = (len(vp_hits) / vp_total) if vp_total else 0.0
    ncldv_ratio = (len(ncldv_hits) / ncldv_total) if ncldv_total else 0.0
    mirus_ratio = (len(mirus_hits) / mirus_total) if mirus_total else 0.0
    ppv_ratio = (ppv_hits / ppv_total) if ppv_total else 0.0

    return {
        "vp_completeness": f"{len(vp_hits)}/{vp_total}",
        "vp_completeness_ratio": vp_ratio,
        "ncldv_completeness": f"{len(ncldv_hits)}/{ncldv_total}",
        "ncldv_completeness_ratio": ncldv_ratio,
        "mirus_completeness": f"{len(mirus_hits)}/{mirus_total}",
        "mirus_completeness_ratio": mirus_ratio,
        "ppv_completeness": f"{ppv_hits}/{ppv_total}",
        "ppv_completeness_ratio": ppv_ratio,
    }


def infer_ppv_subtype(hallmark_genes: list[str]) -> str:
    """Return VP or PLV only for unambiguous subtype-specific marker evidence."""
    has_vp = False
    has_plv = False
    for gene in hallmark_genes:
        key = gene.lower()
        for label, prefixes in VP_CORE_MARKER_PREFIXES.items():
            if label != "atpase" and key.startswith(prefixes):
                has_vp = True
        if key.startswith(PLV_CORE_MARKER_PREFIXES["mcp"]):
            has_plv = True

    if has_vp == has_plv:
        return ""
    if has_vp:
        return "VP"
    if has_plv:
        return "PLV"
    return ""


def infer_likely_family(result: "VerificationResult") -> str:
    classification = result.region_classification or ""
    if canonical_family(classification) in {"NCLDV", "MIRUS", "PPV", "CRESS"}:
        return canonical_family(classification)
    if classification == "MIXED":
        marker_family = canonical_family(result.marker_dominant_family)
        if marker_family in {"NCLDV", "MIRUS", "PPV"} and result.marker_dominant_fraction >= 0.70:
            return marker_family
        return classification

    # ppv_subtype keeps the VP/PLV subgroup; the family label is the lineage.
    if result.ppv_subtype:
        return canonical_family(result.ppv_subtype)

    marker_family = canonical_family(result.marker_dominant_family)
    if marker_family in {"NCLDV", "MIRUS", "PPV", "CRESS"}:
        return marker_family

    tax_family = canonical_family(result.gene_taxonomy_dominant_family)
    if tax_family in {"NCLDV", "MIRUS", "PPV", "CRESS"}:
        return tax_family

    return "UNKNOWN"


def infer_likely_group(
    scored_hallmarks: list[tuple[str, float]],
    annotation_index: Optional[dict[str, dict]] = None,
) -> str:
    """Best-hit capscan group (Trimcap/PgVV/Alpenseevirus/...) for an EVE's hallmark
    MCP markers. Returns the Bellas&Sommaruga 2026 group of the highest-scoring hallmark
    marker that carries a capscan group, resolving below the coarse PLV class. The shared
    double-jelly-roll MCP means one region can hit several group profiles, so the
    strongest hit is the reliable label (a majority count fragments across cross-hits).
    Empty when no hallmark marker carries a capscan group.

    Args:
        scored_hallmarks: ``(hallmark_gene_name, hmm_score)`` pairs for the EVE.
        annotation_index: marker-name -> annotation dict, keyed lowercase.
    """
    if not scored_hallmarks or not annotation_index:
        return ""
    best_group = ""
    best_score = float("-inf")
    for gene, score in scored_hallmarks:
        info = annotation_index.get((gene or "").lower())
        if info and info.get("capscan_group") and score > best_score:
            best_score = score
            best_group = info["capscan_group"]
    return best_group


def _confidence_tier_for_score(score: float, *, high: float, low: float) -> str:
    """Return the confidence tier for a score without mutating a result."""
    if score >= high:
        return "HIGH"
    if score >= low:
        return "MEDIUM"
    return "LOW"


def assign_confidence_tier(result: "VerificationResult", high: float = 0.7, low: float = 0.2) -> str:
    """Assign confidence tier based on final confidence score.

    Tiers:
      - HIGH: confidence >= high - high probability of true EVE
      - MEDIUM: low <= confidence < high - moderate evidence, needs validation
      - LOW: confidence < low - weak evidence, likely false positive

    Args:
        result: VerificationResult with final_confidence set
        high: Threshold for HIGH tier (default 0.8)
        low: Threshold for MEDIUM tier (default 0.4)

    Returns:
        Tier string: "HIGH", "MEDIUM", or "LOW"
    """
    return _confidence_tier_for_score(result.final_confidence, high=high, low=low)


def compute_marker_score(result: "VerificationResult") -> float:
    """Compute marker score from hallmark evidence."""
    # seed_sources records provenance, not marker identity. The authoritative
    # has_mcp flag is populated in _process_hallmark_hits via is_mcp_gene.
    has_mcp = result.has_mcp
    total_markers = result.hallmark_count
    diversity_score = min(result.hallmark_diversity / 5.0, 1.0) if total_markers > 0 else 0.0

    if result.marker_category_hits:
        return (
            0.45 * (1.0 if has_mcp else 0.0) +
            0.35 * result.marker_complement_score +
            0.20 * diversity_score
        )
    return (
        0.50 * (1.0 if has_mcp else 0.0) +
        0.30 * diversity_score +
        0.20 * min(total_markers / 8.0, 1.0)
    )


def compute_family_consistency_score(result: "VerificationResult") -> float:
    """Small bonus when marker, taxonomy, and interproscan evidence agree on family."""
    votes: Counter[str] = Counter()
    # Canonicalise before voting: a region labelled PPV and a marker family
    # labelled VP are the same lineage agreeing, not two families disagreeing.
    region_family = canonical_family(result.region_classification)
    taxonomy_family = canonical_family(result.gene_taxonomy_dominant_family)
    marker_family_vote = canonical_family(result.marker_dominant_family)
    if region_family in {"NCLDV", "MIRUS", "PPV"}:
        votes[region_family] += 1
    if taxonomy_family in {"NCLDV", "MIRUS", "PPV"}:
        if result.gene_taxonomy_dominant_fraction >= 0.25:
            votes[taxonomy_family] += 1
    if marker_family_vote in {"NCLDV", "MIRUS", "PPV"}:
        if result.marker_complement_score > 0:
            votes[marker_family_vote] += 1
    if result.interproscan_family_hits:
        # Deduplicate after canonicalising: legacy ["VP", "PLV"] is one lineage,
        # and would otherwise look like two families and suppress the vote.
        family_hits = {
            canonical_family(f)
            for f in result.interproscan_family_hits
            if canonical_family(f) in {"NCLDV", "MIRUS", "PPV"}
        }
        if len(family_hits) == 1:
            votes[next(iter(family_hits))] += 1
    if not votes:
        return 0.0
    top = max(votes.values())
    if top >= 3:
        return 0.08
    if top == 2:
        return 0.04
    return 0.0


def should_accept_mcp_override(result: "VerificationResult") -> bool:
    """Return True when MCP evidence is present for override acceptance."""
    if result.has_mcp and result.hallmark_count >= 1:
        return True
    if result.interproscan_category_hits and "capsid" in result.interproscan_category_hits:
        return True
    for keyword in result.interproscan_keyword_hits:
        key = keyword.lower()
        if "capsid" in key or "mcp" in key:
            return True
    return False


class VerificationStatus(Enum):
    """Confidence tier status of EVE verification."""

    HIGH_CONFIDENCE = "high_confidence"
    MEDIUM_CONFIDENCE = "medium_confidence"
    AMBIGUOUS = "ambiguous"
    LOW_CONFIDENCE_INITIAL = "low_confidence_initial"
    LOW_CONFIDENCE_TIEBREAKER = "low_confidence_tiebreaker"


@dataclass(frozen=True, slots=True)
class CompositionAblationEffect:
    """Per-candidate A5 counterfactual and aggregate-ready counters."""

    opportunities: int = 0
    interventions: int = 0
    changed: int = 0
    composition_score: float = 0.0
    reference_confidence: Optional[float] = None
    selected_confidence: Optional[float] = None
    reference_tier: str = ""
    selected_tier: str = ""

    def to_dict(self) -> dict[str, object]:
        """Return deterministic JSON-compatible metadata."""
        return {
            "opportunities": self.opportunities,
            "interventions": self.interventions,
            "changed": self.changed,
            "composition_score": self.composition_score,
            "reference_confidence": self.reference_confidence,
            "selected_confidence": self.selected_confidence,
            "reference_tier": self.reference_tier,
            "selected_tier": self.selected_tier,
        }


def evaluate_composition_ablation_effect(
    *,
    ablation_id: AblationID,
    composition_score: float,
    reference_confidence: float,
    selected_confidence: float,
    high_threshold: float,
    low_threshold: float,
    reference_tier: Optional[str] = None,
    selected_tier: Optional[str] = None,
) -> CompositionAblationEffect:
    """Compare A5 with A0 without changing a ``VerificationResult``.

    An A5 candidate with a nonzero active composition score is one opportunity
    and one applied intervention. ``changed`` records a difference in its final
    score or tier. Zero-composition candidates and other arms return this
    counter group's canonical zero value.
    """
    if not isinstance(ablation_id, AblationID):
        raise TypeError("ablation_id must be an AblationID")
    if ablation_id is not AblationID.A5:
        return CompositionAblationEffect()

    resolved_reference_tier = reference_tier or _confidence_tier_for_score(
        reference_confidence,
        high=high_threshold,
        low=low_threshold,
    )
    resolved_selected_tier = selected_tier or _confidence_tier_for_score(
        selected_confidence,
        high=high_threshold,
        low=low_threshold,
    )
    opportunity = int(composition_score > 0.0)
    changed = int(
        opportunity > 0
        and (
            reference_confidence != selected_confidence
            or resolved_reference_tier != resolved_selected_tier
        )
    )
    return CompositionAblationEffect(
        opportunities=opportunity,
        interventions=opportunity,
        changed=changed,
        composition_score=float(composition_score),
        reference_confidence=float(reference_confidence),
        selected_confidence=float(selected_confidence),
        reference_tier=resolved_reference_tier,
        selected_tier=resolved_selected_tier,
    )


@dataclass
class VerificationResult:
    """
    Complete verification result for an EVE candidate.

    Contains the verification status, all evidence, and final scores.
    """

    # Identity
    eve_id: str
    scaffold: str
    start: int
    end: int

    # Status
    status: VerificationStatus = VerificationStatus.AMBIGUOUS
    final_confidence: float = 0.0
    ablation_id: AblationID = AblationID.A0
    composition_ablation_effect: CompositionAblationEffect = field(
        default_factory=CompositionAblationEffect
    )

    # Input evidence
    crf_confidence: float = 0.0
    crf_posterior: float = 0.0

    # Tie-breaker results
    coherence_analysis: Optional[CoherenceAnalysis] = None
    structural_results: list[StructuralHomologyResult] = field(default_factory=list)
    phylogenetic_result: Optional[PhylogeneticValidationResult] = None
    tmvec_all_proteins: list[dict] = field(default_factory=list)

    # Aggregated scores
    coherence_score: float = 0.0
    structural_score: float = 0.0
    phylogenetic_score: float = 0.0

    # Inspectable breakdown of how ``final_confidence`` was assembled.
    # Populated by :func:`calculate_eve_confidence` so that reviewers can
    # see which components drove the score and detect double-counting
    # across marker / gene-taxonomy / seed-cluster / family-consistency
    # channels (same Diamond hit can feed several of these).
    score_components: dict = field(default_factory=dict)

    # Evidence summary
    hallmark_count: int = 0
    hallmark_diversity: int = 0
    has_virus_specific_marker: bool = False
    has_structural_support: bool = False
    has_phylogenetic_support: bool = False
    hallmark_genes: list[str] = field(default_factory=list)
    marker_category_hits: list[str] = field(default_factory=list)
    marker_family_hits: list[str] = field(default_factory=list)
    marker_complement_score: float = 0.0
    marker_dominant_family: str = "UNKNOWN"
    marker_dominant_fraction: float = 0.0
    family_consistency_score: float = 0.0
    vp_completeness: str = "0/4"
    vp_completeness_ratio: float = 0.0
    ppv_completeness: str = "0/0"
    ppv_completeness_ratio: float = 0.0
    ncldv_completeness: str = "0/0"
    ncldv_completeness_ratio: float = 0.0
    mirus_completeness: str = "0/0"
    mirus_completeness_ratio: float = 0.0

    # Phylogenetic details
    gvclass_domain: str = ""
    gvclass_percent: float = 0.0
    diamond_domain: str = ""
    diamond_percent: float = 0.0
    has_mcp: bool = False
    mcp_gene_ids: list[str] = field(default_factory=list)
    tier1_bypassed_marker_ids: list[str] = field(default_factory=list)
    tier1_bypassed_marker_models: list[str] = field(default_factory=list)
    is_chimeric: bool = False

    # Taxonomy (if determinable)
    predicted_taxonomy: str = ""
    taxonomy_confidence: float = 0.0

    # Region classification from Phase 1 seed markers (NCLDV, VP, PLV, MIRUS, MIXED, UNKNOWN)
    region_classification: str = ""
    region_classification_ncldv_markers: int = 0
    region_classification_vp_plv_markers: int = 0
    region_classification_mirus_markers: int = 0
    ppv_subtype: str = ""
    likely_family: str = "UNKNOWN"
    likely_group: str = ""  # capscan Bellas2026 group (Trimcap/PgVV/...) below the class
    confidence_tier: str = "LOW"  # HIGH, MEDIUM, or LOW
    candidate_start: Optional[int] = None
    candidate_end: Optional[int] = None

    # Contig-edge detection (partial EVE flag - no penalty, just informational)
    partial_eve: bool = False  # True if EVE is at contig boundary
    partial_eve_at_start: bool = False  # EVE starts near contig start
    partial_eve_at_end: bool = False  # EVE ends near contig end
    scaffold_length: int = 0  # Length of scaffold for context

    # Within-genome clustering (similar EVEs provide additional evidence)
    cluster_id: int = -1  # -1 means singleton
    cluster_size: int = 1  # Number of similar EVEs in cluster
    similar_eve_count: int = 0  # Number of EVEs with ANI >= 80%
    max_cluster_ani: float = 0.0  # Maximum ANI to any similar EVE
    clustering_bonus: float = 0.0  # Confidence bonus from clustering
    candidate_length: int = 0
    candidate_reduction_bp: int = 0
    candidate_reduction_reason: str = ""
    candidate_common_euk_taxonomy: str = ""
    phase1_host_signature_host_prefixes: str = ""
    phase1_host_signature_top_tokens: str = ""
    phase1_host_signature_token_count: int = 0

    # Per-gene taxonomy summary from Phase 2b Diamond
    gene_taxonomy_total: int = 0  # Interior genes only
    gene_taxonomy_ncldv_top10: int = 0
    gene_taxonomy_mirus_top10: int = 0
    gene_taxonomy_vp_plv_top10: int = 0
    gene_taxonomy_phage_top10: int = 0
    gene_taxonomy_viral_top10: int = 0  # Interior viral genes only (from boundary_diamond)
    gene_taxonomy_cellular: int = 0
    gene_taxonomy_unknown: int = 0
    gene_taxonomy_has_ncldv_mirus: bool = False
    gene_taxonomy_has_vp_plv: bool = False
    gene_taxonomy_dominant_family: str = "UNKNOWN"
    gene_taxonomy_dominant_fraction: float = 0.0
    gene_taxonomy_records: list[dict] = field(default_factory=list)

    # NEW: Flanking gene tracking (for taxonomy expansion fix)
    gene_taxonomy_total_with_flanking: int = 0  # Interior + flanking
    gene_taxonomy_flanking_count: int = 0  # Number of flanking genes analyzed
    gene_taxonomy_viral_interior: int = 0  # Viral genes inside boundary
    gene_taxonomy_viral_flanking: int = 0  # Viral genes in flanking regions
    gene_taxonomy_ncldv_mirus_interior: int = 0  # NCLDV/MIRUS inside boundary
    gene_taxonomy_ncldv_mirus_flanking: int = 0  # NCLDV/MIRUS in flanking
    gene_taxonomy_vp_plv_interior: int = 0  # VP/PLV inside boundary
    gene_taxonomy_vp_plv_flanking: int = 0  # VP/PLV in flanking

    # Seed metadata (passed through from Phase 1 via Phase 2)
    seed_sources: list[str] = field(default_factory=list)
    seed_confidence: str = ""
    seed_hhg_score: float = 0.0
    seed_novelty_score: float = 0.0
    seed_compositional_score: float = 0.0

    # Composition evidence features
    kfd: float = 0.0
    gc_deviation: float = 0.0
    cub_deviation: float = 0.0

    # Gene-level taxonomy (Step 9 - final per-gene Diamond)
    genes_with_ncldv_mirus_top10: int = 0
    genes_with_vp_plv_top10: int = 0
    genes_with_high_pident_euk: int = 0
    high_confidence_euk_genes: int = 0  # Top-3 all EUK >=70% pident, no viral in top-10
    gene_count: int = 0
    host_signature_gene_count: int = 0
    host_signature_fraction: float = 0.0
    host_signature_weighted_mean: float = 0.0
    host_signature_evalue_weighted: float = 0.0  # ISSUE 4: E-value weighted penalty

    # Taxonomy distribution analysis (EVE vs host gene discrimination)
    taxonomy_distribution_viral_score: float = 0.0  # Overall viral likelihood (0-1)
    taxonomy_distribution_diversity: float = 0.0  # Genus diversity in hits
    taxonomy_distribution_host_overlap: float = 0.0  # Overlap with host baseline
    taxonomy_distribution_non_euk_fraction: float = 0.0  # Fraction of non-EUK hits
    taxonomy_distribution_genes_analyzed: int = 0  # Genes with distribution analysis
    taxonomy_distribution_likely_viral_genes: int = 0  # Genes with viral_score >= 0.4
    taxonomy_distribution_likely_host_genes: int = 0  # Genes with host pattern
    taxonomy_distribution_baseline_source: str = "default"
    taxonomy_distribution_baseline_markers: int = 0
    taxonomy_distribution_baseline_genera: int = 0
    taxonomy_distribution_baseline_diversity: float = 0.0

    # InterProScan annotation summary
    interproscan_total_hits: int = 0
    interproscan_viral_hits: int = 0
    interproscan_keyword_hits: list[str] = field(default_factory=list)
    interproscan_score: float = 0.0
    interproscan_category_hits: list[str] = field(default_factory=list)
    interproscan_family_hits: list[str] = field(default_factory=list)
    interproscan_category_score: float = 0.0

    # NUMT detection (InterProScan-based mitochondrial markers)
    interproscan_numt_hits: int = 0
    interproscan_numt_markers: list[str] = field(default_factory=list)
    numt_flag: str = "NONE"  # NONE/DETECTED

    # Jelly roll (DJR/SJR) MCP classification
    jelly_roll_djr_count: int = 0  # Number of DJR-classified MCP proteins
    jelly_roll_sjr_count: int = 0  # Number of SJR-classified MCP proteins
    jelly_roll_total_mcp: int = 0  # Total MCP proteins classified
    jelly_roll_avg_confidence: float = 0.0  # Average classification confidence
    jelly_roll_confidence_bonus: float = 0.0  # Confidence bonus from validated DJR
    jelly_roll_mcp_proteins: list[dict] = field(default_factory=list)  # Per-protein results

    @property
    def is_high_or_medium_confidence(self) -> bool:
        """Whether the EVE meets minimum confidence threshold (MEDIUM or HIGH tier)."""
        return self.confidence_tier in {"HIGH", "MEDIUM"}

    @property
    def is_accepted(self) -> bool:
        """Legacy alias for HIGH/MEDIUM confidence, not the v2 output gate."""
        return self.is_high_or_medium_confidence

    @property
    def length(self) -> int:
        return self.end - self.start

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "eve_id": self.eve_id,
            "scaffold": self.scaffold,
            "start": self.start,
            "end": self.end,
            "length": self.length,
            "status": self.status.value,
            "final_confidence": self.final_confidence,
            "ablation_id": self.ablation_id.value,
            "composition_ablation_effect": self.composition_ablation_effect.to_dict(),
            "crf_confidence": self.crf_confidence,
            "crf_posterior": self.crf_posterior,
            "coherence_score": self.coherence_score,
            "structural_score": self.structural_score,
            "phylogenetic_score": self.phylogenetic_score,
            "score_components": dict(self.score_components) if self.score_components else {},
            "hallmark_count": self.hallmark_count,
            "hallmark_diversity": self.hallmark_diversity,
            "has_virus_specific_marker": self.has_virus_specific_marker,
            "has_structural_support": self.has_structural_support,
            "has_phylogenetic_support": self.has_phylogenetic_support,
            "hallmark_genes": self.hallmark_genes,
            "marker_category_hits": self.marker_category_hits,
            "marker_family_hits": self.marker_family_hits,
            "marker_complement_score": self.marker_complement_score,
            "marker_dominant_family": self.marker_dominant_family,
            "marker_dominant_fraction": self.marker_dominant_fraction,
            "family_consistency_score": self.family_consistency_score,
            "vp_completeness": self.vp_completeness,
            "vp_completeness_ratio": self.vp_completeness_ratio,
            "ppv_completeness": self.ppv_completeness,
            "ppv_completeness_ratio": self.ppv_completeness_ratio,
            "ncldv_completeness": self.ncldv_completeness,
            "ncldv_completeness_ratio": self.ncldv_completeness_ratio,
            "mirus_completeness": self.mirus_completeness,
            "mirus_completeness_ratio": self.mirus_completeness_ratio,
            "gvclass_domain": self.gvclass_domain,
            "gvclass_percent": self.gvclass_percent,
            "diamond_domain": self.diamond_domain,
            "diamond_percent": self.diamond_percent,
            "has_mcp": self.has_mcp,
            "mcp_gene_ids": self.mcp_gene_ids,
            "tier1_bypassed_marker_ids": self.tier1_bypassed_marker_ids,
            "tier1_bypassed_marker_models": self.tier1_bypassed_marker_models,
            "is_chimeric": self.is_chimeric,
            "predicted_taxonomy": self.predicted_taxonomy,
            "taxonomy_confidence": self.taxonomy_confidence,
            "gene_taxonomy_total": self.gene_taxonomy_total,
            "gene_taxonomy_ncldv_top10": self.gene_taxonomy_ncldv_top10,
            "gene_taxonomy_mirus_top10": self.gene_taxonomy_mirus_top10,
            "gene_taxonomy_vp_plv_top10": self.gene_taxonomy_vp_plv_top10,
            "gene_taxonomy_phage_top10": self.gene_taxonomy_phage_top10,
            "gene_taxonomy_viral_top10": self.gene_taxonomy_viral_top10,
            "gene_taxonomy_total_with_flanking": self.gene_taxonomy_total_with_flanking,
            "gene_taxonomy_flanking_count": self.gene_taxonomy_flanking_count,
            "gene_taxonomy_viral_interior": self.gene_taxonomy_viral_interior,
            "gene_taxonomy_viral_flanking": self.gene_taxonomy_viral_flanking,
            "gene_taxonomy_cellular": self.gene_taxonomy_cellular,
            "gene_taxonomy_unknown": self.gene_taxonomy_unknown,
            "gene_taxonomy_has_ncldv_mirus": self.gene_taxonomy_has_ncldv_mirus,
            "gene_taxonomy_has_vp_plv": self.gene_taxonomy_has_vp_plv,
            "gene_taxonomy_dominant_family": self.gene_taxonomy_dominant_family,
            "gene_taxonomy_dominant_fraction": self.gene_taxonomy_dominant_fraction,
            "ppv_subtype": self.ppv_subtype,
            "likely_family": self.likely_family,
            "likely_group": self.likely_group,
            "confidence_tier": self.confidence_tier,
            "candidate_start": self.candidate_start,
            "candidate_end": self.candidate_end,
            "candidate_length": self.candidate_length,
            "candidate_reduction_bp": self.candidate_reduction_bp,
            "candidate_reduction_reason": self.candidate_reduction_reason,
            "candidate_common_euk_taxonomy": self.candidate_common_euk_taxonomy,
            "phase1_host_signature_host_prefixes": self.phase1_host_signature_host_prefixes,
            "phase1_host_signature_top_tokens": self.phase1_host_signature_top_tokens,
            "phase1_host_signature_token_count": self.phase1_host_signature_token_count,
            # Seed metadata from Phase 1
            "seed_sources": self.seed_sources,
            "seed_confidence": self.seed_confidence,
            "seed_hhg_score": self.seed_hhg_score,
            "seed_novelty_score": self.seed_novelty_score,
            "seed_compositional_score": self.seed_compositional_score,
            # Composition features
            "kfd": self.kfd,
            "gc_deviation": self.gc_deviation,
            "cub_deviation": self.cub_deviation,
            # Gene-level taxonomy (Step 9)
            "genes_with_ncldv_mirus_top10": self.genes_with_ncldv_mirus_top10,
            "genes_with_vp_plv_top10": self.genes_with_vp_plv_top10,
            "genes_with_high_pident_euk": self.genes_with_high_pident_euk,
            "high_confidence_euk_genes": self.high_confidence_euk_genes,
            "gene_count": self.gene_count,
            "host_signature_gene_count": self.host_signature_gene_count,
            "host_signature_fraction": self.host_signature_fraction,
            "host_signature_weighted_mean": self.host_signature_weighted_mean,
            # Taxonomy distribution analysis
            "taxonomy_distribution_viral_score": self.taxonomy_distribution_viral_score,
            "taxonomy_distribution_diversity": self.taxonomy_distribution_diversity,
            "taxonomy_distribution_host_overlap": self.taxonomy_distribution_host_overlap,
            "taxonomy_distribution_non_euk_fraction": self.taxonomy_distribution_non_euk_fraction,
            "taxonomy_distribution_genes_analyzed": self.taxonomy_distribution_genes_analyzed,
            "taxonomy_distribution_likely_viral_genes": self.taxonomy_distribution_likely_viral_genes,
            "taxonomy_distribution_likely_host_genes": self.taxonomy_distribution_likely_host_genes,
            "taxonomy_distribution_baseline_source": self.taxonomy_distribution_baseline_source,
            "taxonomy_distribution_baseline_markers": self.taxonomy_distribution_baseline_markers,
            "taxonomy_distribution_baseline_genera": self.taxonomy_distribution_baseline_genera,
            "taxonomy_distribution_baseline_diversity": self.taxonomy_distribution_baseline_diversity,
            "interproscan_total_hits": self.interproscan_total_hits,
            "interproscan_viral_hits": self.interproscan_viral_hits,
            "interproscan_keyword_hits": self.interproscan_keyword_hits,
            "interproscan_score": self.interproscan_score,
            "interproscan_category_hits": self.interproscan_category_hits,
            "interproscan_family_hits": self.interproscan_family_hits,
            "interproscan_category_score": self.interproscan_category_score,
            "jelly_roll_djr_count": self.jelly_roll_djr_count,
            "jelly_roll_sjr_count": self.jelly_roll_sjr_count,
            "jelly_roll_total_mcp": self.jelly_roll_total_mcp,
            "jelly_roll_avg_confidence": self.jelly_roll_avg_confidence,
            "jelly_roll_confidence_bonus": self.jelly_roll_confidence_bonus,
            "jelly_roll_mcp_proteins": self.jelly_roll_mcp_proteins,
            "tmvec_per_protein": self.tmvec_all_proteins,
        }


@dataclass
class EvidenceSynthesizerConfig:
    """Configuration for evidence synthesis.

    All EVE candidates receive full analysis. The tier thresholds
    classify final confidence scores into HIGH, MEDIUM, or LOW tiers.
    Defaults are loaded from virosync.config.thresholds centralized configuration.
    """

    ablation_id: AblationID = AblationID.A0

    # Confidence tier thresholds (for output classification)
    # HIGH: confidence >= high_tier_threshold
    # MEDIUM: low_tier_threshold <= confidence < high_tier_threshold
    # LOW: confidence < low_tier_threshold
    high_tier_threshold: float = field(
        default_factory=lambda: get_config().evidence.high_tier_threshold
    )
    low_tier_threshold: float = field(
        default_factory=lambda: get_config().evidence.low_tier_threshold
    )
    use_crf_in_final_score: bool = False
    priority_marker_list: list[str] = field(default_factory=lambda: ["mcp"])
    marker_floor_priority_only: float = 0.55
    marker_floor_priority_plus_family: float = 0.70
    marker_floor_priority_multi_family: float = 0.80
    marker_family_bonus_per_family: float = 0.06
    marker_multi_family_bonus: float = 0.08

    # Tie-breaker weights (normalized based on available evidence)
    # These weights are used when all tie-breakers are enabled
    crf_weight: float = 0.25
    coherence_weight: float = 0.25
    structural_weight: float = 0.20
    phylogenetic_weight: float = 0.30  # GVClass + Diamond carry significant weight
    interproscan_weight: float = 0.10
    marker_weight: float = 0.10
    family_consistency_weight: float = 0.05

    # Optional host signature strings for host-like penalty
    euk_host_signatures: Optional[set[str]] = None
    # Weighted host signature model
    host_signature_model: Optional[dict] = None
    host_signature_score_threshold: float = 0.3
    host_prefixes: Optional[list[str]] = None
    host_label: str = "EUK"

    # Structural analysis settings (Boltz + FoldSeek)
    use_boltz: bool = False
    boltz_use_msa_server: bool = False
    boltz_min_seq_len: int = 100
    boltz_max_seq_len: int = 1000
    boltz_no_kernels: bool = True  # Use --no_kernels for safer Boltz execution
    boltz_mcp_only: bool = True
    structural_max_porfs: int = 50  # Max pORFs to analyze per EVE

    # TMVec database search (fast structural evidence)
    use_tmvec_database: bool = False
    tmvec_databases: Optional[list[str]] = None  # ["bfvd", "cath", "swissprot"]
    tmvec_database_dir: Optional[Path] = None
    tmvec_min_score: float = 0.5
    tmvec_require_gpu: bool = False

    # Phylogenetic validation settings
    use_phylogenetic_validation: bool = False
    gvclass_db: Optional[Path] = None
    diamond_db: Optional[Path] = None

    # InterProScan settings
    use_interproscan: bool = False
    interproscan_dir: Optional[Path] = None
    interproscan_keywords: Optional[list[str]] = None

    # Marker annotation index (for category scoring)
    marker_annotations_path: Optional[Path] = None

    # Phylogenetic rejection override
    # If phylogenetic validation strongly rejects viral, override boundary support.
    phylogenetic_rejection_override: bool = True
    phylogenetic_rejection_threshold: float = field(
        default_factory=lambda: get_config().evidence.phylogenetic_rejection
    )

    # GPU settings
    device: str = "cuda"
    threads: int = 8

    def __post_init__(self) -> None:
        if not isinstance(self.ablation_id, AblationID):
            raise TypeError("ablation_id must be an AblationID")


def load_jelly_roll_data(jelly_roll_path: Path) -> dict[str, list[dict]]:
    """
    Load jelly roll classification data from TSV file.

    Returns a dict mapping pORF base IDs (without domain suffixes) to classification records.

    Args:
        jelly_roll_path: Path to virosync_jelly_roll_proteins.tsv

    Returns:
        Dict mapping pORF ID to list of classification records
    """
    if not jelly_roll_path or not jelly_roll_path.exists():
        return {}

    results: dict[str, list[dict]] = {}
    try:
        with open(jelly_roll_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                # Support both historical and current output schemas:
                # - porf_id/classification
                # - protein_id/type
                porf_id = (row.get("porf_id") or row.get("protein_id") or "").strip()
                if not porf_id:
                    continue
                base_id = porf_id.split("|aa", 1)[0] if "|aa" in porf_id else porf_id

                confidence_raw = row.get("confidence", 0.0)
                try:
                    confidence = float(confidence_raw)
                except (TypeError, ValueError):
                    confidence = 0.0

                classification = (row.get("classification") or row.get("type") or "UNKNOWN").upper()
                record = {
                    "porf_id": porf_id,
                    "base_id": base_id,
                    "marker": row.get("marker", ""),
                    "classification": classification,
                    "confidence": confidence,
                    "evidence": row.get("evidence", ""),
                }
                if base_id not in results:
                    results[base_id] = []
                results[base_id].append(record)
        logger.info(
            "Loaded jelly roll data for %d proteins from %s",
            len(results),
            jelly_roll_path,
        )
    except Exception as e:
        logger.warning(f"Failed to load jelly roll data from {jelly_roll_path}: {e}")

    return results


def _has_priority_marker(result: VerificationResult, priority_markers: list[str]) -> bool:
    marker_set = {m.lower() for m in priority_markers}
    if result.has_mcp and "mcp" in marker_set:
        return True
    # seed_sources intentionally excluded: it only carries
    # {"hhg","novelty","compositional"} source tags, never priority-marker
    # tokens like "mcp"/"polb". See note above compute_marker_score.
    marker_tokens = (
        list(result.hallmark_genes)
        + list(result.interproscan_keyword_hits)
        + list(result.interproscan_category_hits)
    )
    # "mcp" goes through the canonical detector: a substring match promoted any
    # name merely containing the trigram ("ncmcp_pseudoprotein", "Baculovirus
    # mcp-like repeat domain") to the priority-marker confidence floor, which is
    # applied after the host-signature penalties. The remaining markers have no
    # canonical detector and stay on substring matching.
    mcp_is_priority = "mcp" in marker_set
    substring_markers = marker_set - {"mcp"}
    for token in marker_tokens:
        lower = str(token).lower()
        if mcp_is_priority and is_mcp_gene(lower):
            return True
        if any(marker in lower for marker in substring_markers):
            return True
    return False


def calculate_eve_confidence(
    result: VerificationResult,
    crf_confidence: float,
    tmvec_score: Optional[float] = None,
    use_crf_score: bool = False,
    priority_markers: Optional[list[str]] = None,
    marker_floor_priority_only: float = 0.55,
    marker_floor_priority_plus_family: float = 0.70,
    marker_floor_priority_multi_family: float = 0.80,
    marker_family_bonus_per_family: float = 0.06,
    marker_multi_family_bonus: float = 0.08,
    host_signature_score_threshold: float = 0.3,
    ablation_id: AblationID = AblationID.A0,
    high_tier_threshold: float = 0.7,
    low_tier_threshold: float = 0.2,
) -> float:
    """
    Calculate final EVE confidence score integrating all evidence types.

    Implements Step 10 confidence formula from PIPELINE_HMM_GATED_PLAN.md:

    Score Components:
    - Marker score (0.30): Based on marker types and counts
    - Composition score (0.20): Deviation from genome background
    - Gene taxonomy bonus (0.20): Viral (NCLDV/MIRUS/VP/PLV) genes in region
    - Structural score (0.15): TMVec structural similarity
    - Optional boundary confidence (0.15): legacy Phase 2 confidence signal

    Penalties:
    - High-identity EUK (max -0.30): Fraction of genes with >=70% EUK identity

    Args:
        result: VerificationResult with all evidence populated
        crf_confidence: Legacy Phase 2 boundary confidence
        tmvec_score: Optional TMVec structural score
        ablation_id: Mutually exclusive benchmark ablation arm

    Returns:
        Confidence score in [0, 1]
    """
    if not isinstance(ablation_id, AblationID):
        raise TypeError("ablation_id must be an AblationID")

    def _seed_cluster_score(res: VerificationResult) -> float:
        ncldv = res.region_classification_ncldv_markers
        vp_plv = res.region_classification_vp_plv_markers
        mirus = res.region_classification_mirus_markers
        if vp_plv >= 4 or ncldv >= 3 or mirus >= 3:
            return 1.0
        if max(ncldv, vp_plv, mirus) >= 2:
            return 0.5
        return 0.0
    # ─────────────────────────────────────────────────────────────
    # 1. Marker Score (0-1)
    # ─────────────────────────────────────────────────────────────
    marker_score = compute_marker_score(result)

    # ─────────────────────────────────────────────────────────────
    # 2. Composition Score (0-1)
    # ─────────────────────────────────────────────────────────────
    # Higher score = more distinct from host
    kfd_norm = min(result.kfd / 0.3, 1.0)
    gc_dev_norm = min(result.gc_deviation / 0.15, 1.0)
    cub_norm = min(result.cub_deviation / 0.5, 1.0)

    composition_score = 0.4 * kfd_norm + 0.3 * gc_dev_norm + 0.3 * cub_norm

    # ─────────────────────────────────────────────────────────────
    # 3. Gene Taxonomy Bonus (0-1)
    # ─────────────────────────────────────────────────────────────
    # Use interior genes only for viral_fraction: both numerator
    # (gene_taxonomy_viral_top10) and denominator (gene_count) are
    # interior-only counts.  Flanking genes are host context and
    # would unfairly dilute the viral fraction.
    total_genes_interior = result.gene_count  # Interior only

    if total_genes_interior > 0:
        viral_genes = result.gene_taxonomy_viral_top10 or result.genes_with_ncldv_mirus_top10

        viral_fraction = min(viral_genes / total_genes_interior, 1.0)

        dominant_fraction = result.gene_taxonomy_dominant_fraction or 0.0
        gene_bonus = min(1.0, 0.7 * viral_fraction + 0.3 * dominant_fraction)
    else:
        gene_bonus = 0.0

    # ─────────────────────────────────────────────────────────────
    # 4. Structural Score (0-1) - TMVec
    # ─────────────────────────────────────────────────────────────
    if tmvec_score is not None and tmvec_score > 0:
        # Normalize TM-score (0.5+ is significant)
        structural_score = min((tmvec_score - 0.3) / 0.4, 1.0)
        structural_score = max(0.0, structural_score)
    else:
        # Use existing structural score from result
        structural_score = result.structural_score

    # ─────────────────────────────────────────────────────────────
    # 5. Boundary Confidence (0-1; legacy field name)
    # ─────────────────────────────────────────────────────────────
    crf_score = crf_confidence

    # ─────────────────────────────────────────────────────────────
    # 6. High-Confidence EUK Penalty
    # ─────────────────────────────────────────────────────────────
    # A gene is "high-confidence EUK" when its top 3 Diamond hits are
    # ALL eukaryotic with >=70% sequence identity AND it has zero viral
    # hits in the entire top 10.  This catches clear host genes without
    # relying on the host signature model (which can be too strict).
    # Divergent EUK genes (common in EVEs) typically have <70% identity
    # or mixed viral/EUK top-10 hits, so they are NOT penalized.
    high_conf_euk = 0
    for rec in result.gene_taxonomy_records:
        is_flanking = rec.get("is_flanking", False) if isinstance(rec, dict) else getattr(rec, "is_flanking", False)
        has_hit = rec.get("has_hit", False) if isinstance(rec, dict) else getattr(rec, "has_hit", False)
        has_viral = rec.get("has_viral", False) if isinstance(rec, dict) else getattr(rec, "has_viral", False)
        if is_flanking or not has_hit or has_viral:
            continue
        prefixes = rec.get("top10_prefixes", []) if isinstance(rec, dict) else getattr(rec, "top10_prefixes", [])
        pidents = rec.get("top10_pidents", []) if isinstance(rec, dict) else getattr(rec, "top10_pidents", [])
        if not isinstance(prefixes, list) or not isinstance(pidents, list):
            continue
        top3_p = prefixes[:3]
        top3_id = pidents[:3]
        if len(top3_p) < 3 or len(top3_id) < 3:
            continue
        try:
            all_euk = all(str(p).upper().startswith("EUK") for p in top3_p)
            all_high_id = all(float(pid) >= 70.0 for pid in top3_id)
        except (TypeError, ValueError):
            continue
        if all_euk and all_high_id:
            high_conf_euk += 1
    if total_genes_interior > 0 and high_conf_euk > 0:
        euk_fraction = high_conf_euk / total_genes_interior
        euk_penalty = min(euk_fraction * 0.20, 0.12)
    else:
        euk_penalty = 0.0

    # Broad EUK penalty — secondary signal for overwhelmingly eukaryotic
    # regions with negligible viral evidence.  Fires when >=80% of interior
    # genes have high-identity eukaryotic top hits AND <10% have viral hits.
    # This catches FP regions where the strict high-confidence EUK penalty
    # (top-3 all-EUK >=70%, no viral in top-10) is too conservative.
    if total_genes_interior > 0:
        broad_euk_frac = result.genes_with_high_pident_euk / total_genes_interior
        if broad_euk_frac > 0.80 and viral_fraction < 0.10:
            broad_euk_penalty = min(broad_euk_frac * 0.12, 0.12)
            euk_penalty = max(euk_penalty, broad_euk_penalty)

    # Host signature penalty — fires when a significant fraction of
    # interior genes match the host taxonomy model (fraction >= threshold).
    # Real EVEs carry some host-like flanking genes, but their viral
    # core keeps the fraction low.  FP regions that are purely host
    # sequence have fraction well above the threshold.
    host_frac_threshold = host_signature_score_threshold
    if result.host_signature_gene_count > 0 and result.host_signature_fraction >= host_frac_threshold:
        host_excess = result.host_signature_fraction - host_frac_threshold
        host_signature_penalty = min(host_excess * 0.60, 0.18)
    else:
        host_signature_penalty = 0.0

    # ─────────────────────────────────────────────────────────────
    # 7. Taxonomy Distribution Modulation
    # ─────────────────────────────────────────────────────────────
    # High taxonomy distribution viral score indicates EVE despite EUK hits
    # Reduces EUK penalty when genes show divergent taxonomy (viral pattern)
    tax_dist_viral_score = result.taxonomy_distribution_viral_score
    if tax_dist_viral_score >= 0.4:
        # Reduce EUK penalty proportionally to viral score
        reduction_factor = min(tax_dist_viral_score, 0.8)  # Max 80% reduction
        euk_penalty *= (1.0 - reduction_factor)
        host_signature_penalty *= (1.0 - reduction_factor * 0.5)

    # ─────────────────────────────────────────────────────────────
    # InterProScan Bonus (0-1)
    # ─────────────────────────────────────────────────────────────
    interproscan_score = result.interproscan_score or 0.0
    seed_cluster_score = _seed_cluster_score(result)

    # Family consistency bonus (small additive)
    family_consistency_bonus = compute_family_consistency_score(result)

    # ─────────────────────────────────────────────────────────────
    # Size penalty: small regions without MCP are less reliable.
    # Boundary/composition evidence can inflate confidence for 1-2 gene fragments
    # because those signals are not strong region-level support by themselves.
    # ─────────────────────────────────────────────────────────────
    region_genes = result.gene_count or 0
    size_penalty = 0.0
    if not result.has_mcp:
        if region_genes <= 1:
            size_penalty = 0.15
        elif region_genes <= 2:
            size_penalty = 0.08
        elif region_genes <= 3:
            size_penalty = 0.04

    component_scores = {
        "crf": crf_score,
        "composition": composition_score,
        "marker": marker_score,
        "gene_taxonomy": gene_bonus,
        "interpro": interproscan_score,
        "seed_cluster": seed_cluster_score,
    }
    base_component_weights = {
        "crf": 0.18,
        "composition": 0.18,
        "marker": 0.24,
        "gene_taxonomy": 0.24,
        "interpro": 0.08,
        "seed_cluster": 0.08,
    }

    # Boundary-confidence weight reduction for small regions.
    if region_genes <= 2 and use_crf_score:
        base_component_weights["crf"] = 0.06  # 1/3 of normal weight
    elif region_genes <= 4 and use_crf_score:
        base_component_weights["crf"] = 0.12  # 2/3 of normal weight

    # NOTE: Genome-wide composition baseline was evaluated but disabled.
    # When every candidate shares similar composition (e.g., Aurantiochytrium
    # median ~0.71), any across-the-board reduction also eliminates true EVEs
    # that sit close to the MEDIUM threshold.  Future improvement: couple
    # composition normalization with stronger non-composition evidence to
    # provide margin before deflating.
    if not use_crf_score:
        base_component_weights["crf"] = 0.0

    common_bonus = 0.0
    reference_bonus = 0.0

    def _add_common_bonus(value: float) -> None:
        """Add a non-composition term in the same order to both score paths."""

        nonlocal common_bonus, reference_bonus
        common_bonus += value
        reference_bonus += value

    if marker_score >= 0.6 and gene_bonus >= 0.6:
        _add_common_bonus(0.04)
    if seed_cluster_score >= 1.0 and interproscan_score >= 0.5:
        _add_common_bonus(0.03)
    composition_bonus = 0.0
    if crf_score >= 0.9 and composition_score >= 0.6:
        composition_bonus = 0.03
        reference_bonus += composition_bonus
    if family_consistency_bonus > 0:
        _add_common_bonus(min(family_consistency_bonus, 0.05))
    if result.has_mcp:
        _add_common_bonus(0.05)
    completeness_ratio = max(
        result.vp_completeness_ratio,
        result.ppv_completeness_ratio,
        result.ncldv_completeness_ratio,
        result.mirus_completeness_ratio,
    )
    if completeness_ratio >= 0.9:
        _add_common_bonus(0.60)
    elif completeness_ratio >= 0.7:
        _add_common_bonus(0.40)
    elif completeness_ratio >= 0.5:
        _add_common_bonus(0.20)
    ncldv_core_hits = len({g.lower() for g in result.hallmark_genes if g.lower().startswith("gvogm")})
    if ncldv_core_hits >= 4:
        _add_common_bonus(0.20)
    elif ncldv_core_hits >= 3:
        _add_common_bonus(0.10)

    # Within-genome clustering bonus (similar EVEs provide additional evidence)
    if result.clustering_bonus > 0:
        _add_common_bonus(result.clustering_bonus)
        logger.debug(
            f"Clustering bonus for {result.eve_id}: +{result.clustering_bonus:.2f} "
            f"(cluster_size={result.cluster_size}, max_ani={result.max_cluster_ani:.1f}%)"
        )

    # Taxonomy distribution bonus (divergent taxonomy suggests viral origin)
    if tax_dist_viral_score >= 0.6:
        # High viral score = strong evidence of viral origin
        tax_bonus = min(0.15, (tax_dist_viral_score - 0.4) * 0.5)
        _add_common_bonus(tax_bonus)
        logger.debug(
            f"Taxonomy distribution bonus for {result.eve_id}: +{tax_bonus:.2f} "
            f"(viral_score={tax_dist_viral_score:.2f})"
        )

    # Strong bonus when the region diverges from sampled host-control taxonomy fingerprint.
    # Only apply when the baseline and per-gene distribution are actually available.
    if (
        result.taxonomy_distribution_genes_analyzed > 0
        and result.taxonomy_distribution_baseline_markers > 0
    ):
        host_overlap = max(0.0, min(1.0, result.taxonomy_distribution_host_overlap))
        host_divergence = 1.0 - host_overlap
        if host_divergence >= 0.5:
            host_divergence_bonus = min(0.20, 0.05 + (host_divergence - 0.5) * 0.30)
            _add_common_bonus(host_divergence_bonus)
            logger.debug(
                f"Host-divergence bonus for {result.eve_id}: +{host_divergence_bonus:.2f} "
                f"(host_overlap={host_overlap:.2f}, divergence={host_divergence:.2f})"
            )

    # Viral gene enrichment + low host signal bonus
    # Rewards regions where viral genes dominate and host contamination is minimal
    if total_genes_interior > 0:
        viral_fraction_for_bonus = (result.gene_taxonomy_viral_top10 or 0) / total_genes_interior
        if viral_fraction_for_bonus >= 0.15 and result.host_signature_fraction < 0.10:
            viral_host_bonus = min(0.12, viral_fraction_for_bonus * 0.4)
            _add_common_bonus(viral_host_bonus)
            logger.debug(
                f"Viral-enrichment bonus for {result.eve_id}: +{viral_host_bonus:.2f} "
                f"(viral_frac={viral_fraction_for_bonus:.2f}, host_sig_frac={result.host_signature_fraction:.2f})"
            )

    configured_priority_markers = [m.lower() for m in (priority_markers or ["mcp"])]
    has_priority = _has_priority_marker(result, configured_priority_markers)

    families = set()
    if result.region_classification_ncldv_markers > 0:
        families.add("NCLDV")
    if result.region_classification_mirus_markers > 0:
        families.add("MIRUS")
    if result.region_classification_vp_plv_markers > 0:
        # Preplasmiviricota is one family however the vp_/plv_ markers split.
        families.add("PPV")
    family_tokens = (
        list(result.hallmark_genes)
        + list(result.marker_family_hits)
        + list(result.interproscan_family_hits)
        + [result.region_classification, result.gene_taxonomy_dominant_family, result.likely_family]
    )
    for token in family_tokens:
        upper = str(token).upper()
        if "NCLDV" in upper:
            families.add("NCLDV")
        elif "MIRUS" in upper:
            families.add("MIRUS")
        elif "PPV" in upper or "PLV" in upper:
            families.add("PPV")
        elif upper in {"VP", "VIROPHAGE"} or "VP_" in upper or "_VP" in upper:
            families.add("PPV")
    family_support_count = len(families)
    priority_floor = 0.0
    if has_priority:
        if family_support_count >= 2:
            _add_common_bonus(marker_multi_family_bonus)
            priority_floor = marker_floor_priority_multi_family
        elif family_support_count >= 1:
            _add_common_bonus(marker_family_bonus_per_family)
            priority_floor = marker_floor_priority_plus_family
        else:
            priority_floor = marker_floor_priority_only

    # Jelly roll bonus (validated DJR MCP provides structural evidence)
    if result.jelly_roll_confidence_bonus > 0:
        _add_common_bonus(result.jelly_roll_confidence_bonus)
        logger.debug(
            f"Jelly roll bonus for {result.eve_id}: +{result.jelly_roll_confidence_bonus:.2f} "
            f"(DJR={result.jelly_roll_djr_count}, avg_conf={result.jelly_roll_avg_confidence:.2f})"
        )

    viral_evidence = result.gene_taxonomy_viral_interior or result.gene_taxonomy_viral_top10

    def _assemble_confidence(*, include_composition: bool) -> tuple[float, dict[str, object]]:
        """Assemble one score from immutable local evidence."""
        component_weights = dict(base_component_weights)
        if not include_composition or composition_score <= 0.0:
            component_weights["composition"] = 0.0
        weight_total = sum(component_weights.values()) or 1.0
        base = sum(
            (component_weights[name] / weight_total) * component_scores[name]
            for name in component_scores
        )

        selected_composition_bonus = composition_bonus if include_composition else 0.0
        bonus = reference_bonus if include_composition else common_bonus

        cap = 0.95
        if seed_cluster_score >= 1.0 and gene_bonus >= 0.6 and marker_score >= 0.6:
            cap = 0.99
        composition_cap_active = (
            include_composition
            and composition_score >= 0.6
            and seed_cluster_score >= 1.0
            and gene_bonus >= 0.7
            and marker_score >= 0.7
            and (crf_score >= 0.9 or not use_crf_score)
        )
        if composition_cap_active:
            cap = 1.0

        confidence = (
            base
            + bonus
            + (0.15 * structural_score)
            - euk_penalty
            - host_signature_penalty
            - size_penalty
        )
        confidence = min(confidence, cap)
        if priority_floor > 0.0:
            confidence = max(confidence, priority_floor)

        # Regions with zero viral taxonomy hits and no MCP cannot be EVEs.
        if viral_evidence == 0 and not result.has_mcp:
            confidence = min(confidence, 0.05)

        # Regions below 5% viral genes and without MCP stay below MEDIUM.
        if total_genes_interior > 0 and not result.has_mcp and viral_fraction < 0.05:
            confidence = min(confidence, 0.19)

        score_components: dict[str, object] = {
            "scores": dict(component_scores),
            "weights": component_weights,
            "weighted_base": float(base),
            "bonus_total": float(bonus),
            "composition_bonus": float(selected_composition_bonus),
            "composition_cap_active": composition_cap_active,
            "composition_evidence_active": include_composition,
            "priority_floor": float(priority_floor),
            "cap": float(cap),
            "family_consistency_bonus": float(family_consistency_bonus),
            "has_mcp": bool(result.has_mcp),
            "priority_marker_active": bool(has_priority),
            "final_confidence_pre_clamp": float(confidence),
        }
        return max(0.0, min(1.0, confidence)), score_components

    reference_confidence, reference_components = _assemble_confidence(
        include_composition=True
    )
    if ablation_id is AblationID.A5:
        selected_confidence, selected_components = _assemble_confidence(
            include_composition=False
        )
    else:
        selected_confidence = reference_confidence
        selected_components = reference_components

    # Apply the selected score metadata once. The A0 counterfactual above is
    # pure and never overwrites the candidate while A5 is being evaluated.
    result.ablation_id = ablation_id
    result.high_confidence_euk_genes = high_conf_euk
    result.family_consistency_score = family_consistency_bonus
    result.score_components = selected_components
    result.composition_ablation_effect = evaluate_composition_ablation_effect(
        ablation_id=ablation_id,
        composition_score=composition_score,
        reference_confidence=reference_confidence,
        selected_confidence=selected_confidence,
        high_threshold=high_tier_threshold,
        low_threshold=low_tier_threshold,
    )
    return selected_confidence


@dataclass(frozen=True, slots=True)
class _PostScoreOutcome:
    confidence: float
    tier: str
    priority_promoted: bool
    non_hhg_demoted: bool
    non_host_genes: int


def _evaluate_post_score_policy(
    result: VerificationResult,
    confidence: float,
    *,
    high_threshold: float,
    low_threshold: float,
    priority_markers: list[str],
) -> _PostScoreOutcome:
    """Apply final tier policy to a score without mutating ``result``."""
    tier = _confidence_tier_for_score(
        confidence,
        high=high_threshold,
        low=low_threshold,
    )
    priority_promoted = _has_priority_marker(result, priority_markers) and tier == "LOW"
    if priority_promoted:
        tier = "MEDIUM"
        confidence = max(confidence, low_threshold)

    non_host_genes = result.gene_count - result.host_signature_gene_count
    non_hhg_demoted = False
    if tier in {"HIGH", "MEDIUM"} and result.seed_sources and "hhg" not in result.seed_sources:
        non_hhg_quality_pass = (
            result.hallmark_count >= 3
            or (result.hallmark_count >= 1 and result.has_mcp)
            or non_host_genes >= 5
        )
        if not non_hhg_quality_pass:
            non_hhg_demoted = True
            tier = "LOW"
            confidence = min(confidence, low_threshold - 0.001)

    return _PostScoreOutcome(
        confidence=confidence,
        tier=tier,
        priority_promoted=priority_promoted,
        non_hhg_demoted=non_hhg_demoted,
        non_host_genes=non_host_genes,
    )


class EvidenceSynthesizer:
    """
    Main evidence synthesis engine with gated escalation.

    Implements the verification pipeline that combines marker, taxonomy,
    compositional, structural/domain, and phylogenetic evidence.

    Tie-breaker modules:
    1. Evidence Coherence Analysis - checks multi-evidence agreement
    2. Structural Homology - Boltz + FoldSeek against viral proteins
    3. Phylogenetic Validation - GVClass + Diamond BLASTp for final confirmation
    """

    def __init__(
        self,
        config: Optional[EvidenceSynthesizerConfig] = None,
        viral_structure_db: Optional[Path] = None,
        genome_path: Optional[Path] = None,
        work_dir: Optional[Path] = None,
    ):
        """
        Initialize evidence synthesizer.

        Args:
            config: Configuration settings
            viral_structure_db: Path to viral structure database for FoldSeek
            genome_path: Path to input genome FASTA (needed for phylogenetic validation)
            work_dir: Working directory for intermediate files
        """
        self.config = config or EvidenceSynthesizerConfig()
        self.viral_structure_db = viral_structure_db
        self.genome_path = genome_path
        self.work_dir = work_dir

        # Lazy-loaded components
        self._boltz_analyzer = None
        self._phylogenetic_validator = None
        self._tmvec_searcher = None
        self._marker_annotation_index = None

        self._phase1_host_signature_host_prefixes = "."
        self._phase1_host_signature_top_tokens = "."
        self._phase1_host_signature_token_count = 0
        self._initialize_phase1_host_signature_summary()

    def _initialize_phase1_host_signature_summary(self) -> None:
        """Derive a compact Phase 1 host-signature summary from config payload."""
        payload = self.config.host_signature_model or {}
        token_weights = payload.get("token_weights") or {}
        host_prefixes = payload.get("host_prefixes") or self.config.host_prefixes or []

        if host_prefixes:
            self._phase1_host_signature_host_prefixes = ",".join(str(p) for p in host_prefixes)
        self._phase1_host_signature_token_count = len(token_weights)

        if not token_weights:
            return

        top_items = sorted(
            token_weights.items(),
            key=lambda item: (-float(item[1]), str(item[0])),
        )[:10]
        self._phase1_host_signature_top_tokens = "|".join(
            f"{token}:{float(weight):.2f}" for token, weight in top_items
        )

    @property
    def boltz_analyzer(self) -> Optional[BoltzFoldSeekAnalyzer]:
        """Lazy load Boltz + FoldSeek analyzer (optional)."""
        if not self.config.use_boltz:
            return None

        if self._boltz_analyzer is None:
            from .structural_homology import BoltzFoldSeekAnalyzer

            self._boltz_analyzer = BoltzFoldSeekAnalyzer(
                viral_db_path=self.viral_structure_db,
                device=self.config.device,
                threads=self.config.threads,
                use_msa_server=self.config.boltz_use_msa_server,
                min_seq_len=self.config.boltz_min_seq_len,
                max_seq_len=self.config.boltz_max_seq_len,
                no_kernels=self.config.boltz_no_kernels,
            )

        return self._boltz_analyzer

    @property
    def phylogenetic_validator(self) -> Optional[PhylogeneticValidator]:
        """Lazy load phylogenetic validator."""
        if not self.config.use_phylogenetic_validation:
            return None

        if not self.genome_path:
            logger.warning("Genome path not provided, phylogenetic validation disabled")
            return None

        if self._phylogenetic_validator is None:
            phylo_work_dir = self.work_dir / "phylogenetic" if self.work_dir else None
            if phylo_work_dir:
                phylo_work_dir.mkdir(parents=True, exist_ok=True)

            self._phylogenetic_validator = PhylogeneticValidator(
                genome_path=self.genome_path,
                work_dir=phylo_work_dir,
                gvclass_db=self.config.gvclass_db,
                diamond_db=self.config.diamond_db,
                threads=self.config.threads,
            )

        return self._phylogenetic_validator

    @property
    def tmvec_searcher(self) -> Optional[TMVecDatabaseSearch]:
        """Lazy load TMVec database searcher."""
        if not self.config.use_tmvec_database:
            return None
        if self._tmvec_searcher is None:
            from .tmvec_database import TMVecDatabaseSearch

            self._tmvec_searcher = TMVecDatabaseSearch(
                device=self.config.device,
                databases=self.config.tmvec_databases,
                min_tm=0.0,
                database_root=self.config.tmvec_database_dir,
                require_gpu=self.config.tmvec_require_gpu,
            )
        return self._tmvec_searcher

    @property
    def marker_annotation_index(self) -> dict[str, dict]:
        """Lazy load marker annotation index for category scoring."""
        if self._marker_annotation_index is None:
            self._marker_annotation_index = load_marker_annotation_index(
                self.config.marker_annotations_path
            )
        return self._marker_annotation_index or {}

    def _initialize_result(
        self,
        refined_boundary: RefinedBoundary,
    ) -> VerificationResult:
        """Initialize VerificationResult with boundary and seed metadata."""
        eve_id = f"EVE_{refined_boundary.scaffold}_{refined_boundary.start}-{refined_boundary.end}"
        result = VerificationResult(
            eve_id=eve_id,
            scaffold=refined_boundary.scaffold,
            start=refined_boundary.start,
            end=refined_boundary.end,
            ablation_id=self.config.ablation_id,
            crf_confidence=refined_boundary.confidence,
            crf_posterior=refined_boundary.posterior_probability,
            seed_sources=list(getattr(refined_boundary, "seed_sources", [])),
            seed_confidence=getattr(refined_boundary, "seed_confidence", ""),
            seed_hhg_score=getattr(refined_boundary, "seed_hhg_score", 0.0),
            seed_novelty_score=getattr(refined_boundary, "seed_novelty_score", 0.0),
            seed_compositional_score=getattr(refined_boundary, "seed_compositional_score", 0.0),
            region_classification=getattr(refined_boundary, "predicted_family", ""),
            region_classification_ncldv_markers=getattr(refined_boundary, "region_classification_ncldv_markers", 0),
            region_classification_vp_plv_markers=getattr(refined_boundary, "region_classification_vp_plv_markers", 0),
            region_classification_mirus_markers=getattr(refined_boundary, "region_classification_mirus_markers", 0),
            kfd=getattr(refined_boundary, "max_kfd", 0.0),
            gc_deviation=getattr(refined_boundary, "gc_deviation", 0.0),
            cub_deviation=getattr(refined_boundary, "cub_deviation", 0.0),
        )
        candidate_start = getattr(refined_boundary, "candidate_start", None)
        candidate_end = getattr(refined_boundary, "candidate_end", None)
        if candidate_start is not None and candidate_end is not None:
            result.candidate_start = candidate_start
            result.candidate_end = candidate_end
            result.candidate_length = max(0, candidate_end - candidate_start)
            result.candidate_reduction_bp = max(
                0, result.candidate_length - (result.end - result.start)
            )
            result.candidate_reduction_reason = getattr(refined_boundary, "host_trim_reason", "")
        result.phase1_host_signature_host_prefixes = self._phase1_host_signature_host_prefixes
        result.phase1_host_signature_top_tokens = self._phase1_host_signature_top_tokens
        result.phase1_host_signature_token_count = self._phase1_host_signature_token_count
        return result

    def _detect_contig_edge(
        self,
        result: VerificationResult,
        refined_boundary: RefinedBoundary,
        scaffold_lengths: Optional[dict[str, int]],
    ) -> None:
        """Detect if EVE is at contig boundary (partial EVE)."""
        EDGE_BUFFER_BP = 5000
        if not scaffold_lengths:
            return
        scaffold_length = scaffold_lengths.get(refined_boundary.scaffold, 0)
        result.scaffold_length = scaffold_length
        if scaffold_length > 0:
            at_start = refined_boundary.start < EDGE_BUFFER_BP
            at_end = (scaffold_length - refined_boundary.end) < EDGE_BUFFER_BP
            result.partial_eve_at_start = at_start
            result.partial_eve_at_end = at_end
            result.partial_eve = at_start or at_end
            if result.partial_eve:
                logger.debug(
                    f"Partial EVE detected: {result.eve_id} "
                    f"(at_start={at_start}, at_end={at_end}, "
                    f"scaffold_length={scaffold_length})"
                )

    @staticmethod
    def _deduplicate_hallmark_hits(
        hallmark_hits: list,
    ) -> tuple[list, list[str], dict[str, dict]]:
        """Deduplicate hallmark hits by base gene ID.

        When multiple HMM profiles hit the same protein, keep only the
        best-scoring model (highest hmm_score) per gene.

        Returns:
            Tuple of (deduplicated hit list, hallmark gene names, by_gene dict)
        """
        by_gene: dict[str, dict] = {}
        for hit in hallmark_hits:
            if isinstance(hit, dict):
                gene_name = hit.get("hallmark_gene", "")
                porf_id = hit.get("porf_id", "") or ""
                score = float(hit.get("hmm_score", 0.0) or hit.get("score", 0.0) or 0.0)
            else:
                gene_name = getattr(hit, "hallmark_gene", "")
                porf_id = getattr(hit, "porf_id", None) or getattr(hit, "query_porf", "") or ""
                score = float(getattr(hit, "hmm_score", 0.0) or getattr(hit, "score", 0.0) or 0.0)
            if not gene_name:
                continue
            # Extract base gene ID (remove |aaX-Y domain suffix)
            base_gene = porf_id.split("|aa")[0] if porf_id else gene_name
            if base_gene not in by_gene or score > by_gene[base_gene]["score"]:
                by_gene[base_gene] = {
                    "hit": hit,
                    "hallmark_gene": gene_name,
                    "score": score,
                    "base_gene": base_gene,
                }
        deduped_hits = [info["hit"] for info in by_gene.values()]
        hallmark_genes = [info["hallmark_gene"] for info in by_gene.values()]
        return deduped_hits, hallmark_genes, by_gene

    def _process_hallmark_hits(
        self,
        result: VerificationResult,
        hallmark_hits: Optional[list],
    ) -> None:
        """Process hallmark hits and update result with marker evidence.

        Deduplicates by gene: when multiple HMM profiles hit the same protein,
        only the best-scoring model counts. hallmark_count reflects unique
        marker-bearing proteins, not the number of HMM profile hits.
        """
        if not hallmark_hits:
            return
        deduped_hits, hallmark_genes, by_gene = self._deduplicate_hallmark_hits(
            hallmark_hits,
        )
        if not hallmark_genes:
            return
        result.hallmark_genes = hallmark_genes
        result.hallmark_count = len(hallmark_genes)
        result.hallmark_diversity = len(set(hallmark_genes))
        result.has_virus_specific_marker = True
        result.has_mcp = any(is_mcp_gene(g) for g in hallmark_genes)
        result.mcp_gene_ids = [
            info["base_gene"]
            for info in by_gene.values()
            if is_mcp_gene(info["hallmark_gene"])
        ]
        bypassed = [
            info
            for info in by_gene.values()
            if isinstance(info["hit"], dict)
            and bool(info["hit"].get("tier1_bypassed", False))
        ]
        result.tier1_bypassed_marker_ids = [
            info["base_gene"] for info in bypassed
        ]
        result.tier1_bypassed_marker_models = [
            info["hallmark_gene"] for info in bypassed
        ]
        marker_summary = summarize_marker_hits(
            hallmark_genes,
            annotation_index=self.marker_annotation_index,
        )
        (
            identity_qualified_cress_genes,
            specific_top1_cress_genes,
        ) = _cress_gene_support(hallmark_hits)
        canonical_cress_marker_support = (
            len(identity_qualified_cress_genes) >= 2
            or bool(specific_top1_cress_genes)
        )
        compact_cress_boundary = (
            canonical_family(result.region_classification) == "CRESS"
        )
        if compact_cress_boundary and canonical_cress_marker_support:
            marker_summary["families"] = sorted(
                {*marker_summary["families"], "CRESS"}
            )
            if marker_summary["dominant_family"] == "UNKNOWN":
                marker_summary["dominant_family"] = "CRESS"
                marker_summary["dominant_fraction"] = 1.0
        result.marker_category_hits = marker_summary["categories"]
        result.marker_family_hits = marker_summary["families"]
        result.marker_complement_score = marker_summary["complement_score"]
        result.marker_dominant_family = marker_summary["dominant_family"]
        result.marker_dominant_fraction = marker_summary["dominant_fraction"]
        result.has_mcp = result.has_mcp or marker_summary["has_mcp"]
        scored_hallmarks = [(info["hallmark_gene"], info["score"]) for info in by_gene.values()]
        result.likely_group = infer_likely_group(scored_hallmarks, self.marker_annotation_index)
        completeness = compute_marker_completeness(
            hallmark_genes,
            annotation_index=self.marker_annotation_index,
            hallmark_hits=deduped_hits,
        )
        result.vp_completeness = completeness["vp_completeness"]
        result.vp_completeness_ratio = completeness["vp_completeness_ratio"]
        result.ppv_completeness = completeness["ppv_completeness"]
        result.ppv_completeness_ratio = completeness["ppv_completeness_ratio"]
        result.ncldv_completeness = completeness["ncldv_completeness"]
        result.ncldv_completeness_ratio = completeness["ncldv_completeness_ratio"]
        result.mirus_completeness = completeness["mirus_completeness"]
        result.mirus_completeness_ratio = completeness["mirus_completeness_ratio"]
        result.ppv_subtype = infer_ppv_subtype(hallmark_genes)
        if result.region_classification in {"VP", "PLV", "PPV", "UNKNOWN", ""}:
            if result.ppv_subtype in {"VP", "PLV"}:
                # The subgroup stays in ppv_subtype; the class token is the
                # unified Preplasmiviricota lineage.
                result.region_classification = "PPV"

    def _process_gene_taxonomy(
        self,
        result: VerificationResult,
        gene_taxonomy_records: Optional[list],
        gene_taxonomy_summary: Optional[dict],
    ) -> None:
        """Process gene taxonomy and host signature evidence."""
        if gene_taxonomy_summary:
            result.gene_count = gene_taxonomy_summary.get("total", 0)
            result.genes_with_ncldv_mirus_top10 = gene_taxonomy_summary.get("ncldv_mirus", 0)
            result.genes_with_vp_plv_top10 = gene_taxonomy_summary.get("vp_plv", 0)
            result.genes_with_high_pident_euk = gene_taxonomy_summary.get("high_pident_euk", 0)
            result.gene_taxonomy_total = gene_taxonomy_summary.get("total", 0)
            result.gene_taxonomy_ncldv_top10 = gene_taxonomy_summary.get("ncldv_mirus", 0)
            result.gene_taxonomy_vp_plv_top10 = gene_taxonomy_summary.get("vp_plv", 0)
            result.gene_taxonomy_has_ncldv_mirus = gene_taxonomy_summary.get("has_ncldv_mirus", False)
            result.gene_taxonomy_has_vp_plv = gene_taxonomy_summary.get("has_vp_plv", False)
            result.gene_taxonomy_viral_top10 = gene_taxonomy_summary.get("viral_top10", 0)
            result.gene_taxonomy_dominant_family = gene_taxonomy_summary.get(
                "dominant_family", "UNKNOWN"
            )
            result.gene_taxonomy_dominant_fraction = gene_taxonomy_summary.get(
                "dominant_fraction", 0.0
            )

            # NEW: Flanking gene tracking (taxonomy expansion fix)
            result.gene_taxonomy_total_with_flanking = gene_taxonomy_summary.get(
                "total_with_flanking", result.gene_count
            )
            result.gene_taxonomy_flanking_count = gene_taxonomy_summary.get("flanking_genes", 0)
            result.gene_taxonomy_viral_interior = gene_taxonomy_summary.get("viral_interior", 0)
            result.gene_taxonomy_viral_flanking = gene_taxonomy_summary.get("viral_flanking", 0)
            result.gene_taxonomy_ncldv_mirus_interior = gene_taxonomy_summary.get(
                "ncldv_mirus_interior", 0
            )
            result.gene_taxonomy_ncldv_mirus_flanking = gene_taxonomy_summary.get(
                "ncldv_mirus_flanking", 0
            )
            result.gene_taxonomy_vp_plv_interior = gene_taxonomy_summary.get("vp_plv_interior", 0)
            result.gene_taxonomy_vp_plv_flanking = gene_taxonomy_summary.get("vp_plv_flanking", 0)

            # Taxonomy distribution analysis fields
            result.taxonomy_distribution_viral_score = gene_taxonomy_summary.get(
                "taxonomy_distribution_viral_score", 0.0
            )
            result.taxonomy_distribution_diversity = gene_taxonomy_summary.get(
                "taxonomy_distribution_diversity", 0.0
            )
            result.taxonomy_distribution_host_overlap = gene_taxonomy_summary.get(
                "taxonomy_distribution_host_overlap", 0.0
            )
            result.taxonomy_distribution_non_euk_fraction = gene_taxonomy_summary.get(
                "taxonomy_distribution_non_euk_fraction", 0.0
            )
            result.taxonomy_distribution_genes_analyzed = gene_taxonomy_summary.get(
                "taxonomy_distribution_genes_analyzed", 0
            )
            result.taxonomy_distribution_likely_viral_genes = gene_taxonomy_summary.get(
                "taxonomy_distribution_likely_viral_genes", 0
            )
            result.taxonomy_distribution_likely_host_genes = gene_taxonomy_summary.get(
                "taxonomy_distribution_likely_host_genes", 0
            )
            result.taxonomy_distribution_baseline_source = gene_taxonomy_summary.get(
                "taxonomy_distribution_baseline_source", "default"
            )
            result.taxonomy_distribution_baseline_markers = gene_taxonomy_summary.get(
                "taxonomy_distribution_baseline_markers", 0
            )
            result.taxonomy_distribution_baseline_genera = gene_taxonomy_summary.get(
                "taxonomy_distribution_baseline_genera", 0
            )
            result.taxonomy_distribution_baseline_diversity = gene_taxonomy_summary.get(
                "taxonomy_distribution_baseline_diversity", 0.0
            )
        if gene_taxonomy_records:
            result.gene_taxonomy_records = [
                getattr(r, "__dict__", r) for r in gene_taxonomy_records
            ]
            # Filter out flanking genes for host signature calculation
            # Flanking genes should not affect the host signature penalty since
            # they are outside the EVE boundary and gene_count only counts interior genes
            interior_records = [
                r for r in gene_taxonomy_records
                if not (r.get("is_flanking") if isinstance(r, dict) else getattr(r, "is_flanking", False))
            ]
            if self.config.host_signature_model:
                model = HostSignatureModel.from_dict(self.config.host_signature_model)
                host_like, mean_score, evalue_weighted = host_signature_density_evalue_weighted(
                    interior_records,  # Use interior genes only
                    model,
                    score_threshold=self.config.host_signature_score_threshold,
                    host_prefixes=self.config.host_prefixes,
                )
                result.host_signature_gene_count = host_like
                result.host_signature_weighted_mean = mean_score
                result.host_signature_evalue_weighted = evalue_weighted
                if result.gene_count:
                    result.host_signature_fraction = host_like / result.gene_count
                if host_like == 0 and mean_score > 0 and result.genes_with_high_pident_euk > 0:
                    logger.debug(
                        "%s: host_signature_gene_count=0 despite %d high-pident EUK genes "
                        "(mean_score=%.4f, threshold=%.2f) — EUK/host penalties disabled",
                        result.eve_id, result.genes_with_high_pident_euk,
                        mean_score, self.config.host_signature_score_threshold,
                    )
            elif self.config.euk_host_signatures:
                host_like = 0
                for record in interior_records:  # Use interior genes only
                    top1_target = record.get("top1_target") if isinstance(record, dict) else ""
                    host_prefixes = set(self.config.host_prefixes or ["EUK__"])
                    if not top1_target or not any(top1_target.startswith(p) for p in host_prefixes):
                        continue
                    token = top1_target.split("|", 1)[0]
                    if token in self.config.euk_host_signatures:
                        host_like += 1
                result.host_signature_gene_count = host_like
                if result.gene_count:
                    result.host_signature_fraction = host_like / result.gene_count

    def _process_interproscan(
        self,
        result: VerificationResult,
        interproscan_summary: Optional[dict],
    ) -> None:
        """Process InterProScan annotation results."""
        if not interproscan_summary:
            return
        result.interproscan_total_hits = interproscan_summary.get("total_hits", 0)
        result.interproscan_viral_hits = interproscan_summary.get("viral_hits", 0)
        result.interproscan_keyword_hits = interproscan_summary.get("keyword_hits", [])
        result.interproscan_category_hits = interproscan_summary.get("category_hits", [])
        result.interproscan_family_hits = interproscan_summary.get("family_hits", [])
        result.interproscan_category_score = interproscan_summary.get("category_score", 0.0)

        # NUMT detection (mitochondrial markers)
        result.interproscan_numt_hits = interproscan_summary.get("numt_hits", 0)
        result.interproscan_numt_markers = interproscan_summary.get("numt_markers", [])
        result.numt_flag = "DETECTED" if result.interproscan_numt_hits > 0 else "NONE"

        if result.interproscan_total_hits > 0:
            ratio_score = min(
                result.interproscan_viral_hits / result.interproscan_total_hits,
                1.0,
            )
        else:
            ratio_score = 0.0
        result.interproscan_score = max(ratio_score, result.interproscan_category_score)

    def _apply_jelly_roll_summary(
        self,
        result: VerificationResult,
        jelly_roll_summary: Optional[dict],
    ) -> None:
        """Apply per-boundary DJR/SJR MCP summary before confidence scoring."""
        if not jelly_roll_summary:
            return
        result.jelly_roll_djr_count = int(jelly_roll_summary.get("djr_count", 0) or 0)
        result.jelly_roll_sjr_count = int(jelly_roll_summary.get("sjr_count", 0) or 0)
        result.jelly_roll_total_mcp = int(jelly_roll_summary.get("total_mcp", 0) or 0)
        result.jelly_roll_avg_confidence = float(jelly_roll_summary.get("avg_confidence", 0.0) or 0.0)
        result.jelly_roll_confidence_bonus = float(
            jelly_roll_summary.get("confidence_bonus", 0.0) or 0.0
        )
        proteins = jelly_roll_summary.get("mcp_proteins", [])
        result.jelly_roll_mcp_proteins = proteins if isinstance(proteins, list) else []
        if result.jelly_roll_total_mcp > 0:
            result.has_mcp = True

    def _run_tmvec_database_scan(
        self,
        result: VerificationResult,
        porf_sequences: list[tuple[str, str]],
        precomputed_tmvec: Optional[dict[str, dict]] = None,
    ) -> None:
        """
        Run TMVec database scan for proteins in an EVE.

        Args:
            result: VerificationResult to update
            porf_sequences: List of (porf_id, sequence) tuples
            precomputed_tmvec: Optional pre-computed TMVec results from batch processing.
                              If provided, skips embedding and uses these results directly.
        """
        tmvec_records: list[dict] = []
        best_bfvd = 0.0

        # Use pre-computed results if available (from batch processing)
        if precomputed_tmvec is not None:
            for porf_id, sequence in porf_sequences:
                hits = precomputed_tmvec.get(porf_id, {})
                bfvd_hit = hits.get("bfvd")
                cath_hit = hits.get("cath")
                swiss_hit = hits.get("swissprot")

                bfvd_score = bfvd_hit.tm_score if bfvd_hit else 0.0
                cath_score = cath_hit.tm_score if cath_hit else 0.0
                swiss_score = swiss_hit.tm_score if swiss_hit else 0.0
                viral_specificity = bfvd_score - max(cath_score, swiss_score)

                record = {
                    "eve_id": result.eve_id,
                    "porf_id": porf_id,
                    "length": len(sequence),
                    "tmvec_bfvd_score": bfvd_score,
                    "tmvec_bfvd_hit": bfvd_hit.target_id if bfvd_hit else "",
                    "tmvec_bfvd_annotation": bfvd_hit.protein_name if bfvd_hit else "",
                    "tmvec_bfvd_organism": bfvd_hit.organism if bfvd_hit else "",
                    "tmvec_bfvd_lineage": bfvd_hit.lineage if bfvd_hit else "",
                    "tmvec_bfvd_keywords": bfvd_hit.keywords if bfvd_hit else "",
                    "tmvec_cath_score": cath_score,
                    "tmvec_cath_hit": cath_hit.target_id if cath_hit else "",
                    "tmvec_swiss_score": swiss_score,
                    "tmvec_swiss_hit": swiss_hit.target_id if swiss_hit else "",
                    "tmvec_viral_specificity": viral_specificity,
                }
                tmvec_records.append(record)
                best_bfvd = max(best_bfvd, bfvd_score)
        else:
            # Fall back to per-protein search (slower, used when not batch processing)
            searcher = self.tmvec_searcher
            if searcher is None:
                return

            for porf_id, sequence in porf_sequences:
                hits = searcher.search_sequence(sequence)
                bfvd_hit = hits.get("bfvd")
                cath_hit = hits.get("cath")
                swiss_hit = hits.get("swissprot")

                bfvd_score = bfvd_hit.tm_score if bfvd_hit else 0.0
                cath_score = cath_hit.tm_score if cath_hit else 0.0
                swiss_score = swiss_hit.tm_score if swiss_hit else 0.0
                viral_specificity = bfvd_score - max(cath_score, swiss_score)

                record = {
                    "eve_id": result.eve_id,
                    "porf_id": porf_id,
                    "length": len(sequence),
                    "tmvec_bfvd_score": bfvd_score,
                    "tmvec_bfvd_hit": bfvd_hit.target_id if bfvd_hit else "",
                    "tmvec_bfvd_annotation": bfvd_hit.protein_name if bfvd_hit else "",
                    "tmvec_bfvd_organism": bfvd_hit.organism if bfvd_hit else "",
                    "tmvec_bfvd_lineage": bfvd_hit.lineage if bfvd_hit else "",
                    "tmvec_bfvd_keywords": bfvd_hit.keywords if bfvd_hit else "",
                    "tmvec_cath_score": cath_score,
                    "tmvec_cath_hit": cath_hit.target_id if cath_hit else "",
                    "tmvec_swiss_score": swiss_score,
                    "tmvec_swiss_hit": swiss_hit.target_id if swiss_hit else "",
                    "tmvec_viral_specificity": viral_specificity,
                }
                tmvec_records.append(record)
                best_bfvd = max(best_bfvd, bfvd_score)

        result.tmvec_all_proteins = tmvec_records
        if best_bfvd >= self.config.tmvec_min_score:
            result.structural_score = max(result.structural_score, best_bfvd)
            result.has_structural_support = True

    def _filter_boltz_mcp_sequences(
        self,
        hallmark_hits: Optional[list],
        porf_sequences: list[tuple[str, str]],
    ) -> list[tuple[str, str]]:
        if not hallmark_hits:
            return []

        mcp_ids: set[str] = set()
        for hit in hallmark_hits:
            if isinstance(hit, dict):
                gene = (hit.get("hallmark_gene") or "").lower()
                porf_id = hit.get("porf_id") or hit.get("query_porf")
            else:
                gene = (getattr(hit, "hallmark_gene", "") or "").lower()
                porf_id = getattr(hit, "porf_id", None) or getattr(hit, "query_porf", None)
            if not gene or not porf_id:
                continue
            if any(token in gene for token in ("mcp", "jelly", "capsid", "djr", "sjr")):
                mcp_ids.add(str(porf_id))

        if not mcp_ids:
            return []

        return [(pid, seq) for pid, seq in porf_sequences if pid in mcp_ids]

    def _run_tiebreakers(
        self,
        result: VerificationResult,
        refined_boundary: RefinedBoundary,
        window_features: list,
        hallmark_hits: Optional[list],
        novelty_scores: Optional[dict],
        porf_sequences: Optional[list[tuple[str, str]]],
        precomputed_tmvec: Optional[dict[str, dict]] = None,
    ) -> None:
        """Run tie-breaker modules: coherence, structural, phylogenetic.

        Args:
            precomputed_tmvec: Optional pre-computed TMVec results from batch processing.
                              If provided, uses these instead of running TMVec per-protein.
        """
        # Tie-breaker 1: Evidence Coherence Analysis
        crf_states = refined_boundary.state_sequence or []
        crf_posteriors = refined_boundary.state_posteriors
        coherence = analyze_eve_coherence(
            eve_id=result.eve_id,
            scaffold=refined_boundary.scaffold,
            start=refined_boundary.start,
            end=refined_boundary.end,
            window_features=window_features,
            crf_states=crf_states,
            crf_posteriors=crf_posteriors,
            hallmark_hits=hallmark_hits,
            novelty_scores=novelty_scores,
        )
        result.coherence_analysis = coherence
        result.coherence_score = coherence.coherence_score
        if coherence.profile:
            result.hallmark_diversity = max(
                result.hallmark_diversity,
                coherence.profile.hallmark_diversity,
            )
            result.has_virus_specific_marker = (
                result.has_virus_specific_marker or coherence.profile.has_virus_specific_marker
            )
            hallmark_count = sum(
                1 for k in coherence.profile.evidence_counts.keys()
                if k.value.startswith("hallmark_")
            )
            result.hallmark_count = max(result.hallmark_count, hallmark_count)

        # Tie-breaker 2: TMVec database scan across all EVE proteins
        if self.config.use_tmvec_database and porf_sequences:
            self._run_tmvec_database_scan(result, porf_sequences, precomputed_tmvec)

        # Tie-breaker 3: Boltz + FoldSeek structural analysis (optional, MCP-only by default)
        if self.config.use_boltz and self.boltz_analyzer and porf_sequences and self.work_dir:
            sequences_to_analyze = porf_sequences[: self.config.structural_max_porfs]
            if self.config.boltz_mcp_only:
                sequences_to_analyze = self._filter_boltz_mcp_sequences(
                    hallmark_hits,
                    sequences_to_analyze,
                )
            if sequences_to_analyze:
                eve_work_dir = self.work_dir / safe_filename_component(result.eve_id)
                require_strict_child(self.work_dir, eve_work_dir)
                struct_work_dir = eve_work_dir / "structures"
                require_strict_child(self.work_dir, struct_work_dir)
                structural_results = self.boltz_analyzer.analyze_batch(
                    sequences_to_analyze,
                    work_dir=struct_work_dir,
                )
                result.structural_results = structural_results
                if structural_results:
                    scores = [r.structural_evidence_score for r in structural_results]
                    result.structural_score = max(
                        result.structural_score,
                        max(scores) if scores else 0.0,
                    )
                    result.has_structural_support = result.has_structural_support or any(
                        r.supports_viral_origin for r in structural_results
                    )

        # Tie-breaker 4: Phylogenetic Validation
        self._run_phylogenetic_validation(result, refined_boundary)

    def _calculate_final_decision(
        self,
        result: VerificationResult,
        refined_boundary: RefinedBoundary,
        has_hhg_evidence: bool,
    ) -> None:
        """Calculate final score and make accept/reject decision."""
        _ = has_hhg_evidence  # retained for call compatibility; HHG contributes via marker evidence now.
        selected_score = calculate_eve_confidence(
            result=result,
            crf_confidence=refined_boundary.confidence,
            tmvec_score=None,
            use_crf_score=self.config.use_crf_in_final_score,
            priority_markers=self.config.priority_marker_list,
            marker_floor_priority_only=self.config.marker_floor_priority_only,
            marker_floor_priority_plus_family=self.config.marker_floor_priority_plus_family,
            marker_floor_priority_multi_family=self.config.marker_floor_priority_multi_family,
            marker_family_bonus_per_family=self.config.marker_family_bonus_per_family,
            marker_multi_family_bonus=self.config.marker_multi_family_bonus,
            host_signature_score_threshold=self.config.host_signature_score_threshold,
            ablation_id=self.config.ablation_id,
            high_tier_threshold=self.config.high_tier_threshold,
            low_tier_threshold=self.config.low_tier_threshold,
        )
        selected_outcome = _evaluate_post_score_policy(
            result,
            selected_score,
            high_threshold=self.config.high_tier_threshold,
            low_threshold=self.config.low_tier_threshold,
            priority_markers=self.config.priority_marker_list,
        )
        result.final_confidence = selected_outcome.confidence
        result.confidence_tier = selected_outcome.tier

        effect = result.composition_ablation_effect
        if self.config.ablation_id is AblationID.A5 and effect.reference_confidence is not None:
            reference_outcome = _evaluate_post_score_policy(
                result,
                effect.reference_confidence,
                high_threshold=self.config.high_tier_threshold,
                low_threshold=self.config.low_tier_threshold,
                priority_markers=self.config.priority_marker_list,
            )
            result.composition_ablation_effect = evaluate_composition_ablation_effect(
                ablation_id=self.config.ablation_id,
                composition_score=effect.composition_score,
                reference_confidence=reference_outcome.confidence,
                selected_confidence=selected_outcome.confidence,
                high_threshold=self.config.high_tier_threshold,
                low_threshold=self.config.low_tier_threshold,
                reference_tier=reference_outcome.tier,
                selected_tier=selected_outcome.tier,
            )

        if selected_outcome.priority_promoted:
            logger.info(
                "%s: promoted to MEDIUM confidence due to priority-marker evidence",
                result.eve_id,
            )
        if selected_outcome.non_hhg_demoted:
            logger.info(
                "%s: demoted to LOW — non-HHG seed (sources: %s) with insufficient quality "
                "(hallmark=%d, has_mcp=%s, non_host_genes=%d)",
                result.eve_id,
                result.seed_sources,
                result.hallmark_count,
                result.has_mcp,
                selected_outcome.non_host_genes,
            )

        # Set status for backward compatibility with existing code
        if result.confidence_tier == "HIGH":
            logger.info(
                f"{result.eve_id}: HIGH confidence (score: {result.final_confidence:.3f})"
            )
            result.status = VerificationStatus.HIGH_CONFIDENCE
        elif result.confidence_tier == "MEDIUM":
            logger.info(
                f"{result.eve_id}: MEDIUM confidence (score: {result.final_confidence:.3f})"
            )
            result.status = VerificationStatus.MEDIUM_CONFIDENCE
        else:
            logger.info(
                f"{result.eve_id}: LOW confidence (score: {result.final_confidence:.3f})"
            )
            result.status = VerificationStatus.LOW_CONFIDENCE_TIEBREAKER

    def verify_eve(
        self,
        refined_boundary: RefinedBoundary,
        window_features: list,
        hallmark_hits: Optional[list] = None,
        novelty_scores: Optional[dict] = None,
        porf_sequences: Optional[list[tuple[str, str]]] = None,
        gene_taxonomy_records: Optional[list] = None,
        gene_taxonomy_summary: Optional[dict] = None,
        interproscan_summary: Optional[dict] = None,
        jelly_roll_summary: Optional[dict] = None,
        scaffold_lengths: Optional[dict[str, int]] = None,
        precomputed_tmvec: Optional[dict[str, dict]] = None,
    ) -> VerificationResult:
        """
        Verify a single EVE candidate through gated escalation.

        Args:
            refined_boundary: RefinedBoundary from Phase 2
            window_features: Legacy boundary feature payload
            hallmark_hits: HMM hallmark hits in region
            novelty_scores: Legacy compatibility payload
            porf_sequences: List of (porf_id, sequence) for structural analysis
            gene_taxonomy_records: Gene taxonomy records
            gene_taxonomy_summary: Gene taxonomy summary
            interproscan_summary: InterProScan results
            jelly_roll_summary: DJR/SJR MCP summary for this boundary
            scaffold_lengths: Scaffold lengths for contig-edge detection
            precomputed_tmvec: Pre-computed TMVec results from batch processing

        Returns:
            VerificationResult with complete verification
        """
        # Initialize result with boundary and seed metadata
        result = self._initialize_result(refined_boundary)

        # Process evidence: contig-edge, hallmarks, taxonomy, interproscan
        self._detect_contig_edge(result, refined_boundary, scaffold_lengths)
        self._process_hallmark_hits(result, hallmark_hits)
        self._process_gene_taxonomy(result, gene_taxonomy_records, gene_taxonomy_summary)
        self._process_interproscan(result, interproscan_summary)
        self._apply_jelly_roll_summary(result, jelly_roll_summary)

        # Track HHG evidence for decision weighting
        has_hhg_evidence = "hhg" in getattr(refined_boundary, "seed_sources", [])

        # ==================================================
        # Run Full Analysis (no gates - all EVEs get full evaluation)
        # ==================================================
        logger.info(
            f"{result.eve_id}: Running full evidence analysis "
            f"(boundary_confidence={refined_boundary.confidence:.3f})"
        )

        self._run_tiebreakers(
            result, refined_boundary, window_features,
            hallmark_hits, novelty_scores, porf_sequences, precomputed_tmvec
        )

        # Phylogenetic rejection override
        if (
            self.config.phylogenetic_rejection_override
            and result.phylogenetic_result
            and result.phylogenetic_score < self.config.phylogenetic_rejection_threshold
        ):
            if result.phylogenetic_result.has_viral_in_any_neighbor:
                logger.info(
                    f"{result.eve_id}: NOT rejected despite {result.gvclass_domain} domain - "
                    f"viral signal found in top 5 neighbors (score={result.phylogenetic_score:.3f})"
                )
            elif result.phylogenetic_result.rejects_viral:
                logger.info(
                    f"{result.eve_id}: LOW confidence (phylogenetic rejection) "
                    f"(domain={result.gvclass_domain}, score={result.phylogenetic_score:.3f})"
                )
                result.status = VerificationStatus.LOW_CONFIDENCE_TIEBREAKER
                result.confidence_tier = assign_confidence_tier(
                    result,
                    high=self.config.high_tier_threshold,
                    low=self.config.low_tier_threshold,
                )
                result.likely_family = infer_likely_family(result)
                return result

        # Override decisions
        if (
            (result.gene_taxonomy_has_ncldv_mirus or result.gene_taxonomy_has_vp_plv)
            and result.hallmark_count >= 2
            and refined_boundary.confidence >= 0.15
        ):
            # Defensive: ensure gene_count is synced with gene_taxonomy_total
            if result.gene_count == 0 and result.gene_taxonomy_total > 0:
                result.gene_count = result.gene_taxonomy_total

            result.final_confidence = calculate_eve_confidence(
                result=result,
                crf_confidence=refined_boundary.confidence,
                tmvec_score=None,
                host_signature_score_threshold=self.config.host_signature_score_threshold,
                ablation_id=self.config.ablation_id,
                high_tier_threshold=self.config.high_tier_threshold,
                low_tier_threshold=self.config.low_tier_threshold,
            )
            result.confidence_tier = assign_confidence_tier(
                result,
                high=self.config.high_tier_threshold,
                low=self.config.low_tier_threshold,
            )

            # Sync status with confidence_tier
            if result.confidence_tier == "HIGH":
                result.status = VerificationStatus.HIGH_CONFIDENCE
            elif result.confidence_tier == "MEDIUM":
                result.status = VerificationStatus.MEDIUM_CONFIDENCE
            else:
                result.status = VerificationStatus.LOW_CONFIDENCE_TIEBREAKER

            logger.info(
                f"{result.eve_id}: {result.confidence_tier} confidence via marker+taxonomy override "
                f"(markers={result.hallmark_count}, ncldv_mirus={result.gene_taxonomy_ncldv_top10}, genes={result.gene_count}, score={result.final_confidence:.3f})"
            )
            result.likely_family = infer_likely_family(result)
            return result

        if should_accept_mcp_override(result):
            # Defensive: ensure gene_count is synced with gene_taxonomy_total
            if result.gene_count == 0 and result.gene_taxonomy_total > 0:
                result.gene_count = result.gene_taxonomy_total

            result.final_confidence = calculate_eve_confidence(
                result=result,
                crf_confidence=refined_boundary.confidence,
                tmvec_score=None,
                host_signature_score_threshold=self.config.host_signature_score_threshold,
                ablation_id=self.config.ablation_id,
                high_tier_threshold=self.config.high_tier_threshold,
                low_tier_threshold=self.config.low_tier_threshold,
            )
            result.confidence_tier = assign_confidence_tier(
                result,
                high=self.config.high_tier_threshold,
                low=self.config.low_tier_threshold,
            )

            # Sync status with confidence_tier
            if result.confidence_tier == "HIGH":
                result.status = VerificationStatus.HIGH_CONFIDENCE
            elif result.confidence_tier == "MEDIUM":
                result.status = VerificationStatus.MEDIUM_CONFIDENCE
            else:
                result.status = VerificationStatus.LOW_CONFIDENCE_TIEBREAKER

            logger.info(
                f"{result.eve_id}: {result.confidence_tier} confidence via MCP override "
                f"(markers={result.hallmark_count}, mcp={result.has_mcp}, genes={result.gene_count}, score={result.final_confidence:.3f})"
            )
            result.likely_family = infer_likely_family(result)
            return result

        # Final weighted decision
        self._calculate_final_decision(result, refined_boundary, has_hhg_evidence)
        result.confidence_tier = assign_confidence_tier(
            result,
            high=self.config.high_tier_threshold,
            low=self.config.low_tier_threshold,
        )
        result.likely_family = infer_likely_family(result)
        return result

    def _run_phylogenetic_validation(
        self,
        result: VerificationResult,
        refined_boundary: RefinedBoundary,
    ) -> None:
        """
        Run phylogenetic validation (GVClass + Diamond) on a region.

        Updates the result object in place with phylogenetic scores.

        Args:
            result: VerificationResult to update
            refined_boundary: The refined boundary to validate
        """
        if not self.phylogenetic_validator:
            return

        try:
            phylo_result = self.phylogenetic_validator.validate_eve(
                eve_id=result.eve_id,
                scaffold=refined_boundary.scaffold,
                start=refined_boundary.start,
                end=refined_boundary.end,
                run_gvclass=True,
                run_diamond=self.config.diamond_db is not None,
            )

            result.phylogenetic_result = phylo_result
            result.phylogenetic_score = phylo_result.combined_score
            result.has_phylogenetic_support = phylo_result.supports_viral
            # Stage 1A: has_mcp may only be PROMOTED by phylogenetic validation,
            # never demoted. GVClass's MCP marker corpus can disagree with the
            # HMM/Diamond-based detection used upstream; the upstream True must
            # survive GVClass reporting False.
            result.has_mcp = result.has_mcp or phylo_result.has_mcp
            result.is_chimeric = phylo_result.is_chimeric

            # Extract GVClass details
            if phylo_result.gvclass:
                result.gvclass_domain = phylo_result.gvclass.domain
                result.gvclass_percent = phylo_result.gvclass.domain_percent

                # Use GVClass domain for taxonomy prediction
                if phylo_result.gvclass.is_viral:
                    result.predicted_taxonomy = phylo_result.gvclass.domain
                    result.taxonomy_confidence = phylo_result.gvclass.confidence_score

            # Extract Diamond details
            if phylo_result.diamond:
                result.diamond_domain = phylo_result.diamond.best_domain
                result.diamond_percent = phylo_result.diamond.best_domain_percent

            logger.info(
                f"{result.eve_id}: Phylogenetic validation - "
                f"GVClass={result.gvclass_domain} ({result.gvclass_percent:.1f}%), "
                f"Diamond={result.diamond_domain} ({result.diamond_percent:.1f}%), "
                f"score={result.phylogenetic_score:.3f}"
            )

        except Exception as e:
            logger.warning(f"Phylogenetic validation failed for {result.eve_id}: {e}")

    def verify_batch(
        self,
        refined_boundaries: list[RefinedBoundary],
        all_window_features: list[list],
        all_hallmark_hits: Optional[list[list]] = None,
        all_novelty_scores: Optional[list[dict]] = None,
        all_porf_sequences: Optional[list[list[tuple[str, str]]]] = None,
    ) -> list[VerificationResult]:
        """
        Verify multiple EVE candidates.

        This method performs batch TMVec processing ONCE for all proteins across
        all EVEs, then distributes the results to individual EVE verification.
        This is much more efficient than running TMVec per-EVE because:
        1. Model is loaded only once
        2. All proteins are embedded in batches
        3. Database searches are performed sequentially

        Args:
            refined_boundaries: List of RefinedBoundary from Phase 2
            all_window_features: Window features for each EVE
            all_hallmark_hits: Hallmark hits for each EVE
            all_novelty_scores: Novelty scores for each EVE
            all_porf_sequences: pORF sequences for each EVE

        Returns:
            List of VerificationResult objects
        """
        # Run batch TMVec search for all proteins across all EVEs
        precomputed_tmvec: Optional[dict[str, dict]] = None
        if self.config.use_tmvec_database and all_porf_sequences:
            # Collect proteins and deduplicate by pORF ID.
            raw_protein_count = 0
            conflicting_ids = 0
            protein_by_id: dict[str, str] = {}
            for porf_sequences in all_porf_sequences:
                if porf_sequences:
                    raw_protein_count += len(porf_sequences)
                    for porf_id, sequence in porf_sequences:
                        prior = protein_by_id.get(porf_id)
                        if prior is None:
                            protein_by_id[porf_id] = sequence
                        elif prior != sequence:
                            conflicting_ids += 1

            all_proteins = list(protein_by_id.items())

            if all_proteins:
                logger.info(
                    f"TMVec batch: Processing {raw_protein_count} proteins from "
                    f"{len(refined_boundaries)} EVEs"
                )
                if conflicting_ids:
                    logger.warning(
                        "TMVec batch: %d duplicate pORF IDs had conflicting sequences; using first occurrence",
                        conflicting_ids,
                    )
                searcher = self.tmvec_searcher
                if searcher:
                    precomputed_tmvec = searcher.search_batch(all_proteins)
                    logger.info(
                        f"TMVec batch: Completed {len(precomputed_tmvec)} protein searches"
                    )

        results = []

        for i, boundary in enumerate(refined_boundaries):
            window_features = all_window_features[i] if i < len(all_window_features) else []
            hallmark_hits = all_hallmark_hits[i] if all_hallmark_hits and i < len(all_hallmark_hits) else None
            novelty_scores = all_novelty_scores[i] if all_novelty_scores and i < len(all_novelty_scores) else None
            porf_sequences = all_porf_sequences[i] if all_porf_sequences and i < len(all_porf_sequences) else None

            result = self.verify_eve(
                refined_boundary=boundary,
                window_features=window_features,
                hallmark_hits=hallmark_hits,
                novelty_scores=novelty_scores,
                porf_sequences=porf_sequences,
                precomputed_tmvec=precomputed_tmvec,
            )
            results.append(result)

        return results


def synthesize_evidence(
    refined_boundaries: list[RefinedBoundary],
    window_features_list: list[list],
    genome_path: Optional[Path] = None,
    hallmark_hits_list: Optional[list[list]] = None,
    work_dir: Optional[Path] = None,
    viral_structure_db: Optional[Path] = None,
    gvclass_db: Optional[Path] = None,
    diamond_db: Optional[Path] = None,
    enable_phylogenetic: bool = False,
    use_gpu: bool = True,
    threads: int = 8,
    ablation_id: AblationID = AblationID.A0,
) -> list[VerificationResult]:
    """
    Main entry point for Phase 3 evidence synthesis.

    This function runs full evidence synthesis for each Phase 2 boundary,
    including evidence coherence, optional structural homology, and optional
    phylogenetic validation.

    Args:
        refined_boundaries: RefinedBoundary objects from Phase 2
        window_features_list: Window features for each boundary
        genome_path: Path to input genome FASTA (needed for phylogenetic validation)
        hallmark_hits_list: Hallmark hits for each boundary
        work_dir: Working directory
        viral_structure_db: Path to viral structure DB for FoldSeek
        gvclass_db: Path to GVClass database (optional)
        diamond_db: Path to Diamond database for BLASTp validation
        use_gpu: Whether to use GPU for structural analysis
        threads: Number of threads
        ablation_id: Mutually exclusive benchmark ablation arm

    Returns:
        List of VerificationResult objects with final verdicts
    """
    logger.info("=" * 60)
    logger.info("Phase 3: Evidence Synthesis")
    logger.info("=" * 60)

    if not refined_boundaries:
        logger.info("No boundaries to verify")
        return []

    # Configure synthesizer
    config = EvidenceSynthesizerConfig(
        use_boltz=viral_structure_db is not None,
        use_phylogenetic_validation=enable_phylogenetic,
        gvclass_db=gvclass_db,
        diamond_db=diamond_db,
        device="cuda" if use_gpu else "cpu",
        threads=threads,
        ablation_id=ablation_id,
    )

    synthesizer = EvidenceSynthesizer(
        config=config,
        viral_structure_db=viral_structure_db,
        genome_path=genome_path,
        work_dir=work_dir,
    )

    # Run verification
    results = synthesizer.verify_batch(
        refined_boundaries=refined_boundaries,
        all_window_features=window_features_list,
        all_hallmark_hits=hallmark_hits_list,
    )

    # Summary
    accepted = sum(1 for r in results if r.is_accepted)
    high_conf = sum(
        1 for r in results
        if r.status == VerificationStatus.HIGH_CONFIDENCE
    )
    medium_conf = sum(
        1 for r in results
        if r.status == VerificationStatus.MEDIUM_CONFIDENCE
    )
    low_conf = len(results) - accepted

    # Phylogenetic summary
    with_phylo = sum(1 for r in results if r.phylogenetic_result is not None)
    phylo_support = sum(1 for r in results if r.has_phylogenetic_support)

    logger.info("Evidence synthesis complete:")
    logger.info(f"  Total candidates: {len(results)}")
    logger.info(f"  High/Medium confidence: {accepted}")
    logger.info(f"    - High confidence: {high_conf}")
    logger.info(f"    - Medium confidence: {medium_conf}")
    logger.info(f"  Low confidence: {low_conf}")
    logger.info("  Phylogenetic validation:")
    logger.info(f"    - Regions validated: {with_phylo}")
    logger.info(f"    - With phylo support: {phylo_support}")

    return results
