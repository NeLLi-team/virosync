"""Versioned public coordinate, class, and output conventions."""

COORDINATE_SCHEMA_VERSION = 2
OUTPUT_SCHEMA_VERSION = 4
COORDINATE_CONVENTION = "0-based, half-open [start, end)"

EFFECTIVE_EVE_CLASSES = (
    "NCLDV",
    "MIRUS",
    "PPV",
    "CRESS",
    "MIXED",
    "UNKNOWN",
)
CONCRETE_EVE_CLASSES = frozenset({"NCLDV", "MIRUS", "PPV", "CRESS"})
EFFECTIVE_EVE_CLASS_COUNT_KEYS = {
    "NCLDV": "ncldv_count",
    "MIRUS": "mirus_count",
    "PPV": "ppv_count",
    "CRESS": "cress_count",
    "MIXED": "mixed_count",
    "UNKNOWN": "unknown_count",
}

DETAILED_TAXONOMY_PARTITION = (
    "EUK",
    "MITO",
    "PLASTID",
    "BAC",
    "ARC",
    "UNK",
    "NO_HITS",
    "NCLDV",
    "MIRUS",
    "PPV",
    "CRESS",
    "GVMAG",
    "PHAGE",
)

DETAILED_PREDICTION_COLUMNS = (
    # Identity and final calls
    "eve_id",
    "scaffold",
    "start",
    "end",
    "length",
    "confidence_tier",
    "final_confidence",
    "effective_eve_class",
    "likely_family",
    "ppv_subtype",
    "likely_group",
    # Candidate provenance
    "candidate_start",
    "candidate_end",
    "candidate_length",
    "candidate_reduction_bp",
    "candidate_reduction_reason",
    "seed_sources",
    # Marker evidence
    "hallmark_total",
    "hallmark_unique",
    "mcp_gene_ids",
    "tier1_bypassed_marker_count",
    "tier1_bypassed_marker_ids",
    "tier1_bypassed_marker_models",
    "gvogm_count",
    "gvogm_names",
    "og_count",
    "og_names",
    "gvogm_unvalidated_count",
    "gvogm_unvalidated_names",
    "og_unvalidated_count",
    "og_unvalidated_names",
    "marker_complement_score",
    "family_consistency_score",
    "seed_marker_names",
    "other_marker_names",
    "seed_marker_patterns",
    "other_marker_patterns",
    # Gene taxonomy
    "total_proteins",
    "ncldv_top10_proteins",
    "mirus_top10_proteins",
    "ppv_top10_proteins",
    "cress_top10_proteins",
    "taxonomy_best_hits",
    # Composition and host evidence
    "kfd",
    "gc_deviation",
    "region_gc_percent",
    "genome_gc_percent",
    "gc_delta",
    "host_signature_gene_count",
    "host_signature_fraction",
    "host_signature_weighted_mean",
    # InterProScan evidence
    "interproscan_total_hits",
    "interproscan_viral_hits",
    "interproscan_keyword_hits",
    "interproscan_category_score",
    "interproscan_score",
    # Marker-set completeness
    "vp_completeness",
    "ppv_completeness",
    "ncldv_completeness",
    "mirus_completeness",
)

DETAILED_PREDICTION_EXTENDED_COLUMNS = frozenset(
    {
        "seed_sources",
        "marker_complement_score",
        "family_consistency_score",
        "seed_marker_names",
        "other_marker_names",
        "seed_marker_patterns",
        "other_marker_patterns",
        "host_signature_gene_count",
        "host_signature_fraction",
        "host_signature_weighted_mean",
        "interproscan_category_score",
        "vp_completeness",
        "ppv_completeness",
        "ncldv_completeness",
        "mirus_completeness",
    }
)

# GVClass unified Polinton-like viruses and virophages into the single phylum
# Preplasmiviricota, and the v1.0.6 reference bundle followed: it carries 95,947
# ``PPV__`` labels and no ``VP__`` or ``PLV__`` at all. VP and PLV survive here
# only as read aliases so result files written before the migration still parse.
# Nothing produces them any more; canonical_family() is the one place that
# collapses them, so a single Preplasmiviricota region cannot be counted as two
# or three separate families.
PPV_LEGACY_ALIASES = frozenset({"VP", "PLV"})


def canonical_family(value: object) -> str:
    """Return the canonical lineage token for a raw class or family label.

    Legacy ``VP`` and ``PLV`` both resolve to ``PPV``. Any other token is passed
    through upper-cased and stripped, so callers can use this on marker names,
    class fields, and persisted labels alike.
    """
    token = str(value or "").strip().upper()
    return "PPV" if token in PPV_LEGACY_ALIASES else token


def normalize_effective_eve_class(value: object) -> str:
    """Normalize a persisted class into the exhaustive public partition.

    Legacy VP/PLV values are folded onto PPV, so reading an old result file
    yields the same class the current pipeline would emit for that region.
    """
    normalized = canonical_family(value)
    return normalized if normalized in EFFECTIVE_EVE_CLASSES else "UNKNOWN"


def _resolve_non_low_effective_eve_class(
    *,
    region_classification: object = "",
    classification: object = "",
    likely_family: object = "",
) -> str:
    region = canonical_family(region_classification)
    classification_value = canonical_family(classification)
    likely = canonical_family(likely_family)
    if region in CONCRETE_EVE_CLASSES:
        return region
    for value in (classification_value, likely):
        if value in CONCRETE_EVE_CLASSES:
            return value
    if "MIXED" in (region, classification_value, likely):
        return "MIXED"
    return "UNKNOWN"


def resolve_effective_eve_class(
    *,
    confidence_tier: object = "",
    region_classification: object = "",
    classification: object = "",
    likely_family: object = "",
) -> str:
    """Resolve raw labels with the same tier-aware precedence as the v2 gate."""
    tier = str(confidence_tier or "").strip().upper()
    if tier == "LOW":
        family_value = ""
        for value in (likely_family, classification):
            candidate = str(value or "").strip().upper()
            if candidate:
                family_value = candidate
                break
        family = normalize_effective_eve_class(family_value)
        if family != "UNKNOWN":
            return family
        mixed_bridge = _resolve_non_low_effective_eve_class(
            region_classification=region_classification,
            classification=classification,
            likely_family=likely_family,
        )
        return "MIXED" if mixed_bridge == "MIXED" else "UNKNOWN"
    return _resolve_non_low_effective_eve_class(
        region_classification=region_classification,
        classification=classification,
        likely_family=likely_family,
    )


def empty_effective_eve_class_counts() -> dict[str, int]:
    """Return ordered zero counts for the exhaustive class partition."""
    return {key: 0 for key in EFFECTIVE_EVE_CLASS_COUNT_KEYS.values()}


def normalize_effective_eve_class_counts(
    counts: dict[str, int],
) -> dict[str, int]:
    """Fold legacy VP/PLV count keys into the current public partition."""
    normalized = {eve_class: 0 for eve_class in EFFECTIVE_EVE_CLASSES}
    current_keys = set(EFFECTIVE_EVE_CLASSES)
    legacy_keys = {
        "NCLDV",
        "VP",
        "PLV",
        "MIRUS",
        "MIXED",
        "PPV",
        "UNKNOWN",
    }
    if set(counts) not in (current_keys, legacy_keys):
        raise ValueError(
            "class_counts must contain the current or legacy complete "
            "effective EVE class partition"
        )
    for eve_class, count in counts.items():
        normalized[normalize_effective_eve_class(eve_class)] += int(count)
    return normalized


def effective_eve_class_count_total(summary: dict) -> int:
    """Sum the six exclusive public class counts in a result mapping."""
    return sum(int(summary.get(key, 0) or 0) for key in EFFECTIVE_EVE_CLASS_COUNT_KEYS.values())


def coordinate_contract_metadata() -> dict[str, int | str]:
    """Return the exact coordinate contract embedded in output metadata."""
    return {
        "coordinate_schema_version": COORDINATE_SCHEMA_VERSION,
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "coordinate_convention": COORDINATE_CONVENTION,
    }
