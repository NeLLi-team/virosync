"""Versioned public coordinate, class, and output conventions."""

COORDINATE_SCHEMA_VERSION = 2
OUTPUT_SCHEMA_VERSION = 6
COORDINATE_CONVENTION = "0-based, half-open [start, end)"

# The PUBLISHED class partition. Every surface that reports an EVE's class to a
# reader uses these tokens. MIXED was retired here in schema 5: a region whose
# markers disagree is VIRAL_UNKNOWN, and MIXED survives only as a read alias.
EFFECTIVE_EVE_CLASSES = (
    "NCLDV",
    "MIRUS",
    "PPV",
    "CRESS",
    "PHAGE",
    "VIRAL_UNKNOWN",
    "UNKNOWN",
)
CONCRETE_EVE_CLASSES = frozenset({"NCLDV", "MIRUS", "PPV", "CRESS"})
# Published classes that name a lineage. Distinct from CONCRETE_EVE_CLASSES,
# which is the acceptance gate's own set and has no PHAGE: the gate never sees
# a PHAGE label, but taxonomy consensus and ANI propagation both publish one.
# VIRAL_UNKNOWN and UNKNOWN name the absence of a decision and are excluded.
LINEAGE_EVE_CLASSES = CONCRETE_EVE_CLASSES | frozenset({"PHAGE"})
EFFECTIVE_EVE_CLASS_COUNT_KEYS = {
    "NCLDV": "ncldv_count",
    "MIRUS": "mirus_count",
    "PPV": "ppv_count",
    "CRESS": "cress_count",
    "PHAGE": "phage_count",
    "VIRAL_UNKNOWN": "viral_unknown_count",
    "UNKNOWN": "unknown_count",
}
# Count keys that are no longer written but must still parse when a persisted
# summary is read back, each folded onto its current key.
LEGACY_EVE_CLASS_COUNT_KEYS = {
    "mixed_count": "viral_unknown_count",
    "vp_count": "ppv_count",
    "plv_count": "ppv_count",
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
    "canonical_selection_outcome",
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
    # Per-genome ANI clustering of accepted EVEs, and the class an MCP-bearing
    # cluster member propagated to this region
    "ani_cluster_id",
    "ani_cluster_size",
    "ani_max_percent",
    "taxonomy_class_before_ani",
    "taxonomy_class_propagated_from",
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
    """Normalize a persisted class into the exhaustive published partition.

    Legacy VP/PLV values are folded onto PPV and legacy MIXED onto
    VIRAL_UNKNOWN, so reading an old result file yields the same class the
    current pipeline would emit for that region.
    """
    normalized = canonical_family(value)
    if normalized == "MIXED":
        return "VIRAL_UNKNOWN"
    return normalized if normalized in EFFECTIVE_EVE_CLASSES else "UNKNOWN"


# The v2 acceptance gate has its own class vocabulary, deliberately separate
# from the published partition: it branches on MIXED, and folding MIXED onto
# VIRAL_UNKNOWN here would drop every MIXED region out of the gate's accepting
# branches. PHAGE and VIRAL_UNKNOWN are publication-only tokens and never reach
# the gate, so they normalize to UNKNOWN.
_GATE_EVE_CLASSES = frozenset(CONCRETE_EVE_CLASSES | {"MIXED", "UNKNOWN"})


def _normalize_gate_eve_class(value: object) -> str:
    """Normalize a raw label into the v2 acceptance gate's own vocabulary."""
    normalized = canonical_family(value)
    return normalized if normalized in _GATE_EVE_CLASSES else "UNKNOWN"


def _resolve_non_low_effective_eve_class(
    *,
    region_classification: object = "",
    classification: object = "",
    likely_family: object = "",
) -> str:
    region = _normalize_gate_eve_class(region_classification)
    classification_value = _normalize_gate_eve_class(classification)
    likely = _normalize_gate_eve_class(likely_family)
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
    """Resolve raw labels with the same tier-aware precedence as the v2 gate.

    Returns a GATE class, which can be MIXED. Publication surfaces must pass the
    result through :func:`normalize_effective_eve_class`.
    """
    tier = str(confidence_tier or "").strip().upper()
    if tier == "LOW":
        family_value = ""
        for value in (likely_family, classification):
            candidate = str(value or "").strip().upper()
            if candidate:
                family_value = candidate
                break
        family = _normalize_gate_eve_class(family_value)
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
    """Fold legacy VP/PLV/MIXED count keys into the current public partition."""
    normalized = {eve_class: 0 for eve_class in EFFECTIVE_EVE_CLASSES}
    current_keys = set(EFFECTIVE_EVE_CLASSES)
    # Schema 4 published MIXED; the partition before the Preplasmiviricota
    # migration split VP and PLV and had no CRESS.
    schema4_keys = {"NCLDV", "MIRUS", "PPV", "CRESS", "MIXED", "UNKNOWN"}
    legacy_keys = {
        "NCLDV",
        "VP",
        "PLV",
        "MIRUS",
        "MIXED",
        "PPV",
        "UNKNOWN",
    }
    if set(counts) not in (current_keys, schema4_keys, legacy_keys):
        raise ValueError(
            "class_counts must contain the current or legacy complete "
            "effective EVE class partition"
        )
    for eve_class, count in counts.items():
        normalized[normalize_effective_eve_class(eve_class)] += int(count)
    return normalized


def effective_eve_class_count_total(summary: dict) -> int:
    """Sum the exclusive published class counts in a result mapping."""
    return sum(int(summary.get(key, 0) or 0) for key in EFFECTIVE_EVE_CLASS_COUNT_KEYS.values())


def coordinate_contract_metadata() -> dict[str, int | str]:
    """Return the exact coordinate contract embedded in output metadata."""
    return {
        "coordinate_schema_version": COORDINATE_SCHEMA_VERSION,
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "coordinate_convention": COORDINATE_CONVENTION,
    }
