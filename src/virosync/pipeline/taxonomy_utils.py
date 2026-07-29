"""
Shared taxonomy utilities for fingerprinting and aggregation.

Used by both phase1 (marker validation) and phase2 (boundary refinement).
"""

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class TaxonomyFingerprint:
    """Parsed taxonomy fingerprint for efficient comparison."""

    weighted_tokens: dict[str, float] = field(default_factory=dict)
    raw_tokens: dict[str, int] = field(default_factory=dict)

    @classmethod
    def from_string(cls, weighted_str: str, raw_str: str) -> "TaxonomyFingerprint":
        """Parse from comma-separated format."""
        weighted = {}
        for entry in weighted_str.split(",") if weighted_str else []:
            if ":" in entry:
                token, weight_str = entry.split(":", 1)
                try:
                    weighted[token.strip()] = float(weight_str)
                except ValueError:
                    continue

        raw = {}
        for entry in raw_str.split(",") if raw_str else []:
            if ":" in entry:
                token, count_str = entry.split(":", 1)
                try:
                    raw[token.strip()] = int(count_str)
                except ValueError:
                    continue

        return cls(weighted_tokens=weighted, raw_tokens=raw)

    def to_string(self) -> tuple[str, str]:
        """Serialize to comma-separated format."""
        sorted_weighted = sorted(
            self.weighted_tokens.items(), key=lambda x: x[1], reverse=True
        )[:20]
        weighted_str = ",".join(f"{t}:{w:.2f}" for t, w in sorted_weighted)

        sorted_raw = sorted(self.raw_tokens.items(), key=lambda x: x[1], reverse=True)[
            :20
        ]
        raw_str = ",".join(f"{t}:{c}" for t, c in sorted_raw)

        return weighted_str, raw_str


PLACEHOLDER_TOKENS = {
    "unclassified",
    "uncultured",
    "environmental",
    "sp",
    "sp.",
    "unknown",
    "unidentified",
    "candidatus",
    "incertae_sedis",
    "incertae",
    "sedis",
    "cf.",
    "aff.",
    "cf",
    "aff",
}


def iter_taxonomy_tokens(taxonomy_string: str, min_token_length: int) -> list[str]:
    """Split taxonomy lineage and filter placeholder/short tokens."""
    if not taxonomy_string:
        return []
    tokens = []
    for level in taxonomy_string.split("|"):
        token = level.strip()
        if len(token) < min_token_length:
            continue
        token_lower = token.lower()
        if token_lower in PLACEHOLDER_TOKENS or token_lower.endswith("_x"):
            continue
        tokens.append(token)
    return tokens


def _fallback_prefix_tokens(org_id: str, min_token_length: int) -> list[str]:
    if "__" not in org_id:
        return []
    prefix = org_id.split("__", 1)[0]
    if len(prefix) < min_token_length:
        return []
    return [prefix]


def resolve_org_id(target: str, taxonomy_lookup: dict) -> str:
    """Extract organism ID from a Diamond target and resolve against taxonomy lookup.

    Handles two header formats:
      - combined_proteome: EUK__EP00224|protein_id  -> org_id = EUK__EP00224
      - marker.faa:        EUK__EP00224_Organism_Name|protein_id -> needs stripping

    When the full first-field doesn't match the lookup, progressively strips
    trailing _Word segments until a match is found or the prefix__ base remains.
    """
    org_id = target.split("|", 1)[0]
    if not taxonomy_lookup or org_id in taxonomy_lookup:
        return org_id
    # Progressive strip: remove trailing _Word segments
    if "__" not in org_id:
        return org_id
    prefix, suffix = org_id.split("__", 1)
    parts = suffix.split("_")
    # Try removing trailing parts one at a time
    for i in range(len(parts) - 1, 0, -1):
        candidate = prefix + "__" + "_".join(parts[:i])
        if candidate in taxonomy_lookup:
            return candidate
    return org_id


def compute_hit_weight(rank: int, bits: float, weight_mode: str) -> float:
    """
    Compute per-hit weight based on rank or bitscore.

    Args:
        rank: 0-based rank in the top-10 list.
        bits: Bitscore for the hit.
        weight_mode: "rank" or "bitscore".

    Returns:
        Weight for this hit.
    """
    mode = (weight_mode or "rank").lower()
    if mode == "rank":
        return max(0.0, float(10 - rank))
    return float(bits or 0.0)


def aggregate_taxonomy_substrings(
    top10_hits: list[tuple[str, float, float, float]],
    taxonomy_lookup: Optional[dict],
    min_token_length: int = 3,
    weight_mode: str = "rank",
) -> TaxonomyFingerprint:
    """
    Aggregate taxonomy substring counts across all top-10 Diamond hits.

    Algorithm:
    1. For each hit: lookup full taxonomy lineage
    2. Split lineage by '|' and extract all levels (including kingdom)
    3. Accumulate weighted counts (rank-based or bitscore)
    4. Track raw hit counts per token
    5. Filter placeholders and short tokens

    Args:
        top10_hits: List of (target, bits, pident, evalue) tuples
        taxonomy_lookup: Dict mapping organism IDs to taxonomy lineages
        min_token_length: Minimum substring length

    Returns:
        TaxonomyFingerprint with weighted and raw token counts
    """
    if not taxonomy_lookup or not top10_hits:
        return TaxonomyFingerprint(weighted_tokens={}, raw_tokens={})

    weighted_counts = {}
    raw_counts = {}

    for rank, (target, bits, _, _) in enumerate(top10_hits):
        if not target:
            continue
        weight = compute_hit_weight(rank, bits, weight_mode)
        if weight <= 0:
            continue

        # Extract organism ID (handles marker.faa extended headers)
        org_id = resolve_org_id(target, taxonomy_lookup)

        # Lookup taxonomy
        taxonomy_string = taxonomy_lookup.get(org_id, "")
        tokens = iter_taxonomy_tokens(taxonomy_string, min_token_length)
        if not tokens:
            tokens = _fallback_prefix_tokens(org_id, min_token_length)
        if not tokens:
            continue

        for token in tokens:
            token_lower = token.lower()
            weighted_counts[token_lower] = weighted_counts.get(token_lower, 0.0) + weight
            raw_counts[token_lower] = raw_counts.get(token_lower, 0) + 1

    return TaxonomyFingerprint(weighted_tokens=weighted_counts, raw_tokens=raw_counts)


def calculate_fingerprint_overlap(
    gene_fingerprint: TaxonomyFingerprint,
    baseline_fingerprint: dict[str, float],
    min_overlap_score: float = 0.40,
) -> tuple[bool, float]:
    """
    Calculate weighted Jaccard-like overlap between gene and baseline.

    Uses min/max normalization for each token across gene and baseline weights,
    then computes overlap_sum / union_sum.

    Args:
        gene_fingerprint: TaxonomyFingerprint for gene
        baseline_fingerprint: Dict of {token: weight} from control genes
        min_overlap_score: Minimum overlap (0-1) to classify as host

    Returns:
        Tuple of (is_host_like, overlap_score)
    """
    if not gene_fingerprint.weighted_tokens or not baseline_fingerprint:
        return False, 0.0

    gene_total = sum(gene_fingerprint.weighted_tokens.values())
    host_total = sum(baseline_fingerprint.values())
    if gene_total <= 0 or host_total <= 0:
        return False, 0.0

    gene_norm = {
        k: v / gene_total for k, v in gene_fingerprint.weighted_tokens.items()
    }
    host_norm = {k: v / host_total for k, v in baseline_fingerprint.items()}

    all_tokens = set(gene_norm.keys()) | set(host_norm.keys())

    overlap_sum = 0.0
    union_sum = 0.0

    for token in all_tokens:
        gene_weight = gene_norm.get(token, 0.0)
        host_weight = host_norm.get(token, 0.0)

        overlap_sum += min(gene_weight, host_weight)
        union_sum += max(gene_weight, host_weight)

    overlap_score = overlap_sum / union_sum if union_sum > 0 else 0.0
    is_host = overlap_score >= min_overlap_score

    return is_host, overlap_score
