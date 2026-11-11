"""Host signature modeling for EUK-like marker hits.

Builds a weighted token model from eligible marker records using their first ten
DIAMOND matches. Each full-protein query contributes once. The caller excludes
proteins with viral marker-search evidence or reported capsids before host search.
Used for host-like scoring and outward extension stopping.
"""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from math import log, sqrt
from pathlib import Path

from virosync.pipeline.taxonomy_utils import (
    EUK_PEVE_PREFIX,
    TaxonomyFingerprint,
    compute_hit_weight,
    is_assigned_level,
    resolve_org_id,
)

_TAXONOMY_LOOKUP: TaxonomyLabelLookup | None = None
HOST_TOKEN_CACHE_ATTRIBUTE = "_token_cache"

# Process-level cache so the large taxonomy-label TSV (≈168k rows) is parsed
# once per (path, size, mtime) instead of on every phase/task that needs it.
_LABEL_CACHE: dict = {}

_CELLULAR_DOMAINS = {
    "euk": "eukaryota",
    "eukaryota": "eukaryota",
    "bac": "bacteria",
    "bacteria": "bacteria",
    "arc": "archaea",
    "archaea": "archaea",
}
_RANK_FILLERS = frozenset({"na", "singleton"})


class TaxonomyLabelLookup(dict):
    """Lookup table for taxonomy labels keyed by target prefix."""

    @classmethod
    def load(cls, path) -> TaxonomyLabelLookup:
        path = Path(path)
        cache_key = None
        try:
            st = path.stat()
            cache_key = (str(path.resolve()), st.st_size, int(st.st_mtime))
        except OSError:
            cache_key = None
        if cache_key is not None and cache_key in _LABEL_CACHE:
            return _LABEL_CACHE[cache_key]
        data = cls()
        with open(path) as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                key = parts[0].strip()
                label = parts[1].strip()
                if key:
                    data[key] = label
        if cache_key is not None:
            _LABEL_CACHE[cache_key] = data
        return data


def set_taxonomy_lookup(lookup: TaxonomyLabelLookup | None) -> None:
    """Set global taxonomy lookup used by downstream helpers."""
    global _TAXONOMY_LOOKUP
    _TAXONOMY_LOOKUP = lookup


def get_taxonomy_lookup() -> TaxonomyLabelLookup | None:
    """Return the taxonomy lookup set for this process, if any."""
    return _TAXONOMY_LOOKUP


@dataclass
class HostSignatureModel:
    token_weights: dict[str, float] = field(default_factory=dict)
    token_counts: dict[str, int] = field(default_factory=dict)
    token_bits: dict[str, list[float]] = field(default_factory=dict)
    max_weight: float = 0.0
    min_token_length: int = 3
    host_prefixes: list[str] = field(default_factory=list)
    weight_mode: str = "rank"
    # Full labels and hit weights retained as model-construction provenance.
    lineage_weights: dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> dict:
        payload = {
            "token_weights": self.token_weights,
            "token_counts": self.token_counts,
            "lineage_weights": self.lineage_weights,
            "max_weight": self.max_weight,
            "min_token_length": self.min_token_length,
        }
        if self.token_bits:
            payload["token_bits"] = self.token_bits
        if self.host_prefixes:
            payload["host_prefixes"] = list(self.host_prefixes)
        if self.weight_mode:
            payload["weight_mode"] = self.weight_mode
        return payload

    @classmethod
    def from_dict(cls, payload: dict | None) -> HostSignatureModel:
        if not payload:
            return cls()
        return cls(
            token_weights=dict(payload.get("token_weights") or {}),
            token_counts=dict(payload.get("token_counts") or {}),
            token_bits=dict(payload.get("token_bits") or {}),
            lineage_weights=dict(payload.get("lineage_weights") or {}),
            max_weight=float(payload.get("max_weight") or 0.0),
            min_token_length=int(payload.get("min_token_length") or 3),
            host_prefixes=list(payload.get("host_prefixes") or []),
            weight_mode=str(payload.get("weight_mode") or "rank"),
        )


def is_host_marker(hmm_target: str) -> bool:
    """Whether an HMM model is a cellular housekeeping marker (COG or BUSCO eukaryote)."""
    if not hmm_target:
        return False
    return hmm_target.startswith("COG") or "at2759" in hmm_target


def _split_csv_field(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [str(item) for item in value]
    text = str(value).strip()
    if not text:
        return []
    return text.split(",")


def _parse_float_list(values: list[str]) -> list[float]:
    parsed = []
    for item in values:
        try:
            parsed.append(float(item))
        except (TypeError, ValueError):
            parsed.append(0.0)
    return parsed


def _parse_top10_hits(hit) -> list[tuple[str, float, float, float]]:
    if isinstance(hit, dict):
        targets = _split_csv_field(hit.get("top10_targets"))
        bits = _parse_float_list(_split_csv_field(hit.get("top10_bits", hit.get("top10_bitscores"))))
        pidents = _parse_float_list(_split_csv_field(hit.get("top10_pidents")))
        evalues = _parse_float_list(_split_csv_field(hit.get("top10_evalues")))
        best_target = hit.get("best_hit_target", "")
        best_bits = hit.get("best_hit_bits", 0.0)
    else:
        targets = _split_csv_field(getattr(hit, "top10_targets", ""))
        # GeneTaxonomy stores top10_bits; ValidatedMarkerHit stores top10_bitscores.
        bits = _parse_float_list(_split_csv_field(getattr(hit, "top10_bits", getattr(hit, "top10_bitscores", ""))))
        pidents = _parse_float_list(_split_csv_field(getattr(hit, "top10_pidents", "")))
        evalues = _parse_float_list(_split_csv_field(getattr(hit, "top10_evalues", "")))
        best_target = getattr(hit, "best_hit_target", "")
        best_bits = getattr(hit, "best_hit_bits", 0.0)
    if not targets and best_target:
        return [(str(best_target), float(best_bits or 0.0), 0.0, 1.0)]
    max_len = max(len(targets), len(bits), len(pidents), len(evalues), 0)
    while len(bits) < max_len:
        bits.append(0.0)
    while len(pidents) < max_len:
        pidents.append(0.0)
    while len(evalues) < max_len:
        evalues.append(1.0)
    return [(targets[i], bits[i], pidents[i], evalues[i]) for i in range(min(max_len, len(targets)))]


def _parse_top10_prefixes(hit) -> list[str]:
    if isinstance(hit, dict):
        prefixes = _split_csv_field(hit.get("top10_prefixes"))
    else:
        prefixes = _split_csv_field(getattr(hit, "top10_prefixes", ""))
    return prefixes


def _lineage_tokens(lineage: str, min_token_length: int) -> list[str]:
    levels = [level.strip().lower() for level in lineage.split("|")]
    domains = {_CELLULAR_DOMAINS[level] for level in levels if level in _CELLULAR_DOMAINS}
    if len(domains) != 1 or any(level.endswith("-peve") for level in levels):
        return []
    if any(
        level in {"mito", "plastid", "vir", "ncldv", "mirus", "phage", "ppv", "plv", "vp", "cress"} for level in levels
    ):
        return []
    return list(
        dict.fromkeys(
            [
                *domains,
                *(
                    level
                    for level in levels
                    if level not in _CELLULAR_DOMAINS
                    and is_assigned_level(level, min_token_length)
                    and level not in _RANK_FILLERS
                ),
            ]
        )
    )


def _tokens_for_target(target: str, taxonomy_lookup: dict | None, min_token_length: int) -> list[str]:
    # Both identities are required: a bad label must not turn a viral target into host.
    namespace = target.split("__", 1)[0].lower()
    if namespace not in _CELLULAR_DOMAINS or not taxonomy_lookup:
        return []
    lineage = taxonomy_lookup.get(resolve_org_id(target, taxonomy_lookup), "")
    return _lineage_tokens(lineage, min_token_length)


def _append_weights(
    token_bits: dict[str, list[float]],
    token: str,
    weight: float,
    max_bits_per_token: int,
) -> None:
    if max_bits_per_token <= 0:
        return
    token_bits.setdefault(token, []).append(weight)
    if len(token_bits[token]) <= max_bits_per_token:
        return
    token_bits[token].sort(reverse=True)
    del token_bits[token][max_bits_per_token:]


def build_host_signature_model(
    marker_hits: Iterable[Mapping],
    validated_prefixes: set[str] | None = None,
    supporting_prefixes: set[str] | None = None,
    min_token_length: int = 3,
    host_prefixes: set[str] | None = None,
    max_bits_per_token: int = 10,
    weight_mode: str = "rank",
) -> HostSignatureModel:
    """Build a weighted host-signature model from marker hits.

    Uses each record with no excluded viral/supporting top-ten prefix.
    Token weights are weighted sums of taxonomy substrings (split at |).
    """
    validated = validated_prefixes or {
        # PPV__ is the unified Preplasmiviricota prefix and carries 95,947 of the
        # v1.0.6 labels; PLV__/VP__ are kept for pre-migration bundles.
        "NCLDV__",
        "MIRUS__",
        "PPV__",
        "PLV__",
        "VP__",
        "CRESS__",
        "PHAGE__",
        EUK_PEVE_PREFIX,
    }
    supporting = supporting_prefixes or set()
    host_prefixes = host_prefixes or {"EUK__"}

    taxonomy_lookup = _TAXONOMY_LOOKUP
    token_weights: dict[str, float] = {}
    token_counts: dict[str, int] = {}
    token_bits: dict[str, list[float]] = {}
    lineage_weights: dict[str, float] = {}
    seen_queries: set[tuple[str, str]] = set()

    for hit in marker_hits:
        top10_prefixes = _parse_top10_prefixes(hit)[:10]
        if any(p in validated for p in top10_prefixes):
            continue
        if supporting and any(p in supporting for p in top10_prefixes):
            continue
        if host_prefixes and not any(p in host_prefixes for p in top10_prefixes):
            continue
        query = str(_get_record_attr(hit, "query_porf", "")).split("|aa", 1)[0]
        query_key = (str(_get_record_attr(hit, "scaffold", "")), query)
        if query and query_key in seen_queries:
            continue
        if query:
            seen_queries.add(query_key)
        # Eligibility belongs to the record; do not drop cellular noise or rerank hits.
        for rank, (target, bits, _, _) in enumerate(_parse_top10_hits(hit)[:10]):
            weight = compute_hit_weight(rank, bits, weight_mode)
            tokens = _tokens_for_target(str(target), taxonomy_lookup, min_token_length)
            if weight <= 0 or not tokens:
                continue
            lineage = taxonomy_lookup[resolve_org_id(str(target), taxonomy_lookup)]
            key = lineage_key(lineage, min_token_length)
            lineage_weights[key] = lineage_weights.get(key, 0.0) + weight
            for token in tokens:
                token_weights[token] = token_weights.get(token, 0.0) + weight
                token_counts[token] = token_counts.get(token, 0) + 1
                _append_weights(token_bits, token, weight, max_bits_per_token)

    max_weight = max(token_weights.values()) if token_weights else 0.0

    return HostSignatureModel(
        token_weights=token_weights,
        token_counts=token_counts,
        token_bits=token_bits,
        max_weight=max_weight,
        min_token_length=min_token_length,
        host_prefixes=sorted(host_prefixes),
        weight_mode=weight_mode,
        lineage_weights=lineage_weights,
    )


def _get_record_attr(record, name: str, default=None):
    if isinstance(record, dict):
        return record.get(name, default)
    return getattr(record, name, default)


def score_host_signature(fingerprint: TaxonomyFingerprint, model: HostSignatureModel) -> float:
    """Compare filtered cellular token vectors by cosine similarity.

    Call host_lineage_share for records: a generic fingerprint cannot recover
    namespace exclusions or the retained-hit denominator.
    """
    gene = _canonical_token_weights(fingerprint.weighted_tokens)
    host, host_norm = _host_distribution(model)
    norm = sqrt(sum(weight * weight for weight in gene.values())) * host_norm
    return sum(weight * host.get(token, 0.0) for token, weight in gene.items()) / norm if norm else 0.0


def lineage_key(lineage: str, min_token_length: int) -> str:
    """Record the retained cellular taxonomy names for construction provenance."""
    return "|".join(_lineage_tokens(lineage, min_token_length))


def host_token_weights(model: HostSignatureModel) -> dict[str, float]:
    """Return learned tokens, or configured tokens when no cellular model exists."""
    return _host_distribution(model)[0]


def _canonical_token_weights(tokens: Mapping[str, float]) -> dict[str, float]:
    """Use the same named-token and cellular-domain aliases in both vectors."""
    weights: dict[str, float] = {}
    for token, weight in tokens.items():
        token = token.strip().lower()
        if weight <= 0 or token in _RANK_FILLERS or token in {"mito", "plastid"} or token.endswith("-peve"):
            continue
        if token not in _CELLULAR_DOMAINS and not is_assigned_level(token, 0):
            continue
        token = _CELLULAR_DOMAINS.get(token, token)
        weights[token] = weights.get(token, 0.0) + weight
    return weights


def _host_distribution(model: HostSignatureModel) -> tuple[dict[str, float], float]:
    configured = os.environ.get("VIROSYNC_HOST_LINEAGE", "")
    fallback_key = (configured, model.min_token_length) if configured else ("", 0)
    cache = getattr(model, HOST_TOKEN_CACHE_ATTRIBUTE, None)
    # Snapshot equality catches in-place edits without rebuilding the vector and
    # its norm in Python for every boundary gene. The model is fixed during a run.
    if cache is not None and cache[0] == model.token_weights and cache[1] == fallback_key:
        return cache[2], cache[3]
    weights = _canonical_token_weights(model.token_weights)
    if not weights and configured:
        weights = dict.fromkeys(_lineage_tokens(configured, model.min_token_length), 1.0)
    norm = sqrt(sum(weight * weight for weight in weights.values()))
    setattr(model, HOST_TOKEN_CACHE_ATTRIBUTE, (dict(model.token_weights), fallback_key, weights, norm))
    return weights, norm


def host_stop_available(model: HostSignatureModel | None) -> bool:
    """Whether learned or configured cellular tokens can support a host stop."""
    return model is not None and bool(host_token_weights(model))


def host_lineage_share(record, model: HostSignatureModel) -> float:
    """Return token cosine similarity scaled by retained cellular hit mass.

    Keep original hit ranks, including viral, pEVE, organelle and unknown weight
    in the denominator. This score is a similarity, not a host probability.
    """
    weights: dict[str, float] = {}
    total = assigned = 0.0
    for rank, (target, bits, _, _) in enumerate(_parse_top10_hits(record)[:10]):
        # Empty padding keeps its rank position but is not a retained hit.
        # An unidentified hit with a recorded bitscore still contributes mass.
        if not target.strip() and bits <= 0:
            continue
        weight = compute_hit_weight(rank, bits, model.weight_mode)
        if weight <= 0:
            continue
        total += weight
        tokens = _tokens_for_target(str(target), _TAXONOMY_LOOKUP, model.min_token_length)
        if not tokens:
            continue
        assigned += weight
        for token in tokens:
            weights[token] = weights.get(token, 0.0) + weight
    if not total or not assigned:
        return 0.0
    return score_host_signature(TaxonomyFingerprint(weighted_tokens=weights), model) * assigned / total


def host_signature_density_evalue_weighted(
    gene_records: Iterable[Mapping],
    model: HostSignatureModel,
    score_threshold: float = 0.5,
    host_prefixes: list[str] | None = None,
) -> tuple[int, float, float]:
    """Compute host-like counts, mean score, and evalue-weighted mean score."""
    scores = []
    weighted_scores = []
    total_weight = 0.0
    host_like = 0
    host_prefix_set = set(host_prefixes or [])

    for record in gene_records:
        score = host_lineage_share(record, model)
        scores.append(score)
        has_host_prefix = True
        if host_prefix_set:
            top1_target = _get_record_attr(record, "top1_target", "")
            has_host_prefix = any(str(top1_target).startswith(prefix) for prefix in host_prefix_set)
            if not has_host_prefix:
                for target, _bits, _pident, _evalue in _parse_top10_hits(record):
                    if any(str(target).startswith(prefix) for prefix in host_prefix_set):
                        has_host_prefix = True
                        break
        if score >= score_threshold and has_host_prefix:
            host_like += 1

        evalue = _get_record_attr(record, "top1_evalue", 1.0)
        try:
            evalue = float(evalue) if evalue is not None else 1.0
        except (TypeError, ValueError):
            evalue = 1.0
        if evalue <= 0:
            evalue = 1e-300
        # Use -log(evalue) for weighting: smaller evalues → larger weights
        weight = max(0.0, -log(evalue))
        total_weight += weight
        weighted_scores.append(score * weight)

    mean_score = sum(scores) / len(scores) if scores else 0.0
    weighted_mean = sum(weighted_scores) / total_weight if total_weight > 0 else 0.0
    return host_like, mean_score, weighted_mean


def summarize_host_signature_model(
    model: HostSignatureModel | None,
    top_k: int = 20,
) -> list[tuple[str, float, int]]:
    """Return top-k host signature tokens with weights and counts.

    Args:
        model: HostSignatureModel instance
        top_k: Number of top tokens to return

    Returns:
        List of (token, weight, count) sorted by weight desc.
    """
    if not model or not model.token_weights:
        return []

    sorted_items = sorted(
        model.token_weights.items(),
        key=lambda item: item[1],
        reverse=True,
    )[:top_k]

    return [(token, float(weight), int(model.token_counts.get(token, 0))) for token, weight in sorted_items]


def summarize_host_signature_bits(
    model: HostSignatureModel | None,
    top_k: int = 20,
    max_bits: int = 10,
) -> list[tuple[str, float, int, list[float]]]:
    """Return top-k tokens with weights, counts, and top weight samples."""
    if not model or not model.token_weights:
        return []

    sorted_items = sorted(
        model.token_weights.items(),
        key=lambda item: item[1],
        reverse=True,
    )[:top_k]

    summary = []
    for token, weight in sorted_items:
        bits = list(model.token_bits.get(token, []))
        bits.sort(reverse=True)
        summary.append((token, float(weight), int(model.token_counts.get(token, 0)), bits[:max_bits]))
    return summary
