"""
Host signature modeling for EUK-like marker hits.

Builds a weighted token model from unvalidated marker hits whose top-10
Diamond prefixes are purely cellular. Used for host-like scoring and
boundary trimming.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import log
from pathlib import Path
from typing import Iterable, Mapping, Optional

from virosync.pipeline.taxonomy_utils import (
    TaxonomyFingerprint,
    aggregate_taxonomy_substrings,
    compute_hit_weight,
    iter_taxonomy_tokens,
    resolve_org_id,
)


_TAXONOMY_LOOKUP: Optional["TaxonomyLabelLookup"] = None

# Process-level cache so the large taxonomy-label TSV (≈168k rows) is parsed
# once per (path, size, mtime) instead of on every phase/task that needs it.
_LABEL_CACHE: dict = {}

# Domain-level tokens too generic for host signature comparison.
# Every eukaryotic gene shares "euk", so it dominates the model
# without providing discriminative signal between host and viral.
# Viral roots (ncldv, mirus, phage) are excluded to avoid
# artificially inflating overlap between host and viral genes.
DOMAIN_TOKENS = {"euk", "bac", "arc", "vir", "ncldv", "mirus", "phage"}


class TaxonomyLabelLookup(dict):
    """Lookup table for taxonomy labels keyed by target prefix."""

    @classmethod
    def load(cls, path) -> "TaxonomyLabelLookup":
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


def set_taxonomy_lookup(lookup: Optional[TaxonomyLabelLookup]) -> None:
    """Set global taxonomy lookup used by downstream helpers."""
    global _TAXONOMY_LOOKUP
    _TAXONOMY_LOOKUP = lookup


def get_taxonomy_lookup() -> Optional[TaxonomyLabelLookup]:
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

    def to_dict(self) -> dict:
        payload = {
            "token_weights": self.token_weights,
            "token_counts": self.token_counts,
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
    def from_dict(cls, payload: Optional[dict]) -> "HostSignatureModel":
        if not payload:
            return cls()
        return cls(
            token_weights=dict(payload.get("token_weights") or {}),
            token_counts=dict(payload.get("token_counts") or {}),
            token_bits=dict(payload.get("token_bits") or {}),
            max_weight=float(payload.get("max_weight") or 0.0),
            min_token_length=int(payload.get("min_token_length") or 3),
            host_prefixes=list(payload.get("host_prefixes") or []),
            weight_mode=str(payload.get("weight_mode") or "rank"),
        )


def _is_host_marker(hmm_target: str) -> bool:
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
    return [item for item in text.split(",") if item]


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
        bits = _parse_float_list(_split_csv_field(hit.get("top10_bitscores")))
        pidents = _parse_float_list(_split_csv_field(hit.get("top10_pidents")))
        evalues = _parse_float_list(_split_csv_field(hit.get("top10_evalues")))
        best_target = hit.get("best_hit_target", "")
        best_bits = hit.get("best_hit_bits", 0.0)
    else:
        targets = _split_csv_field(getattr(hit, "top10_targets", ""))
        bits = _parse_float_list(_split_csv_field(getattr(hit, "top10_bitscores", "")))
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
    return [
        (targets[i], bits[i], pidents[i], evalues[i])
        for i in range(min(max_len, len(targets)))
    ]


def _parse_top10_prefixes(hit) -> list[str]:
    if isinstance(hit, dict):
        prefixes = _split_csv_field(hit.get("top10_prefixes"))
    else:
        prefixes = _split_csv_field(getattr(hit, "top10_prefixes", ""))
    return prefixes


def _tokens_for_target(
    target: str,
    taxonomy_lookup: Optional[dict],
    min_token_length: int,
) -> list[str]:
    if not target or not taxonomy_lookup:
        return []
    org_id = resolve_org_id(target, taxonomy_lookup)
    taxonomy_string = taxonomy_lookup.get(org_id, "")
    tokens = iter_taxonomy_tokens(taxonomy_string, min_token_length)
    if tokens:
        return [t.lower() for t in tokens]
    if "__" in org_id:
        prefix = org_id.split("__", 1)[0]
        if len(prefix) >= min_token_length:
            return [prefix.lower()]
    return []

def _target_has_prefix(target: str, prefixes: Optional[set[str]]) -> bool:
    if not prefixes:
        return True
    if not target:
        return False
    org_id = target.split("|", 1)[0]
    return any(org_id.startswith(prefix) for prefix in prefixes)


def _filter_hits_by_prefix(
    top10_hits: list[tuple[str, float, float, float]],
    prefixes: Optional[set[str]],
) -> list[tuple[str, float, float, float]]:
    if not prefixes:
        return top10_hits
    return [hit for hit in top10_hits if _target_has_prefix(hit[0], prefixes)]


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
    validated_prefixes: Optional[set[str]] = None,
    supporting_prefixes: Optional[set[str]] = None,
    min_token_length: int = 3,
    host_prefixes: Optional[set[str]] = None,
    max_bits_per_token: int = 10,
    weight_mode: str = "rank",
) -> HostSignatureModel:
    """
    Build a weighted host-signature model from marker hits.

    Uses BUSCO/COG marker hits with NO NCLDV/MIRUS/VP/PLV/CRESS/PHAGE in top-10 prefixes.
    Token weights are weighted sums of taxonomy substrings (split at |).
    """
    validated = validated_prefixes or {
        # PPV__ is the unified Preplasmiviricota prefix and carries 95,947 of the
        # v1.0.6 labels; PLV__/VP__ are kept for pre-migration bundles.
        "NCLDV__", "MIRUS__", "PPV__", "PLV__", "VP__", "CRESS__", "PHAGE__",
    }
    supporting = supporting_prefixes or set()
    host_prefixes = host_prefixes or {"EUK__"}

    taxonomy_lookup = _TAXONOMY_LOOKUP
    token_weights: dict[str, float] = {}
    token_counts: dict[str, int] = {}
    token_bits: dict[str, list[float]] = {}

    for hit in marker_hits:
        hmm_target = hit.get("hmm_target") if isinstance(hit, dict) else getattr(hit, "hmm_target", "")
        if not _is_host_marker(str(hmm_target)):
            continue
        top10_prefixes = _parse_top10_prefixes(hit)
        if any(p in validated for p in top10_prefixes):
            continue
        if supporting and any(p in supporting for p in top10_prefixes):
            continue
        if host_prefixes and not any(p in host_prefixes for p in top10_prefixes):
            continue
        top10_hits = _parse_top10_hits(hit)
        top10_hits = _filter_hits_by_prefix(top10_hits, host_prefixes)
        if not top10_hits or not taxonomy_lookup:
            continue

        fingerprint = aggregate_taxonomy_substrings(
            top10_hits,
            taxonomy_lookup,
            min_token_length=min_token_length,
            weight_mode=weight_mode,
        )
        for token, weight in fingerprint.weighted_tokens.items():
            if token in DOMAIN_TOKENS:
                continue
            token_weights[token] = token_weights.get(token, 0.0) + weight
        for token, count in fingerprint.raw_tokens.items():
            if token in DOMAIN_TOKENS:
                continue
            token_counts[token] = token_counts.get(token, 0) + count

        for rank, (target, bits, _, _) in enumerate(top10_hits):
            weight = compute_hit_weight(rank, bits, weight_mode)
            if weight <= 0:
                continue
            tokens = _tokens_for_target(
                str(target),
                taxonomy_lookup,
                min_token_length,
            )
            for token in tokens:
                if token in DOMAIN_TOKENS:
                    continue
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
    )


def _get_record_attr(record, name: str, default=None):
    if isinstance(record, dict):
        return record.get(name, default)
    return getattr(record, name, default)


def _normalize_weights(weights: Mapping[str, float]) -> dict[str, float]:
    total = sum(weights.values())
    if total <= 0:
        return {}
    return {k: v / total for k, v in weights.items()}


def _weighted_jaccard(
    gene_weights: Mapping[str, float],
    host_weights: Mapping[str, float],
) -> float:
    gene_norm = _normalize_weights(gene_weights)
    host_norm = _normalize_weights(host_weights)
    if not gene_norm or not host_norm:
        return 0.0
    all_tokens = set(gene_norm.keys()) | set(host_norm.keys())
    overlap = 0.0
    union = 0.0
    for token in all_tokens:
        g = gene_norm.get(token, 0.0)
        h = host_norm.get(token, 0.0)
        overlap += min(g, h)
        union += max(g, h)
    return overlap / union if union > 0 else 0.0


def score_host_signature(
    fingerprint: TaxonomyFingerprint,
    model: HostSignatureModel,
) -> float:
    """
    Score a gene fingerprint against the host signature distribution.

    Domain-level tokens (euk, bac, etc.) are excluded from both sides
    before comparison because they are too generic to separate host from viral.
    """
    if not fingerprint or not model or not model.token_weights:
        return 0.0
    gene_weights = {k: v for k, v in fingerprint.weighted_tokens.items() if k not in DOMAIN_TOKENS}
    host_weights = {k: v for k, v in model.token_weights.items() if k not in DOMAIN_TOKENS}
    return _weighted_jaccard(gene_weights, host_weights)




def _fingerprint_from_record(
    record,
    model: HostSignatureModel,
) -> Optional[TaxonomyFingerprint]:
    if not record or not model:
        return None
    fp = _get_record_attr(record, "taxonomy_fingerprint", None)
    if isinstance(fp, TaxonomyFingerprint):
        return fp
    if isinstance(fp, dict) and "weighted_tokens" in fp:
        return TaxonomyFingerprint(
            weighted_tokens=fp.get("weighted_tokens", {}),
            raw_tokens=fp.get("raw_tokens", {}),
        )
    top10_hits = _parse_top10_hits(record)
    taxonomy_lookup = _TAXONOMY_LOOKUP
    if not taxonomy_lookup or not top10_hits:
        return None
    return aggregate_taxonomy_substrings(
        top10_hits,
        taxonomy_lookup,
        min_token_length=model.min_token_length,
        weight_mode=model.weight_mode,
    )


def score_host_signature_record(record, model: HostSignatureModel) -> float:
    """Score a gene record using its top-10 taxonomy distribution."""
    fingerprint = _fingerprint_from_record(record, model)
    if not fingerprint:
        return 0.0
    return score_host_signature(fingerprint, model)




def host_signature_density_evalue_weighted(
    gene_records: Iterable[Mapping],
    model: HostSignatureModel,
    score_threshold: float = 0.5,
    host_prefixes: Optional[list[str]] = None,
) -> tuple[int, float, float]:
    """
    Compute host-like counts, mean score, and evalue-weighted mean score.
    """
    scores = []
    weighted_scores = []
    total_weight = 0.0
    host_like = 0
    host_prefix_set = set(host_prefixes or [])

    for record in gene_records:
        fingerprint = _fingerprint_from_record(record, model)
        if not fingerprint:
            scores.append(0.0)
            continue

        score = score_host_signature(fingerprint, model)
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
    model: Optional[HostSignatureModel],
    top_k: int = 20,
) -> list[tuple[str, float, int]]:
    """
    Return top-k host signature tokens with weights and counts.

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

    return [
        (token, float(weight), int(model.token_counts.get(token, 0)))
        for token, weight in sorted_items
    ]


def summarize_host_signature_bits(
    model: Optional[HostSignatureModel],
    top_k: int = 20,
    max_bits: int = 10,
) -> list[tuple[str, float, int, list[float]]]:
    """
    Return top-k tokens with weights, counts, and top weight samples.
    """
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
