"""
Phase 2a: Trim boundaries using host-signature density.

Uses gene taxonomy records + weighted host-signature model to identify
windows with strong host-like signal and trims regions inward.

Supports two modes of operation:
1. Traditional: Receive gene_taxonomy_map with pre-formatted records
2. Pre-computed: Receive taxonomy directly from boundary_diamond.py GeneTaxonomy objects

When pre-computed taxonomy is provided, it's used instead of requiring
the gene_taxonomy_map format, avoiding duplicate Diamond runs.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

from virosync.ablation import AblationID, InterventionCounts
from virosync.pipeline.host_signatures import (
    HostSignatureModel,
    score_host_signature_record,
)
from virosync.pipeline.phase1.seed_merger import MergedSeed

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy


@dataclass
class HostTrimParams:
    window_bp: int = 5000
    step_bp: int = 1000
    max_host_fraction: float = 0.3
    min_viral_fraction: float = 0.05
    host_score_threshold: float = 0.3
    buffer_kb: int = 5
    use_control_baseline: bool = True


def _convert_gene_taxonomy_to_record(tax: "GeneTaxonomy") -> dict:
    """
    Convert a GeneTaxonomy object to record dict format.

    This bridges the GeneTaxonomy dataclass from boundary_diamond.py
    to the dict format expected by the trimming functions.

    Args:
        tax: GeneTaxonomy object from boundary_diamond.py

    Returns:
        Dict with fields expected by _get_attr and window analysis
    """
    return {
        "porf_id": tax.porf_id,
        "porf_start": tax.start,
        "porf_end": tax.end,
        "top1_target": tax.top1_target,
        "top1_prefix": tax.top1_prefix,
        "top1_pident": tax.top1_pident,
        "top1_evalue": getattr(tax, "top1_evalue", 1.0),
        "top10_targets": getattr(tax, "top10_targets", []),
        "top10_bitscores": getattr(tax, "top10_bits", []),
        "top10_pidents": getattr(tax, "top10_pidents", []),
        "top10_evalues": getattr(tax, "top10_evalues", []),
        "taxonomy_fingerprint": getattr(tax, "taxonomy_fingerprint", None),
        "has_ncldv_mirus": tax.has_ncldv_mirus,
        "has_vp_plv": tax.has_vp_plv,
        "has_viral": tax.has_viral,
        "has_hit": tax.has_hit,
    }


def _get_attr(record, name, default=None):
    if isinstance(record, dict):
        return record.get(name, default)
    return getattr(record, name, default)


def _iter_windows(start: int, end: int, window_bp: int, step_bp: int):
    if window_bp <= 0 or step_bp <= 0:
        return
    pos = start
    while pos < end:
        win_end = min(end, pos + window_bp)
        yield pos, win_end
        if pos + step_bp == pos:
            break
        pos += step_bp


def _window_gene_records(records: list, win_start: int, win_end: int) -> list:
    window_records = []
    for record in records:
        r_start = _get_attr(record, "porf_start", 0)
        r_end = _get_attr(record, "porf_end", 0)
        if r_start < win_end and r_end > win_start:
            window_records.append(record)
    return window_records


def _trim_seed_by_host_signature_normal(
    seed: MergedSeed,
    gene_records: list,
    host_model: HostSignatureModel,
    validated_markers: Optional[list] = None,
    params: Optional[HostTrimParams] = None,
    precomputed_taxonomy: Optional[dict[str, "GeneTaxonomy"]] = None,
) -> tuple[MergedSeed, dict]:
    """
    Trim a single seed based on host-like window density.

    Args:
        seed: MergedSeed to trim
        gene_records: List of gene taxonomy records (legacy format)
        host_model: Host signature model for scoring
        validated_markers: Optional list of validated markers
        params: Trimming parameters
        precomputed_taxonomy: Optional dict mapping pORF ID to GeneTaxonomy
            from boundary_diamond.py. When provided, gene_records from this
            taxonomy are used instead of the gene_records parameter.

    Returns:
        Tuple of (trimmed_seed, summary_dict)
    """
    params = params or HostTrimParams()
    buffer_bp = max(0, params.buffer_kb * 1000)

    # If precomputed_taxonomy provided, extract gene records for this seed's region
    if precomputed_taxonomy is not None:
        gene_records = [
            _convert_gene_taxonomy_to_record(tax)
            for tax in precomputed_taxonomy.values()
            if (
                tax.scaffold == seed.scaffold
                and tax.start < seed.end
                and tax.end > seed.start
            )
        ]

    if not gene_records:
        return seed, {
            "reason": "no_gene_taxonomy",
            "trimmed_start": seed.start,
            "trimmed_end": seed.end,
            "host_consensus_taxonomy": ".",
        }

    record_scores: dict[str, float] = {}
    for record in gene_records:
        porf_id = _get_attr(record, "porf_id", "")
        record_scores[porf_id] = score_host_signature_record(record, host_model)

    marker_positions = []
    if validated_markers:
        for marker in validated_markers:
            if getattr(marker, "validation_status", "") not in ("validated", "validated_novel"):
                continue
            if marker.scaffold != seed.scaffold:
                continue
            marker_positions.append((marker.start + marker.end) // 2)

    windows = []
    for win_start, win_end in _iter_windows(seed.start, seed.end, params.window_bp, params.step_bp):
        win_records = _window_gene_records(gene_records, win_start, win_end)
        if not win_records:
            continue
        host_like = 0
        viral_like = 0
        scores = []
        for record in win_records:
            porf_id = _get_attr(record, "porf_id", "")
            score = record_scores.get(porf_id, 0.0)
            scores.append(score)
            if score >= params.host_score_threshold:
                host_like += 1
            has_viral = bool(
                _get_attr(record, "has_ncldv_mirus", False)
                or _get_attr(record, "has_vp_plv", False)
                or _get_attr(record, "has_viral", False)
            )
            if has_viral:
                viral_like += 1

        gene_count = len(win_records)
        host_fraction = host_like / gene_count if gene_count else 0.0
        viral_fraction = viral_like / gene_count if gene_count else 0.0
        mean_score = sum(scores) / gene_count if gene_count else 0.0
        has_marker = any(win_start <= m < win_end for m in marker_positions)
        good = host_fraction <= params.max_host_fraction and viral_fraction >= params.min_viral_fraction
        windows.append(
            {
                "start": win_start,
                "end": win_end,
                "host_fraction": host_fraction,
                "viral_fraction": viral_fraction,
                "mean_score": mean_score,
                "has_marker": has_marker,
                "good": good,
            }
        )

    if not windows:
        return seed, {
            "reason": "no_windows",
            "trimmed_start": seed.start,
            "trimmed_end": seed.end,
            "host_consensus_taxonomy": ".",
        }

    good_marker = [w for w in windows if w["good"] and w["has_marker"]]
    good_any = [w for w in windows if w["good"]]
    reason = ""

    if good_marker:
        starts = [w["start"] for w in good_marker]
        ends = [w["end"] for w in good_marker]
        reason = "marker_windows"
    elif good_any:
        midpoint = (seed.start + seed.end) // 2
        mid_windows = [w for w in good_any if w["start"] <= midpoint < w["end"]]
        if mid_windows:
            starts = [w["start"] for w in mid_windows]
            ends = [w["end"] for w in mid_windows]
            reason = "midpoint_windows"
        else:
            # fallback: select window with lowest host_fraction, then highest viral_fraction
            best = sorted(
                good_any,
                key=lambda w: (w["host_fraction"], -w["viral_fraction"]),
            )[0]
            starts = [best["start"]]
            ends = [best["end"]]
            reason = "best_window"
    else:
        return seed, {
            "reason": "no_good_windows",
            "trimmed_start": seed.start,
            "trimmed_end": seed.end,
            "host_consensus_taxonomy": ".",
        }

    trimmed_start = max(seed.start, min(starts) - buffer_bp)
    trimmed_end = min(seed.end, max(ends) + buffer_bp)
    if trimmed_end <= trimmed_start:
        trimmed_start, trimmed_end = seed.start, seed.end
        reason = "invalid_trim"

    trimmed_seed = MergedSeed(**{field: getattr(seed, field) for field in seed.__dataclass_fields__})
    trimmed_seed.start = trimmed_start
    trimmed_seed.end = trimmed_end

    summary = {
        "reason": reason,
        "trimmed_start": trimmed_start,
        "trimmed_end": trimmed_end,
        "window_count": len(windows),
        "good_windows": len(good_any),
        "good_marker_windows": len(good_marker),
        "host_consensus_taxonomy": ".",
    }
    return trimmed_seed, summary


def trim_seed_by_host_signature(
    seed: MergedSeed,
    gene_records: list,
    host_model: HostSignatureModel,
    validated_markers: Optional[list] = None,
    params: Optional[HostTrimParams] = None,
    precomputed_taxonomy: Optional[dict[str, "GeneTaxonomy"]] = None,
    ablation_id: AblationID = AblationID.A0,
) -> tuple[MergedSeed, dict]:
    """Evaluate normal host trimming and select the configured arm's coordinates.

    A4 retains the input seed coordinates while recording the normal trim as a
    counterfactual. Counter units are seeds: every A4 call is one opportunity,
    and a normal coordinate change is both an intervention and a changed output.
    Other arms execute normal production behavior and expose zero A4 counters.
    """

    if not isinstance(ablation_id, AblationID):
        raise TypeError("ablation_id must be an AblationID")

    counterfactual_seed, counterfactual_summary = _trim_seed_by_host_signature_normal(
        seed=seed,
        gene_records=gene_records,
        host_model=host_model,
        validated_markers=validated_markers,
        params=params,
        precomputed_taxonomy=precomputed_taxonomy,
    )
    counterfactual_changed = (
        counterfactual_seed.start != seed.start
        or counterfactual_seed.end != seed.end
    )
    bypassed = ablation_id is AblationID.A4 and counterfactual_changed
    counts = (
        InterventionCounts(
            opportunities=1,
            interventions=int(counterfactual_changed),
            changed=int(counterfactual_changed),
        )
        if ablation_id is AblationID.A4
        else InterventionCounts()
    )
    selected_seed = seed if ablation_id is AblationID.A4 else counterfactual_seed

    summary = dict(counterfactual_summary)
    summary.update(
        {
            "ablation_id": ablation_id.value,
            "counterfactual_trimmed_start": counterfactual_seed.start,
            "counterfactual_trimmed_end": counterfactual_seed.end,
            "counterfactual_reason": counterfactual_summary.get("reason", ""),
            "trimmed_start": selected_seed.start,
            "trimmed_end": selected_seed.end,
            "host_coordinate_change_opportunities": counts.opportunities,
            "host_coordinate_change_interventions": counts.interventions,
            "host_coordinate_change_changed": counts.changed,
        }
    )
    if bypassed:
        summary["reason"] = "a4_host_coordinate_change_bypass"
    return selected_seed, summary


def trim_seeds_by_host_signature(
    seeds: list[MergedSeed],
    gene_taxonomy_map: Optional[dict[str, tuple[list, dict]]] = None,
    host_model: Optional[HostSignatureModel] = None,
    validated_markers: Optional[list] = None,
    params: Optional[HostTrimParams] = None,
    precomputed_taxonomy: Optional[dict[str, "GeneTaxonomy"]] = None,
    ablation_id: AblationID = AblationID.A0,
) -> tuple[list[MergedSeed], list[dict]]:
    """
    Trim multiple seeds based on host-like window density.

    Supports two modes:
    1. Traditional: Use gene_taxonomy_map with pre-formatted records per EVE
    2. Pre-computed: Use precomputed_taxonomy dict from boundary_diamond.py

    When precomputed_taxonomy is provided, gene_taxonomy_map is ignored and
    the taxonomy is filtered per-seed based on coordinate overlap.

    Args:
        seeds: List of MergedSeed objects to trim
        gene_taxonomy_map: Legacy format mapping EVE ID to (records, summary)
        host_model: Host signature model for scoring
        validated_markers: Optional list of validated markers
        params: Trimming parameters
        precomputed_taxonomy: Optional dict mapping pORF ID to GeneTaxonomy
            from boundary_diamond.py. When provided, this is used instead
            of gene_taxonomy_map.

    Returns:
        Tuple of (list of trimmed seeds, list of summary dicts)
    """
    if not isinstance(ablation_id, AblationID):
        raise TypeError("ablation_id must be an AblationID")

    trimmed: list[MergedSeed] = []
    summaries: list[dict] = []
    for seed in seeds:
        eve_id = f"EVE_{seed.scaffold}_{seed.start}-{seed.end}"
        gene_records: list = []

        # Use precomputed taxonomy if provided, otherwise fall back to legacy map
        if precomputed_taxonomy is None and gene_taxonomy_map:
            if eve_id in gene_taxonomy_map:
                gene_records = gene_taxonomy_map[eve_id][0] or []

        trimmed_seed, summary = trim_seed_by_host_signature(
            seed=seed,
            gene_records=gene_records,
            host_model=host_model,
            validated_markers=validated_markers,
            params=params,
            precomputed_taxonomy=precomputed_taxonomy,
            ablation_id=ablation_id,
        )
        summary.update(
            {
                "scaffold": seed.scaffold,
                "start": seed.start,
                "end": seed.end,
            }
        )
        trimmed.append(trimmed_seed)
        summaries.append(summary)
    return trimmed, summaries
