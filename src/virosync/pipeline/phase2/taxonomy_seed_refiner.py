"""
Taxonomy-based seed refinement.

Refines seed boundaries using Diamond taxonomy calls from seed interior and
flanking genes collected in Phase 2b.
"""

from dataclasses import dataclass, replace
import logging
from typing import Optional

from virosync.ablation import AblationID, InterventionCounts
from virosync.pipeline.host_signatures import HostSignatureModel, score_host_signature_record
from virosync.pipeline.phase1.seed_merger import MergedSeed
from virosync.pipeline.phase2.boundary_diamond import has_identity_qualified_viral_hit

logger = logging.getLogger(__name__)


class A4HostAwareTaxonomyMLError(RuntimeError):
    """Raised when A4 would otherwise use the host-aware taxonomy ML path."""


@dataclass(frozen=True, slots=True)
class TaxonomySeedRefinement:
    """Selected and normal seed coordinates plus contract-ready A4 counts."""

    ablation_id: AblationID
    selected_seeds: tuple[MergedSeed, ...]
    counterfactual_seeds: tuple[MergedSeed, ...]
    intervention_counts: InterventionCounts


def validate_taxonomy_refinement_mode(
    *,
    ablation_id: AblationID,
    taxonomy_ml_enabled: bool,
) -> None:
    """Fail closed when A4 would enter the host-feature-dependent ML refiner."""

    if not isinstance(ablation_id, AblationID):
        raise TypeError("ablation_id must be an AblationID")
    if type(taxonomy_ml_enabled) is not bool:
        raise TypeError("taxonomy_ml_enabled must be a bool")
    if ablation_id is AblationID.A4 and taxonomy_ml_enabled:
        raise A4HostAwareTaxonomyMLError(
            "A4 cannot use taxonomy ML refinement because its model contains "
            "host-derived features; disable boundary_taxonomy_ml_enabled"
        )


def _has_strong_viral_signal(tax) -> bool:
    """True if gene has an identity-qualified viral top-10 Diamond hit."""
    return has_identity_qualified_viral_hit(
        getattr(tax, "top10_prefixes", []) or [],
        getattr(tax, "top10_pidents", []) or [],
    )


def _is_host_like_gene(
    tax,
    host_prefix: str,
    host_signature_model: Optional[HostSignatureModel],
    host_signature_threshold: float,
) -> bool:
    """Classify a gene as host-like using host prefix + optional Phase 1 model."""
    if not tax:
        return False
    if _has_strong_viral_signal(tax):
        return False

    top1_prefix = str(getattr(tax, "top1_prefix", "") or "")
    top10_prefixes = [str(p) for p in (getattr(tax, "top10_prefixes", []) or [])]
    has_host_prefix = bool(top1_prefix == host_prefix or any(p == host_prefix for p in top10_prefixes))
    if not has_host_prefix:
        return False

    if host_signature_model is None:
        return True

    host_score = score_host_signature_record(tax, host_signature_model)
    return host_score >= host_signature_threshold


def _trim_seed_host_edges(
    seed_start: int,
    seed_end: int,
    eve_porf_ids: list[str],
    taxonomy_map: dict,
    host_prefix: str,
    host_signature_model: Optional[HostSignatureModel],
    host_signature_threshold: float,
) -> tuple[int, int, bool, bool]:
    """
    Trim host-like genes from both edges of the seed interior until non-host is reached.
    """
    eve_genes = []
    for porf_id in eve_porf_ids:
        tax = taxonomy_map.get(porf_id)
        if not tax:
            continue
        eve_genes.append((int(getattr(tax, "start", 0) or 0), int(getattr(tax, "end", 0) or 0), tax))
    if not eve_genes:
        return seed_start, seed_end, False, False

    eve_genes.sort(key=lambda item: (item[0], item[1]))
    trimmed_start = int(seed_start)
    trimmed_end = int(seed_end)
    left_shrunk = False
    right_shrunk = False

    for start, end, tax in eve_genes:
        if end <= trimmed_start or start >= trimmed_end:
            continue
        if _is_host_like_gene(tax, host_prefix, host_signature_model, host_signature_threshold):
            next_start = max(trimmed_start, end)
            if next_start > trimmed_start:
                trimmed_start = next_start
                left_shrunk = True
        else:
            break

    for start, end, tax in reversed(eve_genes):
        if end <= trimmed_start or start >= trimmed_end:
            continue
        if _is_host_like_gene(tax, host_prefix, host_signature_model, host_signature_threshold):
            next_end = min(trimmed_end, start)
            if next_end < trimmed_end:
                trimmed_end = next_end
                right_shrunk = True
        else:
            break

    if trimmed_start >= trimmed_end:
        return seed_start, seed_end, False, False

    return trimmed_start, trimmed_end, left_shrunk, right_shrunk


def _refine_seeds_by_taxonomy_mode(
    merged_seeds: list[MergedSeed],
    taxonomy_map: dict,  # porf_id -> GeneTaxonomy
    seed_gene_mappings: dict,  # seed_id -> SeedGeneMapping
    host_prefix: str = "EUK__",
    expansion_kb: int = 5,
    host_signature_model: Optional[HostSignatureModel] = None,
    host_signature_threshold: float = 0.5,
    *,
    host_coordinate_paths_enabled: bool,
    emit_log: bool,
) -> list[MergedSeed]:
    """
    Refine seed boundaries around viral-positive genes while stopping at host-like edges.

    Strategy:
    - trim host-like genes from seed interior edges until a non-host gene is reached
    - expand only to viral-positive genes in flanks
    - stop expansion on each side once a host-like flanking gene is encountered
    - clamp final bounds to Phase 2b flanking-gene envelope
    """
    refined = []
    expansion_bp = max(0, expansion_kb * 1000)
    changed = 0
    extended = 0
    contracted = 0

    for seed in merged_seeds:
        mapping = seed_gene_mappings.get(seed.seed_id)
        if not mapping:
            refined.append(seed)
            continue

        if host_coordinate_paths_enabled:
            trimmed_start, trimmed_end, left_shrunk, right_shrunk = _trim_seed_host_edges(
                seed_start=seed.start,
                seed_end=seed.end,
                eve_porf_ids=list(mapping.eve_porf_ids),
                taxonomy_map=taxonomy_map,
                host_prefix=host_prefix,
                host_signature_model=host_signature_model,
                host_signature_threshold=host_signature_threshold,
            )
        else:
            trimmed_start, trimmed_end = seed.start, seed.end
            left_shrunk = right_shrunk = False

        left_anchor = trimmed_start
        right_anchor = trimmed_end
        has_viral_anchor = False

        for porf_id in mapping.eve_porf_ids:
            tax = taxonomy_map.get(porf_id)
            if not tax or not _has_strong_viral_signal(tax):
                continue
            start = int(getattr(tax, "start", 0) or 0)
            end = int(getattr(tax, "end", 0) or 0)
            if end <= trimmed_start or start >= trimmed_end:
                continue
            if _is_host_like_gene(tax, host_prefix, host_signature_model, host_signature_threshold):
                continue
            left_anchor = min(left_anchor, start)
            right_anchor = max(right_anchor, end)
            has_viral_anchor = True

        upstream_host_barrier_end = None
        for porf_id in mapping.upstream_porf_ids:
            tax = taxonomy_map.get(porf_id)
            if not tax:
                continue
            if _is_host_like_gene(tax, host_prefix, host_signature_model, host_signature_threshold):
                if host_coordinate_paths_enabled:
                    upstream_host_barrier_end = int(getattr(tax, "end", 0) or 0)
                    break
                continue
            if _has_strong_viral_signal(tax):
                left_anchor = min(left_anchor, int(getattr(tax, "start", 0) or 0))
                right_anchor = max(right_anchor, int(getattr(tax, "end", 0) or 0))
                has_viral_anchor = True

        downstream_host_barrier_start = None
        for porf_id in mapping.downstream_porf_ids:
            tax = taxonomy_map.get(porf_id)
            if not tax:
                continue
            if _is_host_like_gene(tax, host_prefix, host_signature_model, host_signature_threshold):
                if host_coordinate_paths_enabled:
                    downstream_host_barrier_start = int(getattr(tax, "start", 0) or 0)
                    break
                continue
            if _has_strong_viral_signal(tax):
                left_anchor = min(left_anchor, int(getattr(tax, "start", 0) or 0))
                right_anchor = max(right_anchor, int(getattr(tax, "end", 0) or 0))
                has_viral_anchor = True

        if has_viral_anchor:
            raw_start = max(0, left_anchor - expansion_bp)
            raw_end = right_anchor + expansion_bp
        else:
            raw_start = trimmed_start
            raw_end = trimmed_end

        # Keep refinement inside known flanking-gene range from Phase 2b.
        new_start = max(int(mapping.flank_start_bp), raw_start)
        new_end = min(int(mapping.flank_end_bp), raw_end)

        min_start_allowed = trimmed_start if left_shrunk else int(mapping.flank_start_bp)
        if upstream_host_barrier_end is not None:
            min_start_allowed = max(min_start_allowed, upstream_host_barrier_end)
        new_start = max(new_start, min_start_allowed)

        max_end_allowed = trimmed_end if right_shrunk else int(mapping.flank_end_bp)
        if downstream_host_barrier_start is not None:
            max_end_allowed = min(max_end_allowed, downstream_host_barrier_start)
        new_end = min(new_end, max_end_allowed)

        if new_end <= new_start:
            # Keep the host-trimmed interior span if expansion constraints became invalid.
            new_start, new_end = trimmed_start, trimmed_end
        if new_end <= new_start:
            refined.append(seed)
            continue

        if new_start != seed.start or new_end != seed.end:
            changed += 1
            if new_start < seed.start or new_end > seed.end:
                extended += 1
            if new_start > seed.start or new_end < seed.end:
                contracted += 1
        refined.append(replace(seed, start=new_start, end=new_end))

    if changed and emit_log:
        logger.info(
            "Taxonomy seed refinement updated %d/%d seeds (extended=%d, contracted=%d, expansion_kb=%d)",
            changed,
            len(merged_seeds),
            extended,
            contracted,
            expansion_kb,
        )

    return refined


def evaluate_taxonomy_seed_refinement(
    merged_seeds: list[MergedSeed],
    taxonomy_map: dict,
    seed_gene_mappings: dict,
    host_prefix: str = "EUK__",
    expansion_kb: int = 5,
    host_signature_model: Optional[HostSignatureModel] = None,
    host_signature_threshold: float = 0.5,
    *,
    ablation_id: AblationID = AblationID.A0,
    taxonomy_ml_enabled: bool = False,
) -> TaxonomySeedRefinement:
    """Return normal counterfactual and selected heuristic-refinement coordinates.

    A4 disables only host-derived edge contraction and flank barriers. Viral
    positive expansion and the normal flanking-gene envelope remain active.
    Counter units are input seeds; a seed counts as intervened/changed when its
    selected A4 coordinates differ from the normal host-aware coordinates.
    """

    validate_taxonomy_refinement_mode(
        ablation_id=ablation_id,
        taxonomy_ml_enabled=taxonomy_ml_enabled,
    )
    counterfactual = _refine_seeds_by_taxonomy_mode(
        merged_seeds=merged_seeds,
        taxonomy_map=taxonomy_map,
        seed_gene_mappings=seed_gene_mappings,
        host_prefix=host_prefix,
        expansion_kb=expansion_kb,
        host_signature_model=host_signature_model,
        host_signature_threshold=host_signature_threshold,
        host_coordinate_paths_enabled=True,
        emit_log=ablation_id is not AblationID.A4,
    )

    if ablation_id is AblationID.A4:
        selected = _refine_seeds_by_taxonomy_mode(
            merged_seeds=merged_seeds,
            taxonomy_map=taxonomy_map,
            seed_gene_mappings=seed_gene_mappings,
            host_prefix=host_prefix,
            expansion_kb=expansion_kb,
            host_signature_model=host_signature_model,
            host_signature_threshold=host_signature_threshold,
            host_coordinate_paths_enabled=False,
            emit_log=True,
        )
        changed = sum(
            normal.start != ablated.start or normal.end != ablated.end
            for normal, ablated in zip(counterfactual, selected, strict=True)
        )
        counts = InterventionCounts(
            opportunities=len(merged_seeds),
            interventions=changed,
            changed=changed,
        )
    else:
        selected = counterfactual
        counts = InterventionCounts()

    return TaxonomySeedRefinement(
        ablation_id=ablation_id,
        selected_seeds=tuple(selected),
        counterfactual_seeds=tuple(counterfactual),
        intervention_counts=counts,
    )


def refine_seeds_by_taxonomy(
    merged_seeds: list[MergedSeed],
    taxonomy_map: dict,
    seed_gene_mappings: dict,
    host_prefix: str = "EUK__",
    expansion_kb: int = 5,
    host_signature_model: Optional[HostSignatureModel] = None,
    host_signature_threshold: float = 0.5,
    *,
    ablation_id: AblationID = AblationID.A0,
    taxonomy_ml_enabled: bool = False,
) -> list[MergedSeed]:
    """Return the selected arm's seed coordinates with production-compatible API."""

    result = evaluate_taxonomy_seed_refinement(
        merged_seeds=merged_seeds,
        taxonomy_map=taxonomy_map,
        seed_gene_mappings=seed_gene_mappings,
        host_prefix=host_prefix,
        expansion_kb=expansion_kb,
        host_signature_model=host_signature_model,
        host_signature_threshold=host_signature_threshold,
        ablation_id=ablation_id,
        taxonomy_ml_enabled=taxonomy_ml_enabled,
    )
    return list(result.selected_seeds)
