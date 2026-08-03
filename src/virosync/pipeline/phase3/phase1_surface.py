"""Pure A1 adapter from Phase-1 seeds to canonical output records."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from virosync.ablation import AblationID, InterventionCounts
from virosync.pipeline.phase1.hhg_seeding import Anchor
from virosync.pipeline.phase1.seed_merger import MergedSeed
from virosync.output_contract import (
    canonical_family,
    normalize_effective_eve_class,
)
from virosync.pipeline.phase3.evidence_synthesizer import (
    VerificationResult,
    VerificationStatus,
    infer_ppv_subtype,
    summarize_marker_hits,
)
from virosync.pipeline.phase3.mcp_detection import is_mcp_gene


@dataclass(frozen=True, slots=True)
class Phase1SeedSurface:
    """Identity-preserving detailed and canonical surfaces for A1 output."""

    results: tuple[VerificationResult, ...]
    intervention_counts: InterventionCounts

    @property
    def detailed_results(self) -> tuple[VerificationResult, ...]:
        """Return every exported seed for the all-candidate output surface."""

        return self.results

    @property
    def canonical_results(self) -> tuple[VerificationResult, ...]:
        """Return every exported seed as an explicitly preselected prediction."""

        return self.results


def _base_gene_id(anchor: Anchor) -> str:
    porf_id = anchor.porf_id.split("|aa", 1)[0]
    if porf_id:
        return porf_id
    return (
        f"{anchor.scaffold}:{anchor.start}-{anchor.end}:"
        f"{anchor.strand}:{anchor.hallmark_gene}"
    )


def _eligible_anchors(seed: MergedSeed) -> tuple[tuple[str, Anchor], ...]:
    """Return one deterministic best marker hit per seed-overlapping protein."""

    grouped: dict[str, list[Anchor]] = {}
    for anchor in (*seed.anchors, *seed.hhg_anchors):
        if not anchor.hallmark_gene:
            continue
        if anchor.scaffold != seed.scaffold:
            continue
        if anchor.end <= seed.start or anchor.start >= seed.end:
            continue
        grouped.setdefault(_base_gene_id(anchor), []).append(anchor)

    selected: list[tuple[str, Anchor]] = []
    for gene_id in sorted(grouped):
        anchor = min(
            grouped[gene_id],
            key=lambda item: (
                -float(item.score),
                float(item.evalue),
                item.hallmark_gene,
                item.start,
                item.end,
                item.strand,
            ),
        )
        selected.append((gene_id, anchor))
    return tuple(selected)


def _seed_to_result(seed: MergedSeed) -> VerificationResult:
    eligible_anchors = _eligible_anchors(seed)
    hallmark_genes = [anchor.hallmark_gene for _, anchor in eligible_anchors]
    marker_summary = summarize_marker_hits(hallmark_genes)
    mcp_gene_ids = [
        gene_id
        for gene_id, anchor in eligible_anchors
        if is_mcp_gene(anchor.hallmark_gene)
    ]
    region_classification = canonical_family(seed.predicted_family)
    # Phase 1 has no CRESS marker family. Identity-qualified gene taxonomy may
    # assign CRESS later during evidence synthesis.
    likely_family = (
        region_classification
        if region_classification in {"NCLDV", "MIRUS", "PPV", "MIXED"}
        else "UNKNOWN"
    )

    return VerificationResult(
        eve_id=f"EVE_{seed.scaffold}_{seed.start}-{seed.end}",
        scaffold=seed.scaffold,
        start=seed.start,
        end=seed.end,
        status=VerificationStatus.AMBIGUOUS,
        final_confidence=0.0,
        confidence_tier="LOW",
        ablation_id=AblationID.A1,
        score_components={
            "prediction_stage": "phase1_seed_surface",
            "confidence_kind": "not_scored",
        },
        hallmark_count=len(hallmark_genes),
        hallmark_diversity=len(set(hallmark_genes)),
        has_virus_specific_marker=bool(hallmark_genes),
        hallmark_genes=hallmark_genes,
        marker_category_hits=marker_summary["categories"],
        marker_family_hits=marker_summary["families"],
        marker_complement_score=marker_summary["complement_score"],
        marker_dominant_family=marker_summary["dominant_family"],
        marker_dominant_fraction=marker_summary["dominant_fraction"],
        has_mcp=bool(mcp_gene_ids),
        mcp_gene_ids=mcp_gene_ids,
        region_classification=region_classification,
        region_classification_ncldv_markers=seed.region_classification_ncldv_markers,
        region_classification_vp_plv_markers=seed.region_classification_vp_plv_markers,
        region_classification_mirus_markers=seed.region_classification_mirus_markers,
        ppv_subtype=infer_ppv_subtype(hallmark_genes),
        likely_family=likely_family,
        # A1 stops before the gene taxonomy, so there is no weighted vote to
        # publish. The seed marker family is the only taxonomy this arm has, and
        # without this the published class would be UNKNOWN for every A1 seed.
        taxonomy_class=normalize_effective_eve_class(likely_family),
        seed_sources=sorted(set(seed.sources)),
        seed_confidence=seed.confidence,
        seed_hhg_score=seed.hhg_score,
        seed_novelty_score=seed.novelty_score,
        seed_compositional_score=seed.compositional_score,
        kfd=seed.max_kfd,
        gc_deviation=seed.gc_deviation,
        cub_deviation=seed.cub_deviation,
    )


def build_phase1_seed_surface(seeds: Sequence[MergedSeed]) -> Phase1SeedSurface:
    """Adapt Tier-1-vetted Phase-1 seeds without running Phase 2 or Phase 3."""

    results = tuple(_seed_to_result(seed) for seed in seeds)
    exported = len(results)
    return Phase1SeedSurface(
        results=results,
        intervention_counts=InterventionCounts(
            opportunities=len(seeds),
            interventions=exported,
            # A1 cannot know which seed outputs A0 would later reject or alter
            # without running the deliberately skipped phases. The paired
            # A0/A1 output comparison owns that changed-outcome measurement.
            changed=0,
        ),
    )
