"""Pure A1 adapter from Phase-1 seeds to canonical output records."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from virosync.ablation import AblationID, InterventionCounts
from virosync.output_contract import (
    canonical_family,
    normalize_effective_eve_class,
)
from virosync.pipeline.host_signatures import get_taxonomy_lookup
from virosync.pipeline.phase1.hhg_seeding import Anchor
from virosync.pipeline.phase1.seed_merger import MergedSeed
from virosync.pipeline.phase2.boundary_refiner import (
    RefinedBoundary,
    assign_boundary_candidate_ids,
    boundary_candidate_id,
)
from virosync.pipeline.phase2.seed_extension import confirmed_seed_core
from virosync.pipeline.phase3.evidence_synthesizer import (
    VerificationResult,
    VerificationStatus,
    consensus_taxonomy_class,
    deduplicate_preferring_mcp,
    has_validated_capsid_hit,
    infer_likely_family,
    infer_ppv_subtype,
    is_validated_marker,
    summarize_marker_hits,
    validated_capsid_gene_ids,
    with_profile_class_families,
)
from virosync.pipeline.phase3.gene_taxonomy import has_only_peve_viral_hits, peve_profile_class
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
    return f"{anchor.scaffold}:{anchor.start}-{anchor.end}:{anchor.strand}:{anchor.hallmark_gene}"


def _overlapping_anchors(seed: MergedSeed) -> tuple[Anchor, ...]:
    """Return every marker anchor that sits on the seed, before deduplication."""
    return tuple(
        anchor
        for anchor in (*seed.anchors, *seed.hhg_anchors)
        if anchor.hallmark_gene
        and anchor.scaffold == seed.scaffold
        and anchor.end > seed.start
        and anchor.start < seed.end
    )


def _eligible_anchors(anchors: Sequence[Anchor]) -> tuple[tuple[str, Anchor], ...]:
    """Return one deterministic best marker hit per protein."""
    grouped: dict[str, list[Anchor]] = {}
    for anchor in anchors:
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


def _top10_records(anchors: Sequence[Anchor]) -> list[dict]:
    """The top-10 taxonomy fields of anchors, as Phase 3 keeps them for its hallmark hits."""
    return [{"top10_prefixes": anchor.top10_prefixes, "top10_pidents": anchor.top10_pidents} for anchor in anchors]


def _profile_vote_class(validated_anchors: Sequence[Anchor]) -> str | None:
    """Return the Phase 3 vote class when the EUK-pEVE rules decide it, else ``None``.

    A1 stops before the gene taxonomy, so the seed marker family is the only
    taxonomy it has, and without it the published class would be UNKNOWN for every
    seed. The EUK-pEVE rules of Phase 3 apply on top of it, marker by marker. They
    decide the class when a structure-derived profile has only EUK-pEVE references
    (it then votes its class), or when every qualified viral hit of the validated
    markers is EUK-pEVE. The class is then the Phase 3 vote
    (:func:`consensus_taxonomy_class`, with the same taxonomy lookup) over the
    validated markers, whatever the other markers show. Otherwise the class and the
    likely family come from the seed family.
    """
    peve_only = has_only_peve_viral_hits(_top10_records(validated_anchors))
    has_profile_class = any(peve_profile_class(anchor.hallmark_gene, anchor) for anchor in validated_anchors)
    if not (peve_only or has_profile_class):
        return None
    return consensus_taxonomy_class(validated_anchors, taxonomy_lookup=get_taxonomy_lookup())


def _seed_to_result(seed: MergedSeed, *, eve_id: str) -> VerificationResult:
    overlapping_anchors = _overlapping_anchors(seed)
    eligible_anchors = _eligible_anchors(overlapping_anchors)
    hallmark_genes = [anchor.hallmark_gene for _, anchor in eligible_anchors]
    validated_anchors = [anchor for anchor in overlapping_anchors if is_validated_marker(anchor)]
    # The class family comes from the profile that wins each protein in the Phase 3
    # vote, so family and class never disagree on a tie.
    marker_summary = with_profile_class_families(
        summarize_marker_hits(hallmark_genes),
        deduplicate_preferring_mcp(validated_anchors),
    )
    mcp_gene_ids = [gene_id for gene_id, anchor in eligible_anchors if is_mcp_gene(anchor.hallmark_gene)]
    region_classification = canonical_family(seed.predicted_family)
    # Phase 1 has no CRESS marker family. Identity-qualified gene taxonomy may
    # assign CRESS later during evidence synthesis.
    likely_family = region_classification if region_classification in {"NCLDV", "MIRUS", "PPV", "MIXED"} else "UNKNOWN"
    vote_class = _profile_vote_class(validated_anchors)

    result = VerificationResult(
        eve_id=eve_id,
        scaffold=seed.scaffold,
        start=seed.start,
        end=seed.end,
        confirmed_seeds=seed.confirmed_seeds,
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
        # Every overlapping anchor decides, not only the best profile per protein.
        has_capsid_evidence=bool(mcp_gene_ids) or has_validated_capsid_hit(overlapping_anchors),
        mcp_gene_ids=validated_capsid_gene_ids(
            anchor for anchor in overlapping_anchors if is_mcp_gene(anchor.hallmark_gene)
        ),
        capsid_gene_ids=validated_capsid_gene_ids(overlapping_anchors),
        region_classification=region_classification,
        region_classification_ncldv_markers=seed.region_classification_ncldv_markers,
        region_classification_vp_plv_markers=seed.region_classification_vp_plv_markers,
        region_classification_mirus_markers=seed.region_classification_mirus_markers,
        ppv_subtype=infer_ppv_subtype(hallmark_genes),
        likely_family=likely_family,
        # The seed marker family is the published class, unless the EUK-pEVE rules below decide.
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
    if vote_class is not None:
        # The likely family follows the vote by the same rule as Phase 3, so the two
        # agree on the marker dominance of a mixed region and on EUK-pEVE-only ones.
        result.taxonomy_class = vote_class
        result.hallmark_taxonomy_records = _top10_records(validated_anchors)
        result.likely_family = infer_likely_family(result)
    return result


def build_phase1_seed_surface(seeds: Sequence[MergedSeed]) -> Phase1SeedSurface:
    """Adapt Tier-1-vetted Phase-1 seeds without running Phase 2 or Phase 3."""
    seeds = [confirmed_seed_core(seed) for seed in seeds]
    boundaries = assign_boundary_candidate_ids(
        [
            RefinedBoundary(
                scaffold=seed.scaffold,
                start=seed.start,
                end=seed.end,
                seed_id=seed.seed_id,
                confirmed_seeds=seed.confirmed_seeds,
            )
            for seed in seeds
        ]
    )
    results = tuple(
        _seed_to_result(seed, eve_id=boundary_candidate_id(boundary))
        for seed, boundary in zip(seeds, boundaries, strict=True)
    )
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
