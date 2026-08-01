from __future__ import annotations

from virosync.pipeline.phase0.prodigal import GenePrediction
from virosync.pipeline.phase1.region_assembly import (
    ValidatedMarkerHit,
    assemble_compact_cress_regions,
)
from virosync.pipeline.phase1.seed_merger import MergedSeed
from virosync.pipeline.phase1.viral_markers import (
    is_cress_specific_top1_marker,
    is_identity_qualified_cress_marker,
)
from virosync.pipeline.phase2.boundary_diamond import pORF
from virosync.pipeline.phase2.boundary_refiner import (
    RefinedBoundary,
    extend_seeds_by_genes,
    merge_adjacent_viral_boundaries,
)
from virosync.pipeline.phase2.host_signature_trim import (
    trim_seed_by_host_signature,
)
from virosync.pipeline.host_signatures import HostSignatureModel


def _marker(
    gene: int,
    start: int,
    end: int,
    *,
    model: str = "VS000803",
    score: float = 100.0,
    prefixes: str = "CRESS__,EUK__",
    pidents: str = "45.0,70.0",
    status: str = "validated",
    targets: str = "CRESS__reference,EUK__host",
) -> ValidatedMarkerHit:
    return ValidatedMarkerHit(
        query_porf=f"ctg_{gene}|aa1-100",
        scaffold="ctg",
        start=start,
        end=end,
        strand="+",
        hmm_target=model,
        hmm_score=score,
        hmm_evalue=1e-30,
        validation_status=status,
        top10_prefixes=prefixes,
        best_hit_target="CRESS__reference",
        best_hit_pident=45.0,
        best_hit_bits=200.0,
        has_ncldv=0,
        has_mirus=0,
        has_plv=0,
        has_vp=0,
        has_viral=1,
        top10_targets=targets,
        top10_pidents=pidents,
    )


def _gene(gene: int, start: int, end: int) -> GenePrediction:
    return GenePrediction(
        gene_id=f"ctg_{gene}",
        scaffold="ctg",
        start=start,
        end=end,
        strand="+",
        protein="M",
    )


def test_cress_marker_requires_paired_identity_support() -> None:
    assert is_identity_qualified_cress_marker(_marker(1, 100, 400))
    assert not is_identity_qualified_cress_marker(
        _marker(1, 100, 400, pidents="24.96,70.0")
    )
    assert not is_identity_qualified_cress_marker(
        _marker(
            1,
            100,
            400,
            prefixes="EUK__,CRESS__",
            pidents="70.0,24.9",
        )
    )
    assert not is_identity_qualified_cress_marker(
        _marker(1, 100, 400, model="VS000791")
    )


def test_single_gene_specificity_uses_the_top_hit_target() -> None:
    assert is_cress_specific_top1_marker(_marker(1, 100, 400))
    assert is_cress_specific_top1_marker(
        _marker(
            1,
            100,
            400,
            prefixes="PHAGE__,CRESS__",
            pidents="52.9,43.3",
            targets="PHAGE__MONDNA__reference,CRESS__reference",
        )
    )
    assert not is_cress_specific_top1_marker(
        _marker(
            1,
            100,
            400,
            prefixes="EUK__,CRESS__",
            pidents="98.4,27.7",
            targets="EUK__host,CRESS__reference",
        )
    )


def test_single_cress_gene_keeps_exact_gene_bounds() -> None:
    marker = _marker(1, 100, 400)

    regions = assemble_compact_cress_regions(
        [marker],
        {"ctg": [_gene(1, 100, 400)]},
    )

    assert len(regions) == 1
    assert (regions[0].start, regions[0].end) == (100, 400)
    assert regions[0].predicted_family == "CRESS"


def test_nearby_rep_and_capsid_are_grouped_without_padding() -> None:
    rep = _marker(1, 100, 400, model="VS000803", score=200.0)
    duplicate_rep_domain = _marker(
        1,
        100,
        400,
        model="VS000808",
        score=20.0,
    )
    capsid = _marker(3, 900, 1200, model="VS000798", score=80.0)
    gene_order = {
        "ctg": [
            _gene(1, 100, 400),
            _gene(2, 500, 800),
            _gene(3, 900, 1200),
        ]
    }

    regions = assemble_compact_cress_regions(
        [rep, duplicate_rep_domain, capsid],
        gene_order,
    )

    assert len(regions) == 1
    assert (regions[0].start, regions[0].end) == (100, 1200)
    assert [marker.query_porf for marker in regions[0].markers] == [
        rep.query_porf,
        capsid.query_porf,
    ]


def test_distant_cress_genes_remain_independent_insertions() -> None:
    left = _marker(1, 100, 400)
    right = _marker(4, 2000, 2300)
    gene_order = {
        "ctg": [
            _gene(1, 100, 400),
            _gene(2, 500, 800),
            _gene(3, 900, 1200),
            _gene(4, 2000, 2300),
        ]
    }

    regions = assemble_compact_cress_regions([left, right], gene_order)

    assert [(region.start, region.end) for region in regions] == [
        (100, 400),
        (2000, 2300),
    ]


def test_phase2_does_not_pad_or_cross_merge_cress_seed() -> None:
    generic = MergedSeed(
        scaffold="ctg",
        start=0,
        end=100,
        predicted_family="NCLDV",
    )
    cress = MergedSeed(
        scaffold="ctg",
        start=100,
        end=400,
        predicted_family="CRESS",
    )
    proteome_index = {
        "ctg": [
            pORF(id="gene_0", scaffold="ctg", start=0, end=100),
            pORF(id="gene_1", scaffold="ctg", start=100, end=400),
            pORF(id="gene_2", scaffold="ctg", start=400, end=700),
        ]
    }

    observed = extend_seeds_by_genes(
        [generic, cress],
        proteome_index,
        extension_genes=5,
    )

    cress_observed = next(
        seed for seed in observed if seed.predicted_family == "CRESS"
    )
    assert (cress_observed.start, cress_observed.end) == (100, 400)
    assert len(observed) == 2


def test_phase2_does_not_recompute_cress_bounds_from_overlapping_porfs() -> None:
    cress = MergedSeed(
        scaffold="ctg",
        start=100,
        end=400,
        predicted_family="CRESS",
    )
    proteome_index = {
        "ctg": [
            pORF(id="gene_1", scaffold="ctg", start=50, end=450),
        ]
    }

    observed = extend_seeds_by_genes(
        [cress],
        proteome_index,
        extension_genes=5,
    )

    assert [(seed.start, seed.end) for seed in observed] == [(100, 400)]


def test_phase2_does_not_merge_overlapping_cress_seeds() -> None:
    seeds = [
        MergedSeed(
            scaffold="ctg",
            start=100,
            end=400,
            predicted_family="CRESS",
        ),
        MergedSeed(
            scaffold="ctg",
            start=350,
            end=700,
            predicted_family="CRESS",
        ),
    ]
    proteome_index = {
        "ctg": [
            pORF(id="gene_1", scaffold="ctg", start=100, end=400),
            pORF(id="gene_2", scaffold="ctg", start=350, end=700),
        ]
    }

    observed = extend_seeds_by_genes(
        seeds,
        proteome_index,
        extension_genes=5,
    )

    assert [(seed.start, seed.end) for seed in observed] == [
        (100, 400),
        (350, 700),
    ]


def test_host_signature_trim_preserves_cress_bounds() -> None:
    cress = MergedSeed(
        scaffold="ctg",
        start=100,
        end=400,
        predicted_family="CRESS",
    )

    observed, summary = trim_seed_by_host_signature(
        seed=cress,
        gene_records=[],
        host_model=HostSignatureModel(),
    )

    assert observed is cress
    assert (observed.start, observed.end) == (100, 400)
    assert summary["reason"] == "cress_exact_boundary"
    assert summary["host_coordinate_change_opportunities"] == 0


def test_post_taxonomy_merge_preserves_cress_boundary() -> None:
    boundaries = [
        RefinedBoundary(
            scaffold="ctg",
            start=0,
            end=250,
            predicted_family="NCLDV",
        ),
        RefinedBoundary(
            scaffold="ctg",
            start=100,
            end=400,
            predicted_family="CRESS",
        ),
    ]

    observed = merge_adjacent_viral_boundaries(
        boundaries,
        taxonomy_map={},
    )

    assert [
        (boundary.start, boundary.end, boundary.predicted_family)
        for boundary in observed
    ] == [
        (0, 250, "NCLDV"),
        (100, 400, "CRESS"),
    ]
