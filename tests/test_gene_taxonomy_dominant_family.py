"""Regression tests for the region dominant-family label.

A region whose genes carry no viral top-10 prefix has no dominant family. Three
summarizers in `gene_taxonomy` previously took `max()` over an all-zero family
count, which returns the first key rather than nothing, so such regions were
labelled NCLDV. That label is exported in
`virosync_predictions_detailed.tsv` and, worse, counts as family support in
`calculate_eve_confidence`: it lifted the priority floor from 0.55 to 0.70 and
promoted MEDIUM calls to HIGH.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    VerificationResult,
    calculate_eve_confidence,
)
from virosync.pipeline.phase3.gene_taxonomy import (
    GeneDiamondHit,
    GeneTaxonomy,
    materialize_gene_taxonomy_batch_from_cached_hits,
    summarize_dominant_family,
)


def _gene(porf_id: str, prefixes: list[str]) -> GeneTaxonomy:
    return GeneTaxonomy(
        porf_id=porf_id,
        porf_start=0,
        porf_end=300,
        top10_prefixes=prefixes,
        top10_pidents=[50.0] * len(prefixes),
    )


def test_no_viral_prefix_has_no_dominant_family() -> None:
    # The two DIAMOND-backed summarizers cannot be driven without a database,
    # so the shared helper is where their fix is provable.
    genes = [_gene("g1", []), _gene("g2", ["EUK"]), _gene("g3", [])]
    assert summarize_dominant_family(genes) == ("UNKNOWN", 0.0)


def test_dominant_family_reports_majority_and_fraction() -> None:
    genes = [_gene("g1", ["PLV"]), _gene("g2", ["PLV"]), _gene("g3", ["NCLDV"])]
    family, fraction = summarize_dominant_family(genes)
    assert family == "PPV"
    assert fraction == 2 / 3


def test_cress_can_be_the_dominant_family() -> None:
    genes = [_gene("g1", ["CRESS"]), _gene("g2", ["CRESS"]), _gene("g3", ["EUK"])]
    assert summarize_dominant_family(genes) == ("CRESS", 2 / 3)


def test_batch_summary_reports_unknown_without_hits(tmp_path: Path) -> None:
    """End-to-end through the cached-hit path: one gene, zero DIAMOND hits."""
    proteome = tmp_path / "proteome.fasta"
    proteome.write_text(">scaf1_1 # 100 # 400 # 1 # ID=1_1;partial=00\nMKV\n")

    results = materialize_gene_taxonomy_batch_from_cached_hits(
        regions=[{"eve_id": "EVE_1", "scaffold": "scaf1", "start": 0, "end": 1000}],
        proteome_fasta=proteome,
        diamond_hits={},
        output_dir=tmp_path / "out",
    )
    taxonomies, summary = results["EVE_1"]
    assert len(taxonomies) == 1
    assert summary["dominant_family"] == "UNKNOWN"
    assert summary["dominant_fraction"] == 0.0


@pytest.mark.parametrize("pident,qualified_count,component", [(24.999, 0, 0.0), (25.0, 1, 1.0)])
def test_cached_gene_taxonomy_keeps_raw_flags_out_of_confidence(
    tmp_path: Path, pident: float, qualified_count: int, component: float
) -> None:
    proteome = tmp_path / "proteome.fasta"
    proteome.write_text(">scaf1_1 # 100 # 400 # 1 # ID=1_1;partial=00\nMKV\n")
    taxonomies, summary = materialize_gene_taxonomy_batch_from_cached_hits(
        regions=[{"eve_id": "EVE_1", "scaffold": "scaf1", "start": 0, "end": 1000}],
        proteome_fasta=proteome,
        diamond_hits={"scaf1_1": [GeneDiamondHit("scaf1_1", "NCLDV__reference|protein", 1e-5, 80, pident, 90)]},
        output_dir=tmp_path / "out",
    )["EVE_1"]
    result = VerificationResult("EVE_1", "scaf1", 0, 1000)
    EvidenceSynthesizer()._process_gene_taxonomy(result, taxonomies, summary)
    calculate_eve_confidence(result, 0)
    assert result.genes_with_ncldv_mirus_top10 == 1
    assert result.gene_taxonomy_viral_top10 == qualified_count
    assert result.score_components["scores"]["gene_taxonomy"] == component
