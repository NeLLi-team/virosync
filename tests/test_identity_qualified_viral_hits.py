from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from virosync.pipeline.phase1.taxonomy_expansion import (
    is_gene_viral_positive,
    parse_diamond_top10,
)
from virosync.pipeline.phase2.boundary_diamond import has_identity_qualified_viral_hit
from virosync.pipeline.phase3.gene_taxonomy import GeneDiamondHit, classify_gene_taxonomy


def _tax(**overrides):
    values = {
        "has_hit": True,
        "top1_prefix": "EUK__",
        "top10_prefixes": [],
        "top10_pidents": [],
        "top10_bits": [],
        "taxonomy_fingerprint": None,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def test_single_identity_qualified_cress_hit_counts_as_strong_viral() -> None:
    tax = _tax(
        top10_prefixes=["EUK__", "CRESS__"],
        top10_pidents=[94.0, 25.1],
    )

    assert has_identity_qualified_viral_hit(tax.top10_prefixes, tax.top10_pidents)


def test_single_low_identity_viral_hit_does_not_count() -> None:
    tax = _tax(
        top10_prefixes=["CRESS__"],
        top10_pidents=[24.9],
    )

    assert not has_identity_qualified_viral_hit(tax.top10_prefixes, tax.top10_pidents)


def test_taxonomy_expansion_viral_positive_uses_pident_threshold() -> None:
    assert is_gene_viral_positive([("EUK__host|p1", 95.0), ("CRESS__rep|p2", 25.0)])
    assert not is_gene_viral_positive([("CRESS__rep|p2", 24.9)])
    assert not is_gene_viral_positive(["CRESS__legacy_without_pident"])


def test_parse_diamond_top10_preserves_percent_identity(tmp_path: Path) -> None:
    diamond_tsv = tmp_path / "hits.tsv"
    diamond_tsv.write_text("gene1\tEUK__host|p1\t1e-50\t200\t91.5\t98\ngene1\tCRESS__rep|p2\t1e-20\t120\t35.2\t80\n")

    hits = parse_diamond_top10(diamond_tsv)

    assert hits["gene1"] == [("EUK__host|p1", 91.5), ("CRESS__rep|p2", 35.2)]
    assert is_gene_viral_positive(hits["gene1"])


def test_phase3_gene_taxonomy_uses_identity_qualified_viral_hit() -> None:
    high_identity = classify_gene_taxonomy(
        "porf1",
        [GeneDiamondHit("porf1", "CRESS__rep|p1", 1e-20, 100.0, 25.0, 80.0)],
        start=1,
        end=100,
    )
    low_identity = classify_gene_taxonomy(
        "porf2",
        [GeneDiamondHit("porf2", "CRESS__rep|p2", 1e-20, 100.0, 24.9, 80.0)],
        start=1,
        end=100,
    )

    assert high_identity.has_viral is True
    assert low_identity.has_viral is False
