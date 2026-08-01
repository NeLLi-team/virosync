from __future__ import annotations

from pathlib import Path

import pytest

from virosync.pipeline.phase1.hhg_seeding import HMMHit
from virosync.pipeline.phase1.marker_validation import (
    NovelMarkerCriteria,
    ValidationStatus,
    filter_validated_markers,
    validate_hmm_hit,
)
from virosync.pipeline.phase3.gene_taxonomy import extract_prefix


def _write_proteome(path: Path, query: str = "contig_1_1") -> None:
    path.write_text(
        f">{query} # 1 # 300 # + # ID=1_1;partial=00\n"
        f"{'M' * 100}\n"
    )


def _hmm_hit(query: str = "contig_1_1") -> HMMHit:
    return HMMHit(
        query_name=query,
        target_name="VS000803",
        score=120.0,
        evalue=1e-40,
        domain_score=120.0,
        query_start=1,
        query_end=100,
    )


def _write_diamond(path: Path, query: str, targets: list[tuple[str, float, float]]) -> None:
    lines = [
        f"{query}\t{target}\t1e-40\t{bits}\t{pident}\t95.0\n"
        for target, bits, pident in targets
    ]
    path.write_text("".join(lines))


@pytest.mark.parametrize("prefix", ["NCLDV__", "MIRUS__", "PLV__", "VP__", "CRESS__"])
def test_single_validated_prefix_top10_hit_above_identity_threshold_validates_marker(
    tmp_path: Path,
    prefix: str,
) -> None:
    proteome = tmp_path / "proteome.faa"
    diamond = tmp_path / "diamond.tsv"
    _write_proteome(proteome)
    _write_diamond(
        diamond,
        "contig_1_1",
        [
            (f"{prefix}IMGVR_UViG_1|gene_1", 200.0, 85.0),
            ("EUK__host_1", 190.0, 70.0),
            ("EUK__host_2", 180.0, 70.0),
        ],
    )

    markers = filter_validated_markers([_hmm_hit()], diamond, proteome)

    assert len(markers) == 1
    assert markers[0].validation_status == "validated"
    assert markers[0].has_viral == 1


def test_single_cress_top10_hit_below_identity_threshold_does_not_validate_marker(tmp_path: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    diamond = tmp_path / "diamond.tsv"
    _write_proteome(proteome)
    _write_diamond(
        diamond,
        "contig_1_1",
        [
            ("CRESS__IMGVR_UViG_1|gene_1", 200.0, 24.9),
            ("EUK__host_1", 190.0, 70.0),
            ("EUK__host_2", 180.0, 70.0),
        ],
    )

    markers = filter_validated_markers([_hmm_hit()], diamond, proteome)

    assert len(markers) == 1
    assert markers[0].validation_status == "unvalidated"
    assert markers[0].has_viral == 1


def test_cress_identity_is_not_rounded_across_the_validation_threshold(
    tmp_path: Path,
) -> None:
    proteome = tmp_path / "proteome.faa"
    diamond = tmp_path / "diamond.tsv"
    _write_proteome(proteome)
    _write_diamond(
        diamond,
        "contig_1_1",
        [
            ("CRESS__IMGVR_UViG_1|gene_1", 210.0, 24.96),
            ("NCLDV__virus_1", 200.0, 30.0),
            ("EUK__host_1", 190.0, 70.0),
        ],
    )

    markers = filter_validated_markers([_hmm_hit()], diamond, proteome)

    assert markers[0].validation_status == "validated"
    assert markers[0].top10_pidents.split(",")[0] == "24.96"


def test_single_phage_top10_hit_does_not_validate_marker(tmp_path: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    diamond = tmp_path / "diamond.tsv"
    _write_proteome(proteome)
    _write_diamond(
        diamond,
        "contig_1_1",
        [
            ("PHAGE__IMGVR_UViG_1|gene_1", 200.0, 85.0),
            ("EUK__host_1", 190.0, 70.0),
            ("EUK__host_2", 180.0, 70.0),
        ],
    )

    markers = filter_validated_markers([_hmm_hit()], diamond, proteome)

    assert len(markers) == 1
    assert markers[0].validation_status == "unvalidated"
    assert markers[0].has_viral == 1


def test_cress_prefix_is_supported_by_taxonomy_helpers() -> None:
    assert extract_prefix("CRESS__IMGVR_UViG_1|gene_1") == "CRESS"
    assert (
        validate_hmm_hit(
            hmm_score=100.0,
            hmm_coverage=0.9,
            diamond_hits=[
                ("CRESS__IMGVR_UViG_1|gene_1", 200.0, 85.0, 1e-40),
            ],
            novel_criteria=NovelMarkerCriteria(),
            has_nearby_markers=False,
        )
        == ValidationStatus.VALIDATED
    )
