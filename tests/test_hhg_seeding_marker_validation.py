"""Regression test for NameError on `genome_fasta` in hhg_seeding.

Past bug: `validate_hmm_hits_with_combined_db` referenced `genome_fasta`
inside its body, but that name was not a parameter of the function. Any
happy-path run that reached the `filter_validated_markers` call raised
`NameError: name 'genome_fasta' is not defined`. The only caller
(`hhg_seeding_pipeline`) already carried `genome_fasta` in its own scope,
so the fix is to thread it through as an optional parameter.

This test exercises the full happy path (marker_db exists, Diamond
succeeds) with stubbed Diamond and validator calls, and asserts the
function returns cleanly instead of NameError-ing.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from virosync.pipeline.phase1.hhg_seeding import (
    HMMHit,
    validate_hmm_hits_with_combined_db,
)
from virosync.pipeline.phase1.marker_validation import (
    NovelMarkerCriteria,
    filter_validated_markers,
)


def _fake_extract(hits, proteome_fasta, output_fasta):
    """Stub for extract_hmm_hit_sequences: create an empty FASTA and claim 1 seq."""
    Path(output_fasta).write_text(">q1\nM\n")
    return 1


def _fake_run_diamond(hmm_hit_fasta, diamond_db, output_tsv, threads, evalue, max_seqs):
    """Stub for run_diamond_on_hmm_hits: create an empty output TSV."""
    Path(output_tsv).write_text("")


def test_validate_hmm_hits_does_not_nameerror_on_happy_path(tmp_path: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">contig_1_1\nMGKPSALVR\n")

    marker_db = tmp_path / "marker.dmnd"
    marker_db.write_bytes(b"")  # presence is enough; Diamond is mocked

    hits = [
        HMMHit(
            query_name="contig_1_1",
            target_name="mcp",
            score=120.0,
            evalue=1e-40,
            domain_score=120.0,
            query_start=1,
            query_end=9,
        )
    ]

    with (
        patch(
            "virosync.pipeline.phase1.hhg_seeding.extract_hmm_hit_sequences",
            side_effect=_fake_extract,
        ),
        patch(
            "virosync.pipeline.phase1.hhg_seeding.run_diamond_on_hmm_hits",
            side_effect=_fake_run_diamond,
        ),
        patch(
            "virosync.pipeline.phase1.hhg_seeding.filter_validated_markers",
            return_value=[],
        ),
    ):
        validated_hits, validated_markers = validate_hmm_hits_with_combined_db(
            hits=hits,
            proteome_fasta=proteome,
            marker_db=marker_db,
            threads=1,
            top_k=10,
            output_dir=tmp_path / "out",
        )

    # No validated markers were returned (our mock returned []),
    # so validated_hits should be empty, not a NameError traceback.
    assert validated_hits == []
    assert validated_markers == []


def test_validate_hmm_hits_forwards_optional_genome_fasta(tmp_path: Path) -> None:
    """After the fix, `genome_fasta` is an optional kwarg and is forwarded to
    `filter_validated_markers`. Verify the forwarding without assuming the
    caller ever sets it."""
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">contig_1_1\nMGKPSALVR\n")
    genome = tmp_path / "genome.fna"
    genome.write_text(">contig_1\nACGT\n")
    marker_db = tmp_path / "marker.dmnd"
    marker_db.write_bytes(b"")

    hits = [
        HMMHit(
            query_name="contig_1_1",
            target_name="mcp",
            score=120.0,
            evalue=1e-40,
            domain_score=120.0,
            query_start=1,
            query_end=9,
        )
    ]

    seen_kwargs: dict = {}

    def _capture_filter(**kwargs):
        seen_kwargs.update(kwargs)
        return []

    with (
        patch(
            "virosync.pipeline.phase1.hhg_seeding.extract_hmm_hit_sequences",
            side_effect=_fake_extract,
        ),
        patch(
            "virosync.pipeline.phase1.hhg_seeding.run_diamond_on_hmm_hits",
            side_effect=_fake_run_diamond,
        ),
        patch(
            "virosync.pipeline.phase1.hhg_seeding.filter_validated_markers",
            side_effect=_capture_filter,
        ),
    ):
        validate_hmm_hits_with_combined_db(
            hits=hits,
            proteome_fasta=proteome,
            marker_db=marker_db,
            threads=1,
            top_k=10,
            output_dir=tmp_path / "out",
            genome_fasta=genome,
        )

    assert seen_kwargs.get("genome_fasta") == genome


def test_filter_uses_prodigal_coordinates_without_reading_genome(tmp_path: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">contig_1_1 # 1 # 300 # + # ID=1_1;partial=00\n"
        f"{'M' * 100}\n"
        ">contig_1_2 # 301 # 600 # + # ID=1_2;partial=00\n"
        f"{'M' * 100}\n"
        ">contig_2_1 # 1 # 300 # + # ID=2_1;partial=00\n"
        f"{'M' * 100}\n"
    )
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("")
    missing_genome = tmp_path / "not-read.fna"
    hits = [HMMHit(query, "mcp", 120.0, 1e-40, 120.0, 1, 100) for query in ("contig_1_1", "contig_1_2", "contig_2_1")]

    markers = filter_validated_markers(
        hits,
        diamond,
        proteome,
        genome_fasta=missing_genome,
        novel_criteria=NovelMarkerCriteria(
            min_hmm_coverage=0.3,
            require_cluster=False,
        ),
    )

    assert [marker.validation_status for marker in markers] == [
        "validated_novel",
        "validated_novel",
        "validated_novel",
    ]
    assert [(marker.scaffold, marker.start, marker.end) for marker in markers] == [
        ("contig_1", 0, 300),
        ("contig_1", 300, 600),
        ("contig_2", 0, 300),
    ]


def test_novel_coverage_uses_full_protein_length(tmp_path: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">contig_1_1 # 1 # 300 # + # ID=1_1;partial=00\n"
        f"{'M' * 100}\n"
        ">contig_1_2 # 301 # 3300 # + # ID=1_2;partial=00\n"
        f"{'M' * 1000}\n"
    )
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("")
    hits = [
        HMMHit("contig_1_1", "mcp", 80.0, 1e-20, 80.0, 1, 100),
        HMMHit("contig_1_2", "mcp", 80.0, 1e-20, 80.0, 1, 150),
    ]

    markers = filter_validated_markers(
        hits,
        diamond,
        proteome,
        novel_criteria=NovelMarkerCriteria(require_cluster=False),
    )

    assert [marker.validation_status for marker in markers] == [
        "validated_novel",
        "unvalidated",
    ]


def test_novel_marker_with_unknown_protein_length_is_unvalidated(tmp_path: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">other\nMMMM\n")
    genome = tmp_path / "genome.fna"
    genome.write_text(f">contig_1\n{'A' * 300}\n")
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("")
    hit = HMMHit("contig_1_frame=1", "mcp", 80.0, 1e-20, 80.0, 1, 100)

    markers = filter_validated_markers(
        [hit],
        diamond,
        proteome,
        genome_fasta=genome,
        novel_criteria=NovelMarkerCriteria(
            min_hmm_coverage=0.0,
            require_cluster=False,
        ),
    )

    assert markers[0].validation_status == "unvalidated"


def test_nearby_diamond_validated_protein_supports_novel_marker_in_any_order(
    tmp_path: Path,
) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">contig_1_1 # 1 # 300 # + # ID=1_1;partial=00\n"
        f"{'M' * 100}\n"
        ">contig_1_2 # 10001 # 10300 # + # ID=1_2;partial=00\n"
        f"{'M' * 100}\n"
    )
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text(
        "".join(f"contig_1_2\tEUK__host_{index}\t1e-40\t{300 - index}\t80\t100\n" for index in range(10))
        + "contig_1_2\tNCLDV__reference\t1e-40\t200\t35\t100\n"
    )
    novel = HMMHit("contig_1_1", "mcp", 80.0, 1e-20, 80.0, 1, 100)
    validated = HMMHit("contig_1_2", "polb", 80.0, 1e-20, 80.0, 1, 100)

    forward = filter_validated_markers([novel, validated], diamond, proteome, max_seqs=11)
    reverse = filter_validated_markers([validated, novel], diamond, proteome, max_seqs=11)

    expected = {
        "contig_1_1|aa1-100": "validated_novel",
        "contig_1_2|aa1-100": "validated",
    }
    assert {marker.query_porf: marker.validation_status for marker in forward} == expected
    assert {marker.query_porf: marker.validation_status for marker in reverse} == expected


@pytest.mark.parametrize(
    ("support_query", "support_start", "support_pident"),
    [
        ("contig_1_2", 10302, 35.0),
        ("contig_1_2", 10001, 24.9),
        ("contig_2_1", 10001, 35.0),
    ],
    ids=["distant", "weak", "cross-scaffold"],
)
def test_invalid_cluster_support_does_not_validate_novel_marker(
    tmp_path: Path,
    support_query: str,
    support_start: int,
    support_pident: float,
) -> None:
    support_index = support_query.rsplit("_", 1)[1]
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">contig_1_1 # 1 # 300 # + # ID=1_1;partial=00\n"
        f"{'M' * 100}\n"
        f">{support_query} # {support_start} # {support_start + 299} # + # "
        f"ID=2_{support_index};partial=00\n"
        f"{'M' * 100}\n"
    )
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text(f"{support_query}\tNCLDV__reference\t1e-40\t200\t{support_pident}\t100\n")
    hits = [
        HMMHit("contig_1_1", "mcp", 80.0, 1e-20, 80.0, 1, 100),
        HMMHit(support_query, "polb", 80.0, 1e-20, 80.0, 1, 100),
    ]

    markers = filter_validated_markers(hits, diamond, proteome)

    assert markers[0].validation_status == "unvalidated"


def test_nearby_hmm_only_hits_do_not_support_each_other(tmp_path: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">contig_1_1 # 1 # 300 # + # ID=1_1;partial=00\n"
        f"{'M' * 100}\n"
        ">contig_1_2 # 301 # 600 # + # ID=1_2;partial=00\n"
        f"{'M' * 100}\n"
    )
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("")
    hits = [HMMHit(query, "mcp", 80.0, 1e-20, 80.0, 1, 100) for query in ("contig_1_1", "contig_1_2")]

    markers = filter_validated_markers(hits, diamond, proteome)

    assert [marker.validation_status for marker in markers] == [
        "unvalidated",
        "unvalidated",
    ]


def test_second_hmm_segment_on_same_protein_cannot_support_novel_segment(
    tmp_path: Path,
) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(f">contig_1_1 # 1 # 300 # + # ID=1_1;partial=00\n{'M' * 100}\n")
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("contig_1_1|aa61-100\tNCLDV__reference\t1e-40\t200\t35\t100\n")
    hits = [
        HMMHit("contig_1_1", "mcp", 80.0, 1e-20, 80.0, 1, 60),
        HMMHit("contig_1_1", "polb", 80.0, 1e-20, 80.0, 61, 100),
    ]

    markers = filter_validated_markers(hits, diamond, proteome)

    assert [marker.validation_status for marker in markers] == [
        "unvalidated",
        "validated",
    ]


def test_filter_reads_genome_for_frame_coordinates(tmp_path: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">contig_1_frame=1\nMMMMMMMMMM\n")
    genome = tmp_path / "genome.fna"
    genome.write_text(f">contig_1\n{'A' * 300}\n")
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("")
    hit = HMMHit("contig_1_frame=1", "mcp", 120.0, 1e-40, 120.0, 1, 10)

    markers = filter_validated_markers(
        [hit],
        diamond,
        proteome,
        genome_fasta=genome,
        novel_criteria=NovelMarkerCriteria(require_cluster=False),
    )

    assert len(markers) == 1
    assert markers[0].validation_status == "validated_novel"
    assert (markers[0].scaffold, markers[0].start, markers[0].end, markers[0].strand) == (
        "contig_1",
        0,
        30,
        "+",
    )
