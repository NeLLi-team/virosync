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

from virosync.pipeline.phase1.hhg_seeding import (
    HMMHit,
    validate_hmm_hits_with_combined_db,
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

    with patch(
        "virosync.pipeline.phase1.hhg_seeding.extract_hmm_hit_sequences",
        side_effect=_fake_extract,
    ), patch(
        "virosync.pipeline.phase1.hhg_seeding.run_diamond_on_hmm_hits",
        side_effect=_fake_run_diamond,
    ), patch(
        "virosync.pipeline.phase1.hhg_seeding.filter_validated_markers",
        return_value=[],
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

    with patch(
        "virosync.pipeline.phase1.hhg_seeding.extract_hmm_hit_sequences",
        side_effect=_fake_extract,
    ), patch(
        "virosync.pipeline.phase1.hhg_seeding.run_diamond_on_hmm_hits",
        side_effect=_fake_run_diamond,
    ), patch(
        "virosync.pipeline.phase1.hhg_seeding.filter_validated_markers",
        side_effect=_capture_filter,
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
