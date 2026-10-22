"""Each DIAMOND caller uses the search that matches its database."""

from pathlib import Path

import pytest

from virosync.pipeline import search_backend
from virosync.pipeline.phase1.hhg_seeding import _run_diamond_marker_search
from virosync.pipeline.phase1.marker_validation import run_diamond_on_hmm_hits
from virosync.pipeline.phase2.boundary_diamond import run_diamond_blastp as run_boundary_blastp
from virosync.pipeline.phase3.gene_taxonomy import _run_single_diamond_blastp
from virosync.pipeline.phase3.phylogenetic_validation import PhylogeneticValidator


@pytest.fixture
def searches(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Record which search function each caller picks."""
    calls: list[str] = []

    def make_fake(name: str):
        def fake(*, output_tsv: Path, **_kwargs: object) -> Path:
            calls.append(name)
            output_tsv.write_text("", encoding="utf-8")
            return output_tsv

        return fake

    monkeypatch.setattr(search_backend, "run_sequence_search", make_fake("plain"))
    monkeypatch.setattr(search_backend, "run_genome_filtered_search", make_fake("filtered"))
    return calls


def test_marker_validation_search_is_unfiltered(tmp_path: Path, searches: list[str]) -> None:
    run_diamond_on_hmm_hits(tmp_path / "q.faa", tmp_path / "marker.dmnd", tmp_path / "out.tsv")

    assert searches == ["plain"]


def test_hhg_seeding_marker_search_is_unfiltered(tmp_path: Path, searches: list[str]) -> None:
    output = tmp_path / "work" / "out.tsv"
    (output.parent / "tmp").mkdir(parents=True)
    (output.parent / "tmp" / "marker.dmnd").write_text("", encoding="utf-8")

    _run_diamond_marker_search(tmp_path / "q.faa", tmp_path / "marker.faa", output, threads=1)

    assert searches == ["plain"]


def test_boundary_search_is_genome_filtered(tmp_path: Path, searches: list[str]) -> None:
    run_boundary_blastp(tmp_path / "q.faa", tmp_path / "proteome.dmnd", tmp_path / "out.tsv")

    assert searches == ["filtered"]


def test_gene_taxonomy_search_is_genome_filtered(tmp_path: Path, searches: list[str]) -> None:
    _run_single_diamond_blastp(
        query_fasta=tmp_path / "q.faa",
        db_prefix=tmp_path / "proteome",
        output_file=tmp_path / "out.tsv",
        threads=1,
        evalue=1e-5,
        max_seqs=10,
        fast=False,
    )

    assert searches == ["filtered"]


def test_phase3_diamond_check_is_genome_filtered(tmp_path: Path, searches: list[str]) -> None:
    validator = PhylogeneticValidator(
        genome_path=tmp_path / "genome.fna",
        work_dir=tmp_path / "work",
        diamond_db=tmp_path / "db.dmnd",
    )

    validator._run_diamond("EVE1", "ATG" + "GCT" * 100 + "TAA")

    assert searches == ["filtered"]
