"""Every search that feeds taxonomy or validation receives the run's reference mask."""

from pathlib import Path
from types import SimpleNamespace

import pytest

import virosync.pipeline.phase1 as phase1_package
import virosync.pipeline.phase3 as phase3_package
from virosync.orchestration import resource_monitor, tasks, utils
from virosync.pipeline import search_backend
from virosync.pipeline.phase1.hhg_seeding import _run_diamond_marker_search
from virosync.pipeline.phase1.marker_validation import run_diamond_on_hmm_hits
from virosync.pipeline.phase2 import boundary_diamond
from virosync.pipeline.phase2.boundary_diamond import (
    BoundaryDiamondConfig,
    GenomeDiamondQuery,
    pORF,
)
from virosync.pipeline.phase3 import gene_taxonomy
from virosync.pipeline.phase3.evidence_synthesizer import EvidenceSynthesizer, EvidenceSynthesizerConfig
from virosync.pipeline.phase3.phylogenetic_validation import PhylogeneticValidator

MASK = frozenset({"GCA-000001-1", "EP00001"})


@pytest.fixture
def search_calls(monkeypatch: pytest.MonkeyPatch) -> list[tuple[str, object]]:
    """Record the name and the mask of every call to a search function."""
    calls: list[tuple[str, object]] = []

    def make_fake(name: str):
        def fake(*, output_tsv: Path, **kwargs: object) -> Path:
            calls.append((name, kwargs.get("masked_genomes")))
            output_tsv.write_text("", encoding="utf-8")
            return output_tsv

        return fake

    monkeypatch.setattr(search_backend, "run_sequence_search", make_fake("plain"))
    monkeypatch.setattr(search_backend, "run_genome_filtered_search", make_fake("filtered"))
    return calls


def test_marker_validation_search_receives_the_mask(tmp_path: Path, search_calls: list) -> None:
    run_diamond_on_hmm_hits(tmp_path / "q.faa", tmp_path / "marker.dmnd", tmp_path / "out.tsv", masked_genomes=MASK)

    assert search_calls == [("filtered", MASK)]


def test_hhg_seeding_marker_search_receives_the_mask(tmp_path: Path, search_calls: list) -> None:
    output = tmp_path / "work" / "out.tsv"
    (output.parent / "tmp").mkdir(parents=True)
    (output.parent / "tmp" / "marker.dmnd").write_text("", encoding="utf-8")

    _run_diamond_marker_search(tmp_path / "q.faa", tmp_path / "marker.faa", output, threads=1, masked_genomes=MASK)

    assert search_calls == [("filtered", MASK)]


def test_boundary_search_receives_the_mask(tmp_path: Path, search_calls: list) -> None:
    boundary_diamond.run_diamond_blastp(
        tmp_path / "q.faa", tmp_path / "proteome.dmnd", tmp_path / "out.tsv", masked_genomes=MASK
    )

    assert search_calls == [("filtered", MASK)]


def test_gene_taxonomy_search_receives_the_mask(tmp_path: Path, search_calls: list) -> None:
    gene_taxonomy._run_single_diamond_blastp(
        query_fasta=tmp_path / "q.faa",
        db_prefix=tmp_path / "proteome",
        output_file=tmp_path / "out.tsv",
        threads=1,
        evalue=1e-5,
        max_seqs=10,
        fast=False,
        masked_genomes=MASK,
    )

    assert search_calls == [("filtered", MASK)]


def test_phase3_diamond_check_receives_the_mask(tmp_path: Path, search_calls: list) -> None:
    validator = PhylogeneticValidator(
        genome_path=tmp_path / "genome.fna",
        work_dir=tmp_path / "work",
        diamond_db=tmp_path / "db.dmnd",
        masked_genomes=MASK,
    )

    validator._run_diamond("EVE1", "ATG" + "GCT" * 100 + "TAA")

    assert search_calls == [("filtered", MASK)]


def test_searches_without_a_mask_pass_an_empty_one(tmp_path: Path, search_calls: list) -> None:
    run_diamond_on_hmm_hits(tmp_path / "q.faa", tmp_path / "marker.dmnd", tmp_path / "out.tsv")
    boundary_diamond.run_diamond_blastp(tmp_path / "q.faa", tmp_path / "proteome.dmnd", tmp_path / "out.tsv")

    assert search_calls == [("filtered", frozenset()), ("filtered", frozenset())]


@pytest.fixture
def boundary_calls(monkeypatch: pytest.MonkeyPatch) -> list[dict]:
    """Record the keyword arguments of every boundary-module BLASTP call."""
    calls: list[dict] = []

    def fake_blastp(**kwargs: object) -> None:
        calls.append(kwargs)
        Path(kwargs["output"]).write_text("", encoding="utf-8")

    monkeypatch.setattr(boundary_diamond, "run_diamond_blastp", fake_blastp)
    return calls


def test_full_proteome_search_receives_the_mask(tmp_path: Path, boundary_calls: list[dict]) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">p1\nMPEPTIDE\n", encoding="utf-8")

    boundary_diamond.run_full_proteome_diamond(
        proteome_fasta=proteome,
        diamond_db=tmp_path / "combined.dmnd",
        output_dir=tmp_path / "superset",
        max_target_seqs=10,
        threads=1,
        masked_genomes=MASK,
    )

    assert [call["masked_genomes"] for call in boundary_calls] == [MASK]


def test_chunked_search_gives_every_chunk_the_mask(tmp_path: Path, boundary_calls: list[dict]) -> None:
    query = tmp_path / "q.faa"
    query.write_text(">a\nMPEPTIDE\n>b\nMPEPTIDE\n", encoding="utf-8")

    boundary_diamond.run_diamond_chunked(
        query_fasta=query,
        diamond_db=tmp_path / "combined.dmnd",
        output_dir=tmp_path / "chunks",
        chunk_size=1,
        max_target_seqs=10,
        threads=1,
        masked_genomes=MASK,
    )

    assert [call["masked_genomes"] for call in boundary_calls] == [MASK, MASK]


@pytest.mark.parametrize("chunk_size", [100, 1], ids=["single-run", "chunked"])
def test_batched_boundary_search_uses_the_mask_of_its_config(
    tmp_path: Path, boundary_calls: list[dict], chunk_size: int
) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">p1 # 1 # 90 # + # ID=contig_1;\nMPEPTIDE\n>p2 # 101 # 190 # + # ID=contig_2;\nMPEPTIDE\n",
        encoding="utf-8",
    )
    query = GenomeDiamondQuery(eve_porf_ids={"seed": ["p1", "p2"]}, all_porf_ids=["p1", "p2"])

    boundary_diamond.run_batched_diamond(
        query=query,
        proteome_fasta=proteome,
        diamond_db=tmp_path / "combined.dmnd",
        output_dir=tmp_path / "out",
        proteome_index={"contig": [pORF("p1", "contig", 0, 90), pORF("p2", "contig", 100, 190)]},
        config=BoundaryDiamondConfig(chunk_size=chunk_size, masked_genomes=MASK),
    )

    assert boundary_calls
    assert {call["masked_genomes"] for call in boundary_calls} == {MASK}


def test_gene_taxonomy_batch_search_receives_the_mask(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">contig_1_1 # 1 # 90 # + # ID=1_1;\nMPEPTIDE\n", encoding="utf-8")
    calls: list[dict] = []

    def fake_blastp(**kwargs: object) -> None:
        calls.append(kwargs)
        Path(kwargs["output_file"]).write_text("", encoding="utf-8")

    monkeypatch.setattr(gene_taxonomy, "run_diamond_blastp", fake_blastp)

    gene_taxonomy.run_gene_taxonomy_diamond_batch(
        regions=[{"eve_id": "EVE1", "scaffold": "contig_1", "start": 0, "end": 100}],
        proteome_fasta=proteome,
        combined_faa_db=tmp_path / "combined.dmnd",
        output_dir=tmp_path / "taxonomy",
        masked_genomes=MASK,
    )

    assert [call["masked_genomes"] for call in calls] == [MASK]


def test_the_synthesizer_hands_its_mask_to_the_phylogenetic_validator(tmp_path: Path) -> None:
    config = EvidenceSynthesizerConfig(
        use_phylogenetic_validation=True,
        diamond_db=tmp_path / "db.dmnd",
        masked_genomes=MASK,
    )
    synthesizer = EvidenceSynthesizer(config=config, genome_path=tmp_path / "genome.fna", work_dir=tmp_path / "work")

    assert synthesizer.phylogenetic_validator.masked_genomes == MASK


class _NullResourceMonitor:
    def __init__(self, **_kwargs: object) -> None:
        pass

    def __enter__(self) -> "_NullResourceMonitor":
        return self

    def __exit__(self, *_exc: object) -> bool:
        return False


@pytest.fixture
def null_resource_monitor(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resource_monitor, "ResourceMonitor", _NullResourceMonitor)


def test_marker_validation_task_forwards_the_mask(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, null_resource_monitor: None
) -> None:
    from virosync.pipeline.phase1.hhg_seeding import HMMHit

    calls: list[dict] = []
    monkeypatch.setattr(phase1_package, "extract_hmm_hit_sequences", lambda **_kwargs: 1)

    def marker_search(**kwargs):
        calls.append(kwargs)
        kwargs["output_tsv"].write_text("p1\tEUK__host|protein\t1e-5\t100\t80\t100\n")

    monkeypatch.setattr(phase1_package, "run_diamond_on_hmm_hits", marker_search)
    monkeypatch.setattr(phase1_package, "filter_validated_markers", lambda **_kwargs: [])

    monkeypatch.setattr(
        "virosync.pipeline.phase1.marker_validation.run_host_diamond_on_hmm_hits",
        lambda *args, **kwargs: calls.append(kwargs),
    )
    tasks.marker_validation_task(
        hmm_hits=[HMMHit("p1", "model", 10, 1e-5, 10, 1, 1)],
        proteome_path=tmp_path / "phase0" / "run" / "proteome.faa",
        marker_db=tmp_path / "marker.dmnd",
        output_dir=tmp_path / "validation",
        masked_genomes=MASK,
    )

    assert [call["masked_genomes"] for call in calls] == [MASK, MASK]


def test_hhg_seeding_task_forwards_the_mask(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, null_resource_monitor: None
) -> None:
    calls: list[dict] = []
    monkeypatch.setattr(phase1_package, "hhg_seeding_pipeline", lambda **kwargs: calls.append(kwargs) or ([], []))

    tasks.hhg_seeding_task(
        proteome_path=tmp_path / "phase0" / "run" / "proteome.faa",
        hmm_database=tmp_path / "markers.hmm",
        output_dir=tmp_path / "hhg",
        masked_genomes=MASK,
    )

    assert [call["masked_genomes"] for call in calls] == [MASK]


def test_gene_taxonomy_batch_task_forwards_the_mask(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[dict] = []
    monkeypatch.setattr(gene_taxonomy, "run_gene_taxonomy_diamond_batch", lambda **kwargs: calls.append(kwargs) or {})
    monkeypatch.setattr(
        "virosync.orchestration.utils.run_with_monitor",
        lambda task_name, func, monitor_output_dir, **kwargs: func(**kwargs),
    )

    tasks.gene_taxonomy_batch_task(
        regions=[],
        proteome_path=tmp_path / "proteome.faa",
        gene_taxonomy_db=tmp_path / "combined.dmnd",
        output_dir=tmp_path / "taxonomy",
        masked_genomes=MASK,
    )

    assert [call["masked_genomes"] for call in calls] == [MASK]


class _RecordingConfig:
    """Stand-in for EvidenceSynthesizerConfig that records its keyword arguments."""

    instances: list[dict] = []

    def __init__(self, **kwargs: object) -> None:
        type(self).instances.append(kwargs)


class _StopSynthesizer:
    def __init__(self, **_kwargs: object) -> None:
        raise RuntimeError("stop after the config is built")


@pytest.fixture
def recording_config(monkeypatch: pytest.MonkeyPatch, null_resource_monitor: None) -> type[_RecordingConfig]:
    _RecordingConfig.instances = []
    monkeypatch.setattr(utils, "get_genes_for_boundary", lambda **_kwargs: [])
    monkeypatch.setattr(phase3_package, "EvidenceSynthesizerConfig", _RecordingConfig)
    monkeypatch.setattr(phase3_package, "EvidenceSynthesizer", _StopSynthesizer)
    return _RecordingConfig


def test_verify_eve_task_puts_the_mask_in_the_synthesizer_config(
    tmp_path: Path, recording_config: type[_RecordingConfig]
) -> None:
    boundary = SimpleNamespace(scaffold="scaffold", start=10, end=20)

    with pytest.raises(RuntimeError, match="stop after the config"):
        tasks.verify_eve_task(
            hallmark_evidence=[],
            hallmark_witnesses=[],
            hallmark_support=[],
            boundary=boundary,
            genome_path=tmp_path / "genome.fna",
            work_dir=tmp_path / "work",
            proteome_path=tmp_path / "proteome.faa",
            hallmark_hits=[],
            masked_genomes=MASK,
        )

    assert [kwargs["masked_genomes"] for kwargs in recording_config.instances] == [MASK]


def test_batched_verification_puts_the_mask_in_the_synthesizer_config(
    tmp_path: Path, recording_config: type[_RecordingConfig]
) -> None:
    boundary = SimpleNamespace(scaffold="scaffold", start=10, end=20, candidate_id="EVE_scaffold_10-20")

    with pytest.raises(RuntimeError, match="Failed to verify boundary"):
        tasks.verify_eve_candidates_batched_task(
            boundaries=[boundary],
            genome_path=tmp_path / "genome.fna",
            work_dir=tmp_path / "work",
            proteome_path=tmp_path / "proteome.faa",
            hallmark_hits_map={},
            hallmark_evidence_map={},
            hallmark_witnesses_map={},
            candidate_annotations={},
            max_workers=1,
            masked_genomes=MASK,
        )

    assert [kwargs["masked_genomes"] for kwargs in recording_config.instances] == [MASK]
