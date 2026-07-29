from __future__ import annotations

import concurrent.futures
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

import virosync.pipeline.phase3 as phase3_package
from virosync.orchestration import (
    resource_monitor,
    tasks as orchestration_tasks,
    utils as orchestration_utils,
)
from virosync.orchestration.resource_monitor import ResourceMonitor
from virosync.pipeline import search_backend
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.phase3 import gene_taxonomy, phylogenetic_validation
from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    EvidenceSynthesizerConfig,
    VerificationResult,
)
from virosync.pipeline.phase3.gvclass_runner import (
    load_gvclass_id_map,
    parse_gvclass_results,
    write_gvclass_results_tsv,
)
from virosync.pipeline.phase3.phylogenetic_validation import PhylogeneticValidator
from virosync.pipeline.phase3.structural_homology import BoltzFoldSeekAnalyzer
from virosync.utils.path_safety import safe_filename_component, safe_filename_components


def _write_proteome(path: Path) -> None:
    path.write_text(">contig_1_1 # 1 # 90 # + # ID=1_1;\nMPEPTIDE\n")


def test_single_gene_taxonomy_encodes_raw_eve_paths(
    tmp_path: Path,
    monkeypatch,
) -> None:
    raw_id = "../EVE_NODE/1"
    component = safe_filename_component(raw_id)
    proteome = tmp_path / "proteome.faa"
    _write_proteome(proteome)
    output_dir = tmp_path / "taxonomy"
    parent_sentinel = tmp_path / "sentinel.txt"
    parent_sentinel.write_bytes(b"parent sentinel\n")
    diamond_paths: list[Path] = []

    def _fake_diamond(**kwargs) -> None:
        output_file = Path(kwargs["output_file"])
        diamond_paths.append(output_file)
        output_file.write_text("")

    monkeypatch.setattr(gene_taxonomy, "run_diamond_blastp", _fake_diamond)

    taxonomies, summary = gene_taxonomy.run_gene_taxonomy_diamond(
        eve_id=raw_id,
        scaffold="contig_1",
        start=0,
        end=100,
        proteome_fasta=proteome,
        combined_faa_db=tmp_path / "database.dmnd",
        output_dir=output_dir,
    )

    assert len(taxonomies) == 1
    assert summary["total"] == 1
    assert diamond_paths == [output_dir / f"{component}_diamond.tsv"]
    assert (output_dir / f"{component}.tsv").is_file()
    assert not (tmp_path / "EVE_NODE").exists()
    assert parent_sentinel.read_bytes() == b"parent sentinel\n"


def test_batch_gene_taxonomy_encodes_paths_but_preserves_raw_result_key(
    tmp_path: Path,
    monkeypatch,
) -> None:
    raw_id = "EVE_../../escape"
    component = safe_filename_component(raw_id)
    proteome = tmp_path / "proteome.faa"
    _write_proteome(proteome)
    output_dir = tmp_path / "taxonomy"

    def _fake_diamond(**kwargs) -> None:
        Path(kwargs["output_file"]).write_text("")

    monkeypatch.setattr(gene_taxonomy, "run_diamond_blastp", _fake_diamond)

    results = gene_taxonomy.run_gene_taxonomy_diamond_batch(
        regions=[
            {
                "eve_id": raw_id,
                "scaffold": "contig_1",
                "start": 0,
                "end": 100,
            }
        ],
        proteome_fasta=proteome,
        combined_faa_db=tmp_path / "database.dmnd",
        output_dir=output_dir,
    )

    assert list(results) == [raw_id]
    assert (output_dir / f"{component}.tsv").is_file()
    assert not (tmp_path / "escape").exists()


def test_batch_gene_taxonomy_rejects_duplicate_raw_ids_before_writing(
    tmp_path: Path,
) -> None:
    raw_id = "EVE_duplicate"
    output_dir = tmp_path / "taxonomy"
    regions = [
        {"eve_id": raw_id, "scaffold": "contig_1", "start": index, "end": index + 10}
        for index in range(2)
    ]

    with pytest.raises(ValueError, match="EVE_duplicate.*indices 0, 1"):
        gene_taxonomy.run_gene_taxonomy_diamond_batch(
            regions=regions,
            proteome_fasta=tmp_path / "proteome.faa",
            combined_faa_db=tmp_path / "database.dmnd",
            output_dir=output_dir,
        )

    assert output_dir.exists() is False


def test_boundary_work_dir_preflight_encodes_scaffold(
    tmp_path: Path,
) -> None:
    boundary = SimpleNamespace(scaffold="../NODE/1", start=10, end=20)
    raw_boundary_id = f"{boundary.scaffold}_{boundary.start}_{boundary.end}"

    work_dirs = orchestration_tasks._preflight_boundary_work_dirs(
        [boundary],
        tmp_path / "work",
    )

    expected_component = safe_filename_component(f"eve_{raw_boundary_id}")
    assert work_dirs[raw_boundary_id] == tmp_path / "work" / expected_component
    assert (tmp_path / "work").resolve() in work_dirs[raw_boundary_id].resolve().parents


def test_boundary_batch_rejects_duplicates_before_output_or_executor(
    tmp_path: Path,
    monkeypatch,
) -> None:
    boundaries = [
        SimpleNamespace(scaffold="NODE/1", start=10, end=20),
        SimpleNamespace(scaffold="NODE/1", start=10, end=20),
    ]
    work_dir = tmp_path / "work"
    sentinel = tmp_path / "sentinel.txt"
    sentinel.write_bytes(b"parent sentinel\n")
    executor_created = False

    def _unexpected_executor(*args, **kwargs):
        nonlocal executor_created
        executor_created = True
        raise AssertionError("executor must not be created for duplicate boundaries")

    monkeypatch.setattr(concurrent.futures, "ThreadPoolExecutor", _unexpected_executor)

    with pytest.raises(ValueError, match="duplicate boundary work IDs.*indices 0, 1"):
        orchestration_tasks.verify_eve_candidates_batched_task(
            boundaries=boundaries,
            genome_path=tmp_path / "genome.fna",
            work_dir=work_dir,
            proteome_path=tmp_path / "proteome.faa",
            hallmark_hits_map={},
        )

    assert executor_created is False
    assert work_dir.exists() is False
    assert sentinel.read_bytes() == b"parent sentinel\n"


def test_single_boundary_task_executes_without_batch_preflight_name_error(
    tmp_path: Path,
    monkeypatch,
) -> None:
    class _FakeResourceMonitor:
        def __init__(self, **kwargs) -> None:
            pass

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, traceback) -> bool:
            return False

    class _FakeConfig:
        def __init__(self, **kwargs) -> None:
            pass

    class _FakeSynthesizer:
        def __init__(self, **kwargs) -> None:
            pass

        def verify_eve(self, **kwargs) -> VerificationResult:
            boundary = kwargs["refined_boundary"]
            return VerificationResult(
                eve_id=f"EVE_{boundary.scaffold}_{boundary.start}-{boundary.end}",
                scaffold=boundary.scaffold,
                start=boundary.start,
                end=boundary.end,
                confidence_tier="MEDIUM",
                final_confidence=0.5,
            )

    monkeypatch.setattr(resource_monitor, "ResourceMonitor", _FakeResourceMonitor)
    monkeypatch.setattr(orchestration_utils, "get_genes_for_boundary", lambda **kwargs: [])
    monkeypatch.setattr(phase3_package, "EvidenceSynthesizerConfig", _FakeConfig)
    monkeypatch.setattr(phase3_package, "EvidenceSynthesizer", _FakeSynthesizer)
    boundary = RefinedBoundary(scaffold="scaffold", start=10, end=20)

    result = orchestration_tasks.verify_eve_task(
        boundary=boundary,
        genome_path=tmp_path / "genome.fna",
        work_dir=tmp_path / "single_work",
        proteome_path=tmp_path / "proteome.faa",
        hallmark_hits=[],
    )

    assert result.eve_id == "EVE_scaffold_10-20"
    assert (tmp_path / "single_work").is_dir()


def test_resource_metrics_filenames_encode_and_disambiguate_raw_task_ids(
    tmp_path: Path,
    monkeypatch,
) -> None:
    raw_task_ids = ["sample_NODE/a/b_10_20", "sample_NODE/a|b_10_20"]
    sentinel = tmp_path / "sentinel.txt"
    sentinel.write_bytes(b"parent sentinel\n")
    monkeypatch.setattr(resource_monitor.time, "time", lambda: 1234.0)

    for raw_task_id in raw_task_ids:
        monitor = ResourceMonitor(
            task_name="verify_eve",
            genome_id="sample",
            phase="phase3",
            output_dir=tmp_path / "output",
            threads=1,
            task_id=raw_task_id,
        )
        monitor._start_time = 0.0
        monitor._write_metrics()

    metric_files = list((tmp_path / "output" / "resource_metrics").glob("*.json"))
    assert len(metric_files) == 2
    assert {json.loads(path.read_text())["task_id"] for path in metric_files} == set(
        raw_task_ids
    )
    assert sentinel.read_bytes() == b"parent sentinel\n"


def test_optional_boltz_work_dir_encodes_raw_eve_id(tmp_path: Path) -> None:
    class _FakeBoltzAnalyzer:
        def __init__(self) -> None:
            self.work_dirs: list[Path] = []

        def analyze_batch(self, sequences, work_dir):
            self.work_dirs.append(Path(work_dir))
            return []

    raw_id = "EVE_../../NODE/1"
    analyzer = _FakeBoltzAnalyzer()
    synthesizer = EvidenceSynthesizer(
        config=EvidenceSynthesizerConfig(
            use_boltz=True,
            use_phylogenetic_validation=False,
        ),
        work_dir=tmp_path / "work",
    )
    synthesizer._boltz_analyzer = analyzer
    result = VerificationResult(
        eve_id=raw_id,
        scaffold="NODE/1",
        start=10,
        end=20,
    )

    synthesizer._run_tiebreakers(
        result=result,
        refined_boundary=RefinedBoundary(scaffold="NODE/1", start=10, end=20),
        window_features=[],
        hallmark_hits=[{"hallmark_gene": "mcp", "porf_id": "p1"}],
        novelty_scores=None,
        porf_sequences=[("p1", "M" * 120)],
    )

    expected = tmp_path / "work" / safe_filename_component(raw_id) / "structures"
    assert analyzer.work_dirs == [expected]
    assert (tmp_path / "work").resolve() in expected.resolve().parents
    assert not (tmp_path / "NODE").exists()


def test_boltz_yaml_batch_encodes_aliases_and_rejects_duplicates(
    tmp_path: Path,
    monkeypatch,
) -> None:
    raw_ids = ["a/b", "a|b", "../a", "a\\b", "a\nb"]
    components = safe_filename_components(raw_ids, label="protein ID")
    analyzer = BoltzFoldSeekAnalyzer(
        viral_db_path=tmp_path / "structures.db",
        min_seq_len=1,
        max_seq_len=100,
    )
    monkeypatch.setattr(analyzer, "available", lambda: True)
    monkeypatch.setattr(analyzer, "_run_boltz", lambda yaml_dir, output_dir: False)
    work_dir = tmp_path / "boltz"

    assert analyzer.analyze_batch([(raw_id, "M" * 10) for raw_id in raw_ids], work_dir) == []

    yaml_dir = work_dir / "boltz_yaml"
    assert {path.name for path in yaml_dir.glob("*.yaml")} == {
        f"{component}.yaml" for component in components.values()
    }
    assert not (tmp_path / "a").exists()

    duplicate_work_dir = tmp_path / "duplicate_boltz"
    with pytest.raises(ValueError, match="duplicate protein IDs.*indices 0, 1"):
        analyzer.analyze_batch(
            [("duplicate", "M" * 10), ("duplicate", "M" * 10)],
            duplicate_work_dir,
        )
    assert duplicate_work_dir.exists() is False


def test_phylogenetic_gvclass_paths_are_encoded_and_contained(
    tmp_path: Path,
    monkeypatch,
) -> None:
    raw_id = "../EVE_NODE/1"
    component = safe_filename_component(raw_id)
    genome = tmp_path / "genome.fna"
    genome.write_text(">contig_1\n" + "ATG" * 100 + "\n")
    validator = PhylogeneticValidator(genome_path=genome, work_dir=tmp_path / "work")
    commands: list[list[str]] = []

    def _fake_run(command, **kwargs):
        commands.append(command)
        return SimpleNamespace(returncode=0, stderr="")

    monkeypatch.setattr(phylogenetic_validation.subprocess, "run", _fake_run)

    assert validator._run_gvclass(raw_id, "ATG" * 100, 0, 300) is None

    input_dir = validator.gvclass_dir / component
    assert commands[0][1] == str(input_dir)
    assert commands[0][3] == str(validator.gvclass_dir / f"{component}_output")
    assert (input_dir / f"{component}.fna").is_file()
    assert not (tmp_path / "EVE_NODE").exists()


def test_gvclass_manifest_round_trip_restores_raw_ids(tmp_path: Path) -> None:
    raw_ids = ["EVE_NODE/1", "EVE_λ"]
    components = safe_filename_components(raw_ids, label="EVE ID")
    manifest = tmp_path / "gvclass_input" / "manifest.tsv"
    manifest.parent.mkdir()
    with manifest.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["eve_id", "nucleotide_fasta", "protein_fasta", "confidence"])
        for raw_id in raw_ids:
            writer.writerow(
                [
                    raw_id,
                    f"nucleotide/{components[raw_id]}.fna",
                    "",
                    "0.9000",
                ]
            )

    summary = tmp_path / "gvclass_summary.tsv"
    summary.write_text(
        "file\tdomain\tgvog_count\tmcp_count\tmirus_count\n"
        + "\n".join(
            f"/tmp/gvclass/{components[raw_id]}.fna\tNCLDV\t3\t2\t1"
            for raw_id in raw_ids
        )
        + "\n"
    )

    id_map = load_gvclass_id_map(manifest)
    results = parse_gvclass_results(summary, id_map=id_map)
    output_tsv = write_gvclass_results_tsv(results, tmp_path / "gvclass_results.tsv")

    assert list(results) == raw_ids
    with output_tsv.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["eve_id"] for row in rows] == raw_ids
    assert list(parse_gvclass_results(summary)) == [
        f"/tmp/gvclass/{components[raw_id]}" for raw_id in raw_ids
    ]


@pytest.mark.filterwarnings("ignore:Partial codon")
def test_phylogenetic_diamond_paths_are_encoded_and_raw_id_is_preserved(
    tmp_path: Path,
    monkeypatch,
) -> None:
    raw_id = "EVE_../../escape"
    component = safe_filename_component(raw_id)
    genome = tmp_path / "genome.fna"
    genome.write_text(">contig_1\n" + "ATG" * 100 + "\n")
    validator = PhylogeneticValidator(
        genome_path=genome,
        work_dir=tmp_path / "work",
        diamond_db=tmp_path / "database.dmnd",
    )
    calls: list[dict] = []

    def _fake_search(**kwargs) -> None:
        calls.append(kwargs)
        Path(kwargs["output_tsv"]).write_text("")

    monkeypatch.setattr(search_backend, "run_sequence_search", _fake_search)

    result = validator._run_diamond(raw_id, "ATG" * 100)

    assert result is not None
    assert result.eve_id == raw_id
    assert calls[0]["query_fasta"] == validator.diamond_dir / f"{component}.faa"
    assert calls[0]["output_tsv"] == validator.diamond_dir / f"{component}.diamond.tsv"
    assert not (tmp_path / "escape").exists()
