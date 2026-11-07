"""Terminal-repeat planning in forkserver workers returns exactly the serial result."""

from __future__ import annotations

import os
import pickle
import random
import subprocess
import sys
import threading
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from types import SimpleNamespace
from typing import ClassVar

import pytest
from Bio import SeqIO
from test_terminal_repeats import _adjacent_repeat_sequence, _Marker, _planted_repeat

from virosync.pipeline.phase1.seed_merger import ConfirmedSeed
from virosync.pipeline.phase2 import terminal_repeats
from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy, GenomeDiamondQuery, SeedGeneMapping, pORF
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.phase2.terminal_repeats import refine_boundaries_with_terminal_repeats


@dataclass(frozen=True, slots=True)
class _Fixture:
    genome_path: Path
    boundaries: list[RefinedBoundary]
    allowed_intervals: list[tuple[int, int]]
    markers: list[_Marker]
    query: GenomeDiamondQuery
    taxonomy_map: dict[str, GeneTaxonomy]
    proteome_index: dict[str, list[pORF]]

    def refine(self, *, workers: int) -> list[RefinedBoundary]:
        return refine_boundaries_with_terminal_repeats(
            self.boundaries,
            allowed_intervals=self.allowed_intervals,
            raw_genome_path=self.genome_path,
            validated_markers=self.markers,
            boundary_diamond_query=self.query,
            boundary_taxonomy_map=self.taxonomy_map,
            proteome_index=self.proteome_index,
            extension_bp=500,
            workers=workers,
        )


def _boundary(scaffold: str, seed_id: str, start: int, end: int, core: tuple[int, int]) -> RefinedBoundary:
    return RefinedBoundary(
        scaffold=scaffold,
        start=start,
        end=end,
        seed_id=seed_id,
        original_start=core[0],
        original_end=core[1],
        confirmed_seeds=(ConfirmedSeed(seed_id, scaffold, *core),),
    )


def _mapping(scaffold: str, seed_id: str, start: int, end: int, length: int) -> SeedGeneMapping:
    return SeedGeneMapping(
        seed_id=seed_id,
        scaffold=scaffold,
        seed_start=start,
        seed_end=end,
        flank_start_bp=0,
        flank_end_bp=length,
    )


def _viral_gene(porf_id: str, scaffold: str, start: int, end: int) -> tuple[pORF, GeneTaxonomy]:
    taxonomy = GeneTaxonomy(
        porf_id=porf_id,
        scaffold=scaffold,
        start=start,
        end=end,
        top1_target="NCLDV__ref|p1",
        top1_prefix="NCLDV__",
        top1_pident=60.0,
        top10_prefixes=["NCLDV__"],
        has_ncldv_mirus=True,
        has_viral=True,
        has_hit=True,
    )
    return pORF(porf_id, scaffold, start, end), taxonomy


def _build_fixture(tmp_path: Path) -> _Fixture:
    """Four parent boundaries on three scaffolds: one partition, a shared pair, one trim.

    Viral-covered genes sit in the flanks that the TIR pairs trim away.
    """
    adjacent, _ = _adjacent_repeat_sequence()
    planted, _ = _planted_repeat()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">adj\n{adjacent}\n>shared\n{planted}\n>trim\n{planted}\n", encoding="utf-8")
    partition_parent = RefinedBoundary(
        scaffold="adj",
        start=0,
        end=len(adjacent),
        seed_id="adj-1",
        original_start=200,
        original_end=1100,
        confirmed_seeds=(
            ConfirmedSeed("adj-1", "adj", 200, 260),
            ConfirmedSeed("adj-2", "adj", 620, 680),
            ConfirmedSeed("adj-3", "adj", 980, 1100),
        ),
        seed_sources=["hhg", "marker_validation"],
    )
    boundaries = [
        _boundary("shared", "shared-1", 180, 245, (180, 245)),
        partition_parent,
        _boundary("shared", "shared-2", 260, 340, (260, 340)),
        _boundary("trim", "trim-1", 50, 470, (150, 350)),
    ]
    allowed_intervals = [(0, 520), (0, len(adjacent)), (0, 520), (50, 470)]
    markers = [
        _Marker("adj", 220, 235, "first|aa1-15", hmm_target="gvogm0003"),
        _Marker("adj", 235, 250, "first|aa16-30", hmm_target="gvogm0004"),
        _Marker("adj", 640, 660, "second|aa1-20", hmm_target="plv_mcp"),
        _Marker("adj", 1000, 1020, "neighbor|aa1-20", hmm_target="mirus_core"),
        _Marker("shared", 220, 230, "marker-1"),
        _Marker("shared", 280, 290, "marker-2"),
        _Marker("trim", 250, 270, "marker-3"),
    ]
    mappings = {
        mapping.seed_id: mapping
        for mapping in (
            _mapping("adj", "adj-1", 200, 260, len(adjacent)),
            _mapping("shared", "shared-1", 180, 245, len(planted)),
            _mapping("shared", "shared-2", 260, 340, len(planted)),
            _mapping("trim", "trim-1", 150, 350, len(planted)),
        )
    }
    genes = [
        _viral_gene("adj_1", "adj", 20, 80),
        _viral_gene("adj_2", "adj", 880, 940),
        _viral_gene("shared_1", "shared", 190, 240),
        _viral_gene("trim_1", "trim", 60, 90),
        _viral_gene("trim_2", "trim", 430, 460),
    ]
    proteome_index = {scaffold: [] for scaffold in ("adj", "shared", "trim")}
    for gene, _ in genes:
        proteome_index[gene.scaffold].append(gene)
    return _Fixture(
        genome_path=genome_path,
        boundaries=boundaries,
        allowed_intervals=allowed_intervals,
        markers=markers,
        query=GenomeDiamondQuery(seed_gene_mappings=mappings),
        taxonomy_map={gene.id: taxonomy for gene, taxonomy in genes},
        proteome_index=proteome_index,
    )


@pytest.fixture
def fixture(tmp_path: Path) -> _Fixture:
    return _build_fixture(tmp_path)


def test_parallel_planning_equals_serial_planning(fixture: _Fixture) -> None:
    serial = fixture.refine(workers=1)
    parallel = fixture.refine(workers=3)

    assert parallel == serial


def test_parallel_planning_preserves_populated_query_role_ids(fixture: _Fixture) -> None:
    query = fixture.query
    query.eve_porf_ids = {"adj-1": ["adj_2"], "shared-1": ["shared_1"]}
    query.boundary_porf_ids = {"adj-1": ["adj_1"], "trim-1": ["trim_1"]}
    query.control_porf_ids = ["trim_2"]
    query.all_porf_ids = ["adj_1", "adj_2", "shared_1", "trim_1", "trim_2"]
    original_roles = pickle.dumps(
        (query.eve_porf_ids, query.boundary_porf_ids, query.control_porf_ids, query.all_porf_ids)
    )

    serial = fixture.refine(workers=1)
    parallel = fixture.refine(workers=3)

    assert parallel == serial
    assert (
        pickle.dumps((query.eve_porf_ids, query.boundary_porf_ids, query.control_porf_ids, query.all_porf_ids))
        == original_roles
    )


def test_parallel_fixture_changes_boundaries_by_tir(fixture: _Fixture) -> None:
    refined = fixture.refine(workers=3)

    assert len(refined) == 6
    assert sum(boundary.tir_boundary_override for boundary in refined) >= 3
    assert {boundary.tir_status for boundary in refined} >= {"detected", "shared_pair"}
    assert [(boundary.start, boundary.end) for boundary in refined if boundary.scaffold == "trim"] == [(100, 420)]


def test_parallel_planning_raises_same_validation_error(fixture: _Fixture) -> None:
    fixture.allowed_intervals[2] = (300, 340)

    with pytest.raises(ValueError, match="allowed interval must contain parent") as serial_error:
        fixture.refine(workers=1)
    with pytest.raises(ValueError, match="allowed interval must contain parent") as parallel_error:
        fixture.refine(workers=3)

    assert str(parallel_error.value) == str(serial_error.value)


def test_parallel_planning_raises_for_absent_scaffold(fixture: _Fixture) -> None:
    fixture.boundaries[3] = _boundary("missing", "missing-1", 50, 470, (150, 350))

    with pytest.raises(ValueError, match="boundary scaffold absent from raw genome: missing"):
        fixture.refine(workers=3)


def test_concurrent_genomes_each_get_their_own_serial_result(tmp_path: Path) -> None:
    (tmp_path / "first").mkdir()
    (tmp_path / "second").mkdir()
    first = _build_fixture(tmp_path / "first")
    second = _build_fixture(tmp_path / "second")
    del second.boundaries[3], second.allowed_intervals[3]
    serial = {"first": first.refine(workers=1), "second": second.refine(workers=1)}
    parallel: dict[str, list[RefinedBoundary]] = {}
    errors: list[BaseException] = []
    barrier = threading.Barrier(2)

    def run(name: str, genome: _Fixture) -> None:
        try:
            barrier.wait()
            parallel[name] = genome.refine(workers=2)
        except BaseException as error:  # noqa: BLE001  # surface any thread failure in the main thread
            errors.append(error)

    threads = [threading.Thread(target=run, args=item) for item in (("first", first), ("second", second))]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    assert not errors
    assert parallel == serial
    assert serial["first"] != serial["second"]


def test_worker_exception_propagates_to_the_caller(fixture: _Fixture) -> None:
    # A marker without a protein ID passes the scaffold and span filters, then fails
    # inside planning. Nothing is patched, so the failure happens in the real worker.
    broken = SimpleNamespace(
        scaffold="shared", validation_status="validated", start=220, end=230, query_porf=None, hmm_target="x"
    )
    fixture.markers.append(broken)

    with pytest.raises(AttributeError):
        fixture.refine(workers=1)
    with pytest.raises(AttributeError):
        fixture.refine(workers=3)


class _SpyPool(ProcessPoolExecutor):
    """Record the start method and tasks, and the process that ran each task."""

    start_methods: ClassVar[list[str]] = []
    tasks: ClassVar[list[object]] = []
    worker_pids: ClassVar[list[int]] = []

    def __init__(self, *args: object, mp_context: object, **kwargs: object) -> None:
        self.start_methods.append(mp_context.get_start_method())
        super().__init__(*args, mp_context=mp_context, **kwargs)

    def map(self, fn: object, *iterables: object, **kwargs: object):
        tasks = list(iterables[0])
        self.tasks.extend(tasks)
        for worker_pid, result in super().map(partial(_run_with_pid, fn), tasks, **kwargs):
            self.worker_pids.append(worker_pid)
            yield result


def _run_with_pid(fn: object, task: object) -> tuple[int, object]:
    return os.getpid(), fn(task)


def test_planning_runs_in_forkserver_child_processes(fixture: _Fixture, monkeypatch: pytest.MonkeyPatch) -> None:
    for attribute in ("start_methods", "tasks", "worker_pids"):
        monkeypatch.setattr(_SpyPool, attribute, [])
    monkeypatch.setattr(terminal_repeats, "ProcessPoolExecutor", _SpyPool)
    serial = fixture.refine(workers=1)

    parallel = fixture.refine(workers=3)

    # The unpatched task function ran in the children: the spy only adds the PID.
    assert parallel == serial
    assert _SpyPool.start_methods == ["forkserver"]
    assert sorted(len(task.boundaries) for task in _SpyPool.tasks) == [1, 1, 2]
    assert len(_SpyPool.worker_pids) == 3
    assert os.getpid() not in _SpyPool.worker_pids


def test_scaffold_tasks_run_largest_first(fixture: _Fixture) -> None:
    tasks = _tasks(fixture)

    assert [(task.scaffold, len(task.boundaries)) for task in tasks] == [("shared", 2), ("adj", 1), ("trim", 1)]


def _tasks(fixture: _Fixture, **overrides: object) -> list[terminal_repeats._ScaffoldPlanningTask]:
    inputs = {
        "validated_markers": fixture.markers,
        "boundary_diamond_query": fixture.query,
        "boundary_taxonomy_map": fixture.taxonomy_map,
        "proteome_index": fixture.proteome_index,
        "extension_bp": 500,
    } | overrides
    return terminal_repeats._build_scaffold_tasks(
        fixture.boundaries, fixture.allowed_intervals, fixture.genome_path, **inputs
    )


def test_scaffold_task_holds_only_its_scaffolds_data(fixture: _Fixture) -> None:
    tasks = {task.scaffold: task for task in _tasks(fixture)}

    shared = tasks["shared"]
    assert [index for index, _, _ in shared.boundaries] == [0, 2]
    assert {marker.query_porf for marker in shared.validated_markers} == {"marker-1", "marker-2"}
    assert set(shared.boundary_diamond_query.seed_gene_mappings) == {"shared-1", "shared-2"}
    assert set(shared.boundary_taxonomy_map) == {"shared_1"}
    assert shared.proteome_index == {"shared": fixture.proteome_index["shared"]}
    assert tasks["adj"].extension_bp == 500


def test_scaffold_task_keeps_taxonomy_read_by_gene_id_or_by_scaffold(fixture: _Fixture) -> None:
    # Planning reads the map by entry scaffold and by the ID of this scaffold's genes.
    # An entry matching either way must survive the filter.
    by_attribute = GeneTaxonomy(porf_id="loose", scaffold="trim", start=1, end=2)
    by_gene_id = GeneTaxonomy(porf_id="trim_1", scaffold="elsewhere", start=60, end=90)
    taxonomy_map = fixture.taxonomy_map | {"loose": by_attribute, "trim_1": by_gene_id}

    tasks = {task.scaffold: task for task in _tasks(fixture, boundary_taxonomy_map=taxonomy_map)}

    assert tasks["trim"].boundary_taxonomy_map == {
        "loose": by_attribute,
        "trim_1": by_gene_id,
        "trim_2": fixture.taxonomy_map["trim_2"],
    }
    assert "trim_1" not in tasks["adj"].boundary_taxonomy_map


def test_scaffold_task_keeps_an_absent_proteome_key_absent(fixture: _Fixture) -> None:
    proteome_index = {scaffold: genes for scaffold, genes in fixture.proteome_index.items() if scaffold != "adj"}

    tasks = {task.scaffold: task for task in _tasks(fixture, proteome_index=proteome_index)}

    assert tasks["adj"].proteome_index == {}
    assert tasks["trim"].proteome_index == {"trim": fixture.proteome_index["trim"]}


def test_scaffold_task_names_the_scaffold_and_genome_but_holds_no_sequence(fixture: _Fixture) -> None:
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(fixture.genome_path, "fasta")}

    tasks = _tasks(fixture)

    assert {task.scaffold for task in tasks} == {"adj", "shared", "trim"}
    assert {task.raw_genome_path for task in tasks} == {str(fixture.genome_path)}
    assert not any(hasattr(task, "sequence") for task in tasks)
    assert all(sequences[task.scaffold].encode() not in pickle.dumps(task) for task in tasks)


def test_each_worker_reads_the_sequence_of_its_own_scaffold(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # Both scaffolds carry identical boundaries and markers; only the sequence differs.
    planted, _ = _planted_repeat()
    plain = "".join(random.Random(5).choice("acgt") for _ in planted)
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">plain\n{plain}\n>repeat\n{planted.lower()}\n>plain2\n{plain}\n", encoding="utf-8")
    scaffolds = ["plain", "repeat", "plain2"]
    fixture = _Fixture(
        genome_path=genome_path,
        boundaries=[_boundary(name, f"{name}-1", 50, 470, (150, 350)) for name in scaffolds],
        allowed_intervals=[(50, 470)] * 3,
        markers=[_Marker(name, 250, 270, f"{name}-marker") for name in scaffolds],
        query=GenomeDiamondQuery(
            seed_gene_mappings={f"{name}-1": _mapping(name, f"{name}-1", 150, 350, len(planted)) for name in scaffolds}
        ),
        taxonomy_map={},
        proteome_index={name: [] for name in scaffolds},
    )
    for attribute in ("start_methods", "tasks", "worker_pids"):
        monkeypatch.setattr(_SpyPool, attribute, [])
    monkeypatch.setattr(terminal_repeats, "ProcessPoolExecutor", _SpyPool)

    refined = {boundary.scaffold: boundary for boundary in fixture.refine(workers=3)}

    assert refined == {boundary.scaffold: boundary for boundary in fixture.refine(workers=1)}
    assert (refined["repeat"].start, refined["repeat"].end, refined["repeat"].tir_status) == (100, 420, "detected")
    assert (refined["plain"].start, refined["plain"].end, refined["plain"].tir_status) == (50, 470, "not_detected")
    assert (refined["plain2"].start, refined["plain2"].end) == (50, 470)
    assert os.getpid() not in _SpyPool.worker_pids


_ADJACENT, _ = _adjacent_repeat_sequence()


def _viral_entry(porf_id: str, scaffold: str, start: int, end: int, *, pident: float) -> GeneTaxonomy:
    """A gene whose only top-10 hit is NCLDV; it counts as viral evidence from 25% identity."""
    return GeneTaxonomy(
        porf_id=porf_id,
        scaffold=scaffold,
        start=start,
        end=end,
        top10_prefixes=["NCLDV__"],
        top10_pidents=[pident],
    )


def _taxonomy_fixture(tmp_path: Path) -> _Fixture:
    """One parent per scaffold on the two-pair sequence; taxonomy decides each outcome.

    cand_*: the left pair holds only a cellular marker, so it needs a viral gene to count.
    res_*: the tail slice holds only a cellular marker, so it needs a viral gene to stay.
    alias: the map key differs from the entry's ``porf_id``.
    dup: two entries share one ``porf_id``; only the second one is viral.
    res_dup: like dup for the tail slice, but the viral entry is the EARLIER duplicate and
        is the gene's own entry, so no lookup by ``porf_id`` can recover it.
    orphan: the map key is not the ID of the scaffold's gene, so that gene has no taxonomy.
    """
    length = len(_ADJACENT)
    names = ["cand_viral", "cand_weak", "res_viral", "res_weak", "alias", "dup", "res_dup", "orphan"]
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text("".join(f">{name}\n{_ADJACENT}\n" for name in names), encoding="utf-8")
    boundaries, markers, mappings, proteome_index, taxonomy_map = [], [], {}, {}, {}
    for name in names:
        spans = [(200, 260), (620, 680)] + ([(980, 1100)] if name.startswith("res") else [])
        seeds = tuple(ConfirmedSeed(f"{name}-{number}", name, *span) for number, span in enumerate(spans))
        boundaries.append(
            RefinedBoundary(
                scaffold=name,
                start=0,
                end=length,
                seed_id=seeds[0].seed_id,
                original_start=200,
                original_end=spans[-1][1],
                confirmed_seeds=seeds,
                seed_sources=["hhg", "marker_validation"],
            )
        )
        mappings[seeds[0].seed_id] = _mapping(name, seeds[0].seed_id, 200, 260, length)
        if name.startswith("res"):
            markers += [
                _Marker(name, 220, 250, "a|aa1-30", hmm_target="gvogm0003"),
                _Marker(name, 640, 660, "b|aa1-20", hmm_target="plv_mcp"),
                _Marker(name, 1000, 1020, "c|aa1-20", hmm_target="VS000087"),
            ]
            gene_start, gene_end = 1000, 1020
        else:
            markers += [
                _Marker(name, 220, 250, "a|aa1-30", hmm_target="VS000087"),
                _Marker(name, 640, 660, "b|aa1-20", hmm_target="plv_mcp"),
            ]
            gene_start, gene_end = 220, 250
        gene_id = f"{name}_gene"
        proteome_index[name] = [pORF(gene_id, name, gene_start, gene_end)]
        pident = 10.0 if name.endswith("weak") else 60.0
        taxonomy_map[gene_id] = _viral_entry(gene_id, name, gene_start, gene_end, pident=pident)
    taxonomy_map["alias_key"] = taxonomy_map.pop("alias_gene")
    proteome_index["alias"] = [pORF("alias_key", "alias", 220, 250)]
    taxonomy_map["dup_gene"] = _viral_entry("dup_shared", "dup", 220, 250, pident=10.0)
    taxonomy_map["dup_second"] = _viral_entry("dup_shared", "dup", 220, 250, pident=60.0)
    taxonomy_map["res_dup_gene"] = _viral_entry("dup_shared", "res_dup", 1000, 1020, pident=60.0)
    taxonomy_map["res_dup_late"] = _viral_entry("dup_shared", "res_dup", 1000, 1020, pident=10.0)
    taxonomy_map["orphan_key"] = taxonomy_map.pop("orphan_gene")
    return _Fixture(
        genome_path=genome_path,
        boundaries=boundaries,
        allowed_intervals=[(0, length)] * len(boundaries),
        markers=markers,
        query=GenomeDiamondQuery(seed_gene_mappings=mappings),
        taxonomy_map=taxonomy_map,
        proteome_index=proteome_index,
    )


def test_viral_taxonomy_decides_candidates_and_residuals_the_same_in_workers(tmp_path: Path) -> None:
    fixture = _taxonomy_fixture(tmp_path)

    serial = fixture.refine(workers=1)
    parallel = fixture.refine(workers=3)

    spans: dict[str, list[tuple[int, int, str]]] = {}
    for boundary in serial:
        spans.setdefault(boundary.scaffold, []).append((boundary.start, boundary.end, boundary.tir_status))
    assert parallel == serial
    assert spans == {
        # A viral gene in the cellular-only pair makes it a candidate; a weak hit does not.
        "cand_viral": [(100, 420, "detected"), (520, 840, "detected")],
        "cand_weak": [(0, 1200, "ambiguous")],
        # A viral gene in the tail keeps that residual slice; a weak hit drops it.
        "res_viral": [(100, 420, "detected"), (520, 840, "detected"), (840, 1200, "not_detected")],
        "res_weak": [(0, 1200, "ambiguous")],
        # Map key unlike porf_id, and two entries with one porf_id, give the viral result.
        "alias": [(100, 420, "detected"), (520, 840, "detected")],
        "dup": [(100, 420, "detected"), (520, 840, "detected")],
        # The viral entry is the earlier duplicate; a weak entry with the same porf_id follows it.
        "res_dup": [(100, 420, "detected"), (520, 840, "detected"), (840, 1200, "not_detected")],
        # A key that is not the gene's ID leaves that gene uncovered, so the split is refused
        # although the entry itself is viral.
        "orphan": [(0, 1200, "ambiguous")],
    }


def test_pool_is_skipped_when_the_socket_directory_is_too_long(
    fixture: _Fixture, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    class _NoPool:
        def __init__(self, *args: object, **kwargs: object) -> None:
            raise AssertionError("no worker pool may start")

    serial = fixture.refine(workers=1)
    monkeypatch.setattr(terminal_repeats, "_MAX_SOCKET_DIR_BYTES", 1)
    monkeypatch.setattr(terminal_repeats, "ProcessPoolExecutor", _NoPool)

    with caplog.at_level("WARNING", logger=terminal_repeats.logger.name):
        parallel = fixture.refine(workers=3)

    assert parallel == serial
    assert [record.getMessage().split(":")[0] for record in caplog.records] == [
        "Terminal-repeat planning runs serially"
    ]


def test_workers_start_under_a_deep_tmpdir(tmp_path: Path) -> None:
    deep_tmpdir = tmp_path / ("a" * 40) / ("b" * 40) / ("c" * 40)
    deep_tmpdir.mkdir(parents=True)
    assert len(os.fsencode(deep_tmpdir)) > 100
    tests_dir = Path(__file__).resolve().parent
    environment = os.environ | {
        "TMPDIR": str(deep_tmpdir),
        "PYTHONPATH": os.pathsep.join(sys.path),
    }

    # A fresh interpreter has no multiprocessing temp dir yet. The forkserver test asserts
    # the start method and that worker PIDs differ from the parent's.
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "-p",
            "no:cacheprovider",
            f"--basetemp={tmp_path / 'inner'}",
            f"{Path(__file__).name}::test_planning_runs_in_forkserver_child_processes",
        ],
        cwd=tests_dir,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "1 passed" in completed.stdout


def test_concurrent_socket_guards_restore_the_original_tempdir(monkeypatch: pytest.MonkeyPatch) -> None:
    deep_tempdir = "/clusterfs/" + "deep/" * 30
    monkeypatch.setattr(terminal_repeats.tempfile, "tempdir", deep_tempdir)
    first_entered = threading.Event()
    second_started = threading.Event()
    second_entered = threading.Event()
    release_first = threading.Event()
    first_done = threading.Event()
    results: dict[str, bool] = {}
    errors: list[BaseException] = []

    def socket_directory() -> str:
        if threading.current_thread().name == "first-guard":
            first_entered.set()
            assert release_first.wait(3), "first guard was not released"
        else:
            second_entered.set()
            assert first_done.wait(3), "first guard did not finish"
        return "/tmp/socket-guard-test"

    def run_guard() -> None:
        name = threading.current_thread().name
        if name == "second-guard":
            second_started.set()
        try:
            results[name] = terminal_repeats._worker_sockets_fit()
        except BaseException as error:  # noqa: BLE001  # report worker failures in the test thread
            errors.append(error)
        finally:
            if name == "first-guard":
                first_done.set()

    monkeypatch.setattr(terminal_repeats.multiprocessing.util, "get_temp_dir", socket_directory)
    first = threading.Thread(target=run_guard, name="first-guard")
    second = threading.Thread(target=run_guard, name="second-guard")
    first.start()
    try:
        assert first_entered.wait(3), "first guard did not enter"
        second.start()
        assert second_started.wait(3), "second guard did not start"
        second_entered.wait(0.5)
    finally:
        release_first.set()
        first.join(timeout=5)
        if second.ident is not None:
            second.join(timeout=5)

    assert not first.is_alive()
    assert not second.is_alive()
    assert not errors
    assert results == {"first-guard": True, "second-guard": True}
    assert terminal_repeats.tempfile.tempdir == deep_tempdir
