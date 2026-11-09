from __future__ import annotations

import errno
import logging
import os
import pickle
import random
import signal
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import dataclass, replace
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio import SeqIO

from virosync.config import PipelineConfig
from virosync.orchestration._flows.single_genome import phase2
from virosync.orchestration._flows.single_genome.phase_state import load_phase2_state, write_phase2_state
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
from virosync.pipeline.phase2 import terminal_repeats
from virosync.pipeline.phase2.boundary_diamond import (
    GeneTaxonomy,
    GenomeDiamondQuery,
    SeedGeneMapping,
    pORF,
)
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.phase2.terminal_repeats import (
    MAX_TIR_ARM_BP,
    find_target_site_duplication,
    find_terminal_inverted_repeats,
    refine_boundaries_with_terminal_repeats,
)

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _run_fixture_pool(tasks, n_boundaries, workers):
    """Exercise worker lifecycle directly, independent of production admission."""
    from virosync.pipeline.host_signatures import get_taxonomy_lookup

    with tempfile.TemporaryDirectory(prefix="virosync-test-pool-") as directory:
        root = Path(directory)
        with (root / "taxonomy.pkl").open("wb") as stream:
            pickle.dump(get_taxonomy_lookup(), stream)
        count = 0
        for count, task in enumerate(tasks, 1):
            with (root / f"task-{count - 1}.pkl").open("wb") as stream:
                pickle.dump(task, stream)
        return terminal_repeats._dispatch_scaffold_files(root, count, n_boundaries, workers)


def _fixture_wait(path):
    """Bound setup separately from the short worker-failure assertion."""
    deadline = time.monotonic() + 45
    while time.monotonic() < deadline:
        if path.exists() and (text := path.read_text()):
            return text
        time.sleep(0.01)
    raise RuntimeError("fixture setup deadline")


def _runtime_fixture(task):
    """Exercise real spawned-worker lifecycle without biological work."""
    kind, path = task
    if kind in {"stall", "error"} and path:
        root = Path(path)
        (root / f"{kind}.pid").write_text(str(os.getpid()))
        _fixture_wait(root / ("stall.pid" if kind == "error" else "error.pid"))
        if kind == "error":
            _fixture_wait(root / "stall.started")
            (root / "error.at").write_text(str(time.monotonic()))
    if kind == "search-kill":
        from virosync.pipeline.phase2.terminal_repeat_search import find_seed_positions

        find_seed_positions("ACGTACGTACGT", 4, ("ACGT",))
    if kind == "error":
        raise RuntimeError("fixture planning failed")
    if kind == "stall":
        if path:
            (Path(path) / "stall.started").write_text("ready")
        time.sleep(10)
        if path:
            (Path(path) / "stall.completed").touch()
        return []
    if kind == "kill":
        from virosync.pipeline.phase2.terminal_repeat_search import _run_genometools

        code = (
            "import os,pathlib,signal,time; pathlib.Path("
            + repr(path)
            + ").write_text(str(os.getpid())); os.kill(os.getppid(),signal.SIGKILL); time.sleep(60)"
        )
        _run_genometools([sys.executable, "-c", code], stdout=subprocess.PIPE)
    from virosync.pipeline.host_signatures import get_taxonomy_lookup

    return [(0, [get_taxonomy_lookup(), os.environ.get("VIROSYNC_HOST_LINEAGE")])]


def test_pool_failure_cancels_stalled_peer(tmp_path, monkeypatch):
    monkeypatch.setattr(terminal_repeats, "_plan_scaffold_task", _runtime_fixture)

    def deadline(signum, frame):
        raise RuntimeError("test supervisor deadline")

    previous = signal.signal(signal.SIGALRM, deadline)
    signal.alarm(45)
    try:
        with pytest.raises(RuntimeError, match="fixture planning failed"):
            _run_fixture_pool([("stall", str(tmp_path)), ("error", str(tmp_path))], 0, 2)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)
    assert time.monotonic() - float((tmp_path / "error.at").read_text()) < 8
    assert not (tmp_path / "stall.completed").exists()
    for kind in ("stall", "error"):
        assert not Path(f"/proc/{int((tmp_path / f'{kind}.pid').read_text())}").exists()


def test_pool_worker_death_kills_genometools_child(tmp_path, monkeypatch):
    monkeypatch.setattr(terminal_repeats, "_plan_scaffold_task", _runtime_fixture)
    child_file = tmp_path / "child.pid"
    with pytest.raises(RuntimeError, match="worker"):
        _run_fixture_pool([("kill", str(child_file))], 0, 1)
    status = Path(f"/proc/{int(child_file.read_text())}/status")

    def child_is_live():
        try:
            return "State:\tZ" not in status.read_text()
        except (FileNotFoundError, ProcessLookupError):
            return False

    deadline = time.monotonic() + 1
    while child_is_live() and time.monotonic() < deadline:
        time.sleep(0.01)
    assert not child_is_live()


def test_pool_receives_explicit_taxonomy_and_host_context(monkeypatch):
    from virosync.pipeline import host_signatures

    monkeypatch.setattr(terminal_repeats, "_plan_scaffold_task", _runtime_fixture)
    monkeypatch.setattr(host_signatures, "_TAXONOMY_LOOKUP", {"fixture": "taxon"})
    monkeypatch.setenv("VIROSYNC_HOST_LINEAGE", "fixture lineage")
    assert _run_fixture_pool([("context", "")], 1, 1) == [[{"fixture": "taxon"}, "fixture lineage"]]


def _pid_plan(task):
    """Record real execution without running repeat searches."""
    pids = Path(os.environ["VIROSYNC_TEST_PIDS"])
    (pids / str(os.getpid())).touch()
    deadline = time.monotonic() + 8
    while len(list(pids.iterdir())) < 2:
        if time.monotonic() >= deadline:
            raise RuntimeError("fixture did not receive two distinct worker PIDs")
        time.sleep(0.02)
    return [(index, [boundary]) for index, boundary, _ in task.boundaries]


def _skewed_planning_worker(connection, *args):
    marker = Path(os.environ["VIROSYNC_TEST_PIDS"]).parent / "first-worker"
    try:
        descriptor = os.open(marker, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        time.sleep(3)
    else:
        os.close(descriptor)
    terminal_repeats._planning_worker(connection, *args)


@pytest.mark.parametrize("startup_skew", [False, True])
def test_public_pool_uses_two_processes_without_parent_record_cache(tmp_path, monkeypatch, startup_skew):
    fixture = _taxonomy_fixture(tmp_path)
    pids = tmp_path / "pids"
    pids.mkdir()
    monkeypatch.setenv("VIROSYNC_TEST_PIDS", str(pids))
    monkeypatch.setattr(terminal_repeats, "_plan_scaffold_task", _pid_plan)
    if startup_skew:
        monkeypatch.setattr(terminal_repeats, "_planning_worker", _skewed_planning_worker)
    monkeypatch.setattr(terminal_repeats, "_available_pool_resources", lambda: (2, 4 * 1024**3))
    original_index = SeqIO.index

    class RawOnlyIndex:
        def __init__(self, *args):
            self.index = original_index(*args)

        def keys(self):
            return self.index.keys()

        def get_raw(self, key):
            return self.index.get_raw(key)

        def close(self):
            self.index.close()

        def __getitem__(self, key):
            raise AssertionError("parent materialized a cached BioPython record")

        def __contains__(self, key):
            raise AssertionError("membership materialized a BioPython record")

    monkeypatch.setattr(SeqIO, "index", RawOnlyIndex)
    fixture.refine(workers=2)
    assert len(list(pids.iterdir())) == 2
    assert not (pids / str(os.getpid())).exists()


def test_public_admission_rejects_pair_heavy_context_before_spawn(tmp_path, monkeypatch):
    length = 512020
    names = [f"pair-heavy-{index}" for index in range(4)]
    genome = tmp_path / "raw.fna"
    genome.write_text("".join(f">{name}\n" + "A" * length + "\n" for name in names))
    boundaries = [_boundary(name, name, 0, length, (256000, 256100)) for name in names]
    markers = [_Marker(name, 256000, 256100, f"{name}-mcp") for name in names]
    mappings = {name: _mapping(name, name, 256000, 256100, length) for name in names}
    fixture = _Fixture(
        genome, boundaries, [(0, length)] * 4, markers, GenomeDiamondQuery(seed_gene_mappings=mappings), {}, {}
    )
    admitted = []

    def dispatch(root, count, n_boundaries, workers):
        admitted.append(workers)
        return [[boundary] for boundary in boundaries]

    monkeypatch.setattr(terminal_repeats, "_available_pool_resources", lambda: (4, 2 * 1024**3))
    monkeypatch.setattr(terminal_repeats, "_dispatch_scaffold_files", dispatch)
    monkeypatch.setattr(terminal_repeats, "_plan_boundary_partition", lambda boundary, **kwargs: [boundary])
    assert fixture.refine(workers=4) == boundaries
    assert admitted == [], "pair-heavy 512 kb contexts must select serial before worker creation"


def test_marker_pair_union_exceeds_single_rectangle_cap_but_fits_metadata_bound():
    seed = "ACGTGCA"
    reverse = _reverse_complement(seed)
    length = 1000
    index = {seed: tuple(range(80)), reverse: tuple(length - 7 - i for i in reversed(range(80)))}
    limits = [(70, 70), (86, 56)]
    pairs, allowed = terminal_repeats._eligible_region_seed_pairs(index, length, limits)
    assert len(pairs) == 4896 > terminal_repeats._MAX_SEED_PAIRINGS
    assert allowed[:, 0].all()
    marker_bounds = [min((left - 6) * (right - 6), 4**7 * 4096, 32 * (left + right - 12)) for left, right in limits]
    assert sum(marker_bounds) == 8096
    bound, intervals, marker_pairs, chunks = terminal_repeats._seed_pair_cardinality_bounds(length, limits)
    assert bound == intervals == 8096
    assert marker_pairs == 4096
    assert len(pairs) <= bound
    assert terminal_repeats._pair_allocation_bytes(bound, intervals, marker_pairs, chunks, 2) == (
        max(96 * bound, 32 * bound + 160 * marker_pairs + 24 * intervals, 24 * bound + 808 * intervals)
        + 256 * chunks
        + 2 * 4**7
    )


def _error_storage_failure_plan(task):
    original_write = Path.write_text

    def fail_error_file(path, *args, **kwargs):
        if path.name.startswith("error-"):
            raise OSError(errno.ENOSPC, "fixture error storage exhausted")
        return original_write(path, *args, **kwargs)

    Path.write_text = fail_error_file
    raise ValueError("original distinctive planning error " + "🙂" * 10000)


def test_worker_error_survives_unwritable_error_storage(monkeypatch):
    monkeypatch.setattr(terminal_repeats, "_plan_scaffold_task", _error_storage_failure_plan)
    with pytest.raises(RuntimeError, match="original distinctive planning error") as error:
        _run_fixture_pool([("error", "")], 0, 1)
    assert len(pickle.dumps(("error", str(error.value)))) <= 1024


@pytest.mark.parametrize("resources", [(1, 4 * 1024**3), (8, 128 * 1024**2)])
def test_resource_limited_public_pool_uses_in_process_reference(tmp_path, monkeypatch, resources):
    fixture = _taxonomy_fixture(tmp_path)
    monkeypatch.setattr(terminal_repeats, "_available_pool_resources", lambda: resources)
    processes = []

    def reference(boundary, **kwargs):
        processes.append(os.getpid())
        return [boundary]

    monkeypatch.setattr(terminal_repeats, "_plan_boundary_partition", reference)
    monkeypatch.setattr(
        terminal_repeats, "_dispatch_scaffold_files", lambda *args: pytest.fail("spawned under serial cap")
    )
    fixture.refine(workers=8)
    assert processes and set(processes) == {os.getpid()}


def test_worker_sequence_file_is_exact_and_never_indexed(tmp_path, monkeypatch):
    fixture = _taxonomy_fixture(tmp_path)
    task = next(
        iter(
            terminal_repeats._build_scaffold_tasks(
                fixture.boundaries,
                fixture.allowed_intervals,
                fixture.genome_path,
                validated_markers=fixture.markers,
                boundary_diamond_query=fixture.query,
                boundary_taxonomy_map=fixture.taxonomy_map,
                proteome_index=fixture.proteome_index,
                extension_bp=500,
                metadata_families={},
            )
        )
    )
    wrapped = tmp_path / "wrapped.fna"
    wrapped.write_text(f">{task.scaffold}\nacgT\nNn\n")
    task = replace(task, raw_genome_path=str(wrapped))
    spool = tmp_path / "spool"
    spool.mkdir()
    terminal_repeats._prepare_worker_files([task], spool)
    terminal_repeats._write_worker_sequences(spool, 1)
    with (spool / "task-0.pkl").open("rb") as stream:
        task = pickle.load(stream)
    assert Path(task.sequence_path).read_bytes() == b"ACGTNN"
    sequence = " ACGTN\n"
    path = tmp_path / "sequence.txt"
    path.write_bytes(sequence.encode("ascii"))
    task = replace(task, sequence_path=str(path))
    monkeypatch.setattr(SeqIO, "index", lambda *args: pytest.fail("worker indexed FASTA"))
    seen = []

    def capture(boundary, **kwargs):
        seen.append(kwargs["sequence"])
        return [boundary]

    monkeypatch.setattr(terminal_repeats, "_plan_boundary_partition", capture)
    terminal_repeats._plan_scaffold_task(task)
    assert seen and all(value == sequence for value in seen)


def test_forced_worker_failure_removes_actual_repfind_files(tmp_path, monkeypatch):
    receipt = tmp_path / "workdir"
    executable = tmp_path / "gt"
    executable.write_text(
        "#!" + sys.executable + "\nimport os,pathlib,signal,sys,time\n"
        "work=pathlib.Path(sys.argv[sys.argv.index('-db')+1]).parent\n"
        "pathlib.Path(" + repr(str(receipt)) + ").write_text(str(work))\n"
        "os.kill(os.getppid(),signal.SIGKILL)\ntime.sleep(60)\n"
    )
    executable.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path) + os.pathsep + os.environ["PATH"])
    monkeypatch.setattr(terminal_repeats, "_plan_scaffold_task", _runtime_fixture)
    with pytest.raises(RuntimeError, match="worker"):
        _run_fixture_pool([("search-kill", "")], 0, 1)
    work = Path(receipt.read_text())
    assert not work.exists()
    assert not work.parent.exists()


def test_resource_limits_include_finite_cgroup_ancestors(tmp_path, monkeypatch):
    parent = tmp_path / "parent"
    leaf = parent / "leaf"
    leaf.mkdir(parents=True)
    (leaf / "cpu.max").write_text("max 100000")
    (leaf / "memory.max").write_text("max")
    (parent / "cpu.max").write_text("250000 100000")
    (parent / "memory.max").write_text("1000000")
    (parent / "memory.current").write_text("400000")
    monkeypatch.setattr(terminal_repeats, "_cgroup_directories", lambda: [leaf, parent])
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: set(range(8)))
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "4")
    assert terminal_repeats._available_pool_resources() == (2, 600000)
    (parent / "cpu.max").write_text("50000 100000")
    assert terminal_repeats._available_pool_resources()[0] == 1


def test_cgroup_reclaims_only_clean_inactive_file_cache(tmp_path, monkeypatch):
    (tmp_path / "memory.max").write_text("1000000")
    (tmp_path / "memory.current").write_text("900000")
    (tmp_path / "memory.stat").write_text(
        "inactive_file 700000\nfile_dirty 100000\nfile_writeback 50000\nfile 850000\nshmem 150000\n"
    )
    monkeypatch.setattr(terminal_repeats, "_cgroup_directories", lambda: [tmp_path])
    monkeypatch.setattr(terminal_repeats.psutil, "virtual_memory", lambda: SimpleNamespace(available=2000000))
    assert terminal_repeats._available_pool_resources()[1] == 650000
    (tmp_path / "memory.stat").write_text("inactive_file 0\nfile 850000\nshmem 850000\n")
    assert terminal_repeats._available_pool_resources()[1] == 100000


def test_known_finite_cgroup_without_usage_does_not_use_host_budget(tmp_path, monkeypatch):
    for key in ("SLURM_JOB_ID", "SLURM_JOBID", "SLURM_CPUS_PER_TASK"):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setattr(terminal_repeats.psutil, "virtual_memory", lambda: SimpleNamespace(available=16 * 1024**3))
    monkeypatch.setattr(terminal_repeats, "_cgroup_directories", lambda: [tmp_path])
    (tmp_path / "memory.max").write_text(str(128 * 1024**2))
    assert terminal_repeats._available_pool_resources()[1] is None
    monkeypatch.setattr(terminal_repeats, "_cgroup_directories", lambda: [])
    assert terminal_repeats._available_pool_resources()[1] == 16 * 1024**3


def test_unknown_slurm_memory_budget_uses_logged_serial(tmp_path, monkeypatch, caplog):
    fixture = _taxonomy_fixture(tmp_path)
    monkeypatch.setenv("SLURM_JOB_ID", "fixture-job")
    monkeypatch.setattr(terminal_repeats, "_cgroup_directories", lambda: [])
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {2, 3})
    monkeypatch.setattr(
        terminal_repeats, "_prepare_worker_files", lambda *args: pytest.fail("spooled without Slurm budget")
    )
    monkeypatch.setattr(terminal_repeats, "_plan_boundary_partition", lambda boundary, **kwargs: [boundary])
    with caplog.at_level(logging.INFO, logger=terminal_repeats.__name__):
        fixture.refine(workers=2)
    assert "Slurm" in caplog.text and "serial" in caplog.text


@pytest.mark.parametrize("per_cpu", [False, True])
def test_slurm_tighter_job_budget_rejects_looser_v2_ancestor(tmp_path, monkeypatch, per_cpu):
    monkeypatch.setenv("SLURM_JOB_ID", "fixture-job")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.delenv("SLURM_MEM_PER_NODE", raising=False)
    monkeypatch.delenv("SLURM_MEM_PER_CPU", raising=False)
    monkeypatch.setenv("SLURM_MEM_PER_CPU" if per_cpu else "SLURM_MEM_PER_NODE", "8" if per_cpu else "16")
    monkeypatch.setattr(terminal_repeats, "_cgroup_directories", lambda: [tmp_path])
    (tmp_path / "memory.max").write_text(str(64 * 1024**2))
    (tmp_path / "memory.current").write_text("0")
    assert terminal_repeats._available_pool_resources()[1] is None
    (tmp_path / "memory.max").write_text(str(16 * 1024**2))
    assert terminal_repeats._available_pool_resources()[1] == 16 * 1024**2


@pytest.mark.parametrize("shortage", ["memory", "disk"])
def test_spool_preflight_precedes_metadata_or_sequence_reads(tmp_path, monkeypatch, shortage, caplog):
    fixture = _taxonomy_fixture(tmp_path)
    monkeypatch.setattr(
        terminal_repeats, "_available_pool_resources", lambda: (2, 1 if shortage == "memory" else 4 * 1024**3)
    )
    monkeypatch.setattr(
        terminal_repeats, "_prepare_worker_files", lambda *args: pytest.fail("spooled before preflight")
    )
    monkeypatch.setattr(terminal_repeats, "_plan_boundary_partition", lambda boundary, **kwargs: [boundary])
    import shutil

    monkeypatch.setattr(
        shutil, "disk_usage", lambda path: SimpleNamespace(free=0 if shortage == "disk" else 4 * 1024**3)
    )
    with caplog.at_level(logging.INFO, logger=terminal_repeats.__name__):
        fixture.refine(workers=2)
    assert "serial" in caplog.text


@pytest.mark.parametrize("error_number", [errno.ENOSPC, errno.EDQUOT])
@pytest.mark.parametrize("stage", ["metadata", "sequences"])
def test_production_spool_storage_failure_cleans_up_before_serial(tmp_path, monkeypatch, error_number, stage, caplog):
    fixture = _taxonomy_fixture(tmp_path)
    monkeypatch.setattr(terminal_repeats, "_available_pool_resources", lambda: (2, 4 * 1024**3))
    roots = []

    def exhaust(tasks, root):
        roots.append(root)
        (root / "partial-sequence").write_bytes(b"ACGT")
        raise OSError(error_number, "fixture storage exhausted")

    def serial(boundary, **kwargs):
        assert roots and not roots[0].exists()
        return [boundary]

    if stage == "metadata":
        monkeypatch.setattr(terminal_repeats, "_prepare_worker_files", exhaust)
    else:
        monkeypatch.setattr(terminal_repeats, "_write_worker_sequences", lambda root, count: exhaust(None, root))
    monkeypatch.setattr(terminal_repeats, "_plan_boundary_partition", serial)
    with caplog.at_level(logging.INFO, logger=terminal_repeats.__name__):
        fixture.refine(workers=2)
    assert "serial" in caplog.text and "storage" in caplog.text


def test_native_lifetime_wrapper_ignores_python_path_shadowing(tmp_path, monkeypatch):
    from virosync.pipeline.phase2.terminal_repeat_search import _run_genometools

    (tmp_path / "ctypes.py").write_text("raise RuntimeError('shadowed ctypes')\n")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("PYTHONPATH", str(tmp_path))
    _run_genometools(["/bin/true"], stdout=subprocess.PIPE)


def test_pool_memory_cap_accounts_for_context_and_files(monkeypatch):
    monkeypatch.setattr(terminal_repeats, "_available_pool_resources", lambda: (8, 2 * 1024**3))
    limit = terminal_repeats._pool_worker_limit
    small = limit(
        8, sequence_bytes=1000, scan_bytes=1000, task_bytes=100, context_bytes=100, file_bytes=1000, pair_bytes=0
    )
    context = limit(
        8,
        sequence_bytes=1000,
        scan_bytes=1000,
        task_bytes=100,
        context_bytes=256 * 1024**2,
        file_bytes=1000,
        pair_bytes=0,
    )
    files = limit(
        8, sequence_bytes=1000, scan_bytes=1000, task_bytes=100, context_bytes=100, file_bytes=1024**3, pair_bytes=0
    )
    assert small > context == 1
    assert small > files


def test_pool_cpu_limit_is_independent_of_memory(monkeypatch):
    monkeypatch.setattr(terminal_repeats, "_available_pool_resources", lambda: (1, 128 * 1024**3))
    assert (
        terminal_repeats._pool_worker_limit(
            16, sequence_bytes=1000, scan_bytes=1000, task_bytes=100, context_bytes=100, file_bytes=1000, pair_bytes=0
        )
        == 1
    )


def _nonreceiving_worker(connection, *args):
    """Advertise readiness, then leave the parent's task unread."""
    os.setsid()
    marker = Path(os.environ["VIROSYNC_TEST_FIRST_WORKER"])
    try:
        descriptor = os.open(marker, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        marker.with_suffix(".peer.pid").write_text(str(os.getpid()))
        _fixture_wait(marker.with_suffix(".ready"))
        time.sleep(0.5)
        marker.with_suffix(".error.at").write_text(str(time.monotonic()))
        connection.send(("error", "fixture peer error"))
        time.sleep(10)
    else:
        os.write(descriptor, str(os.getpid()).encode())
        os.close(descriptor)
        _fixture_wait(marker.with_suffix(".peer.pid"))
        marker.with_suffix(".ready").write_text(str(time.monotonic()))
        connection.send(("ready", None))
        time.sleep(10)
        marker.with_suffix(".completed").touch()


def test_nonreceiving_worker_does_not_hide_peer_error(tmp_path, monkeypatch):
    monkeypatch.setenv("VIROSYNC_TEST_FIRST_WORKER", str(tmp_path / "first"))
    monkeypatch.setattr(terminal_repeats, "_planning_worker", _nonreceiving_worker)

    def deadline(signum, frame):
        raise RuntimeError("test supervisor deadline")

    previous = signal.signal(signal.SIGALRM, deadline)
    signal.alarm(45)
    try:
        with pytest.raises(RuntimeError, match="fixture peer error"):
            _run_fixture_pool([("x" * 1024**2, ""), ("other", "")], 0, 2)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)
    marker = tmp_path / "first"
    assert time.monotonic() - float(marker.with_suffix(".ready").read_text()) < 5
    assert time.monotonic() - float(marker.with_suffix(".error.at").read_text()) < 5
    assert not marker.with_suffix(".completed").exists()
    for path in (marker, marker.with_suffix(".peer.pid")):
        assert not Path(f"/proc/{int(path.read_text())}").exists()


def _parent_death_plan(task):
    from virosync.pipeline.phase2.terminal_repeat_search import _run_genometools

    code = (
        "import os,pathlib,time; pathlib.Path("
        + repr(task[1])
        + ").write_text(str(os.getpid())+' '+str(os.getppid())); time.sleep(60)"
    )
    ready = Path(task[1]).with_suffix(".ready")
    pending = ready.with_suffix(".pending")
    pending.write_text(f"{os.getpid()} {os.getppid()}")
    pending.replace(ready)
    _run_genometools([sys.executable, "-c", code], stdout=subprocess.PIPE)


def _parent_death_probe(path, serial):
    if serial:
        _parent_death_plan(("", path))
    else:
        terminal_repeats._plan_scaffold_task = _parent_death_plan
        _run_fixture_pool([("", path)], 0, 1)


@pytest.mark.parametrize("serial", [False, True])
def test_parent_death_terminates_native_child(tmp_path, serial):
    receipt = tmp_path / "native.pid"
    ready = receipt.with_suffix(".ready")
    code = (
        "from test_terminal_repeats import _parent_death_probe; _parent_death_probe("
        + repr(str(receipt))
        + ", "
        + repr(serial)
        + ")"
    )
    environment = dict(os.environ, PYTHONPATH=os.pathsep.join(["tests", "src"]))
    parent = subprocess.Popen([sys.executable, "-c", code], env=environment)
    child_ids = []
    try:
        # Imports/taxonomy/spawn took ~9s; bound setup apart from native lifetime checks.
        deadline = time.monotonic() + 45
        while not ready.exists() and parent.poll() is None and time.monotonic() < deadline:
            time.sleep(0.02)
        assert parent.poll() is None
        assert ready.exists()
        launcher_pid, launcher_parent = map(int, ready.read_text().split())
        if serial:
            assert launcher_pid == parent.pid
            assert launcher_parent == os.getpid()
        else:
            assert launcher_pid != parent.pid
            assert launcher_parent == parent.pid
            child_ids = [launcher_pid]
        deadline = time.monotonic() + 8
        while not receipt.exists() and time.monotonic() < deadline:
            time.sleep(0.02)
        assert receipt.exists()
        child_ids = [int(value) for value in receipt.read_text().split()]
        assert (child_ids[1] == parent.pid) is serial
        assert child_ids[1] == launcher_pid
        parent.kill()
        parent.wait(timeout=2)

        def live(pid):
            try:
                return "State:\tZ" not in Path(f"/proc/{pid}/status").read_text()
            except (FileNotFoundError, ProcessLookupError):
                return False

        deadline = time.monotonic() + 2
        while any(live(pid) for pid in child_ids) and time.monotonic() < deadline:
            time.sleep(0.02)
        assert not any(live(pid) for pid in child_ids)
    finally:
        if parent.poll() is None:
            parent.kill()
            parent.wait()
        for pid in child_ids:
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass


@dataclass(frozen=True, slots=True)
class _Marker:
    scaffold: str
    start: int
    end: int
    query_porf: str
    validation_status: str = "validated"
    hmm_target: str = "ncldv_mcp"
    hmm_score: float = 100.0


def _random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _at_rich_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choices("ACGT", weights=(35, 15, 15, 35), k=length))


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def _planted_repeat(
    *,
    arm_length: int = 60,
    mutated_positions: tuple[int, ...] = (),
    seed: int = 17,
) -> tuple[str, tuple[int, int, int, int]]:
    rng = random.Random(seed)
    left_flank = _random_dna(99, rng) + "A"
    arm = _random_dna(arm_length, rng)
    paired_arm = list(arm)
    substitutions = {"A": "C", "C": "G", "G": "T", "T": "A"}
    for position in mutated_positions:
        paired_arm[position] = substitutions[paired_arm[position]]
    interior = "A" + _random_dna(198, rng) + "A"
    right_flank = "A" + _random_dna(99, rng)
    sequence = left_flank + arm + interior + _reverse_complement("".join(paired_arm)) + right_flank
    left_start = len(left_flank)
    left_end = left_start + arm_length
    right_start = left_end + len(interior)
    right_end = right_start + arm_length
    return sequence, (left_start, left_end, right_start, right_end)


def _adjacent_repeat_sequence(
    *,
    identical_arms: bool = False,
) -> tuple[str, tuple[tuple[int, int], tuple[int, int]]]:
    rng = random.Random(221)
    prefix = _random_dna(90, rng) + "A" * 10
    first_arm = _random_dna(60, rng)
    second_arm = first_arm if identical_arms else _random_dna(60, rng)
    first_interior = "A" * 10 + _random_dna(180, rng) + "A" * 10
    between = "A" * 10 + _random_dna(80, rng) + "A" * 10
    second_interior = "A" * 10 + _random_dna(180, rng) + "A" * 10
    tail = "A" * 10 + _random_dna(350, rng)
    sequence = (
        prefix
        + first_arm
        + first_interior
        + _reverse_complement(first_arm)
        + between
        + second_arm
        + second_interior
        + _reverse_complement(second_arm)
        + tail
    )
    return sequence, ((100, 420), (520, 840))


@pytest.mark.parametrize(
    ("mutated_positions", "expected_identity"),
    [
        ((), 1.0),
        ((10, 20, 30, 40, 50, 55), 0.9),
    ],
)
def test_scanner_recovers_exact_planted_coordinates(
    mutated_positions: tuple[int, ...],
    expected_identity: float,
) -> None:
    sequence, expected = _planted_repeat(mutated_positions=mutated_positions)

    candidates = find_terminal_inverted_repeats(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        marker_start=250,
        marker_end=270,
    ).candidates

    assert len(candidates) == 1
    candidate = candidates[0]
    assert (
        candidate.left_start,
        candidate.left_end,
        candidate.right_start,
        candidate.right_end,
    ) == expected
    assert candidate.identity == pytest.approx(expected_identity)


@pytest.mark.parametrize(("arm_length", "detected"), [(49, False), (50, True)])
def test_scanner_enforces_minimum_arm_length(arm_length: int, detected: bool) -> None:
    rng = random.Random(730 + arm_length)
    arm = _random_dna(arm_length, rng)
    prefix = _random_dna(90, rng) + "A" * 10
    interior = "A" * 10 + _random_dna(180, rng) + "G" * 10
    suffix = "G" * 10 + _random_dna(90, rng)
    sequence = prefix + arm + interior + _reverse_complement(arm) + suffix
    expected = (len(prefix), len(prefix) + arm_length + len(interior) + arm_length)
    marker_start = len(prefix) + arm_length + 90

    candidates = find_terminal_inverted_repeats(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        marker_start=marker_start,
        marker_end=marker_start + 20,
    ).candidates

    assert bool(candidates) is detected
    if detected:
        assert (candidates[0].left_start, candidates[0].right_end) == expected


@pytest.mark.parametrize(
    "mutated_positions",
    [(), (10,), (10, 20, 30)],
)
def test_long_repeat_cap_preserves_outer_boundary_coordinates(
    mutated_positions: tuple[int, ...],
) -> None:
    rng = random.Random(8)
    left_flank = _random_dna(99, rng) + "A"
    right_flank = "A" + _random_dna(99, rng)
    arm = _random_dna(MAX_TIR_ARM_BP + 115, rng)
    paired_arm = list(arm)
    substitutions = {"A": "C", "C": "G", "G": "T", "T": "A"}
    for position in mutated_positions:
        paired_arm[position] = substitutions[paired_arm[position]]
    interior = _random_dna(200, rng)
    sequence = left_flank + arm + interior + _reverse_complement("".join(paired_arm)) + right_flank
    marker_start = len(left_flank) + len(arm) + 80

    candidates = find_terminal_inverted_repeats(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        marker_start=marker_start,
        marker_end=marker_start + 20,
    ).candidates

    assert len(candidates) == 1
    assert candidates[0].left_start == len(left_flank)
    assert candidates[0].right_end == len(sequence) - len(right_flank)
    assert candidates[0].arm_length == MAX_TIR_ARM_BP
    assert candidates[0].alignment_length == MAX_TIR_ARM_BP + 115
    assert candidates[0].alignment_capped is True
    assert candidates[0].identity == pytest.approx((MAX_TIR_ARM_BP - len(mutated_positions)) / MAX_TIR_ARM_BP)


def test_capped_repeat_rejects_subthreshold_reported_arm() -> None:
    rng = random.Random(8)
    left_flank = _random_dna(99, rng) + "A"
    right_flank = "A" + _random_dna(99, rng)
    arm = _random_dna(MAX_TIR_ARM_BP + 115, rng)
    paired_arm = list(arm)
    substitutions = {"A": "C", "C": "G", "G": "T", "T": "A"}
    for position in range(5, MAX_TIR_ARM_BP, 9):
        paired_arm[position] = substitutions[paired_arm[position]]
    interior = _random_dna(200, rng)
    sequence = left_flank + arm + interior + _reverse_complement("".join(paired_arm)) + right_flank
    marker_start = len(left_flank) + len(arm) + 80

    assert (
        find_terminal_inverted_repeats(
            sequence,
            scan_start=0,
            scan_end=len(sequence),
            marker_start=marker_start,
            marker_end=marker_start + 20,
        ).candidates
        == []
    )


def test_target_site_duplication_is_descriptive_only() -> None:
    sequence, expected = _planted_repeat()
    bases = list(sequence)
    bases[expected[0] - 6 : expected[0]] = "GAGGCT"
    bases[expected[3] : expected[3] + 6] = "GAGGCT"
    sequence = "".join(bases)
    candidate = find_terminal_inverted_repeats(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        marker_start=250,
        marker_end=270,
    ).candidates[0]

    assert find_target_site_duplication(sequence, candidate) == "GAGGCT"


def test_mavirus_terminal_repeats_recover_published_outer_boundary() -> None:
    """Regress the KU052222 integration site from DOI 10.1038/nature20593."""
    record = SeqIO.read(Path(__file__).parent / "fixtures" / "KU052222.fna", "fasta")
    sequence = str(record.seq)

    candidates = find_terminal_inverted_repeats(
        sequence,
        scan_start=5000,
        scan_end=35190,
        marker_start=15000,
        marker_end=25000,
    ).candidates

    assert len(candidates) == 1
    candidate = candidates[0]
    assert (candidate.left_start, candidate.right_end) == (10000, 30190)
    assert candidate.arm_length == MAX_TIR_ARM_BP
    assert candidate.alignment_length == 600
    assert candidate.alignment_capped is True
    assert find_target_site_duplication(sequence, candidate) == "GAGGCT"


def test_scanner_rejects_random_and_short_period_sequences() -> None:
    rng = random.Random(91)
    for _ in range(100):
        sequence = _random_dna(800, rng)
        assert (
            find_terminal_inverted_repeats(
                sequence,
                scan_start=0,
                scan_end=len(sequence),
                marker_start=390,
                marker_end=410,
            ).candidates
            == []
        )

    low_complexity_arm = "ACGT" * 15
    sequence = (
        _random_dna(100, rng)
        + low_complexity_arm
        + _random_dna(200, rng)
        + _reverse_complement(low_complexity_arm)
        + _random_dna(100, rng)
    )
    assert (
        find_terminal_inverted_repeats(
            sequence,
            scan_start=0,
            scan_end=len(sequence),
            marker_start=250,
            marker_end=270,
        ).candidates
        == []
    )


def test_scanner_rejects_at_rich_composition_matched_arms() -> None:
    rng = random.Random(108)
    left_arm = _at_rich_dna(60, rng)
    composition_matched_arm = list(left_arm)
    rng.shuffle(composition_matched_arm)
    sequence = (
        _at_rich_dna(100, rng)
        + left_arm
        + _at_rich_dna(200, rng)
        + _reverse_complement("".join(composition_matched_arm))
        + _at_rich_dna(100, rng)
    )

    assert (
        find_terminal_inverted_repeats(
            sequence,
            scan_start=0,
            scan_end=len(sequence),
            marker_start=250,
            marker_end=270,
        ).candidates
        == []
    )


def test_marker_anchors_deduplicate_same_outer_pair_by_longest_evidence() -> None:
    sequence, expected = _planted_repeat(arm_length=MAX_TIR_ARM_BP + 115, seed=8)
    groups = [
        terminal_repeats._MarkerGroup(
            protein_id="arm-marker",
            start=250,
            end=270,
            classification_markers=tuple((name, "") for name in frozenset({"marker-a"})),
        ),
        terminal_repeats._MarkerGroup(
            protein_id="core-marker",
            start=795,
            end=815,
            classification_markers=tuple((name, "") for name in frozenset({"marker-b"})),
        ),
    ]

    candidates = terminal_repeats._discover_marker_anchored_candidates(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        groups=groups,
    )

    assert len(candidates) == 1
    candidate = candidates[0]
    assert (candidate.left_start, candidate.right_end) == (expected[0], expected[3])
    assert candidate.arm_length == MAX_TIR_ARM_BP
    assert candidate.alignment_length == MAX_TIR_ARM_BP + 115


def test_shared_scan_index_preserves_candidates() -> None:
    sequence, _ = _planted_repeat(arm_length=MAX_TIR_ARM_BP + 115, seed=8)
    scan_start = 40
    scan_end = len(sequence) - 40
    arguments = {
        "scan_start": scan_start,
        "scan_end": scan_end,
        "marker_start": 795,
        "marker_end": 815,
    }

    uncached = find_terminal_inverted_repeats(sequence, **arguments).candidates
    cached = find_terminal_inverted_repeats(
        sequence,
        **arguments,
        _scan_seed_index=terminal_repeats._build_scan_seed_index(sequence[scan_start:scan_end]),
    ).candidates

    assert cached == uncached


def test_duplicate_repeat_pairs_are_ambiguous() -> None:
    rng = random.Random(32)
    arm = _random_dna(60, rng)
    sequence = (
        _random_dna(99, rng)
        + "A"
        + arm
        + "A"
        + _random_dna(198, rng)
        + "A"
        + _reverse_complement(arm)
        + "A"
        + _random_dna(28, rng)
        + "A"
        + _reverse_complement(arm)
        + "A"
        + _random_dna(99, rng)
    )

    candidates = find_terminal_inverted_repeats(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        marker_start=250,
        marker_end=270,
    ).candidates

    assert len(candidates) == 2
    assert {(candidate.left_start, candidate.right_end) for candidate in candidates} == {
        (100, 420),
        (100, 510),
    }


def _query(sequence_length: int) -> GenomeDiamondQuery:
    mapping = SeedGeneMapping(
        seed_id="seed-1",
        scaffold="contig",
        seed_start=150,
        seed_end=350,
        flank_start_bp=0,
        flank_end_bp=sequence_length,
    )
    return GenomeDiamondQuery(seed_gene_mappings={mapping.seed_id: mapping})


@pytest.mark.parametrize(
    ("boundary_start", "boundary_end", "expected_bounds", "expected_status", "expected_override", "core_bounds"),
    [
        (120, 400, (120, 400), "ambiguous", False, (150, 350)),
        (50, 470, (100, 420), "detected", True, (150, 350)),
        (50, 470, (50, 470), "ambiguous", False, (50, 470)),
    ],
)
def test_unique_repeat_changes_only_reachable_noncore_flanks(
    tmp_path: Path,
    boundary_start: int,
    boundary_end: int,
    expected_bounds: tuple[int, int],
    expected_status: str,
    expected_override: bool,
    core_bounds: tuple[int, int],
) -> None:
    sequence, _ = _planted_repeat()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    boundary = RefinedBoundary(
        scaffold="contig",
        start=boundary_start,
        end=boundary_end,
        seed_id="seed-1",
        original_start=core_bounds[0],
        original_end=core_bounds[1],
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", *core_bounds),),
    )

    refined = refine_boundaries_with_terminal_repeats(
        [boundary],
        allowed_intervals=[(boundary.start, boundary.end)],
        raw_genome_path=genome_path,
        validated_markers=[_Marker("contig", 250, 270, "marker-1")],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
        metadata_families={},
    )[0]

    assert (refined.start, refined.end) == expected_bounds
    assert (refined.pre_tir_start, refined.pre_tir_end) == (boundary_start, boundary_end)
    assert refined.tir_status == expected_status
    assert refined.tir_present is True
    assert refined.tir_boundary_override is expected_override
    assert refined.tir_alignment_capped is False
    assert refined.tir_candidate_count == 1


def test_coverage_veto_retains_candidate_without_changing_boundary(tmp_path: Path) -> None:
    sequence, expected = _planted_repeat()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    boundary = RefinedBoundary(
        scaffold="contig",
        start=120,
        end=400,
        seed_id="seed-1",
        original_start=150,
        original_end=350,
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 150, 350),),
    )

    refined = refine_boundaries_with_terminal_repeats(
        [boundary],
        allowed_intervals=[(boundary.start, boundary.end)],
        raw_genome_path=genome_path,
        validated_markers=[_Marker("contig", 250, 270, "marker-1")],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": [pORF(id="unclassified", scaffold="contig", start=110, end=130, strand="+")]},
        extension_bp=500,
        metadata_families={},
    )[0]

    assert (refined.start, refined.end) == (120, 400)
    assert refined.tir_status == "taxonomy_incomplete"
    assert refined.tir_present is True
    assert refined.tir_boundary_override is False
    assert (refined.tir_left_start, refined.tir_right_end) == (expected[0], expected[3])


def test_no_marker_anchor_is_explicit_and_coordinate_neutral(tmp_path: Path) -> None:
    sequence, _ = _planted_repeat()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    boundary = RefinedBoundary(
        scaffold="contig",
        start=120,
        end=400,
        seed_id="seed-1",
        original_start=150,
        original_end=350,
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 150, 350),),
    )

    refined = refine_boundaries_with_terminal_repeats(
        [boundary],
        allowed_intervals=[(boundary.start, boundary.end)],
        raw_genome_path=genome_path,
        validated_markers=[],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
        metadata_families={},
    )[0]

    assert (refined.start, refined.end) == (120, 400)
    assert refined.tir_status == "no_marker_anchor"
    assert refined.tir_present is False
    assert refined.tir_candidate_count == 0


@pytest.mark.parametrize("parent_start", [0, 120])
def test_disjoint_tir_pairs_partition_parent_and_keep_marker_neighbor(
    tmp_path: Path,
    parent_start: int,
) -> None:
    sequence, expected_pairs = _adjacent_repeat_sequence()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    parent = RefinedBoundary(
        scaffold="contig",
        start=parent_start,
        end=len(sequence),
        seed_id="seed-1",
        original_start=200,
        original_end=1100,
        confirmed_seeds=(
            ConfirmedSeed("seed-1", "contig", 200, 260),
            ConfirmedSeed("seed-2", "contig", 620, 680),
            ConfirmedSeed("seed-3", "contig", 980, 1100),
        ),
        seed_sources=["hhg", "marker_validation"],
        seed_confidence="high",
        seed_hhg_score=0.99,
        predicted_family="MIXED",
        region_classification_ncldv_markers=99,
        region_classification_vp_plv_markers=99,
        region_classification_mirus_markers=99,
        confidence=0.91,
        posterior_probability=0.92,
        cub_deviation=0.4,
    )
    markers = [
        _Marker("contig", 220, 235, "first|aa1-15", hmm_target="gvogm0003"),
        _Marker("contig", 235, 250, "first|aa16-30", hmm_target="gvogm0004"),
        _Marker("contig", 640, 660, "second|aa1-20", hmm_target="plv_mcp"),
        _Marker("contig", 1000, 1020, "neighbor|aa1-20", hmm_target="mirus_core"),
    ]

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        allowed_intervals=[(0, len(sequence))],
        raw_genome_path=genome_path,
        validated_markers=markers,
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
        metadata_families={},
    )

    assert [(boundary.start, boundary.end) for boundary in refined] == [
        *expected_pairs,
        (840, len(sequence)),
    ]
    assert [boundary.tir_status for boundary in refined] == [
        "detected",
        "detected",
        "not_detected",
    ]
    assert [boundary.predicted_family for boundary in refined] == [
        "NCLDV",
        "PPV",
        "MIRUS",
    ]
    # One best HMM model per protein, not all alternative names, contributes family credit.
    assert [boundary.region_classification_ncldv_markers for boundary in refined] == [1, 0, 0]
    assert [(boundary.original_start, boundary.original_end) for boundary in refined] == [
        (200, 260),
        (620, 680),
        (980, 1100),
    ]
    assert [boundary.confirmed_seeds[0].seed_id for boundary in refined] == ["seed-1", "seed-2", "seed-3"]
    assert all((boundary.pre_tir_start, boundary.pre_tir_end) == (parent_start, len(sequence)) for boundary in refined)
    assert [boundary.seed_id for boundary in refined] == ["seed-1", "seed-2", "seed-3"]
    assert all(boundary.confidence == 0.0 for boundary in refined)
    assert all(boundary.posterior_probability == 0.0 for boundary in refined)
    assert all(boundary.cub_deviation == 0.0 for boundary in refined)


def test_partition_cannot_split_whole_seed_for_viral_taxonomy_neighbor(tmp_path: Path) -> None:
    sequence, expected_pairs = _adjacent_repeat_sequence()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    parent = RefinedBoundary(
        scaffold="contig",
        start=0,
        end=len(sequence),
        seed_id="seed-1",
        original_start=0,
        original_end=len(sequence),
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 0, len(sequence)),),
    )
    markers = [
        _Marker("contig", 220, 240, "first", hmm_target="gvogm0003"),
        _Marker("contig", 640, 660, "second", hmm_target="plv_mcp"),
    ]
    viral_gene = pORF(
        id="viral-neighbor",
        scaffold="contig",
        start=1000,
        end=1080,
        strand="+",
    )
    taxonomy = GeneTaxonomy(
        porf_id=viral_gene.id,
        scaffold="contig",
        start=viral_gene.start,
        end=viral_gene.end,
        top10_prefixes=["NCLDV__"],
        top10_pidents=[70.0],
        has_viral=True,
        has_hit=True,
    )

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        allowed_intervals=[(parent.start, parent.end)],
        raw_genome_path=genome_path,
        validated_markers=markers,
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={viral_gene.id: taxonomy},
        proteome_index={"contig": [viral_gene]},
        extension_bp=500,
        metadata_families={},
    )

    assert len(refined) == 1
    assert (refined[0].start, refined[0].end) == (parent.start, parent.end)
    assert refined[0].confirmed_seeds == parent.confirmed_seeds
    assert refined[0].tir_status == "ambiguous"
    assert refined[0].tir_present is True
    assert refined[0].tir_boundary_override is False


def test_partition_retains_parent_when_marker_neighbor_has_no_confirmed_seed(
    tmp_path: Path,
) -> None:
    sequence, expected = _planted_repeat()
    sequence += _random_dna(300, random.Random(407))
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    parent = RefinedBoundary(
        scaffold="contig",
        start=0,
        end=len(sequence),
        seed_id="seed-1",
        original_start=200,
        original_end=300,
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 200, 300),),
    )

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        allowed_intervals=[(parent.start, parent.end)],
        raw_genome_path=genome_path,
        validated_markers=[
            _Marker("contig", 220, 240, "tir-marker"),
            _Marker("contig", 640, 660, "neighbor-marker"),
        ],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
        metadata_families={},
    )

    assert len(refined) == 1
    assert (refined[0].start, refined[0].end) == (parent.start, parent.end)
    assert refined[0].confirmed_seeds == parent.confirmed_seeds
    assert refined[0].tir_status == "ambiguous"
    assert refined[0].tir_present is True
    assert refined[0].tir_boundary_override is False


def test_partition_cannot_emit_unseeded_call_beyond_parent(tmp_path: Path) -> None:
    planted, _ = _planted_repeat()
    sequence = _random_dna(400, random.Random(503)) + planted
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    parent = RefinedBoundary(
        scaffold="contig",
        start=100,
        end=800,
        seed_id="seed-1",
        original_start=100,
        original_end=800,
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 100, 800),),
    )
    viral_gene = pORF(
        id="viral-residual",
        scaffold="contig",
        start=120,
        end=150,
        strand="+",
    )
    taxonomy = GeneTaxonomy(
        porf_id=viral_gene.id,
        scaffold="contig",
        start=viral_gene.start,
        end=viral_gene.end,
        top10_prefixes=["NCLDV__"],
        top10_pidents=[70.0],
        has_viral=True,
        has_hit=True,
    )

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        allowed_intervals=[(parent.start, parent.end)],
        raw_genome_path=genome_path,
        validated_markers=[_Marker("contig", 620, 640, "tir-marker")],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={viral_gene.id: taxonomy},
        proteome_index={"contig": [viral_gene]},
        extension_bp=500,
        metadata_families={},
    )

    assert len(refined) == 1
    assert (refined[0].start, refined[0].end) == (parent.start, parent.end)
    assert refined[0].confirmed_seeds == parent.confirmed_seeds
    assert refined[0].tir_status == "ambiguous"
    assert refined[0].tir_present is True
    assert refined[0].tir_boundary_override is False


def test_identical_tir_cross_pair_keeps_parent_ambiguous(tmp_path: Path) -> None:
    sequence, _ = _adjacent_repeat_sequence(identical_arms=True)
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    parent = RefinedBoundary(
        scaffold="contig",
        start=0,
        end=len(sequence),
        seed_id="seed-1",
        original_start=0,
        original_end=len(sequence),
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 0, len(sequence)),),
    )

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        allowed_intervals=[(parent.start, parent.end)],
        raw_genome_path=genome_path,
        validated_markers=[
            _Marker("contig", 220, 240, "first"),
            _Marker("contig", 640, 660, "second"),
        ],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
        metadata_families={},
    )

    assert len(refined) == 1
    assert (refined[0].start, refined[0].end) == (parent.start, parent.end)
    assert refined[0].tir_status == "ambiguous"
    assert refined[0].tir_candidate_count >= 3
    assert refined[0].tir_boundary_override is False


def test_crossing_porf_does_not_override_tir_partition(tmp_path: Path) -> None:
    sequence, _ = _adjacent_repeat_sequence()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    parent = RefinedBoundary(
        scaffold="contig",
        start=0,
        end=len(sequence),
        seed_id="seed-1",
        original_start=200,
        original_end=1100,
        confirmed_seeds=(
            ConfirmedSeed("seed-1", "contig", 200, 260),
            ConfirmedSeed("seed-2", "contig", 620, 680),
            ConfirmedSeed("seed-3", "contig", 980, 1100),
        ),
    )
    crossing_gene = pORF(
        id="crossing",
        scaffold="contig",
        start=830,
        end=850,
        strand="+",
    )
    taxonomy = GeneTaxonomy(
        porf_id=crossing_gene.id,
        scaffold="contig",
        start=crossing_gene.start,
        end=crossing_gene.end,
        has_viral=False,
        has_hit=False,
    )
    viral_neighbor = pORF(
        id="viral-neighbor",
        scaffold="contig",
        start=1000,
        end=1080,
        strand="+",
    )
    neighbor_taxonomy = GeneTaxonomy(
        porf_id=viral_neighbor.id,
        scaffold="contig",
        start=viral_neighbor.start,
        end=viral_neighbor.end,
        top10_prefixes=["NCLDV__"],
        top10_pidents=[70.0],
        has_viral=True,
        has_hit=True,
    )

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        allowed_intervals=[(parent.start, parent.end)],
        raw_genome_path=genome_path,
        validated_markers=[
            _Marker("contig", 220, 240, "first"),
            _Marker("contig", 640, 660, "second"),
        ],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={
            crossing_gene.id: taxonomy,
            viral_neighbor.id: neighbor_taxonomy,
        },
        proteome_index={"contig": [crossing_gene, viral_neighbor]},
        extension_bp=500,
        metadata_families={},
    )

    assert [(boundary.start, boundary.end) for boundary in refined] == [
        (100, 420),
        (520, 840),
        (840, len(sequence)),
    ]
    assert [boundary.tir_status for boundary in refined] == [
        "detected",
        "detected",
        "not_detected",
    ]


def test_partition_vetoes_complete_marker_group_crossing_child_boundary(
    tmp_path: Path,
) -> None:
    sequence, _ = _adjacent_repeat_sequence()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    parent = RefinedBoundary(
        scaffold="contig",
        start=0,
        end=len(sequence),
        seed_id="seed-1",
        original_start=0,
        original_end=len(sequence),
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 0, len(sequence)),),
    )

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        allowed_intervals=[(parent.start, parent.end)],
        raw_genome_path=genome_path,
        validated_markers=[
            _Marker("contig", 220, 240, "first"),
            _Marker("contig", 640, 660, "second"),
            _Marker("contig", 830, 838, "crossing|aa1-8"),
            _Marker("contig", 842, 850, "crossing|aa9-16"),
        ],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
        metadata_families={},
    )

    assert len(refined) == 1
    assert (refined[0].start, refined[0].end) == (parent.start, parent.end)
    assert refined[0].tir_status == "ambiguous"
    assert refined[0].tir_candidate_count == 2


@pytest.mark.parametrize(
    ("allowed", "status"), [([(180, 245), (260, 340)], "ambiguous"), ([(0, 520), (0, 520)], "shared_pair")]
)
def test_shared_pair_does_not_collapse_distinct_boundaries(
    tmp_path: Path, allowed: list[tuple[int, int]], status: str
) -> None:
    sequence, expected = _planted_repeat()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    boundaries = [
        RefinedBoundary(
            scaffold="contig",
            start=start,
            end=end,
            seed_id=seed_id,
            original_start=start,
            original_end=end,
            confirmed_seeds=(ConfirmedSeed(seed_id, "contig", start, end),),
        )
        for seed_id, start, end in (
            ("seed-1", 180, 245),
            ("seed-2", 260, 340),
        )
    ]
    mappings = {
        boundary.seed_id: SeedGeneMapping(
            seed_id=boundary.seed_id,
            scaffold="contig",
            seed_start=boundary.start,
            seed_end=boundary.end,
            flank_start_bp=0,
            flank_end_bp=len(sequence),
        )
        for boundary in boundaries
    }

    refined = refine_boundaries_with_terminal_repeats(
        boundaries,
        allowed_intervals=allowed,
        raw_genome_path=genome_path,
        validated_markers=[
            _Marker("contig", 220, 230, "marker-1"),
            _Marker("contig", 280, 290, "marker-2"),
        ],
        boundary_diamond_query=GenomeDiamondQuery(seed_gene_mappings=mappings),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
        metadata_families={},
    )

    assert [(boundary.start, boundary.end) for boundary in refined] == [(180, 245), (260, 340)]
    assert {boundary.tir_status for boundary in refined} == {status}
    assert all(boundary.tir_present for boundary in refined)
    assert not any(boundary.tir_boundary_override for boundary in refined)
    assert {(boundary.tir_left_start, boundary.tir_right_end) for boundary in refined} == {(expected[0], expected[3])}


@pytest.mark.parametrize(("allowed_first", "status"), [((180, 245), "ambiguous"), ((0, 520), "shared_pair")])
def test_outward_tir_proposal_does_not_replace_existing_candidate(
    tmp_path: Path, allowed_first: tuple[int, int], status: str
) -> None:
    sequence, expected = _planted_repeat()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    proposing = RefinedBoundary(
        scaffold="contig",
        start=180,
        end=245,
        seed_id="seed-1",
        original_start=180,
        original_end=245,
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 180, 245),),
    )
    untouched = RefinedBoundary(
        scaffold="contig",
        start=expected[0],
        end=expected[3],
        seed_id="seed-2",
        original_start=150,
        original_end=200,
        confirmed_seeds=(ConfirmedSeed("seed-2", "contig", 150, 200),),
    )
    mappings = {
        boundary.seed_id: SeedGeneMapping(
            seed_id=boundary.seed_id,
            scaffold="contig",
            seed_start=boundary.original_start,
            seed_end=boundary.original_end,
            flank_start_bp=0,
            flank_end_bp=len(sequence),
        )
        for boundary in (proposing, untouched)
    }

    refined = refine_boundaries_with_terminal_repeats(
        [proposing, untouched],
        allowed_intervals=[allowed_first, (untouched.start, untouched.end)],
        raw_genome_path=genome_path,
        validated_markers=[_Marker("contig", 220, 230, "marker-1")],
        boundary_diamond_query=GenomeDiamondQuery(seed_gene_mappings=mappings),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
        metadata_families={},
    )

    assert [(boundary.start, boundary.end) for boundary in refined] == [
        (untouched.start, untouched.end),
        (proposing.start, proposing.end),
    ]
    by_seed = {boundary.seed_id: boundary for boundary in refined}
    assert by_seed["seed-1"].tir_status == status
    assert by_seed["seed-1"].tir_boundary_override is False
    assert by_seed["seed-2"].tir_status == "detected"
    assert by_seed["seed-2"].tir_boundary_override is False

    state_path = tmp_path / "refined_state.json"
    write_phase2_state(state_path, refined)
    assert load_phase2_state(state_path) == refined


def test_shared_pair_resolution_handles_restoration_collision() -> None:
    originals = [
        RefinedBoundary(
            scaffold="contig",
            start=start,
            end=end,
            seed_id=f"source-{start}",
            confirmed_seeds=(ConfirmedSeed(f"source-{start}", "contig", start + 10, end - 10),),
        )
        for start, end in ((100, 200), (250, 350), (400, 500))
    ]

    def proposal(boundary: RefinedBoundary, start: int, end: int) -> RefinedBoundary:
        return replace(
            boundary,
            start=start,
            end=end,
            pre_tir_start=boundary.start,
            pre_tir_end=boundary.end,
            tir_present=True,
            tir_status="detected",
            tir_left_start=start,
            tir_left_end=start + 30,
            tir_right_start=end - 30,
            tir_right_end=end,
            tir_identity=1.0,
            tir_boundary_override=True,
            tir_alignment_length=30,
            tir_candidate_count=1,
        )

    proposed = [
        proposal(originals[0], 50, 550),
        proposal(originals[1], 50, 550),
        proposal(originals[2], 100, 200),
    ]

    resolved = terminal_repeats._resolve_partition_conflicts(
        originals,
        [[boundary] for boundary in proposed],
    )

    assert sorted((boundary.start, boundary.end) for boundary in resolved) == [
        (100, 200),
        (250, 350),
        (400, 500),
    ]
    assert {boundary.tir_status for boundary in resolved} == {"shared_pair"}
    assert not any(boundary.tir_boundary_override for boundary in resolved)


@pytest.mark.parametrize(("stops", "expected_allowed"), [((None, None), (0, 700)), ((50, 650), (50, 650))])
def test_phase2_applies_tir_after_merge_before_downstream_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    stops: tuple[int | None, int | None],
    expected_allowed: tuple[int, int],
) -> None:
    events: list[tuple[str, object]] = []
    retained_marker = SimpleNamespace(
        scaffold="contig",
        start=150,
        end=160,
        query_porf="tax",
        validation_status="validated",
        hmm_target="GVOGm0003",
        hmm_score=100.0,
        hmm_evalue=1e-20,
        top10_prefixes="NCLDV__",
        top10_pidents="50.0",
        top10_targets="",
    )
    raw_path = tmp_path / "raw.fna"
    raw_path.write_text(">contig\n" + "ACGT" * 200 + "\n")
    mappings = {
        seed_id: SeedGeneMapping(
            seed_id=seed_id,
            scaffold="contig",
            seed_start=start,
            seed_end=end,
            flank_start_bp=0 if seed_id == "seed-1" else 200,
            flank_end_bp=600 if seed_id == "seed-1" else 700,
        )
        for seed_id, start, end in (("seed-1", 100, 250), ("seed-2", 300, 450))
    }
    mappings["separate"] = SeedGeneMapping("separate", "contig", 200, 300, flank_start_bp=0, flank_end_bp=800)
    query = GenomeDiamondQuery(seed_gene_mappings=mappings)
    taxonomy = GeneTaxonomy(porf_id="tax", scaffold="contig", start=150, end=160)

    monkeypatch.setattr(phase2, "build_proteome_index", lambda _path: {"contig": []})
    monkeypatch.setattr(phase2, "_require_phase2b_gene_taxonomy_db", lambda *_args, **_kwargs: tmp_path / "db")
    monkeypatch.setattr(
        phase2,
        "_run_boundary_diamond",
        lambda **_kwargs: phase2._BoundaryEvidence(
            taxonomy_map={"tax": taxonomy},
            control_stats=None,
            query=query,
        ),
    )

    def fake_merge(boundaries: list[RefinedBoundary], **_kwargs: object) -> list[RefinedBoundary]:
        events.append(("merge", [(boundary.start, boundary.end) for boundary in boundaries]))
        return [
            replace(
                boundaries[0],
                start=80,
                end=470,
                original_start=100,
                original_end=450,
                confirmed_seeds=boundaries[0].confirmed_seeds + boundaries[1].confirmed_seeds,
                left_stop=stops[0],
                right_stop=stops[1],
            )
        ]

    def fake_tir(boundaries: list[RefinedBoundary], **kwargs: object) -> list[RefinedBoundary]:
        assert kwargs["validated_markers"][0] is retained_marker
        assert kwargs["allowed_intervals"] == [expected_allowed]
        events.append(("tir", (boundaries[0].start, boundaries[0].end, kwargs["raw_genome_path"])))
        boundary = boundaries[0]
        return [
            replace(
                boundary,
                start=90,
                end=460,
                pre_tir_start=boundary.start,
                pre_tir_end=boundary.end,
                tir_present=True,
                tir_status="detected",
                tir_left_start=90,
                tir_left_end=120,
                tir_right_start=430,
                tir_right_end=460,
                tir_identity=1.0,
                tir_boundary_override=True,
                tir_alignment_length=30,
                tir_candidate_count=1,
                tir_scan_start=0,
                tir_scan_end=800,
            )
        ]

    def fake_composition(boundaries: list[RefinedBoundary], **_kwargs: object) -> None:
        events.append(("composition", (boundaries[0].start, boundaries[0].end)))

    monkeypatch.setattr(phase2, "merge_reachable_boundaries", fake_merge)
    monkeypatch.setattr(terminal_repeats, "refine_boundaries_with_terminal_repeats", fake_tir)
    monkeypatch.setattr(phase2, "_recalculate_boundary_composition", fake_composition)

    config = PipelineConfig().with_overrides(
        resume=False,
        gene_taxonomy_faa_db=tmp_path / "db",
        boundary_host_trim_enabled=False,
        boundary_taxonomy_ml_enabled=False,
        threads=1,
        gene_taxonomy_threads=1,
    )
    result = phase2._run_phase2_subflow(
        masked_path=tmp_path / "masked.fna",
        proteome_path=tmp_path / "proteome.faa",
        merged_seeds=[
            MergedSeed(
                scaffold="contig",
                start=100,
                end=250,
                seed_id="seed-1",
                confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 100, 250),),
            ),
            MergedSeed(
                scaffold="contig",
                start=300,
                end=450,
                seed_id="seed-2",
                confirmed_seeds=(ConfirmedSeed("seed-2", "contig", 300, 450),),
            ),
        ],
        validated_markers=[retained_marker],
        host_signature_model=None,
        output_dir=tmp_path,
        genome_id="genome",
        config=config,
        genome_start_time=0.0,
        logger=logging.getLogger(__name__),
        raw_genome_path=raw_path,
        metadata_families={},
    )

    assert [(name, value) for name, value in events] == [
        ("merge", [(100, 250), (300, 450)]),
        ("tir", (80, 470, raw_path)),
        ("composition", (90, 460)),
    ]
    assert (result.refined_boundaries[0].start, result.refined_boundaries[0].end) == (90, 460)
    assert (tmp_path / "phase2" / "refined_boundaries.bed").read_text().startswith("contig\t90\t460\t")


@pytest.mark.parametrize(
    ("parent_bounds", "allowed", "stop", "expected"),
    [
        ((120, 400), (0, 520), None, (100, 420)),
        ((100, 400), (0, 520), None, (100, 420)),
        ((50, 400), (0, 520), None, (100, 420)),
        ((120, 400), (101, 520), None, (120, 400)),
        ((120, 400), (100, 420), None, (100, 420)),
        ((120, 400), (0, 520), 410, (120, 400)),
        ((120, 400), (0, 520), 420, (100, 420)),
    ],
)
def test_outward_tir_respects_walked_interval_and_stop(
    tmp_path: Path,
    parent_bounds: tuple[int, int],
    allowed: tuple[int, int],
    stop: int | None,
    expected: tuple[int, int],
) -> None:
    sequence, _ = _planted_repeat()
    genome = tmp_path / "raw.fna"
    genome.write_text(f">contig\n{sequence}\n", encoding="utf-8")
    parent = RefinedBoundary(
        "contig",
        *parent_bounds,
        seed_id="seed-1",
        original_start=150,
        original_end=350,
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 150, 350),),
        right_stop=stop,
    )
    query = _query(len(sequence))

    result = refine_boundaries_with_terminal_repeats(
        [parent],
        allowed_intervals=[allowed],
        raw_genome_path=genome,
        validated_markers=[_Marker("contig", 250, 270, "marker-1")],
        boundary_diamond_query=query,
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
        metadata_families={},
    )
    phase2._annotate_boundary_repeats(result, raw_genome_path=genome, boundary_diamond_query=query, extension_bp=500)
    state = tmp_path / "refined.json"
    write_phase2_state(state, result)
    loaded = load_phase2_state(state)

    assert (result[0].start, result[0].end) == expected
    assert result[0].confirmed_seeds == parent.confirmed_seeds
    assert (result[0].pre_tir_start, result[0].pre_tir_end) == parent_bounds
    assert result[0].repeat_evidence is not None
    assert (result[0].repeat_evidence.assessed_start, result[0].repeat_evidence.assessed_end) == expected
    assert (result[0].repeat_evidence.parent_start, result[0].repeat_evidence.parent_end) == parent_bounds
    assert loaded == result


@pytest.mark.parametrize(
    ("parent_bounds", "marker_span", "expected"),
    [
        ((120, 400), (410, 450), (120, 400, "ambiguous")),
        ((120, 400), (410, 700), (120, 400, "ambiguous")),
        ((120, 400), (0, 110), (120, 400, "ambiguous")),
        ((50, 470), (410, 700), (100, 420, "detected")),
        ((100, 420), (410, 700), (100, 420, "detected")),
    ],
)
def test_flank_marker_veto_applies_only_to_new_outward_endpoints(
    tmp_path: Path,
    parent_bounds: tuple[int, int],
    marker_span: tuple[int, int],
    expected: tuple[int, int, str],
) -> None:
    sequence, _ = _planted_repeat()
    sequence += "A" * 300
    genome = tmp_path / "raw.fna"
    genome.write_text(f">contig\n{sequence}\n", encoding="utf-8")
    parent = RefinedBoundary(
        "contig",
        *parent_bounds,
        seed_id="seed-1",
        original_start=150,
        original_end=350,
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 150, 350),),
    )
    result = refine_boundaries_with_terminal_repeats(
        [parent],
        allowed_intervals=[(0, 520)],
        raw_genome_path=genome,
        validated_markers=[_Marker("contig", 250, 270, "anchor"), _Marker("contig", *marker_span, "flank")],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
        metadata_families={},
    )
    assert (result[0].start, result[0].end, result[0].tir_status) == expected


def _outward_proposal(parent: RefinedBoundary, start: int, end: int) -> RefinedBoundary:
    return replace(
        parent,
        start=start,
        end=end,
        pre_tir_start=parent.start,
        pre_tir_end=parent.end,
        tir_present=True,
        tir_status="detected",
        tir_left_start=start,
        tir_left_end=start + 30,
        tir_right_start=end - 30,
        tir_right_end=end,
        tir_identity=1.0,
        tir_boundary_override=True,
        tir_alignment_length=30,
        tir_candidate_count=1,
    )


@pytest.mark.parametrize(("end", "expected_end", "status"), [(300, 200, "ambiguous"), (250, 250, "detected")])
@pytest.mark.parametrize("neighbor_bounds", [(250, 400), (260, 395)])
def test_new_partial_overlap_reverts_but_touching_neighbor_is_allowed(
    end: int, expected_end: int, status: str, neighbor_bounds: tuple[int, int]
) -> None:
    left = RefinedBoundary("s", 100, 200, seed_id="left", confirmed_seeds=(ConfirmedSeed("left", "s", 120, 180),))
    right = RefinedBoundary("s", 250, 400, seed_id="right", confirmed_seeds=(ConfirmedSeed("right", "s", 260, 390),))
    neighbor = right if neighbor_bounds == (250, 400) else _outward_proposal(right, *neighbor_bounds)
    result = terminal_repeats._resolve_partition_conflicts(
        [left, right], [[_outward_proposal(left, 100, end)], [neighbor]]
    )
    assert (result[0].end, result[0].tir_status) == (expected_end, status)
    assert result[1] == neighbor
    assert result[0].confirmed_seeds == left.confirmed_seeds


def test_partial_overlap_restoration_exposes_another_conflict() -> None:
    parents = [
        RefinedBoundary(
            "s", start, end, seed_id=str(start), confirmed_seeds=(ConfirmedSeed(str(start), "s", start + 10, end - 10),)
        )
        for start, end in [(100, 200), (250, 350), (20, 80)]
    ]
    parents[0].confirmed_seeds = (ConfirmedSeed("100", "s", 160, 180),)
    plans = [[_outward_proposal(parents[0], 150, 280)], [parents[1]], [_outward_proposal(parents[2], 20, 140)]]
    result = terminal_repeats._resolve_partition_conflicts(parents, plans)
    assert [(boundary.start, boundary.end) for boundary in result] == [(20, 80), (100, 200), (250, 350)]
    assert [boundary.tir_status for boundary in result] == ["ambiguous", "ambiguous", "not_assessed"]


@pytest.mark.parametrize(
    ("mode", "expected_allowed", "expected_bounds"),
    [
        ("no_host", (0, 520), (100, 420)),
        ("CRESS", (150, 350), (150, 350)),
        ("A4", (0, 520), (100, 420)),
        ("host", (0, 360), (150, 350)),
    ],
)
def test_phase2_real_tir_uses_walked_authority_and_final_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mode: str,
    expected_allowed: tuple[int, int],
    expected_bounds: tuple[int, int],
) -> None:
    from virosync.pipeline.host_signatures import HostSignatureModel
    from virosync.pipeline.phase2 import seed_extension

    sequence, _ = _planted_repeat()
    genome = tmp_path / "raw.fna"
    genome.write_text(f">contig\n{sequence}\n", encoding="utf-8")
    genes = [pORF(f"host-{start}", "contig", start, start + 10) for start in (360, 380, 400)]
    taxonomy = {
        gene.id: GeneTaxonomy(
            porf_id=gene.id,
            scaffold="contig",
            start=gene.start,
            end=gene.end,
            has_hit=True,
            top10_prefixes=["EUK__"],
            top10_pidents=[80.0],
        )
        for gene in genes
    }
    query = _query(len(sequence))
    marker = SimpleNamespace(
        scaffold="contig",
        start=250,
        end=270,
        query_porf="anchor",
        validation_status="validated",
        hmm_target="GVOGm0003",
        hmm_score=100.0,
        hmm_evalue=1e-20,
        top10_prefixes="NCLDV__",
        top10_pidents="50.0",
        top10_targets="",
    )
    monkeypatch.setattr(phase2, "build_proteome_index", lambda _path: {"contig": genes})
    monkeypatch.setattr(phase2, "_require_phase2b_gene_taxonomy_db", lambda *_args, **_kwargs: tmp_path / "db")
    monkeypatch.setattr(
        phase2,
        "_run_boundary_diamond",
        lambda **_kwargs: phase2._BoundaryEvidence(taxonomy_map=taxonomy, control_stats=None, query=query),
    )
    monkeypatch.setattr(seed_extension, "host_lineage_share", lambda *_args: 1.0)
    native_refiner = terminal_repeats.refine_boundaries_with_terminal_repeats
    native_composition = phase2._recalculate_boundary_composition
    seen = []

    def check_refiner(boundaries, **kwargs):
        assert kwargs["allowed_intervals"] == [expected_allowed]
        return native_refiner(boundaries, **kwargs)

    def check_composition(boundaries, **kwargs):
        seen.extend((boundary.start, boundary.end) for boundary in boundaries)
        native_composition(boundaries, **kwargs)

    monkeypatch.setattr(terminal_repeats, "refine_boundaries_with_terminal_repeats", check_refiner)
    monkeypatch.setattr(phase2, "_recalculate_boundary_composition", check_composition)
    config = PipelineConfig().with_overrides(
        resume=False,
        gene_taxonomy_faa_db=tmp_path / "db",
        boundary_host_trim_enabled=mode != "no_host",
        boundary_taxonomy_ml_enabled=False,
        ablation_id="A4" if mode == "A4" else "A0",
        threads=1,
        gene_taxonomy_threads=1,
    )
    core = ConfirmedSeed("seed-1", "contig", 150, 350)
    result = phase2._run_phase2_subflow(
        masked_path=genome,
        raw_genome_path=genome,
        proteome_path=tmp_path / "proteome.faa",
        merged_seeds=[
            MergedSeed(
                "contig",
                150,
                350,
                seed_id="seed-1",
                confirmed_seeds=(core,),
                predicted_family="CRESS" if mode == "CRESS" else "NCLDV",
            )
        ],
        validated_markers=[marker],
        host_signature_model=HostSignatureModel(token_weights={"host": 1.0}, lineage_weights={"euk|host": 1.0}),
        output_dir=tmp_path,
        genome_id="genome",
        config=config,
        genome_start_time=0.0,
        logger=logging.getLogger(__name__),
        metadata_families={},
    )
    final = result.refined_boundaries[0]
    assert (final.start, final.end) == expected_bounds
    assert seen == [expected_bounds]
    assert final.confirmed_seeds == (core,)
    assert final.region_classification_ncldv_markers == 1
    assert final.candidate_id == f"EVE_contig_{expected_bounds[0]}-{expected_bounds[1]}"
    assert final.repeat_evidence is not None
    assert (final.repeat_evidence.assessed_start, final.repeat_evidence.assessed_end) == expected_bounds
    assert (final.repeat_evidence.parent_start, final.repeat_evidence.parent_end) == (150, 350)
    assert load_phase2_state(tmp_path / "phase2/refined_state.json") == result.refined_boundaries


def test_outward_tir_missing_taxonomy_in_added_flank_vetoes_change(tmp_path: Path) -> None:
    sequence, _ = _planted_repeat()
    genome = tmp_path / "raw.fna"
    genome.write_text(f">contig\n{sequence}\n", encoding="utf-8")
    parent = RefinedBoundary(
        "contig",
        120,
        400,
        seed_id="seed-1",
        original_start=150,
        original_end=350,
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 150, 350),),
    )
    result = refine_boundaries_with_terminal_repeats(
        [parent],
        allowed_intervals=[(0, 520)],
        raw_genome_path=genome,
        validated_markers=[_Marker("contig", 250, 270, "anchor")],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": [pORF("missing", "contig", 405, 415)]},
        extension_bp=500,
        metadata_families={},
    )
    assert (result[0].start, result[0].end, result[0].tir_status) == (120, 400, "taxonomy_incomplete")


def test_shared_pair_restoration_retains_identical_not_partial_collision_evidence() -> None:
    first = ConfirmedSeed("a1", "s", 220, 260)
    second = ConfirmedSeed("a2", "s", 520, 560)
    parent = RefinedBoundary("s", 200, 600, seed_id="a1", confirmed_seeds=(first, second))
    neighbor = RefinedBoundary("s", 50, 150, seed_id="b", confirmed_seeds=(ConfirmedSeed("b", "s", 70, 90),))
    peer = RefinedBoundary("s", 620, 680, seed_id="c", confirmed_seeds=(ConfirmedSeed("c", "s", 630, 650),))
    partial = replace(
        _outward_proposal(parent, 100, 300), confirmed_seeds=(first,), tir_identity=0.93, tsd_sequence="TGCA"
    )
    identical = replace(
        _outward_proposal(parent, 500, 700),
        seed_id="a2",
        confirmed_seeds=(second,),
        tir_identity=0.97,
        tsd_sequence="ACGT",
    )

    resolved = terminal_repeats._resolve_partition_conflicts(
        [parent, neighbor, peer],
        [[partial, identical], [neighbor], [_outward_proposal(peer, 500, 700)]],
    )
    restored = next(boundary for boundary in resolved if boundary.seed_id == "a1")

    assert (restored.start, restored.end, restored.tir_status) == (200, 600, "shared_pair")
    assert restored.confirmed_seeds == (first, second)
    assert (
        restored.tir_left_start,
        restored.tir_left_end,
        restored.tir_right_start,
        restored.tir_right_end,
        restored.tir_identity,
        restored.tsd_sequence,
    ) == (500, 530, 670, 700, 0.97, "ACGT")


@dataclass(frozen=True, slots=True)
class _Fixture:
    genome_path: Path
    boundaries: list[RefinedBoundary]
    allowed_intervals: list[tuple[int, int]]
    markers: list[_Marker]
    query: GenomeDiamondQuery
    taxonomy_map: dict[str, GeneTaxonomy]
    proteome_index: dict[str, list[pORF]]

    def refine(self, workers: int = 1) -> list[RefinedBoundary]:
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
            metadata_families={},
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


def test_serial_planning_preserves_populated_query_role_ids(fixture: _Fixture) -> None:
    query = fixture.query
    query.eve_porf_ids = {"adj-1": ["adj_2"], "shared-1": ["shared_1"]}
    query.boundary_porf_ids = {"adj-1": ["adj_1"], "trim-1": ["trim_1"]}
    query.control_porf_ids = ["trim_2"]
    query.all_porf_ids = ["adj_1", "adj_2", "shared_1", "trim_1", "trim_2"]
    original_roles = pickle.dumps(
        (query.eve_porf_ids, query.boundary_porf_ids, query.control_porf_ids, query.all_porf_ids)
    )

    fixture.refine()

    assert (
        pickle.dumps((query.eve_porf_ids, query.boundary_porf_ids, query.control_porf_ids, query.all_porf_ids))
        == original_roles
    )


def test_serial_fixture_changes_boundaries_by_tir(fixture: _Fixture) -> None:
    refined = fixture.refine()

    assert len(refined) == 6
    assert sum(boundary.tir_boundary_override for boundary in refined) >= 3
    assert {boundary.tir_status for boundary in refined} >= {"detected", "shared_pair"}
    assert [(boundary.start, boundary.end) for boundary in refined if boundary.scaffold == "trim"] == [(100, 420)]
    assert {core for boundary in refined for core in boundary.confirmed_seeds} == {
        core for boundary in fixture.boundaries for core in boundary.confirmed_seeds
    }
    assert all(
        boundary.start <= core.start < core.end <= boundary.end
        for boundary in refined
        for core in boundary.confirmed_seeds
    )


def test_serial_planning_raises_same_validation_error(fixture: _Fixture) -> None:
    fixture.allowed_intervals[2] = (300, 340)

    with pytest.raises(ValueError, match="allowed interval must contain parent"):
        fixture.refine()


def test_serial_planning_raises_for_absent_scaffold(fixture: _Fixture) -> None:
    fixture.boundaries[3] = _boundary("missing", "missing-1", 50, 470, (150, 350))

    with pytest.raises(ValueError, match="boundary scaffold absent from raw genome: missing"):
        fixture.refine()


def test_concurrent_pools_preserve_serial_results_and_deep_tempdir(tmp_path: Path, monkeypatch) -> None:
    deep_tempdir = tmp_path / ("deep-" * 25)
    deep_tempdir.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(deep_tempdir))
    (tmp_path / "first").mkdir()
    (tmp_path / "second").mkdir()
    first = _build_fixture(tmp_path / "first")
    second = _build_fixture(tmp_path / "second")
    del second.boundaries[3], second.allowed_intervals[3]
    serial = {"first": first.refine(), "second": second.refine()}
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
    assert tempfile.gettempdir() == str(deep_tempdir)


def test_planning_exception_propagates_to_the_caller(fixture: _Fixture) -> None:
    # A marker without a protein ID passes the scaffold and span filters, then fails
    # inside serial planning. Nothing is patched.
    broken = SimpleNamespace(
        scaffold="shared", validation_status="validated", start=220, end=230, query_porf=None, hmm_target="x"
    )
    fixture.markers.append(broken)

    with pytest.raises(AttributeError):
        fixture.refine()


def test_serial_planning_reads_each_scaffold_sequence(tmp_path: Path) -> None:
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

    refined = {boundary.scaffold: boundary for boundary in fixture.refine()}

    assert (refined["repeat"].start, refined["repeat"].end, refined["repeat"].tir_status) == (100, 420, "detected")
    assert (refined["plain"].start, refined["plain"].end, refined["plain"].tir_status) == (50, 470, "not_detected")
    assert (refined["plain2"].start, refined["plain2"].end) == (50, 470)


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
    _serial_adjacent_sequence, _ = _adjacent_repeat_sequence()
    length = len(_serial_adjacent_sequence)
    names = ["cand_viral", "cand_weak", "res_viral", "res_weak", "alias", "dup", "res_dup", "orphan"]
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text("".join(f">{name}\n{_serial_adjacent_sequence}\n" for name in names), encoding="utf-8")
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


def test_serial_taxonomy_decides_candidates_and_residuals(tmp_path: Path) -> None:
    fixture = _taxonomy_fixture(tmp_path)

    serial = fixture.refine()

    spans: dict[str, list[tuple[int, int, str]]] = {}
    for boundary in serial:
        spans.setdefault(boundary.scaffold, []).append((boundary.start, boundary.end, boundary.tir_status))
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


def test_spawned_scaffold_plans_match_serial_taxonomy_fixture(tmp_path: Path) -> None:
    """Real repeat planning preserves every boundary field and input order."""
    fixture = _taxonomy_fixture(tmp_path)
    assert fixture.refine(workers=2) == fixture.refine(workers=1)
