from __future__ import annotations

from dataclasses import replace
import subprocess
from types import SimpleNamespace

import pytest

from scripts.ci import check_production_ready
from virosync.report import graphviz_runtime


def test_resource_release_surfaces_parse_and_cross_check() -> None:
    failures: list[str] = []
    check_production_ready.check_resource_version(failures)
    assert failures == []


def test_resource_release_guard_rejects_database_source_drift(monkeypatch) -> None:
    source = dict(
        check_production_ready.ViroSyncDatabaseManager.DATABASE_SOURCES[0]
    )
    source["manifest_sha256"] = "0" * 64
    monkeypatch.setattr(
        check_production_ready.ViroSyncDatabaseManager,
        "DATABASE_SOURCES",
        [source],
    )
    failures: list[str] = []
    check_production_ready.check_resource_version(failures)
    assert any("DATABASE_SOURCES first source differs" in item for item in failures)


def test_resource_release_guard_rejects_source_manifest(monkeypatch) -> None:
    manifest = check_production_ready.load_resource_manifest(
        check_production_ready.ROOT / check_production_ready.RELEASE_MANIFEST_PATH,
        expected_version=check_production_ready.DATABASE_VERSION,
        expected_manifest_sha256=check_production_ready.RESOURCE_MANIFEST_SHA256,
    )
    monkeypatch.setattr(
        check_production_ready,
        "load_resource_manifest",
        lambda *_args, **_kwargs: replace(manifest, bundle_kind="source"),
    )
    failures: list[str] = []
    check_production_ready.check_resource_version(failures)
    assert "tracked release manifest must describe a runtime bundle" in failures


def test_complete_production_guard_passes_tracked_release_surface() -> None:
    # The benchmark harness moved to the virosync-bench repo, so the guard no
    # longer needs a stubbed benchmark surface to pass.
    assert check_production_ready.main() == 0


def test_production_guard_reports_graphviz_render_failure(monkeypatch) -> None:
    monkeypatch.setattr(
        check_production_ready,
        "graphviz_runtime_error",
        lambda: "Graphviz cannot render the required sfdp PNG report: missing PNG plugin",
    )
    failures: list[str] = []

    check_production_ready.check_graphviz_runtime(failures)

    assert failures == [
        "Graphviz cannot render the required sfdp PNG report: missing PNG plugin"
    ]


@pytest.mark.parametrize(
    ("result", "expected"),
    [
        (
            subprocess.CompletedProcess(
                ["dot"], 1, stdout=b"", stderr=b"missing PNG plugin"
            ),
            "Graphviz cannot render the required sfdp PNG report: missing PNG plugin",
        ),
        (
            subprocess.CompletedProcess(
                ["dot"], 0, stdout=b"not a PNG", stderr=b""
            ),
            "Graphviz cannot render the required sfdp PNG report",
        ),
    ],
)
def test_graphviz_runtime_reports_bad_process_results(
    monkeypatch,
    result: subprocess.CompletedProcess,
    expected: str,
) -> None:
    stub = SimpleNamespace(
        run=lambda *args, **kwargs: result,
        TimeoutExpired=subprocess.TimeoutExpired,
    )
    monkeypatch.setattr(graphviz_runtime, "subprocess", stub)

    assert graphviz_runtime.graphviz_runtime_error() == expected


def test_graphviz_runtime_reports_missing_executable(monkeypatch) -> None:
    def missing(*args, **kwargs):
        raise FileNotFoundError("dot is missing")

    stub = SimpleNamespace(
        run=missing,
        TimeoutExpired=subprocess.TimeoutExpired,
    )
    monkeypatch.setattr(graphviz_runtime, "subprocess", stub)

    assert graphviz_runtime.graphviz_runtime_error() == (
        "Graphviz runtime is unavailable: dot is missing"
    )


def test_graphviz_runtime_times_out(monkeypatch) -> None:
    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(args[0], kwargs["timeout"])

    stub = SimpleNamespace(
        run=timeout,
        TimeoutExpired=subprocess.TimeoutExpired,
    )
    monkeypatch.setattr(graphviz_runtime, "subprocess", stub)

    assert graphviz_runtime.graphviz_runtime_error(timeout_seconds=2) == (
        "Graphviz sfdp PNG render timed out after 2 seconds"
    )


@pytest.mark.parametrize(
    "path",
    [
        "benchmarking/harness/results/private.tsv",
        "benchmarking/v2/.pixi/state.json",
        "benchmarking/v2/results/run/native/output.tsv",
        "benchmarking/v2/manifests/sealed/truth.tsv",
        "benchmarking/v2/publication/artifacts/full.tar.gz",
        "benchmarking/v2/data/private.fna",
        "benchmarking/v2/genome.fasta",
        "benchmarking/v2/contigs.fa",
        "benchmarking/v2/big.csv",
        "benchmarking/v2/matrix.parquet",
        "benchmarking/v2/notes.ipynb",
        "benchmarking/v2/extra/foo.tsv",
        "benchmarking/v2/config/.env.local",
        "docs/ms/virosync_methods_benchmark_manuscript.qmd",
        "docs/ms/v2/manuscript.pdf",
        "docs/ms/v2/result.png",
        "docs/ms/v2/figures/result.png",
        "docs/ms/v2/.memd/data/store.bin",
        "docs/ms/v2/memory.md",
    ],
)
def test_generated_private_and_legacy_benchmark_paths_are_forbidden(path: str) -> None:
    assert check_production_ready.is_forbidden_tracked_path(path) is True
