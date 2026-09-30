"""The run command resolves its native tools once before resources or genomes are touched."""

from __future__ import annotations

import tempfile
from collections.abc import Iterator
from pathlib import Path

import pytest
from click.testing import CliRunner

from virosync.cli.main import cli
from virosync.orchestration import cli as orchestration_cli
from virosync.pipeline.phase2 import terminal_repeat_search


@pytest.fixture
def missing_runtime(monkeypatch: pytest.MonkeyPatch) -> Iterator[Path]:
    """Point the native runtime root at an empty directory short enough for the prefix limit."""
    with tempfile.TemporaryDirectory(prefix="vs-preflight-", dir="/tmp") as directory:
        root = Path(directory) / "runtime"
        monkeypatch.setenv("VIROSYNC_PRODIGAL_RUNTIME", str(root))
        yield root


@pytest.fixture
def genome_and_output(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    """Provide one genome and guard the steps that must not run after a failed preflight."""
    genome = tmp_path / "genome.fna"
    genome.write_text(">scaffold_1\nATGAAATAA\n")
    monkeypatch.setattr(
        orchestration_cli,
        "_resolve_pipeline_resources",
        lambda *_args, **_kwargs: pytest.fail("resources were resolved after a failed preflight"),
    )
    monkeypatch.setattr(
        orchestration_cli,
        "run_batch_python",
        lambda **_kwargs: pytest.fail("genomes were dispatched after a failed preflight"),
    )
    return genome, tmp_path / "results"


def test_run_stops_with_setup_hint_when_gene_caller_is_missing(
    missing_runtime: Path, genome_and_output: tuple[Path, Path]
) -> None:
    genome, output = genome_and_output

    result = CliRunner().invoke(cli, ["run", "-i", str(genome), "-o", str(output)])

    assert result.exit_code == 1
    assert "Corrected Prodigal-GV is not installed" in result.output
    assert "python -m virosync.utils.prodigal_runtime setup" in result.output
    assert not output.exists()
    assert not missing_runtime.exists()


def test_run_stops_with_clear_message_when_genometools_is_missing(
    genome_and_output: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    genome, output = genome_and_output
    monkeypatch.setattr(orchestration_cli, "resolve_prodigal_executable", lambda: tmp_path / "prodigal-gv")
    monkeypatch.setattr(terminal_repeat_search.shutil, "which", lambda _name: None)

    result = CliRunner().invoke(cli, ["run", "-i", str(genome), "-o", str(output)])

    assert result.exit_code == 1
    assert "GenomeTools executable 'gt' was not found on PATH" in result.output
    assert not output.exists()


@pytest.mark.parametrize(
    "arguments",
    [
        ["--version"],
        ["run", "--help"],
        ["orchestrate", "resources", "--help"],
        ["orchestrate", "info"],
    ],
)
def test_help_version_and_resource_commands_work_without_native_tools(
    missing_runtime: Path, monkeypatch: pytest.MonkeyPatch, arguments: list[str]
) -> None:
    monkeypatch.setattr(
        orchestration_cli,
        "resolve_prodigal_executable",
        lambda: pytest.fail("native preflight ran for a command that does not analyse genomes"),
    )

    result = CliRunner().invoke(cli, arguments)

    assert result.exit_code == 0, result.output
