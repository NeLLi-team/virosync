"""Run the corrected Prodigal-GV on inputs that broke the previous caller or its validator.

These tests need the installed corrected runtime. They skip only when it is
absent; the local native gate runs them.
"""

from __future__ import annotations

import functools
import subprocess
from collections.abc import Iterator
from pathlib import Path

import pytest
from Bio import SeqIO

from virosync.pipeline.phase0 import prodigal
from virosync.utils import prodigal_runtime
from virosync.utils.prodigal_runtime import resolve_prodigal_executable

_EXAMPLE_GENOME = Path(__file__).resolve().parents[1] / "example/test-1.fna"
# The corrected caller finishes the 400 kb record in about 10 s on a workstation.
_CALLER_TIMEOUT_S = 300
_START_DENSE_RECORD_BP = 400_002


@pytest.fixture(scope="module")
def corrected_runtime() -> Path:
    """Select the verified corrected caller; skip only when no installation exists.

    A present but corrupted or incomplete installation fails instead of skipping.
    """
    runtime = prodigal_runtime._runtime_directory(prodigal_runtime._load_recipe())
    if not runtime.exists():
        pytest.skip(f"corrected Prodigal-GV runtime is not installed at {runtime}")
    return resolve_prodigal_executable()


@pytest.fixture
def bounded_caller(corrected_runtime: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[Path]:
    """Give every caller subprocess a timeout so a regression cannot hang the gate."""
    monkeypatch.setattr(prodigal.subprocess, "run", functools.partial(subprocess.run, timeout=_CALLER_TIMEOUT_S))
    yield corrected_runtime


def test_start_dense_record_beyond_the_old_node_capacity_is_called(tmp_path: Path, bounded_caller: Path) -> None:
    """An (ATG)n record of 400 kb aborted the unpatched caller with a heap overflow in add_nodes."""
    sequence = "ATG" * (_START_DENSE_RECORD_BP // 3)
    genome = tmp_path / "start_dense.fna"
    genome.write_text(">start_dense\n" + "\n".join(sequence[i : i + 60] for i in range(0, len(sequence), 60)) + "\n")

    proteins, genes = prodigal.run_prodigal_genome(genome, tmp_path / "phase0", threads=1)

    assert genes
    assert proteins.stat().st_size > 0
    assert not (tmp_path / "prodigal_diagnostics").exists()


def test_crlf_input_gives_the_same_predictions_as_lf(tmp_path: Path, bounded_caller: Path) -> None:
    """The single-process path accepts a CRLF genome and reports identical proteins and coordinates."""
    record = next(SeqIO.parse(_EXAMPLE_GENOME, "fasta"))
    sequence = str(record.seq[:20_000])
    lines = [">crlf_probe description", *(sequence[i : i + 60] for i in range(0, len(sequence), 60))]
    lf_genome = tmp_path / "lf.fna"
    crlf_genome = tmp_path / "crlf.fna"
    lf_genome.write_bytes(("\n".join(lines) + "\n").encode("ascii"))
    crlf_genome.write_bytes(("\r\n".join(lines) + "\r\n").encode("ascii"))

    _, lf_genes = prodigal.run_prodigal_genome(lf_genome, tmp_path / "lf_phase0", threads=1)
    _, crlf_genes = prodigal.run_prodigal_genome(crlf_genome, tmp_path / "crlf_phase0", threads=1)

    assert lf_genes
    assert crlf_genes == lf_genes
    assert not (tmp_path / "prodigal_diagnostics").exists()
