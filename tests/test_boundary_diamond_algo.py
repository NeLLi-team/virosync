"""Boundary searches preserve DIAMOND's automatic algorithm selection."""

import os
import shlex
import stat
from pathlib import Path

import pytest

from virosync.pipeline.phase2 import boundary_diamond as bd


@pytest.fixture
def recorded_argv(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Put a fake `diamond` on PATH that appends its argv to a file and writes an empty --out."""
    log = tmp_path / "argv.log"
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    script = bin_dir / "diamond"
    script.write_text(
        f'#!/bin/sh\necho "$*" >> {shlex.quote(str(log))}\nwhile [ "$1" != "--out" ]; do shift; done\nshift\n: > "$1"\n'
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    return log


def _calls(log: Path) -> list[str]:
    return log.read_text().splitlines()


@pytest.mark.parametrize("chunk_size", [100, 1], ids=["single-call", "chunked"])
def test_batched_boundary_search_does_not_pin_algo(tmp_path: Path, recorded_argv: Path, chunk_size: int) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">p1 # 1 # 90 # + # ID=contig_1;\nMPEPTIDE\n>p2 # 101 # 190 # + # ID=contig_2;\nMPEPTIDE\n",
        encoding="utf-8",
    )
    query = bd.GenomeDiamondQuery(eve_porf_ids={"seed": ["p1", "p2"]}, all_porf_ids=["p1", "p2"])

    bd.run_batched_diamond(
        query=query,
        proteome_fasta=proteome,
        diamond_db=tmp_path / "db.dmnd",
        output_dir=tmp_path / "out",
        proteome_index={"contig": [bd.pORF("p1", "contig", 0, 90), bd.pORF("p2", "contig", 100, 190)]},
        config=bd.BoundaryDiamondConfig(chunk_size=chunk_size, threads=1),
    )

    calls = _calls(recorded_argv)
    assert len(calls) == (2 if chunk_size == 1 else 1)
    assert all("--algo" not in call for call in calls)


def test_full_proteome_search_does_not_pin_algo(tmp_path: Path, recorded_argv: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">p1\nMPEPTIDE\n", encoding="utf-8")

    bd.run_full_proteome_diamond(proteome, tmp_path / "db.dmnd", tmp_path / "out", max_target_seqs=10, threads=1)

    calls = _calls(recorded_argv)
    assert len(calls) == 1
    assert "--algo" not in calls[0]
