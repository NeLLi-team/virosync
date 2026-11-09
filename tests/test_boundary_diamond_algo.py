"""Boundary searches pin query-indexed DIAMOND for shared taxonomy and viral evidence."""

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


@pytest.mark.parametrize("chunk_size, expected_calls", [(100, 1), (1, 2)], ids=["single-call", "chunked"])
@pytest.mark.parametrize("top_k, viral_hit_depth, expected_depth", [(10, 50, 50), (60, 50, 60)])
def test_batched_boundary_search_pins_algo_for_shared_hits(
    tmp_path: Path,
    recorded_argv: Path,
    chunk_size: int,
    expected_calls: int,
    top_k: int,
    viral_hit_depth: int,
    expected_depth: int,
) -> None:
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
        config=bd.BoundaryDiamondConfig(chunk_size=chunk_size, threads=1, top_k=top_k, viral_hit_depth=viral_hit_depth),
    )

    calls = _calls(recorded_argv)
    assert len(calls) == expected_calls
    assert all("--algo 1" in call for call in calls)
    assert all(f"--max-target-seqs {expected_depth}" in call for call in calls)


def test_full_proteome_search_does_not_pin_algo(tmp_path: Path, recorded_argv: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">p1\nMPEPTIDE\n", encoding="utf-8")

    bd.run_full_proteome_diamond(proteome, tmp_path / "db.dmnd", tmp_path / "out", max_target_seqs=10, threads=1)

    calls = _calls(recorded_argv)
    assert len(calls) == 1
    assert "--algo" not in calls[0]


def test_boundary_search_default_batch_size() -> None:
    """Direct boundary callers use the same batch size as pipeline defaults."""
    assert bd.BoundaryDiamondConfig().chunk_size == 100000


def test_boundary_algo_survives_reduced_thread_retry(
    tmp_path: Path, recorded_argv: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A timed-out boundary call keeps query-indexed selection on its retry."""
    from virosync.pipeline import search_backend

    script = tmp_path / "bin" / "diamond"
    first = tmp_path / "first-attempt"
    script.write_text(
        script.read_text().replace(
            'while [ "$1" != "--out" ];',
            f"if [ ! -e {shlex.quote(str(first))} ]; then touch {shlex.quote(str(first))}; sleep 2; fi\n"
            'while [ "$1" != "--out" ];',
        )
    )
    monkeypatch.setattr(search_backend, "DIAMOND_DEFAULT_TIMEOUT", 0.2)
    monkeypatch.setattr(search_backend, "DIAMOND_RETRY_THREADS", 4)
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">p1\nMPEPTIDE\n", encoding="utf-8")
    bd.run_batched_diamond(
        query=bd.GenomeDiamondQuery(all_porf_ids=["p1"]),
        proteome_fasta=proteome,
        diamond_db=tmp_path / "db.dmnd",
        output_dir=tmp_path / "out",
        proteome_index={"contig": [bd.pORF("p1", "contig", 0, 90)]},
        config=bd.BoundaryDiamondConfig(threads=8),
    )
    first_call, retry_call = _calls(recorded_argv)
    assert "--algo 1" in first_call
    assert "--algo 1" in retry_call
    assert "--threads 8" in first_call
    assert "--threads 4" in retry_call
