"""Index selected exact k-mer occurrences with GenomeTools."""

from __future__ import annotations

import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

__all__ = ["find_seed_positions"]

_DNA_BASES = frozenset("ACGT")


def find_seed_positions(
    sequence: str,
    seed_length: int,
    kmers: tuple[str, ...],
) -> dict[str, tuple[int, ...]]:
    """Return sorted region positions for each selected k-mer that occurs."""
    native_sequence = _normalize_sequence(sequence)
    with tempfile.TemporaryDirectory(prefix="virosync-repfind-") as directory:
        work_dir = Path(directory)
        input_path = work_dir / "input.fna"
        query_path = work_dir / "queries.fna"
        index_path = work_dir / "sequence"
        output_path = work_dir / "matches.tsv"
        input_path.write_text(f">sequence\n{native_sequence}\n", encoding="utf-8")
        _write_queries(query_path, kmers)
        _build_suffix_index(input_path, index_path)
        _run_repfind(index_path, query_path, output_path, seed_length)
        return _read_positions(output_path, kmers)


def _normalize_sequence(sequence: str) -> str:
    """Map ambiguity symbols to the nonmatching GenomeTools symbol."""
    return "".join(base if base in _DNA_BASES else "N" for base in sequence.upper())


def _write_queries(query_path: Path, kmers: tuple[str, ...]) -> None:
    """Write unique selected k-mers in query-number order."""
    with query_path.open("w", encoding="utf-8") as queries:
        for query_number, kmer in enumerate(kmers):
            queries.write(f">{query_number}\n{kmer}\n")


def _build_suffix_index(input_path: Path, index_path: Path) -> None:
    """Build the enhanced suffix array required by repfind."""
    subprocess.run(
        [
            "gt",
            "suffixerator",
            "-db",
            str(input_path),
            "-indexname",
            str(index_path),
            "-dna",
            "-suf",
            "-lcp",
            "-des",
            "-ssp",
            "-sds",
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def _run_repfind(
    index_path: Path,
    query_path: Path,
    output_path: Path,
    seed_length: int,
) -> None:
    """Write exact forward query occurrences to a temporary file."""
    with output_path.open("w", encoding="utf-8") as output:
        subprocess.run(
            [
                "gt",
                "repfind",
                "-ii",
                str(index_path),
                "-q",
                str(query_path),
                "-l",
                str(seed_length),
                "-f",
                "yes",
                "-outfmt",
                "custom",
                "s.start",
                "q.seqnum",
            ],
            check=True,
            stdout=output,
            stderr=subprocess.PIPE,
            text=True,
        )


def _read_positions(
    output_path: Path,
    kmers: tuple[str, ...],
) -> dict[str, tuple[int, ...]]:
    """Map repfind query numbers back to requested k-mers."""
    positions_by_query: dict[int, list[int]] = defaultdict(list)
    with output_path.open(encoding="utf-8") as rows:
        for row in rows:
            if row.startswith("#") or not row.strip():
                continue
            subject_start, query_number = map(int, row.split())
            positions_by_query[query_number].append(subject_start)
    return {
        kmer: tuple(sorted(positions_by_query[query_number]))
        for query_number, kmer in enumerate(kmers)
        if query_number in positions_by_query
    }
