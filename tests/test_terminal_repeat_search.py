"""Tests for bounded native exact k-mer occurrence discovery."""

from __future__ import annotations

from itertools import product

import pytest

from virosync.pipeline.phase2 import terminal_repeat_search
from virosync.pipeline.phase2.terminal_repeat_search import find_seed_positions

_DNA_BASES = frozenset("ACGT")
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def _reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a normalized DNA sequence."""
    return sequence.translate(_COMPLEMENT)[::-1]


def _normalized(sequence: str) -> str:
    """Match the ambiguity semantics of the native adapter."""
    return "".join(base if base in _DNA_BASES else "N" for base in sequence.upper())


def _brute_positions(
    sequence: str,
    seed_length: int,
    kmers: tuple[str, ...],
) -> dict[str, tuple[int, ...]]:
    """Return exact selected k-mer positions by direct comparison."""
    normalized = _normalized(sequence)
    positions: dict[str, tuple[int, ...]] = {}
    for kmer in kmers:
        starts = tuple(
            start
            for start in range(len(normalized) - seed_length + 1)
            if normalized[start : start + seed_length] == kmer
        )
        if starts:
            positions[kmer] = starts
    return positions


def test_native_positions_equal_brute_force_with_ambiguity() -> None:
    sequence = "ttACGTCAGTNNryACTGACGTaa"
    seed_length = 4
    kmers = tuple("".join(bases) for bases in product("ACGT", repeat=seed_length))

    positions = find_seed_positions(sequence, seed_length, kmers)

    assert positions == _brute_positions(sequence, seed_length, kmers)


def test_native_positions_equal_brute_force_for_full_palindrome() -> None:
    left = (
        "GTAGTTCGAAAACTAAGCTAATGGCAATAACTGCCGGCGTAAATACTTAATGCTGTATGCTCGAAGTGGGCACCGCTAATATAGGGCGTTCAGACT"
        "AGCTCGTCATAACAGTGGTGAAGCAACAGTCATACGTAGAAAGACCACGTCTCGAGCATACTCTGATCGGAGTTCGAGCTATCCCGGGAGACCCC"
        "CTCATTGTTGTTGGTTGATTAACAATCGGAATACTTGCTAACTGTTAGGTCGGAGCGGGGCAAGCCATCGACCCTGACATTGCTCCGGGCAAATTC"
        "TCAGCTGACATAC"
    )
    sequence = left + _reverse_complement(left)
    seed_length = 7
    kmers = tuple(sorted({sequence[start : start + seed_length] for start in range(len(sequence) - seed_length + 1)}))

    positions = find_seed_positions(sequence, seed_length, kmers)

    assert len(sequence) == 600
    assert positions == _brute_positions(sequence, seed_length, kmers)


def test_native_positions_remain_linear_for_frequent_motif() -> None:
    motif = "ACGTTGC"
    copies = 800
    sequence = motif * copies

    positions = find_seed_positions(sequence, len(motif), (motif,))

    assert positions == {motif: tuple(range(0, len(sequence), len(motif)))}


@pytest.mark.parametrize("seed_length", [0, -1])
def test_seed_length_below_one_is_rejected_before_any_command(
    monkeypatch: pytest.MonkeyPatch, seed_length: int
) -> None:
    monkeypatch.setattr(
        terminal_repeat_search.subprocess,
        "run",
        lambda *_args, **_kwargs: pytest.fail("gt must not run for an invalid seed length"),
    )

    with pytest.raises(ValueError, match="seed_length must be >= 1"):
        find_seed_positions("ACGTACGTACGT", seed_length, ("ACGT",))


def test_missing_genometools_raises_a_clear_error(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(terminal_repeat_search.shutil, "which", lambda _name: None)

    with pytest.raises(RuntimeError, match="GenomeTools executable 'gt' was not found on PATH"):
        find_seed_positions("ACGTACGTACGT", 4, ("ACGT",))


def test_genometools_failure_reports_its_stderr() -> None:
    with pytest.raises(RuntimeError, match="gt suffixerator failed with exit status 1: .*sequence must not be empty"):
        find_seed_positions("", 7, ("ACGTACG",))
