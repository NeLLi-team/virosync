"""Exactness checks for vectorized marker-interval merging and deduplication."""

from __future__ import annotations

import numpy as np
import pytest
from numpy.typing import NDArray

from virosync.pipeline.phase2 import terminal_repeats as repeats


def _scalar_oracle(
    pairs: NDArray[np.int64],
    left_length: int,
    right_length: int,
) -> list[tuple[int, int, int]]:
    """Return intervals with the frozen scalar merge implementation."""
    earlier, later = pairs[:-1], pairs[1:]
    diagonal = later[:, 0]
    starts = np.maximum(
        np.maximum(0, -diagonal),
        later[:, 1] + repeats._KMER_BP - repeats.MAX_TIR_ARM_BP,
    )
    ends = np.minimum(
        np.minimum(left_length, right_length - diagonal),
        earlier[:, 1] + repeats.MAX_TIR_ARM_BP,
    )
    supported = (
        (earlier[:, 0] == diagonal)
        & (later[:, 1] - earlier[:, 1] <= repeats.MAX_TIR_ARM_BP - repeats._KMER_BP)
        & (ends - starts >= repeats.MIN_TIR_ARM_BP)
    )
    intervals: list[tuple[int, int, int]] = []
    for current_diagonal, start, end in np.column_stack(
        (diagonal[supported], starts[supported], ends[supported])
    ).tolist():
        if intervals and current_diagonal == intervals[-1][0] and start <= intervals[-1][2]:
            previous = intervals[-1]
            intervals[-1] = (
                current_diagonal,
                previous[1],
                max(previous[2], end),
            )
        else:
            intervals.append((current_diagonal, start, end))
    return intervals


def _pairs(*rows: tuple[int, int]) -> NDArray[np.int64]:
    """Build sorted diagonal-position pairs with one seed ID."""
    return np.asarray([(diagonal, position, 0) for diagonal, position in rows], dtype=np.int64).reshape(-1, 3)


@pytest.mark.parametrize(
    ("pairs", "left_length", "right_length", "expected"),
    [
        (_pairs(), 600, 600, []),
        (_pairs((0, 0)), 600, 600, []),
        (
            _pairs((0, 0), (0, 100), (0, 900), (0, 993)),
            1_500,
            1_500,
            [(0, 0, 1_400)],
        ),
        (
            _pairs((0, 0), (0, 100), (0, 900), (0, 994)),
            1_500,
            1_500,
            [(0, 0, 500), (0, 501, 1_400)],
        ),
        (
            _pairs((-200, 200), (-200, 250), (150, 0), (150, 80)),
            700,
            650,
            [(-200, 200, 700), (150, 0, 500)],
        ),
        (_pairs((0, 0), (0, 494), (0, 988)), 1_500, 1_500, []),
    ],
)
def test_marker_intervals_match_deterministic_oracle(
    pairs: NDArray[np.int64],
    left_length: int,
    right_length: int,
    expected: list[tuple[int, int, int]],
) -> None:
    """Preserve empty, clipped, touching, split, diagonal, and support cases."""
    assert _scalar_oracle(pairs, left_length, right_length) == expected
    actual = repeats._marker_search_intervals(pairs, left_length, right_length)

    assert actual.dtype == np.int64
    assert actual.shape == (len(expected), 3)
    assert [tuple(interval) for interval in actual.tolist()] == expected


def test_marker_intervals_match_seeded_random_oracle() -> None:
    """Match the scalar merge on bounded eligible pairs with varied clipping."""
    random = np.random.default_rng(20_260_919)
    for _ in range(500):
        left_length = int(random.integers(50, 2_001))
        right_length = int(random.integers(50, 2_001))
        rows: list[tuple[int, int]] = []
        diagonals = set(random.integers(-900, 901, size=int(random.integers(0, 9))).tolist())
        for diagonal in sorted(diagonals):
            minimum = max(0, -diagonal)
            maximum = min(
                left_length - repeats._KMER_BP,
                right_length - repeats._KMER_BP - diagonal,
            )
            if maximum < minimum:
                continue
            positions = random.integers(
                minimum,
                maximum + 1,
                size=int(random.integers(1, 31)),
            )
            rows.extend((diagonal, int(position)) for position in sorted(positions.tolist()))
        pairs = _pairs(*rows)

        expected = _scalar_oracle(pairs, left_length, right_length)
        actual = repeats._marker_search_intervals(pairs, left_length, right_length)

        assert [tuple(interval) for interval in actual.tolist()] == expected


def test_region_deduplicates_only_exact_marker_intervals(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep unequal overlapping marker intervals while removing exact triples."""
    pairs = _pairs((0, 0), (0, 100), (0, 200), (0, 300))
    allowed = np.ones((3, 1), dtype=np.bool_)
    groups = [
        repeats._MarkerGroup("first", 600, 1_400, frozenset()),
        repeats._MarkerGroup("duplicate", 600, 1_400, frozenset()),
        repeats._MarkerGroup("wider", 700, 1_300, frozenset()),
    ]
    monkeypatch.setattr(repeats, "_informative_seeds", lambda: ("ACGTACG",))
    monkeypatch.setattr(repeats, "find_seed_positions", lambda *_args: {})
    monkeypatch.setattr(
        repeats,
        "_eligible_region_seed_pairs",
        lambda *_args: (pairs, allowed),
    )

    intervals = repeats._region_search_intervals("A" * 2_000, groups, scan_start=0)

    assert intervals == {0: {(0, 600), (0, 700)}}
