from __future__ import annotations

import random

import pytest

from virosync.pipeline.phase2 import terminal_repeats


def _naive_best_group_segment(
    matches: list[bool],
    group_start: int,
    group_end: int,
) -> tuple[int, int, float] | None:
    """Run the frozen quadratic scorer used before the optimized sweep."""
    mismatch_prefix = [0]
    for is_match in matches:
        mismatch_prefix.append(mismatch_prefix[-1] + (not is_match))

    last_start = min(
        group_end,
        len(matches) - terminal_repeats.MIN_TIR_ARM_BP,
    )
    best: tuple[int, int, float] | None = None
    best_key: tuple[int, float, int, int] | None = None
    for start in range(last_start + 1):
        max_length = len(matches) - start
        for length in range(terminal_repeats.MIN_TIR_ARM_BP, max_length + 1):
            contained_window_start = max(group_start, start)
            contained_window_end = min(
                group_end,
                start + length - terminal_repeats.MIN_TIR_ARM_BP,
            )
            if contained_window_start > contained_window_end:
                continue
            if not matches[start] or not matches[start + length - 1]:
                continue
            mismatches = mismatch_prefix[start + length] - mismatch_prefix[start]
            if not terminal_repeats._passes_identity(mismatches, length):
                continue
            match_count = length - mismatches
            identity = match_count / length
            candidate = (start, length, identity)
            score = match_count * terminal_repeats._MATCH_SCORE + mismatches * terminal_repeats._MISMATCH_SCORE
            candidate_key = (score, identity, length, -start)
            if best_key is None or candidate_key > best_key:
                best = candidate
                best_key = candidate_key
    return best


def _assert_group_parity(matches: list[bool]) -> int:
    """Compare optimized and naive results for every qualifying group."""
    groups = terminal_repeats._valid_window_groups(matches)
    for group_start, group_end in groups:
        expected = _naive_best_group_segment(matches, group_start, group_end)

        observed = terminal_repeats._best_group_segment(
            matches,
            group_start,
            group_end,
        )

        assert observed == expected
    return len(groups)


def test_best_group_segment_matches_exhaustive_bounded_patterns() -> None:
    variable_positions = (0, 1, 6, 11, 16, 21, 31, 40, 49, 50, 56, 61)
    sequence_length = terminal_repeats.MIN_TIR_ARM_BP + 12
    group_count = 0
    for mismatch_mask in range(1 << len(variable_positions)):
        matches = [True] * sequence_length
        for bit, position in enumerate(variable_positions):
            if mismatch_mask & (1 << bit):
                matches[position] = False
        group_count += _assert_group_parity(matches)

    assert group_count == 3819


def test_best_group_segment_matches_randomized_naive_oracle() -> None:
    rng = random.Random(20260919)
    group_count = 0
    for _ in range(400):
        sequence_length = rng.randint(terminal_repeats.MIN_TIR_ARM_BP, 350)
        match_probability = rng.uniform(0.86, 0.99)
        matches = [rng.random() < match_probability for _ in range(sequence_length)]
        group_count += _assert_group_parity(matches)

    assert group_count == 438


@pytest.mark.parametrize(
    ("mismatch_count", "expected"),
    [
        (5, (0, 50, 0.9)),
        (6, None),
    ],
)
def test_best_group_segment_enforces_exact_identity_threshold(
    mismatch_count: int,
    expected: tuple[int, int, float] | None,
) -> None:
    matches = [True] * terminal_repeats.MIN_TIR_ARM_BP
    for position in range(5, 5 + mismatch_count):
        matches[position] = False

    assert terminal_repeats._best_group_segment(matches, 0, 0) == expected


def test_best_group_segment_requires_matching_outer_bases() -> None:
    matches = [False, *([True] * terminal_repeats.MIN_TIR_ARM_BP), False]

    observed = terminal_repeats._best_group_segment(matches, 1, 1)

    assert observed == (1, terminal_repeats.MIN_TIR_ARM_BP, 1.0)


def test_best_group_segment_must_contain_a_window_from_its_group() -> None:
    matches = [
        *([True] * 60),
        *([False] * 60),
        *([True] * terminal_repeats.MIN_TIR_ARM_BP),
    ]

    observed = terminal_repeats._best_group_segment(matches, 120, 120)

    assert observed == (120, terminal_repeats.MIN_TIR_ARM_BP, 1.0)


def test_best_group_segment_score_tie_prefers_higher_identity() -> None:
    matches = [
        base == "1"
        for base in (
            "1111110110111111111111111100111111101111111111111111010111111111111111"
            "1111101101111111101011110111111111111010110011111111111111011110101111"
            "1"
        )
    ]

    observed = terminal_repeats._best_group_segment(matches, 0, 57)

    assert observed == (10, 77, 70 / 77)


def test_best_group_segment_full_tie_prefers_earliest_start() -> None:
    matches = [
        base == "1"
        for base in (
            "1111111010011111110011011111111110111010111111111111110111001111010111"
            "1101111110010111111110111101111111111101111101111111111111111010101101"
            "1111011"
        )
    ]

    observed = terminal_repeats._best_group_segment(matches, 81, 83)

    assert observed == (81, 50, 0.9)
