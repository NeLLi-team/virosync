"""Differential tests for region-wide terminal-repeat discovery."""

from __future__ import annotations

import random

import pytest
from test_terminal_repeats import (
    _adjacent_repeat_sequence,
    _planted_repeat,
    _random_dna,
    _reverse_complement,
)

from virosync.pipeline.phase2 import terminal_repeats


def _marker_group(
    protein_id: str,
    start: int,
    end: int,
) -> terminal_repeats._MarkerGroup:
    return terminal_repeats._MarkerGroup(
        protein_id=protein_id,
        start=start,
        end=end,
        marker_names=frozenset(),
    )


def _public_candidate_union(
    sequence: str,
    *,
    scan_start: int,
    scan_end: int,
    groups: list[terminal_repeats._MarkerGroup],
) -> list[terminal_repeats.TerminalRepeat]:
    candidates: dict[tuple[int, int], terminal_repeats.TerminalRepeat] = {}
    for group in groups:
        for candidate in terminal_repeats.find_terminal_inverted_repeats(
            sequence,
            scan_start=scan_start,
            scan_end=scan_end,
            marker_start=group.start,
            marker_end=group.end,
        ):
            terminal_repeats._retain_best_outer_pair(candidates, candidate)
    return sorted(
        candidates.values(),
        key=lambda candidate: (
            candidate.left_start,
            candidate.right_end,
            -candidate.arm_length,
            -candidate.identity,
        ),
    )


def _assert_matches_public_union(
    sequence: str,
    *,
    scan_start: int,
    scan_end: int,
    groups: list[terminal_repeats._MarkerGroup],
) -> list[terminal_repeats.TerminalRepeat]:
    expected = _public_candidate_union(
        sequence,
        scan_start=scan_start,
        scan_end=scan_end,
        groups=groups,
    )

    discovered = terminal_repeats._discover_marker_anchored_candidates(
        sequence,
        scan_start=scan_start,
        scan_end=scan_end,
        groups=groups,
    )

    assert discovered == expected
    return discovered


@pytest.mark.parametrize(
    ("seed", "mutated_positions"),
    [
        (17, ()),
        (31, (4, 20, 45)),
        (53, (10, 20, 30, 40, 50, 55)),
    ],
)
def test_region_discovery_matches_public_scans_for_planted_random_arms(
    seed: int,
    mutated_positions: tuple[int, ...],
) -> None:
    sequence, expected_pair = _planted_repeat(
        seed=seed,
        mutated_positions=mutated_positions,
    )

    discovered = _assert_matches_public_union(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        groups=[_marker_group("core", 250, 270)],
    )

    assert [
        (candidate.left_start, candidate.left_end, candidate.right_start, candidate.right_end)
        for candidate in discovered
    ] == [expected_pair]


def test_region_discovery_matches_clipped_pair_when_marker_crosses_repeat_arm() -> None:
    sequence, _ = _planted_repeat(arm_length=120, seed=18)

    discovered = _assert_matches_public_union(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        groups=[_marker_group("arm", 160, 180)],
    )

    assert [(candidate.left_start, candidate.right_end) for candidate in discovered] == [(100, 540)]
    assert discovered[0].alignment_length == 60


def test_region_discovery_retains_longest_evidence_for_shared_pair() -> None:
    sequence, expected_pair = _planted_repeat(arm_length=615, seed=8)

    discovered = _assert_matches_public_union(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        groups=[
            _marker_group("arm", 250, 270),
            _marker_group("core", 795, 815),
        ],
    )

    assert [(candidate.left_start, candidate.right_end) for candidate in discovered] == [
        (expected_pair[0], expected_pair[3])
    ]
    assert discovered[0].alignment_length == 615
    assert discovered[0].alignment_capped is True


@pytest.mark.parametrize(("identical_arms", "expected_count"), [(False, 2), (True, 3)])
def test_region_discovery_matches_public_scans_for_adjacent_pairs(
    identical_arms: bool,
    expected_count: int,
) -> None:
    sequence, _ = _adjacent_repeat_sequence(identical_arms=identical_arms)

    discovered = _assert_matches_public_union(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        groups=[
            _marker_group("first", 220, 240),
            _marker_group("second", 640, 660),
        ],
    )

    assert len(discovered) == expected_count


def test_region_discovery_matches_public_scan_with_nonzero_offset() -> None:
    planted, expected_pair = _planted_repeat(seed=29)
    rng = random.Random(99)
    prefix = _random_dna(150, rng)
    sequence = prefix + planted + _random_dna(130, rng)
    shifted_pair = tuple(coordinate + len(prefix) for coordinate in expected_pair)

    discovered = _assert_matches_public_union(
        sequence,
        scan_start=70,
        scan_end=len(sequence) - 60,
        groups=[_marker_group("shifted", 400, 420)],
    )

    assert [
        (candidate.left_start, candidate.left_end, candidate.right_start, candidate.right_end)
        for candidate in discovered
    ] == [shifted_pair]


@pytest.mark.parametrize(
    ("left_repeat_count", "right_repeat_count", "exceeds_cap"),
    [(64, 64, False), (65, 64, True)],
)
def test_region_discovery_preserves_rare_seeds_at_repetitive_seed_pairing_cap(
    left_repeat_count: int,
    right_repeat_count: int,
    exceeds_cap: bool,
) -> None:
    planted, expected_pair = _planted_repeat(seed=47)
    motif = "ACGTGCA"
    left_flank = motif * left_repeat_count
    right_flank = motif * right_repeat_count
    sequence = left_flank + "TTT" + planted + "CCC" + _reverse_complement(right_flank)
    shift = len(left_flank) + 3
    marker = _marker_group("rare-seeded", shift + 250, shift + 270)
    left_seed_index = terminal_repeats._kmer_index(sequence[: marker.start])
    reverse_right_seed_index = terminal_repeats._kmer_index(_reverse_complement(sequence[marker.end :]))

    seed_pairings = len(left_seed_index[motif]) * len(reverse_right_seed_index[motif])
    assert (seed_pairings > terminal_repeats._MAX_SEED_PAIRINGS) is exceeds_cap

    discovered = _assert_matches_public_union(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        groups=[marker],
    )

    assert [(candidate.left_start, candidate.right_end) for candidate in discovered] == [
        (expected_pair[0] + shift, expected_pair[3] + shift)
    ]


def test_region_seed_pairs_preserve_eligible_edges_without_cross_product() -> None:
    seed = "ACGTCAG"
    index = {
        seed: tuple(range(0, 800, 8)),
        _reverse_complement(seed): tuple(range(1000, 1800, 8)),
    }

    pairs, _ = terminal_repeats._eligible_region_seed_pairs(index, 1800, [(7, 800), (800, 8)])

    assert len(pairs) == 199
    assert set(pairs[:, 1]) == set(range(0, 800, 8))
    assert set(pairs[:, 0] + pairs[:, 1]) == set(range(1, 800, 8))
    assert not ((pairs[:, 1] > 0) & (pairs[:, 0] + pairs[:, 1] > 1)).any()


def test_region_discovery_skips_fragmented_overfrequent_motif_before_pairing() -> None:
    planted, expected_pair = _planted_repeat(seed=47)
    sequence = "ACGTCAGN" * 800 + planted + "NCTGACGT" * 800

    discovered = _assert_matches_public_union(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        groups=[_marker_group("middle", 6650, 6670)],
    )

    assert [(candidate.left_start, candidate.right_end) for candidate in discovered] == [
        (expected_pair[0] + 6400, expected_pair[3] + 6400)
    ]
