"""Tests for the conditional GC edge trim of V9g (virosync.pipeline.phase2.gc_boundary)."""

import random

from virosync.pipeline.phase2 import gc_boundary


def sequence(gc: float, length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("GC") if rng.random() < gc else rng.choice("AT") for _ in range(length))


HOST, EVE = 0.62, 0.40


def genome(eve_gc: float) -> str:
    """Return 60 kb host, a 20 kb EVE at 60,000-80,000, then 60 kb host."""
    return sequence(HOST, 60_000, 1) + sequence(eve_gc, 20_000, 2) + sequence(HOST, 60_000, 3)


def test_host_like_overextension_is_trimmed_to_the_eve() -> None:
    start, end = gc_boundary.trimmed_span(genome(EVE), 52_000, 90_000, core=(66_000, 74_000))

    assert (start, end) == (60_000, 80_000)


def test_low_contrast_region_is_left_unchanged() -> None:
    assert gc_boundary.trimmed_span(genome(0.61), 52_000, 90_000, core=(66_000, 74_000)) == (52_000, 90_000)


def test_trim_never_enters_the_seed_core() -> None:
    start, end = gc_boundary.trimmed_span(genome(EVE), 52_000, 90_000, core=(56_000, 74_000))

    assert start <= 56_000
