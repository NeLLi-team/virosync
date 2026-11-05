"""Trim EVE region edges that look like host by GC content (experimental V9e).

The central half of a region sets the EVE GC and up to 50 kb on each side set the host GC.
When the two differ by at least 3 percentage points, in either direction, each edge moves
inward in 2 kb windows while a window is closer to the host GC than to the EVE GC. The trim
stops at the first EVE-like window and never enters the confirmed seed cores. Regions with a
smaller contrast keep their span, because 2 kb windows cannot resolve it.
"""

from collections.abc import Iterable
from pathlib import Path

WINDOW_BP = 2_000
FLANK_BP = 50_000
MIN_CONTRAST = 0.03


def trim_boundaries_by_gc(boundaries: Iterable[object], genome_path: Path) -> list[object]:
    """Return the boundaries with host-like edges trimmed in place."""
    from Bio import SeqIO  # noqa: PLC0415  # Matches the Phase 2 composition step.

    trimmed = list(boundaries)
    index = SeqIO.index(str(genome_path), "fasta")
    try:
        for boundary in trimmed:
            if boundary.scaffold not in index:
                continue
            cores = boundary.confirmed_seeds
            core = (min(c.start for c in cores), max(c.end for c in cores)) if cores else None
            sequence = str(index[boundary.scaffold].seq).upper()
            boundary.start, boundary.end = trimmed_span(sequence, boundary.start, boundary.end, core=core)
    finally:
        index.close()
    return trimmed


def trimmed_span(sequence: str, start: int, end: int, *, core: tuple[int, int] | None) -> tuple[int, int]:
    """Return start-end after trimming host-like 2 kb windows from both edges."""
    quarter = (end - start) // 4
    eve_gc = _gc(sequence[start + quarter : end - quarter])
    host_gc = _gc(sequence[max(0, start - FLANK_BP) : start] + sequence[end : end + FLANK_BP])
    if eve_gc is None or host_gc is None or abs(eve_gc - host_gc) < MIN_CONTRAST:
        return start, end
    core_start, core_end = core or (end, start)

    def host_like(left: int, right: int) -> bool:
        window = _gc(sequence[left:right])
        return window is not None and abs(window - host_gc) < abs(window - eve_gc)

    while end - start > 2 * WINDOW_BP and start + WINDOW_BP <= core_start and host_like(start, start + WINDOW_BP):
        start += WINDOW_BP
    while end - start > 2 * WINDOW_BP and end - WINDOW_BP >= core_end and host_like(end - WINDOW_BP, end):
        end -= WINDOW_BP
    return start, end


def _gc(sequence: str) -> float | None:
    """Return the GC fraction over A, C, G and T, or None without such bases."""
    gc = sequence.count("G") + sequence.count("C")
    total = gc + sequence.count("A") + sequence.count("T")
    return gc / total if total else None
