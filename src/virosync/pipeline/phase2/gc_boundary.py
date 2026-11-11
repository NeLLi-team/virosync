"""Trim EVE region edges that look like host by GC content (experimental V9e).

The central half of a region sets the EVE GC and up to 50 kb on each side set the host GC.
When the two differ by at least 3 percentage points, in either direction, each edge moves
inward in 2 kb windows while a window is closer to the host GC than to the EVE GC. The trim
stops at the first EVE-like window and never enters the protected span. That span is the hull
of the confirmed seed cores and of every gene that overlaps the region and has a viral match.
Regions with a smaller contrast keep their span, because 2 kb windows cannot resolve it.

A trimmed end that lands strictly inside a gene moves to the nearest edge of a gene that
contains it when an edge lies in the allowed interval. Overlapping genes or the lack of
a valid edge can leave the end inside a gene. The snap stays between the pre-trim end
and the protected span.
"""

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Literal, Protocol

from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy, has_viral_match
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary

GC_TRIM_REASON = "gc_trim"

WINDOW_BP = 2_000
FLANK_BP = 50_000
MIN_CONTRAST = 0.03


class GeneSpan(Protocol):
    """A gene with 0-based half-open coordinates, as in the Phase 2 proteome index."""

    @property
    def id(self) -> str: ...

    @property
    def start(self) -> int: ...

    @property
    def end(self) -> int: ...


def trim_boundaries_by_gc(
    boundaries: Iterable[RefinedBoundary],
    genome_path: Path,
    genes: Mapping[str, Sequence[GeneSpan]],
    taxonomy: Mapping[str, GeneTaxonomy],
) -> list[RefinedBoundary]:
    """Return the boundaries with host-like edges trimmed in place.

    The trim and the snap keep out of the protected span of each boundary: its confirmed
    seed cores and every overlapping gene with a viral match in ``taxonomy``. A trimmed end
    inside a gene snaps to a gene edge (see ``snap_to_gene_edge``). A side that the trim
    moved gets the stop reason ``gc_trim``; ``left_stop`` and ``right_stop`` stay.
    """
    from Bio import SeqIO  # noqa: PLC0415  # Matches the Phase 2 composition step.

    trimmed = list(boundaries)
    index = SeqIO.index(str(genome_path), "fasta")
    try:
        for boundary in trimmed:
            if boundary.scaffold not in index:
                continue
            sequence = str(index[boundary.scaffold].seq).upper()
            _trim_boundary(boundary, sequence, genes.get(boundary.scaffold, ()), taxonomy)
    finally:
        index.close()
    return trimmed


def _protected_span(
    boundary: RefinedBoundary, genes: Sequence[GeneSpan], taxonomy: Mapping[str, GeneTaxonomy]
) -> tuple[int, int] | None:
    """Return the hull of the seed cores and of the overlapping genes with a viral match."""
    spans = [(core.start, core.end) for core in boundary.confirmed_seeds]
    spans += [
        (gene.start, gene.end)
        for gene in genes
        if gene.start < boundary.end
        and gene.end > boundary.start
        and gene.id in taxonomy
        and has_viral_match(taxonomy[gene.id])
    ]
    return (min(start for start, _ in spans), max(end for _, end in spans)) if spans else None


def _trim_boundary(
    boundary: RefinedBoundary, sequence: str, genes: Sequence[GeneSpan], taxonomy: Mapping[str, GeneTaxonomy]
) -> None:
    """Trim one boundary by GC, snap the moved ends to gene edges and mark the moved sides."""
    protected = _protected_span(boundary, genes, taxonomy)
    pre_start, pre_end = boundary.start, boundary.end
    start, end = trimmed_span(sequence, pre_start, pre_end, core=protected)
    if start != pre_start:
        protected_start = protected[0] if protected else end
        start = snap_to_gene_edge(start, genes, side="left", low=pre_start, high=min(protected_start, end - 1))
    if end != pre_end:
        protected_end = protected[1] if protected else start
        end = snap_to_gene_edge(end, genes, side="right", low=max(protected_end, start + 1), high=pre_end)
    boundary.start, boundary.end = start, end
    if start != pre_start:
        boundary.left_stop_reason = GC_TRIM_REASON
    if end != pre_end:
        boundary.right_stop_reason = GC_TRIM_REASON


def snap_to_gene_edge(
    position: int, genes: Iterable[GeneSpan], *, side: Literal["left", "right"], low: int, high: int
) -> int:
    """Move an end inside a gene to the nearest edge of a gene that contains it.

    Candidates are the start and end of every gene with ``start < position < end`` that lie
    in ``[low, high]``. On a tie the outer edge wins (lower for the left end, higher for the
    right end), so the gene stays inside the region. Without a candidate, or when the end
    lies on a gene edge, ``position`` is returned.
    """
    candidates = [
        edge
        for gene in genes
        if gene.start < position < gene.end
        for edge in (gene.start, gene.end)
        if low <= edge <= high
    ]
    if not candidates:
        return position
    outer = 1 if side == "left" else -1
    return min(candidates, key=lambda edge: (abs(edge - position), outer * edge))


def trimmed_span(sequence: str, start: int, end: int, *, core: tuple[int, int] | None) -> tuple[int, int]:
    """Return start-end after trimming host-like 2 kb windows from both edges, outside ``core``."""
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
