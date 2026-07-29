"""Tests for the Phase-2 validated-marker floor METADATA computation.

``annotate_boundaries_with_marker_floor`` records the min..max span of a
boundary's own validated viral markers as ``marker_floor_start`` /
``marker_floor_end`` WITHOUT mutating the boundary's ``start`` / ``end``. The
Phase-3 re-admit pass consumes this metadata on REJECTED boundaries only, so an
accepted boundary is never modified and the accepted set cannot regress (the
earlier mutating ``clamp`` approach extended already-accepted regions and dropped
them MEDIUM->LOW -> NCLDV loss).
"""

from types import SimpleNamespace

from virosync.pipeline.phase2.boundary_refiner import (
    RefinedBoundary,
    annotate_boundaries_with_marker_floor,
)


def _marker(scaffold, start, end, status="validated", porf=None):
    return SimpleNamespace(
        scaffold=scaffold,
        start=start,
        end=end,
        validation_status=status,
        # Real ValidatedMarkerHit carries the pORF id; the floor counts distinct
        # proteins, so tests must distinguish two genes from two hits on one gene.
        query_porf=porf if porf is not None else f"{scaffold}_{start}_{end}",
    )


def _boundary(scaffold, start, end, orig_start, orig_end, seed_id="seed"):
    return RefinedBoundary(
        scaffold=scaffold,
        start=start,
        end=end,
        seed_id=seed_id,
        original_start=orig_start,
        original_end=orig_end,
    )


def test_floor_metadata_recorded_without_mutating_boundary():
    # Region 1 analogue: 13 kb seed collapsed to 660 bp, 3 validated markers span wider.
    b = _boundary("S1", 59404, 60064, orig_start=54513, orig_end=67844)
    markers = [
        _marker("S1", 59513, 60013),
        _marker("S1", 62449, 65844),
        _marker("S1", 67054, 67842),
    ]
    n = annotate_boundaries_with_marker_floor([b], markers)
    assert n == 1
    # Boundary coordinates are NOT mutated -- only metadata is recorded.
    assert (b.start, b.end) == (59404, 60064)
    # Floor metadata records the widened validated-marker span.
    assert b.marker_floor_start == 59404  # min(59404, 59513) -> existing start kept
    assert b.marker_floor_end == 67842  # max(60064, 67842) -> extended to marker span
    assert b.marker_floor_end - b.marker_floor_start > 5000


def test_floor_requires_two_validated_markers():
    # A single validated marker must NOT record a floor (specificity).
    b = _boundary("S1", 3682, 10647, orig_start=0, orig_end=21500)
    n = annotate_boundaries_with_marker_floor([b], [_marker("S1", 600, 1200)])
    assert n == 0
    assert b.marker_floor_start is None and b.marker_floor_end is None
    assert (b.start, b.end) == (3682, 10647)


def test_floor_ignores_unvalidated_markers():
    # Host-like / unvalidated hits are never protected even if there are many.
    b = _boundary("S1", 3000, 4000, orig_start=0, orig_end=20000)
    markers = [
        _marker("S1", 500, 900, status="unvalidated"),
        _marker("S1", 1000, 1400, status="supported"),
        _marker("S1", 1500, 1900, status="unvalidated"),
    ]
    n = annotate_boundaries_with_marker_floor([b], markers)
    assert n == 0
    assert b.marker_floor_start is None and b.marker_floor_end is None
    assert (b.start, b.end) == (3000, 4000)


def test_floor_is_scoped_to_boundary_own_seed_span():
    # Markers from a DIFFERENT region on the same scaffold (outside this
    # boundary's seed span) must not widen the recorded floor -- otherwise a small
    # region would claim the whole scaffold's viral content.
    b = _boundary("S1", 57601, 60140, orig_start=57000, orig_end=73285)
    markers = [
        _marker("S1", 65859, 68285, porf="S1_60"),  # in-span
        _marker("S1", 69100, 70200, porf="S1_61"),  # in-span, second gene
        _marker("S1", 116748, 118247, porf="S1_120"),  # far away, different region
        _marker("S1", 188295, 189641, porf="S1_190"),  # far away, different region
    ]
    n = annotate_boundaries_with_marker_floor([b], markers)
    assert n == 1
    assert (b.start, b.end) == (57601, 60140)  # not mutated
    # Floor reaches the furthest IN-SPAN marker (70200) and stops there; the
    # markers at 118247 and 189641 belong to other regions and are excluded.
    assert b.marker_floor_start == 57601
    assert b.marker_floor_end == 70200


def test_validated_novel_counts_as_validated():
    b = _boundary("S1", 1000, 1500, orig_start=0, orig_end=20000)
    markers = [
        _marker("S1", 2000, 2500, status="validated_novel"),
        _marker("S1", 8000, 8500, status="validated_novel"),
    ]
    n = annotate_boundaries_with_marker_floor([b], markers)
    assert n == 1
    assert b.marker_floor_start == 1000
    assert b.marker_floor_end == 8500


def test_no_floor_when_boundary_already_contains_span():
    # A boundary already covering its markers yields no re-admit alternative.
    b = _boundary("S1", 1000, 20000, orig_start=1000, orig_end=20000)
    markers = [_marker("S1", 5000, 6000), _marker("S1", 8000, 9000)]
    n = annotate_boundaries_with_marker_floor([b], markers)
    assert n == 0
    assert b.marker_floor_start is None and b.marker_floor_end is None
    assert (b.start, b.end) == (1000, 20000)


def test_empty_inputs_are_safe():
    assert annotate_boundaries_with_marker_floor([], [_marker("S1", 1, 2)]) == 0
    b = _boundary("S1", 1, 2, 0, 100)
    assert annotate_boundaries_with_marker_floor([b], []) == 0
    assert (b.start, b.end) == (1, 2)
    assert b.marker_floor_start is None and b.marker_floor_end is None


def test_two_hits_on_one_gene_do_not_make_a_floor():
    # The floor exists to recover marker-DENSE regions. A single marker-bearing
    # protein that produced two HMM hits is not marker-dense, and letting it
    # widen the span would feed the Phase-3 re-admit a region carrying only one
    # genuine marker gene.
    b = _boundary("S1", 57601, 60140, orig_start=57000, orig_end=73285)
    markers = [
        _marker("S1", 65859, 68285, porf="S1_60"),
        _marker("S1", 65859, 68285, porf="S1_60"),
    ]
    assert annotate_boundaries_with_marker_floor([b], markers) == 0
    assert b.marker_floor_start is None
    assert b.marker_floor_end is None
