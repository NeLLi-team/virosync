"""Blocked kinase markers never make a repeat-bounded or residual candidate."""

from pathlib import Path

from test_terminal_repeats import _adjacent_repeat_sequence, _Marker, _query

from virosync.pipeline.phase1.seed_merger import ConfirmedSeed
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.phase2.terminal_repeats import refine_boundaries_with_terminal_repeats


def _partition(tmp_path: Path, markers: list[_Marker], cores: tuple[ConfirmedSeed, ...]) -> list[RefinedBoundary]:
    sequence, _ = _adjacent_repeat_sequence()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n", encoding="utf-8")
    parent = RefinedBoundary(
        scaffold="contig",
        start=0,
        end=len(sequence),
        seed_id="seed-1",
        original_start=min(core.start for core in cores),
        original_end=max(core.end for core in cores),
        confirmed_seeds=cores,
        seed_sources=["hhg", "marker_validation"],
    )
    return refine_boundaries_with_terminal_repeats(
        [parent],
        raw_genome_path=genome_path,
        validated_markers=markers,
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
    )


def _spans(boundaries: list[RefinedBoundary]) -> list[tuple[int, int]]:
    return [(boundary.start, boundary.end) for boundary in boundaries]


def test_repeat_pair_holding_only_a_blocked_marker_is_not_a_candidate(tmp_path: Path) -> None:
    markers = [
        _Marker("contig", 220, 250, "first|aa1-30", hmm_target="GVOGm0003"),
        _Marker("contig", 640, 660, "second|aa1-20", hmm_target="VS000087"),
        _Marker("contig", 1000, 1020, "neighbor|aa1-20", hmm_target="VS000084"),
    ]

    refined = _partition(tmp_path, markers, (ConfirmedSeed("first", "contig", 200, 260),))

    assert _spans(refined) == [(100, 420)]


def test_blocked_only_residual_is_dropped_and_real_repeat_pairs_stay(tmp_path: Path) -> None:
    markers = [
        _Marker("contig", 220, 250, "first|aa1-30", hmm_target="GVOGm0003"),
        _Marker("contig", 640, 660, "second|aa1-20", hmm_target="plv_mcp"),
        _Marker("contig", 1000, 1020, "neighbor|aa1-20", hmm_target="VS000087_Pkinase"),
    ]

    refined = _partition(
        tmp_path, markers, (ConfirmedSeed("first", "contig", 200, 260), ConfirmedSeed("second", "contig", 620, 680))
    )

    assert _spans(refined) == [(100, 420), (520, 840)]


def test_blocked_marker_inside_a_real_repeat_pair_stays_in_that_child(tmp_path: Path) -> None:
    markers = [
        _Marker("contig", 220, 235, "first|aa1-15", hmm_target="GVOGm0003"),
        _Marker("contig", 235, 250, "kinase|aa1-15", hmm_target="VS000087"),
        _Marker("contig", 640, 660, "second|aa1-20", hmm_target="plv_mcp"),
    ]

    refined = _partition(
        tmp_path, markers, (ConfirmedSeed("first", "contig", 200, 260), ConfirmedSeed("second", "contig", 620, 680))
    )

    assert _spans(refined) == [(100, 420), (520, 840)]
    assert refined[0].seed_sources == ["hhg", "marker_validation"]


def test_blocked_marker_in_a_residual_with_a_real_marker_stays_in_that_residual(tmp_path: Path) -> None:
    markers = [
        _Marker("contig", 220, 250, "first|aa1-30", hmm_target="GVOGm0003"),
        _Marker("contig", 640, 660, "second|aa1-20", hmm_target="plv_mcp"),
        _Marker("contig", 1000, 1020, "neighbor|aa1-20", hmm_target="mirus_core"),
        _Marker("contig", 1040, 1060, "kinase|aa1-20", hmm_target="VS000087"),
    ]

    refined = _partition(
        tmp_path,
        markers,
        (
            ConfirmedSeed("first", "contig", 200, 260),
            ConfirmedSeed("second", "contig", 620, 680),
            ConfirmedSeed("neighbor", "contig", 980, 1080),
        ),
    )

    assert [span[0] for span in _spans(refined)] == [100, 520, 840]


def test_parent_stays_whole_when_every_repeat_pair_holds_only_blocked_markers(tmp_path: Path) -> None:
    markers = [
        _Marker("contig", 220, 250, "first|aa1-30", hmm_target="VS000087"),
        _Marker("contig", 640, 660, "second|aa1-20", hmm_target="VS000084"),
    ]

    refined = _partition(tmp_path, markers, (ConfirmedSeed("independent-source", "contig", 430, 500),))

    assert [(boundary.start, boundary.tir_status) for boundary in refined] == [(0, "not_detected")]
