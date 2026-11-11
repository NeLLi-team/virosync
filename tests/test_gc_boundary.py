"""Tests for the conditional GC edge trim of V9g (virosync.pipeline.phase2.gc_boundary)."""

import random
from dataclasses import dataclass
from pathlib import Path

from virosync.pipeline.phase1.seed_merger import ConfirmedSeed
from virosync.pipeline.phase2 import gc_boundary
from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary


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


@dataclass(frozen=True)
class Gene:
    start: int
    end: int
    id: str = ""


def test_left_end_inside_a_gene_moves_to_the_nearer_gene_edge() -> None:
    genes = [Gene(5_000, 6_000)]

    assert gc_boundary.snap_to_gene_edge(5_300, genes, side="left", low=4_000, high=9_000) == 5_000


def test_right_end_inside_a_gene_moves_to_the_nearer_gene_edge() -> None:
    genes = [Gene(8_000, 9_000)]

    assert gc_boundary.snap_to_gene_edge(8_700, genes, side="right", low=6_000, high=10_000) == 9_000


def test_left_tie_keeps_the_gene_by_taking_the_outer_edge() -> None:
    genes = [Gene(5_000, 6_000)]

    assert gc_boundary.snap_to_gene_edge(5_500, genes, side="left", low=4_000, high=9_000) == 5_000


def test_right_tie_keeps_the_gene_by_taking_the_outer_edge() -> None:
    genes = [Gene(5_000, 6_000)]

    assert gc_boundary.snap_to_gene_edge(5_500, genes, side="right", low=4_000, high=9_000) == 6_000


def test_left_snap_never_moves_outward_past_the_pre_trim_start() -> None:
    genes = [Gene(5_000, 6_000)]

    assert gc_boundary.snap_to_gene_edge(5_300, genes, side="left", low=5_100, high=9_000) == 6_000


def test_right_snap_never_moves_outward_past_the_pre_trim_end() -> None:
    genes = [Gene(8_000, 9_000)]

    assert gc_boundary.snap_to_gene_edge(8_700, genes, side="right", low=6_000, high=8_900) == 8_000


def test_left_snap_never_enters_the_seed_core() -> None:
    genes = [Gene(5_000, 6_000)]

    assert gc_boundary.snap_to_gene_edge(5_700, genes, side="left", low=4_000, high=5_800) == 5_000


def test_right_snap_never_enters_the_seed_core() -> None:
    genes = [Gene(8_000, 9_000)]

    assert gc_boundary.snap_to_gene_edge(8_300, genes, side="right", low=8_200, high=10_000) == 9_000


def test_end_stays_without_a_valid_gene_edge() -> None:
    genes = [Gene(5_000, 6_000)]

    assert gc_boundary.snap_to_gene_edge(5_300, genes, side="left", low=5_100, high=5_800) == 5_300


def test_end_on_a_gene_edge_is_not_snapped() -> None:
    genes = [Gene(5_000, 6_000)]

    assert gc_boundary.snap_to_gene_edge(5_000, genes, side="left", low=4_000, high=9_000) == 5_000


def test_end_inside_nested_genes_takes_the_nearest_edge_of_any_of_them() -> None:
    genes = [Gene(5_000, 7_000), Gene(5_200, 5_600)]

    assert gc_boundary.snap_to_gene_edge(5_300, genes, side="left", low=4_000, high=9_000) == 5_200


def write_genome(tmp_path: Path) -> Path:
    path = tmp_path / "genome.fa"
    path.write_text(f">ctg\n{genome(EVE)}\n", encoding="ascii")
    return path


def boundary(start: int, end: int) -> RefinedBoundary:
    return RefinedBoundary(
        "ctg",
        start,
        end,
        confirmed_seeds=(ConfirmedSeed("core", "ctg", 66_000, 74_000),),
        left_stop_reason="max_genes",
        right_stop_reason="seed_end",
    )


def test_trim_snaps_a_moved_start_to_the_gene_edge_and_marks_its_reason(tmp_path: Path) -> None:
    genes = {"ctg": [Gene(59_000, 61_500)]}

    (trimmed,) = gc_boundary.trim_boundaries_by_gc([boundary(52_000, 80_000)], write_genome(tmp_path), genes, {})

    assert (trimmed.start, trimmed.left_stop_reason) == (59_000, "gc_trim")


def test_trim_leaves_the_reason_of_a_side_it_did_not_move(tmp_path: Path) -> None:
    genes = {"ctg": []}

    (trimmed,) = gc_boundary.trim_boundaries_by_gc([boundary(52_000, 80_000)], write_genome(tmp_path), genes, {})

    assert (trimmed.end, trimmed.right_stop_reason) == (80_000, "seed_end")


def test_trim_keeps_an_unmoved_end_inside_its_gene(tmp_path: Path) -> None:
    genes = {"ctg": [Gene(79_000, 81_000)]}

    (trimmed,) = gc_boundary.trim_boundaries_by_gc([boundary(52_000, 80_000)], write_genome(tmp_path), genes, {})

    assert trimmed.end == 80_000


def test_trim_keeps_the_reason_when_the_snap_returns_to_the_pre_trim_start(tmp_path: Path) -> None:
    genes = {"ctg": [Gene(58_000, 62_500)]}

    (trimmed,) = gc_boundary.trim_boundaries_by_gc([boundary(58_000, 80_000)], write_genome(tmp_path), genes, {})

    assert (trimmed.start, trimmed.left_stop_reason) == (58_000, "max_genes")


def test_trim_marks_the_reason_of_a_moved_start_without_a_gene(tmp_path: Path) -> None:
    (trimmed,) = gc_boundary.trim_boundaries_by_gc([boundary(52_000, 80_000)], write_genome(tmp_path), {}, {})

    assert (trimmed.start, trimmed.left_stop_reason) == (60_000, "gc_trim")


def taxonomy_record(name: str, start: int, end: int, prefix: str) -> GeneTaxonomy:
    return GeneTaxonomy(name, "ctg", start, end, top10_prefixes=[prefix], top10_pidents=[50.0])


def test_trim_does_not_move_a_start_past_a_gene_with_a_viral_match(tmp_path: Path) -> None:
    genes = {"ctg": [Gene(59_000, 61_500, "viral")]}
    taxonomy = {"viral": taxonomy_record("viral", 59_000, 61_500, "NCLDV__")}

    (trimmed,) = gc_boundary.trim_boundaries_by_gc([boundary(52_000, 80_000)], write_genome(tmp_path), genes, taxonomy)

    assert (trimmed.start, trimmed.left_stop_reason) == (58_000, "gc_trim")


def test_trim_does_not_move_an_end_past_a_gene_with_a_viral_match(tmp_path: Path) -> None:
    genes = {"ctg": [Gene(81_000, 84_000, "viral")]}
    taxonomy = {"viral": taxonomy_record("viral", 81_000, 84_000, "NCLDV__")}

    (trimmed,) = gc_boundary.trim_boundaries_by_gc([boundary(60_000, 88_000)], write_genome(tmp_path), genes, taxonomy)

    assert (trimmed.end, trimmed.right_stop_reason) == (84_000, "gc_trim")


def test_trim_still_moves_a_start_past_a_gene_without_a_viral_match(tmp_path: Path) -> None:
    genes = {"ctg": [Gene(59_000, 61_500, "host")]}
    taxonomy = {"host": taxonomy_record("host", 59_000, 61_500, "EUK__")}

    (trimmed,) = gc_boundary.trim_boundaries_by_gc([boundary(52_000, 80_000)], write_genome(tmp_path), genes, taxonomy)

    assert trimmed.start == 59_000


def test_trim_ignores_a_gene_with_a_viral_match_outside_the_region(tmp_path: Path) -> None:
    genes = {"ctg": [Gene(40_000, 45_000, "viral")]}
    taxonomy = {"viral": taxonomy_record("viral", 40_000, 45_000, "NCLDV__")}

    (trimmed,) = gc_boundary.trim_boundaries_by_gc([boundary(52_000, 80_000)], write_genome(tmp_path), genes, taxonomy)

    assert trimmed.start == 60_000


def test_trim_protects_a_validated_viral_marker_without_a_top10_hit(tmp_path: Path) -> None:
    genes = {"ctg": [Gene(59_000, 61_500, "marker")]}
    marker = GeneTaxonomy("marker", "ctg", 59_000, 61_500, phase1_validated_viral_marker=True)

    (trimmed,) = gc_boundary.trim_boundaries_by_gc(
        [boundary(52_000, 80_000)], write_genome(tmp_path), genes, {"marker": marker}
    )

    assert trimmed.start == 58_000
