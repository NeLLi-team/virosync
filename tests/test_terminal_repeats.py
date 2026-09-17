from __future__ import annotations

import logging
import random
from dataclasses import dataclass, replace
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio import SeqIO

from virosync.ablation import InterventionCounts
from virosync.config import PipelineConfig
from virosync.orchestration._flows.single_genome import phase2
from virosync.pipeline.phase1.seed_merger import MergedSeed
from virosync.pipeline.phase2 import boundary_refiner, taxonomy_seed_refiner, terminal_repeats
from virosync.pipeline.phase2.boundary_diamond import (
    GeneTaxonomy,
    GenomeDiamondQuery,
    SeedGeneMapping,
    pORF,
)
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.phase2.terminal_repeats import (
    MAX_TIR_ARM_BP,
    find_target_site_duplication,
    find_terminal_inverted_repeats,
    refine_boundaries_with_terminal_repeats,
)

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True, slots=True)
class _Marker:
    scaffold: str
    start: int
    end: int
    query_porf: str
    validation_status: str = "validated"
    hmm_target: str = "ncldv_mcp"


def _random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _at_rich_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choices("ACGT", weights=(35, 15, 15, 35), k=length))


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def _planted_repeat(
    *,
    arm_length: int = 60,
    mutated_positions: tuple[int, ...] = (),
    seed: int = 17,
) -> tuple[str, tuple[int, int, int, int]]:
    rng = random.Random(seed)
    left_flank = _random_dna(99, rng) + "A"
    arm = _random_dna(arm_length, rng)
    paired_arm = list(arm)
    substitutions = {"A": "C", "C": "G", "G": "T", "T": "A"}
    for position in mutated_positions:
        paired_arm[position] = substitutions[paired_arm[position]]
    interior = "A" + _random_dna(198, rng) + "A"
    right_flank = "A" + _random_dna(99, rng)
    sequence = left_flank + arm + interior + _reverse_complement("".join(paired_arm)) + right_flank
    left_start = len(left_flank)
    left_end = left_start + arm_length
    right_start = left_end + len(interior)
    right_end = right_start + arm_length
    return sequence, (left_start, left_end, right_start, right_end)


def _adjacent_repeat_sequence(
    *,
    identical_arms: bool = False,
) -> tuple[str, tuple[tuple[int, int], tuple[int, int]]]:
    rng = random.Random(221)
    prefix = _random_dna(90, rng) + "A" * 10
    first_arm = _random_dna(60, rng)
    second_arm = first_arm if identical_arms else _random_dna(60, rng)
    first_interior = "A" * 10 + _random_dna(180, rng) + "A" * 10
    between = "A" * 10 + _random_dna(80, rng) + "A" * 10
    second_interior = "A" * 10 + _random_dna(180, rng) + "A" * 10
    tail = "A" * 10 + _random_dna(350, rng)
    sequence = (
        prefix
        + first_arm
        + first_interior
        + _reverse_complement(first_arm)
        + between
        + second_arm
        + second_interior
        + _reverse_complement(second_arm)
        + tail
    )
    return sequence, ((100, 420), (520, 840))


@pytest.mark.parametrize(
    ("mutated_positions", "expected_identity"),
    [
        ((), 1.0),
        ((10, 20, 30, 40, 50, 55), 0.9),
    ],
)
def test_scanner_recovers_exact_planted_coordinates(
    mutated_positions: tuple[int, ...],
    expected_identity: float,
) -> None:
    sequence, expected = _planted_repeat(mutated_positions=mutated_positions)

    candidates = find_terminal_inverted_repeats(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        marker_start=250,
        marker_end=270,
    )

    assert len(candidates) == 1
    candidate = candidates[0]
    assert (
        candidate.left_start,
        candidate.left_end,
        candidate.right_start,
        candidate.right_end,
    ) == expected
    assert candidate.identity == pytest.approx(expected_identity)


@pytest.mark.parametrize(("arm_length", "detected"), [(49, False), (50, True)])
def test_scanner_enforces_minimum_arm_length(arm_length: int, detected: bool) -> None:
    rng = random.Random(730 + arm_length)
    arm = _random_dna(arm_length, rng)
    prefix = _random_dna(90, rng) + "A" * 10
    interior = "A" * 10 + _random_dna(180, rng) + "G" * 10
    suffix = "G" * 10 + _random_dna(90, rng)
    sequence = prefix + arm + interior + _reverse_complement(arm) + suffix
    expected = (len(prefix), len(prefix) + arm_length + len(interior) + arm_length)
    marker_start = len(prefix) + arm_length + 90

    candidates = find_terminal_inverted_repeats(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        marker_start=marker_start,
        marker_end=marker_start + 20,
    )

    assert bool(candidates) is detected
    if detected:
        assert (candidates[0].left_start, candidates[0].right_end) == expected


@pytest.mark.parametrize(
    "mutated_positions",
    [(), (10,), (10, 20, 30)],
)
def test_long_repeat_cap_preserves_outer_boundary_coordinates(
    mutated_positions: tuple[int, ...],
) -> None:
    rng = random.Random(8)
    left_flank = _random_dna(99, rng) + "A"
    right_flank = "A" + _random_dna(99, rng)
    arm = _random_dna(MAX_TIR_ARM_BP + 115, rng)
    paired_arm = list(arm)
    substitutions = {"A": "C", "C": "G", "G": "T", "T": "A"}
    for position in mutated_positions:
        paired_arm[position] = substitutions[paired_arm[position]]
    interior = _random_dna(200, rng)
    sequence = left_flank + arm + interior + _reverse_complement("".join(paired_arm)) + right_flank
    marker_start = len(left_flank) + len(arm) + 80

    candidates = find_terminal_inverted_repeats(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        marker_start=marker_start,
        marker_end=marker_start + 20,
    )

    assert len(candidates) == 1
    assert candidates[0].left_start == len(left_flank)
    assert candidates[0].right_end == len(sequence) - len(right_flank)
    assert candidates[0].arm_length == MAX_TIR_ARM_BP
    assert candidates[0].alignment_length == MAX_TIR_ARM_BP + 115
    assert candidates[0].alignment_capped is True
    assert candidates[0].identity == pytest.approx((MAX_TIR_ARM_BP - len(mutated_positions)) / MAX_TIR_ARM_BP)


def test_capped_repeat_rejects_subthreshold_reported_arm() -> None:
    rng = random.Random(8)
    left_flank = _random_dna(99, rng) + "A"
    right_flank = "A" + _random_dna(99, rng)
    arm = _random_dna(MAX_TIR_ARM_BP + 115, rng)
    paired_arm = list(arm)
    substitutions = {"A": "C", "C": "G", "G": "T", "T": "A"}
    for position in range(5, MAX_TIR_ARM_BP, 9):
        paired_arm[position] = substitutions[paired_arm[position]]
    interior = _random_dna(200, rng)
    sequence = left_flank + arm + interior + _reverse_complement("".join(paired_arm)) + right_flank
    marker_start = len(left_flank) + len(arm) + 80

    assert (
        find_terminal_inverted_repeats(
            sequence,
            scan_start=0,
            scan_end=len(sequence),
            marker_start=marker_start,
            marker_end=marker_start + 20,
        )
        == []
    )


def test_target_site_duplication_is_descriptive_only() -> None:
    sequence, expected = _planted_repeat()
    bases = list(sequence)
    bases[expected[0] - 6 : expected[0]] = "GAGGCT"
    bases[expected[3] : expected[3] + 6] = "GAGGCT"
    sequence = "".join(bases)
    candidate = find_terminal_inverted_repeats(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        marker_start=250,
        marker_end=270,
    )[0]

    assert find_target_site_duplication(sequence, candidate) == "GAGGCT"


def test_mavirus_terminal_repeats_recover_published_outer_boundary() -> None:
    """Regress the KU052222 integration site from DOI 10.1038/nature20593."""
    record = SeqIO.read(Path(__file__).parent / "fixtures" / "KU052222.fna", "fasta")
    sequence = str(record.seq)

    candidates = find_terminal_inverted_repeats(
        sequence,
        scan_start=5000,
        scan_end=35190,
        marker_start=15000,
        marker_end=25000,
    )

    assert len(candidates) == 1
    candidate = candidates[0]
    assert (candidate.left_start, candidate.right_end) == (10000, 30190)
    assert candidate.arm_length == MAX_TIR_ARM_BP
    assert candidate.alignment_length == 600
    assert candidate.alignment_capped is True
    assert find_target_site_duplication(sequence, candidate) == "GAGGCT"


def test_scanner_rejects_random_and_short_period_sequences() -> None:
    rng = random.Random(91)
    for _ in range(100):
        sequence = _random_dna(800, rng)
        assert (
            find_terminal_inverted_repeats(
                sequence,
                scan_start=0,
                scan_end=len(sequence),
                marker_start=390,
                marker_end=410,
            )
            == []
        )

    low_complexity_arm = "ACGT" * 15
    sequence = (
        _random_dna(100, rng)
        + low_complexity_arm
        + _random_dna(200, rng)
        + _reverse_complement(low_complexity_arm)
        + _random_dna(100, rng)
    )
    assert (
        find_terminal_inverted_repeats(
            sequence,
            scan_start=0,
            scan_end=len(sequence),
            marker_start=250,
            marker_end=270,
        )
        == []
    )


def test_scanner_rejects_at_rich_composition_matched_arms() -> None:
    rng = random.Random(108)
    left_arm = _at_rich_dna(60, rng)
    composition_matched_arm = list(left_arm)
    rng.shuffle(composition_matched_arm)
    sequence = (
        _at_rich_dna(100, rng)
        + left_arm
        + _at_rich_dna(200, rng)
        + _reverse_complement("".join(composition_matched_arm))
        + _at_rich_dna(100, rng)
    )

    assert (
        find_terminal_inverted_repeats(
            sequence,
            scan_start=0,
            scan_end=len(sequence),
            marker_start=250,
            marker_end=270,
        )
        == []
    )


def test_marker_anchors_deduplicate_same_outer_pair_by_longest_evidence() -> None:
    sequence, expected = _planted_repeat(arm_length=MAX_TIR_ARM_BP + 115, seed=8)
    groups = [
        terminal_repeats._MarkerGroup(
            protein_id="arm-marker",
            start=250,
            end=270,
            marker_names=frozenset({"marker-a"}),
        ),
        terminal_repeats._MarkerGroup(
            protein_id="core-marker",
            start=795,
            end=815,
            marker_names=frozenset({"marker-b"}),
        ),
    ]

    candidates = terminal_repeats._discover_marker_anchored_candidates(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        groups=groups,
    )

    assert len(candidates) == 1
    candidate = candidates[0]
    assert (candidate.left_start, candidate.right_end) == (expected[0], expected[3])
    assert candidate.arm_length == MAX_TIR_ARM_BP
    assert candidate.alignment_length == MAX_TIR_ARM_BP + 115


def test_shared_scan_index_preserves_candidates() -> None:
    sequence, _ = _planted_repeat(arm_length=MAX_TIR_ARM_BP + 115, seed=8)
    scan_start = 40
    scan_end = len(sequence) - 40
    arguments = {
        "scan_start": scan_start,
        "scan_end": scan_end,
        "marker_start": 795,
        "marker_end": 815,
    }

    uncached = find_terminal_inverted_repeats(sequence, **arguments)
    cached = find_terminal_inverted_repeats(
        sequence,
        **arguments,
        _scan_seed_index=terminal_repeats._build_scan_seed_index(sequence[scan_start:scan_end]),
    )

    assert cached == uncached


def test_duplicate_repeat_pairs_are_ambiguous() -> None:
    rng = random.Random(32)
    arm = _random_dna(60, rng)
    sequence = (
        _random_dna(99, rng)
        + "A"
        + arm
        + "A"
        + _random_dna(198, rng)
        + "A"
        + _reverse_complement(arm)
        + "A"
        + _random_dna(28, rng)
        + "A"
        + _reverse_complement(arm)
        + "A"
        + _random_dna(99, rng)
    )

    candidates = find_terminal_inverted_repeats(
        sequence,
        scan_start=0,
        scan_end=len(sequence),
        marker_start=250,
        marker_end=270,
    )

    assert len(candidates) == 2
    assert {(candidate.left_start, candidate.right_end) for candidate in candidates} == {
        (100, 420),
        (100, 510),
    }


def _query(sequence_length: int) -> GenomeDiamondQuery:
    mapping = SeedGeneMapping(
        seed_id="seed-1",
        scaffold="contig",
        seed_start=150,
        seed_end=350,
        flank_start_bp=0,
        flank_end_bp=sequence_length,
    )
    return GenomeDiamondQuery(seed_gene_mappings={mapping.seed_id: mapping})


@pytest.mark.parametrize(
    ("boundary_start", "boundary_end"),
    [
        (120, 400),
        (50, 470),
    ],
)
def test_unique_repeat_overrides_inward_and_outward_boundaries(
    tmp_path: Path,
    boundary_start: int,
    boundary_end: int,
) -> None:
    sequence, expected = _planted_repeat()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    boundary = RefinedBoundary(
        scaffold="contig",
        start=boundary_start,
        end=boundary_end,
        seed_id="seed-1",
        original_start=150,
        original_end=350,
        marker_floor_start=140,
        marker_floor_end=360,
    )

    refined = refine_boundaries_with_terminal_repeats(
        [boundary],
        raw_genome_path=genome_path,
        validated_markers=[_Marker("contig", 250, 270, "marker-1")],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
    )[0]

    assert (refined.start, refined.end) == (expected[0], expected[3])
    assert (refined.pre_tir_start, refined.pre_tir_end) == (boundary_start, boundary_end)
    assert refined.tir_status == "detected"
    assert refined.tir_present is True
    assert refined.tir_boundary_override is True
    assert refined.tir_alignment_capped is False
    assert refined.tir_candidate_count == 1
    assert refined.marker_floor_start is None
    assert refined.marker_floor_end is None


def test_coverage_veto_retains_candidate_without_changing_boundary(tmp_path: Path) -> None:
    sequence, expected = _planted_repeat()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    boundary = RefinedBoundary(
        scaffold="contig",
        start=120,
        end=400,
        seed_id="seed-1",
        original_start=150,
        original_end=350,
    )

    refined = refine_boundaries_with_terminal_repeats(
        [boundary],
        raw_genome_path=genome_path,
        validated_markers=[_Marker("contig", 250, 270, "marker-1")],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": [pORF(id="unclassified", scaffold="contig", start=110, end=130, strand="+")]},
        extension_bp=500,
    )[0]

    assert (refined.start, refined.end) == (120, 400)
    assert refined.tir_status == "taxonomy_incomplete"
    assert refined.tir_present is True
    assert refined.tir_boundary_override is False
    assert (refined.tir_left_start, refined.tir_right_end) == (expected[0], expected[3])


def test_no_marker_anchor_is_explicit_and_coordinate_neutral(tmp_path: Path) -> None:
    sequence, _ = _planted_repeat()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    boundary = RefinedBoundary(
        scaffold="contig",
        start=120,
        end=400,
        seed_id="seed-1",
        original_start=150,
        original_end=350,
    )

    refined = refine_boundaries_with_terminal_repeats(
        [boundary],
        raw_genome_path=genome_path,
        validated_markers=[],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
    )[0]

    assert (refined.start, refined.end) == (120, 400)
    assert refined.tir_status == "no_marker_anchor"
    assert refined.tir_present is False
    assert refined.tir_candidate_count == 0


def test_disjoint_tir_pairs_partition_parent_and_keep_marker_neighbor(
    tmp_path: Path,
) -> None:
    sequence, expected_pairs = _adjacent_repeat_sequence()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    parent = RefinedBoundary(
        scaffold="contig",
        start=0,
        end=len(sequence),
        seed_id="seed-1",
        original_start=0,
        original_end=len(sequence),
        seed_sources=["hhg", "marker_validation"],
        seed_confidence="high",
        seed_hhg_score=0.99,
        predicted_family="MIXED",
        region_classification_ncldv_markers=99,
        region_classification_vp_plv_markers=99,
        region_classification_mirus_markers=99,
        confidence=0.91,
        posterior_probability=0.92,
        cub_deviation=0.4,
    )
    markers = [
        _Marker("contig", 220, 235, "first|aa1-15", hmm_target="gvogm0003"),
        _Marker("contig", 235, 250, "first|aa16-30", hmm_target="gvogm0004"),
        _Marker("contig", 640, 660, "second|aa1-20", hmm_target="plv_mcp"),
        _Marker("contig", 1000, 1020, "neighbor|aa1-20", hmm_target="mirus_core"),
    ]

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        raw_genome_path=genome_path,
        validated_markers=markers,
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
    )

    assert [(boundary.start, boundary.end) for boundary in refined] == [
        *expected_pairs,
        (840, len(sequence)),
    ]
    assert [boundary.tir_status for boundary in refined] == [
        "detected",
        "detected",
        "not_detected",
    ]
    assert [boundary.predicted_family for boundary in refined] == [
        "NCLDV",
        "PPV",
        "MIRUS",
    ]
    assert [boundary.region_classification_ncldv_markers for boundary in refined] == [2, 0, 0]
    assert all(boundary.original_start == boundary.start for boundary in refined)
    assert all(boundary.original_end == boundary.end for boundary in refined)
    assert all((boundary.pre_tir_start, boundary.pre_tir_end) == (0, len(sequence)) for boundary in refined)
    assert all(boundary.seed_id == parent.seed_id for boundary in refined)
    assert all(boundary.confidence == 0.0 for boundary in refined)
    assert all(boundary.posterior_probability == 0.0 for boundary in refined)
    assert all(boundary.cub_deviation == 0.0 for boundary in refined)


def test_partition_keeps_marker_free_viral_taxonomy_neighbor(tmp_path: Path) -> None:
    sequence, expected_pairs = _adjacent_repeat_sequence()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    parent = RefinedBoundary(
        scaffold="contig",
        start=0,
        end=len(sequence),
        seed_id="seed-1",
        original_start=0,
        original_end=len(sequence),
    )
    markers = [
        _Marker("contig", 220, 240, "first", hmm_target="gvogm0003"),
        _Marker("contig", 640, 660, "second", hmm_target="plv_mcp"),
    ]
    viral_gene = pORF(
        id="viral-neighbor",
        scaffold="contig",
        start=1000,
        end=1080,
        strand="+",
    )
    taxonomy = GeneTaxonomy(
        porf_id=viral_gene.id,
        scaffold="contig",
        start=viral_gene.start,
        end=viral_gene.end,
        top10_prefixes=["NCLDV__"],
        top10_pidents=[70.0],
        has_viral=True,
        has_hit=True,
    )

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        raw_genome_path=genome_path,
        validated_markers=markers,
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={viral_gene.id: taxonomy},
        proteome_index={"contig": [viral_gene]},
        extension_bp=500,
    )

    assert [(boundary.start, boundary.end) for boundary in refined] == [
        *expected_pairs,
        (840, len(sequence)),
    ]
    assert refined[-1].tir_status == "not_detected"
    assert refined[-1].predicted_family == "UNKNOWN"
    assert refined[-1].seed_sources == []


def test_partition_preserves_marker_in_current_parent_beyond_original_span(
    tmp_path: Path,
) -> None:
    sequence, expected = _planted_repeat()
    sequence += _random_dna(300, random.Random(407))
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    parent = RefinedBoundary(
        scaffold="contig",
        start=0,
        end=len(sequence),
        seed_id="seed-1",
        original_start=200,
        original_end=300,
    )

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        raw_genome_path=genome_path,
        validated_markers=[
            _Marker("contig", 220, 240, "tir-marker"),
            _Marker("contig", 640, 660, "neighbor-marker"),
        ],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
    )

    assert [(boundary.start, boundary.end) for boundary in refined] == [
        (expected[0], expected[3]),
        (expected[3], len(sequence)),
    ]
    assert [boundary.tir_status for boundary in refined] == [
        "detected",
        "not_detected",
    ]


def test_partition_clips_residual_to_current_parent(tmp_path: Path) -> None:
    planted, expected = _planted_repeat()
    sequence = _random_dna(400, random.Random(503)) + planted
    expected = tuple(coordinate + 400 for coordinate in expected)
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    parent = RefinedBoundary(
        scaffold="contig",
        start=100,
        end=200,
        seed_id="seed-1",
        original_start=100,
        original_end=800,
    )
    viral_gene = pORF(
        id="viral-residual",
        scaffold="contig",
        start=120,
        end=150,
        strand="+",
    )
    taxonomy = GeneTaxonomy(
        porf_id=viral_gene.id,
        scaffold="contig",
        start=viral_gene.start,
        end=viral_gene.end,
        top10_prefixes=["NCLDV__"],
        top10_pidents=[70.0],
        has_viral=True,
        has_hit=True,
    )

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        raw_genome_path=genome_path,
        validated_markers=[_Marker("contig", 620, 640, "tir-marker")],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={viral_gene.id: taxonomy},
        proteome_index={"contig": [viral_gene]},
        extension_bp=500,
    )

    assert [(boundary.start, boundary.end) for boundary in refined] == [
        (parent.start, parent.end),
        (expected[0], expected[3]),
    ]
    assert [boundary.tir_status for boundary in refined] == [
        "not_detected",
        "detected",
    ]


def test_identical_tir_cross_pair_keeps_parent_ambiguous(tmp_path: Path) -> None:
    sequence, _ = _adjacent_repeat_sequence(identical_arms=True)
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    parent = RefinedBoundary(
        scaffold="contig",
        start=0,
        end=len(sequence),
        seed_id="seed-1",
        original_start=0,
        original_end=len(sequence),
    )

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        raw_genome_path=genome_path,
        validated_markers=[
            _Marker("contig", 220, 240, "first"),
            _Marker("contig", 640, 660, "second"),
        ],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
    )

    assert len(refined) == 1
    assert (refined[0].start, refined[0].end) == (parent.start, parent.end)
    assert refined[0].tir_status == "ambiguous"
    assert refined[0].tir_candidate_count >= 3
    assert refined[0].tir_boundary_override is False


def test_crossing_porf_does_not_override_tir_partition(tmp_path: Path) -> None:
    sequence, _ = _adjacent_repeat_sequence()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    parent = RefinedBoundary(
        scaffold="contig",
        start=0,
        end=len(sequence),
        seed_id="seed-1",
        original_start=0,
        original_end=len(sequence),
    )
    crossing_gene = pORF(
        id="crossing",
        scaffold="contig",
        start=830,
        end=850,
        strand="+",
    )
    taxonomy = GeneTaxonomy(
        porf_id=crossing_gene.id,
        scaffold="contig",
        start=crossing_gene.start,
        end=crossing_gene.end,
        has_viral=False,
        has_hit=False,
    )
    viral_neighbor = pORF(
        id="viral-neighbor",
        scaffold="contig",
        start=1000,
        end=1080,
        strand="+",
    )
    neighbor_taxonomy = GeneTaxonomy(
        porf_id=viral_neighbor.id,
        scaffold="contig",
        start=viral_neighbor.start,
        end=viral_neighbor.end,
        top10_prefixes=["NCLDV__"],
        top10_pidents=[70.0],
        has_viral=True,
        has_hit=True,
    )

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        raw_genome_path=genome_path,
        validated_markers=[
            _Marker("contig", 220, 240, "first"),
            _Marker("contig", 640, 660, "second"),
        ],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={
            crossing_gene.id: taxonomy,
            viral_neighbor.id: neighbor_taxonomy,
        },
        proteome_index={"contig": [crossing_gene, viral_neighbor]},
        extension_bp=500,
    )

    assert [(boundary.start, boundary.end) for boundary in refined] == [
        (100, 420),
        (520, 840),
        (840, len(sequence)),
    ]
    assert [boundary.tir_status for boundary in refined] == [
        "detected",
        "detected",
        "not_detected",
    ]


def test_partition_vetoes_complete_marker_group_crossing_child_boundary(
    tmp_path: Path,
) -> None:
    sequence, _ = _adjacent_repeat_sequence()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    parent = RefinedBoundary(
        scaffold="contig",
        start=0,
        end=len(sequence),
        seed_id="seed-1",
        original_start=0,
        original_end=len(sequence),
    )

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        raw_genome_path=genome_path,
        validated_markers=[
            _Marker("contig", 220, 240, "first"),
            _Marker("contig", 640, 660, "second"),
            _Marker("contig", 830, 838, "crossing|aa1-8"),
            _Marker("contig", 842, 850, "crossing|aa9-16"),
        ],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
    )

    assert len(refined) == 1
    assert (refined[0].start, refined[0].end) == (parent.start, parent.end)
    assert refined[0].tir_status == "ambiguous"
    assert refined[0].tir_candidate_count == 2


def test_shared_pair_does_not_collapse_distinct_boundaries(tmp_path: Path) -> None:
    sequence, expected = _planted_repeat()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    boundaries = [
        RefinedBoundary(
            scaffold="contig",
            start=start,
            end=end,
            seed_id=seed_id,
            original_start=start,
            original_end=end,
        )
        for seed_id, start, end in (
            ("seed-1", 180, 245),
            ("seed-2", 260, 340),
        )
    ]
    mappings = {
        boundary.seed_id: SeedGeneMapping(
            seed_id=boundary.seed_id,
            scaffold="contig",
            seed_start=boundary.start,
            seed_end=boundary.end,
            flank_start_bp=0,
            flank_end_bp=len(sequence),
        )
        for boundary in boundaries
    }

    refined = refine_boundaries_with_terminal_repeats(
        boundaries,
        raw_genome_path=genome_path,
        validated_markers=[
            _Marker("contig", 220, 230, "marker-1"),
            _Marker("contig", 280, 290, "marker-2"),
        ],
        boundary_diamond_query=GenomeDiamondQuery(seed_gene_mappings=mappings),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
    )

    assert [(boundary.start, boundary.end) for boundary in refined] == [(180, 245), (260, 340)]
    assert {boundary.tir_status for boundary in refined} == {"shared_pair"}
    assert all(boundary.tir_present for boundary in refined)
    assert not any(boundary.tir_boundary_override for boundary in refined)
    assert {(boundary.tir_left_start, boundary.tir_right_end) for boundary in refined} == {(expected[0], expected[3])}


def test_tir_proposals_equal_to_existing_boundary_are_reverted(tmp_path: Path) -> None:
    sequence, expected = _planted_repeat()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    proposing = RefinedBoundary(
        scaffold="contig",
        start=180,
        end=245,
        seed_id="seed-1",
        original_start=180,
        original_end=245,
    )
    untouched = RefinedBoundary(
        scaffold="contig",
        start=expected[0],
        end=expected[3],
        seed_id="seed-2",
        original_start=430,
        original_end=480,
    )
    mappings = {
        boundary.seed_id: SeedGeneMapping(
            seed_id=boundary.seed_id,
            scaffold="contig",
            seed_start=boundary.original_start,
            seed_end=boundary.original_end,
            flank_start_bp=0,
            flank_end_bp=len(sequence),
        )
        for boundary in (proposing, untouched)
    }

    refined = refine_boundaries_with_terminal_repeats(
        [proposing, untouched],
        raw_genome_path=genome_path,
        validated_markers=[_Marker("contig", 220, 230, "marker-1")],
        boundary_diamond_query=GenomeDiamondQuery(seed_gene_mappings=mappings),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
    )

    assert [(boundary.start, boundary.end) for boundary in refined] == [
        (untouched.start, untouched.end),
        (proposing.start, proposing.end),
    ]
    by_seed = {boundary.seed_id: boundary for boundary in refined}
    assert by_seed["seed-1"].tir_status == "shared_pair"
    assert by_seed["seed-1"].tir_boundary_override is False
    assert by_seed["seed-2"].tir_status == "shared_pair"
    assert by_seed["seed-2"].tir_boundary_override is False


def test_shared_pair_resolution_handles_restoration_collision() -> None:
    originals = [
        RefinedBoundary(scaffold="contig", start=start, end=end) for start, end in ((100, 200), (250, 350), (400, 500))
    ]

    def proposal(boundary: RefinedBoundary, start: int, end: int) -> RefinedBoundary:
        return replace(
            boundary,
            start=start,
            end=end,
            pre_tir_start=boundary.start,
            pre_tir_end=boundary.end,
            tir_present=True,
            tir_status="detected",
            tir_left_start=start,
            tir_left_end=start + 30,
            tir_right_start=end - 30,
            tir_right_end=end,
            tir_identity=1.0,
            tir_boundary_override=True,
            tir_alignment_length=30,
            tir_candidate_count=1,
        )

    proposed = [
        proposal(originals[0], 50, 550),
        proposal(originals[1], 50, 550),
        proposal(originals[2], 100, 200),
    ]

    resolved = terminal_repeats._resolve_partition_conflicts(
        originals,
        [[boundary] for boundary in proposed],
    )

    assert sorted((boundary.start, boundary.end) for boundary in resolved) == [
        (100, 200),
        (250, 350),
        (400, 500),
    ]
    assert {boundary.tir_status for boundary in resolved} == {"shared_pair"}
    assert not any(boundary.tir_boundary_override for boundary in resolved)


def test_phase2_applies_tir_after_merge_before_downstream_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    events: list[tuple[str, object]] = []
    raw_path = tmp_path / "raw.fna"
    raw_path.write_text(">contig\n" + "ACGT" * 200 + "\n")
    mappings = {
        seed_id: SeedGeneMapping(
            seed_id=seed_id,
            scaffold="contig",
            seed_start=start,
            seed_end=end,
            flank_start_bp=0,
            flank_end_bp=800,
        )
        for seed_id, start, end in (("seed-1", 100, 250), ("seed-2", 300, 450))
    }
    query = GenomeDiamondQuery(seed_gene_mappings=mappings)
    taxonomy = GeneTaxonomy(porf_id="tax", scaffold="contig", start=150, end=160)

    monkeypatch.setattr(phase2, "build_proteome_index", lambda _path: {"contig": []})
    monkeypatch.setattr(phase2, "_require_phase2b_gene_taxonomy_db", lambda *_args, **_kwargs: tmp_path / "db")
    monkeypatch.setattr(
        phase2,
        "_run_boundary_diamond",
        lambda **_kwargs: phase2._BoundaryEvidence(
            taxonomy_map={"tax": taxonomy},
            control_stats=None,
            query=query,
        ),
    )
    monkeypatch.setattr(boundary_refiner, "extend_seeds_by_genes", lambda seeds, *_args, **_kwargs: seeds)
    monkeypatch.setattr(taxonomy_seed_refiner, "validate_taxonomy_refinement_mode", lambda **_kwargs: None)
    monkeypatch.setattr(
        taxonomy_seed_refiner,
        "evaluate_taxonomy_seed_refinement",
        lambda **kwargs: SimpleNamespace(
            selected_seeds=tuple(kwargs["merged_seeds"]),
            intervention_counts=InterventionCounts(),
        ),
    )

    def fake_merge(boundaries: list[RefinedBoundary], **_kwargs: object) -> list[RefinedBoundary]:
        events.append(("merge", [(boundary.start, boundary.end) for boundary in boundaries]))
        return [
            replace(
                boundaries[0],
                start=80,
                end=470,
                original_start=100,
                original_end=450,
            )
        ]

    def fake_tir(boundaries: list[RefinedBoundary], **kwargs: object) -> list[RefinedBoundary]:
        events.append(("tir", (boundaries[0].start, boundaries[0].end, kwargs["raw_genome_path"])))
        boundary = boundaries[0]
        return [
            replace(
                boundary,
                start=90,
                end=460,
                pre_tir_start=boundary.start,
                pre_tir_end=boundary.end,
                tir_present=True,
                tir_status="detected",
                tir_left_start=90,
                tir_left_end=120,
                tir_right_start=430,
                tir_right_end=460,
                tir_identity=1.0,
                tir_boundary_override=True,
                tir_alignment_length=30,
                tir_candidate_count=1,
                tir_scan_start=0,
                tir_scan_end=800,
            )
        ]

    def fake_marker_floor(boundaries: list[RefinedBoundary], _markers: list[object]) -> int:
        events.append(("marker_floor", (boundaries[0].start, boundaries[0].end)))
        return 0

    def fake_composition(boundaries: list[RefinedBoundary], **_kwargs: object) -> None:
        events.append(("composition", (boundaries[0].start, boundaries[0].end)))

    monkeypatch.setattr(boundary_refiner, "merge_adjacent_viral_boundaries", fake_merge)
    monkeypatch.setattr(terminal_repeats, "refine_boundaries_with_terminal_repeats", fake_tir)
    monkeypatch.setattr(boundary_refiner, "annotate_boundaries_with_marker_floor", fake_marker_floor)
    monkeypatch.setattr(phase2, "_recalculate_boundary_composition", fake_composition)

    config = PipelineConfig().with_overrides(
        resume=False,
        gene_taxonomy_faa_db=tmp_path / "db",
        boundary_host_trim_enabled=False,
        boundary_taxonomy_ml_enabled=False,
        threads=1,
        gene_taxonomy_threads=1,
    )
    result = phase2._run_phase2_subflow(
        masked_path=tmp_path / "masked.fna",
        proteome_path=tmp_path / "proteome.faa",
        merged_seeds=[
            MergedSeed(scaffold="contig", start=100, end=250, seed_id="seed-1"),
            MergedSeed(scaffold="contig", start=300, end=450, seed_id="seed-2"),
        ],
        validated_markers=[object()],
        host_signature_model=None,
        output_dir=tmp_path,
        genome_id="genome",
        config=config,
        genome_start_time=0.0,
        logger=logging.getLogger(__name__),
        raw_genome_path=raw_path,
    )

    assert [(name, value) for name, value in events] == [
        ("merge", [(100, 250), (300, 450)]),
        ("tir", (80, 470, raw_path)),
        ("marker_floor", (90, 460)),
        ("composition", (90, 460)),
    ]
    assert (result.refined_boundaries[0].start, result.refined_boundaries[0].end) == (90, 460)
    assert (tmp_path / "phase2" / "refined_boundaries.bed").read_text().startswith("contig\t90\t460\t")
