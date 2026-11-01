from __future__ import annotations

import logging
import random
from dataclasses import dataclass, replace
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio import SeqIO

from virosync.config import PipelineConfig
from virosync.orchestration._flows.single_genome import phase2
from virosync.orchestration._flows.single_genome.phase_state import load_phase2_state, write_phase2_state
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
from virosync.pipeline.phase2 import terminal_repeats
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
    ).candidates

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
    ).candidates

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
    ).candidates

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
        ).candidates
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
    ).candidates[0]

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
    ).candidates

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
            ).candidates
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
        ).candidates
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
        ).candidates
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

    uncached = find_terminal_inverted_repeats(sequence, **arguments).candidates
    cached = find_terminal_inverted_repeats(
        sequence,
        **arguments,
        _scan_seed_index=terminal_repeats._build_scan_seed_index(sequence[scan_start:scan_end]),
    ).candidates

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
    ).candidates

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
    ("boundary_start", "boundary_end", "expected_bounds", "expected_status", "expected_override", "core_bounds"),
    [
        (120, 400, (120, 400), "ambiguous", False, (150, 350)),
        (50, 470, (100, 420), "detected", True, (150, 350)),
        (50, 470, (50, 470), "ambiguous", False, (50, 470)),
    ],
)
def test_unique_repeat_changes_only_reachable_noncore_flanks(
    tmp_path: Path,
    boundary_start: int,
    boundary_end: int,
    expected_bounds: tuple[int, int],
    expected_status: str,
    expected_override: bool,
    core_bounds: tuple[int, int],
) -> None:
    sequence, _ = _planted_repeat()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    boundary = RefinedBoundary(
        scaffold="contig",
        start=boundary_start,
        end=boundary_end,
        seed_id="seed-1",
        original_start=core_bounds[0],
        original_end=core_bounds[1],
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", *core_bounds),),
    )

    refined = refine_boundaries_with_terminal_repeats(
        [boundary],
        allowed_intervals=[(boundary.start, boundary.end)],
        raw_genome_path=genome_path,
        validated_markers=[_Marker("contig", 250, 270, "marker-1")],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
    )[0]

    assert (refined.start, refined.end) == expected_bounds
    assert (refined.pre_tir_start, refined.pre_tir_end) == (boundary_start, boundary_end)
    assert refined.tir_status == expected_status
    assert refined.tir_present is True
    assert refined.tir_boundary_override is expected_override
    assert refined.tir_alignment_capped is False
    assert refined.tir_candidate_count == 1


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
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 150, 350),),
    )

    refined = refine_boundaries_with_terminal_repeats(
        [boundary],
        allowed_intervals=[(boundary.start, boundary.end)],
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
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 150, 350),),
    )

    refined = refine_boundaries_with_terminal_repeats(
        [boundary],
        allowed_intervals=[(boundary.start, boundary.end)],
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


@pytest.mark.parametrize("parent_start", [0, 120])
def test_disjoint_tir_pairs_partition_parent_and_keep_marker_neighbor(
    tmp_path: Path,
    parent_start: int,
) -> None:
    sequence, expected_pairs = _adjacent_repeat_sequence()
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    parent = RefinedBoundary(
        scaffold="contig",
        start=parent_start,
        end=len(sequence),
        seed_id="seed-1",
        original_start=200,
        original_end=1100,
        confirmed_seeds=(
            ConfirmedSeed("seed-1", "contig", 200, 260),
            ConfirmedSeed("seed-2", "contig", 620, 680),
            ConfirmedSeed("seed-3", "contig", 980, 1100),
        ),
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
        allowed_intervals=[(0, len(sequence))],
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
    assert [(boundary.original_start, boundary.original_end) for boundary in refined] == [
        (200, 260),
        (620, 680),
        (980, 1100),
    ]
    assert [boundary.confirmed_seeds[0].seed_id for boundary in refined] == ["seed-1", "seed-2", "seed-3"]
    assert all((boundary.pre_tir_start, boundary.pre_tir_end) == (parent_start, len(sequence)) for boundary in refined)
    assert [boundary.seed_id for boundary in refined] == ["seed-1", "seed-2", "seed-3"]
    assert all(boundary.confidence == 0.0 for boundary in refined)
    assert all(boundary.posterior_probability == 0.0 for boundary in refined)
    assert all(boundary.cub_deviation == 0.0 for boundary in refined)


def test_partition_cannot_split_whole_seed_for_viral_taxonomy_neighbor(tmp_path: Path) -> None:
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
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 0, len(sequence)),),
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
        allowed_intervals=[(parent.start, parent.end)],
        raw_genome_path=genome_path,
        validated_markers=markers,
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={viral_gene.id: taxonomy},
        proteome_index={"contig": [viral_gene]},
        extension_bp=500,
    )

    assert len(refined) == 1
    assert (refined[0].start, refined[0].end) == (parent.start, parent.end)
    assert refined[0].confirmed_seeds == parent.confirmed_seeds
    assert refined[0].tir_status == "ambiguous"
    assert refined[0].tir_present is True
    assert refined[0].tir_boundary_override is False


def test_partition_retains_parent_when_marker_neighbor_has_no_confirmed_seed(
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
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 200, 300),),
    )

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        allowed_intervals=[(parent.start, parent.end)],
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

    assert len(refined) == 1
    assert (refined[0].start, refined[0].end) == (parent.start, parent.end)
    assert refined[0].confirmed_seeds == parent.confirmed_seeds
    assert refined[0].tir_status == "ambiguous"
    assert refined[0].tir_present is True
    assert refined[0].tir_boundary_override is False


def test_partition_cannot_emit_unseeded_call_beyond_parent(tmp_path: Path) -> None:
    planted, _ = _planted_repeat()
    sequence = _random_dna(400, random.Random(503)) + planted
    genome_path = tmp_path / "raw.fna"
    genome_path.write_text(f">contig\n{sequence}\n")
    parent = RefinedBoundary(
        scaffold="contig",
        start=100,
        end=800,
        seed_id="seed-1",
        original_start=100,
        original_end=800,
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 100, 800),),
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
        allowed_intervals=[(parent.start, parent.end)],
        raw_genome_path=genome_path,
        validated_markers=[_Marker("contig", 620, 640, "tir-marker")],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={viral_gene.id: taxonomy},
        proteome_index={"contig": [viral_gene]},
        extension_bp=500,
    )

    assert len(refined) == 1
    assert (refined[0].start, refined[0].end) == (parent.start, parent.end)
    assert refined[0].confirmed_seeds == parent.confirmed_seeds
    assert refined[0].tir_status == "ambiguous"
    assert refined[0].tir_present is True
    assert refined[0].tir_boundary_override is False


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
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 0, len(sequence)),),
    )

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        allowed_intervals=[(parent.start, parent.end)],
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
        original_start=200,
        original_end=1100,
        confirmed_seeds=(
            ConfirmedSeed("seed-1", "contig", 200, 260),
            ConfirmedSeed("seed-2", "contig", 620, 680),
            ConfirmedSeed("seed-3", "contig", 980, 1100),
        ),
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
        allowed_intervals=[(parent.start, parent.end)],
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
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 0, len(sequence)),),
    )

    refined = refine_boundaries_with_terminal_repeats(
        [parent],
        allowed_intervals=[(parent.start, parent.end)],
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


@pytest.mark.parametrize(
    ("allowed", "status"), [([(180, 245), (260, 340)], "ambiguous"), ([(0, 520), (0, 520)], "shared_pair")]
)
def test_shared_pair_does_not_collapse_distinct_boundaries(
    tmp_path: Path, allowed: list[tuple[int, int]], status: str
) -> None:
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
            confirmed_seeds=(ConfirmedSeed(seed_id, "contig", start, end),),
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
        allowed_intervals=allowed,
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
    assert {boundary.tir_status for boundary in refined} == {status}
    assert all(boundary.tir_present for boundary in refined)
    assert not any(boundary.tir_boundary_override for boundary in refined)
    assert {(boundary.tir_left_start, boundary.tir_right_end) for boundary in refined} == {(expected[0], expected[3])}


@pytest.mark.parametrize(("allowed_first", "status"), [((180, 245), "ambiguous"), ((0, 520), "shared_pair")])
def test_outward_tir_proposal_does_not_replace_existing_candidate(
    tmp_path: Path, allowed_first: tuple[int, int], status: str
) -> None:
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
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 180, 245),),
    )
    untouched = RefinedBoundary(
        scaffold="contig",
        start=expected[0],
        end=expected[3],
        seed_id="seed-2",
        original_start=150,
        original_end=200,
        confirmed_seeds=(ConfirmedSeed("seed-2", "contig", 150, 200),),
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
        allowed_intervals=[allowed_first, (untouched.start, untouched.end)],
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
    assert by_seed["seed-1"].tir_status == status
    assert by_seed["seed-1"].tir_boundary_override is False
    assert by_seed["seed-2"].tir_status == "detected"
    assert by_seed["seed-2"].tir_boundary_override is False

    state_path = tmp_path / "refined_state.json"
    write_phase2_state(state_path, refined)
    assert load_phase2_state(state_path) == refined


def test_shared_pair_resolution_handles_restoration_collision() -> None:
    originals = [
        RefinedBoundary(
            scaffold="contig",
            start=start,
            end=end,
            seed_id=f"source-{start}",
            confirmed_seeds=(ConfirmedSeed(f"source-{start}", "contig", start + 10, end - 10),),
        )
        for start, end in ((100, 200), (250, 350), (400, 500))
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


@pytest.mark.parametrize(("stops", "expected_allowed"), [((None, None), (0, 700)), ((50, 650), (50, 650))])
def test_phase2_applies_tir_after_merge_before_downstream_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    stops: tuple[int | None, int | None],
    expected_allowed: tuple[int, int],
) -> None:
    events: list[tuple[str, object]] = []
    retained_marker = SimpleNamespace(
        scaffold="contig",
        start=150,
        end=160,
        query_porf="tax",
        validation_status="validated",
        hmm_target="GVOGm0003",
        hmm_score=100.0,
        hmm_evalue=1e-20,
        top10_prefixes="NCLDV__",
        top10_pidents="50.0",
        top10_targets="",
    )
    raw_path = tmp_path / "raw.fna"
    raw_path.write_text(">contig\n" + "ACGT" * 200 + "\n")
    mappings = {
        seed_id: SeedGeneMapping(
            seed_id=seed_id,
            scaffold="contig",
            seed_start=start,
            seed_end=end,
            flank_start_bp=0 if seed_id == "seed-1" else 200,
            flank_end_bp=600 if seed_id == "seed-1" else 700,
        )
        for seed_id, start, end in (("seed-1", 100, 250), ("seed-2", 300, 450))
    }
    mappings["separate"] = SeedGeneMapping("separate", "contig", 200, 300, flank_start_bp=0, flank_end_bp=800)
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

    def fake_merge(boundaries: list[RefinedBoundary], **_kwargs: object) -> list[RefinedBoundary]:
        events.append(("merge", [(boundary.start, boundary.end) for boundary in boundaries]))
        return [
            replace(
                boundaries[0],
                start=80,
                end=470,
                original_start=100,
                original_end=450,
                confirmed_seeds=boundaries[0].confirmed_seeds + boundaries[1].confirmed_seeds,
                left_stop=stops[0],
                right_stop=stops[1],
            )
        ]

    def fake_tir(boundaries: list[RefinedBoundary], **kwargs: object) -> list[RefinedBoundary]:
        assert kwargs["validated_markers"][0] is retained_marker
        assert kwargs["allowed_intervals"] == [expected_allowed]
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

    def fake_composition(boundaries: list[RefinedBoundary], **_kwargs: object) -> None:
        events.append(("composition", (boundaries[0].start, boundaries[0].end)))

    monkeypatch.setattr(phase2, "merge_reachable_boundaries", fake_merge)
    monkeypatch.setattr(terminal_repeats, "refine_boundaries_with_terminal_repeats", fake_tir)
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
            MergedSeed(
                scaffold="contig",
                start=100,
                end=250,
                seed_id="seed-1",
                confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 100, 250),),
            ),
            MergedSeed(
                scaffold="contig",
                start=300,
                end=450,
                seed_id="seed-2",
                confirmed_seeds=(ConfirmedSeed("seed-2", "contig", 300, 450),),
            ),
        ],
        validated_markers=[retained_marker],
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
        ("composition", (90, 460)),
    ]
    assert (result.refined_boundaries[0].start, result.refined_boundaries[0].end) == (90, 460)
    assert (tmp_path / "phase2" / "refined_boundaries.bed").read_text().startswith("contig\t90\t460\t")


@pytest.mark.parametrize(
    ("parent_bounds", "allowed", "stop", "expected"),
    [
        ((120, 400), (0, 520), None, (100, 420)),
        ((100, 400), (0, 520), None, (100, 420)),
        ((50, 400), (0, 520), None, (100, 420)),
        ((120, 400), (101, 520), None, (120, 400)),
        ((120, 400), (100, 420), None, (100, 420)),
        ((120, 400), (0, 520), 410, (120, 400)),
        ((120, 400), (0, 520), 420, (100, 420)),
    ],
)
def test_outward_tir_respects_walked_interval_and_stop(
    tmp_path: Path,
    parent_bounds: tuple[int, int],
    allowed: tuple[int, int],
    stop: int | None,
    expected: tuple[int, int],
) -> None:
    sequence, _ = _planted_repeat()
    genome = tmp_path / "raw.fna"
    genome.write_text(f">contig\n{sequence}\n", encoding="utf-8")
    parent = RefinedBoundary(
        "contig",
        *parent_bounds,
        seed_id="seed-1",
        original_start=150,
        original_end=350,
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 150, 350),),
        right_stop=stop,
    )
    query = _query(len(sequence))

    result = refine_boundaries_with_terminal_repeats(
        [parent],
        allowed_intervals=[allowed],
        raw_genome_path=genome,
        validated_markers=[_Marker("contig", 250, 270, "marker-1")],
        boundary_diamond_query=query,
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
    )
    phase2._annotate_boundary_repeats(result, raw_genome_path=genome, boundary_diamond_query=query, extension_bp=500)
    state = tmp_path / "refined.json"
    write_phase2_state(state, result)
    loaded = load_phase2_state(state)

    assert (result[0].start, result[0].end) == expected
    assert result[0].confirmed_seeds == parent.confirmed_seeds
    assert (result[0].pre_tir_start, result[0].pre_tir_end) == parent_bounds
    assert result[0].repeat_evidence is not None
    assert (result[0].repeat_evidence.assessed_start, result[0].repeat_evidence.assessed_end) == expected
    assert (result[0].repeat_evidence.parent_start, result[0].repeat_evidence.parent_end) == parent_bounds
    assert loaded == result


@pytest.mark.parametrize(
    ("parent_bounds", "marker_span", "expected"),
    [
        ((120, 400), (410, 450), (120, 400, "ambiguous")),
        ((120, 400), (410, 700), (120, 400, "ambiguous")),
        ((120, 400), (0, 110), (120, 400, "ambiguous")),
        ((50, 470), (410, 700), (100, 420, "detected")),
        ((100, 420), (410, 700), (100, 420, "detected")),
    ],
)
def test_flank_marker_veto_applies_only_to_new_outward_endpoints(
    tmp_path: Path,
    parent_bounds: tuple[int, int],
    marker_span: tuple[int, int],
    expected: tuple[int, int, str],
) -> None:
    sequence, _ = _planted_repeat()
    sequence += "A" * 300
    genome = tmp_path / "raw.fna"
    genome.write_text(f">contig\n{sequence}\n", encoding="utf-8")
    parent = RefinedBoundary(
        "contig",
        *parent_bounds,
        seed_id="seed-1",
        original_start=150,
        original_end=350,
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 150, 350),),
    )
    result = refine_boundaries_with_terminal_repeats(
        [parent],
        allowed_intervals=[(0, 520)],
        raw_genome_path=genome,
        validated_markers=[_Marker("contig", 250, 270, "anchor"), _Marker("contig", *marker_span, "flank")],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
    )
    assert (result[0].start, result[0].end, result[0].tir_status) == expected


def _outward_proposal(parent: RefinedBoundary, start: int, end: int) -> RefinedBoundary:
    return replace(
        parent,
        start=start,
        end=end,
        pre_tir_start=parent.start,
        pre_tir_end=parent.end,
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


@pytest.mark.parametrize(("end", "expected_end", "status"), [(300, 200, "ambiguous"), (250, 250, "detected")])
@pytest.mark.parametrize("neighbor_bounds", [(250, 400), (260, 395)])
def test_new_partial_overlap_reverts_but_touching_neighbor_is_allowed(
    end: int, expected_end: int, status: str, neighbor_bounds: tuple[int, int]
) -> None:
    left = RefinedBoundary("s", 100, 200, seed_id="left", confirmed_seeds=(ConfirmedSeed("left", "s", 120, 180),))
    right = RefinedBoundary("s", 250, 400, seed_id="right", confirmed_seeds=(ConfirmedSeed("right", "s", 260, 390),))
    neighbor = right if neighbor_bounds == (250, 400) else _outward_proposal(right, *neighbor_bounds)
    result = terminal_repeats._resolve_partition_conflicts(
        [left, right], [[_outward_proposal(left, 100, end)], [neighbor]]
    )
    assert (result[0].end, result[0].tir_status) == (expected_end, status)
    assert result[1] == neighbor
    assert result[0].confirmed_seeds == left.confirmed_seeds


def test_partial_overlap_restoration_exposes_another_conflict() -> None:
    parents = [
        RefinedBoundary(
            "s", start, end, seed_id=str(start), confirmed_seeds=(ConfirmedSeed(str(start), "s", start + 10, end - 10),)
        )
        for start, end in [(100, 200), (250, 350), (20, 80)]
    ]
    parents[0].confirmed_seeds = (ConfirmedSeed("100", "s", 160, 180),)
    plans = [[_outward_proposal(parents[0], 150, 280)], [parents[1]], [_outward_proposal(parents[2], 20, 140)]]
    result = terminal_repeats._resolve_partition_conflicts(parents, plans)
    assert [(boundary.start, boundary.end) for boundary in result] == [(20, 80), (100, 200), (250, 350)]
    assert [boundary.tir_status for boundary in result] == ["ambiguous", "ambiguous", "not_assessed"]


@pytest.mark.parametrize(
    ("mode", "expected_allowed", "expected_bounds"),
    [
        ("no_host", (0, 520), (100, 420)),
        ("CRESS", (150, 350), (150, 350)),
        ("A4", (0, 520), (100, 420)),
        ("host", (0, 360), (150, 350)),
    ],
)
def test_phase2_real_tir_uses_walked_authority_and_final_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mode: str,
    expected_allowed: tuple[int, int],
    expected_bounds: tuple[int, int],
) -> None:
    from virosync.pipeline.host_signatures import HostSignatureModel
    from virosync.pipeline.phase2 import seed_extension

    sequence, _ = _planted_repeat()
    genome = tmp_path / "raw.fna"
    genome.write_text(f">contig\n{sequence}\n", encoding="utf-8")
    genes = [pORF(f"host-{start}", "contig", start, start + 10) for start in (360, 380, 400)]
    taxonomy = {
        gene.id: GeneTaxonomy(porf_id=gene.id, scaffold="contig", start=gene.start, end=gene.end) for gene in genes
    }
    query = _query(len(sequence))
    marker = SimpleNamespace(
        scaffold="contig",
        start=250,
        end=270,
        query_porf="anchor",
        validation_status="validated",
        hmm_target="GVOGm0003",
        hmm_score=100.0,
        hmm_evalue=1e-20,
        top10_prefixes="NCLDV__",
        top10_pidents="50.0",
        top10_targets="",
    )
    monkeypatch.setattr(phase2, "build_proteome_index", lambda _path: {"contig": genes})
    monkeypatch.setattr(phase2, "_require_phase2b_gene_taxonomy_db", lambda *_args, **_kwargs: tmp_path / "db")
    monkeypatch.setattr(
        phase2,
        "_run_boundary_diamond",
        lambda **_kwargs: phase2._BoundaryEvidence(taxonomy_map=taxonomy, control_stats=None, query=query),
    )
    monkeypatch.setattr(seed_extension, "host_lineage_share", lambda *_args: 1.0)
    native_refiner = terminal_repeats.refine_boundaries_with_terminal_repeats
    native_composition = phase2._recalculate_boundary_composition
    seen = []

    def check_refiner(boundaries, **kwargs):
        assert kwargs["allowed_intervals"] == [expected_allowed]
        return native_refiner(boundaries, **kwargs)

    def check_composition(boundaries, **kwargs):
        seen.extend((boundary.start, boundary.end) for boundary in boundaries)
        native_composition(boundaries, **kwargs)

    monkeypatch.setattr(terminal_repeats, "refine_boundaries_with_terminal_repeats", check_refiner)
    monkeypatch.setattr(phase2, "_recalculate_boundary_composition", check_composition)
    config = PipelineConfig().with_overrides(
        resume=False,
        gene_taxonomy_faa_db=tmp_path / "db",
        boundary_host_trim_enabled=mode != "no_host",
        boundary_taxonomy_ml_enabled=False,
        ablation_id="A4" if mode == "A4" else "A0",
        threads=1,
        gene_taxonomy_threads=1,
    )
    core = ConfirmedSeed("seed-1", "contig", 150, 350)
    result = phase2._run_phase2_subflow(
        masked_path=genome,
        raw_genome_path=genome,
        proteome_path=tmp_path / "proteome.faa",
        merged_seeds=[
            MergedSeed(
                "contig",
                150,
                350,
                seed_id="seed-1",
                confirmed_seeds=(core,),
                predicted_family="CRESS" if mode == "CRESS" else "NCLDV",
            )
        ],
        validated_markers=[marker],
        host_signature_model=HostSignatureModel(token_weights={"host": 1.0}),
        output_dir=tmp_path,
        genome_id="genome",
        config=config,
        genome_start_time=0.0,
        logger=logging.getLogger(__name__),
    )
    final = result.refined_boundaries[0]
    assert (final.start, final.end) == expected_bounds
    assert seen == [expected_bounds]
    assert final.confirmed_seeds == (core,)
    assert final.region_classification_ncldv_markers == 1
    assert final.candidate_id == f"EVE_contig_{expected_bounds[0]}-{expected_bounds[1]}"
    assert final.repeat_evidence is not None
    assert (final.repeat_evidence.assessed_start, final.repeat_evidence.assessed_end) == expected_bounds
    assert (final.repeat_evidence.parent_start, final.repeat_evidence.parent_end) == (150, 350)
    assert load_phase2_state(tmp_path / "phase2/refined_state.json") == result.refined_boundaries


def test_outward_tir_missing_taxonomy_in_added_flank_vetoes_change(tmp_path: Path) -> None:
    sequence, _ = _planted_repeat()
    genome = tmp_path / "raw.fna"
    genome.write_text(f">contig\n{sequence}\n", encoding="utf-8")
    parent = RefinedBoundary(
        "contig",
        120,
        400,
        seed_id="seed-1",
        original_start=150,
        original_end=350,
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 150, 350),),
    )
    result = refine_boundaries_with_terminal_repeats(
        [parent],
        allowed_intervals=[(0, 520)],
        raw_genome_path=genome,
        validated_markers=[_Marker("contig", 250, 270, "anchor")],
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": [pORF("missing", "contig", 405, 415)]},
        extension_bp=500,
    )
    assert (result[0].start, result[0].end, result[0].tir_status) == (120, 400, "taxonomy_incomplete")


def test_shared_pair_restoration_retains_identical_not_partial_collision_evidence() -> None:
    first = ConfirmedSeed("a1", "s", 220, 260)
    second = ConfirmedSeed("a2", "s", 520, 560)
    parent = RefinedBoundary("s", 200, 600, seed_id="a1", confirmed_seeds=(first, second))
    neighbor = RefinedBoundary("s", 50, 150, seed_id="b", confirmed_seeds=(ConfirmedSeed("b", "s", 70, 90),))
    peer = RefinedBoundary("s", 620, 680, seed_id="c", confirmed_seeds=(ConfirmedSeed("c", "s", 630, 650),))
    partial = replace(
        _outward_proposal(parent, 100, 300), confirmed_seeds=(first,), tir_identity=0.93, tsd_sequence="TGCA"
    )
    identical = replace(
        _outward_proposal(parent, 500, 700),
        seed_id="a2",
        confirmed_seeds=(second,),
        tir_identity=0.97,
        tsd_sequence="ACGT",
    )

    resolved = terminal_repeats._resolve_partition_conflicts(
        [parent, neighbor, peer],
        [[partial, identical], [neighbor], [_outward_proposal(peer, 500, 700)]],
    )
    restored = next(boundary for boundary in resolved if boundary.seed_id == "a1")

    assert (restored.start, restored.end, restored.tir_status) == (200, 600, "shared_pair")
    assert restored.confirmed_seeds == (first, second)
    assert (
        restored.tir_left_start,
        restored.tir_left_end,
        restored.tir_right_start,
        restored.tir_right_end,
        restored.tir_identity,
        restored.tsd_sequence,
    ) == (500, 530, 670, 700, 0.97, "ACGT")
