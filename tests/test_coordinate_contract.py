from __future__ import annotations

import ast
import importlib
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Lock
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq

from virosync.pipeline.phase0.prodigal import parse_prodigal_header
from virosync.pipeline.phase1.seed_merger import MergedSeed
from virosync.pipeline.phase2.boundary_diamond import (
    filter_taxonomy_to_boundary,
    pORF,
)
from virosync.pipeline.phase2.boundary_refiner import (
    RefinedBoundary,
    extend_seeds_by_genes,
    merge_adjacent_viral_boundaries,
)
from virosync.orchestration._flows.single_genome.phase3 import (
    _build_scaffold_start_index,
    _query_boundary_coordinate_records,
)


CODING_DNA = "ATG" + ("GCT" * 143)
EXPECTED_PROTEIN = "M" + ("A" * 143)
EXPECTED_DIRECT_CONSUMERS = {
    "orchestration/utils.py": 1,
    "pipeline/phase0/prodigal.py": 7,
    "pipeline/phase1/hhg_seeding.py": 1,
    "pipeline/phase1/marker_validation.py": 1,
    "pipeline/phase2/boundary_diamond.py": 1,
    "pipeline/phase3/gene_taxonomy.py": 3,
    "pipeline/phase3/interproscan.py": 1,
    "pipeline/phase3/output_generator.py": 1,
}


def test_phase3_boundary_indexes_match_ordered_half_open_scans() -> None:
    taxonomy_records = [
        SimpleNamespace(scaffold="S1", start=20, end=30, name="first"),
        SimpleNamespace(scaffold="S2", start=10, end=40, name="other"),
        SimpleNamespace(scaffold="S1", start=10, end=20, name="left_touch"),
        SimpleNamespace(scaffold="S1", start=20, end=20, name="zero"),
        SimpleNamespace(scaffold="S1", start=20, end=25, name="second_tie"),
        SimpleNamespace(scaffold="S1", start=30, end=40, name="right_touch"),
        SimpleNamespace(scaffold="S1", start=15, end=35, name="spanning"),
    ]
    markers = [
        SimpleNamespace(scaffold="S1", start=25, end=28, name="marker_first"),
        SimpleNamespace(scaffold="S2", start=21, end=29, name="marker_other"),
        SimpleNamespace(scaffold="S1", start=10, end=20, name="marker_touch"),
        SimpleNamespace(scaffold="S1", start=18, end=23, name="marker_second"),
    ]
    boundary = SimpleNamespace(scaffold="S1", start=20, end=30)
    taxonomy_map = {
        record.name: record for record in taxonomy_records
    }
    expected_taxonomy = filter_taxonomy_to_boundary(taxonomy_map, boundary)
    expected_markers = [
        marker
        for marker in markers
        if marker.scaffold == boundary.scaffold
        and marker.start < boundary.end
        and marker.end > boundary.start
    ]

    observed_taxonomy, observed_markers = _query_boundary_coordinate_records(
        boundary=boundary,
        taxonomy_index=_build_scaffold_start_index(taxonomy_map.values()),
        marker_index=_build_scaffold_start_index(markers),
    )

    assert observed_taxonomy == expected_taxonomy
    assert observed_markers == expected_markers


def _header(start: object, end: object, strand: object) -> str:
    return f"scaffold # {start} # {end} # {strand} # ID=scaffold_1;partial=00"


def test_plus_prodigal_interval_extracts_and_translates_exactly() -> None:
    genome = ("C" * 117) + CODING_DNA + ("G" * 31)

    parsed = parse_prodigal_header(_header(118, 549, 1), "scaffold_1")

    assert parsed == ("scaffold", 117, 549, "+")
    _scaffold, start, end, strand = parsed
    extracted = Seq(genome[start:end])
    assert strand == "+"
    assert len(extracted) == end - start == 432
    assert str(extracted.translate()) == EXPECTED_PROTEIN


def test_minus_prodigal_interval_extracts_and_translates_exactly() -> None:
    reverse_encoded = str(Seq(CODING_DNA).reverse_complement())
    genome = ("C" * 117) + reverse_encoded + ("G" * 31)

    parsed = parse_prodigal_header(_header(118, 549, -1), "scaffold_1")

    assert parsed == ("scaffold", 117, 549, "-")
    _scaffold, start, end, strand = parsed
    extracted = Seq(genome[start:end]).reverse_complement()
    assert strand == "-"
    assert len(extracted) == end - start == 432
    assert str(extracted.translate()) == EXPECTED_PROTEIN


@pytest.mark.parametrize(("token", "expected"), [("+", "+"), ("-", "-")])
def test_symbolic_prodigal_strands_are_normalized(
    token: str,
    expected: str,
) -> None:
    parsed = parse_prodigal_header(_header(1, 3, token), "scaffold_1")

    assert parsed == ("scaffold", 0, 3, expected)


def test_prodigal_start_below_one_is_rejected() -> None:
    with pytest.raises(ValueError, match="start"):
        parse_prodigal_header(_header(0, 3, 1), "scaffold_1")


def test_prodigal_end_before_start_is_rejected() -> None:
    with pytest.raises(ValueError, match="end"):
        parse_prodigal_header(_header(4, 3, 1), "scaffold_1")


@pytest.mark.parametrize(
    "header",
    [_header("not-an-int", 3, 1), _header(1, "not-an-int", 1)],
)
def test_noninteger_prodigal_coordinate_is_rejected(header: str) -> None:
    with pytest.raises(ValueError, match="coordinate"):
        parse_prodigal_header(header, "scaffold_1")


@pytest.mark.parametrize("strand", [0, "forward", ""])
def test_malformed_prodigal_strand_is_rejected(strand: object) -> None:
    with pytest.raises(ValueError, match="strand"):
        parse_prodigal_header(_header(1, 3, strand), "scaffold_1")


def test_structurally_unrelated_header_is_not_treated_as_prodigal() -> None:
    assert (
        parse_prodigal_header(
            "ordinary_protein description",
            "ordinary_protein",
        )
        is None
    )


@pytest.mark.parametrize(
    "header",
    ["scaffold # 1 # 3", " # 1 # 3 # 1"],
)
def test_structurally_malformed_prodigal_header_is_rejected(
    header: str,
) -> None:
    with pytest.raises(ValueError, match="Prodigal"):
        parse_prodigal_header(header, "scaffold_1")


def _is_parser_call(node: ast.AST) -> bool:
    return (
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "parse_prodigal_header"
    )


def _contains_parser_call(node: ast.AST) -> bool:
    return any(_is_parser_call(child) for child in ast.walk(node))


def _coordinate_names_from_parser(function: ast.AST) -> set[str]:
    parsed_names: set[str] = set()
    coordinate_names: set[str] = set()
    for node in ast.walk(function):
        if not isinstance(node, ast.Assign) or not _contains_parser_call(
            node.value
        ):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name):
                parsed_names.add(target.id)
            elif (
                isinstance(target, (ast.Tuple, ast.List))
                and len(target.elts) >= 3
            ):
                coordinate_names.update(
                    item.id
                    for item in target.elts[1:3]
                    if isinstance(item, ast.Name)
                )

    for node in ast.walk(function):
        if not isinstance(node, ast.Assign) or not isinstance(
            node.value, ast.Name
        ):
            continue
        if node.value.id not in parsed_names:
            continue
        for target in node.targets:
            if (
                isinstance(target, (ast.Tuple, ast.List))
                and len(target.elts) >= 3
            ):
                coordinate_names.update(
                    item.id
                    for item in target.elts[1:3]
                    if isinstance(item, ast.Name)
                )
    return coordinate_names


def _double_conversion_lines(function: ast.AST, names: set[str]) -> list[int]:
    lines: set[int] = set()
    for node in ast.walk(function):
        if isinstance(node, ast.AugAssign):
            if (
                isinstance(node.target, ast.Name)
                and node.target.id in names
                and isinstance(node.op, (ast.Add, ast.Sub))
            ):
                lines.add(node.lineno)
        if not isinstance(node, ast.BinOp):
            continue
        if not isinstance(node.left, ast.Name) or node.left.id not in names:
            continue
        if not isinstance(node.op, (ast.Add, ast.Sub)):
            continue
        if isinstance(node.right, ast.Constant) and node.right.value == 1:
            lines.add(node.lineno)
    return sorted(lines)


def _direct_consumer_inventory() -> tuple[
    dict[str, int],
    dict[str, list[int]],
]:
    source_root = Path(__file__).resolve().parents[1] / "src" / "virosync"
    inventory: dict[str, int] = {}
    second_conversions: dict[str, list[int]] = {}
    for path in sorted(source_root.rglob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        count = sum(1 for node in ast.walk(tree) if _is_parser_call(node))
        if not count:
            continue
        relative = path.relative_to(source_root).as_posix()
        inventory[relative] = count
        lines: list[int] = []
        for function in ast.walk(tree):
            if not isinstance(
                function, (ast.FunctionDef, ast.AsyncFunctionDef)
            ):
                continue
            if not any(_is_parser_call(node) for node in ast.walk(function)):
                continue
            names = _coordinate_names_from_parser(function)
            lines.extend(_double_conversion_lines(function, names))
        if lines:
            second_conversions[relative] = sorted(set(lines))
    return inventory, second_conversions


def test_direct_consumers_share_the_normalized_parser_contract() -> None:
    inventory, second_conversions = _direct_consumer_inventory()

    assert inventory == EXPECTED_DIRECT_CONSUMERS
    assert second_conversions == {}


def test_boundary_touching_genes_do_not_overlap(
    tmp_path: Path,
) -> None:
    import virosync.orchestration.utils as orchestration_utils

    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">pORF_scaffold_25_29_late\nLATE\n"
        ">pORF_scaffold_10_20_left\nLEFT\n"
        ">pORF_scaffold_20_25_inside\nINSIDE\n"
        ">pORF_scaffold_30_40_right\nRIGHT\n"
    )
    orchestration_utils._load_proteome_data.cache_clear()

    observed = orchestration_utils.get_overlapping_genes(
        proteome,
        boundary_scaffold="scaffold",
        boundary_start=20,
        boundary_end=30,
    )

    assert observed == {
        "scaffold": [
            ("pORF_scaffold_25_29_late", "LATE"),
            ("pORF_scaffold_20_25_inside", "INSIDE"),
        ]
    }
    orchestration_utils._load_proteome_data.cache_clear()


def test_proteome_cache_coalesces_fill_and_retains_four_genomes(
    monkeypatch,
    tmp_path: Path,
) -> None:
    import virosync.orchestration.utils as orchestration_utils

    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">pORF_scaffold_1_4_gene\nM\n")
    orchestration_utils._load_proteome_data.cache_clear()
    original_parse = orchestration_utils._parse_proteome_records
    call_lock = Lock()
    call_count = 0

    def counted_parse(*args):
        nonlocal call_count
        with call_lock:
            call_count += 1
        time.sleep(0.05)
        return original_parse(*args)

    monkeypatch.setattr(
        orchestration_utils,
        "_parse_proteome_records",
        counted_parse,
    )
    with ThreadPoolExecutor(max_workers=8) as executor:
        list(
            executor.map(
                lambda path: orchestration_utils._cached_proteome_data(path)[0],
                [proteome] * 8,
            )
        )

    assert call_count == 1
    other_proteomes = []
    for index in range(3):
        other = tmp_path / f"proteome_{index}.faa"
        other.write_text(f">pORF_scaffold_{index + 5}_{index + 8}_gene\nM\n")
        other_proteomes.append(other)
    for other in other_proteomes:
        orchestration_utils._cached_proteome_data(other)
    orchestration_utils._cached_proteome_data(proteome)
    assert call_count == 4
    orchestration_utils._load_proteome_data.cache_clear()


def test_one_base_internal_bed_gff_fasta_round_trip() -> None:
    genome = "ACGT"
    parsed = parse_prodigal_header(_header(1, 1, 1), "scaffold_1")
    assert parsed is not None
    scaffold, start, end, _strand = parsed
    boundary = RefinedBoundary(
        scaffold=scaffold,
        start=start,
        end=end,
        core_viral_start=start,
        core_viral_end=end,
    )

    bed = boundary.to_bed_line().split("\t")
    gff = boundary.to_gff_line().split("\t")
    bed_interval = (int(bed[1]), int(bed[2]))
    gff_native = (int(gff[3]), int(gff[4]))
    gff_internal = (gff_native[0] - 1, gff_native[1])
    attrs = dict(
        pair.split("=", 1) for pair in gff[8].split(";") if "=" in pair
    )
    extracted = genome[start:end]

    # This is ViroSync's normalized output contract. Phase 0 genes.gff remains
    # an upstream-native Prodigal GFF and is intentionally not rewritten here.
    assert (start, end) == (0, 1)
    assert bed_interval == (0, 1)
    assert gff_native == (1, 1)
    assert gff_internal == (0, 1)
    # GFF3 attributes carry the same 1-based inclusive coordinates as the
    # columns of their own record, not the internal 0-based half-open values.
    assert (int(attrs["core_start"]), int(attrs["core_end"])) == gff_native
    assert extracted == "A"
    assert len(extracted) == end - start == 1


def test_coordinate_convention_is_versioned_in_output_contract() -> None:
    contract = importlib.import_module("virosync.output_contract")

    assert contract.COORDINATE_SCHEMA_VERSION == 2
    assert contract.OUTPUT_SCHEMA_VERSION == 6
    assert "0-based" in contract.COORDINATE_CONVENTION.lower()
    assert "half-open" in contract.COORDINATE_CONVENTION.lower()


def test_touching_extended_seeds_remain_distinct() -> None:
    seeds = [
        MergedSeed(scaffold="scaffold", start=0, end=10, seed_id="left"),
        MergedSeed(scaffold="scaffold", start=10, end=20, seed_id="right"),
    ]
    proteome_index = {
        "scaffold": [
            pORF(id="gene_left", scaffold="scaffold", start=0, end=10),
            pORF(id="gene_right", scaffold="scaffold", start=10, end=20),
        ]
    }

    observed = extend_seeds_by_genes(seeds, proteome_index, extension_genes=0)

    assert [(seed.start, seed.end) for seed in observed] == [(0, 10), (10, 20)]


def test_gene_extension_does_not_merge_mixed_rescue_and_ordinary_seeds() -> None:
    seeds = [
        MergedSeed(
            scaffold="scaffold",
            start=0,
            end=20,
            sources=["hhg", "marker_validation"],
        ),
        MergedSeed(
            scaffold="scaffold",
            start=10,
            end=30,
            sources=["hhg", "marker_validation", "frameshift_rescue"],
        ),
    ]
    proteome_index = {
        "scaffold": [
            pORF(id="gene_left", scaffold="scaffold", start=0, end=20),
            pORF(id="gene_right", scaffold="scaffold", start=10, end=30),
        ]
    }

    observed = extend_seeds_by_genes(seeds, proteome_index, extension_genes=0)

    assert [(seed.start, seed.end) for seed in observed] == [(0, 30), (0, 30)]
    assert ["frameshift_rescue" in seed.sources for seed in observed] == [False, True]


def test_touching_refined_boundaries_are_not_unconditional_overlaps() -> None:
    boundaries = [
        RefinedBoundary(scaffold="scaffold", start=0, end=10),
        RefinedBoundary(scaffold="scaffold", start=10, end=20),
    ]

    observed = merge_adjacent_viral_boundaries(boundaries, taxonomy_map={})

    assert [(boundary.start, boundary.end) for boundary in observed] == [
        (0, 10),
        (10, 20),
    ]


def test_post_taxonomy_merge_keeps_overlapping_mixed_rescue_boundaries_separate() -> None:
    boundaries = [
        RefinedBoundary(
            scaffold="scaffold",
            start=0,
            end=20,
            seed_sources=["hhg", "marker_validation"],
        ),
        RefinedBoundary(
            scaffold="scaffold",
            start=10,
            end=30,
            seed_sources=["hhg", "marker_validation", "frameshift_rescue"],
        ),
    ]

    observed = merge_adjacent_viral_boundaries(boundaries, taxonomy_map={})

    assert [(boundary.start, boundary.end) for boundary in observed] == [
        (0, 20),
        (10, 30),
    ]


def test_post_taxonomy_merge_keeps_viral_gap_mixed_rescue_boundaries_separate() -> None:
    boundaries = [
        RefinedBoundary(
            scaffold="scaffold",
            start=0,
            end=100,
            seed_sources=["hhg", "marker_validation"],
        ),
        RefinedBoundary(
            scaffold="scaffold",
            start=200,
            end=300,
            seed_sources=["hhg", "marker_validation", "frameshift_rescue"],
        ),
    ]
    taxonomy_map = {
        "gap_gene": SimpleNamespace(
            scaffold="scaffold",
            start=120,
            end=180,
            has_ncldv_mirus=True,
            has_vp_plv=False,
        )
    }

    observed = merge_adjacent_viral_boundaries(boundaries, taxonomy_map)

    assert [(boundary.start, boundary.end) for boundary in observed] == [
        (0, 100),
        (200, 300),
    ]
