"""Tests for the hallmark-composition merge rules of V9g (virosync.pipeline.eve_type_rules)."""

import json
from dataclasses import dataclass

import pytest

from virosync.pipeline import eve_type_rules as rules

CATALOGUE = rules.CATALOGUE


def model(panel: str, family: str) -> str:
    """Return a core model for this panel and family whose other panels share its group."""
    links = json.loads(CATALOGUE.read_text())["links"]
    for link in links:
        if link["stage"] == "phase1_marker" and (link["panel"], link["marker"]) == (panel, family):
            groups = {rules.GROUPS[p] for p, _ in rules.core_roles().get(link["model"], ())}
            if groups == {rules.GROUPS[panel]}:
                return link["model"]
    raise LookupError((panel, family))


def region(*pairs: tuple[str, str]) -> rules.Composition:
    return rules.composition((f"g{i}", model(*pair)) for i, pair in enumerate(pairs))


NCLDV, MIRUS, PLV, CRESS = "Nucleocytoviricota", "Mirusviricota", "Polinton_like_viruses", "CRESS"
POLINTON = "Maverick_Polintons"


def test_complementary_hallmarks_of_one_type_merge() -> None:
    assert rules.may_merge(region((NCLDV, "MCP")), region((NCLDV, "A32")))


def test_repeated_hallmark_without_full_sets_merges() -> None:
    assert rules.may_merge(region((NCLDV, "MCP"), (NCLDV, "A32")), region((NCLDV, "MCP"), (NCLDV, "PolB")))


def test_duplicated_mcp_genes_of_one_degraded_eve_merge() -> None:
    assert rules.may_merge(
        region((NCLDV, "MCP"), (NCLDV, "MCP"), (NCLDV, "SFII")), region((NCLDV, "MCP"), (NCLDV, "SFII"))
    )


def test_two_full_ncldv_sets_stay_separate() -> None:
    left = region((NCLDV, "MCP"), (NCLDV, "PolB"), (NCLDV, "RNAPL"))
    right = region((NCLDV, "MCP"), (NCLDV, "A32"), (NCLDV, "VLTF3"))

    assert not rules.may_merge(left, right)


def test_two_full_plv_sets_stay_separate() -> None:
    assert not rules.may_merge(region((PLV, "MCP"), (PLV, "ATPase")), region((PLV, "MCP"), (PLV, "ATPase")))


def test_plv_mcp_repeat_without_full_set_merges() -> None:
    assert rules.may_merge(region((PLV, "MCP")), region((PLV, "MCP"), (PLV, "ATPase")))


def test_complementary_ncldv_and_mirus_categories_merge() -> None:
    assert rules.may_merge(region((NCLDV, "A32"), (NCLDV, "PolB")), region((MIRUS, "MCP"), (MIRUS, "Portal")))


def test_lone_repeated_mcp_across_compatible_types_merges() -> None:
    assert rules.may_merge(region((NCLDV, "MCP")), region((MIRUS, "MCP")))


def test_full_ncldv_and_full_mirus_sets_stay_separate() -> None:
    left = region((NCLDV, "MCP"), (NCLDV, "PolB"), (NCLDV, "RNAPL"))
    right = region((MIRUS, "Portal"), (MIRUS, "Terminase"), (MIRUS, "RNAP"))

    assert not rules.may_merge(left, right)


def test_complementary_polinton_and_plv_categories_merge() -> None:
    assert rules.may_merge(region((POLINTON, "MCP")), region((PLV, "ATPase")))


def test_polinton_mcp_and_plv_mcp_merge() -> None:
    assert rules.may_merge(region((POLINTON, "MCP")), region((PLV, "MCP")))


@pytest.mark.parametrize("other", [(PLV, "ATPase"), (CRESS, "Rep")])
def test_incompatible_types_never_merge(other: tuple[str, str]) -> None:
    assert not rules.may_merge(region((NCLDV, "A32")), region(other))


def test_mcp_only_region_cannot_set_the_group() -> None:
    assert rules.may_merge(region((PLV, "MCP")), region((NCLDV, "A32"), (NCLDV, "PolB")))


def test_region_without_hallmarks_merges() -> None:
    assert rules.may_merge(region(), region((PLV, "MCP")))


def test_any_preplasmiviricota_types_with_complementary_categories_merge() -> None:
    assert rules.may_merge(region(("Virophages", "Penton")), region((PLV, "MCP")))
    assert rules.may_merge(region(("Virophages", "MCP")), region((POLINTON, "MCP")))


@dataclass
class Core:
    start: int
    end: int


@dataclass
class Boundary:
    scaffold: str
    start: int
    end: int
    confirmed_seeds: list[Core]


@dataclass
class Marker:
    scaffold: str
    start: int
    end: int
    query_porf: str
    hmm_target: str
    validation_status: str = "validated"


def marker(gene: str, start: int, family: str) -> Marker:
    return Marker("s", start, start + 900, gene, model(NCLDV, family))


def full_set(prefix: str, starts: tuple[int, int, int], families: tuple[str, str, str]) -> list[Marker]:
    return [
        marker(f"{prefix}{i}", start, family) for i, (start, family) in enumerate(zip(starts, families, strict=True))
    ]


def test_interleaved_seed_cores_merge_despite_two_full_sets() -> None:
    markers = full_set("a", (1000, 1500, 5000), ("MCP", "PolB", "RNAPL")) + full_set(
        "b", (3000, 7000, 7500), ("MCP", "A32", "VLTF3")
    )
    left = Boundary("s", 0, 9000, [Core(900, 2500), Core(4900, 6000)])
    right = Boundary("s", 2000, 9000, [Core(2900, 4000), Core(6900, 8500)])

    assert rules.boundary_compatibility(markers)(left, right)


def test_side_by_side_full_sets_stay_separate() -> None:
    markers = full_set("a", (1000, 1500, 2000), ("MCP", "PolB", "RNAPL")) + full_set(
        "b", (7000, 7500, 8000), ("MCP", "A32", "VLTF3")
    )
    left = Boundary("s", 0, 5000, [Core(900, 3000)])
    right = Boundary("s", 4000, 9500, [Core(6900, 9000)])

    assert not rules.boundary_compatibility(markers)(left, right)


def test_separate_overlap_cuts_between_seed_cores() -> None:
    left = Boundary("s", 0, 5000, [Core(900, 2000)])
    right = Boundary("s", 4000, 9000, [Core(6900, 8000)])

    cut_left, cut_right = rules.separate_overlap(left, right)

    assert (cut_left.start, cut_left.end, cut_right.start, cut_right.end) == (0, 4450, 4450, 9000)
