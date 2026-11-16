"""Decide from hallmark composition whether two EVE regions may merge (experimental V9).

Each validated Phase 1 marker is mapped through the hallmark panel catalogue to the
(panel, core family) pairs of its model. Regions in different groups never merge
(giant viruses: NCLDV, Mriya, MIRUS; Preplasmiviricota: PLV, Maverick/Polinton, virophage;
CRESS), judged from non-MCP hallmarks. Within a group, regions stay separate only when each
carries a full set of functional categories, because one degraded EVE often carries
duplicated hallmark genes. Regions whose seed cores interleave are one EVE and always
merge; refused regions that overlap are cut apart midway between their seed cores.
"""

import json
from collections.abc import Callable, Iterable
from dataclasses import replace
from functools import cache
from pathlib import Path
from typing import TypeAlias

from virosync.pipeline.phase1.marker_roles import is_ordinary_marker_evidence

CATALOGUE = Path(__file__).resolve().parents[1] / "data" / "hallmark_panels.json"
GROUPS = {
    "Nucleocytoviricota": "giant",
    "Mriyaviruses": "giant",
    "Mirusviricota": "giant",
    "Polinton_like_viruses": "preplasmiviricota",
    "Maverick_Polintons": "preplasmiviricota",
    "Virophages": "preplasmiviricota",
    "CRESS": "cress",
}

Composition: TypeAlias = dict[str, set[tuple[str, str]]]


@cache
def core_roles() -> dict[str, tuple[tuple[str, str], ...]]:
    """Return the (panel, core family) pairs of every Phase 1 marker model."""
    roles: dict[str, list[tuple[str, str]]] = {}
    for link in json.loads(CATALOGUE.read_text(encoding="utf-8"))["links"]:
        if (
            link["stage"] == "phase1_marker"
            and link["role"] == "core_marker"
            and not link["excluded_from_scoring"]
            and link["panel"] in GROUPS
        ):
            roles.setdefault(link["model"], []).append((link["panel"], link["marker"]))
    return {model: tuple(pairs) for model, pairs in roles.items()}


def composition(markers: Iterable[tuple[str, str]]) -> Composition:
    """Map each hallmark gene to its (panel, family) pairs from (gene id, model) pairs."""
    genes: Composition = {}
    for gene, model in markers:
        for pair in core_roles().get(model, ()):
            genes.setdefault(gene, set()).add(pair)
    return genes


def may_merge(left: Composition, right: Composition) -> bool:
    """Return whether two regions with these hallmark compositions may merge.

    Regions whose non-MCP hallmarks place them in different groups never merge; MCP genes
    alone cannot set the group, because double-jelly-roll MCP models cross-hit between giant
    viruses and Polinton-like viruses. Within a group, regions stay separate only when each
    carries a full set of functional categories (replication, structure, packaging,
    transcription): at least 3, or every category the group's core families cover when that
    is fewer.
    """
    if not left or not right:
        return True
    left_groups, right_groups = _groups(left), _groups(right)
    if left_groups and right_groups and not left_groups & right_groups:
        return False
    return not (_is_full_set(left) and _is_full_set(right))


def boundary_compatibility(validated_markers: Iterable[object]) -> Callable[[object, object], bool]:
    """Return a merge predicate over refined boundaries from the validated Phase 1 markers."""
    by_scaffold: dict[str, list[tuple[int, int, str, str]]] = {}
    for marker in validated_markers:
        if not is_ordinary_marker_evidence(marker) or getattr(marker, "validation_status", "validated") != "validated":
            continue
        gene = str(getattr(marker, "query_porf", "")).split("|", 1)[0]
        by_scaffold.setdefault(marker.scaffold, []).append(
            (int(marker.start), int(marker.end), gene, str(marker.hmm_target))
        )

    def region(boundary: object) -> Composition:
        cores = [(core.start, core.end) for core in boundary.confirmed_seeds]
        return composition(
            (gene, model)
            for start, end, gene, model in by_scaffold.get(boundary.scaffold, [])
            if any(core_start <= (start + end) // 2 < core_end for core_start, core_end in cores)
        )

    def compatible(previous: object, boundary: object) -> bool:
        return _interleaved(previous, boundary) or may_merge(region(previous), region(boundary))

    return compatible


def separate_overlap(left: object, right: object) -> tuple[object, object]:
    """Return two refused regions cut apart midway between their seed cores, in input order."""
    if left.end <= right.start or right.end <= left.start:
        return left, right
    first, second = sorted((left, right), key=lambda item: _core_span(item)[0])
    cut = (_core_span(first)[1] + _core_span(second)[0]) // 2
    first = replace(first, end=max(first.start + 1, min(first.end, cut)))
    second = replace(second, start=min(second.end - 1, max(second.start, cut)))
    return (first, second) if _core_span(left) <= _core_span(right) else (second, first)


def _interleaved(left: object, right: object) -> bool:
    """Return whether the seed cores of two regions interleave, which marks one EVE."""
    if not left.confirmed_seeds or not right.confirmed_seeds:
        return False
    left_start, left_end = _core_span(left)
    right_start, right_end = _core_span(right)
    return not (left_end <= right_start or right_end <= left_start)


def _core_span(boundary: object) -> tuple[int, int]:
    """Return the first core start and last core end, or the region span without cores."""
    cores = boundary.confirmed_seeds or [boundary]
    return min(core.start for core in cores), max(core.end for core in cores)


FUNCTIONS = {
    "MCP": "structure",
    "Portal": "structure",
    "Triplex": "structure",
    "Capsid": "structure",
    "Penton": "structure",
    "A32": "packaging",
    "ATPase": "packaging",
    "Terminase": "packaging",
    "Protease": "packaging",
    "PolB": "replication",
    "TopoII": "replication",
    "SFII": "replication",
    "HUH": "replication",
    "Rep": "replication",
    "RNAPL": "transcription",
    "RNAPS": "transcription",
    "RNAP": "transcription",
    "TFIIB": "transcription",
    "VLTF2": "transcription",
    "VLTF3": "transcription",
}
FULL_SET_FUNCTIONS = 3


def _groups(genes: Composition) -> set[str]:
    """Return the groups supported by a region's non-MCP hallmark genes."""
    return {GROUPS[panel] for pairs in genes.values() for panel, family in pairs if family != "MCP"}


def _is_full_set(genes: Composition) -> bool:
    """Return whether a region covers the functional categories that make a full EVE set."""
    groups = _groups(genes) or {GROUPS[panel] for pairs in genes.values() for panel, _ in pairs}
    available = {
        FUNCTIONS[family]
        for pairs in core_roles().values()
        for panel, family in pairs
        if GROUPS[panel] in groups and family in FUNCTIONS
    }
    present = {FUNCTIONS[family] for pairs in genes.values() for _, family in pairs if family in FUNCTIONS}
    return len(present) >= min(FULL_SET_FUNCTIONS, len(available))
