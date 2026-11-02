"""Check that the scoring reference states the shipped catalogue and scoring rules."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import pytest

from virosync.pipeline import hallmark_evidence as hallmarks

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "docs" / "reference" / "scoring.md"
LINKS = json.loads(hallmarks.CATALOGUE_PATH.read_text())["links"]


def _families(panel: str, role: str) -> str:
    """Return one panel's families and exact model counts as the table writes them."""
    models = defaultdict(set)
    capsid = set()
    for link in LINKS:
        if link["panel"] != panel or link["role"] != role or link["excluded_from_scoring"]:
            continue
        if link["marker"] in hallmarks.EXCLUDED_FAMILIES:
            continue
        models[link["marker"]].add(link["model"])
        if link["capsid"]:
            capsid.add(link["marker"])
    cells = [
        f"**{name} ({len(models[name])})**" if name in capsid else f"{name} ({len(models[name])})"
        for name in sorted(models, key=str.lower)
    ]
    return ", ".join(cells) or "none"


def _scored_panels() -> list[str]:
    """Return panels with at least one core model that can score, in catalogue order."""
    scored = (link["panel"] for link in LINKS if link["role"] == "core_marker" and not link["excluded_from_scoring"])
    return list(dict.fromkeys(scored))


@pytest.mark.parametrize("panel", _scored_panels())
def test_scoring_page_lists_each_panel_as_the_catalogue_assigns_it(panel: str) -> None:
    row = f"| {panel} | {_families(panel, 'core_marker')} | {_families(panel, 'support_only')} |"

    assert row in PAGE.read_text(encoding="utf-8")


def test_scoring_page_states_the_catalogue_size() -> None:
    panels = {link["panel"] for link in LINKS}

    assert f"{len(LINKS)} links in {len(panels)}\npanels" in PAGE.read_text(encoding="utf-8")


@pytest.mark.parametrize(
    "statement",
    [
        f"Strong core | E ≤ {hallmarks.STRONG_EVALUE:g}".replace("1e-05", "1e-5"),
        f"E ≤ {hallmarks.WEAK_EVALUE:g}, not excluded",
    ],
)
def test_scoring_page_states_the_evalue_cutoffs(statement: str) -> None:
    assert statement in PAGE.read_text(encoding="utf-8")
