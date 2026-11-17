"""Capsid class rules: a Mirus capsid names Mirusviricota, and capsids name unclassified calls."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from virosync.pipeline.phase3.class_rules import apply_capsid_class_rules, capsid_class
from virosync.pipeline.phase3.evidence_synthesizer import VerificationResult
from virosync.pipeline.phase3.mcp_detection import capsid_profile_class


def _marker(gene: str, model: str, *, start: int = 1_100, score: float = 100.0) -> SimpleNamespace:
    return SimpleNamespace(
        query_porf=f"{gene}|aa1-300",
        scaffold="scaffold",
        start=start,
        end=start + 900,
        hmm_target=model,
        hmm_score=score,
        validation_status="validated",
    )


def _call(eve_class: str, families: list[str] | None = None) -> VerificationResult:
    return VerificationResult(
        eve_id="eve",
        scaffold="scaffold",
        start=1_000,
        end=20_000,
        taxonomy_class=eve_class,
        marker_family_hits=families or [],
    )


@pytest.mark.parametrize(
    ("model", "expected"),
    [
        ("Mirus_MCP_3", "MIRUS"),
        ("MCP_mirus", "MIRUS"),
        ("EVE_CAPSID_NOVEL_0026", "MIRUS"),
        ("PLV_MCP_2", "PPV"),
        ("vp_penton_1", "PPV"),
        ("mcp_poli", "PPV"),
        ("GVOGm0003", "NCLDV"),
        ("OG1352", "NCLDV"),
        ("EVE_CAPSID_NOVEL_0033", "CRESS"),
        ("EVE_CAPSID_NOVEL_0018", "PHAGE"),
        ("VS000005", ""),
        ("GVOGm0013", ""),
    ],
)
def test_capsid_profile_class(model: str, expected: str) -> None:
    assert capsid_profile_class(model) == expected


def test_mirus_capsid_names_mirus_whatever_other_families_say() -> None:
    assert capsid_class(["Mirus_MCP_1", "GVOGm0013"], ["NCLDV", "MIRUS"]) == "MIRUS"


def test_other_capsid_needs_no_conflicting_family() -> None:
    assert capsid_class(["PLV_MCP_1"], ["PPV"]) == "PPV"
    assert capsid_class(["PLV_MCP_1"], ["NCLDV"]) == ""


def test_two_capsid_classes_name_nothing() -> None:
    assert capsid_class(["Mirus_MCP_1", "PLV_MCP_1"], []) == ""


@pytest.mark.parametrize(
    ("eve_class", "families", "models", "expected"),
    [
        ("NCLDV", ["NCLDV"], ["Mirus_MCP_1", "GVOGm0013"], "MIRUS"),
        ("PPV", [], ["Mirus_MCP_1", "PLV_MCP_1"], "PPV"),
        ("VIRAL_UNKNOWN", [], ["PLV_MCP_1"], "PPV"),
        ("UNKNOWN", ["NCLDV"], ["GVOGm0003"], "NCLDV"),
        ("VIRAL_UNKNOWN", ["NCLDV"], ["PLV_MCP_1"], "VIRAL_UNKNOWN"),
        ("VIRAL_UNKNOWN", [], ["VS000005"], "VIRAL_UNKNOWN"),
        ("NCLDV", [], ["PLV_MCP_1"], "NCLDV"),
    ],
)
def test_class_rules_on_calls(eve_class: str, families: list[str], models: list[str], expected: str) -> None:
    call = _call(eve_class, families)
    markers = [_marker(f"g{index}", model, start=1_100 + 1_000 * index) for index, model in enumerate(models)]

    changed = apply_capsid_class_rules([call], markers)

    assert call.taxonomy_class == expected
    assert changed == int(expected != eve_class)


def test_only_markers_overlapping_the_call_count() -> None:
    call = _call("VIRAL_UNKNOWN")
    outside = _marker("g1", "PLV_MCP_1", start=30_000)

    assert apply_capsid_class_rules([call], [outside]) == 0
    assert call.taxonomy_class == "VIRAL_UNKNOWN"


def test_best_model_per_gene_decides() -> None:
    call = _call("VIRAL_UNKNOWN")
    weaker = _marker("g1", "Mirus_MCP_1", score=50.0)
    stronger = _marker("g1", "PLV_MCP_1", score=90.0)

    apply_capsid_class_rules([call], [weaker, stronger])

    assert call.taxonomy_class == "PPV"
