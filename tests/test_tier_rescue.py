"""Tier rescue of NCLDV and Mirus regions without a core hallmark, and the integration bonus."""

from __future__ import annotations

from dataclasses import replace
from types import SimpleNamespace

import pytest

from virosync.ablation import AblationID
from virosync.pipeline.phase3 import acceptance_selection as selection_module
from virosync.pipeline.phase3.acceptance_selection import apply_tier_rescue, select_phase3_acceptance
from virosync.pipeline.phase3.evidence_synthesizer import (
    INTEGRATION_BONUS,
    VerificationResult,
    calculate_eve_confidence,
)
from virosync.pipeline.phase3.output_generator import QualityGateDecision


def _gene(index: int, prefixes: list[str], *, is_flanking: bool = False) -> dict[str, object]:
    return {"porf_id": f"g{index}", "is_flanking": is_flanking, "top10_prefixes": prefixes}


def _region(
    eve_id: str,
    *,
    start: int = 1_000,
    length: int = 6_000,
    eve_class: str = "NCLDV",
    viral_genes: int = 2,
    host_signature_fraction: float = 0.3,
) -> VerificationResult:
    """Return a LOW region without hallmarks that the normal gate rejects."""
    result = VerificationResult(
        eve_id=eve_id,
        scaffold="scaffold",
        start=start,
        end=start + length,
        confidence_tier="LOW",
        region_classification=eve_class,
        taxonomy_class=eve_class,
        host_signature_fraction=host_signature_fraction,
    )
    result.gene_taxonomy_records = [_gene(i, ["EUK__", f"{eve_class}__"]) for i in range(viral_genes)]
    return result


def _rescue(results: list[VerificationResult], ablation_id: AblationID = AblationID.A0):
    return apply_tier_rescue(select_phase3_acceptance(results, ablation_id))


def _integrase(mechanism: str, location: str = "interior") -> dict[str, object]:
    return {"protein_id": "p1", "mechanism": mechanism, "location": location}


def test_gate_rejected_region_meeting_the_rule_is_rescued_as_low() -> None:
    region = _region("eve")

    selection = _rescue([region])

    assert selection.normal_gate_decisions[0].kept is False
    assert selection.canonical_results == (region,)
    assert region.confidence_tier == "LOW"
    assert region.canonical_selection_outcome == "tier_rescue_low"
    assert selection.promoted_low_results == ()


@pytest.mark.parametrize(
    "mark",
    [
        {"tir_status": "detected"},
        {"repeat_evidence": SimpleNamespace(tsd=SimpleNamespace(status="assessed_match"))},
        {"integration_gene_hits": [_integrase("dde_integrase")]},
        {"integration_gene_hits": [_integrase("tyrosine_recombinase", "upstream")]},
        {"integration_gene_hits": [_integrase("serine_recombinase")]},
        {"gc_deviation": 0.10},
    ],
)
def test_rescued_region_with_an_integration_mark_becomes_medium(mark: dict[str, object]) -> None:
    region = replace(_region("eve"), **mark)

    selection = _rescue([region])

    assert selection.canonical_results == (region,)
    assert region.confidence_tier == "MEDIUM"
    assert region.canonical_selection_outcome == "tier_rescue_medium"


@pytest.mark.parametrize(
    "non_mark",
    [
        {"tir_status": "ambiguous"},
        {"repeat_evidence": SimpleNamespace(tsd=SimpleNamespace(status="assessed_no_match"))},
        {"integration_gene_hits": [_integrase("repair_recombinase")]},
        {"gc_deviation": 0.099},
    ],
)
def test_weak_or_non_integration_signals_keep_a_rescue_low(non_mark: dict[str, object]) -> None:
    region = replace(_region("eve"), **non_mark)

    _rescue([region])

    assert region.confidence_tier == "LOW"
    assert region.canonical_selection_outcome == "tier_rescue_low"


@pytest.mark.parametrize(
    "region",
    [
        _region("short", length=4_999),
        _region("one_gene", viral_genes=1),
        _region("host_like", host_signature_fraction=0.51),
        _region("ppv", eve_class="PPV"),
    ],
    ids=["short", "one_gene", "host_like", "ppv"],
)
def test_regions_outside_the_rule_stay_rejected(region: VerificationResult) -> None:
    selection = _rescue([region])

    assert selection.canonical_results == ()
    assert region.canonical_selection_outcome == "normal_gate_rejected"


def test_flanking_genes_do_not_count_as_viral_genes() -> None:
    region = _region("eve", viral_genes=1)
    region.gene_taxonomy_records.append(_gene(9, ["NCLDV__"], is_flanking=True))

    assert _rescue([region]).canonical_results == ()


def test_mirus_genes_count_for_an_ncldv_region_but_are_not_added_to_ncldv_genes() -> None:
    region = _region("eve", viral_genes=0)
    region.gene_taxonomy_records = [_gene(0, ["NCLDV__"]), _gene(1, ["MIRUS__"])]

    assert _rescue([region]).canonical_results == ()


def test_rescue_never_overlaps_a_kept_call() -> None:
    kept = VerificationResult(
        eve_id="kept",
        scaffold="scaffold",
        start=6_000,
        end=12_001,
        confidence_tier="HIGH",
        region_classification="NCLDV",
    )
    region = _region("eve", start=1_000, length=6_000)

    selection = _rescue([kept, region])

    assert selection.canonical_results == (kept,)
    assert region.canonical_selection_outcome == "normal_gate_rejected"


def test_overlapping_rescues_keep_more_viral_genes_then_the_longer_region() -> None:
    fewer_genes = _region("a_fewer_genes", start=1_000, length=9_000, viral_genes=2)
    more_genes = _region("b_more_genes", start=2_000, length=6_000, viral_genes=3)
    longer = _region("c_longer", start=20_000, length=8_000)
    shorter = _region("d_shorter", start=21_000, length=6_000)
    separate = _region("e_separate", start=40_000)

    selection = _rescue([fewer_genes, more_genes, longer, shorter, separate])

    assert selection.canonical_results == (more_genes, longer, separate)
    assert fewer_genes.canonical_selection_outcome == "normal_gate_rejected"


def test_kept_low_call_with_a_mark_is_promoted_to_medium(monkeypatch: pytest.MonkeyPatch) -> None:
    marked = replace(_region("marked"), tir_status="detected")
    unmarked = _region("unmarked", start=20_000)
    monkeypatch.setattr(
        selection_module,
        "evaluate_v2_quality_gate",
        lambda result: QualityGateDecision(True, "NCLDV", "ncldv_mirus_low_promoted", promoted_low=True),
    )

    selection = _rescue([marked, unmarked])

    assert selection.canonical_results == (marked, unmarked)
    assert marked.confidence_tier == "MEDIUM"
    assert marked.canonical_selection_outcome == "tier_promoted_medium"
    assert unmarked.confidence_tier == "LOW"
    assert unmarked.canonical_selection_outcome == "kept"
    assert selection.promoted_low_results == (unmarked,)


def test_a6_selection_is_not_changed() -> None:
    region = replace(_region("eve"), tir_status="detected")

    selection = _rescue([region], AblationID.A6)

    assert selection.canonical_results == (region,)
    assert region.confidence_tier == "LOW"
    assert region.canonical_selection_outcome == "kept"


@pytest.mark.parametrize(
    ("evidence", "has_bonus"),
    [
        ({"tir_status": "detected"}, True),
        ({"repeat_evidence": SimpleNamespace(tsd=SimpleNamespace(status="assessed_match"))}, True),
        ({"integration_gene_hits": [_integrase("dde_integrase")]}, True),
        ({"integration_gene_hits": [_integrase("tyrosine_recombinase", "downstream")]}, False),
        ({"integration_gene_hits": [_integrase("repair_recombinase")]}, False),
        ({"tir_status": "ambiguous"}, False),
    ],
)
def test_integration_bonus_counts_tir_tsd_and_interior_integrase_once(
    evidence: dict[str, object], has_bonus: bool
) -> None:
    bare = VerificationResult(eve_id="eve", scaffold="scaffold", start=0, end=6_000)
    marked = replace(bare, **evidence)
    both = replace(marked, tir_status="detected", integration_gene_hits=[_integrase("dde_integrase")])

    calculate_eve_confidence(bare, crf_confidence=0.0)
    calculate_eve_confidence(marked, crf_confidence=0.0)
    calculate_eve_confidence(both, crf_confidence=0.0)

    expected = INTEGRATION_BONUS if has_bonus else 0.0
    assert marked.score_components["bonus_total"] == pytest.approx(bare.score_components["bonus_total"] + expected)
    assert both.score_components["bonus_total"] == pytest.approx(
        bare.score_components["bonus_total"] + INTEGRATION_BONUS
    )
