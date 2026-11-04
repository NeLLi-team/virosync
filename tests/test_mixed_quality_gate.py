"""MIXED regions are a first-class accepted category in the v2 gate.

A MIXED region is one where the Phase 1 classifier saw markers from more than
one viral family and no family won the tie-break (``classify_region_by_markers``
-> ``"MIXED"``) -- the expected signature of NCLDV-adjacent capscan PLVs
(Aquintoviricetes "Near-" groups, "NCV-like" groups) that carry NCLDV-family
hallmark hits alongside their own capsid (capscan benchmark ds27 = Near-PgVV,
ds29 = VC40). These were dropped as ``unsupported_class``.

MIXED is now scored under the SAME rule as PLV/VP/PPV: an MCP is strong evidence
but NOT required -- >=2 hallmarks with >=1 non-ATPase also qualifies (HIGH/MED),
and a LOW MIXED region is promotable like a LOW PLV/VP region. A high-scoring
MIXED region is kept, not disqualified; only genuinely weak ones are gated out.
"""

from __future__ import annotations

import logging
from types import SimpleNamespace

import pytest

from virosync.pipeline.hallmark_evidence import HallmarkAssessment, HallmarkSummary
from virosync.pipeline.phase3.evidence_synthesizer import VerificationResult
from virosync.pipeline.phase3.output_generator import (
    _resolve_eve_class,
    evaluate_v2_quality_gate,
)
from virosync.pipeline.taxonomy_utils import EUK_PEVE_LABEL


def _result(**overrides) -> SimpleNamespace:
    base = dict(
        confidence_tier="HIGH",
        start=0,
        end=21077,  # ds27 length; well above the 2 kb floor
        hallmark_count=3,
        has_capsid_evidence=True,
        hallmark_genes=["plv_mcp_caps_SP_Aquinto", "COG0532", "VS000001"],
        region_classification="MIXED",
        likely_family="UNKNOWN",
    )
    base.update(overrides)
    return SimpleNamespace(**base)


# --- HIGH/MEDIUM acceptance under the normal viral rule ---


def test_mixed_with_mcp_high_passes() -> None:
    # ds27 / ds29 scenario: HIGH tier, MIXED class, capsid present.
    decision = evaluate_v2_quality_gate(_result())
    assert decision.kept
    assert decision.effective_class == "MIXED"
    assert decision.reason == "mixed_high_medium_pass"


def test_mixed_with_mcp_medium_passes() -> None:
    decision = evaluate_v2_quality_gate(_result(confidence_tier="MEDIUM"))
    assert decision.kept
    assert decision.reason == "mixed_high_medium_pass"


def test_mixed_no_mcp_but_strong_hallmarks_passes() -> None:
    # MCP is NOT required: >=2 hallmarks with >=1 non-ATPase qualifies.
    decision = evaluate_v2_quality_gate(
        _result(has_capsid_evidence=False, hallmark_count=2, hallmark_genes=["COG0532", "VS000001"])
    )
    assert decision.kept
    assert decision.reason == "mixed_high_medium_pass"


def test_mixed_no_mcp_atpase_only_is_gated_out() -> None:
    # ATPase-only support (no MCP) is too weak -> gated out (non_atpase == 0).
    decision = evaluate_v2_quality_gate(
        _result(has_capsid_evidence=False, hallmark_count=2, hallmark_genes=["PLV_PC_054", "GVOGm0760"])
    )
    assert not decision.kept
    assert decision.reason == "mixed_high_medium_gate"


def test_mixed_no_mcp_single_hallmark_is_gated_out() -> None:
    decision = evaluate_v2_quality_gate(
        _result(has_capsid_evidence=False, hallmark_count=1, hallmark_genes=["COG0532"])
    )
    assert not decision.kept
    assert decision.reason == "mixed_high_medium_gate"


def test_mixed_short_is_gated_out() -> None:
    # Below the 2 kb floor even with an MCP.
    decision = evaluate_v2_quality_gate(_result(end=2000))
    assert not decision.kept
    assert decision.reason == "mixed_high_medium_gate"


# --- LOW promotion under the normal viral low rule ---


def test_mixed_low_with_mcp_is_promoted() -> None:
    decision = evaluate_v2_quality_gate(_result(confidence_tier="LOW"))
    assert decision.kept
    assert decision.promoted_low
    assert decision.reason == "mixed_low_promoted"


def test_mixed_low_no_mcp_one_non_atpase_is_gated_out() -> None:
    # This vector used to be promoted at LOW while being rejected at MEDIUM and
    # HIGH, so raising a region's confidence removed it from the output. LOW now
    # applies the same rule as the tiers above it. See
    # tests/test_v2_gate_monotonicity.py for the general property.
    decision = evaluate_v2_quality_gate(
        _result(confidence_tier="LOW", has_capsid_evidence=False, hallmark_count=1, hallmark_genes=["COG0532"])
    )
    assert not decision.kept
    assert decision.reason == "mixed_low_gate"


def test_mixed_low_atpase_only_is_gated_out() -> None:
    decision = evaluate_v2_quality_gate(
        _result(confidence_tier="LOW", has_capsid_evidence=False, hallmark_count=1, hallmark_genes=["GVOGm0760"])
    )
    assert not decision.kept
    assert decision.reason == "mixed_low_gate"


def test_mixed_low_bridge_promotes_when_only_region_is_mixed() -> None:
    # LOW, likely_family=UNKNOWN but a MIXED label present: the bridge resolves
    # MIXED and promotes it under the normal low rule (does not silently drop).
    decision = evaluate_v2_quality_gate(
        _result(
            confidence_tier="LOW",
            region_classification="",
            classification="MIXED",
            likely_family="UNKNOWN",
            has_capsid_evidence=True,
        )
    )
    assert decision.kept
    assert decision.effective_class == "MIXED"
    assert decision.reason == "mixed_low_promoted"


def test_mixed_low_bridge_does_not_override_concrete_classification() -> None:
    # LOW, region=MIXED but a concrete classification=PPV: _resolve_eve_class
    # returns PPV, so the bridge must NOT rewrite the family to MIXED.
    decision = evaluate_v2_quality_gate(
        _result(
            confidence_tier="LOW",
            region_classification="MIXED",
            classification="PPV",
            likely_family="UNKNOWN",
        )
    )
    assert decision.effective_class != "MIXED"


# --- eve_class resolution / precedence ---


def test_mixed_resolved_from_classification_field_only() -> None:
    # ds27/ds29 output shape: region_classification empty, classification=MIXED.
    assert (
        _resolve_eve_class(SimpleNamespace(region_classification="", classification="MIXED", likely_family=""))
        == "MIXED"
    )


def test_concrete_family_wins_over_mixed_region_label() -> None:
    # A concrete fallback family must take precedence over a MIXED region label.
    assert (
        _resolve_eve_class(SimpleNamespace(region_classification="MIXED", classification="PPV", likely_family=""))
        == "PPV"
    )


def test_concrete_likely_family_wins_over_mixed_classification() -> None:
    # classification="MIXED" must NOT shadow a concrete likely_family.
    assert (
        _resolve_eve_class(SimpleNamespace(region_classification="", classification="MIXED", likely_family="PPV"))
        == "PPV"
    )


def test_concrete_likely_family_wins_over_unknown_classification() -> None:
    # classification="UNKNOWN" must NOT shadow a concrete likely_family.
    assert (
        _resolve_eve_class(SimpleNamespace(region_classification="", classification="UNKNOWN", likely_family="PPV"))
        == "PPV"
    )


def test_unknown_still_unsupported() -> None:
    # Guard: a genuinely unknown region is still rejected (no MIXED leakage).
    decision = evaluate_v2_quality_gate(
        _result(region_classification="UNKNOWN", classification="UNKNOWN", likely_family="UNKNOWN")
    )
    assert not decision.kept
    assert decision.reason == "unsupported_class"


def _repeated_peve_low(**overrides) -> VerificationResult:
    """Two separate loci with one model and only pEVE reference support."""
    values = dict(
        eve_id="eve",
        scaffold="sc",
        start=0,
        end=6000,
        confidence_tier="LOW",
        region_classification="MIXED",
        likely_family="MIXED",
        taxonomy_class="VIRAL_UNKNOWN",
        hallmark_count=2,
        hallmark_genes=["VS000254", "VS000254"],
        hallmark_assessment=HallmarkAssessment("sc", 0, 6000, HallmarkSummary(), ()),
        hallmark_taxonomy_records=[{"top10_prefixes": EUK_PEVE_LABEL, "top10_pidents": "60"}],
        gene_taxonomy_records=[{"top10_prefixes": [EUK_PEVE_LABEL], "top10_pidents": [60.0]}],
    )
    values.update(overrides)
    return VerificationResult(**values)


@pytest.mark.parametrize("family", ["NCLDV", "MIRUS", "PPV", "MIXED"])
def test_low_repeated_peve_only_markers_do_not_establish_acceptance(family: str) -> None:
    result = _repeated_peve_low(region_classification=family, likely_family=family)
    decision = evaluate_v2_quality_gate(result)
    assert not decision.kept
    assert not decision.promoted_low
    assert decision.reason == "low_repeated_peve_only_gate"


@pytest.mark.parametrize(
    ("overrides", "kept"),
    [
        pytest.param({"hallmark_genes": ["GVOGm9999", "GVOGM9999"]}, False, id="normalized-same-model"),
        pytest.param({"hallmark_genes": ["VS000254", "VS000151"]}, True, id="distinct-models"),
        pytest.param({"hallmark_genes": ["VS000254"]}, True, id="incomplete-model-list"),
        pytest.param({"hallmark_genes": ["", ""]}, True, id="missing-model-names"),
        pytest.param({"hallmark_assessment": None}, True, id="unknown-core-assessment"),
        pytest.param(
            {"hallmark_assessment": HallmarkAssessment("sc", 0, 6000, HallmarkSummary(core_count=1), ())},
            True,
            id="qualifying-core",
        ),
        pytest.param({"has_capsid_evidence": True}, True, id="capsid-evidence"),
        pytest.param(
            {
                "supplementary_viral_evidence": SimpleNamespace(
                    scaffold="sc", start=0, end=6000, eligible=True, group="NCLDV"
                )
            },
            True,
            id="eligible-module",
        ),
        pytest.param({"taxonomy_class": "NCLDV"}, True, id="lineage-taxonomy"),
        pytest.param({"confidence_tier": "MEDIUM"}, True, id="medium-unchanged"),
        pytest.param(
            {"region_classification": "CRESS", "likely_family": "CRESS", "marker_family_hits": ["CRESS"]},
            True,
            id="cress-identity-route",
        ),
        pytest.param({"hallmark_taxonomy_records": [], "gene_taxonomy_records": []}, True, id="no-viral-evidence"),
        pytest.param(
            {"gene_taxonomy_records": [{"top10_prefixes": ["NCLDV"], "top10_pidents": [25.0]}]},
            True,
            id="interior-direct-viral",
        ),
        pytest.param(
            {"hallmark_taxonomy_records": [{"top10_prefixes": "MIRUS", "top10_pidents": "25"}]},
            True,
            id="marker-direct-viral",
        ),
        pytest.param(
            {"gene_taxonomy_records": [{"top10_prefixes": ["NCLDV"], "top10_pidents": [24.9]}]},
            False,
            id="subthreshold-direct-hit",
        ),
        pytest.param(
            {"gene_taxonomy_records": [{"top10_prefixes": ["NCLDV"], "top10_pidents": [99.0], "is_flanking": True}]},
            False,
            id="flanking-direct-hit",
        ),
        pytest.param(
            {
                "gene_taxonomy_records": [
                    {
                        "top10_prefixes": [EUK_PEVE_LABEL],
                        "top10_pidents": [60.0],
                        "viral_hit_within_depth": True,
                    }
                ]
            },
            False,
            id="deeper-hit-is-not-positive-evidence",
        ),
    ],
)
def test_repeated_peve_low_gate_preserves_other_evidence(overrides: dict, kept: bool) -> None:
    result = _repeated_peve_low(**overrides)
    original = (result.start, result.end, result.final_confidence, result.confidence_tier)
    decision = evaluate_v2_quality_gate(result)
    assert decision.kept is kept
    assert decision.promoted_low is (kept and result.confidence_tier == "LOW")
    assert (decision.reason == "low_repeated_peve_only_gate") is (not kept)
    if kept:
        expected_reason = "mixed_low_promoted"
        if result.confidence_tier == "MEDIUM":
            expected_reason = "mixed_high_medium_pass"
        elif result.likely_family == "CRESS":
            expected_reason = "cress_identity_low_promoted"
        assert decision.reason == expected_reason
    assert (result.start, result.end, result.final_confidence, result.confidence_tier) == original


@pytest.mark.parametrize(
    ("family", "reason"),
    [("NCLDV", "ncldv_mirus_low_promoted"), ("MIRUS", "ncldv_mirus_low_promoted"), ("PPV", "small_dna_low_promoted")],
)
def test_direct_viral_support_preserves_the_low_gate_decision(family: str, reason: str) -> None:
    result = _repeated_peve_low(
        region_classification=family,
        likely_family=family,
        taxonomy_class=family,
        gene_taxonomy_records=[{"top10_prefixes": [family], "top10_pidents": [25.0]}],
    )
    decision = evaluate_v2_quality_gate(result)
    assert decision.kept and decision.promoted_low
    assert decision.reason == reason


def test_supplementary_module_preserves_its_low_promotion_reason() -> None:
    result = _repeated_peve_low(
        region_classification="NCLDV",
        likely_family="NCLDV",
        hallmark_count=1,
        hallmark_genes=["VS000254"],
        supplementary_viral_evidence=SimpleNamespace(scaffold="sc", start=0, end=6000, eligible=True, group="NCLDV"),
    )
    decision = evaluate_v2_quality_gate(result)
    assert decision.kept and decision.promoted_low
    assert decision.reason == "ncldv_mirus_low_supplementary_module"


def test_repeated_peve_low_rejection_has_a_candidate_trace(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.INFO):
        decision = evaluate_v2_quality_gate(_repeated_peve_low())
    assert not decision.kept
    assert "eve (sc:0-6000): low_repeated_peve_only_gate" in caplog.text
