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

from types import SimpleNamespace

from virosync.pipeline.phase3.output_generator import (
    _resolve_eve_class,
    evaluate_v2_quality_gate,
)


def _result(**overrides) -> SimpleNamespace:
    base = dict(
        confidence_tier="HIGH",
        start=0,
        end=21077,  # ds27 length; well above the 2 kb floor
        hallmark_count=3,
        has_mcp=True,
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
        _result(has_mcp=False, hallmark_count=2, hallmark_genes=["COG0532", "VS000001"])
    )
    assert decision.kept
    assert decision.reason == "mixed_high_medium_pass"


def test_mixed_no_mcp_atpase_only_is_gated_out() -> None:
    # ATPase-only support (no MCP) is too weak -> gated out (non_atpase == 0).
    decision = evaluate_v2_quality_gate(
        _result(has_mcp=False, hallmark_count=2, hallmark_genes=["PLV_PC_054", "GVOGm0760"])
    )
    assert not decision.kept
    assert decision.reason == "mixed_high_medium_gate"


def test_mixed_no_mcp_single_hallmark_is_gated_out() -> None:
    decision = evaluate_v2_quality_gate(
        _result(has_mcp=False, hallmark_count=1, hallmark_genes=["COG0532"])
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
        _result(
            confidence_tier="LOW", has_mcp=False, hallmark_count=1, hallmark_genes=["COG0532"]
        )
    )
    assert not decision.kept
    assert decision.reason == "mixed_low_gate"


def test_mixed_low_atpase_only_is_gated_out() -> None:
    decision = evaluate_v2_quality_gate(
        _result(
            confidence_tier="LOW", has_mcp=False, hallmark_count=1, hallmark_genes=["GVOGm0760"]
        )
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
            has_mcp=True,
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
        _resolve_eve_class(
            SimpleNamespace(region_classification="", classification="MIXED", likely_family="")
        )
        == "MIXED"
    )


def test_concrete_family_wins_over_mixed_region_label() -> None:
    # A concrete fallback family must take precedence over a MIXED region label.
    assert (
        _resolve_eve_class(
            SimpleNamespace(region_classification="MIXED", classification="PPV", likely_family="")
        )
        == "PPV"
    )


def test_concrete_likely_family_wins_over_mixed_classification() -> None:
    # classification="MIXED" must NOT shadow a concrete likely_family.
    assert (
        _resolve_eve_class(
            SimpleNamespace(region_classification="", classification="MIXED", likely_family="PPV")
        )
        == "PPV"
    )


def test_concrete_likely_family_wins_over_unknown_classification() -> None:
    # classification="UNKNOWN" must NOT shadow a concrete likely_family.
    assert (
        _resolve_eve_class(
            SimpleNamespace(region_classification="", classification="UNKNOWN", likely_family="PPV")
        )
        == "PPV"
    )


def test_unknown_still_unsupported() -> None:
    # Guard: a genuinely unknown region is still rejected (no MIXED leakage).
    decision = evaluate_v2_quality_gate(
        _result(region_classification="UNKNOWN", classification="UNKNOWN", likely_family="UNKNOWN")
    )
    assert not decision.kept
    assert decision.reason == "unsupported_class"
