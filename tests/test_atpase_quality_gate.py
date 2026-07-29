"""The v2 acceptance gate must hold out packaging-ATPase-only PLV/VP regions.

The broad PLV/VP packaging ATPase (PLV_PC_054 and the VP_ATPase models) cross-hits
ubiquitous cellular P-loop NTPases, so a region whose only hallmark support is an
ATPase is unreliable and must not, on its own, produce an accepted call.
"""

from __future__ import annotations

from types import SimpleNamespace

from virosync.pipeline.phase3.output_generator import evaluate_v2_quality_gate


def _result(**overrides) -> SimpleNamespace:
    base = dict(
        confidence_tier="HIGH",
        start=0,
        end=3000,  # 3 kb, above the 2 kb PLV/VP floor
        hallmark_count=2,
        has_mcp=False,
        hallmark_genes=[],
        region_classification="PLV",
        classification="PLV",
        likely_family="PLV",
    )
    base.update(overrides)
    return SimpleNamespace(**base)


def test_atpase_only_high_medium_plv_is_gated_out() -> None:
    decision = evaluate_v2_quality_gate(
        _result(hallmark_genes=["PLV_PC_054", "VP_ATPase_1"], hallmark_count=2)
    )
    assert not decision.kept
    assert decision.reason == "small_dna_high_medium_gate"


def test_non_atpase_hallmark_high_medium_plv_passes() -> None:
    decision = evaluate_v2_quality_gate(
        _result(hallmark_genes=["PLV_PC_054", "VP_Penton_1"], hallmark_count=2)
    )
    assert decision.kept
    assert decision.reason == "small_dna_high_medium_pass"


def test_mcp_only_high_medium_plv_passes() -> None:
    # An MCP alone is sufficient even without a second hallmark.
    decision = evaluate_v2_quality_gate(
        _result(hallmark_genes=["plv_mcp_caps_PgVV_Aquinto"], hallmark_count=1, has_mcp=True)
    )
    assert decision.kept


def test_atpase_only_low_plv_is_gated_out() -> None:
    decision = evaluate_v2_quality_gate(
        _result(confidence_tier="LOW", hallmark_genes=["PLV_PC_054"], hallmark_count=1)
    )
    assert not decision.kept
    assert decision.reason == "small_dna_low_gate"


def test_mcp_low_plv_is_promoted() -> None:
    decision = evaluate_v2_quality_gate(
        _result(
            confidence_tier="LOW",
            hallmark_genes=["plv_mcp_caps_PgVV_Aquinto"],
            hallmark_count=1,
            has_mcp=True,
        )
    )
    assert decision.kept
    assert decision.reason == "small_dna_low_promoted"
