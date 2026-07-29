"""PPV (Preplasmiviricota) gate and conservative subtype behavior.

After the GVClass PPV unification, a region labeled ``PPV`` must be accepted by the v2
quality gate under the same length/marker rules as the former PLV and VP classes.
VP and PLV are optional PPV subtypes, not separate result classes.
"""

from __future__ import annotations

from types import SimpleNamespace

from virosync.pipeline.phase3.gene_taxonomy import extract_prefix
from virosync.pipeline.phase3.evidence_synthesizer import infer_ppv_subtype
from virosync.pipeline.phase3.output_generator import evaluate_v2_quality_gate


def _result(**overrides) -> SimpleNamespace:
    base = dict(
        confidence_tier="HIGH",
        start=0,
        end=3000,  # 3 kb, above the 2 kb PLV/VP/PPV floor
        hallmark_count=2,
        has_mcp=False,
        hallmark_genes=[],
        region_classification="PPV",
        classification="PPV",
        likely_family="PPV",
    )
    base.update(overrides)
    return SimpleNamespace(**base)


def test_ppv_mcp_high_medium_passes() -> None:
    decision = evaluate_v2_quality_gate(
        _result(hallmark_genes=["plv_mcp_caps_PgVV_Aquinto"], hallmark_count=1, has_mcp=True)
    )
    assert decision.kept
    assert decision.reason == "small_dna_high_medium_pass"


def test_ppv_atpase_only_high_medium_is_gated_out() -> None:
    # ATPase-only PPV region must be held out, exactly like PLV/VP.
    decision = evaluate_v2_quality_gate(
        _result(hallmark_genes=["PLV_PC_054", "VP_ATPase_1"], hallmark_count=2)
    )
    assert not decision.kept
    assert decision.reason == "small_dna_high_medium_gate"


def test_ppv_non_atpase_hallmark_high_medium_passes() -> None:
    decision = evaluate_v2_quality_gate(
        _result(hallmark_genes=["PLV_PC_054", "VP_Penton_1"], hallmark_count=2)
    )
    assert decision.kept


def test_ppv_low_mcp_is_promoted() -> None:
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


def test_ppv_subtype_requires_unambiguous_non_atpase_markers() -> None:
    assert infer_ppv_subtype(["VP_MCP_1", "VP_Penton_1"]) == "VP"
    assert infer_ppv_subtype(["PLV_MCP_1"]) == "PLV"
    assert infer_ppv_subtype(["VP_MCP_1", "PLV_MCP_1"]) == ""
    assert infer_ppv_subtype(["VP_ATPase_1", "PLV_PC_054"]) == ""


def test_extract_prefix_resolves_ppv_and_legacy_vp_plv() -> None:
    assert extract_prefix("PPV__IMGVR_UViG_3300066519|000001_8") == "PPV"
    # transitional: legacy VP/PLV targets still resolve to their fine class
    assert extract_prefix("VP__foo|1") == "VP"
    assert extract_prefix("PLV__bar|2") == "PLV"
