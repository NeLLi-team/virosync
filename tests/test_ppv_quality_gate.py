"""PPV (Preplasmiviricota) regions gate identically to the legacy PLV/VP classes.

After the GVClass PPV unification, a region labeled ``PPV`` must be accepted by the v2
quality gate under the same length/marker rules as the former PLV and VP classes
(virophage + PLV are now class-rank subcategories of one PPV domain). VP/PLV remain
accepted transitionally so a pre-relabel bundle still gates correctly.
"""

from __future__ import annotations

from types import SimpleNamespace

from virosync.pipeline.phase3.gene_taxonomy import extract_prefix
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
    assert decision.reason == "plv_vp_high_medium_pass"


def test_ppv_atpase_only_high_medium_is_gated_out() -> None:
    # ATPase-only PPV region must be held out, exactly like PLV/VP.
    decision = evaluate_v2_quality_gate(
        _result(hallmark_genes=["PLV_PC_054", "VP_ATPase_1"], hallmark_count=2)
    )
    assert not decision.kept
    assert decision.reason == "plv_vp_high_medium_gate"


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
    assert decision.reason == "plv_vp_low_promoted"


def test_extract_prefix_resolves_ppv_and_legacy_vp_plv() -> None:
    assert extract_prefix("PPV__IMGVR_UViG_3300066519|000001_8") == "PPV"
    # transitional: legacy VP/PLV targets still resolve to their fine class
    assert extract_prefix("VP__foo|1") == "VP"
    assert extract_prefix("PLV__bar|2") == "PLV"
