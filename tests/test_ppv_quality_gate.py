"""PPV (Preplasmiviricota) gate and conservative subtype behavior.

After the GVClass PPV unification, a region labeled ``PPV`` must be accepted by the v2
quality gate under the same length/marker rules as the former PLV and VP classes.
VP and PLV are optional PPV subtypes, not separate result classes.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from virosync.pipeline.phase3.evidence_synthesizer import infer_ppv_subtype
from virosync.pipeline.phase3.gene_taxonomy import extract_prefix
from virosync.pipeline.phase3.output_generator import _evaluate_class_quality_gate


def _result(**overrides) -> SimpleNamespace:
    base = dict(
        confidence_tier="HIGH",
        start=0,
        end=3000,  # 3 kb, above the 2 kb PLV/VP/PPV floor
        hallmark_count=2,
        has_capsid_evidence=False,
        hallmark_genes=[],
        region_classification="PPV",
        classification="PPV",
        likely_family="PPV",
    )
    base.update(overrides)
    return SimpleNamespace(**base)


def test_ppv_mcp_high_medium_passes() -> None:
    decision = _evaluate_class_quality_gate(
        _result(hallmark_genes=["plv_mcp_caps_PgVV_Aquinto"], hallmark_count=1, has_capsid_evidence=True)
    )
    assert decision.kept
    assert decision.reason == "small_dna_high_medium_pass"


def test_ppv_atpase_only_high_medium_is_gated_out() -> None:
    # ATPase-only PPV region must be held out, exactly like PLV/VP.
    decision = _evaluate_class_quality_gate(_result(hallmark_genes=["PLV_PC_054", "VP_ATPase_1"], hallmark_count=2))
    assert not decision.kept
    assert decision.reason == "small_dna_high_medium_gate"


def test_ppv_non_atpase_hallmark_high_medium_passes() -> None:
    decision = _evaluate_class_quality_gate(_result(hallmark_genes=["PLV_PC_054", "VP_Penton_1"], hallmark_count=2))
    assert decision.kept


def test_ppv_low_mcp_is_promoted() -> None:
    decision = _evaluate_class_quality_gate(
        _result(
            confidence_tier="LOW",
            hallmark_genes=["plv_mcp_caps_PgVV_Aquinto"],
            hallmark_count=1,
            has_capsid_evidence=True,
        )
    )
    assert decision.kept
    assert decision.reason == "small_dna_low_promoted"


def test_ppv_subtype_requires_unambiguous_non_atpase_markers() -> None:
    assert infer_ppv_subtype(["VP_MCP_1", "VP_Penton_1"]) == "VP"
    assert infer_ppv_subtype(["PLV_MCP_1"]) == "PLV"
    assert infer_ppv_subtype(["VP_MCP_1", "PLV_MCP_1"]) == ""
    assert infer_ppv_subtype(["VP_ATPase_1", "PLV_PC_054"]) == ""


@pytest.mark.parametrize(
    "model",
    [
        "plv_mcp_caps_GKS2_Alpenseevirus_NCV_like",
        "plv_mcp_caps_Zephyrvirus4_NCV_like",
        "plv_mcp_caps_Zephyrvirus5_NCV_like",
        "plv_mcp_caps_VC124_VC134_NCV_like",
        "plv_mcp_caps_VC40_NCV_like",
        "plv_mcp_caps_Zoravirus_VC102_NCV_like",
    ],
)
@pytest.mark.parametrize("other_markers,expected", [([], ""), (["VP_MCP_1"], "VP")])
def test_ppv_subtype_does_not_treat_ncldv_capsid_names_as_plv(model, other_markers, expected):
    assert infer_ppv_subtype([*other_markers, model]) == expected


def test_extract_prefix_resolves_ppv_and_legacy_vp_plv() -> None:
    assert extract_prefix("PPV__IMGVR_UViG_3300066519|000001_8") == "PPV"
    # transitional: legacy VP/PLV targets still resolve to their fine class
    assert extract_prefix("VP__foo|1") == "VP"
    assert extract_prefix("PLV__bar|2") == "PLV"


@pytest.mark.parametrize(
    "markers,prefixes,identities,lineage,expected",
    [
        (["plv_mcp_caps_Mavpol2"], "PPV__", "90", "Preplasmiviricota", ""),
        (["plv_mcp_caps_Mavpol2"], "VP__", "90", "", "VP"),
        (["PLV_PC_054"], "PLV__", "90", "", "PLV"),
        (["PLV_PC_054"], "VP__", "24", "", ""),
        (["PLV_MCP_1"], "VP__", "90", "", ""),
        (["PLV_PC_054"], "PHAGE__", "90", "Preplasmiviricota;Virophages", "VP"),
        (["PLV_PC_054"], "PHAGE__", "90", "Preplasmiviricota", ""),
        (["PLV_PC_054"], "EUK-pEVE", "90", "Virophages", ""),
    ],
)
def test_ppv_subtype_uses_qualified_raw_reference_evidence(markers, prefixes, identities, lineage, expected):
    records = [{"top10_prefixes": prefixes, "top10_pidents": identities, "top10_targets": "PHAGE__ref|p"}]
    assert infer_ppv_subtype(markers, records, {"PHAGE__ref": lineage}) == expected


@pytest.mark.parametrize("prefix,expected", [("VP__", "VP"), ("PLV__", "PLV"), ("PPV__", "")])
def test_phase1_surface_and_phase3_keep_raw_subtype_evidence(tmp_path, monkeypatch, prefix, expected):
    from virosync.pipeline.phase1.hhg_seeding import Anchor
    from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
    from virosync.pipeline.phase3 import evidence_synthesizer, phase1_surface

    monkeypatch.setattr(evidence_synthesizer, "get_taxonomy_lookup", lambda: {})
    monkeypatch.setattr(phase1_surface, "get_taxonomy_lookup", lambda: {})
    anchor = Anchor(
        porf_id="s_1",
        scaffold="s",
        start=100,
        end=900,
        strand="+",
        hallmark_gene="PLV_PC_054",
        score=100,
        evalue=1e-8,
        validation_status="validated",
        top10_prefixes=prefix,
        top10_pidents="90",
        top10_targets=f"{prefix}ref|p",
    )
    seed = MergedSeed("s", 0, 3000, anchors=[anchor], confirmed_seeds=(ConfirmedSeed("seed", "s", 100, 900),))
    surface = phase1_surface.build_phase1_seed_surface(
        [seed], metadata_families={}, marker_annotation_index={}
    ).results[0]
    result = evidence_synthesizer.VerificationResult("eve", "s", 0, 3000, hallmark_genes=["PLV_PC_054"])
    synthesizer = evidence_synthesizer.EvidenceSynthesizer(
        config=evidence_synthesizer.EvidenceSynthesizerConfig(
            use_boltz=False, use_phylogenetic_validation=False, use_tmvec_database=False
        ),
        work_dir=tmp_path,
    )
    synthesizer._assign_taxonomy_class(result, [anchor], [])

    assert surface.ppv_subtype == result.ppv_subtype == expected
    assert surface.taxonomy_class == result.taxonomy_class == "PPV"


def test_ppv_subtype_reads_preserved_raw_prefixes_and_rejects_conflict():
    record = {
        "top10_prefixes": ["PPV", "PPV"],
        "top10_raw_prefixes": ["VP__", "PLV__"],
        "top10_pidents": [90, 24],
        "top10_targets": ["VP__a|p", "PLV__b|p"],
    }
    assert infer_ppv_subtype([], [record]) == "VP"
    record["top10_pidents"] = [90, 90]
    assert infer_ppv_subtype([], [record]) == ""


@pytest.mark.parametrize(
    "lineage,expected",
    [
        ("Preplasmiviricota|Aquintoviricetes", "PLV"),
        ("Preplasmiviricota|Polintoviricetes", ""),
        ("Preplasmiviricota|Maverick", ""),
        ("Preplasmiviricota|Lavidaviridae", ""),
        ("Preplasmiviricota|Maveriviricetes", "VP"),
        ("Preplasmiviricota|Virophaviricetes", "VP"),
        ("Preplasmiviricota|Phypoliviridae", "PLV"),
        ("Preplasmiviricota|Archintovirales", "PLV"),
        ("Preplasmiviricota|Eupolintoviridae", ""),
        ("Preplasmiviricota|Adintoviridae", ""),
        ("Preplasmiviricota|polinton-like viruses", "PLV"),
    ],
)
def test_ppv_subtype_resolves_only_informative_reference_ancestry(lineage, expected):
    record = {"top10_prefixes": "PHAGE__", "top10_targets": "PHAGE__ref|p", "top10_pidents": "90"}
    assert infer_ppv_subtype([], [record], {"PHAGE__ref": lineage}) == expected


def test_ppv_subtype_keeps_empty_prefix_position_aligned_with_identity():
    record = {"top10_prefixes": ",VP__", "top10_targets": ",VP__ref|p", "top10_pidents": "90,24"}
    assert infer_ppv_subtype([], [record]) == ""


@pytest.mark.parametrize("model", ["plv_mcp_caps_Mavpol2", "plv_mcp_caps_VC224_NA"])
def test_curated_cross_panel_mcp_is_ambiguous_until_reference_resolves_it(model):
    assert infer_ppv_subtype([model]) == ""
    record = {"top10_prefixes": "PLV__", "top10_pidents": "90", "top10_targets": "PLV__ref|p"}
    assert infer_ppv_subtype([model], [record]) == "PLV"


def test_ppv_composition_ambiguity_follows_curated_roles_not_model_spelling(monkeypatch):
    from virosync.pipeline import hallmark_evidence

    catalogue = dict(hallmark_evidence._catalogue())
    catalogue["phase1_marker", "PLV_MCP_TEST"] = catalogue["phase1_marker", "plv_mcp_caps_VC224_NA"]
    monkeypatch.setattr(hallmark_evidence, "_catalogue", lambda: catalogue)

    assert infer_ppv_subtype(["PLV_MCP_TEST"]) == ""
