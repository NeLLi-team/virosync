"""Tests for the canonical MCP detection helper.

The pre-Stage-1B code used a lenient ``"mcp" in name.lower()`` substring
match at 10+ sites. This test set pins the new semantics: exact tokens
from a curated corpus + word-char-suffixed family prefixes are positives;
everything else is a negative, including names that contain the trigram
``mcp`` but aren't canonical MCP markers.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from virosync.pipeline.phase3.evidence_synthesizer import (
    load_marker_annotation_index,
    mcp_vote_classes,
    summarize_marker_hits,
)
from virosync.pipeline.phase3.mcp_detection import (
    MCP_HMM_EXACT,
    MCP_HMM_PREFIXES,
    NCLDV_MCP_MODELS,
    is_mcp_gene,
    is_mcp_model,
)

POSITIVES = [
    # NCLDV MCP models (exact)
    "og1352",
    "OG1352",
    "Og1352",
    "og484",
    "OG484",
    "vs000086",
    "VS000086",
    "vs000309",
    "VS000309",
    "gamadvirusmcp",
    "gamadvirusMCP",
    "GamadvirusMCP",
    "gvogm0003",
    "GVOGm0003",
    "GVOGM0003",
    # Curated structure-derived MCP model (exact, taxonomically neutral)
    "eve_capsid_novel_0029",
    "EVE_CAPSID_NOVEL_0029",
    # Bare symbolic MCP token
    "mcp",
    "MCP",
    # Family-scoped prefixes with numeric suffixes
    "plv_mcp_1",
    "PLV_MCP_1",
    "plv_mcp_10",
    "vp_mcp_3",
    "VP_MCP_3",
    "vp_mcp_7",
    "mirus_mcp_2",
    "Mirus_MCP_2",
    # Bare prefix form
    "plv_mcp",
    "vp_mcp",
    "mirus_mcp",
    # Alternate orderings (mcp_mirus vs mirus_mcp, mcp_poli)
    "mcp_mirus",
    "mcp_Mirus",
    "mcp_poli",
    "MCP_POLI",
]


NEGATIVES = [
    "",
    None,
    # Substring matches on "mcp" that are NOT canonical
    "dmcp",
    "mcp_lookalike",
    "NotMirus_MCP",
    "NCLDV_hypothetical",
    "ncmcp_pseudoprotein",
    "somcptest",
    "mcpfoo",  # no word boundary after mcp; and not a registered prefix
    # Real NCLDV markers that aren't MCP
    "polb",
    "polb_mirus",
    "vp_atpase_1",
    "vp_penton_3",
    "vp_pro_1",
    "EVE_CAPSID_NOVEL_0028",
    "EVE_CAPSID_NOVEL_0030",
    "OG1590",
    "OG2068",
    "Mirus_Terminase_merged",
    "gvogm0054",
    "GVOGM0023",
    "GVOGm0461",
    # Pure nonsense
    "hypothetical_protein",
    "unknown",
]


@pytest.mark.parametrize("name", POSITIVES)
def test_is_mcp_gene_positive(name: str) -> None:
    assert is_mcp_gene(name) is True, f"expected True for {name!r}"


@pytest.mark.parametrize("name", NEGATIVES)
def test_is_mcp_gene_negative(name: str | None) -> None:
    assert is_mcp_gene(name) is False, f"expected False for {name!r}"


def test_is_mcp_model_is_alias_of_is_mcp_gene() -> None:
    assert is_mcp_model is is_mcp_gene


def test_canonical_constants_contents() -> None:
    assert "og1352" in NCLDV_MCP_MODELS
    assert "og484" in NCLDV_MCP_MODELS
    assert "vs000086" in NCLDV_MCP_MODELS
    assert "vs000309" in NCLDV_MCP_MODELS
    assert "gamadvirusmcp" in NCLDV_MCP_MODELS
    assert "gvogm0003" in NCLDV_MCP_MODELS
    assert "eve_capsid_novel_0029" not in NCLDV_MCP_MODELS
    assert MCP_HMM_EXACT == NCLDV_MCP_MODELS | frozenset({"eve_capsid_novel_0029", "mcp"})
    # Family prefixes registered
    assert "plv_mcp" in MCP_HMM_PREFIXES
    assert "vp_mcp" in MCP_HMM_PREFIXES
    assert "mirus_mcp" in MCP_HMM_PREFIXES
    assert "mcp_mirus" in MCP_HMM_PREFIXES
    assert "mcp_poli" in MCP_HMM_PREFIXES
    assert "gamadvirusmcp" in MCP_HMM_PREFIXES


@pytest.mark.parametrize(
    "annotation_index",
    [
        None,
        {
            "eve_capsid_novel_0019": {
                "categories": {"capsid"},
                "family": "UNKNOWN",
            },
            "vp_penton_1": {
                "categories": {"capsid", "penton"},
                "family": "PPV",
            },
        },
    ],
)
def test_minor_capsid_and_penton_do_not_promote_mcp(
    annotation_index: dict[str, dict[str, object]] | None,
) -> None:
    summary = summarize_marker_hits(
        ["EVE_CAPSID_NOVEL_0019", "VP_Penton_1"],
        annotation_index,
    )

    assert summary["has_mcp"] is False
    assert "capsid" in summary["categories"]


def test_structure_derived_mcp_has_no_marker_family_vote() -> None:
    summary = summarize_marker_hits(
        ["EVE_CAPSID_NOVEL_0029"],
        {
            "eve_capsid_novel_0029": {
                "categories": {"capsid"},
                "family": "UNKNOWN",
            }
        },
    )

    assert summary["has_mcp"] is True
    assert summary["families"] == []
    assert summary["dominant_family"] == "UNKNOWN"


def test_non_mcp_structure_candidates_keep_annotated_functions(
    tmp_path: Path,
) -> None:
    annotations = tmp_path / "annotations.tsv"
    annotations.write_text(
        "model_name\tsource\tdescription\tmajority_annotation\n"
        "EVE_CAPSID_NOVEL_0001\tEVE_structure_20260811\t"
        "Putative adenain-like protease\t\n"
        "EVE_CAPSID_NOVEL_0006\tEVE_structure_20260811\t"
        "Putative adenain-like protease\t\n"
        "EVE_CAPSID_NOVEL_0015\tEVE_structure_20260811\t"
        "Putative adenain-like protease\t\n"
        "EVE_CAPSID_NOVEL_0018\tEVE_structure_20260811\t"
        "Putative tail tube protein\t\n"
        "EVE_CAPSID_NOVEL_0020\tEVE_structure_20260811\t"
        "Putative adenain-like protease\t\n"
        "EVE_CAPSID_NOVEL_0026\tEVE_structure_20260811\t"
        "Unresolved protease-like protein\t\n"
        "VP_Penton_1\tVP_PLV_update_Dec25\tUnannotated\t\n",
        encoding="utf-8",
    )
    annotation_index = load_marker_annotation_index(annotations)
    expected_categories = {
        "eve_capsid_novel_0001": {"protease"},
        "eve_capsid_novel_0006": {"protease"},
        "eve_capsid_novel_0015": {"protease"},
        "eve_capsid_novel_0018": {"tail_tube"},
        "eve_capsid_novel_0020": {"protease"},
        "eve_capsid_novel_0026": {"protease"},
    }

    assert {model: annotation_index[model]["categories"] for model in expected_categories} == expected_categories
    assert annotation_index["vp_penton_1"]["categories"] == {"penton"}

    summary = summarize_marker_hits(
        list(expected_categories),
        annotation_index,
    )
    assert summary["categories"] == ["protease", "tail_tube"]
    assert summary["complement_score"] == 0.25
    assert summary["has_mcp"] is False
    assert summary["families"] == []
    assert summary["dominant_family"] == "UNKNOWN"


def test_structure_derived_mcp_votes_only_with_sequence_taxonomy() -> None:
    hit = {
        "hallmark_gene": "EVE_CAPSID_NOVEL_0029",
        "porf_id": "protein_1",
        "top10_prefixes": "",
        "top10_pidents": "",
        "top10_targets": "",
        "validation_status": "validated_novel",
    }

    assert mcp_vote_classes([hit]) == set()

    hit["top10_prefixes"] = "PPV"
    hit["top10_pidents"] = "80"
    assert mcp_vote_classes([hit]) == {"PPV"}
