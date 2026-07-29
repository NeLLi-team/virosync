"""Tests for the canonical MCP detection helper.

The pre-Stage-1B code used a lenient ``"mcp" in name.lower()`` substring
match at 10+ sites. This test set pins the new semantics: exact tokens
from a curated corpus + word-char-suffixed family prefixes are positives;
everything else is a negative, including names that contain the trigram
``mcp`` but aren't canonical MCP markers.
"""

from __future__ import annotations

import pytest

from virosync.pipeline.phase3.mcp_detection import (
    MCP_HMM_EXACT,
    MCP_HMM_PREFIXES,
    NCLDV_MCP_MODELS,
    is_mcp_gene,
    is_mcp_model,
)


POSITIVES = [
    # NCLDV MCP models (exact)
    "og1352", "OG1352", "Og1352",
    "og484", "OG484",
    "vs000086", "VS000086",
    "vs000309", "VS000309",
    "gamadvirusmcp", "gamadvirusMCP", "GamadvirusMCP",
    "gvogm0003", "GVOGm0003", "GVOGM0003",
    # Bare symbolic MCP token
    "mcp", "MCP",
    # Family-scoped prefixes with numeric suffixes
    "plv_mcp_1", "PLV_MCP_1", "plv_mcp_10",
    "vp_mcp_3", "VP_MCP_3", "vp_mcp_7",
    "mirus_mcp_2", "Mirus_MCP_2",
    # Bare prefix form
    "plv_mcp", "vp_mcp", "mirus_mcp",
    # Alternate orderings (mcp_mirus vs mirus_mcp, mcp_poli)
    "mcp_mirus", "mcp_Mirus",
    "mcp_poli", "MCP_POLI",
]


NEGATIVES = [
    "",
    None,
    # Substring matches on "mcp" that are NOT canonical
    "dmcp",
    "mcp_lookalike",
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
    # MCP_HMM_EXACT adds the bare "mcp" token on top of the NCLDV set
    assert MCP_HMM_EXACT == NCLDV_MCP_MODELS | frozenset({"mcp"})
    # Family prefixes registered
    assert "plv_mcp" in MCP_HMM_PREFIXES
    assert "vp_mcp" in MCP_HMM_PREFIXES
    assert "mirus_mcp" in MCP_HMM_PREFIXES
    assert "mcp_mirus" in MCP_HMM_PREFIXES
    assert "mcp_poli" in MCP_HMM_PREFIXES
    assert "gamadvirusmcp" in MCP_HMM_PREFIXES
