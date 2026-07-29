"""Canonical MCP (Major Capsid Protein) gene detection.

Single source of truth for deciding whether an HMM model ID or hallmark
gene name refers to a Major Capsid Protein.

Before Stage 1B, this decision was spread across ~10 call sites with a
lenient ``"mcp" in name.lower()`` substring match, which produced false
positives on any name containing the trigram (e.g. ``"dmcp"``,
``"ncmcp_pseudoprotein"``, ``"mcp_lookalike"``). Those false positives
flowed into the priority-marker promotion, the v2 quality gate, and the
scoring MCP bonus, silently inflating confidence.

The canonical corpus is:

- **Exact matches** (``MCP_HMM_EXACT``): the NCLDV MCP models
  (``OG1352``/``VS000086``, ``OG484``/``VS000309``, ``gamadvirusMCP``,
  ``GVOGm0003``) plus the bare symbolic token ``"mcp"`` emitted by some
  pipeline metadata.
- **Prefix matches** (``MCP_HMM_PREFIXES``): family-scoped MCP HMM names
  such as ``plv_mcp_1``, ``vp_mcp_3``, ``mirus_mcp_2``, ``mcp_mirus``,
  ``mcp_poli``. The suffix after the prefix must be alphanumeric /
  underscore only.

Any other name — even one containing the substring ``"mcp"`` — is
rejected.
"""

from __future__ import annotations

import re


NCLDV_MCP_MODELS: frozenset[str] = frozenset({
    "og1352",
    "og484",
    "vs000086",
    "vs000309",
    "gamadvirusmcp",
    "gvogm0003",
})


MCP_HMM_PREFIXES: tuple[str, ...] = (
    "gamadvirusmcp",
    "mcp_mirus",
    "mcp_poli",
    "mirus_mcp",
    "plv_mcp",
    "vp_mcp",
)


MCP_HMM_EXACT: frozenset[str] = NCLDV_MCP_MODELS | frozenset({"mcp"})


_WORD_CHARS_ONLY = re.compile(r"^[a-z0-9_]*$")


def is_mcp_gene(name: str | None) -> bool:
    """Return True when ``name`` identifies an MCP HMM model / hallmark gene.

    Matching rules (input is lower-cased first):

    1. Exact match against :data:`MCP_HMM_EXACT`.
    2. Prefix match against one of :data:`MCP_HMM_PREFIXES`, with the
       suffix required to be word-chars only (``[a-z0-9_]*``).

    Any other substring hit on ``"mcp"`` is rejected.
    """
    if not name:
        return False
    lower = name.lower().strip()
    if not lower:
        return False
    if lower in MCP_HMM_EXACT:
        return True
    for prefix in MCP_HMM_PREFIXES:
        if lower.startswith(prefix):
            suffix = lower[len(prefix):]
            if _WORD_CHARS_ONLY.match(suffix):
                return True
    return False


# Semantic alias: same decision, different caller intent.
is_mcp_model = is_mcp_gene
