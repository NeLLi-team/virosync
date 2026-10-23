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
  ``GVOGm0003``), the curated structure-derived model
  ``EVE_CAPSID_NOVEL_0029``, and the bare symbolic token ``"mcp"`` emitted
  by some pipeline metadata. The structure-derived model carries no lineage
  assignment.
- **Prefix matches** (``MCP_HMM_PREFIXES``): family-scoped MCP HMM names
  such as ``plv_mcp_1``, ``vp_mcp_3``, ``mirus_mcp_2``, ``mcp_mirus``,
  ``mcp_poli``. The suffix after the prefix must be alphanumeric /
  underscore only.

Any other name — even one containing the substring ``"mcp"`` — is
rejected.

Capsid evidence is a wider decision. A hit to one of the capsid-type
structure-derived profiles (:data:`CAPSID_EVIDENCE_MODELS`) is capsid evidence
without being an MCP: scoring and acceptance treat it like an MCP, while
everything that names a major capsid protein still asks :func:`is_mcp_gene`.
"""

from __future__ import annotations

import re

NCLDV_MCP_MODELS: frozenset[str] = frozenset(
    {
        "og1352",
        "og484",
        "vs000086",
        "vs000309",
        "gamadvirusmcp",
        "gvogm0003",
    }
)


MCP_HMM_PREFIXES: tuple[str, ...] = (
    "gamadvirusmcp",
    "mcp_mirus",
    "mcp_poli",
    "mirus_mcp",
    "plv_mcp",
    "vp_mcp",
)


MCP_HMM_EXACT: frozenset[str] = NCLDV_MCP_MODELS | frozenset({"eve_capsid_novel_0029", "mcp"})


# The 27 capsid-type structure-derived profiles: those whose curated function in
# release-manifests/resources_v1_1_0/capsid_model_curation.tsv is a penton, capsid
# or coat protein. The adenain-like and unresolved protease profiles (0001, 0006,
# 0015, 0020, 0026) and the tail-tube profile (0018) stay viral markers without
# capsid evidence. A test checks this set against the curation table.
CAPSID_EVIDENCE_MODELS: frozenset[str] = frozenset(
    {
        "eve_capsid_novel_0002",
        "eve_capsid_novel_0003",
        "eve_capsid_novel_0004",
        "eve_capsid_novel_0005",
        "eve_capsid_novel_0007",
        "eve_capsid_novel_0008",
        "eve_capsid_novel_0009",
        "eve_capsid_novel_0010",
        "eve_capsid_novel_0011",
        "eve_capsid_novel_0012",
        "eve_capsid_novel_0013",
        "eve_capsid_novel_0014",
        "eve_capsid_novel_0016",
        "eve_capsid_novel_0017",
        "eve_capsid_novel_0019",
        "eve_capsid_novel_0021",
        "eve_capsid_novel_0022",
        "eve_capsid_novel_0023",
        "eve_capsid_novel_0024",
        "eve_capsid_novel_0025",
        "eve_capsid_novel_0027",
        "eve_capsid_novel_0028",
        "eve_capsid_novel_0029",
        "eve_capsid_novel_0030",
        "eve_capsid_novel_0031",
        "eve_capsid_novel_0032",
        "eve_capsid_novel_0033",
    }
)


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
            suffix = lower[len(prefix) :]
            if _WORD_CHARS_ONLY.match(suffix):
                return True
    return False


# Semantic alias: same decision, different caller intent.
is_mcp_model = is_mcp_gene


def is_capsid_evidence_gene(name: str | None) -> bool:
    """Return True when ``name`` is an MCP or a capsid-type structure-derived profile.

    The profile name is lower-cased and stripped, as in :func:`is_mcp_gene`.
    """
    if not name:
        return False
    return is_mcp_gene(name) or name.lower().strip() in CAPSID_EVIDENCE_MODELS
