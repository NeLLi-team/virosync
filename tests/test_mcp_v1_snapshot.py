"""v1 snapshot test for the Stage 1B MCP semantics change.

Before Stage 1B, MCP detection used a lenient ``"mcp" in name.lower()``
substring rule at roughly 10 call sites. Stage 1B replaces every site with
the canonical ``is_mcp_gene`` helper. This test pins the semantic diff by:

1. Running the legacy substring logic (``v1_legacy_substring``) against a
   representative corpus of real HMM names scraped from
   ``src/virosync/pipeline/phase1/viral_markers.py`` plus known NCLDV MCP
   model IDs.
2. Running the new ``is_mcp_gene`` helper against the same corpus.
3. Asserting that the set of names whose classification changes is
   exactly the expected **false-positive drop set** — names that only
   matched because of the trigram ``mcp`` and were never real MCP markers.
4. Asserting the **promotion set** (names the new helper calls True but
   v1 didn't) is empty on canonical inputs — the new logic never grows
   the positive set over legitimate markers.

Any future change to ``is_mcp_gene`` should update the expected drop set
explicitly, making the semantic shift auditable in code review.
"""

from __future__ import annotations

from virosync.pipeline.phase3.mcp_detection import is_mcp_gene


def v1_legacy_substring(name: str) -> bool:
    """Reproduce the pre-Stage-1B heuristic (worst-case permissive form).

    Most pre-Stage-1B call sites were variants of
    ``"mcp" in name.lower() or name.lower() == "gvogm0003"``. Some also
    substring-matched the NCLDV model IDs (``og1352``, ``og484``), but the
    result is dominated by the cheap ``"mcp" in ...`` branch. This helper
    reproduces the broader of the two so the snapshot captures every
    name the legacy code would have classified as MCP.
    """
    if not name:
        return False
    lower = name.lower()
    return (
        "mcp" in lower
        or lower == "gvogm0003"
        or lower == "og1352"
        or lower == "og484"
    )


# Canonical corpus drawn from real marker definitions in
# ``pipeline/phase1/viral_markers.py`` + NCLDV MCP model IDs from
# ``pipeline/phase3/mcp_detection.py``. Adding a marker here requires
# updating EXPECTED_FLIPS below if its classification changes.
CANONICAL_CORPUS: list[str] = [
    # NCLDV MCP models (positives under both v1 and v2)
    "og1352",
    "og484",
    "gamadvirusmcp",
    "gvogm0003",
    # NCLDV hallmark genes (short symbolic tokens from NCLDV_PROFILE)
    "mcp",  # generic NCLDV MCP token
    "a32",
    "d5",
    "vltf3",
    "mrnac",
    "polb",
    "rnapl",
    "rnaps",
    "rnr",
    "sfii",
    # Mriyavirus diagnostic + supporting markers
    "vltf2",
    "atpase_pkg",
    "huh_endo",
    "ruvc",
    "pddexk",
    "sf3_hel",
    "sf2_hel",
    "ssb",
    # Mirusvirus
    "mcp_mirus",
    "polb_mirus",
    "hel_mirus",
    # Polintovirus
    "mcp_poli",
    "ppolb",
    "pro_c1",
    "atpase",
    "int_tyr",
    # VP/PLV MCP prefixes (positives under both)
    "vp_mcp_1", "vp_mcp_2", "vp_mcp_3", "vp_mcp_4",
    "vp_mcp_5", "vp_mcp_6", "vp_mcp_7",
    "plv_mcp",
    # VP/PLV non-MCP (negatives under both)
    "vp_atpase_1", "vp_atpase_2", "vp_atpase_3", "vp_atpase_4",
    "vp_penton_1", "vp_penton_2", "vp_penton_3", "vp_penton_4",
    "vp_penton_5", "vp_penton_6", "vp_penton_7",
    "vp_pro_1", "vp_pro_2",
    "plv_pc_054",
    # Additional NCLDV OG/GVOGm IDs seen in real outputs (negatives)
    "OG1590", "OG2068", "OG516",
    "GVOGm0023", "GVOGm0054", "GVOGm0461",
    "Mirus_Terminase_merged",
    "Mirus_Portal",
    "Mirus_Triplex2",
    # Adversarial names that v1 would have wrongly classified as MCP
    # (these are the ones Stage 1B is explicitly meant to drop)
    "dmcp",
    "mcp_lookalike",
    "ncmcp_pseudoprotein",
    "somcptest",
    "mcpfoo",
    "hypothetical_mcp_containing_protein",
]


# Names whose classification changes under Stage 1B. All of these are
# v1=True -> v2=False because the new helper is strictly tighter on
# canonical inputs. If this list shrinks, something about is_mcp_gene
# became more permissive; if it grows, a new false-positive was
# discovered.
EXPECTED_FLIPS_V1_TO_V2 = {
    "dmcp",
    "mcp_lookalike",
    "ncmcp_pseudoprotein",
    "somcptest",
    "mcpfoo",
    "hypothetical_mcp_containing_protein",
}


def test_v1_to_v2_drops_match_expected() -> None:
    drops = {
        name
        for name in CANONICAL_CORPUS
        if name and v1_legacy_substring(name) and not is_mcp_gene(name)
    }
    assert drops == EXPECTED_FLIPS_V1_TO_V2, (
        "Stage 1B drop set diverged from expectation.\n"
        f"Expected: {sorted(EXPECTED_FLIPS_V1_TO_V2)}\n"
        f"Actual:   {sorted(drops)}"
    )


def test_v1_to_v2_does_not_grow_positive_set_on_canonical_inputs() -> None:
    """The new helper must not introduce NEW positives on the canonical
    corpus. (Adversarial negatives like "dmcp" are handled above.)"""
    promotions = {
        name
        for name in CANONICAL_CORPUS
        if name and not v1_legacy_substring(name) and is_mcp_gene(name)
    }
    # v2 may match some names (e.g. "mcp_mirus") that v1's weaker form
    # also matched (it contains "mcp"), so this set should be empty on
    # this corpus. If a future canonical marker is added that v1 missed
    # but v2 should catch, that name must be listed here.
    assert promotions == set(), (
        f"Unexpected new positives under is_mcp_gene: {sorted(promotions)}"
    )


def test_v2_rejects_every_adversarial_name() -> None:
    adversarial = EXPECTED_FLIPS_V1_TO_V2
    for name in adversarial:
        assert is_mcp_gene(name) is False, (
            f"Adversarial name {name!r} must be rejected by is_mcp_gene"
        )


def test_v2_accepts_every_canonical_mcp_name() -> None:
    canonical_mcp = {
        "og1352", "og484", "gamadvirusmcp", "gvogm0003", "mcp",
        "mcp_mirus", "mcp_poli", "plv_mcp",
        "vp_mcp_1", "vp_mcp_2", "vp_mcp_3", "vp_mcp_4",
        "vp_mcp_5", "vp_mcp_6", "vp_mcp_7",
    }
    for name in canonical_mcp:
        assert is_mcp_gene(name) is True, (
            f"Canonical MCP name {name!r} must be accepted by is_mcp_gene"
        )
