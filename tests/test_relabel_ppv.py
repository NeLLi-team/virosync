"""Lock the PPV relabel + 1730-collision __vpdup dedup logic (scripts/relabel_ppv.py).

The dedup must be deterministic, idempotent, prefix-safe, and IDENTICAL between
labels.tsv genome-ids and FASTA headers (so the rebuilt .dmnd target-id <-> labels
lineage lookup stays bijective).
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "relabel_ppv", Path(__file__).resolve().parent.parent / "scripts" / "relabel_ppv.py"
)
rp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rp)


def test_remap_id_collision_dedup() -> None:
    collisions = {"GCA-1"}  # GCA-1 exists under BOTH VP__ and PLV__
    assert rp.remap_id("PLV__GCA-1", collisions) == "PPV__GCA-1"  # PLV wins the bare key
    assert rp.remap_id("VP__GCA-1", collisions) == "PPV__GCA-1__vpdup"  # VP collision -> __vpdup
    assert rp.remap_id("VP__GCA-2", collisions) == "PPV__GCA-2"  # VP, no collision
    assert rp.remap_id("PLV__GCA-3", collisions) == "PPV__GCA-3"


def test_remap_id_preserves_other_domains() -> None:
    c: set[str] = set()
    for gid in ("CRESS__x", "NCLDV__y", "MIRUS__z", "EUK__a", "GVMAG__b", "PHAGE__c", "BAC__d"):
        assert rp.remap_id(gid, c) == gid


def test_remap_never_touches_single_underscore_markers() -> None:
    c: set[str] = set()
    for tok in ("PLV_MCP_caps_1", "VP_MCP_x", "plv_mcp_caps_3", "PLV_unclassified", "VP_unclassified"):
        assert rp.remap_id(tok, c) == tok  # not a double-underscore domain prefix
        assert rp.remap_lineage(tok, c) == tok  # not a leading VP/PLV domain token


def test_remap_lineage_remaps_every_field() -> None:
    # leading bare domain token VP/PLV -> PPV
    assert rp.remap_lineage("VP|Virophaviricetes|Foo", set()) == "PPV|Virophaviricetes|Foo"
    assert rp.remap_lineage("PLV|Polintoviricetes|Bar", set()) == "PPV|Polintoviricetes|Bar"
    assert rp.remap_lineage("PLV|Aquintoviricetes|Baz", set()) == "PPV|Aquintoviricetes|Baz"
    # the trailing field repeats the genome-id with its VP__/PLV__ prefix -> must remap too
    assert (
        rp.remap_lineage("VP|Virophaviricetes|VP_unclassified|VP__Ace", set())
        == "PPV|Virophaviricetes|VP_unclassified|PPV__Ace"
    )
    # ...and the embedded id gets the SAME __vpdup dedup as col0 on a collision
    assert (
        rp.remap_lineage("VP|Virophaviricetes|VP__Dup", {"Dup"})
        == "PPV|Virophaviricetes|PPV__Dup__vpdup"
    )
    # single-underscore lifestyle tag preserved; a non-domain VP-ish token untouched
    assert rp.remap_lineage("EUK|Stramenopiles|VP_like", set()) == "EUK|Stramenopiles|VP_like"


def test_remap_header_matches_remap_id_consistently() -> None:
    c = {"GCA-1"}
    # FASTA header genome-id remap MUST equal the labels col0 remap for the same id
    assert rp.remap_header("VP__GCA-1|000001_8", c) == "PPV__GCA-1__vpdup|000001_8"
    assert rp.remap_header("PLV__GCA-1|000001_8", c) == "PPV__GCA-1|000001_8"
    assert rp.remap_header("EUK__Org|prot7", c) == "EUK__Org|prot7"
    # the id part (pre-'|') equals remap_id on the same id
    assert rp.remap_header("VP__GCA-1|p", c).split("|")[0] == rp.remap_id("VP__GCA-1", c)


def test_idempotent_no_double_relabel() -> None:
    c = {"GCA-1"}
    once = rp.remap_id("VP__GCA-1", c)
    assert once == "PPV__GCA-1__vpdup"
    assert rp.remap_id(once, c) == once  # PPV__ does not start with VP__/PLV__ -> no PPP corruption
    assert rp.remap_lineage(rp.remap_lineage("VP|Virophaviricetes", set()), set()) == "PPV|Virophaviricetes"
