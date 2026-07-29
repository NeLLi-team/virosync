#!/usr/bin/env python3
"""ViroSync PPV relabel: unify the Preplasmiviricota domain prefix VP__/PLV__ -> PPV__.

GVClass adopted the unified ``PPV`` (Preplasmiviricota) domain prefix; ViroSync must
follow so its TIER-1/TIER-2 resources and the v2 gate speak the same vocabulary. Unlike
GVClass (which had already folded VP__ -> PLV__), ViroSync carries BOTH live:

    VP__   virophage references   (lineage domain token "VP",  class Virophaviricetes)
    PLV__  PLV references         (lineage domain token "PLV", classes Polintoviricetes /
                                   Aquintoviricetes)

so a naive ``VP__/PLV__ -> PPV__`` collapse would create duplicate genome-id keys for the
1730 accessions that exist under BOTH prefixes. Those collisions are disambiguated with a
``__vpdup`` suffix on the VP-derived key (GVClass precedent), applied CONSISTENTLY to both
``labels.tsv`` (the key) and the proteome FASTA headers (so the faa<->labels bijection and
the .dmnd target-id <-> lineage lookup stay intact).

Rename rule (double-underscore IDs + leading lineage domain token ONLY):

    genome-id  VP__<acc>   -> PPV__<acc>__vpdup   if <acc> collides with a PLV__<acc>
    genome-id  VP__<acc>   -> PPV__<acc>          otherwise
    genome-id  PLV__<acc>  -> PPV__<acc>
    lineage token "VP"/"PLV" (leading domain)     -> "PPV"

PRESERVED untouched: CRESS__/GVMAG__/PHAGE__/MITO__/PLASTID__/NCLDV__/MIRUS__/BAC__/ARC__/
EUK__ ids and every taxon name (class Virophaviricetes / Polintoviricetes / Aquintoviricetes
survive in the lineage -> they carry the virophage-vs-PLV subcategory). NEVER touched:
single-underscore marker / lifestyle tokens (PLV_MCP_*, VP_MCP_*, plv_mcp_caps_*,
PLV_unclassified, VP_unclassified) -- they do not start with the double-underscore prefix.

Targets (resources/virosync/):
    taxonomy/labels.tsv                col0 id, col1 lineage         (in place, backed up)
    <proteome FASTA>                   ">" headers, streaming        (--proteome IN --proteome-out OUT)

Writes are atomic (*.tmp_ppv then os.replace) and back up to *.bak_ppv once. Idempotent:
re-running when PPV is already present is a no-op. Dry-run by default.

Usage:
    python scripts/relabel_ppv.py                          # dry-run on labels.tsv
    python scripts/relabel_ppv.py --apply                  # relabel labels.tsv (backed up)
    python scripts/relabel_ppv.py --proteome IN.faa --proteome-out OUT.faa --apply
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent / "resources" / "virosync"
LABELS = ROOT / "taxonomy" / "labels.tsv"
BAK = ".bak_ppv"
TMP = ".tmp_ppv"


def collision_accessions(labels_path: Path) -> set[str]:
    """Accessions that appear under BOTH VP__ and PLV__ in labels col0 (the dedup set)."""
    vp: set[str] = set()
    plv: set[str] = set()
    with labels_path.open() as fh:
        for line in fh:
            gid = line.split("\t", 1)[0]
            if gid.startswith("VP__"):
                vp.add(gid[4:])
            elif gid.startswith("PLV__"):
                plv.add(gid[5:])
    return vp & plv


def remap_id(gid: str, collisions: set[str]) -> str:
    """Genome-id (FASTA pre-'|' part / labels col0): VP__/PLV__ -> PPV__ with __vpdup dedup."""
    if gid.startswith("PLV__"):
        return "PPV__" + gid[5:]
    if gid.startswith("VP__"):
        acc = gid[4:]
        return f"PPV__{acc}__vpdup" if acc in collisions else f"PPV__{acc}"
    return gid


def remap_lineage_token(tok: str, collisions: set[str]) -> str:
    """A single lineage field: bare VP/PLV domain -> PPV; an embedded VP__/PLV__ genome-id
    token (the species field repeats the id with its prefix) -> PPV__ with the SAME __vpdup
    dedup as col0; everything else untouched (class names, single-underscore VP_unclassified
    lifestyle tags, other taxa)."""
    if tok in ("VP", "PLV"):
        return "PPV"
    return remap_id(tok, collisions)


def remap_lineage(lin: str, collisions: set[str]) -> str:
    """Remap EVERY '|'-delimited lineage field (not just the leading domain token): the
    lineage's trailing field repeats the genome-id with its VP__/PLV__ prefix and MUST be
    relabeled in lockstep with col0, or the migration leaves residual VP__/PLV__ tokens."""
    return "|".join(remap_lineage_token(tok, collisions) for tok in lin.split("|"))


def remap_header(hdr: str, collisions: set[str]) -> str:
    """FASTA header body (no leading '>'): id part before first '|', protein suffix untouched."""
    if "|" in hdr:
        idpart, prot = hdr.split("|", 1)
        return remap_id(idpart, collisions) + "|" + prot
    return remap_id(hdr, collisions)


def _atomic_write_text(path: Path, text: str) -> None:
    tmp = path.with_name(path.name + TMP)
    tmp.write_text(text)
    os.replace(tmp, path)


def _backup_once(path: Path) -> None:
    bak = path.with_name(path.name + BAK)
    if not bak.exists():
        bak.write_bytes(path.read_bytes())


def render_relabelled_labels(labels_path: Path, collisions: set[str]) -> tuple[str, dict]:
    """Compute the relabeled labels.tsv text + stats WITHOUT writing.

    Detection is separated from the write so the caller can refuse to overwrite when
    the remap would produce duplicate genome-id keys (which would make lineage lookups
    ambiguous) -- the dedup guard must gate the write, not run after it.
    """
    raw = labels_path.read_bytes()
    term = "\r\n" if b"\r\n" in raw else "\n"
    text = raw.decode()
    had_trailing = text.endswith(term)
    body = text[: -len(term)] if had_trailing else text
    out, changed = [], 0
    new_keys: dict[str, int] = {}
    for line in body.split(term):
        cols = line.split("\t")
        cols[0] = remap_id(cols[0], collisions)
        if len(cols) > 1:
            cols[1] = remap_lineage(cols[1], collisions)
        new = "\t".join(cols)
        if new != line:
            changed += 1
        new_keys[cols[0]] = new_keys.get(cols[0], 0) + 1
        out.append(new)
    rendered = term.join(out) + (term if had_trailing else "")
    stats = {
        "rows": len(out),
        "changed": changed,
        "duplicate_keys": sum(1 for c in new_keys.values() if c > 1),
        "ppv_keys": sum(1 for k in new_keys if k.startswith("PPV__")),
    }
    return rendered, stats


def relabel_proteome(in_path: Path, out_path: Path, collisions: set[str], apply: bool) -> dict:
    """Stream the (multi-GB) FASTA, rewriting only '>'-header genome-id prefixes."""
    headers, changed = 0, 0
    if not apply:
        # dry-run: count without writing
        with in_path.open() as fin:
            for line in fin:
                if line.startswith(">"):
                    headers += 1
                    body = line[1:].rstrip("\n")
                    if remap_header(body, collisions) != body:
                        changed += 1
        return {"headers": headers, "changed": changed, "written": False}
    tmp = out_path.with_name(out_path.name + TMP)
    with in_path.open() as fin, tmp.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                headers += 1
                ending = "\n" if line.endswith("\n") else ""
                body = line[1:].rstrip("\n")
                new_body = remap_header(body, collisions)
                if new_body != body:
                    changed += 1
                fout.write(">" + new_body + ending)
            else:
                fout.write(line)
    os.replace(tmp, out_path)
    return {"headers": headers, "changed": changed, "written": True}


def main() -> int:
    ap = argparse.ArgumentParser(description="ViroSync PPV relabel (VP__/PLV__ -> PPV__).")
    ap.add_argument("--apply", action="store_true", help="write changes (atomic; backs up to *.bak_ppv)")
    ap.add_argument("--labels", type=Path, default=LABELS, help="labels.tsv path (read)")
    ap.add_argument(
        "--labels-out",
        type=Path,
        default=None,
        help="write relabeled labels here instead of in place (keeps live labels<->.dmnd "
        "consistent until one atomic swap); when omitted, labels are rewritten in place + backed up",
    )
    ap.add_argument("--proteome", type=Path, default=None, help="input proteome FASTA to relabel")
    ap.add_argument("--proteome-out", type=Path, default=None, help="output relabeled FASTA path")
    args = ap.parse_args()

    mode = "APPLY" if args.apply else "DRY-RUN"
    print(f"=== ViroSync PPV relabel ({mode}) ===")
    if not args.labels.exists():
        print(f"ERROR: labels not found: {args.labels}")
        return 2

    collisions = collision_accessions(args.labels)
    print(f"VP__/PLV__ collision accessions (-> __vpdup on VP-derived): {len(collisions)}")

    rendered, lab = render_relabelled_labels(args.labels, collisions)
    print(
        f"labels.tsv: {lab['changed']}/{lab['rows']} rows changed; "
        f"PPV__ keys={lab['ppv_keys']}; duplicate keys AFTER dedup={lab['duplicate_keys']}"
    )
    # Guard the WRITE: never overwrite labels.tsv if the remap would collide keys.
    if lab["duplicate_keys"]:
        print("ERROR: duplicate genome-id keys remain after dedup -> aborting (lineage lookups would be ambiguous).")
        return 3
    if args.apply and lab["changed"]:
        if args.labels_out is not None:
            _atomic_write_text(args.labels_out, rendered)  # temp/staged output; live labels untouched
        else:
            _backup_once(args.labels)
            _atomic_write_text(args.labels, rendered)

    if args.proteome:
        if not args.proteome.exists():
            print(f"ERROR: proteome not found: {args.proteome}")
            return 2
        if args.apply and not args.proteome_out:
            print("ERROR: --proteome-out required with --apply --proteome")
            return 2
        prot = relabel_proteome(args.proteome, args.proteome_out or args.proteome, collisions, args.apply)
        print(f"proteome: {prot['changed']}/{prot['headers']} headers changed (written={prot['written']})")

    print("APPLIED." if args.apply else "DRY-RUN (re-run with --apply to write).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
