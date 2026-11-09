"""Expand named eukaryotic ancestry without changing supplied identities.

All named ancestors are retained regardless of rank. Source-only tokens survive.
Ambiguous names remain unchanged while known parents retain their ancestry.
Conflicting known taxa gain only shared ancestors; resolution never adds descendants
below the deepest supplied named taxon. Other namespaces
are copied byte for byte. Output is a new label table, never an in-place edit.
"""

from __future__ import annotations

import argparse
import logging
from collections.abc import Iterator
from pathlib import Path
from tempfile import NamedTemporaryFile

logger = logging.getLogger(__name__)

MISSING = {"", "na", "n/a", "unknown", "unclassified", "none"}
EUK_PREFIX = "EUK__"


class NcbiTaxonomy:
    """Scientific names and parent relationships from an NCBI taxonomy dump."""

    def __init__(self, names: dict[str, int | None], name_of: dict[int, str], parent: dict[int, int]) -> None:
        self._names = names
        self._name_of = name_of
        self._parent = parent

    @classmethod
    def load(cls, taxdump: Path) -> NcbiTaxonomy:
        """Read ``names.dmp`` and ``nodes.dmp`` from a taxonomy dump."""
        names: dict[str, int | None] = {}
        name_of: dict[int, str] = {}
        with (taxdump / "names.dmp").open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = [part.strip() for part in line.split("|")]
                if len(parts) > 3 and parts[3] == "scientific name":
                    taxid = int(parts[0])
                    key = parts[1].casefold()
                    if key in names and names[key] != taxid:
                        names[key] = None
                    else:
                        names[key] = taxid
                    name_of[taxid] = parts[1]
        parent: dict[int, int] = {}
        with (taxdump / "nodes.dmp").open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = [part.strip() for part in line.split("|")]
                taxid = int(parts[0])
                parent[taxid] = int(parts[1])
        logger.info("taxonomy: %d names, %d nodes", len(names), len(parent))
        return cls(names, name_of, parent)

    def is_ambiguous(self, name: str) -> bool:
        """Return whether a scientific name identifies multiple taxonomy nodes."""
        key = name.strip().casefold()
        return key in self._names and self._names[key] is None

    def lineage(self, name: str) -> list[str] | None:
        """Return every named ancestor, or None for unknown or ambiguous names."""
        taxid = self._names.get(name.strip().casefold())
        if taxid is None:
            return None
        return self.lineage_for_taxid(taxid)

    def lineage_for_taxid(self, taxid: int) -> list[str] | None:
        """Return ancestors of an explicitly adjudicated taxid without name lookup."""
        ancestors: list[str] = []
        seen: set[int] = set()
        while taxid and taxid != 1:
            if taxid in seen:
                raise ValueError(f"taxonomy cycle at taxid {taxid}")
            seen.add(taxid)
            name = self._name_of.get(taxid, "")
            if name and name.casefold() != "cellular organisms":
                ancestors.append(name)
            taxid = self._parent.get(taxid, 1)
        return list(reversed(ancestors)) or None


def harmonized_line(line: str, taxonomy: NcbiTaxonomy, counts: dict[str, int]) -> str:
    """Add ancestors only when all resolved supplied taxa lie on one lineage."""
    organism, _, lineage = line.rstrip("\n").partition("\t")
    if not lineage or not organism.startswith(EUK_PREFIX):
        counts["other domain"] += 1
        return line
    tokens = [level.strip() for level in lineage.split("|")]
    if not any(token.casefold() in {"euk", "eukaryota"} for token in tokens):
        counts["unresolved"] += 1
        return line
    tokens = [token for token in tokens if token.casefold() not in MISSING | {"euk"}]
    resolved = [ancestry for token in tokens if (ancestry := taxonomy.lineage(token))]
    ancestors = max(resolved, key=len) if resolved else []
    ancestor_names = {name.casefold() for name in ancestors}
    if any(taxonomy.is_ambiguous(token) and token.casefold() not in ancestor_names for token in tokens):
        counts["ambiguous"] += 1
        if not resolved:
            return line
    if not resolved:
        counts["unresolved"] += 1
        return line
    if any(not {name.casefold() for name in lineage} <= ancestor_names for lineage in resolved):
        counts["conflicting"] += 1
        shared = set.intersection(*({name.casefold() for name in lineage} for lineage in resolved))
        ancestors = [name for name in ancestors if name.casefold() in shared]
        ancestor_names = shared
    # Preserve source spelling and source-only clades, not fixed rank slots.
    source_names = {token.casefold(): token for token in tokens}
    retained = [source_names.get(name.casefold(), name) for name in ancestors]
    retained.extend(token for token in tokens if token.casefold() not in ancestor_names)
    counts["expanded"] += 1
    return f"{organism}\tEUK|{'|'.join(retained)}\n"


def harmonize(lines: Iterator[str], taxonomy: NcbiTaxonomy) -> tuple[list[str], dict[str, int]]:
    """Expand lineages and report unresolved, ambiguous and conflicting rows."""
    counts = dict.fromkeys(("expanded", "unresolved", "ambiguous", "conflicting", "other domain"), 0)
    return [harmonized_line(line, taxonomy, counts) for line in lines], counts


def main(argv: list[str] | None = None) -> int:
    """Harmonise a labels table and report the outcome per entry."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels_in", type=Path, help="distributed labels.tsv")
    parser.add_argument("labels_out", type=Path, help="where to write the result")
    parser.add_argument(
        "--taxdump",
        type=Path,
        default=Path.home() / ".taxonkit",
        help="directory holding names.dmp and nodes.dmp",
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    if args.labels_in.resolve() == args.labels_out.resolve() or (
        args.labels_out.exists() and args.labels_in.samefile(args.labels_out)
    ):
        raise ValueError("input and output labels must be distinct files")

    taxonomy = NcbiTaxonomy.load(args.taxdump)
    with args.labels_in.open(encoding="utf-8", errors="replace") as handle:
        harmonized, counts = harmonize(handle, taxonomy)
    with NamedTemporaryFile(mode="w", dir=args.labels_out.parent, encoding="utf-8", delete=False) as handle:
        handle.write("".join(harmonized))
        temporary = Path(handle.name)
    temporary.replace(args.labels_out)
    for key, value in counts.items():
        logger.info("%-16s %d", key, value)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
