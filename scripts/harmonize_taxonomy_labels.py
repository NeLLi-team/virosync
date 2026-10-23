"""Rewrite the eukaryotic lineages of ``taxonomy/labels.tsv`` from NCBI ranks.

The reference set is assembled from several sources, and each keeps the
vocabulary of its own taxonomy. *Aurantiochytrium limacinum* is

    EUK|Sagenista|Labyrinthulomycetes|NA|NA|Aurantiochytrium|...   (EukProt)
    EUK|Sar|Stramenopiles|NA|Thraustochytriidae|Aurantiochytrium|... (NCBI)

so two entries for the same organism share no token above genus. The host
signature of a run is built from whichever entries its marker genes hit, while
the genes it scores hit whichever entries are closest, so a gene and its own
host model can be written in different vocabularies. Any rule that compares
lineage tokens then scores the same gene differently according to which copy of
its relatives the database happens to hold.

This script gives every eukaryotic entry one lineage, taken from the NCBI
taxonomy dump: phylum, class, order, family, genus and species, with ``NA``
where NCBI has no name at that rank. Entries of other domains are copied
unchanged. Nothing is written in place: the harmonised table goes to a new
path, so the distributed resource stays as it is until a release replaces it.

Usage:
    python scripts/harmonize_taxonomy_labels.py LABELS_IN LABELS_OUT [--taxdump DIR]
"""

from __future__ import annotations

import argparse
import logging
from collections.abc import Iterator
from pathlib import Path

logger = logging.getLogger(__name__)

RANKS = ("phylum", "class", "order", "family", "genus", "species")
MISSING = "NA"
EUK_PREFIX = "EUK__"
LINEAGE_FIELDS = 7


class NcbiTaxonomy:
    """Scientific names and parent ranks from an NCBI taxonomy dump."""

    def __init__(
        self, names: dict[str, int], name_of: dict[int, str], parent: dict[int, int], rank: dict[int, str]
    ) -> None:
        self._names = names
        self._name_of = name_of
        self._parent = parent
        self._rank = rank

    @classmethod
    def load(cls, taxdump: Path) -> NcbiTaxonomy:
        """Read ``names.dmp`` and ``nodes.dmp`` from a taxonomy dump."""
        names: dict[str, int] = {}
        name_of: dict[int, str] = {}
        with (taxdump / "names.dmp").open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = [part.strip() for part in line.split("|")]
                if len(parts) > 3 and parts[3] == "scientific name":
                    taxid = int(parts[0])
                    names.setdefault(parts[1].lower(), taxid)
                    name_of[taxid] = parts[1]
        parent: dict[int, int] = {}
        rank: dict[int, str] = {}
        with (taxdump / "nodes.dmp").open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = [part.strip() for part in line.split("|")]
                taxid = int(parts[0])
                parent[taxid] = int(parts[1])
                rank[taxid] = parts[2]
        logger.info("taxonomy: %d names, %d nodes", len(names), len(parent))
        return cls(names, name_of, parent, rank)

    def lineage(self, name: str) -> list[str] | None:
        """Ranked lineage of one scientific name, or None when unknown."""
        taxid = self._names.get(name.strip().lower())
        if taxid is None:
            return None
        found: dict[str, str] = {}
        steps = 0
        while taxid and taxid != 1 and steps < 64:
            node_rank = self._rank.get(taxid, "")
            if node_rank in RANKS and node_rank not in found:
                found[node_rank] = self._name_of.get(taxid, "")
            taxid = self._parent.get(taxid, 1)
            steps += 1
        if not found:
            return None
        return [found.get(level) or MISSING for level in RANKS]


def harmonized_line(line: str, taxonomy: NcbiTaxonomy, counts: dict[str, int]) -> str:
    """Return one table line with a canonical eukaryotic lineage."""
    organism, _, lineage = line.rstrip("\n").partition("\t")
    if not lineage or not organism.startswith(EUK_PREFIX):
        counts["other domain"] += 1
        return line
    levels = [level.strip() for level in lineage.split("|")]
    if len(levels) != LINEAGE_FIELDS:
        counts["unexpected shape"] += 1
        return line
    for source, key in ((levels[-1], "by species"), (levels[-2], "by genus")):
        if source and source.upper() != MISSING:
            resolved = taxonomy.lineage(source)
            if resolved:
                counts[key] += 1
                return f"{organism}\t{levels[0]}|{'|'.join(resolved)}\n"
    counts["unresolved"] += 1
    return line


def harmonize(lines: Iterator[str], taxonomy: NcbiTaxonomy) -> tuple[list[str], dict]:
    """Harmonise every line of the table, reporting what each one took."""
    counts = dict.fromkeys(("by species", "by genus", "unresolved", "other domain", "unexpected shape"), 0)
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

    taxonomy = NcbiTaxonomy.load(args.taxdump)
    with args.labels_in.open(encoding="utf-8", errors="replace") as handle:
        harmonized, counts = harmonize(handle, taxonomy)
    tmp = args.labels_out.with_name(args.labels_out.name + ".tmp")
    tmp.write_text("".join(harmonized), encoding="utf-8")
    tmp.replace(args.labels_out)
    for key, value in counts.items():
        logger.info("%-16s %d", key, value)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
