"""Tests for the taxonomy-label harmoniser."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from harmonize_taxonomy_labels import NcbiTaxonomy, harmonize

# A minimal tree: two thraustochytrid genera under one family, written here the
# way NCBI writes them, plus one green alga.
NAMES = {
    "bigyra": 11,
    "thraustochytrida": 12,
    "thraustochytriidae": 13,
    "aurantiochytrium": 14,
    "aurantiochytrium limacinum": 15,
    "hondaea": 16,
    "hondaea fermentalgiana": 17,
    "chlorophyta": 18,
    "chlamydomonas reinhardtii": 19,
}
NAME_OF = {taxid: name for name, taxid in NAMES.items()}
# 1 is the root, as in a real dump, so no named node may use it.
PARENT = {11: 1, 12: 11, 13: 12, 14: 13, 15: 14, 16: 13, 17: 16, 18: 1, 19: 18}
RANK = {
    11: "phylum",
    12: "order",
    13: "family",
    14: "genus",
    15: "species",
    16: "genus",
    17: "species",
    18: "phylum",
    19: "species",
}


def _taxonomy() -> NcbiTaxonomy:
    return NcbiTaxonomy(NAMES, NAME_OF, PARENT, RANK)


def test_two_sources_of_one_clade_end_up_with_the_same_upper_lineage() -> None:
    lines = [
        "EUK__EP00512\tEUK|Sagenista|Labyrinthulomycetes|NA|NA|Aurantiochytrium|Aurantiochytrium limacinum\n",
        "EUK__GCA-1\tEUK|Sar|Stramenopiles|NA|Thraustochytriidae|Hondaea|Hondaea fermentalgiana\n",
    ]

    harmonized, counts = harmonize(iter(lines), _taxonomy())

    upper = [line.split("\t")[1].split("|")[1:5] for line in harmonized]
    assert upper[0] == upper[1]
    assert upper[0] == ["bigyra", "NA", "thraustochytrida", "thraustochytriidae"]
    assert counts["by species"] == 2


def test_the_species_and_domain_of_each_entry_survive() -> None:
    line = "EUK__X\tEUK|Sar|Stramenopiles|NA|NA|Chlamydomonas|Chlamydomonas reinhardtii\n"

    harmonized, _ = harmonize(iter([line]), _taxonomy())

    fields = harmonized[0].split("\t")[1].strip().split("|")
    assert fields[0] == "EUK"
    assert fields[-1] == "chlamydomonas reinhardtii"
    assert len(fields) == 7


def test_entries_of_other_domains_are_copied_unchanged() -> None:
    line = "BAC__GCA-2\tBAC|Actinomycetota|Actinomycetes|Actinomycetales|X|Y|Y z\n"

    harmonized, counts = harmonize(iter([line]), _taxonomy())

    assert harmonized == [line]
    assert counts["other domain"] == 1


def test_an_unknown_organism_keeps_its_original_lineage() -> None:
    line = "EUK__Z\tEUK|Sar|Stramenopiles|NA|NA|Nothing|Nothing known\n"

    harmonized, counts = harmonize(iter([line]), _taxonomy())

    assert harmonized == [line]
    assert counts["unresolved"] == 1


def test_a_genus_resolves_when_the_species_name_is_unknown() -> None:
    line = "EUK__G\tEUK|Sar|Stramenopiles|NA|NA|Hondaea|Hondaea sp. strain 1\n"

    harmonized, counts = harmonize(iter([line]), _taxonomy())

    assert "thraustochytriidae" in harmonized[0]
    assert counts["by genus"] == 1
