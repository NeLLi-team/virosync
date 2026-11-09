"""Named ancestry survives harmonisation without guessing identities."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from harmonize_taxonomy_labels import NcbiTaxonomy, harmonize, main


def _taxonomy() -> NcbiTaxonomy:
    names = {"eukaryota": 2, "amoebozoa": 3, "discosea": 4, "dactylopodida": 5, "vannella": 6}
    return NcbiTaxonomy(
        names,
        {taxid: name for name, taxid in names.items()},
        {2: 1, 3: 2, 4: 3, 5: 4, 6: 4},
    )


@pytest.mark.parametrize("lineage", ["EUK|NA|Dactylopodida", "EUK|Dactylopodida|NA|Amoebozoa|NA|NA|NA|NA"])
def test_nonstandard_ancestors_survive_variable_shape_and_order(lineage: str) -> None:
    rewritten, counts = harmonize(iter([f"EUK__x\t{lineage}\n"]), _taxonomy())

    assert set(rewritten[0].strip().split("\t")[1].lower().split("|")) == {
        "euk",
        "eukaryota",
        "amoebozoa",
        "discosea",
        "dactylopodida",
    }
    assert counts["expanded"] == 1


def test_source_only_tokens_survive() -> None:
    rewritten, _ = harmonize(iter(["EUK__x\tEUK|CR clade|Dactylopodida|source strain\n"]), _taxonomy())

    assert "CR clade" in rewritten[0]
    assert "source strain" in rewritten[0]


def test_depth_restriction_does_not_gain_descendants() -> None:
    rewritten, _ = harmonize(iter(["EUK__x\tEUK|Discosea|NA\n"]), _taxonomy())

    assert rewritten == ["EUK__x\tEUK|eukaryota|amoebozoa|Discosea\n"]


@pytest.mark.parametrize("line", ["EUK__x\tEUK|NA|NA\n", "BAC__x\tBAC|something\n"])
def test_unresolved_root_and_other_domains_are_unchanged(line: str) -> None:
    rewritten, _ = harmonize(iter([line]), _taxonomy())

    assert rewritten == [line]


def test_conflicting_named_taxa_are_not_reassigned() -> None:
    line = "EUK__x\tEUK|Dactylopodida|Vannella\n"

    rewritten, counts = harmonize(iter([line]), _taxonomy())

    assert rewritten == ["EUK__x\tEUK|eukaryota|amoebozoa|discosea|Dactylopodida|Vannella\n"]
    assert counts["conflicting"] == 1


def test_ambiguous_scientific_name_is_not_guessed(tmp_path: Path) -> None:
    (tmp_path / "names.dmp").write_text(
        "2\t|\tDuplicate\t|\t\t|\tscientific name\t|\n3\t|\tDuplicate\t|\t\t|\tscientific name\t|\n", encoding="utf-8"
    )
    (tmp_path / "nodes.dmp").write_text("2\t|\t1\t|\tgenus\t|\n3\t|\t1\t|\tgenus\t|\n", encoding="utf-8")
    taxonomy = NcbiTaxonomy.load(tmp_path)
    line = "EUK__x\tEUK|Duplicate\n"

    rewritten, counts = harmonize(iter([line]), taxonomy)

    assert rewritten == [line]
    assert counts["ambiguous"] == 1


def test_unique_species_disambiguates_its_own_ancestor_name() -> None:
    taxonomy = NcbiTaxonomy(
        {"shared": None, "species one": 3},
        {2: "Shared", 3: "Species one"},
        {2: 1, 3: 2},
    )

    rewritten, counts = harmonize(iter(["EUK__x\tEUK|Shared|Species one\n"]), taxonomy)

    assert rewritten == ["EUK__x\tEUK|Shared|Species one\n"]
    assert counts["expanded"] == 1


@pytest.mark.parametrize("alias", ["same", "symlink", "hardlink"])
def test_main_refuses_overwriting_input_file(tmp_path: Path, alias: str) -> None:
    source = tmp_path / "labels.tsv"
    source.write_text("EUK__x\tEUK|Dactylopodida\n", encoding="utf-8")
    target = source if alias == "same" else tmp_path / "output.tsv"
    if alias == "symlink":
        target.symlink_to(source)
    elif alias == "hardlink":
        target.hardlink_to(source)

    with pytest.raises(ValueError, match="distinct"):
        main([str(source), str(target), "--taxdump", str(tmp_path / "absent")])
    assert source.read_text(encoding="utf-8") == "EUK__x\tEUK|Dactylopodida\n"


def test_domain_token_last_does_not_discard_first_named_taxon() -> None:
    rewritten, _ = harmonize(iter(["EUK__x\tCR clade|Dactylopodida|EUK\n"]), _taxonomy())

    assert rewritten[0].startswith("EUK__x\tEUK|")
    assert set(rewritten[0].strip().split("\t")[1].split("|")) == {
        "EUK",
        "eukaryota",
        "amoebozoa",
        "discosea",
        "Dactylopodida",
        "CR clade",
    }


def test_output_temporary_name_cannot_overwrite_input(tmp_path: Path) -> None:
    source = tmp_path / "output.tsv.tmp"
    original = "EUK__x\tEUK|Leaf\n"
    source.write_text(original, encoding="utf-8")
    (tmp_path / "names.dmp").write_text("2\t|\tLeaf\t|\t\t|\tscientific name\t|\n", encoding="utf-8")
    (tmp_path / "nodes.dmp").write_text("2\t|\t1\t|\tclade\t|\n", encoding="utf-8")
    assert main([str(source), str(tmp_path / "output.tsv"), "--taxdump", str(tmp_path)]) == 0
    assert source.read_text(encoding="utf-8") == original


def test_ambiguous_leaf_does_not_remove_known_parent_ancestry() -> None:
    taxonomy = NcbiTaxonomy(
        {"parent": 3, "shared": None},
        {2: "Named clade", 3: "Parent"},
        {2: 1, 3: 2},
    )

    rewritten, counts = harmonize(iter(["EUK__x\tEUK|Parent|Shared\n"]), taxonomy)

    assert rewritten == ["EUK__x\tEUK|Named clade|Parent|Shared\n"]
    assert counts["ambiguous"] == 1
    assert counts["expanded"] == 1
