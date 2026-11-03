"""Tests for the reference mask built once per genome run."""

import logging

import pytest

from virosync.pipeline.reference_mask import ReferenceMask, build_reference_mask
from virosync.pipeline.search_backend import select_target_hits

LABELS = {
    "EUK__GCA-000001-1": "EUK|Oomycota|Peronosporomycetes|Peronosporales|Peronosporaceae|Phytophthora|Phytophthora ramorum",
    "EUK-pEVE__GCA-000001-1": "EUK-pEVE|Oomycota-pEVE|Peronosporomycetes-pEVE|Peronosporales-pEVE|Peronosporaceae-pEVE|Phytophthora-pEVE|Phytophthora ramorum-pEVE",
    "PPV__GCA-000001-1": "PPV|MCP in Oomycota|MCP in Peronosporomycetes|MCP in Peronosporales|MCP in Peronosporaceae|MCP in Phytophthora|MCP in Phytophthora ramorum",
    "EUK__GCA-000002-1": "EUK|Oomycota|Peronosporomycetes|Peronosporales|Peronosporaceae|Phytophthora|Phytophthora infestans",
    "EUK__EP00001": "EUK|NA|NA|NA|Collodictyonidae|Collodictyon|Collodictyon triciliatum",
    "NCLDV__V1": "NCLDV|Nucleocytoviricota|Megaviricetes|Imitervirales|Mimiviridae|Mimivirus|Mimivirus bradfordmassiliense",
}


def _mask(genome_id: str = "sample", **kwargs: object) -> ReferenceMask:
    settings: dict[str, object] = {"mask_own_assembly": True, "taxa": (), "labels": None, "self_keys": ()}
    settings.update(kwargs)
    return build_reference_mask(genome_id, **settings)


def test_own_assembly_key_uses_the_label_form_of_a_genbank_accession() -> None:
    mask = _mask("GCA_052324655.1_ASM5232465v1_genomic")

    assert "GCA-052324655-1" in mask.genome_keys
    assert mask.own_assembly_keys == ("GCA-052324655-1", "GCA_052324655.1")


def test_own_assembly_masks_the_raw_accession_that_organelle_rows_use() -> None:
    mask = _mask("GCF_000001215.4_Release_6_genomic")

    assert mask.genome_keys == {"GCF-000001215-4", "GCF_000001215.4"}


@pytest.mark.parametrize("genome_id", ["sample_01", "GCA_12345.1", "NC_000846.1", "GCA-052324655-1"])
def test_own_assembly_is_empty_without_a_full_accession(genome_id: str) -> None:
    mask = _mask(genome_id)

    assert mask.genome_keys == frozenset()
    assert mask.own_assembly_keys == ()


def test_own_assembly_can_be_switched_off() -> None:
    mask = _mask("GCA_052324655.1_genomic", mask_own_assembly=False)

    assert mask == ReferenceMask()


def test_own_assembly_takes_the_first_accession_of_the_genome_id() -> None:
    mask = _mask("GCA_000000001.1_vs_GCF_000000002.1")

    assert mask.own_assembly_keys == ("GCA-000000001-1", "GCA_000000001.1")


def test_species_mask_collects_every_label_of_the_matching_genome() -> None:
    mask = _mask(taxa=("Phytophthora ramorum",), labels=LABELS)

    assert mask.genome_keys == {"GCA-000001-1"}


def test_genus_mask_collects_every_genome_of_the_genus() -> None:
    mask = _mask(taxa=("Phytophthora",), labels=LABELS)

    assert mask.genome_keys == {"GCA-000001-1", "GCA-000002-1"}


def test_taxon_names_match_at_any_rank_ignoring_case_and_padding() -> None:
    mask = _mask(taxa=("  OOMYCOTA ", "collodictyon triciliatum"), labels=LABELS)

    assert mask.genome_keys == {"GCA-000001-1", "GCA-000002-1", "EP00001"}


def test_taxon_mask_does_not_match_part_of_a_rank_name() -> None:
    mask = _mask(taxa=("Phytophthora rama",), labels=LABELS)

    assert mask.genome_keys == frozenset()


def test_masks_combine_the_own_assembly_and_the_taxa() -> None:
    mask = _mask("GCA_052324655.1", taxa=("Collodictyon",), labels=LABELS)

    assert mask.genome_keys == {"EP00001", "GCA-052324655-1", "GCA_052324655.1"}
    assert mask.taxa == ("Collodictyon",)


def test_taxon_without_a_match_is_logged(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.WARNING, logger="virosync.pipeline.reference_mask"):
        mask = _mask(taxa=("Phytophthora ramorum", "Nowhere nothing"), labels=LABELS)

    assert mask.genome_keys == {"GCA-000001-1"}
    assert "Nowhere nothing" in caplog.text
    assert "Phytophthora ramorum" not in caplog.text


def test_taxa_need_the_taxonomy_labels() -> None:
    with pytest.raises(ValueError, match="taxonomy labels"):
        _mask(taxa=("Phytophthora",), labels=None)


def test_record_names_the_taxa_the_keys_and_their_number() -> None:
    mask = _mask("GCA_052324655.1", taxa=("Phytophthora ramorum",), labels=LABELS)

    assert mask.to_record() == {
        "n_genome_keys": 3,
        "taxa": ["Phytophthora ramorum"],
        "own_assembly_keys": ["GCA-052324655-1", "GCA_052324655.1"],
    }


def test_empty_mask_record_is_explicit() -> None:
    assert ReferenceMask().to_record() == {"n_genome_keys": 0, "taxa": [], "own_assembly_keys": []}


def test_a_name_suffixed_marker_id_is_masked_by_its_label_key() -> None:
    labels = {
        "EUK__EP00151": "EUK|Glomeromycota|Glomeromycetes|Glomerales|Glomeraceae|Rhizophagus|Rhizophagus irregularis"
    }
    mask = build_reference_mask("rhizophagus", mask_own_assembly=True, taxa=["Rhizophagus irregularis"], labels=labels)
    subjects = ["EUK__EP00151_Rhizophagus_irregularis|p1", "EUK__EP001512_Other_species|p2", "EUK__EP00999_Other|p3"]

    assert select_target_hits(subjects, 10, mask.genome_keys) == [1, 2]


def test_explicit_self_keys_replace_the_accession_in_the_genome_id() -> None:
    mask = _mask("GCA_052324655.1", self_keys=("EUK__JGI-12345", "IMG-678"))

    assert mask.own_assembly_keys == ("JGI-12345", "IMG-678")
    assert "GCA-052324655-1" not in mask.genome_keys


def test_an_accession_version_written_with_an_underscore_is_parsed() -> None:
    assert _mask("GCA_000002455_1").own_assembly_keys == ("GCA-000002455-1", "GCA_000002455.1")


def test_self_keys_follow_the_own_assembly_switch() -> None:
    assert _mask("x", mask_own_assembly=False, self_keys=("JGI-12345",)).genome_keys == frozenset()


@pytest.mark.parametrize(
    ("key", "other_key"),
    [
        ("sample", "sample_independent_assembly"),
        ("foo", "foo_bar"),
        ("JGI_123", "JGI_123_other"),
        ("IMG_456", "IMG_456_other"),
        ("GWH_ABC", "GWH_ABC_other"),
        ("GCA_000143045.1", "GCA_000143045.1_other"),
    ],
)
def test_explicit_keys_do_not_mask_a_different_underscore_key(key: str, other_key: str) -> None:
    mask = _mask(self_keys=(key,))
    subjects = [f"EUK__{key}|p1", f"EUK__{other_key}|p2"]

    assert select_target_hits(subjects, 10, mask.genome_keys) == [1]


@pytest.mark.parametrize("key", ["GCA-000002455-1", "GCA_000002455.1", "EUK__GCA-000002455-1", "GCF_000001215.4"])
def test_explicit_ncbi_keys_include_only_the_exact_accession_twin(key: str) -> None:
    mask = _mask(self_keys=(key,))
    expected = {"GCF-000001215-4", "GCF_000001215.4"} if "GCF" in key else {"GCA-000002455-1", "GCA_000002455.1"}
    assert mask.genome_keys == expected


def test_taxon_mask_remains_active_when_own_assembly_mask_is_off() -> None:
    mask = _mask("GCA_000000001.1", mask_own_assembly=False, taxa=("Collodictyon",), labels=LABELS)
    assert mask.genome_keys == {"EP00001"}
    assert mask.own_assembly_keys == ()
