from __future__ import annotations

from types import SimpleNamespace

from virosync.pipeline.phase2.boundary_diamond import (
    GeneTaxonomy,
    SeedGeneMapping,
    classify_all_porfs,
    extract_sequences,
    get_flanking_taxonomy,
    pORF,
)
from virosync.pipeline.taxonomy_utils import resolve_org_id


def _taxonomy(porf_id: str, start: int, end: int) -> GeneTaxonomy:
    return GeneTaxonomy(
        porf_id=porf_id,
        scaffold="ctg",
        start=start,
        end=end,
    )


class _IterableIds:
    """Iterable that fails if code repeatedly probes the source container."""

    def __init__(self, values: list[str]) -> None:
        self.values = values

    def __iter__(self):
        return iter(self.values)

    def __contains__(self, value: object) -> bool:
        raise AssertionError(f"unexpected source membership test: {value}")


def test_boundary_query_ids_are_indexed_for_proteome_scans(tmp_path) -> None:
    proteome = tmp_path / "proteome.faa"
    extracted = tmp_path / "query.faa"
    proteome.write_text(
        ">keep # 1 # 9 # 1 # ID=1_1\nMKK\n"
        ">drop # 10 # 18 # 1 # ID=1_2\nMNN\n"
    )

    assert extract_sequences(
        proteome,
        _IterableIds(["keep"]),
        extracted,
    ) == 1
    assert extracted.read_text() == ">keep\nMKK\n"

    taxonomy = classify_all_porfs(
        all_porf_ids=_IterableIds(["keep"]),
        diamond_hits={},
        proteome_index={"ctg": [pORF("keep", "ctg", 0, 9)]},
        host_prefix="EUK__",
    )
    assert list(taxonomy) == ["keep"]


def test_flanks_follow_refined_boundary_after_seed_contraction() -> None:
    porfs = [
        pORF("old-upstream", "ctg", 0, 90),
        pORF("new-upstream-1", "ctg", 100, 190),
        pORF("new-upstream-2", "ctg", 200, 290),
        pORF("core", "ctg", 300, 390),
        pORF("new-downstream", "ctg", 400, 490),
        pORF("old-downstream", "ctg", 500, 590),
    ]
    taxonomy = {
        porf.id: _taxonomy(porf.id, porf.start, porf.end)
        for porf in porfs
    }
    original_seed_mapping = SeedGeneMapping(
        seed_id="seed",
        scaffold="ctg",
        seed_start=100,
        seed_end=500,
        upstream_porf_ids=["old-upstream"],
        downstream_porf_ids=["old-downstream"],
    )

    upstream, downstream = get_flanking_taxonomy(
        taxonomy_map=taxonomy,
        proteome_index={"ctg": porfs},
        refined_boundary=SimpleNamespace(scaffold="ctg", start=300, end=390),
        flank_genes=2,
        seed_mapping=original_seed_mapping,
    )

    assert [gene.porf_id for gene in upstream] == [
        "new-upstream-1",
        "new-upstream-2",
    ]
    assert [gene.porf_id for gene in downstream] == [
        "new-downstream",
        "old-downstream",
    ]


def test_flanks_fall_back_to_seed_mapping_and_exclude_refined_region() -> None:
    taxonomy = {
        "inside-upstream": _taxonomy("inside-upstream", 320, 350),
        "upstream": _taxonomy("upstream", 200, 290),
        "inside-downstream": _taxonomy("inside-downstream", 340, 370),
        "downstream": _taxonomy("downstream", 400, 490),
    }
    seed_mapping = SeedGeneMapping(
        seed_id="seed",
        scaffold="ctg",
        seed_start=100,
        seed_end=500,
        upstream_porf_ids=["inside-upstream", "upstream"],
        downstream_porf_ids=["inside-downstream", "downstream"],
    )

    upstream, downstream = get_flanking_taxonomy(
        taxonomy_map=taxonomy,
        proteome_index={"ctg": []},
        refined_boundary=SimpleNamespace(scaffold="ctg", start=300, end=390),
        seed_mapping=seed_mapping,
    )

    assert [gene.porf_id for gene in upstream] == ["upstream"]
    assert [gene.porf_id for gene in downstream] == ["downstream"]


def test_flanks_are_empty_without_final_boundary_or_seed_mapping() -> None:
    upstream, downstream = get_flanking_taxonomy(
        taxonomy_map={},
        proteome_index={},
        refined_boundary=SimpleNamespace(scaffold="ctg", start=300, end=390),
    )

    assert upstream == []
    assert downstream == []


def test_resolve_org_id_normalizes_legacy_vardna_namespace() -> None:
    lookup = {
        "PHAGE__GCA-000906975-1": (
            "PHAGE|Varidnaviria|Bamfordvirae|Preplasmiviricota"
        )
    }

    resolved = resolve_org_id(
        "PHAGE__VARDNA__GCA-000906975-1_23|protein",
        lookup,
    )

    assert resolved == "PHAGE__GCA-000906975-1"
