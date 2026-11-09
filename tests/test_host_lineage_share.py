"""Native extension consumes the rank-independent host token score."""

import pytest

from virosync.pipeline.host_signatures import HostSignatureModel, TaxonomyLabelLookup, set_taxonomy_lookup
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy, pORF
from virosync.pipeline.phase2.seed_extension import extend_confirmed_seed

LABELS = {"EUK__SELF": "EUK|Amoebozoa|Entamoeba", "EUK__FAR": "BAC|Bacteria", "NCLDV__V": "NCLDV|Virus"}


@pytest.fixture(autouse=True)
def lookup():
    set_taxonomy_lookup(TaxonomyLabelLookup(LABELS))
    yield
    set_taxonomy_lookup(None)


def _model():
    return HostSignatureModel(token_weights={"eukaryota": 10, "amoebozoa": 10, "entamoeba": 10})


@pytest.mark.parametrize("weight_mode", ["rank", "bitscore"])
@pytest.mark.parametrize(
    "target,viral_identity,expected_end,expected_stop",
    [
        ("EUK__SELF", None, 2000, 2100),
        ("EUK__FAR", None, 2800, None),
        ("EUK__SELF", 25.0, 2800, None),
        ("EUK__SELF", 24.9, 2000, 2100),
    ],
)
def test_extension_uses_real_lineage_share(target, viral_identity, expected_end, expected_stop, weight_mode):
    model = _model()
    model.weight_mode = weight_mode
    seed = MergedSeed(
        "s", 1000, 2000, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", 1000, 2000),)
    )
    genes = [
        pORF(id=f"h{index}", scaffold="s", start=start, end=start + 100, strand="+")
        for index, start in enumerate((2100, 2300, 2500, 2700))
    ]
    taxonomy = {}
    for index, gene in enumerate(genes):
        targets = ["NCLDV__V"] if index == 3 else [target]
        pidents = [90.0]
        if index == 1 and viral_identity is not None:
            targets.append("NCLDV__V")
            pidents.append(viral_identity)
        taxonomy[gene.id] = GeneTaxonomy(
            gene.id,
            "s",
            gene.start,
            gene.end,
            top10_targets=targets,
            top10_prefixes=["NCLDV__" if t.startswith("NCLDV") else "EUK__" for t in targets],
            top10_pidents=pidents,
            top10_bits=[100.0] * len(targets),
        )
    boundary = extend_confirmed_seed(seed, genes, taxonomy, (0, 10000), model, 0.35)
    assert (boundary.start, boundary.end) == (1000, expected_end)
    assert boundary.right_stop == expected_stop
    assert boundary.right_stop_reason == ("host_stop" if expected_stop is not None else "search_envelope")
    assert boundary.confirmed_seeds == seed.confirmed_seeds
