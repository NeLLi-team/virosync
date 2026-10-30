"""Tests for the host-lineage share that replaces the aggregate Jaccard.

The old test compared a gene's ten hits with the aggregate of the whole host
model by weighted Jaccard. That asks whether a gene's hit distribution has the
same shape as the pooled distribution of every host marker gene, so in a host
whose clade is densely sequenced no single gene matches and real host genes are
rejected, including genes whose best hit is the host's own species.

The share asks the question the trimming step needs: how much of this gene's
hit weight sits inside the host's own lineage.
"""

from __future__ import annotations

import pytest

from virosync.pipeline.host_signatures import (
    HostSignatureModel,
    TaxonomyLabelLookup,
    dominant_lineage_tokens,
    host_lineage_share,
    score_host_signature,
    set_taxonomy_lookup,
)
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy, pORF
from virosync.pipeline.phase2.seed_extension import extend_confirmed_seed
from virosync.pipeline.taxonomy_utils import aggregate_taxonomy_substrings

HOST = "EUK|NA|Bigyra|Thraustochytrida|Thraustochytriidae|Aurantiochytrium|Aurantiochytrium limacinum"
RELATIVE = "EUK|NA|Bigyra|Thraustochytrida|Thraustochytriidae|Hondaea|Hondaea fermentalgiana"
DISTANT = "EUK|Chlorophyta|Chlorophyceae|Chlamydomonadales|Chlamydomonadaceae|Chlamydomonas|C. reinhardtii"
VIRUS = "MIRUS|Mirusviricota|NA|NA|NA|NA|Mirusvirus sp."

LABELS = {
    "EUK__SELF": HOST,
    "EUK__REL": RELATIVE,
    "EUK__FAR": DISTANT,
    "MIRUS__V": VIRUS,
}


def _model() -> HostSignatureModel:
    """A model whose heaviest tokens are the host's own lineage."""
    return HostSignatureModel(
        token_weights={
            "bigyra": 90.0,
            "thraustochytrida": 80.0,
            "thraustochytriidae": 70.0,
            "aurantiochytrium": 60.0,
            "aurantiochytrium limacinum": 50.0,
            "chlorophyta": 5.0,
        },
        host_prefixes=["EUK__"],
    )


def _record(targets: list[str]) -> dict:
    return {
        "top10_targets": targets,
        "top10_bitscores": [100.0] * len(targets),
        "top10_pidents": [90.0] * len(targets),
        "top10_evalues": [1e-40] * len(targets),
    }


def setup_module() -> None:
    """Register the small taxonomy these tests score against."""
    set_taxonomy_lookup(TaxonomyLabelLookup(LABELS))


def teardown_module() -> None:
    """Leave no lookup registered for other test modules."""
    set_taxonomy_lookup(None)


def test_dominant_tokens_take_the_heaviest_token_at_each_depth() -> None:
    tokens = dominant_lineage_tokens(_model(), LABELS.values())

    assert "aurantiochytrium" in tokens
    assert "chlorophyta" not in tokens  # lighter than bigyra at its depth
    assert "euk" not in tokens  # the domain is too generic to separate


def test_a_gene_hitting_its_own_species_scores_the_maximum() -> None:
    score = host_lineage_share(_record(["EUK__SELF"] * 3), _model())

    assert score == 1.0


def test_a_gene_hitting_only_relatives_still_scores_high() -> None:
    score = host_lineage_share(_record(["EUK__REL"] * 3), _model())

    assert score >= 0.5


def test_a_distant_eukaryote_scores_low() -> None:
    assert host_lineage_share(_record(["EUK__FAR"] * 3), _model()) < 0.2


def test_viral_hits_do_not_count_towards_the_host_share() -> None:
    viral = host_lineage_share(_record(["MIRUS__V"] * 3), _model())
    mixed = host_lineage_share(_record(["EUK__SELF", "MIRUS__V", "MIRUS__V"]), _model())

    assert viral == 0.0
    assert mixed < host_lineage_share(_record(["EUK__SELF"] * 3), _model())


def test_a_model_without_tokens_never_calls_a_gene_host() -> None:
    assert host_lineage_share(_record(["EUK__SELF"]), HostSignatureModel()) == 0.0


def test_gene_taxonomy_bitscore_weighting_uses_stored_hit_bits():
    model = _model()
    model.weight_mode = "bitscore"
    gene = GeneTaxonomy(
        "g",
        "s",
        100,
        200,
        top10_targets=["EUK__SELF", "EUK__FAR"],
        top10_bits=[75.0, 25.0],
        top10_pidents=[90.0, 90.0],
    )
    assert host_lineage_share(gene, model) == 0.75


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
    model.token_weights.pop("chlorophyta")
    model.token_weights.update({f"other_clade_{index}": 50.0 for index in range(30)})
    fingerprint = aggregate_taxonomy_substrings([("EUK__SELF", 100.0, 90.0, 1e-40)], LABELS)
    assert score_host_signature(fingerprint, model) == pytest.approx(0.10447761194029859)
    assert host_lineage_share(_record(["EUK__SELF"]), model) == 1.0
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
