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
    build_host_signature_model,
    dominant_lineage_tokens,
    host_lineage,
    host_lineage_share,
    host_rank_credit,
    host_stop_available,
    lineage_key,
    score_host_signature,
    set_taxonomy_lookup,
)
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy, pORF
from virosync.pipeline.phase2.seed_extension import extend_confirmed_seed
from virosync.pipeline.taxonomy_utils import aggregate_taxonomy_substrings, iter_taxonomy_tokens

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
        lineage_weights={
            lineage_key(HOST, 3): 60.0,
            lineage_key(RELATIVE, 3): 20.0,
            lineage_key(DISTANT, 3): 5.0,
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


def test_dominant_tokens_follow_the_signature_from_the_top() -> None:
    tokens = dominant_lineage_tokens(_model())

    assert "aurantiochytrium" in tokens
    assert "chlorophyta" not in tokens  # a minority of the signature
    assert "euk" not in tokens  # the domain is too generic to separate


CHLAMYDOMONAS = "EUK|Chlorophyta|Chlorophyceae|Chlamydomonadales|Chlamydomonadaceae|Chlamydomonas|C. reinhardtii"
VOLVOX = "EUK|Chlorophyta|Chlorophyceae|Chlamydomonadales|Volvocaceae|Volvox|Volvox carteri"
TETRABAENA = "EUK|Chlorophyta|Chlorophyceae|Chlamydomonadales|Tetrabaenaceae|Tetrabaena|T. socialis"
TREBOUXIA = "EUK|Chlorophyta|Trebouxiophyceae|Chlorellales|Chlorellaceae|Chlorella|C. variabilis"
PLASTID = "PLASTID|PLASTID|PLASTID|PLASTID|PLASTID|PLASTID|PLASTID"
UNLABELLED_MITO = "MITO|unclassified|unclassified|unclassified|unclassified|unclassified|MITO__A"
OOMYCETE = "EUK|Oomycota|Peronosporea|Peronosporales|Pythiaceae|Pythium|Pythium oligandrum"
BACTERIUM = "BAC|Pseudomonadota|Gammaproteobacteria|Enterobacterales|Enterobacteriaceae|Escherichia|E. coli"
BODO = "EUK|Euglenozoa|Kinetoplastea|Bodonida|Bodonidae|Bodo|Bodo saltans"


def _tree(weights: dict[str, float]) -> HostSignatureModel:
    return HostSignatureModel(
        token_weights={"present": 1.0},
        lineage_weights={lineage_key(lineage, 3): weight for lineage, weight in weights.items()},
    )


def test_a_pooled_organelle_label_cannot_take_over_the_host_lineage() -> None:
    model = _tree({PLASTID: 3741.0, UNLABELLED_MITO: 500.0, CHLAMYDOMONAS: 2806.0, VOLVOX: 400.0})

    tokens = dominant_lineage_tokens(model)

    assert "plastid" not in tokens
    assert {"chlorophyta", "chlorophyceae", "chlamydomonadales"} <= tokens


def test_the_lineage_stops_at_the_last_level_with_a_majority() -> None:
    model = _tree({VOLVOX: 1233.0, CHLAMYDOMONAS: 977.0, TETRABAENA: 776.0})

    tokens = dominant_lineage_tokens(model)

    assert tokens == {"chlorophyta", "chlorophyceae", "chlamydomonadales"}


def test_an_even_split_is_not_a_majority() -> None:
    tokens = dominant_lineage_tokens(_tree({CHLAMYDOMONAS: 300.0, TREBOUXIA: 300.0}))

    assert tokens == {"chlorophyta"}


def test_unassigned_ranks_do_not_hide_a_clade_from_its_competitors() -> None:
    model = _tree({HOST: 2400.0, RELATIVE: 2400.0, OOMYCETE: 900.0})

    tokens = dominant_lineage_tokens(model)

    assert "bigyra" in tokens
    assert "oomycota" not in tokens


def test_each_hit_counts_once_when_the_domain_is_chosen() -> None:
    tokens = dominant_lineage_tokens(_tree({BODO: 100.0, BACTERIUM: 150.0}))

    assert "pseudomonadota" in tokens
    assert "euglenozoa" not in tokens


def test_a_child_cannot_compete_with_its_own_ancestor() -> None:
    tokens = dominant_lineage_tokens(_tree({BODO: 100.0}))

    assert {"euglenozoa", "kinetoplastea", "bodo saltans"} <= tokens


def test_a_shared_name_under_another_parent_adds_no_support() -> None:
    model = _tree(
        {
            "EUK|Alpha|Xeno|NA|NA|NA|Xeno one": 40.0,
            "EUK|Alpha|Yankee|NA|NA|NA|Yankee one": 35.0,
            "EUK|Alpha|Zulu|NA|NA|NA|Zulu one": 25.0,
            "EUK|Beta|Xeno|NA|NA|NA|Xeno two": 30.0,
        }
    )

    assert dominant_lineage_tokens(model) == {"alpha"}


def test_weight_under_another_parent_cannot_hide_a_majority() -> None:
    model = _tree(
        {
            "EUK|Alpha|Xeno|NA|NA|NA|Xeno one": 60.0,
            "EUK|Alpha|Yankee|NA|NA|NA|Yankee one": 25.0,
            "EUK|Alpha|Zulu|NA|NA|NA|Zulu one": 15.0,
            "EUK|Beta|Yankee|NA|NA|NA|Yankee two": 30.0,
        }
    )

    assert {"alpha", "xeno"} <= dominant_lineage_tokens(model)


def test_filtered_ranks_are_skipped_on_the_way_down() -> None:
    model = _tree({"BAC|Bacillota|Clostridia|Eubacteriales|Clostridiaceae|Clostridium_X|C. one": 100.0})

    tokens = dominant_lineage_tokens(model)

    assert "c. one" in tokens
    assert "clostridium_x" not in tokens


def test_viral_lineages_never_enter_the_host_lineage() -> None:
    model = _tree({VIRUS: 5000.0, CHLAMYDOMONAS: 10.0})

    assert "mirusviricota" not in dominant_lineage_tokens(model)


def test_an_organelle_only_signature_gives_no_lineage() -> None:
    assert dominant_lineage_tokens(_tree({PLASTID: 3741.0})) == frozenset()


def test_a_signature_without_positive_weight_gives_no_lineage() -> None:
    assert dominant_lineage_tokens(_tree({CHLAMYDOMONAS: 0.0})) == frozenset()


def test_the_model_records_each_hit_under_its_full_lineage() -> None:
    marker = {
        "hmm_target": "COG0085",
        "top10_prefixes": "EUK__,EUK__",
        "top10_targets": "EUK__SELF|p1,EUK__REL|p2",
        "top10_bitscores": "100,90",
        "top10_pidents": "90,80",
        "top10_evalues": "1e-40,1e-30",
    }

    model = build_host_signature_model([marker], host_prefixes={"EUK__"}, weight_mode="rank")

    assert set(model.lineage_weights) == {lineage_key(HOST, 3), lineage_key(RELATIVE, 3)}
    assert model.lineage_weights[lineage_key(HOST, 3)] > model.lineage_weights[lineage_key(RELATIVE, 3)]


def test_the_lineage_survives_a_serialization_round_trip() -> None:
    model = _tree({CHLAMYDOMONAS: 10.0})

    assert HostSignatureModel.from_dict(model.to_dict()).lineage_weights == model.lineage_weights


def test_a_name_at_the_wrong_rank_earns_no_host_share() -> None:
    labels = {"EUK__OTHER": "EUK|Rhodophyta|Florideophyceae|Corallinales|Corallinaceae|Alpha|Alpha one"}
    set_taxonomy_lookup(TaxonomyLabelLookup({**LABELS, **labels}))
    model = _tree({"EUK|Alpha|Xeno|NA|NA|NA|Xeno one": 40.0, "EUK|Alpha|Yankee|NA|NA|NA|Yankee one": 40.0})

    score = host_lineage_share(_record(["EUK__OTHER"] * 3), model)

    set_taxonomy_lookup(TaxonomyLabelLookup(LABELS))
    assert score == 0.0


def test_an_organelle_hit_earns_no_host_share() -> None:
    set_taxonomy_lookup(TaxonomyLabelLookup({**LABELS, "PLASTID__A": PLASTID}))
    model = _tree({CHLAMYDOMONAS: 10.0})

    score = host_lineage_share(_record(["PLASTID__A"] * 3), model)

    set_taxonomy_lookup(TaxonomyLabelLookup(LABELS))
    assert score == 0.0


def test_a_short_lineage_gives_a_full_share_to_a_hit_in_that_clade() -> None:
    set_taxonomy_lookup(TaxonomyLabelLookup({**LABELS, "EUK__GREEN": TREBOUXIA}))
    model = _tree({CHLAMYDOMONAS: 300.0, TREBOUXIA: 300.0})

    score = host_lineage_share(_record(["EUK__GREEN"] * 3), model)

    set_taxonomy_lookup(TaxonomyLabelLookup(LABELS))
    assert score == 1.0


SPLIT_PHYLA = {
    "EUK|Alpha|Xeno|NA|NA|NA|Xeno one": 40.0,
    "EUK|Beta|Yankee|NA|NA|NA|Yankee one": 35.0,
    "EUK|Gamma|Zulu|NA|NA|NA|Zulu one": 25.0,
}


@pytest.mark.parametrize(
    "lineage,expected",
    [
        ("EUK|Alpha|Other|NA|NA|NA|Other one", 1.0),
        ("EUK|Beta|Other|NA|NA|NA|Other one", 35 / 40),
        ("EUK|Gamma|Other|NA|NA|NA|Other one", 25 / 40),
        ("EUK|Delta|Other|NA|NA|NA|Other one", 0.0),
    ],
)
def test_a_split_level_stays_in_the_signature_with_weighted_credit(lineage: str, expected: float) -> None:
    set_taxonomy_lookup(TaxonomyLabelLookup({**LABELS, "EUK__HIT": lineage}))
    model = _tree(SPLIT_PHYLA)

    score = host_lineage_share(_record(["EUK__HIT"]), model)

    set_taxonomy_lookup(TaxonomyLabelLookup(LABELS))
    assert host_stop_available(model)
    assert score == pytest.approx(expected)


def test_a_credit_of_exactly_the_threshold_passes_it() -> None:
    set_taxonomy_lookup(TaxonomyLabelLookup({**LABELS, "EUK__BETA": "EUK|Beta|Far|NA|NA|NA|Far one"}))
    model = _tree(
        {
            "EUK|Alpha|Xeno|NA|NA|NA|Xeno one": 40.0,
            "EUK|Beta|Yankee|NA|NA|NA|Yankee one": 14.0,
            "EUK|Gamma|Zulu|NA|NA|NA|Zulu one": 40.0,
            "EUK|Delta|Other|NA|NA|NA|Other one": 6.0,
        }
    )

    score = host_lineage_share(_record(["EUK__BETA"]), model)

    set_taxonomy_lookup(TaxonomyLabelLookup(LABELS))
    assert score >= 0.35


def test_a_split_level_counts_by_the_support_of_its_leader() -> None:
    labels = {
        "EUK__UNSEEN": "EUK|Chlorophyta|Unseenclass|NA|NA|NA|Unseen one",
        "EUK__SECOND": "EUK|Chlorophyta|Bclass|NA|NA|NA|Second one",
    }
    set_taxonomy_lookup(TaxonomyLabelLookup({**LABELS, **labels}))
    model = _tree(
        {
            "EUK|Chlorophyta|Aclass|NA|NA|NA|A one": 40.0,
            "EUK|Chlorophyta|Bclass|NA|NA|NA|B one": 35.0,
            "EUK|Chlorophyta|Cclass|NA|NA|NA|C one": 25.0,
        }
    )

    unseen, second = (host_lineage_share(_record([target]), model) for target in labels)

    set_taxonomy_lookup(TaxonomyLabelLookup(LABELS))
    assert unseen == pytest.approx(1.0 / 1.4)
    assert second == pytest.approx((1.0 + 0.4 * 35 / 40) / 1.4)


def test_a_minority_clade_at_a_majority_rank_earns_partial_credit() -> None:
    set_taxonomy_lookup(TaxonomyLabelLookup({**LABELS, "EUK__MINOR": "EUK|Beta|Yankee|NA|NA|NA|Yankee two"}))
    model = _tree({"EUK|Alpha|Xeno|NA|NA|NA|Xeno one": 60.0, "EUK|Beta|Yankee|NA|NA|NA|Yankee one": 40.0})

    score = host_lineage_share(_record(["EUK__MINOR"]), model)

    set_taxonomy_lookup(TaxonomyLabelLookup(LABELS))
    assert dominant_lineage_tokens(model) == {"alpha", "xeno", "xeno one"}
    assert score == pytest.approx(0.6 * 40 / 60 / 2.6)


def test_a_missing_parent_name_does_not_split_a_clade() -> None:
    model = _tree(
        {
            "EUK|NA|Xeno|NA|NA|NA|Xeno one": 30.0,
            "EUK|Alpha|Xeno|NA|NA|NA|Xeno two": 40.0,
            "EUK|NA|Yankee|NA|NA|NA|Yankee one": 30.0,
        }
    )

    tokens = dominant_lineage_tokens(model)

    assert "xeno" in tokens
    assert "yankee" not in tokens


def test_a_hit_with_no_lower_names_does_not_vote() -> None:
    model = _tree({"EUK|NA|NA|NA|NA|NA|NA": 100.0, "EUK|Alpha|Xeno|NA|NA|NA|NA": 90.0})

    assert dominant_lineage_tokens(model) == {"alpha", "xeno"}


def test_na_is_not_a_clade_at_any_token_length() -> None:
    model = HostSignatureModel(
        token_weights={"present": 1.0},
        lineage_weights={
            lineage_key("EUK|NA|Alpha", 2): 50.0,
            lineage_key("EUK|NA|Beta", 2): 50.0,
        },
    )

    assert "na" not in dominant_lineage_tokens(model)


def test_a_model_with_tokens_but_no_lineage_tree_is_rejected() -> None:
    payload = _tree({CHLAMYDOMONAS: 10.0}).to_dict()
    del payload["lineage_weights"]

    with pytest.raises(ValueError, match="lineage_weights"):
        HostSignatureModel.from_dict(payload)


def test_an_empty_model_still_loads() -> None:
    assert HostSignatureModel.from_dict(HostSignatureModel().to_dict()).lineage_weights == {}


def test_a_contradicted_parent_stops_the_host_credit() -> None:
    labels = {"EUK__OTHER": "EUK|Beta|Sharedclass|Sharedorder|NA|NA|Other one"}
    set_taxonomy_lookup(TaxonomyLabelLookup({**LABELS, **labels}))
    model = _tree({"EUK|Alpha|Sharedclass|Sharedorder|NA|NA|Host one": 100.0})

    score = host_lineage_share(_record(["EUK__OTHER"] * 3), model)

    set_taxonomy_lookup(TaxonomyLabelLookup(LABELS))
    assert score == 0.0


def test_a_hit_with_an_unassigned_parent_keeps_its_lower_credit() -> None:
    labels = {"EUK__PARTIAL": "EUK|NA|Sharedclass|NA|NA|NA|Other one"}
    set_taxonomy_lookup(TaxonomyLabelLookup({**LABELS, **labels}))
    model = _tree({"EUK|Alpha|Sharedclass|NA|NA|NA|Host one": 100.0})

    score = host_lineage_share(_record(["EUK__PARTIAL"] * 3), model)

    set_taxonomy_lookup(TaxonomyLabelLookup(LABELS))
    assert score == pytest.approx(1 / 3)


def test_token_filtering_keeps_measuring_the_original_text() -> None:
    assert iter_taxonomy_tokens("EUK|\u0130X|Alpha", 3) == ["EUK", "Alpha"]


@pytest.mark.parametrize(
    "weights",
    [
        {CHLAMYDOMONAS: 100.0, BACTERIUM: 100.0},
        {CHLAMYDOMONAS: 40.0, BACTERIUM: 35.0, "ARC|Asgard|Lokiarchaeia|NA|NA|NA|Loki one": 25.0},
    ],
)
def test_the_heaviest_domain_is_the_host_domain(weights: dict[str, float]) -> None:
    tokens = dominant_lineage_tokens(_tree(weights))

    assert "chlorophyta" in tokens
    assert not tokens & {"pseudomonadota", "asgard"}


def test_unassigned_dominance_descends_before_a_split_rank_can_stop_the_walk() -> None:
    model = _tree(
        {
            "EUK|NA|Knownclass|NA|NA|NA|Known one": 900.0,
            "EUK|Alpha|Xeno|NA|NA|NA|Xeno one": 40.0,
            "EUK|Beta|Yankee|NA|NA|NA|Yankee one": 35.0,
            "EUK|Gamma|Zulu|NA|NA|NA|Zulu one": 25.0,
        }
    )

    tokens = dominant_lineage_tokens(model)

    assert "knownclass" in tokens
    assert not tokens & {"alpha", "beta", "gamma"}


def test_singleton_is_not_a_clade() -> None:
    model = _tree(
        {
            "EUK|Alpha|Xeno|singleton|singleton|singleton|singleton": 60.0,
            "EUK|Alpha|Yankee|singleton|singleton|singleton|singleton": 40.0,
        }
    )

    assert "singleton" not in dominant_lineage_tokens(model)


@pytest.mark.parametrize("target", ["EUK__MISSING", "EUK__ROOTLESS"])
def test_a_hit_without_a_usable_lineage_only_dilutes_the_share(target: str) -> None:
    set_taxonomy_lookup(TaxonomyLabelLookup({**LABELS, "EUK__ROOTLESS": "Alpha"}))
    model = _tree({HOST: 100.0})

    score = host_lineage_share(_record(["EUK__SELF", target]), model)

    set_taxonomy_lookup(TaxonomyLabelLookup(LABELS))
    assert score == pytest.approx(10 / 19)


def test_the_threshold_separates_one_of_three_ranks_from_two() -> None:
    labels = {
        "EUK__ONE": "EUK|Alpha|Other|Other|NA|NA|Other one",
        "EUK__TWO": "EUK|Alpha|Beta|Other|NA|NA|Other two",
    }
    set_taxonomy_lookup(TaxonomyLabelLookup({**LABELS, **labels}))
    model = _tree({"EUK|Alpha|Beta|Gamma|NA|NA|NA": 100.0})

    one, two = (host_lineage_share(_record([target]), model) for target in labels)

    set_taxonomy_lookup(TaxonomyLabelLookup(LABELS))
    assert one == pytest.approx(1 / 3)
    assert two == pytest.approx(2 / 3)
    assert one < 0.35 < two


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
    far = host_rank_credit(DISTANT, host_lineage(model), model.min_token_length)

    assert 0 < far < 0.1  # a clade the signature holds only a little of
    assert host_lineage_share(gene, model) == pytest.approx(0.75 + 0.25 * far)


@pytest.mark.parametrize(
    "lineage,available",
    [
        (CHLAMYDOMONAS, True),
        (PLASTID, False),
        ("EUK|NA|NA|NA|NA|NA|NA", False),
    ],
)
def test_the_host_stop_is_available_only_with_a_host_lineage(lineage: str, available: bool) -> None:
    seed = MergedSeed("s", 1000, 2000, seed_id="seed", confirmed_seeds=(ConfirmedSeed("seed", "s", 1000, 2000),))
    gene = pORF(id="g", scaffold="s", start=2100, end=2200, strand="+")
    taxonomy = {gene.id: GeneTaxonomy(gene.id, "s", gene.start, gene.end)}

    boundary = extend_confirmed_seed(seed, [gene], taxonomy, (0, 10000), _tree({lineage: 10.0}), 0.35)

    assert boundary.host_stop_available is available


@pytest.mark.parametrize("weight_mode", ["rank", "bitscore"])
@pytest.mark.parametrize(
    "target,viral_identity,expected_end,expected_stop",
    [
        ("EUK__SELF", None, 2000, 2100),
        ("EUK__FAR", None, 2800, None),
        ("EUK__SELF", 25.0, 2000, 2100),
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


def test_fractional_tie_does_not_create_a_majority_rank() -> None:
    model = _tree({"EUK|Alpha|First": 0.01, "EUK|Alpha|Second": 0.05, "EUK|Beta|Third": 0.06})
    lineage = host_lineage(model)
    assert len(lineage.ranks) == 1
    assert lineage.ranks[0].support == 0.5
    assert dominant_lineage_tokens(model) == frozenset()


def test_no_taxonomy_lookup_produces_no_host_model_lineage() -> None:
    marker = {
        "hmm_target": "COG0085",
        "top10_prefixes": "EUK__",
        "top10_targets": "EUK__SELF|p1",
        "top10_bitscores": "100",
        "top10_pidents": "90",
        "top10_evalues": "1e-40",
    }
    set_taxonomy_lookup(None)
    try:
        model = build_host_signature_model([marker])
    finally:
        set_taxonomy_lookup(TaxonomyLabelLookup(LABELS))
    assert model.lineage_weights == {}
    assert not host_stop_available(model)
