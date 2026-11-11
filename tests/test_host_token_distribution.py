"""Host support uses complete weighted cellular token distributions."""

from math import sqrt
from types import SimpleNamespace

import pytest

from virosync.pipeline import host_signatures as hs


def record(targets, *, query="", bits=None):
    return SimpleNamespace(
        query_porf=query,
        top10_targets=targets,
        top10_prefixes=[target.split("__")[0] + "__" for target in targets],
        top10_bits=bits or [100] * len(targets),
    )


@pytest.fixture
def labels(monkeypatch):
    lookup = hs.TaxonomyLabelLookup(
        {
            "EUK__host": "EUK|Amoebozoa|Entamoeba",
            "EUK__second": "Eukaryota|Amoebozoa|Entamoeba",
            "EUK__relative": "EUK|Amoebozoa|Acanthamoeba",
            "EUK__reorder": "EUK|Entamoeba|NA|singleton|Amoebozoa|unknown",
            "EUK__root": "Eukaryota",
            "BAC__root": "Bacteria",
            "ARC__root": "Archaea",
            "EUK__suffix": "EUK|Amoebozoa|Entamoeba-pEVE",
            "EUK__badlabel": "NCLDV|Amoebozoa|Entamoeba",
            **{
                f"{prefix}__bad": "EUK|Amoebozoa|Entamoeba"
                for prefix in ("NCLDV", "MITO", "PLASTID", "EUK-pEVE", "UNKNOWN")
            },
        }
    )
    monkeypatch.setattr(hs, "_TAXONOMY_LOOKUP", lookup)
    return lookup


def model():
    return hs.HostSignatureModel(token_weights={"eukaryota": 10, "amoebozoa": 10, "entamoeba": 10})


def test_root_only_domain_can_support_host(labels):
    assert hs.host_lineage_share(record(["EUK__root"]), model()) == pytest.approx(1 / sqrt(3))
    assert hs.host_lineage_share(record(["BAC__root"]), model()) == 0
    assert hs.host_lineage_share(record(["ARC__root"]), model()) == 0


def test_full_distribution_and_order_invariance(labels):
    assert hs.host_lineage_share(record(["EUK__host"]), model()) == pytest.approx(1)
    assert hs.host_lineage_share(record(["EUK__second"]), model()) == pytest.approx(1)
    assert hs.host_lineage_share(record(["EUK__reorder"]), model()) == pytest.approx(1)
    assert hs.host_lineage_share(record(["EUK__relative"]), model()) == pytest.approx(2 / 3)


@pytest.mark.parametrize(
    "target",
    [
        "NCLDV__bad",
        "MITO__bad",
        "PLASTID__bad",
        "EUK-pEVE__bad",
        "UNKNOWN__bad",
        "EUK__suffix",
        "EUK__badlabel",
        "EUK__missing",
    ],
)
def test_whole_hit_exclusions_preserve_mass(labels, target):
    assert hs.host_lineage_share(record([target]), model()) == 0
    assert hs.host_lineage_share(record(["EUK__host", target]), model()) == pytest.approx(10 / 19)


def test_original_weights_pool_noise_and_deduplicate_same_query(labels):
    first = record(["BAC__root", "EUK__host", "EUK__second"], query="gene1|aa1-20")
    duplicate_projection = record(["BAC__root", "EUK__host", "EUK__second"], query="gene1|aa30-50")
    second = record(["EUK__relative", "ARC__root"], query="gene2")
    pooled = hs.build_host_signature_model([first, duplicate_projection, second])
    assert pooled.token_weights == {
        "bacteria": 10,
        "eukaryota": 27,
        "amoebozoa": 27,
        "entamoeba": 17,
        "acanthamoeba": 10,
        "archaea": 9,
    }
    assert pooled.token_counts["eukaryota"] == 3
    assert pooled.token_bits["entamoeba"] == [9, 8]


def test_repeated_tokens_sum_to_nineteen(labels):
    pooled = hs.build_host_signature_model([record(["EUK__host", "EUK__second"])])
    assert pooled.token_weights == {"eukaryota": 19, "amoebozoa": 19, "entamoeba": 19}


@pytest.mark.parametrize("targets", [["", "EUK__host"], ",EUK__host"])
def test_empty_padding_preserves_original_rank_without_inventing_mass(labels, targets):
    gene = {"top10_targets": targets, "top10_prefixes": ["", "EUK__"]}

    pooled = hs.build_host_signature_model([gene])

    assert pooled.token_weights == {"eukaryota": 9, "amoebozoa": 9, "entamoeba": 9}
    assert hs.host_lineage_share(gene, model()) == pytest.approx(1)


@pytest.mark.parametrize("target", ["", "EUK__missing", "UNKNOWN__bad"])
def test_unidentified_hit_with_bitscore_retains_original_mass(labels, target):
    gene = record([target, "EUK__host"], bits=[100, 100])

    assert hs.host_lineage_share(gene, model()) == pytest.approx(9 / 19)


@pytest.mark.parametrize("bits", [["", 20, 80], ",20,80"])
def test_empty_bitscore_slot_does_not_shift_later_hit_weights(labels, bits):
    host = model()
    host.weight_mode = "bitscore"
    gene = {
        "top10_targets": "UNKNOWN__bad,EUK__host,NCLDV__bad",
        "top10_bits": bits,
    }

    assert hs.host_lineage_share(gene, host) == pytest.approx(0.2)


def test_bitscore_mass_and_density_use_complete_hits(labels):
    host = model()
    host.weight_mode = "bitscore"
    gene = record(["EUK__host", "NCLDV__bad"], bits=[20, 80])
    assert hs.host_lineage_share(gene, host) == pytest.approx(0.2)
    assert hs.host_signature_density_evalue_weighted([gene], host)[1] == pytest.approx(0.2)


def test_token_fallback_only_for_no_learned_cellular_signature(labels, monkeypatch):
    monkeypatch.setenv("VIROSYNC_HOST_LINEAGE", "EUK|Entamoeba|Amoebozoa")
    assert hs.host_lineage_share(record(["EUK__host"]), hs.HostSignatureModel()) == pytest.approx(1)
    assert hs.host_lineage_share(
        record(["BAC__root"]), hs.HostSignatureModel(token_weights={"bacteria": 10})
    ) == pytest.approx(1)


def test_token_model_without_lineage_schema_is_valid(labels):
    restored = hs.HostSignatureModel.from_dict({"token_weights": {"eukaryota": 1}})
    assert hs.host_stop_available(restored)
    assert hs.host_lineage_share(record(["EUK__root"]), restored) == pytest.approx(1)


def test_several_markers_pool_host_and_other_eukaryote_noise(labels):
    labels["EUK__outlier"] = "EUK|Fungi|Ascomycota|Saccharomycetes|Saccharomyces"
    pooled = hs.build_host_signature_model(
        [
            record(["EUK__host"], query="h1"),
            record(["EUK__host"], query="h2"),
            record(["EUK__host"], query="h3"),
            record(["EUK__outlier"], query="outlier"),
        ]
    )
    assert pooled.token_weights == {
        "eukaryota": 40,
        "amoebozoa": 30,
        "entamoeba": 30,
        "fungi": 10,
        "ascomycota": 10,
        "saccharomycetes": 10,
        "saccharomyces": 10,
    }
    assert hs.host_lineage_share(record(["EUK__host"]), pooled) > hs.host_lineage_share(
        record(["EUK__outlier"]), pooled
    )
    assert hs.host_lineage_share(record(["EUK__outlier"]), model()) < 0.35
    assert hs.host_lineage_share(record(["EUK__root"]), model()) > 0.35


def test_missing_lookup_empty_and_organelle_only_models(labels, monkeypatch):
    monkeypatch.delenv("VIROSYNC_HOST_LINEAGE", raising=False)
    assert hs.host_lineage_share(record(["EUK__host"]), hs.HostSignatureModel()) == 0
    assert not hs.host_stop_available(hs.HostSignatureModel(token_weights={"plastid": 50}))
    monkeypatch.setattr(hs, "_TAXONOMY_LOOKUP", None)
    assert hs.build_host_signature_model([record(["EUK__host"])]).token_weights == {}
    assert hs.host_lineage_share(record(["EUK__host"]), model()) == 0


def test_domain_position_and_alias_duplication_do_not_change_weight(labels):
    labels["EUK__reverse"] = "Entamoeba|Amoebozoa|Eukaryota"
    labels["EUK__leading_filler"] = "NA|Entamoeba|EUK|Amoebozoa|Eukaryota"
    assert hs.host_lineage_share(record(["EUK__reverse"]), model()) == pytest.approx(1)
    pooled = hs.build_host_signature_model([record(["EUK__leading_filler"])])
    assert pooled.token_weights == {"eukaryota": 10, "amoebozoa": 10, "entamoeba": 10}


def test_cellular_namespace_and_label_mismatch_uses_label_distribution(labels):
    labels["EUK__bacterial_label"] = "Bacteria"
    assert hs.host_lineage_share(record(["EUK__bacterial_label"]), model()) == 0
    bacterial = hs.HostSignatureModel(token_weights={"bacteria": 10})
    assert hs.host_lineage_share(record(["EUK__bacterial_label"]), bacterial) == pytest.approx(1)


def test_model_score_tracks_inplace_weight_changes(labels):
    host = model()
    assert hs.host_lineage_share(record(["EUK__host"]), host) == pytest.approx(1)
    host.token_weights.clear()
    host.token_weights["bacteria"] = 10
    assert hs.host_lineage_share(record(["EUK__host"]), host) == 0
    assert hs.host_lineage_share(record(["BAC__root"]), host) == pytest.approx(1)


def test_repeated_names_within_one_lineage_count_once(labels):
    labels["EUK__repeat"] = "EUK|Amoebozoa|Amoebozoa|Entamoeba|Eukaryota"

    pooled = hs.build_host_signature_model([record(["EUK__repeat"])])

    assert pooled.token_weights == {"eukaryota": 10, "amoebozoa": 10, "entamoeba": 10}


def test_scoring_uses_same_canonical_token_space_on_both_sides():
    gene = hs.TaxonomyFingerprint(weighted_tokens={" EUK ": 10, "Amoebozoa": 10, "NA": 90, "MITO": 90})
    host = hs.HostSignatureModel(token_weights={"eukaryota": 10, "amoebozoa": 10, "singleton": 90})

    assert hs.score_host_signature(gene, host) == pytest.approx(1)


def test_learned_root_only_signature_can_enable_host_stopping(labels):
    pooled = hs.build_host_signature_model([record(["EUK__root"], query="host1")])

    assert pooled.token_weights == {"eukaryota": 10}
    assert hs.host_stop_available(pooled)
    assert hs.host_lineage_share(record(["EUK__root"]), pooled) == pytest.approx(1)
