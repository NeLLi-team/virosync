"""Marker alternatives and deeper validation hits must not reweight the host."""

from types import SimpleNamespace

import pytest

from virosync.pipeline import host_signatures as hs


@pytest.fixture(autouse=True)
def labels(monkeypatch):
    monkeypatch.setattr(
        hs,
        "_TAXONOMY_LOOKUP",
        {
            "EUK__a": "EUK|Metamonada|NA|Trichomonadida",
            "EUK__b": "EUK|Discoba|NA|Euglenida",
        },
    )


def marker(query, target="EUK__a", scaffold="ctg"):
    return {
        "query_porf": query,
        "scaffold": scaffold,
        "hmm_target": "COG0001",
        "top10_targets": [target],
        "top10_prefixes": [target.split("__")[0] + "__"],
        "top10_bitscores": [100.0],
        "top10_pidents": [90.0],
        "top10_evalues": [1e-20],
    }


@pytest.mark.parametrize("mode", ["rank", "bitscore"])
@pytest.mark.parametrize("as_object", [False, True])
def test_alternatives_and_domains_cannot_flip_host(mode, as_object):
    original = [marker("p1|aa1-100"), marker("p2", "EUK__b"), marker("p3", "EUK__b")]
    alternatives = [dict(original[0], hmm_target="GVOGm0003"), dict(original[0], query_porf="p1|aa20-80")]
    convert = (lambda rows: [SimpleNamespace(**row) for row in rows]) if as_object else list
    baseline = hs.build_host_signature_model(convert(original), weight_mode=mode)
    expanded = hs.build_host_signature_model(convert(original + alternatives), weight_mode=mode)
    assert expanded.to_dict() == baseline.to_dict()
    assert hs.host_lineage(expanded) == hs.host_lineage(baseline)


@pytest.mark.parametrize("mode", ["rank", "bitscore"])
@pytest.mark.parametrize("late_prefix", ["EUK__", "NCLDV__", "SUPPORT__"])
def test_validation_tail_does_not_change_host_estimate(mode, late_prefix):
    record = marker("p1")
    for field in ("top10_targets", "top10_prefixes", "top10_bitscores", "top10_pidents", "top10_evalues"):
        record[field] *= 10
    expanded = {key: list(value) if isinstance(value, list) else value for key, value in record.items()}
    expanded["top10_targets"] += [late_prefix + "b"] * 10
    expanded["top10_prefixes"] += [late_prefix] * 10
    expanded["top10_bitscores"] += [1000.0] * 10
    expanded["top10_pidents"] += [10.0] * 10
    expanded["top10_evalues"] += [1.0] * 10
    kwargs = {"weight_mode": mode, "supporting_prefixes": {"SUPPORT__"}}
    assert (
        hs.build_host_signature_model([expanded], **kwargs).to_dict()
        == hs.build_host_signature_model([record], **kwargs).to_dict()
    )


@pytest.mark.parametrize("mode", ["rank", "bitscore"])
def test_host_only_below_tenth_hit_cannot_enter_estimate(mode):
    record = marker("p1")
    record.update(
        top10_prefixes=["OTHER__"] * 10 + ["EUK__"],
        top10_targets=["OTHER__x"] * 10 + ["EUK__a"],
        top10_bitscores=[100.0] * 11,
    )
    assert not hs.build_host_signature_model([record], weight_mode=mode).lineage_weights


def test_independent_proteins_and_scaffolds_remain_independent():
    one = hs.build_host_signature_model([marker("p1")])
    three = hs.build_host_signature_model([marker("p1"), marker("p2"), marker("p1", scaffold="other")])
    assert three.lineage_weights == {key: value * 3 for key, value in one.lineage_weights.items()}


def test_first_protein_record_retains_original_evidence():
    first = marker("p1|aa1-100", "EUK__b")
    later = marker("p1|aa20-80", "EUK__a")
    assert hs.build_host_signature_model([first, later]).to_dict() == hs.build_host_signature_model([first]).to_dict()
