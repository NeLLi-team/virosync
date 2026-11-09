"""Viral marker evidence gates host training; each eligible protein keeps its original hit weights."""

from types import SimpleNamespace

import pytest

from virosync.ablation import AblationID
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
def test_domain_records_retain_original_host_weights(mode, as_object):
    original = [marker("p1|aa1-100"), marker("p2", "EUK__b"), marker("p3", "EUK__b")]
    alternatives = [dict(original[0], hmm_target="GVOGm0003"), dict(original[0], query_porf="p1|aa20-80")]
    convert = (lambda rows: [SimpleNamespace(**row) for row in rows]) if as_object else list
    baseline = hs.build_host_signature_model(convert(original), weight_mode=mode)
    expanded = hs.build_host_signature_model(convert(original + alternatives), weight_mode=mode)
    assert (
        expanded.lineage_weights["metamonada|trichomonadida"] == baseline.lineage_weights["metamonada|trichomonadida"]
    )
    assert expanded.lineage_weights["discoba|euglenida"] == baseline.lineage_weights["discoba|euglenida"]


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


def test_repeated_full_protein_search_contributes_once():
    first = marker("p1|aa1-100", "EUK__b")
    later = marker("p1|aa20-80", "EUK__a")
    assert hs.build_host_signature_model([first, later]).to_dict() == hs.build_host_signature_model([first]).to_dict()


@pytest.mark.parametrize("masked", [frozenset(), frozenset({"self"})])
def test_host_search_keeps_prior_candidate_limit_and_duplicate_order(tmp_path, monkeypatch, masked):
    from virosync.pipeline import search_backend as backend
    from virosync.pipeline.phase1 import marker_validation as mv

    observed = []

    def fake_diamond(**kwargs):
        observed.append(kwargs)
        rows = [f"q\tEUK__same|p{i}\t1e-10\t{100 - i}\t90\t100\n" for i in range(12)]
        rows.append("q\tPLV__virus|p\t1e-10\t80\t30\t100\n")
        targets = [row.split("\t")[1] for row in rows]
        assert backend.select_target_hits(targets, 20) == [0, 12]
        kwargs["output_tsv"].write_text("".join(rows[: kwargs["max_target_seqs"]]))

    monkeypatch.setattr(backend, "_run_diamond", fake_diamond)
    path = tmp_path / "host.tsv"
    mv.run_host_diamond_on_hmm_hits(tmp_path / "q.faa", tmp_path / "db.dmnd", path, threads=1, masked_genomes=masked)
    assert observed[0]["max_target_seqs"] == (100 if masked else 10)
    assert observed[0]["extra_flags"] == ["--sensitive", "--ext-chunk-size", "400"]
    assert len(path.read_text().splitlines()) == 10
    assert all("EUK__same|" in line for line in path.read_text().splitlines())


def test_host_task_keeps_rank13_viral_marker_out_of_host_votes(tmp_path, monkeypatch):
    from virosync.orchestration import tasks
    from virosync.pipeline import search_backend as backend
    from virosync.pipeline.phase1.hhg_seeding import HMMHit

    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">ctg_1 # 1 # 300 # 1 # ID=1_1\n" + "M" * 100 + "\n")
    calls = []

    def fake_filtered(**kwargs):
        calls.append(("marker", kwargs["max_seqs"]))
        kwargs["output_tsv"].write_text("ctg_1\tPLV__v|p\t1e-20\t200\t50\t100\n")

    def fake_raw(**kwargs):
        calls.append(("host", kwargs["max_target_seqs"]))
        kwargs["output_tsv"].write_text("ctg_1\tEUK__a|p\t1e-20\t200\t80\t100\n")

    monkeypatch.setattr("virosync.pipeline.phase1.run_diamond_on_hmm_hits", fake_filtered)
    monkeypatch.setattr(backend, "run_sequence_search", fake_raw)
    hit = HMMHit("ctg_1", "PLV_MCP", 100, 1e-20, 90, 1, 100)
    marker, host_records = tasks.marker_validation_task(
        [hit], proteome, tmp_path / "db", tmp_path / "out", max_seqs=27, threads=1
    )
    assert calls == [("marker", 27)]
    assert marker[0].top10_prefixes.split(",")[0] == "PLV__"
    assert host_records == []


def test_host_sampling_excludes_pfam_overridden_records(tmp_path, monkeypatch):
    from virosync.orchestration import tasks
    from virosync.pipeline.phase1.hhg_seeding import HMMHit

    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">ctg_1 # 1 # 300 # 1 # ID=1_1\n" + "M" * 100 + "\n")

    def search(*args, **kwargs):
        target = args[2] if args else kwargs["output_tsv"]
        target.write_text("ctg_1\tEUK__a|p\t1e-20\t200\t80\t100\n")

    monkeypatch.setattr("virosync.pipeline.phase1.run_diamond_on_hmm_hits", search)
    monkeypatch.setattr("virosync.pipeline.phase1.marker_validation.run_host_diamond_on_hmm_hits", search)
    hit = HMMHit("ctg_1", "VS000753", 100, 1e-20, 90, 1, 100, pfam_override="contradicted")
    markers, host_records = tasks.marker_validation_task([hit], proteome, tmp_path / "db", tmp_path / "out", threads=1)
    assert markers == []
    assert hit.pfam_override == "contradicted"
    assert host_records == []
    assert not hs.build_host_signature_model(host_records).lineage_weights


@pytest.mark.parametrize("no_hits", [True, False])
def test_no_query_host_evidence_is_explicit_and_does_not_search(tmp_path, monkeypatch, no_hits):
    from virosync.orchestration import tasks

    monkeypatch.setattr("virosync.pipeline.phase1.extract_hmm_hit_sequences", lambda **kwargs: 0)
    monkeypatch.setattr("virosync.pipeline.phase1.run_diamond_on_hmm_hits", lambda **kwargs: pytest.fail("No query"))
    monkeypatch.setattr(
        "virosync.pipeline.phase1.marker_validation.run_host_diamond_on_hmm_hits",
        lambda *args, **kwargs: pytest.fail("No query"),
    )
    out = tmp_path / "validation"
    assert tasks.marker_validation_task([] if no_hits else [object()], tmp_path / "p", tmp_path / "db", out) == ([], [])
    assert (out / "host_diamond_top10.tsv").read_bytes() == b""
    assert (out / "host_hmm_hit_porfs.faa").read_bytes() == b""


@pytest.mark.parametrize(
    "domain,queries_expected,weights_expected",
    [
        ("Wrong", ["ctg_3"], {"discoba|euglenida": 55.0}),
        (
            "Polymerase",
            ["ctg_3"],
            {"discoba|euglenida": 55.0},
        ),
    ],
)
def test_pfam_capsid_retention_does_not_expand_host_queries_or_votes(
    tmp_path, monkeypatch, domain, queries_expected, weights_expected
):
    from virosync.orchestration import tasks
    from virosync.pipeline.phase1.hhg_seeding import HMMHit
    from virosync.pipeline.phase1.pfam_arbitration import ModelPfamAnnotation, arbitrate_hits

    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        "".join(
            f">ctg_{i} # {1000 * i + 1} # {1000 * i + 300} # 1 # ID=1_{i}\n" + "M" * 100 + "\n" for i in range(1, 4)
        )
    )
    hits = [
        HMMHit(f"ctg_{i}", model, score, 1e-20, score, 1, 100)
        for i in (1, 2)
        for model, score in [("GVOGm0003", 200), ("GVOGm0054", 100)]
    ]
    hits.append(HMMHit("ctg_3", "COG0001", 100, 1e-20, 100, 1, 100))
    selected, _ = arbitrate_hits(
        hits,
        {"ctg_1": {domain}, "ctg_2": {domain}},
        {
            "GVOGm0003": ModelPfamAnnotation(frozenset({"CapsidDomain"}), ""),
            "GVOGm0054": ModelPfamAnnotation(frozenset({"Polymerase"}), ""),
        },
    )
    queries = []

    def diamond(**kwargs):
        ids = [line[1:] for line in kwargs["query_fasta"].read_text().splitlines() if line.startswith(">")]
        queries.append(ids)
        kwargs["output_tsv"].write_text(
            "".join(
                f"{q}\tEUK__{'b' if q == 'ctg_3' else 'a'}|p{i}\t1e-20\t{300 - i}\t80\t100\n"
                for q in ids
                for i in range(10)
            )
        )

    monkeypatch.setattr("virosync.pipeline.search_backend._run_diamond", diamond)
    markers, host = tasks.marker_validation_task(selected, proteome, tmp_path / "db", tmp_path / "out", threads=1)
    assert {m.query_porf.split("|aa", 1)[0] for m in markers if m.hmm_target == "GVOGm0003"} == {"ctg_1", "ctg_2"}
    from virosync.pipeline.phase1.marker_roles import MarkerRole, decide_marker_hit_role

    assert all(
        decide_marker_hit_role(m, ablation_id=AblationID.A0).role is MarkerRole.REPORTED_CAPSID
        for m in markers
        if m.hmm_target == "GVOGm0003"
    )
    assert queries[-1] == queries_expected
    assert {m.query_porf.split("|aa", 1)[0] for m in host} == set(queries_expected)
    model = hs.build_host_signature_model(host)
    assert model.lineage_weights == weights_expected
    assert (tmp_path / "out/host_hmm_hit_porfs.faa").read_text().count(">") == len(queries_expected)


@pytest.mark.parametrize("prefix", ["NCLDV__", "GVMAG__", "PHAGE__"])
def test_retained_rank20_viral_match_excludes_non_capsid_host(tmp_path, monkeypatch, prefix):
    from virosync.orchestration import tasks
    from virosync.pipeline.phase1.hhg_seeding import HMMHit

    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">ctg_1 # 1 # 300 # 1 # ID=1_1\n" + "M" * 100 + "\n")
    calls = []

    def diamond(**kwargs):
        calls.append(kwargs)
        rows = [f"ctg_1\tEUK__g{i}|p\t1e-20\t{300 - i}\t80\t100\n" for i in range(19)]
        rows.append(f"ctg_1\t{prefix}virus|p\t1e-20\t100\t40\t100\n")
        kwargs["output_tsv"].write_text("".join(rows))

    monkeypatch.setattr("virosync.pipeline.search_backend._run_diamond", diamond)
    markers, host = tasks.marker_validation_task(
        [HMMHit("ctg_1", "GVOGm0054", 100, 1e-20, 100, 1, 100)],
        proteome,
        tmp_path / "db",
        tmp_path / "out",
        threads=1,
    )
    assert len(markers[0].top10_prefixes.split(",")) == 20
    assert host == []
    assert len(calls) == 1
    assert (tmp_path / "out/host_hmm_hit_porfs.faa").read_text() == ""
