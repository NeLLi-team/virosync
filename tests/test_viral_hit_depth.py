"""A viral hit deeper than the top 10 keeps a gene from counting as host."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from virosync.config import ApplicationConfig, ConfigError, PipelineConfig
from virosync.orchestration._flows.single_genome.manifest import _compute_config_fingerprint
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
from virosync.pipeline.phase1.viral_markers import viral_marker_gene_ids
from virosync.pipeline.phase2 import boundary_diamond, seed_extension
from virosync.pipeline.phase2.boundary_diamond import (
    BoundaryDiamondConfig,
    DiamondHit,
    GeneTaxonomy,
    GenomeDiamondQuery,
    build_gene_taxonomy_record,
    classify_cached_diamond_query,
    has_viral_hit_veto,
    pORF,
    write_taxonomy_map,
)

VIRAL_TARGET = "NCLDV__virus|mcp"


def _hit(rank: int, target: str, *, pident: float = 80.0) -> DiamondHit:
    return DiamondHit(query="g1", target=target, evalue=1e-20, bits=1000.0 - rank, pident=pident, qcov=90.0)


def _host_hits_with_viral_at(viral_rank: int, *, total: int = 60, pident: float = 80.0) -> list[DiamondHit]:
    """Return host-labelled hits whose only viral hit sits at the 1-based ``viral_rank``."""
    return [
        _hit(rank, VIRAL_TARGET, pident=pident) if rank == viral_rank else _hit(rank, f"EUK__HOST{rank}|protein")
        for rank in range(1, total + 1)
    ]


def _classify(hits: list[DiamondHit], *, viral_hit_depth: int = 10) -> GeneTaxonomy:
    proteome_index = {"ctg": [pORF("g1", "ctg", 0, 300)]}
    query = GenomeDiamondQuery(all_porf_ids=["g1"])
    return classify_cached_diamond_query(
        query,
        {"g1": hits[:10]},
        proteome_index,
        BoundaryDiamondConfig(viral_hit_depth=viral_hit_depth),
        viral_veto_hits={"g1": hits} if viral_hit_depth > 10 else None,
    )["g1"]


@pytest.mark.parametrize(("depth", "expected"), [(10, False), (14, False), (15, True), (50, True)])
def test_depth_flag_is_true_when_the_viral_hit_rank_is_within_the_depth(depth: int, expected: bool) -> None:
    taxonomy = _classify(_host_hits_with_viral_at(15), viral_hit_depth=depth)

    assert taxonomy.viral_hit_within_depth is expected


def test_top10_fields_hold_the_first_ten_hits_of_a_deeper_list() -> None:
    hits = _host_hits_with_viral_at(15)

    deep = _classify(hits, viral_hit_depth=50)
    today = _classify(hits[:10], viral_hit_depth=10)

    assert len(deep.top10_targets) == 10
    assert (
        replace(
            deep,
            viral_hit_within_depth=today.viral_hit_within_depth,
            host_competition_viral_bits=today.host_competition_viral_bits,
            host_competition_viral_target=today.host_competition_viral_target,
        )
        == today
    )
    assert deep.host_competition_assessed and today.host_competition_assessed
    assert deep.host_competition_viral_bits == 985.0
    assert today.host_competition_viral_bits is None


def test_has_viral_keeps_the_top10_rule() -> None:
    taxonomy = _classify(_host_hits_with_viral_at(15), viral_hit_depth=50)

    assert taxonomy.has_viral is False
    assert taxonomy.viral_hit_within_depth is True


def test_a_deep_viral_hit_below_the_identity_floor_does_not_set_the_flag() -> None:
    taxonomy = _classify(_host_hits_with_viral_at(15, pident=24.9), viral_hit_depth=50)

    assert taxonomy.viral_hit_within_depth is False


def test_a_gene_without_hits_has_no_depth_flag() -> None:
    assert _classify([], viral_hit_depth=50).viral_hit_within_depth is False


@pytest.mark.parametrize("viral_rank", [1, 5, 10, 11, 20])
def test_at_depth_10_the_flag_equals_the_top10_viral_call(viral_rank: int) -> None:
    taxonomy = _classify(_host_hits_with_viral_at(viral_rank, total=10), viral_hit_depth=10)

    assert taxonomy.viral_hit_within_depth is taxonomy.has_viral


def test_a_top10_viral_hit_sets_both_flags() -> None:
    taxonomy = _classify(_host_hits_with_viral_at(3), viral_hit_depth=50)

    assert (taxonomy.has_viral, taxonomy.viral_hit_within_depth) == (True, True)


@pytest.mark.parametrize(
    ("top10_viral", "deep_viral", "expected"),
    [(False, False, False), (True, False, True), (False, True, True), (True, True, True)],
)
def test_veto_is_a_top10_qualified_hit_or_the_depth_flag(top10_viral: bool, deep_viral: bool, expected: bool) -> None:
    taxonomy = GeneTaxonomy(
        "g1",
        "ctg",
        0,
        300,
        top10_prefixes=["NCLDV__" if top10_viral else "EUK__"],
        top10_pidents=[80.0],
        viral_hit_within_depth=deep_viral,
    )

    assert has_viral_hit_veto(taxonomy) is expected


@pytest.mark.parametrize(("depth", "expected_flag"), [(10, False), (50, True)])
def test_the_phase3_gene_record_carries_the_depth_flag(depth: int, expected_flag: bool) -> None:
    taxonomy = _classify(_host_hits_with_viral_at(15), viral_hit_depth=depth)

    assert build_gene_taxonomy_record(taxonomy)["viral_hit_within_depth"] is expected_flag


def _write_query_inputs(tmp_path: Path) -> tuple[GenomeDiamondQuery, Path, dict[str, list[pORF]]]:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">g1 # 1 # 90 # + # ID=ctg_1;\nMPEPTIDE\n>g2 # 101 # 190 # + # ID=ctg_2;\nMPEPTIDE\n")
    query = GenomeDiamondQuery(
        eve_porf_ids={"seed": ["g1"]},
        boundary_porf_ids={"seed": ["g2"]},
        all_porf_ids=["g1", "g2"],
    )
    proteome_index = {"ctg": [pORF("g1", "ctg", 0, 90), pORF("g2", "ctg", 100, 190)]}
    return query, proteome, proteome_index


class _FakeSearch:
    """Write the first ``max_target_seqs`` hits of every query and remember the limits."""

    def __init__(self, hits_by_porf: dict[str, list[DiamondHit]]) -> None:
        self.hits_by_porf = hits_by_porf
        self.limits: list[int] = []

    def __call__(self, **kwargs: object) -> None:
        limit = kwargs["max_target_seqs"]
        self.limits.append(limit)
        query_ids = [
            line[1:].strip() for line in Path(kwargs["query"]).read_text().splitlines() if line.startswith(">")
        ]
        lines = [
            f"{porf_id}\t{hit.target}\t{hit.evalue}\t{hit.bits}\t{hit.pident}\t{hit.qcov}\n"
            for porf_id in query_ids
            for hit in self.hits_by_porf.get(porf_id, [])[:limit]
        ]
        Path(kwargs["output"]).write_text("".join(lines))


@pytest.mark.parametrize("chunk_size", [10, 1], ids=["single_run", "chunked"])
@pytest.mark.parametrize(("depth", "expected_flag"), [(10, False), (50, True)])
def test_batched_search_keeps_viral_hit_depth_hits_per_query(
    tmp_path: Path, monkeypatch, chunk_size: int, depth: int, expected_flag: bool
) -> None:
    query, proteome, proteome_index = _write_query_inputs(tmp_path)
    search = _FakeSearch({"g1": _host_hits_with_viral_at(15)})
    monkeypatch.setattr(boundary_diamond, "run_diamond_blastp", search)
    config = BoundaryDiamondConfig(chunk_size=chunk_size, viral_hit_depth=depth)

    taxonomy = boundary_diamond.run_batched_diamond(
        query=query,
        proteome_fasta=proteome,
        diamond_db=tmp_path / "combined.dmnd",
        output_dir=tmp_path / "out",
        proteome_index=proteome_index,
        config=config,
    )

    assert set(search.limits) == ({10, depth} if depth > 10 else {10})
    assert len(taxonomy["g1"].top10_targets) == 10
    assert taxonomy["g1"].viral_hit_within_depth is expected_flag


@pytest.mark.parametrize(("depth", "expected_flag"), [(10, False), (50, True)])
def test_cached_search_slices_to_viral_hit_depth(tmp_path: Path, depth: int, expected_flag: bool) -> None:
    query, _proteome, proteome_index = _write_query_inputs(tmp_path)

    taxonomy = classify_cached_diamond_query(
        query=query,
        diamond_hits={"g1": _host_hits_with_viral_at(15)[:10]},
        viral_veto_hits={"g1": _host_hits_with_viral_at(15)} if depth > 10 else None,
        proteome_index=proteome_index,
        config=BoundaryDiamondConfig(viral_hit_depth=depth),
    )

    assert len(taxonomy["g1"].top10_targets) == 10
    assert taxonomy["g1"].viral_hit_within_depth is expected_flag


def test_cached_and_batched_classification_agree_at_every_depth(tmp_path: Path, monkeypatch) -> None:
    query, proteome, proteome_index = _write_query_inputs(tmp_path)
    hits = {"g1": _host_hits_with_viral_at(15), "g2": _host_hits_with_viral_at(40)}
    monkeypatch.setattr(boundary_diamond, "run_diamond_blastp", _FakeSearch(hits))
    config = BoundaryDiamondConfig(viral_hit_depth=30)

    batched = boundary_diamond.run_batched_diamond(
        query=query,
        proteome_fasta=proteome,
        diamond_db=tmp_path / "combined.dmnd",
        output_dir=tmp_path / "out",
        proteome_index=proteome_index,
        config=config,
    )
    cached = classify_cached_diamond_query(query, hits, proteome_index, config, viral_veto_hits=hits)

    assert batched == cached
    assert (batched["g1"].viral_hit_within_depth, batched["g2"].viral_hit_within_depth) == (True, False)


def test_taxonomy_map_preserves_original_fields_and_exports_measured_competition(tmp_path: Path) -> None:
    deep = _classify(_host_hits_with_viral_at(15), viral_hit_depth=50)
    today = _classify(_host_hits_with_viral_at(15)[:10], viral_hit_depth=10)

    write_taxonomy_map({"g1": deep}, tmp_path / "deep.tsv")
    write_taxonomy_map({"g1": today}, tmp_path / "today.tsv")
    deep_header, deep_row = (tmp_path / "deep.tsv").read_text().splitlines()
    today_header, today_row = (tmp_path / "today.tsv").read_text().splitlines()

    deep_fields = dict(zip(deep_header.split("\t"), deep_row.split("\t"), strict=True))
    today_fields = dict(zip(today_header.split("\t"), today_row.split("\t"), strict=True))
    assert deep_header == today_header
    expected_changes = {"viral_hit_within_depth", "host_competition_viral_bits", "host_competition_viral_target"}
    assert {key for key in deep_fields if deep_fields[key] != today_fields[key]} == expected_changes
    assert deep_fields["viral_hit_within_depth"] == "1"
    assert today_fields["viral_hit_within_depth"] == "0"
    assert deep_fields["host_competition_assessed"] == "1"
    assert today_fields["host_competition_viral_bits"] == "."
    assert float(deep_fields["host_competition_viral_bits"]) == 985.0
    assert deep_fields["host_competition_viral_target"] == VIRAL_TARGET


def _extend(taxonomy: dict[str, GeneTaxonomy], genes: list[pORF], monkeypatch) -> object:
    monkeypatch.setattr(seed_extension, "host_lineage_share", lambda tax, model: 1.0)
    seed = MergedSeed("ctg", 1000, 2000, seed_id="s", confirmed_seeds=(ConfirmedSeed("s", "ctg", 1000, 2000),))
    model = SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0})
    return seed_extension.extend_confirmed_seed(seed, genes, taxonomy, (0, 10000), model, 0.35)


@pytest.mark.parametrize(
    ("depth", "end", "reason"), [(10, 2200, "host_stop"), (50, 2200, "host_stop")], ids=["depth_10", "depth_50"]
)
def test_three_cellular_dominant_hosts_stop_even_with_weak_deep_viral_hits(
    monkeypatch, depth: int, end: int, reason: str
) -> None:
    genes = [
        pORF("viral", "ctg", 2100, 2200),
        pORF("h1", "ctg", 2300, 2400),
        pORF("h2", "ctg", 2500, 2600),
        pORF("h3", "ctg", 2700, 2800),
        pORF("remote", "ctg", 2900, 3100),
    ]
    hits = {
        "viral": [_hit(1, VIRAL_TARGET)],
        "h1": _host_hits_with_viral_at(15),
        "h2": _host_hits_with_viral_at(30),
        "h3": _host_hits_with_viral_at(45),
        "remote": [_hit(1, VIRAL_TARGET)],
    }
    taxonomy = classify_cached_diamond_query(
        GenomeDiamondQuery(all_porf_ids=[gene.id for gene in genes]),
        {gene_id: rows[:10] for gene_id, rows in hits.items()},
        {"ctg": genes},
        BoundaryDiamondConfig(viral_hit_depth=depth),
        viral_veto_hits=hits if depth > 10 else None,
    )

    for name in ("h1", "h2", "h3"):
        assert taxonomy[name].has_viral is False
        assert taxonomy[name].viral_hit_within_depth is (depth == 50)
        if depth == 50:
            assert taxonomy[name].host_competition_cellular_bits > taxonomy[name].host_competition_viral_bits
    boundary = _extend(taxonomy, genes, monkeypatch)

    assert (boundary.end, boundary.right_stop_reason) == (end, reason)


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"viral_hit_depth": 9}, "phase2.viral_hit_depth must be >= 10"),
        ({"viral_hit_depth": 0}, "phase2.viral_hit_depth must be >= 10"),
    ],
)
def test_viral_hit_depth_below_ten_is_rejected(overrides: dict[str, int], message: str) -> None:
    with pytest.raises(ConfigError, match=message):
        PipelineConfig.from_dict({"phase2": overrides})


def test_viral_hit_depth_of_ten_is_accepted() -> None:
    assert PipelineConfig.from_dict({"phase2": {"viral_hit_depth": 10}}).phase2.viral_hit_depth == 10


@pytest.mark.parametrize("path", ["config/orchestration.yaml", "config/orchestration_archaeal.yaml"])
def test_shipped_yaml_declares_the_default_depth(path: str) -> None:
    shipped = ApplicationConfig.from_yaml(Path(path)).pipeline

    assert "viral_hit_depth:" in Path(path).read_text()
    assert shipped.phase2.viral_hit_depth == PipelineConfig().phase2.viral_hit_depth


def test_viral_hit_depth_reaches_the_flow_and_changes_the_run_fingerprint() -> None:
    base = PipelineConfig()
    deeper = base.with_overrides(boundary_viral_hit_depth=base.phase2.viral_hit_depth + 5)

    assert deeper.to_flow_kwargs()["boundary_viral_hit_depth"] == base.phase2.viral_hit_depth + 5
    assert _compute_config_fingerprint(base.to_flow_kwargs()) != _compute_config_fingerprint(deeper.to_flow_kwargs())


def _marker(porf: str, *, status: str = "validated", target: str = "GVOGm0003") -> SimpleNamespace:
    return SimpleNamespace(query_porf=porf, validation_status=status, hmm_target=target)


def test_viral_marker_genes_are_the_validated_markers_by_base_gene_id() -> None:
    markers = [_marker("g1|aa1"), _marker("g2"), _marker("g3", target="vs000012_vp_mcp_1")]

    assert viral_marker_gene_ids(markers) == frozenset({"g1", "g2", "g3"})


@pytest.mark.parametrize(
    ("status", "target"),
    [
        ("validated_novel", "GVOGm0003"),
        ("supported", "GVOGm0003"),
        ("unvalidated", "GVOGm0003"),
        ("validated", "VS000087_Pkinase"),
        ("validated", "COG0085"),
        ("validated", "123at2759"),
    ],
    ids=["novel", "supported", "unvalidated", "host_kinase", "cog", "busco"],
)
def test_viral_marker_genes_exclude_unvalidated_and_host_markers(status: str, target: str) -> None:
    assert viral_marker_gene_ids([_marker("g1", status=status, target=target)]) == frozenset()


def test_serialized_phase1_marker_veto_is_not_top10_viral_evidence():
    from virosync.pipeline.phase3.evidence_synthesizer import _high_confidence_euk_gene_count

    taxonomy = _classify([_hit(i, f"EUK__host{i}|p") for i in range(3)])
    taxonomy.phase1_validated_viral_marker = True
    assert has_viral_hit_veto(taxonomy)
    assert taxonomy.has_viral is False
    record = build_gene_taxonomy_record(taxonomy)
    assert record["phase1_validated_viral_marker"] is True
    assert _high_confidence_euk_gene_count([record]) == 0


def test_marker_ids_preserve_frameshift_identity_and_drop_empty_base():
    assert viral_marker_gene_ids([_marker("g|frameshift_1|aa2"), _marker("|aa2")]) == frozenset({"g|frameshift_1"})


def test_deep_only_hit_does_not_promote_repeat_slice_support():
    from virosync.pipeline.phase2.terminal_repeats import _has_viral_gene

    deep = _classify(_host_hits_with_viral_at(15), viral_hit_depth=50)
    assert has_viral_hit_veto(deep)
    assert not _has_viral_gene("ctg", 0, 300, {"g1": deep})
    top10 = _classify(_host_hits_with_viral_at(1), viral_hit_depth=50)
    assert _has_viral_gene("ctg", 0, 300, {"g1": top10})


def test_phase3_preserves_serialized_marker_flag_without_fabricating_other_gene_evidence(tmp_path, monkeypatch):
    import logging

    from virosync.orchestration._flows.single_genome import phase3

    records = [
        {"porf_id": "g1", "phase1_validated_viral_marker": True},
        {"porf_id": "g2", "phase1_validated_viral_marker": False},
    ]
    evidence = SimpleNamespace(boundary_id="b", hallmarks=[], gene_taxonomy=(records, {}), interproscan=None)
    monkeypatch.setattr(phase3, "_build_boundary_evidence", lambda **kwargs: evidence)
    monkeypatch.setattr(phase3, "_prepare_candidate_annotations", lambda *args, **kwargs: {})
    monkeypatch.setattr(phase3, "_run_verification", lambda **kwargs: kwargs["gene_taxonomy_map"])
    result = phase3._verify_native_boundaries(
        [object()],
        boundary_taxonomy_map={},
        validated_markers=[],
        proteome_index={},
        boundary_diamond_query=None,
        interproscan_map={},
        use_precomputed_taxonomy=True,
        proteome_path=tmp_path / "proteome.faa",
        config=PipelineConfig(),
        logger=logging.getLogger(__name__),
        hallmark_evidence=[],
        hallmark_witnesses=[],
        output_dir=tmp_path,
        jelly_roll_map={},
        masked_path=tmp_path / "genome.fna",
        host_signatures=set(),
        host_signature_model_payload=None,
        precomputed_tmvec=None,
        tmvec_device="cpu",
    )
    assert [row["phase1_validated_viral_marker"] for row in result["b"][0]] == [True, False]
