"""Review regressions for source identity, connected growth, and safe partitions."""

import logging
from dataclasses import replace
from types import SimpleNamespace

import pytest

from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.phase2.seed_extension import (
    accept_seed_partition,
    extend_confirmed_seed,
    merge_reachable_boundaries,
)


def boundary(identity, start, end, *, rescue=False):
    return RefinedBoundary(
        "s",
        start,
        end,
        seed_id=identity,
        confirmed_seeds=(ConfirmedSeed(identity, "s", start, end),),
        seed_sources=["frameshift_rescue"] if rescue else ["hhg"],
    )


@pytest.mark.parametrize(
    "intervals",
    [
        ((100, 5000), (2500, 3000), (4000, 6000)),
        ((100, 2000), (100, 1000), (1500, 2500)),
        ((100, 3000), (500, 1000), (1000, 2000)),
    ],
)
@pytest.mark.parametrize("separate_kind", ["rescue", "CRESS"])
def test_nested_incompatible_seed_does_not_hide_compatible_overlap(intervals, separate_kind):
    left, nested, right = intervals
    independent = boundary("b", *nested, rescue=separate_kind == "rescue")
    independent.predicted_family = "CRESS" if separate_kind == "CRESS" else ""
    merged = merge_reachable_boundaries([boundary("a", *left), independent, boundary("c", *right)])
    assert len(merged) == 2
    assert {frozenset(core.seed_id for core in row.confirmed_seeds) for row in merged} == {
        frozenset({"a", "c"}),
        frozenset({"b"}),
    }


@pytest.mark.parametrize("direction", [-1, 1])
@pytest.mark.parametrize("pending", [1, 4, 5])
def test_support_gap_requires_five_pending_genes(direction, pending):
    seed = MergedSeed("s", 1000, 2000, seed_id="a", confirmed_seeds=(ConfirmedSeed("a", "s", 1000, 2000),))
    starts = [2000 + 100 * i if direction == 1 else 900 - 100 * i for i in range(pending)]
    genes = [SimpleNamespace(id=str(i), start=start, end=start + 50) for i, start in enumerate(starts)]
    tax = {gene.id: SimpleNamespace(top10_prefixes=[], top10_pidents=[]) for gene in genes}
    result = extend_confirmed_seed(seed, genes, tax, (0, 3000), None, 0.35)
    reason = result.right_stop_reason if direction == 1 else result.left_stop_reason
    assert reason == ("support_gap" if pending == 5 else "search_envelope")


@pytest.mark.parametrize("direction", [-1, 1])
@pytest.mark.parametrize("positive_position", [5, 6])
def test_support_on_fifth_gene_connects_but_sixth_does_not(direction, positive_position):
    seed = MergedSeed("s", 1000, 2000, seed_id="a", confirmed_seeds=(ConfirmedSeed("a", "s", 1000, 2000),))
    starts = [2000 + 100 * i if direction == 1 else 900 - 100 * i for i in range(positive_position)]
    genes = [SimpleNamespace(id=str(i), start=start, end=start + 50) for i, start in enumerate(starts)]
    tax = {
        gene.id: SimpleNamespace(top10_prefixes=["NCLDV__"] if i == positive_position - 1 else [], top10_pidents=[30.0])
        for i, gene in enumerate(genes)
    }
    result = extend_confirmed_seed(seed, genes, tax, (0, 3000), None, 0.35)
    assert (result.end > 2000 if direction == 1 else result.start < 1000) == (positive_position == 5)


def test_partition_updates_container_identity_and_rejects_overlaps():
    parent = boundary("a", 100, 900)
    parent.confirmed_seeds = (ConfirmedSeed("a", "s", 100, 200), ConfirmedSeed("b", "s", 700, 800))
    children = accept_seed_partition(parent, [replace(parent, start=100, end=400), replace(parent, start=500, end=900)])
    assert [child.seed_id for child in children] == ["a", "b"]
    assert (
        accept_seed_partition(parent, [replace(parent, start=100, end=600), replace(parent, start=500, end=900)])
        is None
    )


@pytest.mark.parametrize(
    "observed,has_source",
    [
        ((), True),
        ((ConfirmedSeed("extra", "s", 100, 200),), True),
        ((ConfirmedSeed("a", "s", 100, 201),), True),
        ((ConfirmedSeed("extra", "s", 100, 200),), False),
    ],
)
def test_authenticated_resume_rejects_lost_or_changed_source_set(tmp_path, monkeypatch, observed, has_source):
    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome import phase2
    from virosync.orchestration._flows.single_genome.phase2_resume_state import PHASE2_RESUME_STATE_FILENAME
    from virosync.orchestration._flows.single_genome.phase_state import PHASE2_STATE_FILENAME

    directory = tmp_path / "phase2"
    directory.mkdir()
    (directory / PHASE2_STATE_FILENAME).touch()
    (directory / PHASE2_RESUME_STATE_FILENAME).touch()
    rows = [RefinedBoundary("s", 100, 300, seed_id=core.seed_id, confirmed_seeds=(core,)) for core in observed]
    state = SimpleNamespace(
        refined_boundaries=rows, boundary_taxonomy_map={}, boundary_control_stats=None, boundary_diamond_query=None
    )
    monkeypatch.setattr(phase2, "load_phase2_resume_state", lambda _: state)
    monkeypatch.setattr(phase2, "load_phase2_state", lambda _: rows)
    monkeypatch.setattr(phase2, "build_proteome_index", lambda _: {})
    monkeypatch.setattr(
        phase2,
        "call_task",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("unexpected output task before source validation")
        ),
    )
    seed = MergedSeed("s", 100, 200, seed_id="a", confirmed_seeds=(ConfirmedSeed("a", "s", 100, 200),))
    with pytest.raises(ValueError, match="source membership"):
        phase2._run_phase2_subflow(
            masked_path=tmp_path / "genome",
            proteome_path=tmp_path / "proteome",
            merged_seeds=[seed] if has_source else [],
            validated_markers=[],
            host_signature_model=None,
            output_dir=tmp_path,
            genome_id="g",
            config=PipelineConfig().with_overrides(resume=True),
            genome_start_time=0,
            logger=logging.getLogger(__name__),
            resume_authorized=True,
        )


@pytest.mark.parametrize(
    "change",
    [{"seed_id": "renamed"}, {"start": 99}, {"end": 301}, {"confirmed_seeds": (ConfirmedSeed("a", "s", 110, 290),)}],
)
def test_phase1_checkpoint_rejects_source_container_disagreement(change):
    from virosync.orchestration._flows.single_genome.phase1_state import Phase1StateError, _seed_to_document

    seed = MergedSeed("s", 100, 300, seed_id="a", confirmed_seeds=(ConfirmedSeed("a", "s", 100, 300),))
    with pytest.raises(Phase1StateError, match="source"):
        _seed_to_document(replace(seed, **change), 0)


@pytest.mark.parametrize("field,value", [("seed_id", "renamed"), ("start", 99), ("end", 301)])
def test_phase1_checkpoint_read_rejects_source_container_disagreement(field, value):
    from virosync.orchestration._flows.single_genome.phase1_state import (
        Phase1StateError,
        _seed_from_document,
        _seed_to_document,
    )

    seed = MergedSeed("s", 100, 300, seed_id="a", confirmed_seeds=(ConfirmedSeed("a", "s", 100, 300),))
    document = _seed_to_document(seed, 0)
    document[field] = value
    with pytest.raises(Phase1StateError, match="source"):
        _seed_from_document(document, 0)


def test_child_uses_full_parent_taxonomy_coverage_after_identity_update():
    from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy, filter_taxonomy_to_boundary

    parent = boundary("a", 100, 900)
    parent.confirmed_seeds = (ConfirmedSeed("a", "s", 100, 200), ConfirmedSeed("b", "s", 700, 800))
    children = accept_seed_partition(parent, [replace(parent, start=100, end=400), replace(parent, start=500, end=900)])
    expanded = GeneTaxonomy("expanded", "s", 550, 600)
    core = GeneTaxonomy("core", "s", 750, 780)
    assert children[1].seed_id == "b"
    assert filter_taxonomy_to_boundary({"expanded": expanded, "core": core}, children[1]) == [expanded, core]


def test_partition_keeps_parent_when_original_cores_overlap():
    parent = boundary("a", 100, 900)
    parent.confirmed_seeds = (ConfirmedSeed("a", "s", 100, 400), ConfirmedSeed("b", "s", 300, 800))
    assert (
        accept_seed_partition(parent, [replace(parent, start=100, end=400), replace(parent, start=300, end=900)])
        is None
    )


def test_nested_separate_seed_does_not_hide_host_stop():
    first = boundary("a", 100, 5000)
    first.right_stop = 5500
    first.right_stop_reason = "host_stop"
    merged = merge_reachable_boundaries([first, boundary("b", 2500, 3000, rescue=True), boundary("c", 4000, 6000)])
    assert len(merged) == 3


def test_successful_ml_support_uses_shared_controller_and_host_stop(tmp_path, monkeypatch):
    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome import phase2
    from virosync.pipeline import host_signatures
    from virosync.pipeline.host_signatures import HostSignatureModel, TaxonomyLabelLookup
    from virosync.pipeline.phase2 import taxonomy_ml_refiner
    from virosync.pipeline.phase2.boundary_diamond import GenomeDiamondQuery, SeedGeneMapping

    monkeypatch.setattr(host_signatures, "_TAXONOMY_LOOKUP", TaxonomyLabelLookup({"EUK__HOST": "EUK|host_lineage"}))

    genes = [
        SimpleNamespace(id="near", start=2100, end=2200),
        SimpleNamespace(id="h1", start=2300, end=2400),
        SimpleNamespace(id="h2", start=2500, end=2600),
        SimpleNamespace(id="h3", start=2700, end=2800),
        SimpleNamespace(id="distant", start=2900, end=3000),
    ]
    taxonomy = {
        gene.id: SimpleNamespace(
            top10_prefixes=[],
            top10_pidents=[],
            top10_targets=["EUK__HOST"] if gene.id.startswith("h") else [],
        )
        for gene in genes
    }
    query = GenomeDiamondQuery(
        seed_gene_mappings={"a": SeedGeneMapping("a", "s", 1000, 2000, flank_start_bp=0, flank_end_bp=5000)}
    )
    monkeypatch.setattr(phase2, "build_proteome_index", lambda _: {"s": genes})
    monkeypatch.setattr(phase2, "_require_phase2b_gene_taxonomy_db", lambda *args, **kwargs: tmp_path / "db")
    monkeypatch.setattr(
        phase2, "_run_boundary_diamond", lambda **kwargs: phase2._BoundaryEvidence(taxonomy, None, query)
    )
    monkeypatch.setattr(phase2, "_recalculate_boundary_composition", lambda *args, **kwargs: None)
    monkeypatch.setattr(phase2, "_annotate_boundary_repeats", lambda *args, **kwargs: None)
    monkeypatch.setattr(phase2, "_write_phase2_checkpoints", lambda **kwargs: tmp_path / "boundaries.bed")
    monkeypatch.setattr(
        taxonomy_ml_refiner, "predict_taxonomy_support", lambda **kwargs: {"a": frozenset({"near", "distant"})}
    )
    seed = MergedSeed("s", 1000, 2000, seed_id="a", confirmed_seeds=(ConfirmedSeed("a", "s", 1000, 2000),))
    result = phase2._run_phase2_subflow(
        masked_path=tmp_path / "genome",
        proteome_path=tmp_path / "proteome",
        merged_seeds=[seed],
        validated_markers=[],
        host_signature_model=HostSignatureModel(token_weights={"host_lineage": 1.0}),
        output_dir=tmp_path,
        genome_id="g",
        config=PipelineConfig().with_overrides(resume=False, boundary_taxonomy_ml_enabled=True, threads=1),
        genome_start_time=0,
        logger=logging.getLogger(__name__),
    )
    assert (result.refined_boundaries[0].start, result.refined_boundaries[0].end) == (1000, 2200)
    assert result.refined_boundaries[0].right_stop == 2300
    assert result.refined_boundaries[0].right_stop_reason == "host_stop"
    assert result.refined_boundaries[0].confirmed_seeds == seed.confirmed_seeds
