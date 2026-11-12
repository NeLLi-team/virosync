"""Behavioral contract for seed-preserving outward growth."""

from types import SimpleNamespace
from unittest.mock import patch

import pytest

from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
from virosync.pipeline.phase2.boundary_diamond import has_viral_hit_veto
from virosync.pipeline.phase2.seed_extension import extend_confirmed_seed


def gene(name, start, end, *, viral=False, host=False, deep_viral=False):
    return SimpleNamespace(
        id=name,
        porf_id=name,
        start=start,
        end=end,
        top10_prefixes=["NCLDV__"] if viral else ["EUK__"],
        top10_pidents=[30.0],
        viral_hit_within_depth=viral or deep_viral,
        phase1_validated_viral_marker=False,
        host_competition_assessed=False,
        host_competition_cellular_bits=None,
        host_competition_viral_bits=None,
        host_competition_cellular_target="",
        host_competition_viral_target="",
        host=host,
    )


def grow(genes, *, start=1000, end=2000, viral_marker_genes=frozenset()):
    seed = MergedSeed(
        "s", start, end, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", start, end),)
    )
    for item in genes:
        item.phase1_validated_viral_marker = item.id in viral_marker_genes
    taxonomy = {g.id: g for g in genes}
    model = SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0})
    with patch(
        "virosync.pipeline.phase2.seed_extension.host_lineage_share", side_effect=lambda g, m: 1.0 if g.host else 0.2
    ):
        return extend_confirmed_seed(
            seed, genes, taxonomy, (0, 10000), model, 0.35, viral_marker_genes=viral_marker_genes
        )


def test_complete_source_core_and_original_identity_survive():
    boundary = grow([gene("inside", 1300, 1700, host=True)])
    assert (boundary.start, boundary.end) == (1000, 2000)
    assert boundary.confirmed_seeds[0].seed_id == "original"
    assert (boundary.confirmed_seeds[0].start, boundary.confirmed_seeds[0].end) == (1000, 2000)


def test_host_stop_can_exclude_search_padding_without_cutting_seed_gene():
    seed = MergedSeed("s", 0, 10000, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", 4000, 5000),))
    genes = [
        gene("core", 4000, 5000, host=True),
        gene("near", 5100, 5200, viral=True),
        gene("h1", 5300, 5400, host=True),
        gene("h2", 5500, 5600, host=True),
        gene("h3", 5700, 5800, host=True),
        gene("remote", 5900, 6000, viral=True),
    ]
    model = SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0})
    with patch(
        "virosync.pipeline.phase2.seed_extension.host_lineage_share", side_effect=lambda g, m: 1.0 if g.host else 0.2
    ):
        boundary = extend_confirmed_seed(seed, genes, {g.id: g for g in genes}, (0, 10000), model, 0.35)

    assert (boundary.start, boundary.end) == (4000, 5200)
    assert boundary.right_stop == 5300
    assert boundary.confirmed_seeds == seed.confirmed_seeds
    assert (boundary.original_start, boundary.original_end) == (0, 10000)


def test_three_host_genes_stop_remote_viral_support():
    boundary = grow(
        [
            gene("v", 2100, 2200, viral=True),
            gene("h1", 2300, 2400, host=True),
            gene("h2", 2500, 2600, host=True),
            gene("h3", 2700, 2800, host=True),
            gene("remote", 2900, 3100, viral=True),
        ]
    )
    assert boundary.end == 2200
    assert boundary.right_stop_reason == "host_stop"
    assert boundary.right_stop == 2300


def test_deep_viral_hits_veto_host_stop_without_competition_assessment():
    boundary = grow(
        [
            gene("v", 2100, 2200, viral=True),
            gene("h1", 2300, 2400, host=True, deep_viral=True),
            gene("h2", 2500, 2600, host=True, deep_viral=True),
            gene("h3", 2700, 2800, host=True, deep_viral=True),
            gene("remote", 2900, 3100, viral=True),
        ]
    )
    assert boundary.end == 3100
    assert boundary.right_stop is None


def test_one_deep_viral_gene_breaks_host_run():
    boundary = grow(
        [
            gene("v", 2100, 2200, viral=True),
            gene("h1", 2300, 2400, host=True),
            gene("h2", 2500, 2600, host=True, deep_viral=True),
            gene("h3", 2700, 2800, host=True),
            gene("remote", 2900, 3100, viral=True),
        ]
    )
    assert boundary.end == 3100
    assert boundary.right_stop is None


def test_three_host_like_genes_that_are_validated_viral_markers_do_not_stop_extension():
    boundary = grow(
        [
            gene("v", 2100, 2200, viral=True),
            gene("h1", 2300, 2400, host=True),
            gene("h2", 2500, 2600, host=True),
            gene("h3", 2700, 2800, host=True),
            gene("remote", 2900, 3100, viral=True),
        ],
        viral_marker_genes=frozenset({"h1", "h2", "h3"}),
    )
    assert boundary.end == 3100
    assert boundary.right_stop is None


def test_a_validated_viral_marker_gene_counts_as_connected_support():
    boundary = grow([gene("m", 2100, 2200, host=True)], viral_marker_genes=frozenset({"m"}))
    assert boundary.end == 2200


def test_a_viral_marker_in_the_middle_of_a_host_run_breaks_it():
    boundary = grow(
        [
            gene("h1", 2300, 2400, host=True),
            gene("m", 2500, 2600, host=True),
            gene("h3", 2700, 2800, host=True),
            gene("remote", 2900, 3100, viral=True),
        ],
        viral_marker_genes=frozenset({"m"}),
    )
    assert boundary.end == 3100


def test_isolated_host_and_euk_best_hit_can_be_connected():
    boundary = grow([gene("h", 2100, 2200, host=True), gene("e", 2300, 2400), gene("v", 2500, 2600, viral=True)])
    assert boundary.end == 2600


def test_long_supported_gene_cannot_jump_over_overlapping_host_block():
    boundary = grow(
        [
            gene("long", 1900, 9000, viral=True),
            gene("h1", 2100, 2300, host=True),
            gene("h2", 2200, 2400, host=True),
            gene("h3", 2300, 2500, host=True),
        ]
    )
    assert boundary.end == 2000
    assert boundary.right_stop_reason == "host_stop"
    assert boundary.confirmed_seeds == (ConfirmedSeed("original", "s", 1000, 2000),)


def test_host_block_crossing_core_stops_at_core_edge():
    boundary = grow(
        [
            gene("h1", 1900, 2300, host=True),
            gene("h2", 1950, 2400, host=True),
            gene("h3", 1990, 2500, host=True),
            gene("v", 2600, 2700, viral=True),
        ]
    )
    assert boundary.end == 2000


def test_support_gap_does_not_admit_remote_positive():
    boundary = grow(
        [
            gene("a", 2100, 2200),
            gene("b", 2300, 2400),
            gene("c", 2500, 2600),
            gene("d", 2700, 2800),
            gene("e", 2900, 3000),
            gene("v", 3100, 3200, viral=True),
        ]
    )
    assert boundary.end == 2000
    assert boundary.right_stop_reason == "support_gap"


def test_missing_taxonomy_stops_growth_but_one_searched_no_hit_can_be_bridged():
    seed = MergedSeed(
        "s", 1000, 2000, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", 1000, 2000),)
    )
    unknown = gene("u", 2100, 2200)
    unknown.top10_prefixes = []
    unknown.top10_pidents = []
    positive = gene("v", 2300, 2400, viral=True)
    missing = extend_confirmed_seed(seed, [unknown, positive], {"v": positive}, (0, 10000), None, 0.35)
    searched = extend_confirmed_seed(seed, [unknown, positive], {"u": unknown, "v": positive}, (0, 10000), None, 0.35)
    assert missing.end == 2000
    assert missing.right_stop_reason == "taxonomy_unavailable"
    assert searched.end == 2400

    unknown.start, unknown.end = 2500, 2600
    hosts = [gene(f"h{i}", 2100 + i * 200, 2200 + i * 200, host=True) for i in (0, 1, 3, 4, 5)]
    with patch("virosync.pipeline.phase2.seed_extension.host_lineage_share", side_effect=lambda record, model: 1.0):
        interrupted = extend_confirmed_seed(
            seed,
            [*hosts, unknown],
            {item.id: item for item in hosts},
            (0, 10000),
            SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0}),
            0.35,
        )
    assert interrupted.end == 2000
    assert interrupted.right_stop == 2100
    assert interrupted.right_stop_reason == "host_stop"
    assert interrupted.host_stop_available is False
    assert interrupted.confirmed_seeds == seed.confirmed_seeds


def test_repeat_partition_cannot_cut_core_or_create_unseeded_residual():
    from dataclasses import replace

    from virosync.pipeline.phase2.seed_extension import accept_seed_partition

    parent = grow([gene("v", 2100, 2500, viral=True)])
    assert (
        accept_seed_partition(parent, [replace(parent, start=1100)], allowed_interval=(parent.start, parent.end))
        is None
    )
    assert (
        accept_seed_partition(
            parent,
            [replace(parent, end=2000), replace(parent, start=2000)],
            allowed_interval=(parent.start, parent.end),
        )
        is None
    )
    assert (
        accept_seed_partition(parent, [replace(parent, end=2200)], allowed_interval=(parent.start, parent.end))[
            0
        ].confirmed_seeds
        == parent.confirmed_seeds
    )


def test_independent_seeds_sharing_search_envelopes_remain_independent():
    from virosync.pipeline.phase2.seed_extension import merge_reachable_boundaries

    left = grow([], start=1000, end=2000)
    right = grow([], start=3000, end=4000)
    assert len(merge_reachable_boundaries([left, right])) == 2


def test_ablation_disables_host_stop_without_losing_core():
    seed = MergedSeed(
        "s", 1000, 2000, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", 1000, 2000),)
    )
    genes = [
        gene("h1", 2100, 2200, host=True),
        gene("h2", 2300, 2400, host=True),
        gene("h3", 2500, 2600, host=True),
        gene("v", 2700, 2800, viral=True),
    ]
    model = SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0})
    with patch(
        "virosync.pipeline.phase2.seed_extension.host_lineage_share", side_effect=lambda g, m: 1.0 if g.host else 0.2
    ):
        boundary = extend_confirmed_seed(
            seed, genes, {g.id: g for g in genes}, (0, 10000), model, 0.35, disable_host_stop=True
        )
    assert (boundary.start, boundary.end) == (1000, 2800)
    assert boundary.right_stop is None


def _stub_phase2_boundary_search(
    tmp_path, monkeypatch, *, mapped, genes=None, flank_end_bp=10000, unsearched=frozenset()
):
    """Stub the Phase 2b search so only the taxonomy-ML step varies; return the shared subflow arguments."""
    import logging

    from virosync.orchestration._flows.single_genome import phase2
    from virosync.pipeline.phase2.boundary_diamond import GenomeDiamondQuery, SeedGeneMapping

    genes = [gene("v", 8500, 9000, viral=True)] if genes is None else genes
    query = GenomeDiamondQuery(
        seed_gene_mappings={
            "original": SeedGeneMapping("original", "s", 1000, 2000, flank_start_bp=0, flank_end_bp=flank_end_bp)
        }
        if mapped
        else {}
    )
    monkeypatch.setattr(phase2, "build_proteome_index", lambda _: {"s": genes})
    monkeypatch.setattr(phase2, "_require_phase2b_gene_taxonomy_db", lambda *args, **kwargs: tmp_path / "db")

    def searched_evidence(**kwargs):
        assert [(seed.start, seed.end) for seed in kwargs["merged_seeds"]] == [(1000, 2000)]
        return phase2._BoundaryEvidence({item.id: item for item in genes if item.id not in unsearched}, None, query)

    monkeypatch.setattr(phase2, "_run_boundary_diamond", searched_evidence)
    monkeypatch.setattr(phase2, "_recalculate_boundary_composition", lambda *args, **kwargs: None)
    monkeypatch.setattr(phase2, "_annotate_boundary_repeats", lambda *args, **kwargs: None)
    monkeypatch.setattr(phase2, "_write_phase2_checkpoints", lambda **kwargs: tmp_path / "boundaries.bed")
    return dict(
        masked_path=tmp_path / "genome.fna",
        proteome_path=tmp_path / "proteome.faa",
        validated_markers=[],
        host_signature_model=None,
        output_dir=tmp_path,
        genome_id="g",
        genome_start_time=0.0,
        logger=logging.getLogger(__name__),
        metadata_families={},
    )


@pytest.mark.parametrize(
    ("family", "mapped", "expected_end"),
    [("NCLDV", True, 9000), ("CRESS", True, 1700), ("NCLDV", False, 1700)],
)
def test_phase2_searched_envelope_preserves_core_with_ml_off_or_failed(
    tmp_path, monkeypatch, family, mapped, expected_end
):
    import logging

    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome import phase2
    from virosync.pipeline.phase2 import taxonomy_ml_refiner
    from virosync.pipeline.phase2.boundary_diamond import GenomeDiamondQuery, SeedGeneMapping

    genes = [gene("v", 8500, 9000, viral=True)]
    query = GenomeDiamondQuery(
        seed_gene_mappings={
            "original": SeedGeneMapping("original", "s", 1000, 2000, flank_start_bp=0, flank_end_bp=10000)
        }
        if mapped
        else {}
    )
    monkeypatch.setattr(phase2, "build_proteome_index", lambda _: {"s": genes})
    monkeypatch.setattr(phase2, "_require_phase2b_gene_taxonomy_db", lambda *args, **kwargs: tmp_path / "db")

    def searched_evidence(**kwargs):
        assert [(seed.start, seed.end) for seed in kwargs["merged_seeds"]] == [(1000, 2000)]
        return phase2._BoundaryEvidence({"v": genes[0]}, None, query)

    monkeypatch.setattr(phase2, "_run_boundary_diamond", searched_evidence)
    monkeypatch.setattr(phase2, "_recalculate_boundary_composition", lambda *args, **kwargs: None)
    monkeypatch.setattr(phase2, "_annotate_boundary_repeats", lambda *args, **kwargs: None)
    monkeypatch.setattr(phase2, "_write_phase2_checkpoints", lambda **kwargs: tmp_path / "boundaries.bed")
    monkeypatch.setattr(
        taxonomy_ml_refiner,
        "predict_taxonomy_support",
        lambda **kwargs: (_ for _ in ()).throw(RuntimeError("fit failed")),
    )
    config = PipelineConfig().with_overrides(
        resume=False, gene_taxonomy_faa_db=tmp_path / "db", boundary_taxonomy_ml_enabled=False, threads=1
    )
    arguments = dict(
        masked_path=tmp_path / "genome.fna",
        proteome_path=tmp_path / "proteome.faa",
        validated_markers=[],
        host_signature_model=None,
        output_dir=tmp_path,
        genome_id="g",
        genome_start_time=0.0,
        logger=logging.getLogger(__name__),
    )
    off = phase2._run_phase2_subflow(
        merged_seeds=[
            MergedSeed(
                "s",
                1000,
                2000,
                seed_id="original",
                predicted_family=family,
                confirmed_seeds=(ConfirmedSeed("original", "s", 1300, 1700),),
            )
        ],
        config=config,
        **arguments,
        metadata_families={},
    )
    failed = phase2._run_phase2_subflow(
        merged_seeds=[
            MergedSeed(
                "s",
                1000,
                2000,
                seed_id="original",
                predicted_family=family,
                confirmed_seeds=(ConfirmedSeed("original", "s", 1300, 1700),),
            )
        ],
        config=config.with_overrides(boundary_taxonomy_ml_enabled=True),
        **arguments,
        metadata_families={},
    )
    assert off.refined_boundaries == failed.refined_boundaries
    assert (off.refined_boundaries[0].start, off.refined_boundaries[0].end) == (1300, expected_end)
    assert off.refined_boundaries[0].confirmed_seeds[0].end == 1700


def test_phase2_passthrough_emits_protected_gene_hull_without_search_padding():
    from virosync.orchestration._flows.single_genome.phase2 import _seeds_to_refined_boundaries

    seed = MergedSeed("s", 0, 10000, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", 4000, 5000),))

    (boundary,) = _seeds_to_refined_boundaries([seed])

    assert (boundary.start, boundary.end) == (4000, 5000)
    assert boundary.confirmed_seeds == seed.confirmed_seeds


def test_displaced_mutable_coordinates_cannot_relocate_confirmed_seed():
    from virosync.pipeline.phase1.seed_merger import ConfirmedSeed

    seed = MergedSeed(
        "s", 8000, 9000, seed_id="renumbered", confirmed_seeds=(ConfirmedSeed("phase1", "s", 1000, 2000),)
    )
    boundary = extend_confirmed_seed(seed, [], {}, (0, 10000), None, 0.35)
    assert (boundary.start, boundary.end) == (1000, 2000)
    assert boundary.confirmed_seeds[0].seed_id == "phase1"


def test_ml_positive_host_genes_without_viral_identity_still_form_host_stop():
    seed = MergedSeed(
        "s", 1000, 2000, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", 1000, 2000),)
    )
    genes = [
        gene("h1", 2100, 2200, host=True),
        gene("h2", 2300, 2400, host=True),
        gene("h3", 2500, 2600, host=True),
        gene("ml", 2700, 2800),
    ]
    with patch(
        "virosync.pipeline.phase2.seed_extension.host_lineage_share", side_effect=lambda g, m: 1.0 if g.host else 0.2
    ):
        boundary = extend_confirmed_seed(
            seed,
            genes,
            {g.id: g for g in genes},
            (0, 10000),
            SimpleNamespace(token_weights={"h": 1}, lineage_weights={"euk|h": 1.0}),
            0.35,
            ml_support=frozenset({"ml", "h1", "h2", "h3"}),
        )
    assert boundary.end == 2000
    assert boundary.right_stop == 2100


def test_overlap_merge_preserves_all_original_membership():
    from virosync.pipeline.phase1.seed_merger import ConfirmedSeed
    from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
    from virosync.pipeline.phase2.seed_extension import merge_reachable_boundaries

    left = RefinedBoundary("s", 1000, 3000, confirmed_seeds=(ConfirmedSeed("a", "s", 1000, 2000),))
    right = RefinedBoundary("s", 2500, 4000, confirmed_seeds=(ConfirmedSeed("b", "s", 3000, 4000),))
    merged = merge_reachable_boundaries([left, right])
    assert len(merged) == 1
    assert merged[0].confirmed_seeds == left.confirmed_seeds + right.confirmed_seeds
    right.left_stop = 2500
    assert len(merge_reachable_boundaries([left, right])) == 2


def test_ml_features_do_not_contain_host_model_outputs():
    from virosync.pipeline.phase2.taxonomy_ml_refiner import FEATURE_COLUMNS

    assert "host_signature_score" not in FEATURE_COLUMNS
    assert "neighbor_host_score_mean" not in FEATURE_COLUMNS


def test_ml_fit_returns_gene_support_without_coordinate_changes(tmp_path):
    from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy, GenomeDiamondQuery, SeedGeneMapping
    from virosync.pipeline.phase2.taxonomy_ml_refiner import predict_taxonomy_support

    positive = GeneTaxonomy(
        "core",
        "s",
        1200,
        1300,
        top1_prefix="NCLDV__",
        top10_prefixes=["NCLDV__"],
        top10_pidents=[50.0],
        top10_bits=[100.0],
        has_hit=True,
        has_viral=True,
    )
    flank = GeneTaxonomy(
        "flank",
        "s",
        2300,
        2400,
        top1_prefix="NCLDV__",
        top10_prefixes=["NCLDV__"],
        top10_pidents=[50.0],
        top10_bits=[100.0],
        has_hit=True,
        has_viral=True,
    )
    control = GeneTaxonomy(
        "control",
        "s",
        8000,
        8100,
        top1_prefix="EUK__",
        top10_prefixes=["EUK__"],
        top10_pidents=[90.0],
        top10_bits=[100.0],
        has_hit=True,
    )
    mapping = SeedGeneMapping("original", "s", 1000, 2000, eve_porf_ids=["core"], downstream_porf_ids=["flank"])
    query = GenomeDiamondQuery(seed_gene_mappings={"original": mapping}, control_porf_ids=["control"])
    support = predict_taxonomy_support(
        {"core": positive, "flank": flank, "control": control}, query, output_dir=tmp_path
    )
    assert "flank" in support["original"]
    assert (mapping.seed_start, mapping.seed_end) == (1000, 2000)
    assert "host_signature_score" not in (tmp_path / "feature_table.tsv").read_text()


def test_controller_rejects_missing_confirmed_membership():
    import pytest

    with pytest.raises(ValueError, match="no confirmed seed"):
        extend_confirmed_seed(MergedSeed("s", 1000, 2000, seed_id="late"), [], {}, (0, 10000), None, 0.35)


def test_left_three_host_run_stops_remote_viral_support():
    result = grow(
        [
            gene("near", 850, 950, viral=True),
            gene("h1", 700, 800, host=True),
            gene("h2", 500, 600, host=True),
            gene("h3", 300, 400, host=True),
            gene("remote", 100, 200, viral=True),
        ]
    )
    assert result.start == 850
    assert result.left_stop == 800
    assert result.left_stop_reason == "host_stop"


def test_left_crossing_supported_gene_cannot_jump_host_block():
    result = grow(
        [
            gene("long", 100, 1100, viral=True),
            gene("h1", 800, 900, host=True),
            gene("h2", 700, 800, host=True),
            gene("h3", 600, 700, host=True),
        ]
    )
    assert result.start == 1000
    assert result.left_stop == 900
    assert result.left_stop_reason == "host_stop"


def test_independent_seed_outputs_can_overlap_without_loss_or_stop_crossing():
    from virosync.pipeline.phase2.seed_extension import merge_reachable_boundaries

    left_seed = MergedSeed("s", 1000, 2000, seed_id="left", confirmed_seeds=(ConfirmedSeed("left", "s", 1000, 2000),))
    right_seed = MergedSeed(
        "s", 3000, 5000, seed_id="right", confirmed_seeds=(ConfirmedSeed("right", "s", 3000, 5000),)
    )
    genes = [
        gene("v", 2100, 3200, viral=True),
        gene("h1", 3300, 3400, host=True),
        gene("h2", 3500, 3600, host=True),
        gene("h3", 3700, 3800, host=True),
    ]
    taxonomy = {item.id: item for item in genes}
    with patch(
        "virosync.pipeline.phase2.seed_extension.host_lineage_share",
        side_effect=lambda record, model: 1.0 if record.host else 0.2,
    ):
        left = extend_confirmed_seed(
            left_seed,
            genes,
            taxonomy,
            (0, 6000),
            SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0}),
            0.35,
        )
        right = extend_confirmed_seed(
            right_seed,
            genes,
            taxonomy,
            (0, 6000),
            SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0}),
            0.35,
        )
    outputs = merge_reachable_boundaries([left, right])
    assert [(item.start, item.end) for item in outputs] == [(1000, 3200), (2100, 5000)]
    assert outputs[0].end <= outputs[0].right_stop == 3300
    assert outputs[0].confirmed_seeds == left_seed.confirmed_seeds
    assert outputs[1].confirmed_seeds == right_seed.confirmed_seeds


@pytest.mark.parametrize("missing_taxonomy", [False, True])
def test_supported_gene_crossing_limit_does_not_commit_partial_growth(missing_taxonomy):
    seed = MergedSeed(
        "s", 1000, 2000, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", 1000, 2000),)
    )
    near = gene("near", 2100, 2300, viral=True)
    long = gene("long", 2200, 6000, viral=True)
    unknown = gene("unknown", 3000, 3100)
    genes = [near, long, unknown] if missing_taxonomy else [near, long]
    result = extend_confirmed_seed(seed, genes, {"near": near, "long": long}, (0, 5000), None, 0.35)
    assert (result.start, result.end) == (1000, 2300)
    assert result.right_stop_reason == ("taxonomy_unavailable" if missing_taxonomy else "search_envelope")
    assert result.confirmed_seeds == seed.confirmed_seeds


def test_three_seed_merge_keeps_stop_inside_independent_core_and_exact_ownership():
    from virosync.pipeline.phase2.seed_extension import merge_reachable_boundaries, validate_source_membership

    left = MergedSeed("s", 1000, 2000, seed_id="left", confirmed_seeds=(ConfirmedSeed("left", "s", 1000, 2000),))
    middle = MergedSeed("s", 3000, 5000, seed_id="middle", confirmed_seeds=(ConfirmedSeed("middle", "s", 3000, 5000),))
    right = MergedSeed("s", 5500, 6500, seed_id="right", confirmed_seeds=(ConfirmedSeed("right", "s", 5500, 6500),))
    genes = [
        gene("left", 800, 900, viral=True),
        gene("shared", 2100, 3200, viral=True),
        gene("h1", 3300, 3400, host=True),
        gene("h2", 3500, 3600, host=True),
        gene("h3", 3700, 3800, host=True),
        gene("bridge", 4800, 5600, viral=True),
        gene("right", 6600, 6700, viral=True),
    ]
    with patch(
        "virosync.pipeline.phase2.seed_extension.host_lineage_share", side_effect=lambda g, m: 1.0 if g.host else 0.2
    ):
        boundaries = [
            extend_confirmed_seed(
                seed,
                genes,
                {g.id: g for g in genes},
                envelope,
                SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0}),
                0.35,
            )
            for seed, envelope in zip((left, middle, right), ((0, 4000), (2100, 6700), (4800, 7000)), strict=True)
        ]
    outputs = merge_reachable_boundaries(boundaries)
    assert [(item.start, item.end) for item in outputs] == [(800, 3200), (2100, 6700)]
    assert outputs[0].right_stop == 3300
    assert middle.start < outputs[0].right_stop < middle.end
    assert outputs[0].confirmed_seeds == left.confirmed_seeds
    assert outputs[1].confirmed_seeds == middle.confirmed_seeds + right.confirmed_seeds
    validate_source_membership([left, middle, right], outputs)


def test_deep_only_hit_is_not_positive_extension_support():
    boundary = grow([gene("deep", 2100, 2200, host=True, deep_viral=True)])
    assert boundary.end == 2000
    assert boundary.right_stop is None


@pytest.mark.parametrize("direction", [1, -1])
@pytest.mark.parametrize("middle_kind", ["neutral", "nohit", "deep_viral"])
@pytest.mark.parametrize("tail_hosts", [2, 3])
def test_neutral_genes_preserve_host_count_and_deep_viral_resets_it(direction, middle_kind, tail_hosts):
    from virosync.pipeline.phase3.joint_seed_rescue import _gap_has_host_barrier

    genes = [gene(f"h{i}", 2100 + i * 200, 2200 + i * 200, host=True) for i in range(3 + tail_hosts)]
    middle = genes[2]
    middle.host = middle_kind == "deep_viral"
    middle.viral_hit_within_depth = middle_kind == "deep_viral"
    if middle_kind == "nohit":
        middle.top10_prefixes = []
        middle.top10_pidents = []
    if direction == -1:
        for item in genes:
            item.start, item.end = 10000 - item.end, 10000 - item.start
    start, end = (1000, 2000) if direction == 1 else (8000, 9000)
    boundary = grow(genes, start=start, end=end)
    stop = boundary.right_stop if direction == 1 else boundary.left_stop
    expected_stop = (
        (2100 if direction == 1 else 7900)
        if middle_kind != "deep_viral"
        else ((2700 if direction == 1 else 7300) if tail_hosts == 3 else None)
    )
    assert stop == expected_stop
    assert (boundary.start, boundary.end) == (start, end)
    assert boundary.confirmed_seeds == (ConfirmedSeed("original", "s", start, end),)
    children = (SimpleNamespace(end=0), SimpleNamespace(start=10000))
    with patch(
        "virosync.pipeline.phase3.joint_seed_rescue.host_lineage_share",
        side_effect=lambda record, model: 1.0 if record.host else 0.2,
    ):
        blocked = _gap_has_host_barrier(children, genes, {g.id: g for g in genes}, None, 0.35, frozenset())
    assert blocked is (middle_kind != "deep_viral" or tail_hosts == 3)


@pytest.mark.parametrize("pident,expected_veto,expected_end", [(24.9, False, 2000), (25.0, True, 2200)])
def test_host_like_top10_viral_hit_vetoes_and_supports_extension_at_the_identity_floor(
    pident, expected_veto, expected_end
):
    hit = gene("mixed", 2100, 2200, host=True, viral=True)
    hit.top10_prefixes = ["EUK__"] * 9 + ["NCLDV__"]
    hit.top10_pidents = [90.0] * 9 + [pident]
    hit.viral_hit_within_depth = False
    assert has_viral_hit_veto(hit) is expected_veto
    boundary = grow([hit])
    assert boundary.end == expected_end
    assert boundary.right_stop is None


@pytest.mark.parametrize("direction", [1, -1])
def test_five_searched_no_hit_fragments_do_not_discard_connected_viral_span(direction):
    genes = [gene("near", 2100, 2200, viral=True)]
    genes.extend(gene(f"n{i}", 2300 + i * 200, 2400 + i * 200) for i in range(5))
    for item in genes[1:]:
        item.top10_prefixes = []
        item.top10_pidents = []
    genes.append(gene("remote", 3300, 3400, viral=True))
    if direction == -1:
        for item in genes:
            item.start, item.end = 10000 - item.end, 10000 - item.start
    start, end = (1000, 2000) if direction == 1 else (8000, 9000)
    boundary = grow(genes, start=start, end=end)
    assert (boundary.start, boundary.end) == ((1000, 3400) if direction == 1 else (6600, 9000))
    assert (boundary.right_stop if direction == 1 else boundary.left_stop) is None
    assert boundary.confirmed_seeds == (ConfirmedSeed("original", "s", start, end),)


@pytest.mark.parametrize("support", ["marker", "ml"])
def test_searched_no_hit_marker_and_ml_support_are_not_skipped(support):
    seed = MergedSeed(
        "s", 1000, 2000, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", 1000, 2000),)
    )
    hit = gene("supported", 2100, 2200)
    hit.top10_prefixes = []
    hit.top10_pidents = []
    hit.phase1_validated_viral_marker = support == "marker"
    boundary = extend_confirmed_seed(
        seed,
        [hit],
        {hit.id: hit},
        (0, 10000),
        None,
        0.35,
        ml_support=frozenset({hit.id}) if support == "ml" else frozenset(),
    )
    assert boundary.end == 2200
    assert boundary.right_stop is None
    assert boundary.confirmed_seeds == seed.confirmed_seeds


@pytest.mark.parametrize("direction", [1, -1])
@pytest.mark.parametrize("reset_kind", ["deep_viral", "marker_member"])
def test_empty_top10_independent_evidence_resets_host_run_without_growth(direction, reset_kind):
    genes = [gene(f"h{i}", 2100 + i * 200, 2200 + i * 200, host=True) for i in range(5)]
    middle = genes[2]
    middle.top10_prefixes = []
    middle.top10_pidents = []
    middle.viral_hit_within_depth = reset_kind == "deep_viral"
    assert middle.phase1_validated_viral_marker is False
    marker_ids = frozenset({middle.id}) if reset_kind == "marker_member" else frozenset()
    if direction == -1:
        for item in genes:
            item.start, item.end = 10000 - item.end, 10000 - item.start
    start, end = (1000, 2000) if direction == 1 else (8000, 9000)
    seed = MergedSeed(
        "s", start, end, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", start, end),)
    )
    model = SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0})
    with patch("virosync.pipeline.phase2.seed_extension.host_lineage_share", return_value=1.0):
        boundary = extend_confirmed_seed(
            seed,
            genes,
            {item.id: item for item in genes},
            (0, 10000),
            model,
            0.35,
            viral_marker_genes=marker_ids,
        )
    assert (boundary.right_stop if direction == 1 else boundary.left_stop) is None
    assert (boundary.start, boundary.end) == (start, end)
    assert boundary.confirmed_seeds == seed.confirmed_seeds
    assert boundary.host_stop_available is True


@pytest.mark.parametrize("direction", [1, -1])
@pytest.mark.parametrize("support_kind", ["foreign", "ml"])
def test_nonviral_growth_support_does_not_reset_host_count(direction, support_kind):
    genes = [gene(f"h{i}", 2100 + i * 200, 2200 + i * 200, host=True) for i in range(4)]
    middle = genes[2]
    middle.host = False
    if direction == -1:
        for item in genes:
            item.start, item.end = 10000 - item.end, 10000 - item.start
    start, end = (1000, 2000) if direction == 1 else (8000, 9000)
    seed = MergedSeed(
        "s", start, end, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", start, end),)
    )
    model = SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0})
    support = frozenset({middle.id}) if support_kind == "ml" else frozenset()
    with patch(
        "virosync.pipeline.phase2.seed_extension.host_lineage_share",
        side_effect=lambda record, model: 1.0 if record.host else (0.0 if support_kind == "foreign" else 0.2),
    ):
        unblocked = extend_confirmed_seed(
            seed, [middle], {middle.id: middle}, (0, 10000), model, 0.35, ml_support=support
        )
        blocked = extend_confirmed_seed(
            seed, genes, {item.id: item for item in genes}, (0, 10000), model, 0.35, ml_support=support
        )

    assert (unblocked.start, unblocked.end) == ((1000, 2600) if direction == 1 else (7400, 9000))
    assert (blocked.right_stop if direction == 1 else blocked.left_stop) == (2100 if direction == 1 else 7900)
    assert (blocked.right_stop_reason if direction == 1 else blocked.left_stop_reason) == "host_stop"
    assert (blocked.start, blocked.end) == (start, end)
    assert blocked.confirmed_seeds == seed.confirmed_seeds


def test_nonviral_foreign_support_does_not_reset_joint_host_count():
    from virosync.pipeline.phase3.joint_seed_rescue import _gap_has_host_barrier

    genes = [
        gene("h1", 2100, 2200, host=True),
        gene("h2", 2300, 2400, host=True),
        gene("foreign", 2500, 2600),
        gene("h3", 2700, 2800, host=True),
    ]
    children = (SimpleNamespace(end=2000), SimpleNamespace(start=3000))
    with patch(
        "virosync.pipeline.phase3.joint_seed_rescue.host_lineage_share",
        side_effect=lambda record, model: 1.0 if record.host else 0.0,
    ):
        blocked = _gap_has_host_barrier(children, genes, {g.id: g for g in genes}, None, 0.35, frozenset())

    assert blocked is True


def foreign_gene(name, start, end):
    """Return a gene with cellular hits and a host share of 0.0, so it counts as foreign."""
    item = gene(name, start, end)
    item.share = 0.0
    return item


def no_hit_gene(name, start, end):
    """Return a searched gene without any hit."""
    item = gene(name, start, end)
    item.top10_prefixes = []
    item.top10_pidents = []
    return item


def extend(genes, *, start=1000, end=2000, taxonomy=None, envelope=(0, 10000), **kwargs):
    """Extend a seed; a gene flagged host has share 1.0, a gene with a share attribute that share, others 0.2."""
    seed = MergedSeed(
        "s", start, end, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", start, end),)
    )
    model = SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0})
    records = {item.id: item for item in genes} if taxonomy is None else taxonomy
    with patch(
        "virosync.pipeline.phase2.seed_extension.host_lineage_share",
        side_effect=lambda record, _model: getattr(record, "share", 1.0 if record.host else 0.2),
    ):
        return extend_confirmed_seed(seed, genes, records, envelope, model, 0.35, **kwargs)


def test_host_like_gene_with_a_viral_top10_hit_supports_growth():
    boundary = extend([gene("m", 2100, 2200, host=True, viral=True)])
    assert boundary.end == 2200


def test_five_host_like_genes_with_viral_top10_hits_do_not_end_the_walk():
    boundary = extend(
        [
            gene("m1", 2100, 2200, host=True, viral=True),
            gene("m2", 2300, 2400, host=True, viral=True),
            gene("m3", 2500, 2600, host=True, viral=True),
            gene("m4", 2700, 2800, host=True, viral=True),
            gene("m5", 2900, 3000, host=True, viral=True),
        ]
    )
    assert boundary.end == 3000
    assert boundary.right_stop_reason == "search_envelope"


def test_marker_set_membership_alone_does_not_support_growth_at_a_host_like_share():
    boundary = extend([gene("m", 2100, 2200, host=True)], viral_marker_genes=frozenset({"m"}))
    assert boundary.end == 2000


def test_marker_set_membership_supports_growth_below_the_host_share_threshold():
    boundary = extend([gene("m", 2100, 2200)], viral_marker_genes=frozenset({"m"}))
    assert boundary.end == 2200


def real_stop_genes():
    """Return a gene that ends inside the first host gene of a stop, then three host genes."""
    return [
        gene("v", 72678, 73992, viral=True),
        gene("h1", 73980, 75960, host=True),
        gene("h2", 76000, 78466, host=True),
        gene("h3", 78488, 78890, host=True),
    ]


def test_supporting_gene_ending_inside_the_first_host_gene_is_accepted():
    boundary = extend(real_stop_genes(), start=70000, end=72000, envelope=(0, 100000))
    assert boundary.end == 73992


def test_right_stop_after_a_gene_ending_inside_the_first_host_gene_equals_the_region_end():
    boundary = extend(real_stop_genes(), start=70000, end=72000, envelope=(0, 100000))
    assert boundary.right_stop == boundary.end
    assert boundary.right_stop_reason == "host_stop"


def mirrored_real_stop_genes():
    """Return the genes of real_stop_genes mirrored to the left of the seed (x becomes 200000 - x)."""
    return [
        gene("v", 126008, 127322, viral=True),
        gene("h1", 124040, 126020, host=True),
        gene("h2", 121534, 124000, host=True),
        gene("h3", 121110, 121512, host=True),
    ]


def test_supporting_gene_ending_inside_the_first_host_gene_on_the_left_is_accepted():
    boundary = extend(mirrored_real_stop_genes(), start=128000, end=130000, envelope=(0, 200000))
    assert boundary.start == 126008


def test_left_stop_after_a_gene_ending_inside_the_first_host_gene_equals_the_region_start():
    boundary = extend(mirrored_real_stop_genes(), start=128000, end=130000, envelope=(0, 200000))
    assert boundary.left_stop == boundary.start
    assert boundary.left_stop_reason == "host_stop"


def test_supporting_gene_ending_beyond_the_first_host_gene_is_not_accepted():
    genes = [gene("v", 72678, 76500, viral=True), *real_stop_genes()[1:]]
    boundary = extend(genes, start=70000, end=72000, envelope=(0, 100000))
    assert boundary.end == 72000
    assert boundary.right_stop == 73980


def test_supporting_gene_ending_inside_the_first_missing_record_is_rejected():
    supported = gene("v", 2100, 2500, viral=True)
    boundary = extend([supported, gene("u", 2400, 2600)], taxonomy={"v": supported})
    assert boundary.end == 2000
    assert boundary.right_stop_reason == "taxonomy_unavailable"


def test_supporting_gene_ending_beyond_the_first_gene_without_a_record_is_not_accepted():
    supported = gene("v", 2100, 2700, viral=True)
    boundary = extend([supported, gene("u", 2400, 2600)], taxonomy={"v": supported})
    assert boundary.end == 2000


def test_supporting_gene_crossing_the_envelope_edge_is_not_accepted():
    boundary = extend([gene("v", 2400, 3200, viral=True)], envelope=(0, 3000))
    assert boundary.end == 2000
    assert boundary.right_stop_reason == "search_envelope"


def test_supporting_gene_crossing_the_walk_limit_is_not_accepted():
    boundary = extend([gene("v", 2400, 3200, viral=True)], envelope=(0, 2500), walk_limit=(0, 3000))
    assert boundary.end == 2000
    assert boundary.right_stop_reason == "search_envelope"


def test_three_host_genes_followed_by_two_supporting_genes_in_five_are_bridged():
    boundary = extend(
        [
            gene("h1", 2100, 2200, host=True),
            gene("h2", 2300, 2400, host=True),
            gene("h3", 2500, 2600, host=True),
            gene("s1", 2700, 2800, viral=True),
            gene("n1", 2900, 3000),
            gene("s2", 3100, 3200, viral=True),
        ]
    )
    assert boundary.end == 3200
    assert boundary.right_stop is None


def test_three_host_genes_followed_by_two_supporting_genes_on_the_left_are_bridged():
    boundary = extend(
        [
            gene("h1", 7800, 7900, host=True),
            gene("h2", 7600, 7700, host=True),
            gene("h3", 7400, 7500, host=True),
            gene("s1", 7200, 7300, viral=True),
            gene("n1", 7000, 7100),
            gene("s2", 6800, 6900, viral=True),
        ],
        start=8000,
        end=9000,
    )
    assert boundary.start == 6800
    assert boundary.left_stop is None


def test_three_host_genes_followed_by_one_supporting_gene_stop_the_walk():
    boundary = extend(
        [
            gene("h1", 2100, 2200, host=True),
            gene("h2", 2300, 2400, host=True),
            gene("h3", 2500, 2600, host=True),
            gene("s1", 2700, 2800, viral=True),
            gene("n1", 2900, 3000),
        ]
    )
    assert boundary.end == 2000
    assert boundary.right_stop == 2100
    assert boundary.right_stop_reason == "host_stop"


def marked_marker_gene(name, start, end):
    """Return a gene without a viral hit that Phase 1 validated as a viral marker."""
    item = gene(name, start, end)
    item.phase1_validated_viral_marker = True
    return item


def test_two_foreign_supporters_do_not_bridge_a_host_run():
    boundary = extend(
        [
            gene("h1", 2100, 2200, host=True),
            gene("h2", 2300, 2400, host=True),
            gene("h3", 2500, 2600, host=True),
            foreign_gene("f1", 2700, 2800),
            foreign_gene("f2", 2900, 3000),
        ]
    )
    assert boundary.end == 2000
    assert boundary.right_stop == 2100


def test_two_foreign_supporters_do_not_bridge_a_host_run_on_the_left():
    boundary = extend(
        [
            gene("h1", 7800, 7900, host=True),
            gene("h2", 7600, 7700, host=True),
            gene("h3", 7400, 7500, host=True),
            foreign_gene("f1", 7200, 7300),
            foreign_gene("f2", 7000, 7100),
        ],
        start=8000,
        end=9000,
    )
    assert boundary.start == 8000
    assert boundary.left_stop == 7900


def test_two_ml_supported_genes_do_not_bridge_a_host_run():
    boundary = extend(
        [
            gene("h1", 2100, 2200, host=True),
            gene("h2", 2300, 2400, host=True),
            gene("h3", 2500, 2600, host=True),
            gene("m1", 2700, 2800),
            gene("m2", 2900, 3000),
        ],
        ml_support=frozenset({"m1", "m2"}),
    )
    assert boundary.end == 2000
    assert boundary.right_stop == 2100


def test_two_ml_supported_genes_do_not_bridge_a_host_run_on_the_left():
    boundary = extend(
        [
            gene("h1", 7800, 7900, host=True),
            gene("h2", 7600, 7700, host=True),
            gene("h3", 7400, 7500, host=True),
            gene("m1", 7200, 7300),
            gene("m2", 7000, 7100),
        ],
        start=8000,
        end=9000,
        ml_support=frozenset({"m1", "m2"}),
    )
    assert boundary.start == 8000
    assert boundary.left_stop == 7900


def test_two_marker_set_members_without_a_viral_match_do_not_bridge_a_host_run():
    boundary = extend(
        [
            gene("h1", 2100, 2200, host=True),
            gene("h2", 2300, 2400, host=True),
            gene("h3", 2500, 2600, host=True),
            gene("k1", 2700, 2800),
            gene("k2", 2900, 3000),
        ],
        viral_marker_genes=frozenset({"k1", "k2"}),
    )
    assert boundary.end == 2000
    assert boundary.right_stop == 2100


def test_two_marker_set_members_without_a_viral_match_do_not_bridge_a_host_run_on_the_left():
    boundary = extend(
        [
            gene("h1", 7800, 7900, host=True),
            gene("h2", 7600, 7700, host=True),
            gene("h3", 7400, 7500, host=True),
            gene("k1", 7200, 7300),
            gene("k2", 7000, 7100),
        ],
        start=8000,
        end=9000,
        viral_marker_genes=frozenset({"k1", "k2"}),
    )
    assert boundary.start == 8000
    assert boundary.left_stop == 7900


def test_one_viral_match_and_one_foreign_supporter_do_not_bridge_a_host_run():
    boundary = extend(
        [
            gene("h1", 2100, 2200, host=True),
            gene("h2", 2300, 2400, host=True),
            gene("h3", 2500, 2600, host=True),
            gene("v1", 2700, 2800, viral=True),
            foreign_gene("f1", 2900, 3000),
        ]
    )
    assert boundary.end == 2000
    assert boundary.right_stop == 2100


def test_one_viral_match_and_one_foreign_supporter_do_not_bridge_a_host_run_on_the_left():
    boundary = extend(
        [
            gene("h1", 7800, 7900, host=True),
            gene("h2", 7600, 7700, host=True),
            gene("h3", 7400, 7500, host=True),
            gene("v1", 7200, 7300, viral=True),
            foreign_gene("f1", 7000, 7100),
        ],
        start=8000,
        end=9000,
    )
    assert boundary.start == 8000
    assert boundary.left_stop == 7900


def test_two_validated_viral_markers_without_a_top10_hit_bridge_a_host_run():
    boundary = extend(
        [
            gene("h1", 2100, 2200, host=True),
            gene("h2", 2300, 2400, host=True),
            gene("h3", 2500, 2600, host=True),
            marked_marker_gene("k1", 2700, 2800),
            marked_marker_gene("k2", 2900, 3000),
        ]
    )
    assert boundary.end == 3000
    assert boundary.right_stop is None


def test_second_supporting_gene_at_the_fifth_entry_bridges_the_host_run():
    boundary = extend(
        [
            gene("h1", 2100, 2200, host=True),
            gene("h2", 2300, 2400, host=True),
            gene("h3", 2500, 2600, host=True),
            gene("s1", 2700, 2800, viral=True),
            gene("n1", 2900, 3000),
            gene("n2", 3100, 3200),
            gene("n3", 3300, 3400),
            gene("s2", 3500, 3600, viral=True),
        ]
    )
    assert boundary.end == 3600
    assert boundary.right_stop is None


def test_second_supporting_gene_at_the_sixth_entry_does_not_bridge_the_host_run():
    boundary = extend(
        [
            gene("h1", 2100, 2200, host=True),
            gene("h2", 2300, 2400, host=True),
            gene("h3", 2500, 2600, host=True),
            gene("s1", 2700, 2800, viral=True),
            gene("n1", 2900, 3000),
            gene("n2", 3100, 3200),
            gene("n3", 3300, 3400),
            gene("n4", 3500, 3600),
            gene("s2", 3700, 3800, viral=True),
        ]
    )
    assert boundary.end == 2000
    assert boundary.right_stop == 2100


def test_bridged_host_genes_count_toward_the_five_gene_support_gap():
    boundary = extend(
        [
            gene("h1", 2100, 2200, host=True),
            gene("h2", 2300, 2400, host=True),
            gene("h3", 2500, 2600, host=True),
            gene("n1", 2700, 2800),
            gene("n2", 2900, 3000),
            gene("s1", 3100, 3200, viral=True),
            gene("s2", 3300, 3400, viral=True),
        ]
    )
    assert boundary.end == 2000
    assert boundary.right_stop == 2100
    assert boundary.right_stop_reason == "support_gap"


def test_supporters_beyond_a_gene_without_a_record_do_not_bridge_a_host_run():
    genes = [
        gene("h1", 2100, 2200, host=True),
        gene("h2", 2300, 2400, host=True),
        gene("h3", 2500, 2600, host=True),
        gene("u", 2700, 2800),
        gene("s1", 2900, 3000, viral=True),
        gene("s2", 3100, 3200, viral=True),
    ]
    boundary = extend(genes, taxonomy={item.id: item for item in genes if item.id != "u"})
    assert boundary.end == 2000
    assert boundary.right_stop == 2100
    assert boundary.right_stop_reason == "host_stop"


def test_supporters_beyond_a_gene_without_a_record_do_not_bridge_a_host_run_on_the_left():
    genes = [
        gene("h1", 7800, 7900, host=True),
        gene("h2", 7600, 7700, host=True),
        gene("h3", 7400, 7500, host=True),
        gene("u", 7200, 7300),
        gene("s1", 7000, 7100, viral=True),
        gene("s2", 6800, 6900, viral=True),
    ]
    boundary = extend(genes, taxonomy={item.id: item for item in genes if item.id != "u"}, start=8000, end=9000)
    assert boundary.start == 8000
    assert boundary.left_stop == 7900


def test_supporters_before_a_gene_without_a_record_still_bridge_a_host_run():
    genes = [
        gene("h1", 2100, 2200, host=True),
        gene("h2", 2300, 2400, host=True),
        gene("h3", 2500, 2600, host=True),
        gene("s1", 2700, 2800, viral=True),
        gene("s2", 2900, 3000, viral=True),
        gene("u", 3100, 3200),
    ]
    boundary = extend(genes, taxonomy={item.id: item for item in genes if item.id != "u"})
    assert boundary.end == 3000
    assert boundary.right_stop is None


def test_tied_host_stop_and_missing_gene_stop_before_the_gap_on_the_right():
    genes = [
        gene("v", 2100, 2600, viral=True),
        gene("u", 2300, 2400),
        gene("h1", 2300, 2800, host=True),
        gene("h2", 2900, 3000, host=True),
        gene("h3", 3100, 3200, host=True),
    ]
    boundary = extend(genes, taxonomy={item.id: item for item in genes if item.id != "u"})
    assert boundary.end == 2000
    assert boundary.right_stop == 2300


def test_tied_host_stop_and_missing_gene_stop_before_the_gap_on_the_left():
    genes = [
        gene("v", 7400, 7900, viral=True),
        gene("u", 7600, 7700),
        gene("h1", 7200, 7700, host=True),
        gene("h2", 7000, 7100, host=True),
        gene("h3", 6800, 6900, host=True),
    ]
    boundary = extend(genes, taxonomy={item.id: item for item in genes if item.id != "u"}, start=8000, end=9000)
    assert boundary.start == 8000
    assert boundary.left_stop == 7700


def test_tied_host_stop_and_missing_gene_reject_a_gene_ending_inside_both():
    genes = [
        gene("v", 2100, 2350, viral=True),
        gene("u", 2300, 2400),
        gene("h1", 2300, 2800, host=True),
        gene("h2", 2900, 3000, host=True),
        gene("h3", 3100, 3200, host=True),
    ]
    boundary = extend(genes, taxonomy={item.id: item for item in genes if item.id != "u"})
    assert boundary.end == 2000


def test_tied_host_stop_and_missing_gene_keep_the_gap_edge_before_host_end():
    genes = [
        gene("v", 2100, 2600, viral=True),
        gene("u", 2300, 2800),
        gene("h1", 2300, 2400, host=True),
        gene("h2", 2900, 3000, host=True),
        gene("h3", 3100, 3200, host=True),
    ]
    boundary = extend(genes, taxonomy={item.id: item for item in genes if item.id != "u"})
    assert boundary.end == 2000


def test_disabled_host_stop_keeps_the_span_of_the_walk_with_the_host_stop():
    genes = [
        gene("v", 2100, 2700, viral=True),
        gene("h1", 2300, 2800, host=True),
        gene("h2", 2350, 2400, host=True),
        gene("h3", 2400, 2450, host=True),
    ]
    with_stop = extend(genes, start=1300, end=2000, envelope=(1000, 2500))
    without_stop = extend(genes, start=1300, end=2000, envelope=(1000, 2500), disable_host_stop=True)
    assert with_stop.end == 2700
    assert without_stop.end == 2700
    assert without_stop.right_stop is None


def test_uncrossed_bridged_host_run_sets_the_stop_on_the_left():
    boundary = extend(
        [
            gene("h1", 7800, 7900, host=True),
            gene("h2", 7600, 7700, host=True),
            gene("h3", 7400, 7500, host=True),
            gene("n1", 7200, 7300),
            gene("n2", 7000, 7100),
            gene("s1", 6800, 6900, viral=True),
            gene("s2", 6600, 6700, viral=True),
        ],
        start=8000,
        end=9000,
    )
    assert boundary.start == 8000
    assert boundary.left_stop == 7900
    assert boundary.left_stop_reason == "support_gap"


def test_stop_of_an_uncrossed_bridged_run_is_the_accepted_end_when_that_is_further_out():
    boundary = extend(
        [
            gene("s0", 2100, 2450, viral=True),
            gene("h1", 2300, 2400, host=True),
            gene("h2", 2500, 2600, host=True),
            gene("h3", 2700, 2800, host=True),
            gene("n1", 2900, 3000),
            gene("n2", 3100, 3200),
            gene("s1", 3300, 3400, viral=True),
            gene("s2", 3500, 3600, viral=True),
        ]
    )
    assert boundary.end == 2450
    assert boundary.right_stop == 2450


def test_crossed_first_run_and_uncrossed_second_run_stop_at_the_second_on_the_right():
    boundary = extend(
        [
            gene("h1", 2100, 2200, host=True),
            gene("h2", 2300, 2400, host=True),
            gene("h3", 2500, 2600, host=True),
            gene("s1", 2700, 2800, viral=True),
            gene("s2", 2900, 3000, viral=True),
            gene("h4", 3100, 3200, host=True),
            gene("h5", 3300, 3400, host=True),
            gene("h6", 3500, 3600, host=True),
        ]
    )
    assert boundary.end == 3000
    assert boundary.right_stop == 3100
    assert boundary.right_stop_reason == "host_stop"


def test_crossed_first_run_and_uncrossed_second_run_stop_at_the_second_on_the_left():
    boundary = extend(
        [
            gene("h1", 7800, 7900, host=True),
            gene("h2", 7600, 7700, host=True),
            gene("h3", 7400, 7500, host=True),
            gene("s1", 7200, 7300, viral=True),
            gene("s2", 7000, 7100, viral=True),
            gene("h4", 6800, 6900, host=True),
            gene("h5", 6600, 6700, host=True),
            gene("h6", 6400, 6500, host=True),
        ],
        start=8000,
        end=9000,
    )
    assert boundary.start == 7000
    assert boundary.left_stop == 6900
    assert boundary.left_stop_reason == "host_stop"


def test_disabled_host_stop_returns_no_stop_for_an_uncrossed_bridged_run():
    boundary = extend(
        [
            gene("h1", 2100, 2200, host=True),
            gene("h2", 2300, 2400, host=True),
            gene("h3", 2500, 2600, host=True),
            gene("n1", 2700, 2800),
            gene("n2", 2900, 3000),
            gene("s1", 3100, 3200, viral=True),
            gene("s2", 3300, 3400, viral=True),
        ],
        disable_host_stop=True,
    )
    assert boundary.right_stop is None


def test_unbridged_run_that_limits_the_walk_sets_a_stop_when_hosts_start_before_the_limit():
    boundary = extend(
        [
            gene("h1", 1700, 3000, host=True),
            gene("h2", 1800, 2100, host=True),
            gene("h3", 1850, 2150, host=True),
            gene("s", 1900, 2200, viral=True),
            gene("h4", 2300, 2400, host=True),
            gene("h5", 2500, 2600, host=True),
            gene("h6", 2700, 2800, host=True),
        ]
    )
    assert boundary.end == 2200
    assert boundary.right_stop == 2200


def test_unbridged_run_that_limits_the_walk_sets_a_stop_when_hosts_start_before_the_limit_on_the_left():
    boundary = extend(
        [
            gene("h1", 7000, 8300, host=True),
            gene("h2", 7900, 8200, host=True),
            gene("h3", 7850, 8150, host=True),
            gene("s", 7800, 8100, viral=True),
            gene("h4", 7600, 7700, host=True),
            gene("h5", 7400, 7500, host=True),
            gene("h6", 7200, 7300, host=True),
        ],
        start=8000,
        end=9000,
    )
    assert boundary.start == 7800
    assert boundary.left_stop == 7800


def test_disabled_host_stop_does_not_take_the_reach_of_the_removed_host_gene():
    genes = [
        gene("h1", 2100, 4000, host=True),
        gene("h2", 2200, 2250, host=True),
        gene("h3", 2300, 2350, host=True),
        gene("s", 2400, 3500, viral=True),
        gene("u", 2500, 2600),
    ]
    boundary = extend(genes, taxonomy={item.id: item for item in genes if item.id != "u"}, disable_host_stop=True)
    assert boundary.end == 2000


def test_disabled_host_stop_rejects_gap_overlap_after_a_short_first_host():
    genes = [
        gene("h1", 2100, 2150, host=True),
        gene("h2", 2200, 2250, host=True),
        gene("h3", 2300, 2350, host=True),
        gene("s", 2400, 3000, viral=True),
        gene("u", 2500, 3100),
    ]
    boundary = extend(genes, taxonomy={item.id: item for item in genes if item.id != "u"}, disable_host_stop=True)
    assert boundary.end == 2000


def test_walk_limit_stops_before_the_first_gene_without_a_record():
    from virosync.pipeline.phase2.seed_extension import seed_walk_limit

    genes = [gene("a", 2100, 2200), gene("b", 2300, 2400), gene("c", 2500, 2600), gene("d", 2700, 2800)]
    records = {item.id: item for item in genes if item.id != "c"}
    assert seed_walk_limit((1000, 2000), genes, records) == (1000, 2400)


def test_walk_limit_extends_to_the_left_over_genes_with_records():
    from virosync.pipeline.phase2.seed_extension import seed_walk_limit

    genes = [gene("a", 800, 900), gene("b", 600, 700), gene("c", 400, 500), gene("d", 200, 300)]
    records = {item.id: item for item in genes if item.id != "c"}
    assert seed_walk_limit((1000, 2000), genes, records) == (600, 2000)


def test_walk_limit_equals_the_envelope_when_the_first_outer_gene_has_no_record():
    from virosync.pipeline.phase2.seed_extension import seed_walk_limit

    genes = [gene("a", 2100, 2200), gene("b", 2300, 2400)]
    assert seed_walk_limit((1000, 2000), genes, {"b": genes[1]}) == (1000, 2000)


def test_walk_limit_covers_a_gene_with_a_record_that_crosses_the_envelope_edge():
    from virosync.pipeline.phase2.seed_extension import seed_walk_limit

    genes = [gene("a", 1900, 2300), gene("b", 2400, 2500)]
    assert seed_walk_limit((1000, 2000), genes, {"a": genes[0]}) == (1000, 2300)


def test_walk_limit_ignores_genes_inside_the_envelope():
    from virosync.pipeline.phase2.seed_extension import seed_walk_limit

    genes = [gene("a", 1100, 1200), gene("b", 2100, 2200)]
    assert seed_walk_limit((1000, 2000), genes, {"b": genes[1]}) == (1000, 2200)


def test_genes_outside_the_envelope_are_ignored_without_a_walk_limit():
    boundary = extend([gene("out", 2700, 2800, viral=True)], envelope=(0, 2500))
    assert boundary.end == 2000


def test_viral_match_outside_the_envelope_supports_growth():
    boundary = extend([gene("out", 2700, 2800, viral=True)], envelope=(0, 2500), walk_limit=(0, 4000))
    assert boundary.end == 2800


def test_host_like_gene_with_a_viral_match_outside_the_envelope_supports_growth():
    boundary = extend([gene("out", 2700, 2800, host=True, viral=True)], envelope=(0, 2500), walk_limit=(0, 4000))
    assert boundary.end == 2800


def test_walk_to_the_limit_reports_the_search_envelope():
    boundary = extend([gene("out", 2700, 2800, viral=True)], envelope=(0, 2500), walk_limit=(0, 2800))
    assert boundary.right_stop_reason == "search_envelope"


def test_foreign_gene_inside_the_envelope_supports_growth():
    boundary = extend([foreign_gene("in", 2100, 2200)], envelope=(0, 2500), walk_limit=(0, 4000))
    assert boundary.end == 2200


def test_foreign_gene_outside_the_envelope_does_not_support_growth():
    boundary = extend([foreign_gene("out", 2700, 2800)], envelope=(0, 2500), walk_limit=(0, 4000))
    assert boundary.end == 2000


def test_marker_set_membership_outside_the_envelope_does_not_support_growth():
    boundary = extend(
        [gene("out", 2700, 2800)],
        envelope=(0, 2500),
        walk_limit=(0, 4000),
        viral_marker_genes=frozenset({"out"}),
    )
    assert boundary.end == 2000


def test_ml_support_outside_the_envelope_does_not_support_growth():
    boundary = extend(
        [gene("out", 2700, 2800)], envelope=(0, 2500), walk_limit=(0, 4000), ml_support=frozenset({"out"})
    )
    assert boundary.end == 2000


def test_ml_support_inside_the_envelope_supports_growth():
    boundary = extend([gene("in", 2100, 2200)], envelope=(0, 2500), walk_limit=(0, 4000), ml_support=frozenset({"in"}))
    assert boundary.end == 2200


def test_searched_genes_without_hits_outside_the_envelope_count_toward_the_gap():
    boundary = extend(
        [
            no_hit_gene("n1", 2700, 2800),
            no_hit_gene("n2", 2900, 3000),
            no_hit_gene("n3", 3100, 3200),
            no_hit_gene("n4", 3300, 3400),
            no_hit_gene("n5", 3500, 3600),
            gene("far", 3700, 3800, viral=True),
        ],
        envelope=(0, 2500),
        walk_limit=(0, 4000),
    )
    assert boundary.end == 2000
    assert boundary.right_stop_reason == "support_gap"


def test_searched_genes_without_hits_inside_the_envelope_are_skipped():
    boundary = extend(
        [
            no_hit_gene("n1", 2100, 2200),
            no_hit_gene("n2", 2300, 2400),
            no_hit_gene("n3", 2500, 2600),
            no_hit_gene("n4", 2700, 2800),
            no_hit_gene("n5", 2900, 3000),
            gene("far", 3100, 3200, viral=True),
        ],
        envelope=(0, 4000),
        walk_limit=(0, 4000),
    )
    assert boundary.end == 3200


def test_searched_genes_without_hits_inside_the_envelope_do_not_use_up_the_bridge_window():
    boundary = extend(
        [
            gene("h1", 2100, 2200, host=True),
            gene("h2", 2300, 2400, host=True),
            gene("h3", 2500, 2600, host=True),
            no_hit_gene("n1", 2700, 2800),
            no_hit_gene("n2", 2900, 3000),
            no_hit_gene("n3", 3100, 3200),
            no_hit_gene("n4", 3300, 3400),
            no_hit_gene("n5", 3500, 3600),
            gene("s1", 3700, 3800, viral=True),
            gene("s2", 3900, 4000, viral=True),
        ],
        envelope=(0, 5000),
        walk_limit=(0, 5000),
    )
    assert boundary.end == 4000
    assert boundary.right_stop is None


def test_searched_genes_without_hits_outside_the_envelope_use_up_the_bridge_window():
    boundary = extend(
        [
            gene("h1", 2100, 2200, host=True),
            gene("h2", 2300, 2400, host=True),
            gene("h3", 2500, 2600, host=True),
            no_hit_gene("n1", 2700, 2800),
            no_hit_gene("n2", 2900, 3000),
            no_hit_gene("n3", 3100, 3200),
            no_hit_gene("n4", 3300, 3400),
            no_hit_gene("n5", 3500, 3600),
            gene("s1", 3700, 3800, viral=True),
            gene("s2", 3900, 4000, viral=True),
        ],
        envelope=(0, 2050),
        walk_limit=(0, 5000),
    )
    assert boundary.end == 2000
    assert boundary.right_stop == 2100


def test_three_host_genes_outside_the_envelope_stop_the_walk():
    boundary = extend(
        [
            gene("near", 2100, 2200, viral=True),
            gene("h1", 2700, 2800, host=True),
            gene("h2", 2900, 3000, host=True),
            gene("h3", 3100, 3200, host=True),
            gene("far", 3300, 3400, viral=True),
        ],
        envelope=(0, 2500),
        walk_limit=(0, 4000),
    )
    assert boundary.end == 2200
    assert boundary.right_stop == 2700


def test_viral_gene_outside_the_envelope_resets_the_host_run():
    boundary = extend(
        [
            gene("h1", 2700, 2800, host=True),
            gene("h2", 2900, 3000, host=True),
            gene("v", 3100, 3200, viral=True),
            gene("h3", 3300, 3400, host=True),
        ],
        envelope=(0, 2500),
        walk_limit=(0, 4000),
    )
    assert boundary.right_stop is None


def test_gene_without_a_record_outside_the_envelope_leaves_the_host_stop_available():
    boundary = extend([gene("out", 2700, 2800)], taxonomy={}, envelope=(0, 2500), walk_limit=(0, 4000))
    assert boundary.host_stop_available is True
    assert boundary.right_stop_reason == "taxonomy_unavailable"


def test_gene_without_a_record_inside_the_envelope_removes_the_host_stop():
    boundary = extend([gene("in", 2100, 2200)], taxonomy={}, envelope=(0, 2500), walk_limit=(0, 4000))
    assert boundary.host_stop_available is False


def test_phase2_rejects_support_crossing_an_unsearched_gene_at_the_envelope_edge(tmp_path, monkeypatch):
    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome import phase2

    arguments = _stub_phase2_boundary_search(
        tmp_path,
        monkeypatch,
        mapped=True,
        genes=[gene("v", 1500, 2050, viral=True), gene("u", 1950, 2100)],
        unsearched=frozenset({"u"}),
    )
    config = PipelineConfig().with_overrides(
        resume=False, gene_taxonomy_faa_db=tmp_path / "db", boundary_taxonomy_ml_enabled=False, threads=1
    )
    seed = MergedSeed(
        "s",
        1000,
        2000,
        seed_id="original",
        predicted_family="CRESS",
        confirmed_seeds=(ConfirmedSeed("original", "s", 1300, 1700),),
    )

    result = phase2._run_phase2_subflow(merged_seeds=[seed], config=config, **arguments)

    assert (result.refined_boundaries[0].start, result.refined_boundaries[0].end) == (1300, 1700)


def test_phase2_walks_past_the_searched_envelope_over_genes_with_records(tmp_path, monkeypatch):
    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome import phase2

    arguments = _stub_phase2_boundary_search(
        tmp_path,
        monkeypatch,
        mapped=True,
        flank_end_bp=3000,
        genes=[
            gene("v1", 2500, 2600, viral=True),
            gene("v2", 3500, 3600, viral=True),
            gene("v3", 4500, 4600, viral=True),
        ],
    )
    config = PipelineConfig().with_overrides(
        resume=False, gene_taxonomy_faa_db=tmp_path / "db", boundary_taxonomy_ml_enabled=False, threads=1
    )
    seed = MergedSeed(
        "s",
        1000,
        2000,
        seed_id="original",
        predicted_family="NCLDV",
        confirmed_seeds=(ConfirmedSeed("original", "s", 1300, 1700),),
    )

    result = phase2._run_phase2_subflow(merged_seeds=[seed], config=config, **arguments)

    assert (result.refined_boundaries[0].start, result.refined_boundaries[0].end) == (1300, 4600)


@pytest.mark.parametrize("direction,expected_span,expected_stop", [(1, (1000, 2000), 2100), (-1, (8000, 9000), 7900)])
def test_missing_record_inside_host_run_blocks_later_bridge(direction, expected_span, expected_stop):
    genes = [
        gene("h1", 2100, 2200, host=True),
        gene("f", 2250, 2280),
        gene("h2", 2300, 2400, host=True),
        gene("u", 2500, 2600),
        gene("h3", 2700, 2800, host=True),
        gene("s1", 2900, 3000, viral=True),
        gene("s2", 3100, 3200, viral=True),
    ]
    genes[1].share = 0.0
    if direction == -1:
        for item in genes:
            item.start, item.end = 10000 - item.end, 10000 - item.start
    boundary = extend(
        genes,
        start=expected_span[0],
        end=expected_span[1],
        taxonomy={item.id: item for item in genes if item.id != "u"},
    )
    assert (boundary.start, boundary.end) == expected_span
    assert (boundary.right_stop if direction == 1 else boundary.left_stop) == expected_stop


@pytest.mark.parametrize("direction,core", [(1, (1000, 2000)), (-1, (8000, 9000))])
@pytest.mark.parametrize("disable_host_stop", [False, True])
@pytest.mark.parametrize(
    "support_end,expected_right,expected_left", [(2900, 2000, 8000), (2600, 2000, 8000), (2500, 2500, 7500)]
)
def test_nested_missing_record_caps_whole_gene_acceptance(
    direction, core, disable_host_stop, support_end, expected_right, expected_left
):
    genes = [
        gene("v", 2050, support_end, viral=True),
        gene("h1", 2100, 3000, host=True),
        gene("h2", 2150, 2200, host=True),
        gene("h3", 2250, 2300, host=True),
        gene("u", 2500, 2600),
    ]
    if direction == -1:
        for item in genes:
            item.start, item.end = 10000 - item.end, 10000 - item.start
    boundary = extend(
        genes,
        start=core[0],
        end=core[1],
        taxonomy={item.id: item for item in genes if item.id != "u"},
        disable_host_stop=disable_host_stop,
    )
    assert (boundary.start, boundary.end) == ((1000, expected_right) if direction == 1 else (expected_left, 9000))


@pytest.mark.parametrize("direction,core", [(1, (1000, 2000)), (-1, (8000, 9000))])
@pytest.mark.parametrize("disable_host_stop,expected_right,expected_left", [(False, 2000, 8000), (True, 2000, 8000)])
def test_missing_record_before_nested_host_caps_both_walks(
    direction, core, disable_host_stop, expected_right, expected_left
):
    genes = [
        gene("v", 2050, 2900, viral=True),
        gene("u", 2100, 3000),
        gene("h1", 2500, 2600, host=True),
        gene("h2", 2550, 2600, host=True),
        gene("h3", 2570, 2650, host=True),
    ]
    if direction == -1:
        for item in genes:
            item.start, item.end = 10000 - item.end, 10000 - item.start
    boundary = extend(
        genes,
        start=core[0],
        end=core[1],
        taxonomy={item.id: item for item in genes if item.id != "u"},
        disable_host_stop=disable_host_stop,
    )
    assert (boundary.start, boundary.end) == ((1000, expected_right) if direction == 1 else (expected_left, 9000))


def test_phase2_asymmetric_bridge_reports_overlap_without_changing_seed_ownership(tmp_path, monkeypatch):
    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome import phase2
    from virosync.pipeline.phase2.boundary_diamond import GenomeDiamondQuery, SeedGeneMapping

    genes = [
        gene("a", 1800, 1900, viral=True),
        gene("h1", 2100, 2200, host=True),
        gene("h2", 2300, 2400, host=True),
        gene("h3", 2500, 2600, host=True),
        gene("b1", 3500, 3600, viral=True),
        gene("b2", 3700, 3800, viral=True),
    ]
    seeds = [
        MergedSeed(
            "s",
            1000,
            2000,
            seed_id="a",
            predicted_family="NCLDV",
            confirmed_seeds=(ConfirmedSeed("a", "s", 1000, 2000),),
        ),
        MergedSeed(
            "s",
            3500,
            3800,
            seed_id="b",
            predicted_family="NCLDV",
            confirmed_seeds=(ConfirmedSeed("b", "s", 3500, 3800),),
        ),
    ]
    arguments = _stub_phase2_boundary_search(tmp_path, monkeypatch, mapped=True, genes=genes)
    query = GenomeDiamondQuery(
        seed_gene_mappings={
            "a": SeedGeneMapping("a", "s", 1000, 2000, flank_start_bp=0, flank_end_bp=2050),
            "b": SeedGeneMapping("b", "s", 3500, 3800, flank_start_bp=3400, flank_end_bp=4000),
        }
    )
    monkeypatch.setattr(
        phase2,
        "_run_boundary_diamond",
        lambda **kwargs: phase2._BoundaryEvidence({item.id: item for item in genes}, None, query),
    )
    monkeypatch.setattr(
        "virosync.pipeline.phase2.seed_extension.host_lineage_share", lambda record, model: 1.0 if record.host else 0.2
    )
    arguments["host_signature_model"] = SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0})
    config = PipelineConfig().with_overrides(
        resume=False, gene_taxonomy_faa_db=tmp_path / "db", boundary_taxonomy_ml_enabled=False, threads=1
    )
    result = phase2._run_phase2_subflow(merged_seeds=seeds, config=config, **arguments)

    assert len(result.refined_boundaries) == 2
    first, second = result.refined_boundaries
    assert (first.start, first.end, second.start, second.end) == (1000, 3800, 3500, 3800)
    assert first.confirmed_seeds == seeds[0].confirmed_seeds
    assert second.confirmed_seeds == seeds[1].confirmed_seeds
    assert first.right_stop is None
    assert second.left_stop == 2600
    assert max(first.start, second.start) < min(first.end, second.end)


@pytest.mark.parametrize("direction", [1, -1])
@pytest.mark.parametrize("host_before_gap", [False, True])
def test_four_base_gap_overlap_rejects_whole_support_in_normal_and_a4(direction, host_before_gap):
    genes = [gene("safe", 2010, 2040, viral=True), gene("support", 2050, 2504, viral=True), gene("missing", 2500, 3000)]
    if host_before_gap:
        genes += [
            gene("h1", 2100, 3100, host=True),
            gene("h2", 2150, 2200, host=True),
            gene("h3", 2250, 2300, host=True),
        ]
    if direction == -1:
        for item in genes:
            item.start, item.end = 10000 - item.end, 10000 - item.start
    core = (1000, 2000) if direction == 1 else (8000, 9000)
    taxonomy = {item.id: item for item in genes if item.id != "missing"}
    normal = extend(genes, start=core[0], end=core[1], taxonomy=taxonomy)
    a4 = extend(genes, start=core[0], end=core[1], taxonomy=taxonomy, disable_host_stop=True)
    expected = (1000, 2040) if direction == 1 else (7960, 9000)
    assert (normal.start, normal.end) == expected
    assert (a4.start, a4.end) == expected
    assert a4.start <= normal.start and a4.end >= normal.end
    assert normal.confirmed_seeds == a4.confirmed_seeds == (ConfirmedSeed("original", "s", *core),)


@pytest.mark.parametrize("direction", [1, -1])
def test_missing_gene_straddling_core_preserves_seed_and_phase3_refusal(direction):
    from virosync.orchestration._flows.single_genome.phase3 import _query_boundary_coordinate_records

    genes = [gene("missing", 1900, 2100), gene("support", 1950, 2104, viral=True)]
    if direction == -1:
        for item in genes:
            item.start, item.end = 10000 - item.end, 10000 - item.start
    core = (1000, 2000) if direction == 1 else (8000, 9000)
    taxonomy = {"support": genes[1]}
    boundary = extend(genes, start=core[0], end=core[1], taxonomy=taxonomy)
    assert (boundary.start, boundary.end) == core
    assert boundary.confirmed_seeds == (ConfirmedSeed("original", "s", *core),)
    with pytest.raises(RuntimeError, match="1 overlapping genes were not searched"):
        _query_boundary_coordinate_records(
            boundary=boundary, taxonomy_index={}, marker_index={}, taxonomy_map=taxonomy, proteome_index={"s": genes}
        )
