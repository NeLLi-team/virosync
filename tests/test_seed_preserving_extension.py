"""Behavioral contract for seed-preserving outward growth."""

from types import SimpleNamespace
from unittest.mock import patch

from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
from virosync.pipeline.phase2.seed_extension import extend_confirmed_seed


def gene(name, start, end, *, viral=False, host=False):
    return SimpleNamespace(
        id=name,
        porf_id=name,
        start=start,
        end=end,
        top10_prefixes=["NCLDV__"] if viral else ["EUK__"],
        top10_pidents=[30.0],
        host=host,
    )


def grow(genes, *, start=1000, end=2000):
    seed = MergedSeed(
        "s", start, end, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", start, end),)
    )
    taxonomy = {g.id: g for g in genes}
    model = SimpleNamespace(token_weights={"host": 1})
    with patch("virosync.pipeline.phase2.seed_extension.host_lineage_share", side_effect=lambda g, m: float(g.host)):
        return extend_confirmed_seed(seed, genes, taxonomy, (0, 10000), model, 0.35)


def test_gene_free_padding_and_original_identity_survive():
    boundary = grow([gene("inside", 1300, 1700, host=True)])
    assert (boundary.start, boundary.end) == (1000, 2000)
    assert boundary.confirmed_seeds[0].seed_id == "original"
    assert (boundary.confirmed_seeds[0].start, boundary.confirmed_seeds[0].end) == (1000, 2000)


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
    assert boundary.end == 2100
    assert boundary.right_stop_reason == "host_stop"


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


def test_missing_taxonomy_stops_growth_but_searched_no_hit_is_pending():
    seed = MergedSeed(
        "s", 1000, 2000, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", 1000, 2000),)
    )
    unknown = gene("u", 2100, 2200)
    positive = gene("v", 2300, 2400, viral=True)
    missing = extend_confirmed_seed(seed, [unknown, positive], {"v": positive}, (0, 10000), None, 0.35)
    searched = extend_confirmed_seed(seed, [unknown, positive], {"u": unknown, "v": positive}, (0, 10000), None, 0.35)
    assert missing.end == 2000
    assert missing.right_stop_reason == "taxonomy_unavailable"
    assert searched.end == 2400


def test_repeat_partition_cannot_cut_core_or_create_unseeded_residual():
    from dataclasses import replace

    from virosync.pipeline.phase2.seed_extension import accept_seed_partition

    parent = grow([gene("v", 2100, 2500, viral=True)])
    assert accept_seed_partition(parent, [replace(parent, start=1100)]) is None
    assert accept_seed_partition(parent, [replace(parent, end=2000), replace(parent, start=2000)]) is None
    assert accept_seed_partition(parent, [replace(parent, end=2200)])[0].confirmed_seeds == parent.confirmed_seeds


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
    model = SimpleNamespace(token_weights={"host": 1})
    with patch("virosync.pipeline.phase2.seed_extension.host_lineage_share", side_effect=lambda g, m: float(g.host)):
        boundary = extend_confirmed_seed(
            seed, genes, {g.id: g for g in genes}, (0, 10000), model, 0.35, disable_host_stop=True
        )
    assert (boundary.start, boundary.end) == (1000, 2800)
    assert boundary.right_stop is None


def test_phase2_ml_failure_and_off_keep_complete_original_seed(tmp_path, monkeypatch):
    import logging

    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome import phase2
    from virosync.pipeline.phase2 import taxonomy_ml_refiner
    from virosync.pipeline.phase2.boundary_diamond import GenomeDiamondQuery, SeedGeneMapping

    genes = [gene("v", 2300, 2500, viral=True)]
    query = GenomeDiamondQuery(
        seed_gene_mappings={
            "original": SeedGeneMapping("original", "s", 1000, 2000, flank_start_bp=0, flank_end_bp=10000)
        }
    )
    monkeypatch.setattr(phase2, "build_proteome_index", lambda _: {"s": genes})
    monkeypatch.setattr(phase2, "_require_phase2b_gene_taxonomy_db", lambda *args, **kwargs: tmp_path / "db")
    monkeypatch.setattr(
        phase2, "_run_boundary_diamond", lambda **kwargs: phase2._BoundaryEvidence({"v": genes[0]}, None, query)
    )
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
                "s", 1000, 2000, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", 1000, 2000),)
            )
        ],
        config=config,
        **arguments,
    )
    failed = phase2._run_phase2_subflow(
        merged_seeds=[
            MergedSeed(
                "s", 1000, 2000, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", 1000, 2000),)
            )
        ],
        config=config.with_overrides(boundary_taxonomy_ml_enabled=True),
        **arguments,
    )
    assert off.refined_boundaries == failed.refined_boundaries
    assert (off.refined_boundaries[0].start, off.refined_boundaries[0].end) == (1000, 2500)
    assert off.refined_boundaries[0].confirmed_seeds[0].end == 2000


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
    with patch("virosync.pipeline.phase2.seed_extension.host_lineage_share", side_effect=lambda g, m: float(g.host)):
        boundary = extend_confirmed_seed(
            seed,
            genes,
            {g.id: g for g in genes},
            (0, 10000),
            SimpleNamespace(token_weights={"h": 1}),
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
    assert result.start == 900
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
        side_effect=lambda record, model: float(record.host),
    ):
        left = extend_confirmed_seed(
            left_seed, genes, taxonomy, (0, 6000), SimpleNamespace(token_weights={"host": 1}), 0.35
        )
        right = extend_confirmed_seed(
            right_seed, genes, taxonomy, (0, 6000), SimpleNamespace(token_weights={"host": 1}), 0.35
        )
    outputs = merge_reachable_boundaries([left, right])
    assert [(item.start, item.end) for item in outputs] == [(1000, 3200), (2100, 5000)]
    assert outputs[0].end <= outputs[0].right_stop == 3300
    assert outputs[0].confirmed_seeds == left_seed.confirmed_seeds
    assert outputs[1].confirmed_seeds == right_seed.confirmed_seeds
