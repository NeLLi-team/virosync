"""Family association does not manufacture roles or count model alternatives."""

from types import SimpleNamespace

import pytest

from virosync.pipeline.phase1 import viral_markers as families
from virosync.pipeline.phase1.seed_merger import classify_marker_hits
from virosync.pipeline.phase3.evidence_synthesizer import (
    load_marker_annotation_index,
    summarize_marker_hits,
)


def hit(gene, model, score=100, **fields):
    """Make an ordinary saved-marker-shaped fixture."""
    return SimpleNamespace(
        query_porf=gene,
        hmm_target=model,
        hmm_score=score,
        validation_status=fields.pop("validation_status", "validated"),
        pfam_override=fields.pop("pfam_override", ""),
        **fields,
    )


def test_guarded_family_map_excludes_cellular_credit():
    index = {"vs000084": {"family": "NCLDV"}, "cog0504": {"family": "NCLDV"}, "vs000492": {"family": "NCLDV"}}
    assert families.metadata_marker_families(index) == {"vs000492": "NCLDV"}


def test_best_ordinary_hit_then_distinct_models_and_allowlist():
    hits = [
        hit("one", "VS000492", 100),
        hit("one", "GVOGm1", 200, pfam_override="PF1"),
        hit("two", "VS000492", 100),
        hit("three", "other", 200),
        hit("three", "VS000492", 100),
    ]
    summary = classify_marker_hits(hits, metadata_families={"vs000492": "NCLDV"})
    assert summary["ncldv_markers"] == 1
    assert classify_marker_hits(hits[3:], ["VS000492"], metadata_families={"vs000492": "NCLDV"})["ncldv_markers"] == 0


def test_blocked_winner_does_not_choose_lower_alternative_for_credit():
    hits = [hit("one", "VS000084", 200), hit("one", "VS000492", 100)]
    assert classify_marker_hits(hits, metadata_families={"vs000492": "NCLDV"})["ncldv_markers"] == 0


def test_first_seen_tie_is_preserved():
    hits = [hit("one", "unknown", 100), hit("one", "VS000492", 100)]
    assert classify_marker_hits(hits, metadata_families={"vs000492": "NCLDV"})["ncldv_markers"] == 0


def test_named_family_precedes_conflicting_metadata_in_both_consumers():
    summary = classify_marker_hits([hit("one", "GVOGm1")], metadata_families={"gvogm1": "MIRUS"})
    assert summary["classification"] == "NCLDV"
    assert summarize_marker_hits(["GVOGm1"], {"gvogm1": {"family": "MIRUS"}})["families"] == ["NCLDV"]


@pytest.mark.parametrize("model", ["EVE_CAPSID_NOVEL_0001", "EVE_CAPSID_NOVEL_0002", "EVE_CAPSID_NOVEL_0033"])
def test_structure_metadata_never_supplies_training_class(tmp_path, model):
    path = tmp_path / "annotations.tsv"
    path.write_text(
        "model_name\tsource\tdescription\n" + model + "\tEVE_structure_20260811\tPolinton-like NCLDV protein\n"
    )
    assert load_marker_annotation_index(path)[model.lower()]["family"] == "UNKNOWN"


@pytest.mark.parametrize("prefix,expected", [("NCLDV__", 0), ("", 0), ("EUK-pEVE__", 1)])
def test_structure_class_remains_reference_qualified(prefix, expected):
    marker = hit("one", "EVE_CAPSID_NOVEL_0008", top10_prefixes=prefix, top10_pidents="80")
    summary = classify_marker_hits([marker], metadata_families={})
    assert summary["vp_plv_markers"] == expected


def test_cress_profile_does_not_fall_into_metadata_family():
    marker = hit("one", "EVE_CAPSID_NOVEL_0033", top10_prefixes="EUK-pEVE__", top10_pidents="80")
    summary = classify_marker_hits([marker], metadata_families={"eve_capsid_novel_0033": "NCLDV"})
    assert summary["ncldv_markers"] == 0


def test_seven_associated_models_do_not_count_blocked_models_or_copies():
    models = ["VS000361", "VS000050", "VS000508", "VS000767", "VS000417", "VS000532", "VS000492"]
    index = {name.lower(): {"family": "NCLDV"} for name in models + ["VS000084", "VS000087"]}
    markers = [hit(str(i), model) for i, model in enumerate(models)]
    markers += [hit("blocked1", "VS000084"), hit("blocked2", "VS000087"), hit("copy", "VS000492")]
    summary = classify_marker_hits(markers, metadata_families=families.metadata_marker_families(index))
    assert summary["ncldv_markers"] == 7
    assert families.is_seed_blocked_model("VS000084")


def test_repeat_child_and_refresh_keep_one_noncellular_type():
    from virosync.pipeline.phase1.seed_merger import ConfirmedSeed
    from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
    from virosync.pipeline.phase2.seed_extension import refresh_boundary_classification
    from virosync.pipeline.phase2.terminal_repeats import _child_boundary, _marker_groups

    parent = RefinedBoundary(
        "s", 100, 400, original_start=100, original_end=400, confirmed_seeds=(ConfirmedSeed("seed", "s", 200, 240),)
    )
    markers = [
        hit("blocked", "VS000084", scaffold="s", start=150, end=180, hmm_evalue=1e-30),
        hit("viral", "VS000492", scaffold="s", start=200, end=240, hmm_evalue=1e-30),
    ]
    mapping = families.metadata_marker_families({"vs000084": {"family": "NCLDV"}, "vs000492": {"family": "NCLDV"}})
    groups = _marker_groups(parent, markers)
    child = _child_boundary(
        parent, start=120, end=300, groups=groups, scan_start=0, scan_end=500, metadata_families=mapping
    )
    assert child.region_classification_ncldv_markers == 1
    assert groups[0].is_seed_eligible is False
    assert (groups[0].start, groups[0].end) == (150, 180)
    refresh_boundary_classification([parent], markers, metadata_families=mapping)
    assert parent.region_classification_ncldv_markers == 1
    assert parent.confirmed_seeds == (ConfirmedSeed("seed", "s", 200, 240),)


def test_repeat_child_selects_one_best_model_without_shrinking_group():
    from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
    from virosync.pipeline.phase2.terminal_repeats import _child_boundary, _marker_groups

    parent = RefinedBoundary("s", 100, 400, original_start=100, original_end=400)
    markers = [
        hit("same|aa1-50", "GVOGm1", 100, scaffold="s", start=150, end=200),
        hit("same|aa51-100", "GVOGm2", 200, scaffold="s", start=200, end=250),
    ]
    groups = _marker_groups(parent, markers)
    child = _child_boundary(parent, start=100, end=400, groups=groups, scan_start=0, scan_end=500, metadata_families={})
    assert child.region_classification_ncldv_markers == 1
    assert (groups[0].start, groups[0].end) == (150, 250)


def test_reporting_retains_profile_association_without_seed_vote():
    marker = hit("one", "EVE_CAPSID_NOVEL_0033", top10_prefixes="EUK-pEVE__", top10_pidents="80")
    summary = summarize_marker_hits(
        [marker.hmm_target], {marker.hmm_target.lower(): {"family": "NCLDV"}}, hallmark_hits=[marker]
    )
    assert summary["family_counts"] == {"CRESS": 1}
    assert classify_marker_hits([marker], metadata_families={})["ncldv_markers"] == 0


def test_a1_surface_preserves_metadata_counts_and_report_association():
    from virosync.pipeline.phase1.hhg_seeding import Anchor
    from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
    from virosync.pipeline.phase3.phase1_surface import build_phase1_seed_surface

    anchor = Anchor("gene", "s", 100, 200, "+", "VS000492", 100, 1e-30, validation_status="validated")
    seed = MergedSeed("s", 50, 250, anchors=[anchor], confirmed_seeds=(ConfirmedSeed("seed", "s", 100, 200),))
    surface = build_phase1_seed_surface(
        [seed], metadata_families={"vs000492": "NCLDV"}, marker_annotation_index={"vs000492": {"family": "NCLDV"}}
    )
    assert surface.results[0].region_classification_ncldv_markers == 1
    assert surface.results[0].marker_family_hits == ["NCLDV"]
    assert surface.results[0].score_components["confidence_kind"] == "not_scored"


def test_repeat_worker_forwards_metadata_map(tmp_path):
    from dataclasses import replace

    from test_terminal_repeats import _adjacent_repeat_sequence, _query

    from virosync.pipeline.phase1.seed_merger import ConfirmedSeed
    from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
    from virosync.pipeline.phase2.terminal_repeats import (
        _build_scaffold_tasks,
        _plan_boundary_partition,
        _plan_scaffold_task,
    )

    sequence, _ = _adjacent_repeat_sequence()
    parent = RefinedBoundary(
        "contig",
        0,
        len(sequence),
        seed_id="seed-1",
        original_start=200,
        original_end=680,
        confirmed_seeds=(ConfirmedSeed("first", "contig", 200, 260), ConfirmedSeed("second", "contig", 620, 680)),
    )
    markers = [
        hit("first", "VS000492", scaffold="contig", start=220, end=250),
        hit("second", "VS000361", scaffold="contig", start=640, end=660),
    ]
    kwargs = dict(
        validated_markers=markers,
        boundary_diamond_query=_query(len(sequence)),
        boundary_taxonomy_map={},
        proteome_index={"contig": []},
        extension_bp=500,
        metadata_families={"vs000492": "NCLDV", "vs000361": "NCLDV"},
    )
    direct = _plan_boundary_partition(parent, allowed_interval=(0, len(sequence)), sequence=sequence, **kwargs)
    task = next(iter(_build_scaffold_tasks([parent], [(0, len(sequence))], tmp_path / "raw.fna", **kwargs)))
    path = tmp_path / "sequence.txt"
    path.write_text(sequence)
    worker = _plan_scaffold_task(replace(task, sequence_path=str(path)))[0][1]
    assert worker == direct
    assert [child.region_classification_ncldv_markers for child in worker] == [1, 1]


def test_blocked_model_keeps_reporting_provenance_without_classification_credit():
    index = {"vs000084": {"family": "NCLDV"}}
    assert summarize_marker_hits(["VS000084"], index)["family_counts"] == {"NCLDV": 1}
    assert families.metadata_marker_families(index) == {}


@pytest.mark.parametrize("rejected", [{"validation_status": "supported"}, {"pfam_override": "PF1"}])
def test_ordinary_filter_precedes_best_hit_for_the_only_protein(rejected):
    markers = [hit("only", "GVOGm1", 200, **rejected), hit("only", "VS000492", 100)]
    summary = classify_marker_hits(markers, metadata_families={"vs000492": "NCLDV"})
    assert summary["ncldv_markers"] == 1
    assert summary["classification"] == "NCLDV"
