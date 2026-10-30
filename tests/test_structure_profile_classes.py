"""A structure-derived profile gives its training-EVE class to hits that name no lineage.

The 33 ``EVE_CAPSID_NOVEL`` profiles come from EVEs of known classes. A profile
whose training EVEs hold one class at 80% or more (``STRUCTURE_PROFILE_CLASSES``)
lends that class to a marker hit whose qualified viral references are all
EUK-pEVE, which names no lineage. With any lineage hit, the lineage decides.
"""

from __future__ import annotations

import csv
from pathlib import Path
from types import SimpleNamespace

import pytest

from virosync.pipeline import host_signatures
from virosync.pipeline.phase1.hhg_seeding import Anchor
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
from virosync.pipeline.phase1.viral_markers import (
    classify_region_by_markers,
    get_region_classification_summary,
)
from virosync.pipeline.phase2 import terminal_repeats
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    EvidenceSynthesizerConfig,
    VerificationResult,
    consensus_taxonomy_class,
    infer_likely_family,
    mcp_vote_classes,
    taxonomy_class_votes,
)
from virosync.pipeline.phase3.mcp_detection import (
    CAPSID_EVIDENCE_MODELS,
    STRUCTURE_PROFILE_CLASSES,
    structure_profile_class,
)
from virosync.pipeline.phase3.output_generator import evaluate_v2_quality_gate
from virosync.pipeline.phase3.phase1_surface import build_phase1_seed_surface

CURATION_TABLE = (
    Path(__file__).resolve().parents[1] / "release-manifests" / "resources_v1_1_0" / "capsid_model_curation.tsv"
)
CLASS_COLUMNS = ("PPV", "NCLDV", "MIRUS", "PHAGE", "CRESS", "VIRAL_UNKNOWN")
PPV_PROFILE = "EVE_CAPSID_NOVEL_0008"
MCP_PROFILE = "EVE_CAPSID_NOVEL_0029"
UNCLASSED_PROFILE = "EVE_CAPSID_NOVEL_0002"
PEVE_ONLY = ("EUK-pEVE__", "80")


def _curated_profile_classes() -> dict[str, str]:
    """Class held by 80% of a profile's training EVEs, VIRAL_UNKNOWN in the denominator."""
    with CURATION_TABLE.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    classes = {}
    for row in rows:
        counts = {name: int(row[f"context_{name}"]) for name in CLASS_COLUMNS}
        total = sum(counts.values())
        name, count = max(counts.items(), key=lambda item: item[1])
        if name != "VIRAL_UNKNOWN" and count * 5 >= total * 4:
            classes[row["source_model"].lower()] = name
    return classes


def _name(number: int) -> str:
    return f"eve_capsid_novel_{number:04d}"


# --- The table ---


def test_profile_classes_match_the_curation_table() -> None:
    assert STRUCTURE_PROFILE_CLASSES == _curated_profile_classes()
    assert len(STRUCTURE_PROFILE_CLASSES) == 32
    assert all(name == name.lower() for name in STRUCTURE_PROFILE_CLASSES)


def test_profile_classes_hold_29_ppv_and_one_profile_each_of_mirus_cress_and_phage() -> None:
    by_class = {
        label: sorted(name for name, value in STRUCTURE_PROFILE_CLASSES.items() if value == label)
        for label in ("PPV", "MIRUS", "CRESS", "PHAGE", "NCLDV")
    }

    assert len(by_class["PPV"]) == 29
    assert by_class["MIRUS"] == [_name(26)]
    assert by_class["CRESS"] == [_name(33)]
    assert by_class["PHAGE"] == [_name(18)]
    assert by_class["NCLDV"] == []


def test_the_profile_without_a_dominant_class_has_none() -> None:
    assert _name(2) not in STRUCTURE_PROFILE_CLASSES
    assert structure_profile_class(UNCLASSED_PROFILE) == ""


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("EVE_CAPSID_NOVEL_0008", "PPV"),
        ("  eve_capsid_novel_0008 ", "PPV"),
        ("EVE_CAPSID_NOVEL_0026", "MIRUS"),
        ("EVE_CAPSID_NOVEL_0033", "CRESS"),
        ("EVE_CAPSID_NOVEL_0018", "PHAGE"),
        ("EVE_CAPSID_NOVEL_0029", "PPV"),
        ("EVE_CAPSID_NOVEL_0002", ""),
        ("GVOGm0003", ""),
        ("", ""),
        (None, ""),
    ],
)
def test_structure_profile_class_normalizes_the_name(name: str | None, expected: str) -> None:
    assert structure_profile_class(name) == expected


def test_the_non_capsid_profiles_with_a_class_stay_without_capsid_evidence() -> None:
    assert _name(26) in STRUCTURE_PROFILE_CLASSES
    assert _name(26) not in CAPSID_EVIDENCE_MODELS


# --- The vote ---


def _hit(
    gene: str,
    prefixes: str,
    pidents: str,
    porf_id: str = "g1",
    status: str = "validated",
    score: float = 100.0,
) -> dict:
    return {
        "hallmark_gene": gene,
        "porf_id": porf_id,
        "top10_prefixes": prefixes,
        "top10_pidents": pidents,
        "top10_targets": "",
        "validation_status": status,
        "hmm_score": score,
    }


def _gene(porf_id: str, prefixes: str, pidents: str) -> dict:
    return {
        "porf_id": porf_id,
        "top10_prefixes": prefixes,
        "top10_pidents": pidents,
        "top10_targets": "",
        "is_flanking": False,
    }


@pytest.mark.parametrize(
    ("number", "profile_class"),
    [(8, "PPV"), (26, "MIRUS"), (33, "CRESS"), (18, "PHAGE"), (29, "PPV")],
)
def test_peve_only_marker_of_a_classed_profile_votes_its_class_as_an_ordinary_marker(
    number: int, profile_class: str
) -> None:
    markers = [_hit(_name(number), *PEVE_ONLY)]
    genes = [_gene("g1", *PEVE_ONLY)]

    assert taxonomy_class_votes(markers, gene_taxonomy_records=genes) == {profile_class: 2}
    assert mcp_vote_classes(markers) == set()


def test_peve_only_marker_votes_the_marker_weight_without_a_gene_record() -> None:
    markers = [_hit(PPV_PROFILE, *PEVE_ONLY)]

    assert taxonomy_class_votes(markers) == {"PPV": 2}
    assert consensus_taxonomy_class(markers) == "PPV"


def test_peve_only_marker_confirmed_by_a_lineage_gene_record_gets_the_split_weights() -> None:
    markers = [_hit(PPV_PROFILE, *PEVE_ONLY)]
    genes = [_gene("g1", "PPV", "60")]

    assert taxonomy_class_votes(markers, gene_taxonomy_records=genes) == {"PPV": 3}


def test_the_mcp_profile_with_only_peve_hits_votes_ppv_without_the_mcp_weight_or_override() -> None:
    mcp = _hit(MCP_PROFILE, *PEVE_ONLY, porf_id="gm")
    genes = [_gene("a", "NCLDV", "80"), _gene("b", "NCLDV", "70"), _gene("c", "NCLDV", "60")]

    assert taxonomy_class_votes([mcp], gene_taxonomy_records=genes) == {"PPV": 2, "NCLDV": 3}
    assert consensus_taxonomy_class([mcp], gene_taxonomy_records=genes) == "NCLDV"


def test_a_classed_marker_with_a_lineage_hit_votes_its_lineage_as_before() -> None:
    markers = [_hit(PPV_PROFILE, "EUK-pEVE__,NCLDV__", "80,70")]

    assert taxonomy_class_votes(markers) == {"NCLDV": 2}
    assert consensus_taxonomy_class(markers) == "NCLDV"


def test_the_mcp_profile_with_a_lineage_hit_keeps_the_flat_mcp_vote() -> None:
    markers = [_hit(MCP_PROFILE, "EUK-pEVE__,NCLDV__", "80,70")]

    assert taxonomy_class_votes(markers) == {"NCLDV": 5}
    assert mcp_vote_classes(markers) == {"NCLDV"}


def test_an_unclassed_profile_with_only_peve_hits_still_casts_no_class() -> None:
    markers = [_hit(UNCLASSED_PROFILE, *PEVE_ONLY)]

    assert taxonomy_class_votes(markers) == {"VIRAL_UNKNOWN": 2}
    assert consensus_taxonomy_class(markers) == "VIRAL_UNKNOWN"


def test_a_classed_profile_without_any_qualified_hit_votes_nothing() -> None:
    markers = [_hit(PPV_PROFILE, "", "", status="validated_novel"), _hit(PPV_PROFILE, "EUK-pEVE__", "20", "g2")]

    assert taxonomy_class_votes(markers) == {}


def test_peve_only_gene_records_still_cast_no_class() -> None:
    genes = [_gene("g1", *PEVE_ONLY), _gene("g2", *PEVE_ONLY)]

    assert taxonomy_class_votes([], gene_taxonomy_records=genes) == {"VIRAL_UNKNOWN": 2}
    assert consensus_taxonomy_class([], gene_taxonomy_records=genes) == "VIRAL_UNKNOWN"


def test_unclassed_peve_only_markers_add_no_vote_beside_a_classed_one() -> None:
    markers = [_hit(PPV_PROFILE, *PEVE_ONLY, porf_id="g1"), _hit("GVOGm0003", *PEVE_ONLY, porf_id="g2")]

    assert taxonomy_class_votes(markers) == {"PPV": 2}


def test_classed_markers_of_two_classes_tie_into_viral_unknown() -> None:
    markers = [_hit(PPV_PROFILE, *PEVE_ONLY, porf_id="g1"), _hit(_name(26), *PEVE_ONLY, porf_id="g2")]

    assert taxonomy_class_votes(markers) == {"PPV": 2, "MIRUS": 2}
    assert consensus_taxonomy_class(markers) == "VIRAL_UNKNOWN"


def test_peve_only_mcp_does_not_win_deduplication_over_a_higher_scoring_classed_profile() -> None:
    markers = [
        _hit("GVOGm0003", *PEVE_ONLY, score=50.0),
        _hit(PPV_PROFILE, *PEVE_ONLY, score=100.0),
    ]

    assert taxonomy_class_votes(markers) == {"PPV": 2}
    assert consensus_taxonomy_class(markers) == "PPV"


def test_peve_only_mcp_profile_does_not_beat_a_higher_scoring_profile_of_another_class() -> None:
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    markers = [
        _hit(MCP_PROFILE, *PEVE_ONLY, porf_id="sc|1", score=90.0),
        _hit(_name(26), *PEVE_ONLY, porf_id="sc|1", score=100.0),
    ]
    result = VerificationResult(eve_id="e", scaffold="sc", start=0, end=3000)

    synthesizer._process_hallmark_hits(result, markers)
    synthesizer._assign_taxonomy_class(result, markers, [])

    assert taxonomy_class_votes(markers) == {"MIRUS": 2}
    assert result.taxonomy_class == "MIRUS"
    assert result.marker_dominant_family == "MIRUS"
    assert result.marker_family_hits == ["MIRUS"]


def test_mcp_with_a_lineage_hit_still_wins_deduplication_over_a_higher_scoring_profile() -> None:
    markers = [
        _hit("GVOGm0003", "NCLDV__", "80", score=50.0),
        _hit(UNCLASSED_PROFILE, "NCLDV__", "80", score=100.0),
    ]

    assert taxonomy_class_votes(markers) == {"NCLDV": 5}
    assert mcp_vote_classes(markers) == {"NCLDV"}


# --- Phase 3 family and gate ---


def _verify_one_marker_region(
    monkeypatch: pytest.MonkeyPatch,
    *,
    gene: str,
    end: int,
    region_class: str = "",
    prefixes: str = "EUK-pEVE",
    pidents: str = "80",
    status: str = "validated",
) -> VerificationResult:
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    monkeypatch.setattr(synthesizer, "_run_tiebreakers", lambda *args, **kwargs: None)
    hit = _hit(gene, prefixes, pidents, porf_id="sc|1", status=status)
    boundary = RefinedBoundary(
        scaffold="sc",
        start=0,
        end=end,
        confidence=0.95,
        predicted_family=region_class,
        confirmed_seeds=(ConfirmedSeed("marker-core", "sc", 0, 900),),
    )
    return synthesizer.verify_eve(
        boundary,
        window_features=[],
        hallmark_hits=[hit],
        gene_taxonomy_records=[_gene("sc|1", prefixes, pidents)],
        gene_taxonomy_summary={"total": 1},
    )


@pytest.mark.parametrize("region_class", ["", "UNKNOWN", "PPV"])
def test_a_classed_peve_only_region_publishes_its_class_and_takes_that_gate_path(
    monkeypatch: pytest.MonkeyPatch, region_class: str
) -> None:
    result = _verify_one_marker_region(monkeypatch, gene=PPV_PROFILE, end=3000, region_class=region_class)
    decision = evaluate_v2_quality_gate(result)

    assert result.taxonomy_class == "PPV"
    assert result.likely_family == "PPV"
    assert result.marker_dominant_family == "PPV"
    assert result.marker_family_hits == ["PPV"]
    assert (decision.kept, decision.effective_class, decision.reason) == (True, "PPV", "small_dna_high_medium_pass")


def test_a_short_classed_peve_only_region_fails_the_ppv_size_gate(monkeypatch: pytest.MonkeyPatch) -> None:
    result = _verify_one_marker_region(monkeypatch, gene=PPV_PROFILE, end=1000)
    decision = evaluate_v2_quality_gate(result)

    assert result.likely_family == "PPV"
    assert (decision.kept, decision.reason) == (False, "small_dna_high_medium_gate")


@pytest.mark.parametrize(
    ("number", "profile_class", "family"),
    [(26, "MIRUS", "MIRUS"), (33, "CRESS", "CRESS"), (18, "PHAGE", "UNKNOWN")],
)
def test_the_other_classes_publish_their_class_and_family(
    monkeypatch: pytest.MonkeyPatch, number: int, profile_class: str, family: str
) -> None:
    result = _verify_one_marker_region(monkeypatch, gene=_name(number), end=3000)

    assert result.taxonomy_class == profile_class
    assert result.likely_family == family


@pytest.mark.parametrize("gene", [UNCLASSED_PROFILE, "GVOGm0003"])
def test_peve_only_region_without_a_classed_profile_stays_viral_unknown_and_mixed(
    monkeypatch: pytest.MonkeyPatch, gene: str
) -> None:
    result = _verify_one_marker_region(monkeypatch, gene=gene, end=3000)
    decision = evaluate_v2_quality_gate(result)

    assert result.taxonomy_class == "VIRAL_UNKNOWN"
    assert result.likely_family == "MIXED"
    assert (decision.effective_class, decision.reason) == ("MIXED", "mixed_high_medium_pass")


def test_a_classed_profile_with_a_lineage_hit_keeps_the_lineage_and_gives_no_class_family(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    result = _verify_one_marker_region(
        monkeypatch, gene=PPV_PROFILE, end=3000, prefixes="EUK-pEVE,NCLDV", pidents="80,70"
    )

    assert result.taxonomy_class == "NCLDV"
    assert result.marker_family_hits == []
    assert result.likely_family == "UNKNOWN"


def test_a_classed_profile_without_references_gives_no_class_or_family(monkeypatch: pytest.MonkeyPatch) -> None:
    result = _verify_one_marker_region(
        monkeypatch, gene=PPV_PROFILE, end=3000, prefixes="", pidents="", status="validated_novel"
    )
    decision = evaluate_v2_quality_gate(result)

    assert result.taxonomy_class == "VIRAL_UNKNOWN"
    assert result.marker_family_hits == []
    assert result.likely_family == "UNKNOWN"
    assert (decision.kept, decision.reason) == (False, "unsupported_class")


# --- Phase 1 region family ---


@pytest.mark.parametrize(
    ("class_by_marker", "expected"),
    [
        ({_name(8): "PPV"}, "PPV"),
        ({_name(26): "MIRUS"}, "MIRUS"),
        ({_name(33): "CRESS"}, "UNKNOWN"),
        ({_name(18): "PHAGE"}, "UNKNOWN"),
        ({}, "UNKNOWN"),
    ],
)
def test_region_family_counts_a_marker_with_a_supplied_class(class_by_marker: dict[str, str], expected: str) -> None:
    marker_names = {"EVE_CAPSID_NOVEL_0008", "EVE_CAPSID_NOVEL_0026", "EVE_CAPSID_NOVEL_0033", "EVE_CAPSID_NOVEL_0018"}
    names = {name for name in marker_names if name.lower() in class_by_marker}

    assert classify_region_by_markers(names, class_by_marker=class_by_marker) == expected


def test_region_family_summary_counts_the_class_markers() -> None:
    summary = get_region_classification_summary(
        {"EVE_CAPSID_NOVEL_0008", "GVOGm0003"}, class_by_marker={_name(8): "PPV"}
    )

    assert summary["classification"] == "MIXED"
    assert (summary["ncldv_markers"], summary["vp_plv_markers"], summary["mirus_markers"]) == (1, 1, 0)


def _anchor(
    gene: str,
    prefixes: str,
    pidents: str,
    porf_id: str = "gene-a",
    *,
    targets: str = "",
    score: float = 50.0,
) -> Anchor:
    return Anchor(
        porf_id=porf_id,
        scaffold="ctg",
        start=10,
        end=40,
        strand="+",
        hallmark_gene=gene,
        score=score,
        evalue=1e-20,
        validation_status="validated",
        top10_prefixes=prefixes,
        top10_pidents=pidents,
        top10_targets=targets,
    )


def _seed(*anchors: Anchor, predicted_family: str = "") -> MergedSeed:
    return MergedSeed(
        scaffold="ctg",
        start=0,
        end=3000,
        seed_id="seed",
        confirmed_seeds=(ConfirmedSeed("seed", "ctg", 0, 3000),),
        sources=["marker_validation"],
        hhg_score=0.0,
        novelty_score=0.0,
        compositional_score=0.0,
        max_kfd=0.0,
        gc_deviation=0.0,
        cub_deviation=0.0,
        confidence="medium",
        predicted_family=predicted_family,
        region_classification_ncldv_markers=0,
        region_classification_vp_plv_markers=0,
        region_classification_mirus_markers=0,
        anchors=list(anchors),
        hhg_anchors=[],
    )


@pytest.mark.parametrize(
    ("gene", "prefixes", "pidents", "family", "vp_plv"),
    [
        (PPV_PROFILE, "EUK-pEVE__", "80", "PPV", 1),
        (_name(26), "EUK-pEVE__", "80", "MIRUS", 0),
        (PPV_PROFILE, "EUK-pEVE__,NCLDV__", "80,70", "UNKNOWN", 0),
        (PPV_PROFILE, "", "", "UNKNOWN", 0),
        (UNCLASSED_PROFILE, "EUK-pEVE__", "80", "UNKNOWN", 0),
    ],
)
def test_seed_classification_uses_the_class_of_a_peve_only_profile(
    gene: str, prefixes: str, pidents: str, family: str, vp_plv: int
) -> None:
    seed = _seed(_anchor(gene, prefixes, pidents))

    seed.compute_classification(None)

    assert seed.predicted_family == family
    assert seed.region_classification_vp_plv_markers == vp_plv


def test_seed_classification_mixes_a_class_marker_with_a_name_marker_of_another_family() -> None:
    seed = _seed(
        _anchor(PPV_PROFILE, "EUK-pEVE__", "80", porf_id="gene-a"),
        _anchor("GVOGm0003", "NCLDV__", "80", porf_id="gene-b"),
    )

    seed.compute_classification(None)

    assert seed.predicted_family == "MIXED"


@pytest.mark.parametrize(
    ("prefixes", "pidents", "family"),
    [("EUK-pEVE__", "80", "PPV"), ("NCLDV__", "80", "UNKNOWN")],
)
def test_terminal_repeat_children_keep_the_class_of_a_peve_only_profile(
    prefixes: str, pidents: str, family: str
) -> None:
    parent = RefinedBoundary(scaffold="ctg", start=0, end=3000, original_start=0, original_end=3000)
    marker = type(
        "Marker",
        (),
        {
            "scaffold": "ctg",
            "start": 100,
            "end": 900,
            "query_porf": "p1",
            "validation_status": "validated",
            "hmm_target": PPV_PROFILE,
            "top10_prefixes": prefixes,
            "top10_pidents": pidents,
        },
    )()
    groups = terminal_repeats._marker_groups(parent, [marker])

    child = terminal_repeats._child_boundary(parent, start=0, end=1000, groups=groups, scan_start=0, scan_end=1000)

    assert child.predicted_family == family


def _repeat_marker(gene: str, porf_id: str, *, score: float = 100.0, start: int = 100) -> SimpleNamespace:
    return SimpleNamespace(
        scaffold="ctg",
        start=start,
        end=start + 300,
        query_porf=porf_id,
        validation_status="validated",
        hmm_target=gene,
        hmm_score=score,
        top10_prefixes="EUK-pEVE__",
        top10_pidents="80",
    )


def test_terminal_repeat_children_count_one_profile_per_protein() -> None:
    parent = RefinedBoundary(scaffold="ctg", start=0, end=3000, original_start=0, original_end=3000)
    markers = [
        _repeat_marker("mirus_mcp_1", "g1", start=100),
        _repeat_marker(_name(26), "g2", score=100.0, start=500),
        _repeat_marker(PPV_PROFILE, "g2", score=90.0, start=500),
        _repeat_marker("plv_pc_054", "g3", start=900),
    ]
    groups = terminal_repeats._marker_groups(parent, markers)

    child = terminal_repeats._child_boundary(parent, start=0, end=1300, groups=groups, scan_start=0, scan_end=1300)

    assert child.predicted_family == "MIRUS"
    assert (child.region_classification_mirus_markers, child.region_classification_vp_plv_markers) == (2, 1)


def test_terminal_repeat_children_keep_the_higher_scoring_class_of_one_protein() -> None:
    parent = RefinedBoundary(scaffold="ctg", start=0, end=3000, original_start=0, original_end=3000)
    markers = [
        _repeat_marker(_name(26), "g1", score=90.0),
        _repeat_marker(PPV_PROFILE, "g1", score=100.0),
    ]
    groups = terminal_repeats._marker_groups(parent, markers)

    child = terminal_repeats._child_boundary(parent, start=0, end=1000, groups=groups, scan_start=0, scan_end=1000)

    assert child.predicted_family == "PPV"


# --- A1 seed export ---


def _a1_result(*anchors: Anchor, predicted_family: str = "") -> VerificationResult:
    return build_phase1_seed_surface([_seed(*anchors, predicted_family=predicted_family)]).results[0]


@pytest.mark.parametrize(
    ("number", "taxonomy_class", "likely_family"),
    [(8, "PPV", "PPV"), (26, "MIRUS", "MIRUS"), (33, "CRESS", "CRESS"), (18, "PHAGE", "UNKNOWN")],
)
def test_a1_publishes_the_class_of_a_peve_only_classed_profile(
    number: int, taxonomy_class: str, likely_family: str
) -> None:
    result = _a1_result(_anchor(_name(number), "EUK-pEVE__", "80"))

    assert (result.taxonomy_class, result.likely_family) == (taxonomy_class, likely_family)


@pytest.mark.parametrize("gene", [UNCLASSED_PROFILE, "GVOGm0003"])
def test_a1_keeps_viral_unknown_and_mixed_without_a_classed_profile(gene: str) -> None:
    result = _a1_result(_anchor(gene, "EUK-pEVE__", "80"))

    assert (result.taxonomy_class, result.likely_family) == ("VIRAL_UNKNOWN", "MIXED")


def test_a1_ties_two_classed_profiles_of_different_classes_into_viral_unknown() -> None:
    result = _a1_result(
        _anchor(PPV_PROFILE, "EUK-pEVE__", "80", porf_id="gene-a"),
        _anchor(_name(26), "EUK-pEVE__", "80", porf_id="gene-b"),
    )

    assert (result.taxonomy_class, result.likely_family) == ("VIRAL_UNKNOWN", "MIXED")


def test_a1_keeps_the_seed_family_when_a_lineage_hit_joins_a_classed_profile() -> None:
    result = _a1_result(_anchor(PPV_PROFILE, "EUK-pEVE__,NCLDV__", "80,70"))

    assert (result.taxonomy_class, result.likely_family) == ("UNKNOWN", "UNKNOWN")


def test_a1_votes_a_classed_peve_only_marker_beside_a_marker_with_lineage_evidence() -> None:
    classed = _anchor(PPV_PROFILE, "EUK-pEVE__", "80", porf_id="gene-a")
    lineage = _anchor(UNCLASSED_PROFILE, "NCLDV__", "60", porf_id="gene-b")
    phase3_markers = [
        _hit(PPV_PROFILE, "EUK-pEVE__", "80", porf_id="gene-a"),
        _hit(UNCLASSED_PROFILE, "NCLDV__", "60", porf_id="gene-b"),
    ]

    result = _a1_result(classed, lineage, predicted_family="PPV")

    assert taxonomy_class_votes(phase3_markers) == {"PPV": 2, "NCLDV": 2}
    assert result.taxonomy_class == "VIRAL_UNKNOWN"
    assert result.likely_family == "PPV"


def test_a1_gives_the_mcp_override_of_a_lineage_marker_beside_a_classed_peve_only_marker() -> None:
    classed = _anchor(PPV_PROFILE, "EUK-pEVE__", "80", porf_id="gene-a")
    mcp = _anchor("GVOGm0003", "NCLDV__", "60", porf_id="gene-b")

    result = _a1_result(classed, mcp, predicted_family="MIXED")

    assert (result.taxonomy_class, result.likely_family) == ("NCLDV", "MIXED")


def test_a1_without_a_classed_peve_only_marker_keeps_the_seed_family_beside_lineage_evidence() -> None:
    peve_only = _anchor("GVOGm0003", "EUK-pEVE__", "80", porf_id="gene-a")
    lineage = _anchor(UNCLASSED_PROFILE, "NCLDV__", "60", porf_id="gene-b")

    result = _a1_result(peve_only, lineage, predicted_family="NCLDV")

    assert (result.taxonomy_class, result.likely_family) == ("NCLDV", "NCLDV")


# --- Phase 1 has no CRESS or PHAGE region family ---


@pytest.mark.parametrize(("number", "profile_class"), [(33, "CRESS"), (18, "PHAGE")])
def test_cress_and_phage_profile_classes_apply_from_phase_3_on_not_in_phase_1(
    monkeypatch: pytest.MonkeyPatch, number: int, profile_class: str
) -> None:
    gene = _name(number)
    seed = _seed(_anchor(gene, "EUK-pEVE__", "80"))
    parent = RefinedBoundary(scaffold="ctg", start=0, end=3000, original_start=0, original_end=3000)
    marker = type(
        "Marker",
        (),
        {
            "scaffold": "ctg",
            "start": 100,
            "end": 900,
            "query_porf": "p1",
            "validation_status": "validated",
            "hmm_target": gene,
            "top10_prefixes": "EUK-pEVE__",
            "top10_pidents": "80",
        },
    )()
    groups = terminal_repeats._marker_groups(parent, [marker])

    seed.compute_classification(None)
    child = terminal_repeats._child_boundary(parent, start=0, end=1000, groups=groups, scan_start=0, scan_end=1000)
    phase3 = _verify_one_marker_region(monkeypatch, gene=gene, end=3000)

    assert seed.predicted_family == "UNKNOWN"
    assert child.predicted_family == "UNKNOWN"
    assert phase3.taxonomy_class == profile_class


def test_a1_resolves_the_lineage_of_a_reference_with_the_phase_3_taxonomy_lookup(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    lookup = {"PHAGE__VARDNA__Sputnik": "PHAGE|Preplasmiviricota|Lavidaviridae|Sputnikvirus"}
    monkeypatch.setattr(host_signatures, "_TAXONOMY_LOOKUP", lookup)
    target = "PHAGE__VARDNA__Sputnik|p1"
    phase3_markers = [
        {**_hit(PPV_PROFILE, "EUK-pEVE__", "80", porf_id="gene-a")},
        {**_hit(UNCLASSED_PROFILE, "PHAGE__", "60", porf_id="gene-b"), "top10_targets": target},
    ]
    classed = _anchor(PPV_PROFILE, "EUK-pEVE__", "80", porf_id="gene-a")
    lineage = _anchor(UNCLASSED_PROFILE, "PHAGE__", "60", porf_id="gene-b", targets=target)

    result = _a1_result(classed, lineage, predicted_family="PPV")

    assert taxonomy_class_votes(phase3_markers, taxonomy_lookup=lookup) == {"PPV": 4}
    assert result.taxonomy_class == "PPV"


def test_a1_breaks_an_equal_score_tie_between_profiles_of_one_protein_as_the_vote_does() -> None:
    first = _anchor(_name(26), "EUK-pEVE__", "80", porf_id="gene-a", score=100.0)
    second = _anchor(PPV_PROFILE, "EUK-pEVE__", "80", porf_id="gene-a", score=100.0)
    phase3_markers = [_hit(_name(26), "EUK-pEVE__", "80"), _hit(PPV_PROFILE, "EUK-pEVE__", "80")]

    result = _a1_result(first, second)

    assert consensus_taxonomy_class(phase3_markers) == "MIRUS"
    assert result.taxonomy_class == "MIRUS"
    assert result.likely_family == "MIRUS"
    assert result.marker_dominant_family == "MIRUS"
    assert result.marker_family_hits == ["MIRUS"]


@pytest.mark.parametrize(
    ("numbers", "taxonomy_class", "likely_family"),
    [
        ((8, 9, 26), "PPV", "MIXED"),
        ((8, 9, 10, 26), "PPV", "PPV"),
    ],
)
def test_a1_and_phase_3_publish_the_same_likely_family_for_mixed_profile_classes(
    numbers: tuple[int, ...], taxonomy_class: str, likely_family: str
) -> None:
    anchors = [_anchor(_name(number), "EUK-pEVE__", "80", porf_id=f"gene-{number}") for number in numbers]
    markers = [_hit(_name(number), "EUK-pEVE__", "80", porf_id=f"gene-{number}") for number in numbers]
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    phase3 = VerificationResult(eve_id="e", scaffold="sc", start=0, end=3000, region_classification="MIXED")
    synthesizer._process_hallmark_hits(phase3, markers)
    synthesizer._assign_taxonomy_class(phase3, markers, [])
    phase3.likely_family = infer_likely_family(phase3)

    a1 = _a1_result(*anchors, predicted_family="MIXED")

    assert (phase3.taxonomy_class, phase3.likely_family) == (taxonomy_class, likely_family)
    assert (a1.taxonomy_class, a1.likely_family) == (taxonomy_class, likely_family)
