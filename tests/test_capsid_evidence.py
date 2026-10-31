"""Capsid-type structure-derived profiles count as capsid evidence, not as MCPs.

``has_capsid_evidence`` is true when ``has_mcp`` is true or any hallmark gene is
one of the 27 capsid-type ``EVE_CAPSID_NOVEL`` profiles. Raw annotation and
strict validation remain separate from catalogue-qualified
scoring, which requires retained HMM measurements.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest

from virosync.orchestration._flows.single_genome.phase1 import _candidate_regions_to_seeds
from virosync.pipeline.hallmark_evidence import HallmarkEvidence
from virosync.pipeline.phase1.hhg_seeding import Anchor, HMMHit
from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit as MarkerHit
from virosync.pipeline.phase1.marker_validation import filter_validated_markers
from virosync.pipeline.phase1.region_assembly import (
    ValidatedMarkerHit as RegionMarkerHit,
)
from virosync.pipeline.phase1.region_assembly import assemble_candidate_regions
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    EvidenceSynthesizerConfig,
    VerificationResult,
    calculate_eve_confidence,
    infer_likely_family,
)
from virosync.pipeline.phase3.output_generator import evaluate_v2_quality_gate
from virosync.pipeline.phase3.phase1_surface import build_phase1_seed_surface

CAPSID_MODEL = "EVE_CAPSID_NOVEL_0008"
NON_CAPSID_MODEL = "EVE_CAPSID_NOVEL_0001"


def _result(**fields) -> VerificationResult:
    return VerificationResult(eve_id="eve", scaffold="scaf", start=0, end=3000, **fields)


# --- Flag basics ---


def test_capsid_evidence_flag_defaults_to_false_and_is_reported() -> None:
    result = _result(has_capsid_evidence=True)

    assert _result().has_capsid_evidence is False
    assert result.to_dict()["has_capsid_evidence"] is True
    assert result.to_dict()["has_mcp"] is False


def _hit(gene: str, porf_id: str = "p1", score: float = 50.0, status: str = "validated") -> dict:
    return {"hallmark_gene": gene, "porf_id": porf_id, "hmm_score": score, "validation_status": status}


def test_public_capsid_gene_ids_preserve_validated_annotations_without_scoring_them() -> None:
    hits = [
        _hit(CAPSID_MODEL, "capsid|aa1-20"),
        _hit(NON_CAPSID_MODEL, "capsid|aa1-20", score=90),
        _hit(CAPSID_MODEL, "novel|aa1-20", status="validated_novel"),
        _hit(CAPSID_MODEL, "unvalidated", status="unvalidated"),
        _hit(NON_CAPSID_MODEL, "protease"),
        _hit("GVOGm0003", "mcp"),
    ]
    original_hits = deepcopy(hits)
    result = _result()
    EvidenceSynthesizer()._process_hallmark_hits(result, hits)
    exported = result.to_dict()
    assert exported["capsid_gene_ids"] == ["capsid", "mcp", "novel"]
    assert exported["mcp_gene_ids"] == ["mcp"]
    assert exported["hallmark_evidence"] == []
    assert hits == original_hits
    score = calculate_eve_confidence(result, 0)
    result.capsid_gene_ids = []
    assert calculate_eve_confidence(result, 0) == score


def test_public_mcp_gene_ids_preserve_outscored_validated_mcp_annotations() -> None:
    hits = [
        _hit("EVE_CAPSID_NOVEL_0029", "mcp|aa1-20", status="validated_novel"),
        _hit(NON_CAPSID_MODEL, "mcp|aa1-20", score=90),
        _hit("EVE_CAPSID_NOVEL_0029", "unvalidated", status="unvalidated"),
    ]
    result = _result()
    EvidenceSynthesizer()._process_hallmark_hits(result, hits)

    assert result.to_dict()["mcp_gene_ids"] == ["mcp"]
    assert result.to_dict()["capsid_gene_ids"] == ["mcp"]
    assert result.hallmark_genes == [NON_CAPSID_MODEL, "EVE_CAPSID_NOVEL_0029"]
    assert result.hallmark_count == 2
    score = calculate_eve_confidence(result, 0)
    result.mcp_gene_ids = []
    assert calculate_eve_confidence(result, 0) == score


@pytest.mark.parametrize(
    ("gene", "has_mcp", "has_capsid_evidence", "mcp_gene_ids"),
    [
        (CAPSID_MODEL, False, True, []),
        ("EVE_CAPSID_NOVEL_0029", True, True, ["p1"]),
        ("GVOGm0003", True, True, ["p1"]),
        (NON_CAPSID_MODEL, False, False, []),
        ("EVE_CAPSID_NOVEL_0018", False, False, []),
        ("GVOGm0100", False, False, []),
    ],
)
def test_hallmark_processing_sets_capsid_evidence_apart_from_mcp(
    gene: str,
    has_mcp: bool,
    has_capsid_evidence: bool,
    mcp_gene_ids: list[str],
) -> None:
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    result = _result()

    synthesizer._process_hallmark_hits(result, [_hit(gene)])

    assert result.has_mcp is has_mcp
    assert result.has_capsid_evidence is has_capsid_evidence
    assert result.mcp_gene_ids == mcp_gene_ids


def test_one_capsid_profile_among_other_hallmarks_is_enough() -> None:
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    result = _result()

    synthesizer._process_hallmark_hits(
        result,
        [_hit("GVOGm0100", "p1"), _hit(CAPSID_MODEL, "p2"), _hit(NON_CAPSID_MODEL, "p3")],
    )

    assert result.has_capsid_evidence is True
    assert result.has_mcp is False


@pytest.mark.parametrize("status", ["validated", "validated_novel"])
def test_stronger_non_capsid_profile_on_the_same_protein_does_not_hide_capsid_evidence(status: str) -> None:
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    result = _result()

    synthesizer._process_hallmark_hits(
        result,
        [_hit(NON_CAPSID_MODEL, "p1", 90.0), _hit(CAPSID_MODEL, "p1", 50.0, status)],
    )

    assert result.hallmark_genes == [NON_CAPSID_MODEL]
    assert result.hallmark_count == 1
    assert result.has_capsid_evidence is True
    assert result.has_mcp is False
    assert result.mcp_gene_ids == []


def test_hidden_unvalidated_capsid_hit_gives_no_capsid_evidence() -> None:
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    result = _result()

    synthesizer._process_hallmark_hits(
        result,
        [_hit(NON_CAPSID_MODEL, "p1", 90.0), _hit(CAPSID_MODEL, "p1", 50.0, "supported")],
    )

    assert result.has_capsid_evidence is False


@pytest.mark.parametrize("status", ["supported", "unvalidated", ""])
def test_unvalidated_capsid_profile_that_wins_deduplication_gives_no_capsid_evidence(status: str) -> None:
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    result = _result()

    synthesizer._process_hallmark_hits(
        result,
        [_hit(NON_CAPSID_MODEL, "p1", 90.0), _hit(CAPSID_MODEL, "p1", 100.0, status)],
    )

    assert result.hallmark_genes == [CAPSID_MODEL]
    assert result.has_capsid_evidence is False
    assert result.has_mcp is False


@pytest.mark.parametrize("status", ["supported", "unvalidated"])
def test_lone_unvalidated_capsid_hit_gives_no_capsid_evidence(status: str) -> None:
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    result = _result()

    synthesizer._process_hallmark_hits(result, [_hit(CAPSID_MODEL, "p1", 100.0, status)])

    assert result.has_capsid_evidence is False


def test_stronger_non_capsid_profile_does_not_hide_a_validated_mcp_hit_from_capsid_evidence() -> None:
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    result = _result()

    synthesizer._process_hallmark_hits(
        result,
        [_hit(NON_CAPSID_MODEL, "p1", 90.0), _hit("EVE_CAPSID_NOVEL_0029", "p1", 50.0)],
    )

    assert result.has_capsid_evidence is True


def _seed(*anchors: Anchor) -> MergedSeed:
    return MergedSeed(
        scaffold="ctg",
        start=0,
        end=100,
        seed_id="seed",
        confirmed_seeds=(ConfirmedSeed("seed", "ctg", 0, 100),),
        sources=["marker_validation"],
        hhg_score=0.0,
        novelty_score=0.0,
        compositional_score=0.0,
        max_kfd=0.0,
        gc_deviation=0.0,
        cub_deviation=0.0,
        confidence="medium",
        predicted_family="UNKNOWN",
        region_classification_ncldv_markers=0,
        region_classification_vp_plv_markers=0,
        region_classification_mirus_markers=0,
        anchors=list(anchors),
        hhg_anchors=[],
    )


def _anchor(gene: str, porf_id: str = "gene-a", score: float = 50.0, status: str = "validated") -> Anchor:
    return Anchor(
        porf_id=porf_id,
        scaffold="ctg",
        start=10,
        end=40,
        strand="+",
        hallmark_gene=gene,
        score=score,
        evalue=1e-20,
        validation_status=status,
    )


def test_a1_seed_surface_marks_capsid_evidence_without_an_mcp() -> None:
    capsid = build_phase1_seed_surface([_seed(_anchor(CAPSID_MODEL))]).results[0]
    mcp = build_phase1_seed_surface([_seed(_anchor("EVE_CAPSID_NOVEL_0029"))]).results[0]
    other = build_phase1_seed_surface([_seed(_anchor(NON_CAPSID_MODEL))]).results[0]

    assert (capsid.has_capsid_evidence, capsid.has_mcp, capsid.mcp_gene_ids) == (True, False, [])
    assert (mcp.has_capsid_evidence, mcp.has_mcp, mcp.mcp_gene_ids) == (True, True, ["gene-a"])
    assert (other.has_capsid_evidence, other.has_mcp) == (False, False)
    assert capsid.to_dict()["capsid_gene_ids"] == ["gene-a"]
    assert mcp.to_dict()["capsid_gene_ids"] == ["gene-a"]
    assert other.to_dict()["capsid_gene_ids"] == []


@pytest.mark.parametrize("status", ["validated", "validated_novel"])
def test_a1_public_mcp_annotation_survives_a_stronger_non_mcp_profile(status: str) -> None:
    seed = _seed(
        _anchor(NON_CAPSID_MODEL, score=90.0),
        _anchor("EVE_CAPSID_NOVEL_0029", score=40.0, status=status),
    )
    result = build_phase1_seed_surface([seed]).results[0]

    assert result.to_dict()["mcp_gene_ids"] == ["gene-a"]
    assert result.has_mcp is False
    assert result.hallmark_genes == [NON_CAPSID_MODEL]
    assert result.hallmark_count == 1
    assert result.final_confidence == 0.0


@pytest.mark.parametrize("status", ["validated", "validated_novel"])
def test_a1_stronger_non_capsid_profile_on_one_gene_does_not_hide_capsid_evidence(status: str) -> None:
    seed = _seed(_anchor(NON_CAPSID_MODEL, score=90.0), _anchor(CAPSID_MODEL, score=50.0, status=status))

    result = build_phase1_seed_surface([seed]).results[0]

    assert result.hallmark_genes == [NON_CAPSID_MODEL]
    assert result.has_capsid_evidence is True
    assert result.has_mcp is False
    assert result.mcp_gene_ids == []


@pytest.mark.parametrize("status", ["supported", "unvalidated", ""])
@pytest.mark.parametrize("capsid_score", [100.0, 50.0, None])
def test_a1_unvalidated_capsid_anchor_gives_no_capsid_evidence(status: str, capsid_score: float | None) -> None:
    capsid = _anchor(CAPSID_MODEL, score=capsid_score or 100.0, status=status)
    anchors = [capsid] if capsid_score is None else [_anchor(NON_CAPSID_MODEL, score=90.0), capsid]

    result = build_phase1_seed_surface([_seed(*anchors)]).results[0]

    assert result.has_capsid_evidence is False
    assert result.has_mcp is False


def test_a1_capsid_anchor_outside_the_seed_gives_no_capsid_evidence() -> None:
    outside = Anchor(
        porf_id="gene-b",
        scaffold="ctg",
        start=200,
        end=240,
        strand="+",
        hallmark_gene=CAPSID_MODEL,
        score=50.0,
        evalue=1e-20,
        validation_status="validated",
    )

    result = build_phase1_seed_surface([_seed(_anchor(NON_CAPSID_MODEL), outside)]).results[0]

    assert result.has_capsid_evidence is False


# --- Phase 1 seeding ---


def _region_marker(
    model: str,
    status: str,
    prefixes: str = "NCLDV__",
    pidents: str = "",
    targets: str = "",
) -> RegionMarkerHit:
    return RegionMarkerHit(
        query_porf="ctgA_1",
        scaffold="ctgA",
        start=1000,
        end=1300,
        strand="+",
        hmm_target=model,
        hmm_score=120.0,
        hmm_evalue=1e-30,
        validation_status=status,
        top10_prefixes=prefixes,
        best_hit_target="",
        best_hit_pident=0.0,
        best_hit_bits=0.0,
        has_ncldv=0,
        has_mirus=0,
        has_plv=0,
        has_vp=0,
        has_viral=0,
        top10_pidents=pidents,
        top10_targets=targets,
    )


@pytest.mark.parametrize("status", ["validated", "validated_novel", "supported"])
def test_phase1_seeds_carry_the_marker_validation_status_on_their_anchors(status: str) -> None:
    marker = _region_marker(CAPSID_MODEL, status)
    region = SimpleNamespace(scaffold="ctgA", start=0, end=2000, markers=[marker], predicted_family="UNKNOWN")

    (seed,) = _candidate_regions_to_seeds([region], [marker], [], [])

    assert [anchor.validation_status for anchor in seed.anchors] == [status]
    assert [anchor.validation_status for anchor in seed.hhg_anchors] == [status]


def test_phase1_seeds_carry_the_marker_top10_taxonomy_on_their_anchors() -> None:
    marker = _region_marker(
        CAPSID_MODEL,
        "validated",
        prefixes="EUK-pEVE__,NCLDV__",
        pidents="80.0,45.5",
        targets="EUK-pEVE__a|p,NCLDV__b|q",
    )
    region = SimpleNamespace(scaffold="ctgA", start=0, end=2000, markers=[marker], predicted_family="UNKNOWN")

    (seed,) = _candidate_regions_to_seeds([region], [marker], [], [])

    for anchor in (*seed.anchors, *seed.hhg_anchors):
        assert (anchor.top10_prefixes, anchor.top10_pidents, anchor.top10_targets) == (
            "EUK-pEVE__,NCLDV__",
            "80.0,45.5",
            "EUK-pEVE__a|p,NCLDV__b|q",
        )


@pytest.mark.parametrize(
    ("model", "can_seed"),
    [
        (CAPSID_MODEL, True),
        ("EVE_CAPSID_NOVEL_0029", True),
        ("OG1352", True),
        (NON_CAPSID_MODEL, False),
        ("EVE_CAPSID_NOVEL_0018", False),
        ("GVOGm0100", False),
    ],
)
def test_validated_novel_hit_seeds_only_for_mcp_or_capsid_evidence(model: str, can_seed: bool) -> None:
    assert _region_marker(model, "validated_novel", prefixes="").is_valid_seed_marker is can_seed


@pytest.mark.parametrize("model", [CAPSID_MODEL, NON_CAPSID_MODEL, "EVE_CAPSID_NOVEL_0018"])
def test_capsid_type_hit_validated_on_peve_references_seeds_like_any_validated_marker(model: str) -> None:
    marker = _region_marker(model, "validated", prefixes="EUK-pEVE__")

    assert marker.is_valid_seed_marker is True
    assert marker.is_mcp is False


@pytest.mark.parametrize("model", [CAPSID_MODEL, NON_CAPSID_MODEL])
def test_marker_validated_on_peve_references_alone_seeds_a_region(tmp_path: Path, model: str) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(f">contig_1_1 # 1000 # 1300 # + # ID=1_1;partial=00\n{'M' * 100}\n")
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("contig_1_1\tEUK-pEVE__EP00022|prot1\t1e-40\t200\t80.0\t95.0\n")
    hit = HMMHit(
        query_name="contig_1_1",
        target_name=model,
        score=120.0,
        evalue=1e-40,
        domain_score=120.0,
        query_start=1,
        query_end=100,
    )
    genome = tmp_path / "genome.fna"
    genome.write_text(f">contig_1\n{'ACGT' * 5000}\n", encoding="utf-8")

    (marker,) = filter_validated_markers([hit], diamond, proteome)
    regions = assemble_candidate_regions(
        [marker],
        genome,
        proteome,
        tmp_path / "out",
        initial_window_bp=100,
        initial_window_genes=1,
        min_markers_initial=1,
        extension_kb=5,
        merge_distance=0,
    )

    assert marker.validation_status == "validated"
    assert len(regions) == 1


def test_unvalidated_and_supported_capsid_hits_do_not_seed() -> None:
    assert _region_marker(CAPSID_MODEL, "unvalidated").is_valid_seed_marker is False
    assert _region_marker(CAPSID_MODEL, "supported").is_valid_seed_marker is False


def test_validated_novel_capsid_hit_alone_seeds_a_region_that_is_not_an_mcp_region(tmp_path: Path) -> None:
    genome = tmp_path / "genome.fna"
    genome.write_text(f">ctgA\n{'ACGT' * 5000}\n", encoding="utf-8")
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">ctgA_1 # 1000 # 1300 # 1 # ID=1_1;partial=00\nMKK\n", encoding="utf-8")

    def assemble(model: str) -> list:
        marker = MarkerHit(
            query_porf="ctgA_1",
            scaffold="ctgA",
            start=1000,
            end=1300,
            strand="+",
            hmm_target=model,
            hmm_score=120.0,
            hmm_evalue=1e-30,
            validation_status="validated_novel",
            top10_prefixes="",
            best_hit_target="",
            best_hit_pident=0.0,
            best_hit_bits=0.0,
            has_ncldv=0,
            has_mirus=0,
            has_plv=0,
            has_vp=0,
            has_viral=0,
        )
        return assemble_candidate_regions(
            [marker],
            genome,
            proteome,
            tmp_path / model,
            initial_window_bp=100,
            initial_window_genes=1,
            min_markers_initial=1,
            extension_kb=5,
            merge_distance=0,
        )

    capsid_regions = assemble(CAPSID_MODEL)

    assert len(capsid_regions) == 1
    assert capsid_regions[0].has_mcp is False
    assert assemble(NON_CAPSID_MODEL) == []


# --- Scoring ---


@pytest.mark.parametrize("evalue,expected_penalty", [(1e-8, 0.0), (0.001, 0.15), (0.1, 0.15)])
def test_small_region_exemption_requires_a_strong_catalogue_capsid(evalue: float, expected_penalty: float) -> None:
    result = _result(
        gene_count=1,
        hallmark_evidence=[
            HallmarkEvidence(
                query_porf="p1|aa1-20",
                gene="p1",
                scaffold="scaf",
                start=100,
                end=200,
                model="GVOGm0003",
                evalue=evalue,
                score=70.0,
                domain_score=60.0,
                query_start=1,
                query_end=20,
                validation_status="validated_novel",
                viral_top10_max_pident=0.0,
                source="retained_hmm",
            )
        ],
    )

    calculate_eve_confidence(result, crf_confidence=0.0)

    assert result.score_components["size_penalty"] == expected_penalty


@pytest.mark.parametrize("gene_count,expected_penalty", [(1, 0.15), (2, 0.08), (3, 0.04), (4, 0.0)])
def test_raw_capsid_flag_cannot_remove_the_size_penalty(gene_count: int, expected_penalty: float) -> None:
    result = _result(gene_count=gene_count, has_capsid_evidence=True, has_mcp=True)

    score = calculate_eve_confidence(result, crf_confidence=0.0)

    assert result.score_components["size_penalty"] == expected_penalty
    assert result.score_components["scores"]["marker"] == 0.0
    assert result.score_components["hallmark_floor"] == 0.0
    assert score < 0.2


# --- Acceptance gate ---


def _gate_result(**overrides) -> SimpleNamespace:
    base = {
        "confidence_tier": "HIGH",
        "start": 0,
        "end": 3000,
        "hallmark_count": 1,
        "has_mcp": False,
        "has_capsid_evidence": False,
        "hallmark_genes": [CAPSID_MODEL],
        "region_classification": "UNKNOWN",
        "likely_family": "MIXED",
    }
    base.update(overrides)
    return SimpleNamespace(**base)


@pytest.mark.parametrize("tier", ["HIGH", "MEDIUM", "LOW"])
@pytest.mark.parametrize("family", ["MIXED", "PPV"])
def test_gate_accepts_a_short_one_hallmark_region_with_capsid_evidence(tier: str, family: str) -> None:
    kept = evaluate_v2_quality_gate(
        _gate_result(confidence_tier=tier, likely_family=family, has_capsid_evidence=True)
    ).kept
    rejected = evaluate_v2_quality_gate(_gate_result(confidence_tier=tier, likely_family=family)).kept

    assert kept is True
    assert rejected is False


@pytest.mark.parametrize("family", ["NCLDV", "MIRUS"])
def test_gate_accepts_a_short_large_virus_region_with_capsid_evidence(family: str) -> None:
    short = {"end": 3000, "likely_family": family}

    assert evaluate_v2_quality_gate(_gate_result(**short, has_capsid_evidence=True)).kept is True
    assert evaluate_v2_quality_gate(_gate_result(**short)).kept is False


def test_peve_only_region_with_capsid_evidence_is_accepted_through_the_mixed_gate() -> None:
    result = _result(
        confidence_tier="MEDIUM",
        taxonomy_class="VIRAL_UNKNOWN",
        hallmark_count=1,
        hallmark_genes=[CAPSID_MODEL],
        has_capsid_evidence=True,
        gene_taxonomy_records=[
            {
                "porf_id": "g1",
                "top10_prefixes": ["EUK-pEVE"],
                "top10_pidents": [80.0],
                "is_flanking": False,
            }
        ],
    )

    result.likely_family = infer_likely_family(result)
    decision = evaluate_v2_quality_gate(result)

    assert result.likely_family == "MIXED"
    assert decision.kept is True
    assert decision.effective_class == "MIXED"
    assert decision.reason == "mixed_high_medium_pass"
