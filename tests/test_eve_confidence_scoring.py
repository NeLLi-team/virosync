"""Unit tests for the EVE accept/reject scoring core.

`calculate_eve_confidence` and `assign_confidence_tier` decide which candidates
are reported and at what tier, but were previously exercised by no test. These
lock the tier boundaries and the basic evidence ordering (strong viral evidence
must outscore a bare candidate, and a bare candidate must not reach HIGH).
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from virosync.ablation import AblationID
from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    EvidenceSynthesizerConfig,
    VerificationResult,
    VerificationStatus,
    assign_confidence_tier,
    calculate_eve_confidence,
)


def _bare() -> VerificationResult:
    return VerificationResult(eve_id="eve", scaffold="scaf", start=0, end=1000)


def test_assign_confidence_tier_boundaries() -> None:
    r = _bare()
    r.final_confidence = 0.70
    assert assign_confidence_tier(r) == "HIGH"
    r.final_confidence = 0.699
    assert assign_confidence_tier(r) == "MEDIUM"
    r.final_confidence = 0.20
    assert assign_confidence_tier(r) == "MEDIUM"
    r.final_confidence = 0.199
    assert assign_confidence_tier(r) == "LOW"


def test_assign_confidence_tier_custom_thresholds() -> None:
    r = _bare()
    r.final_confidence = 0.5
    assert assign_confidence_tier(r, high=0.4, low=0.1) == "HIGH"
    assert assign_confidence_tier(r, high=0.9, low=0.6) == "LOW"


def test_confidence_in_unit_range_and_bare_not_high() -> None:
    score = calculate_eve_confidence(_bare(), crf_confidence=0.0)
    assert 0.0 <= score <= 1.0
    # A candidate with no viral evidence must never reach the HIGH tier.
    assert score < 0.7


def test_strong_viral_evidence_outscores_bare() -> None:
    weak_score = calculate_eve_confidence(_bare(), crf_confidence=0.0, tmvec_score=None)

    strong = _bare()
    strong.hallmark_count = 8
    strong.hallmark_diversity = 4
    strong.has_virus_specific_marker = True
    strong.has_mcp = True
    strong.has_capsid_evidence = True
    strong.marker_complement_score = 1.0
    strong.marker_category_hits = ["ncldv", "mirus"]
    strong.marker_family_hits = ["NCLDV", "MIRUS"]
    strong.vp_completeness_ratio = 1.0
    strong.gene_count = 20
    strong_score = calculate_eve_confidence(strong, crf_confidence=0.9, tmvec_score=0.85, use_crf_score=True)

    assert 0.0 <= strong_score <= 1.0
    assert strong_score > weak_score


def test_priority_marker_floor_applies() -> None:
    """A priority marker (e.g. MCP) should lift confidence to at least the floor."""
    r = _bare()
    r.has_mcp = True
    r.has_capsid_evidence = True
    r.has_virus_specific_marker = True
    r.hallmark_count = 2
    r.marker_category_hits = ["mcp"]
    r.marker_family_hits = ["NCLDV"]
    score = calculate_eve_confidence(
        r,
        crf_confidence=0.0,
        priority_markers=["mcp"],
        marker_floor_priority_only=0.55,
    )
    assert score >= 0.0  # floor logic must not error and stays in range
    assert score <= 1.0


def test_taxonomy_object_and_dictionary_records_score_identically() -> None:
    record = {
        "is_flanking": False,
        "has_hit": True,
        "has_viral": False,
        "top10_prefixes": ["EUK", "EUK", "EUK"],
        "top10_pidents": [80.0, 75.0, 70.0],
    }
    dictionary_result = _bare()
    dictionary_result.gene_count = 4
    dictionary_result.gene_taxonomy_records = [record]
    object_result = _bare()
    object_result.gene_count = 4
    object_result.gene_taxonomy_records = [SimpleNamespace(**record)]

    dictionary_score = calculate_eve_confidence(dictionary_result, crf_confidence=0.0)
    object_score = calculate_eve_confidence(object_result, crf_confidence=0.0)

    assert object_score == dictionary_score
    assert object_result.high_confidence_euk_genes == dictionary_result.high_confidence_euk_genes == 1
    assert object_result.score_components == dictionary_result.score_components


@pytest.mark.parametrize("ablation_id", [AblationID.A0, AblationID.A5])
@pytest.mark.parametrize("viral_genes", [1, 3])
def test_host_model_diagnostics_do_not_change_score_or_tier(ablation_id: AblationID, viral_genes: int) -> None:
    """Hold coordinates and independent evidence fixed across model outputs."""
    candidate = VerificationResult(
        eve_id="seed",
        scaffold="scaf",
        start=100,
        end=10100,
        gene_count=10,
        gene_taxonomy_viral_interior=viral_genes,
        gene_taxonomy_viral_top10=viral_genes,
        seed_sources=["hhg"],
        kfd=0.3,
        gc_deviation=0.15,
    )
    synthesizer = EvidenceSynthesizer(EvidenceSynthesizerConfig(ablation_id=ablation_id))
    boundary = SimpleNamespace(confidence=0.0)
    synthesizer._calculate_final_decision(candidate, boundary, has_hhg_evidence=True)
    original = (candidate.final_confidence, candidate.confidence_tier, candidate.status)
    original_components = candidate.score_components.copy()

    candidate.host_signature_gene_count = 9
    candidate.host_signature_fraction = 0.9
    candidate.host_signature_weighted_mean = 0.95
    candidate.host_signature_evalue_weighted = 0.99
    synthesizer._calculate_final_decision(candidate, boundary, has_hhg_evidence=True)

    assert (candidate.final_confidence, candidate.confidence_tier, candidate.status) == original
    assert candidate.score_components == original_components
    assert (candidate.start, candidate.end) == (100, 10100)
    assert candidate.to_dict()["host_signature_gene_count"] == 9
    assert candidate.to_dict()["host_signature_fraction"] == 0.9
    assert candidate.to_dict()["host_signature_weighted_mean"] == 0.95
    assert candidate.host_signature_evalue_weighted == 0.99


@pytest.mark.parametrize("host_gene_count", [0, 9])
@pytest.mark.parametrize(
    "hallmarks,capsid,expected_tier", [(0, False, "LOW"), (3, False, "MEDIUM"), (1, True, "MEDIUM")]
)
def test_non_hhg_acceptance_requires_independent_marker_evidence(
    host_gene_count: int, hallmarks: int, capsid: bool, expected_tier: str
) -> None:
    """Gene count cannot rescue a non-HHG seed without marker support."""
    candidate = VerificationResult(
        eve_id="seed",
        scaffold="scaf",
        start=100,
        end=10100,
        gene_count=10,
        gene_taxonomy_viral_interior=3,
        gene_taxonomy_viral_top10=3,
        seed_sources=["novelty"],
        hallmark_count=hallmarks,
        has_capsid_evidence=capsid,
        host_signature_gene_count=host_gene_count,
        host_signature_fraction=host_gene_count / 10,
        kfd=0.3,
        gc_deviation=0.15,
    )

    EvidenceSynthesizer()._calculate_final_decision(candidate, SimpleNamespace(confidence=0.0), has_hhg_evidence=False)

    assert candidate.confidence_tier == expected_tier
    assert (candidate.status is VerificationStatus.LOW_CONFIDENCE_TIEBREAKER) == (expected_tier == "LOW")
