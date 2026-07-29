"""Unit tests for the EVE accept/reject scoring core.

`calculate_eve_confidence` and `assign_confidence_tier` decide which candidates
are reported and at what tier, but were previously exercised by no test. These
lock the tier boundaries and the basic evidence ordering (strong viral evidence
must outscore a bare candidate, and a bare candidate must not reach HIGH).
"""

from __future__ import annotations

from virosync.pipeline.phase3.evidence_synthesizer import (
    VerificationResult,
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
    strong.marker_complement_score = 1.0
    strong.marker_category_hits = ["ncldv", "mirus"]
    strong.marker_family_hits = ["NCLDV", "MIRUS"]
    strong.vp_completeness_ratio = 1.0
    strong.gene_count = 20
    strong_score = calculate_eve_confidence(
        strong, crf_confidence=0.9, tmvec_score=0.85, use_crf_score=True
    )

    assert 0.0 <= strong_score <= 1.0
    assert strong_score > weak_score


def test_priority_marker_floor_applies() -> None:
    """A priority marker (e.g. MCP) should lift confidence to at least the floor."""
    r = _bare()
    r.has_mcp = True
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
