from __future__ import annotations

from types import SimpleNamespace

import pytest

from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    EvidenceSynthesizerConfig,
    VerificationResult,
    infer_likely_family,
)
from virosync.pipeline.phase3.output_generator import evaluate_v2_quality_gate


def _result(region_classification: str = "UNKNOWN") -> VerificationResult:
    return VerificationResult(
        eve_id="EVE_ctg_100-1180",
        scaffold="ctg",
        start=100,
        end=1180,
        region_classification=region_classification,
    )


def _cress_hit(**overrides: object) -> dict[str, object]:
    hit: dict[str, object] = {
        "porf_id": "ctg_1|aa3-330",
        "hallmark_gene": "VS000804",
        "hmm_score": 61.431,
        "validation_status": "validated",
        "top10_prefixes": "CRESS__,EUK__,EUK__",
        "top10_targets": "CRESS__reference,EUK__host_1,EUK__host_2",
        "top10_pidents": "45.5,40.9,37.8",
    }
    hit.update(overrides)
    return hit


def test_identity_qualified_cress_rep_survives_phase3_family_synthesis() -> None:
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    result = _result("CRESS")

    synthesizer._process_hallmark_hits(result, [_cress_hit()])

    assert result.region_classification == "CRESS"
    assert result.marker_family_hits == ["CRESS"]
    assert result.marker_dominant_family == "CRESS"
    assert infer_likely_family(result) == "CRESS"


def test_cress_capsid_is_not_reported_as_a_large_virus_mcp(
    tmp_path,
) -> None:
    annotations = tmp_path / "annotations.tsv"
    annotations.write_text(
        "model_name\tsource\tdescription\tmajority_annotation\n"
        "VS000798\tCRESS_ssDNA_Oliver2026\t"
        "CRESS ssDNA capsid profile\tCRESS ssDNA virus capsid protein\n"
    )
    synthesizer = EvidenceSynthesizer(
        config=EvidenceSynthesizerConfig(
            marker_annotations_path=annotations,
        )
    )
    result = _result("CRESS")

    synthesizer._process_hallmark_hits(
        result,
        [_cress_hit(hallmark_gene="VS000798")],
    )

    assert "capsid" in result.marker_category_hits
    assert result.has_mcp is False


@pytest.mark.parametrize(
    "overrides",
    [
        {"top10_pidents": "24.9,40.9,37.8"},
        {"top10_pidents": "24.96,40.9,37.8"},
        {"top10_prefixes": "EUK__,PHAGE__", "top10_pidents": "45.5,40.9"},
        {"validation_status": "validated_novel"},
        {"hallmark_gene": "VS000791"},
    ],
)
def test_cress_family_requires_cress_model_and_identity_qualified_taxonomy(
    overrides: dict[str, object],
) -> None:
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    result = _result()

    synthesizer._process_hallmark_hits(result, [_cress_hit(**overrides)])

    assert result.region_classification == "UNKNOWN"
    assert "CRESS" not in result.marker_family_hits


@pytest.mark.parametrize(
    ("tier", "expected_reason"),
    [
        ("HIGH", "cress_identity_high_medium_pass"),
        ("MEDIUM", "cress_identity_high_medium_pass"),
        ("LOW", "cress_identity_low_promoted"),
    ],
)
def test_single_identity_qualified_cress_gene_passes_without_two_kb_floor(
    tier: str,
    expected_reason: str,
) -> None:
    decision = evaluate_v2_quality_gate(
        SimpleNamespace(
            confidence_tier=tier,
            start=100,
            end=1180,
            hallmark_count=1,
            hallmark_genes=["VS000804"],
            marker_family_hits=["CRESS"],
            has_mcp=False,
            region_classification="CRESS",
            classification="CRESS",
            likely_family="CRESS",
        )
    )

    assert decision.kept
    assert decision.effective_class == "CRESS"
    assert decision.reason == expected_reason
    assert decision.promoted_low is (tier == "LOW")


def test_euk_top_hit_single_cress_gene_remains_discovery_only() -> None:
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    result = _result("CRESS")

    synthesizer._process_hallmark_hits(
        result,
        [
            _cress_hit(
                top10_prefixes="EUK__,CRESS__",
                top10_targets="EUK__host,CRESS__reference",
                top10_pidents="98.4,27.7",
            )
        ],
    )
    result.confidence_tier = "MEDIUM"

    assert result.region_classification == "CRESS"
    assert "CRESS" not in result.marker_family_hits
    assert not evaluate_v2_quality_gate(result).kept


def test_two_qualified_cress_genes_pass_with_euk_top_hits() -> None:
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    result = _result("CRESS")
    weak_hit = {
        "top10_prefixes": "EUK__,CRESS__",
        "top10_targets": "EUK__host,CRESS__reference",
        "top10_pidents": "70.0,30.0",
    }

    synthesizer._process_hallmark_hits(
        result,
        [
            _cress_hit(porf_id="ctg_1|aa1-100", **weak_hit),
            _cress_hit(porf_id="ctg_2|aa1-100", **weak_hit),
        ],
    )
    result.confidence_tier = "MEDIUM"

    assert result.marker_family_hits == ["CRESS"]
    assert evaluate_v2_quality_gate(result).kept


def test_generic_boundary_overlapping_cress_hit_is_not_promoted() -> None:
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    result = _result("UNKNOWN")

    synthesizer._process_hallmark_hits(result, [_cress_hit()])
    result.confidence_tier = "MEDIUM"
    result.likely_family = "CRESS"

    assert result.region_classification == "UNKNOWN"
    assert "CRESS" not in result.marker_family_hits
    decision = evaluate_v2_quality_gate(result)
    assert not decision.kept
    assert decision.reason == "cress_identity_required"


@pytest.mark.parametrize("tier", ["HIGH", "LOW"])
def test_unqualified_single_cress_label_remains_gated_out(tier: str) -> None:
    decision = evaluate_v2_quality_gate(
        SimpleNamespace(
            confidence_tier=tier,
            start=100,
            end=1180,
            hallmark_count=1,
            hallmark_genes=["VS000804"],
            marker_family_hits=[],
            has_mcp=False,
            region_classification="UNKNOWN",
            classification="",
            likely_family="CRESS",
        )
    )

    assert not decision.kept
    assert decision.reason == "cress_identity_required"
