"""Unit tests for the EVE accept/reject scoring core.

`calculate_eve_confidence` and `assign_confidence_tier` decide which candidates
are reported and at what tier, but were previously exercised by no test. These
lock the tier boundaries and the basic evidence ordering (strong viral evidence
must outscore a bare candidate, and a bare candidate must not reach HIGH).
"""

from __future__ import annotations

import logging
from types import SimpleNamespace

import pytest

from virosync.ablation import AblationID
from virosync.config import PipelineConfig
from virosync.orchestration._flows.single_genome.phase3 import _summarize_boundary_taxonomy
from virosync.pipeline.hallmark_evidence import HallmarkEvidence
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed
from virosync.pipeline.phase2.boundary_diamond import DiamondHit, classify_all_porfs, pORF
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    EvidenceSynthesizerConfig,
    VerificationResult,
    assign_confidence_tier,
    calculate_eve_confidence,
)


def _bare() -> VerificationResult:
    return VerificationResult(eve_id="eve", scaffold="scaf", start=0, end=1000)


@pytest.mark.parametrize("pident,qualified_count,taxonomy_score", [(24.999, 0, 0.0), (25.0, 3, 1.0)])
@pytest.mark.parametrize("seed_marker_count,seed_score", [(0, 0.0), (3, 1.0)])
def test_native_boundary_taxonomy_scores_only_identity_qualified_viral_genes(
    pident: float, qualified_count: int, taxonomy_score: float, seed_marker_count: int, seed_score: float
) -> None:
    genes = [pORF("g1", "scaf", 100, 200), pORF("g2", "scaf", 300, 400), pORF("g3", "scaf", 500, 600)]
    hits = {gene.id: [DiamondHit(gene.id, "NCLDV__reference|protein", 1e-5, 80, pident, 90)] for gene in genes}
    taxonomy = classify_all_porfs([gene.id for gene in genes], hits, {"scaf": genes}, host_prefix="EUK__")
    boundary = RefinedBoundary(
        "scaf",
        0,
        1000,
        confidence=0,
        region_classification_ncldv_markers=seed_marker_count,
        confirmed_seeds=(ConfirmedSeed("seed", "scaf", 0, 1000),),
    )
    records, summary = _summarize_boundary_taxonomy(
        boundary=boundary,
        filtered_taxonomy=list(taxonomy.values()),
        boundary_taxonomy_map=taxonomy,
        proteome_index={"scaf": genes},
        boundary_diamond_query=None,
        config=PipelineConfig(),
        logger=logging.getLogger(__name__),
    )
    result = EvidenceSynthesizer(EvidenceSynthesizerConfig(use_tmvec_database=False)).verify_eve(
        boundary,
        [],
        gene_taxonomy_records=records,
        gene_taxonomy_summary=summary,
    )
    assert summary["ncldv_mirus"] == result.genes_with_ncldv_mirus_top10 == 3
    assert summary["viral_top10"] == result.gene_taxonomy_viral_top10 == qualified_count
    assert result.score_components["scores"]["seed_cluster"] == seed_score
    assert result.score_components["scores"]["gene_taxonomy"] == taxonomy_score
    if qualified_count == 0:
        assert result.final_confidence == pytest.approx(0.085 if seed_marker_count else 0.0)


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
    strong.hallmark_evidence = [
        HallmarkEvidence(
            query_porf=f"{gene}|aa1-20",
            gene=gene,
            scaffold="scaf",
            start=start,
            end=start + 100,
            model=model,
            evalue=1e-8,
            score=70.0,
            domain_score=60.0,
            query_start=1,
            query_end=20,
            validation_status="validated_novel",
            viral_top10_max_pident=0.0,
            source="retained_hmm",
        )
        for gene, model, start in [("capsid", "GVOGm0003", 100), ("atpase", "GVOGm0760", 300)]
    ]
    strong.gene_count = 20
    strong_score = calculate_eve_confidence(strong, crf_confidence=0.9, tmvec_score=0.85, use_crf_score=True)

    assert 0.0 <= strong_score <= 1.0
    assert strong_score >= 0.7
    assert strong_score > weak_score


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
