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
        "top10_targets": ["EUK__a", "EUK__b", "EUK__c"],
        "top10_bitscores": [100.0, 90.0, 80.0],
        "top10_evalues": [1e-8] * 3,
        "top10_qcovhsp": [80.0] * 3,
        "porf_id": "host",
        "scaffold": "scaf",
        "start": 100,
        "end": 200,
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


def _supplementary_gene(name, start, end, *, group="NCLDV__", qcov=80.0, viral_bits=120.0, cellular_bits=100.0):
    from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy, build_gene_taxonomy_record

    targets = [f"{group}reference|{name}"]
    prefixes, identities, evalues, bits, coverage = [group], [30.0], [1e-8], [viral_bits], [qcov]
    if cellular_bits is not None:
        targets.append(f"EUK__host|{name}")
        prefixes.append("EUK__")
        identities.append(10.0)
        evalues.append(0.1)
        bits.append(cellular_bits)
        coverage.append(10.0)
    return build_gene_taxonomy_record(
        GeneTaxonomy(
            name,
            "scaf",
            start,
            end,
            top10_targets=targets,
            top10_prefixes=prefixes,
            top10_pidents=identities,
            top10_evalues=evalues,
            top10_bits=bits,
            top10_qcovhsp=coverage,
        )
    )


def _supplementary_candidate(*, end=6000, total=3, records=None):
    return VerificationResult(
        "module",
        "scaf",
        0,
        end,
        likely_family="NCLDV",
        region_classification="NCLDV",
        gene_count=total,
        gene_taxonomy_total=total,
        gene_taxonomy_viral_top10=3,
        gene_taxonomy_records=records
        if records is not None
        else [
            _supplementary_gene("a", 100, 200),
            _supplementary_gene("b", 300, 400),
            _supplementary_gene("c", 500, 600),
        ],
    )


def test_supplementary_module_removes_only_zero_core_cap_without_changing_weighted_evidence():
    eligible = _supplementary_candidate()
    baseline = _supplementary_candidate(
        records=[
            _supplementary_gene("a", 100, 200, qcov=49.9),
            _supplementary_gene("b", 300, 400),
            _supplementary_gene("c", 500, 600),
        ]
    )
    eligible_score = calculate_eve_confidence(eligible, 0.0)
    baseline_score = calculate_eve_confidence(baseline, 0.0)
    assert eligible_score > 0.2
    assert baseline_score == pytest.approx(0.199)
    assert eligible.supplementary_viral_evidence.eligible
    assert not baseline.supplementary_viral_evidence.eligible
    assert eligible.score_components["scores"] == baseline.score_components["scores"]
    assert eligible.score_components["weighted_base"] == baseline.score_components["weighted_base"]
    assert eligible.score_components["bonus_total"] == baseline.score_components["bonus_total"]
    assert eligible.to_dict()["supplementary_viral_evidence"]["counted_gene_ids"] == ("a", "b", "c")


def test_supplementary_zero_core_still_obeys_high_cap():
    result = _supplementary_candidate()
    result.structural_score = 1.0
    result.interproscan_score = 1.0
    result.kfd = 1.0
    result.gc_deviation = 1.0
    result.region_classification_ncldv_markers = 10
    assert calculate_eve_confidence(result, 1.0, use_crf_score=True) == pytest.approx(0.699)


@pytest.mark.parametrize("end,kept", [(6000, True), (5000, False)])
def test_supplementary_module_low_gate_uses_same_record_and_keeps_length_rule(end, kept):
    from virosync.pipeline.phase3.output_generator import evaluate_v2_quality_gate

    result = _supplementary_candidate(end=end, total=100)
    result.final_confidence = calculate_eve_confidence(result, 0.0)
    result.confidence_tier = assign_confidence_tier(result)
    assert result.confidence_tier == "LOW"
    assessment = result.supplementary_viral_evidence
    decision = evaluate_v2_quality_gate(result)
    assert decision.kept is kept
    assert result.supplementary_viral_evidence is assessment
    result.end += 1
    with pytest.raises(ValueError, match="supplementary.*interval"):
        evaluate_v2_quality_gate(result)


@pytest.mark.parametrize(
    "qcov,cellular_bits,eligible", [(50.0, 100, True), (49.9, 100, False), (80.0, 120, False), (80.0, 130, False)]
)
def test_supplementary_coverage_and_unfiltered_cellular_competition(qcov, cellular_bits, eligible):
    result = _supplementary_candidate(
        records=[
            _supplementary_gene("a", 100, 200, qcov=qcov, cellular_bits=cellular_bits),
            _supplementary_gene("b", 300, 400),
            _supplementary_gene("c", 500, 600),
        ]
    )
    calculate_eve_confidence(result, 0.0)
    assert result.supplementary_viral_evidence.eligible is eligible


def test_supplementary_missing_coverage_is_error_but_empty_no_hit_record_is_valid():
    result = _supplementary_candidate()
    del result.gene_taxonomy_records[0]["top10_qcovhsp"]
    with pytest.raises(ValueError, match="query coverage"):
        calculate_eve_confidence(result, 0.0)
    result.gene_taxonomy_records = [{"porf_id": "empty", "scaffold": "scaf", "start": 100, "end": 200}]
    calculate_eve_confidence(result, 0.0)
    assert not result.supplementary_viral_evidence.eligible


@pytest.mark.parametrize("start,end,eligible", [(500, 6000, True), (500, 6001, False), (-1, 50, False)])
def test_supplementary_requires_full_containment_without_buffer(start, end, eligible):
    result = _supplementary_candidate(
        records=[
            _supplementary_gene("a", 100, 200),
            _supplementary_gene("b", 300, 400),
            _supplementary_gene("c", start, end),
        ]
    )
    calculate_eve_confidence(result, 0.0)
    assert result.supplementary_viral_evidence.eligible is eligible


def test_supplementary_deduplicates_physical_genes_and_excludes_opposite_strand_overlap():
    duplicate = _supplementary_gene("alias", 100, 200)
    overlap = _supplementary_gene("opposite", 399, 600)
    overlap["strand"] = "-"
    result = _supplementary_candidate(
        records=[
            _supplementary_gene("a", 100, 200),
            duplicate,
            _supplementary_gene("b", 300, 400),
            overlap,
        ]
    )
    calculate_eve_confidence(result, 0.0)
    assert result.supplementary_viral_evidence.counted_gene_ids == ("a", "b")
    assert not result.supplementary_viral_evidence.eligible


def test_supplementary_measured_only_selection_is_independent_and_abutting_genes_count():
    result = _supplementary_candidate(
        records=[
            _supplementary_gene("a", 0, 10),
            _supplementary_gene("no_cellular", 10, 20, cellular_bits=None),
            _supplementary_gene("b", 10, 30),
            _supplementary_gene("c", 30, 40),
        ]
    )
    calculate_eve_confidence(result, 0.0)
    evidence = result.supplementary_viral_evidence
    assert evidence.counted_gene_ids == ("a", "no_cellular", "c")
    assert evidence.no_cellular_gene_ids == ("no_cellular",)
    assert evidence.measured_only_gene_ids == ("a", "b", "c")
    assert evidence.eligible_without_no_cellular


@pytest.mark.parametrize(
    "region,family,groups",
    [
        ("MIRUS", "NCLDV", ("NCLDV__",) * 3),
        ("NCLDV", "NCLDV", ("NCLDV__", "NCLDV__", "MIRUS__")),
        ("UNKNOWN", "UNKNOWN", ("NCLDV__",) * 3),
    ],
)
def test_supplementary_does_not_reclassify_or_sum_groups(region, family, groups):
    from virosync.pipeline.phase3.output_generator import evaluate_v2_quality_gate

    result = _supplementary_candidate(
        records=[
            _supplementary_gene("a", 100, 200, group=groups[0]),
            _supplementary_gene("b", 300, 400, group=groups[1]),
            _supplementary_gene("c", 500, 600, group=groups[2]),
        ]
    )
    result.region_classification, result.likely_family = region, family
    result.final_confidence = calculate_eve_confidence(result, 0.0)
    result.confidence_tier = assign_confidence_tier(result)
    assert not result.supplementary_viral_evidence.eligible
    assert not evaluate_v2_quality_gate(result).kept
    assert (result.region_classification, result.likely_family) == (region, family)


@pytest.mark.parametrize("identity,evalue,eligible", [(25, 1e-5, True), (24.99, 1e-5, False), (25, 1.01e-5, False)])
def test_supplementary_viral_identity_and_evalue_boundaries(identity, evalue, eligible):
    result = _supplementary_candidate()
    result.gene_taxonomy_records[0]["top10_pidents"][0] = identity
    result.gene_taxonomy_records[0]["top10_evalues"][0] = evalue
    calculate_eve_confidence(result, 0.0)
    assert result.supplementary_viral_evidence.eligible is eligible


def test_supplementary_preserves_existing_capsid_floor_and_original_status():
    result = _supplementary_candidate(total=100)
    result.hallmark_evidence = [
        HallmarkEvidence(
            query_porf="capsid|aa1-20",
            gene="capsid",
            scaffold="scaf",
            start=100,
            end=200,
            model="GVOGm0003",
            evalue=0.001,
            score=70,
            domain_score=60,
            query_start=1,
            query_end=20,
            validation_status="unvalidated",
            viral_top10_max_pident=0,
            source="retained_hmm",
        )
    ]
    assert calculate_eve_confidence(result, 0.0) == pytest.approx(0.2)
    assert result.supplementary_viral_evidence.eligible
    assert result.hallmark_evidence[0].validation_status == "unvalidated"


def test_supplementary_alias_cellular_competitor_blocks_duplicate_viral_record():
    result = _supplementary_candidate()
    alias = _supplementary_gene("alias", 100, 200, cellular_bits=200)
    result.gene_taxonomy_records.append(alias)
    calculate_eve_confidence(result, 0.0)
    assert not result.supplementary_viral_evidence.eligible
    assert result.supplementary_viral_evidence.counted_gene_ids == ("b", "c")


@pytest.mark.parametrize("total,tier", [(100, "LOW"), (3, "MEDIUM")])
def test_supplementary_native_verify_infers_existing_family_before_assessment(total, tier):
    from virosync.pipeline.phase3.acceptance_selection import select_phase3_acceptance
    from virosync.pipeline.phase3.eve_ani_clustering import unsupported_eve_ids

    boundary = RefinedBoundary(
        "scaf",
        0,
        6000,
        predicted_family="NCLDV",
        confirmed_seeds=(ConfirmedSeed("core", "scaf", 100, 200),),
    )
    records = _supplementary_candidate().gene_taxonomy_records
    result = EvidenceSynthesizer(EvidenceSynthesizerConfig(use_tmvec_database=False)).verify_eve(
        boundary,
        [],
        gene_taxonomy_records=records,
        gene_taxonomy_summary={"total": total, "viral_top10": 3, "ncldv_mirus": 3},
    )
    assert result.likely_family == "NCLDV"
    assert result.supplementary_viral_evidence.eligible
    assert result.supplementary_viral_evidence.low_class == "NCLDV"
    assert result.confidence_tier == tier
    assert select_phase3_acceptance([result], AblationID.A0).canonical_results == (result,)
    assert unsupported_eve_ids([result]) == set()


def _native_organelle_records(target_prefix, competitor_bits, weak_euk):
    from virosync.pipeline.phase2.boundary_diamond import (
        DiamondHit,
        build_gene_taxonomy_record,
        classify_all_porfs,
        extract_prefix,
        pORF,
    )

    genes = [pORF(f"native{i}", "scaf", i * 100, i * 100 + 90) for i in range(1, 4)]
    hits = {}
    for gene in genes:
        target = f"{target_prefix}reference|{gene.id}"
        assert extract_prefix(target) == "UNKNOWN"
        rows = [
            DiamondHit(gene.id, f"NCLDV__viral|{gene.id}", 1e-8, 100, 30, 80),
            DiamondHit(gene.id, target, 0.5, competitor_bits, 5, 2),
        ]
        if weak_euk:
            rows.append(DiamondHit(gene.id, f"EUK__host|{gene.id}", 0.1, 50, 10, 5))
        hits[gene.id] = sorted(rows, key=lambda row: row.bits, reverse=True)
    taxonomy = classify_all_porfs([gene.id for gene in genes], hits, {"scaf": genes}, host_prefix="EUK__")
    records = [build_gene_taxonomy_record(taxonomy[gene.id]) for gene in genes]
    assert all("UNKNOWN" in record["top10_prefixes"] for record in records)
    return records


@pytest.mark.parametrize("prefix", ["MITO__", "PLASTID__"])
@pytest.mark.parametrize("bits", [120, 100, 80])
@pytest.mark.parametrize("weak_euk", [False, True])
def test_native_unknown_organelle_competes_without_viral_filters(prefix, bits, weak_euk):
    records = _native_organelle_records(prefix, bits, weak_euk)
    candidate = _supplementary_candidate(records=records)
    calculate_eve_confidence(candidate, 0.0)
    evidence = candidate.supplementary_viral_evidence
    assert evidence.eligible is (bits < 100)
    assert len(evidence.counted_gene_ids) == (3 if bits < 100 else 0)
    assert len(evidence.measured_gene_ids) == (3 if bits < 100 else 0)
    assert evidence.no_cellular_gene_ids == ()
    assert evidence.eligible_without_no_cellular is (bits < 100)


@pytest.mark.parametrize("prefix", ["MITOLOGY__", "PLASTIDLIKE__", "XMITO__"])
def test_unrelated_unknown_target_is_not_an_organelle_competitor(prefix):
    records = _native_organelle_records(prefix, 200, False)
    candidate = _supplementary_candidate(records=records)
    calculate_eve_confidence(candidate, 0.0)
    evidence = candidate.supplementary_viral_evidence
    assert evidence.eligible
    assert len(evidence.no_cellular_gene_ids) == 3
    assert not evidence.eligible_without_no_cellular
