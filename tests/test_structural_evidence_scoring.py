from __future__ import annotations

from pathlib import Path

from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    EvidenceSynthesizerConfig,
    VerificationResult,
)
from virosync.pipeline.phase3.structural_homology import StructuralHomologyResult
from virosync.pipeline.phase3.tmvec_database import TMVecHit


def _result() -> VerificationResult:
    return VerificationResult(
        eve_id="EVE_scaffold_1-100",
        scaffold="scaffold",
        start=1,
        end=100,
    )


def test_tmvec_structural_support_requires_configured_score_threshold(
    tmp_path: Path,
) -> None:
    synthesizer = EvidenceSynthesizer(
        config=EvidenceSynthesizerConfig(
            use_tmvec_database=True,
            tmvec_min_score=0.5,
        ),
        work_dir=tmp_path,
    )
    result = _result()

    synthesizer._run_tmvec_database_scan(
        result,
        [("p1", "M" * 50)],
        precomputed_tmvec={
            "p1": {
                "bfvd": TMVecHit(
                    target_id="bfvd_low",
                    tm_score=0.49,
                    database="bfvd",
                ),
            },
        },
    )

    assert result.structural_score == 0.0
    assert result.has_structural_support is False
    assert result.tmvec_all_proteins[0]["eve_id"] == result.eve_id
    assert result.tmvec_all_proteins[0]["tmvec_bfvd_score"] == 0.49
    assert not {
        "tmvec_cath_score",
        "tmvec_swiss_score",
        "tmvec_viral_specificity",
    } & result.tmvec_all_proteins[0].keys()

    synthesizer._run_tmvec_database_scan(
        result,
        [("p1", "M" * 50)],
        precomputed_tmvec={
            "p1": {
                "bfvd": TMVecHit(
                    target_id="bfvd_hit",
                    tm_score=0.50,
                    database="bfvd",
                ),
            },
        },
    )

    assert result.structural_score == 0.50
    assert result.has_structural_support is True


def test_boltz_results_do_not_overwrite_existing_tmvec_support(
    tmp_path: Path,
) -> None:
    class FakeBoltzAnalyzer:
        def analyze_batch(self, sequences, work_dir):
            return [
                StructuralHomologyResult(
                    porf_id="p1",
                    prediction=None,
                    structural_evidence_score=0.2,
                )
            ]

    synthesizer = EvidenceSynthesizer(
        config=EvidenceSynthesizerConfig(
            use_boltz=True,
            use_tmvec_database=False,
            use_phylogenetic_validation=False,
        ),
        work_dir=tmp_path,
    )
    synthesizer._boltz_analyzer = FakeBoltzAnalyzer()
    result = _result()
    result.structural_score = 0.8
    result.has_structural_support = True
    boundary = RefinedBoundary(scaffold="scaffold", start=1, end=100)

    synthesizer._run_tiebreakers(
        result=result,
        refined_boundary=boundary,
        window_features=[],
        hallmark_hits=[
            {
                "hallmark_gene": "mcp",
                "porf_id": "p1",
            }
        ],
        novelty_scores=None,
        porf_sequences=[("p1", "M" * 120)],
    )

    assert result.structural_score == 0.8
    assert result.has_structural_support is True
