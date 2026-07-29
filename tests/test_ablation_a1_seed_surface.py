from __future__ import annotations

import pytest

from virosync.ablation import AblationID, InterventionCounts
from virosync.pipeline.phase1.hhg_seeding import Anchor
from virosync.pipeline.phase1.seed_merger import MergedSeed
from virosync.pipeline.phase3.evidence_synthesizer import VerificationStatus
from virosync.pipeline.phase3.phase1_surface import build_phase1_seed_surface
from virosync.orchestration._flows.single_genome.orchestrator import (
    _write_combined_eve_fasta,
)


def _anchor(
    porf_id: str,
    hallmark_gene: str,
    score: float,
    *,
    start: int = 10,
    end: int = 40,
    scaffold: str = "ctg",
    evalue: float = 1e-20,
) -> Anchor:
    return Anchor(
        porf_id=porf_id,
        scaffold=scaffold,
        start=start,
        end=end,
        strand="+",
        hallmark_gene=hallmark_gene,
        score=score,
        evalue=evalue,
    )


def test_a1_exports_exact_phase1_coordinates_and_unscored_provenance() -> None:
    lower_scoring_alias = _anchor("gene-a|aa1-90", "GVOGm0003", 30.0)
    selected_gene_a = _anchor("gene-a", "GVOGm0004", 50.0, evalue=1e-30)
    mcp = _anchor("gene-b", "plv_mcp_1", 80.0, start=60, end=95)
    outside_seed = _anchor(
        "gene-c",
        "mirus_mcp_1",
        90.0,
        start=110,
        end=140,
    )
    other_scaffold = _anchor(
        "gene-d",
        "vp_mcp_1",
        100.0,
        scaffold="other",
    )
    seed = MergedSeed(
        scaffold="ctg",
        start=0,
        end=100,
        seed_id="order-dependent-seed-id",
        sources=["marker_validation", "hhg", "hhg"],
        hhg_score=20.0,
        novelty_score=0.25,
        compositional_score=0.5,
        max_kfd=0.125,
        gc_deviation=0.2,
        cub_deviation=0.0,
        confidence="high",
        predicted_family="PLV",
        region_classification_ncldv_markers=1,
        region_classification_vp_plv_markers=2,
        region_classification_mirus_markers=0,
        anchors=[lower_scoring_alias, selected_gene_a, mcp, outside_seed, other_scaffold],
        hhg_anchors=[selected_gene_a, mcp],
    )

    surface = build_phase1_seed_surface([seed])
    result = surface.results[0]

    assert surface.detailed_results is surface.results
    assert surface.canonical_results is surface.results
    assert surface.intervention_counts == InterventionCounts(1, 1, 0)
    assert result.eve_id == "EVE_ctg_0-100"
    assert (result.scaffold, result.start, result.end, result.length) == (
        "ctg",
        0,
        100,
        100,
    )
    assert result.ablation_id is AblationID.A1
    assert result.status is VerificationStatus.AMBIGUOUS
    assert result.final_confidence == 0.0
    assert result.confidence_tier == "LOW"
    assert result.score_components == {
        "prediction_stage": "phase1_seed_surface",
        "confidence_kind": "not_scored",
    }

    assert result.hallmark_genes == ["GVOGm0004", "plv_mcp_1"]
    assert result.hallmark_count == 2
    assert result.hallmark_diversity == 2
    assert result.has_virus_specific_marker is True
    assert result.has_mcp is True
    assert result.mcp_gene_ids == ["gene-b"]
    # plv_ markers report under the unified Preplasmiviricota family.
    assert result.marker_family_hits == ["NCLDV", "PPV"]

    assert result.seed_sources == ["hhg", "marker_validation"]
    assert result.seed_confidence == "high"
    assert result.seed_hhg_score == 20.0
    assert result.seed_novelty_score == 0.25
    assert result.seed_compositional_score == 0.5
    # plv_ seeds surface under the unified Preplasmiviricota class.
    assert result.region_classification == "PPV"
    assert result.likely_family == "PPV"
    assert result.region_classification_ncldv_markers == 1
    assert result.region_classification_vp_plv_markers == 2
    assert result.region_classification_mirus_markers == 0
    assert result.kfd == 0.125
    assert result.gc_deviation == 0.2
    assert result.cub_deviation == 0.0


def test_a1_identity_and_hallmarks_are_stable_across_input_order() -> None:
    first = _anchor("gene-a", "GVOGm0004", 50.0)
    second = _anchor("gene-b", "plv_mcp_1", 80.0, start=60, end=95)
    seed_a = MergedSeed(
        scaffold="ctg",
        start=5,
        end=100,
        seed_id="seed_0_ctg_5",
        anchors=[second, first],
        hhg_anchors=[first],
    )
    seed_b = MergedSeed(
        scaffold="ctg",
        start=5,
        end=100,
        seed_id="seed_99_ctg_5",
        anchors=[first, second],
        hhg_anchors=[second],
    )

    result_a = build_phase1_seed_surface([seed_a]).results[0]
    result_b = build_phase1_seed_surface([seed_b]).results[0]

    assert result_a.eve_id == result_b.eve_id == "EVE_ctg_5-100"
    assert result_a.hallmark_genes == result_b.hallmark_genes == [
        "GVOGm0004",
        "plv_mcp_1",
    ]
    assert result_a.mcp_gene_ids == result_b.mcp_gene_ids == ["gene-b"]


def test_a1_empty_surface_has_zero_intervention_counts() -> None:
    surface = build_phase1_seed_surface([])

    assert surface.results == ()
    assert surface.canonical_results == ()
    assert surface.intervention_counts == InterventionCounts()


def test_a1_combined_fasta_is_available_before_report_generation(tmp_path) -> None:
    genome_path = tmp_path / "masked.fna"
    genome_path.write_text(">ctg\nAACCGGTTAACC\n", encoding="utf-8")
    seed = MergedSeed(
        scaffold="ctg",
        start=2,
        end=8,
        seed_id="seed_0_ctg_2",
        anchors=[_anchor("gene-a", "GVOGm0004", 50.0, start=2, end=8)],
    )
    surface = build_phase1_seed_surface([seed])

    output_path = _write_combined_eve_fasta(
        output_dir=tmp_path,
        genome_id="test-genome",
        genome_path=genome_path,
        results=surface.canonical_results,
    )

    assert output_path == tmp_path / "test-genome_eves.fna"
    assert output_path.read_text(encoding="utf-8").splitlines() == [
        ">EVE_ctg_2-8 scaffold=ctg start=2 end=8 tier=LOW confidence=0.0000",
        "CCGGTT",
    ]


def test_combined_fasta_is_not_created_for_an_empty_surface(tmp_path) -> None:
    genome_path = tmp_path / "masked.fna"
    genome_path.write_text(">ctg\nAACCGGTT\n", encoding="utf-8")

    output_path = _write_combined_eve_fasta(
        output_dir=tmp_path,
        genome_id="test-genome",
        genome_path=genome_path,
        results=(),
    )

    assert output_path is None
    assert not (tmp_path / "test-genome_eves.fna").exists()


def test_combined_fasta_fails_closed_when_results_do_not_match_the_genome(
    tmp_path,
) -> None:
    genome_path = tmp_path / "masked.fna"
    genome_path.write_text(">other\nAACCGGTT\n", encoding="utf-8")
    seed = MergedSeed(
        scaffold="ctg",
        start=2,
        end=8,
        seed_id="seed_0_ctg_2",
        anchors=[_anchor("gene-a", "GVOGm0004", 50.0, start=2, end=8)],
    )
    surface = build_phase1_seed_surface([seed])

    with pytest.raises(RuntimeError, match="produced no file"):
        _write_combined_eve_fasta(
            output_dir=tmp_path,
            genome_id="test-genome",
            genome_path=genome_path,
            results=surface.canonical_results,
        )
