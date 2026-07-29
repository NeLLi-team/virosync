from __future__ import annotations

from pathlib import Path

from virosync.ablation import AblationID
from virosync.orchestration._flows.single_genome.phase1 import (
    _region_coordinate_surface,
)
from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit
from virosync.pipeline.phase1.region_assembly import assemble_candidate_regions
from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    EvidenceSynthesizerConfig,
    VerificationResult,
)


def _strong_supported_marker() -> ValidatedMarkerHit:
    return ValidatedMarkerHit(
        query_porf="ctg_1",
        scaffold="ctg",
        start=99,
        end=300,
        strand="+",
        hmm_target="GVOGm0003",
        hmm_score=100.0,
        hmm_evalue=1e-20,
        validation_status="supported",
        top10_prefixes="EUK__",
        best_hit_target="EUK__host",
        best_hit_pident=80.0,
        best_hit_bits=200.0,
        has_ncldv=0,
        has_mirus=0,
        has_plv=0,
        has_vp=0,
        has_viral=0,
    )


def _inputs(tmp_path: Path) -> tuple[Path, Path]:
    genome = tmp_path / "genome.fna"
    genome.write_text(">ctg\n" + "ACGT" * 250 + "\n")
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">ctg_1 # 100 # 300 # 1 # ID=1_1;partial=00\n"
        "M" * 67
        + "\n"
    )
    return genome, proteome


def test_a2_strong_supported_marker_creates_explicit_bypass_seed(
    tmp_path: Path,
) -> None:
    genome, proteome = _inputs(tmp_path)
    marker = _strong_supported_marker()

    normal = assemble_candidate_regions(
        [marker],
        genome,
        proteome,
        tmp_path / "a0",
        min_markers_initial=1,
        extension_kb=0,
        merge_distance=0,
        ablation_id=AblationID.A0,
        single_marker_min_score=100.0,
    )
    bypassed = assemble_candidate_regions(
        [marker],
        genome,
        proteome,
        tmp_path / "a2",
        min_markers_initial=1,
        extension_kb=0,
        merge_distance=0,
        ablation_id=AblationID.A2,
        single_marker_min_score=100.0,
    )

    assert normal == []
    assert len(bypassed) == 1
    assert (bypassed[0].scaffold, bypassed[0].start, bypassed[0].end) == (
        "ctg",
        99,
        300,
    )
    assert len(bypassed[0].markers) == 1
    retained = bypassed[0].markers[0]
    assert retained.validation_status == "supported"
    assert retained.tier1_bypassed is True
    assert retained.is_validated is True
    assert retained.is_valid_seed_marker is True
    assert _region_coordinate_surface(normal) == set()
    assert _region_coordinate_surface(bypassed) == {("ctg", 99, 300)}


def test_counterfactual_region_assembly_does_not_write_artifacts(
    tmp_path: Path,
) -> None:
    genome, proteome = _inputs(tmp_path)
    output_dir = tmp_path / "counterfactual"

    regions = assemble_candidate_regions(
        [_strong_supported_marker()],
        genome,
        proteome,
        output_dir,
        min_markers_initial=1,
        extension_kb=0,
        merge_distance=0,
        ablation_id=AblationID.A2,
        single_marker_min_score=100.0,
        write_outputs=False,
    )

    assert len(regions) == 1
    assert not output_dir.exists()


def test_a2_bypass_label_survives_hallmark_evidence_synthesis() -> None:
    synthesizer = EvidenceSynthesizer(
        config=EvidenceSynthesizerConfig(ablation_id=AblationID.A2)
    )
    result = VerificationResult(
        eve_id="EVE_ctg_99-300",
        scaffold="ctg",
        start=99,
        end=300,
        ablation_id=AblationID.A2,
    )

    synthesizer._process_hallmark_hits(
        result,
        [
            {
                "porf_id": "ctg_1",
                "hallmark_gene": "GVOGm0003",
                "hmm_score": 100.0,
                "tier1_bypassed": True,
                "validation_status": "supported",
            }
        ],
    )

    assert result.hallmark_count == 1
    assert result.tier1_bypassed_marker_ids == ["ctg_1"]
    assert result.tier1_bypassed_marker_models == ["GVOGm0003"]
    document = result.to_dict()
    assert document["tier1_bypassed_marker_ids"] == ["ctg_1"]
    assert document["tier1_bypassed_marker_models"] == ["GVOGm0003"]
