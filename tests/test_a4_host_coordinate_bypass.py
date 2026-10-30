"""The Phase 2 A4 intervention bypasses exterior stops, never confirmed cores."""

import logging
from pathlib import Path

import pytest

from virosync.ablation import AblationID, InterventionCounts
from virosync.config import PipelineConfig
from virosync.orchestration._flows.single_genome import phase2
from virosync.pipeline import host_signatures
from virosync.pipeline.host_signatures import HostSignatureModel, TaxonomyLabelLookup
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy, GenomeDiamondQuery, SeedGeneMapping, pORF


@pytest.mark.parametrize(
    "ablation_id,expected_end,expected_counts",
    [(AblationID.A0, 100, InterventionCounts()), (AblationID.A4, 500, InterventionCounts(1, 1, 1))],
)
def test_phase2_a4_bypasses_three_gene_stop_and_records_counterfactual(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    ablation_id: AblationID,
    expected_end: int,
    expected_counts: InterventionCounts,
) -> None:
    monkeypatch.setattr(host_signatures, "_TAXONOMY_LOOKUP", TaxonomyLabelLookup({"EUK__host": "EUK|host"}))
    genome = tmp_path / "genome.fna"
    genome.write_text(">ctg\n" + "ACGT" * 200 + "\n", encoding="utf-8")
    genes = [pORF(id=f"gene{i}", scaffold="ctg", start=100 * i, end=100 * i + 100) for i in range(1, 5)]
    taxonomy = {
        gene.id: GeneTaxonomy(
            porf_id=gene.id,
            scaffold="ctg",
            start=gene.start,
            end=gene.end,
            has_hit=True,
            top1_prefix="EUK__" if i < 3 else "NCLDV__",
            top1_target="EUK__host" if i < 3 else "NCLDV__virus",
            top10_prefixes=["EUK__"] if i < 3 else ["NCLDV__"],
            top10_pidents=[80.0],
            top10_targets=["EUK__host"] if i < 3 else ["NCLDV__virus"],
            top10_bits=[100.0],
        )
        for i, gene in enumerate(genes)
    }
    query = GenomeDiamondQuery(
        seed_gene_mappings={
            "seed": SeedGeneMapping(
                seed_id="seed", scaffold="ctg", seed_start=0, seed_end=100, flank_start_bp=0, flank_end_bp=600
            ),
        }
    )
    monkeypatch.setattr(phase2, "build_proteome_index", lambda _path: {"ctg": genes})
    monkeypatch.setattr(phase2, "_require_phase2b_gene_taxonomy_db", lambda *_args, **_kwargs: tmp_path / "db")
    monkeypatch.setattr(
        phase2,
        "_run_boundary_diamond",
        lambda **_kwargs: phase2._BoundaryEvidence(
            taxonomy_map=taxonomy,
            control_stats=None,
            query=query,
        ),
    )
    config = PipelineConfig().with_overrides(resume=False, threads=1, gene_taxonomy_threads=1)
    config.ablation.id = ablation_id
    result = phase2._run_phase2_subflow(
        masked_path=genome,
        proteome_path=tmp_path / "proteome.faa",
        merged_seeds=[
            MergedSeed(
                scaffold="ctg",
                start=0,
                end=100,
                seed_id="seed",
                confirmed_seeds=(ConfirmedSeed("seed", "ctg", 0, 100),),
            )
        ],
        validated_markers=[],
        host_signature_model=HostSignatureModel(token_weights={"host": 1.0}),
        output_dir=tmp_path,
        genome_id="genome",
        config=config,
        genome_start_time=0,
        logger=logging.getLogger(__name__),
    )

    assert isinstance(result, phase2.Phase2Result)
    assert (result.refined_boundaries[0].start, result.refined_boundaries[0].end) == (0, expected_end)
    assert result.refined_boundaries[0].confirmed_seeds == (ConfirmedSeed("seed", "ctg", 0, 100),)
    assert result.ablation_counts == expected_counts
