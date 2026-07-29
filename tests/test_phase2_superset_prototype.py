from __future__ import annotations

import logging
import time
from pathlib import Path

import pytest

from virosync.config import ApplicationConfig, ConfigError, PipelineConfig
from virosync.orchestration._flows.single_genome import orchestrator, phase2
from virosync.orchestration._flows.single_genome.manifest import (
    _compute_config_fingerprint,
)
from virosync.orchestration._flows.single_genome.phase_state import (
    phase2_state_to_document,
)
from virosync.pipeline.host_signatures import HostSignatureModel
from virosync.pipeline.phase1.seed_merger import MergedSeed
from virosync.pipeline.phase2 import boundary_diamond
from virosync.pipeline.phase2.boundary_diamond import (
    BoundaryDiamondConfig,
    DiamondHit,
    GenomeDiamondQuery,
    SeedGeneMapping,
    classify_cached_diamond_query,
    pORF,
)
from virosync.pipeline.phase3 import gene_taxonomy


def _hit(query: str, target: str, bits: float) -> DiamondHit:
    return DiamondHit(
        query=query,
        target=target,
        evalue=1e-20,
        bits=bits,
        pident=80.0,
        qcov=90.0,
    )


def _write_proteome(path: Path) -> None:
    records = []
    for index, start in enumerate(range(1, 602, 100), start=1):
        end = start + 89
        records.append(
            f">contig_{index} # {start} # {end} # + # "
            f"ID=contig_{index};\nMPEPTIDE\n"
        )
    path.write_text("".join(records))


def test_full_proteome_search_is_one_unchunked_raw_id_pass(
    tmp_path: Path,
    monkeypatch,
) -> None:
    proteome = tmp_path / "proteome.faa"
    _write_proteome(proteome)
    calls: list[dict] = []

    def fake_search(**kwargs) -> None:
        calls.append(kwargs)
        Path(kwargs["output"]).write_text(
            "contig_1\tNCLDV__virus\t1e-20\t100\t80\t90\n"
        )

    monkeypatch.setattr(boundary_diamond, "run_diamond_blastp", fake_search)

    hits = boundary_diamond.run_full_proteome_diamond(
        proteome_fasta=proteome,
        diamond_db=tmp_path / "combined.dmnd",
        output_dir=tmp_path / "superset",
        max_target_seqs=17,
        threads=3,
    )

    assert len(calls) == 1
    assert calls[0]["query"] == proteome
    assert calls[0]["max_target_seqs"] == 17
    assert list(hits) == ["contig_1"]


def test_cached_phase2a_reconstructs_overlaps_no_hits_and_top10(
    tmp_path: Path,
    monkeypatch,
) -> None:
    proteome = tmp_path / "proteome.faa"
    _write_proteome(proteome)
    hits = {
        "contig_1": [
            _hit("contig_1", f"NCLDV__virus_{index}", 100.0 - index)
            for index in range(12)
        ],
        "contig_3": [_hit("contig_3", "EUK__host", 80.0)],
    }
    regions = [
        {"eve_id": "left", "scaffold": "contig", "start": 0, "end": 350},
        {"eve_id": "overlap", "scaffold": "contig", "start": 50, "end": 350},
    ]

    def fake_region_search(**kwargs) -> None:
        query_ids = [
            line[1:].strip()
            for line in Path(kwargs["query_fasta"]).read_text().splitlines()
            if line.startswith(">")
        ]
        lines = []
        for query_id in query_ids:
            raw_porf_id = query_id.split("|", 1)[1]
            for hit in hits.get(raw_porf_id, []):
                lines.append(
                    f"{query_id}\t{hit.target}\t{hit.evalue}\t{hit.bits}\t"
                    f"{hit.pident}\t{hit.qcov}\n"
                )
        Path(kwargs["output_file"]).write_text("".join(lines))

    monkeypatch.setattr(
        gene_taxonomy,
        "run_diamond_blastp",
        fake_region_search,
    )
    legacy = gene_taxonomy.run_gene_taxonomy_diamond_batch(
        regions=regions,
        proteome_fasta=proteome,
        combined_faa_db=tmp_path / "combined.dmnd",
        output_dir=tmp_path / "legacy",
    )
    result = gene_taxonomy.materialize_gene_taxonomy_batch_from_cached_hits(
        regions=regions,
        proteome_fasta=proteome,
        diamond_hits=hits,
        output_dir=tmp_path / "cached",
    )

    assert result == legacy
    left_records, left_summary = result["left"]
    overlap_records, overlap_summary = result["overlap"]
    assert [record.porf_id for record in left_records] == [
        "contig_1",
        "contig_3",
        "contig_2",
        "contig_4",
    ]
    assert [record.porf_id for record in overlap_records] == [
        "contig_1",
        "contig_3",
        "contig_2",
        "contig_4",
    ]
    assert len(left_records[0].top10_targets) == 10
    assert left_records[2].top1_prefix == "UNKNOWN"
    assert left_records[2].porf_start == 100
    assert left_summary["total"] == overlap_summary["total"] == 4
    assert (tmp_path / "cached" / "left.tsv").read_text() == (
        tmp_path / "legacy" / "left.tsv"
    ).read_text()
    assert (tmp_path / "cached" / "overlap.tsv").read_text() == (
        tmp_path / "legacy" / "overlap.tsv"
    ).read_text()


def test_cached_phase2b_slices_exact_query_and_preserves_controls(
    tmp_path: Path,
    monkeypatch,
) -> None:
    query = GenomeDiamondQuery(
        eve_porf_ids={"seed": ["p1"]},
        boundary_porf_ids={"seed": ["p2"]},
        control_porf_ids=["p3"],
        all_porf_ids=["p1", "p2", "p3"],
        seed_gene_mappings={
            "seed": SeedGeneMapping(
                seed_id="seed",
                scaffold="contig",
                seed_start=0,
                seed_end=100,
                eve_porf_ids=["p1"],
                downstream_porf_ids=["p2"],
            )
        },
    )
    proteome_index = {
        "contig": [
            pORF("p1", "contig", 0, 90),
            pORF("p2", "contig", 100, 190),
            pORF("p3", "contig", 200, 290),
        ]
    }
    hits = {
        "p1": [
            _hit("p1", f"NCLDV__virus_{index}", 100.0 - index)
            for index in range(12)
        ],
        "p2": [_hit("p2", "EUK__host", 90.0)],
    }

    config = BoundaryDiamondConfig(top_k=3, chunk_size=100)
    taxonomy = classify_cached_diamond_query(
        query=query,
        diamond_hits=hits,
        proteome_index=proteome_index,
        config=config,
    )

    assert list(taxonomy) == query.all_porf_ids
    assert len(taxonomy["p1"].top10_targets) == 3
    assert taxonomy["p3"].has_hit is False
    assert query.control_porf_ids == ["p3"]
    assert query.seed_gene_mappings["seed"].downstream_porf_ids == ["p2"]

    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">p1 # 1 # 90 # + # ID=contig_1;\nMPEPTIDE\n"
        ">p2 # 101 # 190 # + # ID=contig_2;\nMPEPTIDE\n"
        ">p3 # 201 # 290 # + # ID=contig_3;\nMPEPTIDE\n"
    )

    def fake_boundary_search(**kwargs) -> None:
        lines = []
        for porf_id in query.all_porf_ids:
            for hit in hits.get(porf_id, [])[: kwargs["max_target_seqs"]]:
                lines.append(
                    f"{porf_id}\t{hit.target}\t{hit.evalue}\t{hit.bits}\t"
                    f"{hit.pident}\t{hit.qcov}\n"
                )
        Path(kwargs["output"]).write_text("".join(lines))

    monkeypatch.setattr(
        boundary_diamond,
        "run_diamond_blastp",
        fake_boundary_search,
    )
    legacy = boundary_diamond.run_batched_diamond(
        query=query,
        proteome_fasta=proteome,
        diamond_db=tmp_path / "combined.dmnd",
        output_dir=tmp_path / "legacy_boundary",
        proteome_index=proteome_index,
        config=config,
    )
    assert taxonomy == legacy


def _run_phase2(
    *,
    output_dir: Path,
    masked: Path,
    proteome: Path,
    database: Path,
    prototype_enabled: bool,
) -> dict:
    return phase2._run_phase2_subflow(
        masked_path=masked,
        proteome_path=proteome,
        merged_seeds=[
            MergedSeed(
                scaffold="contig",
                start=100,
                end=350,
                seed_id="seed-1",
                sources=["hhg", "marker_validation"],
                confidence="high",
            )
        ],
        validated_markers=[],
        host_signature_model=HostSignatureModel(token_weights={"host": 1.0}),
        output_dir=output_dir,
        genome_id="genome",
        resume=False,
        refined_bed=output_dir / "phase2" / "refined_boundaries.bed",
        gene_taxonomy_faa_db=database,
        marker_db=None,
        taxonomy_labels_file=None,
        host_prefixes=["EUK__"],
        host_label="EUK",
        high_pident_host_threshold=70.0,
        boundary_host_trim_enabled=True,
        boundary_host_trim_window_bp=100,
        boundary_host_trim_step_bp=50,
        boundary_host_trim_max_host_fraction=0.5,
        boundary_host_trim_min_viral_fraction=0.0,
        boundary_host_trim_score_threshold=0.5,
        boundary_host_trim_buffer_kb=0,
        boundary_host_trim_min_overlap_score=0.2,
        boundary_host_signature_min_token_len=3,
        taxonomy_weight_mode="rank",
        boundary_taxonomy_ml_enabled=False,
        boundary_taxonomy_ml_model="logreg",
        boundary_taxonomy_ml_threshold=0.5,
        boundary_taxonomy_ml_neighbor_window=1,
        boundary_diamond_flank_genes=1,
        boundary_diamond_control_sample_size=1,
        boundary_diamond_control_min_distance=1,
        boundary_diamond_top_k=10,
        boundary_diamond_chunk_size=2,
        boundary_diamond_random_seed=42,
        threads=1,
        gene_taxonomy_threads=1,
        extended_output=False,
        search_backend="diamond",
        genome_start_time=time.time(),
        logger=logging.getLogger(__name__),
        boundary_diamond_superset_prototype_enabled=prototype_enabled,
    )


def test_superset_flow_uses_one_search_and_matches_legacy_phase2_surfaces(
    tmp_path: Path,
    monkeypatch,
) -> None:
    masked = tmp_path / "masked.fna"
    masked.write_text(">contig\n" + "ACGT" * 200 + "\n")
    proteome = tmp_path / "proteome.faa"
    _write_proteome(proteome)
    database = tmp_path / "combined.dmnd"
    database.write_bytes(b"database")
    raw_hits = {
        "contig_2": [_hit("contig_2", "EUK__host", 100.0)],
        "contig_3": [
            _hit("contig_3", f"NCLDV__virus_{index}", 90.0 - index)
            for index in range(12)
        ],
    }
    search_counts = {"legacy": 0, "prototype": 0}

    def fake_call_task(_task, **kwargs):
        search_counts["legacy"] += 1
        return gene_taxonomy.materialize_gene_taxonomy_batch_from_cached_hits(
            regions=kwargs["regions"],
            proteome_fasta=kwargs["proteome_path"],
            diamond_hits=raw_hits,
            output_dir=kwargs["output_dir"],
            high_pident_euk_threshold=kwargs["high_pident_host_threshold"],
        )

    def fake_legacy_boundary(**kwargs):
        search_counts["legacy"] += 1
        return classify_cached_diamond_query(
            query=kwargs["query"],
            diamond_hits=raw_hits,
            proteome_index=kwargs["proteome_index"],
            config=kwargs["config"],
            taxonomy_lookup=kwargs["taxonomy_lookup"],
        )

    def fake_superset(**kwargs):
        search_counts["prototype"] += 1
        output_dir = Path(kwargs["output_dir"])
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "full_proteome.tsv").write_text(
            "contig_2\tEUK__host\t1e-20\t100\t80\t90\n"
        )
        return raw_hits

    monkeypatch.setattr(phase2, "call_task", fake_call_task)
    monkeypatch.setattr(phase2, "run_batched_diamond", fake_legacy_boundary)
    monkeypatch.setattr(phase2, "run_full_proteome_diamond", fake_superset)

    legacy = _run_phase2(
        output_dir=tmp_path / "legacy",
        masked=masked,
        proteome=proteome,
        database=database,
        prototype_enabled=False,
    )
    prototype = _run_phase2(
        output_dir=tmp_path / "prototype",
        masked=masked,
        proteome=proteome,
        database=database,
        prototype_enabled=True,
    )

    assert search_counts == {"legacy": 2, "prototype": 1}
    assert phase2_state_to_document(legacy["refined_boundaries"]) == (
        phase2_state_to_document(prototype["refined_boundaries"])
    )
    assert legacy["boundary_taxonomy_map"] == prototype["boundary_taxonomy_map"]
    assert legacy["boundary_control_stats"] == prototype["boundary_control_stats"]
    assert legacy["boundary_diamond_query"] == prototype["boundary_diamond_query"]
    prototype_artifacts = orchestrator._phase_artifacts(
        tmp_path / "prototype",
        2,
    )
    assert "phase2/superset_diamond/full_proteome.tsv" in {
        artifact.relative_path for artifact in prototype_artifacts
    }


def test_superset_opt_in_round_trips_and_changes_provenance_fingerprint() -> None:
    default = PipelineConfig.from_dict({})
    enabled = ApplicationConfig.from_dict(
        {
            "phase2": {
                "boundary_diamond": {
                    "superset_prototype_enabled": True,
                }
            }
        }
    ).pipeline

    assert default.phase2.diamond_superset_prototype_enabled is False
    assert enabled.phase2.diamond_superset_prototype_enabled is True
    assert (
        enabled.to_flow_kwargs()[
            "boundary_diamond_superset_prototype_enabled"
        ]
        is True
    )
    assert _compute_config_fingerprint(default.to_flow_kwargs()) != (
        _compute_config_fingerprint(enabled.to_flow_kwargs())
    )


def test_superset_prototype_rejects_a_different_boundary_hit_limit() -> None:
    with pytest.raises(
        ConfigError,
        match="diamond_superset_prototype_enabled requires phase2.diamond_top_k=10",
    ):
        PipelineConfig.from_dict(
            {
                "phase2": {
                    "diamond_superset_prototype_enabled": True,
                    "diamond_top_k": 3,
                }
            }
        )
