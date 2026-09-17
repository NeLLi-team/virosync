"""Phase 3 subflow: evidence synthesis, verification, tiering."""

import json
from bisect import bisect_left
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

from virosync.ablation import AblationID, InterventionCounts
from virosync.config.pipeline_config import PipelineConfig
from virosync.orchestration.runtime import call_task
from virosync.orchestration.utils import get_genes_for_boundary
from virosync.output_contract import (
    EFFECTIVE_EVE_CLASSES,
    canonical_family,
    normalize_effective_eve_class,
)
from virosync.pipeline.phase1.marker_roles import decide_marker_hit_role
from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit
from virosync.pipeline.phase1.viral_markers import get_assembly_mode
from virosync.pipeline.phase2.boundary_diamond import (
    MIN_VIRAL_HIT_PIDENT,
    build_gene_taxonomy_record,
    get_flanking_taxonomy,
    missing_boundary_taxonomy_ids,
)
from virosync.pipeline.phase3.acceptance_selection import (
    select_phase3_acceptance,
)
from virosync.pipeline.phase3.eve_ani_clustering import (
    cluster_accepted_eves,
    recluster_survivors,
    unsupported_eve_ids,
)
from virosync.pipeline.phase3.evidence_synthesizer import VerificationResult
from virosync.pipeline.phase3.output_generator import (
    _is_atpase_marker,
    evaluate_v2_quality_gate,
)

from .loaders import (
    _load_interproscan_summary,
    _load_tmvec_cache,
    _serialize_tmvec_cache,
)

_ScaffoldStartIndex = dict[
    str,
    tuple[
        tuple[int, ...],
        tuple[tuple[int, int, int, object], ...],
        int,
    ],
]


@dataclass(frozen=True, slots=True)
class Phase3Result:
    """Outputs needed to publish a completed Phase 3 run."""

    verification_results: list
    accepted_results: list
    promoted_low_results: list
    classification_stats: dict[str, int]
    accepted: int
    elapsed: float
    ablation_counts: InterventionCounts = field(default_factory=InterventionCounts)


def _build_scaffold_start_index(records: Iterable[object]) -> _ScaffoldStartIndex:
    """Index coordinate records by scaffold and start while retaining input order."""
    grouped: dict[str, list[tuple[int, int, int, object]]] = defaultdict(list)
    for ordinal, record in enumerate(records):
        grouped[record.scaffold].append((record.start, record.end, ordinal, record))

    index: _ScaffoldStartIndex = {}
    for scaffold, entries in grouped.items():
        entries.sort(key=lambda item: item[0])
        index[scaffold] = (
            tuple(item[0] for item in entries),
            tuple(entries),
            max(max(0, item[1] - item[0]) for item in entries),
        )
    return index


def _query_scaffold_index(
    index: _ScaffoldStartIndex,
    *,
    scaffold: str,
    start: int,
    end: int,
) -> list[object]:
    """Return half-open overlaps in the records' original insertion order."""
    scaffold_index = index.get(scaffold)
    if scaffold_index is None:
        return []
    starts, entries, max_length = scaffold_index
    lower = bisect_left(starts, start - max_length)
    upper = bisect_left(starts, end)
    matches = [item for item in entries[lower:upper] if item[0] < end and item[1] > start]
    matches.sort(key=lambda item: item[2])
    return [item[3] for item in matches]


def _query_boundary_coordinate_records(
    *,
    boundary,
    taxonomy_index: _ScaffoldStartIndex,
    marker_index: _ScaffoldStartIndex,
    taxonomy_map: dict | None = None,
    proteome_index: dict | None = None,
) -> tuple[list[object], list[object]]:
    """Return ordered taxonomy and marker overlaps for one boundary."""
    query = {
        "scaffold": boundary.scaffold,
        "start": boundary.start,
        "end": boundary.end,
    }
    if taxonomy_map is not None and proteome_index is not None:
        missing_ids = missing_boundary_taxonomy_ids(
            **query,
            taxonomy_map=taxonomy_map,
            proteome_index=proteome_index,
        )
        if missing_ids:
            raise RuntimeError(
                "Phase 3 taxonomy coverage is incomplete for "
                f"{boundary.scaffold}:{boundary.start}-{boundary.end}: "
                f"{len(missing_ids)} overlapping genes were not searched"
            )

    return (
        _query_scaffold_index(taxonomy_index, **query),
        _query_scaffold_index(marker_index, **query),
    )


def _is_marker_floor_recovery_candidate(result) -> bool:
    """Return whether selection excluded a candidate for a recoverable reason."""
    return result.canonical_selection_outcome in {
        "normal_gate_rejected",
        "rescue_marker_excluded",
    }


def _resolve_tmvec_device(device: str) -> str:
    """Honor the configured TMVec device without implicit promotion/demotion."""
    if device == "cpu":
        return "cpu"
    if device != "cuda":
        raise RuntimeError(f"Unsupported TMVec device: {device}")
    try:
        import torch
    except Exception as exc:
        raise RuntimeError("TMVec CUDA validation requires PyTorch") from exc
    if not torch.cuda.is_available():
        raise RuntimeError("TMVec requested CUDA but CUDA is not available")
    return "cuda"


def _run_interproscan(
    *,
    regions: list[dict[str, object]],
    proteome_path: Path,
    output_dir: Path,
    config: PipelineConfig,
    resume_authorized: bool,
    logger: object,
) -> object | None:
    """Load or compute optional InterProScan evidence for all candidates."""
    phase3 = config.phase3
    if not phase3.interproscan_enabled:
        return None

    summary_path = output_dir / "phase3" / "interproscan" / "interproscan_summary.tsv"
    interproscan_map = None
    if config.execution.resume and resume_authorized and summary_path.exists():
        try:
            interproscan_map = _load_interproscan_summary(summary_path)
            logger.info(
                "Phase 3 resume: loaded InterProScan summary for %d candidates",
                len(interproscan_map),
            )
        except Exception as exc:
            logger.warning(
                "Phase 3 resume: failed loading InterProScan summary (%s); rerunning",
                exc,
            )

    if interproscan_map is not None:
        return interproscan_map
    if not phase3.interproscan_dir:
        logger.warning("Phase 3: InterProScan enabled but interproscan_dir not set; skipping")
        return None
    if not phase3.interproscan_dir.exists():
        logger.warning("Phase 3: InterProScan dir not found: %s; skipping", phase3.interproscan_dir)
        return None

    from virosync.orchestration.tasks import interproscan_batch_task

    threads = config.compute.effective_threads()
    configured_threads = config.compute.interproscan_threads
    interpro_threads = configured_threads if configured_threads is not None else threads
    interpro_threads = max(1, min(interpro_threads, threads))
    logger.info("Phase 3: InterProScan threads=%s (batch)", interpro_threads)
    try:
        return call_task(
            interproscan_batch_task,
            regions=regions,
            proteome_path=proteome_path,
            interproscan_dir=phase3.interproscan_dir,
            output_dir=output_dir / "phase3" / "interproscan",
            threads=interpro_threads,
            keywords=phase3.interproscan_keywords,
            applications=phase3.interproscan_applications,
        )
    except Exception as exc:
        logger.error(
            "Phase 3: InterProScan failed after it was enabled; refusing silent fallback: %s",
            exc,
        )
        raise


def _precompute_tmvec(
    *,
    refined_boundaries: list,
    proteome_path: Path,
    output_dir: Path,
    config: PipelineConfig,
    resume_authorized: bool,
    logger: object,
) -> tuple[object | None, str]:
    """Load or compute the shared TMVec result set for Phase 3."""
    phase3 = config.phase3
    requested_device = config.compute.device.value
    tmvec_device = _resolve_tmvec_device(requested_device) if phase3.use_tmvec_database else requested_device
    if not phase3.use_tmvec_database:
        return None, tmvec_device

    cache_path = output_dir / "phase3" / "tmvec" / "precomputed_tmvec.json"
    if config.execution.resume and resume_authorized and cache_path.exists():
        try:
            cached = _load_tmvec_cache(cache_path)
            logger.info(
                "Phase 3 resume: loaded cached TMVec hits for %d proteins",
                len(cached),
            )
            return cached, tmvec_device
        except Exception as exc:
            logger.warning(
                "Phase 3 resume: failed loading TMVec cache (%s); recomputing",
                exc,
            )

    if not refined_boundaries:
        return None, tmvec_device

    raw_protein_count = 0
    conflicting_ids = 0
    protein_by_id: dict[str, str] = {}
    for boundary in refined_boundaries:
        porf_sequences = get_genes_for_boundary(
            proteome_path=proteome_path,
            scaffold=boundary.scaffold,
            start=boundary.start,
            end=boundary.end,
            max_porfs=10000,
        )
        if porf_sequences:
            raw_protein_count += len(porf_sequences)
            for porf_id, sequence in porf_sequences:
                prior = protein_by_id.get(porf_id)
                if prior is None:
                    protein_by_id[porf_id] = sequence
                elif prior != sequence:
                    conflicting_ids += 1

    all_proteins = list(protein_by_id.items())
    if not all_proteins:
        return None, tmvec_device

    logger.info(
        "Phase 3: TMVec batch - collecting %d proteins from %d EVEs (unique pORFs=%d)",
        raw_protein_count,
        len(refined_boundaries),
        len(all_proteins),
    )
    if conflicting_ids:
        logger.warning(
            "Phase 3: TMVec batch - %d duplicate pORF IDs had conflicting sequences; using first occurrence",
            conflicting_ids,
        )
    try:
        from virosync.pipeline.phase3.tmvec_database import TMVecDatabaseSearch

        searcher = TMVecDatabaseSearch(
            device=tmvec_device,
            databases=phase3.tmvec_databases or ["bfvd"],
            min_tm=0.0,
            database_root=phase3.tmvec_database_dir,
            require_gpu=phase3.tmvec_require_gpu,
            fail_on_unavailable=True,
        )
        precomputed_tmvec = searcher.search_batch(all_proteins)
        logger.info(
            "Phase 3: TMVec batch - completed %d protein searches",
            len(precomputed_tmvec) if precomputed_tmvec else 0,
        )
        predictor = getattr(searcher, "_predictor", None)
        if predictor is not None:
            oom = getattr(predictor, "_batch_oom_fallbacks", 0)
            per_seq_failures = getattr(predictor, "_per_seq_fallback_failures", 0)
            if oom or per_seq_failures:
                logger.warning(
                    "Phase 3: TMVec GPU degradation detected — batch_oom_fallbacks=%d per_seq_fallback_failures=%d",
                    oom,
                    per_seq_failures,
                )
        if precomputed_tmvec:
            _serialize_tmvec_cache(precomputed_tmvec, cache_path)
            logger.info("Phase 3: TMVec cache written to %s", cache_path)
        return precomputed_tmvec, tmvec_device
    except Exception as exc:
        logger.error(
            "Phase 3: TMVec failed after preflight enabled it; refusing silent fallback: %s",
            exc,
        )
        raise


def _load_or_classify_jelly_roll(
    *,
    validated_hits_tsv: Path,
    output_dir: Path,
    config: PipelineConfig,
    resume_authorized: bool,
    logger: object,
) -> dict[str, list[dict]]:
    """Load or compute MCP DJR/SJR classifications used during scoring."""
    from virosync.pipeline.phase3.evidence_synthesizer import load_jelly_roll_data

    jelly_roll_output = output_dir / "phase3_synthesis" / "virosync_jelly_roll_proteins.tsv"
    marker_sequences_path = output_dir / "phase1" / "marker_validation" / "hmm_hit_porfs.faa"
    interproscan_batch_path = output_dir / "phase3" / "interproscan" / "interproscan_batch.tsv"
    foldseek_results_path = output_dir / "structural_analysis" / "foldseek_pdb_results.tsv"

    jelly_roll_map: dict[str, list[dict]] = {}
    if config.execution.resume and resume_authorized and jelly_roll_output.exists():
        jelly_roll_map = load_jelly_roll_data(jelly_roll_output)
        if jelly_roll_map:
            logger.info(
                "Phase 3 resume: loaded jelly-roll classifications for %d MCP proteins",
                len(jelly_roll_map),
            )

    if jelly_roll_map:
        return jelly_roll_map
    if not validated_hits_tsv.exists() or not marker_sequences_path.exists():
        logger.info("Phase 3: skipping jelly-roll classification before scoring (missing marker validation files)")
        return jelly_roll_map

    from virosync.orchestration.tasks import classify_jelly_roll_task

    logger.info("Phase 3: classifying MCP proteins as DJR/SJR before confidence scoring")
    call_task(
        classify_jelly_roll_task,
        marker_hits_path=validated_hits_tsv,
        sequences_path=marker_sequences_path,
        output_path=jelly_roll_output,
        interproscan_path=interproscan_batch_path if interproscan_batch_path.exists() else None,
        tmvec_results_path=None,
        foldseek_results_path=foldseek_results_path if foldseek_results_path.exists() else None,
    )
    jelly_roll_map = load_jelly_roll_data(jelly_roll_output)
    if jelly_roll_map:
        logger.info(
            "Phase 3: jelly-roll classifications ready for %d MCP proteins",
            len(jelly_roll_map),
        )
    return jelly_roll_map


@dataclass(frozen=True, slots=True)
class _BoundaryEvidence:
    """Evidence maps produced for one candidate boundary."""

    boundary_id: str
    hallmarks: list[dict[str, object]]
    gene_taxonomy: object | None
    interproscan: object | None


def _summarize_boundary_taxonomy(
    *,
    boundary: object,
    filtered_taxonomy: list[object],
    boundary_taxonomy_map: dict,
    proteome_index: dict,
    boundary_diamond_query: object | None,
    config: PipelineConfig,
    logger: object,
) -> tuple[list[dict], dict[str, object]] | None:
    """Build the interior and flank taxonomy evidence for one boundary."""
    if not filtered_taxonomy:
        return None

    host_prefix = f"{config.host.label}__"
    gene_tax_records = []
    n_viral_interior = 0
    n_ncldv_mirus_interior = 0
    n_vp_plv_interior = 0
    n_host = 0
    for taxonomy in filtered_taxonomy:
        gene_tax_records.append(build_gene_taxonomy_record(taxonomy, is_flanking=False))
        n_viral_interior += int(taxonomy.has_viral)
        n_ncldv_mirus_interior += int(taxonomy.has_ncldv_mirus)
        n_vp_plv_interior += int(taxonomy.has_vp_plv)
        n_host += int(taxonomy.top1_prefix == host_prefix)

    seed_mapping = None
    if boundary_diamond_query and boundary.seed_id:
        seed_mapping = boundary_diamond_query.seed_gene_mappings.get(boundary.seed_id)
    upstream_taxonomy, downstream_taxonomy = get_flanking_taxonomy(
        taxonomy_map=boundary_taxonomy_map,
        proteome_index=proteome_index,
        refined_boundary=boundary,
        flank_genes=config.phase2.diamond_flank_genes,
        seed_mapping=seed_mapping,
    )
    n_flanking = 0
    n_viral_flanking = 0
    n_ncldv_mirus_flanking = 0
    n_vp_plv_flanking = 0
    for flank_position, records in (
        ("upstream", upstream_taxonomy),
        ("downstream", downstream_taxonomy),
    ):
        for taxonomy in records:
            gene_tax_records.append(
                build_gene_taxonomy_record(
                    taxonomy,
                    is_flanking=True,
                    flank_position=flank_position,
                )
            )
            n_flanking += 1
            n_viral_flanking += int(taxonomy.has_viral)
            n_ncldv_mirus_flanking += int(taxonomy.has_ncldv_mirus)
            n_vp_plv_flanking += int(taxonomy.has_vp_plv)

    family_counts = {
        family: sum(
            1
            for taxonomy in filtered_taxonomy
            if family
            in {
                canonical_family(prefix.rstrip("_"))
                for prefix, pident in zip(
                    taxonomy.top10_prefixes or [],
                    taxonomy.top10_pidents or [],
                )
                if pident >= MIN_VIRAL_HIT_PIDENT
            }
        )
        for family in ("NCLDV", "MIRUS", "PPV", "CRESS")
    }
    dominant_family = "UNKNOWN"
    dominant_fraction = 0.0
    if filtered_taxonomy:
        dominant_family = max(family_counts, key=family_counts.get)
        max_count = family_counts[dominant_family]
        if max_count > 0:
            dominant_fraction = max_count / len(filtered_taxonomy)
        else:
            dominant_family = "UNKNOWN"

    n_viral_total = n_viral_interior + n_viral_flanking
    n_ncldv_mirus_total = n_ncldv_mirus_interior + n_ncldv_mirus_flanking
    summary = {
        "total": len(filtered_taxonomy),
        "total_with_flanking": len(gene_tax_records),
        "flanking_genes": n_flanking,
        "ncldv_mirus": n_ncldv_mirus_interior,
        "vp_plv": n_vp_plv_interior,
        "viral_top10": n_viral_interior,
        "high_pident_euk": n_host,
        "has_ncldv_mirus": n_ncldv_mirus_interior > 0,
        "has_vp_plv": n_vp_plv_interior > 0,
        "dominant_family": dominant_family,
        "dominant_fraction": dominant_fraction,
        "viral_interior": n_viral_interior,
        "viral_flanking": n_viral_flanking,
        "ncldv_mirus_interior": n_ncldv_mirus_interior,
        "ncldv_mirus_flanking": n_ncldv_mirus_flanking,
        "vp_plv_interior": n_vp_plv_interior,
        "vp_plv_flanking": n_vp_plv_flanking,
    }
    boundary_id = f"{boundary.scaffold}_{boundary.start}_{boundary.end}"
    logger.info(
        "%s: Reused Phase 2b taxonomy - %d interior + %d flanking genes, "
        "%d viral (%.1f%%), %d NCLDV/MIRUS (interior: %d viral, flanking: %d viral)",
        boundary_id,
        len(filtered_taxonomy),
        n_flanking,
        n_viral_total,
        100.0 * n_viral_total / max(1, len(gene_tax_records)),
        n_ncldv_mirus_total,
        n_viral_interior,
        n_viral_flanking,
    )
    return gene_tax_records, summary


def _build_boundary_hallmarks(
    boundary_markers: list[object],
    config: PipelineConfig,
) -> list[dict[str, object]]:
    """Build retained marker evidence for one boundary."""
    hallmarks = []
    single_marker_min_score = get_assembly_mode(config.phase1.assembly_mode.value).single_marker_min_score
    for marker in boundary_markers:
        marker_role = decide_marker_hit_role(
            marker,
            ablation_id=config.ablation.id,
            single_marker_min_score=single_marker_min_score,
        )
        if not marker_role.is_retained_evidence:
            continue
        porf_id = getattr(marker, "porf_id", None) or getattr(marker, "query_porf", None)
        hallmarks.append(
            {
                "start": marker.start,
                "end": marker.end,
                "porf_id": porf_id,
                "hallmark_gene": marker.hmm_target,
                "hmm_score": marker.hmm_score,
                "score": marker.hmm_score,
                "top10_prefixes": getattr(marker, "top10_prefixes", ""),
                "top10_targets": getattr(marker, "top10_targets", ""),
                "top10_pidents": getattr(marker, "top10_pidents", ""),
                "best_hit_target": getattr(marker, "best_hit_target", ""),
                "validation_status": marker_role.original_validation_status,
                "tier1_bypassed": marker_role.is_tier1_bypassed,
            }
        )
    return hallmarks


def _build_boundary_evidence(
    *,
    boundary: object,
    taxonomy_index: _ScaffoldStartIndex,
    marker_index: _ScaffoldStartIndex,
    boundary_taxonomy_map: dict,
    proteome_index: dict,
    boundary_diamond_query: object | None,
    interproscan_map: object | None,
    use_precomputed_taxonomy: bool,
    proteome_path: Path,
    config: PipelineConfig,
    logger: object,
) -> _BoundaryEvidence:
    """Build the complete verification evidence for one candidate boundary."""
    filtered_taxonomy, boundary_markers = _query_boundary_coordinate_records(
        boundary=boundary,
        taxonomy_index=taxonomy_index,
        marker_index=marker_index,
        taxonomy_map=boundary_taxonomy_map if use_precomputed_taxonomy else None,
        proteome_index=proteome_index if use_precomputed_taxonomy else None,
    )
    gene_taxonomy = None
    if use_precomputed_taxonomy:
        gene_taxonomy = _summarize_boundary_taxonomy(
            boundary=boundary,
            filtered_taxonomy=filtered_taxonomy,
            boundary_taxonomy_map=boundary_taxonomy_map,
            proteome_index=proteome_index,
            boundary_diamond_query=boundary_diamond_query,
            config=config,
            logger=logger,
        )
    if gene_taxonomy is None and config.ablation.id is AblationID.A3:
        porf_sequences = get_genes_for_boundary(
            proteome_path=proteome_path,
            scaffold=boundary.scaffold,
            start=boundary.start,
            end=boundary.end,
            max_porfs=10000,
        )
        gene_count = len(porf_sequences)
        gene_taxonomy = (
            [],
            {
                "total": gene_count,
                "total_with_flanking": gene_count,
                "flanking_genes": 0,
            },
        )

    eve_id = f"EVE_{boundary.scaffold}_{boundary.start}-{boundary.end}"
    interproscan = interproscan_map.get(eve_id) if interproscan_map else None
    return _BoundaryEvidence(
        boundary_id=f"{boundary.scaffold}_{boundary.start}_{boundary.end}",
        hallmarks=_build_boundary_hallmarks(boundary_markers, config),
        gene_taxonomy=gene_taxonomy,
        interproscan=interproscan,
    )


def _run_verification(
    *,
    boundaries: list,
    hallmark_hits_map: dict,
    gene_taxonomy_map: dict,
    interproscan_map: dict,
    jelly_roll_map: dict,
    masked_path: Path,
    proteome_path: Path,
    output_dir: Path,
    host_signatures: set,
    host_signature_model_payload: dict | None,
    precomputed_tmvec: object | None,
    tmvec_device: str,
    config: PipelineConfig,
) -> list:
    """Run the shared batched verification path for candidate boundaries."""
    from virosync.orchestration.tasks import verify_eve_candidates_batched_task

    databases = config.databases
    phase3 = config.phase3
    return call_task(
        verify_eve_candidates_batched_task,
        boundaries=boundaries,
        genome_path=masked_path,
        work_dir=output_dir / "phase3",
        proteome_path=proteome_path,
        hallmark_hits_map=hallmark_hits_map,
        novelty_scores={},
        gene_taxonomy_map=gene_taxonomy_map or None,
        interproscan_map=interproscan_map or None,
        jelly_roll_map=jelly_roll_map or None,
        euk_host_signatures=host_signatures,
        host_signature_model=host_signature_model_payload,
        host_signature_score_threshold=phase3.host_signature_evidence_threshold,
        host_prefixes=config.host.prefixes,
        host_label=config.host.label,
        high_tier_threshold=phase3.high_tier_threshold,
        low_tier_threshold=phase3.low_tier_threshold,
        use_crf_in_final_score=phase3.use_crf_in_final_score,
        priority_marker_list=phase3.priority_marker_list,
        marker_floor_priority_only=phase3.marker_floor_priority_only,
        marker_floor_priority_plus_family=phase3.marker_floor_priority_plus_family,
        marker_floor_priority_multi_family=phase3.marker_floor_priority_multi_family,
        marker_family_bonus_per_family=phase3.marker_family_bonus_per_family,
        marker_multi_family_bonus=phase3.marker_multi_family_bonus,
        skip_structural=phase3.skip_structural,
        use_boltz=phase3.use_boltz,
        boltz_mcp_only=phase3.boltz_mcp_only,
        boltz_use_msa_server=phase3.boltz_use_msa_server,
        boltz_min_seq_len=phase3.boltz_min_seq_len,
        boltz_max_seq_len=phase3.boltz_max_seq_len,
        boltz_no_kernels=phase3.boltz_no_kernels,
        use_tmvec_database=phase3.use_tmvec_database,
        tmvec_databases=phase3.tmvec_databases,
        tmvec_database_dir=phase3.tmvec_database_dir,
        tmvec_min_score=phase3.tmvec_min_score,
        tmvec_require_gpu=phase3.tmvec_require_gpu,
        device=tmvec_device,
        viral_structure_db=phase3.viral_structure_db,
        gvclass_db=databases.gvclass_db,
        diamond_db=databases.diamond_db,
        enable_phylogenetic=phase3.enable_phylogenetic,
        taxonomy_labels_file=databases.taxonomy_labels_file,
        hmm_database=databases.hmm_database,
        precomputed_tmvec=precomputed_tmvec,
        max_workers=config.compute.effective_threads(),
        ablation_id=config.ablation.id,
    )


def _run_phase3_subflow(
    masked_path: Path,
    proteome_path: Path,
    validated_markers: list,
    host_signatures: set,
    host_signature_model_payload: dict | None,
    refined_boundaries: list,
    boundary_taxonomy_map: dict,
    boundary_diamond_query,
    proteome_index: dict,
    merged_seeds: list,
    output_dir: Path,
    config: PipelineConfig,
    validated_hits_tsv: Path,
    logger: object,
    resume_authorized: bool = False,
) -> Phase3Result:
    """Phase 3: Evidence synthesis and verification.

    This phase verifies candidate regions through:
    1. Loading validated markers on resume
    2. InterProScan batch analysis (optional)
    3. TMVec structural similarity batch (optional)
    4. Building verification maps
    5. Batched verification of all candidates
    6. Computing confidence tiers and classification stats

    Args:
        masked_path: Path to masked genome FASTA
        proteome_path: Path to protein FASTA
        validated_markers: List of validated marker hits from Phase 1
        config: Resolved nested pipeline configuration.
        logger: Logger instance

    Returns:
        Explicit Phase 3 outputs needed by final publication.
    """
    import time

    resume = config.execution.resume
    ablation_id = config.ablation.id
    threads = config.compute.effective_threads()
    host_prefixes = config.host.prefixes

    phase3_start = time.time()
    logger.info("-" * 60)
    logger.info(f"Phase 3: Verifying {len(refined_boundaries)} candidates")
    logger.info("  Steps: gene taxonomy -> evidence synthesis -> structural (if enabled)")

    # Load validated markers on resume if needed
    if (not validated_markers) and resume and resume_authorized and validated_hits_tsv.exists():
        from virosync.pipeline.host_signatures import HostSignatureModel
        from virosync.pipeline.phase1.marker_validation import (
            collect_host_signatures,
            load_validated_marker_hits,
        )

        validated_markers = load_validated_marker_hits(validated_hits_tsv)
        host_signatures = collect_host_signatures(
            validated_markers,
            host_prefixes=set(host_prefixes),
        )
        model_path = output_dir / "phase1" / "marker_validation" / "host_signature_model.json"
        if model_path.exists():
            with model_path.open() as handle:
                host_signature_model_payload = json.load(handle)
                # Parse eagerly and discard: the payload is consumed downstream
                # as a raw dict, so this is the only place a resumed run rejects
                # a well-formed but wrongly-typed model file before Phase 3 work.
                HostSignatureModel.from_dict(host_signature_model_payload)
        logger.info("Resume: loaded %d validated marker hits", len(validated_markers))

    # NOTE: Phase 3 gene taxonomy is ALWAYS pre-computed in Phase 2b
    # No fallback Diamond batch in Phase 3 (Phase 2b is mandatory)

    # Prepare regions payload for gene taxonomy batch task
    regions_payload = [
        {
            "eve_id": f"EVE_{b.scaffold}_{b.start}-{b.end}",
            "scaffold": b.scaffold,
            "start": b.start,
            "end": b.end,
        }
        for b in refined_boundaries
    ]
    # Use pre-computed taxonomy from Phase 2b (MANDATORY - no Phase 3 Diamond fallback)
    # Phase 2b runs Diamond ONCE per genome on all seeds + +/-20 flanking genes
    # Boundary constraints ensure refined boundaries never exceed Phase 2b coverage
    use_precomputed_taxonomy = bool(boundary_taxonomy_map)
    if use_precomputed_taxonomy:
        logger.info(
            "Phase 3: Using pre-computed taxonomy from Phase 2b (%d pORFs)",
            len(boundary_taxonomy_map),
        )
    elif merged_seeds and ablation_id is AblationID.A3:
        logger.info("A3: gene taxonomy unavailable because Tier-2 search is bypassed")
    elif merged_seeds:
        # Outside A3, Phase 2b is mandatory when seeds exist.
        raise RuntimeError("Phase 3 requires complete Phase 2b gene taxonomy when seeds exist")

    interproscan_map = _run_interproscan(
        regions=regions_payload,
        proteome_path=proteome_path,
        output_dir=output_dir,
        config=config,
        resume_authorized=resume_authorized,
        logger=logger,
    )
    precomputed_tmvec, tmvec_device = _precompute_tmvec(
        refined_boundaries=refined_boundaries,
        proteome_path=proteome_path,
        output_dir=output_dir,
        config=config,
        resume_authorized=resume_authorized,
        logger=logger,
    )
    jelly_roll_map = _load_or_classify_jelly_roll(
        validated_hits_tsv=validated_hits_tsv,
        output_dir=output_dir,
        config=config,
        resume_authorized=resume_authorized,
        logger=logger,
    )

    # ==================================================
    # Build maps for batched verification
    # ==================================================
    hallmark_hits_map = {}
    gene_taxonomy_map_batched = {}
    interproscan_map_batched = {}
    boundary_taxonomy_index = _build_scaffold_start_index(boundary_taxonomy_map.values())
    validated_marker_index = _build_scaffold_start_index(validated_markers)

    logger.info(
        "Phase 3: Preparing data for batched verification of %d boundaries",
        len(refined_boundaries),
    )
    for boundary in refined_boundaries:
        evidence = _build_boundary_evidence(
            boundary=boundary,
            taxonomy_index=boundary_taxonomy_index,
            marker_index=validated_marker_index,
            boundary_taxonomy_map=boundary_taxonomy_map,
            proteome_index=proteome_index,
            boundary_diamond_query=boundary_diamond_query,
            interproscan_map=interproscan_map,
            use_precomputed_taxonomy=use_precomputed_taxonomy,
            proteome_path=proteome_path,
            config=config,
            logger=logger,
        )
        hallmark_hits_map[evidence.boundary_id] = evidence.hallmarks
        if evidence.gene_taxonomy:
            gene_taxonomy_map_batched[evidence.boundary_id] = evidence.gene_taxonomy
        if evidence.interproscan:
            interproscan_map_batched[evidence.boundary_id] = evidence.interproscan

    # ==================================================
    # Run batched verification
    # ==================================================
    logger.info("Phase 3: Running batched verification with %d threads", threads)
    verification_results = _run_verification(
        boundaries=refined_boundaries,
        hallmark_hits_map=hallmark_hits_map,
        gene_taxonomy_map=gene_taxonomy_map_batched,
        interproscan_map=interproscan_map_batched,
        jelly_roll_map=jelly_roll_map,
        masked_path=masked_path,
        proteome_path=proteome_path,
        output_dir=output_dir,
        host_signatures=host_signatures,
        host_signature_model_payload=host_signature_model_payload,
        precomputed_tmvec=precomputed_tmvec,
        tmvec_device=tmvec_device,
        config=config,
    )

    # First acceptance pass. It identifies which candidates the v2 gate
    # rejects; the marker-floor re-admit below can only add to that set.
    acceptance_selection = select_phase3_acceptance(
        verification_results,
        ablation_id,
    )

    # === Phase-3 marker-floor re-admit (recall recovery) ===
    # Phase-2 host trimming can collapse a marker-dense NCLDV/MIRUS seed below 5 kb,
    # stripping validated hallmark markers out of the boundary so the v2 gate
    # rejects a genuine EVE. Rather than MUTATE the boundary in Phase 2 (extending an
    # ALREADY-ACCEPTED region pulls host genes in, lowers confidence, drops the tier
    # MEDIUM->LOW, and the stricter LOW NCLDV gate then rejects it -> NCLDV genes lost
    # on rhizophagus/tstriata), recover ONLY here and ONLY additively: for each
    # REJECTED region whose own seed span carried >=2 validated markers
    # (annotate_boundaries_with_marker_floor stored marker_floor_start/end),
    # synthesize an ALTERNATIVE at the floored span via the same evidence path and
    # re-admit it iff it (a) independently passes the v2 gate, (b) carries a genuine
    # non-ATPase viral hallmark inside the floored span, (c) has TIER-2 interior
    # viral fraction >= READMIT_MIN_VIRAL_FRACTION (so host-with-a-marker regions are
    # rejected), and (d) does not overlap any accepted region. No boundary already
    # accepted is modified, and the alternatives only ever ADD candidates, so the
    # accepted set cannot shrink. The selector is re-run afterwards over the
    # extended list. Direct-overlap arbitration prevents a non-overlapping
    # alternative from changing an original winner through a transitive bridge.
    #
    READMIT_MIN_VIRAL_FRACTION = 0.10
    # Regions accepted so far, seeded from the first pass and extended as
    # alternatives qualify, so no two re-admits can overlap each other.
    readmit_accepted = list(acceptance_selection.canonical_results)
    boundary_by_region = {(b.scaffold, b.start, b.end): b for b in refined_boundaries}
    readmit_boundaries = []
    for r in verification_results:
        if not _is_marker_floor_recovery_candidate(r):
            continue
        boundary = boundary_by_region.get((r.scaffold, r.start, r.end))
        if boundary is None:
            continue
        if getattr(boundary, "tir_boundary_override", False):
            continue
        floor_start = getattr(boundary, "marker_floor_start", None)
        floor_end = getattr(boundary, "marker_floor_end", None)
        if floor_start is None or floor_end is None:
            continue
        # Only a strictly wider floored span gives the gate something new to see.
        if not (floor_start < boundary.start or floor_end > boundary.end):
            continue
        readmit_boundaries.append(
            replace(
                boundary,
                start=min(boundary.start, floor_start),
                end=max(boundary.end, floor_end),
            )
        )

    n_readmitted = 0
    if readmit_boundaries:
        from .phase2 import _recalculate_boundary_composition

        _recalculate_boundary_composition(
            readmit_boundaries,
            masked_path=masked_path,
        )
        # Evidence parity is partial by construction. Hallmarks and Phase-2b gene
        # taxonomy are recomputed for the floored span, but InterProScan and TMVec
        # were precomputed for the ORIGINAL boundaries only, so proteins that the
        # floor newly includes carry no InterPro or structural evidence here.
        # That can only lower an alternative's score, so it costs recoveries
        # rather than admitting anything extra. Re-running those two tools per
        # alternative is the fix if the recovery rate proves too low.
        alt_hh_map, alt_gt_map, alt_ip_map = {}, {}, {}
        for alt_boundary in readmit_boundaries:
            evidence = _build_boundary_evidence(
                boundary=alt_boundary,
                taxonomy_index=boundary_taxonomy_index,
                marker_index=validated_marker_index,
                boundary_taxonomy_map=boundary_taxonomy_map,
                proteome_index=proteome_index,
                boundary_diamond_query=boundary_diamond_query,
                interproscan_map=interproscan_map,
                use_precomputed_taxonomy=use_precomputed_taxonomy,
                proteome_path=proteome_path,
                config=config,
                logger=logger,
            )
            alt_hh_map[evidence.boundary_id] = evidence.hallmarks
            if evidence.gene_taxonomy:
                alt_gt_map[evidence.boundary_id] = evidence.gene_taxonomy
            if evidence.interproscan:
                alt_ip_map[evidence.boundary_id] = evidence.interproscan
        alt_results = _run_verification(
            boundaries=readmit_boundaries,
            hallmark_hits_map=alt_hh_map,
            gene_taxonomy_map=alt_gt_map,
            interproscan_map=alt_ip_map,
            jelly_roll_map=jelly_roll_map,
            masked_path=masked_path,
            proteome_path=proteome_path,
            output_dir=output_dir,
            host_signatures=host_signatures,
            host_signature_model_payload=host_signature_model_payload,
            precomputed_tmvec=precomputed_tmvec,
            tmvec_device=tmvec_device,
            config=config,
        )
        # Deterministic order so the highest-confidence alternative wins any overlap.
        alt_results.sort(key=lambda a: (-getattr(a, "final_confidence", 0.0), a.scaffold, a.start))
        for alt in alt_results:
            # (a) the floored alternative must independently pass the v2 gate.
            if not evaluate_v2_quality_gate(alt).kept:
                continue
            # (b) viral-content criterion: a genuine non-ATPase viral hallmark must
            # lie inside the floored span (blocks host regions whose >=2 validated
            # markers are sparse / host-like by TIER-2).
            non_atpase_hallmark = sum(1 for g in (getattr(alt, "hallmark_genes", []) or []) if not _is_atpase_marker(g))
            if non_atpase_hallmark < 1:
                continue
            # (c) TIER-2 viral-content floor: the floored span's interior gene
            # content must be genuinely viral, not merely carry a sparse marker.
            # A validated marker is a TIER-1 (HMM/Diamond) hit; host regions whose
            # >=2 validated markers are host-like under the TIER-2 proteome
            # classifier have near-zero interior viral fraction. Requiring
            # >=READMIT_MIN_VIRAL_FRACTION admits genuine marker-dense NCLDV EVEs
            # while rejecting host-with-a-marker regions. Measured on three real
            # genomes with v1.0.6 resources: every region this pass admitted
            # carried viral best-hit proteins (none were zero-viral), and no
            # previously accepted region was lost. The 0.10 value itself has not
            # been swept, so it remains a working default rather than a tuned
            # optimum.
            total_interior = getattr(alt, "gene_taxonomy_total", 0) or 0
            viral_interior = getattr(alt, "gene_taxonomy_viral_top10", 0) or 0
            alt_viral_fraction = viral_interior / total_interior if total_interior else 0.0
            if alt_viral_fraction < READMIT_MIN_VIRAL_FRACTION:
                continue
            # (d) dedup: never re-admit something overlapping an already-accepted
            # region; readmit_accepted grows as we re-admit, so two alternatives
            # cannot overlap each other either.
            if any(a.scaffold == alt.scaffold and a.start < alt.end and a.end > alt.start for a in readmit_accepted):
                continue
            readmit_accepted.append(alt)
            # Also surface the alternative as a verified candidate. The second
            # selection pass below produces the explicit canonical subset used
            # by output generation. Only fully qualified alternatives reach it.
            verification_results.append(alt)
            n_readmitted += 1
        if n_readmitted:
            # The rejected originals stay in verification_results; the qualifying
            # alternatives are appended alongside them. Every count below is
            # recomputed from the second selection, so no manual bookkeeping.
            logger.info(
                "Phase 3 marker-floor re-admit: recovered %d marker-dense EVE(s) "
                "excluded after host-trim boundary collapse",
                n_readmitted,
            )

        if n_readmitted:
            # Re-run selection over the extended candidate list so the counts,
            # classification statistics, promoted-LOW list, and ablation
            # bookkeeping below all derive from one consistent selection.
            # Appending to accepted_results after the fact would leave
            # normal_gate_decisions and the classified-vs-accepted invariant
            # out of step with the accepted set.
            acceptance_selection = select_phase3_acceptance(
                verification_results,
                ablation_id,
            )

    ablation_counts = acceptance_selection.intervention_counts
    if ablation_id is AblationID.A5:
        effects = [result.composition_ablation_effect for result in verification_results]
        ablation_counts = InterventionCounts(
            opportunities=sum(effect.opportunities for effect in effects),
            interventions=sum(effect.interventions for effect in effects),
            changed=sum(effect.changed for effect in effects),
        )
    accepted_results = list(acceptance_selection.canonical_results)
    # Cluster the accepted EVEs before any class is counted or persisted: the
    # manifest normalizes the persisted effective_eve_class through the contract
    # and the orchestrator raises when persisted and in-memory counts disagree,
    # so a class rewrite after this point would fail the run.
    _edges_path, ani_pairs = cluster_accepted_eves(
        accepted_results,
        genome_fasta=masked_path,
        output_dir=output_dir,
        threads=threads,
    )
    # Then drop the EVEs with no viral evidence at all that no marker-bearing
    # relative vouches for. This needs the clusters, and it is the one place
    # Phase 3 changes the accepted set after the gate has run.
    gate_kept = len(accepted_results)
    unsupported = unsupported_eve_ids(accepted_results)
    if unsupported:
        for result in accepted_results:
            if result.eve_id in unsupported:
                result.canonical_selection_outcome = "unsupported_no_viral_evidence"
        accepted_results = [result for result in accepted_results if result.eve_id not in unsupported]
        # Cluster sizes counted the dropped members; recount over the survivors
        # so no published row claims a relative that is not published.
        recluster_survivors(accepted_results, ani_pairs)
    # Count only what the gate itself rejected. The unsupported-EVE removal is
    # reported separately; folding it in here would blame the gate for it.
    quality_gate_dropped = len(verification_results) - gate_kept
    counterfactual_quality_gate_dropped = sum(
        not decision.kept for decision in acceptance_selection.normal_gate_decisions
    )
    accepted = len(accepted_results)
    _annotate_integration_genes(
        verification_results,
        proteome_path=proteome_path,
        validated_markers=validated_markers,
        output_dir=output_dir,
        config=config,
    )
    phase3_elapsed = time.time() - phase3_start

    # Compute confidence tier distributions for all candidates and the
    # selected canonical surface. Under A6, canonical LOW rows can include
    # candidates that the normal gate rejected.
    candidate_tier_counts = {"HIGH": 0, "MEDIUM": 0, "LOW": 0}
    for r in verification_results:
        tier = getattr(r, "confidence_tier", "LOW")
        if tier in candidate_tier_counts:
            candidate_tier_counts[tier] += 1
    tier_counts = {"HIGH": 0, "MEDIUM": 0, "LOW": 0}
    for r in accepted_results:
        tier = getattr(r, "confidence_tier", "LOW")
        if tier in tier_counts:
            tier_counts[tier] += 1

    # Compute classification statistics for canonical accepted regions.
    classification_stats = {eve_class: 0 for eve_class in EFFECTIVE_EVE_CLASSES}
    accepted_bp = 0
    total_genes = 0
    total_hallmarks = 0

    # Iterate the published list, not acceptance_selection.candidates: the
    # unsupported-EVE drop above can remove a candidate the gate kept, and these
    # counts must match what is written out.
    for r in accepted_results:
        accepted_bp += r.end - r.start
        total_genes += getattr(r, "gene_count", 0)
        total_hallmarks += getattr(r, "hallmark_count", 0)
        # Published class, not the gate decision: acceptance is the gate's job,
        # the label is the taxonomy consensus over the region's own markers.
        cls = normalize_effective_eve_class(getattr(r, "taxonomy_class", "UNKNOWN"))
        classification_stats[cls] += 1

    classified = sum(classification_stats[key] for key in EFFECTIVE_EVE_CLASSES)
    if classified != accepted:
        raise RuntimeError(
            "exclusive effective-class counts do not sum to accepted predictions: "
            f"accepted={accepted} classified={classified}"
        )

    # Log detailed results
    if verification_results:
        logger.info(
            f"Phase 3 complete: {phase3_elapsed:.1f}s, "
            f"verified {len(verification_results)} candidates: "
            f"HIGH={candidate_tier_counts['HIGH']}, "
            f"MEDIUM={candidate_tier_counts['MEDIUM']}, "
            f"LOW={candidate_tier_counts['LOW']}"
        )
        logger.info(
            "  %s kept %d/%d predictions "
            "(HIGH=%d, MEDIUM=%d, LOW=%d; selection-dropped=%d; "
            "no-viral-evidence-dropped=%d; normal-gate-rejected=%d)",
            "A6 acceptance bypass" if ablation_id is AblationID.A6 else "Canonical v2 gate",
            accepted,
            len(verification_results),
            tier_counts["HIGH"],
            tier_counts["MEDIUM"],
            tier_counts["LOW"],
            quality_gate_dropped,
            len(unsupported),
            counterfactual_quality_gate_dropped,
        )
        if accepted_results:
            logger.info(f"  Total regions: {accepted_bp:,} bp, {total_genes} genes")
            # Log classification breakdown
            cls_parts = [f"{k}={v}" for k, v in classification_stats.items() if v > 0]
            if cls_parts:
                logger.info(f"  Classifications: {', '.join(cls_parts)}")
    else:
        logger.info(f"Phase 3 complete: {phase3_elapsed:.1f}s, no candidates to verify")

    return Phase3Result(
        verification_results=verification_results,
        accepted_results=accepted_results,
        # Must stay a subset of accepted_results: a dropped EVE is no longer
        # canonical, so it cannot remain in the promoted-LOW subset either.
        promoted_low_results=[
            result for result in acceptance_selection.promoted_low_results if result.eve_id not in unsupported
        ],
        classification_stats=classification_stats,
        accepted=accepted,
        ablation_counts=ablation_counts,
        elapsed=phase3_elapsed,
    )


def _annotate_integration_genes(
    results: list[VerificationResult],
    *,
    proteome_path: Path,
    validated_markers: list[ValidatedMarkerHit],
    output_dir: Path,
    config: PipelineConfig,
) -> None:
    """Attach integration annotations without modifying coordinates or confidence."""
    from virosync.pipeline.phase3.integration_genes import (
        IntegrationRegion,
        collect_annotated_integration_genes,
        scan_integration_genes,
    )

    if not results:
        return
    flank_bp = config.phase1.extension_kb * 1000
    regions = [
        IntegrationRegion(
            region_id=result.eve_id,
            scaffold=result.scaffold,
            scan_start=max(
                0,
                min(result.start, result.candidate_start if result.candidate_start is not None else result.start)
                - flank_bp,
            ),
            scan_end=max(result.end, result.candidate_end if result.candidate_end is not None else result.end)
            + flank_bp,
            interior_start=result.start,
            interior_end=result.end,
        )
        for result in results
    ]
    hits = scan_integration_genes(proteome_path, regions, threads=config.compute.effective_threads())
    model_annotations = (
        Path(config.databases.hmm_database).parent / "model_annotations_with_interpro.tsv"
        if config.databases.hmm_database
        else None
    )
    interproscan_path = output_dir / "phase3" / "interproscan" / "interproscan_batch.tsv"
    hits.extend(
        collect_annotated_integration_genes(
            proteome_path,
            validated_markers,
            regions,
            model_annotations_path=model_annotations,
            interproscan_path=interproscan_path if config.phase3.interproscan_enabled else None,
        )
    )
    by_region: dict[str, list[dict[str, object]]] = defaultdict(list)
    for hit in hits:
        by_region[hit.region_id].append(asdict(hit))
    for result in results:
        result.integration_gene_hits = by_region[result.eve_id]
