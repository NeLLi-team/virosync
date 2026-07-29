"""Phase 3 subflow: evidence synthesis, verification, tiering."""

import json
from dataclasses import replace
from pathlib import Path
from typing import Optional

from virosync.ablation import AblationID, InterventionCounts
from virosync.output_contract import (
    EFFECTIVE_EVE_CLASSES,
    normalize_effective_eve_class,
)
from virosync.pipeline.phase3.acceptance_selection import (
    select_phase3_acceptance,
)
from virosync.pipeline.phase3.output_generator import (
    _is_atpase_marker,
    evaluate_v2_quality_gate,
)
from virosync.pipeline.phase1.marker_roles import decide_marker_hit_role
from virosync.pipeline.phase1.viral_markers import get_assembly_mode
from virosync.orchestration.runtime import call_task
from virosync.pipeline.phase2.boundary_diamond import (
    filter_taxonomy_to_boundary,
    get_flanking_taxonomy,
    build_gene_taxonomy_record,
)
from virosync.orchestration.utils import get_genes_for_boundary

from .loaders import (
    _load_interproscan_summary,
    _load_tmvec_cache,
    _serialize_tmvec_cache,
)


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


def _run_phase3_subflow(
    # Core inputs from Phase 0
    masked_path: Path,
    proteome_path: Path,
    # Core inputs from Phase 1
    validated_markers: list,
    host_signatures: set,
    host_signature_model,
    host_signature_model_payload: Optional[dict],
    # Core inputs from Phase 2
    refined_boundaries: list,
    boundary_taxonomy_map: dict,
    boundary_control_stats,
    boundary_diamond_query,
    proteome_index: dict,
    boundaries_bed: Path,
    merged_seeds: list,
    # Core identifiers
    output_dir: Path,
    genome_id: str,
    # Resume configuration
    resume: bool,
    validated_hits_tsv: Path,
    # Database parameters
    gene_taxonomy_faa_db: Optional[Path],
    marker_db: Optional[Path],
    marker_faa_db: Optional[Path],
    marker_faa_dir: Optional[Path],
    faa_dir: Optional[Path],
    diamond_db: Optional[Path],
    gvclass_db: Optional[Path],
    hmm_database: Optional[Path],
    viral_structure_db: Optional[Path],
    taxonomy_labels_file: Optional[Path],
    # Host configuration
    host_prefixes: list[str],
    host_label: str,
    high_pident_host_threshold: float,
    boundary_host_trim_score_threshold: float,
    host_signature_evidence_threshold: float,
    boundary_diamond_flank_genes: int,
    # Verification parameters
    high_tier_threshold: float,
    low_tier_threshold: float,
    use_crf_in_final_score: bool,
    priority_marker_list: Optional[list[str]],
    marker_floor_priority_only: float,
    marker_floor_priority_plus_family: float,
    marker_floor_priority_multi_family: float,
    marker_family_bonus_per_family: float,
    marker_multi_family_bonus: float,
    enable_phylogenetic: bool,
    # Structural analysis parameters
    skip_structural: bool,
    use_boltz: bool,
    boltz_mcp_only: bool,
    boltz_use_msa_server: bool,
    boltz_min_seq_len: int,
    boltz_max_seq_len: int,
    boltz_no_kernels: bool,
    use_tmvec_database: bool,
    tmvec_require_gpu: bool,
    tmvec_databases: Optional[list[str]],
    tmvec_database_dir: Optional[Path],
    tmvec_min_score: float,
    device: str,
    # InterProScan parameters
    interproscan_enabled: bool,
    interproscan_dir: Optional[Path],
    interproscan_keywords: Optional[list[str]],
    interproscan_threads: Optional[int],
    interproscan_applications: Optional[list[str]],
    # Database rebuild (for Phase 3 fallback)
    rebuild_db: bool,
    # Threading
    threads: int,
    gene_taxonomy_threads: Optional[int],
    # Logger
    logger,
    # Set only after schema-v3 marker validation by the orchestrator.
    resume_authorized: bool = False,
    ablation_id: AblationID = AblationID.A0,
    assembly_mode: str = "default",
) -> dict:
    """
    Phase 3: Evidence synthesis and verification.

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
        ... (see function signature for all parameters)
        logger: Logger instance

    Returns:
        dict with keys:
            - verification_results: List of VerificationResult objects
            - gene_taxonomy_map_batched: Dict of gene taxonomy by boundary
            - hallmark_hits_map: Dict of hallmark hits by boundary
            - interproscan_map_batched: Dict of InterProScan results by boundary
            - precomputed_tmvec: Dict of TMVec results
            - tier_counts: Dict with HIGH/MEDIUM/LOW counts
            - classification_stats: Dict with classification counts
            - accepted: Count of accepted predictions
            - accepted_bp: Total basepairs of accepted regions
            - total_genes: Total gene count
            - total_hallmarks: Total hallmark count
            - elapsed: Phase 3 elapsed time in seconds
    """
    import time

    phase3_start = time.time()
    logger.info("-" * 60)
    logger.info(f"Phase 3: Verifying {len(refined_boundaries)} candidates")
    logger.info("  Steps: gene taxonomy -> evidence synthesis -> structural (if enabled)")

    # Load validated markers on resume if needed
    if (
        (not validated_markers)
        and resume
        and resume_authorized
        and validated_hits_tsv.exists()
    ):
        from virosync.pipeline.phase1.marker_validation import (
            load_validated_marker_hits,
            collect_host_signatures,
        )
        from virosync.pipeline.host_signatures import HostSignatureModel

        validated_markers = load_validated_marker_hits(validated_hits_tsv)
        host_signatures = collect_host_signatures(
            validated_markers,
            host_prefixes=set(host_prefixes),
        )
        model_path = output_dir / "phase1" / "marker_validation" / "host_signature_model.json"
        if model_path.exists():
            with model_path.open() as handle:
                host_signature_model_payload = json.load(handle)
                host_signature_model = HostSignatureModel.from_dict(host_signature_model_payload)
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
    gene_taxonomy_map = None
    interproscan_map = None
    interproscan_summary_path = output_dir / "phase3" / "interproscan" / "interproscan_summary.tsv"

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
        logger.info(
            "A3: gene taxonomy unavailable because Tier-2 search is bypassed"
        )
    elif merged_seeds:
        # Outside A3, Phase 2b is mandatory when seeds exist.
        logger.error(
            "Phase 3: No taxonomy data available! Phase 2b should have run. "
            "Gene taxonomy scores will be unavailable."
        )

    if interproscan_enabled:
        if resume and resume_authorized and interproscan_summary_path.exists():
            try:
                interproscan_map = _load_interproscan_summary(interproscan_summary_path)
                logger.info(
                    "Phase 3 resume: loaded InterProScan summary for %d candidates",
                    len(interproscan_map),
                )
            except Exception as exc:
                logger.warning(
                    "Phase 3 resume: failed loading InterProScan summary (%s); rerunning",
                    exc,
                )
                interproscan_map = None

        if interproscan_map is not None:
            pass
        elif not interproscan_dir:
            logger.warning("Phase 3: InterProScan enabled but interproscan_dir not set; skipping")
        elif not Path(interproscan_dir).exists():
            logger.warning("Phase 3: InterProScan dir not found: %s; skipping", interproscan_dir)
        else:
            from virosync.orchestration.tasks import interproscan_batch_task

            interpro_threads = interproscan_threads if interproscan_threads is not None else threads
            interpro_threads = max(1, min(interpro_threads, threads))
            logger.info(
                "Phase 3: InterProScan threads=%s (batch)",
                interpro_threads,
            )
            try:
                interproscan_map = call_task(
                    interproscan_batch_task,
                    regions=regions_payload,
                    proteome_path=proteome_path,
                    interproscan_dir=Path(interproscan_dir),
                    output_dir=output_dir / "phase3" / "interproscan",
                    threads=interpro_threads,
                    keywords=interproscan_keywords,
                    applications=interproscan_applications,
                )
            except Exception as exc:
                logger.warning("Phase 3: batch InterProScan failed: %s", exc)
                interproscan_map = None

    # ==================================================
    # Batch TMVec precomputation for ALL EVE proteins
    # ==================================================
    tmvec_device = _resolve_tmvec_device(device) if use_tmvec_database else device

    precomputed_tmvec = None
    tmvec_cache_path = output_dir / "phase3" / "tmvec" / "precomputed_tmvec.json"
    if (
        use_tmvec_database
        and resume
        and resume_authorized
        and tmvec_cache_path.exists()
    ):
        try:
            precomputed_tmvec = _load_tmvec_cache(tmvec_cache_path)
            logger.info(
                "Phase 3 resume: loaded cached TMVec hits for %d proteins",
                len(precomputed_tmvec),
            )
        except Exception as exc:
            logger.warning(
                "Phase 3 resume: failed loading TMVec cache (%s); recomputing",
                exc,
            )
            precomputed_tmvec = None

    if use_tmvec_database and refined_boundaries and precomputed_tmvec is None:
        # Collect all proteins from all EVEs and deduplicate by pORF ID.
        raw_protein_count = 0
        conflicting_ids = 0
        protein_by_id: dict[str, str] = {}
        for boundary in refined_boundaries:
            porf_sequences = get_genes_for_boundary(
                proteome_path=Path(proteome_path),
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

        if all_proteins:
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
                    databases=tmvec_databases or ["bfvd"],
                    min_tm=0.0,  # Record all hits; scoring threshold applied later
                    database_root=tmvec_database_dir,
                    require_gpu=tmvec_require_gpu,
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
                    per_seq_failures = getattr(
                        predictor, "_per_seq_fallback_failures", 0
                    )
                    if oom or per_seq_failures:
                        logger.warning(
                            "Phase 3: TMVec GPU degradation detected — "
                            "batch_oom_fallbacks=%d per_seq_fallback_failures=%d",
                            oom,
                            per_seq_failures,
                        )
                if precomputed_tmvec:
                    _serialize_tmvec_cache(precomputed_tmvec, tmvec_cache_path)
                    logger.info("Phase 3: TMVec cache written to %s", tmvec_cache_path)
            except Exception as exc:
                logger.error(
                    "Phase 3: TMVec failed after preflight enabled it; "
                    "refusing silent fallback: %s",
                    exc,
                )
                raise

    # ==================================================
    # Pre-score MCP DJR/SJR classification (feeds confidence bonus)
    # ==================================================
    from virosync.pipeline.phase3.evidence_synthesizer import load_jelly_roll_data

    jelly_roll_map: dict[str, list[dict]] = {}
    jelly_roll_output = output_dir / "phase3_synthesis" / "virosync_jelly_roll_proteins.tsv"
    marker_hits_path = Path(validated_hits_tsv)
    marker_sequences_path = output_dir / "phase1" / "marker_validation" / "hmm_hit_porfs.faa"
    interproscan_batch_path = output_dir / "phase3" / "interproscan" / "interproscan_batch.tsv"
    foldseek_results_path = output_dir / "structural_analysis" / "foldseek_pdb_results.tsv"

    if resume and resume_authorized and jelly_roll_output.exists():
        jelly_roll_map = load_jelly_roll_data(jelly_roll_output)
        if jelly_roll_map:
            logger.info(
                "Phase 3 resume: loaded jelly-roll classifications for %d MCP proteins",
                len(jelly_roll_map),
            )

    if not jelly_roll_map and marker_hits_path.exists() and marker_sequences_path.exists():
        from virosync.orchestration.tasks import classify_jelly_roll_task

        logger.info("Phase 3: classifying MCP proteins as DJR/SJR before confidence scoring")
        call_task(
            classify_jelly_roll_task,
            marker_hits_path=marker_hits_path,
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
    elif not jelly_roll_map:
        logger.info(
            "Phase 3: skipping jelly-roll classification before scoring "
            "(missing marker validation files)"
        )

    # ==================================================
    # Build maps for batched verification
    # ==================================================
    hallmark_hits_map = {}
    gene_taxonomy_map_batched = {}
    interproscan_map_batched = {}

    def build_boundary_evidence(boundary):
        """Build (boundary_id, hallmarks, gene_tax_result, interproscan_result)
        for one boundary from Phase-2b precomputed taxonomy + validated markers.

        Factored out of the main verification loop so the Phase-3 marker-floor
        re-admit pass can synthesize a REJECTED boundary's floored alternative the
        SAME way the normal pipeline does (.
        """
        boundary_id = f"{boundary.scaffold}_{boundary.start}_{boundary.end}"
        gene_tax_result = None
        interproscan_result = None

        # Get gene taxonomy for this boundary
        if use_precomputed_taxonomy:
            # Use pre-computed taxonomy from Phase 2b (filter to refined boundary)
            filtered_taxonomy = filter_taxonomy_to_boundary(boundary_taxonomy_map, boundary)
            if filtered_taxonomy:
                # Convert GeneTaxonomy objects to the format expected by verify_eve_task
                # The expected format is: (list of records, summary dict)
                gene_tax_records = []
                n_viral_interior = 0
                n_ncldv_mirus_interior = 0
                n_vp_plv_interior = 0
                n_host = 0
                host_prefix = f"{host_label}__"
                for tax in filtered_taxonomy:
                    gene_tax_records.append(
                        build_gene_taxonomy_record(tax, is_flanking=False)
                    )
                    if tax.has_viral:
                        n_viral_interior += 1
                    if tax.has_ncldv_mirus:
                        n_ncldv_mirus_interior += 1
                    if tax.has_vp_plv:
                        n_vp_plv_interior += 1
                    if tax.top1_prefix == host_prefix:
                        n_host += 1

                # Get flanking gene taxonomy (+/-20 genes outside EVE boundary)
                # Use pre-computed seed mapping if available for reliable retrieval
                seed_mapping = None
                if boundary_diamond_query and boundary.seed_id:
                    seed_mapping = boundary_diamond_query.seed_gene_mappings.get(boundary.seed_id)
                upstream_tax, downstream_tax = get_flanking_taxonomy(
                    taxonomy_map=boundary_taxonomy_map,
                    proteome_index=proteome_index,
                    refined_boundary=boundary,
                    flank_genes=boundary_diamond_flank_genes,
                    seed_mapping=seed_mapping,
                )
                n_flanking = 0
                n_viral_flanking = 0
                n_ncldv_mirus_flanking = 0
                n_vp_plv_flanking = 0
                for tax in upstream_tax:
                    gene_tax_records.append(
                        build_gene_taxonomy_record(
                            tax, is_flanking=True, flank_position="upstream"
                        )
                    )
                    n_flanking += 1
                    # Count viral genes in flanking regions
                    if tax.has_viral:
                        n_viral_flanking += 1
                    if tax.has_ncldv_mirus:
                        n_ncldv_mirus_flanking += 1
                    if tax.has_vp_plv:
                        n_vp_plv_flanking += 1
                for tax in downstream_tax:
                    gene_tax_records.append(
                        build_gene_taxonomy_record(
                            tax, is_flanking=True, flank_position="downstream"
                        )
                    )
                    n_flanking += 1
                    # Count viral genes in flanking regions
                    if tax.has_viral:
                        n_viral_flanking += 1
                    if tax.has_ncldv_mirus:
                        n_ncldv_mirus_flanking += 1
                    if tax.has_vp_plv:
                        n_vp_plv_flanking += 1

                # Combine interior and flanking counts for diagnostics only.
                # Core confidence/classification fields below use interior genes only.
                n_viral_total = n_viral_interior + n_viral_flanking
                n_ncldv_mirus_total = n_ncldv_mirus_interior + n_ncldv_mirus_flanking
                n_vp_plv_total = n_vp_plv_interior + n_vp_plv_flanking

                # Calculate dominant_family based on NCLDV/MIRUS/VP/PLV counts in top10
                # Use interior genes only to avoid flanking contamination.
                # boundary_diamond stores prefixes with trailing underscores.
                all_genes_for_family = list(filtered_taxonomy)
                family_counts = {
                    "NCLDV": sum(
                        1 for t in all_genes_for_family
                        if any(p in {"NCLDV__", "NCLDV"} for p in (t.top10_prefixes or []))
                    ),
                    "MIRUS": sum(
                        1 for t in all_genes_for_family
                        if any(p in {"MIRUS__", "MIRUS"} for p in (t.top10_prefixes or []))
                    ),
                    "VP": sum(
                        1 for t in all_genes_for_family
                        if any(p in {"VP__", "VP"} for p in (t.top10_prefixes or []))
                    ),
                    "PLV": sum(
                        1 for t in all_genes_for_family
                        if any(p in {"PLV__", "PLV"} for p in (t.top10_prefixes or []))
                    ),
                    "PPV": sum(
                        1 for t in all_genes_for_family
                        if any(p in {"PPV__", "PPV"} for p in (t.top10_prefixes or []))
                    ),
                }
                dominant_family = "UNKNOWN"
                dominant_fraction = 0.0
                total_genes_for_family = len(all_genes_for_family)
                if total_genes_for_family > 0:
                    dominant_family = max(family_counts, key=family_counts.get)
                    max_count = family_counts[dominant_family]
                    if max_count > 0:
                        dominant_fraction = max_count / total_genes_for_family
                    else:
                        dominant_family = "UNKNOWN"

                gene_tax_summary = {
                    "total": len(filtered_taxonomy),  # Interior only for gene_count
                    "total_with_flanking": len(gene_tax_records),  # Interior + flanking
                    "flanking_genes": n_flanking,
                    "ncldv_mirus": n_ncldv_mirus_interior,  # Interior only
                    "vp_plv": n_vp_plv_interior,  # Interior only
                    "viral_top10": n_viral_interior,  # Interior only
                    "high_pident_euk": n_host,  # Interior only (host penalty)
                    "has_ncldv_mirus": n_ncldv_mirus_interior > 0,
                    "has_vp_plv": n_vp_plv_interior > 0,
                    "dominant_family": dominant_family,
                    "dominant_fraction": dominant_fraction,
                    # NEW: Debugging fields to track interior vs flanking contributions
                    "viral_interior": n_viral_interior,
                    "viral_flanking": n_viral_flanking,
                    "ncldv_mirus_interior": n_ncldv_mirus_interior,
                    "ncldv_mirus_flanking": n_ncldv_mirus_flanking,
                    "vp_plv_interior": n_vp_plv_interior,
                    "vp_plv_flanking": n_vp_plv_flanking,
                }
                gene_tax_result = (gene_tax_records, gene_tax_summary)
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
        # Fall back to Phase 3 batch Diamond if Phase 2b didn't cover this boundary
        if gene_tax_result is None and gene_taxonomy_map:
            # Use Phase 3 batch Diamond result (fallback or primary)
            gene_tax_result = gene_taxonomy_map.get(
                f"EVE_{boundary.scaffold}_{boundary.start}-{boundary.end}"
            )
        if gene_tax_result is None and ablation_id is AblationID.A3:
            porf_sequences = get_genes_for_boundary(
                proteome_path=Path(proteome_path),
                scaffold=boundary.scaffold,
                start=boundary.start,
                end=boundary.end,
                max_porfs=10000,
            )
            gene_count = len(porf_sequences)
            gene_tax_result = (
                [],
                {
                    "total": gene_count,
                    "total_with_flanking": gene_count,
                    "flanking_genes": 0,
                },
            )
        if interproscan_map:
            interproscan_result = interproscan_map.get(
                f"EVE_{boundary.scaffold}_{boundary.start}-{boundary.end}"
            )

        boundary_hallmarks = []
        single_marker_min_score = get_assembly_mode(
            assembly_mode
        ).single_marker_min_score
        for marker in validated_markers:
            if marker.scaffold != boundary.scaffold:
                continue
            if marker.start >= boundary.end or marker.end <= boundary.start:
                continue
            marker_role = decide_marker_hit_role(
                marker,
                ablation_id=ablation_id,
                single_marker_min_score=single_marker_min_score,
            )
            if not marker_role.is_retained_evidence:
                continue
            porf_id = getattr(marker, "porf_id", None) or getattr(marker, "query_porf", None)
            boundary_hallmarks.append(
                {
                    "start": marker.start,
                    "end": marker.end,
                    "porf_id": porf_id,
                    "hallmark_gene": marker.hmm_target,
                    "hmm_score": marker.hmm_score,
                    "score": marker.hmm_score,
                    "top10_prefixes": getattr(marker, "top10_prefixes", ""),
                    "validation_status": marker_role.original_validation_status,
                    "tier1_bypassed": marker_role.is_tier1_bypassed,
                }
            )

        return boundary_id, boundary_hallmarks, gene_tax_result, interproscan_result

    from virosync.orchestration.tasks import verify_eve_candidates_batched_task

    def run_verification(boundary_list, hh_map, gt_map, ip_map):
        """Run batched Phase-3 verification for a set of boundaries + their
        evidence maps. Shared by the main pass and the marker-floor re-admit pass,
        so an alternative boundary is scored through the same synthesis path
        (hallmark_count, non_atpase_hallmark, has_mcp, viral_fraction, confidence,
        tier, eve_class).
        """
        return call_task(
            verify_eve_candidates_batched_task,
            boundaries=boundary_list,
            genome_path=masked_path,
            work_dir=output_dir / "phase3",
            proteome_path=proteome_path,
            hallmark_hits_map=hh_map,
            novelty_scores={},
            gene_taxonomy_map=gt_map if gt_map else None,
            interproscan_map=ip_map if ip_map else None,
            jelly_roll_map=jelly_roll_map if jelly_roll_map else None,
            euk_host_signatures=host_signatures,
            host_signature_model=host_signature_model_payload,
            host_signature_score_threshold=host_signature_evidence_threshold,
            host_prefixes=host_prefixes,
            host_label=host_label,
            high_tier_threshold=high_tier_threshold,
            low_tier_threshold=low_tier_threshold,
            use_crf_in_final_score=use_crf_in_final_score,
            priority_marker_list=priority_marker_list,
            marker_floor_priority_only=marker_floor_priority_only,
            marker_floor_priority_plus_family=marker_floor_priority_plus_family,
            marker_floor_priority_multi_family=marker_floor_priority_multi_family,
            marker_family_bonus_per_family=marker_family_bonus_per_family,
            marker_multi_family_bonus=marker_multi_family_bonus,
            skip_structural=skip_structural,
            use_boltz=use_boltz,
            boltz_mcp_only=boltz_mcp_only,
            boltz_use_msa_server=boltz_use_msa_server,
            boltz_min_seq_len=boltz_min_seq_len,
            boltz_max_seq_len=boltz_max_seq_len,
            boltz_no_kernels=boltz_no_kernels,
            use_tmvec_database=use_tmvec_database,
            tmvec_databases=tmvec_databases,
            tmvec_database_dir=tmvec_database_dir,
            tmvec_min_score=tmvec_min_score,
            tmvec_require_gpu=tmvec_require_gpu,
            device=tmvec_device,
            viral_structure_db=viral_structure_db,
            gvclass_db=gvclass_db,
            diamond_db=diamond_db,
            enable_phylogenetic=enable_phylogenetic,
            taxonomy_labels_file=Path(taxonomy_labels_file) if taxonomy_labels_file else None,
            hmm_database=hmm_database,
            precomputed_tmvec=precomputed_tmvec,
            max_workers=threads,
            ablation_id=ablation_id,
        )

    logger.info(
        "Phase 3: Preparing data for batched verification of %d boundaries",
        len(refined_boundaries),
    )
    for boundary in refined_boundaries:
        boundary_id, boundary_hallmarks, gene_tax_result, interproscan_result = (
            build_boundary_evidence(boundary)
        )
        hallmark_hits_map[boundary_id] = boundary_hallmarks
        if gene_tax_result:
            gene_taxonomy_map_batched[boundary_id] = gene_tax_result
        if interproscan_result:
            interproscan_map_batched[boundary_id] = interproscan_result

    # ==================================================
    # Run batched verification
    # ==================================================
    logger.info("Phase 3: Running batched verification with %d threads", threads)
    verification_results = run_verification(
        refined_boundaries,
        hallmark_hits_map,
        gene_taxonomy_map_batched,
        interproscan_map_batched,
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
    # extended list; it evaluates each candidate independently, so the originals
    # receive the same decisions they got in the first pass.
    #
    READMIT_MIN_VIRAL_FRACTION = 0.10
    accepted_ids = {id(r) for r in acceptance_selection.canonical_results}
    # Regions accepted so far, seeded from the first pass and extended as
    # alternatives qualify, so no two re-admits can overlap each other.
    readmit_accepted = list(acceptance_selection.canonical_results)
    boundary_by_region = {(b.scaffold, b.start, b.end): b for b in refined_boundaries}
    readmit_boundaries = []
    for r in verification_results:
        if id(r) in accepted_ids:
            continue  # never touch / re-evaluate an accepted region
        boundary = boundary_by_region.get((r.scaffold, r.start, r.end))
        if boundary is None:
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
        # Evidence parity is partial by construction. Hallmarks and Phase-2b gene
        # taxonomy are recomputed for the floored span, but InterProScan and TMVec
        # were precomputed for the ORIGINAL boundaries only, so proteins that the
        # floor newly includes carry no InterPro or structural evidence here.
        # That can only lower an alternative's score, so it costs recoveries
        # rather than admitting anything extra. Re-running those two tools per
        # alternative is the fix if the recovery rate proves too low.
        alt_hh_map, alt_gt_map, alt_ip_map = {}, {}, {}
        for alt_boundary in readmit_boundaries:
            bid, b_hallmarks, gt_result, ip_result = build_boundary_evidence(alt_boundary)
            alt_hh_map[bid] = b_hallmarks
            if gt_result:
                alt_gt_map[bid] = gt_result
            if ip_result:
                alt_ip_map[bid] = ip_result
        alt_results = run_verification(
            readmit_boundaries, alt_hh_map, alt_gt_map, alt_ip_map
        )
        # Deterministic order so the highest-confidence alternative wins any overlap.
        alt_results.sort(
            key=lambda a: (-getattr(a, "final_confidence", 0.0), a.scaffold, a.start)
        )
        for alt in alt_results:
            # (a) the floored alternative must independently pass the v2 gate.
            if not evaluate_v2_quality_gate(alt).kept:
                continue
            # (b) viral-content criterion: a genuine non-ATPase viral hallmark must
            # lie inside the floored span (blocks host regions whose >=2 validated
            # markers are sparse / host-like by TIER-2).
            non_atpase_hallmark = sum(
                1 for g in (getattr(alt, "hallmark_genes", []) or [])
                if not _is_atpase_marker(g)
            )
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
            if any(
                a.scaffold == alt.scaffold and a.start < alt.end and a.end > alt.start
                for a in readmit_accepted
            ):
                continue
            readmit_accepted.append(alt)
            # Also surface the alternative as a verified candidate so it flows
            # through the canonical output path (generate_outputs_task re-applies
            # the v2 gate to verification_results; only gate-passers are written to
            # virosync_predictions.* / gene_taxonomy/). Only fully-qualified alts
            # (gate-pass AND non-ATPase hallmark AND no overlap with an accepted
            # region) reach this point, so the re-gated output set still equals the
            # de-duplicated accepted set.
            verification_results.append(alt)
            n_readmitted += 1
        if n_readmitted:
            # The rejected originals stay in verification_results; the qualifying
            # alternatives are appended alongside them. Every count below is
            # recomputed from the second selection, so no manual bookkeeping.
            logger.info(
                "Phase 3 marker-floor re-admit: recovered %d marker-dense EVE(s) "
                "rejected by the v2 gate after host-trim boundary collapse",
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
    quality_gate_dropped = len(verification_results) - len(accepted_results)
    counterfactual_quality_gate_dropped = sum(
        not decision.kept
        for decision in acceptance_selection.normal_gate_decisions
    )
    accepted = len(accepted_results)
    phase3_elapsed = time.time() - phase3_start

    # Compute confidence tier distributions for all candidates and the
    # selected canonical surface. Under A6, canonical LOW rows can include
    # candidates that the normal gate rejected.
    candidate_tier_counts = {"HIGH": 0, "MEDIUM": 0, "LOW": 0}
    for r in verification_results:
        tier = getattr(r, 'confidence_tier', 'LOW')
        if tier in candidate_tier_counts:
            candidate_tier_counts[tier] += 1
    tier_counts = {"HIGH": 0, "MEDIUM": 0, "LOW": 0}
    for r in accepted_results:
        tier = getattr(r, 'confidence_tier', 'LOW')
        if tier in tier_counts:
            tier_counts[tier] += 1

    # Compute classification statistics for canonical accepted regions.
    classification_stats = {
        "NCLDV": 0,
        "VP": 0,
        "PLV": 0,
        "MIRUS": 0,
        "MIXED": 0,
        "PPV": 0,
        "UNKNOWN": 0,
    }
    accepted_bp = 0
    total_genes = 0
    total_hallmarks = 0

    for candidate in acceptance_selection.candidates:
        if not candidate.canonical_kept:
            continue
        r = candidate.result
        accepted_bp += r.end - r.start
        total_genes += getattr(r, 'gene_count', 0)
        total_hallmarks += getattr(r, 'hallmark_count', 0)
        cls = normalize_effective_eve_class(
            candidate.normal_gate_decision.effective_class
        )
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
            "(HIGH=%d, MEDIUM=%d, LOW=%d; dropped=%d; "
            "normal-gate-rejected=%d)",
            "A6 acceptance bypass" if ablation_id is AblationID.A6 else "Canonical v2 gate",
            accepted,
            len(verification_results),
            tier_counts["HIGH"],
            tier_counts["MEDIUM"],
            tier_counts["LOW"],
            quality_gate_dropped,
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

    return {
        "verification_results": verification_results,
        "accepted_results": accepted_results,
        "promoted_low_results": list(acceptance_selection.promoted_low_results),
        "gene_taxonomy_map_batched": gene_taxonomy_map_batched,
        "hallmark_hits_map": hallmark_hits_map,
        "interproscan_map_batched": interproscan_map_batched,
        "precomputed_tmvec": precomputed_tmvec,
        "tier_counts": tier_counts,
        "candidate_tier_counts": candidate_tier_counts,
        "classification_stats": classification_stats,
        "accepted": accepted,
        "accepted_bp": accepted_bp,
        "total_genes": total_genes,
        "total_hallmarks": total_hallmarks,
        "quality_gate_dropped": quality_gate_dropped,
        "counterfactual_quality_gate_dropped": counterfactual_quality_gate_dropped,
        "ablation_counts": ablation_counts,
        "elapsed": phase3_elapsed,
    }
