"""Phase 2 subflow: boundary refinement (host-trim, taxonomy, Diamond)."""

from dataclasses import dataclass
from pathlib import Path

from virosync.ablation import AblationID, InterventionCounts
from virosync.config import PipelineConfig
from virosync.orchestration.runtime import call_task
from virosync.orchestration.tasks import generate_outputs_task
from virosync.pipeline.phase2.boundary_diamond import (
    BoundaryDiamondConfig,
    ControlStats,
    GeneTaxonomy,
    GenomeDiamondQuery,
    build_proteome_index,
    classify_cached_diamond_query,
    collect_query_proteins,
    compute_control_stats,
    run_batched_diamond,
    run_full_proteome_diamond,
    write_control_stats,
    write_taxonomy_map,
)
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.utils.atomic_write import atomic_write_context

from .manifest import (
    _empty_prediction_summary,
    _write_empty_run_log,
)
from .phase2_resume_state import (
    PHASE2_RESUME_STATE_FILENAME,
    load_phase2_resume_state,
    write_phase2_resume_state,
)
from .phase_state import (
    PHASE2_STATE_FILENAME,
    load_phase2_state,
    phase2_state_to_document,
    write_phase2_state,
)
from .reports import _generate_required_reports
from .resume import _require_phase2b_gene_taxonomy_db


@dataclass(frozen=True, slots=True)
class Phase2Result:
    """Inputs produced for evidence synthesis."""

    refined_boundaries: list[RefinedBoundary]
    boundary_taxonomy_map: dict[str, GeneTaxonomy]
    boundary_control_stats: ControlStats | None
    boundary_diamond_query: GenomeDiamondQuery | None
    proteome_index: dict
    boundaries_bed: Path
    elapsed: float
    goto_phase3: bool = False
    phase_outcome: str = "complete"
    ablation_counts: InterventionCounts = InterventionCounts()


@dataclass(frozen=True, slots=True)
class Phase2Terminal:
    """A complete public flow result produced before Phase 3."""

    result: dict[str, object]


@dataclass(frozen=True, slots=True)
class _BoundaryEvidence:
    """Genome-wide taxonomy evidence used for boundary refinement."""

    taxonomy_map: dict[str, GeneTaxonomy]
    control_stats: ControlStats
    query: GenomeDiamondQuery


def _sum_intervention_counts(
    *counts: InterventionCounts,
) -> InterventionCounts:
    """Sum candidate-level counts from independent hooks in one phase."""
    return InterventionCounts(
        opportunities=sum(item.opportunities for item in counts),
        interventions=sum(item.interventions for item in counts),
        changed=sum(item.changed for item in counts),
    )


def _seeds_to_refined_boundaries(
    merged_seeds: list,
) -> list[RefinedBoundary]:
    """Convert exact seed intervals to Phase-2 boundaries."""
    return [
        RefinedBoundary(
            scaffold=seed.scaffold,
            start=seed.start,
            end=seed.end,
            seed_id=seed.seed_id,
            original_start=seed.start,
            original_end=seed.end,
            confidence=0.0,
            posterior_probability=0.0,
            seed_sources=list(seed.sources),
            seed_confidence=seed.confidence,
            seed_hhg_score=seed.hhg_score,
            seed_novelty_score=seed.novelty_score,
            seed_compositional_score=seed.compositional_score,
            seed_has_mcp=seed.has_mcp,
            gc_deviation=0.0,
            cub_deviation=getattr(seed, "cub_deviation", 0.0),
            max_kfd=0.0,
            predicted_family=getattr(seed, "predicted_family", ""),
            region_classification_ncldv_markers=getattr(
                seed,
                "region_classification_ncldv_markers",
                0,
            ),
            region_classification_vp_plv_markers=getattr(
                seed,
                "region_classification_vp_plv_markers",
                0,
            ),
            region_classification_mirus_markers=getattr(
                seed,
                "region_classification_mirus_markers",
                0,
            ),
            candidate_start=getattr(seed, "host_trim_original_start", None),
            candidate_end=getattr(seed, "host_trim_original_end", None),
            host_trim_reason=getattr(seed, "host_trim_reason", ""),
            host_trim_common_euk_taxonomy=getattr(
                seed,
                "host_trim_common_euk_taxonomy",
                "",
            ),
        )
        for seed in merged_seeds
    ]


def _recalculate_boundary_composition(
    boundaries: list[RefinedBoundary],
    *,
    masked_path: Path,
) -> None:
    """Set GC and k-mer deviations from each boundary's current coordinates."""
    from Bio import SeqIO

    from virosync.features.compositional import (
        BackgroundModel,
        calculate_gc_deviation,
        calculate_kfd,
    )

    if not boundaries:
        return

    scaffold_index = SeqIO.index(str(masked_path), "fasta")
    try:
        background_chunks: list[str] = []
        sampled_bases = 0
        for scaffold_id in scaffold_index:
            seq = str(scaffold_index[scaffold_id].seq)
            if not seq:
                continue
            remaining = 10_000_000 - sampled_bases
            if remaining <= 0:
                break
            background_chunks.append(seq[:remaining])
            sampled_bases += min(len(seq), remaining)
        background_seq = "".join(background_chunks)
        bg_model = BackgroundModel.from_sequence(background_seq, k=4) if len(background_seq) >= 100 else None

        for boundary in boundaries:
            boundary.gc_deviation = 0.0
            boundary.max_kfd = 0.0
            if bg_model is None or boundary.scaffold not in scaffold_index:
                continue
            scaffold_sequence = str(scaffold_index[boundary.scaffold].seq)
            region_sequence = scaffold_sequence[boundary.start : boundary.end]
            if len(region_sequence) < 100:
                continue
            boundary.gc_deviation = calculate_gc_deviation(
                region_sequence,
                bg_model.gc_content,
            )
            boundary.max_kfd = calculate_kfd(
                region_sequence,
                bg_model.kmer_freqs,
                k=4,
            )
    finally:
        scaffold_index.close()


def _write_phase2_checkpoints(
    *,
    output_dir: Path,
    refined_boundaries: list[RefinedBoundary],
    boundary_taxonomy_map: dict,
    boundary_control_stats,
    boundary_diamond_query,
) -> Path:
    """Write the lossless Phase-2 state and its BED report."""
    phase2_dir = output_dir / "phase2"
    phase2_dir.mkdir(parents=True, exist_ok=True)
    boundaries_bed_path = phase2_dir / "refined_boundaries.bed"
    with atomic_write_context(boundaries_bed_path, "w") as handle:
        for boundary in refined_boundaries:
            eve_id = f"EVE_{boundary.scaffold}_{boundary.start}-{boundary.end}"
            score = int(boundary.confidence * 1000)
            handle.write(f"{boundary.scaffold}\t{boundary.start}\t{boundary.end}\t{eve_id}\t{score}\t.\n")
    write_phase2_state(
        phase2_dir / PHASE2_STATE_FILENAME,
        refined_boundaries,
    )
    write_phase2_resume_state(
        phase2_dir / PHASE2_RESUME_STATE_FILENAME,
        refined_boundaries=refined_boundaries,
        boundary_taxonomy_map=boundary_taxonomy_map,
        boundary_control_stats=boundary_control_stats,
        boundary_diamond_query=boundary_diamond_query,
    )
    return boundaries_bed_path


def _run_boundary_diamond(
    *,
    merged_seeds: list,
    proteome_index: dict,
    proteome_path: Path,
    gene_taxonomy_db: Path,
    output_dir: Path,
    config: PipelineConfig,
    superset_diamond_hits: dict | None,
    logger: object,
) -> _BoundaryEvidence:
    """Classify the shared seed, flank, and control protein query."""
    phase2 = config.phase2
    threads = config.compute.gene_taxonomy_threads or config.compute.effective_threads()
    target_control_genes = max(1, phase2.diamond_control_sample_size)
    control_min_distance = max(1, phase2.diamond_control_min_distance)
    logger.info(
        "Phase 2b: Control sampling target: %d controls (min_distance=%d genes)",
        target_control_genes,
        control_min_distance,
    )
    diamond_config = BoundaryDiamondConfig(
        flank_genes=phase2.diamond_flank_genes,
        control_sample_size=target_control_genes,
        control_min_distance=control_min_distance,
        control_region_genes=11,
        top_k=phase2.diamond_top_k,
        chunk_size=phase2.diamond_chunk_size,
        threads=threads,
        random_seed=phase2.diamond_random_seed,
        host_prefix=f"{config.host.label}__",
        taxonomy_weight_mode=phase2.taxonomy_weight_mode,
        search_backend=config.compute.search_backend.value,
    )
    query = collect_query_proteins(
        merged_seeds=merged_seeds,
        proteome_index=proteome_index,
        config=diamond_config,
    )
    logger.info(
        "Phase 2b: Collected %d query proteins (%d EVE regions, %d controls)",
        len(query.all_porf_ids),
        len(query.eve_porf_ids),
        len(query.control_porf_ids),
    )

    diamond_dir = output_dir / "phase2" / "boundary_diamond"
    diamond_dir.mkdir(parents=True, exist_ok=True)
    tax_lookup = None
    taxonomy_labels_file = config.databases.taxonomy_labels_file
    if taxonomy_labels_file and taxonomy_labels_file.exists():
        from virosync.pipeline.host_signatures import TaxonomyLabelLookup

        tax_lookup = TaxonomyLabelLookup.load(taxonomy_labels_file)
        logger.info(
            "Phase 2b: Loaded taxonomy labels for fingerprinting (%d entries)",
            len(tax_lookup),
        )

    if superset_diamond_hits is not None:
        taxonomy_map = classify_cached_diamond_query(
            query=query,
            diamond_hits=superset_diamond_hits,
            proteome_index=proteome_index,
            config=diamond_config,
            taxonomy_lookup=tax_lookup,
        )
    else:
        taxonomy_map = run_batched_diamond(
            query=query,
            proteome_fasta=proteome_path,
            diamond_db=gene_taxonomy_db,
            output_dir=diamond_dir,
            proteome_index=proteome_index,
            config=diamond_config,
            taxonomy_lookup=tax_lookup,
        )
    logger.info("Phase 2b: Diamond complete (%d pORFs classified)", len(taxonomy_map))

    control_taxonomy = [taxonomy_map[pid] for pid in query.control_porf_ids if pid in taxonomy_map]
    control_stats = compute_control_stats(control_taxonomy, diamond_config.host_prefix)
    logger.info(
        "Phase 2b: Control stats: n_genes=%d, host_freq=%.2f, no_hit_freq=%.2f, dominant=%s",
        control_stats.n_genes,
        control_stats.host_frequency,
        control_stats.no_hit_frequency,
        control_stats.dominant_organism,
    )
    write_taxonomy_map(taxonomy_map, diamond_dir / "taxonomy_map.tsv")
    write_control_stats(control_stats, diamond_dir / "control_stats.json")
    return _BoundaryEvidence(
        taxonomy_map=taxonomy_map,
        control_stats=control_stats,
        query=query,
    )


def _run_phase2_subflow(
    masked_path: Path,
    proteome_path: Path,
    merged_seeds: list,
    validated_markers: list,
    host_signature_model,
    output_dir: Path,
    genome_id: str,
    config: PipelineConfig,
    genome_start_time: float,
    logger,
    config_fingerprint: str | None = None,
    resume_authorized: bool = False,
    threads: int | None = None,
    raw_genome_path: Path | None = None,
) -> Phase2Result | Phase2Terminal:
    """Phase 2: Boundary refinement (gene extension, Diamond taxonomy, host trimming).

    This phase refines seed boundaries through:
    - Gene-based seed extension (±5 genes, merge overlapping)
    - Phase 2a: Host-signature trimming (optional)
    - Phase 2b: Batched Diamond for gene taxonomy
    - Phase 2c: Taxonomy-based seed refinement
    - Phase 2f: Host taxonomy trimming (density-aware walk)
    - Post-processing: boundary constraints, adjacent region merging

    Args:
        masked_path: Path to masked genome FASTA
        raw_genome_path: Optional unmasked genome FASTA used for repeat evidence.
        proteome_path: Path to protein FASTA
        merged_seeds: List of MergedSeed objects from Phase 1
        validated_markers: List of validated marker hits from Phase 1
        host_signature_model: HostSignatureModel from Phase 1
        config: Resolved pipeline configuration.
        logger: Logger instance

    Returns:
        Inputs for Phase 3, or a terminal public flow result.
    """
    import time

    databases = config.databases
    compute = config.compute
    host = config.host
    phase2 = config.phase2
    phase3 = config.phase3
    resume = config.execution.resume
    ablation_id = config.ablation.id
    threads = threads if threads is not None else compute.effective_threads()
    gene_taxonomy_faa_db = databases.gene_taxonomy_faa_db
    taxonomy_labels_file = databases.taxonomy_labels_file
    host_label = host.label
    high_pident_host_threshold = host.high_pident_threshold
    boundary_host_trim_enabled = phase2.host_trim_enabled
    boundary_host_trim_window_bp = phase2.host_trim_window_bp
    boundary_host_trim_step_bp = phase2.host_trim_step_bp
    boundary_host_trim_max_host_fraction = phase2.host_trim_max_host_fraction
    boundary_host_trim_min_viral_fraction = phase2.host_trim_min_viral_fraction
    boundary_host_trim_score_threshold = phase2.host_trim_score_threshold
    boundary_host_trim_buffer_kb = phase2.host_trim_buffer_kb
    boundary_host_trim_min_overlap_score = phase2.host_trim_min_overlap_score
    boundary_host_signature_min_token_len = phase2.host_signature_min_token_len
    taxonomy_weight_mode = phase2.taxonomy_weight_mode
    boundary_taxonomy_ml_enabled = phase2.taxonomy_ml_enabled
    boundary_taxonomy_ml_model = phase2.taxonomy_ml_model
    boundary_taxonomy_ml_threshold = phase2.taxonomy_ml_threshold
    boundary_taxonomy_ml_neighbor_window = phase2.taxonomy_ml_neighbor_window
    boundary_diamond_flank_genes = phase2.diamond_flank_genes
    boundary_diamond_top_k = phase2.diamond_top_k
    boundary_diamond_random_seed = phase2.diamond_random_seed
    boundary_diamond_superset_prototype_enabled = phase2.diamond_superset_prototype_enabled
    gene_taxonomy_threads = compute.gene_taxonomy_threads
    extended_output = phase3.extended_output
    search_backend = compute.search_backend.value
    refined_bed = output_dir / "phase2" / "refined_boundaries.bed"

    phase2_start = time.time()
    host_coordinate_counts = InterventionCounts()

    # === PHASE 2a: Host-signature trimming (optional) ===
    phase2_state_path = output_dir / "phase2" / PHASE2_STATE_FILENAME
    phase2_resume_state_path = output_dir / "phase2" / PHASE2_RESUME_STATE_FILENAME
    if ablation_id is AblationID.A3 and merged_seeds and not (resume and resume_authorized):
        logger.info(
            "A3: forwarding %d exact Phase-1 seed intervals to Phase 3",
            len(merged_seeds),
        )
        proteome_index = build_proteome_index(proteome_path)
        refined_boundaries = _seeds_to_refined_boundaries(merged_seeds)
        _recalculate_boundary_composition(
            refined_boundaries,
            masked_path=masked_path,
        )
        boundaries_bed_path = _write_phase2_checkpoints(
            output_dir=output_dir,
            refined_boundaries=refined_boundaries,
            boundary_taxonomy_map={},
            boundary_control_stats=None,
            boundary_diamond_query=None,
        )
        return Phase2Result(
            refined_boundaries=refined_boundaries,
            boundary_taxonomy_map={},
            boundary_control_stats=None,
            boundary_diamond_query=None,
            proteome_index=proteome_index,
            boundaries_bed=boundaries_bed_path,
            elapsed=time.time() - phase2_start,
            phase_outcome="passthrough",
            ablation_counts=InterventionCounts(
                opportunities=len(merged_seeds),
                interventions=len(refined_boundaries),
                changed=0,
            ),
        )

    superset_diamond_hits = None
    superset_proteome_index = None
    if boundary_diamond_superset_prototype_enabled and not (resume and resume_authorized) and merged_seeds:
        superset_db = _require_phase2b_gene_taxonomy_db(
            gene_taxonomy_faa_db,
            has_seeds=True,
        )
        superset_proteome_index = build_proteome_index(proteome_path)
        superset_top_k = max(10, boundary_diamond_top_k)
        logger.info(
            "Phase 2 superset prototype: searching the full proteome once (top_k=%d)",
            superset_top_k,
        )
        superset_diamond_hits = run_full_proteome_diamond(
            proteome_fasta=proteome_path,
            diamond_db=superset_db,
            output_dir=output_dir / "phase2" / "superset_diamond",
            max_target_seqs=superset_top_k,
            threads=gene_taxonomy_threads or threads,
            search_backend=search_backend,
        )
    if resume and resume_authorized:
        missing_state = [path for path in (phase2_state_path, phase2_resume_state_path) if not path.is_file()]
        if missing_state:
            raise ValueError(
                "authenticated Phase 2 is missing lossless resume state: "
                + ", ".join(str(path) for path in missing_state)
            )
        logger.info("Phase 2a: host trim skipped (verified schema-v3 resume)")
    elif boundary_host_trim_enabled and host_signature_model and host_signature_model.token_weights:
        gene_taxonomy_db = Path(gene_taxonomy_faa_db) if gene_taxonomy_faa_db else None
        if gene_taxonomy_db and merged_seeds:
            from virosync.orchestration.tasks import gene_taxonomy_batch_task
            from virosync.pipeline.phase2.host_signature_trim import (
                HostTrimParams,
                trim_seeds_by_host_signature,
            )

            host_trim_dir = output_dir / "phase2" / "host_trim"
            host_trim_dir.mkdir(parents=True, exist_ok=True)
            regions_payload = [
                {
                    "eve_id": f"EVE_{s.scaffold}_{s.start}-{s.end}",
                    "scaffold": s.scaffold,
                    "start": s.start,
                    "end": s.end,
                }
                for s in merged_seeds
            ]
            host_trim_threads = gene_taxonomy_threads or threads
            logger.info(
                "Phase 2a: Running gene taxonomy for host-trim (%d regions)",
                len(regions_payload),
            )
            # Pass Phase 1 marker_validation dir for taxonomy distribution analysis
            marker_validation_dir = output_dir / "phase1" / "marker_validation"
            # Pass HMM hits file for building genome-specific host baseline from COG/BUSCO markers
            hmm_hits_file = output_dir / "phase1" / "hmm" / "hmm_hits.tsv"
            if superset_diamond_hits is not None:
                from virosync.pipeline.phase3.gene_taxonomy import (
                    materialize_gene_taxonomy_batch_from_cached_hits,
                )

                gene_tax_map = materialize_gene_taxonomy_batch_from_cached_hits(
                    regions=regions_payload,
                    proteome_fasta=proteome_path,
                    diamond_hits=superset_diamond_hits,
                    output_dir=host_trim_dir / "gene_taxonomy",
                    high_pident_euk_threshold=high_pident_host_threshold,
                )
            else:
                gene_tax_map = call_task(
                    gene_taxonomy_batch_task,
                    regions=regions_payload,
                    proteome_path=proteome_path,
                    gene_taxonomy_db=gene_taxonomy_db,
                    output_dir=host_trim_dir / "gene_taxonomy",
                    threads=host_trim_threads,
                    high_pident_host_threshold=high_pident_host_threshold,
                    host_label=host_label,
                    marker_validation_dir=marker_validation_dir if marker_validation_dir.exists() else None,
                    hmm_hits_file=hmm_hits_file if hmm_hits_file.exists() else None,
                    search_backend=search_backend,
                )
            trim_params = HostTrimParams(
                window_bp=boundary_host_trim_window_bp,
                step_bp=boundary_host_trim_step_bp,
                max_host_fraction=boundary_host_trim_max_host_fraction,
                min_viral_fraction=boundary_host_trim_min_viral_fraction,
                host_score_threshold=boundary_host_trim_score_threshold,
                buffer_kb=boundary_host_trim_buffer_kb,
            )
            merged_seeds, trim_summaries = trim_seeds_by_host_signature(
                seeds=merged_seeds,
                gene_taxonomy_map=gene_tax_map or {},
                host_model=host_signature_model,
                validated_markers=validated_markers,
                params=trim_params,
                ablation_id=ablation_id,
            )
            host_coordinate_counts = _sum_intervention_counts(
                host_coordinate_counts,
                InterventionCounts(
                    opportunities=sum(int(row["host_coordinate_change_opportunities"]) for row in trim_summaries),
                    interventions=sum(int(row["host_coordinate_change_interventions"]) for row in trim_summaries),
                    changed=sum(int(row["host_coordinate_change_changed"]) for row in trim_summaries),
                ),
            )
            if trim_summaries and len(trim_summaries) == len(merged_seeds):
                for seed, summary in zip(merged_seeds, trim_summaries):
                    seed.host_trim_original_start = summary.get("start")
                    seed.host_trim_original_end = summary.get("end")
                    seed.host_trimmed_start = summary.get("trimmed_start")
                    seed.host_trimmed_end = summary.get("trimmed_end")
                    seed.host_trim_reason = summary.get("reason", "")
                    seed.host_trim_common_euk_taxonomy = ""
            summary_path = host_trim_dir / "host_signature_trim.tsv"
            with atomic_write_context(summary_path, "w") as handle:
                handle.write(
                    "scaffold\tstart\tend\ttrimmed_start\ttrimmed_end\treason\t"
                    "window_count\tgood_windows\tgood_marker_windows\thost_consensus_taxonomy\n"
                )
                for row in trim_summaries:
                    handle.write(
                        f"{row['scaffold']}\t{row['start']}\t{row['end']}\t"
                        f"{row['trimmed_start']}\t{row['trimmed_end']}\t{row['reason']}\t"
                        f"{row.get('window_count', 0)}\t{row.get('good_windows', 0)}\t"
                        f"{row.get('good_marker_windows', 0)}\t"
                        f"{row.get('host_consensus_taxonomy', '.')}\n"
                    )
            logger.info("Phase 2a: wrote host trim summary %s", summary_path)
        else:
            logger.info("Phase 2a: host trim skipped (no gene taxonomy DB or seeds)")

    if not merged_seeds:
        logger.warning("No seeds found - pipeline complete with 0 predictions")
        output_files_empty = call_task(
            generate_outputs_task,
            verification_results=[],
            output_dir=output_dir,
            genome_path=masked_path,
            proteome_path=proteome_path,
            accepted_only=True,
            extended_output=extended_output,
        )
        total_elapsed = time.time() - genome_start_time
        logger.info(f"Genome {genome_id} complete: {total_elapsed:.1f}s (0 predictions)")
        output_files = {
            "accepted": output_files_empty,
            **_generate_required_reports(
                output_dir=output_dir,
                genome_id=genome_id,
                taxonomy_labels_file=taxonomy_labels_file,
                logger=logger,
            ),
        }
        _write_empty_run_log(
            output_dir=output_dir,
            genome_id=genome_id,
            input_path=masked_path,
            reason="no seeds",
            elapsed_sec=total_elapsed,
            output_files=output_files,
            fingerprint=config_fingerprint,
        )
        return Phase2Terminal(
            {
                "genome_id": genome_id,
                "success": True,
                **_empty_prediction_summary(),
                "output_files": output_files,
                "elapsed_sec": total_elapsed,
                "ablation_counts": host_coordinate_counts,
            }
        )

    if resume and resume_authorized:
        resume_state = load_phase2_resume_state(phase2_resume_state_path)
        boundary_report_state = load_phase2_state(phase2_state_path)
        if phase2_state_to_document(resume_state.refined_boundaries) != phase2_state_to_document(boundary_report_state):
            raise ValueError("authenticated Phase 2 boundary and resume checkpoints disagree")

        refined_boundaries = resume_state.refined_boundaries
        boundary_taxonomy_map = resume_state.boundary_taxonomy_map
        boundary_control_stats = resume_state.boundary_control_stats
        boundary_diamond_query = resume_state.boundary_diamond_query
        proteome_index = build_proteome_index(proteome_path)
        logger.info(
            "Phase 2 resume: loaded %d boundaries, %d taxonomy records, and the exact Diamond query from %s",
            len(refined_boundaries),
            len(boundary_taxonomy_map),
            phase2_resume_state_path,
        )

        phase2_elapsed = time.time() - phase2_start
        return Phase2Result(
            refined_boundaries=refined_boundaries,
            boundary_taxonomy_map=boundary_taxonomy_map,
            boundary_control_stats=boundary_control_stats,
            boundary_diamond_query=boundary_diamond_query,
            proteome_index=proteome_index,
            goto_phase3=True,
            boundaries_bed=refined_bed,
            elapsed=phase2_elapsed,
        )

    # === Gene-based seed extension ===
    # Extend each seed by ±5 genes and merge overlapping seeds BEFORE Diamond.
    # Seeds are now gene-anchored before taxonomy refinement.
    if merged_seeds:
        from virosync.pipeline.phase2.boundary_refiner import extend_seeds_by_genes

        proteome_index = (
            superset_proteome_index if superset_proteome_index is not None else build_proteome_index(proteome_path)
        )
        pre_extend_count = len(merged_seeds)
        merged_seeds = extend_seeds_by_genes(
            merged_seeds,
            proteome_index,
            extension_genes=5,
        )
        logger.info(
            "Phase 2: Extended seeds by ±5 genes and merged overlapping: %d -> %d seeds",
            pre_extend_count,
            len(merged_seeds),
        )

    # === PHASE 2b: Batched Diamond for boundary refinement ===
    # Run Diamond ONCE per genome on ALL seeds + flanking genes + control genes
    # Phase 2b is MANDATORY for flanking gene integration - provides taxonomy for all genes
    # within the max boundary extension range (+/-flank_genes from each seed)
    boundary_taxonomy_map: dict[str, GeneTaxonomy] = {}
    boundary_control_stats: ControlStats | None = None
    boundary_diamond_query: GenomeDiamondQuery | None = None

    phase2b_db = _require_phase2b_gene_taxonomy_db(
        gene_taxonomy_faa_db,
        has_seeds=bool(merged_seeds),
    )

    if phase2b_db and merged_seeds:
        logger.info("-" * 60)
        logger.info("Phase 2b: Running batched Diamond for boundary refinement")
        logger.info(
            "Phase 2b: Proteome index: %d scaffolds, %d total pORFs",
            len(proteome_index),
            sum(len(porfs) for porfs in proteome_index.values()),
        )
        boundary_evidence = _run_boundary_diamond(
            merged_seeds=merged_seeds,
            proteome_index=proteome_index,
            proteome_path=proteome_path,
            gene_taxonomy_db=phase2b_db,
            output_dir=output_dir,
            config=config,
            superset_diamond_hits=superset_diamond_hits,
            logger=logger,
        )
        boundary_taxonomy_map = boundary_evidence.taxonomy_map
        boundary_control_stats = boundary_evidence.control_stats
        boundary_diamond_query = boundary_evidence.query

    # === PHASE 2c: Taxonomy-based seed refinement ===
    if merged_seeds and boundary_taxonomy_map and boundary_diamond_query:
        from virosync.pipeline.phase2 import taxonomy_seed_refiner
        from virosync.pipeline.phase2.taxonomy_ml_refiner import (
            refine_seeds_by_taxonomy_ml,
        )

        logger.info("-" * 60)
        logger.info("Phase 2c: Taxonomy-based seed refinement")
        logger.info("-" * 60)

        def _count_boundary_updates(original, refined) -> int:
            if not isinstance(refined, list) or len(refined) != len(original):
                return 0
            updated = 0
            for old_seed, new_seed in zip(original, refined):
                if old_seed.start != new_seed.start or old_seed.end != new_seed.end:
                    updated += 1
            return updated

        heuristic_kwargs = {
            "merged_seeds": merged_seeds,
            "taxonomy_map": boundary_taxonomy_map,
            "seed_gene_mappings": boundary_diamond_query.seed_gene_mappings,
            "host_prefix": f"{host_label}__",
            "expansion_kb": 5,
            "host_signature_model": host_signature_model,
            "host_signature_threshold": boundary_host_trim_score_threshold,
        }
        taxonomy_refined_seeds = None

        taxonomy_seed_refiner.validate_taxonomy_refinement_mode(
            ablation_id=ablation_id,
            taxonomy_ml_enabled=boundary_taxonomy_ml_enabled,
        )
        if boundary_taxonomy_ml_enabled:
            try:
                ml_result = refine_seeds_by_taxonomy_ml(
                    merged_seeds=merged_seeds,
                    taxonomy_map=boundary_taxonomy_map,
                    boundary_query=boundary_diamond_query,
                    host_signature_model=host_signature_model,
                    model_type=boundary_taxonomy_ml_model,
                    host_prefix=f"{host_label}__",
                    probability_threshold=boundary_taxonomy_ml_threshold,
                    neighbor_window=boundary_taxonomy_ml_neighbor_window,
                    output_dir=output_dir / "phase2" / "taxonomy_ml",
                    random_state=boundary_diamond_random_seed,
                    taxonomy_weight_mode=taxonomy_weight_mode,
                    host_signature_threshold=boundary_host_trim_score_threshold,
                )
            except Exception as exc:
                logger.warning(
                    "Phase 2c: Taxonomy ML refinement failed (%s); using heuristic fallback",
                    exc,
                )
            else:
                ml_updates = _count_boundary_updates(merged_seeds, ml_result)
                if ml_updates > 0:
                    taxonomy_refined_seeds = ml_result
                    logger.info(
                        "Phase 2c: Taxonomy ML refinement updated %d/%d seeds "
                        "(model=%s, threshold=%.3f, neighbor_window=%d)",
                        ml_updates,
                        len(merged_seeds),
                        boundary_taxonomy_ml_model,
                        boundary_taxonomy_ml_threshold,
                        boundary_taxonomy_ml_neighbor_window,
                    )
                else:
                    logger.info(
                        "Phase 2c: Taxonomy ML refinement produced no boundary updates; using heuristic fallback"
                    )

        if taxonomy_refined_seeds is None:
            taxonomy_evaluation = taxonomy_seed_refiner.evaluate_taxonomy_seed_refinement(
                **heuristic_kwargs,
                ablation_id=ablation_id,
            )
            taxonomy_refined_seeds = list(taxonomy_evaluation.selected_seeds)
            host_coordinate_counts = _sum_intervention_counts(
                host_coordinate_counts,
                taxonomy_evaluation.intervention_counts,
            )
            heuristic_updates = _count_boundary_updates(merged_seeds, taxonomy_refined_seeds)
            logger.info(
                "Phase 2c: Heuristic taxonomy refinement updated %d/%d seeds",
                heuristic_updates,
                len(merged_seeds),
            )

        # Use taxonomy-refined seeds for downstream steps
        merged_seeds = taxonomy_refined_seeds
    else:
        logger.warning("Skipping taxonomy-based seed refinement (missing Diamond data)")

    # Initialize Phase 2 output variables
    boundaries_bed_path = refined_bed
    refined_boundaries = []

    # === Convert seeds directly to RefinedBoundary ===
    # Boundaries are gene-anchored from extension above.
    # Boundary confidence fields remain available for output compatibility.

    if merged_seeds:
        refined_boundaries = _seeds_to_refined_boundaries(merged_seeds)
        logger.info(
            "Phase 2: Converted %d gene-extended seeds to RefinedBoundary",
            len(refined_boundaries),
        )

    # Phase 2f: Host Taxonomy Trimming
    if refined_boundaries and boundary_taxonomy_map and boundary_control_stats and boundary_diamond_query:
        from virosync.pipeline.phase2.boundary_diamond import build_taxonomy_consensus
        from virosync.pipeline.phase2.boundary_refiner import trim_boundary_by_host_taxonomy

        logger.info("-" * 60)
        logger.info("Phase 2c: Host Taxonomy Trimming")
        logger.info("-" * 60)

        # Build taxonomic consensus from control genes using full taxonomy lookup
        control_taxonomy = [
            boundary_taxonomy_map[pid]
            for pid in boundary_diamond_query.control_porf_ids
            if pid in boundary_taxonomy_map
        ]

        # Build host baseline fingerprint from control genes
        from virosync.pipeline.host_signatures import TaxonomyLabelLookup
        from virosync.pipeline.phase2.boundary_diamond import build_host_baseline_fingerprint

        host_baseline_fingerprint = {}
        taxonomy_consensus = "unknown"
        if taxonomy_labels_file and Path(taxonomy_labels_file).exists():
            tax_lookup = TaxonomyLabelLookup.load(Path(taxonomy_labels_file))
            taxonomy_consensus = build_taxonomy_consensus(
                control_taxonomy,
                taxonomy_lookup=tax_lookup,
                host_prefix=f"{host_label}__",
                min_token_length=boundary_host_signature_min_token_len,
                weight_mode=taxonomy_weight_mode,
            )

            # Build host baseline fingerprint
            host_baseline_fingerprint = build_host_baseline_fingerprint(
                control_taxonomy,
                host_prefix=f"{host_label}__",
                min_token_count=3,
                min_weight_fraction=0.10,
            )

            if host_baseline_fingerprint:
                logger.info(
                    "Phase 2c: Host baseline fingerprint built (%d tokens)",
                    len(host_baseline_fingerprint),
                )
            else:
                logger.warning("Phase 2c: Host baseline fingerprint empty - using fallback prefix matching")

        # Log host taxonomy consensus with full lineage
        logger.info(f"Host Taxonomy Consensus (from {boundary_control_stats.n_genes} control genes):")
        logger.info(f"  Lineage: {taxonomy_consensus}")
        logger.info(f"  Host frequency: {boundary_control_stats.host_frequency:.2%}")
        logger.info(f"  No-hit frequency: {boundary_control_stats.no_hit_frequency:.2%}")
        logger.info(f"  Mean pident (host hits): {boundary_control_stats.mean_pident:.1f}%")
        logger.info(f"  Dominant organism: {boundary_control_stats.dominant_organism}")

        # Apply trimming to each boundary
        trimmed_boundaries = []
        n_trimmed = 0
        total_upstream_bp = 0
        total_downstream_bp = 0

        # Get min_overlap_score from config
        min_overlap = boundary_host_trim_min_overlap_score

        _n_no_seed_id = 0
        _n_no_mapping = 0
        for boundary in refined_boundaries:
            if boundary.predicted_family == "CRESS":
                trimmed_boundaries.append(boundary)
                continue
            seed_mapping = None
            if boundary.seed_id:
                seed_mapping = boundary_diamond_query.seed_gene_mappings.get(boundary.seed_id)
            if not boundary.seed_id:
                _n_no_seed_id += 1
            elif not seed_mapping:
                _n_no_mapping += 1

            if seed_mapping:
                new_start, new_end, stats = trim_boundary_by_host_taxonomy(
                    boundary=boundary,
                    seed_mapping=seed_mapping,
                    taxonomy_map=boundary_taxonomy_map,
                    control_stats=boundary_control_stats,
                    host_prefix=f"{host_label}__",
                    host_baseline_fingerprint=host_baseline_fingerprint,
                    host_signature_model=host_signature_model,
                    host_signature_threshold=boundary_host_trim_score_threshold,
                    min_overlap_score=min_overlap,
                    taxonomy_weight_mode=taxonomy_weight_mode,
                    unknown_neighbor_window=1,
                    unknown_host_penalty=2.0,
                    unknown_viral_bonus=2.0,
                    ablation_id=ablation_id,
                )
                host_coordinate_counts = _sum_intervention_counts(
                    host_coordinate_counts,
                    InterventionCounts(
                        opportunities=int(stats["host_coordinate_change_opportunities"]),
                        interventions=int(stats["host_coordinate_change_interventions"]),
                        changed=int(stats["host_coordinate_change_changed"]),
                    ),
                )

                if stats.get("trimmed"):
                    n_trimmed += 1
                    total_upstream_bp += stats.get("upstream_bp", 0)
                    total_downstream_bp += stats.get("downstream_bp", 0)

                    logger.info(
                        f"Host trim: {boundary.seed_id} "
                        f"{boundary.start}-{boundary.end} -> {new_start}-{new_end} "
                        f"(upstream={stats['upstream_bp']}bp, downstream={stats['downstream_bp']}bp)"
                    )
                    if stats.get("upstream_stopped_by"):
                        logger.info(f"  Upstream stopped by: {stats['upstream_stopped_by']}")
                    if stats.get("downstream_stopped_by"):
                        logger.info(f"  Downstream stopped by: {stats['downstream_stopped_by']}")

                    # Update boundary
                    boundary.start = new_start
                    boundary.end = new_end

            trimmed_boundaries.append(boundary)

        refined_boundaries = trimmed_boundaries

        logger.info(
            f"Host taxonomy trimming: {n_trimmed}/{len(trimmed_boundaries)} boundaries trimmed "
            f"(upstream={total_upstream_bp:,}bp, downstream={total_downstream_bp:,}bp total) "
            f"[no_seed_id={_n_no_seed_id}, no_mapping={_n_no_mapping}]"
        )
        logger.info("-" * 60)

    # Apply boundary constraints to ensure boundaries don't exceed +/-N genes from original seed
    # This is critical: the flanking genes in Phase 2b only cover +/-flank_genes from the seed
    if refined_boundaries and boundary_diamond_query and boundary_diamond_query.seed_gene_mappings:
        from virosync.pipeline.phase2.boundary_refiner import constrain_to_seed_bounds

        constrained_boundaries = []
        n_constrained = 0
        for boundary in refined_boundaries:
            seed_mapping = boundary_diamond_query.seed_gene_mappings.get(boundary.seed_id)
            if seed_mapping:
                original_start, original_end = boundary.start, boundary.end
                constrained = constrain_to_seed_bounds(boundary, seed_mapping)
                constrained_boundaries.append(constrained)
                if constrained.start != original_start or constrained.end != original_end:
                    n_constrained += 1
            else:
                # No seed mapping available (e.g., resume from BED) - keep original
                constrained_boundaries.append(boundary)
                if boundary.seed_id:
                    logger.debug(
                        "No seed mapping for %s - boundary not constrained",
                        boundary.seed_id,
                    )
        refined_boundaries = constrained_boundaries
        if n_constrained > 0:
            logger.info(
                "Boundary constraint applied to %d/%d regions (max +/-%d genes from seed)",
                n_constrained,
                len(refined_boundaries),
                boundary_diamond_flank_genes,
            )

    # Check for adjacent EVE regions that should be merged based on taxonomy
    # If the gap between two EVEs contains viral genes, merge them
    if refined_boundaries and boundary_taxonomy_map and len(refined_boundaries) > 1:
        from virosync.pipeline.phase2.boundary_refiner import merge_adjacent_viral_boundaries

        pre_merge_count = len(refined_boundaries)
        refined_boundaries = merge_adjacent_viral_boundaries(
            boundaries=refined_boundaries,
            taxonomy_map=boundary_taxonomy_map,
            proteome_index=proteome_index,
            max_gap_bp=10000,  # Max 10kb gap to consider merging
            min_viral_fraction=0.3,  # At least 30% of gap genes must be viral
        )
        if len(refined_boundaries) < pre_merge_count:
            logger.info(
                "Post-taxonomy merge: %d -> %d EVE regions",
                pre_merge_count,
                len(refined_boundaries),
            )

    # Terminal repeats are read from the raw genome because Phase 0 hard-masks
    # repeats. Run this after every heuristic coordinate change and merge, but
    # before marker-floor metadata and coordinate-derived composition evidence.
    if refined_boundaries and raw_genome_path is not None and boundary_diamond_query is not None:
        from virosync.pipeline.phase2.terminal_repeats import (
            refine_boundaries_with_terminal_repeats,
        )

        refined_boundaries = refine_boundaries_with_terminal_repeats(
            refined_boundaries,
            raw_genome_path=raw_genome_path,
            validated_markers=validated_markers,
            boundary_diamond_query=boundary_diamond_query,
            boundary_taxonomy_map=boundary_taxonomy_map,
            proteome_index=proteome_index,
            extension_bp=config.phase1.extension_kb * 1000,
        )
        n_tir_overrides = sum(boundary.tir_boundary_override for boundary in refined_boundaries)
        n_tir_present = sum(boundary.tir_present for boundary in refined_boundaries)
        logger.info(
            "Terminal inverted repeats: %d candidate-bearing boundaries, %d authoritative overrides",
            n_tir_present,
            n_tir_overrides,
        )

    # === Validated-marker floor metadata (Phase-3 re-admit input) ===
    # Phase-2 host trimming can collapse a marker-dense NCLDV/MIRUS boundary below
    # its validated-marker span, stripping hallmark genes out so the v2 gate
    # (length > 5000 OR has_mcp) rejects a genuine EVE. We do NOT mutate the
    # boundary here (extending an accepted region pulls host genes in and regresses
    # its confidence/tier -> NCLDV loss;. Instead we
    # only record the floored-alternative span as metadata on each boundary; the
    # Phase-3 re-admit pass consumes it on REJECTED boundaries only. Scoped to each
    # boundary's own seed span, validated markers only -> host regions get no floor.
    if refined_boundaries and validated_markers:
        from virosync.pipeline.phase2.boundary_refiner import (
            annotate_boundaries_with_marker_floor,
        )

        n_marker_floored = annotate_boundaries_with_marker_floor(refined_boundaries, validated_markers)
        if n_marker_floored:
            logger.info(
                "Validated-marker floor recorded on %d/%d boundaries (Phase-3 re-admit candidates)",
                n_marker_floored,
                len(refined_boundaries),
            )

    _recalculate_boundary_composition(
        refined_boundaries,
        masked_path=masked_path,
    )

    # Log region statistics with before/after comparison
    if refined_boundaries and merged_seeds:
        logger.info("-" * 60)
        logger.info("Phase 2: Boundary Refinement Summary")
        logger.info("-" * 60)

        # Before (seeds)
        length_before = sum(s.length for s in merged_seeds)
        avg_before = length_before // len(merged_seeds)

        # After (refined boundaries)
        length_after = sum(b.end - b.start for b in refined_boundaries)
        avg_after = length_after // len(refined_boundaries)
        min_bp = min(b.end - b.start for b in refined_boundaries)
        max_bp = max(b.end - b.start for b in refined_boundaries)

        # Shrinkage
        shrinkage_pct = (1 - length_after / length_before) * 100 if length_before > 0 else 0

        logger.info(f"  Boundaries refined: {len(refined_boundaries)}")
        logger.info(f"  Total length: {length_before:,} bp -> {length_after:,} bp ({shrinkage_pct:+.1f}%)")
        logger.info(f"  Average length: {avg_before:,} bp -> {avg_after:,} bp")
        logger.info(f"  Range: {min_bp:,}-{max_bp:,} bp")
        logger.info("-" * 60)
    elif refined_boundaries:
        # No merged_seeds for comparison (e.g., resume)
        total_bp = sum(b.end - b.start for b in refined_boundaries)
        avg_bp = total_bp // len(refined_boundaries) if refined_boundaries else 0
        logger.info(f"Refined boundaries: {len(refined_boundaries)} regions, total={total_bp:,}bp, avg={avg_bp:,}bp")
    else:
        logger.info("Refined boundaries: 0")

    boundaries_bed_path = _write_phase2_checkpoints(
        output_dir=output_dir,
        refined_boundaries=refined_boundaries,
        boundary_taxonomy_map=boundary_taxonomy_map,
        boundary_control_stats=boundary_control_stats,
        boundary_diamond_query=boundary_diamond_query,
    )

    phase2_elapsed = time.time() - phase2_start
    logger.info(
        "Phase 2 complete: %.1fs, wrote %s, %s, and %s",
        phase2_elapsed,
        boundaries_bed_path,
        phase2_state_path,
        phase2_resume_state_path,
    )

    if not refined_boundaries:
        logger.warning("No boundaries refined - pipeline complete with 0 predictions")
        output_files_empty = call_task(
            generate_outputs_task,
            verification_results=[],
            output_dir=output_dir,
            genome_path=masked_path,
            proteome_path=proteome_path,
            accepted_only=True,
            extended_output=extended_output,
        )
        total_elapsed = time.time() - genome_start_time
        logger.info(f"Genome {genome_id} complete: {total_elapsed:.1f}s (0 predictions)")
        output_files = {
            "accepted": output_files_empty,
            "phase2_boundaries": str(boundaries_bed_path),
            **_generate_required_reports(
                output_dir=output_dir,
                genome_id=genome_id,
                taxonomy_labels_file=taxonomy_labels_file,
                logger=logger,
            ),
        }
        _write_empty_run_log(
            output_dir=output_dir,
            genome_id=genome_id,
            input_path=masked_path,
            reason="no refined boundaries",
            elapsed_sec=total_elapsed,
            output_files=output_files,
            fingerprint=config_fingerprint,
        )
        return Phase2Terminal(
            {
                "genome_id": genome_id,
                "success": True,
                **_empty_prediction_summary(),
                "output_files": output_files,
                "elapsed_sec": total_elapsed,
                "ablation_counts": host_coordinate_counts,
            }
        )

    return Phase2Result(
        refined_boundaries=refined_boundaries,
        boundary_taxonomy_map=boundary_taxonomy_map,
        boundary_control_stats=boundary_control_stats,
        boundary_diamond_query=boundary_diamond_query,
        proteome_index=proteome_index,
        boundaries_bed=boundaries_bed_path,
        elapsed=phase2_elapsed,
        ablation_counts=host_coordinate_counts,
    )
