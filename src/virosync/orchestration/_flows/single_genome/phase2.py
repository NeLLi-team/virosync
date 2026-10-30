"""Phase 2 subflow: bounded, seed-preserving growth from gene evidence."""

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from virosync.ablation import AblationID, InterventionCounts
from virosync.config import PipelineConfig
from virosync.orchestration.runtime import call_task
from virosync.orchestration.tasks import generate_outputs_task
from virosync.pipeline.phase1.viral_markers import get_assembly_mode
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
from virosync.pipeline.phase2.boundary_refiner import (
    RefinedBoundary,
    assign_boundary_candidate_ids,
    boundary_candidate_id,
)
from virosync.pipeline.phase2.repeat_evidence import assess_repeat_evidence
from virosync.pipeline.phase2.seed_extension import (
    confirmed_seed_core,
    extend_confirmed_seed,
    merge_reachable_boundaries,
    refresh_boundary_classification,
    validate_source_membership,
)
from virosync.pipeline.phase2.terminal_repeats import TerminalRepeat
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
            confirmed_seeds=seed.confirmed_seeds,
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
        )
        for seed in merged_seeds
    ]


def _annotate_boundary_repeats(
    boundaries: list[RefinedBoundary],
    *,
    raw_genome_path: Path | None,
    boundary_diamond_query: GenomeDiamondQuery | None,
    extension_bp: int,
) -> None:
    """Attach raw-sequence repeat evidence after assigning public candidate IDs."""
    from Bio import SeqIO

    if not boundaries:
        return
    by_scaffold: dict[str, list[RefinedBoundary]] = defaultdict(list)
    coverage_by_scaffold: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for boundary in boundaries:
        by_scaffold[boundary.scaffold].append(boundary)
    if boundary_diamond_query is not None:
        for mapping in boundary_diamond_query.seed_gene_mappings.values():
            coverage_by_scaffold[mapping.scaffold].append((mapping.flank_start_bp, mapping.flank_end_bp))

    scaffold_index = SeqIO.index(str(raw_genome_path), "fasta") if raw_genome_path is not None else None
    try:
        for scaffold, scaffold_boundaries in by_scaffold.items():
            sequence = (
                str(scaffold_index[scaffold].seq) if scaffold_index is not None and scaffold in scaffold_index else None
            )
            for boundary in scaffold_boundaries:
                boundary.repeat_evidence = assess_repeat_evidence(
                    sequence,
                    scaffold=scaffold,
                    input_id=boundary.seed_id or boundary_candidate_id(boundary),
                    start=boundary.start,
                    end=boundary.end,
                    parent_start=boundary.pre_tir_start if boundary.pre_tir_start is not None else boundary.start,
                    parent_end=boundary.pre_tir_end if boundary.pre_tir_end is not None else boundary.end,
                    coverage_intervals=tuple(coverage_by_scaffold[scaffold]),
                    extension_bp=extension_bp,
                    tir=_retained_terminal_repeat(boundary),
                )
    finally:
        if scaffold_index is not None:
            scaffold_index.close()


def _retained_terminal_repeat(boundary: RefinedBoundary) -> TerminalRepeat | None:
    """Recover the retained legacy TIR pair without changing its interpretation."""
    if not boundary.tir_present:
        return None
    assert boundary.tir_left_start is not None
    assert boundary.tir_left_end is not None
    assert boundary.tir_right_start is not None
    assert boundary.tir_right_end is not None
    return TerminalRepeat(
        left_start=boundary.tir_left_start,
        left_end=boundary.tir_left_end,
        right_start=boundary.tir_right_start,
        right_end=boundary.tir_right_end,
        identity=boundary.tir_identity,
        alignment_length=boundary.tir_alignment_length,
        alignment_capped=boundary.tir_alignment_capped,
    )


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
            score = int(boundary.confidence * 1000)
            handle.write(
                f"{boundary.scaffold}\t{boundary.start}\t{boundary.end}\t"
                f"{boundary_candidate_id(boundary)}\t{score}\t.\n"
            )
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
    """Search bounded gene context and extend complete confirmed seeds outward.

    The shared controller applies connected support and permanent host stops.
    Repeat partitions must retain whole source seeds within the reachable span.

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
    boundary_host_trim_enabled = phase2.host_trim_enabled
    boundary_host_trim_score_threshold = phase2.host_trim_score_threshold
    taxonomy_weight_mode = phase2.taxonomy_weight_mode
    boundary_taxonomy_ml_enabled = phase2.taxonomy_ml_enabled
    boundary_taxonomy_ml_model = phase2.taxonomy_ml_model
    boundary_taxonomy_ml_threshold = phase2.taxonomy_ml_threshold
    boundary_taxonomy_ml_neighbor_window = phase2.taxonomy_ml_neighbor_window
    boundary_diamond_top_k = phase2.diamond_top_k
    boundary_diamond_random_seed = phase2.diamond_random_seed
    boundary_diamond_superset_prototype_enabled = phase2.diamond_superset_prototype_enabled
    gene_taxonomy_threads = compute.gene_taxonomy_threads
    extended_output = phase3.extended_output
    search_backend = compute.search_backend.value
    refined_bed = output_dir / "phase2" / "refined_boundaries.bed"

    phase2_start = time.time()
    host_coordinate_counts = InterventionCounts()
    merged_seeds = [confirmed_seed_core(seed) for seed in merged_seeds]

    # Authenticated state is the only resume source.
    phase2_state_path = output_dir / "phase2" / PHASE2_STATE_FILENAME
    phase2_resume_state_path = output_dir / "phase2" / PHASE2_RESUME_STATE_FILENAME
    if ablation_id is AblationID.A3 and merged_seeds and not (resume and resume_authorized):
        logger.info(
            "A3: forwarding %d exact Phase-1 seed intervals to Phase 3",
            len(merged_seeds),
        )
        proteome_index = build_proteome_index(proteome_path)
        refined_boundaries = _seeds_to_refined_boundaries(merged_seeds)
        validate_source_membership(merged_seeds, refined_boundaries)
        refresh_boundary_classification(
            refined_boundaries,
            validated_markers,
            ablation_id=ablation_id,
            single_marker_min_score=get_assembly_mode(config.phase1.assembly_mode.value).single_marker_min_score,
        )
        _recalculate_boundary_composition(
            refined_boundaries,
            masked_path=masked_path,
        )
        refined_boundaries = assign_boundary_candidate_ids(refined_boundaries)
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
        logger.info("Phase 2: using authenticated seed-preserving resume state")

    if not merged_seeds and not (resume and resume_authorized):
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
        validate_source_membership(merged_seeds, refined_boundaries)
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

    proteome_index = (
        superset_proteome_index if superset_proteome_index is not None else build_proteome_index(proteome_path)
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

    ml_support: dict[str, frozenset[str]] = {}
    if boundary_taxonomy_ml_enabled and boundary_diamond_query is not None:
        from virosync.pipeline.phase2.taxonomy_ml_refiner import predict_taxonomy_support

        try:
            ml_support = predict_taxonomy_support(
                taxonomy_map=boundary_taxonomy_map,
                boundary_query=boundary_diamond_query,
                model_type=boundary_taxonomy_ml_model,
                host_prefix=f"{host_label}__",
                probability_threshold=boundary_taxonomy_ml_threshold,
                neighbor_window=boundary_taxonomy_ml_neighbor_window,
                output_dir=output_dir / "phase2" / "taxonomy_ml",
                random_state=boundary_diamond_random_seed,
                taxonomy_weight_mode=taxonomy_weight_mode,
            )
        except Exception as exc:
            logger.warning("Taxonomy ML failed (%s); using deterministic extension", exc)

    refined_boundaries = []
    extension_bp = config.phase1.extension_kb * 1000
    for seed in merged_seeds:
        mapping = boundary_diamond_query.seed_gene_mappings.get(seed.seed_id) if boundary_diamond_query else None
        envelope = (seed.start, seed.end)
        if mapping is not None and seed.predicted_family != "CRESS":
            envelope = (
                max(0, seed.start - extension_bp, min(seed.start, mapping.flank_start_bp)),
                min(seed.end + extension_bp, max(seed.end, mapping.flank_end_bp)),
            )
        kwargs = dict(
            seed=seed,
            genes=proteome_index.get(seed.scaffold, []),
            taxonomy=boundary_taxonomy_map,
            envelope=envelope,
            host_model=host_signature_model if boundary_host_trim_enabled else None,
            host_threshold=boundary_host_trim_score_threshold,
            ml_support=ml_support.get(seed.seed_id, frozenset()),
        )
        boundary = extend_confirmed_seed(**kwargs, disable_host_stop=ablation_id is AblationID.A4)
        if ablation_id is AblationID.A4:
            normal = extend_confirmed_seed(**kwargs)
            affected = (normal.start, normal.end) != (boundary.start, boundary.end)
            host_coordinate_counts = _sum_intervention_counts(
                host_coordinate_counts,
                InterventionCounts(opportunities=int(affected), interventions=int(affected), changed=int(affected)),
            )
        refined_boundaries.append(boundary)
    refined_boundaries = merge_reachable_boundaries(refined_boundaries)

    # Terminal repeats are read from the raw genome because Phase 0 hard-masks
    # repeats. Run this after every heuristic coordinate change and merge, but
    # before coordinate-derived composition evidence.
    if (
        config.phase2.terminal_repeat_enabled
        and refined_boundaries
        and raw_genome_path is not None
        and boundary_diamond_query is not None
    ):
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

    validate_source_membership(merged_seeds, refined_boundaries)
    refresh_boundary_classification(
        refined_boundaries,
        validated_markers,
        ablation_id=ablation_id,
        single_marker_min_score=get_assembly_mode(config.phase1.assembly_mode.value).single_marker_min_score,
    )
    _recalculate_boundary_composition(
        refined_boundaries,
        masked_path=masked_path,
    )
    refined_boundaries = assign_boundary_candidate_ids(refined_boundaries)
    _annotate_boundary_repeats(
        refined_boundaries,
        raw_genome_path=raw_genome_path,
        boundary_diamond_query=boundary_diamond_query,
        extension_bp=config.phase1.extension_kb * 1000,
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
