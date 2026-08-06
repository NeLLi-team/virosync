"""Phase 1 subflow: HMM scan -> marker validation -> region assembly."""

import json
from pathlib import Path
from typing import Optional

from virosync.ablation import AblationID, InterventionCounts
from virosync.utils.atomic_write import atomic_write_context
from virosync.orchestration.runtime import call_task
from virosync.orchestration.tasks import (
    frameshift_screening_task,
    generate_outputs_task,
    hhg_seeding_task,
)
from virosync.pipeline.phase1.hhg_seeding import Anchor
from virosync.pipeline.phase1.marker_roles import decide_marker_hit_role
from virosync.pipeline.phase1.viral_markers import get_assembly_mode
from virosync.pipeline.phase3.mcp_detection import is_mcp_gene
from virosync.orchestration._flows.utils import (
    build_marker_faa,
    ensure_combined_faa,
    log_region_statistics,
)

from .manifest import (
    _empty_prediction_summary,
    _write_empty_run_log,
)
from .phase1_state import (
    PHASE1_STATE_FILENAME,
    load_phase1_state,
    write_phase1_state,
)
from .reports import _generate_required_reports


def _seed_annotation_markers(
    validated_markers: list,
    *additional_marker_groups: list,
) -> list:
    """Preserve the historical marker surface used to annotate existing seeds."""

    markers = list(validated_markers)
    for group in additional_marker_groups:
        for marker in group:
            if marker not in markers:
                markers.append(marker)
    return markers


def _region_coordinate_surface(regions: list) -> set[tuple[str, int, int]]:
    """Return the normalized seed-coordinate surface for a region collection."""

    return {(region.scaffold, region.start, region.end) for region in regions}


def _run_phase1_subflow(
    # Core inputs (from Phase 0)
    masked_path: Path,
    proteome_path: Path,
    repeat_regions: list,
    # Core identifiers
    output_dir: Path,
    genome_id: str,
    # HMM & Database parameters
    hmm_database: Optional[Path],
    hmm_allowlist: Optional[Path],
    hmm_chunk_size: Optional[int],
    frameshift_screening_enabled: bool,
    marker_faa_db: Optional[Path],
    marker_faa_dir: Optional[Path],
    marker_db: Optional[Path],
    faa_dir: Optional[Path],
    gene_taxonomy_faa_db: Optional[Path],
    # Taxonomy parameters
    taxonomy_labels_file: Optional[Path],
    host_prefixes: list[str],
    host_label: str,
    taxonomy_weight_mode: str,
    # Host taxonomy deviation parameters
    host_taxonomy_deviation_enabled: bool,
    host_taxonomy_deviation_allow_seeds: bool,
    host_taxonomy_deviation_min_token_len: int,
    host_taxonomy_deviation_min_tokens: int,
    host_taxonomy_deviation_overlap_threshold: float,
    host_taxonomy_deviation_max_pident: float,
    host_taxonomy_deviation_max_hits: int,
    host_taxonomy_deviation_window_bp: int,
    host_taxonomy_deviation_window_count: int,
    host_taxonomy_deviation_window_seed: int,
    host_taxonomy_deviation_window_min_markers: int,
    host_taxonomy_deviation_seed_window_bp: int,
    host_taxonomy_deviation_seed_min_markers: int,
    marker_validation_top_k: int,
    novel_marker_min_score: float,
    novel_marker_min_coverage: float,
    novel_marker_require_cluster: bool,
    # Region assembly parameters
    initial_window_bp: int,
    initial_window_genes: int,
    min_markers_initial: int,
    extension_kb: int,
    merge_distance: int,
    # Host signature parameter shared with Phase 2.
    boundary_host_signature_min_token_len: int,
    # Workflow configuration
    rebuild_db: bool,
    assembly_mode: str,
    extended_output: bool,
    resume: bool,
    # Threading
    threads: int,
    # Search backend
    search_backend: str,
    # Logger
    logger,
    # Resume config fingerprint (written to early-exit completion manifests)
    config_fingerprint: Optional[str] = None,
    # Set only after schema-v3 marker validation by the orchestrator.
    resume_authorized: bool = False,
    ablation_id: AblationID = AblationID.A0,
) -> dict:
    """
    Phase 1: Seeding (HMM scan -> marker validation -> region assembly).

    This phase identifies candidate viral regions through:
    1. HMM search for viral markers
    2. Diamond validation of HMM hits
    3. Region assembly from validated markers
    4. Taxonomy expansion for low-marker regions (optional)
    5. Composition-based clustering and expansion

    Args:
        masked_path: Path to masked genome FASTA (from Phase 0)
        proteome_path: Path to protein FASTA (from Phase 0)
        repeat_regions: List of RepeatRegion objects (from Phase 0)
        output_dir: Base output directory
        genome_id: Genome identifier for logging
        hmm_database: Path to HMM database
        hmm_allowlist: Path to HMM allowlist file
        ... (see function signature for all parameters)
        logger: Logger instance

    Returns:
        dict with keys:
            - merged_seeds: List of MergedSeed objects
            - validated_markers: List of validated marker hits
            - host_signature_model: HostSignatureModel for host taxonomy
            - host_signatures: Set of host taxonomy signatures
            - background: Background nucleotide composition model
            - gene_data: Gene prediction data
            - host_deviation_summary: Summary of host taxonomy deviation analysis
            - elapsed: Phase 1 elapsed time in seconds

        Or error dict with keys:
            - genome_id, success=False, error, predictions=0, accepted=0, output_files
    """
    import time

    phase1_start = time.time()
    logger.info("-" * 60)
    logger.info("Phase 1: Seeding (HMM scan -> marker validation -> region assembly)")

    # Use prodigal-gv proteome for HMM query
    hmm_query_fasta = proteome_path

    # Initialize outputs
    validated_markers = []
    host_signatures: set[str] = set()
    host_signature_model = None
    host_signature_model_payload = None
    background = None
    gene_data = None
    host_deviation_summary: dict | None = None
    merged_seeds = []
    frameshift_hits = []

    phase1_state_path = output_dir / "phase1" / PHASE1_STATE_FILENAME

    if resume and resume_authorized:
        from virosync.pipeline.host_signatures import (
            TaxonomyLabelLookup,
            set_taxonomy_lookup,
        )

        # Initialize taxonomy lookup on resume (needed for host signature scoring)
        if taxonomy_labels_file and Path(taxonomy_labels_file).exists():
            tax_lookup = TaxonomyLabelLookup.load(Path(taxonomy_labels_file))
            set_taxonomy_lookup(tax_lookup)

        if not phase1_state_path.is_file():
            raise ValueError(
                "authenticated Phase 1 is missing phase1/resume_state.json"
            )
        confirmed_frameshift_faa = (
            output_dir
            / "phase1"
            / "frameshift_screening"
            / "confirmed_frameshift_proteins.faa"
        )
        if frameshift_screening_enabled and not confirmed_frameshift_faa.is_file():
            raise ValueError(
                "authenticated Phase 1 is missing the confirmed frameshift protein FAA"
            )
        confirmed_frameshift_tsv = confirmed_frameshift_faa.with_name(
            "confirmed_frameshift_markers.tsv"
        )
        if frameshift_screening_enabled and not confirmed_frameshift_tsv.is_file():
            raise ValueError(
                "authenticated Phase 1 is missing the confirmed frameshift marker table"
            )
        logger.info("Phase 1 resume: loading exact cached Phase-1 state")
        phase1_state = load_phase1_state(phase1_state_path)
        validated_markers = phase1_state.validated_markers
        merged_seeds = phase1_state.merged_seeds
        host_signature_model = phase1_state.host_signature_model
        host_signatures = phase1_state.host_signatures
        host_deviation_summary = phase1_state.host_deviation_summary
        logger.info(
            "Phase 1 resume: loaded host signature model (%d tokens)",
            len(host_signature_model.token_weights),
        )

        phase1_elapsed = time.time() - phase1_start
        logger.info(
            "Phase 1 resume: loaded %d markers and %d seed regions",
            len(validated_markers),
            len(merged_seeds),
        )
        logger.info(f"Phase 1 complete: {phase1_elapsed:.1f}s")

        classification_counts = {}
        for seed in merged_seeds:
            cls = seed.predicted_family or "UNKNOWN"
            classification_counts[cls] = classification_counts.get(cls, 0) + 1
        if classification_counts:
            logger.info(f"Seed classifications: {classification_counts}")

        return {
            "merged_seeds": merged_seeds,
            "validated_markers": validated_markers,
            "host_signature_model": host_signature_model,
            "host_signatures": host_signatures,
            "background": background,
            "gene_data": gene_data,
            "host_deviation_summary": host_deviation_summary,
            "elapsed": phase1_elapsed,
        }

    # === HMM-GATED WORKFLOW ===
    logger.info("Using HMM-gated Diamond workflow")

    if not hmm_database or not Path(hmm_database).exists():
        logger.error("HMM-gated workflow requires HMM database")
        return {
            "genome_id": genome_id,
            "success": False,
            "error": "HMM database required for HMM-gated workflow",
            "predictions": 0,
            "accepted": 0,
            "output_files": {},
        }

    if not faa_dir and not marker_db:
        logger.error("HMM-gated workflow requires combined FAA database directory")
        return {
            "genome_id": genome_id,
            "success": False,
            "error": "Combined FAA database required for HMM-gated workflow",
            "predictions": 0,
            "accepted": 0,
            "output_files": {},
        }

    # Get combined FAA database (for Diamond marker validation)
    combined_faa = None
    if marker_db:
        # Use pre-built combined FAA directly
        combined_faa = Path(marker_db)
    else:
        # Build derived databases inside this run, never in fingerprinted sources.
        derived_database_dir = output_dir / "phase1" / "database"
        marker_faa = None
        if marker_faa_db:
            marker_faa = Path(marker_faa_db)
        else:
            marker_faa_dir_path = Path(marker_faa_dir) if marker_faa_dir else None
            if not marker_faa_dir_path:
                logger.error("HMM-gated workflow requires marker FAA DB or directory for validation")
                return {
                    "genome_id": genome_id,
                    "success": False,
                    "error": "Marker FAA DB or directory required for HMM-gated workflow",
                    "predictions": 0,
                    "accepted": 0,
                    "output_files": {},
                }
            marker_faa = build_marker_faa(
                marker_faa_dir_path,
                derived_database_dir / "marker.faa",
                logger,
                rebuild=True,
            )
            if not marker_faa:
                return {
                    "genome_id": genome_id,
                    "success": False,
                    "error": "Failed to build marker.faa for HMM-gated workflow",
                    "predictions": 0,
                    "accepted": 0,
                    "output_files": {},
                }

        combined_faa = ensure_combined_faa(
            Path(faa_dir),
            logger,
            marker_faa=marker_faa,
            rebuild=True,
            output_path=derived_database_dir / "combined.faa",
        )
    if not combined_faa:
        return {
            "genome_id": genome_id,
            "success": False,
            "error": "Failed to build combined.faa for HMM-gated workflow",
            "predictions": 0,
            "accepted": 0,
            "output_files": {},
        }

    if frameshift_screening_enabled:
        from virosync.pipeline.phase1.frameshift_screening import (
            write_confirmed_frameshift_markers,
        )

        frameshift_hits = call_task(
            frameshift_screening_task,
            masked_fasta=masked_path,
            hmm_database=Path(hmm_database),
            output_dir=output_dir / "phase1" / "frameshift_screening",
            threads=threads,
        )
        logger.info("Frameshift screening event-bearing hits: %d", len(frameshift_hits))
        confirmed_frameshift_faa = (
            output_dir
            / "phase1"
            / "frameshift_screening"
            / "confirmed_frameshift_proteins.faa"
        )
        confirmed_frameshift_faa.parent.mkdir(parents=True, exist_ok=True)
        with atomic_write_context(confirmed_frameshift_faa, "w"):
            pass
        write_confirmed_frameshift_markers(
            [],
            confirmed_frameshift_faa.with_name("confirmed_frameshift_markers.tsv"),
        )

    # Step 1: HMM search (raw hits only)
    logger.info("Step 1: Running HMM search...")
    hhg_seeds, hhg_hits = call_task(
        hhg_seeding_task,
        proteome_path=proteome_path,
        hmm_database=hmm_database,
        hmm_query_fasta=hmm_query_fasta,
        genome_fasta=masked_path,
        hmm_allowlist=hmm_allowlist,
        marker_faa_dir=None,
        marker_db=None,
        window_size=20000,
        threads=threads,
        assembly_mode=assembly_mode if assembly_mode != "default" else None,
        hmm_chunk_size=hmm_chunk_size,
        output_dir=output_dir / "phase1" / "hmm",
    )
    logger.info(f"HMM hits: {len(hhg_hits)}")

    if not hhg_hits and not frameshift_hits:
        logger.warning(
            "No protein HMM or frameshift hits found - pipeline complete with 0 predictions"
        )
        elapsed_sec = time.time() - phase1_start
        output_files_empty = call_task(
            generate_outputs_task,
            verification_results=[],
            output_dir=output_dir,
            genome_path=masked_path,
            proteome_path=proteome_path,
            accepted_only=True,
            extended_output=extended_output,
        )
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
            reason="no protein HMM or frameshift hits",
            elapsed_sec=elapsed_sec,
            output_files=output_files,
            fingerprint=config_fingerprint,
        )
        return {
            "genome_id": genome_id,
            "success": True,
            **_empty_prediction_summary(),
            "output_files": output_files,
            "elapsed_sec": elapsed_sec,
            "ablation_counts": InterventionCounts(),
        }

    # Import here to avoid circular dependency
    from virosync.orchestration.tasks import (
        marker_validation_task,
        region_assembly_task,
    )
    from virosync.pipeline.phase1.marker_validation import (
        NovelMarkerCriteria,
        VALIDATED_PREFIXES,
        VALIDATION_MIN_PIDENT,
        collect_host_signatures,
    )
    from virosync.pipeline.host_signatures import (
        build_host_signature_model,
        TaxonomyLabelLookup,
        set_taxonomy_lookup,
        summarize_host_signature_model,
        summarize_host_signature_bits,
    )

    # Initialize taxonomy lookup for full lineage-based host signature comparison
    tax_lookup = None
    if taxonomy_labels_file and Path(taxonomy_labels_file).exists():
        tax_lookup = TaxonomyLabelLookup.load(Path(taxonomy_labels_file))
        set_taxonomy_lookup(tax_lookup)
        logger.info(
            "Loaded taxonomy labels: %d entries from %s",
            len(tax_lookup),
            taxonomy_labels_file,
        )

    # Step 2-3: Marker validation (HMM-gated Diamond)
    logger.info("Steps 2-3: Running marker validation (HMM-gated Diamond)...")
    validated_markers = call_task(
        marker_validation_task,
        hmm_hits=hhg_hits,
        proteome_path=hmm_query_fasta or proteome_path,
        marker_db=combined_faa,
        output_dir=output_dir / "phase1" / "marker_validation",
        genome_path=masked_path,
        threads=threads,
        taxonomy_labels_file=Path(taxonomy_labels_file) if taxonomy_labels_file else None,
        taxonomy_weight_mode=taxonomy_weight_mode,
        search_backend=search_backend,
        max_seqs=marker_validation_top_k,
        novel_criteria=NovelMarkerCriteria(
            min_hmm_score=novel_marker_min_score,
            min_hmm_coverage=novel_marker_min_coverage,
            require_cluster=novel_marker_require_cluster,
        ),
    )
    if frameshift_hits:
        from Bio import SeqIO

        from virosync.pipeline.phase1.frameshift_screening import (
            rescued_protein_id,
            select_confirmed_frameshift_markers,
            write_confirmed_frameshift_faa,
            write_confirmed_frameshift_markers,
        )
        from virosync.pipeline.phase1.hhg_seeding import HMMHit

        frameshift_dir = output_dir / "phase1" / "frameshift_screening"
        candidate_faa = frameshift_dir / "frameshift_candidates.faa"
        candidate_lengths = {
            record.id: len(record.seq) for record in SeqIO.parse(candidate_faa, "fasta")
        }
        frameshift_by_protein_id = {
            rescued_protein_id(hit): hit for hit in frameshift_hits
        }
        rescue_hmm_hits = [
            HMMHit(
                query_name=protein_id,
                target_name=hit.query_name,
                score=hit.score,
                evalue=hit.evalue,
                domain_score=hit.score,
                query_start=1,
                query_end=candidate_lengths[protein_id],
            )
            for protein_id, hit in frameshift_by_protein_id.items()
        ]
        frameshift_validation_dir = frameshift_dir / "validation"
        rescue_validation = call_task(
            marker_validation_task,
            hmm_hits=rescue_hmm_hits,
            proteome_path=candidate_faa,
            marker_db=combined_faa,
            output_dir=frameshift_validation_dir,
            genome_path=masked_path,
            threads=threads,
            taxonomy_labels_file=(
                Path(taxonomy_labels_file) if taxonomy_labels_file else None
            ),
            taxonomy_weight_mode=taxonomy_weight_mode,
            search_backend=search_backend,
            max_seqs=marker_validation_top_k,
            novel_criteria=NovelMarkerCriteria(
                min_hmm_score=novel_marker_min_score,
                min_hmm_coverage=novel_marker_min_coverage,
                require_cluster=novel_marker_require_cluster,
            ),
        )
        confirmed_frameshift_markers = select_confirmed_frameshift_markers(
            rescue_validation,
            frameshift_by_protein_id,
            frameshift_validation_dir / "diamond_top10.tsv",
            validated_prefixes=VALIDATED_PREFIXES,
            min_pident=VALIDATION_MIN_PIDENT,
        )
        confirmed_faa = frameshift_dir / "confirmed_frameshift_proteins.faa"
        confirmed_proteins = write_confirmed_frameshift_faa(
            candidate_faa,
            confirmed_frameshift_markers,
            confirmed_faa,
        )
        write_confirmed_frameshift_markers(
            confirmed_frameshift_markers,
            frameshift_dir / "confirmed_frameshift_markers.tsv",
        )
        validated_markers.extend(confirmed_frameshift_markers)
        logger.info(
            "Frameshift rescue validation: candidates=%d diamond_validated=%d confirmed_loci=%d faa_records=%d",
            len(frameshift_hits),
            sum(
                marker.validation_status == "validated"
                for marker in rescue_validation
            ),
            len(confirmed_frameshift_markers),
            confirmed_proteins,
        )
    validated_counts = {
        "validated": 0,
        "validated_novel": 0,
        "supported": 0,
        "unvalidated": 0,
    }
    for hit in validated_markers:
        status = getattr(hit, "validation_status", "")
        if status in validated_counts:
            validated_counts[status] += 1
    logger.info(
        "Marker validation totals: total=%d validated=%d validated_novel=%d supported=%d unvalidated=%d",
        len(validated_markers),
        validated_counts["validated"],
        validated_counts["validated_novel"],
        validated_counts["supported"],
        validated_counts["unvalidated"],
    )

    # Write extended taxonomy with full lineage strings (if taxonomy lookup available)
    if tax_lookup:
        from virosync.pipeline.phase1.marker_validation import write_extended_taxonomy
        marker_val_dir = output_dir / "phase1" / "marker_validation"
        diamond_top10_tsv = marker_val_dir / "diamond_top10.tsv"
        if diamond_top10_tsv.exists():
            write_extended_taxonomy(diamond_top10_tsv, marker_val_dir, tax_lookup)

    single_marker_min_score = get_assembly_mode(
        assembly_mode
    ).single_marker_min_score
    marker_roles = {
        id(hit): decide_marker_hit_role(
            hit,
            ablation_id=ablation_id,
            single_marker_min_score=single_marker_min_score,
        )
        for hit in validated_markers
    }
    production_validated_markers = [
        hit
        for hit in validated_markers
        if marker_roles[id(hit)].is_production_validated
    ]
    validated_only = [
        hit
        for hit in validated_markers
        if marker_roles[id(hit)].is_retained_evidence
    ]
    tier1_bypassed_markers = [
        hit
        for hit in validated_markers
        if marker_roles[id(hit)].is_tier1_bypassed
    ]
    phase1_ablation_counts = InterventionCounts(
        opportunities=len(tier1_bypassed_markers),
        interventions=len(tier1_bypassed_markers),
        changed=0,
    )
    host_signatures = collect_host_signatures(
        validated_markers,
        host_prefixes=set(host_prefixes),
    )
    host_signature_model = build_host_signature_model(
        validated_markers,
        min_token_length=boundary_host_signature_min_token_len,
        host_prefixes=set(host_prefixes),
        weight_mode=taxonomy_weight_mode,
    )
    host_signature_model_payload = host_signature_model.to_dict() if host_signature_model else None
    if host_signature_model:
        model_path = output_dir / "phase1" / "marker_validation" / "host_signature_model.json"
        model_path.parent.mkdir(parents=True, exist_ok=True)
        with atomic_write_context(model_path, "w") as handle:
            json.dump(host_signature_model_payload, handle, indent=2)
    if host_signature_model and host_signature_model.token_weights:
        logger.info(
            "Host signature model: %d tokens (max_weight=%.2f)",
            len(host_signature_model.token_weights),
            host_signature_model.max_weight,
        )
        summary = summarize_host_signature_model(host_signature_model, top_k=20)
        if summary:
            logger.info(
                "Host signature model details (top %d tokens, min_len=%d):",
                len(summary),
                host_signature_model.min_token_length,
            )
            for token, weight, count in summary:
                logger.info("  %s\t%.4f\t%d", token, weight, count)
        debug_summary = summarize_host_signature_bits(
            host_signature_model, top_k=50, max_bits=10
        )
        if debug_summary:
            debug_path = output_dir / "phase1" / "marker_validation" / "host_signature_model_debug.tsv"
            with debug_path.open("w") as debug_handle:
                debug_handle.write("token\tweight_sum\thit_count\ttop_weights\n")
                for token, weight, count, bits in debug_summary:
                    bits_str = ",".join(f"{b:.1f}" for b in bits)
                    debug_handle.write(f"{token}\t{weight:.2f}\t{count}\t{bits_str}\n")
            logger.info("Host signature debug summary written to %s", debug_path)
    else:
        logger.info("Host signature model: no tokens (skipping)")
    host_deviation_markers: list = []
    seedable_deviation_markers: list = []
    if host_taxonomy_deviation_enabled:
        from virosync.pipeline.phase1.taxonomy_expansion import (
            identify_host_deviation_markers,
            filter_deviation_seed_markers,
        )
        from virosync.pipeline.phase1.region_assembly import assemble_candidate_regions
        deviation_dir = output_dir / "phase1" / "marker_validation"
        deviation_regions = None
        if host_taxonomy_deviation_window_count > 0:
            deviation_regions = assemble_candidate_regions(
                validated_hits=validated_only,
                genome_fasta=masked_path,
                proteome_fasta=proteome_path,
                output_dir=output_dir / "phase1" / "region_assembly_pre",
                initial_window_bp=initial_window_bp,
                initial_window_genes=initial_window_genes,
                min_markers_initial=min_markers_initial,
                extension_kb=extension_kb,
                merge_distance=merge_distance,
                ablation_id=ablation_id,
                single_marker_min_score=single_marker_min_score,
                write_outputs=False,
            )
        deviation_hits = identify_host_deviation_markers(
            marker_hits=validated_markers,
            diamond_top10_tsv=deviation_dir / "diamond_top10.tsv",
            output_dir=deviation_dir,
            taxonomy_lookup=tax_lookup,
            host_prefixes=set(host_prefixes),
            min_token_len=host_taxonomy_deviation_min_token_len,
            max_hits=host_taxonomy_deviation_max_hits,
            min_tokens=host_taxonomy_deviation_min_tokens,
            overlap_threshold=host_taxonomy_deviation_overlap_threshold,
            max_pident=host_taxonomy_deviation_max_pident,
            candidate_regions=deviation_regions,
            genome_fasta=masked_path,
            background_window_bp=host_taxonomy_deviation_window_bp,
            background_window_count=host_taxonomy_deviation_window_count,
            background_window_seed=host_taxonomy_deviation_window_seed,
            background_window_min_markers=host_taxonomy_deviation_window_min_markers,
        )
        host_deviation_markers = deviation_hits or []
        if host_taxonomy_deviation_allow_seeds and host_deviation_markers:
            seedable_deviation_markers = filter_deviation_seed_markers(
                host_deviation_markers,
                window_bp=host_taxonomy_deviation_seed_window_bp,
                min_markers=host_taxonomy_deviation_seed_min_markers,
            )
        baseline_path = deviation_dir / "host_taxonomy_baseline.json"
        baseline_payload = None
        if baseline_path.exists():
            try:
                with baseline_path.open() as handle:
                    baseline_payload = json.load(handle)
            except Exception as e:
                logger.warning(
                    "Failed to load host taxonomy baseline from %s: %s",
                    baseline_path, e
                )
                baseline_payload = None
        host_deviation_summary = {
            "enabled": True,
            "markers_total": len(host_deviation_markers),
            "markers_seedable": len(seedable_deviation_markers),
            "baseline": baseline_payload,
            "report_path": str(deviation_dir / "host_taxonomy_deviation.tsv"),
        }
        logger.info(
            "Host taxonomy deviation markers: %d (allow_seeds=%s)",
            len(host_deviation_markers),
            host_taxonomy_deviation_allow_seeds,
        )

    # Allow all validated markers for seeding (allowlist is informational only)
    logger.info(
        "Markers eligible for seeding (NCLDV/MIRUS in top-5): %d",
        len(validated_only),
    )

    if not validated_only and not (host_taxonomy_deviation_allow_seeds and seedable_deviation_markers):
        logger.warning("No validated markers - pipeline complete with 0 predictions")
        elapsed_sec = time.time() - phase1_start
        output_files_empty = call_task(
            generate_outputs_task,
            verification_results=[],
            output_dir=output_dir,
            genome_path=masked_path,
            proteome_path=proteome_path,
            accepted_only=True,
            extended_output=extended_output,
        )
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
            reason="no validated markers",
            elapsed_sec=elapsed_sec,
            output_files=output_files,
            fingerprint=config_fingerprint,
        )
        return {
            "genome_id": genome_id,
            "success": True,
            **_empty_prediction_summary(),
            "output_files": output_files,
            "elapsed_sec": elapsed_sec,
            "ablation_counts": phase1_ablation_counts,
        }

    seed_markers = validated_only
    extension_markers = host_deviation_markers
    if host_taxonomy_deviation_allow_seeds and host_deviation_markers:
        seed_markers = validated_only + seedable_deviation_markers
        seedable_ids = {m.query_porf for m in seedable_deviation_markers}
        extension_markers = [m for m in host_deviation_markers if m.query_porf not in seedable_ids]

    # Step 4: Region assembly (iterative extension)
    logger.info("Step 4: Assembling candidate regions...")
    candidate_regions = call_task(
        region_assembly_task,
        validated_markers=seed_markers,
        genome_path=masked_path,
        proteome_path=proteome_path,
        output_dir=output_dir / "phase1" / "region_assembly",
        initial_window_bp=initial_window_bp,
        initial_window_genes=initial_window_genes,
        min_markers_initial=min_markers_initial,
        extension_kb=extension_kb,
        merge_distance=merge_distance,
        ablation_id=ablation_id,
        single_marker_min_score=single_marker_min_score,
    )
    log_region_statistics(candidate_regions, logger, "Candidate regions")

    changed_seed_regions = 0
    if ablation_id is AblationID.A2 and tier1_bypassed_markers:
        from virosync.pipeline.phase1.region_assembly import assemble_candidate_regions

        counterfactual_seed_markers = list(production_validated_markers)
        if host_taxonomy_deviation_allow_seeds and host_deviation_markers:
            counterfactual_seed_markers.extend(seedable_deviation_markers)
        counterfactual_regions = assemble_candidate_regions(
            validated_hits=counterfactual_seed_markers,
            genome_fasta=masked_path,
            proteome_fasta=proteome_path,
            output_dir=output_dir / "phase1" / "region_assembly_a0_counterfactual",
            initial_window_bp=initial_window_bp,
            initial_window_genes=initial_window_genes,
            min_markers_initial=min_markers_initial,
            extension_kb=extension_kb,
            merge_distance=merge_distance,
            ablation_id=AblationID.A0,
            single_marker_min_score=single_marker_min_score,
            write_outputs=False,
        )
        changed_seed_regions = len(
            _region_coordinate_surface(candidate_regions)
            - _region_coordinate_surface(counterfactual_regions)
        )

    # Simply convert candidate_regions to MergedSeed format for Phase 2
    # Boundary refinement will handle boundaries (no composition expansion needed)
    from virosync.pipeline.phase1.frameshift_screening import (
        is_rescued_protein_id,
    )
    from virosync.pipeline.phase1.seed_merger import MergedSeed

    logger.info("Converting %d marker-based regions to seeds",
               len(candidate_regions))

    merged_seeds = []
    anchor_markers = _seed_annotation_markers(
        validated_markers,
        seedable_deviation_markers,
        extension_markers,
    )
    for region in candidate_regions:
        anchors = []
        region_anchor_markers = (
            region.markers
            if getattr(region, "predicted_family", "") == "CRESS"
            else anchor_markers
        )
        for marker in region_anchor_markers:
            if marker.scaffold != region.scaffold:
                continue
            if marker.start >= region.end or marker.end <= region.start:
                continue
            porf_id = getattr(marker, "porf_id", None) or getattr(marker, "query_porf", "")
            anchors.append(
                Anchor(
                    porf_id=porf_id,
                    scaffold=marker.scaffold,
                    start=marker.start,
                    end=marker.end,
                    strand=getattr(marker, "strand", "+"),
                    hallmark_gene=marker.hmm_target,
                    score=marker.hmm_score,
                    evalue=getattr(marker, "hmm_evalue", 0.0),
                )
            )
        has_rescued_marker = any(
            is_rescued_protein_id(marker.query_porf) for marker in region.markers
        )
        has_ordinary_marker = any(
            not is_rescued_protein_id(marker.query_porf) for marker in region.markers
        )
        sources = ["hhg", "marker_validation"] if has_ordinary_marker else []
        if has_rescued_marker:
            sources.append("frameshift_rescue")
        merged_seeds.append(MergedSeed(
            scaffold=region.scaffold,
            start=region.start,
            end=region.end,
            sources=sources,
            confidence="high" if any(
                is_mcp_gene(a.hallmark_gene) for a in anchors
            ) else "medium",
            hhg_score=len(anchors) * 10.0,
            novelty_score=0.0,
            compositional_score=0.0,
            mean_kfd=0.0,
            max_kfd=0.0,
            mean_composite=0.0,
            max_composite=0.0,
            gc_deviation=0.0,
            n_windows=0,
            cluster_ids=[],
            anchors=anchors,
            hhg_anchors=anchors,
            predicted_family=getattr(region, "predicted_family", ""),
        ))

    # Assign stable seed_id (same format as seed_merger.py line 328-330)
    for idx, seed in enumerate(merged_seeds):
        seed.seed_id = f"seed_{idx}_{seed.scaffold}_{seed.start}"



    phase1_elapsed = time.time() - phase1_start
    logger.info(f"Phase 1 complete: {phase1_elapsed:.1f}s")

    # Compute region classification for each seed based on marker types
    for seed in merged_seeds:
        if seed.predicted_family != "CRESS":
            seed.compute_classification(None)
    phase1_ablation_counts = InterventionCounts(
        opportunities=phase1_ablation_counts.opportunities,
        interventions=phase1_ablation_counts.interventions,
        changed=changed_seed_regions,
    )
    # Log classification summary
    classification_counts = {}
    for seed in merged_seeds:
        cls = seed.predicted_family or "UNKNOWN"
        classification_counts[cls] = classification_counts.get(cls, 0) + 1
    if classification_counts:
        logger.info(f"Seed classifications: {classification_counts}")

    if host_signature_model is None:
        raise ValueError("Phase 1 did not produce a host signature model")
    write_phase1_state(
        phase1_state_path,
        validated_markers=validated_markers,
        merged_seeds=merged_seeds,
        host_signature_model=host_signature_model,
        host_signatures=host_signatures,
        host_deviation_summary=host_deviation_summary,
    )
    logger.info("Phase 1 resume state written to %s", phase1_state_path)

    return {
        "merged_seeds": merged_seeds,
        "validated_markers": validated_markers,
        "host_signature_model": host_signature_model,
        "host_signatures": host_signatures,
        "background": background,
        "gene_data": gene_data,
        "host_deviation_summary": host_deviation_summary,
        "elapsed": phase1_elapsed,
        "ablation_counts": phase1_ablation_counts,
    }
