"""
ViroSync orchestration task functions.

Each pipeline phase is decomposed into granular Python functions for reuse by
the single-genome orchestrator.

Reviewer Notes Addressed:
- Avoid passing large in-memory objects: tasks load data from disk paths
- result_storage_key interpolation: removed, using default per-run keys
"""

import logging
from pathlib import Path
from typing import Optional

# Import pipeline functions for testing/mocking
from virosync.ablation import AblationID
from virosync.config import MaskingBackend, MaskingConfig
from virosync.pipeline.phase0.masking import (
    MaskingResult,
    mask_genome_pipeline,
    quick_mask as _quick_mask,
)
from virosync.pipeline.phase1.taxonomy_expansion import filter_regions_by_taxonomy_expansion
from virosync.orchestration.runtime import get_orchestration_logger
from virosync.utils.path_safety import require_strict_child, safe_filename_components

# Kept as a public module attribute for callers/tests; it is no longer an automatic
# failure fallback.
quick_mask = _quick_mask


def task(*_args, **_kwargs):
    """No-op decorator preserving the historical task declaration style."""

    def _decorate(func):
        return func

    return _decorate


def _boundary_run_id(boundary) -> str:
    """Return the raw boundary identifier used in maps and logs."""
    return f"{boundary.scaffold}_{boundary.start}_{boundary.end}"


def _preflight_boundary_work_dirs(
    boundaries: list,
    work_dir: Path,
) -> dict[str, Path]:
    """Map boundaries to contained work dirs before any verification worker starts."""
    boundary_ids = [_boundary_run_id(boundary) for boundary in boundaries]
    raw_components = [f"eve_{boundary_id}" for boundary_id in boundary_ids]
    filename_components = safe_filename_components(
        raw_components,
        label="boundary work ID",
    )
    work_dir = Path(work_dir)
    work_dirs: dict[str, Path] = {}
    for boundary_id, raw_component in zip(boundary_ids, raw_components, strict=True):
        candidate = work_dir / filename_components[raw_component]
        require_strict_child(work_dir, candidate)
        work_dirs[boundary_id] = candidate
    return work_dirs


# Helpers

# === Phase 0 Tasks ===


@task(
    name="mask_genome",
    description="Mask repeats in genome using TRF + RepeatMasker",
    retries=2,
    retry_delay_seconds=30,
    persist_result=True,
)
def mask_genome_task(
    genome_path: Path,
    output_dir: Path,
    threads: int = 8,
    masking: Optional[MaskingConfig] = None,
    skip_masking: Optional[bool] = None,
) -> MaskingResult:
    """
    Identify repeats and optionally mask the genome.

    Args:
        genome_path: Path to input genome FASTA
        output_dir: Output directory for masked genome
        threads: Number of threads for masking tools
        masking: Canonical backend, target, and failure policy
        skip_masking: Deprecated compatibility override (True maps to off)

    Returns:
        Immutable masking result with persisted status identity
    """
    from virosync.orchestration.resource_monitor import ResourceMonitor

    logger = get_orchestration_logger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    masking = masking or MaskingConfig()
    if isinstance(masking, dict):
        masking = MaskingConfig(**masking)
    if skip_masking is not None:
        if type(skip_masking) is not bool:
            raise TypeError("skip_masking compatibility override must be a boolean")
        if skip_masking:
            masking = masking.with_backend(MaskingBackend.OFF)
        elif masking.backend is MaskingBackend.OFF:
            masking = masking.with_backend(MaskingBackend.TRF_REPEATMASKER)

    logger.info(
        "mask_genome: input=%s output_dir=%s threads=%s backend=%s policy=%s",
        genome_path,
        output_dir,
        threads,
        masking.backend.value,
        masking.failure_policy.value,
    )

    with ResourceMonitor(
        task_name="mask_genome",
        genome_id=Path(genome_path).stem,
        phase="phase0",
        output_dir=output_dir.parent,
        threads=threads,
        task_id=Path(genome_path).stem,
    ):
        result = mask_genome_pipeline(
            Path(genome_path),
            output_dir / "masking",
            threads=threads,
            config=masking,
        )
        logger.info(
            "mask_genome: status=%s effective_backend=%s masked_bases=%d "
            "repeat_regions=%d output=%s",
            result.status,
            result.effective_backend.value,
            result.masked_bases,
            len(result.repeat_regions),
            result.output_path,
        )
        return result


@task(
    name="generate_proteome",
    description="Run prodigal-gv to generate genome-wide gene predictions",
    retries=1,
    persist_result=True,
)
def generate_proteome_task(
    genome_path: Path,
    output_dir: Path,
    threads: int = 1,
) -> tuple[Path, int]:
    """
    Generate genome-wide gene predictions via prodigal-gv.

    Args:
        genome_path: Path to (masked) genome FASTA
        output_dir: Output directory for proteome
    Returns:
        Tuple of (proteome_path, n_genes)
    """
    from virosync.pipeline.phase0 import run_prodigal_genome
    from virosync.orchestration.resource_monitor import ResourceMonitor

    logger = get_orchestration_logger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    proteome_path = output_dir / "proteome.fasta"
    logger.info(
        "generate_proteome: genome=%s output=%s",
        genome_path,
        proteome_path,
    )
    with ResourceMonitor(
        task_name="generate_proteome",
        genome_id=Path(genome_path).stem,
        phase="phase0",
        output_dir=output_dir.parent,
        threads=threads,
        task_id=Path(genome_path).stem,
    ):
        proteins_path, genes = run_prodigal_genome(
            genome_fasta=Path(genome_path),
            output_dir=output_dir,
            threads=threads,
        )
        n_genes = len(genes)

    logger.info("generate_proteome: n_genes=%s", n_genes)
    return proteome_path, n_genes


# === Phase 1 Tasks (Parallel Seeding Strategies) ===


@task(
    name="frameshift_screening",
    description="Run frameshift-sensitive VS marker rescue screening",
    persist_result=True,
)
def frameshift_screening_task(
    masked_fasta: Path,
    hmm_database: Path,
    output_dir: Path,
    threads: int = 8,
) -> list:
    """Run the optional BATH rescue screen against the masked assembly."""

    from virosync.pipeline.phase1.frameshift_screening import (
        run_frameshift_screening,
    )

    return run_frameshift_screening(
        masked_fasta=Path(masked_fasta),
        hmm_database=Path(hmm_database),
        output_dir=Path(output_dir),
        threads=threads,
    )


@task(
    name="hhg_seeding",
    description="HMM hallmark gene seeding",
    retries=2,
    retry_delay_seconds=60,
    timeout_seconds=3600,  # 1 hour max
    persist_result=True,
)
def hhg_seeding_task(
    proteome_path: Path,
    hmm_database: Path,
    hmm_query_fasta: Optional[Path] = None,
    genome_fasta: Optional[Path] = None,
    hmm_allowlist: Optional[Path] = None,
    marker_faa_dir: Optional[Path] = None,
    marker_db: Optional[Path] = None,
    marker_top_k: int = 10,
    window_size: int = 20000,
    threads: int = 8,
    assembly_mode: Optional[str] = None,
    hmm_chunk_size: Optional[int] = None,
    output_dir: Optional[Path] = None,
) -> tuple[list, list]:
    """
    Run HHG seeding strategy.

    Searches proteome against HMM database to find hallmark genes,
    then creates seed windows around anchors.

    Args:
        proteome_path: Path to proteome FASTA
        hmm_database: Path to HMM database file
        hmm_allowlist: Optional allowlist of HMM names
        marker_faa_dir: Optional directory with marker FAA files
        marker_top_k: Top K markers to use
        window_size: Size of seed windows around anchors
        threads: Number of threads for HMM search
        assembly_mode: Assembly mode (default, fragmented, relaxed, strict)
        hmm_chunk_size: Chunk size for HMM processing

    Returns:
        Tuple of (hhg_seeds, hhg_hits)
    """
    from virosync.pipeline.phase1 import hhg_seeding_pipeline
    from virosync.orchestration.resource_monitor import ResourceMonitor

    logger = get_orchestration_logger(__name__)
    output_dir = Path(output_dir) if output_dir else Path(proteome_path).parent.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    query_fasta = Path(hmm_query_fasta) if hmm_query_fasta else Path(proteome_path)
    logger.info(
        "hhg_seeding: proteome=%s hmm_query=%s hmm_db=%s allowlist=%s marker_faa_dir=%s marker_db=%s output_dir=%s",
        proteome_path,
        query_fasta,
        hmm_database,
        hmm_allowlist,
        marker_faa_dir,
        marker_db,
        output_dir,
    )
    with ResourceMonitor(
        task_name="hhg_seeding",
        genome_id=Path(proteome_path).parent.parent.name,
        phase="phase1",
        output_dir=output_dir,
        threads=threads,
        task_id=Path(proteome_path).parent.parent.name,
    ):
        result = hhg_seeding_pipeline(
            proteome_fasta=query_fasta,
            genome_fasta=Path(genome_fasta) if genome_fasta else None,
            hmm_file=Path(hmm_database),
            hmm_allowlist=Path(hmm_allowlist) if hmm_allowlist else None,
            marker_faa_dir=Path(marker_faa_dir) if marker_faa_dir else None,
            marker_db=Path(marker_db) if marker_db else None,
            marker_top_k=marker_top_k,
            window_size=window_size,
            threads=threads,
            hmm_chunk_size=hmm_chunk_size,
            assembly_mode=assembly_mode,
            return_hits=True,
            output_dir=output_dir,
        )

    # Handle both return formats
    if isinstance(result, tuple):
        return result
    return result, []


# === Phase 1 HMM-Gated Workflow Tasks ===


@task(
    name="marker_validation",
    description="HMM-gated Diamond marker validation",
    retries=2,
    retry_delay_seconds=60,
    persist_result=True,
)
def marker_validation_task(
    hmm_hits: list,
    proteome_path: Path,
    marker_db: Path,
    output_dir: Path,
    genome_path: Optional[Path] = None,
    threads: int = 8,
    evalue: float = 1e-5,
    max_seqs: int = 10,
    novel_criteria=None,
    taxonomy_labels_file: Optional[Path] = None,
    taxonomy_weight_mode: str = "rank",
    search_backend: str = "diamond",
) -> list:
    """
    Run HMM-gated Diamond validation workflow (Steps 2-3).

    This is the key optimization: only run Diamond on HMM-hit sequences,
    not all genes. Validates markers based on top-10 taxonomy (NCLDV/MIRUS in top-10).

    Args:
        hmm_hits: List of HMMHit objects from hhg_seeding_task
        proteome_path: Path to conceptual proteome FASTA (HMM query FASTA)
        marker_db: Path to marker Diamond database (marker.dmnd - TIER 1)
        output_dir: Output directory for validation results
        genome_path: Optional genome FASTA for coordinate mapping
        threads: Number of threads for Diamond
        evalue: E-value cutoff for Diamond
        max_seqs: Top-K hits to retrieve (default 10)
        novel_criteria: Optional HMM-only novel-marker acceptance criteria
        taxonomy_labels_file: Optional path to taxonomy labels TSV for host fingerprinting

    Returns:
        List of ValidatedMarkerHit objects
    """
    from virosync.pipeline.phase1 import (
        extract_hmm_hit_sequences,
        run_diamond_on_hmm_hits,
        filter_validated_markers,
    )
    from virosync.pipeline.host_signatures import TaxonomyLabelLookup
    from virosync.orchestration.resource_monitor import ResourceMonitor

    logger = get_orchestration_logger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load taxonomy lookup if available
    taxonomy_lookup = None
    if taxonomy_labels_file and Path(taxonomy_labels_file).exists():
        logger.info(f"Loading taxonomy lookup from {taxonomy_labels_file}")
        taxonomy_lookup = TaxonomyLabelLookup.load(Path(taxonomy_labels_file))
        logger.info(f"Loaded {len(taxonomy_lookup)} taxonomy labels")

    if not hmm_hits:
        logger.info("marker_validation: skipped (no HMM hits)")
        return []

    logger.info(
        "marker_validation: hits=%s proteome=%s db=%s output_dir=%s threads=%s",
        len(hmm_hits),
        proteome_path,
        marker_db,
        output_dir,
        threads,
    )
    with ResourceMonitor(
        task_name="marker_validation",
        genome_id=Path(proteome_path).parent.parent.name,
        phase="phase1",
        output_dir=output_dir.parent,
        threads=threads,
        task_id=Path(proteome_path).parent.parent.name,
    ):
        # Step 2a: Extract HMM-hit sequences
        tmp_faa_path = output_dir / "hmm_hit_porfs.faa"

        n_extracted = extract_hmm_hit_sequences(
            hmm_hits=hmm_hits,
            proteome_fasta=Path(proteome_path),
            output_fasta=tmp_faa_path,
        )

        if n_extracted == 0:
            return []

        # Step 2b: Run Diamond on HMM-hit sequences only
        # Use --sensitive for Tier-1 marker validation to detect divergent
        # viral homologs in ancient EVE proteins.
        diamond_output = output_dir / "diamond_top10.tsv"
        run_diamond_on_hmm_hits(
            hmm_hit_fasta=tmp_faa_path,
            diamond_db=Path(marker_db),
            output_tsv=diamond_output,
            threads=threads,
            evalue=evalue,
            max_seqs=max_seqs,
            search_backend=search_backend,
            sensitive=True,
        )

        # Step 3: Filter validated markers based on top-10 taxonomy
        validated_markers = filter_validated_markers(
            hmm_hits=hmm_hits,
            diamond_output_tsv=diamond_output,
            proteome_fasta=Path(proteome_path),
            genome_fasta=Path(genome_path) if genome_path else None,
            output_dir=output_dir,
            taxonomy_lookup=taxonomy_lookup,
            taxonomy_weight_mode=taxonomy_weight_mode,
            novel_criteria=novel_criteria,
            max_seqs=max_seqs,
        )

        logger.info("marker_validation: validated=%s", len(validated_markers))
        return validated_markers


@task(
    name="region_assembly",
    description="Iterative marker-driven region assembly",
    retries=1,
    persist_result=True,
)
def region_assembly_task(
    validated_markers: list,
    genome_path: Path,
    proteome_path: Path,
    output_dir: Path,
    initial_window_bp: int = 10000,
    initial_window_genes: int = 5,
    min_markers_initial: int = 2,
    extension_kb: int = 5,
    merge_distance: int = 1000,
    ablation_id: AblationID = AblationID.A0,
    single_marker_min_score: float = 50.0,
) -> list:
    """
    Assemble candidate regions from validated markers (Step 4).

    Uses iterative extension to capture adjacent markers until no more can be found.

    Args:
        validated_markers: List of ValidatedMarkerHit from marker_validation_task
        genome_path: Path to genome FASTA (for scaffold lengths)
        proteome_path: Path to proteome FASTA (for gene-based distances)
        output_dir: Output directory for marker_seed_regions.bed
        initial_window_bp: Max bp between markers for initial clustering
        initial_window_genes: Max genes between markers for initial clustering
        min_markers_initial: Min markers to form initial cluster
        extension_kb: Extension distance from outermost markers (kb)
        merge_distance: Max gap to merge overlapping regions (bp)

    Returns:
        List of CandidateRegion objects
    """
    from virosync.pipeline.phase1 import assemble_candidate_regions
    from virosync.orchestration.resource_monitor import ResourceMonitor

    logger = get_orchestration_logger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not validated_markers:
        logger.info("region_assembly: skipped (no validated markers)")
        return []

    logger.info(
        "region_assembly: validated=%s genome=%s output_dir=%s",
        len(validated_markers),
        genome_path,
        output_dir,
    )
    with ResourceMonitor(
        task_name="region_assembly",
        genome_id=Path(genome_path).parent.parent.name,
        phase="phase1",
        output_dir=output_dir.parent,
        threads=1,
        task_id=Path(genome_path).parent.parent.name,
    ):
        candidate_regions = assemble_candidate_regions(
            validated_hits=validated_markers,
            genome_fasta=Path(genome_path),
            proteome_fasta=Path(proteome_path),
            output_dir=output_dir,
            initial_window_bp=initial_window_bp,
            initial_window_genes=initial_window_genes,
            min_markers_initial=min_markers_initial,
            extension_kb=extension_kb,
            merge_distance=merge_distance,
            ablation_id=ablation_id,
            single_marker_min_score=single_marker_min_score,
        )

        logger.info("region_assembly: candidates=%s", len(candidate_regions))
        return candidate_regions


@task(
    name="taxonomy_expansion",
    description="Validate low-marker regions via flanking gene taxonomy",
    retries=1,
    timeout_seconds=3600,  # 60 min max for large genomes
    persist_result=True,
)
def taxonomy_expansion_task(
    candidate_regions: list,
    proteome_path: Path,
    gene_taxonomy_faa_db: Path,
    marker_count_threshold: int = 1,
    flank_genes: int = 5,
    min_viral_genes_total: int = 3,
    min_viral_genes_non_marker: int = 2,
    short_scaffold_min_fraction: float = 0.20,
    require_multi_family: bool = False,
    batch_diamond: bool = False,
    threads: int = 4,
    output_dir: Optional[Path] = None,
    search_backend: str = "diamond",
) -> list:
    """
    Filter low-marker regions using taxonomy expansion (Step 4.5).

    Extracts ±N genes around validated markers and runs Diamond BLASTP vs
    combined_proteome.dmnd to validate regions with viral flanking gene taxonomy.

    EVE-specific approach:
    - Gene-level: one viral top-10 hit at >=25% identity → viral-positive
    - Region-level: ≥3 viral-positive genes (≥2 non-marker) → MEDIUM/HIGH confidence
    - MCP boost: MCP marker presence increases confidence

    Args:
        candidate_regions: Regions from region_assembly_task
        proteome_path: Path to proteome FASTA
        gene_taxonomy_faa_db: Path to combined_proteome.dmnd (40GB full database)
        marker_count_threshold: Regions with ≤N markers undergo expansion (default 1)
        flank_genes: Number of genes to check on each side (default 5 = 11 total)
        min_viral_genes_total: Min total viral genes (default 3, including marker)
        min_viral_genes_non_marker: Min non-marker viral genes (default 2)
        short_scaffold_min_fraction: Viral fraction for short scaffolds (default 0.20 = 20%)
        require_multi_family: Require ≥2 distinct viral families (default False)
        batch_diamond: Deprecated no-op, retained for signature compatibility.
            The chunked batch implementation was removed as unreachable; DIAMOND
            robustness now comes from search_backend's watchdog and
            reduced-thread retry.
        threads: Number of threads
        output_dir: Output directory for expansion results

    Returns:
        List of filtered CandidateRegion objects (only accepted regions)
    """
    from virosync.pipeline.phase0.prodigal import load_gene_predictions
    from virosync.orchestration.resource_monitor import ResourceMonitor

    logger = get_orchestration_logger(__name__)

    # Fail-fast validation: missing database is a systemic error
    if not gene_taxonomy_faa_db:
        raise ValueError(
            "Gene taxonomy database path not provided. "
            "Cannot perform taxonomy expansion without combined_proteome.dmnd database."
        )

    if not Path(gene_taxonomy_faa_db).exists():
        raise FileNotFoundError(
            f"Gene taxonomy database not found: {gene_taxonomy_faa_db}. "
            "Cannot perform taxonomy expansion without combined_proteome.dmnd database."
        )

    if not candidate_regions:
        logger.info("taxonomy_expansion: skipped (no candidate regions)")
        return []

    if output_dir is None:
        output_dir = Path(proteome_path).parents[1] / "phase1" / "taxonomy_expansion"

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "taxonomy_expansion: regions=%s db=%s flank_genes=%s threshold=%s/%s",
        len(candidate_regions),
        gene_taxonomy_faa_db,
        flank_genes,
        min_viral_genes_total,
        min_viral_genes_non_marker,
    )

    # Load gene predictions for flanking gene extraction
    gene_data = load_gene_predictions(Path(proteome_path))

    with ResourceMonitor(
        task_name="taxonomy_expansion",
        genome_id=Path(proteome_path).parent.parent.name,
        phase="phase1",
        output_dir=output_dir.parent,
        threads=threads,
        task_id=Path(proteome_path).parent.parent.name,
    ):
        filtered_regions, expansion_results = filter_regions_by_taxonomy_expansion(
            candidate_regions=candidate_regions,
            proteome_path=Path(proteome_path),
            gene_taxonomy_db=Path(gene_taxonomy_faa_db),
            gene_data=gene_data,
            marker_count_threshold=marker_count_threshold,
            flank_genes=flank_genes,
            min_viral_genes_total=min_viral_genes_total,
            min_viral_genes_non_marker=min_viral_genes_non_marker,
            short_scaffold_min_fraction=short_scaffold_min_fraction,
            batch_diamond=batch_diamond,
            threads=threads,
            output_dir=output_dir,
            search_backend=search_backend,
        )

        # Log statistics by confidence level
        high_conf = sum(1 for r in expansion_results if r.accepted and r.expansion_confidence == "HIGH")
        medium_conf = sum(1 for r in expansion_results if r.accepted and r.expansion_confidence == "MEDIUM")
        low_conf = sum(1 for r in expansion_results if r.accepted and r.expansion_confidence == "LOW")

        logger.info(
            "taxonomy_expansion: %s total → %s kept (HIGH=%s, MEDIUM=%s, LOW=%s)",
            len(candidate_regions),
            len(filtered_regions),
            high_conf,
            medium_conf,
            low_conf,
        )

        return filtered_regions


@task(
    name="gene_taxonomy_batch",
    description="Batch Diamond taxonomy for all EVE candidates",
)
def gene_taxonomy_batch_task(
    regions: list[dict],
    proteome_path: Path,
    gene_taxonomy_db: Path,
    output_dir: Path,
    threads: int = 4,
    high_pident_host_threshold: float = 70.0,
    host_label: Optional[str] = None,
    marker_validation_dir: Optional[Path] = None,
    hmm_hits_file: Optional[Path] = None,
    search_backend: str = "diamond",
) -> dict[str, tuple[list, dict]]:
    """
    Run Diamond gene taxonomy for all candidates in one batch using genome-wide prodigal genes.

    Uses TIER 2 database (combined_proteome.dmnd) for comprehensive gene taxonomy.
    """
    from virosync.pipeline.phase3.gene_taxonomy import run_gene_taxonomy_diamond_batch
    from virosync.orchestration.utils import run_with_monitor

    return run_with_monitor(
        task_name="gene_taxonomy_batch",
        func=run_gene_taxonomy_diamond_batch,
        monitor_output_dir=output_dir.parent,
        regions=regions,
        proteome_fasta=proteome_path,
        combined_faa_db=gene_taxonomy_db,
        output_dir=output_dir,
        threads=threads,
        high_pident_euk_threshold=high_pident_host_threshold,
        search_backend=search_backend,
    )


@task(
    name="interproscan_batch",
    description="Batch InterProScan annotation for all EVE candidates",
)
def interproscan_batch_task(
    regions: list[dict],
    proteome_path: Path,
    interproscan_dir: Path,
    output_dir: Path,
    threads: int = 4,
    keywords: Optional[list[str]] = None,
    applications: Optional[list[str]] = None,
) -> dict[str, dict]:
    """
    Run InterProScan for all candidates in one batch using genome-wide prodigal genes.
    """
    from virosync.pipeline.phase3.interproscan import run_interproscan_batch
    from virosync.orchestration.utils import run_with_monitor

    return run_with_monitor(
        task_name="interproscan_batch",
        func=run_interproscan_batch,
        monitor_output_dir=output_dir.parent,
        regions=regions,
        proteome_fasta=proteome_path,
        interproscan_dir=interproscan_dir,
        output_dir=output_dir,
        threads=threads,
        keywords=keywords,
        applications=applications,
    )


# === Phase 3 Tasks (Parallel Evidence Synthesis) ===


@task(
    name="verify_eve_candidate",
    description="Evidence synthesis for single EVE candidate",
    retries=1,
    timeout_seconds=600,  # 10 min per candidate
    persist_result=True,
)
def verify_eve_task(
    boundary,
    genome_path: Path,
    work_dir: Path,
    proteome_path: Path,
    hallmark_hits: Optional[list] = None,
    novelty_scores: Optional[dict] = None,
    gene_taxonomy_result: Optional[tuple[list, dict]] = None,
    interproscan_result: Optional[dict] = None,
    euk_host_signatures: Optional[set[str]] = None,
    host_signature_model: Optional[dict] = None,
    host_signature_score_threshold: float = 0.3,
    host_prefixes: Optional[list[str]] = None,
    host_label: str = "EUK",
    high_tier_threshold: float = 0.7,
    low_tier_threshold: float = 0.2,
    use_crf_in_final_score: bool = False,
    priority_marker_list: Optional[list[str]] = None,
    marker_floor_priority_only: float = 0.55,
    marker_floor_priority_plus_family: float = 0.70,
    marker_floor_priority_multi_family: float = 0.80,
    marker_family_bonus_per_family: float = 0.06,
    marker_multi_family_bonus: float = 0.08,
    skip_structural: bool = False,
    use_boltz: bool = False,
    boltz_mcp_only: bool = True,
    boltz_use_msa_server: bool = False,
    boltz_min_seq_len: int = 100,
    boltz_max_seq_len: int = 1000,
    boltz_no_kernels: bool = True,
    use_tmvec_database: bool = False,
    tmvec_databases: Optional[list[str]] = None,
    tmvec_database_dir: Optional[Path] = None,
    tmvec_min_score: float = 0.5,
    tmvec_require_gpu: bool = False,
    device: str = "cuda",
    viral_structure_db: Optional[Path] = None,
    gvclass_db: Optional[Path] = None,
    diamond_db: Optional[Path] = None,
    enable_phylogenetic: bool = False,
    max_porfs: int = 10000,
    hmm_database: Optional[Path] = None,
    precomputed_tmvec: Optional[dict] = None,
    taxonomy_labels_file: Optional[Path] = None,
    ablation_id: AblationID = AblationID.A0,
):
    """
    Verify a single EVE candidate using evidence synthesis.

    Assigns confidence tiers (HIGH/MEDIUM/LOW) from marker, taxonomy,
    composition, optional structural/domain, and phylogenetic evidence.

    Args:
        boundary: RefinedBoundary from Phase 2
        genome_path: Path to genome FASTA
        work_dir: Working directory for intermediate files
        proteome_path: Path to proteome FASTA for gene extraction
        hallmark_hits: HMM hits for evidence
        novelty_scores: Legacy compatibility payload; unused by the active seeding path
        host_prefixes: Prefixes used to label host taxa in gene taxonomy
        host_label: Host label string (EUK/ARC)
        high_tier_threshold: Threshold for HIGH confidence tier (default 0.8)
        low_tier_threshold: Threshold for LOW confidence tier (default 0.4)
        skip_structural: Skip slow structural prediction (Boltz)
        use_boltz: Enable Boltz + FoldSeek structural homology (disabled by default)
        boltz_mcp_only: Only run Boltz on MCP candidates (DJR/SJR)
        boltz_use_msa_server: Use Boltz MSA server (requires internet)
        boltz_min_seq_len: Minimum sequence length for Boltz prediction
        boltz_max_seq_len: Maximum sequence length for Boltz prediction
        boltz_no_kernels: Use --no_kernels flag for safer Boltz execution
        use_tmvec_database: Enable TMVec database search (optional; disabled by default)
        tmvec_databases: TMVec databases to search (bfvd/cath/swissprot)
        tmvec_min_score: TMVec score threshold for structural support
        device: Device for TMVec/Boltz (cuda/cpu)
        viral_structure_db: FoldSeek database path for structural homology
        gvclass_db: Path to GVClass database
        diamond_db: Path to Diamond database
        max_porfs: Maximum genes to return (default 10000, effectively no limit for TMVec)
        hmm_database: Optional HMM database path (used to locate marker annotations)

    Returns:
        VerificationResult with final status
    """
    from virosync.pipeline.phase3 import (
        EvidenceSynthesizer,
        EvidenceSynthesizerConfig,
    )
    from virosync.orchestration.resource_monitor import ResourceMonitor
    from virosync.orchestration.utils import get_genes_for_boundary

    logger = get_orchestration_logger(__name__)
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    if taxonomy_labels_file and Path(taxonomy_labels_file).exists():
        from virosync.pipeline.host_signatures import TaxonomyLabelLookup, set_taxonomy_lookup
        lookup = TaxonomyLabelLookup.load(Path(taxonomy_labels_file))
        set_taxonomy_lookup(lookup)
        logger.info("Loaded taxonomy labels for host signature (%d entries)", len(lookup))

    region_len = boundary.end - boundary.start
    structural_desc = ""
    methods = []
    if use_tmvec_database:
        methods.append("TMVec")
    boltz_enabled = use_boltz and not skip_structural
    if boltz_enabled:
        methods.append("Boltz/FoldSeek")
    if methods:
        structural_desc = f" + structural ({'/'.join(methods)})"
    logger.info(
        "verify_eve: %s:%d-%d (%d bp) [phylo=%s%s]",
        boundary.scaffold,
        boundary.start,
        boundary.end,
        region_len,
        "GVClass/Diamond" if enable_phylogenetic else "disabled",
        structural_desc,
    )
    marker_annotations_path = None
    if hmm_database:
        candidate = Path(hmm_database).parent / "model_annotations_with_interpro.tsv"
        if candidate.exists():
            marker_annotations_path = candidate

    # Taxonomy lookup is already loaded above (cached); no need to reload here.

    config = EvidenceSynthesizerConfig(
        ablation_id=ablation_id,
        high_tier_threshold=high_tier_threshold,
        low_tier_threshold=low_tier_threshold,
        use_crf_in_final_score=use_crf_in_final_score,
        priority_marker_list=priority_marker_list or ["mcp"],
        marker_floor_priority_only=marker_floor_priority_only,
        marker_floor_priority_plus_family=marker_floor_priority_plus_family,
        marker_floor_priority_multi_family=marker_floor_priority_multi_family,
        marker_family_bonus_per_family=marker_family_bonus_per_family,
        marker_multi_family_bonus=marker_multi_family_bonus,
        use_boltz=boltz_enabled,
        boltz_mcp_only=boltz_mcp_only,
        boltz_use_msa_server=boltz_use_msa_server,
        boltz_min_seq_len=boltz_min_seq_len,
        boltz_max_seq_len=boltz_max_seq_len,
        boltz_no_kernels=boltz_no_kernels,
        use_tmvec_database=use_tmvec_database,
        tmvec_databases=tmvec_databases,
        tmvec_database_dir=Path(tmvec_database_dir) if tmvec_database_dir else None,
        tmvec_min_score=tmvec_min_score,
        tmvec_require_gpu=tmvec_require_gpu,
        device=device,
        gvclass_db=Path(gvclass_db) if gvclass_db else None,
        diamond_db=Path(diamond_db) if diamond_db else None,
        use_phylogenetic_validation=enable_phylogenetic,
        marker_annotations_path=marker_annotations_path,
        euk_host_signatures=euk_host_signatures,
        host_signature_model=host_signature_model,
        host_signature_score_threshold=host_signature_score_threshold,
        host_prefixes=host_prefixes,
        host_label=host_label,
    )

    output_dir = work_dir.parent
    with ResourceMonitor(
        task_name="verify_eve",
        genome_id=Path(genome_path).stem,
        phase="phase3",
        output_dir=output_dir,
        threads=1,
        task_id=f"{Path(genome_path).stem}_{boundary.scaffold}_{boundary.start}_{boundary.end}",
    ):
        synthesizer = EvidenceSynthesizer(
            config=config,
            viral_structure_db=Path(viral_structure_db) if viral_structure_db else None,
            genome_path=Path(genome_path),
            work_dir=work_dir,
        )

        # Get overlapping genes for this boundary
        porf_sequences = get_genes_for_boundary(
            proteome_path=Path(proteome_path),
            scaffold=boundary.scaffold,
            start=boundary.start,
            end=boundary.end,
            max_porfs=max_porfs,
        )
        logger.info(
            "  Region has %d genes, %d hallmark markers",
            len(porf_sequences),
            len(hallmark_hits),
        )

        # Get window features from boundary (set in Phase 2)
        window_features = getattr(boundary, "window_features", [])

        gene_taxonomy_records = []
        gene_taxonomy_summary = {}
        if gene_taxonomy_result:
            gene_taxonomy_records, gene_taxonomy_summary = gene_taxonomy_result
        result = synthesizer.verify_eve(
            refined_boundary=boundary,
            window_features=window_features,
            hallmark_hits=hallmark_hits,
            novelty_scores=novelty_scores,
            porf_sequences=porf_sequences,
            gene_taxonomy_records=gene_taxonomy_records,
            gene_taxonomy_summary=gene_taxonomy_summary,
            interproscan_summary=interproscan_result,
            precomputed_tmvec=precomputed_tmvec,
        )

        # Log verification result
        status_str = result.confidence_tier
        logger.info(
            "  → %s (confidence=%.3f, status=%s)",
            status_str,
            result.final_confidence,
            result.status.name if hasattr(result.status, 'name') else result.status,
        )

        result.gene_taxonomy_records = [
            getattr(r, "__dict__", r) for r in gene_taxonomy_records
        ]
        if gene_taxonomy_summary:
            result.gene_taxonomy_total = gene_taxonomy_summary.get("total", 0)
            result.gene_taxonomy_ncldv_top10 = gene_taxonomy_summary.get("ncldv_mirus", 0)
            result.gene_taxonomy_mirus_top10 = 0
            result.gene_taxonomy_vp_plv_top10 = gene_taxonomy_summary.get("vp_plv", 0)
            result.gene_taxonomy_phage_top10 = 0
            result.gene_taxonomy_viral_top10 = gene_taxonomy_summary.get("viral_top10", 0)
            result.gene_taxonomy_cellular = gene_taxonomy_summary.get("high_pident_euk", 0)
            result.gene_taxonomy_unknown = 0
            result.gene_taxonomy_has_ncldv_mirus = gene_taxonomy_summary.get("has_ncldv_mirus", False)
            result.gene_taxonomy_has_vp_plv = gene_taxonomy_summary.get("has_vp_plv", False)
            result.gene_taxonomy_dominant_family = gene_taxonomy_summary.get(
                "dominant_family", "UNKNOWN"
            )
            result.gene_taxonomy_dominant_fraction = gene_taxonomy_summary.get(
                "dominant_fraction", 0.0
            )
            result.gene_count = gene_taxonomy_summary.get("total", 0)
            result.genes_with_ncldv_mirus_top10 = gene_taxonomy_summary.get("ncldv_mirus", 0)
            result.genes_with_vp_plv_top10 = gene_taxonomy_summary.get("vp_plv", 0)
            result.genes_with_high_pident_euk = gene_taxonomy_summary.get("high_pident_euk", 0)
        return result


def _base_porf_id(porf_id: str) -> str:
    """Normalize pORF identifiers by stripping optional domain suffixes."""
    return porf_id.split("|aa", 1)[0] if "|aa" in porf_id else porf_id


def _build_jelly_roll_summary_for_boundary(
    porf_sequences: list[tuple[str, str]],
    jelly_roll_map: Optional[dict[str, list[dict]]],
) -> Optional[dict]:
    """Aggregate DJR/SJR evidence for a single boundary from classified MCP proteins."""
    if not porf_sequences or not jelly_roll_map:
        return None

    best_records: dict[str, dict] = {}
    for porf_id, _sequence in porf_sequences:
        for lookup_id in (porf_id, _base_porf_id(porf_id)):
            for record in jelly_roll_map.get(lookup_id, []):
                protein_id = str(record.get("porf_id") or lookup_id)
                classification = str(record.get("classification") or "UNKNOWN").upper()
                evidence = str(record.get("evidence") or "")
                marker = str(record.get("marker") or "")
                try:
                    confidence = float(record.get("confidence", 0.0))
                except (TypeError, ValueError):
                    confidence = 0.0

                normalized = {
                    "porf_id": protein_id,
                    "classification": classification,
                    "confidence": confidence,
                    "evidence": evidence,
                    "marker": marker,
                }
                prior = best_records.get(protein_id)
                if prior is None or confidence > prior["confidence"]:
                    best_records[protein_id] = normalized

    if not best_records:
        return None

    records = sorted(best_records.values(), key=lambda rec: rec["porf_id"])
    total_mcp = len(records)
    djr_confidences = [r["confidence"] for r in records if r["classification"] == "DJR"]
    sjr_count = sum(1 for r in records if r["classification"] == "SJR")
    djr_count = len(djr_confidences)
    avg_confidence = sum(r["confidence"] for r in records) / total_mcp

    confidence_bonus = 0.0
    if djr_count > 0:
        avg_djr_confidence = sum(djr_confidences) / djr_count
        djr_fraction = djr_count / total_mcp
        # Moderate DJR-only contribution: strongest when all MCPs support DJR.
        confidence_bonus = min(0.15, 0.05 * djr_count * avg_djr_confidence * djr_fraction)

    return {
        "djr_count": djr_count,
        "sjr_count": sjr_count,
        "total_mcp": total_mcp,
        "avg_confidence": avg_confidence,
        "confidence_bonus": confidence_bonus,
        "mcp_proteins": records,
    }


@task(
    name="verify_eve_candidates_batched",
    description="Batch evidence synthesis for all EVE candidates using ThreadPoolExecutor",
    retries=1,
    persist_result=True,
)
def verify_eve_candidates_batched_task(
    boundaries: list,
    genome_path: Path,
    work_dir: Path,
    proteome_path: Path,
    hallmark_hits_map: dict,
    novelty_scores: Optional[dict] = None,
    gene_taxonomy_map: Optional[dict] = None,
    interproscan_map: Optional[dict] = None,
    jelly_roll_map: Optional[dict[str, list[dict]]] = None,
    euk_host_signatures: Optional[set[str]] = None,
    host_signature_model: Optional[dict] = None,
    host_signature_score_threshold: float = 0.3,
    host_prefixes: Optional[list[str]] = None,
    host_label: str = "EUK",
    high_tier_threshold: float = 0.7,
    low_tier_threshold: float = 0.2,
    use_crf_in_final_score: bool = False,
    priority_marker_list: Optional[list[str]] = None,
    marker_floor_priority_only: float = 0.55,
    marker_floor_priority_plus_family: float = 0.70,
    marker_floor_priority_multi_family: float = 0.80,
    marker_family_bonus_per_family: float = 0.06,
    marker_multi_family_bonus: float = 0.08,
    skip_structural: bool = False,
    use_boltz: bool = False,
    boltz_mcp_only: bool = True,
    boltz_use_msa_server: bool = False,
    boltz_min_seq_len: int = 100,
    boltz_max_seq_len: int = 1000,
    boltz_no_kernels: bool = True,
    use_tmvec_database: bool = False,
    tmvec_databases: Optional[list[str]] = None,
    tmvec_database_dir: Optional[Path] = None,
    tmvec_min_score: float = 0.5,
    tmvec_require_gpu: bool = False,
    device: str = "cuda",
    viral_structure_db: Optional[Path] = None,
    gvclass_db: Optional[Path] = None,
    diamond_db: Optional[Path] = None,
    enable_phylogenetic: bool = False,
    max_porfs: int = 10000,
    hmm_database: Optional[Path] = None,
    precomputed_tmvec: Optional[dict] = None,
    taxonomy_labels_file: Optional[Path] = None,
    max_workers: int = 32,
    ablation_id: AblationID = AblationID.A0,
):
    """
    Verify ALL EVE candidates in parallel using ThreadPoolExecutor.

    This keeps candidate verification within one Python task function
    that internally parallelizes with Python threads. I/O-bound operations
    (file reading, gene extraction) release the GIL, allowing true parallelism.

    Args:
        boundaries: List of RefinedBoundary objects from Phase 2
        genome_path: Path to genome FASTA
        work_dir: Base working directory
        proteome_path: Path to proteome FASTA
        hallmark_hits_map: Dict mapping boundary_id → hallmark hits
        novelty_scores: Novelty scores dict
        gene_taxonomy_map: Dict mapping boundary_id → (records, summary)
        interproscan_map: Dict mapping boundary_id → interproscan results
        jelly_roll_map: Dict mapping pORF ID → DJR/SJR classification records
        [remaining args same as verify_eve_task]
        max_workers: Number of parallel threads (default 32)

    Returns:
        List of VerificationResult objects
    """
    from virosync.pipeline.phase3 import (
        EvidenceSynthesizer,
        EvidenceSynthesizerConfig,
    )
    from virosync.orchestration.utils import get_genes_for_boundary
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import logging

    logger = get_orchestration_logger(__name__)
    work_dir = Path(work_dir)
    boundary_work_dirs = _preflight_boundary_work_dirs(boundaries, work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    structural_desc = ""
    methods = []
    if use_tmvec_database:
        methods.append("TMVec")
    boltz_enabled = use_boltz and not skip_structural
    if boltz_enabled:
        methods.append("Boltz/FoldSeek")
    if methods:
        structural_desc = f" + structural ({'/'.join(methods)})"

    logger.info(
        "verify_eve_batched: Processing %d candidates with %d threads [phylo=%s%s]",
        len(boundaries),
        max_workers,
        "GVClass/Diamond" if enable_phylogenetic else "disabled",
        structural_desc,
    )

    # Determine marker annotations path once
    marker_annotations_path = None
    if hmm_database:
        candidate = Path(hmm_database).parent / "model_annotations_with_interpro.tsv"
        if candidate.exists():
            marker_annotations_path = candidate

    # Create shared config (read-only, thread-safe)
    config = EvidenceSynthesizerConfig(
        ablation_id=ablation_id,
        high_tier_threshold=high_tier_threshold,
        low_tier_threshold=low_tier_threshold,
        use_crf_in_final_score=use_crf_in_final_score,
        priority_marker_list=priority_marker_list or ["mcp"],
        marker_floor_priority_only=marker_floor_priority_only,
        marker_floor_priority_plus_family=marker_floor_priority_plus_family,
        marker_floor_priority_multi_family=marker_floor_priority_multi_family,
        marker_family_bonus_per_family=marker_family_bonus_per_family,
        marker_multi_family_bonus=marker_multi_family_bonus,
        use_boltz=boltz_enabled,
        boltz_mcp_only=boltz_mcp_only,
        boltz_use_msa_server=boltz_use_msa_server,
        boltz_min_seq_len=boltz_min_seq_len,
        boltz_max_seq_len=boltz_max_seq_len,
        boltz_no_kernels=boltz_no_kernels,
        use_tmvec_database=use_tmvec_database,
        tmvec_databases=tmvec_databases,
        tmvec_database_dir=Path(tmvec_database_dir) if tmvec_database_dir else None,
        tmvec_min_score=tmvec_min_score,
        tmvec_require_gpu=tmvec_require_gpu,
        device=device,
        gvclass_db=Path(gvclass_db) if gvclass_db else None,
        diamond_db=Path(diamond_db) if diamond_db else None,
        use_phylogenetic_validation=enable_phylogenetic,
        marker_annotations_path=marker_annotations_path,
        euk_host_signatures=euk_host_signatures,
        host_signature_model=host_signature_model,
        host_signature_score_threshold=host_signature_score_threshold,
        host_prefixes=host_prefixes,
        host_label=host_label,
    )

    # Ensure taxonomy lookup is set for host signature fallback scoring
    if taxonomy_labels_file and Path(taxonomy_labels_file).exists():
        from virosync.pipeline.host_signatures import TaxonomyLabelLookup, set_taxonomy_lookup
        lookup = TaxonomyLabelLookup.load(Path(taxonomy_labels_file))
        set_taxonomy_lookup(lookup)
        logger.info("Batched: loaded taxonomy labels for host signature (%d entries)", len(lookup))

    def verify_single_boundary(boundary):
        """Worker function to verify a single boundary."""
        boundary_id = _boundary_run_id(boundary)

        try:
            # Each thread creates its own synthesizer (GPU operations may not be thread-safe)
            boundary_work_dir = boundary_work_dirs[boundary_id]
            boundary_work_dir.mkdir(parents=True, exist_ok=True)

            synthesizer = EvidenceSynthesizer(
                config=config,
                viral_structure_db=Path(viral_structure_db) if viral_structure_db else None,
                genome_path=Path(genome_path),
                work_dir=boundary_work_dir,
            )

            # Get overlapping genes for this boundary
            porf_sequences = get_genes_for_boundary(
                proteome_path=Path(proteome_path),
                scaffold=boundary.scaffold,
                start=boundary.start,
                end=boundary.end,
                max_porfs=max_porfs,
            )

            # Get boundary-specific data from maps
            hallmark_hits = hallmark_hits_map.get(boundary_id, [])
            gene_taxonomy_result = gene_taxonomy_map.get(boundary_id) if gene_taxonomy_map else None
            interproscan_result = interproscan_map.get(boundary_id) if interproscan_map else None
            jelly_roll_summary = _build_jelly_roll_summary_for_boundary(
                porf_sequences,
                jelly_roll_map,
            )

            gene_taxonomy_records = []
            gene_taxonomy_summary = {}
            if gene_taxonomy_result:
                gene_taxonomy_records, gene_taxonomy_summary = gene_taxonomy_result

            # Get window features from boundary
            window_features = getattr(boundary, "window_features", [])

            # Run verification
            result = synthesizer.verify_eve(
                refined_boundary=boundary,
                window_features=window_features,
                hallmark_hits=hallmark_hits,
                novelty_scores=novelty_scores,
                porf_sequences=porf_sequences,
                gene_taxonomy_records=gene_taxonomy_records,
                gene_taxonomy_summary=gene_taxonomy_summary,
                interproscan_summary=interproscan_result,
                jelly_roll_summary=jelly_roll_summary,
                precomputed_tmvec=precomputed_tmvec,
            )

            # Populate gene taxonomy fields
            result.gene_taxonomy_records = [
                getattr(r, "__dict__", r) for r in gene_taxonomy_records
            ]
            if gene_taxonomy_summary:
                result.gene_taxonomy_total = gene_taxonomy_summary.get("total", 0)
                result.gene_taxonomy_ncldv_top10 = gene_taxonomy_summary.get("ncldv_mirus", 0)
                result.gene_taxonomy_mirus_top10 = 0
                result.gene_taxonomy_vp_plv_top10 = gene_taxonomy_summary.get("vp_plv", 0)
                result.gene_taxonomy_phage_top10 = 0
                result.gene_taxonomy_viral_top10 = gene_taxonomy_summary.get("viral_top10", 0)
                result.gene_taxonomy_cellular = gene_taxonomy_summary.get("high_pident_euk", 0)
                result.gene_taxonomy_unknown = 0
                result.gene_taxonomy_has_ncldv_mirus = gene_taxonomy_summary.get("has_ncldv_mirus", False)
                result.gene_taxonomy_has_vp_plv = gene_taxonomy_summary.get("has_vp_plv", False)
                result.gene_taxonomy_dominant_family = gene_taxonomy_summary.get("dominant_family", "UNKNOWN")
                result.gene_taxonomy_dominant_fraction = gene_taxonomy_summary.get("dominant_fraction", 0.0)
                result.gene_count = gene_taxonomy_summary.get("total", 0)
                result.genes_with_ncldv_mirus_top10 = gene_taxonomy_summary.get("ncldv_mirus", 0)
                result.genes_with_vp_plv_top10 = gene_taxonomy_summary.get("vp_plv", 0)
                result.genes_with_high_pident_euk = gene_taxonomy_summary.get("high_pident_euk", 0)

            # Log result
            region_len = boundary.end - boundary.start
            logging.info(
                "verify_eve: %s:%d-%d (%d bp) → %s (confidence=%.3f, genes=%d, markers=%d)",
                boundary.scaffold,
                boundary.start,
                boundary.end,
                region_len,
                result.confidence_tier,
                result.final_confidence,
                len(porf_sequences),
                len(hallmark_hits),
            )

            return result

        except Exception as e:
            logger.exception("Failed to verify boundary %s", boundary_id)
            raise RuntimeError(f"Failed to verify boundary {boundary_id}") from e

    # Process all boundaries in parallel using ThreadPoolExecutor
    results = []
    completed = 0
    total = len(boundaries)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all boundaries
        future_to_boundary = {
            executor.submit(verify_single_boundary, boundary): boundary
            for boundary in boundaries
        }

        # Collect results as they complete
        for future in as_completed(future_to_boundary):
            boundary = future_to_boundary[future]
            result = future.result()
            if result is not None:
                results.append(result)
            completed += 1

            # Progress logging every 10 completed or at end
            if completed % 10 == 0 or completed == total:
                logger.info(
                    "verify_eve_batched: Progress %d/%d (%.1f%%)",
                    completed,
                    total,
                    100.0 * completed / total,
                )

    # Thread completion order is not reproducible between runs. Sort once here
    # so every downstream consumer (BED, GFF3, evidence_profiles.json) inherits
    # the same stable order.
    results.sort(key=lambda r: (r.scaffold, r.start, r.end, r.eve_id))

    logger.info(
        "verify_eve_batched: Completed %d/%d candidates (%d succeeded)",
        completed,
        total,
        len(results),
    )

    return results


@task(
    name="classify_jelly_roll",
    description="Classify MCP proteins as DJR/SJR using multi-signal approach",
    persist_result=True,
)
def classify_jelly_roll_task(
    marker_hits_path: Path,
    sequences_path: Path,
    output_path: Path,
    interproscan_path: Optional[Path] = None,
    tmvec_results_path: Optional[Path] = None,
    foldseek_results_path: Optional[Path] = None,
) -> Path:
    """
    Classify MCP proteins as DJR (Double Jelly Roll) or SJR (Single Jelly Roll).

    Uses multi-signal classification approach:
    1. InterProScan domain counting (highest confidence)
    2. Multiple HMM domain hits detection
    3. TMVec reference similarity
    4. FoldSeek structural hits
    5. Length heuristics (fallback)

    Args:
        marker_hits_path: Path to validated_marker_hits.tsv from phase1
        sequences_path: Path to hmm_hit_porfs.faa containing HMM-hit sequences
        output_path: Output TSV file for jelly roll classifications
        interproscan_path: Optional path to InterProScan batch results
        tmvec_results_path: Optional path to TMVec results TSV
        foldseek_results_path: Optional path to FoldSeek results TSV

    Returns:
        Path to output TSV file
    """
    import sys
    from pathlib import Path as PathLib
    sys.path.insert(0, str(PathLib(__file__).parent.parent.parent.parent / "scripts"))

    from classify_jelly_roll import (
        load_marker_hits,
        load_sequences,
        load_interproscan_domains,
        load_tmvec_results,
        load_foldseek_results,
        count_hmm_domains_per_protein,
        classify_proteins,
        write_results,
    )
    from virosync.orchestration.resource_monitor import ResourceMonitor

    logger = get_orchestration_logger(__name__)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(
        "classify_jelly_roll: marker_hits=%s sequences=%s output=%s",
        marker_hits_path,
        sequences_path,
        output_path,
    )

    with ResourceMonitor(
        task_name="classify_jelly_roll",
        genome_id=Path(marker_hits_path).parent.parent.parent.name,
        phase="phase3",
        output_dir=output_path.parent,
        threads=1,
        task_id=Path(marker_hits_path).parent.parent.parent.name,
    ):
        # Load marker hits (MCP markers only)
        marker_hits = load_marker_hits(Path(marker_hits_path), mcp_only=True)
        logger.info(f"Loaded {len(marker_hits)} MCP marker hits")

        if not marker_hits:
            logger.warning("No MCP marker hits found, writing empty output")
            write_results([], output_path)
            return output_path

        # Count HMM domains per protein
        hmm_domain_counts = count_hmm_domains_per_protein(marker_hits)
        multi_domain = sum(1 for c in hmm_domain_counts.values() if c > 1)
        logger.info(f"Found {multi_domain} proteins with multiple HMM domains")

        # Load optional signal sources
        interproscan_domains = None
        if interproscan_path and Path(interproscan_path).exists():
            interproscan_domains = load_interproscan_domains(Path(interproscan_path))
            logger.info(f"Loaded {len(interproscan_domains)} InterProScan domain annotations")

        tmvec_results = None
        if tmvec_results_path and Path(tmvec_results_path).exists():
            tmvec_results = load_tmvec_results(Path(tmvec_results_path))
            logger.info(f"Loaded TMVec scores for {len(tmvec_results)} proteins")

        foldseek_results = None
        if foldseek_results_path and Path(foldseek_results_path).exists():
            foldseek_results = load_foldseek_results(Path(foldseek_results_path))
            logger.info(f"Loaded FoldSeek results for {len(foldseek_results)} queries")

        # Load sequences
        porf_ids = set(marker_hits.keys())
        sequences = load_sequences(Path(sequences_path), porf_ids)
        logger.info(f"Loaded {len(sequences)} sequences")

        # Classify proteins
        classifications = classify_proteins(
            marker_hits,
            sequences,
            interproscan_domains=interproscan_domains,
            hmm_domain_counts=hmm_domain_counts,
            tmvec_results=tmvec_results,
            foldseek_results=foldseek_results,
        )

        # Log summary
        djr_count = sum(1 for c in classifications if c.jelly_roll_type == "DJR")
        sjr_count = sum(1 for c in classifications if c.jelly_roll_type == "SJR")
        logger.info(f"Classified {len(classifications)} proteins: {djr_count} DJR, {sjr_count} SJR")

        # Write results
        write_results(classifications, output_path)

        return output_path


@task(
    name="generate_outputs",
    description="Generate final output files (BED, GFF3, TSV)",
    persist_result=True,
)
def generate_outputs_task(
    verification_results: list,
    output_dir: Path,
    genome_path: Optional[Path] = None,
    proteome_path: Optional[Path] = None,
    accepted_only: bool = False,
    extended_output: bool = True,
    seed_marker_allowlist: Optional[list[str]] = None,
    export_all_eve_sequences: bool = False,
    canonical_results: Optional[list] = None,
    promoted_low_results: Optional[list] = None,
) -> dict:
    """
    Generate output files from verification results.

    Produces BED, GFF3, TSV, and optionally FASTA files.

    Args:
        verification_results: List of VerificationResult objects
        output_dir: Output directory for files
        genome_path: Path to genome for sequence extraction
        proteome_path: Path to proteome for sequence extraction
        accepted_only: Only include accepted results (default: False, include all)
        extended_output: Include extended TSV columns for new evidence fields
        seed_marker_allowlist: Optional allowlist of seed markers for output reporting
        canonical_results: Explicit canonical subset already selected by the
            Phase-3 acceptance policy.
        promoted_low_results: Identity-preserving subset of canonical results
            that the normal quality gate promoted from LOW confidence.

    Returns:
        Dictionary mapping output type to file path
    """
    from virosync.pipeline.phase3 import OutputGenerator

    logger = get_orchestration_logger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "generate_outputs: output_dir=%s accepted_only=%s",
        output_dir,
        accepted_only,
    )
    output_gen = OutputGenerator(
        output_dir=output_dir,
        genome_fasta=Path(genome_path) if genome_path else None,
        proteome_fasta=Path(proteome_path) if proteome_path else None,
        extended_output=extended_output,
        seed_marker_allowlist=seed_marker_allowlist,
        export_all_eve_sequences=export_all_eve_sequences,
    )

    return output_gen.generate_all(
        verification_results,
        accepted_only=accepted_only,
        canonical_results=canonical_results,
        promoted_low_results=promoted_low_results,
    )


# === Artifact Generation Tasks ===


@task(
    name="create_summary_artifact",
    description="Create markdown summary artifact",
    persist_result=False,
)
def create_summary_artifact_task(
    genome_id: str,
    n_seeds: int,
    n_boundaries: int,
    n_verified: int,
    tier_counts: dict[str, int],
    output_files: dict,
) -> str:
    """
    Create a markdown summary artifact placeholder.

    Args:
        genome_id: Genome identifier
        n_seeds: Number of seeds found
        n_boundaries: Number of boundaries refined
        n_verified: Number of candidates verified
        tier_counts: Dictionary with HIGH/MEDIUM/LOW tier counts
        output_files: Dictionary of output file paths

    Returns:
        Artifact key
    """
    import re

    def _sanitize_key(value: str) -> str:
        sanitized = re.sub(r"[^a-z0-9-]+", "-", value.lower()).strip("-")
        return sanitized or "virosync-results"

    high = tier_counts.get("HIGH", 0)
    medium = tier_counts.get("MEDIUM", 0)
    low = tier_counts.get("LOW", 0)

    summary = f"""
## ViroSync Results: {genome_id}

| Metric | Value |
|--------|-------|
| Seeds found | {n_seeds} |
| Boundaries refined | {n_boundaries} |
| Candidates verified | {n_verified} |
| **HIGH confidence** | **{high}** |
| **MEDIUM confidence** | **{medium}** |
| **LOW confidence** | **{low}** |

### Output Files
"""
    for file_type, file_path in output_files.items():
        summary += f"- {file_type.upper()}: `{file_path}`\n"

    artifact_key = _sanitize_key(f"virosync-{genome_id}")
    logging.getLogger(__name__).debug(
        "Summary artifact key for %s: %s (%d chars)",
        genome_id,
        artifact_key,
        len(summary),
    )

    return artifact_key
