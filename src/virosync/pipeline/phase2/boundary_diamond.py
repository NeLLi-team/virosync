"""
Batched Diamond BLAST for boundary refinement with host control sampling.

Runs gene taxonomy Diamond ONCE per genome on:
1. All pORFs within ALL candidate EVE regions
2. +/-10 genes flanking each boundary (for boundary context)
3. Control genes sampled from regions outside ALL EVEs and boundaries

Key optimizations:
- Single Diamond run for all seeds (Diamond is fastest with many queries)
- Chunked at 10k proteins per query file
- Top-10 hits for robust neighbor analysis
- Deterministic sampling for reproducibility

This module is part of Phase 2 boundary refinement, implementing the batched
Diamond BLAST strategy from DIAMOND_BLAST_REFACTOR_PLAN.md Section 3.

Key design principles:
1. Run Diamond ONCE per genome, not per-seed
2. Collect ALL proteins from ALL candidate EVE regions + boundaries
3. Add control genes (total, not per-EVE) from outside all EVE/boundary regions
4. Single query file -> single Diamond run
5. Chunk at 10k proteins for optimal Diamond performance
6. Top-10 hits per query protein
7. Include +/-10 genes around each EVE (for boundary context)
8. Re-filter taxonomy after trimming to refined ORF set
"""

from __future__ import annotations

import logging
import random
import subprocess
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from Bio import SeqIO

from virosync.pipeline.phase0.prodigal import parse_prodigal_header
from virosync.utils.atomic_write import atomic_write_context
from virosync.pipeline.taxonomy_utils import (
    TaxonomyFingerprint,
    aggregate_taxonomy_substrings,
    compute_hit_weight,
    resolve_org_id,
)

logger = logging.getLogger(__name__)


def _get_logger() -> logging.Logger:
    return logger


# Taxonomy prefix definitions with consistent __ suffix
VIRAL_PREFIXES = {"NCLDV__", "MIRUS__", "VP__", "PLV__", "PPV__", "CRESS__", "GVMAG__", "PHAGE__"}
CELLULAR_PREFIXES = {"EUK__", "BAC__", "ARC__"}

# Minimum identity for a viral top-10 hit to count as viral evidence.
MIN_VIRAL_HIT_PIDENT = 25.0


def _listify_top10_field(value) -> list:
    if value is None:
        return []
    if isinstance(value, str):
        return [item for item in value.split(",") if item]
    return list(value)


def has_identity_qualified_viral_hit(
    prefixes,
    pidents,
    viral_prefixes: set[str] | None = None,
    min_pident: float = MIN_VIRAL_HIT_PIDENT,
) -> bool:
    """Return True for any viral top-10 hit with percent identity >= threshold."""
    prefix_values = _listify_top10_field(prefixes)
    pident_values = _listify_top10_field(pidents)
    allowed = viral_prefixes or VIRAL_PREFIXES
    allowed = {p if str(p).endswith("__") else f"{p}__" for p in allowed}

    for idx, raw_prefix in enumerate(prefix_values):
        prefix = str(raw_prefix).strip()
        if prefix and not prefix.endswith("__"):
            prefix = f"{prefix}__"
        if prefix not in allowed:
            continue
        try:
            pident = float(pident_values[idx])
        except (IndexError, TypeError, ValueError):
            continue
        if pident >= min_pident:
            return True
    return False


@dataclass
class BoundaryDiamondConfig:
    """Configuration for batched boundary Diamond analysis.

    Attributes:
        flank_genes: Number of genes on each side of EVE boundary
            (for boundary context)
        control_sample_size: Total control genes to sample (genome-wide, not per-EVE)
        control_min_distance: Minimum genes away from any EVE boundary for controls
        control_region_genes: Number of consecutive genes per control sample
        top_k: Top-N hits to retrieve from Diamond
        chunk_size: Maximum proteins per Diamond chunk for memory management
        threads: Number of threads for Diamond execution
        random_seed: Seed for deterministic control sampling
        host_prefix: Expected host taxonomy prefix (EUK__ or ARC__ for archaeal hosts)
    """

    flank_genes: int = 10
    control_sample_size: int = 100
    control_min_distance: int = 30
    control_region_genes: int = 11
    top_k: int = 10
    chunk_size: int = 10000
    threads: int = 8
    random_seed: int = 42
    host_prefix: str = "EUK__"
    taxonomy_weight_mode: str = "rank"
    search_backend: str = "diamond"


@dataclass
class SeedGeneMapping:
    """Mapping of seed region to its gene IDs for boundary enforcement.

    Uses ordered lists for deterministic output.
    Includes stable seed_id for boundary-to-seed mapping.

    Attributes:
        seed_id: Stable ID from MergedSeed (format: "seed_{idx}_{scaffold}_{start}")
        scaffold: Scaffold/contig name
        seed_start: Original seed start (pre-host-trim)
        seed_end: Original seed end (pre-host-trim)
        eve_porf_ids: Genes overlapping seed (ordered by position)
        upstream_porf_ids: Upstream flanking genes (ordered, closest to EVE first)
        downstream_porf_ids: Downstream flanking genes (ordered, closest to EVE first)
        flank_start_idx: First gene index in proteome (for boundary constraint)
        flank_end_idx: Last gene index in proteome (for boundary constraint)
        flank_start_bp: Genomic position of first flanking gene
        flank_end_bp: Genomic position of last flanking gene
        flank_genes_config: The flank_genes value used (for audit)
    """

    seed_id: str
    scaffold: str
    seed_start: int
    seed_end: int
    eve_porf_ids: list[str] = field(default_factory=list)
    upstream_porf_ids: list[str] = field(default_factory=list)
    downstream_porf_ids: list[str] = field(default_factory=list)
    flank_start_idx: int = 0
    flank_end_idx: int = 0
    flank_start_bp: int = 0
    flank_end_bp: int = 0
    flank_genes_config: int = 10


@dataclass
class GenomeDiamondQuery:
    """Aggregated query for all EVEs in a genome.

    Contains protein IDs organized by role (EVE interior, boundary, control).
    The all_porf_ids field is the union of all other sets.
    The seed_gene_mappings provides stable mapping for boundary enforcement.

    Attributes:
        eve_porf_ids: Mapping of seed_id to pORF IDs within EVE boundaries
        boundary_porf_ids: Mapping of seed_id to pORF IDs in boundary regions
        control_porf_ids: Control pORF IDs (shared across all EVEs)
        all_porf_ids: Union of all above for Diamond query
        seed_gene_mappings: Mapping of seed_id to SeedGeneMapping for boundary enforcement
    """

    eve_porf_ids: dict[str, list[str]] = field(default_factory=dict)
    boundary_porf_ids: dict[str, list[str]] = field(default_factory=dict)
    control_porf_ids: list[str] = field(default_factory=list)
    all_porf_ids: list[str] = field(default_factory=list)
    seed_gene_mappings: dict[str, SeedGeneMapping] = field(default_factory=dict)


@dataclass
class GeneTaxonomy:
    """Taxonomy classification for a single pORF.

    Contains both the best hit information and aggregate statistics
    across top-10 hits.

    Attributes:
        porf_id: Protein ORF identifier
        scaffold: Scaffold/contig name
        start: Start position (bp)
        end: End position (bp)
        top1_target: Best hit target ID
        top1_prefix: Taxonomy prefix of best hit
        top1_pident: Percent identity of best hit
        top10_prefixes: List of all top-10 prefixes
        taxonomy_fingerprint: Aggregated taxonomy fingerprint from top-10 hits
        has_ncldv_mirus: Whether NCLDV or MIRUS appears in top-10
        has_vp_plv: Whether VP or PLV appears in top-10
        has_viral: Whether any viral prefix appears in top-10
        has_hit: Whether Diamond found any hit (False = no-hit entry)
    """

    porf_id: str
    scaffold: str
    start: int
    end: int
    top1_target: str = ""
    top1_prefix: str = "UNKNOWN"
    top1_pident: float = 0.0
    top1_evalue: float = 1.0
    top10_prefixes: list[str] = field(default_factory=list)
    top10_targets: list[str] = field(default_factory=list)
    top10_bits: list[float] = field(default_factory=list)
    top10_pidents: list[float] = field(default_factory=list)
    top10_evalues: list[float] = field(default_factory=list)
    taxonomy_fingerprint: Optional[TaxonomyFingerprint] = None
    has_ncldv_mirus: bool = False
    has_vp_plv: bool = False
    has_viral: bool = False
    has_hit: bool = False


@dataclass
class ControlStats:
    """Statistics from control (host) region.

    Used to establish baseline host signal for comparison with EVE regions.

    Attributes:
        n_genes: Total number of control genes
        n_no_hits: Number of genes without Diamond hits
        no_hit_frequency: Proportion of genes with no hits
        host_frequency: Proportion with host prefix as best hit
        mean_pident: Mean percent identity for host hits
        dominant_organism: Most common organism in host hits
        host_prefix: Which prefix was used (EUK__ or ARC__)
    """

    n_genes: int = 0
    n_no_hits: int = 0
    no_hit_frequency: float = 0.0
    host_frequency: float = 0.0
    mean_pident: float = 0.0
    dominant_organism: str = "unknown"
    host_prefix: str = "EUK__"


def build_gene_taxonomy_record(
    tax: GeneTaxonomy,
    is_flanking: bool = False,
    flank_position: Optional[str] = None,
) -> dict[str, Any]:
    """
    Build standardized gene taxonomy record for evidence synthesis.

    Converts GeneTaxonomy object to dict format expected by verify_eve_task.
    Handles prefix stripping (EUK__ → EUK) for output compatibility.

    Args:
        tax: GeneTaxonomy object from boundary_diamond
        is_flanking: Whether gene is flanking region (vs. EVE interior)
        flank_position: "upstream" or "downstream" if is_flanking=True

    Returns:
        Dictionary with standardized taxonomy fields
    """

    # Strip trailing underscores from prefixes
    top1_prefix_stripped = (
        tax.top1_prefix.rstrip("_") if tax.top1_prefix else "UNKNOWN"
    )
    top10_prefixes_stripped = (
        [p.rstrip("_") for p in tax.top10_prefixes]
        if tax.top10_prefixes
        else []
    )

    # Serialize taxonomy_fingerprint so host signature scoring can use it
    fp_dict = None
    if tax.taxonomy_fingerprint is not None:
        fp_dict = {
            "weighted_tokens": tax.taxonomy_fingerprint.weighted_tokens,
            "raw_tokens": tax.taxonomy_fingerprint.raw_tokens,
        }

    return {
        "porf_id": tax.porf_id,
        "scaffold": tax.scaffold,
        "start": tax.start,
        "end": tax.end,
        "top1_target": tax.top1_target,
        "top1_prefix": top1_prefix_stripped,
        "top1_pident": tax.top1_pident,
        "top1_evalue": tax.top1_evalue,
        "top10_prefixes": top10_prefixes_stripped,
        "top10_targets": tax.top10_targets,
        "top10_bitscores": tax.top10_bits,
        "top10_pidents": tax.top10_pidents,
        "top10_evalues": tax.top10_evalues,
        "taxonomy_fingerprint": fp_dict,
        "has_ncldv_mirus": tax.has_ncldv_mirus,
        "has_vp_plv": tax.has_vp_plv,
        "has_viral": tax.has_viral,
        "has_hit": tax.has_hit,
        "is_flanking": is_flanking,
        "flank_position": flank_position,
    }


@dataclass
class DiamondHit:
    """Single Diamond hit result."""

    query: str
    target: str
    evalue: float
    bits: float
    pident: float
    qcov: float


@dataclass
class pORF:
    """Simple pORF representation for indexing."""

    id: str
    scaffold: str
    start: int
    end: int
    strand: str = "+"


def extract_prefix(target_id: str) -> str:
    """
    Extract taxonomy prefix from target ID.

    Returns prefix WITH trailing underscores for consistency.

    Args:
        target_id: Target sequence ID (e.g., "NCLDV__GVOGm0003_...")

    Returns:
        Taxonomy prefix (e.g., "NCLDV__") or "UNKNOWN" if not found

    Examples:
        >>> extract_prefix("NCLDV__GVOGm0003_Marseillevirus")
        'NCLDV__'
        >>> extract_prefix("EUK__Arabidopsis_thaliana|protein123")
        'EUK__'
        >>> extract_prefix("unknown_protein")
        'UNKNOWN'
    """
    for prefix in VIRAL_PREFIXES | CELLULAR_PREFIXES:
        if target_id.startswith(prefix):
            return prefix
    return "UNKNOWN"


def extract_organism(target_id: str) -> str:
    """
    Extract organism name from target ID.

    Args:
        target_id: Target sequence ID (e.g., "EUK__Arabidopsis_thaliana|protein123")

    Returns:
        Organism name or "unknown"
    """
    if "__" not in target_id:
        return "unknown"
    # Split on prefix, then extract organism before | separator
    parts = target_id.split("__", 1)
    if len(parts) < 2:
        return "unknown"
    organism_part = parts[1].split("|")[0]
    return organism_part if organism_part else "unknown"


def build_proteome_index(proteome_fasta: Path) -> dict[str, list[pORF]]:
    """
    Build an index of pORFs organized by scaffold.

    Args:
        proteome_fasta: Path to proteome FASTA file

    Returns:
        Dict mapping scaffold name to sorted list of pORF objects
    """
    index: dict[str, list[pORF]] = {}

    for record in SeqIO.parse(proteome_fasta, "fasta"):
        parsed = parse_prodigal_header(record.description, record.id)
        if not parsed:
            continue
        scaffold, start, end, strand = parsed
        porf = pORF(id=record.id, scaffold=scaffold, start=start, end=end, strand=strand)

        if scaffold not in index:
            index[scaffold] = []
        index[scaffold].append(porf)

    # Sort by start position within each scaffold
    for scaffold in index:
        index[scaffold].sort(key=lambda p: p.start)

    return index


def collect_query_proteins(
    merged_seeds: list,
    proteome_index: dict[str, list[pORF]],
    config: BoundaryDiamondConfig,
) -> GenomeDiamondQuery:
    """
    Collect all proteins needed for Diamond from ALL seeds.

    This is called ONCE per genome to build a single query file.

    Args:
        merged_seeds: List of MergedSeed objects from Phase 1
        proteome_index: Scaffold -> sorted list of pORF objects
        config: BoundaryDiamondConfig with sampling parameters

    Returns:
        GenomeDiamondQuery containing all pORF IDs organized by role
    """
    active_logger = _get_logger()

    # Runtime introspection and verification
    active_logger.info("=" * 80)
    active_logger.info("COLLECT_QUERY_PROTEINS START")
    active_logger.info("Module: %s", __file__)
    active_logger.info(
        "Config: flank_genes=%d, control_sample_size=%d, control_min_distance=%d, control_region_genes=%d",
        config.flank_genes,
        config.control_sample_size,
        config.control_min_distance,
        config.control_region_genes,
    )
    active_logger.info("Input: %d merged seeds", len(merged_seeds))

    # Log git SHA if available (for verifying which code is running)
    try:
        git_sha = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
            cwd=Path(__file__).parent
        ).decode().strip()
        active_logger.info("Git SHA: %s", git_sha)
    except Exception as e:
        active_logger.debug("Git SHA unavailable: %s", e)
        active_logger.info("Git SHA: unavailable (not in git repo or git not installed)")

    rng = random.Random(config.random_seed)

    eve_porf_ids: dict[str, list[str]] = {}
    boundary_porf_ids: dict[str, list[str]] = {}
    all_eve_boundary_indices: dict[str, set[int]] = {}

    seed_gene_mappings: dict[str, SeedGeneMapping] = {}

    # Pass 1: Collect EVE and boundary pORFs for all seeds
    for seed in merged_seeds:
        scaffold = seed.scaffold
        scaffold_porfs = proteome_index.get(scaffold, [])
        if not scaffold_porfs:
            continue

        # Find pORFs within EVE
        eve_porfs = []
        eve_indices = []
        for i, p in enumerate(scaffold_porfs):
            # Consider a pORF "within" if it overlaps with the seed region
            if p.start < seed.end and p.end > seed.start:
                eve_porfs.append(p)
                eve_indices.append(i)

        if not eve_indices:
            continue

        # Use stable seed_id from MergedSeed (assigned in seed_merger.py)
        seed_id = getattr(seed, "seed_id", "") or f"EVE_{scaffold}_{seed.start}-{seed.end}"

        # Store EVE pORF IDs (ordered by position)
        eve_porf_ids[seed_id] = [p.id for p in eve_porfs]

        # Find boundary pORFs (+/-flank_genes)
        first_idx = min(eve_indices)
        last_idx = max(eve_indices)
        boundary_start_idx = max(0, first_idx - config.flank_genes)
        boundary_end_idx = min(len(scaffold_porfs) - 1, last_idx + config.flank_genes)

        # Split flanking genes into upstream and downstream (ordered, closest first)
        upstream_porfs = scaffold_porfs[boundary_start_idx:first_idx]
        downstream_porfs = scaffold_porfs[last_idx + 1 : boundary_end_idx + 1]

        # Reverse upstream so closest to EVE is first
        upstream_porf_ids_ordered = [p.id for p in reversed(upstream_porfs)]
        downstream_porf_ids_ordered = [p.id for p in downstream_porfs]

        # Boundary includes genes outside the EVE but within flank range
        boundary_porfs = upstream_porfs + downstream_porfs
        boundary_porf_ids[seed_id] = [p.id for p in boundary_porfs]

        # Log per-EVE region collection
        active_logger.info("EVE %s: collected %d EVE genes, boundary: upstream=%d, downstream=%d",
                          seed_id, len(eve_porfs), len(upstream_porfs), len(downstream_porfs))

        # Calculate flanking bounds for boundary constraint enforcement
        flank_start_bp = scaffold_porfs[boundary_start_idx].start if boundary_start_idx < len(scaffold_porfs) else seed.start
        flank_end_bp = scaffold_porfs[boundary_end_idx].end if boundary_end_idx < len(scaffold_porfs) else seed.end

        # Create SeedGeneMapping for boundary enforcement
        seed_gene_mappings[seed_id] = SeedGeneMapping(
            seed_id=seed_id,
            scaffold=scaffold,
            seed_start=seed.start,
            seed_end=seed.end,
            eve_porf_ids=[p.id for p in eve_porfs],  # Ordered by position
            upstream_porf_ids=upstream_porf_ids_ordered,
            downstream_porf_ids=downstream_porf_ids_ordered,
            flank_start_idx=boundary_start_idx,
            flank_end_idx=boundary_end_idx,
            flank_start_bp=flank_start_bp,
            flank_end_bp=flank_end_bp,
            flank_genes_config=config.flank_genes,
        )

        # Track excluded indices for control sampling
        if scaffold not in all_eve_boundary_indices:
            all_eve_boundary_indices[scaffold] = set()
        all_eve_boundary_indices[scaffold].update(
            range(boundary_start_idx, boundary_end_idx + 1)
        )

    # Pass 2: Sample control pORFs from outside ALL EVE/boundary regions
    control_porf_ids = sample_control_porfs_genome_wide(
        proteome_index=proteome_index,
        excluded_indices=all_eve_boundary_indices,
        sample_size=config.control_sample_size,
        min_distance=config.control_min_distance,
        rng=rng,
        stretch_size=config.control_region_genes,
    )

    # Track deduplication statistics
    total_eve_genes_before_dedup = sum(len(ids) for ids in eve_porf_ids.values())
    total_boundary_genes_before_dedup = sum(len(ids) for ids in boundary_porf_ids.values())
    total_before_dedup = total_eve_genes_before_dedup + total_boundary_genes_before_dedup + len(control_porf_ids)

    # Single source of truth: compute expected counts from actual collected genes
    active_logger.info("Collected %d EVE regions", len(eve_porf_ids))
    active_logger.info("Before dedup - EVE: %d, Boundary: %d, Controls: %d, Total: %d",
                      total_eve_genes_before_dedup, total_boundary_genes_before_dedup,
                      len(control_porf_ids), total_before_dedup)
    active_logger.info("EXPECTED (computed): EVE: %d, Boundary: %d, Controls: %d, Total: %d",
                      total_eve_genes_before_dedup, total_boundary_genes_before_dedup,
                      len(control_porf_ids), total_before_dedup)

    # Build union of all pORF IDs (ordered, deduplicated)
    seen = set()
    all_porf_ids = []

    # Add EVE pORFs first (ordered by seed, then by position within seed)
    eve_duplicates = 0
    for ids in eve_porf_ids.values():
        for pid in ids:
            if pid not in seen:
                seen.add(pid)
                all_porf_ids.append(pid)
            else:
                eve_duplicates += 1

    # Add boundary pORFs second
    boundary_duplicates = 0
    for ids in boundary_porf_ids.values():
        for pid in ids:
            if pid not in seen:
                seen.add(pid)
                all_porf_ids.append(pid)
            else:
                boundary_duplicates += 1

    # Add control pORFs last
    control_duplicates = 0
    for pid in control_porf_ids:
        if pid not in seen:
            seen.add(pid)
            all_porf_ids.append(pid)
        else:
            control_duplicates += 1

    total_duplicates = eve_duplicates + boundary_duplicates + control_duplicates
    dedup_pct = (total_duplicates / total_before_dedup * 100) if total_before_dedup > 0 else 0

    logger.info(
        "Query protein deduplication: %d → %d (removed %d duplicates, %.1f%%)",
        total_before_dedup,
        len(all_porf_ids),
        total_duplicates,
        dedup_pct,
    )
    logger.info(
        "  EVE duplicates: %d/%d, Boundary duplicates: %d/%d, Control duplicates: %d/%d",
        eve_duplicates,
        total_eve_genes_before_dedup,
        boundary_duplicates,
        total_boundary_genes_before_dedup,
        control_duplicates,
        len(control_porf_ids),
    )

    logger.info(
        "Collected query proteins: %d EVE regions, %d EVE pORFs (deduplicated), "
        "%d boundary pORFs (deduplicated), %d controls, %d total",
        len(eve_porf_ids),
        total_eve_genes_before_dedup - eve_duplicates,
        total_boundary_genes_before_dedup - boundary_duplicates,
        len(control_porf_ids) - control_duplicates,
        len(all_porf_ids),
    )

    # Final logging before return
    active_logger.info("Deduplication: %d duplicates found (EVE: %d, Boundary: %d, Control: %d)",
                      total_duplicates, eve_duplicates, boundary_duplicates, control_duplicates)
    active_logger.info("COLLECT_QUERY_PROTEINS RETURN - %d total genes (after dedup)", len(all_porf_ids))
    active_logger.info("=" * 80)

    return GenomeDiamondQuery(
        eve_porf_ids=eve_porf_ids,
        boundary_porf_ids=boundary_porf_ids,
        control_porf_ids=control_porf_ids,
        all_porf_ids=all_porf_ids,
        seed_gene_mappings=seed_gene_mappings,
    )


def sample_control_porfs_genome_wide(
    proteome_index: dict[str, list[pORF]],
    excluded_indices: dict[str, set[int]],
    sample_size: int,
    min_distance: int,
    rng: random.Random,
    stretch_size: int = 11,
) -> list[str]:
    """
    Sample control pORFs from regions distant from ALL EVEs.

    Sampling is:
    - Deterministic (seeded RNG)
    - Genome-wide (across all scaffolds)
    - Excludes all candidate viral regions
    - Samples in stretches of consecutive genes (default 11) to preserve gene neighborhood context

    Args:
        proteome_index: Scaffold -> sorted list of pORF objects
        excluded_indices: Scaffold -> set of excluded pORF indices
        sample_size: Target number of control genes
        min_distance: Minimum genes away from any EVE boundary
        rng: Seeded random number generator
        stretch_size: Number of consecutive genes per control sample (default 11)

    Returns:
        Deterministically ordered control pORF IDs
    """
    # Step 1: Build eligible gene stretches (outside exclusion zones)
    eligible_stretches = []  # List of (scaffold, start_idx, end_idx, porfs)

    for scaffold, porfs in proteome_index.items():
        excluded = excluded_indices.get(scaffold, set())

        # Expand exclusion zone by min_distance
        expanded_excluded = set()
        for idx in excluded:
            expanded_excluded.update(
                range(max(0, idx - min_distance), min(len(porfs), idx + min_distance + 1))
            )

        # Find contiguous eligible regions
        current_stretch_start = None
        for i in range(len(porfs)):
            if i not in expanded_excluded:
                if current_stretch_start is None:
                    current_stretch_start = i
            else:
                # End of stretch
                if current_stretch_start is not None:
                    stretch_length = i - current_stretch_start
                    if stretch_length >= stretch_size:
                        eligible_stretches.append(
                            (scaffold, current_stretch_start, i, porfs[current_stretch_start:i])
                        )
                    current_stretch_start = None

        # Handle stretch at end of scaffold
        if current_stretch_start is not None:
            stretch_length = len(porfs) - current_stretch_start
            if stretch_length >= stretch_size:
                eligible_stretches.append(
                    (scaffold, current_stretch_start, len(porfs), porfs[current_stretch_start:])
                )

    if not eligible_stretches:
        logger.warning("No eligible control stretches found (min_distance=%d, stretch_size=%d)", min_distance, stretch_size)
        return []

    # Step 2: Create sampling windows of stretch_size genes
    sampling_windows = []  # List of (scaffold, start_idx, genes_in_window)

    for scaffold, start_idx, end_idx, stretch_porfs in eligible_stretches:
        # Create non-overlapping windows (chunks) of stretch_size genes
        # Use step=stretch_size to avoid overlap, which would cause the sampling loop
        # to collect many more genes than intended when accumulating to target size
        for window_start in range(0, len(stretch_porfs) - stretch_size + 1, stretch_size):
            window_genes = stretch_porfs[window_start : window_start + stretch_size]
            sampling_windows.append((scaffold, start_idx + window_start, window_genes))

    if not sampling_windows:
        logger.warning(
            "No sampling windows found (eligible stretches: %d, stretch_size: %d)",
            len(eligible_stretches),
            stretch_size,
        )
        return []

    # Step 3: Randomly sample windows until reaching target gene count
    # Calculate the actual maximum unique genes available (not window count × stretch_size,
    # since windows overlap). Sum the length of each eligible stretch.
    max_possible_genes = sum(len(stretch_porfs) for _, _, _, stretch_porfs in eligible_stretches)
    if max_possible_genes < sample_size:
        logger.info(
            "Control sampling: max possible %d unique genes (requested %d), using all %d eligible stretches",
            max_possible_genes,
            sample_size,
            len(eligible_stretches),
        )
        # Use all genes from all eligible stretches
        selected_genes = set()
        for _, _, _, stretch_porfs in eligible_stretches:
            selected_genes.update(p.id for p in stretch_porfs)
        return sorted(selected_genes)

    # Sample windows randomly
    selected_genes = set()
    sampled_windows = []
    rng.shuffle(sampling_windows)  # Shuffle to randomize selection

    for window in sampling_windows:
        if len(selected_genes) >= sample_size:
            break
        scaffold, start_idx, window_genes = window
        # Add all genes from this window
        selected_genes.update(p.id for p in window_genes)
        sampled_windows.append((scaffold, start_idx, len(window_genes)))

    logger.info(
        "Control sampling: selected %d windows (%d genes) from %d available windows",
        len(sampled_windows),
        len(selected_genes),
        len(sampling_windows),
    )

    return sorted(selected_genes)


def extract_sequences(
    proteome_fasta: Path,
    porf_ids: set[str],
    output_fasta: Path,
) -> int:
    """
    Extract sequences for specified pORF IDs to output FASTA.

    Args:
        proteome_fasta: Path to full proteome FASTA
        porf_ids: Set of pORF IDs to extract
        output_fasta: Path to write extracted sequences

    Returns:
        Number of sequences extracted
    """
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    extracted = 0

    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(proteome_fasta, "fasta"):
            if record.id in porf_ids:
                out_handle.write(f">{record.id}\n{record.seq}\n")
                extracted += 1

    return extracted


def run_diamond_blastp(
    query: Path,
    db: Path,
    output: Path,
    max_target_seqs: int = 10,
    threads: int = 8,
    evalue: float = 1e-5,
    search_backend: str = "diamond",
) -> None:
    """
    Run protein sequence search.

    Args:
        query: Query FASTA file
        db: Diamond database path (.dmnd)
        output: Output TSV file
        max_target_seqs: Maximum hits per query
        threads: Number of threads
        evalue: E-value cutoff
        search_backend: "diamond" (the only supported backend)
    """
    from virosync.pipeline.search_backend import run_sequence_search

    run_sequence_search(
        query_fasta=query,
        db_path=db,
        output_tsv=output,
        threads=threads,
        backend=search_backend,
        evalue=evalue,
        max_target_seqs=max_target_seqs,
    )


def parse_diamond_output(output_file: Path) -> dict[str, list[DiamondHit]]:
    """
    Parse Diamond output TSV into hits by query.

    Args:
        output_file: Diamond output TSV

    Returns:
        Dict mapping query ID to list of DiamondHit (sorted by bits, descending)
    """
    hits_by_query: dict[str, list[DiamondHit]] = {}

    if not output_file.exists():
        logger.warning("Diamond output not found: %s", output_file)
        return hits_by_query

    with open(output_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue

            query, target, evalue, bits, pident, qcov = parts[:6]
            try:
                hit = DiamondHit(
                    query=query,
                    target=target,
                    evalue=float(evalue),
                    bits=float(bits),
                    pident=float(pident),
                    qcov=float(qcov),
                )
                hits_by_query.setdefault(query, []).append(hit)
            except ValueError:
                continue

    # Sort by bits score (descending) for each query
    for query in hits_by_query:
        hits_by_query[query].sort(key=lambda h: h.bits, reverse=True)

    return hits_by_query


def run_full_proteome_diamond(
    proteome_fasta: Path,
    diamond_db: Path,
    output_dir: Path,
    *,
    max_target_seqs: int,
    threads: int,
    search_backend: str = "diamond",
) -> dict[str, list[DiamondHit]]:
    """Search the full proteome in one unchunked DIAMOND invocation.

    This is the opt-in Phase 2 superset prototype seam. Query IDs remain the
    raw pORF IDs so one result can be sliced for both host trimming and boundary
    refinement.
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "full_proteome.tsv"
    run_diamond_blastp(
        query=proteome_fasta,
        db=diamond_db,
        output=output_file,
        max_target_seqs=max_target_seqs,
        threads=threads,
        search_backend=search_backend,
    )
    return parse_diamond_output(output_file)


def classify_cached_diamond_query(
    query: GenomeDiamondQuery,
    diamond_hits: dict[str, list[DiamondHit]],
    proteome_index: dict[str, list[pORF]],
    config: BoundaryDiamondConfig,
    taxonomy_lookup: Optional[dict] = None,
) -> dict[str, GeneTaxonomy]:
    """Classify an exact boundary query by slicing full-proteome raw hits."""

    selected_hits = {
        porf_id: diamond_hits.get(porf_id, [])[: config.top_k]
        for porf_id in query.all_porf_ids
    }
    return classify_all_porfs(
        all_porf_ids=query.all_porf_ids,
        diamond_hits=selected_hits,
        proteome_index=proteome_index,
        host_prefix=config.host_prefix,
        taxonomy_lookup=taxonomy_lookup,
        taxonomy_weight_mode=config.taxonomy_weight_mode,
    )


def split_fasta(
    input_fasta: Path,
    chunk_size: int,
    output_dir: Path,
) -> list[Path]:
    """
    Split FASTA file into chunks of specified size.

    Args:
        input_fasta: Input FASTA file
        chunk_size: Maximum sequences per chunk
        output_dir: Directory to write chunk files

    Returns:
        List of paths to chunk FASTA files
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    chunks: list[Path] = []
    current_chunk: list = []
    chunk_num = 0

    for record in SeqIO.parse(input_fasta, "fasta"):
        current_chunk.append(record)

        if len(current_chunk) >= chunk_size:
            chunk_num += 1
            chunk_path = output_dir / f"chunk_{chunk_num}.faa"
            with open(chunk_path, "w") as out_handle:
                for rec in current_chunk:
                    out_handle.write(f">{rec.id}\n{rec.seq}\n")
            chunks.append(chunk_path)
            current_chunk = []

    # Write remaining sequences
    if current_chunk:
        chunk_num += 1
        chunk_path = output_dir / f"chunk_{chunk_num}.faa"
        with open(chunk_path, "w") as out_handle:
            for rec in current_chunk:
                out_handle.write(f">{rec.id}\n{rec.seq}\n")
        chunks.append(chunk_path)

    return chunks


def run_diamond_chunked(
    query_fasta: Path,
    diamond_db: Path,
    output_dir: Path,
    chunk_size: int,
    max_target_seqs: int,
    threads: int,
    search_backend: str = "diamond",
) -> dict[str, list[DiamondHit]]:
    """
    Run Diamond in chunks for large query sets.

    Chunks are processed sequentially to avoid memory issues.
    Results are merged.

    Args:
        query_fasta: Path to query FASTA
        diamond_db: Path to Diamond database
        output_dir: Directory for output files
        chunk_size: Maximum sequences per chunk
        max_target_seqs: Maximum hits per query
        threads: Number of threads

    Returns:
        Dict mapping query ID to list of DiamondHit
    """
    chunks_dir = output_dir / "chunks"
    chunks = split_fasta(query_fasta, chunk_size, chunks_dir)

    logger.info("Chunked %s into %d chunks of max %d sequences", query_fasta, len(chunks), chunk_size)

    all_hits: dict[str, list[DiamondHit]] = {}

    for i, chunk_fasta in enumerate(chunks):
        chunk_output = output_dir / f"diamond_chunk_{i + 1}.tsv"
        logger.info("Processing chunk %d/%d", i + 1, len(chunks))

        run_diamond_blastp(
            query=chunk_fasta,
            db=diamond_db,
            output=chunk_output,
            max_target_seqs=max_target_seqs,
            threads=threads,
            search_backend=search_backend,
        )

        chunk_hits = parse_diamond_output(chunk_output)
        all_hits.update(chunk_hits)

    return all_hits


def classify_all_porfs(
    all_porf_ids: set[str],
    diamond_hits: dict[str, list[DiamondHit]],
    proteome_index: dict[str, list[pORF]],
    host_prefix: str,
    taxonomy_lookup: Optional[dict] = None,
    taxonomy_weight_mode: str = "rank",
) -> dict[str, GeneTaxonomy]:
    """
    Classify all pORFs, including those with no Diamond hits.

    IMPORTANT: Explicitly materializes "no hit" entries so downstream
    code can distinguish between "no hit" and "missing data".

    Args:
        all_porf_ids: Set of all pORF IDs to classify
        diamond_hits: Diamond results keyed by pORF ID
        proteome_index: Scaffold -> sorted list of pORF objects
        host_prefix: Host taxonomy prefix (e.g., "EUK__")
        taxonomy_lookup: Optional taxonomy label lookup for fingerprinting

    Returns:
        Dict mapping pORF ID to GeneTaxonomy
    """
    # Build pORF lookup for coordinates
    porf_lookup: dict[str, pORF] = {}
    for porfs in proteome_index.values():
        for p in porfs:
            if p.id in all_porf_ids:
                porf_lookup[p.id] = p

    taxonomy_map: dict[str, GeneTaxonomy] = {}

    for porf_id in all_porf_ids:
        porf = porf_lookup.get(porf_id)
        if not porf:
            # pORF not found in index - create minimal entry
            taxonomy_map[porf_id] = GeneTaxonomy(
                porf_id=porf_id,
                scaffold="unknown",
                start=0,
                end=0,
                top1_prefix="UNKNOWN",
                has_hit=False,
            )
            continue

        hits = diamond_hits.get(porf_id, [])

        if not hits:
            # Explicitly materialize no-hit entry
            taxonomy_map[porf_id] = GeneTaxonomy(
                porf_id=porf_id,
                scaffold=porf.scaffold,
                start=porf.start,
                end=porf.end,
                top1_target="",
                top1_prefix="UNKNOWN",
                top1_pident=0.0,
                top1_evalue=1.0,
                top10_prefixes=[],
                top10_targets=[],
                top10_bits=[],
                top10_pidents=[],
                top10_evalues=[],
                has_ncldv_mirus=False,
                has_vp_plv=False,
                has_viral=False,
                has_hit=False,
            )
        else:
            # Classify based on hits (take top-10)
            top10 = hits[:10]
            prefixes = [extract_prefix(h.target) for h in top10]
            top10_targets = [h.target for h in top10]
            top10_bits = [h.bits for h in top10]
            top10_pidents = [h.pident for h in top10]
            top10_evalues = [h.evalue for h in top10]
            top1 = top10[0]

            has_ncldv_mirus = any(p in {"NCLDV__", "MIRUS__"} for p in prefixes)
            has_vp_plv = any(p in {"VP__", "PLV__", "PPV__"} for p in prefixes)
            has_viral = has_identity_qualified_viral_hit(prefixes, top10_pidents)

            # Calculate taxonomy fingerprint if lookup available
            fingerprint = None
            if taxonomy_lookup:
                top10_tuples = [
                    (h.target, h.bits, h.pident, h.evalue) for h in top10
                ]
                fingerprint = aggregate_taxonomy_substrings(
                    top10_tuples,
                    taxonomy_lookup,
                    min_token_length=3,
                    weight_mode=taxonomy_weight_mode,
                )

            taxonomy_map[porf_id] = GeneTaxonomy(
                porf_id=porf_id,
                scaffold=porf.scaffold,
                start=porf.start,
                end=porf.end,
                top1_target=top1.target,
                top1_prefix=extract_prefix(top1.target),
                top1_pident=top1.pident,
                top1_evalue=top1.evalue,
                top10_prefixes=prefixes,
                top10_targets=top10_targets,
                top10_bits=top10_bits,
                top10_pidents=top10_pidents,
                top10_evalues=top10_evalues,
                taxonomy_fingerprint=fingerprint,
                has_ncldv_mirus=has_ncldv_mirus,
                has_vp_plv=has_vp_plv,
                has_viral=has_viral,
                has_hit=True,
            )

    return taxonomy_map


def build_host_baseline_fingerprint(
    control_taxonomy: list[GeneTaxonomy],
    host_prefix: str = "EUK__",
    min_token_count: int = 3,
    min_weight_fraction: float = 0.10,
) -> dict[str, float]:
    """
    Build host baseline fingerprint from control genes.

    Filters to host-prefix genes only to avoid viral contamination,
    then aggregates tokens that appear frequently enough.

    Args:
        control_taxonomy: List of GeneTaxonomy for control genes
        host_prefix: Filter to genes with this prefix (e.g., "EUK__")
        min_token_count: Minimum genes token must appear in
        min_weight_fraction: Minimum average weight

    Returns:
        Dict of {token: normalized_weight}
    """
    if not control_taxonomy:
        return {}

    # Filter to host-prefix genes only (avoid viral contamination)
    host_genes = [
        g
        for g in control_taxonomy
        if g.top1_prefix == host_prefix and g.taxonomy_fingerprint
    ]

    if len(host_genes) < min_token_count:
        # Not enough control genes - return empty (will use fallback logic)
        return {}

    token_weights = {}  # {token: [weights]}

    for gene_tax in host_genes:
        if not gene_tax.taxonomy_fingerprint:
            continue

        for token, weight in gene_tax.taxonomy_fingerprint.weighted_tokens.items():
            if token not in token_weights:
                token_weights[token] = []
            token_weights[token].append(weight)

    # Calculate baseline
    baseline = {}
    for token, weights in token_weights.items():
        if len(weights) < min_token_count:
            continue

        frequency = len(weights) / len(host_genes)
        avg_weight = sum(weights) / len(weights)
        combined = frequency * avg_weight

        if combined >= min_weight_fraction:
            baseline[token] = combined

    return baseline


def run_batched_diamond(
    query: GenomeDiamondQuery,
    proteome_fasta: Path,
    diamond_db: Path,
    output_dir: Path,
    proteome_index: dict[str, list[pORF]],
    config: BoundaryDiamondConfig,
    taxonomy_lookup: Optional[dict] = None,
) -> dict[str, GeneTaxonomy]:
    """
    Run Diamond ONCE for all proteins in the genome query.

    Returns taxonomy keyed by pORF ID for fast lookup.

    Args:
        query: GenomeDiamondQuery with all pORF IDs to query
        proteome_fasta: Path to proteome FASTA
        diamond_db: Path to Diamond database
        output_dir: Directory for output files
        proteome_index: Scaffold -> sorted list of pORF objects
        config: BoundaryDiamondConfig with execution parameters
        taxonomy_lookup: Optional taxonomy label lookup for fingerprinting

    Returns:
        Dict mapping pORF ID to GeneTaxonomy
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    active_logger = _get_logger()

    # WRITE PATH VERIFICATION: Log before writing
    query_fasta = output_dir / "genome_boundary_query.faa"
    active_logger.info("WRITE PATH: About to write %d genes to %s", len(query.all_porf_ids), query_fasta)
    active_logger.info("WRITE PATH: EVE regions=%d, Control genes=%d",
                      len(query.eve_porf_ids), len(query.control_porf_ids))

    # Write debug artifact to file for post-run inspection.
    debug_file = output_dir / "collect_query_proteins_debug.txt"
    with open(debug_file, "w") as f:
        f.write(f"Module: {__file__}\n")
        f.write(f"EVE regions: {len(query.eve_porf_ids)}\n")
        f.write(f"All genes to write: {len(query.all_porf_ids)}\n")
        f.write(f"Controls: {len(query.control_porf_ids)}\n")
        f.write(f"Output file: {query_fasta}\n")

    # Extract all sequences to single file
    n_extracted = extract_sequences(proteome_fasta, query.all_porf_ids, query_fasta)

    # WRITE PATH VERIFICATION: Verify after writing
    actual_written = sum(1 for line in query_fasta.read_text().split("\n") if line.startswith(">"))
    active_logger.info("WRITE PATH: Verified %d genes written to file", actual_written)

    # CRITICAL ASSERTION: Catch gene count mismatches immediately
    if actual_written != len(query.all_porf_ids):
        error_msg = f"Gene count mismatch: {actual_written} written vs {len(query.all_porf_ids)} expected"
        active_logger.error("WRITE PATH ERROR: %s", error_msg)
        raise AssertionError(error_msg)

    logger.info(
        "Extracted %d sequences for Diamond (requested %d)",
        n_extracted,
        len(query.all_porf_ids),
    )

    n_sequences = n_extracted
    safe_chunk_size = max(1, config.chunk_size)
    estimated_chunks = max(1, (n_sequences + safe_chunk_size - 1) // safe_chunk_size)
    active_logger.info(
        "Phase 2b: Diamond settings: threads=%d top_k=%d chunk_size=%d queries=%d estimated_chunks=%d",
        config.threads,
        config.top_k,
        config.chunk_size,
        n_sequences,
        estimated_chunks,
    )

    # Chunk if needed
    if n_sequences <= config.chunk_size:
        # Single run
        diamond_output = output_dir / "boundary_diamond.tsv"
        run_diamond_blastp(
            query=query_fasta,
            db=diamond_db,
            output=diamond_output,
            max_target_seqs=config.top_k,
            threads=config.threads,
            search_backend=config.search_backend,
        )
        hits = parse_diamond_output(diamond_output)
    else:
        # Chunked execution
        hits = run_diamond_chunked(
            query_fasta=query_fasta,
            diamond_db=diamond_db,
            output_dir=output_dir,
            chunk_size=config.chunk_size,
            max_target_seqs=config.top_k,
            threads=config.threads,
            search_backend=config.search_backend,
        )

    # Parse and classify all pORFs
    taxonomy_map = classify_all_porfs(
        all_porf_ids=query.all_porf_ids,
        diamond_hits=hits,
        proteome_index=proteome_index,
        host_prefix=config.host_prefix,
        taxonomy_lookup=taxonomy_lookup,
        taxonomy_weight_mode=config.taxonomy_weight_mode,
    )

    # Log summary statistics
    n_with_hits = sum(1 for t in taxonomy_map.values() if t.has_hit)
    n_viral = sum(1 for t in taxonomy_map.values() if t.has_viral)
    n_ncldv_mirus = sum(1 for t in taxonomy_map.values() if t.has_ncldv_mirus)

    logger.info(
        "Diamond classification: %d total, %d with hits, %d viral, %d NCLDV/MIRUS",
        len(taxonomy_map),
        n_with_hits,
        n_viral,
        n_ncldv_mirus,
    )

    return taxonomy_map


def filter_taxonomy_to_boundary(
    taxonomy_map: dict[str, GeneTaxonomy],
    refined_boundary,
) -> list[GeneTaxonomy]:
    """
    Filter pre-computed taxonomy to refined boundary.

    Called AFTER trimming to get taxonomy for final EVE region.
    Does NOT re-run Diamond.

    Args:
        taxonomy_map: Pre-computed taxonomy keyed by pORF ID
        refined_boundary: RefinedBoundary object with final coordinates

    Returns:
        List of GeneTaxonomy for pORFs within refined boundary
    """
    return [
        tax
        for tax in taxonomy_map.values()
        if (
            tax.scaffold == refined_boundary.scaffold
            and tax.start < refined_boundary.end
            and tax.end > refined_boundary.start
        )
    ]


def get_flanking_taxonomy(
    taxonomy_map: dict[str, GeneTaxonomy],
    proteome_index: dict[str, list[pORF]],
    refined_boundary,
    flank_genes: int = 10,
    seed_mapping: Optional[SeedGeneMapping] = None,
) -> tuple[list[GeneTaxonomy], list[GeneTaxonomy]]:
    """
    Get flanking gene taxonomy for genes outside but near an EVE boundary.

    Returns genes that are:
    - On the same scaffold as the boundary
    - Within flank_genes positions from the boundary
    - NOT overlapping with the EVE region

    Flanks are selected relative to the final refined boundary. The seed
    mapping is only a fallback when the final boundary cannot be placed in the
    proteome index. This matters when refinement contracts a seed: genes that
    were inside the seed then become final-boundary flanks.

    Args:
        taxonomy_map: Pre-computed taxonomy keyed by pORF ID
        proteome_index: Scaffold -> sorted list of pORF objects
        refined_boundary: RefinedBoundary object with final coordinates
        flank_genes: Number of genes on each side to include
        seed_mapping: Optional Phase 2b mapping used only as a fallback

    Returns:
        Tuple of (upstream_taxonomy, downstream_taxonomy) lists
    """
    scaffold = refined_boundary.scaffold
    scaffold_porfs = proteome_index.get(scaffold, [])

    # Prefer positions relative to the final boundary. Pre-computed seed flank
    # IDs describe the original seed and can omit genes exposed by contraction.
    eve_indices = [
        index
        for index, porf in enumerate(scaffold_porfs)
        if porf.start < refined_boundary.end and porf.end > refined_boundary.start
    ]
    if eve_indices:
        first_eve_idx = min(eve_indices)
        last_eve_idx = max(eve_indices)
        upstream_start_idx = max(0, first_eve_idx - flank_genes)
        downstream_end_idx = min(
            len(scaffold_porfs),
            last_eve_idx + 1 + flank_genes,
        )
        upstream_porfs = scaffold_porfs[upstream_start_idx:first_eve_idx]
        downstream_porfs = scaffold_porfs[last_eve_idx + 1:downstream_end_idx]

        upstream_taxonomy = [
            taxonomy_map[porf.id]
            for porf in upstream_porfs
            if porf.id in taxonomy_map
        ]
        downstream_taxonomy = [
            taxonomy_map[porf.id]
            for porf in downstream_porfs
            if porf.id in taxonomy_map
        ]
        missing_upstream = len(upstream_porfs) - len(upstream_taxonomy)
        missing_downstream = len(downstream_porfs) - len(downstream_taxonomy)
        if missing_upstream > 0 or missing_downstream > 0:
            logger.warning(
                "Flanking gene taxonomy incomplete for %s:%d-%d: "
                "%d/%d upstream missing, %d/%d downstream missing",
                scaffold,
                refined_boundary.start,
                refined_boundary.end,
                missing_upstream,
                len(upstream_porfs),
                missing_downstream,
                len(downstream_porfs),
            )
        return upstream_taxonomy, downstream_taxonomy

    # Fall back to the original seed mapping only when the final boundary
    # cannot be located in the proteome index.
    if seed_mapping is not None:
        upstream_taxonomy = []
        downstream_taxonomy = []
        missing_upstream = 0
        missing_downstream = 0
        filtered_upstream = 0
        filtered_downstream = 0

        # Use pre-computed upstream pORF IDs (ordered, closest to EVE first)
        # Filter out genes that now fall inside the refined boundary
        for porf_id in seed_mapping.upstream_porf_ids:
            tax = taxonomy_map.get(porf_id)
            if tax:
                # Skip genes that overlap with refined boundary (now inside EVE)
                if tax.start < refined_boundary.end and tax.end > refined_boundary.start:
                    filtered_upstream += 1
                    continue
                upstream_taxonomy.append(tax)
            else:
                missing_upstream += 1

        # Use pre-computed downstream pORF IDs (ordered, closest to EVE first)
        # Filter out genes that now fall inside the refined boundary
        for porf_id in seed_mapping.downstream_porf_ids:
            tax = taxonomy_map.get(porf_id)
            if tax:
                # Skip genes that overlap with refined boundary (now inside EVE)
                if tax.start < refined_boundary.end and tax.end > refined_boundary.start:
                    filtered_downstream += 1
                    continue
                downstream_taxonomy.append(tax)
            else:
                missing_downstream += 1

        # Log if any flanking genes are missing from taxonomy map
        if missing_upstream > 0 or missing_downstream > 0:
            logger.warning(
                "Flanking gene taxonomy incomplete for %s: "
                "%d/%d upstream missing, %d/%d downstream missing",
                seed_mapping.seed_id,
                missing_upstream,
                len(seed_mapping.upstream_porf_ids),
                missing_downstream,
                len(seed_mapping.downstream_porf_ids),
            )

        # Log if any flanking genes were filtered (now inside refined boundary)
        if filtered_upstream > 0 or filtered_downstream > 0:
            logger.debug(
                "Flanking genes filtered for %s (now inside refined boundary): "
                "%d upstream, %d downstream",
                seed_mapping.seed_id,
                filtered_upstream,
                filtered_downstream,
            )

        return upstream_taxonomy, downstream_taxonomy

    return [], []


def compute_control_stats(
    control_taxonomy: list[GeneTaxonomy],
    host_prefix: str,
) -> ControlStats:
    """
    Compute statistics from control region for comparison.

    Uses consistent prefix format (with __).
    Handles no-hit entries explicitly.

    Args:
        control_taxonomy: List of GeneTaxonomy from control pORFs
        host_prefix: Host taxonomy prefix (e.g., "EUK__")

    Returns:
        ControlStats with baseline host signal metrics
    """
    n_genes = len(control_taxonomy)
    if n_genes == 0:
        return ControlStats(
            n_genes=0,
            n_no_hits=0,
            no_hit_frequency=0.0,
            host_frequency=0.0,
            mean_pident=0.0,
            dominant_organism="unknown",
            host_prefix=host_prefix,
        )

    n_no_hits = sum(1 for g in control_taxonomy if not g.has_hit)
    n_host = sum(1 for g in control_taxonomy if g.top1_prefix == host_prefix)

    host_pidents = [
        g.top1_pident for g in control_taxonomy if g.top1_prefix == host_prefix
    ]

    # Count organisms
    organisms = [
        extract_organism(g.top1_target)
        for g in control_taxonomy
        if g.top1_prefix == host_prefix
    ]
    org_counts = Counter(organisms)
    dominant = org_counts.most_common(1)[0][0] if org_counts else "unknown"

    return ControlStats(
        n_genes=n_genes,
        n_no_hits=n_no_hits,
        no_hit_frequency=n_no_hits / n_genes,
        host_frequency=n_host / n_genes,
        mean_pident=sum(host_pidents) / len(host_pidents) if host_pidents else 0.0,
        dominant_organism=dominant,
        host_prefix=host_prefix,
    )


def build_taxonomy_consensus(
    control_taxonomy: list[GeneTaxonomy],
    taxonomy_lookup: Optional[dict[str, str]] = None,
    host_prefix: str = "EUK__",
    min_token_length: int = 3,
    weight_mode: str = "rank",
) -> str:
    """
    Build taxonomic consensus lineage from control genes.

    Uses taxonomy lookup table to get full taxonomy strings for control genes,
    splits at "|" to extract taxonomy levels, weights hits by rank/bitscore,
    and builds a consensus lineage showing the dominant taxonomic path.

    Args:
        control_taxonomy: List of GeneTaxonomy from control pORFs
        taxonomy_lookup: Taxonomy label lookup (organism_id -> full_taxonomy)
        host_prefix: Host taxonomy prefix (e.g., "EUK__")

    Returns:
        Consensus lineage string (e.g., "EUK|Discosea|Longamoebia|Thecamoebida|Stenamoeba")
        or "unknown" if no consensus can be built
    """
    if not control_taxonomy or not taxonomy_lookup:
        return "unknown"

    level_weights: list[dict[str, float]] = []

    for gene in control_taxonomy:
        if not gene.top10_targets:
            continue
        for rank, (target, bits) in enumerate(zip(gene.top10_targets, gene.top10_bits)):
            if not target:
                continue
            if host_prefix and not str(target).startswith(host_prefix):
                continue
            weight = compute_hit_weight(rank, bits, weight_mode)
            if weight <= 0:
                continue
            org_id = resolve_org_id(str(target).split("|", 1)[0], taxonomy_lookup)
            taxonomy_string = taxonomy_lookup.get(org_id, "")
            if not taxonomy_string:
                continue
            levels = [lvl.strip() for lvl in taxonomy_string.split("|") if lvl.strip()]
            for idx, level in enumerate(levels):
                if len(level) < min_token_length:
                    continue
                while len(level_weights) <= idx:
                    level_weights.append({})
                level_weights[idx][level] = level_weights[idx].get(level, 0.0) + weight

    if not level_weights:
        return "unknown"

    consensus_levels = []
    for level_dict in level_weights:
        if not level_dict:
            break
        token = max(level_dict.items(), key=lambda item: item[1])[0]
        consensus_levels.append(token)

    return "|".join(consensus_levels) if consensus_levels else "unknown"


def write_taxonomy_map(
    taxonomy_map: dict[str, GeneTaxonomy],
    output_path: Path,
) -> None:
    """
    Write taxonomy map to TSV file.

    Args:
        taxonomy_map: Dict mapping pORF ID to GeneTaxonomy
        output_path: Path to write TSV file
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with atomic_write_context(output_path, "w") as f:
        # Header
        f.write(
            "porf_id\tscaffold\tstart\tend\ttop1_target\ttop1_prefix\t"
            "top1_pident\ttop1_evalue\ttop10_prefixes\ttop10_targets\t"
            "top10_bits\ttop10_pidents\ttop10_evalues\t"
            "has_ncldv_mirus\thas_vp_plv\thas_viral\thas_hit\n"
        )
        # Data
        for porf_id, tax in sorted(taxonomy_map.items()):
            prefixes_str = ",".join(tax.top10_prefixes) if tax.top10_prefixes else ""
            targets_str = ",".join(tax.top10_targets) if tax.top10_targets else ""
            bits_str = ",".join(f"{b:.1f}" for b in tax.top10_bits) if tax.top10_bits else ""
            pidents_str = ",".join(f"{p:.2f}" for p in tax.top10_pidents) if tax.top10_pidents else ""
            evalues_str = ",".join(f"{e:.2e}" for e in tax.top10_evalues) if tax.top10_evalues else ""
            f.write(
                f"{tax.porf_id}\t{tax.scaffold}\t{tax.start}\t{tax.end}\t"
                f"{tax.top1_target}\t{tax.top1_prefix}\t{tax.top1_pident:.2f}\t"
                f"{tax.top1_evalue:.2e}\t{prefixes_str}\t{targets_str}\t"
                f"{bits_str}\t{pidents_str}\t{evalues_str}\t"
                f"{int(tax.has_ncldv_mirus)}\t{int(tax.has_vp_plv)}\t"
                f"{int(tax.has_viral)}\t{int(tax.has_hit)}\n"
            )

    logger.info("Wrote taxonomy map: %s", output_path)


def write_control_stats(
    control_stats: ControlStats,
    output_path: Path,
) -> None:
    """
    Write control statistics to JSON file.

    Args:
        control_stats: ControlStats object
        output_path: Path to write JSON file
    """
    import json

    output_path.parent.mkdir(parents=True, exist_ok=True)

    stats_dict = {
        "n_genes": control_stats.n_genes,
        "n_no_hits": control_stats.n_no_hits,
        "no_hit_frequency": control_stats.no_hit_frequency,
        "host_frequency": control_stats.host_frequency,
        "mean_pident": control_stats.mean_pident,
        "dominant_organism": control_stats.dominant_organism,
        "host_prefix": control_stats.host_prefix,
    }

    with atomic_write_context(output_path, "w") as f:
        json.dump(stats_dict, f, indent=2)

    logger.info("Wrote control stats: %s", output_path)
