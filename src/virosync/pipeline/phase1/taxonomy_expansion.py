"""
Taxonomy-Based Expansion for Low-Marker EVE Regions.

Validates single-marker or low-marker regions by expanding ±N genes
and checking if flanking genes have viral taxonomy in top-10 Diamond hits.

EVE Biology Context:
-------------------
EVEs (Endogenous Viral Elements) are ancient viral sequences integrated
into eukaryotic chromosomes. Because they are endogenized in host genomes:

1. **Host hits dominate top Diamond ranks** (higher bitscore, same genome)
2. **Viral signal appears at ranks 4-10** (lower bitscore, diverged)
3. **Cannot use bitscore margin or top-1 requirement** (would eliminate all EVEs)
4. **Must analyze top-10 DISTRIBUTION** (viral presence despite host dominance)

Example Diamond top-10 for EVE gene:
  Rank 1-3: EUK__ (host genome, 90-95% identity, 450-500 bits)
  Rank 4-7: NCLDV__ (viral origin, 70-80% identity, 350-400 bits)
  Rank 8-10: MIRUS__/PLV__/VP__ (viral, 65-75% identity, 300-350 bits)

Decision: Gene is VIRAL-POSITIVE (viral prefix present in top-10 at >=25% identity).

Confidence Scoring:
------------------
- **LOW**: <3 viral-positive genes in window (single marker + <2 flanking)
- **MEDIUM**: ≥3 viral-positive genes AND ≥2 non-marker flanking genes
- **HIGH**: MEDIUM criteria + region contains MCP marker
- **MCP Special**: MCP marker + 1 viral-positive flanking gene → MEDIUM

MCP (Major Capsid Protein) markers are highly specific for viral detection
and boost confidence even with lower flanking gene counts.
"""

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from Bio import SeqIO

from virosync.pipeline.phase0.prodigal import GenePrediction
from virosync.pipeline.search_backend import DIAMOND_DEFAULT_TIMEOUT

logger = logging.getLogger(__name__)

# Viral taxonomy prefixes for EVE detection
# Updated to match marker_validation.py validated prefixes
VALIDATED_PREFIXES = {"NCLDV__", "MIRUS__", "PLV__", "VP__", "PPV__", "CRESS__", "GVMAG__"}
MIN_VIRAL_HIT_PIDENT = 25.0


@dataclass
class TaxonomyExpansionResult:
    """Result of taxonomy expansion for a single region.

    Attributes:
        region: The candidate region being evaluated
        genes_tested: Total genes in flanking window (marker + flanking)
        genes_with_viral_taxonomy: Count of viral-positive genes (viral top-10 hit at >=25% identity)
        non_marker_viral_genes: Count of viral-positive genes excluding validated markers
        has_mcp: Whether region contains an MCP marker
        accepted: Whether region passed expansion thresholds
        expansion_confidence: Confidence level (LOW, MEDIUM, HIGH) based on evidence
        viral_gene_ids: Gene IDs classified as viral-positive
        acceptance_reason: Human-readable explanation of decision
        rejection_reason: Reason for rejection (if not accepted)
    """

    region: object  # CandidateRegion
    genes_tested: int
    genes_with_viral_taxonomy: int
    non_marker_viral_genes: int
    has_mcp: bool
    accepted: bool
    expansion_confidence: str  # LOW, MEDIUM, HIGH
    viral_gene_ids: list[str] = field(default_factory=list)
    acceptance_reason: str = ""
    rejection_reason: str = ""


# ============================================================================
# Batch Diamond Processing Helper Functions
# ============================================================================












def extract_all_flanking_genes_batched(
    candidate_regions: list,
    gene_data: dict[str, list[GenePrediction]],
    proteome_path: Path,
    flank_genes: int = 5,
    output_dir: Path = None,
) -> tuple[Path, dict[str, list[str]]]:
    """
    Extract ALL flanking genes from ALL regions into a SINGLE FASTA.

    This is Phase 1 of batch processing: consolidate all genes into one file
    for a single batched Diamond call.

    Args:
        candidate_regions: List of regions needing taxonomy expansion
        gene_data: Gene predictions by scaffold
        proteome_path: Path to proteome FASTA
        flank_genes: Number of genes on each side (default 5 = 11 total)
        output_dir: Output directory for batch FASTA

    Returns:
        (batch_fasta_path, gene_to_regions_map)
        - batch_fasta_path: FASTA with all genes (~4,378 proteins)
        - gene_to_regions_map: {gene_id: [region_id, ...]} for multi-mapping
    """
    if output_dir is None:
        output_dir = Path(proteome_path).parent / "taxonomy_expansion"

    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("-" * 60)
    logger.info("Phase 1: Extracting flanking genes (batched)")
    logger.info(f"  Regions: {len(candidate_regions)}")
    logger.info(f"  Flank genes: ±{flank_genes} per marker")

    # Load proteome once (efficiency)
    logger.info(f"  Loading proteome: {proteome_path.name}")
    proteome_seqs = {}
    try:
        for record in SeqIO.parse(proteome_path, "fasta"):
            proteome_seqs[record.id] = str(record.seq)
    except Exception as e:
        raise FileNotFoundError(
            f"Failed to load proteome {proteome_path}: {e}"
        ) from e

    logger.info(f"  Proteome loaded: {len(proteome_seqs)} sequences")

    # Track gene → regions mapping (handle multi-mapping)
    gene_to_regions = defaultdict(list)
    all_genes = set()
    total_extractions = 0

    for region in candidate_regions:
        flanking_genes = extract_flanking_genes(region, gene_data, flank_genes)
        total_extractions += len(flanking_genes)

        for gene in flanking_genes:
            gene_id = gene.gene_id
            all_genes.add(gene_id)
            gene_to_regions[gene_id].append(region.region_id)

    # Write batch FASTA
    batch_fasta = output_dir / "all_flanking_genes.faa"
    missing_genes = []

    with batch_fasta.open("w") as out:
        for gene_id in sorted(all_genes):
            if gene_id in proteome_seqs:
                out.write(f">{gene_id}\n{proteome_seqs[gene_id]}\n")
            else:
                missing_genes.append(gene_id)

    if missing_genes and len(missing_genes) == len(all_genes):
        raise ValueError(
            f"ALL genes missing from proteome! Check proteome/gene prediction mismatch. "
            f"Missing: {missing_genes[:5]}"
        )

    if missing_genes:
        logger.warning(
            f"  {len(missing_genes)}/{len(all_genes)} genes missing from proteome"
        )

    # Log stats
    unique_genes = len(all_genes) - len(missing_genes)
    multi_mapped = sum(1 for v in gene_to_regions.values() if len(v) > 1)

    logger.info(f"  Total gene extractions: {total_extractions}")
    logger.info(f"  Unique genes: {unique_genes}")
    logger.info(f"  Deduplication saved: {total_extractions - unique_genes} duplicates")
    logger.info(f"  Avg genes per region: {total_extractions / len(candidate_regions):.1f}")

    # MEDIUM FIX (Codex review): Guard division by zero
    if unique_genes > 0:
        logger.info(f"  Multi-mapped genes: {multi_mapped} ({multi_mapped*100.0/unique_genes:.1f}%)")
    else:
        logger.info(f"  Multi-mapped genes: {multi_mapped} (N/A - no genes extracted)")

    logger.info(f"  Batch FASTA: {batch_fasta}")

    return batch_fasta, dict(gene_to_regions)




def parse_diamond_top10(diamond_tsv: Path) -> dict[str, list[tuple[str, float]]]:
    """
    Parse Diamond TSV into gene-level top-10 hits.

    Args:
        diamond_tsv: Merged Diamond results TSV

    Returns:
        {gene_id: [(target1, pident1), ..., (target10, pident10)]} dictionary
    """
    gene_top10 = defaultdict(list)

    if not diamond_tsv.exists() or diamond_tsv.stat().st_size == 0:
        logger.warning(f"Diamond results empty or missing: {diamond_tsv}")
        return dict(gene_top10)

    with diamond_tsv.open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            query, target = parts[0], parts[1]
            try:
                pident = float(parts[4]) if len(parts) > 4 else 0.0
            except ValueError:
                pident = 0.0
            gene_top10[query].append((target, pident))

    # Limit to top-10 (should already be limited by Diamond --max-target-seqs)
    for gene_id in gene_top10:
        gene_top10[gene_id] = gene_top10[gene_id][:10]

    logger.debug(f"Parsed Diamond results: {len(gene_top10)} genes")

    return dict(gene_top10)


def classify_region_from_gene_hits(
    region: object,
    gene_top10: dict[str, list[tuple[str, float]]],
    gene_to_region_map: dict[str, list[str]],
    gene_data: dict[str, list[GenePrediction]],
    min_viral_genes_total: int = 3,
    min_viral_genes_non_marker: int = 2,
    short_scaffold_min_fraction: float = 0.20,
    flank_genes: int = 5,
) -> TaxonomyExpansionResult:
    """
    Classify a single region using pre-computed gene-level Diamond hits.

    This is Phase 3 of batch processing: classify regions using the merged
    Diamond results from Phase 2.

    Args:
        region: Candidate region to classify
        gene_top10: Pre-computed gene-level top-10 hits {gene_id: [(target, pident)]}
        gene_to_region_map: Multi-mapping {gene_id: [region_ids]}
        gene_data: Gene predictions by scaffold
        min_viral_genes_total: Min total viral genes (default 3)
        min_viral_genes_non_marker: Min non-marker viral genes (default 2)
        short_scaffold_min_fraction: Viral fraction for short scaffolds (default 0.20)
        flank_genes: Number of flanking genes (for logging)

    Returns:
        TaxonomyExpansionResult with classification
    """
    # Extract flanking genes for this region
    flanking_genes = extract_flanking_genes(region, gene_data, flank_genes)

    if len(flanking_genes) < 2:
        return TaxonomyExpansionResult(
            region=region,
            genes_tested=len(flanking_genes),
            genes_with_viral_taxonomy=0,
            non_marker_viral_genes=0,
            has_mcp=False,
            accepted=False,
            expansion_confidence="LOW",
            viral_gene_ids=[],
            rejection_reason=f"Insufficient genes ({len(flanking_genes)} < 2)",
        )

    # Classify genes using pre-computed hits
    marker_gene_ids = {m.query_porf for m in region.markers}
    viral_positive_genes = []

    for gene in flanking_genes:
        gene_id = gene.gene_id
        top10_hits = gene_top10.get(gene_id, [])

        if is_gene_viral_positive(top10_hits, VALIDATED_PREFIXES):
            viral_positive_genes.append(gene)

    # Counts
    total_genes = len(flanking_genes)
    total_viral = len(viral_positive_genes)
    non_marker_viral = sum(
        1 for g in viral_positive_genes if g.gene_id not in marker_gene_ids
    )

    # Check for MCP markers
    has_mcp = any(m.is_mcp for m in region.markers)

    # Apply thresholds (same logic as single-region expansion)
    accepted = False
    confidence = "LOW"
    acceptance_reason = ""
    rejection_reason = ""

    expected_window = 2 * flank_genes + 1

    # Default threshold (expected window size)
    if total_genes >= expected_window:
        if total_viral >= min_viral_genes_total and non_marker_viral >= min_viral_genes_non_marker:
            accepted = True
            confidence = "HIGH" if has_mcp else "MEDIUM"
            acceptance_reason = (
                f"{total_viral}/{expected_window} genes viral, "
                f"{non_marker_viral} non-marker"
                + (", MCP present → HIGH" if has_mcp else " → MEDIUM")
            )
        elif has_mcp and total_viral >= 2 and non_marker_viral >= 1:
            accepted = True
            confidence = "MEDIUM"
            acceptance_reason = (
                f"{total_viral}/{expected_window} genes viral, "
                f"{non_marker_viral} non-marker, MCP present → MEDIUM (MCP boost)"
            )
        else:
            rejection_reason = (
                f"Insufficient viral genes ({total_viral}/{total_genes}, "
                f"{non_marker_viral} non-marker < {min_viral_genes_non_marker})"
            )

    # Short scaffold threshold
    elif total_genes >= 2:
        viral_fraction = total_viral / total_genes

        if (
            viral_fraction >= short_scaffold_min_fraction
            and total_viral >= 2
            and non_marker_viral >= 1
        ):
            accepted = True
            confidence = "HIGH" if (has_mcp and total_viral >= 3) else "MEDIUM"
            acceptance_reason = (
                f"{total_viral}/{total_genes} genes viral ({viral_fraction:.1%}), "
                f"{non_marker_viral} non-marker"
                + (" → HIGH" if (has_mcp and total_viral >= 3) else " → MEDIUM")
            )
        else:
            rejection_reason = (
                f"Short scaffold: {total_viral}/{total_genes} genes viral "
                f"({viral_fraction:.1%} < {short_scaffold_min_fraction:.0%})"
            )

    else:
        rejection_reason = f"Insufficient genes ({total_genes} < 2)"

    return TaxonomyExpansionResult(
        region=region,
        genes_tested=total_genes,
        genes_with_viral_taxonomy=total_viral,
        non_marker_viral_genes=non_marker_viral,
        has_mcp=has_mcp,
        accepted=accepted,
        expansion_confidence=confidence,
        viral_gene_ids=[g.gene_id for g in viral_positive_genes],
        acceptance_reason=acceptance_reason,
        rejection_reason=rejection_reason,
    )


def extract_flanking_genes(
    region: object,  # CandidateRegion
    gene_data: dict[str, list[GenePrediction]],
    flank_genes: int = 5,
) -> list[GenePrediction]:
    """
    Extract ±flank_genes around all markers in a region.

    For each validated marker in the region, extracts N genes upstream
    and N genes downstream. Combines into a single deduplicated list.

    Args:
        region: Candidate EVE region with validated markers
        gene_data: Gene predictions by scaffold (from load_gene_predictions)
        flank_genes: Number of genes to extract on each side (default 5)

    Returns:
        List of GenePrediction objects for flanking window
        Typically 11 genes = 1 marker + 5 upstream + 5 downstream
    """
    scaffold_genes = gene_data.get(region.scaffold, [])
    if not scaffold_genes:
        return []

    # Get all marker positions
    marker_positions = set()
    for marker in region.markers:
        marker_positions.add((marker.start, marker.end))

    # Find indices of genes overlapping markers
    marker_indices = set()
    for idx, gene in enumerate(scaffold_genes):
        if (gene.start, gene.end) in marker_positions:
            marker_indices.add(idx)

    if not marker_indices:
        return []

    # Expand ±flank_genes from each marker
    flanking_indices = set()
    for idx in marker_indices:
        start_idx = max(0, idx - flank_genes)
        end_idx = min(len(scaffold_genes), idx + flank_genes + 1)
        flanking_indices.update(range(start_idx, end_idx))

    # Return all genes in flanking window (sorted, deduplicated)
    return [scaffold_genes[idx] for idx in sorted(flanking_indices)]


def is_gene_viral_positive(
    top10_hits: list[str | tuple[str, float]],
    validated_prefixes: set[str] = VALIDATED_PREFIXES,
    min_pident: float = MIN_VIRAL_HIT_PIDENT,
) -> bool:
    """
    Classify gene as viral-positive based on top-10 Diamond hits.

    EVE-specific rule: Gene is viral-positive if any viral prefix appears
    in top-10 with percent identity >= 25.0, regardless of rank or bitscore.

    Rationale:
    - Host hits dominate top ranks for endogenized viral genes
    - Viral signal appears at ranks 4-10 (lower bitscore due to divergence)
    - Bitscore margin or top-1 requirement would systematically eliminate EVEs

    Args:
        top10_hits: List of target IDs or (target ID, pident) tuples from Diamond (max 10)
        validated_prefixes: Set of viral taxonomy prefixes
        min_pident: Minimum percent identity for a viral hit

    Returns:
        True if any identity-qualified viral prefix is found in top-10
    """
    for hit in top10_hits:
        if isinstance(hit, (tuple, list)):
            target = str(hit[0]) if hit else ""
            try:
                pident = float(hit[1])
            except (IndexError, TypeError, ValueError):
                pident = 0.0
        else:
            target = str(hit)
            pident = 0.0
        # Extract prefix (format: "PREFIX__Family|Genus|...")
        if "__" in target:
            prefix = target.split("__")[0] + "__"
            if prefix in validated_prefixes and pident >= min_pident:
                return True
    return False




def filter_regions_by_taxonomy_expansion(
    candidate_regions: list,
    proteome_path: Path,
    gene_taxonomy_db: Path,
    gene_data: dict[str, list[GenePrediction]],
    marker_count_threshold: int = 1,
    flank_genes: int = 5,
    min_viral_genes_total: int = 3,
    min_viral_genes_non_marker: int = 2,
    short_scaffold_min_fraction: float = 0.20,
    batch_diamond: bool = False,
    threads: int = 4,
    output_dir: Optional[Path] = None,
    search_backend: str = "diamond",
) -> tuple[list, list[TaxonomyExpansionResult]]:
    """
    Filter low-marker regions by taxonomy expansion.

    Regions with >marker_count_threshold markers are kept automatically.
    Regions with ≤marker_count_threshold markers undergo taxonomy expansion.

    Example: threshold=1 means single-marker regions (≤1) get expanded,
    multi-marker regions (>1) are auto-accepted.


    Args:
        candidate_regions: Regions from region assembly
        proteome_path: Path to proteome FASTA
        gene_taxonomy_db: Path to combined_proteome.dmnd (40GB)
        gene_data: Gene predictions by scaffold
        marker_count_threshold: Regions with ≤N markers undergo expansion (default 1)
        flank_genes: Number of genes to check on each side (default 5)
        min_viral_genes_total: Minimum total viral genes (default 3)
        min_viral_genes_non_marker: Minimum non-marker viral genes (default 2)
        short_scaffold_min_fraction: Viral fraction for short scaffolds (default 0.20)
        batch_diamond: Deprecated no-op; see taxonomy_expansion_task.
        threads: Number of threads
        output_dir: Optional output directory

    Returns:
        Tuple of (filtered_regions, expansion_results)
    """
    # Fail-fast validation: check database exists before processing
    if not gene_taxonomy_db.exists():
        raise FileNotFoundError(
            f"Gene taxonomy database not found: {gene_taxonomy_db}. "
            "Cannot perform taxonomy expansion without combined_proteome.dmnd database."
        )

    # CRITICAL FIX (Codex review): Ensure output_dir is never None
    if output_dir is None:
        output_dir = Path(proteome_path).parent / "phase1" / "taxonomy_expansion"

    output_dir.mkdir(parents=True, exist_ok=True)

    filtered_regions = []
    expansion_results = []

    # Track statistics
    auto_accepted = 0
    expansion_tested = 0
    expansion_accepted_high = 0
    expansion_accepted_medium = 0
    expansion_accepted_low = 0
    expansion_rejected = 0

    # Separate high-marker (auto-accept) from low-marker (need expansion)
    # Semantics: regions with > threshold markers auto-accept, ≤ threshold undergo expansion
    # Example: threshold=1 means single-marker regions (≤1) get expanded, multi-marker (>1) auto-accept
    high_marker_regions = []
    low_marker_regions = []

    for region in candidate_regions:
        if region.marker_count > marker_count_threshold:
            high_marker_regions.append(region)
        else:
            low_marker_regions.append(region)

    # Auto-accept high-marker regions
    filtered_regions.extend(high_marker_regions)
    auto_accepted = len(high_marker_regions)

    # Process low-marker regions
    if not low_marker_regions:
        logger.info("No low-marker regions to expand")
        return filtered_regions, expansion_results

    logger.info(
        f"Taxonomy expansion: {len(candidate_regions)} regions → "
        f"{len(high_marker_regions)} auto-accepted (>{marker_count_threshold} markers), "
        f"{len(low_marker_regions)} require expansion (≤{marker_count_threshold} markers)"
    )

    # ========== BATCH PROCESSING ==========
    # Replace serial processing with batched Diamond approach

    # Phase 1: Extract ALL flanking genes
    batch_fasta, gene_to_region_map = extract_all_flanking_genes_batched(
        candidate_regions=low_marker_regions,
        gene_data=gene_data,
        proteome_path=proteome_path,
        flank_genes=flank_genes,
        output_dir=output_dir,
    )

    # Phase 2: Run batched Diamond BLASTP
    # SIMPLIFIED: Single Diamond call with all threads (no chunking)
    # Chunking causes I/O contention with large databases - better to use all threads in one process
    logger.info("-" * 60)
    logger.info("Phase 2: Running batched Diamond BLASTP (single process)")

    merged_diamond_tsv = output_dir / "all_flanking_diamond.tsv"

    # Check if already completed (checkpointing)
    if merged_diamond_tsv.exists() and merged_diamond_tsv.stat().st_size > 0:
        logger.info(f"Diamond results already exist (checkpoint): {merged_diamond_tsv}")
    else:
        # Run Diamond with all threads in single process
        logger.info(f"  Query: {batch_fasta.name}")
        logger.info(f"  Database: {gene_taxonomy_db.name} ({gene_taxonomy_db.stat().st_size / 1e9:.1f} GB)")
        logger.info(f"  Threads: {threads}")
        logger.info(f"  Timeout: {DIAMOND_DEFAULT_TIMEOUT / 60:.0f} minutes")

        from virosync.pipeline.phase1.marker_validation import run_diamond_on_hmm_hits

        try:
            run_diamond_on_hmm_hits(
                hmm_hit_fasta=batch_fasta,
                diamond_db=gene_taxonomy_db,
                output_tsv=merged_diamond_tsv,
                threads=threads,
                evalue=1e-5,
                max_seqs=10,
                search_backend=search_backend,
            )
            logger.info(f"  Diamond complete: {merged_diamond_tsv.stat().st_size / 1e6:.1f} MB")
        except Exception as e:
            logger.error(f"Diamond failed: {e}")
            raise

    # Phase 3: Parse results and classify regions
    logger.info("-" * 60)
    logger.info("Phase 3: Classifying regions by flanking gene taxonomy")

    gene_top10 = parse_diamond_top10(merged_diamond_tsv)

    try:
        from tqdm import tqdm
        has_tqdm = True
    except ImportError:
        has_tqdm = False

    iterator = low_marker_regions
    if has_tqdm:
        iterator = tqdm(low_marker_regions, desc="Classifying regions")

    for region in iterator:
        expansion_tested += 1

        result = classify_region_from_gene_hits(
            region=region,
            gene_top10=gene_top10,
            gene_to_region_map=gene_to_region_map,
            gene_data=gene_data,
            min_viral_genes_total=min_viral_genes_total,
            min_viral_genes_non_marker=min_viral_genes_non_marker,
            short_scaffold_min_fraction=short_scaffold_min_fraction,
            flank_genes=flank_genes,
        )

        expansion_results.append(result)

        # CHANGED (Jan 2026): Include ALL regions (accepted and rejected)
        # Rejected regions are marked as LOW confidence but still processed through Phase 2/3
        # This allows host gene comparison and full evidence synthesis for all candidates
        if result.accepted:
            # Set confidence from expansion result
            region.taxonomy_expansion_confidence = result.expansion_confidence
            filtered_regions.append(region)

            if result.expansion_confidence == "HIGH":
                expansion_accepted_high += 1
            elif result.expansion_confidence == "MEDIUM":
                expansion_accepted_medium += 1
            else:
                expansion_accepted_low += 1

            logger.debug(
                f"Taxonomy expansion: {result.expansion_confidence} confidence - "
                f"{region.region_id} ({result.genes_with_viral_taxonomy}/"
                f"{result.genes_tested} genes viral, {result.non_marker_viral_genes} "
                f"non-marker, MCP={result.has_mcp})"
            )
        else:
            # Mark rejected regions as LOW confidence but still include them
            region.taxonomy_expansion_confidence = "LOW"
            filtered_regions.append(region)
            expansion_rejected += 1

            logger.debug(
                f"Taxonomy expansion: LOW confidence (insufficient evidence) - "
                f"{region.region_id} ({result.genes_with_viral_taxonomy}/"
                f"{result.genes_tested} genes viral) - {result.rejection_reason}"
            )

    logger.info(
        f"Taxonomy expansion: {len(candidate_regions)} total → "
        f"{len(filtered_regions)} kept (ALL regions kept, marked with confidence: "
        f"auto_accepted={auto_accepted}, "
        f"expansion_HIGH={expansion_accepted_high}, expansion_MEDIUM={expansion_accepted_medium}, "
        f"expansion_LOW_accepted={expansion_accepted_low}, expansion_LOW_insufficient={expansion_rejected})"
    )

    # Write summary TSV
    if output_dir:
        summary_file = output_dir / "taxonomy_expansion_summary.tsv"
        with summary_file.open("w") as handle:
            handle.write(
                "region_id\tscaffold\tstart\tend\tmarker_count\t"
                "genes_tested\tgenes_viral\tnon_marker_viral\thas_mcp\t"
                "accepted\tconfidence\tacceptance_reason\trejection_reason\n"
            )
            for result in expansion_results:
                handle.write(
                    f"{result.region.region_id}\t{result.region.scaffold}\t"
                    f"{result.region.start}\t{result.region.end}\t"
                    f"{result.region.marker_count}\t"
                    f"{result.genes_tested}\t{result.genes_with_viral_taxonomy}\t"
                    f"{result.non_marker_viral_genes}\t{int(result.has_mcp)}\t"
                    f"{int(result.accepted)}\t{result.expansion_confidence}\t"
                    f"{result.acceptance_reason}\t{result.rejection_reason}\n"
                )
        logger.info(f"Wrote taxonomy expansion summary: {summary_file}")

    return filtered_regions, expansion_results
