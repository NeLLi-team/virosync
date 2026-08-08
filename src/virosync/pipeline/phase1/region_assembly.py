"""
Iterative Marker-Driven Region Assembly.

Step 4 of Phase 1 Rewrite: Assembles candidate EVE regions from validated marker hits
using iterative extension until no more markers can be captured.

This module:
1. Groups validated markers within initial clustering windows
2. Iteratively extends regions to capture adjacent markers
3. Merges overlapping regions
4. Outputs marker seed regions for Phase 2

Key concepts:
- Initial clustering: Group markers within initial_window_bp OR initial_window_genes
- Iterative extension: Extend by extension_kb, check for new markers, repeat until no more markers found
- Validated marker: HMM hit confirmed by Diamond top-5 taxonomy (NCLDV/MIRUS prefix)
- Candidate region: A genomic interval containing ≥2 validated markers after extension

Algorithm:
The iterative extension captures more EVE content by discovering markers that
weren't initially clustered but fall within the extended boundaries. This is crucial
for degraded or fragmented EVEs where markers may be spaced further apart.
"""

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from virosync.ablation import AblationID
from virosync.pipeline.phase0.prodigal import GenePrediction, load_gene_predictions
from virosync.pipeline.phase1.marker_roles import decide_marker_hit_role
from virosync.pipeline.phase1.viral_markers import (
    base_marker_gene_id,
    is_identity_qualified_cress_marker,
)

logger = logging.getLogger(__name__)


@dataclass
class ValidatedMarkerHit:
    """Represents a validated HMM marker hit with Diamond support.

    This dataclass mirrors marker_validation.ValidatedMarkerHit for compatibility.
    """

    query_porf: str  # pORF ID (renamed from porf_id to match marker_validation.py)
    scaffold: str
    start: int  # Nucleotide start
    end: int  # Nucleotide end
    strand: str  # Strand (+ or -)
    hmm_target: str  # HMM profile name (e.g., GVOGm0003, OG1234)
    hmm_score: float  # HMM bit score
    hmm_evalue: float  # HMM E-value
    validation_status: str  # "validated", "validated_novel", "supported", "unvalidated"
    top10_prefixes: str  # Comma-separated taxonomy prefixes (changed from list[str] to match marker_validation.py)
    best_hit_target: str  # Best Diamond hit target ID
    best_hit_pident: float  # Percent identity of best Diamond hit
    best_hit_bits: float  # Bit score of best Diamond hit
    has_ncldv: int
    has_mirus: int
    has_plv: int
    has_vp: int
    has_viral: int
    # Optional fields with defaults
    top10_targets: str = ""
    top10_pidents: str = ""
    top10_bitscores: str = ""
    top10_evalues: str = ""
    taxonomy_substring_counts: str = ""
    taxonomy_raw_counts: str = ""
    tier1_bypassed: bool = False

    @property
    def is_validated(self) -> bool:
        """Check if marker is validated (NCLDV/MIRUS in top-5) or validated_novel (HMM-only)."""
        return self.tier1_bypassed or self.validation_status in (
            "validated",
            "validated_novel",
        )

    @property
    def is_mcp(self) -> bool:
        """Check if marker is a Major Capsid Protein (MCP) marker.

        MCP markers include:
        - GVOGm0003 (NCLDV MCP)
        - VS000086/OG1352, VS000309/OG484 (NCLDV MCP PFAM - Group II dsDNA virus capsid)
        - gamadvirusMCP (Gamadnavirus MCP)
        - PLV_MCP (Polinton-like virus MCP)
        - VP_MCP (Virophage MCP)
        - Mirus_MCP (Mirus virus MCP)
        - Any marker with "mcp" substring in name

        Returns:
            True if marker is an MCP marker
        """
        from virosync.pipeline.phase3.mcp_detection import is_mcp_gene
        return is_mcp_gene(self.hmm_target)

    @property
    def is_gvogm(self) -> bool:
        """Check if marker is GVOGm."""
        return self.hmm_target.lower().startswith("gvogm")

    @property
    def is_valid_seed_marker(self) -> bool:
        """Check if marker is a valid seed for region assembly.

        UPDATED (Jan 2026): Expanded to allow ALL validated markers to seed regions,
        with taxonomy-based expansion filtering false positives in Step 4.5.

        Seeding rules:
        1. validation_status="validated" (NCLDV/MIRUS/PLV/VP in Diamond top-10) → CAN SEED
        2. validation_status="validated_novel" (HMM-only, no Diamond hits) → CAN SEED ONLY IF MCP

        Rationale:
        - All HMM hits are validated via marker_validation.py (Diamond top-10 taxonomy)
        - validated_novel has InterProScan/TMVec support but no Diamond hits
        - MCP markers are highly specific for viral detection
        - Low-marker regions undergo taxonomy expansion validation (Step 4.5)

        Returns:
            True if marker can seed an EVE region
        """
        # Validated markers with Diamond top-10 viral taxonomy → auto-seed
        if self.tier1_bypassed:
            return True

        if self.validation_status == "validated":
            return True

        # validated_novel (HMM-only, no Diamond hits) → ONLY MCP can seed
        if self.validation_status == "validated_novel":
            return self.is_mcp  # Use expanded MCP detection

        # unvalidated or supported → cannot seed
        return False


@dataclass
class CandidateRegion:
    """Represents a candidate EVE region assembled from validated markers."""

    scaffold: str
    start: int  # Nucleotide start (0-based)
    end: int  # Nucleotide end (exclusive)
    markers: list[ValidatedMarkerHit] = field(default_factory=list)
    region_id: str = ""
    has_mcp: bool = False
    # Region classification based on seed markers (NCLDV, VP, PLV, MIRUS, MIXED, UNKNOWN)
    predicted_family: str = ""
    # Confidence from taxonomy expansion (LOW, MEDIUM, HIGH, or None for auto-accepted)
    taxonomy_expansion_confidence: Optional[str] = None

    @property
    def length(self) -> int:
        """Region length in bp."""
        return self.end - self.start

    @property
    def marker_count(self) -> int:
        """Number of markers in region."""
        return len(self.markers)

    @property
    def marker_types(self) -> set[str]:
        """Unique HMM marker types in region."""
        return set(m.hmm_target for m in self.markers)

    @property
    def marker_types_str(self) -> str:
        """Comma-separated list of marker types."""
        return ",".join(sorted(self.marker_types))

    def update_has_mcp(self) -> None:
        """Update MCP flag based on markers."""
        self.has_mcp = any(m.is_mcp for m in self.markers)


def count_genes_between(
    marker1: ValidatedMarkerHit,
    marker2: ValidatedMarkerHit,
    gene_order: dict[str, list[GenePrediction]],
) -> int:
    """
    Count the number of genes between two markers on the same scaffold.

    Args:
        marker1: First marker hit
        marker2: Second marker hit
        gene_order: Gene order by scaffold from load_gene_predictions()

    Returns:
        Number of genes between markers (0 if adjacent or on different scaffolds)
    """
    if marker1.scaffold != marker2.scaffold:
        return float("inf")  # Different scaffolds

    scaffold_genes = gene_order.get(marker1.scaffold, [])
    if not scaffold_genes:
        return 0

    start = min(marker1.start, marker2.start)
    end = max(marker1.end, marker2.end)
    count = 0
    for gene in scaffold_genes:
        if gene.start >= start and gene.end <= end:
            count += 1
    return max(0, count - 2)


def initial_clustering(
    validated_hits: list[ValidatedMarkerHit],
    scaffold: str,
    initial_window_bp: int,
    initial_window_genes: int,
    min_markers_initial: int,
    gene_order: dict[str, list[GenePrediction]],
) -> list[list[ValidatedMarkerHit]]:
    """
    Group validated markers within initial clustering windows.

    Markers are clustered if they are within initial_window_bp (bp distance)
    OR within initial_window_genes (gene distance).

    Args:
        validated_hits: List of validated marker hits for this scaffold
        scaffold: Scaffold name
        initial_window_bp: Maximum bp distance between markers for clustering
        initial_window_genes: Maximum gene distance between markers for clustering
        min_markers_initial: Minimum markers required to form a cluster (must include
            at least one is_valid_seed_marker)
        gene_order: Gene order by scaffold

    Returns:
        List of marker clusters (each cluster is a list of ValidatedMarkerHit)
    """
    if not validated_hits:
        return []

    # Sort by genomic position
    hits = sorted(validated_hits, key=lambda h: h.start)

    clusters = []
    current_cluster = [hits[0]]

    for hit in hits[1:]:
        # Check distance from last marker in current cluster
        last_marker = current_cluster[-1]

        # Distance in bp
        distance_bp = hit.start - last_marker.end

        # Distance in genes
        distance_genes = count_genes_between(last_marker, hit, gene_order)

        # Cluster if within either window
        if distance_bp <= initial_window_bp or distance_genes <= initial_window_genes:
            current_cluster.append(hit)
        else:
            # Finalize current cluster
            if len(current_cluster) >= min_markers_initial:
                if any(m.is_valid_seed_marker for m in current_cluster):
                    clusters.append(current_cluster)
            elif len(current_cluster) == 1:
                marker = current_cluster[0]
                if marker.is_valid_seed_marker:
                    clusters.append(current_cluster)
            # Start new cluster
            current_cluster = [hit]

    # Don't forget last cluster
    if len(current_cluster) >= min_markers_initial:
        if any(m.is_valid_seed_marker for m in current_cluster):
            clusters.append(current_cluster)
    elif len(current_cluster) == 1:
        marker = current_cluster[0]
        if marker.is_valid_seed_marker:
            clusters.append(current_cluster)

    logger.debug(
        f"  {scaffold}: {len(hits)} markers → {len(clusters)} initial clusters "
        f"(window: {initial_window_bp} bp or {initial_window_genes} genes, min: {min_markers_initial} markers)"
    )

    return clusters


def find_validated_markers_in_range(
    all_hits: list[ValidatedMarkerHit],
    scaffold: str,
    start: int,
    end: int,
    exclude_ids: set[str],
) -> list[ValidatedMarkerHit]:
    """
    Find validated markers within a genomic range, excluding already-included markers.

    Args:
        all_hits: All validated marker hits
        scaffold: Scaffold name
        start: Range start (bp)
        end: Range end (bp)
        exclude_ids: Set of query IDs already in the region

    Returns:
        List of new markers within range
    """
    new_markers = []

    for hit in all_hits:
        if hit.scaffold != scaffold:
            continue
        if hit.query_porf in exclude_ids:
            continue
        # Check if marker overlaps range [start, end)
        if hit.start < end and hit.end > start:
            new_markers.append(hit)

    return new_markers


def iterative_extension(
    cluster: list[ValidatedMarkerHit],
    scaffold: str,
    all_validated_hits: list[ValidatedMarkerHit],
    extension_kb: int,
    scaffold_length: int,
) -> CandidateRegion:
    """
    Iteratively extend a marker cluster until no more markers can be captured.

    Algorithm:
    1. Start with initial cluster markers
    2. Extend region by extension_kb from outermost markers
    3. Check if extension captures new validated markers
    4. If yes, add them and repeat from step 2
    5. If no, finalize region boundaries

    Args:
        cluster: Initial marker cluster
        scaffold: Scaffold name
        all_validated_hits: All validated hits on this scaffold
        extension_kb: Extension distance from outermost markers (kb)
        scaffold_length: Length of scaffold (bp) for boundary clamping

    Returns:
        CandidateRegion with final boundaries and all captured markers
    """
    extension_bp = extension_kb * 1000

    # Initialize region from cluster
    region = CandidateRegion(
        scaffold=scaffold,
        start=min(h.start for h in cluster) - extension_bp,
        end=max(h.end for h in cluster) + extension_bp,
        markers=cluster.copy(),
    )

    # Clamp to scaffold boundaries
    region.start = max(0, region.start)
    region.end = min(scaffold_length, region.end)

    # Enforce minimum region length (2 * extension)
    min_len = extension_bp * 2
    if region.length < min_len:
        deficit = min_len - region.length
        grow_left = deficit // 2
        grow_right = deficit - grow_left
        region.start = max(0, region.start - grow_left)
        region.end = min(scaffold_length, region.end + grow_right)
        if region.length < min_len:
            # If still short due to boundaries, extend where possible
            if region.start == 0:
                region.end = min(scaffold_length, region.end + (min_len - region.length))
            elif region.end == scaffold_length:
                region.start = max(0, region.start - (min_len - region.length))

    # Track included marker IDs
    included_ids = {m.query_porf for m in region.markers}

    iteration = 0
    max_iterations = 100  # Safety limit

    while iteration < max_iterations:
        iteration += 1

        # Find new markers within current boundaries
        new_markers = find_validated_markers_in_range(
            all_hits=all_validated_hits,
            scaffold=scaffold,
            start=region.start,
            end=region.end,
            exclude_ids=included_ids,
        )

        if not new_markers:
            # No more markers to capture - stop iteration
            break

        # Add new markers to region
        region.markers.extend(new_markers)
        included_ids.update(m.query_porf for m in new_markers)

        # Update boundaries to outermost marker ± extension
        outermost_start = min(h.start for h in region.markers)
        outermost_end = max(h.end for h in region.markers)

        region.start = outermost_start - extension_bp
        region.end = outermost_end + extension_bp

        # Clamp to scaffold boundaries
        region.start = max(0, region.start)
        region.end = min(scaffold_length, region.end)

        logger.debug(
            f"    Iteration {iteration}: captured {len(new_markers)} new markers, "
            f"region now {region.start}-{region.end} ({region.marker_count} total markers)"
        )

    if iteration >= max_iterations:
        logger.warning(f"Iterative extension hit max iterations ({max_iterations}) for {scaffold}")

    # Update MCP flag
    region.update_has_mcp()

    logger.debug(
        f"  Final region: {scaffold}:{region.start}-{region.end} "
        f"({region.length:,} bp, {region.marker_count} markers, {len(region.marker_types)} types)"
    )

    return region


def assemble_compact_cress_regions(
    cress_hits: list[ValidatedMarkerHit],
    gene_order: dict[str, list[GenePrediction]],
    *,
    max_gap_bp: int = 10_000,
    max_intervening_genes: int = 1,
) -> list[CandidateRegion]:
    """Build exact gene-bounded regions for compact CRESS insertions."""

    best_by_gene: dict[str, ValidatedMarkerHit] = {}
    for hit in cress_hits:
        if not is_identity_qualified_cress_marker(hit):
            continue
        gene_id = base_marker_gene_id(hit.query_porf)
        current = best_by_gene.get(gene_id)
        if current is None or hit.hmm_score > current.hmm_score:
            best_by_gene[gene_id] = hit

    hits_by_scaffold: dict[str, list[ValidatedMarkerHit]] = defaultdict(list)
    for hit in best_by_gene.values():
        hits_by_scaffold[hit.scaffold].append(hit)

    regions: list[CandidateRegion] = []
    for scaffold, scaffold_hits in sorted(hits_by_scaffold.items()):
        ordered = sorted(scaffold_hits, key=lambda hit: (hit.start, hit.end))
        cluster = [ordered[0]]
        for hit in ordered[1:]:
            previous = cluster[-1]
            gap_bp = max(0, hit.start - previous.end)
            intervening_genes = count_genes_between(
                previous,
                hit,
                gene_order,
            )
            if (
                gap_bp <= max_gap_bp
                and intervening_genes <= max_intervening_genes
            ):
                cluster.append(hit)
                continue
            regions.append(
                CandidateRegion(
                    scaffold=scaffold,
                    start=min(marker.start for marker in cluster),
                    end=max(marker.end for marker in cluster),
                    markers=list(cluster),
                    predicted_family="CRESS",
                )
            )
            cluster = [hit]
        regions.append(
            CandidateRegion(
                scaffold=scaffold,
                start=min(marker.start for marker in cluster),
                end=max(marker.end for marker in cluster),
                markers=list(cluster),
                predicted_family="CRESS",
            )
        )

    return regions


def merge_overlapping_regions(
    regions: list[CandidateRegion],
    merge_distance: int = 1000,
) -> list[CandidateRegion]:
    """
    Merge overlapping or adjacent regions on the same scaffold.

    Args:
        regions: List of candidate regions
        merge_distance: Maximum gap between regions to merge (bp)

    Returns:
        List of merged regions
    """
    if not regions:
        return []

    # Group by scaffold
    scaffold_regions = defaultdict(list)
    for region in regions:
        scaffold_regions[region.scaffold].append(region)

    merged = []

    for scaffold, scaffold_region_list in scaffold_regions.items():
        # Sort by start position
        scaffold_region_list.sort(key=lambda r: r.start)

        current = scaffold_region_list[0]

        for region in scaffold_region_list[1:]:
            # Check if regions overlap or are within merge distance
            if region.start <= current.end + merge_distance:
                # Merge regions
                current = CandidateRegion(
                    scaffold=scaffold,
                    start=current.start,
                    end=max(current.end, region.end),
                    markers=current.markers + region.markers,
                )
                current.update_has_mcp()
            else:
                # Finalize current and start new
                merged.append(current)
                current = region

        # Don't forget last region
        merged.append(current)

    logger.info(f"Merged {len(regions)} regions → {len(merged)} final regions")

    return merged


def load_scaffold_lengths(genome_fasta: Path) -> dict[str, int]:
    """
    Load scaffold lengths from genome FASTA.

    Args:
        genome_fasta: Path to genome FASTA file

    Returns:
        Dictionary mapping scaffold name to length (bp)
    """
    from Bio import SeqIO

    scaffold_lengths = {}
    for record in SeqIO.parse(genome_fasta, "fasta"):
        scaffold_lengths[record.id] = len(record.seq)

    return scaffold_lengths


def assemble_candidate_regions(
    validated_hits: list[ValidatedMarkerHit],
    genome_fasta: Path,
    proteome_fasta: Path,
    output_dir: Path,
    initial_window_bp: int = 10000,
    initial_window_genes: int = 5,
    min_markers_initial: int = 2,
    extension_kb: int = 5,
    merge_distance: int = 1000,
    ablation_id: AblationID = AblationID.A0,
    single_marker_min_score: float = 50.0,
    write_outputs: bool = True,
) -> list[CandidateRegion]:
    """
    Assemble candidate EVE regions from validated marker hits using iterative extension.

    This is the main pipeline function implementing Step 4 of Phase 1.

    Algorithm:
    1. Load scaffold lengths and gene order
    2. Group markers by scaffold
    3. For each scaffold:
       a. Initial clustering: group markers within windows
       b. Iterative extension: extend each cluster until no more markers found
    4. Merge overlapping regions
    5. Write output BED file

    Args:
        validated_hits: List of validated marker hits from Step 3
        genome_fasta: Path to genome FASTA (for scaffold lengths)
        proteome_fasta: Path to proteome FASTA (for gene-based distances)
        output_dir: Output directory for marker_seed_regions.bed
        initial_window_bp: Max bp between markers for initial clustering (default: 10kb)
        initial_window_genes: Max genes between markers for initial clustering (default: 5)
        min_markers_initial: Min markers to form initial cluster (default: 2). Cluster
            must contain at least one is_valid_seed_marker.
        extension_kb: Extension distance from outermost markers (default: 5kb)
        merge_distance: Max gap to merge overlapping regions (default: 1kb)
        write_outputs: Write the BED and detailed TSV artifacts when true.

    Returns:
        List of CandidateRegion objects
    """
    logger.info("=" * 60)
    logger.info("Step 4: Marker-Driven Region Assembly")
    logger.info("=" * 60)

    if type(write_outputs) is not bool:
        raise TypeError("write_outputs must be a bool")
    if write_outputs:
        output_dir.mkdir(parents=True, exist_ok=True)

    # Load scaffold lengths
    logger.info("Loading scaffold lengths...")
    scaffold_lengths = load_scaffold_lengths(genome_fasta)
    logger.info(f"  Loaded {len(scaffold_lengths)} scaffolds")

    # Load gene order for gene-based distances
    logger.info("Loading gene genomic order...")
    gene_order = load_gene_predictions(proteome_fasta)
    logger.info(f"  Loaded gene order for {len(gene_order)} scaffolds")

    # Normalize validated hits (supports marker_validation.ValidatedMarkerHit)
    normalized_hits: list[ValidatedMarkerHit] = []
    for hit in validated_hits:
        query_porf = getattr(hit, "query_porf", None)
        porf_id = query_porf or ""
        top10_prefixes = getattr(hit, "top10_prefixes", "")
        # Ensure it's a string (convert list if necessary)
        if isinstance(top10_prefixes, list):
            top10_prefixes = ",".join(top10_prefixes)

        role = decide_marker_hit_role(
            hit,
            ablation_id=ablation_id,
            single_marker_min_score=single_marker_min_score,
        )
        normalized_hits.append(
            ValidatedMarkerHit(
                query_porf=porf_id or "",
                scaffold=hit.scaffold,
                start=hit.start,
                end=hit.end,
                strand=getattr(hit, "strand", "+"),
                hmm_target=hit.hmm_target,
                hmm_score=hit.hmm_score,
                hmm_evalue=getattr(hit, "hmm_evalue", 0.0),
                validation_status=hit.validation_status,
                top10_prefixes=top10_prefixes,
                best_hit_target=getattr(hit, "best_hit_target", ""),
                best_hit_pident=getattr(hit, "best_hit_pident", 0.0),
                best_hit_bits=getattr(hit, "best_hit_bits", 0.0),
                has_ncldv=getattr(hit, "has_ncldv", 0),
                has_mirus=getattr(hit, "has_mirus", 0),
                has_plv=getattr(hit, "has_plv", 0),
                has_vp=getattr(hit, "has_vp", 0),
                has_viral=getattr(hit, "has_viral", 0),
                top10_targets=getattr(hit, "top10_targets", ""),
                top10_pidents=getattr(hit, "top10_pidents", ""),
                top10_bitscores=getattr(hit, "top10_bitscores", ""),
                top10_evalues=getattr(hit, "top10_evalues", ""),
                taxonomy_substring_counts=getattr(hit, "taxonomy_substring_counts", ""),
                taxonomy_raw_counts=getattr(hit, "taxonomy_raw_counts", ""),
                tier1_bypassed=role.is_tier1_bypassed,
            )
        )

    # Group validated hits by scaffold
    hits_by_scaffold = defaultdict(list)
    for hit in normalized_hits:
        if hit.is_validated:  # Only use validated markers
            hits_by_scaffold[hit.scaffold].append(hit)

    logger.info(f"Validated markers on {len(hits_by_scaffold)} scaffolds")

    candidate_regions = []
    compact_cress_regions = assemble_compact_cress_regions(
        [
            hit
            for scaffold_hits in hits_by_scaffold.values()
            for hit in scaffold_hits
        ],
        gene_order,
    )
    cress_gene_keys = {
        (marker.scaffold, base_marker_gene_id(marker.query_porf))
        for region in compact_cress_regions
        for marker in region.markers
    }
    total_clusters = 0
    extended_clusters = 0
    extension_added_markers = 0

    # Process each scaffold
    for scaffold, scaffold_hits in sorted(hits_by_scaffold.items()):
        scaffold_hits = [
            hit
            for hit in scaffold_hits
            if (
                hit.scaffold,
                base_marker_gene_id(hit.query_porf),
            )
            not in cress_gene_keys
        ]
        if not scaffold_hits:
            continue
        logger.info(f"Processing scaffold: {scaffold} ({len(scaffold_hits)} markers)")

        # Step 4a: Initial clustering
        clusters = initial_clustering(
            validated_hits=scaffold_hits,
            scaffold=scaffold,
        initial_window_bp=initial_window_bp,
        initial_window_genes=initial_window_genes,
        min_markers_initial=min_markers_initial,
        gene_order=gene_order,
    )

        cluster_marker_ids = {m.query_porf for cluster in clusters for m in cluster}
        single_marker_clusters = []
        for marker in scaffold_hits:
            if marker.query_porf in cluster_marker_ids:
                continue
            # Accept any valid seed marker with NCLDV/MIRUS/viral taxonomy
            if marker.is_valid_seed_marker:
                single_marker_clusters.append([marker])

        if single_marker_clusters:
            logger.info(
                "  Added %d single-marker clusters (valid seed markers with NCLDV/MIRUS/VP/PLV)",
                len(single_marker_clusters),
            )

        clusters.extend(single_marker_clusters)

        if not clusters:
            logger.info(f"  No clusters formed (need ≥{min_markers_initial} markers per cluster)")
            continue

        # Step 4b: Iterative extension
        scaffold_length = scaffold_lengths.get(scaffold, float("inf"))
        logger.info(f"  Iteratively extending {len(clusters)} clusters...")

        for cluster in clusters:
            total_clusters += 1
            region = iterative_extension(
                cluster=cluster,
                scaffold=scaffold,
                all_validated_hits=scaffold_hits,
                extension_kb=extension_kb,
                scaffold_length=scaffold_length,
            )
            added_markers = region.marker_count - len(cluster)
            if added_markers > 0:
                extended_clusters += 1
                extension_added_markers += added_markers
            candidate_regions.append(region)

    # Step 4c: Merge overlapping regions
    logger.info("Merging overlapping regions...")
    pre_merge_count = len(candidate_regions)
    candidate_regions = merge_overlapping_regions(candidate_regions, merge_distance=merge_distance)
    merged_count = pre_merge_count - len(candidate_regions)
    if merged_count > 0:
        logger.info("  Merged %d overlapping regions", merged_count)
    candidate_regions.extend(compact_cress_regions)
    candidate_regions.sort(key=lambda region: (region.scaffold, region.start, region.end))
    if compact_cress_regions:
        logger.info(
            "Added %d compact, gene-bounded CRESS regions",
            len(compact_cress_regions),
        )

    # Assign region IDs
    for idx, region in enumerate(candidate_regions, start=1):
        region.region_id = f"marker_seed_{idx}"

    if write_outputs:
        # Write output BED file
        bed_path = output_dir / "marker_seed_regions.bed"
        with bed_path.open("w") as handle:
            handle.write("# scaffold\tstart\tend\tregion_id\tmarker_count\tmarker_types\thas_mcp\n")
            for region in candidate_regions:
                handle.write(
                    f"{region.scaffold}\t{region.start}\t{region.end}\t{region.region_id}\t"
                    f"{region.marker_count}\t{region.marker_types_str}\t{1 if region.has_mcp else 0}\n"
                )

        logger.info(f"Wrote marker seed regions: {bed_path}")

        # Write detailed TSV for marker-supported regions
        tsv_path = output_dir / "marker_seed_regions.tsv"
        with tsv_path.open("w") as handle:
            handle.write(
                "region_id\tscaffold\tstart\tend\tlength_bp\tmarker_count\t"
                "ncldv_top10_count\tmirus_top10_count\tplv_top10_count\t"
                "vp_top10_count\thas_mcp\tmarker_targets\tmarker_query_ids\n"
            )
            for region in candidate_regions:
                ncldv_top10 = sum(1 for m in region.markers if getattr(m, "has_ncldv", 0))
                mirus_top10 = sum(1 for m in region.markers if getattr(m, "has_mirus", 0))
                plv_top10 = sum(1 for m in region.markers if getattr(m, "has_plv", 0))
                vp_top10 = sum(1 for m in region.markers if getattr(m, "has_vp", 0))
                marker_targets = ",".join(sorted(m.hmm_target for m in region.markers))
                marker_query_ids = ",".join(sorted(m.query_porf for m in region.markers))
                handle.write(
                    f"{region.region_id}\t{region.scaffold}\t{region.start}\t{region.end}\t"
                    f"{region.length}\t{region.marker_count}\t{ncldv_top10}\t{mirus_top10}\t"
                    f"{plv_top10}\t{vp_top10}\t{1 if region.has_mcp else 0}\t"
                    f"{marker_targets}\t{marker_query_ids}\n"
                )

        logger.info(f"Wrote marker seed detail TSV: {tsv_path}")

    # Log statistics
    if candidate_regions:
        total_length = sum(r.length for r in candidate_regions)
        total_markers = sum(r.marker_count for r in candidate_regions)
        unique_marker_types = set()
        for r in candidate_regions:
            unique_marker_types.update(r.marker_types)
        mcp_regions = sum(1 for r in candidate_regions if r.has_mcp)

        logger.info("=" * 60)
        logger.info("Region Assembly Statistics:")
        logger.info(f"  Initial clusters: {total_clusters}")
        logger.info(
            f"  Extended clusters (±{extension_kb}kb): {extended_clusters} "
            f"(added markers: {extension_added_markers})"
        )
        logger.info(f"  Total regions: {len(candidate_regions)}")
        logger.info(f"  Total coverage: {total_length:,} bp")
        logger.info(f"  Total markers: {total_markers}")
        logger.info(f"  Unique marker types: {len(unique_marker_types)}")
        logger.info(f"  Regions with MCP: {mcp_regions} ({100*mcp_regions/len(candidate_regions):.1f}%)")
        logger.info(f"  Mean markers per region: {total_markers/len(candidate_regions):.1f}")
        logger.info(f"  Mean region length: {total_length/len(candidate_regions):,.0f} bp")
        logger.info("=" * 60)

    return candidate_regions
