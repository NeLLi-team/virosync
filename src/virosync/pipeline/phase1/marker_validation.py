"""
HMM-Gated Diamond Marker Validation.

This module implements the key optimization from PIPELINE_HMM_GATED_PLAN.md:
Only run Diamond on HMM-hit pORFs (not all pORFs), validating viral markers
with top-10 taxonomy analysis using the SMALL marker.dmnd database.

=== ViroSync's Two-Tier Diamond Strategy ===

TIER 1 (THIS MODULE): marker.dmnd validation (Phase 1, Steps 2-3)
  - Database: marker.dmnd (~1,500 curated marker sequences)
  - Input: Only HMM-hit proteins (1-10% of proteome)
  - Purpose: Validate markers (NCLDV/MIRUS/VP/PLV/CRESS top-10 hit at >=25% identity,
             or HMM-only novel marker gates)
  - Output: Validated markers that become EVE region seeds
  - Speed: FAST (small database, small input)

TIER 2 (Phase 2b/3): combined_proteome.dmnd taxonomy
  - Database: combined_proteome.dmnd (millions of sequences)
  - Input: All genes in/near EVE candidate regions
  - Purpose: Full taxonomic classification for boundary refinement
  - Speed: SLOWER but comprehensive (only run on narrowed regions)

Key optimization steps:
1. Run HMM search first (fast, targeted)
2. Extract only HMM-hit pORF sequences (typically 1-10% of total)
3. Run Diamond vs marker.dmnd (TIER 1 - small, fast)
4. Validate markers based on top-10 hit taxonomy
5. Assemble regions from validated markers, extend ±5kb, merge adjacent
6. Run Diamond vs combined_proteome.dmnd ONLY on candidate regions (TIER 2)

This reduces Phase 1 runtime from hours to minutes on large genomes.

Taxonomy prefixes (from marker.dmnd database):
- EUK__: Eukaryote (cellular)
- BAC__: Bacteria (cellular)
- ARC__: Archaea (cellular)
- NCLDV__: Giant viruses (Nucleocytoviricota) - VALIDATES marker
- MIRUS__: Mirusviricota - VALIDATES marker
- VP__: Virophages - VALIDATES marker
- PLV__: Polinton-like viruses - VALIDATES marker
- CRESS__: CRESS/ssDNA viruses - VALIDATES marker
- GVMAG__: Giant virus MAGs - SUPPORTS marker
- PHAGE__: Bacteriophages - viral top-10 evidence

Validation rules in the production filter:
- "validated": at least one validated-prefix top-10 hit with >=25% identity
- "validated_novel": no Diamond match, but HMM-only gates pass
- "supported": GVMAG__ top-10 support below the validation threshold
- "unvalidated": cellular only, unknown, or below thresholds
"""

import logging
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional

from pyhmmer.easel import SequenceFile
from Bio import SeqIO

from virosync.pipeline.phase0.prodigal import parse_prodigal_header
from virosync.pipeline.phase1.coordinates import parse_frame_id, aa_to_nt_coords
from virosync.pipeline.taxonomy_utils import aggregate_taxonomy_substrings, resolve_org_id

logger = logging.getLogger(__name__)


# Taxonomy prefix definitions
VALIDATED_PREFIXES = {"NCLDV__", "MIRUS__", "PLV__", "VP__", "PPV__", "CRESS__"}
SUPPORTING_PREFIXES = {"GVMAG__", "PHAGE__"}
CELLULAR_PREFIXES = {"EUK__", "BAC__", "ARC__"}
VALIDATION_MIN_PIDENT = 25.0


class ValidationStatus(Enum):
    """Marker validation status.

    Attributes:
        VALIDATED: At least one validated-prefix top-10 hit with >=25% identity
        VALIDATED_NOVEL: HMM hit with no Diamond match, but meets gating criteria
        SUPPORTED: GVMAG in top-10 below the validation threshold
        UNVALIDATED: Cellular only or below threshold
    """

    VALIDATED = "validated"
    VALIDATED_NOVEL = "validated_novel"
    SUPPORTED = "supported"
    UNVALIDATED = "unvalidated"


@dataclass
class NovelMarkerCriteria:
    """Criteria for accepting HMM-only markers as validated.

    When an HMM hit has no Diamond matches, these criteria determine
    whether it should be accepted as a novel viral marker.

    Attributes:
        min_hmm_score: Minimum HMM bit score required
        min_hmm_coverage: Minimum query coverage (0.0-1.0)
        require_cluster: If True, marker must be near other markers
        max_confidence_weight: Maximum contribution to confidence scoring
    """

    min_hmm_score: float = 30.0
    min_hmm_coverage: float = 0.5
    require_cluster: bool = True
    max_confidence_weight: float = 0.7




@dataclass
class ValidatedMarkerHit:
    """
    Represents a validated HMM marker hit with Diamond taxonomy support.

    Attributes:
        query_porf: pORF ID
        scaffold: Scaffold name
        start: pORF start position (bp)
        end: pORF end position (bp)
        strand: Strand (+ or -)
        hmm_target: HMM profile name (e.g., GVOGm0003, mcp, polb)
        hmm_score: HMM bit score
        hmm_evalue: HMM E-value
        validation_status: "validated" | "validated_novel" | "supported" | "unvalidated"
        top10_prefixes: Comma-separated list of taxonomy prefixes from top-10 hits
        best_hit_target: Best Diamond hit target ID
        best_hit_pident: Percent identity of best hit
        best_hit_bits: Bit score of best hit
        has_ncldv: 1 if NCLDV__ in top-10, else 0
        has_mirus: 1 if MIRUS__ in top-10, else 0
        has_viral: 1 if any viral prefix in top-10, else 0
        top10_targets: All 10 target IDs (comma-separated)
        top10_pidents: All 10 pident values (comma-separated)
        top10_bitscores: All 10 bitscores (comma-separated)
        top10_evalues: All 10 evalues (comma-separated)
        taxonomy_substring_counts: Weighted taxonomy counts (e.g., "Tubulinea:7.2,Amoebozoa:6.8")
        taxonomy_raw_counts: Raw (unweighted) occurrence counts (e.g., "Tubulinea:8,Amoebozoa:7")
    """
    query_porf: str
    scaffold: str
    start: int
    end: int
    strand: str
    hmm_target: str
    hmm_score: float
    hmm_evalue: float
    validation_status: str
    top10_prefixes: str
    best_hit_target: str
    best_hit_pident: float
    best_hit_bits: float
    has_ncldv: int
    has_mirus: int
    has_plv: int
    has_vp: int
    has_viral: int
    top10_targets: str = ""
    top10_pidents: str = ""
    top10_bitscores: str = ""
    top10_evalues: str = ""
    taxonomy_substring_counts: str = ""
    taxonomy_raw_counts: str = ""

    @property
    def is_validated(self) -> bool:
        """Check if marker is validated by top-10 viral support or validated_novel (HMM-only)."""
        return self.validation_status in ("validated", "validated_novel")

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
    def is_valid_seed_marker(self) -> bool:
        """Check if marker is a valid seed for region assembly.

        UPDATED (Jan 2026): Expanded to allow ALL validated markers to seed regions,
        with taxonomy-based expansion filtering false positives in Step 4.5.

        Seeding rules:
        1. validation_status="validated" (NCLDV/MIRUS/PLV/VP/CRESS Diamond hit >=25% identity) → CAN SEED
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
        if self.validation_status == "validated":
            return True

        # validated_novel (HMM-only, no Diamond hits) → ONLY MCP can seed
        if self.validation_status == "validated_novel":
            return self.is_mcp  # Use expanded MCP detection

        # unvalidated or supported → cannot seed
        return False

    @property
    def is_gvogm(self) -> bool:
        """Check if marker is GVOGm."""
        return self.hmm_target.lower().startswith("gvogm")


def validate_hmm_hit(
    hmm_score: float,
    hmm_coverage: float,
    diamond_hits: list[tuple[str, float, float, float]],
    novel_criteria: NovelMarkerCriteria,
    has_nearby_markers: bool = False,
) -> ValidationStatus:
    """
    Validate an HMM hit using Diamond taxonomy and gating criteria.

    If Diamond hits exist, validation is based on taxonomy prefixes.
    If no Diamond hits exist, applies gating criteria to determine
    if the marker should be accepted as validated_novel.

    Args:
        hmm_score: HMM bit score for this hit
        hmm_coverage: Query coverage from HMM alignment (0.0-1.0)
        diamond_hits: List of (target, bits, pident, evalue) tuples from Diamond
        novel_criteria: Criteria for accepting HMM-only markers
        has_nearby_markers: Whether other validated markers are nearby

    Returns:
        ValidationStatus enum indicating the validation result
    """
    if not diamond_hits:
        # No Diamond hits - check gating criteria for novel markers
        meets_score = hmm_score >= novel_criteria.min_hmm_score
        meets_coverage = hmm_coverage >= novel_criteria.min_hmm_coverage
        meets_cluster = not novel_criteria.require_cluster or has_nearby_markers

        if meets_score and meets_coverage and meets_cluster:
            return ValidationStatus.VALIDATED_NOVEL
        else:
            return ValidationStatus.UNVALIDATED

    # Diamond hits exist - classify based on taxonomy prefixes.
    top10_prefixes = [parse_taxonomy_prefix(target) for target, _, _, _ in diamond_hits[:10]]

    max_validated_pident = max(
        (pident for (_, _, pident, _), p in zip(diamond_hits[:10], top10_prefixes) if p in VALIDATED_PREFIXES),
        default=0.0,
    )
    if max_validated_pident >= VALIDATION_MIN_PIDENT:
        return ValidationStatus.VALIDATED

    # Check for supporting prefixes below the validation threshold.
    if any(p == "GVMAG__" for p in top10_prefixes):
        return ValidationStatus.SUPPORTED

    # All cellular or unknown
    return ValidationStatus.UNVALIDATED


def extract_hmm_hit_sequences(
    hmm_hits: list,
    proteome_fasta: Path,
    output_fasta: Path,
) -> int:
    """
    Extract full protein sequences for HMM-hit pORFs to a FASTA file.

    Each unique protein is written once (deduplicated by base pORF name).
    Diamond searches the full protein rather than just the HMM-aligned
    sub-region, which improves sensitivity for divergent markers.

    Args:
        hmm_hits: List of HMMHit objects from hhg_seeding.py
        proteome_fasta: Path to full conceptual proteome FASTA
        output_fasta: Path to write extracted sequences

    Returns:
        Number of sequences extracted
    """
    # Get unique pORF IDs from HMM hits
    porf_ids = {hit.query_name for hit in hmm_hits}

    if not porf_ids:
        logger.warning("No HMM hits provided for sequence extraction")
        return 0

    # Load sequences for the HMM-hit IDs
    sequences = {}
    with SequenceFile(proteome_fasta) as seq_handle:
        for record in seq_handle:
            name = record.name.decode() if isinstance(record.name, bytes) else record.name
            if name in porf_ids:
                sequence = record.sequence.decode() if isinstance(record.sequence, bytes) else str(record.sequence)
                sequences[name] = sequence

    # Write full protein sequences, deduplicated by base pORF name.
    # Each protein is written once even if multiple HMM models hit it.
    extracted = 0
    total_len = 0
    written_porfs: set[str] = set()
    with output_fasta.open("w") as out_handle:
        for hit in hmm_hits:
            base_name = hit.query_name
            if base_name in written_porfs:
                continue
            seq = sequences.get(base_name)
            if not seq:
                continue
            out_handle.write(f">{base_name}\n{seq}\n")
            written_porfs.add(base_name)
            extracted += 1
            total_len += len(seq)

    avg_len = (total_len / extracted) if extracted else 0.0
    logger.info(
        "Extracted %d full-length HMM-hit proteins (from %d HMM hits, %d unique pORFs); avg_len=%.1f aa",
        extracted,
        len(hmm_hits),
        len(porf_ids),
        avg_len,
    )
    return extracted


def run_diamond_on_hmm_hits(
    hmm_hit_fasta: Path,
    diamond_db: Path,
    output_tsv: Path,
    threads: int = 4,
    evalue: float = 1e-5,
    max_seqs: int = 10,
    search_backend: str = "diamond",
    sensitive: bool = False,
) -> Path:
    """
    Run sequence search on HMM-hit pORFs only (NOT all pORFs).

    This is the core optimization: instead of searching millions of pORFs,
    we only search the ~1-10% that have HMM hits.

    Args:
        hmm_hit_fasta: FASTA with HMM-hit pORFs only
        diamond_db: Diamond database path (.dmnd)
        output_tsv: Output TSV file path
        threads: Number of CPU threads
        evalue: E-value cutoff
        max_seqs: Maximum hits per query (default 10 for top-10 analysis)
        search_backend: "diamond" (the only supported backend)
        sensitive: Use Diamond --sensitive mode for better detection of
            divergent homologs (slower but more sensitive)

    Returns:
        Path to output TSV file
    """
    from virosync.pipeline.search_backend import run_sequence_search

    mode_str = " (sensitive)" if sensitive else ""
    logger.info(
        "Running %s%s on HMM-hit pORFs (max %d hits per query)",
        search_backend, mode_str, max_seqs,
    )
    logger.info("  Query: %s", hmm_hit_fasta.name)
    logger.info("  Database: %s", diamond_db.name)

    extra_flags = ["--sensitive"] if sensitive else None
    run_sequence_search(
        query_fasta=hmm_hit_fasta,
        db_path=diamond_db,
        output_tsv=output_tsv,
        threads=threads,
        backend=search_backend,
        evalue=evalue,
        max_target_seqs=max_seqs,
        extra_flags=extra_flags,
    )

    logger.info(f"Diamond search complete: {output_tsv}")
    return output_tsv


def collect_host_signatures(
    marker_hits: list,
    host_prefixes: Optional[set[str]] = None,
) -> set[str]:
    """
    Collect host-signature tokens from unvalidated marker hits.
    """
    prefixes = host_prefixes or {"EUK__"}
    signatures: set[str] = set()
    for hit in marker_hits:
        status = hit.get("validation_status") if isinstance(hit, dict) else getattr(hit, "validation_status", "")
        if status != "unvalidated":
            continue
        target = hit.get("best_hit_target") if isinstance(hit, dict) else getattr(hit, "best_hit_target", "")
        if not target:
            continue
        if any(str(target).startswith(p) for p in prefixes):
            signatures.add(str(target).split("|", 1)[0])
    return signatures


def write_extended_taxonomy(
    diamond_top10_tsv: Path,
    output_dir: Path,
    tax_lookup: dict,
) -> Path:
    """
    Write diamond_top10_taxonomy.tsv with full lineage labels appended.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "diamond_top10_taxonomy.tsv"

    with diamond_top10_tsv.open() as inp, output_path.open("w") as out:
        out.write(
            "query\ttarget\tevalue\tbitscore\tpident\tqcov\ttaxonomy_label\n"
        )
        for line in inp:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            target = parts[1]
            token = resolve_org_id(target.split("|", 1)[0], tax_lookup)
            label = tax_lookup.get(token, "")
            out.write("\t".join(parts + [label]) + "\n")

    logger.info("Wrote extended taxonomy: %s", output_path)
    return output_path


def parse_taxonomy_prefix(target_id: str) -> str:
    """
    Extract taxonomy prefix from target sequence ID.

    Args:
        target_id: Target sequence ID (e.g., "NCLDV__GVOGm0003_protein123")

    Returns:
        Taxonomy prefix (e.g., "NCLDV__") or "UNKNOWN__" if not found
    """
    for prefix in VALIDATED_PREFIXES | SUPPORTING_PREFIXES | CELLULAR_PREFIXES:
        if target_id.startswith(prefix):
            return prefix
    return "UNKNOWN__"




def filter_validated_markers(
    hmm_hits: list,
    diamond_output_tsv: Path,
    proteome_fasta: Path,
    genome_fasta: Optional[Path] = None,
    output_dir: Optional[Path] = None,
    novel_criteria: Optional[NovelMarkerCriteria] = None,
    taxonomy_lookup: Optional[dict] = None,
    taxonomy_weight_mode: str = "rank",
    max_seqs: int = 10,
) -> list[ValidatedMarkerHit]:
    """
    Filter HMM hits to validated markers based on Diamond top-10 taxonomy.

    Validation rules:
    1. "validated": at least one top-10 hit from a validated viral prefix
       (NCLDV__, MIRUS__, VP__, PLV__, or CRESS__) with >=25% identity
    2. "validated_novel": HMM hit with no Diamond match, meets gating criteria
    3. "supported": GVMAG__ top-10 support below the validation threshold
    4. "unvalidated": cellular only, unknown, or below thresholds

    Args:
        hmm_hits: List of HMMHit objects from hhg_seeding.py
        diamond_output_tsv: Diamond search results TSV
        proteome_fasta: Path to conceptual proteome (for coordinate parsing)
        genome_fasta: Optional path to genome FASTA for coordinate resolution
        output_dir: Optional directory to write output files
        novel_criteria: Criteria for accepting HMM-only markers as validated_novel
        max_seqs: Number of ranked Diamond hits retained per marker

    Returns:
        List of ValidatedMarkerHit objects (includes all hits with metadata)
    """
    def _marker_group(name: str) -> str:
        lower = name.lower()
        if lower.startswith("gvogm"):
            return "GVOGm"
        if lower.startswith("og") or lower.startswith("mog"):
            return "OG"
        return "OTHER"

    # Parse pORF coordinates from headers
    porf_coords = {}
    for record in SeqIO.parse(proteome_fasta, "fasta"):
        name = record.id
        prodigal_parsed = parse_prodigal_header(record.description, name)
        if prodigal_parsed:
            scaffold, start, end, strand = prodigal_parsed
            porf_coords[name] = (scaffold, start, end, strand)

    contig_lengths = {}
    if genome_fasta:
        with SequenceFile(genome_fasta) as seq_handle:
            for record in seq_handle:
                name = record.name.decode() if isinstance(record.name, bytes) else record.name
                seq = record.sequence.decode() if isinstance(record.sequence, bytes) else str(record.sequence)
                contig_lengths[name] = len(seq)

    # Parse Diamond results
    diamond_hits = {}  # query -> list of (rank, target, bits, pident, evalue)

    if diamond_output_tsv.exists() and diamond_output_tsv.stat().st_size > 0:
        with diamond_output_tsv.open() as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                query, target, evalue, bits, pident, qcov = parts

                try:
                    bits_val = float(bits)
                    pident_val = float(pident)
                    evalue_val = float(evalue)
                except ValueError:
                    continue

                if query not in diamond_hits:
                    diamond_hits[query] = []
                diamond_hits[query].append((target, bits_val, pident_val, evalue_val))

    # Build HMM hit lookup keyed by query|aaX-Y (unique per HMM hit).
    # Diamond results are keyed by base protein name (full-protein search),
    # so we also build a mapping from base name → Diamond hits.
    hmm_hit_lookup = {}  # query|aaX-Y -> (target_name, score, evalue, qstart, qend)
    for hit in hmm_hits:
        query_id = f"{hit.query_name}|aa{hit.query_start}-{hit.query_end}"
        if query_id not in hmm_hit_lookup:
            hmm_hit_lookup[query_id] = (
                hit.target_name,
                hit.score,
                hit.evalue,
                hit.query_start,
                hit.query_end,
            )
        elif hit.score > hmm_hit_lookup[query_id][1]:
            # Keep best HMM score if multiple hits for identical segment
            hmm_hit_lookup[query_id] = (
                hit.target_name,
                hit.score,
                hit.evalue,
                hit.query_start,
                hit.query_end,
            )

    hmm_group_counts = {"GVOGm": 0, "OG": 0, "OTHER": 0}
    for hmm_target, _score, _evalue, _qstart, _qend in hmm_hit_lookup.values():
        hmm_group_counts[_marker_group(hmm_target)] += 1

    # Use default novel criteria if not provided
    if novel_criteria is None:
        novel_criteria = NovelMarkerCriteria()

    # Validate markers
    validated_markers = []
    stats = {
        "validated": 0,
        "validated_novel": 0,
        "supported": 0,
        "unvalidated": 0,
        "no_diamond_hits": 0,
    }

    for query_porf, (hmm_target, hmm_score, hmm_evalue, qstart, qend) in hmm_hit_lookup.items():
        # Get coordinates
        coords = porf_coords.get(query_porf)
        base_query = query_porf.split("|aa", 1)[0]
        if isinstance(base_query, bytes):
            base_query = base_query.decode()
        if base_query.startswith("b'") and base_query.endswith("'"):
            base_query = base_query[2:-1]
        if not coords:
            coords = porf_coords.get(base_query)
        if not coords and contig_lengths:
            contig, frame, offset = parse_frame_id(base_query)
            contig_len = contig_lengths.get(contig)
            if frame and contig_len:
                start, end, strand = aa_to_nt_coords(qstart, qend, contig_len, frame, offset=offset)
                scaffold = contig
                coords = (scaffold, start, end, strand)
        if not coords:
            logger.warning(f"Could not parse coordinates for {query_porf}")
            continue
        scaffold, start, end, strand = coords

        # Get the configured number of ranked Diamond hits.
        # Diamond queries use base protein name (full-protein search),
        # so look up by base name; fall back to full query_porf for
        # backward compatibility with old-format result files.
        top10_hits = diamond_hits.get(base_query, []) or diamond_hits.get(query_porf, [])
        if not top10_hits:
            # No Diamond hits - check if qualifies as validated_novel
            stats["no_diamond_hits"] += 1

            # Calculate HMM coverage (query coverage from HMM alignment)
            # Coverage = (query_end - query_start + 1) / query_length
            # Note: We approximate coverage using the aligned region proportion
            hmm_aligned_len = qend - qstart + 1
            # Use a reasonable estimate of typical pORF length (300 aa) if unknown
            # This is a conservative estimate; actual coverage may be higher
            hmm_coverage = min(1.0, hmm_aligned_len / 300.0)

            # Check for nearby markers on same scaffold (simplified check)
            # A marker is "nearby" if another HMM hit exists on the same scaffold
            has_nearby_markers = sum(
                1 for other_porf, (_, _, _, _, _) in hmm_hit_lookup.items()
                if other_porf != query_porf
                and other_porf.split("|aa", 1)[0].rsplit("_", 1)[0] == scaffold
            ) > 0

            # Use validate_hmm_hit to determine status
            validation_result = validate_hmm_hit(
                hmm_score=hmm_score,
                hmm_coverage=hmm_coverage,
                diamond_hits=[],
                novel_criteria=novel_criteria,
                has_nearby_markers=has_nearby_markers,
            )
            validation_status = validation_result.value

            if validation_result == ValidationStatus.VALIDATED_NOVEL:
                stats["validated_novel"] += 1
            else:
                stats["unvalidated"] += 1

            validated_markers.append(ValidatedMarkerHit(
                query_porf=query_porf,
                scaffold=scaffold,
                start=start,
                end=end,
                strand=strand,
                hmm_target=hmm_target,
                hmm_score=hmm_score,
                hmm_evalue=hmm_evalue,
                validation_status=validation_status,
                top10_prefixes="",
                best_hit_target="",
                best_hit_pident=0.0,
                best_hit_bits=0.0,
                has_ncldv=0,
                has_mirus=0,
                has_plv=0,
                has_vp=0,
                has_viral=0,
            ))
            continue

        # Sort by bit score (descending) and take the configured top K.
        top10_hits.sort(key=lambda x: x[1], reverse=True)
        top10_hits = top10_hits[:max_seqs]

        # Pad to exactly K hits with empty entries for stable TSV columns.
        while len(top10_hits) < max_seqs:
            top10_hits.append(("", 0.0, 0.0, 1.0))

        # Extract ALL details from top-10 hits
        top10_targets = ",".join(target for target, _, _, _ in top10_hits)
        top10_pidents = ",".join(f"{pident:.1f}" for _, _, pident, _ in top10_hits)
        top10_bitscores = ",".join(f"{bits:.1f}" for _, bits, _, _ in top10_hits)
        top10_evalues = ",".join(f"{evalue:.2e}" for _, _, _, evalue in top10_hits)

        # Extract taxonomy prefixes
        top10_prefixes = [parse_taxonomy_prefix(target) for target, _, _, _ in top10_hits]
        top10_prefixes_str = ",".join(top10_prefixes)

        # Check for validated/supported status
        has_ncldv = int(any(p == "NCLDV__" for p in top10_prefixes))
        has_mirus = int(any(p == "MIRUS__" for p in top10_prefixes))
        # PPV (Preplasmiviricota) is the unified VP+PLV domain; fold PPV__ into the
        # PLV-class signal so it counts toward has_viral and the VP/PLV region support.
        has_plv = int(any(p in ("PLV__", "PPV__") for p in top10_prefixes))
        has_vp = int(any(p == "VP__" for p in top10_prefixes))
        has_cress = any(p == "CRESS__" for p in top10_prefixes)
        has_gvmag = any(p == "GVMAG__" for p in top10_prefixes)
        has_phage = any(p == "PHAGE__" for p in top10_prefixes)
        has_viral = int(has_ncldv or has_mirus or has_plv or has_vp or has_cress or has_gvmag or has_phage)

        max_validated_pident = max(
            (pident for (_, _, pident, _), p in zip(top10_hits, top10_prefixes) if p in VALIDATED_PREFIXES),
            default=0.0,
        )

        # Determine validation status
        # Require one validated viral reference in the top-10 at >=25% identity.
        if max_validated_pident >= VALIDATION_MIN_PIDENT:
            validation_status = "validated"
            stats["validated"] += 1
        elif has_gvmag:
            validation_status = "supported"
            stats["supported"] += 1
        else:
            validation_status = "unvalidated"
            stats["unvalidated"] += 1

        # Get best hit info
        best_target, best_bits, best_pident, best_evalue = top10_hits[0]

        # Aggregate taxonomy substrings for unvalidated cellular hits
        # This creates a taxonomic fingerprint for novel host lineages
        taxonomy_substring_counts = ""
        taxonomy_raw_counts = ""
        if validation_status == "unvalidated":
            # Require at least one known cellular prefix to avoid viral contamination
            cellular_count = sum(1 for p in top10_prefixes if p in CELLULAR_PREFIXES)
            if cellular_count > 0 and taxonomy_lookup:
                fingerprint = aggregate_taxonomy_substrings(
                    top10_hits,
                    taxonomy_lookup,
                    min_token_length=3,
                    weight_mode=taxonomy_weight_mode,
                )
                taxonomy_substring_counts, taxonomy_raw_counts = fingerprint.to_string()

        validated_markers.append(ValidatedMarkerHit(
            query_porf=query_porf,
            scaffold=scaffold,
            start=start,
            end=end,
            strand=strand,
            hmm_target=hmm_target,
            hmm_score=hmm_score,
            hmm_evalue=hmm_evalue,
            validation_status=validation_status,
            top10_prefixes=top10_prefixes_str,
            best_hit_target=best_target,
            best_hit_pident=best_pident,
            best_hit_bits=best_bits,
            has_ncldv=has_ncldv,
            has_mirus=has_mirus,
            has_plv=has_plv,
            has_vp=has_vp,
            has_viral=has_viral,
            top10_targets=top10_targets,
            top10_pidents=top10_pidents,
            top10_bitscores=top10_bitscores,
            top10_evalues=top10_evalues,
            taxonomy_substring_counts=taxonomy_substring_counts,
            taxonomy_raw_counts=taxonomy_raw_counts,
        ))

    logger.info("Marker validation results:")
    logger.info(
        "  HMM hits by group: GVOGm=%s OG=%s OTHER=%s",
        hmm_group_counts["GVOGm"],
        hmm_group_counts["OG"],
        hmm_group_counts["OTHER"],
    )
    logger.info(f"  Validated (validated-prefix top-10 hit >=25% identity): {stats['validated']}")
    logger.info(f"  Validated novel (HMM-only, meets criteria): {stats['validated_novel']}")
    logger.info(f"  Supported (GVMAG in top-10): {stats['supported']}")
    logger.info(f"  Unvalidated (cellular only): {stats['unvalidated']}")
    logger.info(f"  No Diamond hits: {stats['no_diamond_hits']}")

    group_counts = {
        "GVOGm": {"total": 0, "ncldv": 0, "mirus": 0},
        "OG": {"total": 0, "ncldv": 0, "mirus": 0},
        "OTHER": {"total": 0, "ncldv": 0, "mirus": 0},
    }
    taxonomy_counts = {
        "NCLDV__": 0,
        "MIRUS__": 0,
        "PPV__": 0,
        "PLV__": 0,
        "VP__": 0,
        "CRESS__": 0,
        "GVMAG__": 0,
        "PHAGE__": 0,
        "EUK__": 0,
        "BAC__": 0,
        "ARC__": 0,
        "OTHER": 0,
        "NO_HITS": 0,
    }
    mirus_total = 0
    ncldv_total = 0
    for hit in validated_markers:
        group = _marker_group(hit.hmm_target)
        group_counts[group]["total"] += 1
        if hit.has_ncldv:
            group_counts[group]["ncldv"] += 1
            ncldv_total += 1
        if hit.has_mirus:
            group_counts[group]["mirus"] += 1
            mirus_total += 1

        if not hit.top10_prefixes:
            taxonomy_counts["NO_HITS"] += 1
        else:
            top1_prefix = hit.top10_prefixes.split(",")[0]
            if hit.has_ncldv:
                taxonomy_counts["NCLDV__"] += 1
            elif hit.has_mirus:
                taxonomy_counts["MIRUS__"] += 1
            elif "PPV__" in hit.top10_prefixes.split(","):
                taxonomy_counts["PPV__"] += 1
            elif "PLV__" in hit.top10_prefixes.split(","):
                taxonomy_counts["PLV__"] += 1
            elif "VP__" in hit.top10_prefixes.split(","):
                taxonomy_counts["VP__"] += 1
            elif "CRESS__" in hit.top10_prefixes.split(","):
                taxonomy_counts["CRESS__"] += 1
            elif top1_prefix in taxonomy_counts:
                taxonomy_counts[top1_prefix] += 1
            else:
                taxonomy_counts["OTHER"] += 1

    logger.info("  Diamond top-10 by marker group:")
    for group in ("GVOGm", "OG", "OTHER"):
        logger.info(
            "    %s: total=%s ncldv_top10=%s mirus_top10=%s",
            group,
            group_counts[group]["total"],
            group_counts[group]["ncldv"],
            group_counts[group]["mirus"],
        )
    logger.info(
        "  Top-10 viral prefixes: NCLDV=%s MIRUS=%s (of %s total hits)",
        ncldv_total,
        mirus_total,
        len(validated_markers),
    )
    logger.info("  Taxonomy breakdown (top1 overridden by NCLDV/MIRUS in top-10):")
    for prefix in ("NCLDV__", "MIRUS__", "PPV__", "PLV__", "VP__", "CRESS__", "GVMAG__", "PHAGE__", "EUK__", "BAC__", "ARC__", "OTHER", "NO_HITS"):
        logger.info("    %s: %s", prefix, taxonomy_counts[prefix])

    # Write output files
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)

        # Write validated_marker_hits.tsv
        validated_path = output_dir / "validated_marker_hits.tsv"
        with validated_path.open("w") as f:
            f.write("query_porf\tscaffold\tstart\tend\tstrand\thmm_target\thmm_score\t"
                   "validation_status\ttop10_prefixes\tbest_hit_pident\tbest_hit_target\t"
                   "best_hit_bits\thas_ncldv\thas_mirus\thas_viral\t"
                   "top10_targets\ttop10_pidents\ttop10_bitscores\ttop10_evalues\t"
                   "taxonomy_substring_counts\ttaxonomy_raw_counts\n")
            for vm in validated_markers:
                f.write(f"{vm.query_porf}\t{vm.scaffold}\t{vm.start}\t{vm.end}\t{vm.strand}\t"
                       f"{vm.hmm_target}\t{vm.hmm_score:.3f}\t{vm.validation_status}\t"
                       f"{vm.top10_prefixes}\t{vm.best_hit_pident:.1f}\t{vm.best_hit_target}\t"
                       f"{vm.best_hit_bits:.1f}\t{vm.has_ncldv}\t{vm.has_mirus}\t{vm.has_viral}\t"
                       f"{vm.top10_targets}\t{vm.top10_pidents}\t{vm.top10_bitscores}\t"
                       f"{vm.top10_evalues}\t{vm.taxonomy_substring_counts}\t{vm.taxonomy_raw_counts}\n")
        logger.info(f"Wrote validated marker hits: {validated_path}")

        # Write diamond_top10_taxonomy.tsv
        taxonomy_path = output_dir / "diamond_top10_taxonomy.tsv"
        with taxonomy_path.open("w") as f:
            f.write("query_porf\thmm_target\ttop1_prefix\ttop10_prefixes\t"
                   "has_ncldv\thas_mirus\thas_viral\n")
            for vm in validated_markers:
                top1_prefix = vm.top10_prefixes.split(",")[0] if vm.top10_prefixes else ""
                f.write(f"{vm.query_porf}\t{vm.hmm_target}\t{top1_prefix}\t"
                       f"{vm.top10_prefixes}\t{vm.has_ncldv}\t{vm.has_mirus}\t{vm.has_viral}\n")
        logger.info(f"Wrote top-10 taxonomy summary: {taxonomy_path}")

    return validated_markers


def load_validated_marker_hits(tsv_path: Path) -> list[ValidatedMarkerHit]:
    """
    Load validated marker hits from a TSV produced by filter_validated_markers.
    """
    hits: list[ValidatedMarkerHit] = []
    if not tsv_path.exists():
        return hits
    with tsv_path.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            hits.append(
                ValidatedMarkerHit(
                    query_porf=parts[idx.get("query_porf", 0)],
                    scaffold=parts[idx.get("scaffold", 1)],
                    start=int(parts[idx.get("start", 2)]),
                    end=int(parts[idx.get("end", 3)]),
                    strand=parts[idx.get("strand", 4)],
                    hmm_target=parts[idx.get("hmm_target", 5)],
                    hmm_score=float(parts[idx.get("hmm_score", 6)]),
                    hmm_evalue=0.0,
                    validation_status=parts[idx.get("validation_status", 7)],
                    top10_prefixes=parts[idx.get("top10_prefixes", 8)],
                    best_hit_pident=float(parts[idx.get("best_hit_pident", 9)]),
                    best_hit_target=parts[idx.get("best_hit_target", 10)],
                    best_hit_bits=float(parts[idx.get("best_hit_bits", 11)]),
                    has_ncldv=int(parts[idx.get("has_ncldv", 12)]) if "has_ncldv" in idx else int(parts[12]),
                    has_mirus=int(parts[idx.get("has_mirus", 13)]) if "has_mirus" in idx else int(parts[13]),
                    has_plv=int(parts[idx["has_plv"]]) if "has_plv" in idx else 0,
                    has_vp=int(parts[idx["has_vp"]]) if "has_vp" in idx else 0,
                    has_viral=int(parts[idx.get("has_viral", 14)]) if "has_viral" in idx else int(parts[14]),
                    # New fields with backwards compatibility
                    top10_targets=parts[idx["top10_targets"]] if "top10_targets" in idx else "",
                    top10_pidents=parts[idx["top10_pidents"]] if "top10_pidents" in idx else "",
                    top10_bitscores=parts[idx["top10_bitscores"]] if "top10_bitscores" in idx else "",
                    top10_evalues=parts[idx["top10_evalues"]] if "top10_evalues" in idx else "",
                    taxonomy_substring_counts=parts[idx["taxonomy_substring_counts"]] if "taxonomy_substring_counts" in idx else "",
                    taxonomy_raw_counts=parts[idx["taxonomy_raw_counts"]] if "taxonomy_raw_counts" in idx else "",
                )
            )
    return hits
