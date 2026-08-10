"""
Hierarchical Hallmark Graph (HHG) Seeding.

Path A of Phase 1: Identifies EVE seeds based on viral hallmark gene HMM hits.

This module:
1. Runs HMM search using pyhmmer against viral hallmark HMM profiles
2. Identifies "anchor" pORFs with strong hallmark hits
3. Calculates neighbor density scores in genomic windows
4. Forms seed regions from high-scoring anchor neighborhoods

Key concepts:
- Anchor pORF: A pORF with a significant hit to a viral hallmark gene HMM
- Neighbor score: Density of viral-like pORFs in a genomic window around an anchor
- Seed: A genomic region centered on one or more anchors, candidate for EVE expansion

Assembly Modes:
- default: Standard mode for well-assembled genomes (requires marker diversity)
- fragmented: For MAGs/fragmented assemblies (accepts single high-scoring markers)
- relaxed: Exploratory mode (maximum sensitivity)
- strict: High-confidence only (maximum specificity)
"""

import logging
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pyhmmer
from Bio import SeqIO
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMMFile, HMM

from virosync.pipeline.phase0.translation import PORF
from virosync.pipeline.phase0.prodigal import parse_prodigal_header
from virosync.pipeline.phase1.viral_markers import (
    VIRAL_FAMILIES,
    AssemblyMode,
    get_assembly_mode,
    get_family_for_markers,
    ALL_DIAGNOSTIC_MARKERS,
)
from virosync.pipeline.phase1.marker_validation import (
    extract_hmm_hit_sequences,
    run_diamond_on_hmm_hits,
    filter_validated_markers,
    ValidatedMarkerHit,
)
from virosync.pipeline.phase1.coordinates import parse_frame_id, aa_to_nt_coords

logger = logging.getLogger(__name__)


# Legacy definitions for backward compatibility
# These are now superseded by viral_markers.py
VIRUS_SPECIFIC_MARKERS = {"mcp", "a32", "d5", "vltf3", "mrnac", "sfii"}


@dataclass
class HMMHit:
    """Represents a single HMM hit from hmmsearch."""

    query_name: str  # pORF ID
    target_name: str  # HMM profile name (hallmark gene)
    score: float  # Bit score
    evalue: float  # E-value
    domain_score: float  # Best domain score
    query_start: int  # Hit start in query (aa)
    query_end: int  # Hit end in query (aa)


@dataclass
class Anchor:
    """Represents an anchor pORF with strong hallmark gene hit."""

    porf_id: str
    scaffold: str
    start: int  # Nucleotide start
    end: int  # Nucleotide end
    strand: str
    hallmark_gene: str
    score: float
    evalue: float


@dataclass
class HHGSeed:
    """Represents a seed region from HHG seeding."""

    scaffold: str
    start: int
    end: int
    anchors: list[Anchor] = field(default_factory=list)
    neighbor_score: float = 0.0
    source: str = "hhg"
    predicted_family: Optional[str] = None  # Predicted viral family (ncldv, mriyavirus, etc.)

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def num_anchors(self) -> int:
        return len(self.anchors)

    @property
    def hallmark_genes(self) -> list[str]:
        """List of unique hallmark genes in this seed."""
        return list(set(a.hallmark_gene for a in self.anchors))


def load_hmm_profiles(
    hmm_file: Path,
    allowlist: Optional[set[str]] = None,
) -> list[HMM]:
    """
    Load HMM profiles from a file.

    Args:
        hmm_file: Path to HMM profile file (can contain multiple profiles)

    Returns:
        List of pyhmmer HMM objects
    """
    profiles = []
    with HMMFile(hmm_file) as hmm_handle:
        for hmm in hmm_handle:
            hmm_name = hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
            if allowlist and hmm_name not in allowlist:
                continue
            profiles.append(hmm)

    logger.info(f"Loaded {len(profiles)} HMM profiles from {hmm_file}")
    return profiles


def load_hmm_allowlist(allowlist_path: Optional[Path]) -> Optional[set[str]]:
    """Load allowed HMM names (one per line)."""
    if not allowlist_path:
        return None
    if not allowlist_path.exists():
        logger.warning("HMM allowlist not found: %s", allowlist_path)
        return None
    allowlist = set()
    with allowlist_path.open() as handle:
        for line in handle:
            name = line.strip()
            if name and not name.startswith("#"):
                allowlist.add(name)
    if not allowlist:
        logger.warning("HMM allowlist is empty: %s", allowlist_path)
    return allowlist if allowlist else None


def load_faa_markers(marker_faa_dir: Optional[Path]) -> Optional[set[str]]:
    """Load marker names from FAA filenames (stem)."""
    if not marker_faa_dir:
        return None
    marker_faa_dir = Path(marker_faa_dir)
    if not marker_faa_dir.exists():
        logger.warning("Marker FAA dir not found: %s", marker_faa_dir)
        return None
    markers = {p.stem for p in marker_faa_dir.glob("*.faa")}
    if not markers:
        logger.warning("No FAA markers found in: %s", marker_faa_dir)
        return None
    return markers


def _load_porf_sequences(
    proteome_fasta: Path,
    porf_ids: set[str],
) -> dict[str, str]:
    """Load only sequences needed for hit validation."""
    sequences = {}
    with SequenceFile(proteome_fasta) as handle:
        for record in handle:
            name = record.name.decode() if isinstance(record.name, bytes) else record.name
            if name in porf_ids:
                sequences[name] = str(record.sequence)
    return sequences


def _load_porf_order(
    proteome_fasta: Path,
) -> dict[str, list[str]]:
    """Load pORF order per scaffold for neighbor checks."""
    porfs_by_scaffold: dict[str, list[tuple[int, str]]] = defaultdict(list)
    with SequenceFile(proteome_fasta) as handle:
        for record in handle:
            name = record.name.decode() if isinstance(record.name, bytes) else record.name
            porf = PORF.parse_header(name)
            if not porf:
                continue
            porfs_by_scaffold[porf.scaffold].append((porf.start, name))
    ordered = {}
    for scaffold, porfs in porfs_by_scaffold.items():
        ordered[scaffold] = [p_id for _, p_id in sorted(porfs)]
    return ordered


def _run_diamond_marker_search(
    query_fasta: Path,
    target_faa: Path,
    output_file: Path,
    threads: int,
) -> None:
    tmp_dir = output_file.parent / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = tmp_dir / target_faa.stem
    db_path = db_prefix.with_suffix(".dmnd")
    # Build the small per-marker Diamond DB once (idempotent). 5-min timeout
    # is generous — these are tiny inputs.
    if not db_path.exists():
        subprocess.run(
            ["diamond", "makedb", "--in", str(target_faa), "-d", str(db_prefix)],
            check=True,
            capture_output=True,
            text=True,
            timeout=300,
            start_new_session=True,
        )

    # Route the actual blastp through the central hardened wrapper so this
    # call benefits from the per-call tempdir + retry that prevents the
    # Diamond 2.1.21 finalization-phase futex deadlock.
    from virosync.pipeline.search_backend import run_sequence_search

    run_sequence_search(
        query_fasta=query_fasta,
        db_path=db_path,
        output_tsv=output_file,
        threads=threads,
        backend="diamond",
        evalue=1e-5,
        max_target_seqs=10,
        output_columns=["qseqid", "sseqid", "bitscore"],
        timeout=600,
    )


def validate_hmm_hits_with_markers(
    hits: list[HMMHit],
    proteome_fasta: Path,
    marker_faa_dir: Path,
    threads: int = 4,
    top_k: int = 10,
    neighbor_genes: int = 5,
    output_dir: Optional[Path] = None,
) -> list[HMMHit]:
    """
    Validate HMM hits with marker-specific FAA databases.

    For each HMM hit:
      1) BLAST (Diamond) query pORFs against corresponding marker FAA.
      2) Accept if NCLDV, MIRUS, VP, or PLV appears in top-k hits.
      3) Keep only hits that also have another validated marker within +/- neighbor_genes.
    """
    if not hits:
        return []

    marker_faa_dir = Path(marker_faa_dir)
    if not marker_faa_dir.exists():
        logger.warning("Marker FAA dir not found: %s", marker_faa_dir)
        return hits

    porf_ids = {h.query_name for h in hits}
    porf_sequences = _load_porf_sequences(proteome_fasta, porf_ids)
    if not porf_sequences:
        logger.warning("No pORF sequences loaded for marker validation")
        return hits

    hits_by_marker: dict[str, list[HMMHit]] = defaultdict(list)
    for hit in hits:
        hits_by_marker[hit.target_name].append(hit)

    validated: set[str] = set()
    validated_marker: dict[str, str] = {}

    diamond_top_hits: list[tuple[str, str, int, str, float, str]] = []
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        for marker, marker_hits in hits_by_marker.items():
            marker_faa = marker_faa_dir / f"{marker}.faa"
            if not marker_faa.exists():
                logger.warning("Marker FAA not found for %s: %s", marker, marker_faa)
                continue

            query_fasta = tmpdir_path / f"{marker}.queries.faa"
            output_file = tmpdir_path / f"{marker}.diamond.tsv"

            with query_fasta.open("w") as handle:
                for hit in marker_hits:
                    seq = porf_sequences.get(hit.query_name)
                    if not seq:
                        continue
                    handle.write(f">{hit.query_name}\n{seq}\n")

            if query_fasta.stat().st_size == 0:
                continue

            try:
                _run_diamond_marker_search(
                    query_fasta=query_fasta,
                    target_faa=marker_faa,
                    output_file=output_file,
                    threads=threads,
                )
            except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                logger.warning("Diamond failed for %s: %s", marker, exc)
                continue

            if not output_file.exists() or output_file.stat().st_size == 0:
                continue

            top_hits: dict[str, list[tuple[float, str]]] = defaultdict(list)
            with output_file.open() as handle:
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 3:
                        continue
                    query, target, bits = parts[0], parts[1], parts[2]
                    try:
                        bits_val = float(bits)
                    except ValueError:
                        bits_val = 0.0
                    top_hits[query].append((bits_val, target))

            for query, hits_list in top_hits.items():
                hits_list.sort(key=lambda x: x[0], reverse=True)
                top_targets = [t for _, t in hits_list[:top_k]]
                for rank, (bits_val, target) in enumerate(hits_list[:top_k], start=1):
                    is_validated = any(
                        target.startswith(p)
                        for p in ("NCLDV__", "MIRUS__", "VP__", "PLV__", "PPV__", "CRESS__")
                    )
                    diamond_top_hits.append(
                        (
                            marker,
                            query,
                            rank,
                            target,
                            bits_val,
                            "1" if is_validated else "0",
                        )
                    )
                if any(
                    t.startswith(p)
                    for t in top_targets
                    for p in ("NCLDV__", "MIRUS__", "VP__", "PLV__", "PPV__", "CRESS__")
                ):
                    validated.add(query)
                    validated_marker[query] = marker

    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        top10_path = output_dir / "diamond_top10.tsv"
        with top10_path.open("w") as handle:
            handle.write("marker\tquery\trank\ttarget\tbits\tis_ncldv_or_mirus\n")
            for row in diamond_top_hits:
                handle.write("\t".join(str(x) for x in row) + "\n")
        logger.info("Wrote Diamond top-10 hits: %s", top10_path)

    if not validated:
        logger.info("No HMM hits validated by marker databases")
        return []

    # Neighbor support: another validated marker within +/- neighbor_genes
    porf_order = _load_porf_order(proteome_fasta)
    porf_index = {
        scaffold: {porf_id: idx for idx, porf_id in enumerate(ids)}
        for scaffold, ids in porf_order.items()
    }
    neighbor_supported: set[str] = set()
    for hit in hits:
        if hit.query_name not in validated:
            continue
        coords = parse_porf_coordinates(hit.query_name)
        if not coords:
            continue
        scaffold = coords[0]
        idx = porf_index.get(scaffold, {}).get(hit.query_name)
        if idx is None:
            continue
        start_idx = max(0, idx - neighbor_genes)
        end_idx = min(len(porf_order[scaffold]) - 1, idx + neighbor_genes)
        for neighbor_id in porf_order[scaffold][start_idx : end_idx + 1]:
            if neighbor_id == hit.query_name:
                continue
            if neighbor_id in validated and validated_marker.get(neighbor_id) != hit.target_name:
                neighbor_supported.add(hit.query_name)
                neighbor_supported.add(neighbor_id)

    if neighbor_supported:
        logger.info(
            "Marker validation: %d hits validated, %d had neighbor support; keeping all validated hits",
            len(validated),
            len(neighbor_supported),
        )
    else:
        logger.info(
            "Marker validation: %d hits validated, none had neighbor support; keeping validated hits",
            len(validated),
        )
    return [h for h in hits if h.query_name in validated]


def validate_hmm_hits_with_combined_db(
    hits: list[HMMHit],
    proteome_fasta: Path,
    marker_db: Path,
    threads: int = 4,
    top_k: int = 10,
    evalue: float = 1e-5,
    output_dir: Optional[Path] = None,
    genome_fasta: Optional[Path] = None,
) -> tuple[list[HMMHit], list[ValidatedMarkerHit]]:
    """
    Validate HMM hits using HMM-gated Diamond workflow (PIPELINE_HMM_GATED_PLAN.md Steps 2-3).

    This is the KEY OPTIMIZATION: Only run Diamond on HMM-hit pORFs, not all pORFs.
    This reduces runtime from hours to minutes on large genomes.

    Workflow:
    1. Extract sequences for HMM-hit pORFs only (1-10% of total pORFs)
    2. Run Diamond blastp on these sequences against combined FAA DB
    3. Parse top-10 hits and validate markers based on taxonomy prefixes
    4. Return validated hits and detailed marker metadata

    Args:
        hits: List of HMMHit objects from HMM search
        proteome_fasta: Path to conceptual proteome FASTA
        marker_db: Path to marker Diamond database (marker.dmnd - TIER 1, ~1,500 curated marker sequences)
        threads: Number of CPU threads
        top_k: Number of top hits to analyze (default 10)
        evalue: E-value cutoff for Diamond
        output_dir: Optional directory for output files

    Returns:
        Tuple of:
        - List of validated HMMHit objects (only "validated" and "supported" markers)
        - List of ValidatedMarkerHit objects with full metadata
    """
    if not hits:
        logger.info("No HMM hits to validate")
        return [], []

    if not marker_db.exists():
        logger.warning(f"Combined FAA database not found: {marker_db}")
        logger.warning("Skipping marker validation (returning all HMM hits)")
        return hits, []

    logger.info("=" * 60)
    logger.info("HMM-Gated Diamond Marker Validation")
    logger.info("=" * 60)
    logger.info(f"Validating {len(hits)} HMM hits")
    logger.info(f"Database: {marker_db}")

    # Step 1: Extract HMM-hit pORF sequences to temp FASTA
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        hmm_hit_fasta = tmpdir_path / "hmm_hit_sequences.faa"

        logger.info("Step 1: Extracting HMM-hit pORF sequences...")
        num_extracted = extract_hmm_hit_sequences(hits, proteome_fasta, hmm_hit_fasta)

        if num_extracted == 0:
            logger.warning("No sequences extracted for HMM hits")
            return [], []

        # Step 2: Run Diamond on HMM-hit pORFs only
        logger.info("Step 2: Running Diamond on HMM-hit pORFs (NOT all pORFs)...")
        diamond_output = tmpdir_path / "diamond_top10.tsv"

        try:
            run_diamond_on_hmm_hits(
                hmm_hit_fasta=hmm_hit_fasta,
                diamond_db=marker_db,
                output_tsv=diamond_output,
                threads=threads,
                evalue=evalue,
                max_seqs=top_k,
            )
        except Exception as e:
            logger.warning(f"Diamond search failed: {e}")
            logger.warning("Returning all HMM hits without validation")
            return hits, []

        # Step 3: Filter validated markers based on top-10 taxonomy
        logger.info("Step 3: Validating markers based on top-10 taxonomy...")
        validated_markers = filter_validated_markers(
            hmm_hits=hits,
            diamond_output_tsv=diamond_output,
            proteome_fasta=proteome_fasta,
            genome_fasta=genome_fasta,
            output_dir=output_dir,
        )

        # Copy Diamond results to output_dir if provided
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            diamond_out_path = output_dir / "diamond_top10.tsv"
            if diamond_output.exists():
                import shutil
                shutil.copy(diamond_output, diamond_out_path)
                logger.info(f"Copied Diamond top-10 results to: {diamond_out_path}")

    # Filter to only "validated" and "supported" markers
    validated_porf_ids = {
        vm.query_porf for vm in validated_markers
        if vm.validation_status in ("validated", "supported")
    }

    validated_hits = [h for h in hits if h.query_name in validated_porf_ids]

    logger.info(f"Marker validation complete:")
    logger.info(f"  Total HMM hits: {len(hits)}")
    logger.info(f"  Validated/Supported: {len(validated_hits)}")
    logger.info(f"  Reduction: {len(hits) - len(validated_hits)} hits filtered")
    logger.info("=" * 60)

    return validated_hits, validated_markers


# Broad packaging-ATPase marker(s) (A32 / FtsK-HerA family) that cross-hit ubiquitous
# cellular P-loop NTPases. A per-marker GA gathering floor is applied ONLY to these at
# HMM-search time; other markers keep E-value reporting (their GA cutoffs are not
# calibrated for gathering use, and enforcing them globally over-restricts seeding).
# The virophage VP_ATPase models carry higher GA cutoffs that would also reject divergent
# true-positive ATPases, so VP ATPase-only regions are handled at the Phase-3 acceptance
# gate (which requires a non-ATPase hallmark) rather than filtered at search time.
GA_FLOOR_MARKERS: frozenset = frozenset({"PLV_PC_054"})

# capscan MCP profiles (named *_caps_*) carry no GA cutoff and fire promiscuously at
# E<=1e-5 on giant-virus metagenomes (cross-hitting NCLDV/Mirus capsids). Apply capscan's
# own strong-hit bitscore threshold (>=75, capscan.sh) so they detect their groups without
# seeding false positives. Enforced via the same enforce_ga_cutoffs path.
CAPS_BITSCORE_FLOOR = 75.0


def run_hmmsearch(
    proteome_fasta: Path,
    hmm_profiles: list[HMM],
    evalue_cutoff: Optional[float] = None,
    threads: int = 4,
    chunk_size: Optional[int] = None,
    enforce_ga_cutoffs: bool = False,
) -> list[HMMHit]:
    """
    Run HMM search using pyhmmer.

    Args:
        proteome_fasta: Path to conceptual proteome FASTA
        hmm_profiles: List of pyhmmer HMM profile objects
        evalue_cutoff: Maximum E-value for reporting hits. None reports all hits.
        threads: Number of CPU threads
        chunk_size: Number of sequences per chunk (None to load all at once)
        enforce_ga_cutoffs: If True, apply per-marker bitscore floors (sequence and
            best-domain score) in addition to the E-value, for two marker classes:
            (1) the broad packaging-ATPase markers in GA_FLOOR_MARKERS, gated at their
            own GA gathering cutoff (suppresses A32 / FtsK-HerA cross-hits onto cellular
            P-loop NTPases); and (2) the integrated capscan MCP markers (name contains
            "_caps_"), which carry no GA cutoff, gated at CAPS_BITSCORE_FLOOR (=75, the
            capscan strong-hit threshold) to stop promiscuous cross-hits on giant-virus
            metagenomes. All other markers keep plain E-value reporting (their GA cutoffs
            are not calibrated for gathering use, and enforcing them globally collapses
            seeding).

    Returns:
        List of HMMHit objects for reported hits
    """
    hits = []

    alphabet = pyhmmer.easel.Alphabet.amino()
    search_options = (
        {"E": evalue_cutoff}
        if evalue_cutoff is not None
        else {
            "E": float("inf"),
            "domE": float("inf"),
            "incE": float("inf"),
            "incdomE": float("inf"),
        }
    )

    # Per-model gathering (GA) bitscore floors. When enabled, a hit must also meet
    # the model's own calibrated GA cutoff (sequence and best-domain score), not just
    # the E-value. This suppresses sub-GA cross-hits from broad markers (e.g. the
    # A32 / FtsK-HerA packaging ATPase) onto ubiquitous cellular P-loop NTPases.
    ga_map: dict[str, tuple[float, float]] = {}
    if enforce_ga_cutoffs:
        for _hmm in hmm_profiles:
            _name = _hmm.name.decode() if isinstance(_hmm.name, bytes) else _hmm.name
            if _name in GA_FLOOR_MARKERS and _hmm.cutoffs.gathering_available():
                ga_map[_name] = tuple(_hmm.cutoffs.gathering)
            elif "_caps_" in _name:
                # Use a calibrated GA cutoff if the capscan profile carries one;
                # otherwise fall back to the fixed bitscore floor.
                if _hmm.cutoffs.gathering_available():
                    ga_map[_name] = tuple(_hmm.cutoffs.gathering)
                else:
                    ga_map[_name] = (CAPS_BITSCORE_FLOOR, CAPS_BITSCORE_FLOOR)

    def search_sequences(sequence_batch: list[pyhmmer.easel.DigitalSequence]) -> None:
        for hmm, top_hits in zip(
            hmm_profiles,
            pyhmmer.hmmsearch(
                hmm_profiles,
                sequence_batch,
                cpus=threads,
                parallel="targets",
                **search_options,
            ),
        ):
            hmm_name = hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name

            for hit in top_hits:
                if evalue_cutoff is not None and hit.evalue > evalue_cutoff:
                    continue
                target_name = hit.name.decode() if isinstance(hit.name, bytes) else hit.name

                # Get best domain
                best_domain = None
                for domain in hit.domains:
                    if best_domain is None or domain.score > best_domain.score:
                        best_domain = domain

                if best_domain:
                    if enforce_ga_cutoffs and hmm_name in ga_map:
                        seq_ga, dom_ga = ga_map[hmm_name]
                        if hit.score < seq_ga or best_domain.score < dom_ga:
                            continue
                    hits.append(
                        HMMHit(
                            query_name=target_name,
                            target_name=hmm_name,
                            score=hit.score,
                            evalue=hit.evalue,
                            domain_score=best_domain.score,
                            query_start=best_domain.env_from,
                            query_end=best_domain.env_to,
                        )
                    )

    max_seq_length = 100000
    with SequenceFile(proteome_fasta, digital=True, alphabet=alphabet) as seq_handle:
        if chunk_size and chunk_size > 0:
            logger.info(
                "Searching pORFs in chunks of %d against %d HMM profiles",
                chunk_size,
                len(hmm_profiles),
            )
            batch: list[pyhmmer.easel.DigitalSequence] = []
            total = 0
            for seq in seq_handle:
                if len(seq) > max_seq_length:
                    logger.warning(
                        "Skipping long sequence %s (%d aa) - exceeds pyhmmer limit",
                        seq.name.decode() if isinstance(seq.name, bytes) else seq.name,
                        len(seq),
                    )
                    continue
                batch.append(seq)
                if len(batch) >= chunk_size:
                    total += len(batch)
                    logger.info(
                        "Running HMM chunk (%d-%d sequences)",
                        total - len(batch) + 1,
                        total,
                    )
                    search_sequences(batch)
                    batch = []
            if batch:
                total += len(batch)
                logger.info(
                    "Running HMM chunk (%d-%d sequences)",
                    total - len(batch) + 1,
                    total,
                )
                search_sequences(batch)
            logger.info("Searched %d pORFs across %d HMM profiles", total, len(hmm_profiles))
        else:
            sequences = []
            for seq in seq_handle:
                if len(seq) > max_seq_length:
                    logger.warning(
                        "Skipping long sequence %s (%d aa) - exceeds pyhmmer limit",
                        seq.name.decode() if isinstance(seq.name, bytes) else seq.name,
                        len(seq),
                    )
                    continue
                sequences.append(seq)
            logger.info("Searching %d pORFs against %d HMM profiles", len(sequences), len(hmm_profiles))
            search_sequences(sequences)

    if evalue_cutoff is None:
        logger.info("Found %d HMM hits with no HMM E-value reporting cutoff", len(hits))
    else:
        logger.info("Found %d HMM hits at E <= %s", len(hits), evalue_cutoff)
    return hits




def parse_porf_coordinates(
    porf_id: str,
    coord_lookup: Optional[dict[str, tuple[str, int, int, str]]] = None,
) -> Optional[tuple[str, int, int, str]]:
    """
    Parse pORF ID to extract genomic coordinates.

    Expected format: pORF_1|scaffold:start-end_strand:+_frame:1

    Returns:
        Tuple of (scaffold, start, end, strand) or None if parsing fails
    """
    if coord_lookup and porf_id in coord_lookup:
        return coord_lookup[porf_id]
    porf = PORF.parse_header(porf_id)
    if porf:
        return (porf.scaffold, porf.start, porf.end, porf.strand)
    return None


def identify_anchors(
    hits: list[HMMHit],
    min_score: float = 50.0,
    min_evalue: float = 1e-5,
    genome_lengths: Optional[dict[str, int]] = None,
    coord_lookup: Optional[dict[str, tuple[str, int, int, str]]] = None,
) -> list[Anchor]:
    """
    Identify anchor pORFs from HMM hits.

    Anchors are pORFs with strong hits to viral hallmark genes.
    If a pORF has multiple hits, the best scoring hit is used.

    Args:
        hits: List of HMM hits
        min_score: Minimum bit score for anchor designation
        min_evalue: Maximum E-value for anchor designation

    Returns:
        List of Anchor objects
    """
    # Group hits by pORF, keeping best hit per HMM
    porf_best_hits = defaultdict(dict)

    for hit in hits:
        porf_id = hit.query_name
        hmm_name = hit.target_name

        if hmm_name not in porf_best_hits[porf_id]:
            porf_best_hits[porf_id][hmm_name] = hit
        elif hit.score > porf_best_hits[porf_id][hmm_name].score:
            porf_best_hits[porf_id][hmm_name] = hit

    # Select anchors based on best hit
    anchors = []

    for porf_id, hmm_hits in porf_best_hits.items():
        # Get best overall hit for this pORF
        best_hit = max(hmm_hits.values(), key=lambda h: h.score)

        if best_hit.score >= min_score and best_hit.evalue <= min_evalue:
            coords = parse_porf_coordinates(porf_id, coord_lookup=coord_lookup)
            if not coords and genome_lengths:
                contig, frame, offset = parse_frame_id(porf_id)
                contig_len = genome_lengths.get(contig)
                if frame and contig_len:
                    start, end, strand = aa_to_nt_coords(
                        best_hit.query_start,
                        best_hit.query_end,
                        contig_len,
                        frame,
                        offset=offset,
                    )
                    coords = (contig, start, end, strand)
            if coords:
                scaffold, start, end, strand = coords
                anchors.append(
                    Anchor(
                        porf_id=porf_id,
                        scaffold=scaffold,
                        start=start,
                        end=end,
                        strand=strand,
                        hallmark_gene=best_hit.target_name,
                        score=best_hit.score,
                        evalue=best_hit.evalue,
                    )
                )

    logger.info(f"Identified {len(anchors)} anchor pORFs (score >= {min_score})")
    return anchors


def calculate_neighbor_scores(
    anchors: list[Anchor],
    all_hits: list[HMMHit],
    window_size: int = 50000,
) -> dict[str, float]:
    """
    Calculate neighbor density score for each anchor.

    The neighbor score measures the concentration of other viral-like pORFs
    in the genomic vicinity of an anchor.

    Args:
        anchors: List of anchor pORFs
        all_hits: All HMM hits (including weak ones)
        window_size: Size of window around anchor to scan (bp)

    Returns:
        Dictionary mapping pORF ID to neighbor score
    """
    # Build coordinate index for all hits
    hit_coords = {}
    for hit in all_hits:
        coords = parse_porf_coordinates(hit.query_name)
        if coords:
            scaffold, start, end, strand = coords
            if hit.query_name not in hit_coords:
                hit_coords[hit.query_name] = {
                    "scaffold": scaffold,
                    "start": start,
                    "end": end,
                    "score": hit.score,
                }

    neighbor_scores = {}

    for anchor in anchors:
        neighbor_count = 0
        neighbor_score_sum = 0.0

        for porf_id, coords in hit_coords.items():
            # Skip self
            if porf_id == anchor.porf_id:
                continue

            # Check same scaffold
            if coords["scaffold"] != anchor.scaffold:
                continue

            # Check within window
            distance = min(
                abs(coords["start"] - anchor.start),
                abs(coords["end"] - anchor.end),
            )

            if distance <= window_size:
                neighbor_count += 1
                neighbor_score_sum += coords["score"]

        # Normalize by window size (score per kb)
        normalized_score = neighbor_score_sum / (window_size / 1000)
        neighbor_scores[anchor.porf_id] = normalized_score

    return neighbor_scores


def form_seeds(
    anchors: list[Anchor],
    neighbor_scores: dict[str, float],
    min_neighbor_score: float = 10.0,
    flank_size: int = 50000,
    merge_distance: int = 10000,
    allow_isolated_anchors: bool = True,
    isolated_anchor_min_score: float = 100.0,
) -> list[HHGSeed]:
    """
    Form seed regions from qualified anchors.

    Seeds are formed by:
    1. Filtering anchors by neighbor score threshold
    2. Accepting isolated high-scoring anchors (for fragmented EVEs)
    3. Extending anchor regions by flank_size
    4. Merging overlapping/adjacent regions

    Args:
        anchors: List of anchor pORFs
        neighbor_scores: Dictionary of neighbor scores per anchor
        min_neighbor_score: Minimum neighbor score for seed formation
        flank_size: Distance to extend from anchor edges
        merge_distance: Maximum gap between anchors to merge into one seed
        allow_isolated_anchors: Accept high-scoring anchors even without
                               neighbor support (for fragmented EVEs)
        isolated_anchor_min_score: Minimum anchor HMM score to accept
                                   isolated anchors (default 100.0)

    Returns:
        List of HHGSeed objects
    """
    # Filter anchors by neighbor score OR high individual score
    qualified_anchors = []
    isolated_count = 0

    for a in anchors:
        neighbor_score = neighbor_scores.get(a.porf_id, 0)
        if neighbor_score >= min_neighbor_score:
            qualified_anchors.append(a)
        elif allow_isolated_anchors and a.score >= isolated_anchor_min_score:
            # Accept isolated high-scoring anchors for fragmented EVEs
            qualified_anchors.append(a)
            isolated_count += 1

    if isolated_count > 0:
        logger.info(f"  Accepted {isolated_count} isolated high-scoring anchors (score >= {isolated_anchor_min_score})")

    if not qualified_anchors:
        logger.info("No anchors met neighbor score threshold or isolated anchor criteria")
        return []

    # Group by scaffold
    scaffold_anchors = defaultdict(list)
    for anchor in qualified_anchors:
        scaffold_anchors[anchor.scaffold].append(anchor)

    seeds = []

    for scaffold, scaffold_anchor_list in scaffold_anchors.items():
        # Sort by position
        scaffold_anchor_list.sort(key=lambda a: a.start)

        # Merge overlapping regions
        current_anchors = [scaffold_anchor_list[0]]
        current_start = max(0, scaffold_anchor_list[0].start - flank_size)
        current_end = scaffold_anchor_list[0].end + flank_size

        for anchor in scaffold_anchor_list[1:]:
            anchor_start = max(0, anchor.start - flank_size)
            anchor_end = anchor.end + flank_size

            if anchor_start - current_end <= merge_distance:
                # Merge
                current_anchors.append(anchor)
                current_end = max(current_end, anchor_end)
            else:
                # Emit current seed
                max_neighbor_score = max(
                    neighbor_scores.get(a.porf_id, 0)
                    for a in current_anchors
                )
                seeds.append(
                    HHGSeed(
                        scaffold=scaffold,
                        start=current_start,
                        end=current_end,
                        anchors=current_anchors.copy(),
                        neighbor_score=max_neighbor_score,
                    )
                )
                # Start new seed
                current_anchors = [anchor]
                current_start = anchor_start
                current_end = anchor_end

        # Emit final seed
        max_neighbor_score = max(
            neighbor_scores.get(a.porf_id, 0)
            for a in current_anchors
        )
        seeds.append(
            HHGSeed(
                scaffold=scaffold,
                start=current_start,
                end=current_end,
                anchors=current_anchors,
                neighbor_score=max_neighbor_score,
            )
        )

    logger.info(f"Formed {len(seeds)} HHG seeds from {len(qualified_anchors)} qualified anchors")
    return seeds


def filter_seeds_by_diversity(
    seeds: list[HHGSeed],
    min_marker_types: int = 2,
    high_diversity_threshold: int = 3,
    high_neighbor_score_threshold: float = 5.0,
    assembly_mode: Optional[AssemblyMode] = None,
) -> list[HHGSeed]:
    """
    Filter seeds requiring marker diversity.

    Real NCLDV EVEs contain multiple different hallmark genes. Seeds with only
    one marker type (especially universal markers like PolB, RNAPL) are likely
    false positives from host genes.

    Assembly modes modify filtering behavior:
    - default: Requires marker diversity (min 2 types)
    - fragmented: Accepts single high-scoring diagnostic markers
    - relaxed: Accepts any marker with sufficient score
    - strict: Requires high diversity (min 3 types)

    Filtering logic (default mode):
    - Seeds with ≥high_diversity_threshold marker types: accept (high diversity)
    - Seeds with min_marker_types and virus-specific marker: accept
    - Seeds with min_marker_types and high neighbor score: accept (degraded EVEs)
    - Seeds with <min_marker_types: reject (too ambiguous)

    Args:
        seeds: List of HHGSeed objects
        min_marker_types: Minimum number of different marker types required
        high_diversity_threshold: Number of marker types that bypasses virus-specific requirement
        high_neighbor_score_threshold: Neighbor score that allows universal-only markers
        assembly_mode: AssemblyMode configuration (overrides other params if set)

    Returns:
        Filtered list of HHGSeed objects
    """
    if not seeds:
        return []

    # Apply assembly mode settings if provided
    if assembly_mode is not None:
        min_marker_types = assembly_mode.min_marker_types
        high_neighbor_score_threshold = assembly_mode.neighbor_score_threshold

        # In fragmented/relaxed mode, accept single diagnostic markers
        if assembly_mode.accept_single_diagnostic:
            logger.info(f"Assembly mode '{assembly_mode.mode}': accepting single diagnostic markers")

    filtered = []
    rejected_diversity = 0
    rejected_no_specific = 0
    accepted_high_diversity = 0
    accepted_virus_specific = 0
    accepted_high_neighbor = 0
    accepted_single_diagnostic = 0
    accepted_high_score_single = 0
    accepted_single_gvogm = 0

    for seed in seeds:
        marker_types = set(a.hallmark_gene.lower() for a in seed.anchors)
        num_markers = len(marker_types)
        max_anchor_score = max(a.score for a in seed.anchors) if seed.anchors else 0

        # Check for diagnostic markers using new viral_markers module
        diagnostic_markers = marker_types & ALL_DIAGNOSTIC_MARKERS
        has_diagnostic = len(diagnostic_markers) > 0

        # Predict viral family based on markers
        predicted_family = get_family_for_markers(marker_types)

        # === FRAGMENTED/RELAXED MODE: Accept single high-scoring markers ===
        if assembly_mode is not None and assembly_mode.accept_single_diagnostic:
            # Family-specific overrides: a predicted family may declare a
            # stricter single-marker bit-score floor or a higher minimum-
            # marker requirement than the generic assembly-mode threshold.
            # The Mriyavirus profile for instance is tightened to reject
            # the permissive 60-bit / 1-marker calls that matched TE-borne
            # HUH endonuclease / SF3 helicase homologs.
            family_profile = VIRAL_FAMILIES.get(predicted_family)
            family_single_min_score = assembly_mode.single_marker_min_score
            family_min_markers = 1
            if family_profile is not None:
                family_single_min_score = max(
                    family_single_min_score,
                    float(family_profile.single_marker_min_score),
                )
                family_min_markers = max(
                    family_min_markers,
                    int(family_profile.min_markers_fragmented),
                )

            family_marker_gate_ok = num_markers >= family_min_markers

            # Accept ANY marker with high enough score in fragmented mode,
            # subject to the family-specific floor and marker-count gate.
            if (
                family_marker_gate_ok
                and max_anchor_score >= family_single_min_score
            ):
                logger.debug(
                    f"Accepted seed {seed.scaffold}:{seed.start}-{seed.end}: "
                    f"high-scoring marker (score={max_anchor_score:.1f}, markers={marker_types}, "
                    f"predicted_family={predicted_family}, "
                    f"family_single_min_score={family_single_min_score:.1f}, "
                    f"family_min_markers={family_min_markers})"
                )
                # Annotate seed with predicted family
                seed.predicted_family = predicted_family
                filtered.append(seed)
                accepted_high_score_single += 1
                continue

            # Accept diagnostic markers even with lower scores
            if has_diagnostic and family_marker_gate_ok:
                logger.debug(
                    f"Accepted seed {seed.scaffold}:{seed.start}-{seed.end}: "
                    f"diagnostic marker(s) {diagnostic_markers} in fragmented mode "
                    f"(predicted_family={predicted_family}, "
                    f"family_min_markers={family_min_markers})"
                )
                seed.predicted_family = predicted_family
                filtered.append(seed)
                accepted_single_diagnostic += 1
                continue

        # === STANDARD MODE: Require diversity ===

        if num_markers == 1:
            marker_name = next(iter(marker_types))
            if marker_name.startswith("gvogm"):
                logger.debug(
                    f"Accepted seed {seed.scaffold}:{seed.start}-{seed.end}: "
                    f"single GVOGm marker ({marker_name})"
                )
                seed.predicted_family = predicted_family
                filtered.append(seed)
                accepted_single_gvogm += 1
                continue

        # Require minimum diversity
        if num_markers < min_marker_types:
            logger.debug(
                f"Rejected seed {seed.scaffold}:{seed.start}-{seed.end}: "
                f"only {num_markers} marker type(s): {marker_types}"
            )
            rejected_diversity += 1
            continue

        # High diversity: accept regardless of virus-specific markers
        if num_markers >= high_diversity_threshold:
            logger.debug(
                f"Accepted seed {seed.scaffold}:{seed.start}-{seed.end}: "
                f"high diversity ({num_markers} markers: {marker_types})"
            )
            seed.predicted_family = predicted_family
            filtered.append(seed)
            accepted_high_diversity += 1
            continue

        # Medium diversity: check for virus-specific marker (using both old and new definitions)
        virus_specific = marker_types & (VIRUS_SPECIFIC_MARKERS | ALL_DIAGNOSTIC_MARKERS)
        if virus_specific:
            logger.debug(
                f"Accepted seed {seed.scaffold}:{seed.start}-{seed.end}: "
                f"has virus-specific marker(s): {virus_specific}"
            )
            seed.predicted_family = predicted_family
            filtered.append(seed)
            accepted_virus_specific += 1
            continue

        # Medium diversity without virus-specific: accept if high neighbor score
        # (indicates true viral clustering, likely a degraded EVE)
        if seed.neighbor_score >= high_neighbor_score_threshold:
            logger.debug(
                f"Accepted seed {seed.scaffold}:{seed.start}-{seed.end}: "
                f"high neighbor score ({seed.neighbor_score:.1f}) with {num_markers} markers"
            )
            seed.predicted_family = predicted_family
            filtered.append(seed)
            accepted_high_neighbor += 1
            continue

        # Medium diversity, low neighbor score, no virus-specific: reject
        logger.debug(
            f"Rejected seed {seed.scaffold}:{seed.start}-{seed.end}: "
            f"only {num_markers} markers without virus-specific (has: {marker_types}), "
            f"neighbor_score={seed.neighbor_score:.1f}"
        )
        rejected_no_specific += 1

    # Build summary message
    accepted_parts = []
    if accepted_high_diversity > 0:
        accepted_parts.append(f"{accepted_high_diversity} high-diversity")
    if accepted_virus_specific > 0:
        accepted_parts.append(f"{accepted_virus_specific} virus-specific")
    if accepted_high_neighbor > 0:
        accepted_parts.append(f"{accepted_high_neighbor} high-neighbor")
    if accepted_single_diagnostic > 0:
        accepted_parts.append(f"{accepted_single_diagnostic} single-diagnostic")
    if accepted_high_score_single > 0:
        accepted_parts.append(f"{accepted_high_score_single} high-score-single")
    if accepted_single_gvogm > 0:
        accepted_parts.append(f"{accepted_single_gvogm} single-gvogm")

    rejected_parts = []
    if rejected_diversity > 0:
        rejected_parts.append(f"{rejected_diversity} low-diversity")
    if rejected_no_specific > 0:
        rejected_parts.append(f"{rejected_no_specific} medium-no-specific")

    logger.info(
        f"Marker diversity filter: {len(seeds)} → {len(filtered)} seeds "
        f"(accepted: {', '.join(accepted_parts) or 'none'}; "
        f"rejected: {', '.join(rejected_parts) or 'none'})"
    )

    return filtered


def hhg_seeding_pipeline(
    proteome_fasta: Path,
    hmm_file: Path,
    genome_fasta: Optional[Path] = None,
    hmm_allowlist: Optional[Path] = None,
    marker_faa_dir: Optional[Path] = None,
    marker_db: Optional[Path] = None,
    marker_top_k: int = 10,
    marker_neighbor_genes: int = 5,
    evalue_cutoff: Optional[float] = None,
    min_anchor_score: float = 50.0,
    min_neighbor_score: float = 5.0,
    window_size: int = 50000,
    threads: int = 4,
    hmm_chunk_size: Optional[int] = None,
    enforce_ga_cutoffs: bool = True,
    min_marker_types: int = 2,
    high_diversity_threshold: int = 3,
    allow_isolated_anchors: bool = True,
    isolated_anchor_min_score: float = 100.0,
    assembly_mode: Optional[str] = None,
    return_hits: bool = False,
    output_dir: Optional[Path] = None,
) -> list[HHGSeed] | tuple[list[HHGSeed], list[HMMHit]]:
    """
    Full HHG seeding pipeline.

    This is the main entry point for Path A of Phase 1.

    Marker Validation Options (choose ONE):
    1. marker_db: HMM-gated Diamond workflow (CURRENT, RECOMMENDED)
       - Points to marker.dmnd (SMALL, ~1,500 curated marker sequences)
       - Only runs Diamond on HMM-hit pORFs (1-10% of total proteome)
       - Validates markers with top-10 taxonomy analysis
       - Implements PIPELINE_HMM_GATED_PLAN.md Steps 2-3
    2. marker_faa_dir: Legacy per-marker FAA validation (DEPRECATED, slow)
       - Runs separate Diamond searches for each marker type (MCP, PolB, etc.)
       - Requires neighbor support within genes
       - No longer used (set to null in config)

    Args:
        proteome_fasta: Path to conceptual proteome from Phase 0
        hmm_file: Path to viral hallmark HMM database
        hmm_allowlist: Optional path to file with allowed HMM names (one per line)
        marker_faa_dir: Optional directory containing marker-specific FAA files (LEGACY, not used)
        marker_db: Path to marker.dmnd for HMM-gated validation (TIER 1, ~1,500 curated marker sequences)
        marker_top_k: Top-k hits to check for NCLDV/MIRUS/VP/PLV support (default 10)
        marker_neighbor_genes: Neighbor window (genes) for supporting markers (legacy mode only)
        evalue_cutoff: E-value threshold for HMM search. None reports all hits.
        min_anchor_score: Minimum bit score for anchor designation
        min_neighbor_score: Minimum neighbor score for seed formation
        window_size: Window size for neighbor scoring (bp)
        threads: Number of CPU threads
        hmm_chunk_size: Number of sequences per hmmsearch chunk (None to disable)
        min_marker_types: Minimum number of different marker types required per seed
        high_diversity_threshold: Number of marker types that bypasses virus-specific requirement
        allow_isolated_anchors: Accept high-scoring anchors without neighbor support
                               (for fragmented EVEs). Default True.
        isolated_anchor_min_score: Minimum HMM score to accept isolated anchors.
                                   Default 100.0 (strong hit only).
        assembly_mode: Assembly mode name ("default", "fragmented", "relaxed", "strict").
                      If provided, overrides min_marker_types and related settings.
                      Use "fragmented" for MAGs and highly fragmented assemblies.

    Returns:
        If return_hits=False: List of HHGSeed objects
        If return_hits=True: Tuple of (seeds, HMM hits)

    Example:
        >>> # Standard mode for assembled genomes
        >>> seeds = hhg_seeding_pipeline(
        ...     Path("conceptual_proteome.faa"),
        ...     Path("data/databases/viralrecall_v3/hmm/gvog.hmm"),
        ...     threads=16
        ... )
        >>> print(f"Found {len(seeds)} candidate EVE regions")

        >>> # Fragmented mode for MAGs
        >>> seeds = hhg_seeding_pipeline(
        ...     Path("mag_proteome.faa"),
        ...     Path("data/databases/viralrecall_v3/hmm/gvog.hmm"),
        ...     threads=16,
        ...     assembly_mode="fragmented"
        ... )
        >>> print(f"Found {len(seeds)} candidate regions in MAG")
    """
    # Get assembly mode configuration
    mode_config = None
    if assembly_mode is not None:
        mode_config = get_assembly_mode(assembly_mode)
        logger.info(f"Using assembly mode: {assembly_mode}")
        # Override settings from mode config
        min_anchor_score = mode_config.single_marker_min_score
        min_neighbor_score = mode_config.neighbor_score_threshold
        min_marker_types = mode_config.min_marker_types
        allow_isolated_anchors = mode_config.accept_isolated_anchors

    logger.info("=" * 60)
    logger.info("HHG Seeding Pipeline")
    logger.info("=" * 60)

    # Step 1: Load HMM profiles
    logger.info("Step 1: Loading HMM profiles...")
    allowlist = load_hmm_allowlist(hmm_allowlist)
    faa_markers = load_faa_markers(marker_faa_dir) if marker_faa_dir else None
    if faa_markers is not None:
        allowlist = allowlist & faa_markers if allowlist else faa_markers
        logger.info(
            "Filtered HMMs to %d markers with FAA files", len(allowlist)
        )
    hmm_profiles = load_hmm_profiles(hmm_file, allowlist=allowlist)

    # Step 2: Run HMM search
    logger.info("Step 2: Running HMM search...")
    hits = run_hmmsearch(
        proteome_fasta,
        hmm_profiles,
        evalue_cutoff,
        threads,
        chunk_size=hmm_chunk_size,
        enforce_ga_cutoffs=enforce_ga_cutoffs,
    )

    if not hits:
        logger.warning("No HMM hits found")
        return ([], []) if return_hits else []

    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        hits_path = output_dir / "hmm_hits.tsv"
        with hits_path.open("w") as handle:
            handle.write("query\ttarget\tscore\tevalue\tdomain_score\tquery_start\tquery_end\n")
            for hit in hits:
                handle.write(
                    f"{hit.query_name}\t{hit.target_name}\t{hit.score:.3f}\t{hit.evalue:.3e}\t"
                    f"{hit.domain_score:.3f}\t{hit.query_start}\t{hit.query_end}\n"
                )
        logger.info("Wrote HMM hits: %s", hits_path)

    # Step 2b: Validate markers (choose validation method)
    if marker_db:
        # NEW: HMM-gated Diamond validation (FAST, recommended)
        logger.info("Step 2b: HMM-gated marker validation with combined database...")
        hits, validated_markers = validate_hmm_hits_with_combined_db(
            hits=hits,
            proteome_fasta=proteome_fasta,
            marker_db=marker_db,
            threads=threads,
            top_k=marker_top_k,
            output_dir=output_dir,
            genome_fasta=genome_fasta,
        )
        if not hits:
            logger.warning("No HMM hits after HMM-gated validation")
            return ([], []) if return_hits else []
    elif marker_faa_dir:
        # LEGACY: Per-marker FAA validation (SLOW, deprecated)
        logger.warning("Using legacy per-marker validation (DEPRECATED)")
        logger.warning("Consider using marker_db parameter for 10-100x speedup")
        hits = validate_hmm_hits_with_markers(
            hits=hits,
            proteome_fasta=proteome_fasta,
            marker_faa_dir=marker_faa_dir,
            threads=threads,
            top_k=marker_top_k,
            neighbor_genes=marker_neighbor_genes,
            output_dir=output_dir,
        )
        if not hits:
            logger.warning("No HMM hits after marker validation")
            return ([], []) if return_hits else []

    if output_dir:
        output_dir = Path(output_dir)
        validated_path = output_dir / "hmm_hits_validated.tsv"
        with validated_path.open("w") as handle:
            handle.write("query\ttarget\tscore\tevalue\tdomain_score\tquery_start\tquery_end\n")
            for hit in hits:
                handle.write(
                    f"{hit.query_name}\t{hit.target_name}\t{hit.score:.3f}\t{hit.evalue:.3e}\t"
                    f"{hit.domain_score:.3f}\t{hit.query_start}\t{hit.query_end}\n"
                )
        logger.info("Wrote validated HMM hits: %s", validated_path)

    # Step 3: Identify anchors
    logger.info("Step 3: Identifying anchor sequences from HMM hits...")
    genome_lengths = None
    if genome_fasta:
        genome_lengths = {
            rec.id: len(rec.seq)
            for rec in SeqIO.parse(genome_fasta, "fasta")
        }
    coord_lookup = {}
    try:
        for record in SeqIO.parse(proteome_fasta, "fasta"):
            parsed = parse_prodigal_header(record.description, record.id)
            if parsed:
                scaffold, start, end, strand = parsed
                coord_lookup[record.id] = (scaffold, start, end, strand)
    except FileNotFoundError:
        coord_lookup = {}
    anchors = identify_anchors(
        hits,
        min_anchor_score,
        genome_lengths=genome_lengths,
        coord_lookup=coord_lookup if coord_lookup else None,
    )

    if not anchors:
        logger.warning("No anchors identified from HMM hits")
        return ([], hits) if return_hits else []

    # Step 4: Calculate neighbor scores
    logger.info("Step 4: Calculating neighbor scores...")
    neighbor_scores = calculate_neighbor_scores(anchors, hits, window_size)

    # Step 5: Form seeds
    logger.info("Step 5: Forming seeds...")
    logger.info(f"  Neighbor score threshold: {min_neighbor_score}")
    logger.info(f"  Allow isolated anchors: {allow_isolated_anchors} (min score: {isolated_anchor_min_score})")
    seeds = form_seeds(
        anchors,
        neighbor_scores,
        min_neighbor_score=min_neighbor_score,
        flank_size=window_size,
        allow_isolated_anchors=allow_isolated_anchors,
        isolated_anchor_min_score=isolated_anchor_min_score,
    )

    # Step 6: Filter by marker diversity
    logger.info("Step 6: Filtering by marker diversity...")
    seeds = filter_seeds_by_diversity(
        seeds,
        min_marker_types=min_marker_types,
        high_diversity_threshold=high_diversity_threshold,
        assembly_mode=mode_config,
    )

    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        bed_path = output_dir / "hhg_seeds.bed"
        with bed_path.open("w") as handle:
            for idx, seed in enumerate(seeds, start=1):
                name = f"HHG_{idx}_{seed.scaffold}_{seed.start}_{seed.end}"
                score = min(int(seed.neighbor_score * 100), 1000)
                handle.write(
                    f"{seed.scaffold}\t{seed.start}\t{seed.end}\t{name}\t{score}\t.\n"
                )
        logger.info("Wrote HHG seeds BED: %s", bed_path)

    # Log summary
    if seeds:
        total_length = sum(s.length for s in seeds)
        total_anchors = sum(s.num_anchors for s in seeds)
        unique_hallmarks = set()
        for s in seeds:
            unique_hallmarks.update(s.hallmark_genes)

        logger.info(f"HHG seeding complete:")
        logger.info(f"  Seeds: {len(seeds)}")
        logger.info(f"  Total coverage: {total_length:,} bp")
        logger.info(f"  Total anchors: {total_anchors}")
        logger.info(f"  Unique hallmark genes: {len(unique_hallmarks)}")

    return (seeds, hits) if return_hits else seeds
