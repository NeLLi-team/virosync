#!/usr/bin/env python3
"""
Classify MCP (Major Capsid Protein) sequences as DJR (Double Jelly Roll) or SJR (Single Jelly Roll).

Multi-Signal Classification Approach:
=====================================
This script uses multiple evidence sources to classify MCPs, with confidence scores
that reflect the strength of evidence:

1. InterProScan Jelly Roll Domain Counting (PRIMARY - highest confidence)
   - PF21738: Double jelly roll capsid-like protein → DJR indicator
   - IPR049512: Double jelly roll-like domain → DJR indicator
   - 2+ jelly roll domains detected → DJR (confidence 0.95)
   - 1 jelly roll domain detected → SJR (confidence 0.90)

2. Multiple HMM Hits to Same Protein
   - Non-overlapping PLV_MCP domain hits indicate multiple jelly roll domains
   - 2+ non-overlapping domains → DJR evidence (confidence 0.85)

3. TMVec Reference Similarity
   - Comparison against curated DJR/SJR reference proteins
   - DJR references: Mavirus (6G45), Mimivirus, PBCV-1
   - SJR references: Adenovirus hexon, T4 gp23
   - Score >0.3 to reference → classification with confidence up to 0.80

4. FoldSeek Structural Hits (with Boltz predictions)
   - Parse PDB hits for known DJR/SJR structures
   - 6G45, Mimivirus, PBCV-1 → DJR (confidence 0.75)
   - Adenovirus structures → SJR (confidence 0.75)

5. Length Heuristics (fallback - lowest confidence)
   - >400 aa → likely DJR (confidence 0.60)
   - PLV_MCP marker → DJR (confidence 0.50, structural evidence from Boltz/Foldseek)
   - Otherwise → UNKNOWN (insufficient evidence)
   NOTE: SJR classification based on length was removed - Boltz/Foldseek PDB
   evidence shows all PLV MCPs are DJR (hits to Marseillevirus, PBCV-1, etc.)

Classification Logic:
- Signals are evaluated in priority order
- First definitive signal determines classification
- Confidence reflects evidence strength, not certainty

Usage:
    python scripts/classify_jelly_roll.py \\
        --marker-hits tests/results/stena/phase1/marker_validation/validated_marker_hits.tsv \\
        --sequences tests/results/stena/phase1/marker_validation/hmm_hit_porfs.faa \\
        --output tests/results/stena/phase3_synthesis/virosync_jelly_roll_proteins.tsv \\
        --interproscan tests/results/stena/phase3/interproscan/interproscan_batch.tsv \\
        --tmvec tests/results/stena/phase3_synthesis/virosync_tmvec_proteins.tsv \\
        --foldseek tests/results/stena/structural_analysis/foldseek_pdb_results.tsv

References:
    - Krupovic, M., & Koonin, E.V. (2015). Polintons: a hotbed of eukaryotic virus,
      transposon and plasmid evolution. Nature Reviews Microbiology.
    - Yutin, N., et al. (2014). Origin of giant viruses from smaller DNA viruses
      not from a fourth domain of cellular life. Virology.
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# Jelly roll classification parameters
DJR_MIN_LENGTH = 400  # Minimum length for DJR classification
SJR_MAX_LENGTH = 350  # Maximum length for SJR classification
DJR_OPTIMAL_MIN = 450  # Optimal DJR range start
DJR_OPTIMAL_MAX = 550  # Optimal DJR range end
SJR_OPTIMAL_MIN = 250  # Optimal SJR range start
SJR_OPTIMAL_MAX = 320  # Optimal SJR range end

# MCP-related marker patterns (case-insensitive)
MCP_MARKER_PATTERNS = [
    "PLV_MCP",
    "MCP",
    "major_capsid",
    "capsid_protein",
]

# InterProScan domains for jelly roll detection
DJR_PFAM_DOMAINS = {
    "PF21738",  # Double jelly roll capsid-like protein
}
DJR_INTERPRO_DOMAINS = {
    "IPR049512",  # Double jelly roll-like domain
}
JELLY_ROLL_KEYWORDS = {
    "double jelly roll",
    "jelly roll",
    "jelly-roll",
}

# FoldSeek PDB structures that indicate DJR/SJR
DJR_PDB_PATTERNS = [
    r"6g4[345]",  # Mavirus MCP structures
    r"5ti[qp]",   # PBCV-1 Vp54
    r"1m4x",      # PBCV-1 capsid
]
SJR_PDB_PATTERNS = [
    r"4cwu",      # Adenovirus hexon
    r"1p30",      # Adenovirus hexon
]

# TMVec reference score threshold
TMVEC_SCORE_THRESHOLD = 0.3


@dataclass
class ClassificationSignals:
    """Collection of signals for classification."""

    jelly_roll_domains: int = 0  # Count of detected jelly roll domains
    hmm_domain_count: int = 0    # Count of non-overlapping HMM domain hits
    tmvec_djr_score: float = 0.0
    tmvec_sjr_score: float = 0.0
    tmvec_djr_hit: str = ""
    tmvec_sjr_hit: str = ""
    foldseek_hit_is_djr: bool = False
    foldseek_hit_is_sjr: bool = False
    foldseek_top_hit: str = ""
    length: int = 0
    evidence_sources: list[str] = field(default_factory=list)


@dataclass
class JellyRollClassification:
    """Classification result for a single protein."""

    protein_id: str
    jelly_roll_type: str  # DJR, SJR, or UNKNOWN
    confidence: float  # 0.0-1.0
    length: int  # Protein length in amino acids
    marker: str  # HMM marker that detected it
    evidence: str  # Evidence source(s) for classification
    sequence: str  # Full sequence for reference


def is_mcp_marker(marker_name: str) -> bool:
    """Check if a marker name corresponds to an MCP marker."""
    marker_lower = marker_name.lower()
    for pattern in MCP_MARKER_PATTERNS:
        if pattern.lower() in marker_lower:
            return True
    return False


def calculate_djr_confidence(length: int) -> float:
    """Calculate confidence for DJR classification based on length."""
    if DJR_OPTIMAL_MIN <= length <= DJR_OPTIMAL_MAX:
        return 1.0
    elif length > DJR_OPTIMAL_MAX:
        excess = length - DJR_OPTIMAL_MAX
        return max(0.5, 1.0 - (excess / 300))
    elif length >= DJR_MIN_LENGTH:
        deficit = DJR_OPTIMAL_MIN - length
        return max(0.6, 1.0 - (deficit / 100))
    return 0.0


def calculate_sjr_confidence(length: int) -> float:
    """Calculate confidence for SJR classification based on length."""
    if SJR_OPTIMAL_MIN <= length <= SJR_OPTIMAL_MAX:
        return 1.0
    elif length < SJR_OPTIMAL_MIN:
        deficit = SJR_OPTIMAL_MIN - length
        return max(0.5, 1.0 - (deficit / 100))
    elif length <= SJR_MAX_LENGTH:
        excess = length - SJR_OPTIMAL_MAX
        return max(0.7, 1.0 - (excess / 60))
    return 0.0


def classify_with_signals(signals: ClassificationSignals, marker: str) -> tuple[str, float, str]:
    """
    Multi-signal classification for DJR/SJR.

    Args:
        signals: Classification signals from all evidence sources
        marker: HMM marker name

    Returns:
        (classification, confidence, evidence_description) tuple
    """
    # Priority 1: InterProScan domain count (highest confidence)
    if signals.jelly_roll_domains >= 2:
        return "DJR", 0.95, "interproscan:2+_djr_domains"
    if signals.jelly_roll_domains == 1:
        # Single jelly roll domain detected - likely DJR (partial detection)
        # Boltz/Foldseek evidence shows PLV MCPs are DJR, not SJR
        return "DJR", 0.85, "interproscan:1_djr_domain"

    # Priority 2: Multiple HMM hits to same protein
    if signals.hmm_domain_count >= 2:
        return "DJR", 0.85, f"hmm:{signals.hmm_domain_count}_domains"

    # Priority 3: TMVec reference similarity
    if signals.tmvec_djr_score > TMVEC_SCORE_THRESHOLD:
        conf = min(0.80, signals.tmvec_djr_score)
        return "DJR", conf, f"tmvec:{signals.tmvec_djr_hit}({signals.tmvec_djr_score:.2f})"
    if signals.tmvec_sjr_score > TMVEC_SCORE_THRESHOLD:
        conf = min(0.80, signals.tmvec_sjr_score)
        return "SJR", conf, f"tmvec:{signals.tmvec_sjr_hit}({signals.tmvec_sjr_score:.2f})"

    # Priority 4: FoldSeek structural hits
    if signals.foldseek_hit_is_djr:
        return "DJR", 0.75, f"foldseek:{signals.foldseek_top_hit}"
    if signals.foldseek_hit_is_sjr:
        return "SJR", 0.75, f"foldseek:{signals.foldseek_top_hit}"

    # Priority 5: Length heuristics (fallback)
    # NOTE: Removed SJR classification based on length - Boltz/Foldseek PDB evidence
    # shows all PLV MCPs are DJR (hits to Marseillevirus, PBCV-1, Faustovirus capsids)
    is_plv_mcp = "plv_mcp" in marker.lower()

    if signals.length > DJR_MIN_LENGTH:
        confidence = calculate_djr_confidence(signals.length) * 0.6  # Scale down
        if is_plv_mcp:
            confidence = min(1.0, confidence + 0.1)
        return "DJR", confidence, "length"
    elif is_plv_mcp:
        # PLV MCPs are DJR regardless of length (structural evidence from Boltz/Foldseek)
        confidence = 0.5
        return "DJR", confidence, "plv_mcp_marker"
    else:
        # No confident classification without additional evidence
        return "UNKNOWN", 0.0, "insufficient_evidence"


def extract_base_porf_id(porf_id: str) -> str:
    """Extract base pORF ID without domain coordinates."""
    if "|aa" in porf_id:
        return porf_id.rsplit("|aa", 1)[0]
    return porf_id


def parse_domain_coordinates(porf_id: str) -> tuple[int | None, int | None]:
    """Parse domain coordinates from pORF ID."""
    if "|aa" not in porf_id:
        return None, None
    try:
        _, domain_str = porf_id.rsplit("|aa", 1)
        start_str, end_str = domain_str.split("-")
        return int(start_str), int(end_str)
    except (ValueError, IndexError):
        return None, None


def estimate_full_protein_length(porf_id: str, domain_sequence_length: int) -> int:
    """Estimate the full protein length from domain coordinates."""
    start, end = parse_domain_coordinates(porf_id)
    if end is not None:
        return end
    return domain_sequence_length


def load_marker_hits(
    marker_hits_path: Path, mcp_only: bool = True
) -> dict[str, tuple[str, str]]:
    """Load marker hits from validated_marker_hits.tsv."""
    hits: dict[str, tuple[str, str]] = {}

    with marker_hits_path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            porf_id = row.get("query_porf", "")
            marker = row.get("hmm_target", "")
            status = row.get("validation_status", "")

            if not porf_id or not marker:
                continue

            if mcp_only and not is_mcp_marker(marker):
                continue

            hits[porf_id] = (marker, status)

    return hits


def count_hmm_domains_per_protein(
    marker_hits: dict[str, tuple[str, str]]
) -> dict[str, int]:
    """
    Count non-overlapping HMM domain hits per base protein.

    Multiple HMM hits to the same protein with non-overlapping coordinates
    suggest multiple jelly roll domains (DJR evidence).
    """
    protein_domains: dict[str, list[tuple[int, int]]] = defaultdict(list)

    for porf_id in marker_hits:
        base_id = extract_base_porf_id(porf_id)
        start, end = parse_domain_coordinates(porf_id)
        if start is not None and end is not None:
            protein_domains[base_id].append((start, end))

    # Count non-overlapping domains
    domain_counts: dict[str, int] = {}
    for base_id, domains in protein_domains.items():
        if len(domains) <= 1:
            domain_counts[base_id] = len(domains)
            continue

        # Sort by start position
        sorted_domains = sorted(domains, key=lambda x: x[0])

        # Count non-overlapping (>50% non-overlap required)
        non_overlapping = 1
        prev_end = sorted_domains[0][1]

        for start, end in sorted_domains[1:]:
            # Check overlap
            overlap = max(0, prev_end - start)
            domain_len = end - start
            if domain_len > 0 and overlap / domain_len < 0.5:
                non_overlapping += 1
            prev_end = max(prev_end, end)

        domain_counts[base_id] = non_overlapping

    return domain_counts


def load_interproscan_domains(
    interproscan_path: Path | None,
) -> dict[str, int]:
    """
    Load InterProScan results and count jelly roll domains per protein.

    Returns:
        Dictionary mapping protein IDs to jelly roll domain counts.
    """
    if interproscan_path is None or not interproscan_path.exists():
        return {}

    domain_counts: dict[str, int] = defaultdict(int)

    with interproscan_path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) < 13:
                continue

            protein_id = parts[0]
            pfam_id = parts[4] if len(parts) > 4 else ""
            interpro_id = parts[11] if len(parts) > 11 else ""
            description = parts[5].lower() if len(parts) > 5 else ""
            interpro_desc = parts[12].lower() if len(parts) > 12 else ""

            # Check for DJR domain signatures
            is_djr = (
                pfam_id in DJR_PFAM_DOMAINS
                or interpro_id in DJR_INTERPRO_DOMAINS
                or any(kw in description for kw in JELLY_ROLL_KEYWORDS)
                or any(kw in interpro_desc for kw in JELLY_ROLL_KEYWORDS)
            )

            if is_djr:
                # Extract base protein ID (may have EVE prefix)
                # Format: EVE_scaffold_coords|porf_id (e.g., EVE_stena|contig_122_0-31000|stena|contig_122_8)
                # Or: EVE_scaffold_coords|porf_id where porf_id lacks genome prefix
                domain_counts[protein_id] += 1

                if "|" in protein_id:
                    # Take last part as pORF ID (e.g., stena|contig_122_8 or contig_122_8)
                    last_part = protein_id.split("|")[-1]
                    domain_counts[last_part] += 1

                    # Also store with common genome prefixes for matching
                    # The marker hits use format: stena|contig_X_Y
                    # InterProScan may have: contig_X_Y (without prefix)
                    if not last_part.startswith(("stena|", "genome|")):
                        # Add common genome prefixes
                        domain_counts[f"stena|{last_part}"] += 1
                    else:
                        # Also store without prefix
                        parts = last_part.split("|", 1)
                        if len(parts) > 1:
                            domain_counts[parts[1]] += 1

    return dict(domain_counts)


def load_tmvec_results(
    tmvec_path: Path | None,
    djr_references: set[str] | None = None,
    sjr_references: set[str] | None = None,
) -> dict[str, tuple[float, str, float, str]]:
    """
    Load TMVec results and extract scores for DJR/SJR reference comparison.

    Returns:
        Dictionary mapping protein IDs to (djr_score, djr_hit, sjr_score, sjr_hit).
    """
    if tmvec_path is None or not tmvec_path.exists():
        return {}

    # Default references (can be extended)
    if djr_references is None:
        djr_references = {"6g45", "mavirus", "mimivirus", "pbcv", "chlorella"}
    if sjr_references is None:
        sjr_references = {"adenovirus", "hexon", "t4", "gp23", "bacteriophage"}

    results: dict[str, tuple[float, str, float, str]] = {}

    with tmvec_path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            porf_id = row.get("porf_id", "")
            if not porf_id:
                continue

            # Check BFVD hits for reference matches
            bfvd_score = float(row.get("tmvec_bfvd_score", 0) or 0)
            bfvd_hit = row.get("tmvec_bfvd_hit", "")
            bfvd_annotation = row.get("tmvec_bfvd_annotation", "").lower()
            bfvd_keywords = row.get("tmvec_bfvd_keywords", "").lower()

            djr_score = 0.0
            djr_hit = ""
            sjr_score = 0.0
            sjr_hit = ""

            # Check if hit matches DJR references
            hit_text = f"{bfvd_hit} {bfvd_annotation} {bfvd_keywords}".lower()
            for ref in djr_references:
                if ref in hit_text:
                    if bfvd_score > djr_score:
                        djr_score = bfvd_score
                        djr_hit = bfvd_hit
                    break

            # Check if hit matches SJR references
            for ref in sjr_references:
                if ref in hit_text:
                    if bfvd_score > sjr_score:
                        sjr_score = bfvd_score
                        sjr_hit = bfvd_hit
                    break

            results[porf_id] = (djr_score, djr_hit, sjr_score, sjr_hit)

    return results


def load_foldseek_results(
    foldseek_path: Path | None,
) -> dict[str, tuple[bool, bool, str]]:
    """
    Load FoldSeek results and check for DJR/SJR structural hits.

    Returns:
        Dictionary mapping query IDs to (is_djr, is_sjr, top_hit).
    """
    if foldseek_path is None or not foldseek_path.exists():
        return {}

    results: dict[str, tuple[bool, bool, str]] = {}

    with foldseek_path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue

            query = parts[0]
            target = parts[1].lower()

            # Check against DJR PDB patterns
            is_djr = any(re.search(pattern, target) for pattern in DJR_PDB_PATTERNS)

            # Check against SJR PDB patterns
            is_sjr = any(re.search(pattern, target) for pattern in SJR_PDB_PATTERNS)

            # Only store first (best) hit per query
            if query not in results:
                results[query] = (is_djr, is_sjr, parts[1])

    return results


def load_sequences(
    sequences_path: Path, porf_ids: set[str] | None = None
) -> dict[str, str]:
    """Load protein sequences from FASTA file."""
    sequences: dict[str, str] = {}
    for record in SeqIO.parse(sequences_path, "fasta"):
        seq_id = record.id
        if porf_ids is not None and seq_id not in porf_ids:
            continue
        sequences[seq_id] = str(record.seq)
    return sequences


def classify_proteins(
    marker_hits: dict[str, tuple[str, str]],
    sequences: dict[str, str],
    interproscan_domains: dict[str, int] | None = None,
    hmm_domain_counts: dict[str, int] | None = None,
    tmvec_results: dict[str, tuple[float, str, float, str]] | None = None,
    foldseek_results: dict[str, tuple[bool, bool, str]] | None = None,
) -> list[JellyRollClassification]:
    """
    Classify all MCP proteins using multi-signal approach.
    """
    classifications: list[JellyRollClassification] = []
    processed_ids: set[str] = set()

    interproscan_domains = interproscan_domains or {}
    hmm_domain_counts = hmm_domain_counts or {}
    tmvec_results = tmvec_results or {}
    foldseek_results = foldseek_results or {}

    for porf_id, (marker, status) in marker_hits.items():
        base_id = extract_base_porf_id(porf_id)
        if base_id in processed_ids:
            continue

        sequence = sequences.get(porf_id) or sequences.get(base_id)
        if not sequence:
            logger.debug(f"Sequence not found for {porf_id}")
            continue

        domain_length = len(sequence)
        estimated_length = estimate_full_protein_length(porf_id, domain_length)

        # Build classification signals
        signals = ClassificationSignals(length=estimated_length)

        # Signal 1: InterProScan domains
        for lookup_id in [porf_id, base_id]:
            if lookup_id in interproscan_domains:
                signals.jelly_roll_domains = max(
                    signals.jelly_roll_domains,
                    interproscan_domains[lookup_id]
                )
                signals.evidence_sources.append("interproscan")
                break

        # Signal 2: Multiple HMM domains
        if base_id in hmm_domain_counts:
            signals.hmm_domain_count = hmm_domain_counts[base_id]
            if signals.hmm_domain_count > 1:
                signals.evidence_sources.append("hmm_domains")

        # Signal 3: TMVec results
        for lookup_id in [porf_id, base_id]:
            if lookup_id in tmvec_results:
                djr_score, djr_hit, sjr_score, sjr_hit = tmvec_results[lookup_id]
                signals.tmvec_djr_score = djr_score
                signals.tmvec_djr_hit = djr_hit
                signals.tmvec_sjr_score = sjr_score
                signals.tmvec_sjr_hit = sjr_hit
                if djr_score > TMVEC_SCORE_THRESHOLD or sjr_score > TMVEC_SCORE_THRESHOLD:
                    signals.evidence_sources.append("tmvec")
                break

        # Signal 4: FoldSeek results (look for protein in Boltz model names)
        for query, (is_djr, is_sjr, top_hit) in foldseek_results.items():
            if base_id in query or porf_id in query:
                signals.foldseek_hit_is_djr = is_djr
                signals.foldseek_hit_is_sjr = is_sjr
                signals.foldseek_top_hit = top_hit
                if is_djr or is_sjr:
                    signals.evidence_sources.append("foldseek")
                break

        # Classify using multi-signal approach
        jelly_roll_type, confidence, evidence = classify_with_signals(signals, marker)

        classifications.append(
            JellyRollClassification(
                protein_id=porf_id,
                jelly_roll_type=jelly_roll_type,
                confidence=confidence,
                length=estimated_length,
                marker=marker,
                evidence=evidence,
                sequence=sequence,
            )
        )
        processed_ids.add(base_id)

    return classifications


def write_results(
    classifications: list[JellyRollClassification], output_path: Path
) -> None:
    """Write classification results to TSV file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["protein_id", "type", "confidence", "length", "marker", "evidence"])

        for c in classifications:
            writer.writerow(
                [
                    c.protein_id,
                    c.jelly_roll_type,
                    f"{c.confidence:.3f}",
                    c.length,
                    c.marker,
                    c.evidence,
                ]
            )

    logger.info(f"Wrote {len(classifications)} classifications to {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Classify MCP sequences as DJR or SJR using multi-signal approach.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python scripts/classify_jelly_roll.py \\
        --marker-hits tests/results/stena/phase1/marker_validation/validated_marker_hits.tsv \\
        --sequences tests/results/stena/phase1/marker_validation/hmm_hit_porfs.faa \\
        --output tests/results/stena/phase3_synthesis/virosync_jelly_roll_proteins.tsv \\
        --interproscan tests/results/stena/phase3/interproscan/interproscan_batch.tsv
        """,
    )
    parser.add_argument(
        "--marker-hits",
        type=Path,
        required=True,
        help="Path to validated_marker_hits.tsv from phase1 marker validation",
    )
    parser.add_argument(
        "--sequences",
        type=Path,
        required=True,
        help="Path to hmm_hit_porfs.faa containing HMM-hit pORF sequences",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output TSV file path for jelly roll classifications",
    )
    parser.add_argument(
        "--interproscan",
        type=Path,
        help="Path to InterProScan batch results TSV (optional, improves accuracy)",
    )
    parser.add_argument(
        "--tmvec",
        type=Path,
        help="Path to TMVec results TSV (optional, provides structural similarity)",
    )
    parser.add_argument(
        "--foldseek",
        type=Path,
        help="Path to FoldSeek results TSV (optional, from Boltz predictions)",
    )
    parser.add_argument(
        "--all-markers",
        action="store_true",
        help="Classify all markers, not just MCP-related ones",
    )
    parser.add_argument(
        "--include-sequences",
        action="store_true",
        help="Include full sequences in output (creates larger file)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Validate input files
    if not args.marker_hits.exists():
        raise FileNotFoundError(f"Marker hits file not found: {args.marker_hits}")
    if not args.sequences.exists():
        raise FileNotFoundError(f"Sequences file not found: {args.sequences}")

    # Load marker hits
    logger.info(f"Loading marker hits from {args.marker_hits}")
    mcp_only = not args.all_markers
    marker_hits = load_marker_hits(args.marker_hits, mcp_only=mcp_only)
    logger.info(f"Found {len(marker_hits)} {'MCP' if mcp_only else 'total'} marker hits")

    if not marker_hits:
        logger.warning("No marker hits found. Check if the marker-hits file contains MCP markers.")
        write_results([], args.output)
        return

    # Count HMM domains per protein
    hmm_domain_counts = count_hmm_domains_per_protein(marker_hits)
    multi_domain = sum(1 for c in hmm_domain_counts.values() if c > 1)
    logger.info(f"Found {multi_domain} proteins with multiple HMM domain hits")

    # Load optional signal sources
    interproscan_domains = None
    if args.interproscan:
        logger.info(f"Loading InterProScan results from {args.interproscan}")
        interproscan_domains = load_interproscan_domains(args.interproscan)
        djr_count = sum(1 for c in interproscan_domains.values() if c > 0)
        logger.info(f"Found {djr_count} proteins with jelly roll domain annotations")

    tmvec_results = None
    if args.tmvec:
        logger.info(f"Loading TMVec results from {args.tmvec}")
        tmvec_results = load_tmvec_results(args.tmvec)
        logger.info(f"Loaded TMVec scores for {len(tmvec_results)} proteins")

    foldseek_results = None
    if args.foldseek:
        logger.info(f"Loading FoldSeek results from {args.foldseek}")
        foldseek_results = load_foldseek_results(args.foldseek)
        djr_hits = sum(1 for r in foldseek_results.values() if r[0])
        sjr_hits = sum(1 for r in foldseek_results.values() if r[1])
        logger.info(f"Found {djr_hits} DJR and {sjr_hits} SJR structural hits")

    # Load sequences
    logger.info(f"Loading sequences from {args.sequences}")
    porf_ids = set(marker_hits.keys())
    sequences = load_sequences(args.sequences, porf_ids)
    logger.info(f"Loaded {len(sequences)} sequences")

    # Classify proteins
    logger.info("Classifying proteins with multi-signal approach...")
    classifications = classify_proteins(
        marker_hits,
        sequences,
        interproscan_domains=interproscan_domains,
        hmm_domain_counts=hmm_domain_counts,
        tmvec_results=tmvec_results,
        foldseek_results=foldseek_results,
    )
    logger.info(f"Classified {len(classifications)} proteins")

    # Print summary
    djr_count = sum(1 for c in classifications if c.jelly_roll_type == "DJR")
    sjr_count = sum(1 for c in classifications if c.jelly_roll_type == "SJR")
    unknown_count = sum(1 for c in classifications if c.jelly_roll_type == "UNKNOWN")

    logger.info("Classification summary:")
    logger.info(f"  DJR (Double Jelly Roll): {djr_count}")
    logger.info(f"  SJR (Single Jelly Roll): {sjr_count}")
    logger.info(f"  UNKNOWN: {unknown_count}")

    if classifications:
        avg_confidence = sum(c.confidence for c in classifications) / len(classifications)
        logger.info(f"  Average confidence: {avg_confidence:.3f}")

        # Evidence source breakdown
        evidence_counts: dict[str, int] = defaultdict(int)
        for c in classifications:
            source = c.evidence.split(":")[0] if ":" in c.evidence else c.evidence
            evidence_counts[source] += 1
        logger.info("  Evidence sources:")
        for source, count in sorted(evidence_counts.items(), key=lambda x: -x[1]):
            logger.info(f"    {source}: {count}")

        lengths = [c.length for c in classifications]
        logger.info(f"  Length range: {min(lengths)}-{max(lengths)} aa")
        logger.info(f"  Mean length: {sum(lengths)/len(lengths):.1f} aa")

    # Write results
    write_results(classifications, args.output)

    # Optionally write extended output with sequences
    if args.include_sequences:
        extended_output = args.output.with_suffix(".with_sequences.tsv")
        with extended_output.open("w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(
                ["protein_id", "type", "confidence", "length", "marker", "evidence", "sequence"]
            )
            for c in classifications:
                writer.writerow(
                    [
                        c.protein_id,
                        c.jelly_roll_type,
                        f"{c.confidence:.3f}",
                        c.length,
                        c.marker,
                        c.evidence,
                        c.sequence,
                    ]
                )
        logger.info(f"Wrote extended output with sequences to {extended_output}")


if __name__ == "__main__":
    main()
