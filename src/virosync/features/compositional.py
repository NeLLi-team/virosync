"""
Compositional features for EVE detection.

This module provides sequence-based features that can identify viral regions
WITHOUT relying on HMM homology searches. These features detect compositional
anomalies that distinguish viral DNA from host DNA.

Key features:
- KFD (K-mer Frequency Deviation): Jensen-Shannon divergence of k-mer frequencies
- CUB (Codon Usage Bias): Deviation from host codon usage patterns
- GC Content: Local GC percentage
- pORF Density: Gene density per window

These features are critical for:
1. Detecting highly divergent/novel viruses with no HMM hits
2. CRF boundary detection (Tier 1 and Tier 2)
3. Evidence coherence scoring
"""

import logging
from collections import Counter
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.spatial.distance import jensenshannon

logger = logging.getLogger(__name__)


# Standard genetic code codon table
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# All possible codons (excluding stop codons for CUB)
ALL_CODONS = [c for c in CODON_TABLE.keys() if CODON_TABLE[c] != '*']


@dataclass
class BackgroundModel:
    """Pre-computed background statistics for a genome."""

    kmer_freqs: dict[str, float]  # k-mer frequencies
    codon_freqs: dict[str, float]  # Codon frequencies
    gc_content: float  # Overall GC content
    k: int  # k-mer size used

    @classmethod
    def from_sequence(cls, sequence: str, k: int = 4) -> "BackgroundModel":
        """
        Compute background model from a genome sequence.

        Args:
            sequence: Full genome sequence (uppercase, may contain N)
            k: K-mer size for KFD calculation

        Returns:
            BackgroundModel with pre-computed frequencies
        """
        sequence = sequence.upper()

        # Compute k-mer frequencies
        kmer_counts = Counter()
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if 'N' not in kmer:
                kmer_counts[kmer] += 1

        total_kmers = sum(kmer_counts.values())
        kmer_freqs = {kmer: count / total_kmers for kmer, count in kmer_counts.items()} if total_kmers > 0 else {}

        # Compute codon frequencies from all frames
        codon_counts = Counter()
        for frame in range(3):
            for i in range(frame, len(sequence) - 2, 3):
                codon = sequence[i:i+3]
                if 'N' not in codon and codon in CODON_TABLE:
                    codon_counts[codon] += 1

        total_codons = sum(codon_counts.values())
        codon_freqs = {codon: count / total_codons for codon, count in codon_counts.items()} if total_codons > 0 else {}

        # Compute GC content
        gc_count = sequence.count('G') + sequence.count('C')
        valid_bases = len(sequence) - sequence.count('N')
        gc_content = gc_count / valid_bases if valid_bases > 0 else 0.5

        return cls(
            kmer_freqs=kmer_freqs,
            codon_freqs=codon_freqs,
            gc_content=gc_content,
            k=k,
        )

    def save(self, path: Path) -> None:
        """Save background model to JSON."""
        import json
        data = {
            "kmer_freqs": self.kmer_freqs,
            "codon_freqs": self.codon_freqs,
            "gc_content": self.gc_content,
            "k": self.k,
        }
        with open(path, 'w') as f:
            json.dump(data, f)

    @classmethod
    def load(cls, path: Path) -> "BackgroundModel":
        """Load background model from JSON."""
        import json
        with open(path) as f:
            data = json.load(f)
        return cls(
            kmer_freqs=data["kmer_freqs"],
            codon_freqs=data["codon_freqs"],
            gc_content=data["gc_content"],
            k=data["k"],
        )


def calculate_kfd(
    sequence: str,
    background_freqs: dict[str, float],
    k: int = 4,
) -> float:
    """
    Calculate K-mer Frequency Deviation using Jensen-Shannon divergence.

    KFD measures how different the k-mer composition of a local window is
    from the genome-wide background. Viral insertions often have distinct
    k-mer signatures that differ from the host genome.

    Args:
        sequence: Local DNA sequence window (e.g., 1000 bp)
        background_freqs: Pre-computed k-mer frequencies for the host genome
        k: K-mer length (default: 4, giving 256 possible k-mers)

    Returns:
        Jensen-Shannon divergence score (0-1). Higher = more divergent from host.
        Returns 0.0 if sequence is too short or has too many Ns.

    Example:
        >>> bg = BackgroundModel.from_sequence(genome_seq)
        >>> kfd = calculate_kfd(window_seq, bg.kmer_freqs, k=4)
        >>> if kfd > 0.15:
        ...     print("Compositionally divergent region!")
    """
    sequence = sequence.upper()

    # Avoid calculation on sequences with too many unknown bases
    if len(sequence) < k:
        return 0.0

    n_count = sequence.count('N')
    if n_count / len(sequence) > 0.5:
        return 0.0

    # Generate all possible k-mers for the alphabet
    alphabet = 'ATGC'
    all_kmers = [''.join(p) for p in product(alphabet, repeat=k)]
    kmer_map = {kmer: i for i, kmer in enumerate(all_kmers)}

    # Calculate k-mer frequencies for the local sequence window
    window_counts = Counter()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if 'N' not in kmer:
            window_counts[kmer] += 1

    total_window_kmers = sum(window_counts.values())
    if total_window_kmers == 0:
        return 0.0

    window_freqs = {kmer: count / total_window_kmers for kmer, count in window_counts.items()}

    # Create probability distribution vectors in a consistent order
    p_vec = np.zeros(len(all_kmers), dtype=float)
    q_vec = np.zeros(len(all_kmers), dtype=float)

    for kmer, index in kmer_map.items():
        p_vec[index] = window_freqs.get(kmer, 0.0)
        q_vec[index] = background_freqs.get(kmer, 0.0)

    # Add small pseudocount to avoid division by zero
    epsilon = 1e-10
    p_vec = p_vec + epsilon
    q_vec = q_vec + epsilon

    # Normalize to ensure valid probability distributions
    p_vec = p_vec / p_vec.sum()
    q_vec = q_vec / q_vec.sum()

    # Calculate Jensen-Shannon divergence
    js_divergence = jensenshannon(p_vec, q_vec)

    return float(js_divergence) if not np.isnan(js_divergence) else 0.0


def calculate_cub_deviation(
    coding_sequences: list[str],
    host_codon_freqs: dict[str, float],
) -> float:
    """
    Calculate Codon Usage Bias deviation from host reference.

    CUB measures how different the codon usage of genes in a region is
    from the typical host codon usage. Viral genes often have distinct
    codon preferences reflecting their evolutionary history.

    Args:
        coding_sequences: List of coding sequences (nucleotide, in-frame)
        host_codon_freqs: Reference codon frequencies for the host

    Returns:
        Mean absolute deviation of codon frequencies (0-1).
        Higher values indicate more divergent codon usage.
        Returns 0.0 if insufficient coding sequence.

    Example:
        >>> bg = BackgroundModel.from_sequence(genome_seq)
        >>> porf_seqs = [porf.nt_sequence for porf in window_porfs]
        >>> cub = calculate_cub_deviation(porf_seqs, bg.codon_freqs)
        >>> if cub > 0.05:
        ...     print("Unusual codon usage detected!")
    """
    if not coding_sequences:
        return 0.0

    # Combine all coding sequences and count codons
    all_codons_str = ""
    for seq in coding_sequences:
        seq = seq.upper()
        # Only use complete codons
        trim_len = len(seq) - (len(seq) % 3)
        all_codons_str += seq[:trim_len]

    if len(all_codons_str) < 30:  # Minimum threshold for meaningful calculation
        return 0.0

    # Count observed codons (excluding stop codons and Ns)
    observed_counts = Counter()
    for i in range(0, len(all_codons_str), 3):
        codon = all_codons_str[i:i+3]
        if 'N' not in codon and codon in ALL_CODONS:
            observed_counts[codon] += 1

    total_observed_codons = sum(observed_counts.values())
    if total_observed_codons == 0:
        return 0.0

    # Convert observed counts to frequencies
    observed_freqs = {codon: count / total_observed_codons for codon, count in observed_counts.items()}

    # Compare observed frequencies to host reference frequencies
    deviations = []
    for codon in ALL_CODONS:
        observed = observed_freqs.get(codon, 0.0)
        expected = host_codon_freqs.get(codon, 0.0)
        deviations.append(abs(observed - expected))

    # Return the mean absolute deviation
    return sum(deviations) / len(deviations) if deviations else 0.0


def calculate_gc_content(sequence: str) -> float:
    """
    Calculate GC content of a sequence.

    Args:
        sequence: DNA sequence (may contain N)

    Returns:
        GC fraction (0-1), or 0.5 if sequence has no valid bases
    """
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    valid_bases = len(sequence) - sequence.count('N')
    return gc_count / valid_bases if valid_bases > 0 else 0.5


def calculate_gc_deviation(sequence: str, background_gc: float) -> float:
    """
    Calculate deviation of local GC content from background.

    Args:
        sequence: Local DNA sequence
        background_gc: Background GC content (0-1)

    Returns:
        Absolute deviation from background GC (0-0.5)
    """
    local_gc = calculate_gc_content(sequence)
    return abs(local_gc - background_gc)


@dataclass
class WindowFeatures:
    """Compositional features for a genomic window."""

    scaffold: str
    start: int
    end: int
    kfd: float  # K-mer frequency deviation
    cub: float  # Codon usage bias deviation
    gc_content: float  # Local GC content
    gc_deviation: float  # Deviation from background GC
    porf_density: float  # Genes per kb
    porf_count: int  # Number of genes in window

    @property
    def composite_score(self) -> float:
        """
        Calculate composite compositional anomaly score.

        Combines multiple signals into a single score indicating
        how "non-host-like" this window appears.

        Returns:
            Score from 0 to 1, where higher = more anomalous
        """
        # Weighted combination of features
        # KFD is most informative, followed by CUB and GC deviation
        score = (
            0.4 * min(self.kfd / 0.3, 1.0) +  # Normalize KFD (typical range 0-0.3)
            0.3 * min(self.cub / 0.1, 1.0) +  # Normalize CUB (typical range 0-0.1)
            0.2 * min(self.gc_deviation / 0.15, 1.0) +  # Normalize GC dev
            0.1 * min(self.porf_density / 2.0, 1.0)  # High pORF density is suspicious
        )
        return min(score, 1.0)


def calculate_window_features(
    sequence: str,
    scaffold: str,
    start: int,
    end: int,
    background: BackgroundModel,
    porf_sequences: Optional[list[str]] = None,
) -> WindowFeatures:
    """
    Calculate all compositional features for a genomic window.

    Args:
        sequence: DNA sequence of the window
        scaffold: Scaffold/chromosome name
        start: Window start coordinate
        end: Window end coordinate
        background: Pre-computed background model
        porf_sequences: Optional list of gene coding sequences in this window

    Returns:
        WindowFeatures with all computed metrics
    """
    window_size_kb = (end - start) / 1000

    kfd = calculate_kfd(sequence, background.kmer_freqs, k=background.k)
    gc = calculate_gc_content(sequence)
    gc_dev = abs(gc - background.gc_content)

    if porf_sequences:
        cub = calculate_cub_deviation(porf_sequences, background.codon_freqs)
        porf_count = len(porf_sequences)
        porf_density = porf_count / window_size_kb if window_size_kb > 0 else 0
    else:
        cub = 0.0
        porf_count = 0
        porf_density = 0.0

    return WindowFeatures(
        scaffold=scaffold,
        start=start,
        end=end,
        kfd=kfd,
        cub=cub,
        gc_content=gc,
        gc_deviation=gc_dev,
        porf_density=porf_density,
        porf_count=porf_count,
    )


def scan_genome_windows(
    genome_fasta: Path,
    background: BackgroundModel,
    window_size: int = 1000,
    step_size: int = 500,
) -> list[WindowFeatures]:
    """
    Scan genome with sliding windows and compute compositional features.

    This provides a genome-wide compositional profile that can be used
    for seeding (identifying anomalous regions) without HMM search.

    Args:
        genome_fasta: Path to genome FASTA file
        background: Pre-computed background model
        window_size: Window size in bp (default: 1000)
        step_size: Step size for sliding window (default: 500)

    Returns:
        List of WindowFeatures for all windows
    """
    from Bio import SeqIO

    all_features = []

    for record in SeqIO.parse(genome_fasta, "fasta"):
        sequence = str(record.seq).upper()
        scaffold = record.id

        for start in range(0, len(sequence) - window_size + 1, step_size):
            end = start + window_size
            window_seq = sequence[start:end]

            features = calculate_window_features(
                sequence=window_seq,
                scaffold=scaffold,
                start=start,
                end=end,
                background=background,
            )
            all_features.append(features)

    logger.info(f"Computed compositional features for {len(all_features)} windows")
    return all_features


def identify_compositional_anomalies(
    features: list[WindowFeatures],
    kfd_threshold: float = 0.15,
    composite_threshold: float = 0.5,
    min_consecutive: int = 3,
) -> list[tuple[str, int, int, float]]:
    """
    Identify regions with compositional anomalies that may indicate viral content.

    This can serve as an alternative/complementary seeding strategy when
    HMM search fails to find hallmark genes.

    Args:
        features: List of WindowFeatures from genome scan
        kfd_threshold: Minimum KFD for anomaly detection
        composite_threshold: Minimum composite score for anomaly
        min_consecutive: Minimum consecutive anomalous windows to call a region

    Returns:
        List of (scaffold, start, end, max_score) for anomalous regions
    """
    # Group by scaffold
    by_scaffold: dict[str, list[WindowFeatures]] = {}
    for f in features:
        if f.scaffold not in by_scaffold:
            by_scaffold[f.scaffold] = []
        by_scaffold[f.scaffold].append(f)

    anomalous_regions = []

    for scaffold, scaffold_features in by_scaffold.items():
        # Sort by position
        scaffold_features.sort(key=lambda x: x.start)

        # Find runs of anomalous windows
        current_run = []

        for feat in scaffold_features:
            is_anomalous = (
                feat.kfd >= kfd_threshold or
                feat.composite_score >= composite_threshold
            )

            if is_anomalous:
                current_run.append(feat)
            else:
                if len(current_run) >= min_consecutive:
                    start = current_run[0].start
                    end = current_run[-1].end
                    max_score = max(f.composite_score for f in current_run)
                    anomalous_regions.append((scaffold, start, end, max_score))
                current_run = []

        # Handle run at end
        if len(current_run) >= min_consecutive:
            start = current_run[0].start
            end = current_run[-1].end
            max_score = max(f.composite_score for f in current_run)
            anomalous_regions.append((scaffold, start, end, max_score))

    logger.info(f"Identified {len(anomalous_regions)} compositionally anomalous regions")
    return anomalous_regions
