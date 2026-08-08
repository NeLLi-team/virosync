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

import numpy as np
from scipy.spatial.distance import jensenshannon

logger = logging.getLogger(__name__)


@dataclass
class BackgroundModel:
    """Pre-computed background statistics for a genome."""

    kmer_freqs: dict[str, float]  # k-mer frequencies
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

        # Compute GC content
        gc_count = sequence.count('G') + sequence.count('C')
        valid_bases = len(sequence) - sequence.count('N')
        gc_content = gc_count / valid_bases if valid_bases > 0 else 0.5

        return cls(
            kmer_freqs=kmer_freqs,
            gc_content=gc_content,
            k=k,
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
