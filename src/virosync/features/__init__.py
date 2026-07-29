"""
ViroSync Features Module.

Provides sequence-based features for EVE detection and characterization.

Modules:
- compositional: K-mer frequency deviation, codon usage bias, GC content
"""

from .compositional import (
    BackgroundModel,
    WindowFeatures,
    calculate_kfd,
    calculate_cub_deviation,
    calculate_gc_content,
    calculate_gc_deviation,
    calculate_window_features,
    scan_genome_windows,
    identify_compositional_anomalies,
)

__all__ = [
    "BackgroundModel",
    "WindowFeatures",
    "calculate_kfd",
    "calculate_cub_deviation",
    "calculate_gc_content",
    "calculate_gc_deviation",
    "calculate_window_features",
    "scan_genome_windows",
    "identify_compositional_anomalies",
]
