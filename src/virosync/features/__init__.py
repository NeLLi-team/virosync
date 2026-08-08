"""
ViroSync Features Module.

Provides sequence-based features for EVE detection and characterization.

Modules:
- compositional: K-mer frequency deviation, GC content
"""

from .compositional import (
    BackgroundModel,
    WindowFeatures,
    calculate_kfd,
    calculate_gc_content,
    calculate_gc_deviation,
)

__all__ = [
    "BackgroundModel",
    "WindowFeatures",
    "calculate_kfd",
    "calculate_gc_content",
    "calculate_gc_deviation",
]
