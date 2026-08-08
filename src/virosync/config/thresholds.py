"""
Centralized Threshold Configuration for ViroSync.

This module provides a single source of truth for all pipeline thresholds,
enabling easy tuning without modifying source code across multiple files.

Usage:
    from virosync.config.thresholds import get_config, ViroSyncConfig

    # Get default config
    config = get_config()

    # Use thresholds
    if plddt > config.structural.plddt_high_confidence:
        ...
"""

from dataclasses import dataclass, field
from typing import Optional
import logging

logger = logging.getLogger(__name__)


@dataclass
class StructuralThresholds:
    """Thresholds for structural homology analysis (Phase 3).

    These control structure prediction quality filtering
    and FoldSeek/TMvec structural similarity cutoffs.
    """

    # pLDDT (predicted Local Distance Difference Test) thresholds
    # Higher = more confident structure prediction
    plddt_high_confidence: float = 70.0      # Accept structure as confident
    plddt_very_high_confidence: float = 85.0  # Very high confidence

    # TM-score thresholds (0-1, higher = more similar structures)
    tm_score_significant: float = 0.5        # Significant structural match
    tm_score_highly_significant: float = 0.7  # Highly significant (likely homolog)

    # E-value thresholds for FoldSeek hits
    evalue_significant: float = 1e-3
    evalue_highly_significant: float = 1e-10


@dataclass
class EvidenceThresholds:
    """Thresholds for evidence synthesis (Phase 3).

    All EVE candidates receive full analysis. The tier thresholds
    classify final confidence scores into HIGH, MEDIUM, or LOW tiers.
    """

    # Confidence tier thresholds (for output classification)
    # HIGH: confidence >= high_tier_threshold
    # MEDIUM: low_tier_threshold <= confidence < high_tier_threshold
    # LOW: confidence < low_tier_threshold
    high_tier_threshold: float = 0.7
    low_tier_threshold: float = 0.2

    # Phylogenetic evidence thresholds
    phylogenetic_rejection: float = 0.3  # Below this = reject as contamination


@dataclass
class DatabaseConfig:
    """Configuration for database naming conventions (Issue #7).

    Maps database prefixes to taxonomic classifications.
    Allows customization for different database sources.
    """

    # Prefixes that indicate viral sequences
    viral_prefixes: list[str] = field(default_factory=lambda: [
        "NCLDV__",
        "MIRUS__",
        "GVMAG__",
        "PHAGE__",
        "CRESS__",
        "VIRUS__",
    ])


@dataclass
class ViroSyncConfig:
    """Master configuration containing all thresholds.

    This is the main configuration class that aggregates all threshold
    categories. Use get_config() to get the default singleton instance.

    Example:
        config = ViroSyncConfig()
        config.structural.plddt_high_confidence = 75.0  # Customize
    """

    structural: StructuralThresholds = field(default_factory=StructuralThresholds)
    evidence: EvidenceThresholds = field(default_factory=EvidenceThresholds)
    database: DatabaseConfig = field(default_factory=DatabaseConfig)


# Global singleton instance
_global_config: Optional[ViroSyncConfig] = None


def get_config() -> ViroSyncConfig:
    """Get the global configuration singleton.

    Returns:
        ViroSyncConfig instance (creates default if not set)

    Example:
        from virosync.config.thresholds import get_config

        config = get_config()
        if score > config.structural.tm_score_significant:
            print("Significant structural match!")
    """
    global _global_config
    if _global_config is None:
        _global_config = ViroSyncConfig()
    return _global_config
