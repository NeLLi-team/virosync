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

    # Load custom config from YAML
    config = ViroSyncConfig.from_yaml("custom_thresholds.yaml")

    # Export current config
    config.to_yaml("current_thresholds.yaml")
"""

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional
import json
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
    plddt_minimum: float = 50.0              # Minimum to consider

    # TM-score thresholds (0-1, higher = more similar structures)
    tm_score_significant: float = 0.5        # Significant structural match
    tm_score_highly_significant: float = 0.7  # Highly significant (likely homolog)
    tm_score_minimum: float = 0.3            # Minimum to report

    # E-value thresholds for FoldSeek hits
    evalue_significant: float = 1e-3
    evalue_highly_significant: float = 1e-10

    # ESM-2 sequence length limit (truncates longer sequences)
    esm2_max_length: int = 1022


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
    phylogenetic_viral_strong: float = 0.8   # Strong viral evidence
    phylogenetic_viral_likely: float = 0.6   # Likely viral
    phylogenetic_nonviral: float = 0.4       # Likely non-viral

    # Contamination detection
    contamination_score_threshold: float = 0.5


@dataclass
class EvidenceGraphThresholds:
    """Thresholds for evidence graph coherence scoring (Phase 3)."""

    # Coherence score thresholds
    coherence_strong: float = 0.7     # Strong support
    coherence_moderate: float = 0.5   # Moderate support
    coherence_weak: float = 0.3       # Weak support

    # Feature thresholds for evidence nodes
    kfd_significant: float = 0.3
    cub_deviation_significant: float = 0.1
    gc_deviation_significant: float = 0.1
    novelty_significant: float = 0.7
    viral_prob_strong: float = 0.9


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

    # Prefixes that indicate host/cellular sequences
    host_prefixes: list[str] = field(default_factory=lambda: [
        "EUK__",
        "BAC__",
        "ARC__",
        "CELLULAR__",
    ])

    # Default for unknown prefixes
    unknown_classification: str = "UNKNOWN"


@dataclass
class ViroSyncConfig:
    """Master configuration containing all thresholds.

    This is the main configuration class that aggregates all threshold
    categories. Use get_config() to get the default singleton instance.

    Example:
        config = ViroSyncConfig()
        config.structural.plddt_high_confidence = 75.0  # Customize
        config.to_yaml("my_config.yaml")
    """

    structural: StructuralThresholds = field(default_factory=StructuralThresholds)
    evidence: EvidenceThresholds = field(default_factory=EvidenceThresholds)
    evidence_graph: EvidenceGraphThresholds = field(default_factory=EvidenceGraphThresholds)
    database: DatabaseConfig = field(default_factory=DatabaseConfig)

    @classmethod
    def from_yaml(cls, path: Path) -> "ViroSyncConfig":
        """Load configuration from YAML file.

        Args:
            path: Path to YAML configuration file

        Returns:
            ViroSyncConfig instance with loaded values
        """
        try:
            import yaml
        except ImportError:
            raise ImportError("PyYAML required for YAML config. Install with: pip install pyyaml")

        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Config file not found: {path}")

        with open(path) as f:
            data = yaml.safe_load(f)

        return cls._from_dict(data)

    @classmethod
    def from_json(cls, path: Path) -> "ViroSyncConfig":
        """Load configuration from JSON file.

        Args:
            path: Path to JSON configuration file

        Returns:
            ViroSyncConfig instance with loaded values
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Config file not found: {path}")

        with open(path) as f:
            data = json.load(f)

        return cls._from_dict(data)

    @classmethod
    def _from_dict(cls, data: dict) -> "ViroSyncConfig":
        """Create config from dictionary."""
        config = cls()

        # Map section names to dataclass types
        sections = {
            "structural": StructuralThresholds,
            "evidence": EvidenceThresholds,
            "evidence_graph": EvidenceGraphThresholds,
            "database": DatabaseConfig,
        }

        for section_name, section_cls in sections.items():
            if section_name in data:
                section_data = data[section_name]
                current = getattr(config, section_name)
                for key, value in section_data.items():
                    if hasattr(current, key):
                        setattr(current, key, value)
                    else:
                        logger.warning(f"Unknown config key: {section_name}.{key}")

        return config

    def to_yaml(self, path: Path) -> None:
        """Export configuration to YAML file.

        Args:
            path: Path for output YAML file
        """
        try:
            import yaml
        except ImportError:
            raise ImportError("PyYAML required for YAML config. Install with: pip install pyyaml")

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        with open(path, "w") as f:
            yaml.dump(asdict(self), f, default_flow_style=False, sort_keys=False)

        logger.info(f"Exported config to {path}")

    def to_json(self, path: Path) -> None:
        """Export configuration to JSON file.

        Args:
            path: Path for output JSON file
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        with open(path, "w") as f:
            json.dump(asdict(self), f, indent=2)

        logger.info(f"Exported config to {path}")

    def to_dict(self) -> dict:
        """Convert configuration to dictionary."""
        return asdict(self)


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


def set_config(config: ViroSyncConfig) -> None:
    """Set the global configuration singleton.

    Args:
        config: ViroSyncConfig instance to use globally

    Example:
        from virosync.config.thresholds import set_config, ViroSyncConfig

        custom_config = ViroSyncConfig.from_yaml("my_thresholds.yaml")
        set_config(custom_config)
    """
    global _global_config
    _global_config = config
    logger.info("Global ViroSync config updated")


def reset_config() -> None:
    """Reset global configuration to defaults."""
    global _global_config
    _global_config = ViroSyncConfig()
    logger.info("Global ViroSync config reset to defaults")
