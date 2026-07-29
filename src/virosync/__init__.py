"""
ViroSync: Framework for discovering ancient viral elements in eukaryotic genomes.

A bioinformatics framework for discovering candidate giant endogenous viral elements
(EVEs) in eukaryotic genomes using HMM-gated marker discovery, taxonomy-guided
boundary refinement, and multi-evidence confidence scoring.
"""

__version__ = "1.0.0"
__author__ = "ViroSync Development Team"

from pathlib import Path

# Project root directory
PROJECT_ROOT = Path(__file__).parent.parent.parent

# Default paths
DATA_DIR = PROJECT_ROOT / "data"
CONFIGS_DIR = PROJECT_ROOT / "config"
MODELS_DIR = PROJECT_ROOT / "models"
