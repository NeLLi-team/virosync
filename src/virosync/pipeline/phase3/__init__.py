"""
Phase 3: Evidence Synthesis & Verification.

This module implements the verification pipeline that scores EVE candidates
using marker, taxonomy, composition, optional structural/domain, and optional
phylogenetic evidence.

Tie-breaker modules:
- Evidence Correlation Graph (Coherence Score)
- Structural Homology (Boltz + FoldSeek)
- Phylogenetic Validation (GVClass + Diamond BLASTp) - optional

When enabled, phylogenetic validation provides a final check, running GVClass
classification and Diamond BLASTp against reference proteomes to confirm or
reject candidate regions.
"""

from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    EvidenceSynthesizerConfig,
)

from virosync.pipeline.phase3.output_generator import OutputGenerator

__all__ = [
    "EvidenceSynthesizer",
    "EvidenceSynthesizerConfig",
    "OutputGenerator",
]
