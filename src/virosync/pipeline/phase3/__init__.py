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

import logging

logger = logging.getLogger(__name__)

from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    EvidenceSynthesizerConfig,
    VerificationResult,
    VerificationStatus,
    synthesize_evidence,
    calculate_eve_confidence,
)

from virosync.pipeline.phase3.evidence_graph import (
    EvidenceType,
    WindowEvidence,
    EvidenceProfile,
    EvidenceCorrelationGraph,
    CoherenceAnalysis,
    build_evidence_profile,
    analyze_eve_coherence,
)

_STRUCTURAL_HOMOLOGY_AVAILABLE = True
try:
    from virosync.pipeline.phase3.structural_homology import (
        StructurePrediction,
        FoldSeekHit,
        TMvecHit,
        StructuralHomologyResult,
        TMvec2Predictor,
        ESM2Embedder,
        FoldSeekSearcher,
        BoltzFoldSeekAnalyzer,
    )
except Exception as exc:  # pragma: no cover - optional dependency path
    _STRUCTURAL_HOMOLOGY_AVAILABLE = False
    logger.debug("Structural homology imports unavailable: %s", exc)
    StructurePrediction = None
    FoldSeekHit = None
    TMvecHit = None
    StructuralHomologyResult = None
    TMvec2Predictor = None
    ESM2Embedder = None
    FoldSeekSearcher = None
    BoltzFoldSeekAnalyzer = None

from virosync.pipeline.phase3.phylogenetic_validation import (
    PhylogeneticValidator,
    PhylogeneticValidationResult,
    PhylogeneticVerdict,
    GVClassValidation,
    DiamondValidation,
    validate_eve_regions,
)

from virosync.pipeline.phase3.output_generator import (
    OutputGenerator,
    generate_outputs,
)

from virosync.pipeline.phase3.gene_taxonomy import (
    run_gene_taxonomy_diamond,
    GeneTaxonomy,
    classify_gene_taxonomy,
    extract_prefix,
)

__all__ = [
    # Evidence Synthesizer
    "EvidenceSynthesizer",
    "EvidenceSynthesizerConfig",
    "VerificationResult",
    "VerificationStatus",
    "synthesize_evidence",
    # Evidence Graph
    "EvidenceType",
    "WindowEvidence",
    "EvidenceProfile",
    "EvidenceCorrelationGraph",
    "CoherenceAnalysis",
    "build_evidence_profile",
    "analyze_eve_coherence",
    # Phylogenetic Validation
    "PhylogeneticValidator",
    "PhylogeneticValidationResult",
    "PhylogeneticVerdict",
    "GVClassValidation",
    "DiamondValidation",
    "validate_eve_regions",
    # Output Generator
    "OutputGenerator",
    "generate_outputs",
    # Gene Taxonomy (Step 9)
    "run_gene_taxonomy_diamond",
    "GeneTaxonomy",
    "classify_gene_taxonomy",
    "extract_prefix",
    # Confidence Calculation (Step 10)
    "calculate_eve_confidence",
]

if _STRUCTURAL_HOMOLOGY_AVAILABLE:
    __all__.extend(
        [
            "StructurePrediction",
            "FoldSeekHit",
            "TMvecHit",
            "StructuralHomologyResult",
            "TMvec2Predictor",
            "ESM2Embedder",
            "FoldSeekSearcher",
            "BoltzFoldSeekAnalyzer",
        ]
    )
