"""Phase 0: Pre-processing and Conceptual Proteome Generation.

This phase prepares the input genome for EVE detection:
1. Masks low-complexity and repetitive regions
2. Generates conceptual proteome via 6-frame translation
3. Filters pORFs by minimum length threshold

Output: Multi-FASTA of potential ORFs with encoded coordinates.
"""

from .masking import (
    MaskedRegion,
    MaskingBackendError,
    MaskingResult,
    apply_mask,
    identify_repeats,
    load_masking_result,
    mask_genome_pipeline,
    parse_repeatmasker_output,
    parse_trf_output,
    quick_mask,
    run_repeatmasker,
    run_trf,
    validate_masking_result,
    write_masking_status,
)
from .prodigal import (
    GenePrediction,
    load_gene_predictions,
    parse_prodigal_header,
    run_prodigal_genome,
)
from .translation import (
    PORF,
)

__all__ = [
    # Masking
    "MaskedRegion",
    "MaskingBackendError",
    "MaskingResult",
    "identify_repeats",
    "load_masking_result",
    "mask_genome_pipeline",
    "quick_mask",
    "run_trf",
    "run_repeatmasker",
    "parse_trf_output",
    "parse_repeatmasker_output",
    "apply_mask",
    "validate_masking_result",
    "write_masking_status",
    # Prodigal
    "GenePrediction",
    "run_prodigal_genome",
    "load_gene_predictions",
    "parse_prodigal_header",
    # Translation
    "PORF",
]
