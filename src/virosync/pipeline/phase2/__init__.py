"""Phase 2: Boundary Refinement.

The active workflow uses gene extension, batched Diamond taxonomy, and host
taxonomy trimming.
"""

from virosync.pipeline.phase2.boundary_diamond import (
    CELLULAR_PREFIXES,
    VIRAL_PREFIXES,
    BoundaryDiamondConfig,
    ControlStats,
    GeneTaxonomy,
    GenomeDiamondQuery,
    build_proteome_index,
    classify_all_porfs,
    collect_query_proteins,
    compute_control_stats,
    extract_prefix,
    filter_taxonomy_to_boundary,
    run_batched_diamond,
    run_diamond_chunked,
    sample_control_porfs_genome_wide,
)
from virosync.pipeline.phase2.boundary_refiner import (
    RefinedBoundary,
)

__all__ = [
    # Boundary Refinement
    "RefinedBoundary",
    # Boundary Diamond
    "BoundaryDiamondConfig",
    "GenomeDiamondQuery",
    "GeneTaxonomy",
    "ControlStats",
    "collect_query_proteins",
    "sample_control_porfs_genome_wide",
    "run_batched_diamond",
    "run_diamond_chunked",
    "classify_all_porfs",
    "extract_prefix",
    "filter_taxonomy_to_boundary",
    "compute_control_stats",
    "build_proteome_index",
    "VIRAL_PREFIXES",
    "CELLULAR_PREFIXES",
]
