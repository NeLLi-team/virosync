"""Phase 1: Unified Seeding.

HHG Seeding (the only active path):
    Uses HMM search against viral hallmark genes to find anchor pORFs,
    then scores neighborhoods by viral signal density. Marker validation
    refines the anchor set; region assembly turns validated anchors into
    candidate loci for Phase 2 boundary refinement.

Legacy novelty and compositional seeding paths are retired from the active
workflow; their compatibility modules are not part of the public Phase 1
surface.
"""

from .hhg_seeding import (
    Anchor,
    HHGSeed,
    HMMHit,
    calculate_neighbor_scores,
    form_seeds,
    hhg_seeding_pipeline,
    identify_anchors,
    load_hmm_profiles,
    run_hmmsearch,
)
from .marker_validation import (
    NovelMarkerCriteria,
    ValidatedMarkerHit,
    ValidationStatus,
    extract_hmm_hit_sequences,
    filter_validated_markers,
    run_diamond_on_hmm_hits,
    validate_hmm_hit,
)
from .region_assembly import (
    CandidateRegion,
    assemble_candidate_regions,
    initial_clustering,
    iterative_extension,
    merge_overlapping_regions,
)
from .seed_merger import (
    MergedSeed,
)

__all__ = [
    # HHG Seeding
    "HMMHit",
    "Anchor",
    "HHGSeed",
    "load_hmm_profiles",
    "run_hmmsearch",
    "identify_anchors",
    "calculate_neighbor_scores",
    "form_seeds",
    "hhg_seeding_pipeline",
    # Seed Merger (MergedSeed dataclass — used by Phase 2)
    "MergedSeed",
    # Marker Validation (Pipeline Rewrite)
    "ValidationStatus",
    "NovelMarkerCriteria",
    "ValidatedMarkerHit",
    "validate_hmm_hit",
    "extract_hmm_hit_sequences",
    "run_diamond_on_hmm_hits",
    "filter_validated_markers",
    # Region Assembly (Pipeline Rewrite)
    "CandidateRegion",
    "assemble_candidate_regions",
    "initial_clustering",
    "iterative_extension",
    "merge_overlapping_regions",
]
