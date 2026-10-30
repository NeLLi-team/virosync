"""Boundary records and public candidate identities for seed-preserving growth."""

from __future__ import annotations

import hashlib
import logging
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING

import numpy as np

from virosync.pipeline.phase1.seed_merger import ConfirmedSeed

if TYPE_CHECKING:
    from virosync.pipeline.phase2.repeat_evidence import RepeatEvidence

logger = logging.getLogger(__name__)

TIR_STATUSES = frozenset(
    {
        "ambiguous",
        "detected",
        "no_marker_anchor",
        "not_assessed",
        "not_detected",
        "shared_pair",
        "taxonomy_incomplete",
    }
)


@dataclass
class RefinedBoundary:
    """Refined EVE boundary from Phase 2 boundary processing.

    Contains precise coordinates and confidence metrics.
    """

    scaffold: str
    start: int
    end: int
    candidate_id: str = field(default="", kw_only=True)
    confirmed_seeds: tuple[ConfirmedSeed, ...] = field(default=(), kw_only=True)
    left_stop_reason: str = field(default="", kw_only=True)
    right_stop_reason: str = field(default="", kw_only=True)
    left_stop: int | None = field(default=None, kw_only=True)
    right_stop: int | None = field(default=None, kw_only=True)
    host_stop_available: bool = field(default=False, kw_only=True)

    # Original seed info
    seed_id: str = ""  # Stable ID from MergedSeed for boundary-to-seed mapping
    original_start: int = 0
    original_end: int = 0

    # Terminal inverted-repeat evidence. Coordinates are 0-based, half-open.
    # ``pre_tir_*`` preserves the heuristic boundary that the accepted repeat
    # pair replaced, while ``tir_scan_*`` records the searched interval.
    pre_tir_start: int | None = None
    pre_tir_end: int | None = None
    tir_present: bool = False
    tir_status: str = "not_assessed"
    tir_left_start: int | None = None
    tir_left_end: int | None = None
    tir_right_start: int | None = None
    tir_right_end: int | None = None
    tir_identity: float = 0.0
    tir_boundary_override: bool = False
    tir_alignment_capped: bool = False
    tir_alignment_length: int = 0
    tir_candidate_count: int = 0
    tir_scan_start: int | None = None
    tir_scan_end: int | None = None
    tsd_sequence: str = ""
    repeat_evidence: RepeatEvidence | None = None

    # Seed evidence metadata (passed from Phase 1)
    seed_sources: list[str] = field(default_factory=list)  # ["hhg", "novelty", "compositional"]
    seed_confidence: str = "low"  # "high", "medium", "low"
    seed_hhg_score: float = 0.0
    seed_novelty_score: float = 0.0
    seed_compositional_score: float = 0.0
    seed_has_mcp: bool = False
    # Region classification from Phase 1 (NCLDV, VP, PLV, MIRUS, MIXED, UNKNOWN)
    predicted_family: str = ""
    region_classification_ncldv_markers: int = 0
    region_classification_vp_plv_markers: int = 0
    region_classification_mirus_markers: int = 0

    # Boundary confidence
    confidence: float = 0.0
    posterior_probability: float = 0.0

    # Anatomical breakdown
    core_viral_start: int | None = None
    core_viral_end: int | None = None
    flank_5_start: int | None = None
    flank_5_end: int | None = None
    flank_3_start: int | None = None
    flank_3_end: int | None = None

    # CRF state sequence
    state_sequence: list[int] = field(default_factory=list)
    state_posteriors: np.ndarray = field(default=None)

    # Evidence
    hallmark_genes: list[str] = field(default_factory=list)
    max_kfd: float = 0.0
    gc_deviation: float = 0.0
    cub_deviation: float = 0.0
    mean_novelty: float = 0.0

    # Window features for Phase 3 coherence analysis
    # Stored as list of WindowFeatures from Tier 2 (or Tier 1 if skip_tier2)
    window_features: list = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start

    def to_bed_line(self) -> str:
        """Format as BED line."""
        score = int(min(1000, self.confidence * 1000))
        return f"{self.scaffold}\t{self.start}\t{self.end}\t{boundary_candidate_id(self)}\t{score}\t."

    def to_gff_line(self) -> str:
        """Format as GFF3 line with attributes."""
        score = int(min(1000, self.confidence * 1000))
        attrs = [
            f"ID={boundary_candidate_id(self)}",
            f"confidence={self.confidence:.3f}",
            f"posterior={self.posterior_probability:.3f}",
        ]
        if self.core_viral_start is not None:
            # Same 0-based half-open to 1-based inclusive shift as the
            # coordinate columns below, so the attributes match their record.
            attrs.append(f"core_start={self.core_viral_start + 1}")
            attrs.append(f"core_end={self.core_viral_end}")
        if self.hallmark_genes:
            attrs.append(f"hallmarks={','.join(self.hallmark_genes)}")

        return f"{self.scaffold}\tViroSync\tEVE\t{self.start + 1}\t{self.end}\t{score}\t.\t.\t{';'.join(attrs)}"


def boundary_candidate_id(boundary: RefinedBoundary) -> str:
    """Return the public candidate ID, including the transient singleton default."""
    candidate_id = getattr(boundary, "candidate_id", "")
    if candidate_id:
        return candidate_id
    return f"EVE_{boundary.scaffold}_{boundary.start}-{boundary.end}"


def assign_boundary_candidate_ids(
    boundaries: Sequence[RefinedBoundary],
    *,
    reserved_candidate_ids: Iterable[str] = (),
) -> list[RefinedBoundary]:
    """Assign deterministic IDs without renaming already assigned candidates."""
    boundary_list = list(boundaries)
    existing_ids = [boundary.candidate_id for boundary in boundary_list]
    if any(existing_ids):
        if not all(existing_ids):
            raise ValueError("candidate IDs must be either all assigned or all empty")
        duplicates = sorted(candidate_id for candidate_id, count in Counter(existing_ids).items() if count > 1)
        if duplicates:
            raise ValueError(f"duplicate assigned candidate IDs: {duplicates}")
        reserved_conflicts = sorted(set(existing_ids) & set(reserved_candidate_ids))
        if reserved_conflicts:
            raise ValueError(f"assigned candidate IDs are reserved: {reserved_conflicts}")
        return boundary_list

    reserved_ids = set(reserved_candidate_ids)
    base_ids = [boundary_candidate_id(boundary) for boundary in boundary_list]
    base_counts = Counter(base_ids)
    used_ids = set(reserved_ids)
    assigned: list[RefinedBoundary] = []
    for boundary, base_id in zip(boundary_list, base_ids, strict=True):
        candidate_id = base_id
        if base_counts[base_id] > 1 or base_id in reserved_ids:
            digest = hashlib.sha256(boundary.seed_id.encode("utf-8")).hexdigest()[:16]
            candidate_id = f"{base_id}-c{digest}"
        if candidate_id in used_ids:
            raise ValueError(
                "cannot assign distinct candidate IDs to boundaries with the same "
                f"coordinates and seed_id: {base_id!r}, {boundary.seed_id!r}"
            )
        used_ids.add(candidate_id)
        assigned.append(replace(boundary, candidate_id=candidate_id))
    return assigned
