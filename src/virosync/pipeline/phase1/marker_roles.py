"""Pure marker-role decisions for the Tier-1 taxonomy-gate ablation.

The role is deliberately separate from ``validation_status``.  In particular,
an A2 bypass does not rewrite a Tier-1 ``supported`` or ``unvalidated`` result;
it records why that original evidence is retained downstream.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math
from typing import Final, Protocol

from virosync.ablation import AblationID


DEFAULT_A2_SINGLE_MARKER_MIN_SCORE: Final = 50.0
A2_MAX_HMM_EVALUE: Final = 1e-5

_PRODUCTION_VALIDATION_STATUSES: Final = frozenset({"validated", "validated_novel"})
_A2_BYPASS_VALIDATION_STATUSES: Final = frozenset({"supported", "unvalidated"})


class MarkerRole(str, Enum):
    """Closed downstream role assigned to one Tier-1 marker result."""

    PRODUCTION_VALIDATED = "production_validated"
    TIER1_BYPASSED = "tier1_bypassed"
    REJECTED = "rejected"


class MarkerRoleEvidence(Protocol):
    """Structural fields required from either validated-marker hit class."""

    validation_status: str
    hmm_score: float
    hmm_evalue: float


@dataclass(frozen=True, slots=True)
class MarkerRoleDecision:
    """Role plus the unchanged Tier-1 status on which it was based."""

    role: MarkerRole
    original_validation_status: str

    @property
    def is_production_validated(self) -> bool:
        """Return whether normal Tier-1 validation assigned this role."""

        return self.role is MarkerRole.PRODUCTION_VALIDATED

    @property
    def is_tier1_bypassed(self) -> bool:
        """Return whether A2 retained a marker rejected by normal Tier 1."""

        return self.role is MarkerRole.TIER1_BYPASSED

    @property
    def is_rejected(self) -> bool:
        """Return whether the marker cannot seed a region or become retained evidence.

        Historical A0 behavior still permits rejected overlapping HMM hits to
        annotate a seed created by production-validated markers.
        """

        return self.role is MarkerRole.REJECTED

    @property
    def is_retained_evidence(self) -> bool:
        """Return whether downstream stages may retain this marker as evidence.

        This is not a complete region-seed decision.  Existing production
        restrictions, such as the MCP requirement for ``validated_novel``,
        still apply independently.
        """

        return not self.is_rejected


def _validate_single_marker_min_score(value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError("single_marker_min_score must be a number")
    floor = float(value)
    if not math.isfinite(floor) or floor < 0.0:
        raise ValueError("single_marker_min_score must be finite and nonnegative")
    return floor


def _has_strong_hmm_evidence(
    *,
    hmm_score: float,
    hmm_evalue: float,
    single_marker_min_score: float,
) -> bool:
    if isinstance(hmm_score, bool) or not isinstance(hmm_score, (int, float)):
        return False
    if isinstance(hmm_evalue, bool) or not isinstance(hmm_evalue, (int, float)):
        return False
    score = float(hmm_score)
    evalue = float(hmm_evalue)
    return (
        math.isfinite(score)
        and math.isfinite(evalue)
        and evalue >= 0.0
        and score >= single_marker_min_score
        and evalue <= A2_MAX_HMM_EVALUE
    )


def decide_marker_role(
    *,
    ablation_id: AblationID,
    validation_status: str,
    hmm_score: float,
    hmm_evalue: float,
    single_marker_min_score: float = DEFAULT_A2_SINGLE_MARKER_MIN_SCORE,
) -> MarkerRoleDecision:
    """Assign a deterministic downstream role without altering Tier-1 status.

    ``single_marker_min_score`` should be the selected assembly mode's floor.
    The fallback of 50 is the lowest existing production assembly-mode floor.
    Only A2 can retain ``supported`` or ``unvalidated`` markers, and then only
    when both the score floor and the fixed ``1e-5`` E-value ceiling pass.
    """

    if not isinstance(ablation_id, AblationID):
        raise TypeError("ablation_id must be an AblationID")
    if not isinstance(validation_status, str):
        raise TypeError("validation_status must be a string")

    floor = _validate_single_marker_min_score(single_marker_min_score)

    if validation_status in _PRODUCTION_VALIDATION_STATUSES:
        role = MarkerRole.PRODUCTION_VALIDATED
    elif (
        ablation_id is AblationID.A2
        and validation_status in _A2_BYPASS_VALIDATION_STATUSES
        and _has_strong_hmm_evidence(
            hmm_score=hmm_score,
            hmm_evalue=hmm_evalue,
            single_marker_min_score=floor,
        )
    ):
        role = MarkerRole.TIER1_BYPASSED
    else:
        role = MarkerRole.REJECTED

    return MarkerRoleDecision(
        role=role,
        original_validation_status=validation_status,
    )


def decide_marker_hit_role(
    hit: MarkerRoleEvidence,
    *,
    ablation_id: AblationID,
    single_marker_min_score: float = DEFAULT_A2_SINGLE_MARKER_MIN_SCORE,
) -> MarkerRoleDecision:
    """Assign a role from either existing ``ValidatedMarkerHit`` shape."""

    return decide_marker_role(
        ablation_id=ablation_id,
        validation_status=hit.validation_status,
        hmm_score=hit.hmm_score,
        hmm_evalue=hit.hmm_evalue,
        single_marker_min_score=single_marker_min_score,
    )
