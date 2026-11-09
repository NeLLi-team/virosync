"""Pure marker-role decisions for the Tier-1 taxonomy-gate ablation.

The role is deliberately separate from ``validation_status``.  In particular,
an A2 bypass does not rewrite a Tier-1 ``supported`` or ``unvalidated`` result;
it records why that original evidence is retained downstream.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from typing import Final, Protocol

from virosync.ablation import AblationID

MIN_CAPSID_BITSCORE: Final = 10.0
DEFAULT_A2_SINGLE_MARKER_MIN_SCORE: Final = 50.0
A2_MAX_HMM_EVALUE: Final = 1e-5

PRODUCTION_VALIDATION_STATUSES: Final = frozenset({"validated", "validated_novel"})
_A2_BYPASS_VALIDATION_STATUSES: Final = frozenset({"supported", "unvalidated"})


class MarkerRole(str, Enum):
    """Closed downstream role assigned to one Tier-1 marker result."""

    PRODUCTION_VALIDATED = "production_validated"
    REPORTED_CAPSID = "reported_capsid"
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

    if validation_status in PRODUCTION_VALIDATION_STATUSES:
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


def is_arbitration_selected_marker(hit: object) -> bool:
    """Check Pfam selection before validation has assigned a downstream role."""
    override = hit.get("pfam_override", "") if isinstance(hit, dict) else getattr(hit, "pfam_override", "")
    return not override


def passes_capsid_score_floor(hit: object) -> bool:
    """Require ten bits only for semantic capsid/Penton admission and evidence."""
    from virosync.pipeline.phase3.mcp_detection import is_capsid_evidence_gene

    fields = (
        hit
        if isinstance(hit, dict)
        else {
            key: getattr(hit, key)
            for key in ("hmm_target", "hallmark_gene", "model", "target_name", "hmm_score", "score")
            if hasattr(hit, key)
        }
    )
    model = next(
        (fields.get(key) for key in ("hmm_target", "hallmark_gene", "model", "target_name") if fields.get(key)), ""
    )
    if not is_capsid_evidence_gene(model):
        return True
    score = fields.get("hmm_score", fields.get("score", 0.0))
    if isinstance(score, bool) or not isinstance(score, (int, float, str)):
        return False
    try:
        score = float(score)
    except ValueError:
        return False
    return math.isfinite(score) and score >= MIN_CAPSID_BITSCORE


def is_ordinary_marker_evidence(hit: object) -> bool:
    """Allow production validation or an explicit noncapsid A2 bypass for ordinary credit."""
    from virosync.pipeline.phase3.mcp_detection import is_capsid_evidence_gene

    if not passes_capsid_score_floor(hit) or not is_arbitration_selected_marker(hit):
        return False
    if isinstance(hit, dict):
        status = hit.get("validation_status", "")
        bypassed = hit.get("tier1_bypassed", False)
        model = hit.get("hmm_target") or hit.get("hallmark_gene") or hit.get("model") or ""
    else:
        status = getattr(hit, "validation_status", "")
        bypassed = getattr(hit, "tier1_bypassed", False)
        model = getattr(hit, "hmm_target", "") or getattr(hit, "hallmark_gene", "") or getattr(hit, "model", "")
    return status in PRODUCTION_VALIDATION_STATUSES or (
        bypassed is True and bool(model) and not is_capsid_evidence_gene(model)
    )


def decide_marker_hit_role(
    hit: MarkerRoleEvidence,
    *,
    ablation_id: AblationID,
    single_marker_min_score: float = DEFAULT_A2_SINGLE_MARKER_MIN_SCORE,
) -> MarkerRoleDecision:
    """Assign a role from either existing ``ValidatedMarkerHit`` shape."""
    from virosync.pipeline.phase3.mcp_detection import is_capsid_evidence_gene

    if not passes_capsid_score_floor(hit):
        return MarkerRoleDecision(MarkerRole.REJECTED, hit.validation_status)

    decision = decide_marker_role(
        ablation_id=ablation_id,
        validation_status=hit.validation_status,
        hmm_score=hit.hmm_score,
        hmm_evalue=hit.hmm_evalue,
        single_marker_min_score=single_marker_min_score,
    )
    if is_capsid_evidence_gene(getattr(hit, "hmm_target", "")) and (
        decision.is_rejected or decision.is_tier1_bypassed or not is_arbitration_selected_marker(hit)
    ):
        return MarkerRoleDecision(MarkerRole.REPORTED_CAPSID, hit.validation_status)
    return decision
