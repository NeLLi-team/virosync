from __future__ import annotations

from dataclasses import dataclass
import math

import pytest

from virosync.ablation import AblationID
from virosync.orchestration._flows.single_genome.phase1 import (
    _seed_annotation_markers,
)
from virosync.pipeline.phase1.marker_roles import (
    A2_MAX_HMM_EVALUE,
    DEFAULT_A2_SINGLE_MARKER_MIN_SCORE,
    MarkerRole,
    decide_marker_hit_role,
    decide_marker_role,
)
from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit as ValidationMarkerHit
from virosync.pipeline.phase1.region_assembly import ValidatedMarkerHit as RegionMarkerHit
from virosync.pipeline.phase1.viral_markers import get_assembly_mode


@dataclass
class _Hit:
    validation_status: str
    hmm_score: float
    hmm_evalue: float


def _complete_marker_fields() -> dict[str, object]:
    return {
        "query_porf": "scaffold_1_1",
        "scaffold": "scaffold_1",
        "start": 10,
        "end": 200,
        "strand": "+",
        "hmm_target": "GVOGm0003",
        "hmm_score": 50.0,
        "hmm_evalue": 1e-8,
        "validation_status": "supported",
        "top10_prefixes": "GVMAG__",
        "best_hit_target": "GVMAG__example",
        "best_hit_pident": 20.0,
        "best_hit_bits": 45.0,
        "has_ncldv": 0,
        "has_mirus": 0,
        "has_plv": 0,
        "has_vp": 0,
        "has_viral": 1,
    }


@pytest.mark.parametrize("status", ["validated", "validated_novel"])
@pytest.mark.parametrize("ablation_id", list(AblationID))
def test_production_statuses_keep_production_role(
    status: str,
    ablation_id: AblationID,
) -> None:
    decision = decide_marker_role(
        ablation_id=ablation_id,
        validation_status=status,
        hmm_score=float("nan"),
        hmm_evalue=float("inf"),
    )

    assert decision.role is MarkerRole.PRODUCTION_VALIDATED
    assert decision.original_validation_status == status
    assert decision.is_production_validated
    assert decision.is_retained_evidence
    assert not decision.is_tier1_bypassed
    assert not decision.is_rejected


@pytest.mark.parametrize("status", ["supported", "unvalidated"])
def test_a2_bypasses_eligible_tier1_rejections_at_inclusive_boundaries(status: str) -> None:
    decision = decide_marker_role(
        ablation_id=AblationID.A2,
        validation_status=status,
        hmm_score=DEFAULT_A2_SINGLE_MARKER_MIN_SCORE,
        hmm_evalue=A2_MAX_HMM_EVALUE,
    )

    assert decision.role is MarkerRole.TIER1_BYPASSED
    assert decision.original_validation_status == status
    assert decision.is_tier1_bypassed
    assert decision.is_retained_evidence


@pytest.mark.parametrize("ablation_id", [arm for arm in AblationID if arm is not AblationID.A2])
def test_strong_rejected_marker_is_bypassed_only_by_a2(ablation_id: AblationID) -> None:
    decision = decide_marker_role(
        ablation_id=ablation_id,
        validation_status="supported",
        hmm_score=500.0,
        hmm_evalue=1e-100,
    )

    assert decision.role is MarkerRole.REJECTED
    assert decision.is_rejected
    assert not decision.is_retained_evidence


@pytest.mark.parametrize(
    ("status", "score", "evalue"),
    [
        ("supported", DEFAULT_A2_SINGLE_MARKER_MIN_SCORE - 0.001, A2_MAX_HMM_EVALUE),
        ("supported", DEFAULT_A2_SINGLE_MARKER_MIN_SCORE, A2_MAX_HMM_EVALUE * 1.001),
        ("unvalidated", float("nan"), 1e-20),
        ("unvalidated", float("inf"), 1e-20),
        ("unvalidated", 500.0, float("nan")),
        ("unvalidated", 500.0, float("inf")),
        ("unvalidated", 500.0, -1e-20),
        ("unknown", 500.0, 1e-20),
    ],
)
def test_a2_rejects_ineligible_status_or_hmm_evidence(
    status: str,
    score: float,
    evalue: float,
) -> None:
    decision = decide_marker_role(
        ablation_id=AblationID.A2,
        validation_status=status,
        hmm_score=score,
        hmm_evalue=evalue,
    )

    assert decision.role is MarkerRole.REJECTED
    assert decision.original_validation_status == status


@pytest.mark.parametrize(
    ("mode_name", "floor"),
    [
        ("default", 100.0),
        ("fragmented", 70.0),
        ("relaxed", 50.0),
        ("strict", 150.0),
    ],
)
def test_actual_assembly_mode_floor_can_be_applied(mode_name: str, floor: float) -> None:
    assembly_mode = get_assembly_mode(mode_name)
    assert assembly_mode.single_marker_min_score == floor

    below = decide_marker_role(
        ablation_id=AblationID.A2,
        validation_status="supported",
        hmm_score=math.nextafter(floor, -math.inf),
        hmm_evalue=1e-20,
        single_marker_min_score=assembly_mode.single_marker_min_score,
    )
    at_floor = decide_marker_role(
        ablation_id=AblationID.A2,
        validation_status="supported",
        hmm_score=floor,
        hmm_evalue=1e-20,
        single_marker_min_score=assembly_mode.single_marker_min_score,
    )

    assert below.role is MarkerRole.REJECTED
    assert at_floor.role is MarkerRole.TIER1_BYPASSED


def test_hit_helper_preserves_original_object_and_status() -> None:
    hit = _Hit(
        validation_status="unvalidated",
        hmm_score=50.0,
        hmm_evalue=1e-8,
    )

    decision = decide_marker_hit_role(hit, ablation_id=AblationID.A2)

    assert decision.role is MarkerRole.TIER1_BYPASSED
    assert decision.original_validation_status == "unvalidated"
    assert hit == _Hit(
        validation_status="unvalidated",
        hmm_score=50.0,
        hmm_evalue=1e-8,
    )


@pytest.mark.parametrize("hit_type", [ValidationMarkerHit, RegionMarkerHit])
def test_hit_helper_accepts_both_existing_validated_marker_shapes(hit_type: type[object]) -> None:
    hit = hit_type(**_complete_marker_fields())

    decision = decide_marker_hit_role(hit, ablation_id=AblationID.A2)  # type: ignore[arg-type]

    assert decision.role is MarkerRole.TIER1_BYPASSED
    assert decision.original_validation_status == "supported"
    assert hit.validation_status == "supported"  # type: ignore[attr-defined]


@pytest.mark.parametrize("floor", [True, "50", -1.0, float("nan"), float("inf")])
def test_invalid_assembly_floor_is_rejected(floor: object) -> None:
    error = TypeError if isinstance(floor, (bool, str)) else ValueError
    with pytest.raises(error):
        decide_marker_role(
            ablation_id=AblationID.A2,
            validation_status="supported",
            hmm_score=50.0,
            hmm_evalue=1e-8,
            single_marker_min_score=floor,  # type: ignore[arg-type]
        )


def test_closed_role_values_are_stable() -> None:
    assert [role.value for role in MarkerRole] == [
        "production_validated",
        "tier1_bypassed",
        "rejected",
    ]


def test_raw_ablation_string_is_not_accepted() -> None:
    with pytest.raises(TypeError, match="AblationID"):
        decide_marker_role(
            ablation_id="A2",  # type: ignore[arg-type]
            validation_status="supported",
            hmm_score=50.0,
            hmm_evalue=1e-8,
        )


def test_seed_annotation_preserves_rejected_historical_marker_surface() -> None:
    production_marker = _Hit("validated", 50.0, 1e-8)
    rejected_marker = _Hit("unvalidated", 50.0, 1e-8)
    deviation_marker = _Hit("supported", 60.0, 1e-9)

    markers = _seed_annotation_markers(
        [production_marker, rejected_marker],
        [rejected_marker, deviation_marker],
    )

    assert markers == [production_marker, rejected_marker, deviation_marker]
