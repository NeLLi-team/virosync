from __future__ import annotations

import math
from dataclasses import dataclass

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
    fields = _complete_marker_fields()
    fields["hmm_target"] = "GVOGm0054"
    hit = hit_type(**fields)

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
        "reported_capsid",
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


@pytest.mark.parametrize("status", ["supported", "unvalidated"])
def test_a2_does_not_count_already_retained_capsid_as_bypass(status):
    fields = _complete_marker_fields()
    fields.update(hmm_target="GVOGm0003", validation_status=status)
    hit = ValidationMarkerHit(**fields)
    assert decide_marker_hit_role(hit, ablation_id=AblationID.A0).role is MarkerRole.REPORTED_CAPSID
    decision = decide_marker_hit_role(hit, ablation_id=AblationID.A2)
    assert decision.role is MarkerRole.REPORTED_CAPSID
    assert not decision.is_tier1_bypassed
    assert hit.validation_status == status


def test_override_does_not_erase_eligible_viral_support_on_same_gene():
    from dataclasses import replace

    from virosync.pipeline.phase1.viral_markers import viral_marker_gene_ids

    fields = _complete_marker_fields()
    fields.update(validation_status="validated", hmm_target="GVOGm0003")
    override = ValidationMarkerHit(**fields, pfam_override="reassigned")
    selected = replace(override, hmm_target="GVOGm0054", pfam_override="")
    assert viral_marker_gene_ids([override]) == frozenset()
    assert viral_marker_gene_ids([override, selected]) == frozenset({"scaffold_1_1"})


@pytest.mark.parametrize("status", ["validated", "validated_novel", "unvalidated"])
@pytest.mark.parametrize("score,retained", [(9.999, False), (10.0, True)])
def test_capsid_score_floor_precedes_validation(status, score, retained):
    from types import SimpleNamespace

    marker = SimpleNamespace(hmm_target="VP_Penton_1", hmm_score=score, hmm_evalue=1e-20, validation_status=status)
    assert decide_marker_hit_role(marker, ablation_id=AblationID.A0).is_retained_evidence is retained


def test_ineligible_capsid_cannot_hide_ordinary_hit_on_same_protein():
    from types import SimpleNamespace

    from virosync.pipeline.phase1.seed_merger import best_marker_hits

    weak = SimpleNamespace(
        query_porf="p|aa1-20", hmm_target="VP_Penton_1", hmm_score=9.9, validation_status="validated"
    )
    ordinary = SimpleNamespace(
        query_porf="p|aa1-20", hmm_target="GVOGm0760", hmm_score=9.0, validation_status="validated"
    )
    assert best_marker_hits([weak, ordinary]) == [ordinary]


@pytest.mark.parametrize("sequence_score,domain_score,eligible", [(9.9, 100, False), (10, -5, True)])
def test_capsid_floor_uses_full_sequence_score(sequence_score, domain_score, eligible):
    from virosync.pipeline.phase1.hhg_seeding import HMMHit
    from virosync.pipeline.phase1.marker_roles import passes_capsid_score_floor

    raw = HMMHit("p", "VP_Penton_1", sequence_score, 1e-30, domain_score, 1, 20)
    assert passes_capsid_score_floor(raw) is eligible


@pytest.mark.parametrize("score,eligible", [(9.999, False), (10, True)])
def test_direct_region_assembly_capsid_admission_obeys_floor(score, eligible):
    fields = _complete_marker_fields()
    fields.update(hmm_target="VP_Penton_1", hmm_score=score, validation_status="validated", tier1_bypassed=True)
    assert RegionMarkerHit(**fields).is_valid_seed_marker is eligible


@pytest.mark.parametrize(
    "score,eligible",
    [
        ("9.999", False),
        ("10", True),
        ("75.0", True),
        ("bad", False),
        ("nan", False),
        ("inf", False),
        ("", False),
        (None, False),
        (True, False),
    ],
)
def test_serialized_capsid_scores_preserve_eligibility(score, eligible):
    from virosync.pipeline.phase1.marker_roles import passes_capsid_score_floor

    assert passes_capsid_score_floor({"hmm_target": "GVOGm0003", "hmm_score": score}) is eligible
