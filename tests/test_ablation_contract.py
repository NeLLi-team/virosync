from dataclasses import FrozenInstanceError
import hashlib
import json

import pytest

from virosync.ablation import (
    ABLATION_CONTRACT_SCHEMA,
    ABLATION_CONTRACT_SHA256,
    ABLATION_EVENTS_SCHEMA,
    ABLATION_POLICIES,
    AblationContractError,
    AblationCounters,
    AblationEvents,
    AblationID,
    InterventionCounts,
    InterventionKey,
    ablation_contract_bytes,
    ablation_contract_document,
    ablation_policy,
    validate_ablation_events_bytes,
    validate_ablation_events_document,
)


EXPECTED_POLICIES = {
    AblationID.A0: (
        "Full production workflow",
        "Reference",
        None,
        (),
    ),
    AblationID.A1: (
        "Export Tier-1-vetted Phase-1 merged_seeds as predictions; do not run Phase 2 or Phase 3",
        "At least one A1 seed that A0 later rejects or changes",
        1,
        (InterventionKey.PHASE1_SEED_SURFACE_EXPORT,),
    ),
    AblationID.A2: (
        "Admit strong HMM-qualified anchors rejected by Tier 1 into region assembly and retain them as explicitly bypassed marker evidence downstream",
        "At least one bypassed anchor creates or changes a seed",
        1,
        (InterventionKey.TIER1_TAXONOMY_GATE_BYPASS,),
    ),
    AblationID.A3: (
        "Forward exact Phase-1 merged_seeds into Phase 3; skip Phase-2 extension, Tier-2 search, and boundary refinement",
        "A0 changes a boundary that A3 leaves at the Phase-1 coordinates",
        2,
        (InterventionKey.PHASE2_EXPANSION_REFINEMENT_BYPASS,),
    ),
    AblationID.A4: (
        "Disable all coordinate-changing host paths in Phase 2a, taxonomy-refinement barriers, and Phase 2f; retain Phase-3 host-confidence penalties",
        "A0 trims or stops one boundary that A4 retains",
        2,
        (InterventionKey.HOST_COORDINATE_CHANGE_BYPASS,),
    ),
    AblationID.A5: (
        "Compute GC/KFD fields for audit, but remove their weight, bonuses, and caps from confidence calculation; do not claim a CUB ablation because CUB is inactive",
        "One nonzero-composition candidate changes score or tier",
        3,
        (InterventionKey.COMPOSITION_EVIDENCE_BYPASS,),
    ),
    AblationID.A6: (
        "Accept every scored Phase-3 candidate; compute the normal class/length/marker gate decisions only as counterfactual events",
        "At least one candidate A0 would reject is retained",
        3,
        (InterventionKey.FINAL_ACCEPTANCE_GATE_BYPASS,),
    ),
}


def _event(ablation_id: AblationID = AblationID.A3) -> AblationEvents:
    return AblationEvents(
        ablation_id=ablation_id,
        counters=AblationCounters.for_ablation(
            ablation_id,
            opportunities=8 if ablation_id is not AblationID.A0 else 0,
            interventions=3 if ablation_id is not AblationID.A0 else 0,
            changed=1 if ablation_id is not AblationID.A0 else 0,
        ),
    )


def test_closed_ids_and_frozen_policy_registry_match_decision() -> None:
    assert [item.value for item in AblationID] == [
        "A0",
        "A1",
        "A2",
        "A3",
        "A4",
        "A5",
        "A6",
    ]
    assert tuple(ABLATION_POLICIES) == tuple(AblationID)
    assert {
        ablation_id: (
            policy.frozen_intervention,
            policy.required_counterfactual,
            policy.intervention_phase,
            policy.active_intervention_keys,
        )
        for ablation_id, policy in ABLATION_POLICIES.items()
    } == EXPECTED_POLICIES

    with pytest.raises(TypeError):
        ABLATION_POLICIES[AblationID.A0] = ABLATION_POLICIES[AblationID.A1]
    with pytest.raises(FrozenInstanceError):
        ABLATION_POLICIES[AblationID.A1].frozen_intervention = "changed"


def test_a0_has_no_active_intervention_and_each_other_key_is_unique() -> None:
    assert ablation_policy(AblationID.A0).active_intervention_keys == ()
    assert tuple(
        key for ablation_id in tuple(AblationID)[1:] for key in ablation_policy(ablation_id).active_intervention_keys
    ) == tuple(InterventionKey)
    with pytest.raises(TypeError, match="AblationID"):
        ablation_policy("A0")  # type: ignore[arg-type]


def test_contract_document_is_fresh_and_canonical_digest_is_golden() -> None:
    first = ablation_contract_document()
    first["schema"] = "mutated"
    assert ablation_contract_document()["schema"] == ABLATION_CONTRACT_SCHEMA

    content = ablation_contract_bytes()
    assert content == json.dumps(
        ablation_contract_document(),
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    assert hashlib.sha256(content).hexdigest() == ABLATION_CONTRACT_SHA256
    assert ABLATION_CONTRACT_SHA256 == ("bcde3b3053760d897e1756454038e0ba8aa5bf75148e0c6eabc56cb080653bd2")


@pytest.mark.parametrize("ablation_id", list(AblationID))
def test_selected_counter_group_round_trips_and_others_are_zero(
    ablation_id: AblationID,
) -> None:
    events = _event(ablation_id)
    active = frozenset(ablation_policy(ablation_id).active_intervention_keys)
    for key in InterventionKey:
        counts = events.counters.for_key(key)
        if key in active:
            assert counts == InterventionCounts(8, 3, 1)
        else:
            assert counts.is_zero

    assert events.counters.total_opportunities == (8 if active else 0)
    assert events.counters.total_interventions == (3 if active else 0)
    assert events.counters.total_changed == (1 if active else 0)
    assert AblationEvents.from_bytes(events.to_bytes()) == events


def test_a0_rejects_counter_activity() -> None:
    with pytest.raises(AblationContractError, match="A0 cannot"):
        AblationCounters.for_ablation(
            AblationID.A0,
            opportunities=1,
            interventions=1,
            changed=1,
        )
    with pytest.raises(AblationContractError, match="non-selected"):
        AblationEvents(
            ablation_id=AblationID.A0,
            counters=AblationCounters(phase1_seed_surface_export=InterventionCounts(1, 1, 1)),
        )


def test_non_selected_counter_groups_must_remain_zero() -> None:
    with pytest.raises(AblationContractError, match="non-selected"):
        AblationEvents(
            ablation_id=AblationID.A2,
            counters=AblationCounters(
                tier1_taxonomy_gate_bypass=InterventionCounts(3, 1, 1),
                final_acceptance_gate_bypass=InterventionCounts(2, 1, 0),
            ),
        )


def test_counter_groups_require_typed_immutable_counts() -> None:
    with pytest.raises(TypeError, match="must be InterventionCounts"):
        AblationCounters(
            tier1_taxonomy_gate_bypass={  # type: ignore[arg-type]
                "opportunities": 1,
                "interventions": 1,
                "changed": 1,
            }
        )


@pytest.mark.parametrize(
    ("counts", "message"),
    [
        ((-1, 0, 0), "nonnegative"),
        ((1, -1, 0), "nonnegative"),
        ((1, 1, -1), "nonnegative"),
        ((1, 2, 0), "cannot exceed opportunities"),
        ((2, 1, 2), "cannot exceed interventions"),
        ((True, 0, 0), "must be an integer"),
        ((1, False, 0), "must be an integer"),
        ((1, 1, True), "must be an integer"),
    ],
)
def test_counter_values_are_strict_and_monotonic(
    counts: tuple[object, object, object],
    message: str,
) -> None:
    with pytest.raises(AblationContractError, match=message):
        InterventionCounts(*counts)  # type: ignore[arg-type]


def test_zero_selected_counter_group_is_valid_for_a_biological_run() -> None:
    events = AblationEvents(
        ablation_id=AblationID.A6,
        counters=AblationCounters.for_ablation(AblationID.A6),
    )
    assert events.counters.total_interventions == 0
    assert AblationEvents.from_bytes(events.to_bytes()) == events


def test_event_document_has_exact_schema_and_all_counter_groups() -> None:
    document = _event().to_document()
    assert set(document) == {
        "schema",
        "contract_sha256",
        "ablation_id",
        "counters",
    }
    assert document["schema"] == ABLATION_EVENTS_SCHEMA
    assert document["contract_sha256"] == ABLATION_CONTRACT_SHA256
    assert document["ablation_id"] == "A3"
    assert set(document["counters"]) == {key.value for key in InterventionKey}
    assert validate_ablation_events_document(document) == _event()


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda document: document.update(extra=True), "missing or unknown"),
        (lambda document: document.pop("counters"), "missing or unknown"),
        (
            lambda document: document.__setitem__("schema", "wrong/v1"),
            "unsupported",
        ),
        (
            lambda document: document.__setitem__("contract_sha256", "0" * 64),
            "SHA256 mismatch",
        ),
        (lambda document: document.__setitem__("ablation_id", "A7"), "unknown"),
        (
            lambda document: document.__setitem__("ablation_id", 3),
            "must be a string",
        ),
        (
            lambda document: document["counters"].pop("phase1_seed_surface_export"),
            "registered intervention keys",
        ),
        (
            lambda document: document["counters"].update(unknown={}),
            "registered intervention keys",
        ),
        (
            lambda document: document["counters"]["phase2_expansion_refinement_bypass"].update(extra=0),
            "must contain exactly",
        ),
    ],
)
def test_document_validation_rejects_schema_drift(mutation, message: str) -> None:
    document = _event().to_document()
    mutation(document)
    with pytest.raises(AblationContractError, match=message):
        validate_ablation_events_document(document)


def test_canonical_bytes_are_deterministic_and_strict() -> None:
    events = _event(AblationID.A5)
    content = events.to_bytes()
    assert content == events.to_bytes()
    assert validate_ablation_events_bytes(content) == events

    decoded = json.loads(content)
    noncanonical = json.dumps(decoded, indent=2, sort_keys=False).encode("utf-8")
    with pytest.raises(AblationContractError, match="canonical form"):
        validate_ablation_events_bytes(noncanonical)
    with pytest.raises(AblationContractError, match="canonical form"):
        validate_ablation_events_bytes(content + b"\n")
    with pytest.raises(TypeError, match="must be bytes"):
        validate_ablation_events_bytes(bytearray(content))  # type: ignore[arg-type]


def test_bytes_validation_rejects_invalid_json_and_large_input() -> None:
    with pytest.raises(AblationContractError, match="UTF-8 JSON"):
        validate_ablation_events_bytes(b"{not json}")
    with pytest.raises(AblationContractError, match="UTF-8 JSON"):
        validate_ablation_events_bytes(b"\xff")
    with pytest.raises(AblationContractError, match="too large"):
        validate_ablation_events_bytes(b" " * (64 * 1024 + 1))
    deeply_nested = b"[" * 2_000 + b"0" + b"]" * 2_000
    with pytest.raises(AblationContractError, match="valid UTF-8 JSON"):
        validate_ablation_events_bytes(deeply_nested)


def test_event_types_are_immutable() -> None:
    events = _event()
    with pytest.raises(FrozenInstanceError):
        events.ablation_id = AblationID.A0
    with pytest.raises(FrozenInstanceError):
        events.counters.phase1_seed_surface_export = InterventionCounts()
