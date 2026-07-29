"""Closed, deterministic contract for the ViroSync benchmark ablations.

This module is intentionally leaf-like: it defines policy and serialization,
but performs no I/O and imports no configuration or orchestration code.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
import hashlib
import json
from types import MappingProxyType
from typing import Final


ABLATION_CONTRACT_SCHEMA: Final = "virosync.ablation_contract/v1"
ABLATION_EVENTS_SCHEMA: Final = "virosync.ablation_events/v1"
MAX_ABLATION_EVENTS_BYTES: Final = 64 * 1024


class AblationContractError(ValueError):
    """Raised when an ablation policy, counter, or document is invalid."""


class AblationID(str, Enum):
    """The complete set of mutually exclusive benchmark ablations."""

    A0 = "A0"
    A1 = "A1"
    A2 = "A2"
    A3 = "A3"
    A4 = "A4"
    A5 = "A5"
    A6 = "A6"


class InterventionKey(str, Enum):
    """Stable names for the six intervention-counter groups."""

    PHASE1_SEED_SURFACE_EXPORT = "phase1_seed_surface_export"
    TIER1_TAXONOMY_GATE_BYPASS = "tier1_taxonomy_gate_bypass"
    PHASE2_EXPANSION_REFINEMENT_BYPASS = "phase2_expansion_refinement_bypass"
    HOST_COORDINATE_CHANGE_BYPASS = "host_coordinate_change_bypass"
    COMPOSITION_EVIDENCE_BYPASS = "composition_evidence_bypass"
    FINAL_ACCEPTANCE_GATE_BYPASS = "final_acceptance_gate_bypass"


@dataclass(frozen=True, slots=True)
class AblationPolicy:
    """Immutable scientific definition for one ablation arm."""

    ablation_id: AblationID
    frozen_intervention: str
    required_counterfactual: str
    intervention_phase: int | None
    active_intervention_keys: tuple[InterventionKey, ...]


ABLATION_POLICIES: Final[Mapping[AblationID, AblationPolicy]] = MappingProxyType(
    {
        AblationID.A0: AblationPolicy(
            ablation_id=AblationID.A0,
            frozen_intervention="Full production workflow",
            required_counterfactual="Reference",
            intervention_phase=None,
            active_intervention_keys=(),
        ),
        AblationID.A1: AblationPolicy(
            ablation_id=AblationID.A1,
            frozen_intervention=(
                "Export Tier-1-vetted Phase-1 merged_seeds as predictions; do not run Phase 2 or Phase 3"
            ),
            required_counterfactual=("At least one A1 seed that A0 later rejects or changes"),
            intervention_phase=1,
            active_intervention_keys=(InterventionKey.PHASE1_SEED_SURFACE_EXPORT,),
        ),
        AblationID.A2: AblationPolicy(
            ablation_id=AblationID.A2,
            frozen_intervention=(
                "Admit strong HMM-qualified anchors rejected by Tier 1 into region "
                "assembly and retain them as explicitly bypassed marker evidence "
                "downstream"
            ),
            required_counterfactual=("At least one bypassed anchor creates or changes a seed"),
            intervention_phase=1,
            active_intervention_keys=(InterventionKey.TIER1_TAXONOMY_GATE_BYPASS,),
        ),
        AblationID.A3: AblationPolicy(
            ablation_id=AblationID.A3,
            frozen_intervention=(
                "Forward exact Phase-1 merged_seeds into Phase 3; skip Phase-2 "
                "extension, Tier-2 search, and boundary refinement"
            ),
            required_counterfactual=("A0 changes a boundary that A3 leaves at the Phase-1 coordinates"),
            intervention_phase=2,
            active_intervention_keys=(InterventionKey.PHASE2_EXPANSION_REFINEMENT_BYPASS,),
        ),
        AblationID.A4: AblationPolicy(
            ablation_id=AblationID.A4,
            frozen_intervention=(
                "Disable all coordinate-changing host paths in Phase 2a, "
                "taxonomy-refinement barriers, and Phase 2f; retain Phase-3 "
                "host-confidence penalties"
            ),
            required_counterfactual=("A0 trims or stops one boundary that A4 retains"),
            intervention_phase=2,
            active_intervention_keys=(InterventionKey.HOST_COORDINATE_CHANGE_BYPASS,),
        ),
        AblationID.A5: AblationPolicy(
            ablation_id=AblationID.A5,
            frozen_intervention=(
                "Compute GC/KFD fields for audit, but remove their weight, "
                "bonuses, and caps from confidence calculation; do not claim a "
                "CUB ablation because CUB is inactive"
            ),
            required_counterfactual=("One nonzero-composition candidate changes score or tier"),
            intervention_phase=3,
            active_intervention_keys=(InterventionKey.COMPOSITION_EVIDENCE_BYPASS,),
        ),
        AblationID.A6: AblationPolicy(
            ablation_id=AblationID.A6,
            frozen_intervention=(
                "Accept every scored Phase-3 candidate; compute the normal "
                "class/length/marker gate decisions only as counterfactual events"
            ),
            required_counterfactual=("At least one candidate A0 would reject is retained"),
            intervention_phase=3,
            active_intervention_keys=(InterventionKey.FINAL_ACCEPTANCE_GATE_BYPASS,),
        ),
    }
)


def _validate_policy_registry() -> None:
    if tuple(ABLATION_POLICIES) != tuple(AblationID):
        raise RuntimeError("ablation policy registry must cover A0-A6 in order")
    assigned = tuple(key for policy in ABLATION_POLICIES.values() for key in policy.active_intervention_keys)
    if ABLATION_POLICIES[AblationID.A0].active_intervention_keys:
        raise RuntimeError("A0 must not activate an intervention")
    if ABLATION_POLICIES[AblationID.A0].intervention_phase is not None:
        raise RuntimeError("A0 must not own an intervention phase")
    if assigned != tuple(InterventionKey):
        raise RuntimeError("each non-reference ablation must own one unique counter group")
    for ablation_id, policy in ABLATION_POLICIES.items():
        if policy.ablation_id is not ablation_id:
            raise RuntimeError("ablation policy registry key does not match policy")
        if ablation_id is not AblationID.A0 and policy.intervention_phase not in {1, 2, 3}:
            raise RuntimeError("each non-reference ablation must own Phase 1, 2, or 3")


_validate_policy_registry()


def ablation_policy(ablation_id: AblationID) -> AblationPolicy:
    """Return the frozen policy for an already validated ablation ID."""

    if not isinstance(ablation_id, AblationID):
        raise TypeError("ablation_id must be an AblationID")
    return ABLATION_POLICIES[ablation_id]


def _canonical_json_bytes(document: object) -> bytes:
    return json.dumps(
        document,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def ablation_contract_document() -> dict[str, object]:
    """Return a new JSON-compatible copy of the complete frozen contract."""

    return {
        "schema": ABLATION_CONTRACT_SCHEMA,
        "counter_fields": ["opportunities", "interventions", "changed"],
        "counter_invariants": {
            "a0_all_zero": True,
            "changed_lte_interventions_lte_opportunities": True,
            "non_selected_groups_zero": True,
        },
        "intervention_keys": [key.value for key in InterventionKey],
        "policies": [
            {
                "ablation_id": policy.ablation_id.value,
                "frozen_intervention": policy.frozen_intervention,
                "required_counterfactual": policy.required_counterfactual,
                "intervention_phase": policy.intervention_phase,
                "active_intervention_keys": [key.value for key in policy.active_intervention_keys],
            }
            for policy in ABLATION_POLICIES.values()
        ],
    }


def ablation_contract_bytes() -> bytes:
    """Return the sole canonical byte representation of the policy contract."""

    return _canonical_json_bytes(ablation_contract_document())


ABLATION_CONTRACT_SHA256: Final = hashlib.sha256(ablation_contract_bytes()).hexdigest()


def _strict_nonnegative_int(value: object, field: str) -> int:
    if type(value) is not int:
        raise AblationContractError(f"{field} must be an integer")
    if value < 0:
        raise AblationContractError(f"{field} must be nonnegative")
    return value


@dataclass(frozen=True, slots=True)
class InterventionCounts:
    """Opportunity, intervention, and changed-outcome counts for one hook."""

    opportunities: int = 0
    interventions: int = 0
    changed: int = 0

    def __post_init__(self) -> None:
        opportunities = _strict_nonnegative_int(
            self.opportunities,
            "opportunities",
        )
        interventions = _strict_nonnegative_int(
            self.interventions,
            "interventions",
        )
        changed = _strict_nonnegative_int(self.changed, "changed")
        if interventions > opportunities:
            raise AblationContractError("interventions cannot exceed opportunities")
        if changed > interventions:
            raise AblationContractError("changed outcomes cannot exceed interventions")

    @property
    def is_zero(self) -> bool:
        """Return whether the counter group records no activity."""

        return not (self.opportunities or self.interventions or self.changed)

    def to_document(self) -> dict[str, int]:
        """Return a JSON-compatible counter document."""

        return {
            "opportunities": self.opportunities,
            "interventions": self.interventions,
            "changed": self.changed,
        }

    @classmethod
    def from_document(
        cls,
        document: object,
        *,
        context: str,
    ) -> InterventionCounts:
        """Parse one counter group with an exact, closed field set."""

        if not isinstance(document, Mapping):
            raise AblationContractError(f"{context} must be an object")
        required = {"opportunities", "interventions", "changed"}
        if set(document) != required:
            raise AblationContractError(f"{context} must contain exactly {sorted(required)}")
        return cls(
            opportunities=_strict_nonnegative_int(
                document["opportunities"],
                f"{context}.opportunities",
            ),
            interventions=_strict_nonnegative_int(
                document["interventions"],
                f"{context}.interventions",
            ),
            changed=_strict_nonnegative_int(
                document["changed"],
                f"{context}.changed",
            ),
        )


_ZERO_COUNTS = InterventionCounts()


@dataclass(frozen=True, slots=True)
class AblationCounters:
    """Closed counter schema containing every named intervention group."""

    phase1_seed_surface_export: InterventionCounts = _ZERO_COUNTS
    tier1_taxonomy_gate_bypass: InterventionCounts = _ZERO_COUNTS
    phase2_expansion_refinement_bypass: InterventionCounts = _ZERO_COUNTS
    host_coordinate_change_bypass: InterventionCounts = _ZERO_COUNTS
    composition_evidence_bypass: InterventionCounts = _ZERO_COUNTS
    final_acceptance_gate_bypass: InterventionCounts = _ZERO_COUNTS

    def __post_init__(self) -> None:
        for key in InterventionKey:
            if type(self.for_key(key)) is not InterventionCounts:
                raise TypeError(f"{key.value} must be InterventionCounts")

    def for_key(self, key: InterventionKey) -> InterventionCounts:
        """Return one named counter group."""

        if not isinstance(key, InterventionKey):
            raise TypeError("key must be an InterventionKey")
        return getattr(self, key.value)

    def validate_for(self, ablation_id: AblationID) -> None:
        """Reject activity outside the counter group selected by an arm."""

        active = frozenset(ablation_policy(ablation_id).active_intervention_keys)
        invalid = [key.value for key in InterventionKey if key not in active and not self.for_key(key).is_zero]
        if invalid:
            raise AblationContractError("non-selected ablation counter groups must remain zero: " + ", ".join(invalid))

    @property
    def total_opportunities(self) -> int:
        return sum(self.for_key(key).opportunities for key in InterventionKey)

    @property
    def total_interventions(self) -> int:
        return sum(self.for_key(key).interventions for key in InterventionKey)

    @property
    def total_changed(self) -> int:
        return sum(self.for_key(key).changed for key in InterventionKey)

    def to_document(self) -> dict[str, dict[str, int]]:
        """Return every counter group, including required zero groups."""

        return {key.value: self.for_key(key).to_document() for key in InterventionKey}

    @classmethod
    def from_document(cls, document: object) -> AblationCounters:
        """Parse the complete closed counter schema."""

        if not isinstance(document, Mapping):
            raise AblationContractError("counters must be an object")
        required = {key.value for key in InterventionKey}
        if set(document) != required:
            raise AblationContractError("counters must contain exactly the registered intervention keys")
        return cls(
            **{
                key.value: InterventionCounts.from_document(
                    document[key.value],
                    context=f"counters.{key.value}",
                )
                for key in InterventionKey
            }
        )

    @classmethod
    def for_ablation(
        cls,
        ablation_id: AblationID,
        *,
        opportunities: int = 0,
        interventions: int = 0,
        changed: int = 0,
    ) -> AblationCounters:
        """Build a counter set with activity only in the selected group."""

        counts = InterventionCounts(
            opportunities=opportunities,
            interventions=interventions,
            changed=changed,
        )
        keys = ablation_policy(ablation_id).active_intervention_keys
        if not keys:
            if not counts.is_zero:
                raise AblationContractError("A0 cannot record an intervention")
            return cls()
        if len(keys) != 1:
            raise RuntimeError("one selected counter group is required per arm")
        return cls(**{keys[0].value: counts})


@dataclass(frozen=True, slots=True)
class AblationEvents:
    """Validated content of one future ``ablation_events.json`` artifact."""

    ablation_id: AblationID
    counters: AblationCounters

    def __post_init__(self) -> None:
        if not isinstance(self.ablation_id, AblationID):
            raise TypeError("ablation_id must be an AblationID")
        if type(self.counters) is not AblationCounters:
            raise TypeError("counters must be AblationCounters")
        self.counters.validate_for(self.ablation_id)

    def to_document(self) -> dict[str, object]:
        """Return the exact JSON document authenticated for this run."""

        return {
            "schema": ABLATION_EVENTS_SCHEMA,
            "contract_sha256": ABLATION_CONTRACT_SHA256,
            "ablation_id": self.ablation_id.value,
            "counters": self.counters.to_document(),
        }

    def to_bytes(self) -> bytes:
        """Return the sole canonical UTF-8 JSON representation."""

        return _canonical_json_bytes(self.to_document())

    @classmethod
    def from_document(cls, document: object) -> AblationEvents:
        """Validate a decoded ablation-events document."""

        if not isinstance(document, Mapping):
            raise AblationContractError("ablation events document must be an object")
        required = {"schema", "contract_sha256", "ablation_id", "counters"}
        if set(document) != required:
            raise AblationContractError("ablation events document has missing or unknown fields")
        if document["schema"] != ABLATION_EVENTS_SCHEMA:
            raise AblationContractError("unsupported ablation events schema")
        if document["contract_sha256"] != ABLATION_CONTRACT_SHA256:
            raise AblationContractError("ablation contract SHA256 mismatch")
        raw_id = document["ablation_id"]
        if type(raw_id) is not str:
            raise AblationContractError("ablation_id must be a string")
        try:
            ablation_id = AblationID(raw_id)
        except ValueError as error:
            raise AblationContractError(f"unknown ablation_id: {raw_id!r}") from error
        counters = AblationCounters.from_document(document["counters"])
        return cls(ablation_id=ablation_id, counters=counters)

    @classmethod
    def from_bytes(cls, content: bytes) -> AblationEvents:
        """Validate exact canonical artifact bytes and return typed content."""

        if type(content) is not bytes:
            raise TypeError("ablation events content must be bytes")
        if len(content) > MAX_ABLATION_EVENTS_BYTES:
            raise AblationContractError("ablation events document is too large")
        try:
            document = json.loads(content.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError, RecursionError) as error:
            raise AblationContractError("ablation events content is not valid UTF-8 JSON") from error
        events = cls.from_document(document)
        if events.to_bytes() != content:
            raise AblationContractError("ablation events bytes are not in canonical form")
        return events


def validate_ablation_events_document(document: object) -> AblationEvents:
    """Validate a decoded event document using the closed current contract."""

    return AblationEvents.from_document(document)


def validate_ablation_events_bytes(content: bytes) -> AblationEvents:
    """Validate canonical event bytes using the closed current contract."""

    return AblationEvents.from_bytes(content)
