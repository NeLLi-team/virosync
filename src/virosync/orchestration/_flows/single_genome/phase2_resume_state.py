"""Lossless JSON persistence for the complete Phase-2 resume state.

Phase-2 TSV and BED files are reporting artifacts.  They round numeric values
and omit fields that Phase 3 consumes, so an authenticated resume must use this
closed checkpoint instead.  Mapping values are encoded as ordered entries to
preserve their insertion order as well as their contents.
"""

from __future__ import annotations

import json
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, fields
from numbers import Integral, Real
from pathlib import Path

from virosync.pipeline.phase2.boundary_diamond import (
    ControlStats,
    GeneTaxonomy,
    GenomeDiamondQuery,
    SeedGeneMapping,
)
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.taxonomy_utils import TaxonomyFingerprint
from virosync.utils.atomic_write import atomic_write

from .phase_state import phase2_state_from_document, phase2_state_to_document


PHASE2_RESUME_STATE_FILENAME = "resume_state.json"
PHASE2_RESUME_STATE_SCHEMA_VERSION = 1
PHASE2_RESUME_STATE_ARTIFACT_TYPE = "virosync.phase2.resume_state"
PHASE2_RESUME_STATE_SCHEMA = (
    f"{PHASE2_RESUME_STATE_ARTIFACT_TYPE}"
    f"/v{PHASE2_RESUME_STATE_SCHEMA_VERSION}"
)

_TOP_LEVEL_FIELDS = {
    "artifact_type",
    "schema_version",
    "refined_boundaries",
    "boundary_taxonomy_map",
    "boundary_control_stats",
    "boundary_diamond_query",
}
_ENTRY_FIELDS = {"key", "value"}

_GENE_STRING_FIELDS = (
    "porf_id",
    "scaffold",
    "top1_target",
    "top1_prefix",
)
_GENE_INTEGER_FIELDS = ("start", "end")
_GENE_FLOAT_FIELDS = ("top1_pident", "top1_evalue")
_GENE_STRING_LIST_FIELDS = ("top10_prefixes", "top10_targets")
_GENE_FLOAT_LIST_FIELDS = ("top10_bits", "top10_pidents", "top10_evalues")
_GENE_BOOLEAN_FIELDS = (
    "has_ncldv_mirus",
    "has_vp_plv",
    "has_viral",
    "has_hit",
)
_GENE_SPECIAL_FIELDS = ("taxonomy_fingerprint",)
_GENE_FIELDS = (
    *_GENE_STRING_FIELDS,
    *_GENE_INTEGER_FIELDS,
    *_GENE_FLOAT_FIELDS,
    *_GENE_STRING_LIST_FIELDS,
    *_GENE_FLOAT_LIST_FIELDS,
    *_GENE_BOOLEAN_FIELDS,
    *_GENE_SPECIAL_FIELDS,
)
_GENE_FIELD_SET = set(_GENE_FIELDS)

_FINGERPRINT_FIELDS = ("weighted_tokens", "raw_tokens")
_FINGERPRINT_FIELD_SET = set(_FINGERPRINT_FIELDS)

_CONTROL_INTEGER_FIELDS = ("n_genes", "n_no_hits")
_CONTROL_FLOAT_FIELDS = ("no_hit_frequency", "host_frequency", "mean_pident")
_CONTROL_STRING_FIELDS = ("dominant_organism", "host_prefix")
_CONTROL_FIELDS = (
    *_CONTROL_INTEGER_FIELDS,
    *_CONTROL_FLOAT_FIELDS,
    *_CONTROL_STRING_FIELDS,
)
_CONTROL_FIELD_SET = set(_CONTROL_FIELDS)

_SEED_STRING_FIELDS = ("seed_id", "scaffold")
_SEED_INTEGER_FIELDS = (
    "seed_start",
    "seed_end",
    "flank_start_idx",
    "flank_end_idx",
    "flank_start_bp",
    "flank_end_bp",
    "flank_genes_config",
)
_SEED_STRING_LIST_FIELDS = (
    "eve_porf_ids",
    "upstream_porf_ids",
    "downstream_porf_ids",
)
_SEED_FIELDS = (
    *_SEED_STRING_FIELDS,
    *_SEED_INTEGER_FIELDS,
    *_SEED_STRING_LIST_FIELDS,
)
_SEED_FIELD_SET = set(_SEED_FIELDS)

_QUERY_STRING_MAP_FIELDS = ("eve_porf_ids", "boundary_porf_ids")
_QUERY_STRING_LIST_FIELDS = ("control_porf_ids", "all_porf_ids")
_QUERY_SPECIAL_FIELDS = ("seed_gene_mappings",)
_QUERY_FIELDS = (
    *_QUERY_STRING_MAP_FIELDS,
    *_QUERY_STRING_LIST_FIELDS,
    *_QUERY_SPECIAL_FIELDS,
)
_QUERY_FIELD_SET = set(_QUERY_FIELDS)


class Phase2ResumeStateError(ValueError):
    """Raised when the Phase-2 resume checkpoint violates its schema."""


@dataclass(frozen=True)
class Phase2ResumeState:
    """Exact Phase-2 objects passed from the orchestrator into Phase 3."""

    refined_boundaries: list[RefinedBoundary]
    boundary_taxonomy_map: dict[str, GeneTaxonomy]
    boundary_control_stats: ControlStats | None
    boundary_diamond_query: GenomeDiamondQuery | None


def _require_exact_fields(
    value: object,
    expected: set[str],
    context: str,
) -> dict[str, object]:
    if not isinstance(value, dict):
        raise Phase2ResumeStateError(f"{context} must be a JSON object")
    actual = set(value)
    if actual != expected:
        raise Phase2ResumeStateError(
            f"{context} fields differ from the schema; "
            f"missing={sorted(expected - actual)}, "
            f"extra={sorted(actual - expected)}"
        )
    return value


def _require_string(value: object, context: str) -> str:
    if not isinstance(value, str):
        raise Phase2ResumeStateError(f"{context} must be a string")
    return value


def _require_integer(value: object, context: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise Phase2ResumeStateError(f"{context} must be an integer")
    return int(value)


def _require_finite_float(value: object, context: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise Phase2ResumeStateError(f"{context} must be a number")
    result = float(value)
    if not math.isfinite(result):
        raise Phase2ResumeStateError(f"{context} must be finite")
    return result


def _require_boolean(value: object, context: str) -> bool:
    if not isinstance(value, bool):
        raise Phase2ResumeStateError(f"{context} must be a boolean")
    return value


def _require_string_list(value: object, context: str) -> list[str]:
    if not isinstance(value, list):
        raise Phase2ResumeStateError(f"{context} must be a list")
    return [
        _require_string(item, f"{context}[{index}]")
        for index, item in enumerate(value)
    ]


def _require_float_list(value: object, context: str) -> list[float]:
    if not isinstance(value, list):
        raise Phase2ResumeStateError(f"{context} must be a list")
    return [
        _require_finite_float(item, f"{context}[{index}]")
        for index, item in enumerate(value)
    ]


def _validate_model_fields() -> None:
    expected = (
        (GeneTaxonomy, _GENE_FIELD_SET),
        (TaxonomyFingerprint, _FINGERPRINT_FIELD_SET),
        (ControlStats, _CONTROL_FIELD_SET),
        (SeedGeneMapping, _SEED_FIELD_SET),
        (GenomeDiamondQuery, _QUERY_FIELD_SET),
    )
    for model, expected_fields in expected:
        actual_fields = {item.name for item in fields(model)}
        if actual_fields != expected_fields:
            raise Phase2ResumeStateError(
                f"{model.__name__} fields changed without a Phase-2 resume "
                "state schema update"
            )


def _require_exact_instance(value: object, model: type, context: str) -> None:
    if type(value) is not model:
        raise Phase2ResumeStateError(f"{context} must be {model.__name__}")
    expected = {item.name for item in fields(model)}
    if set(vars(value)) != expected:
        raise Phase2ResumeStateError(
            f"{context} contains unsupported dynamic fields"
        )


def _mapping_to_entries(
    value: object,
    encode_value: Callable[[object, str], object],
    context: str,
) -> list[dict[str, object]]:
    if type(value) is not dict:
        raise Phase2ResumeStateError(f"{context} must be a dictionary")
    entries: list[dict[str, object]] = []
    for index, (key, item) in enumerate(value.items()):
        entry_context = f"{context}[{index}]"
        entries.append(
            {
                "key": _require_string(key, f"{entry_context}.key"),
                "value": encode_value(item, f"{entry_context}.value"),
            }
        )
    return entries


def _mapping_from_entries(
    value: object,
    decode_value: Callable[[object, str], object],
    context: str,
) -> dict[str, object]:
    if not isinstance(value, list):
        raise Phase2ResumeStateError(f"{context} must be a list of entries")
    result: dict[str, object] = {}
    for index, raw_entry in enumerate(value):
        entry_context = f"{context}[{index}]"
        entry = _require_exact_fields(raw_entry, _ENTRY_FIELDS, entry_context)
        key = _require_string(entry["key"], f"{entry_context}.key")
        if key in result:
            raise Phase2ResumeStateError(
                f"{context} contains duplicate key {key!r}"
            )
        result[key] = decode_value(entry["value"], f"{entry_context}.value")
    return result




def _integer_to_document(value: object, context: str) -> int:
    return _require_integer(value, context)


def _float_to_document(value: object, context: str) -> float:
    return _require_finite_float(value, context)


def _string_list_to_document(value: object, context: str) -> list[str]:
    return _require_string_list(value, context)


def _fingerprint_to_document(
    value: object,
    context: str,
) -> dict[str, object] | None:
    if value is None:
        return None
    _require_exact_instance(value, TaxonomyFingerprint, context)
    return {
        "weighted_tokens": _mapping_to_entries(
            value.weighted_tokens,
            _float_to_document,
            f"{context}.weighted_tokens",
        ),
        "raw_tokens": _mapping_to_entries(
            value.raw_tokens,
            _integer_to_document,
            f"{context}.raw_tokens",
        ),
    }


def _fingerprint_from_document(
    value: object,
    context: str,
) -> TaxonomyFingerprint | None:
    if value is None:
        return None
    document = _require_exact_fields(value, _FINGERPRINT_FIELD_SET, context)
    weighted_tokens = _mapping_from_entries(
        document["weighted_tokens"],
        _require_finite_float,
        f"{context}.weighted_tokens",
    )
    raw_tokens = _mapping_from_entries(
        document["raw_tokens"],
        _require_integer,
        f"{context}.raw_tokens",
    )
    return TaxonomyFingerprint(
        weighted_tokens=dict(weighted_tokens),
        raw_tokens=dict(raw_tokens),
    )


def _gene_to_document(value: object, context: str) -> dict[str, object]:
    _require_exact_instance(value, GeneTaxonomy, context)
    document: dict[str, object] = {}
    for name in _GENE_STRING_FIELDS:
        document[name] = _require_string(
            getattr(value, name), f"{context}.{name}"
        )
    for name in _GENE_INTEGER_FIELDS:
        document[name] = _require_integer(
            getattr(value, name), f"{context}.{name}"
        )
    for name in _GENE_FLOAT_FIELDS:
        document[name] = _require_finite_float(
            getattr(value, name), f"{context}.{name}"
        )
    for name in _GENE_STRING_LIST_FIELDS:
        document[name] = _require_string_list(
            getattr(value, name), f"{context}.{name}"
        )
    for name in _GENE_FLOAT_LIST_FIELDS:
        document[name] = _require_float_list(
            getattr(value, name), f"{context}.{name}"
        )
    for name in _GENE_BOOLEAN_FIELDS:
        document[name] = _require_boolean(
            getattr(value, name), f"{context}.{name}"
        )
    document["taxonomy_fingerprint"] = _fingerprint_to_document(
        value.taxonomy_fingerprint,
        f"{context}.taxonomy_fingerprint",
    )
    return document


def _gene_from_document(value: object, context: str) -> GeneTaxonomy:
    document = _require_exact_fields(value, _GENE_FIELD_SET, context)
    kwargs: dict[str, object] = {}
    for name in _GENE_STRING_FIELDS:
        kwargs[name] = _require_string(document[name], f"{context}.{name}")
    for name in _GENE_INTEGER_FIELDS:
        kwargs[name] = _require_integer(document[name], f"{context}.{name}")
    for name in _GENE_FLOAT_FIELDS:
        kwargs[name] = _require_finite_float(
            document[name], f"{context}.{name}"
        )
    for name in _GENE_STRING_LIST_FIELDS:
        kwargs[name] = _require_string_list(
            document[name], f"{context}.{name}"
        )
    for name in _GENE_FLOAT_LIST_FIELDS:
        kwargs[name] = _require_float_list(document[name], f"{context}.{name}")
    for name in _GENE_BOOLEAN_FIELDS:
        kwargs[name] = _require_boolean(document[name], f"{context}.{name}")
    kwargs["taxonomy_fingerprint"] = _fingerprint_from_document(
        document["taxonomy_fingerprint"],
        f"{context}.taxonomy_fingerprint",
    )
    return GeneTaxonomy(**kwargs)


def _control_to_document(
    value: object,
    context: str,
) -> dict[str, object] | None:
    if value is None:
        return None
    _require_exact_instance(value, ControlStats, context)
    document: dict[str, object] = {}
    for name in _CONTROL_INTEGER_FIELDS:
        document[name] = _require_integer(
            getattr(value, name), f"{context}.{name}"
        )
    for name in _CONTROL_FLOAT_FIELDS:
        document[name] = _require_finite_float(
            getattr(value, name), f"{context}.{name}"
        )
    for name in _CONTROL_STRING_FIELDS:
        document[name] = _require_string(
            getattr(value, name), f"{context}.{name}"
        )
    return document


def _control_from_document(value: object, context: str) -> ControlStats | None:
    if value is None:
        return None
    document = _require_exact_fields(value, _CONTROL_FIELD_SET, context)
    kwargs: dict[str, object] = {}
    for name in _CONTROL_INTEGER_FIELDS:
        kwargs[name] = _require_integer(document[name], f"{context}.{name}")
    for name in _CONTROL_FLOAT_FIELDS:
        kwargs[name] = _require_finite_float(
            document[name], f"{context}.{name}"
        )
    for name in _CONTROL_STRING_FIELDS:
        kwargs[name] = _require_string(document[name], f"{context}.{name}")
    return ControlStats(**kwargs)


def _seed_mapping_to_document(
    value: object,
    context: str,
) -> dict[str, object]:
    _require_exact_instance(value, SeedGeneMapping, context)
    document: dict[str, object] = {}
    for name in _SEED_STRING_FIELDS:
        document[name] = _require_string(
            getattr(value, name), f"{context}.{name}"
        )
    for name in _SEED_INTEGER_FIELDS:
        document[name] = _require_integer(
            getattr(value, name), f"{context}.{name}"
        )
    for name in _SEED_STRING_LIST_FIELDS:
        document[name] = _require_string_list(
            getattr(value, name), f"{context}.{name}"
        )
    return document


def _seed_mapping_from_document(
    value: object,
    context: str,
) -> SeedGeneMapping:
    document = _require_exact_fields(value, _SEED_FIELD_SET, context)
    kwargs: dict[str, object] = {}
    for name in _SEED_STRING_FIELDS:
        kwargs[name] = _require_string(document[name], f"{context}.{name}")
    for name in _SEED_INTEGER_FIELDS:
        kwargs[name] = _require_integer(document[name], f"{context}.{name}")
    for name in _SEED_STRING_LIST_FIELDS:
        kwargs[name] = _require_string_list(
            document[name], f"{context}.{name}"
        )
    return SeedGeneMapping(**kwargs)


def _query_to_document(
    value: object,
    context: str,
) -> dict[str, object] | None:
    if value is None:
        return None
    _require_exact_instance(value, GenomeDiamondQuery, context)
    document: dict[str, object] = {}
    for name in _QUERY_STRING_MAP_FIELDS:
        document[name] = _mapping_to_entries(
            getattr(value, name),
            _string_list_to_document,
            f"{context}.{name}",
        )
    for name in _QUERY_STRING_LIST_FIELDS:
        document[name] = _require_string_list(
            getattr(value, name), f"{context}.{name}"
        )
    document["seed_gene_mappings"] = _mapping_to_entries(
        value.seed_gene_mappings,
        _seed_mapping_to_document,
        f"{context}.seed_gene_mappings",
    )
    return document


def _query_from_document(
    value: object,
    context: str,
) -> GenomeDiamondQuery | None:
    if value is None:
        return None
    document = _require_exact_fields(value, _QUERY_FIELD_SET, context)
    kwargs: dict[str, object] = {}
    for name in _QUERY_STRING_MAP_FIELDS:
        kwargs[name] = dict(
            _mapping_from_entries(
                document[name],
                _require_string_list,
                f"{context}.{name}",
            )
        )
    for name in _QUERY_STRING_LIST_FIELDS:
        kwargs[name] = _require_string_list(
            document[name], f"{context}.{name}"
        )
    kwargs["seed_gene_mappings"] = dict(
        _mapping_from_entries(
            document["seed_gene_mappings"],
            _seed_mapping_from_document,
            f"{context}.seed_gene_mappings",
        )
    )
    return GenomeDiamondQuery(**kwargs)


def phase2_resume_state_to_document(
    state: Phase2ResumeState,
) -> dict[str, object]:
    """Convert exact Phase-2 inputs into the closed JSON-compatible schema."""

    _validate_model_fields()
    _require_exact_instance(state, Phase2ResumeState, "phase2 resume state")
    if not isinstance(state.refined_boundaries, Sequence) or isinstance(
        state.refined_boundaries, (str, bytes)
    ):
        raise Phase2ResumeStateError("refined_boundaries must be a sequence")
    return {
        "artifact_type": PHASE2_RESUME_STATE_ARTIFACT_TYPE,
        "schema_version": PHASE2_RESUME_STATE_SCHEMA_VERSION,
        "refined_boundaries": phase2_state_to_document(
            state.refined_boundaries
        ),
        "boundary_taxonomy_map": _mapping_to_entries(
            state.boundary_taxonomy_map,
            _gene_to_document,
            "boundary_taxonomy_map",
        ),
        "boundary_control_stats": _control_to_document(
            state.boundary_control_stats,
            "boundary_control_stats",
        ),
        "boundary_diamond_query": _query_to_document(
            state.boundary_diamond_query,
            "boundary_diamond_query",
        ),
    }


def phase2_resume_state_from_document(document: object) -> Phase2ResumeState:
    """Validate a decoded document and reconstruct all Phase-3 inputs."""

    _validate_model_fields()
    payload = _require_exact_fields(
        document,
        _TOP_LEVEL_FIELDS,
        "phase2 resume state",
    )
    schema_version = payload["schema_version"]
    valid_schema_version = (
        not isinstance(schema_version, bool)
        and isinstance(schema_version, int)
        and schema_version == PHASE2_RESUME_STATE_SCHEMA_VERSION
    )
    if not valid_schema_version:
        raise Phase2ResumeStateError(
            f"unsupported Phase-2 resume schema_version: {schema_version!r}"
        )
    artifact_type = _require_string(
        payload["artifact_type"], "phase2 resume state.artifact_type"
    )
    if artifact_type != PHASE2_RESUME_STATE_ARTIFACT_TYPE:
        raise Phase2ResumeStateError(
            f"unsupported Phase-2 resume artifact_type: {artifact_type!r}"
        )
    taxonomy_map = _mapping_from_entries(
        payload["boundary_taxonomy_map"],
        _gene_from_document,
        "boundary_taxonomy_map",
    )
    return Phase2ResumeState(
        refined_boundaries=phase2_state_from_document(
            payload["refined_boundaries"]
        ),
        boundary_taxonomy_map=dict(taxonomy_map),
        boundary_control_stats=_control_from_document(
            payload["boundary_control_stats"],
            "boundary_control_stats",
        ),
        boundary_diamond_query=_query_from_document(
            payload["boundary_diamond_query"],
            "boundary_diamond_query",
        ),
    )


def _object_without_duplicate_keys(
    pairs: list[tuple[str, object]],
) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise Phase2ResumeStateError(
                f"duplicate JSON key in Phase-2 resume state: {key!r}"
            )
        result[key] = value
    return result


def _reject_nonstandard_constant(value: str) -> object:
    raise Phase2ResumeStateError(
        f"non-finite JSON number in Phase-2 resume state: {value}"
    )


def write_phase2_resume_state(
    path: str | Path,
    *,
    refined_boundaries: Sequence[RefinedBoundary],
    boundary_taxonomy_map: dict[str, GeneTaxonomy],
    boundary_control_stats: ControlStats | None,
    boundary_diamond_query: GenomeDiamondQuery | None,
) -> None:
    """Atomically write the complete Phase-2 resume checkpoint."""

    state = Phase2ResumeState(
        refined_boundaries=list(refined_boundaries),
        boundary_taxonomy_map=boundary_taxonomy_map,
        boundary_control_stats=boundary_control_stats,
        boundary_diamond_query=boundary_diamond_query,
    )
    document = phase2_resume_state_to_document(state)
    try:
        content = json.dumps(
            document,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
    except (TypeError, ValueError) as exc:
        raise Phase2ResumeStateError(
            f"Phase-2 resume state is not JSON-safe: {exc}"
        ) from exc
    atomic_write(path, content + "\n")


def load_phase2_resume_state(path: str | Path) -> Phase2ResumeState:
    """Load and strictly validate the complete Phase-2 resume checkpoint."""

    state_path = Path(path)
    try:
        content = state_path.read_text(encoding="utf-8")
        document = json.loads(
            content,
            object_pairs_hook=_object_without_duplicate_keys,
            parse_constant=_reject_nonstandard_constant,
        )
    except Phase2ResumeStateError:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise Phase2ResumeStateError(
            f"invalid Phase-2 resume state {state_path}: {exc}"
        ) from exc
    return phase2_resume_state_from_document(document)
