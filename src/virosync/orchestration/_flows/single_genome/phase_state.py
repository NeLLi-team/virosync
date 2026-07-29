"""Lossless JSON persistence for Phase-2 refined boundary state.

The BED export is a reporting artifact and cannot reconstruct the evidence that
Phase 3 consumes.  This module stores that evidence with an explicit, closed
schema.  It intentionally supports only JSON primitives, ``RefinedBoundary``,
``WindowFeatures``, and numeric posterior arrays.
"""

from __future__ import annotations

import json
import math
from collections.abc import Sequence
from dataclasses import fields
from numbers import Integral, Real
from pathlib import Path

import numpy as np

from virosync.features.compositional import WindowFeatures
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.utils.atomic_write import atomic_write


PHASE2_STATE_FILENAME = "refined_state.json"
PHASE2_STATE_SCHEMA_VERSION = 2
PHASE2_STATE_ARTIFACT_TYPE = "virosync.phase2.refined_boundaries"
PHASE2_STATE_SCHEMA = (
    f"{PHASE2_STATE_ARTIFACT_TYPE}/v{PHASE2_STATE_SCHEMA_VERSION}"
)

_TOP_LEVEL_FIELDS = {"artifact_type", "schema_version", "boundaries"}
_POSTERIOR_FIELDS = {"dtype", "shape", "data"}
_POSTERIOR_DTYPES = {"float16", "float32", "float64"}

_BOUNDARY_STRING_FIELDS = (
    "scaffold",
    "seed_id",
    "host_trim_reason",
    "host_trim_common_euk_taxonomy",
    "seed_confidence",
    "predicted_family",
)
_BOUNDARY_INTEGER_FIELDS = (
    "start",
    "end",
    "original_start",
    "original_end",
    "region_classification_ncldv_markers",
    "region_classification_vp_plv_markers",
    "region_classification_mirus_markers",
)
_BOUNDARY_OPTIONAL_INTEGER_FIELDS = (
    "candidate_start",
    "candidate_end",
    "core_viral_start",
    "core_viral_end",
    "flank_5_start",
    "flank_5_end",
    "flank_3_start",
    "flank_3_end",
    "marker_floor_start",
    "marker_floor_end",
)
_BOUNDARY_FLOAT_FIELDS = (
    "seed_hhg_score",
    "seed_novelty_score",
    "seed_compositional_score",
    "confidence",
    "posterior_probability",
    "max_kfd",
    "gc_deviation",
    "cub_deviation",
    "mean_novelty",
)
_BOUNDARY_STRING_LIST_FIELDS = ("seed_sources", "hallmark_genes")
_BOUNDARY_INTEGER_LIST_FIELDS = ("state_sequence",)
_BOUNDARY_BOOLEAN_FIELDS = ("seed_has_mcp",)
_BOUNDARY_SPECIAL_FIELDS = ("state_posteriors", "window_features")
_BOUNDARY_FIELDS = (
    *_BOUNDARY_STRING_FIELDS,
    *_BOUNDARY_INTEGER_FIELDS,
    *_BOUNDARY_OPTIONAL_INTEGER_FIELDS,
    *_BOUNDARY_FLOAT_FIELDS,
    *_BOUNDARY_STRING_LIST_FIELDS,
    *_BOUNDARY_INTEGER_LIST_FIELDS,
    *_BOUNDARY_BOOLEAN_FIELDS,
    *_BOUNDARY_SPECIAL_FIELDS,
)
_BOUNDARY_FIELD_SET = set(_BOUNDARY_FIELDS)

_WINDOW_STRING_FIELDS = ("scaffold",)
_WINDOW_INTEGER_FIELDS = ("start", "end", "porf_count")
_WINDOW_FLOAT_FIELDS = (
    "kfd",
    "cub",
    "gc_content",
    "gc_deviation",
    "porf_density",
)
_WINDOW_FIELDS = (
    *_WINDOW_STRING_FIELDS,
    *_WINDOW_INTEGER_FIELDS,
    *_WINDOW_FLOAT_FIELDS,
)
_WINDOW_FIELD_SET = set(_WINDOW_FIELDS)


class Phase2StateError(ValueError):
    """Raised when Phase-2 resume state does not match the closed schema."""


def _require_exact_fields(
    value: object,
    expected: set[str],
    context: str,
) -> dict[str, object]:
    if not isinstance(value, dict):
        raise Phase2StateError(f"{context} must be a JSON object")
    actual = set(value)
    if actual != expected:
        raise Phase2StateError(
            f"{context} fields differ from the schema; "
            f"missing={sorted(expected - actual)}, extra={sorted(actual - expected)}"
        )
    return value


def _require_string(value: object, context: str) -> str:
    if not isinstance(value, str):
        raise Phase2StateError(f"{context} must be a string")
    return value


def _require_integer(value: object, context: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise Phase2StateError(f"{context} must be an integer")
    return int(value)


def _require_optional_integer(value: object, context: str) -> int | None:
    if value is None:
        return None
    return _require_integer(value, context)


def _require_finite_float(value: object, context: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise Phase2StateError(f"{context} must be a number")
    result = float(value)
    if not math.isfinite(result):
        raise Phase2StateError(f"{context} must be finite")
    return result


def _require_boolean(value: object, context: str) -> bool:
    if not isinstance(value, bool):
        raise Phase2StateError(f"{context} must be a boolean")
    return value


def _require_string_list(value: object, context: str) -> list[str]:
    if not isinstance(value, list):
        raise Phase2StateError(f"{context} must be a list")
    return [
        _require_string(item, f"{context}[{index}]")
        for index, item in enumerate(value)
    ]


def _require_integer_list(value: object, context: str) -> list[int]:
    if not isinstance(value, list):
        raise Phase2StateError(f"{context} must be a list")
    return [
        _require_integer(item, f"{context}[{index}]")
        for index, item in enumerate(value)
    ]


def _validate_model_fields() -> None:
    boundary_fields = {item.name for item in fields(RefinedBoundary)}
    if boundary_fields != _BOUNDARY_FIELD_SET:
        raise Phase2StateError(
            "RefinedBoundary fields changed without a Phase-2 state schema update"
        )
    window_fields = {item.name for item in fields(WindowFeatures)}
    if window_fields != _WINDOW_FIELD_SET:
        raise Phase2StateError(
            "WindowFeatures fields changed without a Phase-2 state schema update"
        )


def _posterior_to_document(value: object, context: str) -> dict[str, object] | None:
    if value is None:
        return None
    if type(value) is not np.ndarray:
        raise Phase2StateError(f"{context} must be a NumPy array or null")
    dtype_name = value.dtype.name
    if dtype_name not in _POSTERIOR_DTYPES:
        raise Phase2StateError(
            f"{context} dtype must be one of {sorted(_POSTERIOR_DTYPES)}"
        )
    if value.ndim != 2:
        raise Phase2StateError(f"{context} must be a two-dimensional array")
    if not np.isfinite(value).all():
        raise Phase2StateError(f"{context} must contain only finite values")
    return {
        "dtype": dtype_name,
        "shape": [int(dimension) for dimension in value.shape],
        "data": value.reshape(-1).tolist(),
    }


def _posterior_from_document(value: object, context: str) -> np.ndarray | None:
    if value is None:
        return None
    document = _require_exact_fields(value, _POSTERIOR_FIELDS, context)
    dtype_name = _require_string(document["dtype"], f"{context}.dtype")
    if dtype_name not in _POSTERIOR_DTYPES:
        raise Phase2StateError(
            f"{context}.dtype must be one of {sorted(_POSTERIOR_DTYPES)}"
        )
    raw_shape = document["shape"]
    if not isinstance(raw_shape, list) or len(raw_shape) != 2:
        raise Phase2StateError(f"{context}.shape must contain two dimensions")
    shape = tuple(
        _require_integer(dimension, f"{context}.shape[{index}]")
        for index, dimension in enumerate(raw_shape)
    )
    if any(dimension < 0 for dimension in shape):
        raise Phase2StateError(f"{context}.shape dimensions must be non-negative")

    raw_data = document["data"]
    if not isinstance(raw_data, list):
        raise Phase2StateError(f"{context}.data must be a flat list")
    expected_size = math.prod(shape)
    if len(raw_data) != expected_size:
        raise Phase2StateError(
            f"{context}.data length {len(raw_data)} does not match shape {shape}"
        )
    data = [
        _require_finite_float(item, f"{context}.data[{index}]")
        for index, item in enumerate(raw_data)
    ]
    try:
        with np.errstate(over="ignore", invalid="ignore"):
            result = np.asarray(data, dtype=np.dtype(dtype_name)).reshape(shape)
    except (TypeError, ValueError, OverflowError) as exc:
        raise Phase2StateError(
            f"{context}.data cannot be reshaped to {shape}"
        ) from exc
    if not np.isfinite(result).all():
        raise Phase2StateError(
            f"{context}.data cannot be represented as finite {dtype_name} values"
        )
    return result


def _window_to_document(window: object, context: str) -> dict[str, object]:
    if type(window) is not WindowFeatures:
        raise Phase2StateError(f"{context} must be WindowFeatures")
    if set(vars(window)) != _WINDOW_FIELD_SET:
        raise Phase2StateError(f"{context} contains unsupported dynamic fields")
    document: dict[str, object] = {}
    for name in _WINDOW_STRING_FIELDS:
        document[name] = _require_string(getattr(window, name), f"{context}.{name}")
    for name in _WINDOW_INTEGER_FIELDS:
        document[name] = _require_integer(getattr(window, name), f"{context}.{name}")
    for name in _WINDOW_FLOAT_FIELDS:
        document[name] = _require_finite_float(
            getattr(window, name), f"{context}.{name}"
        )
    return document


def _window_from_document(value: object, context: str) -> WindowFeatures:
    document = _require_exact_fields(value, _WINDOW_FIELD_SET, context)
    kwargs: dict[str, object] = {}
    for name in _WINDOW_STRING_FIELDS:
        kwargs[name] = _require_string(document[name], f"{context}.{name}")
    for name in _WINDOW_INTEGER_FIELDS:
        kwargs[name] = _require_integer(document[name], f"{context}.{name}")
    for name in _WINDOW_FLOAT_FIELDS:
        kwargs[name] = _require_finite_float(document[name], f"{context}.{name}")
    return WindowFeatures(**kwargs)


def _boundary_to_document(boundary: object, index: int) -> dict[str, object]:
    context = f"boundaries[{index}]"
    if type(boundary) is not RefinedBoundary:
        raise Phase2StateError(f"{context} must be RefinedBoundary")
    if set(vars(boundary)) != _BOUNDARY_FIELD_SET:
        raise Phase2StateError(f"{context} contains unsupported dynamic fields")

    document: dict[str, object] = {}
    for name in _BOUNDARY_STRING_FIELDS:
        document[name] = _require_string(
            getattr(boundary, name), f"{context}.{name}"
        )
    for name in _BOUNDARY_INTEGER_FIELDS:
        document[name] = _require_integer(
            getattr(boundary, name), f"{context}.{name}"
        )
    for name in _BOUNDARY_OPTIONAL_INTEGER_FIELDS:
        document[name] = _require_optional_integer(
            getattr(boundary, name), f"{context}.{name}"
        )
    for name in _BOUNDARY_FLOAT_FIELDS:
        document[name] = _require_finite_float(
            getattr(boundary, name), f"{context}.{name}"
        )
    for name in _BOUNDARY_STRING_LIST_FIELDS:
        document[name] = _require_string_list(
            getattr(boundary, name), f"{context}.{name}"
        )
    for name in _BOUNDARY_INTEGER_LIST_FIELDS:
        document[name] = _require_integer_list(
            getattr(boundary, name), f"{context}.{name}"
        )
    for name in _BOUNDARY_BOOLEAN_FIELDS:
        document[name] = _require_boolean(
            getattr(boundary, name), f"{context}.{name}"
        )
    document["state_posteriors"] = _posterior_to_document(
        boundary.state_posteriors,
        f"{context}.state_posteriors",
    )
    if not isinstance(boundary.window_features, list):
        raise Phase2StateError(f"{context}.window_features must be a list")
    document["window_features"] = [
        _window_to_document(window, f"{context}.window_features[{window_index}]")
        for window_index, window in enumerate(boundary.window_features)
    ]
    return document


def _boundary_from_document(value: object, index: int) -> RefinedBoundary:
    context = f"boundaries[{index}]"
    document = _require_exact_fields(value, _BOUNDARY_FIELD_SET, context)
    kwargs: dict[str, object] = {}
    for name in _BOUNDARY_STRING_FIELDS:
        kwargs[name] = _require_string(document[name], f"{context}.{name}")
    for name in _BOUNDARY_INTEGER_FIELDS:
        kwargs[name] = _require_integer(document[name], f"{context}.{name}")
    for name in _BOUNDARY_OPTIONAL_INTEGER_FIELDS:
        kwargs[name] = _require_optional_integer(
            document[name], f"{context}.{name}"
        )
    for name in _BOUNDARY_FLOAT_FIELDS:
        kwargs[name] = _require_finite_float(document[name], f"{context}.{name}")
    for name in _BOUNDARY_STRING_LIST_FIELDS:
        kwargs[name] = _require_string_list(document[name], f"{context}.{name}")
    for name in _BOUNDARY_INTEGER_LIST_FIELDS:
        kwargs[name] = _require_integer_list(document[name], f"{context}.{name}")
    for name in _BOUNDARY_BOOLEAN_FIELDS:
        kwargs[name] = _require_boolean(document[name], f"{context}.{name}")
    kwargs["state_posteriors"] = _posterior_from_document(
        document["state_posteriors"],
        f"{context}.state_posteriors",
    )
    raw_windows = document["window_features"]
    if not isinstance(raw_windows, list):
        raise Phase2StateError(f"{context}.window_features must be a list")
    kwargs["window_features"] = [
        _window_from_document(window, f"{context}.window_features[{window_index}]")
        for window_index, window in enumerate(raw_windows)
    ]
    return RefinedBoundary(**kwargs)


def phase2_state_to_document(
    boundaries: Sequence[RefinedBoundary],
) -> dict[str, object]:
    """Convert refined boundaries to the closed JSON-compatible schema."""

    _validate_model_fields()
    if isinstance(boundaries, (str, bytes)) or not isinstance(boundaries, Sequence):
        raise Phase2StateError("boundaries must be a sequence")
    return {
        "artifact_type": PHASE2_STATE_ARTIFACT_TYPE,
        "schema_version": PHASE2_STATE_SCHEMA_VERSION,
        "boundaries": [
            _boundary_to_document(boundary, index)
            for index, boundary in enumerate(boundaries)
        ],
    }


def phase2_state_from_document(document: object) -> list[RefinedBoundary]:
    """Validate a decoded JSON document and reconstruct refined boundaries."""

    _validate_model_fields()
    payload = _require_exact_fields(document, _TOP_LEVEL_FIELDS, "phase2 state")
    schema_version = payload["schema_version"]
    if (
        isinstance(schema_version, bool)
        or not isinstance(schema_version, int)
        or schema_version != PHASE2_STATE_SCHEMA_VERSION
    ):
        raise Phase2StateError(
            f"unsupported Phase-2 state schema_version: {schema_version!r}"
        )
    artifact_type = _require_string(
        payload["artifact_type"], "phase2 state.artifact_type"
    )
    if artifact_type != PHASE2_STATE_ARTIFACT_TYPE:
        raise Phase2StateError(
            f"unsupported Phase-2 state artifact_type: {artifact_type!r}"
        )
    raw_boundaries = payload["boundaries"]
    if not isinstance(raw_boundaries, list):
        raise Phase2StateError("phase2 state.boundaries must be a list")
    return [
        _boundary_from_document(boundary, index)
        for index, boundary in enumerate(raw_boundaries)
    ]


def _object_without_duplicate_keys(
    pairs: list[tuple[str, object]],
) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise Phase2StateError(f"duplicate JSON key in Phase-2 state: {key!r}")
        result[key] = value
    return result


def _reject_nonstandard_constant(value: str) -> object:
    raise Phase2StateError(f"non-finite JSON number in Phase-2 state: {value}")


def write_phase2_state(
    path: str | Path,
    boundaries: Sequence[RefinedBoundary],
) -> None:
    """Atomically write refined boundaries as canonical schema-v1 JSON."""

    document = phase2_state_to_document(boundaries)
    try:
        content = json.dumps(
            document,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
    except (TypeError, ValueError) as exc:
        raise Phase2StateError(f"Phase-2 state is not JSON-safe: {exc}") from exc
    atomic_write(path, content + "\n")


def load_phase2_state(path: str | Path) -> list[RefinedBoundary]:
    """Load and strictly validate a Phase-2 refined-boundary JSON artifact."""

    state_path = Path(path)
    try:
        content = state_path.read_text(encoding="utf-8")
        document = json.loads(
            content,
            object_pairs_hook=_object_without_duplicate_keys,
            parse_constant=_reject_nonstandard_constant,
        )
    except Phase2StateError:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise Phase2StateError(f"invalid Phase-2 state {state_path}: {exc}") from exc
    return phase2_state_from_document(document)
