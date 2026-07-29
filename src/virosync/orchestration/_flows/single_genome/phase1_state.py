"""Lossless JSON persistence for Phase-1 resume state.

The Phase-1 TSV files are reporting artifacts: numeric values are rounded and
several fields consumed by later phases are omitted.  This module persists the
exact in-memory state behind an explicit, closed JSON schema.
"""

from __future__ import annotations

import json
import math
from collections.abc import Sequence
from dataclasses import dataclass, fields
from numbers import Integral, Real
from pathlib import Path

from virosync.pipeline.host_signatures import HostSignatureModel
from virosync.pipeline.phase1.hhg_seeding import Anchor
from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit
from virosync.pipeline.phase1.seed_merger import MergedSeed
from virosync.utils.atomic_write import atomic_write


PHASE1_STATE_FILENAME = "resume_state.json"
PHASE1_STATE_SCHEMA_VERSION = 1
PHASE1_STATE_ARTIFACT_TYPE = "virosync.phase1.resume_state"
PHASE1_STATE_SCHEMA = (
    f"{PHASE1_STATE_ARTIFACT_TYPE}/v{PHASE1_STATE_SCHEMA_VERSION}"
)

_TOP_LEVEL_FIELDS = {
    "artifact_type",
    "schema_version",
    "validated_markers",
    "merged_seeds",
    "host_signature_model",
    "host_signatures",
    "host_deviation_summary",
}

_MARKER_STRING_FIELDS = (
    "query_porf",
    "scaffold",
    "strand",
    "hmm_target",
    "validation_status",
    "top10_prefixes",
    "best_hit_target",
    "top10_targets",
    "top10_pidents",
    "top10_bitscores",
    "top10_evalues",
    "taxonomy_substring_counts",
    "taxonomy_raw_counts",
)
_MARKER_INTEGER_FIELDS = (
    "start",
    "end",
    "has_ncldv",
    "has_mirus",
    "has_plv",
    "has_vp",
    "has_viral",
)
_MARKER_FLOAT_FIELDS = (
    "hmm_score",
    "hmm_evalue",
    "best_hit_pident",
    "best_hit_bits",
)
_MARKER_FIELDS = (
    *_MARKER_STRING_FIELDS,
    *_MARKER_INTEGER_FIELDS,
    *_MARKER_FLOAT_FIELDS,
)
_MARKER_FIELD_SET = set(_MARKER_FIELDS)

_ANCHOR_STRING_FIELDS = (
    "porf_id",
    "scaffold",
    "strand",
    "hallmark_gene",
)
_ANCHOR_INTEGER_FIELDS = ("start", "end")
_ANCHOR_FLOAT_FIELDS = ("score", "evalue")
_ANCHOR_FIELDS = (
    *_ANCHOR_STRING_FIELDS,
    *_ANCHOR_INTEGER_FIELDS,
    *_ANCHOR_FLOAT_FIELDS,
)
_ANCHOR_FIELD_SET = set(_ANCHOR_FIELDS)

_SEED_STRING_FIELDS = (
    "scaffold",
    "seed_id",
    "confidence",
    "predicted_family",
    "host_trim_reason",
    "host_trim_common_euk_taxonomy",
)
_SEED_INTEGER_FIELDS = (
    "start",
    "end",
    "n_windows",
    "region_classification_ncldv_markers",
    "region_classification_vp_plv_markers",
    "region_classification_mirus_markers",
)
_SEED_OPTIONAL_INTEGER_FIELDS = (
    "host_trim_original_start",
    "host_trim_original_end",
    "host_trimmed_start",
    "host_trimmed_end",
)
_SEED_FLOAT_FIELDS = (
    "hhg_score",
    "novelty_score",
    "compositional_score",
    "mean_kfd",
    "mean_composite",
    "max_kfd",
    "max_composite",
    "gc_deviation",
    "cub_deviation",
    "priority",
    "score",
)
_SEED_STRING_LIST_FIELDS = ("sources",)
_SEED_INTEGER_LIST_FIELDS = ("cluster_ids",)
_SEED_ANCHOR_LIST_FIELDS = ("anchors", "hhg_anchors")
_SEED_FIELDS = (
    *_SEED_STRING_FIELDS,
    *_SEED_INTEGER_FIELDS,
    *_SEED_OPTIONAL_INTEGER_FIELDS,
    *_SEED_FLOAT_FIELDS,
    *_SEED_STRING_LIST_FIELDS,
    *_SEED_INTEGER_LIST_FIELDS,
    *_SEED_ANCHOR_LIST_FIELDS,
)
_SEED_FIELD_SET = set(_SEED_FIELDS)

_HOST_MODEL_FIELDS = {
    "token_weights",
    "token_counts",
    "token_bits",
    "max_weight",
    "min_token_length",
    "host_prefixes",
    "weight_mode",
}
_HOST_DEVIATION_FIELDS = {
    "enabled",
    "markers_total",
    "markers_seedable",
    "baseline",
    "report_path",
}


class Phase1StateError(ValueError):
    """Raised when Phase-1 resume state does not match the closed schema."""


@dataclass(frozen=True)
class Phase1ResumeState:
    """Exact Phase-1 values consumed by Phase 2, Phase 3, and the run log."""

    validated_markers: list[ValidatedMarkerHit]
    merged_seeds: list[MergedSeed]
    host_signature_model: HostSignatureModel
    host_signatures: set[str]
    host_deviation_summary: dict[str, object] | None


def _require_exact_fields(
    value: object,
    expected: set[str],
    context: str,
) -> dict[str, object]:
    if not isinstance(value, dict):
        raise Phase1StateError(f"{context} must be a JSON object")
    actual = set(value)
    if actual != expected:
        raise Phase1StateError(
            f"{context} fields differ from the schema; "
            f"missing={sorted(expected - actual)}, extra={sorted(actual - expected)}"
        )
    return value


def _require_string(value: object, context: str) -> str:
    if not isinstance(value, str):
        raise Phase1StateError(f"{context} must be a string")
    return value


def _require_integer(value: object, context: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise Phase1StateError(f"{context} must be an integer")
    return int(value)


def _require_optional_integer(value: object, context: str) -> int | None:
    if value is None:
        return None
    return _require_integer(value, context)


def _require_finite_float(value: object, context: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise Phase1StateError(f"{context} must be a number")
    result = float(value)
    if not math.isfinite(result):
        raise Phase1StateError(f"{context} must be finite")
    return result


def _require_boolean(value: object, context: str) -> bool:
    if not isinstance(value, bool):
        raise Phase1StateError(f"{context} must be a boolean")
    return value


def _require_string_list(value: object, context: str) -> list[str]:
    if not isinstance(value, list):
        raise Phase1StateError(f"{context} must be a list")
    return [
        _require_string(item, f"{context}[{index}]")
        for index, item in enumerate(value)
    ]


def _require_integer_list(value: object, context: str) -> list[int]:
    if not isinstance(value, list):
        raise Phase1StateError(f"{context} must be a list")
    return [
        _require_integer(item, f"{context}[{index}]")
        for index, item in enumerate(value)
    ]


def _require_string_number_map(
    value: object,
    context: str,
    *,
    integer_values: bool,
) -> dict[str, int] | dict[str, float]:
    if not isinstance(value, dict):
        raise Phase1StateError(f"{context} must be a JSON object")
    result: dict[str, int] | dict[str, float] = {}
    for key, item in value.items():
        token = _require_string(key, f"{context} key")
        if integer_values:
            result[token] = _require_integer(item, f"{context}[{token!r}]")
        else:
            result[token] = _require_finite_float(
                item,
                f"{context}[{token!r}]",
            )
    return result


def _validate_model_fields() -> None:
    model_fields = (
        (ValidatedMarkerHit, _MARKER_FIELD_SET, "ValidatedMarkerHit"),
        (Anchor, _ANCHOR_FIELD_SET, "Anchor"),
        (MergedSeed, _SEED_FIELD_SET, "MergedSeed"),
        (HostSignatureModel, _HOST_MODEL_FIELDS, "HostSignatureModel"),
    )
    for model_type, expected, name in model_fields:
        actual = {item.name for item in fields(model_type)}
        if actual != expected:
            raise Phase1StateError(
                f"{name} fields changed without a Phase-1 state schema update"
            )


def _marker_to_document(marker: object, index: int) -> dict[str, object]:
    context = f"validated_markers[{index}]"
    if type(marker) is not ValidatedMarkerHit:
        raise Phase1StateError(f"{context} must be ValidatedMarkerHit")
    if set(vars(marker)) != _MARKER_FIELD_SET:
        raise Phase1StateError(f"{context} contains unsupported dynamic fields")
    document: dict[str, object] = {}
    for name in _MARKER_STRING_FIELDS:
        document[name] = _require_string(getattr(marker, name), f"{context}.{name}")
    for name in _MARKER_INTEGER_FIELDS:
        document[name] = _require_integer(getattr(marker, name), f"{context}.{name}")
    for name in _MARKER_FLOAT_FIELDS:
        document[name] = _require_finite_float(
            getattr(marker, name),
            f"{context}.{name}",
        )
    return document


def _marker_from_document(value: object, index: int) -> ValidatedMarkerHit:
    context = f"validated_markers[{index}]"
    document = _require_exact_fields(value, _MARKER_FIELD_SET, context)
    kwargs: dict[str, object] = {}
    for name in _MARKER_STRING_FIELDS:
        kwargs[name] = _require_string(document[name], f"{context}.{name}")
    for name in _MARKER_INTEGER_FIELDS:
        kwargs[name] = _require_integer(document[name], f"{context}.{name}")
    for name in _MARKER_FLOAT_FIELDS:
        kwargs[name] = _require_finite_float(document[name], f"{context}.{name}")
    return ValidatedMarkerHit(**kwargs)


def _anchor_to_document(anchor: object, context: str) -> dict[str, object]:
    if type(anchor) is not Anchor:
        raise Phase1StateError(f"{context} must be Anchor")
    if set(vars(anchor)) != _ANCHOR_FIELD_SET:
        raise Phase1StateError(f"{context} contains unsupported dynamic fields")
    document: dict[str, object] = {}
    for name in _ANCHOR_STRING_FIELDS:
        document[name] = _require_string(getattr(anchor, name), f"{context}.{name}")
    for name in _ANCHOR_INTEGER_FIELDS:
        document[name] = _require_integer(getattr(anchor, name), f"{context}.{name}")
    for name in _ANCHOR_FLOAT_FIELDS:
        document[name] = _require_finite_float(
            getattr(anchor, name),
            f"{context}.{name}",
        )
    return document


def _anchor_from_document(value: object, context: str) -> Anchor:
    document = _require_exact_fields(value, _ANCHOR_FIELD_SET, context)
    kwargs: dict[str, object] = {}
    for name in _ANCHOR_STRING_FIELDS:
        kwargs[name] = _require_string(document[name], f"{context}.{name}")
    for name in _ANCHOR_INTEGER_FIELDS:
        kwargs[name] = _require_integer(document[name], f"{context}.{name}")
    for name in _ANCHOR_FLOAT_FIELDS:
        kwargs[name] = _require_finite_float(document[name], f"{context}.{name}")
    return Anchor(**kwargs)


def _seed_to_document(seed: object, index: int) -> dict[str, object]:
    context = f"merged_seeds[{index}]"
    if type(seed) is not MergedSeed:
        raise Phase1StateError(f"{context} must be MergedSeed")
    if set(vars(seed)) != _SEED_FIELD_SET:
        raise Phase1StateError(f"{context} contains unsupported dynamic fields")
    document: dict[str, object] = {}
    for name in _SEED_STRING_FIELDS:
        document[name] = _require_string(getattr(seed, name), f"{context}.{name}")
    for name in _SEED_INTEGER_FIELDS:
        document[name] = _require_integer(getattr(seed, name), f"{context}.{name}")
    for name in _SEED_OPTIONAL_INTEGER_FIELDS:
        document[name] = _require_optional_integer(
            getattr(seed, name),
            f"{context}.{name}",
        )
    for name in _SEED_FLOAT_FIELDS:
        document[name] = _require_finite_float(
            getattr(seed, name),
            f"{context}.{name}",
        )
    for name in _SEED_STRING_LIST_FIELDS:
        document[name] = _require_string_list(
            getattr(seed, name),
            f"{context}.{name}",
        )
    for name in _SEED_INTEGER_LIST_FIELDS:
        document[name] = _require_integer_list(
            getattr(seed, name),
            f"{context}.{name}",
        )
    for name in _SEED_ANCHOR_LIST_FIELDS:
        anchors = getattr(seed, name)
        if not isinstance(anchors, list):
            raise Phase1StateError(f"{context}.{name} must be a list")
        document[name] = [
            _anchor_to_document(anchor, f"{context}.{name}[{anchor_index}]")
            for anchor_index, anchor in enumerate(anchors)
        ]
    return document


def _seed_from_document(value: object, index: int) -> MergedSeed:
    context = f"merged_seeds[{index}]"
    document = _require_exact_fields(value, _SEED_FIELD_SET, context)
    kwargs: dict[str, object] = {}
    for name in _SEED_STRING_FIELDS:
        kwargs[name] = _require_string(document[name], f"{context}.{name}")
    for name in _SEED_INTEGER_FIELDS:
        kwargs[name] = _require_integer(document[name], f"{context}.{name}")
    for name in _SEED_OPTIONAL_INTEGER_FIELDS:
        kwargs[name] = _require_optional_integer(
            document[name],
            f"{context}.{name}",
        )
    for name in _SEED_FLOAT_FIELDS:
        kwargs[name] = _require_finite_float(document[name], f"{context}.{name}")
    for name in _SEED_STRING_LIST_FIELDS:
        kwargs[name] = _require_string_list(document[name], f"{context}.{name}")
    for name in _SEED_INTEGER_LIST_FIELDS:
        kwargs[name] = _require_integer_list(document[name], f"{context}.{name}")
    for name in _SEED_ANCHOR_LIST_FIELDS:
        raw_anchors = document[name]
        if not isinstance(raw_anchors, list):
            raise Phase1StateError(f"{context}.{name} must be a list")
        kwargs[name] = [
            _anchor_from_document(
                anchor,
                f"{context}.{name}[{anchor_index}]",
            )
            for anchor_index, anchor in enumerate(raw_anchors)
        ]
    return MergedSeed(**kwargs)


def _host_model_to_document(model: object) -> dict[str, object]:
    context = "host_signature_model"
    if type(model) is not HostSignatureModel:
        raise Phase1StateError(f"{context} must be HostSignatureModel")
    if set(vars(model)) != _HOST_MODEL_FIELDS:
        raise Phase1StateError(f"{context} contains unsupported dynamic fields")
    token_bits: dict[str, list[float]] = {}
    if not isinstance(model.token_bits, dict):
        raise Phase1StateError(f"{context}.token_bits must be a JSON object")
    for token, raw_bits in model.token_bits.items():
        key = _require_string(token, f"{context}.token_bits key")
        if not isinstance(raw_bits, list):
            raise Phase1StateError(f"{context}.token_bits[{key!r}] must be a list")
        token_bits[key] = [
            _require_finite_float(bit, f"{context}.token_bits[{key!r}][{index}]")
            for index, bit in enumerate(raw_bits)
        ]
    return {
        "token_weights": _require_string_number_map(
            model.token_weights,
            f"{context}.token_weights",
            integer_values=False,
        ),
        "token_counts": _require_string_number_map(
            model.token_counts,
            f"{context}.token_counts",
            integer_values=True,
        ),
        "token_bits": token_bits,
        "max_weight": _require_finite_float(
            model.max_weight,
            f"{context}.max_weight",
        ),
        "min_token_length": _require_integer(
            model.min_token_length,
            f"{context}.min_token_length",
        ),
        "host_prefixes": _require_string_list(
            model.host_prefixes,
            f"{context}.host_prefixes",
        ),
        "weight_mode": _require_string(
            model.weight_mode,
            f"{context}.weight_mode",
        ),
    }


def _host_model_from_document(value: object) -> HostSignatureModel:
    context = "host_signature_model"
    document = _require_exact_fields(value, _HOST_MODEL_FIELDS, context)
    raw_bits = document["token_bits"]
    if not isinstance(raw_bits, dict):
        raise Phase1StateError(f"{context}.token_bits must be a JSON object")
    token_bits: dict[str, list[float]] = {}
    for token, bits in raw_bits.items():
        key = _require_string(token, f"{context}.token_bits key")
        if not isinstance(bits, list):
            raise Phase1StateError(f"{context}.token_bits[{key!r}] must be a list")
        token_bits[key] = [
            _require_finite_float(bit, f"{context}.token_bits[{key!r}][{index}]")
            for index, bit in enumerate(bits)
        ]
    return HostSignatureModel(
        token_weights=dict(
            _require_string_number_map(
                document["token_weights"],
                f"{context}.token_weights",
                integer_values=False,
            )
        ),
        token_counts=dict(
            _require_string_number_map(
                document["token_counts"],
                f"{context}.token_counts",
                integer_values=True,
            )
        ),
        token_bits=token_bits,
        max_weight=_require_finite_float(
            document["max_weight"],
            f"{context}.max_weight",
        ),
        min_token_length=_require_integer(
            document["min_token_length"],
            f"{context}.min_token_length",
        ),
        host_prefixes=_require_string_list(
            document["host_prefixes"],
            f"{context}.host_prefixes",
        ),
        weight_mode=_require_string(
            document["weight_mode"],
            f"{context}.weight_mode",
        ),
    )


def _json_safe_value(value: object, context: str) -> object:
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, Integral):
        return int(value)
    if isinstance(value, Real):
        return _require_finite_float(value, context)
    if isinstance(value, list):
        return [
            _json_safe_value(item, f"{context}[{index}]")
            for index, item in enumerate(value)
        ]
    if isinstance(value, dict):
        result: dict[str, object] = {}
        for key, item in value.items():
            name = _require_string(key, f"{context} key")
            result[name] = _json_safe_value(item, f"{context}[{name!r}]")
        return result
    raise Phase1StateError(f"{context} contains unsupported JSON value {type(value)}")


def _host_deviation_to_document(value: object) -> dict[str, object] | None:
    if value is None:
        return None
    context = "host_deviation_summary"
    document = _require_exact_fields(value, _HOST_DEVIATION_FIELDS, context)
    return {
        "enabled": _require_boolean(document["enabled"], f"{context}.enabled"),
        "markers_total": _require_integer(
            document["markers_total"],
            f"{context}.markers_total",
        ),
        "markers_seedable": _require_integer(
            document["markers_seedable"],
            f"{context}.markers_seedable",
        ),
        "baseline": _json_safe_value(
            document["baseline"],
            f"{context}.baseline",
        ),
        "report_path": _require_string(
            document["report_path"],
            f"{context}.report_path",
        ),
    }


def phase1_state_to_document(
    *,
    validated_markers: Sequence[ValidatedMarkerHit],
    merged_seeds: Sequence[MergedSeed],
    host_signature_model: HostSignatureModel,
    host_signatures: set[str],
    host_deviation_summary: dict[str, object] | None,
) -> dict[str, object]:
    """Convert Phase-1 values to the closed JSON-compatible schema."""

    _validate_model_fields()
    if isinstance(validated_markers, (str, bytes)) or not isinstance(
        validated_markers,
        Sequence,
    ):
        raise Phase1StateError("validated_markers must be a sequence")
    if isinstance(merged_seeds, (str, bytes)) or not isinstance(
        merged_seeds,
        Sequence,
    ):
        raise Phase1StateError("merged_seeds must be a sequence")
    if not isinstance(host_signatures, set):
        raise Phase1StateError("host_signatures must be a set")
    signatures = sorted(
        _require_string(item, "host_signatures item")
        for item in host_signatures
    )
    return {
        "artifact_type": PHASE1_STATE_ARTIFACT_TYPE,
        "schema_version": PHASE1_STATE_SCHEMA_VERSION,
        "validated_markers": [
            _marker_to_document(marker, index)
            for index, marker in enumerate(validated_markers)
        ],
        "merged_seeds": [
            _seed_to_document(seed, index)
            for index, seed in enumerate(merged_seeds)
        ],
        "host_signature_model": _host_model_to_document(host_signature_model),
        "host_signatures": signatures,
        "host_deviation_summary": _host_deviation_to_document(
            host_deviation_summary
        ),
    }


def phase1_state_from_document(document: object) -> Phase1ResumeState:
    """Validate a decoded JSON document and reconstruct Phase-1 values."""

    _validate_model_fields()
    payload = _require_exact_fields(document, _TOP_LEVEL_FIELDS, "phase1 state")
    schema_version = payload["schema_version"]
    if (
        isinstance(schema_version, bool)
        or not isinstance(schema_version, int)
        or schema_version != PHASE1_STATE_SCHEMA_VERSION
    ):
        raise Phase1StateError(
            f"unsupported Phase-1 state schema_version: {schema_version!r}"
        )
    artifact_type = _require_string(
        payload["artifact_type"],
        "phase1 state.artifact_type",
    )
    if artifact_type != PHASE1_STATE_ARTIFACT_TYPE:
        raise Phase1StateError(
            f"unsupported Phase-1 state artifact_type: {artifact_type!r}"
        )
    raw_markers = payload["validated_markers"]
    if not isinstance(raw_markers, list):
        raise Phase1StateError("phase1 state.validated_markers must be a list")
    raw_seeds = payload["merged_seeds"]
    if not isinstance(raw_seeds, list):
        raise Phase1StateError("phase1 state.merged_seeds must be a list")
    raw_signatures = payload["host_signatures"]
    signatures = _require_string_list(
        raw_signatures,
        "phase1 state.host_signatures",
    )
    if signatures != sorted(set(signatures)):
        raise Phase1StateError(
            "phase1 state.host_signatures must be sorted and unique"
        )
    return Phase1ResumeState(
        validated_markers=[
            _marker_from_document(marker, index)
            for index, marker in enumerate(raw_markers)
        ],
        merged_seeds=[
            _seed_from_document(seed, index)
            for index, seed in enumerate(raw_seeds)
        ],
        host_signature_model=_host_model_from_document(
            payload["host_signature_model"]
        ),
        host_signatures=set(signatures),
        host_deviation_summary=_host_deviation_to_document(
            payload["host_deviation_summary"]
        ),
    )


def _object_without_duplicate_keys(
    pairs: list[tuple[str, object]],
) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise Phase1StateError(f"duplicate JSON key in Phase-1 state: {key!r}")
        result[key] = value
    return result


def _reject_nonstandard_constant(value: str) -> object:
    raise Phase1StateError(f"non-finite JSON number in Phase-1 state: {value}")


def write_phase1_state(
    path: str | Path,
    *,
    validated_markers: Sequence[ValidatedMarkerHit],
    merged_seeds: Sequence[MergedSeed],
    host_signature_model: HostSignatureModel,
    host_signatures: set[str],
    host_deviation_summary: dict[str, object] | None,
) -> None:
    """Atomically write exact Phase-1 resume state as canonical JSON."""

    document = phase1_state_to_document(
        validated_markers=validated_markers,
        merged_seeds=merged_seeds,
        host_signature_model=host_signature_model,
        host_signatures=host_signatures,
        host_deviation_summary=host_deviation_summary,
    )
    try:
        content = json.dumps(
            document,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
    except (TypeError, ValueError) as exc:
        raise Phase1StateError(f"Phase-1 state is not JSON-safe: {exc}") from exc
    atomic_write(path, content + "\n")


def load_phase1_state(path: str | Path) -> Phase1ResumeState:
    """Load and strictly validate a Phase-1 resume-state JSON artifact."""

    state_path = Path(path)
    try:
        content = state_path.read_text(encoding="utf-8")
        document = json.loads(
            content,
            object_pairs_hook=_object_without_duplicate_keys,
            parse_constant=_reject_nonstandard_constant,
        )
    except Phase1StateError:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise Phase1StateError(f"invalid Phase-1 state {state_path}: {exc}") from exc
    return phase1_state_from_document(document)
