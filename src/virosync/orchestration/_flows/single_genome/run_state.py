"""Schema-v3 run identity, phase markers, and guarded resume state.

This module is deliberately leaf-like.  It knows how to authenticate files and
state, but it does not import the single-genome orchestrator or any phase code.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from contextvars import ContextVar
import csv
from dataclasses import dataclass, fields, is_dataclass
import fcntl
from functools import wraps
import hashlib
import io
import json
import math
import os
from pathlib import Path, PurePosixPath
import stat
import tempfile
from typing import Any
from urllib.parse import unquote

from virosync.ablation import (
    ABLATION_CONTRACT_SHA256,
    MAX_ABLATION_EVENTS_BYTES,
    AblationEvents,
    AblationID,
    ablation_policy,
    validate_ablation_events_bytes,
)


RUN_STATE_SCHEMA_VERSION = 3
RUN_STATE_FILENAME = "virosync_run_state.json"
PHASE_MARKER_FILENAMES = tuple(
    f"phase{phase}.complete.json" for phase in range(4)
)
RUN_STATUSES = frozenset({"running", "failed", "success"})
PHASE_OUTCOMES = frozenset(
    {"complete", "passthrough", "terminal_zero", "terminal_ablation"}
)
_TERMINAL_PHASE_OUTCOMES = frozenset({"terminal_zero", "terminal_ablation"})

_SHA256_LENGTH = 64
_MAX_STATE_BYTES = 16 * 1024 * 1024
# The largest measured Phase-1 state is 597,742,200 bytes; 1 GiB gives 1.80x room.
_MAX_CHECKPOINT_BYTES = 1024 * 1024 * 1024
_ARTIFACT_OBSERVATION_CACHE: ContextVar[
    dict[tuple[str, str, str], tuple[int, str, int | None]] | None
] = ContextVar("artifact_observation_cache", default=None)
_INPUT_SCAFFOLD_CACHE: ContextVar[
    dict[tuple[str, int, str], dict[str, int]] | None
] = ContextVar("input_scaffold_cache", default=None)
_RUNTIME_ENVIRONMENT_SHA256: str | None = None
_PHASE_DIRECTORIES = {
    0: ("phase0",),
    1: ("phase1",),
    2: ("phase2",),
    3: ("phase3", "phase3_synthesis", "structural_analysis"),
}
_FINAL_OUTPUTS = (
    "ablation_events.json",
    "run.log",
    "virosync_run_complete.json",
    "virosync_predictions.tsv",
    "virosync_predictions.bed",
    "virosync_predictions.gff3",
    "virosync_predictions_detailed.tsv",
    "virosync_summary.json",
    "virosync_tsv_invariant_report.tsv",
    "gvclass_results.tsv",
    "host_signature_model.png",
)


@contextmanager
def _validation_cache_scope() -> Iterator[None]:
    """Read each authenticated path once within one public decision.

    A cache hit is the decision's existing snapshot and does not re-stat the path.
    """

    artifact_token = _ARTIFACT_OBSERVATION_CACHE.set({})
    scaffold_token = _INPUT_SCAFFOLD_CACHE.set({})
    try:
        yield
    finally:
        _INPUT_SCAFFOLD_CACHE.reset(scaffold_token)
        _ARTIFACT_OBSERVATION_CACHE.reset(artifact_token)


def _scoped_validation_cache(function: Any) -> Any:
    @wraps(function)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        with _validation_cache_scope():
            return function(*args, **kwargs)

    return wrapper

# These are the minimum reloadable artifacts for a non-terminal phase.  A phase
# marker may record additional diagnostics, but it cannot make an empty or
# report-only directory authoritative.  Phase 1's JSON path is the lossless
# checkpoint consumed by schema-v3 resume; the public TSVs remain reports.
_PHASE_REQUIRED_PATHS = {
    0: frozenset(
        {
            "phase0/ablation_events.json",
            "phase0/masking/masking_status.json",
            "phase0/proteome.fasta",
            "phase0/genes.gff",
        }
    ),
    1: frozenset(
        {
            "phase1/ablation_events.json",
            "phase1/resume_state.json",
        }
    ),
    2: frozenset(
        {
            "phase2/ablation_events.json",
            "phase2/refined_boundaries.bed",
            "phase2/refined_state.json",
            "phase2/resume_state.json",
        }
    ),
    3: frozenset({"phase3/ablation_events.json"}),
}

# Full and terminal-zero runs place the canonical exports in different
# directories.  Each tuple is one required logical artifact with allowed paths.
_FINAL_REQUIRED_PATH_GROUPS = (
    ("ablation_events.json",),
    (
        "virosync_predictions.tsv",
        "phase3_synthesis/virosync_predictions.tsv",
    ),
    ("virosync_predictions_detailed.tsv",),
    (
        "virosync_predictions.bed",
        "phase3_synthesis/virosync_predictions.bed",
    ),
    (
        "virosync_predictions.gff3",
        "phase3_synthesis/virosync_predictions.gff3",
    ),
    (
        "virosync_summary.json",
        "phase3_synthesis/virosync_summary.json",
    ),
    ("virosync_tsv_invariant_report.tsv",),
    ("run.log",),
    ("virosync_run_complete.json",),
    ("notebooks/jupyter/eve_analysis.ipynb",),
)

_KNOWN_ARTIFACT_SCHEMAS = {
    "ablation_events.json": "virosync.ablation_events/v1",
    "phase0/ablation_events.json": "virosync.ablation_events/v1",
    "phase1/ablation_events.json": "virosync.ablation_events/v1",
    "phase2/ablation_events.json": "virosync.ablation_events/v1",
    "phase3/ablation_events.json": "virosync.ablation_events/v1",
    "phase0/masking/masking_status.json": "masking-status-v1",
    "phase1/resume_state.json": "virosync.phase1.resume_state/v1",
    "phase1/frameshift_screening/frameshift_hits.tsv": "frameshift-hits-v1",
    "phase1/frameshift_screening/confirmed_frameshift_proteins.faa": (
        "frameshift-rescued-proteins-v1"
    ),
    "phase1/frameshift_screening/confirmed_frameshift_markers.tsv": (
        "frameshift-rescued-markers-v1"
    ),
    "phase1/pfam_arbitration.tsv": "pfam-arbitration-v1",
    "phase2/refined_state.json": "virosync.phase2.refined_boundaries/v2",
    "phase2/resume_state.json": "virosync.phase2.resume_state/v1",
    "virosync_predictions.tsv": "canonical-predictions-v6",
    "phase3_synthesis/virosync_predictions.tsv": "canonical-predictions-v6",
    "virosync_predictions_detailed.tsv": "detailed-predictions-v6",
    "phase3_synthesis/virosync_predictions_detailed.tsv": (
        "detailed-predictions-v6"
    ),
    "virosync_predictions.bed": "canonical-predictions-bed-v1",
    "phase3_synthesis/virosync_predictions.bed": (
        "canonical-predictions-bed-v1"
    ),
    "virosync_predictions.gff3": "canonical-predictions-gff3-v1",
    "phase3_synthesis/virosync_predictions.gff3": (
        "canonical-predictions-gff3-v1"
    ),
    "virosync_summary.json": "virosync-summary-v3",
    "phase3_synthesis/virosync_summary.json": "virosync-summary-v3",
    # Pinned but not required: a genome with fewer than two accepted EVEs, or one
    # whose EVEs skani cannot sketch, legitimately publishes a header-only edge
    # table, and a run that never reaches Phase 3 publishes none at all.
    "phase3_synthesis/eve_ani_edges.tsv": "eve-ani-edges-v1",
    "virosync_tsv_invariant_report.tsv": "tsv-invariant-report-v1",
    "run.log": "run-log-v1",
    "virosync_run_complete.json": "completion-manifest-v2",
    "notebooks/jupyter/eve_analysis.ipynb": "eve-analysis-notebook-v1",
}


@dataclass(frozen=True)
class InputIdentity:
    """Whole-file identity for the submitted genome."""

    size: int
    sha256: str


@dataclass(frozen=True)
class ConfigIdentity:
    """Identity of the effective output-determining configuration."""

    sha256: str
    ablation_id: str
    ablation_contract_sha256: str


@dataclass(frozen=True)
class CodeIdentity:
    """Installed ViroSync version and source-tree identity."""

    version: str
    source_sha256: str


@dataclass(frozen=True)
class EnvironmentIdentity:
    """Identity of the output-determining execution environment."""

    lock_sha256: str
    runtime_sha256: str
    requested_device: str
    effective_device: str
    sha256: str


@dataclass(frozen=True)
class ResourceIdentity:
    """Manifest-bound identity for one enabled resource."""

    name: str
    kind: str
    version: str
    manifest_sha256: str


@dataclass(frozen=True)
class ArtifactIdentity:
    """Exact identity of one required phase or final artifact."""

    relative_path: str
    size: int
    sha256: str
    schema: str
    row_count: int | None


@dataclass(frozen=True)
class PhaseRecord:
    """Authenticated completion record for one pipeline phase."""

    schema_version: int
    phase: int
    run_fingerprint: str
    dependency_sha256: str
    artifacts: tuple[ArtifactIdentity, ...]
    outcome: str
    requested_masking: dict[str, object] | None = None
    actual_masking: dict[str, object] | None = None


@dataclass(frozen=True)
class ResumePlan:
    """Sequentially authenticated phase prefix available to a retry."""

    reusable_phases: tuple[int, ...]
    restart_phase: int
    terminal_phase: int | None = None
    completed: bool = False
    reason: str | None = None


@dataclass(frozen=True)
class RunState:
    """Authoritative state for one run fingerprint and attempt."""

    schema_version: int
    run_fingerprint: str
    status: str
    attempt: int
    identities: dict[str, object]
    result: dict[str, object] | None = None
    failure: dict[str, str] | None = None
    artifacts: tuple[ArtifactIdentity, ...] = ()


def _jsonable(value: object) -> object:
    if is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: _jsonable(getattr(value, field.name))
            for field in fields(value)
        }
    if isinstance(value, Mapping):
        result: dict[str, object] = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise TypeError("canonical JSON object keys must be strings")
            result[key] = _jsonable(item)
        return result
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("canonical JSON does not allow NaN or infinity")
        return value
    raise TypeError(f"value is not representable as canonical JSON: {type(value)!r}")


def canonical_json_bytes(value: object) -> bytes:
    """Encode *value* as deterministic UTF-8 JSON without insignificant bytes."""

    return json.dumps(
        _jsonable(value),
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def canonical_sha256(value: object) -> str:
    """Return the SHA-256 of :func:`canonical_json_bytes`."""

    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def _require_sha256(value: object, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != _SHA256_LENGTH
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{label} must be a lowercase 64-character SHA-256")
    return value


def _require_nonnegative_int(value: object, label: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{label} must be a non-negative integer")
    return value


def _object_without_duplicate_keys(
    pairs: list[tuple[str, object]],
) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def _normalized_relative_path(value: object, label: str = "relative path") -> str:
    if not isinstance(value, str) or not value or "\\" in value:
        raise ValueError(f"{label} must be a normalized POSIX relative path")
    relative = PurePosixPath(value)
    if (
        relative.is_absolute()
        or str(relative) != value
        or any(part in {"", ".", ".."} for part in relative.parts)
    ):
        raise ValueError(f"{label} must be a normalized POSIX relative path")
    return value


def _require_directory_no_follow(path: Path, label: str) -> Path:
    candidate = Path(path)
    try:
        metadata = candidate.lstat()
    except OSError as exc:
        raise ValueError(f"{label} is not accessible: {candidate}") from exc
    if not stat.S_ISDIR(metadata.st_mode) or stat.S_ISLNK(metadata.st_mode):
        raise ValueError(f"{label} must be a real directory: {candidate}")
    return candidate


def _open_regular_no_follow(path: Path) -> tuple[int, os.stat_result]:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags)
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise ValueError(f"artifact must be a regular file: {path}")
    except Exception:
        os.close(descriptor)
        raise
    return descriptor, metadata


_DIRECTORY_OPEN_FLAGS = (
    os.O_RDONLY
    | getattr(os, "O_CLOEXEC", 0)
    | getattr(os, "O_DIRECTORY", 0)
    | getattr(os, "O_NOFOLLOW", 0)
)


def _open_artifact_no_follow(
    root: Path,
    relative_path: str,
) -> tuple[int, os.stat_result]:
    relative = _normalized_relative_path(relative_path, "artifact relative_path")
    descriptors: list[int] = []
    try:
        current_fd = os.open(root, _DIRECTORY_OPEN_FLAGS)
        descriptors.append(current_fd)
        for component in PurePosixPath(relative).parts[:-1]:
            current_fd = os.open(
                component,
                _DIRECTORY_OPEN_FLAGS,
                dir_fd=current_fd,
            )
            descriptors.append(current_fd)
        leaf_flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
        leaf_flags |= getattr(os, "O_NOFOLLOW", 0)
        descriptor = os.open(
            PurePosixPath(relative).parts[-1],
            leaf_flags,
            dir_fd=current_fd,
        )
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            os.close(descriptor)
            raise ValueError(f"artifact must be a regular file: {relative}")
        return descriptor, metadata
    finally:
        for directory_fd in reversed(descriptors):
            os.close(directory_fd)




def _sha256_regular_file(path: Path) -> tuple[int, str]:
    descriptor, metadata = _open_regular_no_follow(path)
    digest = hashlib.sha256()
    try:
        with os.fdopen(descriptor, "rb", closefd=True) as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
            after = os.fstat(handle.fileno())
    except Exception:
        try:
            os.close(descriptor)
        except OSError:
            pass
        raise
    if (
        (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
        != (metadata.st_dev, metadata.st_ino, metadata.st_size, metadata.st_mtime_ns)
    ):
        raise ValueError(f"file changed while hashing: {path}")
    return metadata.st_size, digest.hexdigest()


def _observe_artifact(
    root: Path,
    relative_path: str,
    schema: str,
) -> tuple[int, str, int | None]:
    cache_key = (str(root.absolute()), relative_path, schema)
    cache = _ARTIFACT_OBSERVATION_CACHE.get()
    if cache is not None and cache_key in cache:
        return cache[cache_key]
    descriptor, metadata = _open_artifact_no_follow(root, relative_path)
    signature = (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
    )
    digest = hashlib.sha256()
    suffix = PurePosixPath(relative_path).suffix.lower()
    lowered = schema.lower()
    if suffix == ".tsv" or "tsv" in lowered or "table" in lowered:
        row_mode = "table"
    elif suffix == ".csv" or "csv" in lowered:
        row_mode = "table"
    elif suffix in {".faa", ".fna", ".fa", ".fasta"} or "fasta" in lowered:
        row_mode = "fasta"
    elif suffix in {".bed", ".gff", ".gff3"} or any(
        token in lowered for token in ("bed", "gff")
    ):
        row_mode = "records"
    else:
        row_mode = None
    try:
        with os.fdopen(descriptor, "rb", closefd=True) as handle:
            if row_mode is None:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
                row_count = None
            else:
                nonempty_rows = 0
                fasta_rows = 0
                record_rows = 0
                for raw_line in handle:
                    digest.update(raw_line)
                    line = raw_line.decode("utf-8")
                    if line.strip():
                        nonempty_rows += 1
                    if raw_line.startswith(b">"):
                        fasta_rows += 1
                    if line.strip() and not line.lstrip().startswith("#"):
                        record_rows += 1
                if row_mode == "table":
                    if nonempty_rows == 0:
                        if schema in {
                            "canonical-predictions-v6",
                            "detailed-predictions-v6",
                        }:
                            raise ValueError(
                                "final prediction table has no header: "
                                f"{relative_path}"
                            )
                        row_count = 0
                    else:
                        row_count = nonempty_rows - 1
                elif row_mode == "fasta":
                    row_count = fasta_rows
                else:
                    row_count = record_rows
            after = os.fstat(handle.fileno())
    except Exception:
        try:
            os.close(descriptor)
        except OSError:
            pass
        raise
    after_identity = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    )
    if after_identity != signature:
        raise ValueError(f"artifact changed while hashing: {relative_path}")
    observation = (metadata.st_size, digest.hexdigest(), row_count)
    if cache is not None:
        cache[cache_key] = observation
    return observation


def _read_artifact_json(root: Path, relative_path: str) -> dict[str, object]:
    descriptor, metadata = _open_artifact_no_follow(root, relative_path)
    if metadata.st_size <= 0 or metadata.st_size > _MAX_CHECKPOINT_BYTES:
        os.close(descriptor)
        raise ValueError(
            f"artifact JSON has invalid size: {relative_path} "
            f"(observed={metadata.st_size}, allowed=1-{_MAX_CHECKPOINT_BYTES})"
        )
    with os.fdopen(descriptor, "rb", closefd=True) as handle:
        content = handle.read()
        after = os.fstat(handle.fileno())
    if (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
    ) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    ):
        raise ValueError(f"artifact changed while reading: {relative_path}")
    try:
        payload = json.loads(
            content.decode("utf-8"),
            object_pairs_hook=_object_without_duplicate_keys,
        )
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid artifact JSON: {relative_path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"artifact JSON must contain an object: {relative_path}")
    return payload


def build_artifact_identity(
    path: str | Path,
    *,
    root: str | Path,
    schema: str,
    row_count: int | None = None,
) -> ArtifactIdentity:
    """Hash one normalized, regular, non-symlink artifact below *root*."""

    if not isinstance(schema, str) or not schema.strip():
        raise ValueError("artifact schema must be a non-empty string")
    root_path = _require_directory_no_follow(Path(root), "artifact root")
    supplied = Path(path)
    if "\\" in os.fspath(path) or any(part == ".." for part in supplied.parts):
        raise ValueError("artifact path must be normalized")
    if supplied.is_absolute():
        try:
            relative_path = supplied.relative_to(root_path).as_posix()
        except ValueError as exc:
            raise ValueError(f"artifact escapes root: {supplied}") from exc
    else:
        relative_path = supplied.as_posix()
    relative_path = _normalized_relative_path(relative_path, "artifact path")
    size, digest, observed_rows = _observe_artifact(
        root_path,
        relative_path,
        schema,
    )
    if row_count is not None:
        expected_rows = _require_nonnegative_int(row_count, "artifact row_count")
        if observed_rows is None:
            raise ValueError(
                f"artifact row_count is not defined for schema {schema!r}"
            )
        if observed_rows != expected_rows:
            raise ValueError(
                f"artifact row_count mismatch for {relative_path}: "
                f"{observed_rows} != {expected_rows}"
            )
    return ArtifactIdentity(
        relative_path=relative_path,
        size=size,
        sha256=digest,
        schema=schema,
        row_count=observed_rows,
    )


def _coerce_artifact(value: ArtifactIdentity | Mapping[str, object]) -> ArtifactIdentity:
    if isinstance(value, ArtifactIdentity):
        artifact = value
    elif isinstance(value, Mapping):
        required = {"relative_path", "size", "sha256", "schema", "row_count"}
        if set(value) != required:
            raise ValueError("artifact identity fields differ from schema v3")
        artifact = ArtifactIdentity(
            relative_path=value["relative_path"],  # type: ignore[arg-type]
            size=value["size"],  # type: ignore[arg-type]
            sha256=value["sha256"],  # type: ignore[arg-type]
            schema=value["schema"],  # type: ignore[arg-type]
            row_count=value["row_count"],  # type: ignore[arg-type]
        )
    else:
        raise TypeError("artifact identity must be an ArtifactIdentity or mapping")
    _normalized_relative_path(artifact.relative_path, "artifact relative_path")
    _require_nonnegative_int(artifact.size, "artifact size")
    _require_sha256(artifact.sha256, "artifact sha256")
    if not isinstance(artifact.schema, str) or not artifact.schema.strip():
        raise ValueError("artifact schema must be a non-empty string")
    if artifact.row_count is not None:
        _require_nonnegative_int(artifact.row_count, "artifact row_count")
    return artifact


def _artifacts_by_path(
    artifacts: Sequence[ArtifactIdentity | Mapping[str, object]],
) -> dict[str, ArtifactIdentity]:
    normalized = [_coerce_artifact(artifact) for artifact in artifacts]
    by_path = {artifact.relative_path: artifact for artifact in normalized}
    if len(by_path) != len(normalized):
        raise ValueError("artifact set contains duplicate relative paths")
    return by_path


def _require_known_artifact_schemas(
    artifacts: Mapping[str, ArtifactIdentity],
) -> None:
    for relative_path, expected_schema in _KNOWN_ARTIFACT_SCHEMAS.items():
        artifact = artifacts.get(relative_path)
        if artifact is not None and artifact.schema != expected_schema:
            raise ValueError(
                f"artifact {relative_path} must use schema {expected_schema!r}"
            )


def _require_final_artifact_set(
    artifacts: Sequence[ArtifactIdentity | Mapping[str, object]],
) -> dict[str, ArtifactIdentity]:
    by_path = _artifacts_by_path(artifacts)
    _require_known_artifact_schemas(by_path)
    missing = [
        choices
        for choices in _FINAL_REQUIRED_PATH_GROUPS
        if not any(path in by_path for path in choices)
    ]
    if missing:
        rendered = [" or ".join(group) for group in missing]
        raise ValueError(f"final artifact set is incomplete: {rendered!r}")
    return by_path


def _require_phase_artifact_set(record: PhaseRecord) -> None:
    by_path = _artifacts_by_path(record.artifacts)
    _require_known_artifact_schemas(by_path)
    if record.phase == 0:
        if record.requested_masking is None or record.actual_masking is None:
            raise ValueError("Phase 0 requires requested and actual masking state")
        status = by_path.get("phase0/masking/masking_status.json")
        if status is not None and record.actual_masking.get("status_sha256") != status.sha256:
            raise ValueError("Phase 0 masking state does not bind masking_status.json")
    required = set(_PHASE_REQUIRED_PATHS[record.phase])
    if record.outcome in _TERMINAL_PHASE_OUTCOMES or record.phase == 3:
        _require_final_artifact_set(record.artifacts)
    if record.outcome == "terminal_zero" and record.phase not in {0, 3}:
        # Early exits need not have produced the normal phase checkpoint, but
        # their final zero-result contract must be complete.
        required.clear()
    missing = sorted(required - set(by_path))
    if missing:
        raise ValueError(
            f"Phase {record.phase} artifact set is incomplete: {missing!r}"
        )


def _validate_phase0_binding(
    root: Path,
    record: PhaseRecord,
    identities: Mapping[str, object],
) -> None:
    if record.phase != 0:
        return
    if record.requested_masking != identities.get("requested_masking"):
        raise ValueError("Phase 0 requested masking differs from the run identity")
    status_artifact = _artifacts_by_path(record.artifacts)[
        "phase0/masking/masking_status.json"
    ]
    status_payload = _read_artifact_json(root, status_artifact.relative_path)
    required_status_fields = {
        "schema_version",
        "status",
        "requested_backend",
        "effective_backend",
        "failure_policy",
        "input_sha256",
        "output_path",
        "output_sha256",
        "benchmark_eligible",
        "result_fingerprint",
    }
    if (
        status_payload.get("schema_version") != 1
        or not required_status_fields.issubset(status_payload)
        or not isinstance(status_payload.get("benchmark_eligible"), bool)
    ):
        raise ValueError("masking status is not canonical schema v1")
    fingerprint_payload = dict(status_payload)
    result_fingerprint = fingerprint_payload.pop("result_fingerprint")
    masking_fingerprint = hashlib.sha256(
        json.dumps(
            fingerprint_payload,
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()
    if result_fingerprint != masking_fingerprint:
        raise ValueError("masking status result fingerprint is invalid")
    from virosync.pipeline.phase0.masking import (
        masking_result_from_status_payload,
    )

    masking_result_from_status_payload(
        status_payload,
        status_path=root / status_artifact.relative_path,
        status_sha256=status_artifact.sha256,
    )
    if (
        status_payload.get("requested_backend")
        != record.requested_masking.get("backend")
        or status_payload.get("failure_policy")
        != record.requested_masking.get("failure_policy")
    ):
        raise ValueError("masking status differs from the requested masking identity")
    expected_actual = {**status_payload, "status_sha256": status_artifact.sha256}
    if record.actual_masking != expected_actual:
        raise ValueError("Phase 0 actual masking differs from masking_status.json")
    input_identity = identities.get("input")
    if not isinstance(input_identity, Mapping):
        raise ValueError("run identity has no input binding")
    if status_payload.get("input_sha256") != input_identity.get("sha256"):
        raise ValueError("Phase 0 consumed input differs from the run identity")
    output_sha256 = _require_sha256(
        status_payload.get("output_sha256"),
        "masking output_sha256",
    )
    output_path_raw = status_payload.get("output_path")
    if not isinstance(output_path_raw, str) or not output_path_raw:
        raise ValueError("masking status has no output path")
    output_path = Path(output_path_raw)
    if status_payload.get("effective_backend") == "off":
        if output_sha256 != input_identity.get("sha256"):
            raise ValueError("unmasked Phase 0 output differs from the input identity")
        size, observed_sha256 = _sha256_regular_file(output_path)
        if (
            size != input_identity.get("size")
            or observed_sha256 != output_sha256
        ):
            raise ValueError("unmasked Phase 0 output is missing or stale")
        return
    try:
        relative_output = output_path.resolve(strict=True).relative_to(
            root.resolve(strict=True)
        ).as_posix()
    except (OSError, ValueError) as exc:
        raise ValueError("masked Phase 0 output is outside the run directory") from exc
    output_artifact = _artifacts_by_path(record.artifacts).get(relative_output)
    if output_artifact is None or output_artifact.sha256 != output_sha256:
        raise ValueError("Phase 0 marker does not authenticate the masked output")


def _validate_phase_checkpoint(
    root: Path,
    record: PhaseRecord,
    identities: Mapping[str, object] | None = None,
) -> None:
    if identities is None:
        raise ValueError("phase validation requires the run identity")
    phase_events = _validate_ablation_event_binding(
        root,
        f"phase{record.phase}/ablation_events.json",
        identities,
    )
    owner_phase = ablation_policy(phase_events.ablation_id).intervention_phase
    if (
        owner_phase is not None
        and record.phase < owner_phase
        and phase_events.counters.total_opportunities != 0
    ):
        raise ValueError("ablation counters are nonzero before their owner phase")
    if record.phase > 0:
        previous_events = _validate_ablation_event_binding(
            root,
            f"phase{record.phase - 1}/ablation_events.json",
            identities,
        )
        previous_counts = previous_events.counters.to_document()
        current_counts = phase_events.counters.to_document()
        for key, fields in previous_counts.items():
            if any(
                current_counts[key][field] < value
                for field, value in fields.items()
            ):
                raise ValueError(
                    "ablation counters decreased across phase fragments"
                )
        if (
            owner_phase is not None
            and record.phase > owner_phase
            and phase_events.counters != previous_events.counters
        ):
            raise ValueError("ablation counters changed after their owner phase")
    if record.phase == 3 or record.outcome in _TERMINAL_PHASE_OUTCOMES:
        root_events = _validate_ablation_event_binding(
            root,
            "ablation_events.json",
            identities,
        )
        if root_events != phase_events:
            raise ValueError(
                "final ablation events differ from the terminal phase fragment"
            )
    if record.outcome not in {"complete", "passthrough", "terminal_ablation"}:
        return
    artifacts = _artifacts_by_path(record.artifacts)
    if record.phase == 1:
        from .phase1_state import phase1_state_from_document

        phase1_state_from_document(
            _read_artifact_json(root, "phase1/resume_state.json")
        )
    elif record.phase == 2:
        from .phase2_resume_state import phase2_resume_state_from_document
        from .phase_state import phase2_state_from_document, phase2_state_to_document

        boundaries = phase2_state_from_document(
            _read_artifact_json(root, "phase2/refined_state.json")
        )
        resume_state = phase2_resume_state_from_document(
            _read_artifact_json(root, "phase2/resume_state.json")
        )
        if not boundaries:
            raise ValueError("completed Phase 2 checkpoint contains no boundaries")
        if phase2_state_to_document(boundaries) != phase2_state_to_document(
            resume_state.refined_boundaries
        ):
            raise ValueError("Phase 2 boundary checkpoints disagree")
        if identities is None:
            raise ValueError("Phase 2 validation requires the run identity")
        scaffold_lengths = _authenticated_scaffold_lengths(identities)
        for boundary in boundaries:
            scaffold_length = scaffold_lengths.get(boundary.scaffold)
            if scaffold_length is None or boundary.end > scaffold_length:
                raise ValueError(
                    "Phase 2 boundary lies outside the authenticated input FASTA"
                )
        bed = artifacts["phase2/refined_boundaries.bed"]
        if bed.row_count != len(boundaries):
            raise ValueError("Phase 2 BED row count differs from its checkpoint")
        descriptor, metadata = _open_artifact_no_follow(
            root,
            "phase2/refined_boundaries.bed",
        )
        observed_boundaries: list[tuple[str, int, int, str, int, str]] = []
        with os.fdopen(
            descriptor,
            "r",
            encoding="utf-8",
            closefd=True,
        ) as handle:
            for row_number, line in enumerate(handle, start=1):
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 6:
                    raise ValueError(
                        f"Phase 2 BED row {row_number} does not contain six fields"
                    )
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                    score = int(fields[4])
                except ValueError as exc:
                    raise ValueError(
                        f"Phase 2 BED row {row_number} has invalid numerics"
                    ) from exc
                if start < 0 or end <= start or not 0 <= score <= 1000:
                    raise ValueError(
                        f"Phase 2 BED row {row_number} violates BED6 semantics"
                    )
                observed_boundaries.append(
                    (fields[0], start, end, fields[3], score, fields[5])
                )
            after = os.fstat(handle.fileno())
        if (
            metadata.st_dev,
            metadata.st_ino,
            metadata.st_size,
            metadata.st_mtime_ns,
            metadata.st_ctime_ns,
        ) != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
            after.st_ctime_ns,
        ):
            raise ValueError("Phase 2 BED changed while validating")
        expected_boundaries = [
            (
                boundary.scaffold,
                boundary.start,
                boundary.end,
                f"EVE_{boundary.scaffold}_{boundary.start}-{boundary.end}",
                int(boundary.confidence * 1000),
                ".",
            )
            for boundary in boundaries
        ]
        if observed_boundaries != expected_boundaries:
            raise ValueError("Phase 2 BED differs from its lossless checkpoint")


def validate_artifact_identity(
    artifact: ArtifactIdentity | Mapping[str, object],
    *,
    root: str | Path,
) -> bool:
    """Return whether an artifact still has its recorded identity and row count."""

    try:
        expected = _coerce_artifact(artifact)
        size, digest, observed_rows = _observe_artifact(
            Path(root),
            expected.relative_path,
            expected.schema,
        )
        if size != expected.size or digest != expected.sha256:
            return False
        if observed_rows is not None or expected.row_count is not None:
            return observed_rows == expected.row_count
        return True
    except (OSError, UnicodeError, ValueError, TypeError, csv.Error):
        return False


def build_input_identity(path: str | Path) -> InputIdentity:
    """Build the whole-file input identity used by the run fingerprint."""

    candidate = Path(path)
    try:
        metadata = candidate.lstat()
    except OSError as exc:
        raise ValueError(f"input is not accessible: {candidate}") from exc
    if not stat.S_ISREG(metadata.st_mode) or stat.S_ISLNK(metadata.st_mode):
        raise ValueError(f"input must be a regular non-symlink file: {candidate}")
    size, digest = _sha256_regular_file(candidate)
    return InputIdentity(size=size, sha256=digest)


def _authenticated_scaffold_lengths(
    identities: Mapping[str, object],
) -> dict[str, int]:
    """Read the fingerprinted FASTA once and return its sequence lengths."""

    input_path = identities.get("input_path")
    input_identity = identities.get("input")
    if not isinstance(input_path, str) or not isinstance(input_identity, Mapping):
        raise ValueError("run identity has no authenticated input FASTA")
    expected_size = _require_nonnegative_int(
        input_identity.get("size"),
        "input size",
    )
    expected_sha256 = _require_sha256(
        input_identity.get("sha256"),
        "input sha256",
    )
    cache_key = (input_path, expected_size, expected_sha256)
    cache = _INPUT_SCAFFOLD_CACHE.get()
    if cache is not None and cache_key in cache:
        return dict(cache[cache_key])
    descriptor, metadata = _open_regular_no_follow(Path(input_path))
    signature = (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
    )
    digest = hashlib.sha256()
    lengths: dict[str, int] = {}
    current: str | None = None
    try:
        with os.fdopen(descriptor, "rb", closefd=True) as handle:
            for raw_line in handle:
                digest.update(raw_line)
                if raw_line.startswith(b">"):
                    try:
                        header = raw_line[1:].decode("utf-8").strip()
                    except UnicodeDecodeError as exc:
                        raise ValueError("input FASTA header is not UTF-8") from exc
                    current = header.split(maxsplit=1)[0] if header else ""
                    if not current or current in lengths:
                        raise ValueError(
                            "input FASTA contains an empty or duplicate scaffold ID"
                        )
                    lengths[current] = 0
                    continue
                sequence = b"".join(raw_line.split())
                if not sequence:
                    continue
                if current is None:
                    raise ValueError("input FASTA sequence precedes its header")
                lengths[current] += len(sequence)
            after = os.fstat(handle.fileno())
    except Exception:
        try:
            os.close(descriptor)
        except OSError:
            pass
        raise
    before_identity = signature
    after_identity = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    )
    if before_identity != after_identity:
        raise ValueError("input FASTA changed while validating coordinates")
    if not lengths:
        raise ValueError("input FASTA contains no scaffold records")
    if metadata.st_size != expected_size or digest.hexdigest() != expected_sha256:
        raise ValueError("input FASTA differs from the run identity")
    if cache is not None:
        cache[cache_key] = dict(lengths)
    return lengths


def _scan_regular_tree(
    root: Path,
    *,
    label: str,
    python_only: bool = False,
) -> tuple[tuple[str, ...], tuple[object, ...]]:
    """Return stable regular-file paths plus a metadata mutation signature."""

    root = _require_directory_no_follow(root, label)
    files: list[str] = []
    signature: list[object] = []

    def visit(directory: Path) -> None:
        with os.scandir(directory) as entries:
            ordered = sorted(entries, key=lambda entry: entry.name)
        for entry in ordered:
            metadata = entry.stat(follow_symlinks=False)
            path = Path(entry.path)
            relative = path.relative_to(root).as_posix()
            _normalized_relative_path(relative, f"{label} path")
            if stat.S_ISLNK(metadata.st_mode):
                raise ValueError(f"{label} contains a symlink: {relative}")
            if stat.S_ISDIR(metadata.st_mode):
                signature.append(
                    (
                        "directory",
                        relative,
                        metadata.st_dev,
                        metadata.st_ino,
                        metadata.st_mtime_ns,
                        metadata.st_ctime_ns,
                    )
                )
                visit(path)
                continue
            if not stat.S_ISREG(metadata.st_mode):
                raise ValueError(f"{label} contains a special file: {relative}")
            if python_only and path.suffix != ".py":
                continue
            files.append(relative)
            signature.append(
                (
                    "file",
                    relative,
                    metadata.st_dev,
                    metadata.st_ino,
                    metadata.st_size,
                    metadata.st_mtime_ns,
                    metadata.st_ctime_ns,
                )
            )

    visit(root)
    return tuple(files), tuple(signature)


def _hash_relative_file(root: Path, relative_path: str) -> tuple[int, str]:
    descriptor, metadata = _open_artifact_no_follow(root, relative_path)
    digest = hashlib.sha256()
    try:
        with os.fdopen(descriptor, "rb", closefd=True) as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
            after = os.fstat(handle.fileno())
    except Exception:
        try:
            os.close(descriptor)
        except OSError:
            pass
        raise
    before_identity = (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
    )
    after_identity = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    )
    if before_identity != after_identity:
        raise ValueError(f"file changed while hashing: {relative_path}")
    return metadata.st_size, digest.hexdigest()


def _synthetic_resource_inventory(root: Path) -> list[dict[str, object]]:
    root = _require_directory_no_follow(root, "resource root")
    relative_paths, signature = _scan_regular_tree(
        root,
        label="resource tree",
    )
    if not relative_paths:
        raise ValueError(f"resource tree contains no files: {root}")
    inventory = [
        {"path": relative, "size": size, "sha256": digest}
        for relative in relative_paths
        for size, digest in [_hash_relative_file(root, relative)]
    ]
    after_paths, after_signature = _scan_regular_tree(
        root,
        label="resource tree",
    )
    if after_paths != relative_paths or after_signature != signature:
        raise ValueError("resource tree changed while hashing")
    return sorted(inventory, key=lambda item: str(item["path"]))


def build_resource_identity(
    name: str,
    version: str,
    manifest: str | Path,
    *,
    kind: str,
) -> ResourceIdentity:
    """Build an R7-manifest identity or a synthetic whole-tree identity.

    A canonical R7 manifest is cheap to identify because its payload digest is
    already the authenticated bundle identity.  Optional resources without that
    contract receive a deterministic content manifest over every regular file.
    """

    if not isinstance(name, str) or not name.strip():
        raise ValueError("resource name must be a non-empty string")
    if not isinstance(version, str) or not version.strip():
        raise ValueError("resource version must be a non-empty string")
    if kind not in {"core", "optional"}:
        raise ValueError("resource kind must be 'core' or 'optional'")
    candidate = Path(manifest)
    manifest_path = candidate / "RESOURCE_MANIFEST.json" if candidate.is_dir() else candidate
    from virosync.utils.resource_manifest import (
        ResourceManifestError,
        load_resource_manifest,
    )

    if kind == "core":
        parsed = load_resource_manifest(
            manifest_path,
            expected_version=version,
        )
        return ResourceIdentity(
            name=name,
            kind=kind,
            version=parsed.version,
            manifest_sha256=parsed.manifest_sha256,
        )

    try:
        parsed = load_resource_manifest(
            manifest_path,
            expected_version=version,
        )
    except (OSError, ResourceManifestError, ValueError):
        if candidate.is_dir():
            inventory = _synthetic_resource_inventory(candidate)
        else:
            try:
                metadata = candidate.lstat()
            except OSError as exc:
                raise ValueError(f"resource file is not accessible: {candidate}") from exc
            if not stat.S_ISREG(metadata.st_mode) or stat.S_ISLNK(metadata.st_mode):
                raise ValueError(
                    f"resource file must be regular and non-symlink: {candidate}"
                )
            size, file_digest = _sha256_regular_file(candidate)
            inventory = [
                {
                    "path": candidate.name,
                    "size": size,
                    "sha256": file_digest,
                }
            ]
        digest = canonical_sha256(
            {
                "schema_version": 1,
                "resource_name": name,
                "resource_version": version,
                "files": inventory,
            }
        )
        return ResourceIdentity(
            name=name,
            kind=kind,
            version=version,
            manifest_sha256=digest,
        )
    return ResourceIdentity(
        name=name,
        kind=kind,
        version=parsed.version,
        manifest_sha256=parsed.manifest_sha256,
    )


def build_resource_set_identity(
    name: str,
    version: str,
    members: Mapping[str, str | Path],
    *,
    kind: str = "optional",
) -> ResourceIdentity:
    """Hash a declared set of regular resource files as one synthetic manifest."""

    if not isinstance(name, str) or not name.strip():
        raise ValueError("resource name must be a non-empty string")
    if not isinstance(version, str) or not version.strip():
        raise ValueError("resource version must be a non-empty string")
    if kind not in {"core", "optional"}:
        raise ValueError("resource kind must be 'core' or 'optional'")
    if not members:
        raise ValueError("resource dependency set must not be empty")
    inventory: list[dict[str, object]] = []
    for label, member in sorted(members.items()):
        normalized_label = _normalized_relative_path(label, "resource member label")
        candidate = Path(member)
        try:
            metadata = candidate.lstat()
        except OSError as exc:
            raise ValueError(f"resource member is not accessible: {candidate}") from exc
        if not stat.S_ISREG(metadata.st_mode) or stat.S_ISLNK(metadata.st_mode):
            raise ValueError(
                f"resource member must be regular and non-symlink: {candidate}"
            )
        size, digest = _sha256_regular_file(candidate)
        inventory.append(
            {
                "path": normalized_label,
                "size": size,
                "sha256": digest,
            }
        )
    return ResourceIdentity(
        name=name,
        kind=kind,
        version=version,
        manifest_sha256=canonical_sha256(
            {
                "schema_version": 1,
                "resource_name": name,
                "resource_version": version,
                "files": inventory,
            }
        ),
    )


def build_code_identity(
    source_root: str | Path,
    *,
    version: str,
) -> CodeIdentity:
    """Hash the installed Python source inventory without following symlinks."""

    if not isinstance(version, str) or not version:
        raise ValueError("code version must be a non-empty string")
    root = _require_directory_no_follow(Path(source_root), "source root")
    relative_paths, signature = _scan_regular_tree(
        root,
        label="source tree",
        python_only=True,
    )
    if not relative_paths:
        raise ValueError(f"source root contains no Python files: {root}")
    inventory = [
        {"path": relative, "size": size, "sha256": digest}
        for relative in relative_paths
        for size, digest in [_hash_relative_file(root, relative)]
    ]
    after_paths, after_signature = _scan_regular_tree(
        root,
        label="source tree",
        python_only=True,
    )
    if after_paths != relative_paths or after_signature != signature:
        raise ValueError("source tree changed while hashing")
    return CodeIdentity(
        version=version,
        source_sha256=canonical_sha256(
            {"schema_version": 1, "files": inventory}
        ),
    )


def build_environment_identity(
    lock_path: str | Path,
    *,
    requested_device: str,
    effective_device: str,
) -> EnvironmentIdentity:
    """Bind the lockfile plus requested and effective accelerator choices."""

    if not isinstance(requested_device, str) or not requested_device:
        raise ValueError("requested_device must be a non-empty string")
    if not isinstance(effective_device, str) or not effective_device:
        raise ValueError("effective_device must be a non-empty string")
    lock = build_input_identity(lock_path)
    payload = {
        "lock_sha256": lock.sha256,
        "runtime_sha256": runtime_environment_sha256(),
        "requested_device": requested_device,
        "effective_device": effective_device,
    }
    return EnvironmentIdentity(**payload, sha256=canonical_sha256(payload))


def runtime_environment_sha256() -> str:
    """Hash the effective interpreter and installed distribution inventory."""

    global _RUNTIME_ENVIRONMENT_SHA256
    if _RUNTIME_ENVIRONMENT_SHA256 is not None:
        return _RUNTIME_ENVIRONMENT_SHA256
    import platform
    import sys
    from importlib.metadata import distributions

    installed = sorted(
        (
            distribution.metadata.get("Name", "unknown"),
            distribution.version,
        )
        for distribution in distributions()
    )
    uname = os.uname()
    _RUNTIME_ENVIRONMENT_SHA256 = canonical_sha256(
        {
            "schema_version": 1,
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "platform": {
                "sys_platform": sys.platform,
                "machine": platform.machine(),
                "os": {
                    "sysname": uname.sysname,
                    "release": uname.release,
                    "version": uname.version,
                    "machine": uname.machine,
                },
            },
            "installed_distributions": installed,
        }
    )
    return _RUNTIME_ENVIRONMENT_SHA256


def _validate_run_identity(identity: Mapping[str, object]) -> dict[str, object]:
    document = _jsonable(identity)
    if not isinstance(document, dict):
        raise TypeError("run identity must be a JSON object")
    required = {
        "genome_id",
        "input_path",
        "output_dir",
        "input",
        "config",
        "code",
        "environment",
        "coordinate_schema_version",
        "coordinate_convention",
        "output_schema_version",
        "summary_schema_version",
        "requested_masking",
        "resources",
    }
    if set(document) != required:
        raise ValueError(
            "run identity fields differ from schema v3; "
            f"missing={sorted(required - set(document))}, "
            f"extra={sorted(set(document) - required)}"
        )

    input_identity = document["input"]
    config_identity = document["config"]
    code_identity = document["code"]
    environment_identity = document["environment"]
    if not isinstance(document["genome_id"], str) or not document["genome_id"]:
        raise ValueError("run identity genome_id must be a non-empty string")
    for field in ("input_path", "output_dir"):
        value = document[field]
        if not isinstance(value, str) or not Path(value).is_absolute():
            raise ValueError(f"run identity {field} must be an absolute path")
    if not isinstance(input_identity, dict):
        raise ValueError("run identity input must be an object")
    if not isinstance(config_identity, dict):
        raise ValueError("run identity config must be an object")
    if not isinstance(code_identity, dict):
        raise ValueError("run identity code must be an object")
    if not isinstance(environment_identity, dict):
        raise ValueError("run identity environment must be an object")
    nested_fields = (
        (input_identity, {"size", "sha256"}, "input"),
        (
            config_identity,
            {"sha256", "ablation_id", "ablation_contract_sha256"},
            "config",
        ),
        (code_identity, {"version", "source_sha256"}, "code"),
        (
            environment_identity,
            {
                "lock_sha256",
                "runtime_sha256",
                "requested_device",
                "effective_device",
                "sha256",
            },
            "environment",
        ),
    )
    for nested, expected_fields, label in nested_fields:
        if set(nested) != expected_fields:
            raise ValueError(f"run identity {label} fields differ from schema v3")
    _require_nonnegative_int(input_identity.get("size"), "input size")
    _require_sha256(input_identity.get("sha256"), "input sha256")
    _require_sha256(config_identity.get("sha256"), "config sha256")
    try:
        AblationID(config_identity.get("ablation_id"))
    except (TypeError, ValueError) as exc:
        raise ValueError("config ablation_id must be one of A0-A6") from exc
    _require_sha256(
        config_identity.get("ablation_contract_sha256"),
        "config ablation_contract_sha256",
    )
    if config_identity["ablation_contract_sha256"] != ABLATION_CONTRACT_SHA256:
        raise ValueError("config ablation contract does not match this ViroSync build")
    if not isinstance(code_identity.get("version"), str) or not code_identity.get(
        "version"
    ):
        raise ValueError("code version must be a non-empty string")
    _require_sha256(code_identity.get("source_sha256"), "code source_sha256")
    _require_sha256(environment_identity.get("lock_sha256"), "environment lock_sha256")
    _require_sha256(
        environment_identity.get("runtime_sha256"),
        "environment runtime_sha256",
    )
    for field in ("requested_device", "effective_device"):
        if not isinstance(environment_identity.get(field), str) or not environment_identity[
            field
        ]:
            raise ValueError(f"environment {field} must be a non-empty string")
    _require_sha256(environment_identity.get("sha256"), "environment sha256")
    environment_payload = {
        "lock_sha256": environment_identity["lock_sha256"],
        "runtime_sha256": environment_identity["runtime_sha256"],
        "requested_device": environment_identity["requested_device"],
        "effective_device": environment_identity["effective_device"],
    }
    if canonical_sha256(environment_payload) != environment_identity["sha256"]:
        raise ValueError("environment sha256 does not match its identity fields")
    for field in (
        "coordinate_schema_version",
        "output_schema_version",
        "summary_schema_version",
    ):
        value = document[field]
        if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
            raise ValueError(f"{field} must be a positive integer")
    if document["coordinate_convention"] != "0-based, half-open [start, end)":
        raise ValueError("coordinate_convention is not the canonical interval contract")
    requested_masking = document["requested_masking"]
    if not isinstance(requested_masking, dict):
        raise ValueError("requested_masking must be an object")
    masking_fields = {
        "backend",
        "failure_policy",
        "fallback_backend",
        "repeatmasker_species",
        "repeatmasker_library",
        "repeatmasker_library_sha256",
    }
    if set(requested_masking) != masking_fields:
        raise ValueError("requested_masking fields differ from schema v3")
    resources = document["resources"]
    if not isinstance(resources, list):
        raise ValueError("resources must be a list")
    normalized_resources: list[dict[str, object]] = []
    for resource in resources:
        if not isinstance(resource, dict):
            raise ValueError("each resource identity must be an object")
        if set(resource) != {"name", "kind", "version", "manifest_sha256"}:
            raise ValueError("resource identity fields differ from schema v3")
        for field in ("name", "kind", "version"):
            if not isinstance(resource.get(field), str) or not resource[field]:
                raise ValueError(f"resource {field} must be a non-empty string")
        if resource["kind"] not in {"core", "optional"}:
            raise ValueError("resource kind must be 'core' or 'optional'")
        _require_sha256(resource.get("manifest_sha256"), "resource manifest_sha256")
        normalized_resources.append(resource)
    document["resources"] = sorted(
        normalized_resources,
        key=lambda resource: (
            str(resource["name"]),
            str(resource["version"]),
            str(resource["manifest_sha256"]),
        ),
    )
    return document


def compute_run_fingerprint(payload: Mapping[str, object]) -> str:
    """Hash the complete canonical run identity as lowercase SHA-256."""

    return canonical_sha256(_validate_run_identity(payload))


def build_phase_record(
    *,
    phase: int,
    run_fingerprint: str,
    dependency_sha256: str,
    artifacts: Sequence[ArtifactIdentity | Mapping[str, object]],
    outcome: str,
    requested_masking: Mapping[str, object] | None = None,
    actual_masking: Mapping[str, object] | None = None,
) -> PhaseRecord:
    """Build a deterministic phase completion record."""

    if phase not in range(len(PHASE_MARKER_FILENAMES)):
        raise ValueError("phase must be in the range 0..3")
    _require_sha256(run_fingerprint, "run_fingerprint")
    _require_sha256(dependency_sha256, "dependency_sha256")
    if outcome not in PHASE_OUTCOMES:
        raise ValueError(f"invalid phase outcome: {outcome!r}")
    if outcome == "passthrough" and phase != 2:
        raise ValueError("passthrough is only valid for Phase 2")
    if outcome == "terminal_ablation" and phase != 1:
        raise ValueError("terminal_ablation is only valid for Phase 1")
    normalized_artifacts = tuple(
        sorted(
            (_coerce_artifact(artifact) for artifact in artifacts),
            key=lambda artifact: artifact.relative_path,
        )
    )
    paths = [artifact.relative_path for artifact in normalized_artifacts]
    if len(paths) != len(set(paths)):
        raise ValueError("phase record contains duplicate artifact paths")
    if (requested_masking is None) != (actual_masking is None):
        raise ValueError("requested and actual masking identities must be paired")
    if phase != 0 and (requested_masking is not None or actual_masking is not None):
        raise ValueError("masking identities belong only to Phase 0")
    requested = _jsonable(requested_masking) if requested_masking is not None else None
    actual = _jsonable(actual_masking) if actual_masking is not None else None
    if requested is not None and not isinstance(requested, dict):
        raise TypeError("requested_masking must be an object")
    if actual is not None and not isinstance(actual, dict):
        raise TypeError("actual_masking must be an object")
    return PhaseRecord(
        schema_version=RUN_STATE_SCHEMA_VERSION,
        phase=phase,
        run_fingerprint=run_fingerprint,
        dependency_sha256=dependency_sha256,
        artifacts=normalized_artifacts,
        outcome=outcome,
        requested_masking=requested,
        actual_masking=actual,
    )




def marker_sha256(path: str | Path, phase: int | None = None) -> str:
    """Hash the exact published marker bytes used by the next phase."""

    candidate = Path(path)
    if phase is not None:
        if phase not in range(len(PHASE_MARKER_FILENAMES)):
            raise ValueError("phase must be in the range 0..3")
        candidate = candidate / PHASE_MARKER_FILENAMES[phase]
    size, digest = _sha256_regular_file(candidate)
    if size <= 0:
        raise ValueError(f"phase marker is empty: {candidate}")
    return digest


def atomic_write_json(path: str | Path, payload: object) -> Path:
    """Durably publish complete JSON with replace and file/directory fsync."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    parent = _require_directory_no_follow(destination.parent, "JSON parent")
    try:
        existing = destination.lstat()
    except FileNotFoundError:
        existing = None
    if existing is not None and (
        not stat.S_ISREG(existing.st_mode) or stat.S_ISLNK(existing.st_mode)
    ):
        raise ValueError(f"JSON destination must be a regular file: {destination}")
    content = (
        json.dumps(
            _jsonable(payload),
            ensure_ascii=False,
            allow_nan=False,
            indent=2,
            sort_keys=True,
        )
        + "\n"
    ).encode("utf-8")
    if len(content) > _MAX_STATE_BYTES:
        raise ValueError(
            f"JSON state exceeds the {_MAX_STATE_BYTES}-byte reader limit"
        )
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb", closefd=True) as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
            os.fchmod(handle.fileno(), 0o644)
        os.replace(temporary, destination)
        directory_descriptor = os.open(
            parent,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_CLOEXEC", 0),
        )
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return destination


@contextmanager
def sibling_run_lock(output_dir: str | Path) -> Iterator[Path]:
    """Serialize attempts with a non-symlink lock beside the output directory."""

    output = Path(output_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    parent = _require_directory_no_follow(output.parent, "run-lock parent")
    lock_path = parent / f".{output.name}.virosync-run.lock"
    flags = os.O_RDWR | os.O_CREAT | getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(lock_path, flags, 0o600)
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
            raise ValueError(f"run lock must be a single-link regular file: {lock_path}")
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        current = lock_path.lstat()
        if (
            not stat.S_ISREG(current.st_mode)
            or current.st_nlink != 1
            or (current.st_dev, current.st_ino)
            != (metadata.st_dev, metadata.st_ino)
        ):
            raise ValueError(f"run lock changed while acquiring it: {lock_path}")
        yield lock_path
    finally:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            os.close(descriptor)


def _read_json(path: Path) -> dict[str, object]:
    descriptor, metadata = _open_regular_no_follow(path)
    if metadata.st_size <= 0 or metadata.st_size > _MAX_STATE_BYTES:
        os.close(descriptor)
        raise ValueError(f"state JSON has invalid size: {path}")
    with os.fdopen(descriptor, "rb", closefd=True) as handle:
        content = handle.read()
    try:
        payload = json.loads(
            content.decode("utf-8"),
            object_pairs_hook=_object_without_duplicate_keys,
        )
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid state JSON: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"state JSON must contain an object: {path}")
    return payload


def _phase_record_from_payload(payload: Mapping[str, object]) -> PhaseRecord:
    required = {
        "schema_version",
        "phase",
        "run_fingerprint",
        "dependency_sha256",
        "artifacts",
        "outcome",
        "requested_masking",
        "actual_masking",
    }
    if set(payload) != required or payload.get("schema_version") != RUN_STATE_SCHEMA_VERSION:
        raise ValueError("phase marker fields differ from schema v3")
    artifacts = payload["artifacts"]
    if not isinstance(artifacts, list):
        raise ValueError("phase artifacts must be a list")
    phase = payload["phase"]
    if not isinstance(phase, int) or isinstance(phase, bool):
        raise ValueError("phase must be an integer")
    requested = payload["requested_masking"]
    actual = payload["actual_masking"]
    if requested is not None and not isinstance(requested, dict):
        raise ValueError("requested_masking must be an object or null")
    if actual is not None and not isinstance(actual, dict):
        raise ValueError("actual_masking must be an object or null")
    record = build_phase_record(
        phase=phase,
        run_fingerprint=payload["run_fingerprint"],  # type: ignore[arg-type]
        dependency_sha256=payload["dependency_sha256"],  # type: ignore[arg-type]
        artifacts=artifacts,  # type: ignore[arg-type]
        outcome=payload["outcome"],  # type: ignore[arg-type]
        requested_masking=requested,
        actual_masking=actual,
    )
    if _jsonable(record) != dict(payload):
        raise ValueError("phase marker is not canonical schema-v3 state")
    return record


def _load_phase_record(path: Path) -> PhaseRecord:
    return _phase_record_from_payload(_read_json(path))




def _validated_phase_prefix(
    output_dir: Path,
    expected_run_fingerprint: str,
) -> tuple[tuple[PhaseRecord, ...], str | None]:
    try:
        state = _load_run_state(output_dir)
    except (FileNotFoundError, OSError, TypeError, ValueError) as exc:
        return (), f"missing or invalid run state: {exc}"
    if state.run_fingerprint != expected_run_fingerprint:
        return (), "run state has a stale fingerprint"
    records: list[PhaseRecord] = []
    dependency = expected_run_fingerprint
    for phase, filename in enumerate(PHASE_MARKER_FILENAMES):
        marker = output_dir / filename
        try:
            marker.lstat()
        except FileNotFoundError:
            return tuple(records), f"missing {filename}"
        try:
            record = _load_phase_record(marker)
        except (OSError, TypeError, ValueError) as exc:
            return tuple(records), f"invalid {filename}: {exc}"
        if record.phase != phase:
            return tuple(records), f"{filename} records phase {record.phase}"
        if record.run_fingerprint != expected_run_fingerprint:
            return tuple(records), f"{filename} has a stale run fingerprint"
        if record.dependency_sha256 != dependency:
            return tuple(records), f"{filename} has a stale phase dependency"
        try:
            _require_phase_artifact_set(record)
            _validate_phase0_binding(output_dir, record, state.identities)
            _validate_phase_checkpoint(output_dir, record, state.identities)
        except (OSError, TypeError, ValueError) as exc:
            return tuple(records), f"{filename} has an incomplete artifact set: {exc}"
        if any(
            not validate_artifact_identity(artifact, root=output_dir)
            for artifact in record.artifacts
        ):
            return tuple(records), f"{filename} has a stale artifact"
        records.append(record)
        dependency = marker_sha256(marker)
        if record.outcome in _TERMINAL_PHASE_OUTCOMES:
            if any(
                _entry_exists(output_dir / later)
                for later in PHASE_MARKER_FILENAMES[phase + 1:]
            ):
                return tuple(records[:-1]), f"{filename} precedes a downstream marker"
            return tuple(records), None
    return tuple(records), None


def _fsync_artifact(root: Path, artifact: ArtifactIdentity) -> None:
    descriptor, _metadata = _open_artifact_no_follow(
        root,
        artifact.relative_path,
    )
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)

    relative = _normalized_relative_path(
        artifact.relative_path,
        "artifact relative_path",
    )
    directory_descriptors: list[int] = []
    try:
        current_fd = os.open(root, _DIRECTORY_OPEN_FLAGS)
        directory_descriptors.append(current_fd)
        for component in PurePosixPath(relative).parts[:-1]:
            current_fd = os.open(
                component,
                _DIRECTORY_OPEN_FLAGS,
                dir_fd=current_fd,
            )
            directory_descriptors.append(current_fd)
        for directory_fd in reversed(directory_descriptors):
            os.fsync(directory_fd)
    finally:
        for directory_fd in reversed(directory_descriptors):
            os.close(directory_fd)


@_scoped_validation_cache
def publish_phase_completion(
    output_dir: str | Path,
    record: PhaseRecord | None = None,
    **record_fields: Any,
) -> Path:
    """Validate, fsync, and atomically publish one chained phase marker."""

    root = _require_directory_no_follow(Path(output_dir), "output directory")
    if record is None:
        record = build_phase_record(**record_fields)
    elif record_fields:
        raise TypeError("record fields cannot accompany an explicit PhaseRecord")
    if not isinstance(record, PhaseRecord):
        raise TypeError("record must be a PhaseRecord")
    _require_phase_artifact_set(record)
    current_state = _load_run_state(root)
    if current_state.run_fingerprint != record.run_fingerprint:
        raise ValueError("phase marker run fingerprint differs from run state")
    if current_state.status != "running":
        raise ValueError("phase markers can only be published for a running attempt")
    _validate_phase0_binding(root, record, current_state.identities)
    _validate_phase_checkpoint(root, record, current_state.identities)
    if any(
        not validate_artifact_identity(artifact, root=root)
        for artifact in record.artifacts
    ):
        raise ValueError(f"Phase {record.phase} artifacts failed identity validation")
    prefix, reason = _validated_phase_prefix(root, record.run_fingerprint)
    if len(prefix) < record.phase:
        raise ValueError(reason or "prior phase marker is invalid")
    if (
        prefix
        and prefix[-1].outcome in _TERMINAL_PHASE_OUTCOMES
        and record.phase > prefix[-1].phase
    ):
        raise ValueError("no phase may be published after a terminal completion")
    if record.phase == 0:
        expected_dependency = record.run_fingerprint
    else:
        expected_dependency = marker_sha256(
            root / PHASE_MARKER_FILENAMES[record.phase - 1]
        )
    if record.dependency_sha256 != expected_dependency:
        raise ValueError("phase dependency does not match the authenticated prefix")
    if any(
        _entry_exists(root / filename)
        for filename in PHASE_MARKER_FILENAMES[record.phase + 1:]
    ):
        raise ValueError("downstream phase markers must be invalidated before publish")
    for artifact in record.artifacts:
        _fsync_artifact(root, artifact)
    return atomic_write_json(root / PHASE_MARKER_FILENAMES[record.phase], record)


def _state_from_payload(payload: Mapping[str, object]) -> RunState:
    required = {
        "schema_version",
        "run_fingerprint",
        "status",
        "attempt",
        "identities",
        "result",
        "failure",
    }
    allowed = required | {"artifacts"}
    if (
        not required.issubset(payload)
        or not set(payload).issubset(allowed)
        or payload.get("schema_version") != RUN_STATE_SCHEMA_VERSION
    ):
        raise ValueError("run-state fields differ from schema v3")
    run_fingerprint = _require_sha256(payload["run_fingerprint"], "run_fingerprint")
    status = payload["status"]
    if status not in RUN_STATUSES:
        raise ValueError(f"invalid run status: {status!r}")
    attempt = payload["attempt"]
    if not isinstance(attempt, int) or isinstance(attempt, bool) or attempt <= 0:
        raise ValueError("run attempt must be a positive integer")
    identities = payload["identities"]
    if not isinstance(identities, dict):
        raise ValueError("run identities must be an object")
    identities = _validate_run_identity(identities)
    if run_fingerprint != compute_run_fingerprint(identities):
        raise ValueError("run fingerprint does not match the canonical run identity")
    artifacts = payload.get("artifacts", [])
    if not isinstance(artifacts, list):
        raise ValueError("run artifacts must be a list")
    normalized_artifacts = tuple(_coerce_artifact(item) for item in artifacts)
    result = payload["result"]
    if result is not None and not isinstance(result, dict):
        raise ValueError("run result must be an object or null")
    failure = payload["failure"]
    if failure is not None:
        if not isinstance(failure, dict) or set(failure) != {"error_type", "message"}:
            raise ValueError("run failure must contain error_type and message")
        if any(not isinstance(failure[field], str) for field in failure):
            raise ValueError("run failure fields must be strings")
    state = RunState(
        schema_version=RUN_STATE_SCHEMA_VERSION,
        run_fingerprint=run_fingerprint,
        status=status,
        attempt=attempt,
        identities=identities,
        result=result,
        failure=failure,  # type: ignore[arg-type]
        artifacts=normalized_artifacts,
    )
    if status in {"running", "failed"} and (state.artifacts or state.result is not None):
        raise ValueError(f"{status} state cannot carry success artifacts")
    if status == "running" and state.failure is not None:
        raise ValueError("running state cannot carry failure details")
    if status == "failed" and state.failure is None:
        raise ValueError("failed state must carry an error")
    if status == "success" and state.failure is not None:
        raise ValueError("success state cannot carry an error")
    if status == "success" and state.result is None:
        raise ValueError("success state must carry a result")
    return state


def _load_run_state(output_dir: Path) -> RunState:
    return _state_from_payload(_read_json(output_dir / RUN_STATE_FILENAME))


def load_run_state(output_dir: str | Path) -> RunState:
    """Load a structurally valid schema-v3 run state."""

    root = _require_directory_no_follow(Path(output_dir), "output directory")
    return _load_run_state(root)


def _unlink_no_follow(path: Path) -> None:
    try:
        metadata = path.lstat()
    except FileNotFoundError:
        return
    if stat.S_ISDIR(metadata.st_mode) and not stat.S_ISLNK(metadata.st_mode):
        raise ValueError(f"refusing to unlink a directory as a file: {path}")
    path.unlink()


def _entry_exists(path: Path) -> bool:
    """Return true for any directory entry, including dangling symlinks."""

    try:
        path.lstat()
    except FileNotFoundError:
        return False
    return True


@_scoped_validation_cache
def publish_run_started(
    output_dir: str | Path,
    run_fingerprint: str | None = None,
    *,
    identities: Mapping[str, object] | None = None,
    preserve_success_artifacts: bool = False,
) -> RunState:
    """Publish ``running`` before Phase 0 and advance the attempt counter."""

    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    root = _require_directory_no_follow(root, "output directory")
    previous: RunState | None
    try:
        previous = _load_run_state(root)
    except FileNotFoundError:
        previous = None
    except (OSError, TypeError, ValueError):
        _unlink_no_follow(root / RUN_STATE_FILENAME)
        previous = None

    normalized_identities: dict[str, object] | None = None
    if identities is not None:
        normalized_identities = _validate_run_identity(identities)
        identity_fingerprint = compute_run_fingerprint(normalized_identities)
        if run_fingerprint is None:
            run_fingerprint = identity_fingerprint
        elif run_fingerprint != identity_fingerprint:
            raise ValueError(
                "run fingerprint does not match the canonical run identity"
            )
    if run_fingerprint is None:
        raise ValueError("run_fingerprint or identity is required")
    _require_sha256(run_fingerprint, "run_fingerprint")

    if previous is not None and previous.run_fingerprint == run_fingerprint:
        if normalized_identities is None:
            normalized_identities = previous.identities
        elif canonical_json_bytes(normalized_identities) != canonical_json_bytes(
            previous.identities
        ):
            raise ValueError("same fingerprint has a different run identity")
        attempt = previous.attempt + 1
        if previous.status == "success":
            if preserve_success_artifacts:
                records, reason = _validated_phase_prefix(root, run_fingerprint)
                if not records or (
                    records[-1].phase != len(PHASE_MARKER_FILENAMES) - 1
                    and records[-1].outcome not in _TERMINAL_PHASE_OUTCOMES
                ):
                    raise ValueError(
                        reason
                        or "cannot preserve success artifacts without a final marker"
                    )
            else:
                invalidate_from_phase(root, from_phase=0)
    else:
        if preserve_success_artifacts:
            raise ValueError(
                "success artifacts can only be preserved for the same fingerprint"
            )
        if normalized_identities is None:
            raise ValueError("a new run fingerprint requires the complete identity")
        attempt = 1
        if previous is not None:
            _unlink_no_follow(root / RUN_STATE_FILENAME)
        invalidate_from_phase(root, from_phase=0)

    assert normalized_identities is not None
    state = RunState(
        schema_version=RUN_STATE_SCHEMA_VERSION,
        run_fingerprint=run_fingerprint,
        status="running",
        attempt=attempt,
        identities=normalized_identities,
    )
    atomic_write_json(root / RUN_STATE_FILENAME, state)
    return state


def publish_run_failed(
    output_dir: str | Path,
    *,
    run_fingerprint: str,
    error_type: str,
    message: str,
) -> RunState:
    """Atomically replace the current running state with ``failed``."""

    root = _require_directory_no_follow(Path(output_dir), "output directory")
    current = _load_run_state(root)
    if current.status != "running":
        raise ValueError("only a running attempt can transition to failed")
    if current.run_fingerprint != run_fingerprint:
        raise ValueError("failed state run fingerprint mismatch")
    if not isinstance(error_type, str) or not error_type:
        raise ValueError("error_type must be a non-empty string")
    if not isinstance(message, str):
        raise TypeError("message must be a string")
    failed = RunState(
        schema_version=RUN_STATE_SCHEMA_VERSION,
        run_fingerprint=current.run_fingerprint,
        status="failed",
        attempt=current.attempt,
        identities=current.identities,
        failure={"error_type": error_type, "message": message},
    )
    atomic_write_json(root / RUN_STATE_FILENAME, failed)
    return failed


def _validated_result(result: Mapping[str, object]) -> dict[str, object]:
    normalized = _jsonable(result)
    if not isinstance(normalized, dict):
        raise TypeError("success result must be an object")
    required = {
        "terminal_phase",
        "canonical_rows",
        "detailed_rows",
        "accepted_bp",
        "class_counts",
        "tier_counts",
        "benchmark_eligible",
    }
    missing = required - set(normalized)
    if missing:
        raise ValueError(f"success result is missing fields: {sorted(missing)!r}")
    canonical_rows = _require_nonnegative_int(
        normalized["canonical_rows"], "canonical_rows"
    )
    detailed_rows = _require_nonnegative_int(
        normalized["detailed_rows"], "detailed_rows"
    )
    _require_nonnegative_int(normalized["accepted_bp"], "accepted_bp")
    if detailed_rows < canonical_rows:
        raise ValueError("detailed_rows cannot be smaller than canonical_rows")
    for field in ("class_counts", "tier_counts"):
        counts = normalized[field]
        if not isinstance(counts, dict):
            raise ValueError(f"{field} must be an object")
        for name, value in counts.items():
            if not isinstance(name, str) or not name:
                raise ValueError(f"{field} keys must be non-empty strings")
            _require_nonnegative_int(value, f"{field}.{name}")
        if sum(counts.values()) != canonical_rows:  # type: ignore[arg-type]
            raise ValueError(f"{field} must sum to canonical_rows")
    from virosync.output_contract import normalize_effective_eve_class_counts

    normalized["class_counts"] = normalize_effective_eve_class_counts(
        normalized["class_counts"]
    )
    if set(normalized["tier_counts"]) != {"HIGH", "MEDIUM", "LOW"}:
        raise ValueError("tier_counts must contain HIGH, MEDIUM, and LOW")
    promoted_low_rows = normalized.get(
        "promoted_low_rows",
        normalized["tier_counts"]["LOW"],
    )
    promoted_low_rows = _require_nonnegative_int(
        promoted_low_rows,
        "promoted_low_rows",
    )
    if promoted_low_rows > normalized["tier_counts"]["LOW"]:
        raise ValueError("promoted_low_rows cannot exceed canonical LOW rows")
    normalized["promoted_low_rows"] = promoted_low_rows
    terminal_phase = normalized["terminal_phase"]
    if terminal_phase is not None and terminal_phase not in range(
        len(PHASE_MARKER_FILENAMES)
    ):
        raise ValueError("terminal_phase must be null or an integer in the range 0..3")
    if isinstance(terminal_phase, bool):
        raise ValueError("terminal_phase must not be boolean")
    if not isinstance(normalized["benchmark_eligible"], bool):
        raise ValueError("benchmark_eligible must be boolean")
    return normalized


def _validate_invariant_report(path: Path, *, expected_rows: int) -> None:
    descriptor, _metadata = _open_regular_no_follow(path)
    with os.fdopen(
        descriptor,
        "r",
        encoding="utf-8",
        newline="",
        closefd=True,
    ) as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    summary_header = [
        "status",
        "rows_checked",
        "issue_count",
        "error_count",
        "warning_count",
    ]
    if len(rows) < 2 or rows[0] != summary_header:
        raise ValueError("invariant report has an invalid summary schema")
    if len(rows[1]) != len(summary_header):
        raise ValueError("invariant report has an invalid summary row")
    status = rows[1][0]
    if status not in {"PASS", "PASS_WITH_WARNINGS"}:
        raise ValueError("invariant report did not pass")
    try:
        rows_checked = int(rows[1][1])
        issue_count = int(rows[1][2])
        error_count = int(rows[1][3])
        warning_count = int(rows[1][4])
    except ValueError as exc:
        raise ValueError("invariant report counts must be integers") from exc
    if (
        rows_checked != expected_rows
        or issue_count < 0
        or error_count != 0
        or warning_count != issue_count
    ):
        raise ValueError("invariant report counts disagree with final outputs")
    issue_rows = [row for row in rows[2:] if any(field.strip() for field in row)]
    if issue_rows:
        if issue_rows[0] != ["eve_id", "check", "severity", "message"]:
            raise ValueError("invariant report has an invalid issue schema")
        issue_rows = issue_rows[1:]
    if len(issue_rows) != issue_count:
        raise ValueError("invariant report issue_count disagrees with issue rows")
    if any(len(row) != 4 or row[2].strip().lower() != "warning" for row in issue_rows):
        raise ValueError("invariant report contains a fatal or malformed issue")
    if status == "PASS" and issue_count:
        raise ValueError("PASS invariant report cannot contain issues")
    if status == "PASS_WITH_WARNINGS" and not issue_count:
        raise ValueError("PASS_WITH_WARNINGS requires at least one warning")


def _validate_completion_manifest(
    path: Path,
    *,
    run_fingerprint: str,
    benchmark_eligible: bool,
    identities: Mapping[str, object],
) -> None:
    payload = _read_json(path)
    if payload.get("schema_version") != 2 or payload.get("status") != "success":
        raise ValueError("completion manifest is not a schema-v2 success")
    if payload.get("config_fingerprint") != run_fingerprint:
        raise ValueError("completion manifest has a stale run fingerprint")
    for field in (
        "genome_id",
        "coordinate_schema_version",
        "coordinate_convention",
        "output_schema_version",
    ):
        if payload.get(field) != identities.get(field):
            raise ValueError(
                f"completion manifest {field} differs from the run identity"
            )
    if not isinstance(payload.get("generated_at"), str) or not isinstance(
        payload.get("reason"),
        str,
    ):
        raise ValueError("completion manifest timestamps/reason are malformed")
    output_files = payload.get("output_files")
    if not isinstance(output_files, dict):
        raise ValueError("completion manifest has no output_files object")

    root = path.parent.absolute()

    def validate_output_paths(value: object) -> None:
        if isinstance(value, dict):
            for child in value.values():
                validate_output_paths(child)
            return
        if isinstance(value, (list, tuple)):
            for child in value:
                validate_output_paths(child)
            return
        if value is None:
            return
        if not isinstance(value, str) or not value:
            raise ValueError("completion manifest output paths are malformed")
        candidate = Path(value)
        if candidate.is_absolute():
            try:
                relative = candidate.absolute().relative_to(root).as_posix()
            except ValueError as exc:
                raise ValueError(
                    "completion manifest output path escapes the run directory"
                ) from exc
        else:
            relative = candidate.as_posix()
        _normalized_relative_path(relative, "completion output path")

    validate_output_paths(output_files)
    masking = payload.get("masking_status")
    if (
        not isinstance(masking, dict)
        or masking.get("benchmark_eligible") is not benchmark_eligible
    ):
        raise ValueError("completion manifest benchmark eligibility is inconsistent")
    status_path = path.parent / "phase0" / "masking" / "masking_status.json"
    from virosync.pipeline.phase0.masking import load_masking_result

    result = load_masking_result(status_path)
    status_payload = _read_json(status_path)
    expected_masking = {
        "path": "phase0/masking/masking_status.json",
        "sha256": result.status_sha256,
        "result_fingerprint": status_payload.get("result_fingerprint"),
        "benchmark_eligible": status_payload.get("benchmark_eligible"),
        "status": status_payload.get("status"),
    }
    if masking != expected_masking:
        raise ValueError("completion manifest masking identity is inconsistent")
    expected_effective = hashlib.sha256(
        f"{run_fingerprint}|{result.status_sha256}".encode()
    ).hexdigest()
    if payload.get("effective_masking_fingerprint") != expected_effective:
        raise ValueError("completion manifest effective masking identity is stale")


def _load_json_document(path: Path, label: str) -> dict[str, object]:
    try:
        payload = _read_json(path)
    except ValueError as exc:
        raise ValueError(f"{label} is not valid JSON") from exc
    if not payload:
        raise ValueError(f"{label} must not be empty")
    return payload


def _read_delimited_artifact(
    root: Path,
    relative_path: str,
    *,
    delimiter: str,
) -> tuple[list[str], list[dict[str, str]]]:
    descriptor, metadata = _open_artifact_no_follow(root, relative_path)
    with os.fdopen(descriptor, "rb", closefd=True) as binary:
        text = io.TextIOWrapper(binary, encoding="utf-8", newline="")
        try:
            reader = csv.DictReader(text, delimiter=delimiter)
            fields = list(reader.fieldnames or [])
            rows = [dict(row) for row in reader]
            after = os.fstat(binary.fileno())
        finally:
            text.detach()
    if (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
    ) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    ):
        raise ValueError(f"tabular artifact changed while reading: {relative_path}")
    if (
        not fields
        or any(not field.strip() for field in fields)
        or len(fields) != len(set(fields))
    ):
        raise ValueError(f"tabular artifact has an invalid header: {relative_path}")
    if any(None in row for row in rows):
        raise ValueError(f"tabular artifact has an extra-width row: {relative_path}")
    return fields, rows


def _prediction_coordinates(
    root: Path,
    relative_path: str,
) -> tuple[
    dict[str, tuple[str, int, int, int, str, float, str]],
    dict[str, object],
]:
    from virosync.output_contract import (
        EFFECTIVE_EVE_CLASSES,
        PPV_LEGACY_ALIASES,
        normalize_effective_eve_class,
    )

    fields, rows = _read_delimited_artifact(
        root,
        relative_path,
        delimiter="\t",
    )
    required = {
        "eve_id",
        "scaffold",
        "start",
        "end",
        "length",
        "confidence_tier",
        "final_confidence",
        "effective_eve_class",
    }
    if not required.issubset(fields):
        raise ValueError(
            f"prediction table {relative_path} lacks required schema columns"
        )
    coordinates: dict[
        str,
        tuple[str, int, int, int, str, float, str],
    ] = {}
    accepted_bp = 0
    class_counts = {eve_class: 0 for eve_class in EFFECTIVE_EVE_CLASSES}
    tier_counts = {tier: 0 for tier in ("HIGH", "MEDIUM", "LOW")}
    for row_number, row in enumerate(rows, start=2):
        eve_id = str(row.get("eve_id") or "").strip()
        scaffold = str(row.get("scaffold") or "").strip()
        if not eve_id or not scaffold or eve_id in coordinates:
            raise ValueError(
                f"prediction table {relative_path} has an invalid ID at row {row_number}"
            )
        try:
            start = int(row["start"])
            end = int(row["end"])
            length = int(row["length"])
            confidence = float(row["final_confidence"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"prediction table {relative_path} has invalid numeric fields "
                f"at row {row_number}"
            ) from exc
        if start < 0 or end <= start or length != end - start:
            raise ValueError(
                f"prediction table {relative_path} violates coordinates at "
                f"row {row_number}"
            )
        if not math.isfinite(confidence) or not 0.0 <= confidence <= 1.0:
            raise ValueError(
                f"prediction table {relative_path} has invalid confidence"
            )
        persisted_class = str(row["effective_eve_class"]).strip().upper()
        if persisted_class not in (
            set(EFFECTIVE_EVE_CLASSES) | set(PPV_LEGACY_ALIASES) | {"MIXED"}
        ):
            raise ValueError(
                f"prediction table {relative_path} has an invalid effective class"
            )
        # Fold legacy VP/PLV onto PPV and legacy MIXED onto VIRAL_UNKNOWN exactly
        # like the manifest summary does, so a pre-migration result directory
        # still validates instead of forcing a full recompute.
        eve_class = normalize_effective_eve_class(persisted_class)
        tier = str(row["confidence_tier"]).strip().upper()
        if tier not in tier_counts:
            raise ValueError(
                f"prediction table {relative_path} has an invalid confidence tier"
            )
        coordinates[eve_id] = (
            scaffold,
            start,
            end,
            length,
            tier,
            confidence,
            eve_class,
        )
        accepted_bp += length
        class_counts[eve_class] += 1
        tier_counts[tier] += 1
    return coordinates, {
        "accepted_bp": accepted_bp,
        "class_counts": class_counts,
        "tier_counts": tier_counts,
    }


def _validate_bed_export(
    root: Path,
    relative_path: str,
    canonical: Mapping[
        str,
        tuple[str, int, int, int, str, float, str],
    ],
) -> None:
    descriptor, metadata = _open_artifact_no_follow(root, relative_path)
    observed: dict[str, tuple[str, int, int, int, str]] = {}
    with os.fdopen(descriptor, "r", encoding="utf-8", closefd=True) as handle:
        for row_number, line in enumerate(handle, start=1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                raise ValueError(f"BED row {row_number} does not contain six fields")
            try:
                start, end, score = int(fields[1]), int(fields[2]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f"BED row {row_number} has invalid numerics") from exc
            if (
                start < 0
                or end <= start
                or not 0 <= score <= 1000
                or fields[5] != "."
            ):
                raise ValueError(f"BED row {row_number} violates BED6 semantics")
            if fields[3] in observed:
                raise ValueError(f"BED contains duplicate EVE ID {fields[3]!r}")
            observed[fields[3]] = (fields[0], start, end, score, fields[5])
        after = os.fstat(handle.fileno())
    expected = {
        eve_id: (
            values[0],
            values[1],
            values[2],
            int(min(1000, values[5] * 1000)),
            ".",
        )
        for eve_id, values in canonical.items()
    }
    if observed != expected:
        raise ValueError("BED coordinates differ from canonical predictions")
    if (metadata.st_dev, metadata.st_ino, metadata.st_mtime_ns) != (
        after.st_dev,
        after.st_ino,
        after.st_mtime_ns,
    ):
        raise ValueError("BED changed while validating")


def _validate_gff_export(
    root: Path,
    relative_path: str,
    canonical: Mapping[
        str,
        tuple[str, int, int, int, str, float, str],
    ],
) -> None:
    descriptor, metadata = _open_artifact_no_follow(root, relative_path)
    observed: dict[str, tuple[str, int, int, int, str, str, float]] = {}
    has_version = False
    with os.fdopen(descriptor, "r", encoding="utf-8", closefd=True) as handle:
        for row_number, line in enumerate(handle, start=1):
            if line.startswith("##gff-version 3"):
                has_version = True
                continue
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if (
                len(fields) != 9
                or fields[1] != "ViroSync"
                or fields[2] != "EVE"
                or fields[6] != "."
                or fields[7] != "."
            ):
                raise ValueError(f"GFF3 row {row_number} has an invalid schema")
            try:
                start, end = int(fields[3]) - 1, int(fields[4])
                score = int(fields[5])
            except ValueError as exc:
                raise ValueError(f"GFF3 row {row_number} has invalid coordinates") from exc
            attributes = {
                key: unquote(value)
                for token in fields[8].split(";")
                if "=" in token
                for key, value in [token.split("=", 1)]
            }
            eve_id = attributes.get("ID", "")
            try:
                confidence = float(attributes.get("confidence", ""))
            except ValueError as exc:
                raise ValueError(
                    f"GFF3 row {row_number} has invalid confidence"
                ) from exc
            if (
                not eve_id
                or eve_id in observed
                or attributes.get("Name") != eve_id
                or start < 0
                or end <= start
                or not 0 <= score <= 1000
                or not math.isfinite(confidence)
            ):
                raise ValueError(f"GFF3 row {row_number} violates coordinate semantics")
            observed[eve_id] = (
                unquote(fields[0]),
                start,
                end,
                score,
                fields[6],
                fields[7],
                confidence,
            )
        after = os.fstat(handle.fileno())
    expected = {
        eve_id: (
            values[0],
            values[1],
            values[2],
            int(min(1000, values[5] * 1000)),
            ".",
            ".",
            values[5],
        )
        for eve_id, values in canonical.items()
    }
    if not has_version or observed != expected:
        raise ValueError("GFF3 coordinates differ from canonical predictions")
    if (metadata.st_dev, metadata.st_ino, metadata.st_mtime_ns) != (
        after.st_dev,
        after.st_ino,
        after.st_mtime_ns,
    ):
        raise ValueError("GFF3 changed while validating")


def _validate_notebook(root: Path) -> None:
    notebook = _read_artifact_json(
        root,
        "notebooks/jupyter/eve_analysis.ipynb",
    )
    if (
        notebook.get("nbformat") != 4
        or not isinstance(notebook.get("nbformat_minor"), int)
        or not isinstance(notebook.get("cells"), list)
        or not isinstance(notebook.get("metadata"), dict)
    ):
        raise ValueError("EVE analysis notebook is not valid nbformat v4 JSON")


def _validate_run_log(root: Path) -> None:
    descriptor, metadata = _open_artifact_no_follow(root, "run.log")
    if metadata.st_size > _MAX_STATE_BYTES:
        os.close(descriptor)
        raise ValueError("run.log is unexpectedly large")
    with os.fdopen(descriptor, "r", encoding="utf-8", closefd=True) as handle:
        content = handle.read()
    if not content.startswith("# ViroSync Run Log:") or "## Results Summary" not in content:
        raise ValueError("run.log does not contain the ViroSync completion contract")


def _validate_ablation_event_binding(
    root: Path,
    relative_path: str,
    identities: Mapping[str, object],
) -> AblationEvents:
    descriptor, metadata = _open_artifact_no_follow(root, relative_path)
    if metadata.st_size > MAX_ABLATION_EVENTS_BYTES:
        os.close(descriptor)
        raise ValueError("ablation event document is unexpectedly large")
    with os.fdopen(descriptor, "rb", closefd=True) as handle:
        events = validate_ablation_events_bytes(handle.read())
    config_identity = identities.get("config")
    if not isinstance(config_identity, Mapping):
        raise ValueError("run identity has no ablation-bound config identity")
    if events.ablation_id.value != config_identity.get("ablation_id"):
        raise ValueError("ablation event ID differs from the run identity")
    if config_identity.get("ablation_contract_sha256") != ABLATION_CONTRACT_SHA256:
        raise ValueError("ablation event contract differs from the run identity")
    return events


def _validate_success_artifacts(
    root: Path,
    artifacts: Sequence[ArtifactIdentity | Mapping[str, object]],
    result: Mapping[str, object],
    *,
    run_fingerprint: str,
    identities: Mapping[str, object],
) -> tuple[ArtifactIdentity, ...]:
    normalized = tuple(
        sorted(
            (_coerce_artifact(artifact) for artifact in artifacts),
            key=lambda artifact: artifact.relative_path,
        )
    )
    by_path = _require_final_artifact_set(normalized)
    if any(not validate_artifact_identity(artifact, root=root) for artifact in normalized):
        raise ValueError("success artifact identity validation failed")
    _validate_ablation_event_binding(
        root,
        "ablation_events.json",
        identities,
    )
    canonical = [
        artifact
        for artifact in normalized
        if Path(artifact.relative_path).name == "virosync_predictions.tsv"
    ]
    detailed = [
        artifact
        for artifact in normalized
        if Path(artifact.relative_path).name == "virosync_predictions_detailed.tsv"
    ]
    if not canonical or not detailed:
        raise ValueError("success state requires canonical and detailed prediction tables")
    if any(artifact.row_count != result["canonical_rows"] for artifact in canonical):
        raise ValueError("canonical prediction row count disagrees with summary")
    if any(artifact.row_count != result["detailed_rows"] for artifact in detailed):
        raise ValueError("detailed prediction row count disagrees with summary")
    for suffix in (".bed", ".gff3"):
        exports = [
            artifact
            for artifact in normalized
            if artifact.relative_path.endswith(f"virosync_predictions{suffix}")
        ]
        if any(artifact.row_count != result["canonical_rows"] for artifact in exports):
            raise ValueError(
                f"canonical {suffix} row count disagrees with summary"
            )
    for relative_path in (
        "run.log",
        "virosync_run_complete.json",
        "virosync_tsv_invariant_report.tsv",
        "notebooks/jupyter/eve_analysis.ipynb",
    ):
        if by_path[relative_path].size <= 0:
            raise ValueError(f"required final artifact is empty: {relative_path}")
    canonical_paths = {
        artifact.relative_path: _prediction_coordinates(
            root,
            artifact.relative_path,
        )
        for artifact in canonical
    }
    detailed_paths = {
        artifact.relative_path: _prediction_coordinates(
            root,
            artifact.relative_path,
        )
        for artifact in detailed
    }
    if len(
        {
            (artifact.size, artifact.sha256, artifact.row_count)
            for artifact in canonical
        }
    ) != 1:
        raise ValueError("duplicate canonical prediction tables disagree")
    if len(
        {
            (artifact.size, artifact.sha256, artifact.row_count)
            for artifact in detailed
        }
    ) != 1:
        raise ValueError("duplicate detailed prediction tables disagree")
    canonical_path = next(
        path
        for path in (
            "phase3_synthesis/virosync_predictions.tsv",
            "virosync_predictions.tsv",
        )
        if path in canonical_paths
    )
    detailed_path = next(
        path
        for path in (
            "virosync_predictions_detailed.tsv",
            "phase3_synthesis/virosync_predictions_detailed.tsv",
        )
        if path in detailed_paths
    )
    canonical_coordinates, canonical_counts = canonical_paths[canonical_path]
    detailed_coordinates, _detailed_counts = detailed_paths[detailed_path]
    if any(value != canonical_paths[canonical_path] for value in canonical_paths.values()):
        raise ValueError("duplicate canonical prediction tables differ semantically")
    if any(value != detailed_paths[detailed_path] for value in detailed_paths.values()):
        raise ValueError("duplicate detailed prediction tables differ semantically")
    if canonical_counts != {
        "accepted_bp": result["accepted_bp"],
        "class_counts": result["class_counts"],
        "tier_counts": result["tier_counts"],
    }:
        raise ValueError("canonical prediction counts disagree with run result")
    scaffold_lengths = _authenticated_scaffold_lengths(identities)
    for rows in (canonical_coordinates, detailed_coordinates):
        for values in rows.values():
            scaffold_length = scaffold_lengths.get(values[0])
            if scaffold_length is None or values[2] > scaffold_length:
                raise ValueError(
                    "prediction lies outside the authenticated input FASTA"
                )
    if any(
        detailed_coordinates.get(eve_id) != coordinates
        for eve_id, coordinates in canonical_coordinates.items()
    ):
        raise ValueError(
            "canonical prediction coordinates differ from detailed predictions"
        )
    for bed_path in (
        "virosync_predictions.bed",
        "phase3_synthesis/virosync_predictions.bed",
    ):
        if bed_path in by_path:
            _validate_bed_export(root, bed_path, canonical_coordinates)
    for gff_path in (
        "virosync_predictions.gff3",
        "phase3_synthesis/virosync_predictions.gff3",
    ):
        if gff_path in by_path:
            _validate_gff_export(root, gff_path, canonical_coordinates)
    _validate_notebook(root)
    _validate_run_log(root)
    _validate_invariant_report(
        root / "virosync_tsv_invariant_report.tsv",
        expected_rows=result["detailed_rows"],  # type: ignore[arg-type]
    )
    _validate_completion_manifest(
        root / "virosync_run_complete.json",
        run_fingerprint=run_fingerprint,
        benchmark_eligible=result["benchmark_eligible"],  # type: ignore[arg-type]
        identities=identities,
    )
    expected_statistics = {
        "canonical_predictions": result["canonical_rows"],
        "total_candidates": result["detailed_rows"],
        "total_accepted_length_bp": result["accepted_bp"],
        "high_confidence": result["tier_counts"]["HIGH"],  # type: ignore[index]
        "medium_confidence": result["tier_counts"]["MEDIUM"],  # type: ignore[index]
        "low_confidence": result["tier_counts"]["LOW"],  # type: ignore[index]
        "promoted_low_confidence": result["promoted_low_rows"],
    }
    summary_contracts: list[dict[str, object]] = []
    for summary_path in (
        "virosync_summary.json",
        "phase3_synthesis/virosync_summary.json",
    ):
        if summary_path not in by_path:
            continue
        summary = _load_json_document(root / summary_path, "ViroSync summary")
        for field in (
            "coordinate_schema_version",
            "coordinate_convention",
            "output_schema_version",
        ):
            if summary.get(field) != identities.get(field):
                raise ValueError(
                    f"ViroSync summary {field} differs from the run identity"
                )
        statistics = summary.get("statistics")
        if not isinstance(statistics, dict):
            raise ValueError("ViroSync summary has no statistics object")
        for field, expected in expected_statistics.items():
            if statistics.get(field) != expected:
                raise ValueError(
                    f"ViroSync summary {field} disagrees with run result"
                )
        if not isinstance(summary.get("virosync_version"), str) or not isinstance(
            summary.get("per_scaffold"),
            dict,
        ):
            raise ValueError("ViroSync summary lacks version/per-scaffold metadata")
        summary_contracts.append(
            {
                key: summary.get(key)
                for key in (
                    "virosync_version",
                    "coordinate_schema_version",
                    "coordinate_convention",
                    "output_schema_version",
                    "statistics",
                    "per_scaffold",
                )
            }
        )
    if any(contract != summary_contracts[0] for contract in summary_contracts[1:]):
        raise ValueError("duplicate ViroSync summaries disagree")
    return normalized


def _success_state_is_valid(root: Path, state: RunState) -> bool:
    if state.status != "success" or state.result is None:
        return False
    try:
        result = _validated_result(state.result)
        _validate_success_artifacts(
            root,
            state.artifacts,
            result,
            run_fingerprint=state.run_fingerprint,
            identities=state.identities,
        )
    except (OSError, TypeError, ValueError):
        return False
    return True


def _success_artifacts_match_terminal_record(
    artifacts: Sequence[ArtifactIdentity],
    records: Sequence[PhaseRecord],
) -> bool:
    if not records:
        return False
    recorded = {
        artifact.relative_path: artifact for artifact in records[-1].artifacts
    }
    return all(recorded.get(artifact.relative_path) == artifact for artifact in artifacts)


@_scoped_validation_cache
def publish_run_success(
    output_dir: str | Path,
    *,
    artifacts: Sequence[ArtifactIdentity | Mapping[str, object]],
    result: Mapping[str, object],
    run_fingerprint: str,
) -> RunState:
    """Validate final outputs and phase chain before publishing ``success``."""

    root = _require_directory_no_follow(Path(output_dir), "output directory")
    current = _load_run_state(root)
    if current.status != "running":
        raise ValueError("only a running attempt can transition to success")
    if current.run_fingerprint != run_fingerprint:
        raise ValueError("success state run fingerprint mismatch")
    normalized_result = _validated_result(result)
    normalized_artifacts = _validate_success_artifacts(
        root,
        artifacts,
        normalized_result,
        run_fingerprint=run_fingerprint,
        identities=current.identities,
    )
    phase_records, reason = _validated_phase_prefix(root, current.run_fingerprint)
    if not _success_artifacts_match_terminal_record(
        normalized_artifacts,
        phase_records,
    ):
        raise ValueError(
            "success artifacts are not authenticated by the final phase marker"
        )
    terminal_phase = normalized_result["terminal_phase"]
    if terminal_phase is not None:
        if (
            not phase_records
            or phase_records[-1].outcome not in _TERMINAL_PHASE_OUTCOMES
            or phase_records[-1].phase != terminal_phase
        ):
            raise ValueError("terminal success requires a terminal phase marker")
        final_outcome = phase_records[-1].outcome
        if final_outcome == "terminal_zero":
            if any(
                (
                    normalized_result["canonical_rows"],
                    normalized_result["detailed_rows"],
                    normalized_result["accepted_bp"],
                    *normalized_result["class_counts"].values(),
                    *normalized_result["tier_counts"].values(),
                )
            ):
                raise ValueError(
                    "terminal-zero summaries must contain all-zero result counts"
                )
        else:
            config_identity = current.identities.get("config")
            if (
                terminal_phase != 1
                or not isinstance(config_identity, Mapping)
                or config_identity.get("ablation_id") != AblationID.A1.value
            ):
                raise ValueError(
                    "terminal-ablation success is reserved for the A1 Phase-1 surface"
                )
            if (
                normalized_result["canonical_rows"] <= 0
                or normalized_result["detailed_rows"] <= 0
                or normalized_result["accepted_bp"] <= 0
            ):
                raise ValueError(
                    "terminal-ablation success requires nonzero authenticated output"
                )
            if normalized_result["promoted_low_rows"] != 0:
                raise ValueError(
                    "terminal-ablation LOW rows are unscored, not promoted"
                )
    elif len(phase_records) != len(PHASE_MARKER_FILENAMES):
        raise ValueError(reason or "success requires all four phase markers")
    for artifact in normalized_artifacts:
        _fsync_artifact(root, artifact)
    success = RunState(
        schema_version=RUN_STATE_SCHEMA_VERSION,
        run_fingerprint=current.run_fingerprint,
        status="success",
        attempt=current.attempt,
        identities=current.identities,
        result=normalized_result,
        artifacts=normalized_artifacts,
    )
    atomic_write_json(root / RUN_STATE_FILENAME, success)
    return success


@_scoped_validation_cache
def plan_resume(
    output_dir: str | Path,
    *,
    expected_run_fingerprint: str,
) -> ResumePlan:
    """Validate state and return only the reusable sequential phase prefix."""

    _require_sha256(expected_run_fingerprint, "expected_run_fingerprint")
    root = Path(output_dir)
    try:
        root = _require_directory_no_follow(root, "output directory")
        state = _load_run_state(root)
    except (FileNotFoundError, OSError, TypeError, ValueError) as exc:
        return ResumePlan((), 0, reason=f"missing or invalid run state: {exc}")
    if state.run_fingerprint != expected_run_fingerprint:
        return ResumePlan((), 0, reason="run fingerprint changed")
    records, reason = _validated_phase_prefix(root, expected_run_fingerprint)
    phases = tuple(record.phase for record in records)
    terminal_phase = (
        records[-1].phase
        if records and records[-1].outcome in _TERMINAL_PHASE_OUTCOMES
        else None
    )
    if state.status == "success":
        if not _success_state_is_valid(root, state):
            return ResumePlan(
                phases,
                len(records),
                terminal_phase=terminal_phase,
                reason=reason or "final success artifacts are stale",
            )
        if not _success_artifacts_match_terminal_record(state.artifacts, records):
            return ResumePlan(
                phases,
                len(records),
                terminal_phase=terminal_phase,
                reason="final success artifacts are not bound to the last phase",
            )
        assert state.result is not None
        if state.result["terminal_phase"] != terminal_phase:
            return ResumePlan(
                phases,
                len(records),
                terminal_phase=terminal_phase,
                reason="terminal state disagrees with phases",
            )
        if terminal_phase is None and len(records) != len(PHASE_MARKER_FILENAMES):
            return ResumePlan(
                phases,
                len(records),
                reason=reason or "phase chain is incomplete",
            )
        return ResumePlan(
            phases,
            len(PHASE_MARKER_FILENAMES),
            terminal_phase=terminal_phase,
            completed=True,
        )
    restart_phase = len(records)
    return ResumePlan(
        phases,
        restart_phase,
        terminal_phase=terminal_phase,
        reason=reason,
    )


def _open_directory_fd(path: Path) -> int:
    descriptor = os.open(path, _DIRECTORY_OPEN_FLAGS)
    metadata = os.fstat(descriptor)
    if not stat.S_ISDIR(metadata.st_mode):
        os.close(descriptor)
        raise ValueError(f"path is not a real directory: {path}")
    return descriptor


def _remove_entry_at(parent_fd: int, name: str) -> None:
    try:
        metadata = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except FileNotFoundError:
        return
    if not stat.S_ISDIR(metadata.st_mode) or stat.S_ISLNK(metadata.st_mode):
        try:
            os.unlink(name, dir_fd=parent_fd)
        except FileNotFoundError:
            pass
        return

    try:
        directory_fd = os.open(name, _DIRECTORY_OPEN_FLAGS, dir_fd=parent_fd)
    except OSError as exc:
        raise ValueError(
            f"directory changed during guarded invalidation: {name}"
        ) from exc
    try:
        opened = os.fstat(directory_fd)
        if (opened.st_dev, opened.st_ino) != (metadata.st_dev, metadata.st_ino):
            raise ValueError(
                f"directory changed during guarded invalidation: {name}"
            )
        with os.scandir(directory_fd) as entries:
            child_names = [entry.name for entry in entries]
        for child_name in child_names:
            _remove_entry_at(directory_fd, child_name)
    finally:
        os.close(directory_fd)

    try:
        current = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except FileNotFoundError:
        return
    if (
        not stat.S_ISDIR(current.st_mode)
        or (current.st_dev, current.st_ino) != (metadata.st_dev, metadata.st_ino)
    ):
        raise ValueError(f"directory changed during guarded invalidation: {name}")
    os.rmdir(name, dir_fd=parent_fd)


def _remove_tree_no_follow(path: Path) -> None:
    parent_fd = _open_directory_fd(path.parent)
    try:
        _remove_entry_at(parent_fd, path.name)
    finally:
        os.close(parent_fd)


def _remove_relative_no_follow(root: Path, relative_path: str) -> None:
    relative = _normalized_relative_path(relative_path, "recorded artifact path")
    root_fd = _open_directory_fd(root)
    current_fd = root_fd
    opened_fds: list[int] = []
    try:
        parts = PurePosixPath(relative).parts
        for component in parts[:-1]:
            try:
                metadata = os.stat(
                    component,
                    dir_fd=current_fd,
                    follow_symlinks=False,
                )
            except FileNotFoundError:
                return
            if not stat.S_ISDIR(metadata.st_mode) or stat.S_ISLNK(metadata.st_mode):
                _remove_entry_at(current_fd, component)
                return
            try:
                next_fd = os.open(
                    component,
                    _DIRECTORY_OPEN_FLAGS,
                    dir_fd=current_fd,
                )
            except OSError as exc:
                raise ValueError(
                    "artifact parent changed during guarded invalidation: "
                    f"{component}"
                ) from exc
            opened = os.fstat(next_fd)
            if (opened.st_dev, opened.st_ino) != (metadata.st_dev, metadata.st_ino):
                os.close(next_fd)
                raise ValueError(
                    "artifact parent changed during guarded invalidation: "
                    f"{component}"
                )
            opened_fds.append(next_fd)
            current_fd = next_fd
        _remove_entry_at(current_fd, parts[-1])
    finally:
        for descriptor in reversed(opened_fds):
            os.close(descriptor)
        os.close(root_fd)


def _is_owned_final_name(name: str) -> bool:
    return (
        name in {*_FINAL_OUTPUTS, "notebooks", "provenance.json"}
        or name.endswith("_eves.fna")
        or (name.startswith("eve_") and name.endswith(".png"))
    )


def invalidate_from_phase(output_dir: str | Path, *, from_phase: int) -> None:
    """Clear success and remove a phase suffix without following symlinks."""

    if from_phase not in range(len(PHASE_MARKER_FILENAMES)):
        raise ValueError("from_phase must be in the range 0..3")
    root = _require_directory_no_follow(Path(output_dir), "output directory")

    recorded_artifacts: set[str] = set()
    state_path = root / RUN_STATE_FILENAME
    try:
        state = _load_run_state(root)
    except FileNotFoundError:
        state = None
    except (OSError, TypeError, ValueError):
        _unlink_no_follow(state_path)
        state = None
    if state is not None and state.status == "success":
        recorded_artifacts.update(
            artifact.relative_path for artifact in state.artifacts
        )
        cleared = RunState(
            schema_version=RUN_STATE_SCHEMA_VERSION,
            run_fingerprint=state.run_fingerprint,
            status="running",
            attempt=state.attempt,
            identities=state.identities,
        )
        atomic_write_json(state_path, cleared)

    for phase in range(from_phase, len(PHASE_MARKER_FILENAMES)):
        marker = root / PHASE_MARKER_FILENAMES[phase]
        try:
            record = _load_phase_record(marker)
        except (FileNotFoundError, OSError, TypeError, ValueError):
            continue
        recorded_artifacts.update(
            artifact.relative_path for artifact in record.artifacts
        )

    for phase in range(from_phase, len(PHASE_MARKER_FILENAMES)):
        _unlink_no_follow(root / PHASE_MARKER_FILENAMES[phase])
    for phase in range(from_phase, len(PHASE_MARKER_FILENAMES)):
        for name in _PHASE_DIRECTORIES[phase]:
            _remove_tree_no_follow(root / name)
    root_fd = _open_directory_fd(root)
    try:
        with os.scandir(root_fd) as entries:
            owned_names = [
                entry.name for entry in entries if _is_owned_final_name(entry.name)
            ]
        for name in owned_names:
            _remove_entry_at(root_fd, name)
    finally:
        os.close(root_fd)
    for relative_path in sorted(recorded_artifacts):
        _remove_relative_no_follow(root, relative_path)
