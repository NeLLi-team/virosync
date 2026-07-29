"""
Genome masking with TRF and RepeatMasker.

Masks low-complexity regions and host-specific repeats to reduce false positives
in downstream homology searches and statistical analyses.

This module provides:
- TRF (Tandem Repeats Finder) integration for tandem repeat masking
- RepeatMasker integration for known repeat element masking
- Combined masking pipeline that converts masked regions to 'N' characters
"""

import hashlib
import json
import logging
import re
import shutil
import stat
import subprocess
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from virosync.config import (
    ConfigError,
    MaskingBackend,
    MaskingConfig,
    MaskingFailurePolicy,
)
from virosync.utils.atomic_write import atomic_write_context
from virosync.utils.path_safety import require_strict_child

logger = logging.getLogger(__name__)


MASKING_STATUS_SCHEMA_VERSION = 1
_LEGACY_VERSION = "version not probed in legacy adapter"
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}\Z")


@dataclass(frozen=True)
class MaskedRegion:
    """Represents a masked genomic region."""

    seq_id: str
    start: int  # 0-based
    end: int  # exclusive
    mask_type: str  # 'trf', 'repeatmasker', 'combined'
    annotation: str = ""


class MaskingBackendError(RuntimeError):
    """A requested masking backend could not produce a valid result."""

    def __init__(self, backend: MaskingBackend | str, reason: str):
        self.backend = MaskingBackend(backend)
        self.reason = reason
        super().__init__(f"{self.backend.value} masking backend failed: {reason}")


@dataclass(frozen=True)
class MaskingResult:
    """Immutable masking outcome used by Phase 0 and persisted status."""

    output_path: Path
    repeat_regions: tuple[MaskedRegion, ...]
    requested_backend: MaskingBackend
    effective_backend: MaskingBackend
    failure_policy: MaskingFailurePolicy
    status: str
    legacy_adapter: bool
    backend_versions: tuple[tuple[str, str], ...]
    masked_bases: int
    repeatmasker_species: Optional[str]
    repeatmasker_library: Optional[Path]
    repeatmasker_library_sha256: Optional[str]
    configured_fallback_backend: Optional[MaskingBackend]
    fallback_backend: Optional[MaskingBackend]
    fallback_reason: Optional[str]
    input_sha256: str
    output_sha256: str
    status_path: Optional[Path] = None
    status_sha256: Optional[str] = None

    @property
    def benchmark_eligible(self) -> bool:
        """Only off or strict verified success may enter primary benchmarks."""
        if _masking_state_errors(self):
            return False
        if self.status == "off":
            return True
        if (
            self.status != "success"
            or self.legacy_adapter
            or self.failure_policy is not MaskingFailurePolicy.STRICT
        ):
            return False
        versions = dict(self.backend_versions)
        expected = {
            backend.value for backend in _requested_tools(self.effective_backend)
        }
        return bool(expected) and all(
            versions.get(backend)
            and versions[backend] not in {"failed", _LEGACY_VERSION}
            for backend in expected
        )

    def to_status_payload(
        self,
        *,
        repeat_region_count: Optional[int] = None,
    ) -> dict:
        """Return the canonical JSON-safe status payload."""
        if repeat_region_count is None:
            repeat_region_count = len(self.repeat_regions)
        payload = {
            "schema_version": MASKING_STATUS_SCHEMA_VERSION,
            "status": self.status,
            "requested_backend": self.requested_backend.value,
            "effective_backend": self.effective_backend.value,
            "failure_policy": self.failure_policy.value,
            "legacy_adapter": self.legacy_adapter,
            "backend_versions": dict(self.backend_versions),
            "masked_bases": self.masked_bases,
            "repeat_region_count": repeat_region_count,
            "repeatmasker_species": self.repeatmasker_species,
            "repeatmasker_library": (
                str(self.repeatmasker_library.resolve())
                if self.repeatmasker_library is not None
                else None
            ),
            "repeatmasker_library_sha256": self.repeatmasker_library_sha256,
            "configured_fallback_backend": (
                self.configured_fallback_backend.value
                if self.configured_fallback_backend is not None
                else None
            ),
            "fallback_backend": (
                self.fallback_backend.value
                if self.fallback_backend is not None
                else None
            ),
            "fallback_reason": self.fallback_reason,
            "input_sha256": self.input_sha256,
            "output_path": str(self.output_path.resolve()),
            "output_sha256": self.output_sha256,
            "benchmark_eligible": self.benchmark_eligible,
        }
        canonical = json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
        payload["result_fingerprint"] = hashlib.sha256(canonical).hexdigest()
        return payload


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_input_fasta(path: Path) -> None:
    """Reject unreadable, empty, or ID-ambiguous inputs before any backend runs."""
    try:
        records = list(SeqIO.parse(path, "fasta"))
    except (OSError, UnicodeError, ValueError) as exc:
        raise ValueError(f"invalid input FASTA {path}: {exc}") from exc
    if not records:
        raise ValueError(f"invalid input FASTA {path}: contains no records")
    seen: set[str] = set()
    duplicates: set[str] = set()
    for record in records:
        if record.id in seen:
            duplicates.add(record.id)
        seen.add(record.id)
    if duplicates:
        joined = ", ".join(repr(value) for value in sorted(duplicates))
        raise ValueError(f"invalid input FASTA {path}: duplicate record IDs {joined}")


def _library_sha256(config: MaskingConfig) -> Optional[str]:
    if config.repeatmasker_library is None:
        return None
    library = Path(config.repeatmasker_library)
    if not library.is_file():
        raise ConfigError(
            f"execution.masking.repeatmasker_library is not a file: {library}"
        )
    if library.stat().st_size == 0:
        raise ConfigError(
            f"execution.masking.repeatmasker_library is empty: {library}"
        )
    return _file_sha256(library)


def _is_verified_version(value: object) -> bool:
    return (
        isinstance(value, str)
        and bool(value.strip())
        and value not in {"failed", _LEGACY_VERSION}
    )


def _masking_state_errors(
    result: MaskingResult,
    *,
    repeat_region_count: Optional[int] = None,
) -> list[str]:
    """Return violations of the persisted masking state machine."""
    errors: list[str] = []
    if result.status not in {"off", "success", "fallback"}:
        errors.append(f"unknown status {result.status!r}")
    if type(result.legacy_adapter) is not bool:
        errors.append("legacy_adapter must be a boolean")
    if type(result.masked_bases) is not int or result.masked_bases < 0:
        errors.append("masked_bases must be a non-negative integer")
    if repeat_region_count is not None and (
        type(repeat_region_count) is not int or repeat_region_count < 0
    ):
        errors.append("repeat_region_count must be a non-negative integer")
    for label, digest in (
        ("input_sha256", result.input_sha256),
        ("output_sha256", result.output_sha256),
    ):
        if not isinstance(digest, str) or _SHA256_PATTERN.fullmatch(digest) is None:
            errors.append(f"{label} must be a lowercase SHA256")

    try:
        versions = dict(result.backend_versions)
    except (TypeError, ValueError):
        versions = {}
        errors.append("backend_versions must contain key/value pairs")
    if len(versions) != len(result.backend_versions):
        errors.append("backend_versions contains duplicate tools")
    if any(not isinstance(key, str) or not isinstance(value, str) for key, value in versions.items()):
        errors.append("backend_versions keys and values must be strings")

    if result.failure_policy is MaskingFailurePolicy.STRICT:
        if result.configured_fallback_backend is not None:
            errors.append("strict policy cannot configure a fallback")
    elif result.failure_policy is MaskingFailurePolicy.FALLBACK:
        if result.configured_fallback_backend not in {
            MaskingBackend.OFF,
            MaskingBackend.TRF,
        }:
            errors.append("fallback policy requires an off or trf fallback")
        if result.requested_backend is MaskingBackend.OFF:
            errors.append("off requests cannot use fallback policy")
        if result.configured_fallback_backend is result.requested_backend:
            errors.append("configured fallback must differ from requested backend")
    else:
        errors.append("unknown failure policy")

    has_species = result.repeatmasker_species is not None
    has_library = result.repeatmasker_library is not None
    if has_species and has_library:
        errors.append("RepeatMasker species and library are mutually exclusive")
    if result.requested_backend in {
        MaskingBackend.REPEATMASKER,
        MaskingBackend.TRF_REPEATMASKER,
    } and not (has_species or has_library):
        errors.append("RepeatMasker request has no explicit target")
    if has_library:
        if (
            not isinstance(result.repeatmasker_library_sha256, str)
            or _SHA256_PATTERN.fullmatch(result.repeatmasker_library_sha256) is None
        ):
            errors.append("RepeatMasker library must have a lowercase SHA256")
    elif result.repeatmasker_library_sha256 is not None:
        errors.append("RepeatMasker library SHA256 has no library path")

    required_tools = {
        backend.value for backend in _requested_tools(result.requested_backend)
    }
    if result.status == "off":
        if result.requested_backend is not MaskingBackend.OFF:
            errors.append("off status requires requested backend off")
        if result.effective_backend is not MaskingBackend.OFF:
            errors.append("off status requires effective backend off")
        if result.legacy_adapter:
            errors.append("off status cannot be a legacy adapter success")
        if versions:
            errors.append("off status cannot record backend versions")
        if result.masked_bases != 0:
            errors.append("off status must mask zero bases")
        if result.input_sha256 != result.output_sha256:
            errors.append("off status input and output SHA256 must match")
        if result.fallback_backend is not None or result.fallback_reason is not None:
            errors.append("off status cannot select a fallback")
        if result.repeat_regions:
            errors.append("off status cannot contain repeat regions")
        if repeat_region_count not in {None, 0}:
            errors.append("off status repeat_region_count must be zero")
    elif result.status == "success":
        if result.requested_backend is MaskingBackend.OFF:
            errors.append("success status requires an enabled backend")
        if result.effective_backend is not result.requested_backend:
            errors.append("success status requires requested and effective backends to match")
        if result.fallback_backend is not None or result.fallback_reason is not None:
            errors.append("success status cannot select a fallback")
        if set(versions) != required_tools:
            errors.append("success status must record exactly the requested tools")
        elif result.legacy_adapter:
            if any(value != _LEGACY_VERSION for value in versions.values()):
                errors.append("legacy adapter success has an invalid version marker")
            if result.failure_policy is not MaskingFailurePolicy.STRICT:
                errors.append("legacy adapter success must use strict policy")
        elif any(not _is_verified_version(value) for value in versions.values()):
            errors.append("success status requires verified backend versions")
    elif result.status == "fallback":
        if result.legacy_adapter:
            errors.append("fallback status cannot be a legacy adapter success")
        if result.failure_policy is not MaskingFailurePolicy.FALLBACK:
            errors.append("fallback status requires fallback policy")
        if result.fallback_backend is not result.configured_fallback_backend:
            errors.append("selected fallback must equal configured fallback")
        if result.effective_backend is not result.fallback_backend:
            errors.append("effective backend must equal selected fallback")
        if not result.fallback_reason:
            errors.append("fallback status requires a failure reason")
        allowed_versions = set(required_tools)
        if result.fallback_backend is MaskingBackend.TRF:
            allowed_versions.add("fallback_trf")
            if not _is_verified_version(versions.get("fallback_trf")):
                errors.append("trf fallback requires a verified fallback version")
        elif "fallback_trf" in versions:
            errors.append("off fallback cannot record fallback_trf")
        if not set(versions).issubset(allowed_versions):
            errors.append("fallback status records an unexpected backend version")
        if not any(versions.get(tool) == "failed" for tool in required_tools):
            errors.append("fallback status must record a failed requested tool")
        if any(
            value != "failed" and not _is_verified_version(value)
            for key, value in versions.items()
            if key in required_tools
        ):
            errors.append("fallback status has an invalid requested-tool version")
        if result.effective_backend is MaskingBackend.OFF:
            if result.masked_bases != 0:
                errors.append("off fallback must mask zero bases")
            if result.input_sha256 != result.output_sha256:
                errors.append("off fallback input and output SHA256 must match")
            if result.repeat_regions:
                errors.append("off fallback cannot contain repeat regions")
            if repeat_region_count not in {None, 0}:
                errors.append("off fallback repeat_region_count must be zero")
    return errors


def _validate_masking_state(
    result: MaskingResult,
    *,
    repeat_region_count: Optional[int] = None,
) -> None:
    errors = _masking_state_errors(
        result,
        repeat_region_count=repeat_region_count,
    )
    if errors:
        raise ValueError("invalid masking state: " + "; ".join(errors))


def _require_regular_file(path: Path, label: str) -> None:
    try:
        mode = Path(path).lstat().st_mode
    except OSError as exc:
        raise ValueError(f"{label} is not a readable regular file: {path}") from exc
    if not stat.S_ISREG(mode):
        raise ValueError(f"{label} must be a non-symlink regular file: {path}")


def _status_path(output_dir: Path) -> Path:
    return Path(output_dir) / "masking_status.json"


def write_masking_status(result: MaskingResult, output_dir: Path) -> MaskingResult:
    """Atomically persist a successful/off/fallback result and attach its digest."""
    status_path = _status_path(output_dir)
    status_path.parent.mkdir(parents=True, exist_ok=True)
    with atomic_write_context(status_path, "w") as handle:
        json.dump(result.to_status_payload(), handle, indent=2, sort_keys=True)
        handle.write("\n")
    _require_regular_file(status_path, "masking status")
    return replace(
        result,
        status_path=status_path,
        status_sha256=_file_sha256(status_path),
    )


def _write_failed_masking_status(
    *,
    input_fasta: Path,
    output_dir: Path,
    config: MaskingConfig,
    error: MaskingBackendError,
    backend_versions: dict[str, str],
    selected_fallback_backend: Optional[MaskingBackend] = None,
) -> Path:
    """Atomically persist a typed backend failure without claiming output success."""
    library_sha256 = _library_sha256(config)
    payload = {
        "schema_version": MASKING_STATUS_SCHEMA_VERSION,
        "status": "failed",
        "requested_backend": config.backend.value,
        "effective_backend": None,
        "failure_policy": config.failure_policy.value,
        "legacy_adapter": False,
        "backend_versions": dict(sorted(backend_versions.items())),
        "masked_bases": 0,
        "repeat_region_count": 0,
        "repeatmasker_species": config.repeatmasker_species,
        "repeatmasker_library": (
            str(config.repeatmasker_library.resolve())
            if config.repeatmasker_library is not None
            else None
        ),
        "repeatmasker_library_sha256": library_sha256,
        "configured_fallback_backend": (
            config.fallback_backend.value
            if config.fallback_backend is not None
            else None
        ),
        "fallback_backend": (
            selected_fallback_backend.value
            if selected_fallback_backend is not None
            else None
        ),
        "fallback_reason": str(error),
        "input_sha256": _file_sha256(input_fasta),
        "output_path": None,
        "output_sha256": None,
        "benchmark_eligible": False,
    }
    canonical = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    payload["result_fingerprint"] = hashlib.sha256(canonical).hexdigest()
    status_path = _status_path(output_dir)
    status_path.parent.mkdir(parents=True, exist_ok=True)
    with atomic_write_context(status_path, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    _require_regular_file(status_path, "masking status")
    return status_path


def _load_status_payload(status_path: Path) -> dict:
    """Load and validate the stable status schema and semantic fingerprint."""
    status_path = Path(status_path)
    _require_regular_file(status_path, "masking status")
    try:
        payload = json.loads(status_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid masking status {status_path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"invalid masking status {status_path}: root must be an object")
    if payload.get("schema_version") != MASKING_STATUS_SCHEMA_VERSION:
        raise ValueError(
            f"invalid masking status schema in {status_path}: "
            f"{payload.get('schema_version')!r}"
        )
    fingerprint = payload.get("result_fingerprint")
    semantic = dict(payload)
    semantic.pop("result_fingerprint", None)
    expected = hashlib.sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if fingerprint != expected:
        raise ValueError(f"masking status semantic fingerprint mismatch: {status_path}")
    return payload


def masking_result_from_status_payload(
    payload: dict,
    *,
    status_path: Path,
    status_sha256: str,
    repeat_regions: Optional[tuple[MaskedRegion, ...]] = None,
) -> MaskingResult:
    """Reconstruct and semantically validate a decoded masking-status payload."""

    if payload.get("status") not in {"off", "success", "fallback"}:
        raise ValueError(
            f"masking status is not reusable: {payload.get('status')!r}"
        )
    if payload.get("legacy_adapter") is True:
        raise ValueError("legacy adapter masking status is not reusable")
    output_path_raw = payload.get("output_path")
    if not output_path_raw:
        raise ValueError("masking status has no output_path")
    result = MaskingResult(
        output_path=Path(output_path_raw),
        repeat_regions=repeat_regions or (),
        requested_backend=MaskingBackend(payload["requested_backend"]),
        effective_backend=MaskingBackend(payload["effective_backend"]),
        failure_policy=MaskingFailurePolicy(payload["failure_policy"]),
        status=str(payload["status"]),
        legacy_adapter=payload.get("legacy_adapter", False),
        backend_versions=tuple(
            sorted((str(key), str(value)) for key, value in payload["backend_versions"].items())
        ),
        masked_bases=int(payload["masked_bases"]),
        repeatmasker_species=payload.get("repeatmasker_species"),
        repeatmasker_library=(
            Path(payload["repeatmasker_library"])
            if payload.get("repeatmasker_library")
            else None
        ),
        repeatmasker_library_sha256=payload.get("repeatmasker_library_sha256"),
        configured_fallback_backend=(
            MaskingBackend(payload["configured_fallback_backend"])
            if payload.get("configured_fallback_backend")
            else None
        ),
        fallback_backend=(
            MaskingBackend(payload["fallback_backend"])
            if payload.get("fallback_backend")
            else None
        ),
        fallback_reason=payload.get("fallback_reason"),
        input_sha256=str(payload["input_sha256"]),
        output_sha256=str(payload["output_sha256"]),
        status_path=Path(status_path),
        status_sha256=status_sha256,
    )
    persisted_region_count = payload.get("repeat_region_count")
    _validate_masking_state(
        result,
        repeat_region_count=(
            len(result.repeat_regions)
            if repeat_regions is not None
            else persisted_region_count
        ),
    )
    expected_payload = result.to_status_payload(
        repeat_region_count=(
            None if repeat_regions is not None else persisted_region_count
        )
    )
    mismatched = sorted(
        key
        for key in set(payload) | set(expected_payload)
        if payload.get(key) != expected_payload.get(key)
    )
    if mismatched:
        raise ValueError(
            "masking status result mismatch: " + ", ".join(mismatched)
        )
    return result


def load_masking_result(
    status_path: Path,
    *,
    repeat_regions: Optional[tuple[MaskedRegion, ...]] = None,
    expected_config: Optional[MaskingConfig] = None,
    expected_input: Optional[Path] = None,
) -> MaskingResult:
    """Reconstruct and validate a successful/off/fallback result from status."""
    payload = _load_status_payload(status_path)
    result = masking_result_from_status_payload(
        payload,
        status_path=Path(status_path),
        status_sha256=_file_sha256(Path(status_path)),
        repeat_regions=repeat_regions,
    )
    if expected_config is not None:
        expected_library = (
            str(expected_config.repeatmasker_library.resolve())
            if expected_config.repeatmasker_library is not None
            else None
        )
        mismatches = []
        if result.requested_backend is not expected_config.backend:
            mismatches.append("requested backend")
        if result.failure_policy is not expected_config.failure_policy:
            mismatches.append("failure policy")
        if result.repeatmasker_species != expected_config.repeatmasker_species:
            mismatches.append("RepeatMasker species")
        observed_library = (
            str(result.repeatmasker_library.resolve())
            if result.repeatmasker_library is not None
            else None
        )
        if observed_library != expected_library:
            mismatches.append("RepeatMasker library")
        if result.repeatmasker_library_sha256 != _library_sha256(expected_config):
            mismatches.append("RepeatMasker library SHA256")
        if result.configured_fallback_backend is not expected_config.fallback_backend:
            mismatches.append("configured fallback backend")
        if mismatches:
            raise ValueError(
                "masking status request mismatch: " + ", ".join(mismatches)
            )
    validate_masking_result(
        result,
        expected_input=expected_input,
        verify_repeat_region_count=repeat_regions is not None,
    )
    return result


def validate_masking_result(
    result: MaskingResult,
    *,
    expected_input: Optional[Path] = None,
    verify_repeat_region_count: bool = True,
) -> None:
    """Verify status and output identities before passing output downstream."""
    if result.status_path is None or result.status_sha256 is None:
        raise ValueError("masking result is missing status path or status SHA256")
    _require_regular_file(result.status_path, "masking status")
    if _file_sha256(result.status_path) != result.status_sha256:
        raise ValueError("masking status file SHA256 mismatch")
    payload = _load_status_payload(result.status_path)
    persisted_region_count = payload.get("repeat_region_count")
    _validate_masking_state(
        result,
        repeat_region_count=(
            len(result.repeat_regions)
            if verify_repeat_region_count
            else persisted_region_count
        ),
    )
    expected_payload = result.to_status_payload(
        repeat_region_count=(
            None if verify_repeat_region_count else persisted_region_count
        )
    )
    mismatched = sorted(
        key
        for key in set(payload) | set(expected_payload)
        if payload.get(key) != expected_payload.get(key)
    )
    if mismatched:
        raise ValueError(
            "masking status result mismatch: " + ", ".join(mismatched)
        )
    if result.effective_backend is MaskingBackend.OFF:
        if not result.output_path.is_file():
            raise ValueError(f"masking output is not a file: {result.output_path}")
    else:
        _require_regular_file(result.output_path, "enabled masking output")
    if _file_sha256(result.output_path) != result.output_sha256:
        raise ValueError("masking output SHA256 mismatch")
    if expected_input is not None:
        if _file_sha256(expected_input) != result.input_sha256:
            raise ValueError("masking input SHA256 mismatch")
        if (
            result.effective_backend is MaskingBackend.OFF
            and result.output_path.resolve() != Path(expected_input).resolve()
        ):
            raise ValueError("off masking output path does not match the input")


def run_trf(
    input_fasta: Path,
    output_dir: Path,
    match: int = 2,
    mismatch: int = 7,
    delta: int = 7,
    pm: int = 80,
    pi: int = 10,
    minscore: int = 50,
    maxperiod: int = 500,
) -> Path:
    """
    Run Tandem Repeats Finder on input FASTA.

    Args:
        input_fasta: Input genome FASTA file
        output_dir: Directory for TRF output files
        match, mismatch, delta: Alignment scoring parameters
        pm, pi: Match and indel probabilities
        minscore: Minimum alignment score
        maxperiod: Maximum repeat period size

    Returns:
        Path to TRF .dat output file

    Raises:
        subprocess.CalledProcessError: If TRF fails
        FileNotFoundError: If TRF output is not found
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # TRF writes output to current directory, so we run from output_dir
    cmd = [
        "trf",
        str(input_fasta.absolute()),
        str(match),
        str(mismatch),
        str(delta),
        str(pm),
        str(pi),
        str(minscore),
        str(maxperiod),
        "-d",  # Produce data file
        "-h",  # Suppress HTML output
    ]

    logger.info(f"Running TRF: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, cwd=str(output_dir), check=True, capture_output=True)
    except FileNotFoundError as exc:
        raise MaskingBackendError(MaskingBackend.TRF, "trf executable not found") from exc
    except subprocess.CalledProcessError as exc:
        stderr = (exc.stderr or b"")
        if isinstance(stderr, bytes):
            stderr = stderr.decode(errors="replace")
        reason = stderr.strip() or f"trf exited with status {exc.returncode}"
        raise MaskingBackendError(MaskingBackend.TRF, reason) from exc
    except OSError as exc:
        raise MaskingBackendError(
            MaskingBackend.TRF,
            f"cannot execute trf: {exc}",
        ) from exc

    # Find the .dat output file
    dat_files = list(output_dir.glob("*.dat"))
    if not dat_files:
        raise MaskingBackendError(
            MaskingBackend.TRF,
            f"trf did not produce .dat output in {output_dir}",
        )
    if len(dat_files) != 1:
        raise MaskingBackendError(
            MaskingBackend.TRF,
            f"trf produced {len(dat_files)} .dat outputs instead of one",
        )
    dat_file = dat_files[0]
    try:
        _require_regular_file(dat_file, "TRF output")
        if dat_file.stat().st_size == 0:
            raise ValueError(f"TRF output is empty: {dat_file}")
    except (OSError, ValueError) as exc:
        raise MaskingBackendError(MaskingBackend.TRF, str(exc)) from exc
    return dat_file


def parse_trf_output(trf_dat: Path) -> list[MaskedRegion]:
    """
    Parse TRF .dat output file to extract masked regions.

    TRF .dat format:
    - Lines starting with "Sequence:" indicate new sequence
    - Data lines contain: start end period_size copy_number ...

    Args:
        trf_dat: Path to TRF .dat file

    Returns:
        List of MaskedRegion objects for tandem repeats
    """
    regions = []
    current_seq = None

    with open(trf_dat) as f:
        for line in f:
            line = line.strip()

            # New sequence header
            if line.startswith("Sequence:"):
                header = line.split(maxsplit=1)
                if len(header) != 2 or not header[1].strip():
                    raise ValueError(f"malformed TRF sequence header: {line!r}")
                current_seq = header[1].split()[0]

            # Skip headers and empty lines
            elif not line or line.startswith("Parameters") or line.startswith("@"):
                continue

            # Data line: start end period_size copy_number ...
            elif current_seq:
                parts = line.split()
                if parts and re.match(r"[+-]?\d", parts[0]):
                    if len(parts) < 2:
                        raise ValueError(f"malformed TRF data line: {line!r}")
                    try:
                        start = int(parts[0])
                        end = int(parts[1])
                        # TRF uses 1-based coordinates, convert to 0-based
                        regions.append(
                            MaskedRegion(
                                seq_id=current_seq,
                                start=start - 1,
                                end=end,
                                mask_type="trf",
                                annotation=f"period={parts[2]}" if len(parts) > 2 else "",
                            )
                        )
                    except ValueError as exc:
                        raise ValueError(
                            f"malformed TRF coordinates: {line!r}"
                        ) from exc

    logger.info(f"Parsed {len(regions)} tandem repeat regions from TRF output")
    return regions


def run_repeatmasker(
    input_fasta: Path,
    output_dir: Path,
    species: Optional[str] = None,
    library: Optional[Path] = None,
    threads: int = 4,
    engine: str = "rmblast",
) -> Path:
    """
    Run RepeatMasker on input FASTA.

    Args:
        input_fasta: Input genome FASTA file
        output_dir: Directory for RepeatMasker output
        species: Explicit RepeatMasker taxonomy target
        library: Explicit custom RepeatMasker library
        threads: Number of parallel threads
        engine: Search engine (rmblast, hmmer, etc.)

    Returns:
        Path to masked output FASTA

    Note:
        RepeatMasker soft-masks repeats (lowercase). We later convert to hard mask (N).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    if bool(species) == bool(library):
        raise ConfigError(
            "RepeatMasker requires exactly one of species or library"
        )

    cmd = [
        "RepeatMasker",
        "-e",
        engine,
        "-pa",
        str(threads),
        "-dir",
        str(output_dir),
        "-xsmall",  # Soft-mask (lowercase)
        "-nolow",  # Don't mask low complexity
        str(input_fasta),
    ]
    target = ["-species", str(species)] if species else ["-lib", str(library)]
    cmd[3:3] = target

    logger.info(f"Running RepeatMasker: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        reason = (e.stderr or "").strip() or (
            f"RepeatMasker exited with status {e.returncode}"
        )
        raise MaskingBackendError(MaskingBackend.REPEATMASKER, reason) from e
    except FileNotFoundError as exc:
        raise MaskingBackendError(
            MaskingBackend.REPEATMASKER,
            "RepeatMasker executable not found",
        ) from exc
    except OSError as exc:
        raise MaskingBackendError(
            MaskingBackend.REPEATMASKER,
            f"cannot execute RepeatMasker: {exc}",
        ) from exc

    # Find masked output
    masked_files = list(output_dir.glob("*.masked"))
    if len(masked_files) > 1:
        raise MaskingBackendError(
            MaskingBackend.REPEATMASKER,
            f"RepeatMasker produced {len(masked_files)} masked outputs instead of one",
        )
    if not masked_files:
        out_files = list(output_dir.glob("*.out"))
        for out_file in out_files:
            try:
                contents = out_file.read_text()
            except (OSError, UnicodeError) as exc:
                raise MaskingBackendError(
                    MaskingBackend.REPEATMASKER,
                    f"cannot read RepeatMasker output {out_file}: {exc}",
                ) from exc
            if "no repetitive sequences detected" in contents.lower():
                masked_copy = output_dir / f"{input_fasta.name}.masked"
                with Path(input_fasta).open("rb") as source:
                    with atomic_write_context(masked_copy, "wb") as target:
                        shutil.copyfileobj(source, target)
                _require_regular_file(masked_copy, "RepeatMasker no-repeat output")
                logger.info(
                    "RepeatMasker reported no repeats; using unmasked copy at %s",
                    masked_copy,
                )
                return masked_copy

        raise MaskingBackendError(
            MaskingBackend.REPEATMASKER,
            f"RepeatMasker did not produce masked output in {output_dir}",
        )

    return masked_files[0]


def parse_repeatmasker_output(
    rm_fasta: Path,
    input_fasta: Optional[Path] = None,
) -> list[MaskedRegion]:
    """
    Parse RepeatMasker soft-masked FASTA to extract masked regions.

    Soft-masked regions are in lowercase. We identify runs of lowercase
    characters as masked regions.

    Args:
        rm_fasta: Path to RepeatMasker masked FASTA

    Returns:
        List of MaskedRegion objects for identified repeats
    """
    regions = []
    input_sequences = (
        {
            record.id: str(record.seq)
            for record in SeqIO.parse(input_fasta, "fasta")
        }
        if input_fasta is not None
        else {}
    )

    for record in SeqIO.parse(rm_fasta, "fasta"):
        seq = str(record.seq)
        original = input_sequences.get(record.id)
        in_masked = False
        mask_start = 0

        for i, char in enumerate(seq):
            newly_lowercase = char.islower() and (
                original is None or original[i].isupper()
            )
            if newly_lowercase and not in_masked:
                # Start of masked region
                in_masked = True
                mask_start = i
            elif not newly_lowercase and in_masked:
                # End of masked region
                regions.append(
                    MaskedRegion(
                        seq_id=record.id,
                        start=mask_start,
                        end=i,
                        mask_type="repeatmasker",
                    )
                )
                in_masked = False

        # Handle region extending to end
        if in_masked:
            regions.append(
                MaskedRegion(
                    seq_id=record.id,
                    start=mask_start,
                    end=len(seq),
                    mask_type="repeatmasker",
                )
            )

    logger.info(f"Parsed {len(regions)} repeat regions from RepeatMasker output")
    return regions


def merge_overlapping_regions(regions: list[MaskedRegion]) -> list[MaskedRegion]:
    """
    Merge overlapping or adjacent masked regions.

    Args:
        regions: List of MaskedRegion objects

    Returns:
        List of merged MaskedRegion objects
    """
    if not regions:
        return []

    # Group by sequence
    by_seq = {}
    for r in regions:
        if r.seq_id not in by_seq:
            by_seq[r.seq_id] = []
        by_seq[r.seq_id].append(r)

    merged = []
    for seq_id, seq_regions in by_seq.items():
        # Sort by start position
        seq_regions.sort(key=lambda r: r.start)

        current = seq_regions[0]
        for region in seq_regions[1:]:
            if region.start <= current.end:
                # Overlapping or adjacent, merge
                current = MaskedRegion(
                    seq_id=seq_id,
                    start=current.start,
                    end=max(current.end, region.end),
                    mask_type="combined",
                )
            else:
                merged.append(current)
                current = region

        merged.append(current)

    return merged


def apply_mask(
    input_fasta: Path,
    output_fasta: Path,
    regions: list[MaskedRegion],
    mask_char: str = "N",
) -> int:
    """
    Apply masking to a genome FASTA file.

    Converts specified regions to the mask character (default 'N').

    Args:
        input_fasta: Input genome FASTA
        output_fasta: Output masked FASTA
        regions: List of regions to mask
        mask_char: Character to use for masking

    Returns:
        Total number of bases masked
    """
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    # Build region index for efficient lookup
    region_index = {}
    for r in regions:
        if r.seq_id not in region_index:
            region_index[r.seq_id] = []
        region_index[r.seq_id].append((r.start, r.end))

    masked_records = []
    total_masked = 0
    total_bases = 0

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_list = list(str(record.seq))
        total_bases += len(seq_list)

        # Apply masking
        for start, end in region_index.get(record.id, []):
            for i in range(start, min(end, len(seq_list))):
                if seq_list[i].upper() != mask_char.upper():
                    total_masked += 1
                seq_list[i] = mask_char

        masked_record = SeqRecord(
            Seq("".join(seq_list)),
            id=record.id,
            description=record.description,
        )
        masked_records.append(masked_record)

    with atomic_write_context(output_fasta, "w") as handle:
        SeqIO.write(masked_records, handle, "fasta")
    masked_pct = (total_masked / total_bases * 100.0) if total_bases else 0.0
    logger.info(
        "Masked %s bases (%.2f%%) across %s sequences",
        f"{total_masked:,}",
        masked_pct,
        len(masked_records),
    )

    return total_masked


def identify_repeats(
    input_fasta: Path,
    output_dir: Path,
    species: Optional[str] = None,
    library: Optional[Path] = None,
    threads: int = 4,
    skip_repeatmasker: bool = True,
) -> list[MaskedRegion]:
    """
    Identify repeat regions (TRF + RepeatMasker) without masking.

    Returns a list of MaskedRegion objects (unmerged) for downstream features.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    all_regions: list[MaskedRegion] = []

    # Step 1: TRF
    logger.info("Repeat identification: running TRF...")
    trf_dir = output_dir / "trf"
    trf_dat = run_trf(input_fasta, trf_dir)
    trf_regions = parse_trf_output(trf_dat)
    all_regions.extend(trf_regions)
    logger.info("Repeat identification: %d TRF regions", len(trf_regions))

    # Step 2: RepeatMasker (optional)
    if not skip_repeatmasker:
        logger.info("Repeat identification: running RepeatMasker...")
        rm_dir = output_dir / "repeatmasker"
        rm_masked = run_repeatmasker(
            input_fasta,
            rm_dir,
            species=species,
            library=library,
            threads=threads,
        )
        if rm_masked is None:
            raise MaskingBackendError(
                MaskingBackend.REPEATMASKER,
                "RepeatMasker returned no output",
            )
        rm_regions = parse_repeatmasker_output(rm_masked)
        all_regions.extend(rm_regions)
        logger.info("Repeat identification: %d RepeatMasker regions", len(rm_regions))

    return all_regions


def _probe_backend_version(backend: MaskingBackend) -> str:
    """Return a nonempty version string for a backend that just succeeded."""
    executable = "trf" if backend is MaskingBackend.TRF else "RepeatMasker"
    try:
        result = subprocess.run(
            [executable, "-v"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError) as exc:
        raise MaskingBackendError(
            backend,
            f"cannot capture {executable} version: {exc}",
        ) from exc
    version = (result.stdout or result.stderr or "").strip().splitlines()
    if result.returncode != 0 or not version or not version[0].strip():
        raise MaskingBackendError(
            backend,
            f"{executable} -v did not return a version",
        )
    return version[0].strip()


def _fasta_shape(path: Path) -> list[tuple[str, int]]:
    return [(record.id, len(record.seq)) for record in SeqIO.parse(path, "fasta")]


def _validate_fasta_shape(
    input_fasta: Path,
    output_fasta: Path,
    backend: MaskingBackend,
) -> None:
    """Reject missing, reordered, renamed, or resized masking output."""
    expected = _fasta_shape(input_fasta)
    observed = _fasta_shape(output_fasta) if Path(output_fasta).is_file() else []
    if not expected:
        raise MaskingBackendError(backend, "input FASTA contains no records")
    if observed != expected:
        raise MaskingBackendError(
            backend,
            "masked FASTA IDs/order/lengths do not match the input",
        )


def _validate_repeatmasker_output(input_fasta: Path, output_fasta: Path) -> None:
    """Require RepeatMasker to preserve every base except letter case."""
    input_records = list(SeqIO.parse(input_fasta, "fasta"))
    output_records = (
        list(SeqIO.parse(output_fasta, "fasta"))
        if Path(output_fasta).is_file()
        else []
    )
    expected_shape = [(record.id, len(record.seq)) for record in input_records]
    observed_shape = [(record.id, len(record.seq)) for record in output_records]
    if not expected_shape or observed_shape != expected_shape:
        raise MaskingBackendError(
            MaskingBackend.REPEATMASKER,
            "masked FASTA IDs/order/lengths do not match the input",
        )
    for input_record, output_record in zip(input_records, output_records):
        if str(output_record.seq).upper() != str(input_record.seq).upper():
            raise MaskingBackendError(
                MaskingBackend.REPEATMASKER,
                f"masked FASTA changed sequence content for {input_record.id!r}",
            )


def _validate_trf_output(input_fasta: Path, trf_dat: Path) -> None:
    """Require one nonempty regular TRF artifact with exact input header order."""
    try:
        _require_regular_file(trf_dat, "TRF output")
        if Path(trf_dat).stat().st_size == 0:
            raise ValueError(f"TRF output is empty: {trf_dat}")
        observed_ids: list[str] = []
        with Path(trf_dat).open(encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line.startswith("Sequence:"):
                    continue
                header = line.split(maxsplit=1)
                if len(header) != 2 or not header[1].strip():
                    raise ValueError(
                        f"malformed TRF sequence header: {line!r}"
                    )
                observed_ids.append(header[1].split()[0])
    except (OSError, UnicodeError, ValueError) as exc:
        raise MaskingBackendError(
            MaskingBackend.TRF,
            f"invalid TRF output {trf_dat}: {exc}",
        ) from exc
    expected_ids = [record.id for record in SeqIO.parse(input_fasta, "fasta")]
    if observed_ids != expected_ids:
        raise MaskingBackendError(
            MaskingBackend.TRF,
            "TRF Sequence headers do not match input IDs/order",
        )


def _validate_regions(
    input_fasta: Path,
    regions: list[MaskedRegion],
    backend: MaskingBackend,
) -> None:
    """Reject parser regions that cannot address the input sequence exactly."""
    shape = _fasta_shape(input_fasta)
    lengths = {seq_id: length for seq_id, length in shape}
    if len(lengths) != len(shape):
        raise MaskingBackendError(backend, "input FASTA contains duplicate sequence IDs")
    for region in regions:
        length = lengths.get(region.seq_id)
        if length is None:
            raise MaskingBackendError(
                backend,
                f"repeat region references unknown sequence ID {region.seq_id!r}",
            )
        if region.start < 0 or region.end <= region.start or region.end > length:
            raise MaskingBackendError(
                backend,
                f"invalid repeat interval {region.seq_id}:{region.start}-{region.end} "
                f"for sequence length {length}",
            )


def _run_backend_attempt(
    *,
    backend: MaskingBackend,
    input_fasta: Path,
    attempt_dir: Path,
    config: MaskingConfig,
    threads: int,
) -> tuple[list[MaskedRegion], str]:
    """Run one backend in its isolated attempt directory."""
    try:
        if attempt_dir.exists():
            shutil.rmtree(require_strict_child(attempt_dir.parent, attempt_dir))
        attempt_dir.mkdir(parents=True, exist_ok=True)
        if backend is MaskingBackend.TRF:
            trf_dat = run_trf(input_fasta, attempt_dir)
            _validate_trf_output(input_fasta, trf_dat)
            try:
                regions = parse_trf_output(trf_dat)
            except (UnicodeError, ValueError) as exc:
                raise MaskingBackendError(
                    backend,
                    f"malformed trf output {trf_dat}: {exc}",
                ) from exc
        elif backend is MaskingBackend.REPEATMASKER:
            masked = run_repeatmasker(
                input_fasta,
                attempt_dir,
                species=config.repeatmasker_species,
                library=config.repeatmasker_library,
                threads=threads,
            )
            if masked is None:
                raise MaskingBackendError(
                    MaskingBackend.REPEATMASKER,
                    "RepeatMasker returned no output",
                )
            try:
                _validate_repeatmasker_output(input_fasta, masked)
                regions = parse_repeatmasker_output(masked, input_fasta)
            except (UnicodeError, ValueError) as exc:
                raise MaskingBackendError(
                    backend,
                    f"malformed RepeatMasker output {masked}: {exc}",
                ) from exc
        else:
            raise ValueError(f"not an executable masking backend: {backend.value}")
        _validate_regions(input_fasta, regions, backend)
        return regions, _probe_backend_version(backend)
    except MaskingBackendError:
        raise
    except OSError as exc:
        raise MaskingBackendError(
            backend,
            f"cannot read or write backend attempt data: {exc}",
        ) from exc


def _requested_tools(backend: MaskingBackend) -> tuple[MaskingBackend, ...]:
    if backend is MaskingBackend.TRF:
        return (MaskingBackend.TRF,)
    if backend is MaskingBackend.REPEATMASKER:
        return (MaskingBackend.REPEATMASKER,)
    if backend is MaskingBackend.TRF_REPEATMASKER:
        return (MaskingBackend.TRF, MaskingBackend.REPEATMASKER)
    return ()


def mask_genome_pipeline(
    input_fasta: Path,
    output_dir: Path,
    species: Optional[str] = None,
    library: Optional[Path] = None,
    threads: int = 4,
    skip_repeatmasker: bool = True,
    write_mask: bool = True,
    config: Optional[MaskingConfig] = None,
    **legacy_kwargs,
) -> MaskingResult:
    """Execute one explicit masking request and atomically record its outcome."""
    if "apply_mask" in legacy_kwargs:
        write_mask = legacy_kwargs.pop("apply_mask")
    if legacy_kwargs:
        unexpected = ", ".join(sorted(legacy_kwargs))
        raise TypeError(f"unexpected masking arguments: {unexpected}")
    if type(write_mask) is not bool:
        raise TypeError("apply_mask/write_mask must be a boolean")
    input_fasta = Path(input_fasta)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    _status_path(output_dir).unlink(missing_ok=True)
    _validate_input_fasta(input_fasta)

    legacy_mode = config is None
    if config is None:
        backend = (
            MaskingBackend.OFF
            if not write_mask
            else MaskingBackend.TRF
            if skip_repeatmasker or not (species or library)
            else MaskingBackend.TRF_REPEATMASKER
        )
        config = MaskingConfig(
            backend=backend,
            repeatmasker_species=species,
            repeatmasker_library=library,
        )
    errors = config.validate()
    if errors:
        raise ConfigError("Invalid masking configuration: " + "; ".join(errors))
    library_sha256 = _library_sha256(config)
    input_sha256 = _file_sha256(input_fasta)

    if config.backend is MaskingBackend.OFF:
        result = MaskingResult(
            output_path=input_fasta,
            repeat_regions=(),
            requested_backend=MaskingBackend.OFF,
            effective_backend=MaskingBackend.OFF,
            failure_policy=config.failure_policy,
            status="off",
            legacy_adapter=False,
            backend_versions=(),
            masked_bases=0,
            repeatmasker_species=config.repeatmasker_species,
            repeatmasker_library=config.repeatmasker_library,
            repeatmasker_library_sha256=library_sha256,
            configured_fallback_backend=config.fallback_backend,
            fallback_backend=None,
            fallback_reason=None,
            input_sha256=input_sha256,
            output_sha256=input_sha256,
        )
        result = write_masking_status(result, output_dir)
        validate_masking_result(result, expected_input=input_fasta)
        return result

    all_regions: list[MaskedRegion] = []
    backend_versions: dict[str, str] = {}
    failed: Optional[MaskingBackendError] = None
    if legacy_mode:
        all_regions = identify_repeats(
            input_fasta=input_fasta,
            output_dir=output_dir / "attempt_requested",
            species=config.repeatmasker_species,
            library=config.repeatmasker_library,
            threads=threads,
            skip_repeatmasker=config.backend is MaskingBackend.TRF,
        )
        for backend in _requested_tools(config.backend):
            backend_versions[backend.value] = _LEGACY_VERSION
    else:
        for backend in _requested_tools(config.backend):
            try:
                regions, version = _run_backend_attempt(
                    backend=backend,
                    input_fasta=input_fasta,
                    attempt_dir=output_dir / "attempt_requested" / backend.value,
                    config=config,
                    threads=threads,
                )
            except MaskingBackendError as exc:
                backend_versions[backend.value] = "failed"
                failed = exc
                break
            all_regions.extend(regions)
            backend_versions[backend.value] = version

    effective_backend = config.backend
    status = "success"
    fallback_backend: Optional[MaskingBackend] = None
    fallback_reason: Optional[str] = None
    if failed is not None:
        if config.failure_policy is MaskingFailurePolicy.STRICT:
            _write_failed_masking_status(
                input_fasta=input_fasta,
                output_dir=output_dir,
                config=config,
                error=failed,
                backend_versions=backend_versions,
            )
            raise failed
        fallback_backend = config.fallback_backend
        fallback_reason = str(failed)
        status = "fallback"
        if fallback_backend is MaskingBackend.OFF:
            all_regions = []
            effective_backend = MaskingBackend.OFF
        elif fallback_backend is MaskingBackend.TRF:
            try:
                regions, version = _run_backend_attempt(
                    backend=MaskingBackend.TRF,
                    input_fasta=input_fasta,
                    attempt_dir=output_dir / "attempt_fallback" / "trf",
                    config=config,
                    threads=threads,
                )
            except MaskingBackendError as fallback_error:
                backend_versions["fallback_trf"] = "failed"
                _write_failed_masking_status(
                    input_fasta=input_fasta,
                    output_dir=output_dir,
                    config=config,
                    error=fallback_error,
                    backend_versions=backend_versions,
                    selected_fallback_backend=fallback_backend,
                )
                raise fallback_error
            all_regions = list(regions)
            backend_versions["fallback_trf"] = version
            effective_backend = MaskingBackend.TRF
        else:
            raise AssertionError("validated fallback policy has no backend")

    if effective_backend is MaskingBackend.OFF:
        output_fasta = input_fasta
        masked_bases = 0
        output_sha256 = input_sha256
    else:
        _validate_regions(input_fasta, all_regions, effective_backend)
        merged_regions = merge_overlapping_regions(all_regions)
        output_fasta = output_dir / "genome.masked.fna"
        masked_bases = apply_mask(input_fasta, output_fasta, merged_regions)
        _validate_fasta_shape(input_fasta, output_fasta, effective_backend)
        _require_regular_file(output_fasta, "enabled masking output")
        output_sha256 = _file_sha256(output_fasta)

    result = MaskingResult(
        output_path=output_fasta,
        repeat_regions=tuple(all_regions),
        requested_backend=config.backend,
        effective_backend=effective_backend,
        failure_policy=config.failure_policy,
        status=status,
        legacy_adapter=legacy_mode,
        backend_versions=tuple(sorted(backend_versions.items())),
        masked_bases=masked_bases,
        repeatmasker_species=config.repeatmasker_species,
        repeatmasker_library=config.repeatmasker_library,
        repeatmasker_library_sha256=library_sha256,
        configured_fallback_backend=config.fallback_backend,
        fallback_backend=fallback_backend,
        fallback_reason=fallback_reason,
        input_sha256=input_sha256,
        output_sha256=output_sha256,
    )
    result = write_masking_status(result, output_dir)
    validate_masking_result(result, expected_input=input_fasta)
    return result


def quick_mask(input_fasta: Path, output_fasta: Path, threads: int = 4) -> Path:
    """
    Quick masking using only TRF (faster than full pipeline).

    Useful for initial testing or when RepeatMasker libraries are unavailable.

    Args:
        input_fasta: Input genome FASTA
        output_fasta: Output masked FASTA
        threads: Number of threads (used for future parallel TRF)

    Returns:
        Path to masked genome FASTA
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)

        # Run TRF only
        trf_dat = run_trf(input_fasta, output_dir)
        regions = parse_trf_output(trf_dat)

        # Apply masking
        apply_mask(input_fasta, output_fasta, regions)

    return output_fasta
