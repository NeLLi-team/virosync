"""Authenticated manifest contract for ViroSync core resource bundles."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import subprocess
from typing import BinaryIO, TypeAlias

RESOURCE_MANIFEST_NAME = "RESOURCE_MANIFEST.json"

# RESOURCE_MANIFEST.json authenticates these payloads and deliberately excludes
# itself. DB_METADATA.json is generated for each installation and is not part of
# the release bundle.
CORE_RESOURCE_FILES: tuple[str, ...] = (
    "DB_VERSION",
    "DATABASE_README.txt",
    "models/combined.hmm",
    "models/combined.hmm.h3f",
    "models/combined.hmm.h3i",
    "models/combined.hmm.h3m",
    "models/combined.hmm.h3p",
    "models/model_annotations_with_interpro.tsv",
    "models/og_marker_name_map.tsv",
    "marker/marker.faa",
    "marker/marker.dmnd",
    "genomes/combined_proteome.dmnd",
    "taxonomy/labels.tsv",
)

RUNTIME_RESOURCE_FILES: tuple[str, ...] = (
    "DB_VERSION",
    "DATABASE_README.txt",
    "models/combined.hmm",
    "models/model_annotations_with_interpro.tsv",
    "models/og_marker_name_map.tsv",
    "marker/marker.dmnd",
    "genomes/combined_proteome.dmnd",
    "taxonomy/labels.tsv",
)

SOURCE_RESOURCE_FILES: tuple[str, ...] = (
    "models/combined.hmm.h3f",
    "models/combined.hmm.h3i",
    "models/combined.hmm.h3m",
    "models/combined.hmm.h3p",
    "marker/marker.faa",
)

SEMANTIC_COUNT_KEYS: tuple[str, ...] = (
    "hmm_models",
    "hmm_index_files",
    "model_annotations",
    "og_marker_mappings",
    "marker_proteins",
    "marker_diamond_sequences",
    "proteome_diamond_sequences",
    "taxonomy_labels",
)

REQUIRED_SEMANTIC_COUNT_KEYS: tuple[str, ...] = (
    "hmm_models",
    "model_annotations",
    "marker_proteins",
    "marker_diamond_sequences",
    "proteome_diamond_sequences",
    "taxonomy_labels",
)

OPTIONAL_SEMANTIC_COUNT_KEYS: tuple[str, ...] = (
    "hmm_index_files",
    "og_marker_mappings",
)

_CANONICAL_ROLES: Mapping[str, str] = {
    "DB_VERSION": "bundle_version",
    "DATABASE_README.txt": "bundle_documentation",
    "models/combined.hmm": "hmm_models",
    "models/combined.hmm.h3f": "hmm_index",
    "models/combined.hmm.h3i": "hmm_index",
    "models/combined.hmm.h3m": "hmm_index",
    "models/combined.hmm.h3p": "hmm_index",
    "models/model_annotations_with_interpro.tsv": "model_annotations",
    "models/og_marker_name_map.tsv": "og_marker_mappings",
    "marker/marker.faa": "marker_proteins",
    "marker/marker.dmnd": "marker_diamond_database",
    "genomes/combined_proteome.dmnd": "proteome_diamond_database",
    "taxonomy/labels.tsv": "taxonomy_labels",
}

_COMPATIBLE_ROLE_ALIASES: Mapping[str, str] = {
    "models/model_annotations_with_interpro.tsv": "annotation_table",
    "models/og_marker_name_map.tsv": "annotation_table",
    "marker/marker.dmnd": "diamond_database",
    "genomes/combined_proteome.dmnd": "diamond_database",
}

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_VERSION_RE = re.compile(r"^v[0-9]+\.[0-9]+\.[0-9]+$")
_MANIFEST_MAX_BYTES = 1024 * 1024
_HMM_INDEX_FILES = CORE_RESOURCE_FILES[3:7]


class ResourceManifestError(ValueError):
    """A resource manifest or authenticated payload is invalid."""


def _validate_marker_count_parity(semantic_counts: Mapping[str, int]) -> None:
    marker_proteins = semantic_counts.get("marker_proteins")
    marker_diamond = semantic_counts.get("marker_diamond_sequences")
    if (
        marker_proteins is not None
        and marker_diamond is not None
        and marker_proteins != marker_diamond
    ):
        raise ResourceManifestError(
            "marker protein and marker DIAMOND sequence counts differ: "
            f"{marker_proteins} != {marker_diamond}"
        )


@dataclass(frozen=True)
class ManifestFile:
    """One authenticated resource payload."""

    path: str
    size: int
    sha256: str
    role: str


@dataclass(frozen=True)
class ResourceManifest:
    """Parsed core-resource manifest."""

    schema_version: int
    resource_version: str
    version: str
    files: tuple[ManifestFile, ...]
    semantic_counts: dict[str, int]
    manifest_sha256: str
    bundle_kind: str = "legacy"
    runtime_manifest_sha256: str | None = None


@dataclass(frozen=True)
class ResourceValidationResult:
    """Identity returned after a resource tree passes validation."""

    version: str
    manifest_sha256: str
    semantic_counts: dict[str, int]
    files_verified: int
    full: bool


ResourcePayload: TypeAlias = Path | bytes
DiamondSequenceCounter: TypeAlias = Callable[[ResourcePayload], int]


def sha256_file(path: Path) -> str:
    """Return the lowercase SHA-256 digest of *path*."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _validate_sha256(value: object, field: str) -> str:
    if not isinstance(value, str) or _SHA256_RE.fullmatch(value) is None:
        raise ResourceManifestError(f"{field} must be a lowercase 64-character SHA-256")
    return value


def _validate_version(value: object, field: str) -> str:
    if not isinstance(value, str) or _VERSION_RE.fullmatch(value) is None:
        raise ResourceManifestError(f"{field} must have the form vMAJOR.MINOR.PATCH")
    return value


def _validate_relative_path(value: object) -> str:
    if not isinstance(value, str) or not value or "\\" in value:
        raise ResourceManifestError(f"invalid resource path: {value!r}")
    path = PurePosixPath(value)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise ResourceManifestError(
            f"resource path must be relative and normalized: {value!r}"
        )
    if str(path) != value:
        raise ResourceManifestError(f"resource path must be normalized: {value!r}")
    return value


def _object_without_duplicate_keys(
    pairs: list[tuple[str, object]],
) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ResourceManifestError(
                f"duplicate JSON key in resource manifest: {key!r}"
            )
        result[key] = value
    return result


def _parse_manifest(content: bytes, manifest_sha256: str) -> ResourceManifest:
    try:
        document = json.loads(
            content.decode("utf-8"),
            object_pairs_hook=_object_without_duplicate_keys,
        )
    except ResourceManifestError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ResourceManifestError(f"invalid {RESOURCE_MANIFEST_NAME}: {exc}") from exc

    if not isinstance(document, dict):
        raise ResourceManifestError("resource manifest must be a JSON object")
    schema_version = document.get("schema_version")
    if (
        not isinstance(schema_version, int)
        or isinstance(schema_version, bool)
        or schema_version not in {1, 2}
    ):
        raise ResourceManifestError("resource manifest schema_version must be 1 or 2")

    expected_top_level = {
        "schema_version",
        "resource_version",
        "version",
        "files",
        "semantic_counts",
    }
    bundle_kind = "legacy"
    runtime_manifest_sha256: str | None = None
    if schema_version == 2:
        raw_bundle_kind = document.get("bundle_kind")
        if not isinstance(raw_bundle_kind, str) or raw_bundle_kind not in {
            "runtime",
            "source",
        }:
            raise ResourceManifestError(
                "schema_version 2 resource manifest bundle_kind must be 'runtime' or 'source'"
            )
        bundle_kind = raw_bundle_kind
        expected_top_level.add("bundle_kind")
        if bundle_kind == "source":
            expected_top_level.add("runtime_manifest_sha256")
    if set(document) != expected_top_level:
        missing = sorted(expected_top_level - set(document))
        extra = sorted(set(document) - expected_top_level)
        raise ResourceManifestError(
            f"resource manifest fields differ from schema_version {schema_version}; missing={missing}, extra={extra}"
        )
    if bundle_kind == "source":
        runtime_manifest_sha256 = _validate_sha256(
            document["runtime_manifest_sha256"],
            "runtime_manifest_sha256",
        )

    resource_version = _validate_version(
        document["resource_version"], "resource_version"
    )
    version = _validate_version(document["version"], "version")
    if resource_version != version:
        raise ResourceManifestError(
            f"resource manifest version disagreement: {resource_version!r} != {version!r}"
        )

    raw_files = document["files"]
    if not isinstance(raw_files, list):
        raise ResourceManifestError("resource manifest files must be a list")
    files: list[ManifestFile] = []
    seen_paths: set[str] = set()
    for index, raw_file in enumerate(raw_files):
        if not isinstance(raw_file, dict):
            raise ResourceManifestError(
                f"manifest file entry {index} must be an object"
            )
        if set(raw_file) != {"path", "size", "sha256", "role"}:
            raise ResourceManifestError(
                f"manifest file entry {index} must contain path, size, sha256, and role"
            )
        relative = _validate_relative_path(raw_file["path"])
        if relative in seen_paths:
            raise ResourceManifestError(f"duplicate resource manifest path: {relative}")
        seen_paths.add(relative)
        size = raw_file["size"]
        if not isinstance(size, int) or isinstance(size, bool) or size <= 0:
            raise ResourceManifestError(
                f"manifest size for {relative} must be a positive integer"
            )
        digest = _validate_sha256(raw_file["sha256"], f"sha256 for {relative}")
        role = raw_file["role"]
        if not isinstance(role, str) or not role.strip():
            raise ResourceManifestError(
                f"manifest role for {relative} must be non-empty"
            )
        canonical_role = _CANONICAL_ROLES.get(relative)
        compatible_roles = {canonical_role}
        if schema_version == 1:
            compatible_roles.add(_COMPATIBLE_ROLE_ALIASES.get(relative))
        if role not in compatible_roles:
            raise ResourceManifestError(
                f"manifest role mismatch for {relative}: {role!r} is not one of "
                f"{sorted(value for value in compatible_roles if value is not None)!r}"
            )
        files.append(ManifestFile(relative, size, digest, role))

    payload_files = {
        "legacy": CORE_RESOURCE_FILES,
        "runtime": RUNTIME_RESOURCE_FILES,
        "source": SOURCE_RESOURCE_FILES,
    }[bundle_kind]
    required = set(payload_files)
    present = set(seen_paths)
    if present != required or len(files) != len(payload_files):
        raise ResourceManifestError(
            "resource manifest payload set is incomplete or unexpected; "
            f"missing={sorted(required - present)}, unexpected={sorted(present - required)}"
        )
    if tuple(item.path for item in files) != payload_files:
        raise ResourceManifestError(
            "resource manifest files must use canonical payload order"
        )

    raw_counts = document["semantic_counts"]
    if not isinstance(raw_counts, dict):
        raise ResourceManifestError("semantic_counts must be an object")
    required_counts = (
        set(SEMANTIC_COUNT_KEYS)
        if schema_version == 2
        else set(REQUIRED_SEMANTIC_COUNT_KEYS)
    )
    allowed_counts = set(SEMANTIC_COUNT_KEYS)
    present_counts = set(raw_counts)
    if not required_counts.issubset(present_counts) or not present_counts.issubset(
        allowed_counts
    ):
        raise ResourceManifestError(
            f"semantic_counts keys differ from schema_version {schema_version}; "
            f"missing={sorted(required_counts - present_counts)}, "
            f"unexpected={sorted(present_counts - allowed_counts)}"
        )
    semantic_counts: dict[str, int] = {}
    for key in SEMANTIC_COUNT_KEYS:
        if key not in raw_counts:
            continue
        value = raw_counts[key]
        if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
            raise ResourceManifestError(
                f"semantic_counts.{key} must be a positive integer"
            )
        semantic_counts[key] = value
    if "hmm_index_files" in semantic_counts and semantic_counts[
        "hmm_index_files"
    ] != len(_HMM_INDEX_FILES):
        raise ResourceManifestError(
            f"semantic_counts.hmm_index_files must be {len(_HMM_INDEX_FILES)}"
        )
    _validate_marker_count_parity(semantic_counts)

    return ResourceManifest(
        schema_version=schema_version,
        resource_version=resource_version,
        version=version,
        files=tuple(files),
        semantic_counts=semantic_counts,
        manifest_sha256=manifest_sha256,
        bundle_kind=bundle_kind,
        runtime_manifest_sha256=runtime_manifest_sha256,
    )


def load_resource_manifest(
    path_or_root: str | Path,
    expected_version: str | None = None,
    expected_manifest_sha256: str | None = None,
    expected_runtime_manifest_sha256: str | None = None,
) -> ResourceManifest:
    """Load and validate a manifest from a file or resource root."""

    candidate = Path(path_or_root)
    manifest_path = (
        candidate / RESOURCE_MANIFEST_NAME if candidate.is_dir() else candidate
    )
    descriptor = -1
    try:
        descriptor = os.open(
            manifest_path,
            os.O_RDONLY
            | getattr(os, "O_CLOEXEC", 0)
            | getattr(os, "O_NOFOLLOW", 0),
        )
        opened = os.fstat(descriptor)
        current = manifest_path.lstat()
    except OSError as exc:
        raise ResourceManifestError(
            f"missing {RESOURCE_MANIFEST_NAME}: {manifest_path}"
        ) from exc
    try:
        if (
            not stat.S_ISREG(opened.st_mode)
            or opened.st_nlink != 1
            or not stat.S_ISREG(current.st_mode)
            or current.st_nlink != 1
            or (current.st_dev, current.st_ino) != (opened.st_dev, opened.st_ino)
        ):
            raise ResourceManifestError(
                f"{RESOURCE_MANIFEST_NAME} must be a single-link regular file"
            )
        if opened.st_size <= 0 or opened.st_size > _MANIFEST_MAX_BYTES:
            raise ResourceManifestError(
                f"{RESOURCE_MANIFEST_NAME} has invalid size {opened.st_size}"
            )
        with os.fdopen(descriptor, "rb") as handle:
            descriptor = -1
            content = handle.read(_MANIFEST_MAX_BYTES + 1)
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    if len(content) != opened.st_size or len(content) > _MANIFEST_MAX_BYTES:
        raise ResourceManifestError(
            f"{RESOURCE_MANIFEST_NAME} changed while it was being read"
        )
    manifest_sha256 = _sha256_bytes(content)
    if expected_manifest_sha256 is not None:
        expected_digest = _validate_sha256(
            expected_manifest_sha256,
            "expected_manifest_sha256",
        )
        if manifest_sha256 != expected_digest:
            raise ResourceManifestError(
                f"resource manifest SHA-256 mismatch: {manifest_sha256} != {expected_digest}"
            )

    manifest = _parse_manifest(content, manifest_sha256)
    if expected_version is not None:
        expected = _validate_version(expected_version, "expected_version")
        if manifest.version != expected:
            raise ResourceManifestError(
                f"resource version mismatch: {manifest.version!r} != {expected!r}"
            )
    if expected_runtime_manifest_sha256 is not None:
        expected_runtime_digest = _validate_sha256(
            expected_runtime_manifest_sha256,
            "expected_runtime_manifest_sha256",
        )
        if manifest.bundle_kind != "source":
            raise ResourceManifestError(
                "expected_runtime_manifest_sha256 requires a source manifest"
            )
        if manifest.runtime_manifest_sha256 != expected_runtime_digest:
            raise ResourceManifestError(
                f"runtime manifest SHA-256 mismatch: {manifest.runtime_manifest_sha256} != {expected_runtime_digest}"
            )
    return manifest


def _payload_handle(payload: ResourcePayload) -> BinaryIO:
    if isinstance(payload, bytes):
        return io.BytesIO(payload)
    return payload.open("rb")


def _count_prefixed_lines(payload: ResourcePayload, prefix: bytes) -> int:
    with _payload_handle(payload) as handle:
        return sum(1 for line in handle if line.startswith(prefix))


def _count_tsv_rows(payload: ResourcePayload, *, header: bool) -> int:
    with _payload_handle(payload) as handle:
        rows = sum(1 for line in handle if line.strip())
    return max(0, rows - int(header))


def _count_taxonomy_rows(payload: ResourcePayload) -> int:
    rows = 0
    first: bytes | None = None
    with _payload_handle(payload) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if first is None:
                first = stripped
            rows += 1
    if first is not None:
        columns = [part.strip().lower() for part in first.split(b"\t")]
        if (
            len(columns) >= 2
            and columns[0] in {b"genome", b"genome_id", b"id", b"name"}
            and (b"lineage" in columns[1] or b"taxonomy" in columns[1])
        ):
            rows -= 1
    return rows


def _payload_at(
    root: Path,
    relative: str,
    overrides: Mapping[str, ResourcePayload] | None,
) -> ResourcePayload:
    if overrides is not None and relative in overrides:
        return overrides[relative]
    return root / relative


def diamond_sequence_count(
    database: Path,
    *,
    command_runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
) -> int:
    """Return the sequence count reported by ``diamond dbinfo``."""

    database = Path(database)
    try:
        completed = command_runner(
            ["diamond", "dbinfo", "--db", str(database)],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise ResourceManifestError(
            f"diamond dbinfo failed for {database}: {exc}"
        ) from exc
    output = f"{completed.stdout or ''}\n{completed.stderr or ''}"
    match = re.search(r"^\s*Sequences\s+([0-9][0-9,]*)\s*$", output, re.MULTILINE)
    if match is None:
        raise ResourceManifestError(
            f"diamond dbinfo did not report a sequence count for {database}"
        )
    count = int(match.group(1).replace(",", ""))
    if count <= 0:
        raise ResourceManifestError(f"DIAMOND database has no sequences: {database}")
    return count


def compute_semantic_counts(
    root: str | Path,
    *,
    overrides: Mapping[str, ResourcePayload] | None = None,
    diamond_sequence_counter: DiamondSequenceCounter | None = None,
    include_diamond: bool = True,
) -> dict[str, int]:
    """Compute canonical semantic counts from resource payload contents."""

    resource_root = Path(root)
    counts = {
        "hmm_models": _count_prefixed_lines(
            _payload_at(resource_root, "models/combined.hmm", overrides), b"NAME"
        ),
        "hmm_index_files": sum(
            1
            for relative in _HMM_INDEX_FILES
            if _payload_size(_payload_at(resource_root, relative, overrides)) > 0
        ),
        "model_annotations": _count_tsv_rows(
            _payload_at(
                resource_root,
                "models/model_annotations_with_interpro.tsv",
                overrides,
            ),
            header=True,
        ),
        "og_marker_mappings": _count_tsv_rows(
            _payload_at(resource_root, "models/og_marker_name_map.tsv", overrides),
            header=True,
        ),
        "marker_proteins": _count_prefixed_lines(
            _payload_at(resource_root, "marker/marker.faa", overrides), b">"
        ),
        "taxonomy_labels": _count_taxonomy_rows(
            _payload_at(resource_root, "taxonomy/labels.tsv", overrides)
        ),
    }
    if include_diamond:
        if diamond_sequence_counter is None:

            def _count(payload: ResourcePayload) -> int:
                if not isinstance(payload, Path):
                    raise ResourceManifestError(
                        "DIAMOND validation requires a database file path"
                    )
                return diamond_sequence_count(payload)

            diamond_sequence_counter = _count
        counts["marker_diamond_sequences"] = diamond_sequence_counter(
            _payload_at(resource_root, "marker/marker.dmnd", overrides)
        )
        counts["proteome_diamond_sequences"] = diamond_sequence_counter(
            _payload_at(resource_root, "genomes/combined_proteome.dmnd", overrides)
        )
    return {key: counts[key] for key in SEMANTIC_COUNT_KEYS if key in counts}


def _payload_size(payload: ResourcePayload) -> int:
    return len(payload) if isinstance(payload, bytes) else payload.stat().st_size


def _payload_sha256(payload: ResourcePayload) -> str:
    return (
        _sha256_bytes(payload) if isinstance(payload, bytes) else sha256_file(payload)
    )


def _manifest_document(
    version: str,
    files: tuple[ManifestFile, ...],
    semantic_counts: Mapping[str, int],
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "resource_version": version,
        "version": version,
        "files": [
            {
                "path": item.path,
                "size": item.size,
                "sha256": item.sha256,
                "role": item.role,
            }
            for item in files
        ],
        "semantic_counts": {key: semantic_counts[key] for key in SEMANTIC_COUNT_KEYS},
    }


def _manifest_files(
    resource_root: Path,
    payload_files: tuple[str, ...],
    overrides: Mapping[str, ResourcePayload] | None,
) -> tuple[ManifestFile, ...]:
    files: list[ManifestFile] = []
    for relative in payload_files:
        _validate_relative_path(relative)
        payload = _payload_at(resource_root, relative, overrides)
        try:
            if isinstance(payload, Path):
                metadata = payload.lstat()
                if not stat.S_ISREG(metadata.st_mode) or payload.is_symlink():
                    raise ResourceManifestError(
                        f"resource payload must be a regular file: {relative}"
                    )
                if overrides is None or relative not in overrides:
                    payload.resolve(strict=True).relative_to(
                        resource_root.resolve(strict=True)
                    )
            size = _payload_size(payload)
            digest = _payload_sha256(payload)
        except ResourceManifestError:
            raise
        except ValueError as exc:
            raise ResourceManifestError(
                f"resource payload escapes resource root: {relative}"
            ) from exc
        except OSError as exc:
            raise ResourceManifestError(
                f"missing resource payload {relative}: {exc}"
            ) from exc
        if size <= 0:
            raise ResourceManifestError(f"resource payload is empty: {relative}")
        files.append(
            ManifestFile(
                path=relative,
                size=size,
                sha256=digest,
                role=_CANONICAL_ROLES[relative],
            )
        )
    return tuple(files)


def _validated_semantic_counts(
    resource_root: Path,
    overrides: Mapping[str, ResourcePayload] | None,
    diamond_sequence_counter: DiamondSequenceCounter | None,
) -> dict[str, int]:
    counts = compute_semantic_counts(
        resource_root,
        overrides=overrides,
        diamond_sequence_counter=diamond_sequence_counter,
        include_diamond=True,
    )
    for key in SEMANTIC_COUNT_KEYS:
        value = counts.get(key)
        if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
            raise ResourceManifestError(
                f"computed semantic count {key} must be positive"
            )
    if counts["hmm_index_files"] != len(_HMM_INDEX_FILES):
        raise ResourceManifestError("all four HMM index files are required")
    _validate_marker_count_parity(counts)
    return counts


def build_resource_manifest(
    root: str | Path,
    version: str,
    *,
    overrides: Mapping[str, ResourcePayload] | None = None,
    diamond_sequence_counter: DiamondSequenceCounter | None = None,
) -> tuple[ResourceManifest, bytes]:
    """Build deterministic manifest bytes for a complete resource tree."""

    normalized_version = _validate_version(version, "version")
    resource_root = Path(root)
    immutable_files = _manifest_files(
        resource_root,
        CORE_RESOURCE_FILES,
        overrides,
    )
    counts = _validated_semantic_counts(
        resource_root,
        overrides,
        diamond_sequence_counter,
    )

    document = _manifest_document(normalized_version, immutable_files, counts)
    content = (json.dumps(document, indent=2, sort_keys=True) + "\n").encode("utf-8")
    manifest = ResourceManifest(
        schema_version=1,
        resource_version=normalized_version,
        version=normalized_version,
        files=immutable_files,
        semantic_counts=dict(counts),
        manifest_sha256=_sha256_bytes(content),
    )
    return manifest, content


def build_split_resource_manifests(
    root: str | Path,
    version: str,
    *,
    overrides: Mapping[str, ResourcePayload] | None = None,
    diamond_sequence_counter: DiamondSequenceCounter | None = None,
) -> tuple[
    tuple[ResourceManifest, bytes],
    tuple[ResourceManifest, bytes],
]:
    """Build bound schema-v2 runtime and source manifests."""

    normalized_version = _validate_version(version, "version")
    resource_root = Path(root)
    counts = _validated_semantic_counts(
        resource_root,
        overrides,
        diamond_sequence_counter,
    )

    runtime_files = _manifest_files(
        resource_root,
        RUNTIME_RESOURCE_FILES,
        overrides,
    )
    runtime_document = _manifest_document(
        normalized_version,
        runtime_files,
        counts,
    )
    runtime_document["schema_version"] = 2
    runtime_document["bundle_kind"] = "runtime"
    runtime_content = (
        json.dumps(runtime_document, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    runtime_manifest = ResourceManifest(
        schema_version=2,
        resource_version=normalized_version,
        version=normalized_version,
        files=runtime_files,
        semantic_counts=dict(counts),
        manifest_sha256=_sha256_bytes(runtime_content),
        bundle_kind="runtime",
    )

    source_files = _manifest_files(
        resource_root,
        SOURCE_RESOURCE_FILES,
        overrides,
    )
    source_document = _manifest_document(
        normalized_version,
        source_files,
        counts,
    )
    source_document["schema_version"] = 2
    source_document["bundle_kind"] = "source"
    source_document["runtime_manifest_sha256"] = runtime_manifest.manifest_sha256
    source_content = (
        json.dumps(source_document, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    source_manifest = ResourceManifest(
        schema_version=2,
        resource_version=normalized_version,
        version=normalized_version,
        files=source_files,
        semantic_counts=dict(counts),
        manifest_sha256=_sha256_bytes(source_content),
        bundle_kind="source",
        runtime_manifest_sha256=runtime_manifest.manifest_sha256,
    )
    return (
        (runtime_manifest, runtime_content),
        (source_manifest, source_content),
    )


def _regular_payload(root: Path, item: ManifestFile) -> Path:
    candidate = root / item.path
    try:
        metadata = candidate.lstat()
    except OSError as exc:
        raise ResourceManifestError(f"missing resource payload: {item.path}") from exc
    if not stat.S_ISREG(metadata.st_mode) or candidate.is_symlink():
        raise ResourceManifestError(
            f"resource payload must be a regular file: {item.path}"
        )
    try:
        candidate.resolve(strict=True).relative_to(root.resolve(strict=True))
    except (OSError, ValueError) as exc:
        raise ResourceManifestError(
            f"resource payload escapes resource root: {item.path}"
        ) from exc
    if metadata.st_size <= 0:
        raise ResourceManifestError(f"resource payload is empty: {item.path}")
    if metadata.st_size != item.size:
        raise ResourceManifestError(
            f"resource payload size mismatch for {item.path}: {metadata.st_size} != {item.size}"
        )
    return candidate


def _validate_tree_inventory(
    root: Path,
    payload_files: tuple[str, ...] = CORE_RESOURCE_FILES,
) -> None:
    allowed_files = {
        *payload_files,
        RESOURCE_MANIFEST_NAME,
        "DB_METADATA.json",
    }
    allowed_directories = {
        parent.as_posix()
        for relative in allowed_files
        for parent in PurePosixPath(relative).parents
        if parent.as_posix() != "."
    }
    actual_files: set[str] = set()
    actual_directories: set[str] = set()
    for candidate in root.rglob("*"):
        relative = candidate.relative_to(root).as_posix()
        metadata = candidate.lstat()
        if stat.S_ISLNK(metadata.st_mode):
            raise ResourceManifestError(
                f"resource tree contains an unexpected symlink: {relative}"
            )
        if stat.S_ISREG(metadata.st_mode):
            actual_files.add(relative)
        elif stat.S_ISDIR(metadata.st_mode):
            actual_directories.add(relative)
        else:
            raise ResourceManifestError(
                f"resource tree contains an unexpected special file: {relative}"
            )
    unexpected_files = actual_files - allowed_files
    unexpected_directories = actual_directories - allowed_directories
    if unexpected_files or unexpected_directories:
        raise ResourceManifestError(
            "resource tree contains unexpected paths; "
            f"files={sorted(unexpected_files)}, "
            f"directories={sorted(unexpected_directories)}"
        )


def _validate_hmm_annotation_identity(root: Path) -> None:
    hmm_names: list[str] = []
    with (root / "models/combined.hmm").open("rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2 and fields[0] == "NAME":
                hmm_names.append(fields[1])
    annotation_names: list[str] = []
    with (root / "models/model_annotations_with_interpro.tsv").open(
        "rt", encoding="utf-8"
    ) as handle:
        for line_number, line in enumerate(handle):
            if line_number == 0 or not line.strip():
                continue
            annotation_names.append(line.rstrip("\n").split("\t", 1)[0])
    if len(hmm_names) != len(set(hmm_names)):
        raise ResourceManifestError("HMM NAME values must be unique")
    if len(annotation_names) != len(set(annotation_names)):
        raise ResourceManifestError("model annotation identifiers must be unique")
    if set(hmm_names) != set(annotation_names):
        missing = sorted(set(hmm_names) - set(annotation_names))[:5]
        unexpected = sorted(set(annotation_names) - set(hmm_names))[:5]
        raise ResourceManifestError(
            "HMM and model-annotation identifiers disagree; "
            f"missing={missing}, unexpected={unexpected}"
        )


def validate_resource_tree(
    root: str | Path,
    expected_version: str | None = None,
    expected_manifest_sha256: str | None = None,
    expected_runtime_manifest_sha256: str | None = None,
    verify_hashes: bool = True,
    full: bool = False,
    command_runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
) -> ResourceValidationResult:
    """Validate an extracted resource tree and return its immutable identity.

    Fast validation performs no child-process calls. ``full=True`` additionally
    asks DIAMOND to inspect both databases and compares their sequence counts.
    """

    resource_root = Path(root)
    if not resource_root.is_dir():
        raise ResourceManifestError(
            f"resource root is not a directory: {resource_root}"
        )
    manifest = load_resource_manifest(
        resource_root,
        expected_version=expected_version,
        expected_manifest_sha256=expected_manifest_sha256,
        expected_runtime_manifest_sha256=expected_runtime_manifest_sha256,
    )
    manifest_payload_files = tuple(item.path for item in manifest.files)
    _validate_tree_inventory(resource_root, manifest_payload_files)

    payload_paths: dict[str, Path] = {}
    for item in manifest.files:
        payload = _regular_payload(resource_root, item)
        payload_paths[item.path] = payload
        if verify_hashes:
            actual_digest = sha256_file(payload)
            if actual_digest != item.sha256:
                raise ResourceManifestError(
                    f"resource payload SHA-256 mismatch for {item.path}: "
                    f"{actual_digest} != {item.sha256}"
                )

    if "DB_VERSION" in payload_paths:
        try:
            version_file = (
                payload_paths["DB_VERSION"].read_text(encoding="utf-8").strip()
            )
        except (OSError, UnicodeDecodeError) as exc:
            raise ResourceManifestError(f"invalid DB_VERSION: {exc}") from exc
        if version_file != manifest.version:
            raise ResourceManifestError(
                f"DB_VERSION mismatch: {version_file!r} != {manifest.version!r}"
            )

    if verify_hashes or full:
        if manifest.schema_version == 1:
            actual_counts = compute_semantic_counts(
                resource_root,
                include_diamond=False,
            )
        else:
            actual_counts = {}
            if "models/combined.hmm" in payload_paths:
                actual_counts["hmm_models"] = _count_prefixed_lines(
                    payload_paths["models/combined.hmm"],
                    b"NAME",
                )
            hmm_indices = [
                payload_paths[relative]
                for relative in _HMM_INDEX_FILES
                if relative in payload_paths
            ]
            if hmm_indices:
                actual_counts["hmm_index_files"] = sum(
                    1 for payload in hmm_indices if payload.stat().st_size > 0
                )
            if "models/model_annotations_with_interpro.tsv" in payload_paths:
                actual_counts["model_annotations"] = _count_tsv_rows(
                    payload_paths["models/model_annotations_with_interpro.tsv"],
                    header=True,
                )
            if "models/og_marker_name_map.tsv" in payload_paths:
                actual_counts["og_marker_mappings"] = _count_tsv_rows(
                    payload_paths["models/og_marker_name_map.tsv"],
                    header=True,
                )
            if "marker/marker.faa" in payload_paths:
                actual_counts["marker_proteins"] = _count_prefixed_lines(
                    payload_paths["marker/marker.faa"],
                    b">",
                )
            if "taxonomy/labels.tsv" in payload_paths:
                actual_counts["taxonomy_labels"] = _count_taxonomy_rows(
                    payload_paths["taxonomy/labels.tsv"]
                )
        if "hmm_index_files" in actual_counts and actual_counts[
            "hmm_index_files"
        ] != len(_HMM_INDEX_FILES):
            raise ResourceManifestError(
                f"resource tree must contain {len(_HMM_INDEX_FILES)} non-empty HMM index files"
            )
        for key, actual_count in actual_counts.items():
            if actual_count <= 0:
                raise ResourceManifestError(f"semantic count {key} must be positive")
            if key not in manifest.semantic_counts:
                continue
            expected = manifest.semantic_counts[key]
            if actual_count != expected:
                raise ResourceManifestError(
                    f"semantic count mismatch for {key}: {actual_count} != {expected}"
                )
        if {
            "models/combined.hmm",
            "models/model_annotations_with_interpro.tsv",
        }.issubset(payload_paths):
            _validate_hmm_annotation_identity(resource_root)

    if full:
        diamond_counts = {}
        if "marker/marker.dmnd" in payload_paths:
            diamond_counts["marker_diamond_sequences"] = diamond_sequence_count(
                payload_paths["marker/marker.dmnd"],
                command_runner=command_runner,
            )
        if "genomes/combined_proteome.dmnd" in payload_paths:
            diamond_counts["proteome_diamond_sequences"] = diamond_sequence_count(
                payload_paths["genomes/combined_proteome.dmnd"],
                command_runner=command_runner,
            )
        for key, actual_count in diamond_counts.items():
            expected = manifest.semantic_counts[key]
            if actual_count != expected:
                raise ResourceManifestError(
                    f"semantic count mismatch for {key}: {actual_count} != {expected}"
                )

    return ResourceValidationResult(
        version=manifest.version,
        manifest_sha256=manifest.manifest_sha256,
        semantic_counts=dict(manifest.semantic_counts),
        files_verified=len(manifest.files),
        full=full,
    )
