"""
Database management utilities for ViroSync.
Handles downloading and validating the reference database bundle.
"""

from __future__ import annotations

import json
import logging
import math
import os
import csv
import re
import shutil
import subprocess
import tempfile
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from virosync.utils.resource_installer import (
    ResourceInstallError,
    ResourceSource,
    active_installed_candidate,
    copy_or_download_archive,
    install_core_resources,
    recover_pending_install,
    safe_extract_archive,
    safe_extract_optional_archive,
    sha256_file,
    sibling_install_lock,
    verified_install_receipt,
)
from virosync.utils.resource_manifest import (
    LEGACY_RUNTIME_RESOURCE_FILES,
    RUNTIME_RESOURCE_FILES,
    ResourceManifestError,
    ResourceValidationResult,
    load_resource_manifest,
    validate_resource_tree,
)

logger = logging.getLogger(__name__)

# ``resource_installer`` enforces fcntl.flock(..., fcntl.LOCK_EX) on the stable
# sibling ``.resource-install.lock`` before any recovery, download, or activation.
TMVEC_EMBEDDING_WIDTH = 512
TMVEC_MANIFEST_NAME = "TMVEC_MANIFEST.json"
TMVEC2_BASE_MODEL_ID = "asalam91/lobster_24M"
TMVEC2_BASE_MODEL_REVISION = "9c36ae05d277e312ac319cbc41b5759472f5bd90"
TMVEC2_BASE_WEIGHT_SHA256 = (
    "d80ed1022349db63a51a3ee2ea0ea5f71aa78e36f4a5dd4977ae3da49e6b9aa6"
)
TMVEC2_HEAD_MODEL_ID = "scikit-bio/TMVec-2"
TMVEC2_HEAD_MODEL_REVISION = "91fbaaefbacd72ff6bc2f2126e8a0c165b2a9d92"
TMVEC2_HEAD_WEIGHT_SHA256 = (
    "7739dc359b62712061ad79f01269b37d96eae0a9e4c810c4d5fbef58eab85302"
)
TMVEC2_ARCHITECTURE = {
    "base_embedding_dim": 408,
    "output_dim": TMVEC_EMBEDDING_WIDTH,
    "nhead": 8,
    "num_layers": 4,
    "dim_feedforward": 2048,
    "transformer_activation": "gelu",
    "projection_hidden_dim": 1024,
    "projection_activation": "relu",
    "dropout": 0.2,
    "max_sequence_length": 512,
}
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


class ViroSyncDatabaseManager:
    """Manages ViroSync reference databases."""

    REQUIRED_FILES = [
        "DB_VERSION",
        "DATABASE_README.txt",
        "models/model_annotations_with_interpro.tsv",
        "models/og_marker_name_map.tsv",
        "marker/marker.dmnd",
        "genomes/combined_proteome.dmnd",
        "taxonomy/labels.tsv",
    ]

    TMVEC_REQUIRED_FILES = {
        "bfvd": [
            "bfvd/bfvd_embeddings.npy",
            "bfvd/bfvd_annotations.tsv",
        ],
    }
    INTERPROSCAN_REQUIRED_FILES = ["interproscan.sh"]

    DATABASE_SOURCES = [
        {
            "version": "v1.0.7",
            "source": "https://dl.newlineages.com/virosync/resources_v1_0_7_runtime.tar.gz",
            "filename": "resources_v1_0_7_runtime.tar.gz",
            "archive_size_bytes": 5_877_324_818,
            "payload_size_bytes": 13_137_477_318,
            "archive_sha256": (
                "57daed0b39bf2bc4c4f84ec3b612c6034a3d26ea38e7ec5fba4f4469da36e9a2"
            ),
            "manifest_sha256": (
                "f3aeed77045f4728207c6997f5986ed155056e2b4b2a297574d57686982a18b3"
            ),
        },
    ]
    DATABASE_VERSION = DATABASE_SOURCES[0]["version"]

    @classmethod
    def _project_root(cls) -> Path:
        return Path(__file__).resolve().parents[3]

    @staticmethod
    def normalize_path(path_value: str | Path) -> Path:
        """Return an absolute path without resolving symlinks."""
        path = Path(path_value).expanduser()
        if not path.is_absolute():
            path = Path.cwd() / path
        return path

    @classmethod
    def default_database_path(cls) -> Path:
        return cls._project_root() / "resources" / "virosync"

    @classmethod
    def default_tmvec_path(cls, database_path: Optional[Path] = None) -> Path:
        if database_path:
            core_path = Path(database_path)
        else:
            core_path = cls.default_database_path()
        return core_path.parent / f"{core_path.name}-optional" / "tmvec"

    @classmethod
    def default_interproscan_path(cls, database_path: Optional[Path] = None) -> Path:
        if database_path:
            core_path = Path(database_path)
        else:
            core_path = cls.default_database_path()
        return core_path.parent / f"{core_path.name}-optional" / "interproscan"

    @classmethod
    def get_database_version(cls, database_path: Path) -> str:
        version_file = database_path / "DB_VERSION"
        if version_file.exists():
            content = version_file.read_text().strip()
            if not content:
                return "unknown"
            return next((line.strip() for line in content.splitlines() if line.strip()), "unknown")
        return "unknown"

    @classmethod
    def _write_database_metadata(
        cls,
        target_path: Path,
        component: str,
        version: str,
        source: str,
        extra: Optional[dict] = None,
    ) -> None:
        payload = {
            "component": component,
            "version": version,
            "installed_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "source": source,
        }
        if extra:
            payload.update(extra)
        with open(target_path / "DB_METADATA.json", "w") as handle:
            json.dump(payload, handle, indent=2)

    @staticmethod
    def _certifi_ca_bundle() -> Optional[str]:
        try:
            import certifi
        except ImportError:
            return None
        certifi_path = Path(certifi.where())
        if certifi_path.is_file():
            return str(certifi_path)
        return None

    @classmethod
    def _record_for_source(cls, source: str) -> Optional[dict]:
        normalized = source.replace("file://", "", 1) if source.startswith("file://") else source
        for candidate in cls.DATABASE_SOURCES:
            candidate_source = str(candidate["source"])
            candidate_normalized = (
                candidate_source.replace("file://", "", 1)
                if candidate_source.startswith("file://")
                else candidate_source
            )
            if normalized == candidate_normalized:
                return candidate
        return None

    @classmethod
    def _resolve_core_source(
        cls,
        source: Optional[str],
        version: Optional[str],
        archive_sha256: Optional[str],
        manifest_sha256: Optional[str],
    ) -> ResourceSource:
        if source is None:
            if not cls.DATABASE_SOURCES:
                raise ResourceInstallError("No authenticated core-resource source is configured")
            record = cls.DATABASE_SOURCES[0]
            selected_source = str(record["source"])
        else:
            selected_source = source
            record = cls._record_for_source(source)

        record_version = str(record["version"]) if record is not None else None
        record_archive_sha = (
            record.get("archive_sha256") or record.get("sha256")
            if record is not None
            else None
        )
        record_manifest_sha = (
            record.get("manifest_sha256") or record.get("resource_manifest_sha256")
            if record is not None
            else None
        )
        selected_version = version or record_version
        selected_archive_sha = archive_sha256 or record_archive_sha
        selected_manifest_sha = manifest_sha256 or record_manifest_sha

        if record is not None:
            mismatches = []
            if version is not None and version != record_version:
                mismatches.append(f"version {version!r} != {record_version!r}")
            if archive_sha256 is not None and archive_sha256 != record_archive_sha:
                mismatches.append("archive SHA-256 differs from the pinned source record")
            if manifest_sha256 is not None and manifest_sha256 != record_manifest_sha:
                mismatches.append("manifest SHA-256 differs from the pinned source record")
            if mismatches:
                raise ResourceInstallError(
                    "Configured core-resource identity conflicts with its source record: "
                    + "; ".join(mismatches)
                )

        if not selected_version or not selected_archive_sha or not selected_manifest_sha:
            raise ResourceInstallError(
                "Core-resource source is unpinned; source, version, archive SHA-256, "
                "and manifest SHA-256 are all required"
            )
        filename = (
            str(record.get("filename"))
            if record is not None and record.get("filename")
            else Path(selected_source.replace("file://", "", 1)).name
        )
        return ResourceSource(
            version=selected_version,
            source=selected_source,
            filename=filename or "resources.tar.gz",
            archive_sha256=str(selected_archive_sha),
            manifest_sha256=str(selected_manifest_sha),
        )

    @classmethod
    def _copy_or_download_archive(
        cls,
        source: str,
        archive_path: Path,
        *,
        progress_callback=None,
    ) -> None:
        copy_or_download_archive(
            source,
            archive_path,
            ca_bundle=cls._certifi_ca_bundle(),
            command_runner=subprocess.run,
            progress_callback=progress_callback,
        )

    @staticmethod
    def _safe_extract_archive(archive_path: Path, target_dir: Path) -> None:
        safe_extract_archive(archive_path, target_dir)

    @staticmethod
    def _extract_optional_archive(archive_path: Path, target_dir: Path) -> None:
        """Extract a legacy optional bundle that has no core manifest contract."""
        safe_extract_optional_archive(archive_path, target_dir)

    @classmethod
    def _install_archive(
        cls,
        target_dir: Path,
        source: str,
        filename: Optional[str] = None,
        archive_sha256: Optional[str] = None,
        force: bool = False,
        progress_callback=None,
        payload_preparer=None,
        payload_validator=None,
    ) -> None:
        target_dir.parent.mkdir(parents=True, exist_ok=True)
        if archive_sha256 is not None and _SHA256_RE.fullmatch(archive_sha256) is None:
            raise ResourceInstallError(
                "optional resource archive SHA-256 must be a lowercase 64-character digest"
            )
        archive_name = filename or Path(source).name or "resources.tar.gz"
        with tempfile.TemporaryDirectory(
            prefix=f".{target_dir.name}.optional-stage-",
            dir=target_dir.parent,
        ) as stage_name:
            stage_root = Path(stage_name)
            archive_path = stage_root / Path(archive_name).name
            payload_path = stage_root / "payload"
            cls._copy_or_download_archive(
                source,
                archive_path,
                progress_callback=(
                    (
                        lambda percent: progress_callback(
                            percent * 0.70,
                            f"downloading ({int(percent)}%)",
                        )
                    )
                    if progress_callback is not None
                    else None
                ),
            )
            if archive_sha256 is not None:
                actual_sha256 = sha256_file(archive_path)
                if actual_sha256 != archive_sha256:
                    raise ResourceInstallError(
                        "optional resource archive checksum mismatch: "
                        f"expected {archive_sha256}, found {actual_sha256}"
                    )
            if progress_callback is not None:
                progress_callback(70, "extracting")
            logger.info("Extracting archive: %s", archive_path)
            cls._extract_optional_archive(archive_path, payload_path)
            if payload_preparer is not None:
                payload_preparer(payload_path)
            if payload_validator is not None:
                payload_validator(payload_path)
            if progress_callback is not None:
                progress_callback(90, "activating")
            if force:
                if target_dir.is_symlink() or (
                    target_dir.exists() and not target_dir.is_dir()
                ):
                    target_dir.unlink()
                elif target_dir.exists():
                    shutil.rmtree(target_dir)
                os.rename(payload_path, target_dir)
            else:
                if target_dir.is_symlink() or (
                    target_dir.exists() and not target_dir.is_dir()
                ):
                    raise ResourceInstallError(
                        f"Optional resource target is not a directory: {target_dir}"
                    )
                shutil.copytree(payload_path, target_dir, dirs_exist_ok=True)
            if progress_callback is not None:
                progress_callback(100, "ready")

    @classmethod
    def _validate_core_tree(
        cls,
        database_path: Path,
        *,
        expected_version: Optional[str],
        expected_manifest_sha256: Optional[str],
        verify_hashes: bool,
        full: bool,
        semantic_runner=None,
    ) -> ResourceValidationResult:
        if not full or semantic_runner is None:
            return validate_resource_tree(
                database_path,
                expected_version=expected_version,
                expected_manifest_sha256=expected_manifest_sha256,
                verify_hashes=verify_hashes,
                full=full,
            )

        result = validate_resource_tree(
            database_path,
            expected_version=expected_version,
            expected_manifest_sha256=expected_manifest_sha256,
            verify_hashes=verify_hashes,
            full=False,
        )
        manifest = load_resource_manifest(
            database_path,
            expected_version=expected_version,
            expected_manifest_sha256=expected_manifest_sha256,
        )
        actual_counts = semantic_runner(database_path)
        if not isinstance(actual_counts, dict):
            actual_counts = dict(actual_counts)
        unavailable_runtime_counts = (
            {"hmm_index_files", "marker_proteins"}
            if manifest.bundle_kind == "runtime"
            else set()
        )
        for key, expected in manifest.semantic_counts.items():
            if key in unavailable_runtime_counts:
                continue
            actual = actual_counts.get(key)
            if actual != expected:
                raise ResourceManifestError(
                    f"semantic count mismatch for {key}: {actual!r} != {expected}"
                )
        return replace(result, full=True)

    @classmethod
    def verify_database(
        cls,
        database_path: str | Path,
        *,
        expected_version: Optional[str] = None,
        manifest_sha256: Optional[str] = None,
        full: bool = False,
        semantic_runner=None,
    ) -> ResourceValidationResult:
        """Verify the stable core-resource path and return its manifest identity."""
        return cls._validate_core_tree(
            cls.normalize_path(database_path),
            expected_version=expected_version,
            expected_manifest_sha256=manifest_sha256,
            verify_hashes=full,
            full=full,
            semantic_runner=semantic_runner,
        )

    @classmethod
    def setup_database(
        cls,
        database_path: Optional[str] = None,
        source: Optional[str] = None,
        version: Optional[str] = None,
        archive_sha256: Optional[str] = None,
        manifest_sha256: Optional[str] = None,
        force: bool = False,
        full: bool = True,
        semantic_runner=None,
        fault_injector=None,
        progress_callback=None,
    ) -> Path:
        """Install one authenticated bundle behind the stable database pointer."""
        if database_path:
            db_path = cls.normalize_path(database_path)
        else:
            db_path = cls.default_database_path()

        logger.info("Setting up ViroSync database at: %s", db_path)
        implicit_source = source is None
        selected = cls._resolve_core_source(
            source,
            version,
            archive_sha256,
            manifest_sha256,
        )

        def _verify_tree(path: Path, **kwargs):
            runner = kwargs.pop("command_runner", semantic_runner)
            return cls._validate_core_tree(
                path,
                semantic_runner=runner,
                **kwargs,
            )

        def _copy_archive(source_value: str, archive_path: Path) -> None:
            cls._copy_or_download_archive(
                source_value,
                archive_path,
                progress_callback=(
                    (
                        lambda percent: progress_callback(
                            5 + percent * 0.30,
                            f"downloading core resources ({int(percent)}%)",
                        )
                    )
                    if progress_callback is not None
                    else None
                ),
            )

        installed = install_core_resources(
            db_path,
            selected,
            copy_archive=_copy_archive,
            verify_tree=_verify_tree,
            required_files=list(RUNTIME_RESOURCE_FILES),
            full=full,
            reuse_existing=not force,
            reject_invalid_existing=implicit_source and not force,
            semantic_runner=semantic_runner,
            fault_injector=fault_injector,
            progress_callback=progress_callback,
        )
        logger.info("Database setup complete")
        return installed

    @classmethod
    def hmm_db_path(cls, db_path: Path) -> Path:
        return db_path / "models" / "combined.hmm"

    @classmethod
    def required_files_for_path(cls, db_path: Path) -> list[str]:
        try:
            manifest = load_resource_manifest(db_path)
        except (OSError, ResourceManifestError):
            manifest = None
        if manifest is not None and manifest.schema_version == 2:
            return list(RUNTIME_RESOURCE_FILES)
        return list(LEGACY_RUNTIME_RESOURCE_FILES)

    @classmethod
    def _check_missing_files(cls, db_path: Path) -> list[str]:
        missing = []
        for rel_path in cls.required_files_for_path(db_path):
            candidate = db_path / rel_path
            if not candidate.is_file() or candidate.stat().st_size <= 0:
                missing.append(rel_path)
        return missing

    @classmethod
    def setup_optional_archive(
        cls,
        name: str,
        target_path: Path,
        source: Optional[str],
        required_files: list[str],
        version: str = "unknown",
        archive_sha256: Optional[str] = None,
        force: bool = False,
        progress_callback=None,
    ) -> bool:
        def _validate_payload(path: Path) -> None:
            if name == "tmvec":
                cls.load_tmvec_manifest(path, verify_hashes=True)
                return
            missing = [rel for rel in required_files if not (path / rel).is_file()]
            if missing:
                raise ResourceInstallError(
                    f"{name} installation incomplete; missing files: {missing}"
                )
            if name == "interproscan" and not os.access(
                path / "interproscan.sh", os.X_OK
            ):
                raise ResourceInstallError(
                    "InterProScan installation incomplete; interproscan.sh "
                    "is not executable"
                )

        if (
            name in {"tmvec", "interproscan"}
            and source is not None
            and archive_sha256 is None
        ):
            logger.error("%s archive source requires archive_sha256", name)
            return False

        if source is None:
            try:
                if name == "tmvec":
                    cls._download_tmvec_models(target_path)
                _validate_payload(target_path)
            except (OSError, ResourceInstallError, ValueError) as exc:
                logger.warning(
                    "%s not installed and no source URL/path provided: %s",
                    name,
                    exc,
                )
                return False
            installed_version = (
                cls.load_tmvec_manifest(target_path)["bundle_version"]
                if name == "tmvec"
                else version
            )
            metadata_path = target_path / "DB_METADATA.json"
            if not metadata_path.exists():
                cls._write_database_metadata(
                    target_path=target_path,
                    component=name,
                    version=installed_version,
                    source="preexisting",
                    extra={
                        "required_files": required_files,
                        "note": "Metadata created from existing installation",
                    },
                )
            return True

        try:
            cls._install_archive(
                target_path,
                source,
                archive_sha256=archive_sha256,
                force=force,
                progress_callback=progress_callback,
                payload_preparer=(
                    cls._download_tmvec_models if name == "tmvec" else None
                ),
                payload_validator=_validate_payload,
            )
        except Exception as exc:
            logger.error("%s installation failed: %s", name, exc)
            return False

        installed_version = (
            cls.load_tmvec_manifest(target_path)["bundle_version"]
            if name == "tmvec"
            else version
        )
        cls._write_database_metadata(
            target_path=target_path,
            component=name,
            version=installed_version,
            source=source,
            extra={
                "required_files": required_files,
                "archive_sha256": archive_sha256,
            },
        )
        return True

    @staticmethod
    def _tmvec_npy_header(path: Path) -> tuple[tuple[int, ...], np.dtype]:
        with path.open("rb") as handle:
            version = np.lib.format.read_magic(handle)
            shape, _fortran_order, dtype = np.lib.format._read_array_header(
                handle,
                version,
            )
            payload_offset = handle.tell()
        if dtype.hasobject:
            raise ResourceInstallError(f"TMVec NPY asset must not contain objects: {path}")
        expected_bytes = math.prod(shape) * dtype.itemsize
        payload_bytes = path.stat().st_size - payload_offset
        if payload_bytes != expected_bytes:
            raise ResourceInstallError(
                f"TMVec NPY payload size mismatch for {path}: "
                f"expected {expected_bytes}, found {max(0, payload_bytes)}"
            )
        return tuple(shape), dtype

    @staticmethod
    def _tmvec_member(root: Path, value: object, label: str) -> Path:
        if not isinstance(value, str) or not value.strip():
            raise ResourceInstallError(f"{label} must be a non-empty relative path")
        relative = Path(value)
        if relative.is_absolute() or ".." in relative.parts:
            raise ResourceInstallError(f"{label} must stay inside the TMVec resource root")
        resolved_root = root.resolve()
        resolved = (root / relative).resolve()
        if not resolved.is_relative_to(resolved_root):
            raise ResourceInstallError(f"{label} resolves outside the TMVec resource root")
        if not resolved.is_file() or resolved.stat().st_size <= 0:
            raise ResourceInstallError(f"{label} is missing or empty: {resolved}")
        return resolved

    @classmethod
    def _download_tmvec_models(cls, root: Path) -> None:
        """Fetch hash-pinned model files declared by an authenticated DB bundle."""
        manifest_path = root / TMVEC_MANIFEST_NAME
        try:
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise ResourceInstallError(
                f"cannot read {manifest_path} before model download: {exc}"
            ) from exc
        model = payload.get("model") if isinstance(payload, dict) else None
        if not isinstance(model, dict) or model.get("family") != "tmvec2":
            raise ResourceInstallError("TMVec bundle does not declare model.family=tmvec2")

        contracts = {
            "base": (
                TMVEC2_BASE_MODEL_ID,
                TMVEC2_BASE_MODEL_REVISION,
                {
                    "config.json",
                    "pytorch_model.bin",
                    "special_tokens_map.json",
                    "tokenizer_config.json",
                    "vocab.txt",
                },
            ),
            "head": (
                TMVEC2_HEAD_MODEL_ID,
                TMVEC2_HEAD_MODEL_REVISION,
                {"params.json", "tmvec-2.ckpt"},
            ),
        }
        resolved_root = root.resolve()
        for label, (model_id, revision, required_names) in contracts.items():
            section = model.get(label)
            if (
                not isinstance(section, dict)
                or section.get("id") != model_id
                or section.get("revision") != revision
            ):
                raise ResourceInstallError(
                    f"TMVec bundle must pin {model_id}@{revision} before model download"
                )
            files = section.get("files")
            if not isinstance(files, list) or not files:
                raise ResourceInstallError(f"model.{label}.files must be a non-empty list")
            declared_names = [
                Path(item.get("path", "")).name
                for item in files
                if isinstance(item, dict) and isinstance(item.get("path"), str)
            ]
            if len(declared_names) != len(required_names) or set(
                declared_names
            ) != required_names:
                raise ResourceInstallError(
                    f"model.{label}.files must contain exactly: "
                    + ", ".join(sorted(required_names))
                )
            for index, item in enumerate(files):
                if not isinstance(item, dict):
                    raise ResourceInstallError(
                        f"model.{label}.files[{index}] must be a mapping"
                    )
                relative_text = item.get("path")
                sha256 = item.get("sha256")
                if not isinstance(relative_text, str) or not relative_text:
                    raise ResourceInstallError(
                        f"model.{label}.files[{index}].path must be a relative path"
                    )
                relative = Path(relative_text)
                expected_dir = "lobster_24M" if label == "base" else "tmvec-2"
                if (
                    relative.is_absolute()
                    or ".." in relative.parts
                    or relative.parent != Path("models") / expected_dir
                ):
                    raise ResourceInstallError(
                        f"model.{label} file has an invalid bundle path: {relative_text}"
                    )
                target = (root / relative).resolve()
                if not target.is_relative_to(resolved_root):
                    raise ResourceInstallError(
                        f"model.{label} file resolves outside the bundle: {relative_text}"
                    )
                if not isinstance(sha256, str) or _SHA256_RE.fullmatch(sha256) is None:
                    raise ResourceInstallError(
                        f"model.{label}.files[{index}].sha256 must be a lowercase SHA-256"
                    )
                if target.is_file() and sha256_file(target) == sha256:
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                source = (
                    f"https://huggingface.co/{model_id}/resolve/{revision}/"
                    f"{relative.name}?download=true"
                )
                cls._copy_or_download_archive(source, target)
                if sha256_file(target) != sha256:
                    raise ResourceInstallError(
                        f"downloaded TMVec model checksum mismatch: {target}"
                    )

    @classmethod
    def load_tmvec_manifest(
        cls,
        tmvec_root: str | Path,
        *,
        verify_hashes: bool = False,
        databases: Optional[list[str]] = None,
    ) -> dict:
        """Load and validate one model-bound TMVec2 resource manifest."""
        root = Path(tmvec_root).expanduser()
        manifest_path = root / TMVEC_MANIFEST_NAME
        try:
            if not manifest_path.is_file() or manifest_path.stat().st_size <= 0:
                raise ResourceInstallError(f"missing {manifest_path}")
            if manifest_path.stat().st_size > 1024 * 1024:
                raise ResourceInstallError(f"{manifest_path} exceeds 1 MiB")
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        except ResourceInstallError:
            raise
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise ResourceInstallError(
                f"cannot read {manifest_path}: {exc}"
            ) from exc
        if not isinstance(payload, dict):
            raise ResourceInstallError(f"{manifest_path} must contain a JSON object")
        if payload.get("schema_version") != 1:
            raise ResourceInstallError("TMVec manifest schema_version must be 1")
        bundle_version = payload.get("bundle_version")
        if not isinstance(bundle_version, str) or re.fullmatch(
            r"v[0-9]+\.[0-9]+\.[0-9]+",
            bundle_version,
        ) is None:
            raise ResourceInstallError(
                "TMVec manifest bundle_version must have the form vMAJOR.MINOR.PATCH"
            )

        model = payload.get("model")
        if not isinstance(model, dict) or model.get("family") != "tmvec2":
            raise ResourceInstallError("TMVec manifest model.family must be 'tmvec2'")
        architecture = model.get("architecture")
        if not isinstance(architecture, dict):
            raise ResourceInstallError("TMVec manifest model.architecture must be a mapping")
        for name, expected in TMVEC2_ARCHITECTURE.items():
            actual = architecture.get(name)
            if type(actual) is not type(expected) or actual != expected:
                raise ResourceInstallError(
                    f"TMVec2 architecture mismatch for {name}: "
                    f"expected {expected!r}, found {actual!r}"
                )
        if set(architecture) != set(TMVEC2_ARCHITECTURE):
            unexpected = sorted(set(architecture) - set(TMVEC2_ARCHITECTURE))
            missing = sorted(set(TMVEC2_ARCHITECTURE) - set(architecture))
            raise ResourceInstallError(
                "TMVec2 architecture keys differ from the release contract: "
                f"missing={missing}, unexpected={unexpected}"
            )

        model_contracts = {
            "base": (
                TMVEC2_BASE_MODEL_ID,
                TMVEC2_BASE_MODEL_REVISION,
                "models/lobster_24M/pytorch_model.bin",
                TMVEC2_BASE_WEIGHT_SHA256,
                {
                    "config.json",
                    "pytorch_model.bin",
                    "special_tokens_map.json",
                    "tokenizer_config.json",
                    "vocab.txt",
                },
            ),
            "head": (
                TMVEC2_HEAD_MODEL_ID,
                TMVEC2_HEAD_MODEL_REVISION,
                "models/tmvec-2/tmvec-2.ckpt",
                TMVEC2_HEAD_WEIGHT_SHA256,
                {"params.json", "tmvec-2.ckpt"},
            ),
        }
        manifest_members: set[str] = set()
        for label, (model_id, revision, weight_path, weight_sha, required_names) in (
            model_contracts.items()
        ):
            section = model.get(label)
            if not isinstance(section, dict):
                raise ResourceInstallError(f"TMVec manifest model.{label} must be a mapping")
            if section.get("id") != model_id or section.get("revision") != revision:
                raise ResourceInstallError(
                    f"TMVec manifest model.{label} must pin {model_id}@{revision}"
                )
            files = section.get("files")
            if not isinstance(files, list) or not files:
                raise ResourceInstallError(
                    f"TMVec manifest model.{label}.files must be a non-empty list"
                )
            names: set[str] = set()
            weight_verified = False
            for index, item in enumerate(files):
                if not isinstance(item, dict):
                    raise ResourceInstallError(
                        f"TMVec manifest model.{label}.files[{index}] must be a mapping"
                    )
                relative = item.get("path")
                sha256 = item.get("sha256")
                member = cls._tmvec_member(
                    root,
                    relative,
                    f"model.{label}.files[{index}].path",
                )
                relative_text = str(relative)
                prefix = f"models/{'lobster_24M' if label == 'base' else 'tmvec-2'}/"
                if not relative_text.startswith(prefix):
                    raise ResourceInstallError(
                        f"model.{label} file must be under {prefix}: {relative_text}"
                    )
                if relative_text in manifest_members:
                    raise ResourceInstallError(
                        f"TMVec manifest lists a file more than once: {relative_text}"
                    )
                manifest_members.add(relative_text)
                names.add(Path(relative_text).name)
                if not isinstance(sha256, str) or _SHA256_RE.fullmatch(sha256) is None:
                    raise ResourceInstallError(
                        f"model.{label}.files[{index}].sha256 must be a lowercase SHA-256"
                    )
                if relative_text == weight_path:
                    if sha256 != weight_sha:
                        raise ResourceInstallError(
                            f"{weight_path} does not match the pinned upstream SHA-256"
                        )
                    weight_verified = True
                if verify_hashes and sha256_file(member) != sha256:
                    raise ResourceInstallError(f"TMVec file checksum mismatch: {member}")
            if names != required_names:
                raise ResourceInstallError(
                    f"model.{label}.files must contain exactly: "
                    + ", ".join(sorted(required_names))
                )
            if not weight_verified:
                raise ResourceInstallError(
                    f"model.{label}.files must include {weight_path}"
                )

        smoke = payload.get("smoke_query")
        if not isinstance(smoke, dict):
            raise ResourceInstallError("TMVec manifest smoke_query must be a mapping")
        smoke_id = smoke.get("id")
        sequence = smoke.get("sequence")
        expected_target_id = smoke.get("expected_target_id")
        if not isinstance(smoke_id, str) or not smoke_id:
            raise ResourceInstallError("smoke_query.id must be a non-empty string")
        if (
            not isinstance(sequence, str)
            or not sequence
            or len(sequence) > TMVEC2_ARCHITECTURE["max_sequence_length"]
            or re.fullmatch(r"[A-Z]+", sequence) is None
        ):
            raise ResourceInstallError(
                "smoke_query.sequence must contain 1-512 uppercase amino-acid letters"
            )
        if smoke.get("database") != "bfvd":
            raise ResourceInstallError("smoke_query.database must be 'bfvd'")
        if not isinstance(expected_target_id, str) or not expected_target_id:
            raise ResourceInstallError(
                "smoke_query.expected_target_id must be a non-empty string"
            )
        expected_score = smoke.get("expected_score")
        score_tolerance = smoke.get("score_tolerance")
        if (
            isinstance(expected_score, bool)
            or not isinstance(expected_score, (int, float))
            or not -1.0 <= float(expected_score) <= 1.0
        ):
            raise ResourceInstallError("smoke_query.expected_score must be in [-1, 1]")
        if (
            isinstance(score_tolerance, bool)
            or not isinstance(score_tolerance, (int, float))
            or not 0.0 <= float(score_tolerance) <= 2.0
        ):
            raise ResourceInstallError("smoke_query.score_tolerance must be in [0, 2]")

        references = smoke.get("reference_embeddings")
        if not isinstance(references, dict) or set(references) != {"cpu", "cuda"}:
            raise ResourceInstallError(
                "smoke_query.reference_embeddings must contain cpu and cuda mappings"
            )
        reference_paths: set[str] = set()
        for device_name in ("cpu", "cuda"):
            reference = references[device_name]
            label = f"smoke_query.reference_embeddings.{device_name}"
            if not isinstance(reference, dict):
                raise ResourceInstallError(f"{label} must be a mapping")
            if reference.get("dimensions") != TMVEC_EMBEDDING_WIDTH:
                raise ResourceInstallError(
                    f"{label}.dimensions must be {TMVEC_EMBEDDING_WIDTH}"
                )
            for tolerance_name in ("atol", "rtol"):
                tolerance = reference.get(tolerance_name)
                if (
                    isinstance(tolerance, bool)
                    or not isinstance(tolerance, (int, float))
                    or not 0.0 <= float(tolerance) <= 1.0
                ):
                    raise ResourceInstallError(
                        f"{label}.{tolerance_name} must be in [0, 1]"
                    )
            reference_value = reference.get("path")
            reference_path = cls._tmvec_member(
                root,
                reference_value,
                f"{label}.path",
            )
            if reference_value in reference_paths:
                raise ResourceInstallError(
                    "CPU and CUDA upstream references must use separate files"
                )
            reference_paths.add(reference_value)
            reference_sha = reference.get("sha256")
            if (
                not isinstance(reference_sha, str)
                or _SHA256_RE.fullmatch(reference_sha) is None
            ):
                raise ResourceInstallError(
                    f"{label}.sha256 must be a lowercase SHA-256"
                )
            reference_shape, reference_dtype = cls._tmvec_npy_header(reference_path)
            if reference_shape not in {
                (TMVEC_EMBEDDING_WIDTH,),
                (1, TMVEC_EMBEDDING_WIDTH),
            } or reference_dtype.kind not in "f":
                raise ResourceInstallError(
                    f"{label} must be a 512-value floating NPY array"
                )
            if verify_hashes and sha256_file(reference_path) != reference_sha:
                raise ResourceInstallError(
                    f"TMVec file checksum mismatch: {reference_path}"
                )

        database_payload = payload.get("databases")
        if not isinstance(database_payload, dict) or "bfvd" not in database_payload:
            raise ResourceInstallError("TMVec manifest databases must include bfvd")
        requested_databases = (
            list(database_payload) if databases is None else databases
        )
        unsupported = sorted(set(requested_databases) - set(cls.TMVEC_REQUIRED_FILES))
        if unsupported:
            raise ResourceInstallError(
                "unsupported TMVec database key(s): " + ", ".join(unsupported)
            )
        absent = sorted(set(requested_databases) - set(database_payload))
        if absent:
            raise ResourceInstallError(
                "TMVec manifest does not contain requested database(s): "
                + ", ".join(absent)
            )

        for name, database in database_payload.items():
            if name not in cls.TMVEC_REQUIRED_FILES or not isinstance(database, dict):
                raise ResourceInstallError(f"invalid TMVec database entry: {name}")
            if name == "bfvd":
                attribution = database.get("attribution")
                if attribution != {
                    "name": "BFVD",
                    "creator": "Kim, Rachel Seongeun",
                    "source_url": "https://bfvd.steineggerlab.workers.dev/",
                    "doi": "10.5281/zenodo.13993145",
                    "license": "CC BY 4.0",
                    "license_url": "https://creativecommons.org/licenses/by/4.0/",
                    "changes": (
                        "Converted BFVD protein sequences to "
                        "Lobster-24M/TMVec2 embeddings."
                    ),
                }:
                    raise ResourceInstallError(
                        "BFVD attribution must name the official source and CC BY 4.0 license"
                    )
            resolved_files: dict[str, tuple[Path, dict]] = {}
            for file_kind in ("embeddings", "metadata"):
                file_info = database.get(file_kind)
                if not isinstance(file_info, dict):
                    raise ResourceInstallError(
                        f"databases.{name}.{file_kind} must be a mapping"
                    )
                rows = file_info.get("rows")
                if isinstance(rows, bool) or not isinstance(rows, int) or rows < 1:
                    raise ResourceInstallError(
                        f"databases.{name}.{file_kind}.rows must be a positive integer"
                    )
                member = cls._tmvec_member(
                    root,
                    file_info.get("path"),
                    f"databases.{name}.{file_kind}.path",
                )
                sha256 = file_info.get("sha256")
                if not isinstance(sha256, str) or _SHA256_RE.fullmatch(sha256) is None:
                    raise ResourceInstallError(
                        f"databases.{name}.{file_kind}.sha256 must be a lowercase SHA-256"
                    )
                if verify_hashes and sha256_file(member) != sha256:
                    raise ResourceInstallError(f"TMVec file checksum mismatch: {member}")
                resolved_files[file_kind] = (member, file_info)

            embedding_path, embedding_info = resolved_files["embeddings"]
            embedding_shape, embedding_dtype = cls._tmvec_npy_header(embedding_path)
            expected_rows = embedding_info["rows"]
            if (
                embedding_shape != (expected_rows, TMVEC_EMBEDDING_WIDTH)
                or embedding_dtype.kind not in "iuf"
            ):
                raise ResourceInstallError(
                    f"databases.{name}.embeddings must have shape "
                    f"({expected_rows}, {TMVEC_EMBEDDING_WIDTH}) with a real numeric dtype"
                )

            metadata_path, metadata_info = resolved_files["metadata"]
            if metadata_path.suffix not in {".tsv", ".jsonl"}:
                raise ResourceInstallError(
                    f"databases.{name}.metadata must use TSV or JSONL, not pickle"
                )
            metadata_rows = 0
            found_expected_target = False
            seen_target_ids: set[str] = set()
            try:
                with metadata_path.open("r", encoding="utf-8", newline="") as handle:
                    if metadata_path.suffix == ".tsv":
                        records = csv.DictReader(handle, delimiter="\t")
                        if records.fieldnames is None or "id" not in records.fieldnames:
                            raise ResourceInstallError(
                                f"databases.{name}.metadata TSV must have an id column"
                            )
                        record_iter = records
                    else:
                        def _jsonl_records():
                            for line_number, line in enumerate(handle, start=1):
                                if not line.strip():
                                    continue
                                try:
                                    item = json.loads(line)
                                except json.JSONDecodeError as exc:
                                    raise ResourceInstallError(
                                        f"invalid JSONL in {metadata_path} line {line_number}: {exc}"
                                    ) from exc
                                if not isinstance(item, dict):
                                    raise ResourceInstallError(
                                        f"{metadata_path} line {line_number} must be a JSON object"
                                    )
                                yield item
                        record_iter = _jsonl_records()
                    for item in record_iter:
                        target_id = item.get("id")
                        if not isinstance(target_id, str) or not target_id:
                            raise ResourceInstallError(
                                f"databases.{name}.metadata contains an empty id"
                            )
                        if target_id in seen_target_ids:
                            raise ResourceInstallError(
                                f"databases.{name}.metadata contains duplicate id "
                                f"{target_id!r}"
                            )
                        seen_target_ids.add(target_id)
                        metadata_rows += 1
                        if name == "bfvd" and target_id == expected_target_id:
                            found_expected_target = True
            except (OSError, UnicodeError) as exc:
                raise ResourceInstallError(
                    f"cannot read TMVec metadata {metadata_path}: {exc}"
                ) from exc
            if metadata_rows != metadata_info["rows"] or metadata_rows != expected_rows:
                raise ResourceInstallError(
                    f"databases.{name} row count mismatch: embeddings={expected_rows}, "
                    f"metadata manifest={metadata_info['rows']}, metadata file={metadata_rows}"
                )
            if name == "bfvd" and not found_expected_target:
                raise ResourceInstallError(
                    "smoke_query.expected_target_id is absent from BFVD metadata"
                )

        return payload

    @classmethod
    def missing_tmvec_files(
        cls,
        tmvec_root: Optional[str | Path],
        databases: Optional[list[str]] = None,
    ) -> list[str]:
        root = (
            Path(tmvec_root).expanduser()
            if tmvec_root is not None
            else cls.default_tmvec_path()
        )
        try:
            cls.load_tmvec_manifest(
                root,
                verify_hashes=False,
                databases=["bfvd"] if databases is None else databases,
            )
        except (OSError, ResourceInstallError, ValueError) as exc:
            return [str(exc)]
        return []

    @classmethod
    def interproscan_available(cls, interproscan_dir: Optional[str | Path]) -> bool:
        if interproscan_dir is None:
            return False
        path = Path(interproscan_dir).expanduser()
        return all(
            (path / rel).is_file() and os.access(path / rel, os.X_OK)
            for rel in cls.INTERPROSCAN_REQUIRED_FILES
        )

    @classmethod
    def default_paths(cls, db_path: Path) -> dict[str, Optional[Path]]:
        marker_faa = db_path / "marker" / "marker.faa"
        return {
            "hmm_db": cls.hmm_db_path(db_path),
            "marker_faa_db": marker_faa if marker_faa.is_file() else None,
            "marker_db": db_path / "marker" / "marker.dmnd",
            "gene_taxonomy_faa_db": db_path / "genomes" / "combined_proteome.dmnd",
            "taxonomy_labels_file": db_path / "taxonomy" / "labels.tsv",
            "tmvec_database_dir": cls.default_tmvec_path(db_path),
        }

    @classmethod
    def _trusted_database_root(
        cls,
        database_root: str | Path | None,
        *,
        source: ResourceSource,
    ) -> Path | None:
        """Return a quiescent, receipt-authenticated core-resource root."""
        if database_root is None:
            return None
        root = cls.normalize_path(database_root)
        with sibling_install_lock(root):
            recover_pending_install(root)
            try:
                resource_root = active_installed_candidate(root, source)
                if resource_root is None:
                    return None
                manifest = load_resource_manifest(
                    resource_root,
                    expected_version=source.version,
                    expected_manifest_sha256=source.manifest_sha256,
                )
            except (OSError, ResourceManifestError):
                return None
            if not verified_install_receipt(
                resource_root,
                manifest,
                expected_archive_sha256=source.archive_sha256,
            ):
                return None
        return resource_root

    @classmethod
    def resolve_config_paths(
        cls,
        config: dict,
        config_path: Optional[Path] = None,
    ) -> dict:
        """Resolve relative paths and auto-fill missing database entries."""
        resolved = dict(config)
        base_dir = config_path.parent if config_path else Path.cwd()
        env_database_root = os.environ.get("VIROSYNC_DB_ROOT")
        if env_database_root and not resolved.get("database_root"):
            resolved["database_root"] = env_database_root

        def _abs_path(value: Optional[str | Path]) -> Optional[str]:
            if value is None:
                return None
            value_path = Path(value)
            if not value_path.is_absolute():
                value_path = base_dir / value_path
            return str(value_path)

        for key in (
            "database_root",
            "hmm_db",
            "hmm_database",
            "marker_faa_db",
            "marker_db",
            "gene_taxonomy_faa_db",
            "taxonomy_labels_file",
            "tmvec_database_dir",
            "interproscan_dir",
            "gvclass_path",
            "viral_structure_db",
        ):
            if key in resolved and resolved[key] is not None:
                resolved[key] = _abs_path(resolved[key])

        if isinstance(resolved.get("phase3"), dict):
            phase3_cfg = dict(resolved["phase3"])
            for key in (
                "interproscan_dir",
                "tmvec_database_dir",
                "gvclass_path",
                "viral_structure_db",
            ):
                if phase3_cfg.get(key) is not None:
                    phase3_cfg[key] = _abs_path(phase3_cfg[key])
            resolved["phase3"] = phase3_cfg

        # Check for required databases (support both hmm_db and hmm_database keys)
        hmm_key = "hmm_database" if "hmm_database" in resolved else "hmm_db"
        required_checks = [
            (hmm_key, resolved.get(hmm_key)),
            ("marker_db", resolved.get("marker_db")),
            ("gene_taxonomy_faa_db", resolved.get("gene_taxonomy_faa_db")),
        ]
        missing = [
            key
            for key, value in required_checks
            if not value or not Path(str(value)).exists()
        ]

        if missing:
            db_root = resolved.get("database_root")
            core_source = resolved.get("core_resources_url")
            selected_source = cls._resolve_core_source(
                core_source,
                resolved.get("core_resources_version"),
                resolved.get("core_resources_sha256"),
                resolved.get("core_resources_manifest_sha256"),
            )
            db_path = cls._trusted_database_root(
                db_root,
                source=selected_source,
            )
            if db_path is None:
                db_path = cls.setup_database(
                    database_path=db_root,
                    source=core_source,
                    version=resolved.get("core_resources_version"),
                    archive_sha256=resolved.get("core_resources_sha256"),
                    manifest_sha256=resolved.get("core_resources_manifest_sha256"),
                    full=True,
                )
            defaults = cls.default_paths(db_path)
            for key, value in defaults.items():
                if not resolved.get(key) or not Path(str(resolved[key])).exists():
                    resolved[key] = str(value) if value is not None else None

            # Propagate defaults into the nested phase3 dict so config loading
            # picks up resolved paths instead of None.
            if isinstance(resolved.get("phase3"), dict):
                phase3_cfg = resolved["phase3"]
                for key in ("tmvec_database_dir", "interproscan_dir"):
                    if phase3_cfg.get(key) is None and resolved.get(key):
                        phase3_cfg[key] = resolved[key]

        # Taxonomy labels are optional for DB setup but required for host-signature
        # lineage scoring. Backfill them when the file exists in the resolved DB root.
        if not resolved.get("taxonomy_labels_file"):
            candidate_roots: list[Path] = []
            db_root = resolved.get("database_root")
            if db_root:
                candidate_roots.append(Path(str(db_root)).expanduser())
            for key in (hmm_key, "marker_db", "marker_faa_db", "gene_taxonomy_faa_db"):
                value = resolved.get(key)
                if not value:
                    continue
                path = Path(str(value)).expanduser()
                candidate_roots.append(path.parent.parent)

            seen_roots = set()
            for root in candidate_roots:
                root_resolved = root.resolve()
                if root_resolved in seen_roots:
                    continue
                seen_roots.add(root_resolved)
                labels_path = root_resolved / "taxonomy" / "labels.tsv"
                if labels_path.exists():
                    resolved["taxonomy_labels_file"] = str(labels_path)
                    break

        return resolved
