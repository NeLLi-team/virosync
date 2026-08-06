"""
Database management utilities for ViroSync.
Handles downloading and validating the reference database bundle.
"""

from __future__ import annotations

import json
import logging
import math
import os
import pickletools
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


class ViroSyncDatabaseManager:
    """Manages ViroSync reference databases."""

    HMM_REQUIRED_FILES = ["models/combined.hmm"]
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
            "bfvd/bfvd_annotations.npy",
        ],
        "cath": [
            "cath/cath_large.npy",
            "cath/cath_large_metadata.npy",
        ],
        "swissprot": [
            "swissprot/swiss_large.npy",
            "swissprot/swiss_large_metadata.npy",
        ],
        "pdb": [
            "pdb/embeddings.npy",
            "pdb/metadata.npy",
        ],
    }
    INTERPROSCAN_REQUIRED_FILES = ["interproscan.sh"]

    DATABASE_SOURCES = [
        {
            "version": "v1.0.7",
            "source": "https://dl.newlineages.com/virosync/resources_v1_0_7_runtime.tar.gz",
            "filename": "resources_v1_0_7_runtime.tar.gz",
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
    def _extract_archive(archive_path: Path, target_dir: Path) -> None:
        """Compatibility alias for the strict core-resource extractor."""
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
        force: bool = False,
        progress_callback=None,
    ) -> None:
        target_dir.parent.mkdir(parents=True, exist_ok=True)
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
            if progress_callback is not None:
                progress_callback(70, "extracting")
            logger.info("Extracting archive: %s", archive_path)
            cls._extract_optional_archive(archive_path, payload_path)
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
        force: bool = False,
        progress_callback=None,
    ) -> bool:
        if source is None:
            missing = [rel for rel in required_files if not (target_path / rel).exists()]
            if missing:
                logger.warning(
                    "%s not installed and no source URL/path provided. Missing: %s",
                    name,
                    missing,
                )
                return False
            metadata_path = target_path / "DB_METADATA.json"
            if not metadata_path.exists():
                cls._write_database_metadata(
                    target_path=target_path,
                    component=name,
                    version=version,
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
                force=force,
                progress_callback=progress_callback,
            )
        except Exception as exc:
            logger.error("%s installation failed: %s", name, exc)
            return False

        missing = [rel for rel in required_files if not (target_path / rel).exists()]
        if missing:
            logger.warning("%s installation incomplete; missing files: %s", name, missing)
            return False

        cls._write_database_metadata(
            target_path=target_path,
            component=name,
            version=version,
            source=source,
            extra={"required_files": required_files},
        )
        return True

    @classmethod
    def missing_tmvec_files(
        cls,
        tmvec_root: Optional[str | Path],
        databases: Optional[list[str]] = None,
    ) -> list[str]:
        root = Path(tmvec_root).expanduser() if tmvec_root is not None else cls.default_tmvec_path()
        dbs = databases or ["bfvd"]
        missing: list[str] = []

        def _ready_file(path: Path) -> bool:
            try:
                return path.is_file() and path.stat().st_size > 0
            except OSError:
                return False

        def _npy_header(path: Path) -> tuple[tuple[int, ...], np.dtype, int, int]:
            with path.open("rb") as handle:
                version = np.lib.format.read_magic(handle)
                shape, _fortran_order, dtype = np.lib.format._read_array_header(
                    handle,
                    version,
                )
                payload_offset = handle.tell()
                file_size = path.stat().st_size
                if not dtype.hasobject:
                    expected_bytes = math.prod(shape) * dtype.itemsize
                    if file_size - payload_offset < expected_bytes:
                        raise ValueError(
                            f"numeric payload is truncated: expected {expected_bytes} "
                            f"bytes, found {max(0, file_size - payload_offset)}"
                        )
                else:
                    try:
                        for _opcode, _argument, _position in pickletools.genops(
                            handle
                        ):
                            pass
                    except Exception as exc:
                        raise ValueError(
                            f"object pickle payload is invalid: {exc}"
                        ) from exc
                    if handle.tell() != file_size:
                        raise ValueError(
                            "object pickle payload has unexpected trailing data"
                        )
            return tuple(shape), dtype, payload_offset, file_size

        def _group_issues(group: list[Path]) -> list[str]:
            issues: list[str] = []
            if len(group) != 2:
                return ["expected an embeddings/metadata file pair"]
            headers = []
            for path in group:
                if not _ready_file(path):
                    issues.append(f"{path} must be a non-empty regular file")
                    headers.append(None)
                    continue
                try:
                    headers.append(_npy_header(path))
                except Exception as exc:
                    issues.append(f"{path} is not a valid NPY asset: {exc}")
                    headers.append(None)

            embedding_header, metadata_header = headers
            if embedding_header is not None:
                embedding_shape, embedding_dtype, _, _ = embedding_header
                if len(embedding_shape) != 2:
                    issues.append(
                        f"{group[0]} embeddings must be 2-D, got {embedding_shape}"
                    )
                elif any(size < 1 for size in embedding_shape):
                    issues.append(
                        f"{group[0]} embeddings must have a non-empty shape, "
                        f"got {embedding_shape}"
                    )
                elif embedding_shape[1] != TMVEC_EMBEDDING_WIDTH:
                    issues.append(
                        f"{group[0]} embeddings must have width "
                        f"{TMVEC_EMBEDDING_WIDTH}, got {embedding_shape[1]}"
                    )
                if embedding_dtype.kind not in "iuf":
                    issues.append(
                        f"{group[0]} embeddings must use a real numeric dtype, "
                        f"got {embedding_dtype}"
                    )

            if metadata_header is not None:
                metadata_shape, _, _, _ = metadata_header
                if not metadata_shape:
                    issues.append(f"{group[1]} metadata must have a row dimension")
                elif metadata_shape[0] < 1:
                    issues.append(f"{group[1]} metadata must have at least one row")
                elif (
                    embedding_header is not None
                    and len(embedding_header[0]) == 2
                    and metadata_shape[0] != embedding_header[0][0]
                ):
                    issues.append(
                        "TMVec embedding/metadata row count mismatch: "
                        f"{embedding_header[0][0]} != {metadata_shape[0]}"
                    )
            return issues

        def _candidate_groups(name: str) -> Optional[list[list[Path]]]:
            if name == "bfvd":
                return [
                    [root / "bfvd" / "bfvd_embeddings.npy", root / "bfvd" / "bfvd_annotations.npy"],
                    [root / "bfvd_embeddings.npy", root / "bfvd_annotations.npy"],
                    [
                        root / "tmvec_embeddings" / "bfvd_embeddings.npy",
                        root / "tmvec_embeddings" / "bfvd_annotations.npy",
                    ],
                ]
            if name == "cath":
                return [
                    [root / "cath" / "cath_large.npy", root / "cath" / "cath_large_metadata.npy"],
                    [root / "cath_large.npy", root / "cath_large_metadata.npy"],
                ]
            if name == "swissprot":
                return [
                    [root / "swissprot" / "swiss_large.npy", root / "swissprot" / "swiss_large_metadata.npy"],
                    [root / "swiss_large.npy", root / "swiss_large_metadata.npy"],
                ]
            if name == "pdb":
                return [
                    [root / "pdb" / "embeddings.npy", root / "pdb" / "metadata.npy"],
                    [root / "embeddings.npy", root / "metadata.npy"],
                ]
            rel_paths = cls.TMVEC_REQUIRED_FILES.get(name)
            if rel_paths is None:
                return None
            return [[root / rel for rel in rel_paths]]

        for name in dbs:
            path_groups = _candidate_groups(name)
            if path_groups is None:
                missing.append(f"{name}: unsupported database key")
                continue

            candidate_issues = [_group_issues(group) for group in path_groups]
            if any(not issues for issues in candidate_issues):
                continue

            expected = " OR ".join(
                "[" + ", ".join(str(path) for path in group) + "]"
                for group in path_groups
            )
            details = " OR ".join(
                "; ".join(issues) for issues in candidate_issues if issues
            )
            missing.append(f"{name}: no valid TMVec NPY pair {expected}: {details}")
        return missing

    @classmethod
    def interproscan_available(cls, interproscan_dir: Optional[str | Path]) -> bool:
        if interproscan_dir is None:
            return False
        path = Path(interproscan_dir).expanduser()
        return all((path / rel).exists() for rel in cls.INTERPROSCAN_REQUIRED_FILES)

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

            # Propagate defaults into nested phase3 dict so that
            # _from_flat_dict picks up resolved paths instead of None.
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
