"""Hermetic contracts for authenticated, transactional resource installs.

Every archive and resource tree in this module is a tiny synthetic fixture under
``tmp_path``.  Network access and child processes are rejected so the tests can
never reach a published bundle, a project resource tree, or an external tool.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
import hashlib
import inspect
import io
import json
import os
from pathlib import Path
import re
import socket
import stat
import subprocess
import tarfile
import threading
import time
from typing import Any

import pytest
import yaml

import virosync.utils.database_manager as database_manager
from virosync.utils.database_manager import ViroSyncDatabaseManager
import virosync.utils.resource_installer as resource_installer
import virosync.utils.resource_manifest as resource_manifest

_REAL_SETUP_DATABASE = ViroSyncDatabaseManager.setup_database.__func__
_REAL_RESOLVE_CONFIG_PATHS = ViroSyncDatabaseManager.resolve_config_paths.__func__
_REAL_DIAMOND_SEQUENCE_COUNT = resource_manifest.diamond_sequence_count
_REAL_SUBPROCESS_RUN = subprocess.run

VERSION = "v1.0.6"
PRIOR_VERSION = "v1.0.5"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")

CORE_PAYLOADS = (
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

SEMANTIC_COUNTS = {
    "hmm_models": 1,
    "model_annotations": 1,
    "marker_proteins": 1,
    "marker_diamond_sequences": 1,
    "proteome_diamond_sequences": 1,
    "taxonomy_labels": 1,
}

INSTALL_FAULT_PHASES = (
    "after_download",
    "after_archive_verify",
    "after_extract",
    "after_stage_validate",
    "after_candidate_promote",
    "after_journal_write",
    "after_prior_move",
    "after_pointer_prepare",
    "after_pointer_activate",
    "after_activation_validate",
    "after_journal_clear",
)

PRE_ACTIVATION_PHASES = frozenset(
    INSTALL_FAULT_PHASES[: INSTALL_FAULT_PHASES.index("after_pointer_activate")]
)


class UnexpectedExternalActivity(AssertionError):
    """A resource test attempted network or child-process activity."""


class InjectedCrash(BaseException):
    """Simulate abrupt termination without exercising process management."""


@pytest.fixture(autouse=True)
def restore_real_installer(
    monkeypatch: pytest.MonkeyPatch,
    _block_database_auto_download: None,
) -> None:
    """Override the suite-wide download guard for this installer-only module."""

    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "setup_database",
        classmethod(_REAL_SETUP_DATABASE),
    )


@pytest.fixture(autouse=True)
def forbid_network_and_processes(monkeypatch: pytest.MonkeyPatch) -> None:
    """Make accidental real downloads and semantic tools fail immediately."""

    def _reject_network(*_args: object, **_kwargs: object) -> None:
        raise UnexpectedExternalActivity("resource test attempted network access")

    def _reject_process(*_args: object, **_kwargs: object) -> None:
        raise UnexpectedExternalActivity("resource test attempted a child process")

    def _guard_diamond(
        database: Path,
        *,
        command_runner: Callable[..., object] = _REAL_SUBPROCESS_RUN,
    ) -> int:
        if command_runner is _REAL_SUBPROCESS_RUN:
            raise UnexpectedExternalActivity(
                "resource test attempted DIAMOND without its injected runner"
            )
        return _REAL_DIAMOND_SEQUENCE_COUNT(
            database,
            command_runner=command_runner,
        )

    monkeypatch.setattr(socket, "create_connection", _reject_network)
    monkeypatch.setattr(socket.socket, "connect", _reject_network)
    monkeypatch.setattr(subprocess, "run", _reject_process)
    monkeypatch.setattr(subprocess, "Popen", _reject_process)
    monkeypatch.setattr(subprocess, "check_call", _reject_process)
    monkeypatch.setattr(subprocess, "check_output", _reject_process)
    monkeypatch.setattr(resource_manifest, "diamond_sequence_count", _guard_diamond)


def _sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _payloads(version: str = VERSION) -> dict[str, bytes]:
    return {
        "DB_VERSION": f"{version}\n".encode(),
        "DATABASE_README.txt": f"Synthetic ViroSync resources {version}\n".encode(),
        "models/combined.hmm": (b"HMMER3/f [synthetic]\nNAME  VS000001\nLENG  4\n//\n"),
        "models/combined.hmm.h3f": b"synthetic-h3f\n",
        "models/combined.hmm.h3i": b"synthetic-h3i\n",
        "models/combined.hmm.h3m": b"synthetic-h3m\n",
        "models/combined.hmm.h3p": b"synthetic-h3p\n",
        "models/model_annotations_with_interpro.tsv": (
            b"model\tannotation\nVS000001\tsynthetic\n"
        ),
        "models/og_marker_name_map.tsv": b"model\tmarker\nVS000001\tMCP\n",
        "marker/marker.faa": b">synthetic_marker\nMPEP\n",
        "marker/marker.dmnd": b"synthetic-marker-diamond\n",
        "genomes/combined_proteome.dmnd": b"synthetic-proteome-diamond\n",
        "taxonomy/labels.tsv": b"genome\tlineage\nsynthetic\tNCLDV\n",
    }


def _role(path: str) -> str:
    if path == "DB_VERSION":
        return "bundle_version"
    if path == "DATABASE_README.txt":
        return "bundle_documentation"
    if path == "models/combined.hmm":
        return "hmm_models"
    if ".h3" in path:
        return "hmm_index"
    if path.endswith(".dmnd"):
        return "diamond_database"
    if path.endswith(".faa"):
        return "marker_proteins"
    if path.endswith("labels.tsv"):
        return "taxonomy_labels"
    return "annotation_table"


def _manifest_document(
    payloads: Mapping[str, bytes],
    *,
    version: str = VERSION,
    semantic_counts: Mapping[str, int] | None = None,
) -> dict[str, Any]:
    counts = dict(SEMANTIC_COUNTS)
    if semantic_counts is not None:
        counts.update(semantic_counts)
    return {
        "schema_version": 1,
        "resource_version": version,
        "version": version,
        "files": [
            {
                "path": path,
                "size": len(payloads[path]),
                "sha256": _sha256_bytes(payloads[path]),
                "role": _role(path),
            }
            for path in CORE_PAYLOADS
            if path in payloads
        ],
        "semantic_counts": counts,
    }


def _manifest_bytes(document: Mapping[str, Any]) -> bytes:
    return (json.dumps(document, indent=2, sort_keys=True) + "\n").encode()


def _write_tree(
    root: Path,
    payloads: Mapping[str, bytes],
    *,
    manifest: bytes | None = None,
    include_manifest: bool = True,
) -> Path:
    for relative, content in payloads.items():
        destination = root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    if include_manifest:
        document = _manifest_document(payloads)
        (root / "RESOURCE_MANIFEST.json").write_bytes(
            manifest if manifest is not None else _manifest_bytes(document)
        )
    return root


def _add_bytes(handle: tarfile.TarFile, name: str, content: bytes) -> None:
    member = tarfile.TarInfo(name)
    member.size = len(content)
    member.mode = 0o644
    member.mtime = 0
    member.uid = 0
    member.gid = 0
    member.uname = ""
    member.gname = ""
    handle.addfile(member, io.BytesIO(content))


def _add_member(handle: tarfile.TarFile, member: tarfile.TarInfo) -> None:
    member.mtime = 0
    member.uid = 0
    member.gid = 0
    member.uname = ""
    member.gname = ""
    handle.addfile(member)


def _build_archive(
    path: Path,
    payloads: Mapping[str, bytes] | None = None,
    *,
    manifest: bytes | None = None,
    include_manifest: bool = True,
    root: str = "virosync",
) -> tuple[Path, bytes | None]:
    actual_payloads = dict(payloads if payloads is not None else _payloads())
    actual_manifest = manifest
    if include_manifest and actual_manifest is None:
        actual_manifest = _manifest_bytes(_manifest_document(actual_payloads))
    with tarfile.open(path, "w:gz") as handle:
        for relative, content in actual_payloads.items():
            _add_bytes(handle, f"{root}/{relative}", content)
        if include_manifest and actual_manifest is not None:
            _add_bytes(
                handle,
                f"{root}/RESOURCE_MANIFEST.json",
                actual_manifest,
            )
    return path, actual_manifest


def _build_schema2_archive(
    tmp_path: Path,
    bundle_kind: str,
) -> tuple[Path, bytes, resource_manifest.ResourceManifest]:
    payloads = _payloads()
    payloads["models/pfam_virosync_screening.hmm"] = (
        b"HMMER3/f\nNAME  PfamOne\nGA    10.0 10.0;\n//\n"
    )
    resources = _write_tree(
        tmp_path / "schema2-resources",
        payloads,
        include_manifest=False,
    )
    runtime, source = resource_manifest.build_split_resource_manifests(
        resources,
        VERSION,
        diamond_sequence_counter=lambda _payload: 1,
    )
    if bundle_kind == "runtime":
        manifest, manifest_bytes = runtime
        payload_files = resource_manifest.RUNTIME_RESOURCE_FILES
    elif bundle_kind == "source":
        manifest, manifest_bytes = source
        payload_files = resource_manifest.SOURCE_RESOURCE_FILES
    else:  # pragma: no cover - protects the test helper itself
        raise AssertionError(bundle_kind)
    archive, embedded_manifest = _build_archive(
        tmp_path / f"{bundle_kind}.tar.gz",
        {relative: payloads[relative] for relative in payload_files},
        manifest=manifest_bytes,
    )
    assert embedded_manifest == manifest_bytes
    return archive, manifest_bytes, manifest


class SemanticProbe:
    """Return deterministic semantic counts while recording invocations."""

    def __init__(self, counts: Mapping[str, int] | None = None) -> None:
        self.counts = dict(counts if counts is not None else SEMANTIC_COUNTS)
        self.calls: list[tuple[tuple[object, ...], dict[str, object]]] = []
        self._guard = threading.Lock()
        self.active = 0
        self.max_active = 0
        self.delay_seconds = 0.0

    def __call__(self, *args: object, **kwargs: object) -> dict[str, int]:
        with self._guard:
            self.calls.append((args, kwargs))
            self.active += 1
            self.max_active = max(self.max_active, self.active)
        try:
            if self.delay_seconds:
                time.sleep(self.delay_seconds)
            return dict(self.counts)
        finally:
            with self._guard:
                self.active -= 1


def _setup(
    target: Path,
    archive: Path,
    manifest: bytes,
    *,
    version: str = VERSION,
    archive_sha256: str | None = None,
    manifest_sha256: str | None = None,
    force: bool = True,
    full: bool = True,
    semantic_runner: Callable[..., Mapping[str, int]] | None = None,
    fault_injector: Callable[[str], None] | None = None,
) -> Path:
    result = ViroSyncDatabaseManager.setup_database(
        database_path=str(target),
        source=str(archive),
        version=version,
        force=force,
        archive_sha256=archive_sha256 or _sha256_file(archive),
        manifest_sha256=manifest_sha256 or _sha256_bytes(manifest),
        full=full,
        semantic_runner=semantic_runner or SemanticProbe(),
        fault_injector=fault_injector,
    )
    return Path(result)


def _verify(
    target: Path,
    manifest: bytes,
    *,
    version: str = VERSION,
    full: bool = False,
    semantic_runner: Callable[..., Mapping[str, int]] | None = None,
) -> object:
    return ViroSyncDatabaseManager.verify_database(
        target,
        expected_version=version,
        manifest_sha256=_sha256_bytes(manifest),
        full=full,
        semantic_runner=semantic_runner,
    )


def _safe_extract(archive: Path, target: Path) -> None:
    ViroSyncDatabaseManager._safe_extract_archive(archive, target)


def _expect_rejected(
    operation: Callable[[], object],
    *,
    keywords: tuple[str, ...],
) -> BaseException:
    try:
        operation()
    except UnexpectedExternalActivity as exc:
        pytest.fail(f"installer attempted external activity: {exc}")
    except BaseException as exc:
        message = str(exc).lower()
        assert any(
            keyword in message for keyword in keywords
        ), f"wrong rejection {type(exc).__name__}: {exc}"
        return exc
    pytest.fail("invalid resource input was accepted")


def _snapshot(root: Path) -> dict[str, tuple[str, str]]:
    if not root.exists() and not root.is_symlink():
        return {}
    base = root.resolve(strict=True) if root.is_symlink() else root
    snapshot: dict[str, tuple[str, str]] = {}
    for path in sorted(base.rglob("*")):
        relative = path.relative_to(base).as_posix()
        if path.is_symlink():
            snapshot[relative] = ("symlink", os.readlink(path))
        elif path.is_file():
            snapshot[relative] = ("file", _sha256_file(path))
        elif path.is_dir():
            snapshot[relative] = ("directory", "")
        else:
            snapshot[relative] = ("special", oct(path.lstat().st_mode))
    return snapshot


def _manifest_digest(result: object) -> str | None:
    if isinstance(result, str) and SHA256_RE.fullmatch(result.lower()):
        return result.lower()
    if isinstance(result, Mapping):
        for name in ("manifest_sha256", "resource_manifest_sha256"):
            value = result.get(name)
            if isinstance(value, str) and SHA256_RE.fullmatch(value.lower()):
                return value.lower()
    for name in ("manifest_sha256", "resource_manifest_sha256"):
        value = getattr(result, name, None)
        if isinstance(value, str) and SHA256_RE.fullmatch(value.lower()):
            return value.lower()
    return None


def _source_digest(record: Mapping[str, object], *names: str) -> str | None:
    for name in names:
        value = record.get(name)
        if isinstance(value, str) and SHA256_RE.fullmatch(value.lower()):
            return value.lower()
    return None


def _record_reuse_checks(
    monkeypatch: pytest.MonkeyPatch,
    events: list[str],
) -> None:
    real_receipt_check = resource_installer._install_receipt_is_current
    real_validate = ViroSyncDatabaseManager._validate_core_tree.__func__

    def _record_receipt_check(
        root: Path,
        source: resource_installer.ResourceSource,
        *,
        manifest: object | None = None,
    ) -> bool:
        events.append("receipt")
        return real_receipt_check(root, source, manifest=manifest)

    def _record_verification(cls, *args: object, **kwargs: object) -> object:
        events.append("verify")
        return real_validate(cls, *args, **kwargs)

    monkeypatch.setattr(
        resource_installer,
        "_install_receipt_is_current",
        _record_receipt_check,
    )
    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "_validate_core_tree",
        classmethod(_record_verification),
    )


def _v106_source() -> Mapping[str, object] | None:
    return next(
        (
            candidate
            for candidate in ViroSyncDatabaseManager.DATABASE_SOURCES
            if candidate.get("version") == VERSION
        ),
        None,
    )


def _prior_tree(target: Path) -> dict[str, tuple[str, str]]:
    payloads = _payloads(PRIOR_VERSION)
    manifest = _manifest_bytes(_manifest_document(payloads, version=PRIOR_VERSION))
    _write_tree(target, payloads, manifest=manifest)
    (target / "prior-sentinel.txt").write_text("prior\n", encoding="utf-8")
    return _snapshot(target)


def _assert_no_visible_members(target: Path) -> None:
    if target.exists():
        assert not any(target.rglob("*")), _snapshot(target)


def test_installer_api_exposes_authentication_and_verification_controls() -> None:
    setup_parameters = inspect.signature(
        ViroSyncDatabaseManager.setup_database
    ).parameters
    verify = getattr(ViroSyncDatabaseManager, "verify_database", None)
    extractor = getattr(ViroSyncDatabaseManager, "_safe_extract_archive", None)

    assert {
        "archive_sha256",
        "manifest_sha256",
        "full",
        "semantic_runner",
        "fault_injector",
    } <= set(setup_parameters)
    assert callable(verify)
    assert callable(extractor)
    assert {"expected_version", "manifest_sha256", "full"} <= set(
        inspect.signature(verify).parameters
    )


def test_v106_source_and_shipped_configs_pin_the_same_two_digests() -> None:
    record = _v106_source()
    assert record is not None, "no v1.0.6 source metadata"
    archive_sha = _source_digest(record, "archive_sha256", "sha256")
    manifest_sha = _source_digest(
        record,
        "manifest_sha256",
        "resource_manifest_sha256",
    )
    assert archive_sha is not None
    assert manifest_sha is not None

    repo_root = Path(database_manager.__file__).resolve().parents[3]
    for relative in (
        "config/orchestration.yaml",
        "config/orchestration_archaeal.yaml",
    ):
        document = yaml.safe_load((repo_root / relative).read_text(encoding="utf-8"))
        orchestration = document["orchestration"]
        assert "resources_v1_0_6.tar.gz" in orchestration["core_resources_url"]
        assert orchestration["core_resources_sha256"] == archive_sha
        assert orchestration["core_resources_manifest_sha256"] == manifest_sha


def test_archive_digest_mismatch_from_flipped_byte_fails_before_install(
    tmp_path: Path,
) -> None:
    clean, manifest = _build_archive(tmp_path / "clean.tar.gz")
    assert manifest is not None
    expected_archive_sha = _sha256_file(clean)
    content = bytearray(clean.read_bytes())
    content[len(content) // 2] ^= 1
    flipped = tmp_path / "flipped.tar.gz"
    flipped.write_bytes(content)
    target = tmp_path / "resources" / "virosync"

    _expect_rejected(
        lambda: _setup(
            target,
            flipped,
            manifest,
            archive_sha256=expected_archive_sha,
        ),
        keywords=("archive", "sha", "digest", "checksum"),
    )
    assert not target.exists() and not target.is_symlink()


@pytest.mark.parametrize(
    ("case", "keywords"),
    [
        ("schema", ("schema",)),
        ("path", ("path", "relative", "..")),
        ("duplicate_path", ("duplicate", "path")),
        ("missing_entry", ("missing", "incomplete", "payload")),
        ("size", ("size", "nonzero", "positive")),
        ("role", ("role",)),
        ("version", ("version", VERSION, "v9.9.9")),
    ],
)
def test_manifest_schema_path_size_role_and_version_are_validated(
    tmp_path: Path,
    case: str,
    keywords: tuple[str, ...],
) -> None:
    payloads = _payloads()
    document = _manifest_document(payloads)
    if case == "schema":
        document["schema_version"] = 999
    elif case == "path":
        document["files"][0]["path"] = "../DB_VERSION"
    elif case == "duplicate_path":
        document["files"].append(dict(document["files"][0]))
    elif case == "missing_entry":
        document["files"].pop()
    elif case == "size":
        document["files"][0]["size"] = 0
    elif case == "role":
        document["files"][0]["role"] = "wrong_role"
    elif case == "version":
        document["resource_version"] = "v9.9.9"
        document["version"] = "v9.9.9"
    else:  # pragma: no cover - protects the test table itself
        raise AssertionError(case)

    manifest = _manifest_bytes(document)
    target = _write_tree(tmp_path / case, payloads, manifest=manifest)
    _expect_rejected(
        lambda: _verify(target, manifest),
        keywords=tuple(token.lower() for token in keywords),
    )


@pytest.mark.parametrize(
    "case",
    [
        "traversal",
        "absolute",
        "unexpected_root",
        "symlink",
        "hardlink",
        "duplicate",
        "fifo",
        "device",
    ],
)
def test_safe_extractor_rejects_unsafe_archive_shapes_before_visibility(
    tmp_path: Path,
    case: str,
) -> None:
    archive = tmp_path / f"{case}.tar.gz"
    with tarfile.open(archive, "w:gz") as handle:
        if case == "traversal":
            _add_bytes(handle, "virosync/safe-first.txt", b"partial\n")
            _add_bytes(handle, "virosync/../../escape.txt", b"escape\n")
        elif case == "absolute":
            _add_bytes(handle, "/virosync/absolute.txt", b"escape\n")
        elif case == "unexpected_root":
            _add_bytes(handle, "not-virosync/file.txt", b"wrong root\n")
        elif case in {"symlink", "hardlink"}:
            _add_bytes(handle, "virosync/source.txt", b"source\n")
            member = tarfile.TarInfo(f"virosync/{case}.txt")
            member.type = tarfile.SYMTYPE if case == "symlink" else tarfile.LNKTYPE
            member.linkname = (
                "../../outside" if case == "symlink" else "virosync/source.txt"
            )
            member.mode = 0o777
            _add_member(handle, member)
        elif case == "duplicate":
            _add_bytes(handle, "virosync/duplicate.txt", b"first\n")
            _add_bytes(handle, "virosync/duplicate.txt", b"second\n")
        elif case == "fifo":
            member = tarfile.TarInfo("virosync/fifo")
            member.type = tarfile.FIFOTYPE
            member.mode = 0o600
            _add_member(handle, member)
        elif case == "device":
            member = tarfile.TarInfo("virosync/device")
            member.type = tarfile.CHRTYPE
            member.devmajor = 1
            member.devminor = 3
            member.mode = 0o600
            _add_member(handle, member)

    target = tmp_path / "extract" / "virosync"
    with pytest.raises(Exception):
        _safe_extract(archive, target)
    _assert_no_visible_members(target)
    assert not (tmp_path / "escape.txt").exists()
    assert not (tmp_path / "outside").exists()


@pytest.mark.parametrize(
    ("case", "keywords"),
    [
        ("missing_manifest", ("manifest",)),
        ("missing_file", ("missing", "marker.dmnd", "required")),
        ("zero_byte", ("empty", "zero", "size", "nonzero")),
        ("wrong_version", ("version", PRIOR_VERSION, VERSION)),
        ("same_size_corruption", ("sha", "digest", "checksum", "corrupt")),
    ],
)
def test_invalid_bundle_never_becomes_active(
    tmp_path: Path,
    case: str,
    keywords: tuple[str, ...],
) -> None:
    payloads = _payloads()
    manifest_document = _manifest_document(payloads)
    include_manifest = True
    if case == "missing_manifest":
        include_manifest = False
    elif case == "missing_file":
        payloads.pop("marker/marker.dmnd")
    elif case == "zero_byte":
        payloads["marker/marker.dmnd"] = b""
    elif case == "wrong_version":
        payloads = _payloads(PRIOR_VERSION)
        manifest_document = _manifest_document(payloads, version=PRIOR_VERSION)
    elif case == "same_size_corruption":
        original = payloads["models/og_marker_name_map.tsv"]
        payloads["models/og_marker_name_map.tsv"] = (
            bytes([original[0] ^ 1]) + original[1:]
        )

    manifest = _manifest_bytes(manifest_document)
    archive, embedded_manifest = _build_archive(
        tmp_path / f"{case}.tar.gz",
        payloads,
        manifest=manifest,
        include_manifest=include_manifest,
    )
    target = tmp_path / "resources" / "virosync"
    _expect_rejected(
        lambda: _setup(
            target,
            archive,
            embedded_manifest or manifest,
            version=VERSION,
            manifest_sha256=_sha256_bytes(manifest),
        ),
        keywords=tuple(token.lower() for token in keywords),
    )
    assert not target.exists() and not target.is_symlink()


@pytest.mark.parametrize("case", ["zero_byte", "wrong_version"])
def test_fast_verify_rejects_invalid_preexisting_resources(
    tmp_path: Path,
    case: str,
) -> None:
    if case == "wrong_version":
        payloads = _payloads(PRIOR_VERSION)
        manifest = _manifest_bytes(_manifest_document(payloads, version=PRIOR_VERSION))
    else:
        payloads = _payloads()
        manifest = _manifest_bytes(_manifest_document(payloads))
        payloads["marker/marker.dmnd"] = b""
    target = _write_tree(tmp_path / case, payloads, manifest=manifest)

    _expect_rejected(
        lambda: _verify(target, manifest, version=VERSION, full=False),
        keywords=("empty", "zero", "size", "version", PRIOR_VERSION, VERSION),
    )


def test_unpinned_custom_archive_is_not_stamped_with_default_version(
    tmp_path: Path,
) -> None:
    payloads = _payloads("v9.9.9")
    manifest = _manifest_bytes(_manifest_document(payloads, version="v9.9.9"))
    archive, _ = _build_archive(
        tmp_path / "custom.tar.gz",
        payloads,
        manifest=manifest,
    )
    target = tmp_path / "resources" / "virosync"
    _expect_rejected(
        lambda: ViroSyncDatabaseManager.setup_database(
            database_path=str(target),
            source=str(archive),
            force=True,
            archive_sha256=None,
            manifest_sha256=None,
            full=False,
        ),
        keywords=("pin", "sha", "digest", "checksum"),
    )
    assert not (target / "DB_VERSION").exists()


def test_fast_verification_skips_semantic_tools_and_returns_manifest_identity(
    tmp_path: Path,
) -> None:
    payloads = _payloads()
    manifest = _manifest_bytes(_manifest_document(payloads))
    target = _write_tree(tmp_path / "database", payloads, manifest=manifest)

    def _must_not_run(*_args: object, **_kwargs: object) -> Mapping[str, int]:
        raise AssertionError("fast verification invoked semantic tools")

    result = _verify(
        target,
        manifest,
        full=False,
        semantic_runner=_must_not_run,
    )
    assert _manifest_digest(result) == _sha256_bytes(manifest)


def test_full_verification_uses_injected_semantic_runner(
    tmp_path: Path,
) -> None:
    payloads = _payloads()
    manifest = _manifest_bytes(_manifest_document(payloads))
    target = _write_tree(tmp_path / "database", payloads, manifest=manifest)
    runner = SemanticProbe()

    result = _verify(
        target,
        manifest,
        full=True,
        semantic_runner=runner,
    )

    assert len(runner.calls) == 1
    assert _manifest_digest(result) == _sha256_bytes(manifest)


def test_full_verification_rejects_semantic_count_mismatch(
    tmp_path: Path,
) -> None:
    payloads = _payloads()
    manifest = _manifest_bytes(_manifest_document(payloads))
    target = _write_tree(tmp_path / "database", payloads, manifest=manifest)
    counts = dict(SEMANTIC_COUNTS)
    counts["hmm_models"] = 99

    _expect_rejected(
        lambda: _verify(
            target,
            manifest,
            full=True,
            semantic_runner=SemanticProbe(counts),
        ),
        keywords=("semantic", "count", "hmm", "profile"),
    )


def test_valid_install_is_relative_immutable_idempotent_and_retains_prior(
    tmp_path: Path,
) -> None:
    target = tmp_path / "resources" / "virosync"
    _prior_tree(target)
    archive, manifest = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest is not None

    installed = _setup(target, archive, manifest)

    assert installed == target
    assert target.is_symlink()
    link = Path(os.readlink(target))
    assert not link.is_absolute()
    resolved = target.resolve(strict=True)
    assert resolved.parent == target.parent
    assert ViroSyncDatabaseManager.get_database_version(target) == VERSION
    writable = [
        path.relative_to(resolved).as_posix()
        for path in resolved.rglob("*")
        if path.is_file() and path.stat().st_mode & stat.S_IWUSR
    ]
    assert writable == []
    retained = [
        path
        for path in target.parent.rglob("prior-sentinel.txt")
        if not path.is_relative_to(resolved)
    ]
    assert retained, "successful migration discarded the prior real directory"

    before = _snapshot(target)
    second = _setup(
        target,
        archive,
        manifest,
        force=False,
        full=False,
    )
    assert second == target
    assert target.resolve(strict=True) == resolved
    assert _snapshot(target) == before


def test_rehardening_preserves_payload_ctimes_and_verified_receipt(
    tmp_path: Path,
) -> None:
    target = tmp_path / "resources" / "virosync"
    archive, manifest_bytes = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest_bytes is not None
    _setup(target, archive, manifest_bytes)
    candidate = target.resolve(strict=True)
    manifest = resource_manifest.load_resource_manifest(
        target,
        expected_version=VERSION,
        expected_manifest_sha256=_sha256_bytes(manifest_bytes),
    )
    before = {
        relative: (candidate / relative).lstat().st_ctime_ns
        for relative in CORE_PAYLOADS
    }

    resource_installer._make_tree_immutable(candidate)

    after = {
        relative: (candidate / relative).lstat().st_ctime_ns
        for relative in CORE_PAYLOADS
    }
    assert after == before
    assert resource_installer.verified_install_receipt(candidate, manifest)


def test_valid_receipt_reuse_preserves_install_metadata_bytes_and_stats(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "resources" / "virosync"
    archive, manifest_bytes = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest_bytes is not None
    _setup(target, archive, manifest_bytes)
    metadata_path = target.resolve(strict=True) / "DB_METADATA.json"
    before_bytes = metadata_path.read_bytes()
    before_stat = metadata_path.stat()
    events: list[str] = []
    _record_reuse_checks(monkeypatch, events)

    _setup(target, archive, manifest_bytes, force=False, full=True)

    after_stat = metadata_path.stat()
    assert events == ["receipt"]
    assert metadata_path.read_bytes() == before_bytes
    assert after_stat.st_mtime_ns == before_stat.st_mtime_ns
    assert after_stat.st_ctime_ns == before_stat.st_ctime_ns


def test_missing_receipt_is_verified_before_metadata_is_rebuilt(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "resources" / "virosync"
    archive, manifest_bytes = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest_bytes is not None
    _setup(target, archive, manifest_bytes)
    candidate = target.resolve(strict=True)
    manifest = resource_manifest.load_resource_manifest(target)
    candidate.chmod(0o755)
    (candidate / "DB_METADATA.json").unlink()
    candidate.chmod(0o555)
    events: list[str] = []
    _record_reuse_checks(monkeypatch, events)

    _setup(target, archive, manifest_bytes, force=False, full=True)

    assert events == ["receipt", "verify"]
    assert resource_installer.verified_install_receipt(candidate, manifest)


def test_schema2_runtime_bundle_installs_and_reuses_its_full_receipt(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    archive, manifest_bytes, manifest = _build_schema2_archive(tmp_path, "runtime")
    target = tmp_path / "resources" / "virosync"
    semantic_runner = SemanticProbe(manifest.semantic_counts)

    installed = _setup(
        target,
        archive,
        manifest_bytes,
        semantic_runner=semantic_runner,
    )

    candidate = installed.resolve(strict=True)
    expected_payloads = list(resource_manifest.RUNTIME_RESOURCE_FILES)
    installed_files = {
        path.relative_to(candidate).as_posix()
        for path in candidate.rglob("*")
        if path.is_file()
    }
    assert installed_files == {
        *expected_payloads,
        resource_manifest.RESOURCE_MANIFEST_NAME,
        "DB_METADATA.json",
    }
    metadata = json.loads(
        (candidate / "DB_METADATA.json").read_text(encoding="utf-8")
    )
    assert metadata["required_files"] == expected_payloads
    assert set(metadata["verified_files"]) == set(expected_payloads)
    assert resource_installer.verified_install_receipt(
        candidate,
        manifest,
        expected_archive_sha256=_sha256_file(archive),
    )

    events: list[str] = []
    _record_reuse_checks(monkeypatch, events)
    _setup(
        target,
        archive,
        manifest_bytes,
        force=False,
        full=True,
        semantic_runner=semantic_runner,
    )
    assert events == ["receipt"]


@pytest.mark.parametrize("active_candidate", [False, True])
def test_schema2_source_bundle_is_rejected_as_a_core_install(
    tmp_path: Path,
    active_candidate: bool,
) -> None:
    archive, manifest_bytes, manifest = _build_schema2_archive(tmp_path, "source")
    target = tmp_path / "resources" / "virosync"
    if active_candidate:
        source = ViroSyncDatabaseManager._resolve_core_source(
            str(archive),
            VERSION,
            _sha256_file(archive),
            _sha256_bytes(manifest_bytes),
        )
        candidate = resource_installer._candidate_path(target, source)
        resource_installer.safe_extract_archive(archive, candidate)
        resource_installer._make_tree_immutable(candidate)
        target.symlink_to(candidate.name)

    with pytest.raises(
        resource_installer.ResourceInstallError,
        match="source/repair bundle cannot be installed",
    ):
        _setup(
            target,
            archive,
            manifest_bytes,
            force=False,
            semantic_runner=SemanticProbe(manifest.semantic_counts),
        )

    if active_candidate:
        assert target.is_symlink()
        assert target.resolve(strict=True) == candidate
    else:
        assert not target.exists() and not target.is_symlink()


@pytest.mark.parametrize(
    "mutation",
    ["unexpected_file", "writable_directory", "intermediate_symlink", "wrong_component"],
)
def test_verified_receipt_rejects_structural_mutations(
    tmp_path: Path,
    mutation: str,
) -> None:
    target = tmp_path / "resources" / "virosync"
    archive, manifest_bytes = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest_bytes is not None
    _setup(target, archive, manifest_bytes)
    candidate = target.resolve(strict=True)
    manifest = resource_manifest.load_resource_manifest(target)

    if mutation == "unexpected_file":
        candidate.chmod(0o755)
        unexpected = candidate / "unexpected"
        unexpected.write_text("not authenticated\n", encoding="utf-8")
        unexpected.chmod(0o444)
        candidate.chmod(0o555)
    elif mutation == "writable_directory":
        (candidate / "models").chmod(0o755)
    elif mutation == "intermediate_symlink":
        candidate.chmod(0o755)
        moved = candidate / "models-real"
        (candidate / "models").rename(moved)
        (candidate / "models").symlink_to(moved.name, target_is_directory=True)
        candidate.chmod(0o555)
    else:
        metadata_path = candidate / "DB_METADATA.json"
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        metadata["component"] = "not_virosync_core"
        metadata_path.chmod(0o600)
        metadata_path.write_text(
            json.dumps(metadata, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        metadata_path.chmod(0o444)

    assert not resource_installer.verified_install_receipt(candidate, manifest)


@pytest.mark.parametrize(
    "relative",
    ["DB_METADATA.json", "models/og_marker_name_map.tsv"],
)
def test_stale_receipt_reuse_rejects_hardlinked_files_before_metadata_rewrite(
    tmp_path: Path,
    relative: str,
) -> None:
    target = tmp_path / "resources" / "virosync"
    archive, manifest_bytes = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest_bytes is not None
    _setup(target, archive, manifest_bytes)
    candidate = target.resolve(strict=True)
    metadata_path = candidate / "DB_METADATA.json"
    metadata_before = metadata_path.read_bytes()
    linked = candidate / relative
    outside = tmp_path / f"outside-{linked.name}"
    os.link(linked, outside)

    with pytest.raises(
        resource_installer.ResourceInstallError,
        match="multiply linked file",
    ):
        _setup(target, archive, manifest_bytes, force=False, full=False)

    assert candidate.stat().st_mode & 0o222 == 0
    assert metadata_path.read_bytes() == metadata_before
    assert outside.read_bytes() == linked.read_bytes()


def test_schema3_identity_and_reuse_reject_same_size_installed_payload_corruption(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from virosync.orchestration._flows.single_genome.orchestrator import (
        _enabled_resource_identities,
    )

    target = tmp_path / "resources" / "virosync"
    archive, manifest_bytes = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest_bytes is not None
    _setup(target, archive, manifest_bytes)
    candidate = target.resolve(strict=True)
    manifest = resource_manifest.load_resource_manifest(target)
    payload = candidate / "models/og_marker_name_map.tsv"
    assert _enabled_resource_identities({"hmm_database": candidate / "models/combined.hmm"})
    original = payload.read_bytes()
    corrupted = bytes([original[0] ^ 1]) + original[1:]
    assert len(corrupted) == len(original)
    payload.chmod(0o644)
    payload.write_bytes(corrupted)
    payload.chmod(0o444)

    assert not resource_installer.verified_install_receipt(candidate, manifest)
    with pytest.raises(resource_manifest.ResourceManifestError, match="SHA-256 mismatch"):
        _enabled_resource_identities(
            {"hmm_database": candidate / "models/combined.hmm"}
        )
    events: list[str] = []
    _record_reuse_checks(monkeypatch, events)
    with pytest.raises(Exception, match="SHA-256 mismatch"):
        _setup(target, archive, manifest_bytes, force=False, full=False)
    assert events == ["receipt", "verify"]
    assert target.resolve(strict=True) == candidate
    assert payload.read_bytes() == corrupted


def test_fast_identity_and_reuse_reject_same_size_manifest_mutation(
    tmp_path: Path,
) -> None:
    target = tmp_path / "resources" / "virosync"
    archive, manifest_bytes = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest_bytes is not None
    _setup(target, archive, manifest_bytes)
    candidate = target.resolve(strict=True)
    manifest_path = candidate / "RESOURCE_MANIFEST.json"
    original = manifest_path.read_bytes()
    newline = original.index(b"\n")
    mutated = original[:newline] + b" " + original[newline + 1 :]
    assert len(mutated) == len(original)
    json.loads(mutated)
    manifest_path.chmod(0o644)
    manifest_path.write_bytes(mutated)
    manifest_path.chmod(0o444)
    source = ViroSyncDatabaseManager._resolve_core_source(
        str(archive),
        VERSION,
        _sha256_file(archive),
        _sha256_bytes(manifest_bytes),
    )

    assert ViroSyncDatabaseManager._trusted_database_root(
        target,
        source=source,
    ) is None
    with pytest.raises(Exception, match="manifest SHA-256 mismatch"):
        _setup(target, archive, manifest_bytes, force=False, full=False)


def test_resolve_config_paths_uses_only_a_verified_install_receipt(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "resources" / "virosync"
    archive, manifest_bytes = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest_bytes is not None
    _setup(target, archive, manifest_bytes)
    setup_calls: list[dict[str, object]] = []

    def _reject_setup(cls, **kwargs: object) -> Path:
        setup_calls.append(dict(kwargs))
        raise UnexpectedExternalActivity("verified install triggered setup")

    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "setup_database",
        classmethod(_reject_setup),
    )
    resolved = _REAL_RESOLVE_CONFIG_PATHS(
        ViroSyncDatabaseManager,
        {
            "database_root": str(target),
            "core_resources_url": str(archive),
            "core_resources_version": VERSION,
            "core_resources_sha256": _sha256_file(archive),
            "core_resources_manifest_sha256": _sha256_bytes(manifest_bytes),
            "hmm_db": None,
            "marker_db": None,
            "gene_taxonomy_faa_db": None,
        }
    )

    assert setup_calls == []
    candidate = target.resolve(strict=True)
    assert resolved["hmm_db"] == str(candidate / "models/combined.hmm")
    assert resolved["marker_db"] == str(candidate / "marker/marker.dmnd")
    assert resolved["gene_taxonomy_faa_db"] == str(
        candidate / "genomes/combined_proteome.dmnd"
    )


def test_resolve_config_paths_rejects_a_stale_install_receipt(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "resources" / "virosync"
    archive, manifest_bytes = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest_bytes is not None
    _setup(target, archive, manifest_bytes)
    candidate = target.resolve(strict=True)
    candidate.chmod(0o755)
    (candidate / "DB_METADATA.json").unlink()
    candidate.chmod(0o555)
    setup_calls: list[dict[str, object]] = []

    def _record_setup(cls, **kwargs: object) -> Path:
        setup_calls.append(dict(kwargs))
        return target

    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "setup_database",
        classmethod(_record_setup),
    )
    _REAL_RESOLVE_CONFIG_PATHS(
        ViroSyncDatabaseManager,
        {
            "database_root": str(target),
            "core_resources_url": str(archive),
            "core_resources_version": VERSION,
            "core_resources_sha256": _sha256_file(archive),
            "core_resources_manifest_sha256": _sha256_bytes(manifest_bytes),
            "hmm_db": None,
            "marker_db": None,
            "gene_taxonomy_faa_db": None,
        }
    )

    assert len(setup_calls) == 1
    assert setup_calls[0]["full"] is True


def test_resolve_config_paths_rejects_conflicting_pinned_source_before_reuse(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "resources" / "virosync"
    archive, manifest_bytes = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest_bytes is not None
    archive_sha256 = _sha256_file(archive)
    manifest_sha256 = _sha256_bytes(manifest_bytes)
    _setup(target, archive, manifest_bytes)
    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "DATABASE_SOURCES",
        [
            {
                "version": VERSION,
                "source": str(archive),
                "filename": archive.name,
                "archive_sha256": archive_sha256,
                "manifest_sha256": manifest_sha256,
            }
        ],
    )

    with pytest.raises(
        resource_installer.ResourceInstallError,
        match="conflicts with its source record",
    ):
        _REAL_RESOLVE_CONFIG_PATHS(
            ViroSyncDatabaseManager,
            {
                "database_root": str(target),
                "core_resources_url": str(archive),
                "core_resources_version": VERSION,
                "core_resources_sha256": "0" * 64,
                "core_resources_manifest_sha256": manifest_sha256,
                "hmm_db": None,
                "marker_db": None,
                "gene_taxonomy_faa_db": None,
            },
        )


def test_resolve_config_paths_does_not_trust_an_absolute_active_pointer(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "resources" / "virosync"
    archive, manifest_bytes = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest_bytes is not None
    _setup(target, archive, manifest_bytes)
    candidate = target.resolve(strict=True)
    target.unlink()
    target.symlink_to(candidate, target_is_directory=True)
    setup_calls: list[dict[str, object]] = []

    def _record_setup(cls, **kwargs: object) -> Path:
        setup_calls.append(dict(kwargs))
        return target

    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "setup_database",
        classmethod(_record_setup),
    )
    _REAL_RESOLVE_CONFIG_PATHS(
        ViroSyncDatabaseManager,
        {
            "database_root": str(target),
            "core_resources_url": str(archive),
            "core_resources_version": VERSION,
            "core_resources_sha256": _sha256_file(archive),
            "core_resources_manifest_sha256": _sha256_bytes(manifest_bytes),
            "hmm_db": None,
            "marker_db": None,
            "gene_taxonomy_faa_db": None,
        },
    )

    assert len(setup_calls) == 1


def test_existing_default_tree_must_match_the_pinned_manifest(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    archive, pinned_manifest = _build_archive(tmp_path / "pinned.tar.gz")
    assert pinned_manifest is not None
    payloads = _payloads()
    self_signed_manifest = _manifest_bytes(_manifest_document(payloads)) + b"\n"
    target = _write_tree(
        tmp_path / "resources" / "virosync",
        payloads,
        manifest=self_signed_manifest,
    )
    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "DATABASE_SOURCES",
        [
            {
                "version": VERSION,
                "source": str(archive),
                "filename": archive.name,
                "archive_sha256": _sha256_file(archive),
                "manifest_sha256": _sha256_bytes(pinned_manifest),
            }
        ],
    )

    installed = ViroSyncDatabaseManager.setup_database(
        database_path=str(target),
        source=None,
        force=False,
        full=True,
        semantic_runner=SemanticProbe(),
    )

    assert Path(installed).is_symlink()
    result = ViroSyncDatabaseManager.verify_database(
        installed,
        expected_version=VERSION,
        manifest_sha256=_sha256_bytes(pinned_manifest),
        full=False,
    )
    assert _manifest_digest(result) == _sha256_bytes(pinned_manifest)


@pytest.mark.parametrize("case", ["zero_byte", "wrong_version"])
def test_implicit_setup_rejects_invalid_existing_before_download(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    case: str,
) -> None:
    archive, pinned_manifest = _build_archive(tmp_path / "pinned.tar.gz")
    assert pinned_manifest is not None
    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "DATABASE_SOURCES",
        [
            {
                "version": VERSION,
                "source": str(archive),
                "filename": archive.name,
                "archive_sha256": _sha256_file(archive),
                "manifest_sha256": _sha256_bytes(pinned_manifest),
            }
        ],
    )
    if case == "wrong_version":
        payloads = _payloads(PRIOR_VERSION)
        manifest = _manifest_bytes(
            _manifest_document(payloads, version=PRIOR_VERSION)
        )
    else:
        payloads = _payloads()
        manifest = _manifest_bytes(_manifest_document(payloads))
    target = _write_tree(
        tmp_path / "resources" / "virosync",
        payloads,
        manifest=manifest,
    )
    if case == "zero_byte":
        (target / "marker/marker.dmnd").write_bytes(b"")
    before = _snapshot(target)
    copy_calls: list[str] = []

    def _reject_copy(cls, source: str, _destination: Path) -> None:
        copy_calls.append(source)
        raise UnexpectedExternalActivity("invalid existing tree triggered a download")

    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "_copy_or_download_archive",
        classmethod(_reject_copy),
    )
    with pytest.raises(Exception, match="empty|version|size"):
        ViroSyncDatabaseManager.setup_database(
            database_path=str(target),
            source=None,
            force=False,
            full=True,
            semantic_runner=SemanticProbe(),
        )

    assert copy_calls == []
    assert not target.is_symlink()
    assert _snapshot(target) == before


def test_existing_valid_real_tree_is_migrated_to_an_immutable_pointer(
    tmp_path: Path,
) -> None:
    payloads = _payloads()
    manifest = _manifest_bytes(_manifest_document(payloads))
    target = _write_tree(
        tmp_path / "resources" / "virosync",
        payloads,
        manifest=manifest,
    )
    before = _snapshot(target)
    archive, embedded_manifest = _build_archive(
        tmp_path / "valid.tar.gz",
        payloads,
        manifest=manifest,
    )
    assert embedded_manifest is not None

    installed = _setup(
        target,
        archive,
        embedded_manifest,
        force=False,
    )

    assert installed.is_symlink()
    candidate = installed.resolve(strict=True)
    assert candidate.parent == installed.parent
    assert os.readlink(installed) == candidate.name
    assert all(
        path.stat().st_mode & 0o222 == 0
        for path in (candidate, *candidate.rglob("*"))
    )
    retained = [
        path
        for path in target.parent.iterdir()
        if path.name.startswith(f".{target.name}.legacy-") and path.is_dir()
    ]
    assert len(retained) == 1
    assert _snapshot(retained[0]) == before


def test_normal_setup_recovers_pending_journal_before_fast_return(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "resources" / "virosync"
    _prior_tree(target)
    archive, manifest = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest is not None
    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "DATABASE_SOURCES",
        [
            {
                "version": VERSION,
                "source": str(archive),
                "filename": archive.name,
                "archive_sha256": _sha256_file(archive),
                "manifest_sha256": _sha256_bytes(manifest),
            }
        ],
    )

    def _crash_after_activation(phase: str) -> None:
        if phase == "after_pointer_activate":
            raise InjectedCrash("synthetic crash after pointer activation")

    with pytest.raises(InjectedCrash):
        _setup(
            target,
            archive,
            manifest,
            fault_injector=_crash_after_activation,
        )
    journal = target.parent / f".{target.name}.resource-install-journal.json"
    assert journal.is_file()

    installed = ViroSyncDatabaseManager.setup_database(
        database_path=str(target),
        source=None,
        force=False,
        full=True,
        semantic_runner=SemanticProbe(),
    )

    assert Path(installed).is_symlink()
    assert not journal.exists()
    assert ViroSyncDatabaseManager.get_database_version(target) == VERSION


def test_failed_directory_prior_recovery_never_removes_active_candidate(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "resources" / "virosync"
    _prior_tree(target)
    archive, manifest = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest is not None

    def _crash_after_activation(phase: str) -> None:
        if phase == "after_pointer_activate":
            raise InjectedCrash("synthetic crash after pointer activation")

    with pytest.raises(InjectedCrash):
        _setup(
            target,
            archive,
            manifest,
            fault_injector=_crash_after_activation,
        )
    assert target.is_symlink()
    assert ViroSyncDatabaseManager.get_database_version(target) == VERSION

    real_replace = resource_installer.os.replace

    def _fail_restore(source: str | Path, destination: str | Path) -> None:
        if Path(destination) == target:
            raise OSError("synthetic recovery replacement failure")
        real_replace(source, destination)

    monkeypatch.setattr(resource_installer.os, "replace", _fail_restore)
    with pytest.raises(OSError, match="synthetic recovery replacement failure"):
        resource_installer.recover_pending_install(target)

    assert target.is_symlink()
    assert ViroSyncDatabaseManager.get_database_version(target) == VERSION


def test_candidate_data_and_modes_are_fsynced_before_activation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    archive, manifest = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest is not None
    target = tmp_path / "resources" / "virosync"
    synced_modes: dict[Path, int] = {}
    real_fsync = resource_installer.os.fsync

    def _record_fsync(descriptor: int) -> None:
        try:
            path = Path(os.readlink(f"/proc/self/fd/{descriptor}"))
            synced_modes[path] = os.fstat(descriptor).st_mode
        except OSError:
            pass
        real_fsync(descriptor)

    monkeypatch.setattr(resource_installer.os, "fsync", _record_fsync)

    def _assert_durable_before_activation(phase: str) -> None:
        if phase != "after_candidate_promote":
            return
        candidate = target.parent / f"{target.name}-{VERSION}-{_sha256_bytes(manifest)[:16]}"
        files = {
            path
            for path in candidate.rglob("*")
            if path.is_file() and not path.is_symlink()
        }
        assert files
        assert files <= set(synced_modes)
        assert all(synced_modes[path] & 0o222 == 0 for path in files)
        assert not target.exists() and not target.is_symlink()

    _setup(
        target,
        archive,
        manifest,
        fault_injector=_assert_durable_before_activation,
    )

    candidate = target.resolve(strict=True)
    installed_files = {
        path
        for path in candidate.rglob("*")
        if path.is_file() and not path.is_symlink()
    }
    assert installed_files
    assert installed_files <= set(synced_modes)
    assert all(synced_modes[path] & 0o222 == 0 for path in installed_files)


def test_reused_candidate_is_fsynced_after_rehardening(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    archive, manifest = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest is not None
    target = tmp_path / "resources" / "virosync"
    _setup(target, archive, manifest)
    candidate = target.resolve(strict=True)
    payload = candidate / "DB_VERSION"
    candidate.chmod(0o755)
    payload.chmod(0o644)
    synced_modes: dict[Path, int] = {}
    real_fsync = resource_installer.os.fsync

    def _record_fsync(descriptor: int) -> None:
        try:
            path = Path(os.readlink(f"/proc/self/fd/{descriptor}"))
            synced_modes[path] = os.fstat(descriptor).st_mode
        except OSError:
            pass
        real_fsync(descriptor)

    monkeypatch.setattr(resource_installer.os, "fsync", _record_fsync)
    _setup(target, archive, manifest, force=False, full=False)

    assert payload in synced_modes
    assert synced_modes[payload] & 0o222 == 0
    assert payload.stat().st_mode & 0o222 == 0
    assert candidate.stat().st_mode & 0o222 == 0


def test_fast_setup_cannot_activate_a_bundle(tmp_path: Path) -> None:
    archive, manifest = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest is not None
    target = tmp_path / "resources" / "virosync"

    with pytest.raises(
        resource_installer.ResourceInstallError,
        match="Full semantic validation.*activation",
    ):
        _setup(target, archive, manifest, force=True, full=False)

    assert not target.exists() and not target.is_symlink()


@pytest.mark.parametrize(
    "alias_kind",
    ["target", "candidate", "journal", "lock"],
)
def test_recovery_journal_cannot_alias_control_paths(
    tmp_path: Path,
    alias_kind: str,
) -> None:
    parent = tmp_path / "resources"
    candidate = parent / "virosync-v1.0.6-candidate"
    candidate.mkdir(parents=True)
    (candidate / "sentinel").write_text("valid\n", encoding="utf-8")
    target = parent / "virosync"
    target.symlink_to(candidate.name)
    journal = parent / ".virosync.resource-install-journal.json"
    journal.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "state": "activated",
                "target": target.name,
                "candidate": candidate.name,
                "temporary_pointer": {
                    "target": target.name,
                    "candidate": candidate.name,
                    "journal": journal.name,
                    "lock": f".{target.name}.resource-install.lock",
                }[alias_kind],
                "prior_kind": "missing",
                "prior": None,
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(resource_installer.ResourceInstallError, match="temporary"):
        resource_installer.recover_pending_install(target)

    assert target.is_symlink()
    assert target.resolve(strict=True) == candidate


@pytest.mark.parametrize("kind", ["symlink", "hardlink", "fifo"])
def test_sibling_install_lock_rejects_nonregular_or_aliased_paths(
    tmp_path: Path,
    kind: str,
) -> None:
    archive, manifest = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest is not None
    target = tmp_path / "resources" / "virosync"
    target.parent.mkdir(parents=True)
    outside = tmp_path / "outside-lock"
    outside.write_text("do not lock through this path\n", encoding="utf-8")
    lock = target.parent / f".{target.name}.resource-install.lock"
    if kind == "symlink":
        lock.symlink_to(outside)
    elif kind == "hardlink":
        os.link(outside, lock)
    else:
        os.mkfifo(lock)

    with pytest.raises(resource_installer.ResourceInstallError, match="lock"):
        _setup(target, archive, manifest)

    assert lock.exists() or lock.is_symlink()
    assert outside.read_text(encoding="utf-8") == "do not lock through this path\n"


@pytest.mark.parametrize("kind", ["symlink", "hardlink"])
def test_recovery_rejects_aliased_journal_without_touching_active(
    tmp_path: Path,
    kind: str,
) -> None:
    parent = tmp_path / "resources"
    candidate = parent / "virosync-v1.0.6-0123456789abcdef"
    candidate.mkdir(parents=True)
    target = parent / "virosync"
    target.symlink_to(candidate.name)
    outside = tmp_path / "outside-journal.json"
    outside.write_text("{}\n", encoding="utf-8")
    journal = parent / ".virosync.resource-install-journal.json"
    if kind == "symlink":
        journal.symlink_to(outside)
    else:
        os.link(outside, journal)

    with pytest.raises(resource_installer.ResourceInstallError, match="journal"):
        resource_installer.recover_pending_install(target)

    assert target.is_symlink()
    assert target.resolve(strict=True) == candidate
    assert outside.read_text(encoding="utf-8") == "{}\n"


def test_failed_staged_validation_preserves_prior_real_directory(
    tmp_path: Path,
) -> None:
    target = tmp_path / "resources" / "virosync"
    before = _prior_tree(target)
    payloads = _payloads()
    payloads.pop("marker/marker.dmnd")
    archive, manifest = _build_archive(tmp_path / "invalid.tar.gz", payloads)
    assert manifest is not None

    with pytest.raises(Exception):
        _setup(target, archive, manifest)

    assert not target.is_symlink()
    assert _snapshot(target) == before


@pytest.mark.parametrize("phase", INSTALL_FAULT_PHASES)
def test_every_install_fault_phase_preserves_or_recovers_a_valid_active_tree(
    tmp_path: Path,
    phase: str,
) -> None:
    target = tmp_path / phase / "resources" / "virosync"
    prior_snapshot = _prior_tree(target)
    archive, manifest = _build_archive(tmp_path / phase / "valid.tar.gz")
    assert manifest is not None
    seen: list[str] = []

    def _fail_at_phase(current: str) -> None:
        seen.append(current)
        if current == phase:
            raise RuntimeError(f"synthetic installer fault at {phase}")

    with pytest.raises(RuntimeError, match=re.escape(phase)):
        _setup(
            target,
            archive,
            manifest,
            fault_injector=_fail_at_phase,
        )
    assert phase in seen

    if phase in PRE_ACTIVATION_PHASES:
        assert target.exists() or target.is_symlink()
        assert _snapshot(target) == prior_snapshot
        assert ViroSyncDatabaseManager.get_database_version(target) == PRIOR_VERSION
    else:
        assert target.exists() or target.is_symlink()
        visible_version = ViroSyncDatabaseManager.get_database_version(target)
        assert visible_version in {
            PRIOR_VERSION,
            VERSION,
        }
        if visible_version == VERSION:
            _verify(
                target,
                manifest,
                full=True,
                semantic_runner=SemanticProbe(),
            )
        else:
            assert _snapshot(target) == prior_snapshot

    recovered = _setup(target, archive, manifest)
    assert recovered == target
    assert target.is_symlink()
    assert ViroSyncDatabaseManager.get_database_version(target) == VERSION
    retained_prior = [
        path
        for path in target.parent.rglob("prior-sentinel.txt")
        if not path.is_relative_to(target.resolve(strict=True))
    ]
    assert retained_prior, f"fault phase {phase} discarded the prior install"


def test_recovery_journal_repairs_an_abrupt_crash_after_prior_move(
    tmp_path: Path,
) -> None:
    target = tmp_path / "resources" / "virosync"
    _prior_tree(target)
    archive, manifest = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest is not None

    def _crash_after_prior_move(phase: str) -> None:
        if phase == "after_prior_move":
            raise InjectedCrash("synthetic crash after prior move")

    with pytest.raises(InjectedCrash):
        _setup(
            target,
            archive,
            manifest,
            fault_injector=_crash_after_prior_move,
        )

    recovered = _setup(target, archive, manifest)
    assert recovered == target
    assert target.is_symlink()
    assert ViroSyncDatabaseManager.get_database_version(target) == VERSION
    retained_prior = [
        path
        for path in target.parent.rglob("prior-sentinel.txt")
        if not path.is_relative_to(target.resolve(strict=True))
    ]
    assert retained_prior, "journal recovery lost the migrated real directory"


@pytest.mark.parametrize("phase", ["after_journal_write", "after_pointer_prepare"])
def test_recovery_repairs_abrupt_crashes_before_the_prior_move(
    tmp_path: Path,
    phase: str,
) -> None:
    target = tmp_path / phase / "resources" / "virosync"
    prior_snapshot = _prior_tree(target)
    archive, manifest = _build_archive(tmp_path / phase / "valid.tar.gz")
    assert manifest is not None

    def _crash(current: str) -> None:
        if current == phase:
            raise InjectedCrash(f"synthetic crash at {phase}")

    with pytest.raises(InjectedCrash, match=phase):
        _setup(target, archive, manifest, fault_injector=_crash)

    journal = target.parent / f".{target.name}.resource-install-journal.json"
    assert journal.is_file()
    assert _snapshot(target) == prior_snapshot

    recovered = _setup(target, archive, manifest, force=False)

    assert recovered.is_symlink()
    assert not journal.exists()
    assert ViroSyncDatabaseManager.get_database_version(target) == VERSION


@pytest.mark.parametrize(
    "phase",
    [
        "after_journal_write",
        "after_pointer_prepare",
        "after_prior_move",
        "after_pointer_activate",
    ],
)
def test_first_install_recovers_or_rolls_forward_after_abrupt_crash(
    tmp_path: Path,
    phase: str,
) -> None:
    target = tmp_path / phase / "resources" / "virosync"
    archive, manifest = _build_archive(tmp_path / f"{phase}.tar.gz")
    assert manifest is not None

    def _crash(current: str) -> None:
        if current == phase:
            raise InjectedCrash(f"synthetic first-install crash at {phase}")

    with pytest.raises(InjectedCrash, match=phase):
        _setup(target, archive, manifest, fault_injector=_crash)

    journal = target.parent / f".{target.name}.resource-install-journal.json"
    assert journal.is_file()
    if phase == "after_pointer_activate":
        assert target.is_symlink()
        activated = target.resolve(strict=True)
    else:
        assert not target.exists() and not target.is_symlink()
        activated = None

    recovered = _setup(target, archive, manifest, force=False)

    assert recovered.is_symlink()
    assert not journal.exists()
    assert ViroSyncDatabaseManager.get_database_version(target) == VERSION
    if activated is not None:
        assert target.resolve(strict=True) == activated


def test_recovery_leaves_active_candidate_visible_when_retained_prior_is_missing(
    tmp_path: Path,
) -> None:
    target = tmp_path / "resources" / "virosync"
    _prior_tree(target)
    archive, manifest = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest is not None

    def _crash_after_activation(phase: str) -> None:
        if phase == "after_pointer_activate":
            raise InjectedCrash("synthetic crash after pointer activation")

    with pytest.raises(InjectedCrash):
        _setup(
            target,
            archive,
            manifest,
            fault_injector=_crash_after_activation,
        )
    assert target.is_symlink()
    assert ViroSyncDatabaseManager.get_database_version(target) == VERSION
    retained_sentinel = next(
        path
        for path in target.parent.rglob("prior-sentinel.txt")
        if not path.is_relative_to(target.resolve(strict=True))
    )
    retained = retained_sentinel.parent
    quarantine = tmp_path / "quarantined-prior"
    retained.rename(quarantine)

    with pytest.raises(Exception, match="Retained resource directory is unavailable"):
        _setup(target, archive, manifest)

    assert target.is_symlink()
    assert ViroSyncDatabaseManager.get_database_version(target) == VERSION


def test_reused_candidate_is_made_immutable_again_before_activation(
    tmp_path: Path,
) -> None:
    target = tmp_path / "resources" / "virosync"
    archive, manifest = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest is not None
    _setup(target, archive, manifest)
    candidate = target.resolve(strict=True)
    payload = candidate / "DB_VERSION"
    candidate.chmod(0o755)
    payload.chmod(0o644)
    assert payload.stat().st_mode & stat.S_IWUSR
    parsed_manifest = resource_manifest.load_resource_manifest(target)
    assert not resource_installer.verified_install_receipt(candidate, parsed_manifest)

    _setup(target, archive, manifest)

    assert target.resolve(strict=True) == candidate
    assert not payload.stat().st_mode & stat.S_IWUSR
    assert candidate.stat().st_mode & stat.S_IWUSR == 0
    assert resource_installer.verified_install_receipt(candidate, parsed_manifest)


def test_concurrent_setup_serializes_and_leaves_one_valid_active_bundle(
    tmp_path: Path,
) -> None:
    archive, manifest = _build_archive(tmp_path / "valid.tar.gz")
    assert manifest is not None
    target = tmp_path / "resources" / "virosync"
    runner = SemanticProbe()
    runner.delay_seconds = 0.05
    start = threading.Barrier(3)
    results: list[Path] = []
    failures: list[BaseException] = []

    def _worker() -> None:
        start.wait(timeout=5)
        try:
            result = _setup(target, archive, manifest, semantic_runner=runner)
        except BaseException as exc:  # preserve both thread failures for assertion
            failures.append(exc)
        else:
            results.append(result)

    threads = [threading.Thread(target=_worker, daemon=True) for _ in range(2)]
    for thread in threads:
        thread.start()
    start.wait(timeout=5)
    for thread in threads:
        thread.join(timeout=10)

    assert all(not thread.is_alive() for thread in threads)
    assert failures == []
    assert results == [target, target]
    assert runner.max_active == 1
    assert target.is_symlink()
    assert ViroSyncDatabaseManager.get_database_version(target) == VERSION
    _verify(target, manifest, full=True, semantic_runner=SemanticProbe())
    sibling_locks = [
        path
        for path in target.parent.iterdir()
        if path.is_file() and path.name.endswith(".lock")
    ]
    assert sibling_locks, "setup did not leave a stable sibling lock file"
