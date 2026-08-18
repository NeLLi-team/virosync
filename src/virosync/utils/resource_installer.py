"""Authenticated, transactional installation of ViroSync core resources."""

from __future__ import annotations

from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import logging
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import subprocess
import tarfile
import tempfile
from typing import Any
from uuid import uuid4

from virosync.utils.ssl_env import sanitized_ssl_env
from virosync.utils.resource_manifest import RESOURCE_MANIFEST_NAME, load_resource_manifest


logger = logging.getLogger(__name__)

SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
DOWNLOAD_PERCENT_PATTERN = re.compile(rb"(?<!\d)(\d{1,3})(?:\.\d+)?%")
ARCHIVE_ROOT = "virosync"
INSTALL_METADATA_MAX_BYTES = 1024 * 1024
DOWNLOAD_ERROR_TAIL_BYTES = 64 * 1024

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


class ResourceInstallError(RuntimeError):
    """The core-resource transaction could not be completed safely."""


class ArchiveSafetyError(ResourceInstallError):
    """An archive failed member preflight before extraction."""


@dataclass(frozen=True)
class ResourceSource:
    """Externally authenticated identity for one core-resource archive."""

    version: str
    source: str
    filename: str
    archive_sha256: str
    manifest_sha256: str


FaultInjector = Callable[[str], None]
TreeVerifier = Callable[..., Any]
ArchiveCopier = Callable[[str, Path], None]


def _invoke_fault(fault_injector: FaultInjector | None, phase: str) -> None:
    if phase not in INSTALL_FAULT_PHASES:
        raise AssertionError(f"Unknown resource-install fault phase: {phase}")
    if fault_injector is not None:
        fault_injector(phase)


def normalize_sha256(value: str, label: str) -> str:
    """Return a strict lowercase SHA-256 digest or raise."""
    normalized = str(value).strip().lower()
    if SHA256_PATTERN.fullmatch(normalized) is None:
        raise ResourceInstallError(f"{label} must be a 64-character SHA-256 digest")
    return normalized


def sha256_file(path: Path) -> str:
    """Hash a regular file without loading it into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download_error_detail(stderr: object) -> str:
    """Return a short URL-redacted download diagnostic."""
    if isinstance(stderr, bytes):
        text = stderr.decode("utf-8", errors="replace")
    else:
        text = str(stderr or "")
    lines = [line.strip() for line in text.replace("\r", "\n").splitlines() if line.strip()]
    diagnostic_lines = [
        line
        for line in lines
        if re.search(
            r"error|failed|unable|resolv|connect|certificate|timed out|not found|denied|curl:",
            line,
            re.IGNORECASE,
        )
    ]
    selected = diagnostic_lines[-2:] or lines[-1:]
    redacted: list[str] = []
    for line in selected:
        line = re.sub(r"(?:https?|ftp)://\S+", "<url>", line)
        endpoint_match = re.match(r"(Resolving|Connecting to)\b", line, re.IGNORECASE)
        if endpoint_match:
            outcome = re.search(r"\b(failed|connected)\b.*$", line, re.IGNORECASE)
            line = f"{endpoint_match.group(1)} <host>"
            if outcome:
                line += f": {outcome.group(0)}"
        line = re.sub(
            r"(resolve host:)\s+(?!<url>)\S+",
            r"\1 <host>",
            line,
            flags=re.IGNORECASE,
        )
        redacted.append(line)
    detail = " | ".join(redacted)
    return detail[:500]


def _run_streamed_download(
    command: list[str],
    *,
    env: dict[str, str],
    progress_callback: Callable[[float], None],
    popen_factory: Callable[..., subprocess.Popen] = subprocess.Popen,
) -> subprocess.CompletedProcess:
    """Capture downloader stderr while forwarding parsed percentage updates."""
    process = popen_factory(
        command,
        env=env,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    stderr_tail = bytearray()
    overlap = b""
    last_percent = -1
    stream = process.stderr
    if stream is None:
        raise ResourceInstallError("downloader did not expose its diagnostic stream")
    read_chunk = getattr(stream, "read1", stream.read)
    try:
        while True:
            chunk = read_chunk(4096)
            if not chunk:
                break
            if isinstance(chunk, str):
                chunk = chunk.encode()
            stderr_tail.extend(chunk)
            if len(stderr_tail) > DOWNLOAD_ERROR_TAIL_BYTES:
                del stderr_tail[:-DOWNLOAD_ERROR_TAIL_BYTES]
            scan = overlap + chunk
            for match in DOWNLOAD_PERCENT_PATTERN.finditer(scan):
                percent = int(match.group(1))
                if percent > 100:
                    continue
                if percent > last_percent:
                    progress_callback(float(percent))
                    last_percent = percent
            overlap = scan[-16:]
        returncode = process.wait()
    except BaseException:
        try:
            process.terminate()
            process.wait(timeout=5)
        except (OSError, ProcessLookupError):
            pass
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
        raise
    finally:
        stream.close()
    if returncode == 0 and last_percent < 100:
        progress_callback(100.0)
    return subprocess.CompletedProcess(
        command,
        returncode,
        stdout=b"",
        stderr=bytes(stderr_tail),
    )


def copy_or_download_archive(
    source: str,
    archive_path: Path,
    *,
    ca_bundle: str | None = None,
    command_runner: Callable[..., subprocess.CompletedProcess] = subprocess.run,
    popen_factory: Callable[..., subprocess.Popen] = subprocess.Popen,
    progress_callback: Callable[[float], None] | None = None,
) -> None:
    """Copy a local archive or download one with bounded retries and timeouts."""
    normalized_source = (
        source.replace("file://", "", 1) if source.startswith("file://") else source
    )
    local_path = Path(normalized_source).expanduser()
    if local_path.exists():
        if not local_path.is_file():
            raise ResourceInstallError(f"Archive source is not a file: {local_path}")
        logger.info("Using local archive source: %s", local_path)
        shutil.copy2(local_path, archive_path)
        if progress_callback is not None:
            progress_callback(100.0)
        return

    if not source.startswith(("http://", "https://", "ftp://")):
        raise FileNotFoundError(f"Archive source not found and not a URL: {source}")

    wget_command = [
        "wget",
        "--progress=bar:force:noscroll",
        "--tries=3",
        "--timeout=120",
        "-O",
        str(archive_path),
    ]
    curl_command = [
        "curl",
        "-fL",
        "--retry",
        "3",
        "--retry-delay",
        "5",
        "--connect-timeout",
        "30",
        "--progress-bar",
        "--show-error",
        "-o",
        str(archive_path),
    ]
    if ca_bundle:
        wget_command.extend(["--ca-certificate", ca_bundle])
        curl_command.extend(["--cacert", ca_bundle])
    wget_command.append(source)
    curl_command.append(source)

    errors: list[str] = []
    commands = (("wget", wget_command), ("curl", curl_command))
    last_reported_percent = -1.0

    def report_progress(percent: float) -> None:
        nonlocal last_reported_percent
        if progress_callback is not None and percent > last_reported_percent:
            progress_callback(percent)
            last_reported_percent = percent

    for index, (tool, command) in enumerate(commands):
        if index:
            logger.warning("%s failed, trying %s...", commands[index - 1][0], tool)
        try:
            if progress_callback is not None:
                result = _run_streamed_download(
                    command,
                    env=sanitized_ssl_env(),
                    progress_callback=report_progress,
                    popen_factory=popen_factory,
                )
            else:
                result = command_runner(
                    command,
                    env=sanitized_ssl_env(),
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.PIPE,
                )
        except FileNotFoundError:
            result = subprocess.CompletedProcess(command, 127)
        if (
            result.returncode == 0
            and archive_path.is_file()
            and archive_path.stat().st_size > 0
        ):
            return
        archive_path.unlink(missing_ok=True)
        detail = _download_error_detail(result.stderr)
        suffix = f": {detail}" if detail else ""
        errors.append(f"{tool}: exit code {result.returncode}{suffix}")
    raise ResourceInstallError("archive download failed (" + "; ".join(errors) + ")")


def _archive_member_path(
    member: tarfile.TarInfo,
    archive_root: str,
) -> PurePosixPath | None:
    raw_name = member.name
    if not raw_name or "\\" in raw_name or raw_name.startswith("/"):
        raise ArchiveSafetyError(f"Unsafe absolute or malformed archive path: {raw_name!r}")
    if archive_root == ".":
        if raw_name.rstrip("/") == ".":
            if not member.isdir():
                raise ArchiveSafetyError("Archive root must be a directory entry")
            return None
        if not raw_name.startswith("./"):
            raise ArchiveSafetyError(
                f"Unexpected archive root for {raw_name!r}; expected './'"
            )
        raw_name = raw_name[2:]
    raw_parts = raw_name.rstrip("/").split("/")
    if any(part in {"", ".", ".."} for part in raw_parts):
        raise ArchiveSafetyError(f"Unsafe archive traversal path: {member.name!r}")
    path = PurePosixPath(*raw_parts)
    if archive_root == ".":
        return path
    if path.is_absolute() or not path.parts or path.parts[0] != archive_root:
        raise ArchiveSafetyError(
            f"Unexpected archive root for {member.name!r}; expected {archive_root!r}"
        )
    if len(path.parts) == 1:
        if not member.isdir():
            raise ArchiveSafetyError("Archive root must be a directory entry")
        return None
    return PurePosixPath(*path.parts[1:])


def _preflight_archive(
    handle: tarfile.TarFile,
    archive_root: str,
) -> list[tuple[tarfile.TarInfo, PurePosixPath]]:
    members: list[tuple[tarfile.TarInfo, PurePosixPath]] = []
    seen: set[str] = set()
    for member in handle.getmembers():
        relative = _archive_member_path(member, archive_root)
        if relative is None:
            continue
        normalized = relative.as_posix()
        if normalized in seen:
            raise ArchiveSafetyError(f"Duplicate archive member: {normalized}")
        seen.add(normalized)
        if member.issym() or member.islnk():
            raise ArchiveSafetyError(
                f"Archive links are not permitted in core resources: {normalized}"
            )
        if not (member.isdir() or member.isreg()):
            raise ArchiveSafetyError(
                f"Archive special-file member is not permitted: {normalized}"
            )
        members.append((member, relative))
    if not members:
        raise ArchiveSafetyError("Archive contains no core-resource payloads")
    return members


def _discover_single_archive_root(handle: tarfile.TarFile) -> str:
    members = handle.getmembers()
    dot_roots = [member for member in members if member.name.rstrip("/") == "."]
    if dot_roots:
        if len(dot_roots) != 1 or not dot_roots[0].isdir():
            raise ArchiveSafetyError("Archive root must be one directory entry")
        if any(
            member.name.rstrip("/") != "." and not member.name.startswith("./")
            for member in members
        ):
            raise ArchiveSafetyError("Archive mixes './' and named top-level roots")
        return "."

    roots: set[str] = set()
    for member in members:
        raw_name = member.name
        if not raw_name or "\\" in raw_name or raw_name.startswith("/"):
            raise ArchiveSafetyError(
                f"Unsafe absolute or malformed archive path: {raw_name!r}"
            )
        raw_parts = raw_name.rstrip("/").split("/")
        if any(part in {"", ".", ".."} for part in raw_parts):
            raise ArchiveSafetyError(f"Unsafe archive traversal path: {raw_name!r}")
        roots.add(raw_parts[0])
    if len(roots) != 1:
        raise ArchiveSafetyError(
            f"Archive must have exactly one top-level root; found {sorted(roots)}"
        )
    return roots.pop()


def _safe_extract_archive(
    archive_path: Path,
    target_dir: Path,
    *,
    archive_root: str | None,
    preserve_executable: bool = False,
) -> None:
    archive_path = Path(archive_path)
    target_dir = Path(target_dir)
    if target_dir.exists() or target_dir.is_symlink():
        raise ArchiveSafetyError(f"Extraction target already exists: {target_dir}")

    try:
        with tarfile.open(archive_path, mode="r:*") as handle:
            selected_root = archive_root or _discover_single_archive_root(handle)
            members = _preflight_archive(handle, selected_root)
            target_dir.mkdir(parents=True, mode=0o700)
            try:
                for member, relative in members:
                    destination = target_dir.joinpath(*relative.parts)
                    if member.isdir():
                        destination.mkdir(parents=True, exist_ok=True, mode=0o755)
                        continue
                    destination.parent.mkdir(parents=True, exist_ok=True, mode=0o755)
                    source_handle = handle.extractfile(member)
                    if source_handle is None:
                        raise ArchiveSafetyError(
                            f"Could not read archive member: {relative.as_posix()}"
                        )
                    with source_handle, destination.open("xb") as output_handle:
                        shutil.copyfileobj(source_handle, output_handle, 1024 * 1024)
                    executable = preserve_executable and bool(member.mode & 0o111)
                    destination.chmod(0o755 if executable else 0o644)
            except Exception:
                shutil.rmtree(target_dir, ignore_errors=True)
                raise
    except (tarfile.TarError, OSError) as exc:
        if target_dir.exists() and not target_dir.is_symlink():
            shutil.rmtree(target_dir, ignore_errors=True)
        if isinstance(exc, ArchiveSafetyError):
            raise
        raise ArchiveSafetyError(f"Could not safely extract {archive_path}: {exc}") from exc


def safe_extract_archive(archive_path: Path, target_dir: Path) -> None:
    """Safely extract an authenticated core archive rooted at ``virosync/``."""
    _safe_extract_archive(
        archive_path,
        target_dir,
        archive_root=ARCHIVE_ROOT,
    )


def safe_extract_optional_archive(archive_path: Path, target_dir: Path) -> None:
    """Safely extract an optional archive with one named or ``.`` root."""
    _safe_extract_archive(
        archive_path,
        target_dir,
        archive_root=None,
        preserve_executable=True,
    )


def _validate_extracted_inventory(root: Path, required_files: list[str]) -> None:
    expected_files = {*required_files, "RESOURCE_MANIFEST.json"}
    expected_directories = {
        parent.as_posix()
        for relative in expected_files
        for parent in PurePosixPath(relative).parents
        if parent.as_posix() != "."
    }
    actual_files: set[str] = set()
    actual_directories: set[str] = set()
    for path in root.rglob("*"):
        relative = path.relative_to(root).as_posix()
        metadata = path.lstat()
        if stat.S_ISLNK(metadata.st_mode):
            raise ResourceInstallError(f"Extracted resource tree contains a link: {relative}")
        if stat.S_ISREG(metadata.st_mode):
            actual_files.add(relative)
        elif stat.S_ISDIR(metadata.st_mode):
            actual_directories.add(relative)
        else:
            raise ResourceInstallError(
                f"Extracted resource tree contains a special file: {relative}"
            )
    if actual_files != expected_files or not actual_directories.issubset(expected_directories):
        raise ResourceInstallError(
            "Extracted resource payload set is incomplete or unexpected; "
            f"missing={sorted(expected_files - actual_files)}, "
            f"unexpected={sorted(actual_files - expected_files)}, "
            f"unexpected_directories={sorted(actual_directories - expected_directories)}"
        )


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _fsync_tree(root: Path) -> None:
    """Make regular-file data, modes, and directory entries crash-durable."""
    directories = [root]
    for path in sorted(root.rglob("*")):
        metadata = path.lstat()
        if stat.S_ISLNK(metadata.st_mode):
            raise ResourceInstallError(f"Durable core tree contains a symlink: {path}")
        if stat.S_ISDIR(metadata.st_mode):
            directories.append(path)
            continue
        if not stat.S_ISREG(metadata.st_mode):
            raise ResourceInstallError(f"Durable core tree contains a special file: {path}")
        flags = (
            os.O_RDONLY
            | getattr(os, "O_CLOEXEC", 0)
            | getattr(os, "O_NOFOLLOW", 0)
        )
        descriptor = os.open(path, flags)
        try:
            opened = os.fstat(descriptor)
            if not stat.S_ISREG(opened.st_mode):
                raise ResourceInstallError(
                    f"Durable core payload changed type while opening: {path}"
                )
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    for directory in sorted(
        directories,
        key=lambda item: len(item.parts),
        reverse=True,
    ):
        _fsync_directory(directory)


def _write_journal(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_name(f"{path.name}.tmp-{uuid4().hex}")
    data = (json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n").encode()
    descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.rename(temporary, path)
        _fsync_directory(path.parent)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def _remove_journal(
    path: Path,
    expected_identity: tuple[int, int] | None = None,
) -> None:
    try:
        metadata = path.lstat()
    except FileNotFoundError:
        if expected_identity is not None:
            raise ResourceInstallError(f"Resource recovery journal disappeared: {path}")
        return
    if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
        raise ResourceInstallError(
            f"Resource recovery journal must be a single-link regular file: {path}"
        )
    if expected_identity is not None and (metadata.st_dev, metadata.st_ino) != expected_identity:
        raise ResourceInstallError(f"Resource recovery journal changed while reading: {path}")
    path.unlink()
    _fsync_directory(path.parent)


def _journal_path(target: Path) -> Path:
    return target.parent / f".{target.name}.resource-install-journal.json"


def _lock_path(target: Path) -> Path:
    # Stable sibling lock: every installer for this pointer contends on one inode.
    return target.parent / f".{target.name}.resource-install.lock"


@contextmanager
def sibling_install_lock(target: Path) -> Iterator[None]:
    """Serialize all setup/recovery work for one stable resource pointer."""
    target.parent.mkdir(parents=True, exist_ok=True)
    lock_path = _lock_path(target)
    flags = (
        os.O_RDWR
        | os.O_CREAT
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    try:
        descriptor = os.open(lock_path, flags, 0o600)
    except OSError as exc:
        raise ResourceInstallError(
            f"Could not open regular sibling install lock: {lock_path}"
        ) from exc
    try:
        opened = os.fstat(descriptor)
        if not stat.S_ISREG(opened.st_mode) or opened.st_nlink != 1:
            raise ResourceInstallError(
                f"Sibling install lock must be a single-link regular file: {lock_path}"
            )
        os.fchmod(descriptor, 0o600)
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        current = lock_path.lstat()
        if (
            not stat.S_ISREG(current.st_mode)
            or current.st_nlink != 1
            or (current.st_dev, current.st_ino) != (opened.st_dev, opened.st_ino)
        ):
            raise ResourceInstallError(
                f"Sibling install lock changed while acquiring it: {lock_path}"
            )
        yield
    finally:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        except OSError:
            pass
        os.close(descriptor)


def _safe_sibling(parent: Path, name: str | None, label: str) -> Path | None:
    if name is None:
        return None
    if not name or Path(name).name != name or name in {".", ".."}:
        raise ResourceInstallError(f"Unsafe {label} in recovery journal: {name!r}")
    return parent / name


def _generated_name(target: Path, name: str, kind: str) -> bool:
    prefix = re.escape(f".{target.name}.{kind}-")
    return re.fullmatch(prefix + r"[0-9a-f]{32}", name) is not None


def _atomic_relative_pointer(target: Path, link: str, kind: str) -> None:
    temporary = target.parent / f".{target.name}.{kind}-{uuid4().hex}"
    os.symlink(link, temporary)
    try:
        os.replace(temporary, target)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    _fsync_directory(target.parent)


def _restore_from_journal(target: Path, journal: dict[str, Any]) -> None:
    parent = target.parent
    if not isinstance(journal, dict):
        raise ResourceInstallError(f"Invalid resource recovery journal: {_journal_path(target)}")
    if (
        type(journal.get("schema_version")) is not int
        or journal.get("schema_version") != 1
        or journal.get("target") != target.name
        or journal.get("state") not in {"prepared", "prior_moved", "activated"}
    ):
        raise ResourceInstallError(f"Invalid resource recovery journal: {_journal_path(target)}")

    temporary = _safe_sibling(parent, journal.get("temporary_pointer"), "temporary pointer")
    if temporary is None or not _generated_name(target, temporary.name, "activate"):
        raise ResourceInstallError(
            f"Unsafe temporary pointer in recovery journal: {journal.get('temporary_pointer')!r}"
        )
    candidate = _safe_sibling(parent, journal.get("candidate"), "candidate")
    candidate_pattern = re.compile(
        re.escape(f"{target.name}-") + r"[A-Za-z0-9._-]+-[0-9a-f]{16}"
    )
    if candidate is None or candidate_pattern.fullmatch(candidate.name) is None:
        raise ResourceInstallError(
            f"Unsafe candidate in recovery journal: {journal.get('candidate')!r}"
        )
    protected_names = {
        target.name,
        candidate.name,
        _journal_path(target).name,
        _lock_path(target).name,
    }
    if temporary.name in protected_names:
        raise ResourceInstallError(
            f"Unsafe temporary pointer alias in recovery journal: {temporary.name!r}"
        )
    if temporary.is_symlink():
        if os.readlink(temporary) != candidate.name:
            raise ResourceInstallError(
                f"Recovery temporary pointer selects an unexpected target: {temporary}"
            )
        temporary.unlink()
    elif temporary.exists():
        raise ResourceInstallError(
            f"Recovery temporary pointer is not a symlink: {temporary}"
        )
    prior_kind = journal.get("prior_kind")
    if prior_kind == "directory":
        retained = _safe_sibling(parent, journal.get("prior"), "retained directory")
        if (
            retained is None
            or not _generated_name(target, retained.name, "legacy")
            or retained.name in protected_names
            or retained == temporary
        ):
            raise ResourceInstallError("Recovery journal lacks retained directory")
        retained_ready = retained.is_dir() and not retained.is_symlink()
        if target.is_symlink():
            if not retained_ready:
                raise ResourceInstallError(
                    f"Retained resource directory is unavailable: {retained}"
                )
            active_link = os.readlink(target)
            if active_link == retained.name:
                return
            if active_link != candidate.name:
                raise ResourceInstallError(
                    f"Active resource pointer is unrelated to recovery journal: {target}"
                )
        if target.exists():
            if target.is_symlink():
                pass
            elif retained.exists():
                raise ResourceInstallError(
                    "Both active and retained real directories exist during recovery"
                )
            elif target.is_dir():
                return
            else:
                raise ResourceInstallError(
                    f"Refusing to replace non-directory resource path during recovery: {target}"
                )
        if not retained_ready:
            raise ResourceInstallError(f"Retained resource directory is unavailable: {retained}")
        _atomic_relative_pointer(target, retained.name, "recover")
        return

    if prior_kind == "symlink":
        prior_path = _safe_sibling(parent, journal.get("prior"), "prior pointer")
        if prior_path is None or prior_path.name in {
            target.name,
            temporary.name,
            _journal_path(target).name,
            _lock_path(target).name,
        }:
            raise ResourceInstallError("Recovery journal contains an unsafe prior pointer")
        prior_link = prior_path.name
        resolved_prior = prior_path.resolve(strict=False)
        if resolved_prior.parent != parent.resolve():
            raise ResourceInstallError("Recovery journal prior pointer escapes resource parent")
        if not resolved_prior.is_dir() or resolved_prior.is_symlink():
            raise ResourceInstallError(
                f"Recovery journal prior resource directory is unavailable: {resolved_prior}"
            )
        if target.is_symlink():
            active_link = os.readlink(target)
            if active_link == prior_link:
                return
            if active_link != candidate.name:
                raise ResourceInstallError(
                    f"Active resource pointer is unrelated to recovery journal: {target}"
                )
        elif target.exists():
            raise ResourceInstallError(
                f"Refusing to replace non-pointer resource path during recovery: {target}"
            )
        _atomic_relative_pointer(target, prior_link, "recover")
        return

    if prior_kind == "missing":
        if journal.get("prior") is not None:
            raise ResourceInstallError("Recovery journal has a prior path for missing state")
        if target.is_symlink():
            if os.readlink(target) != candidate.name:
                raise ResourceInstallError(
                    f"Active resource pointer is unrelated to recovery journal: {target}"
                )
            if candidate.is_symlink() or not candidate.is_dir():
                raise ResourceInstallError(
                    f"Activated recovery candidate is unavailable: {candidate}"
                )
            _fsync_directory(parent)
            return
        elif target.exists():
            raise ResourceInstallError(
                f"Refusing to remove non-pointer resource path during recovery: {target}"
            )
        return
    raise ResourceInstallError(f"Unknown prior resource state in journal: {prior_kind!r}")


def recover_pending_install(target: Path) -> bool:
    """Roll back an interrupted activation while retaining every resource tree."""
    target = Path(target)
    journal_path = _journal_path(target)
    flags = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    try:
        descriptor = os.open(journal_path, flags)
    except FileNotFoundError:
        return False
    except OSError as exc:
        raise ResourceInstallError(
            f"Could not open regular resource recovery journal: {journal_path}"
        ) from exc
    try:
        opened = os.fstat(descriptor)
        current = journal_path.lstat()
        if (
            not stat.S_ISREG(opened.st_mode)
            or opened.st_nlink != 1
            or opened.st_size > 64 * 1024
            or not stat.S_ISREG(current.st_mode)
            or current.st_nlink != 1
            or (current.st_dev, current.st_ino) != (opened.st_dev, opened.st_ino)
        ):
            raise ResourceInstallError(
                f"Resource recovery journal must be a small single-link regular file: {journal_path}"
            )
        with os.fdopen(descriptor, "r", encoding="utf-8") as handle:
            descriptor = -1
            journal = json.load(handle)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ResourceInstallError(f"Unreadable resource recovery journal: {journal_path}") from exc
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    _restore_from_journal(target, journal)
    _remove_journal(journal_path, (opened.st_dev, opened.st_ino))
    return True


def _unique_retained_path(target: Path) -> Path:
    while True:
        candidate = target.parent / f".{target.name}.legacy-{uuid4().hex}"
        if not candidate.exists() and not candidate.is_symlink():
            return candidate


def _prior_state(target: Path) -> tuple[str, str | None]:
    if target.is_symlink():
        link = os.readlink(target)
        if Path(link).is_absolute():
            raise ResourceInstallError("Active resource pointer must be relative")
        resolved = (target.parent / link).resolve(strict=False)
        if resolved.parent != target.parent.resolve():
            raise ResourceInstallError("Active resource pointer escapes its parent")
        return "symlink", link
    if target.exists():
        if not target.is_dir():
            raise ResourceInstallError(f"Active resource path is not a directory: {target}")
        return "directory", _unique_retained_path(target).name
    return "missing", None


def _activate_candidate(
    target: Path,
    candidate: Path,
    *,
    fault_injector: FaultInjector | None,
) -> None:
    parent = target.parent
    if candidate.parent != parent:
        raise ResourceInstallError("Candidate and active resource path are not siblings")
    prior_kind, prior = _prior_state(target)
    temporary = parent / f".{target.name}.activate-{uuid4().hex}"
    journal = {
        "schema_version": 1,
        "state": "prepared",
        "target": target.name,
        "candidate": candidate.name,
        "temporary_pointer": temporary.name,
        "prior_kind": prior_kind,
        "prior": prior,
    }
    journal_path = _journal_path(target)
    try:
        _write_journal(journal_path, journal)
        _invoke_fault(fault_injector, "after_journal_write")
        os.symlink(candidate.name, temporary)
        _invoke_fault(fault_injector, "after_pointer_prepare")
        if prior_kind == "directory":
            retained = parent / str(prior)
            os.replace(target, retained)
            _fsync_directory(parent)
        journal["state"] = "prior_moved"
        _write_journal(journal_path, journal)
        _invoke_fault(fault_injector, "after_prior_move")
        os.replace(temporary, target)
        _fsync_directory(parent)
        journal["state"] = "activated"
        _write_journal(journal_path, journal)
        _invoke_fault(fault_injector, "after_pointer_activate")
        if not target.is_symlink() or target.resolve(strict=True) != candidate.resolve(strict=True):
            raise ResourceInstallError("Atomic resource activation did not select the candidate")
        _invoke_fault(fault_injector, "after_activation_validate")
        _remove_journal(journal_path)
        _invoke_fault(fault_injector, "after_journal_clear")
    except Exception:
        try:
            try:
                journal_path.lstat()
            except FileNotFoundError:
                journal_present = False
            else:
                journal_present = True
            if journal_present:
                recover_pending_install(target)
            elif temporary.is_symlink():
                temporary.unlink()
        except Exception:
            logger.exception("Resource activation rollback requires recovery on the next setup")
        raise


def _chmod_if_needed(path: Path, mode: int) -> None:
    if stat.S_IMODE(path.lstat().st_mode) != mode:
        path.chmod(mode)


def _make_tree_immutable(root: Path, *, finalize_root: bool = True) -> None:
    paths = sorted(root.rglob("*"), key=lambda item: len(item.parts), reverse=True)
    for path in paths:
        metadata = path.lstat()
        if stat.S_ISREG(metadata.st_mode) and metadata.st_nlink != 1:
            raise ResourceInstallError(
                f"Immutable core tree contains a multiply linked file: {path}"
            )
    _chmod_if_needed(root, 0o755)
    for path in paths:
        if path.is_symlink():
            raise ResourceInstallError(f"Immutable core tree contains a symlink: {path}")
        if path.is_file():
            executable = bool(path.stat().st_mode & stat.S_IXUSR)
            _chmod_if_needed(path, 0o555 if executable else 0o444)
        elif path.is_dir():
            _chmod_if_needed(path, 0o555)
        else:
            raise ResourceInstallError(f"Immutable core tree contains a special file: {path}")
    _chmod_if_needed(root, 0o555 if finalize_root else 0o755)


def _make_tree_removable(root: Path) -> None:
    if not root.exists() or root.is_symlink():
        return
    for path in root.rglob("*"):
        if path.is_dir() and not path.is_symlink():
            path.chmod(0o700)
        elif path.is_file():
            path.chmod(0o600)
    root.chmod(0o700)


def _write_install_metadata(
    root: Path,
    source: ResourceSource,
    required_files: list[str],
) -> None:
    verified_files = {}
    for relative in required_files:
        metadata = (root / relative).lstat()
        verified_files[relative] = {
            "device": metadata.st_dev,
            "inode": metadata.st_ino,
            "size": metadata.st_size,
            "mtime_ns": metadata.st_mtime_ns,
            "ctime_ns": metadata.st_ctime_ns,
        }
    payload = {
        "component": "virosync_core",
        "version": source.version,
        "installed_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": source.source,
        "archive_sha256": source.archive_sha256,
        "manifest_sha256": source.manifest_sha256,
        "required_files": required_files,
        "verified_files": verified_files,
    }
    destination = root / "DB_METADATA.json"
    if destination.exists():
        destination.chmod(0o600)
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    destination.chmod(0o444)


def verified_install_receipt(
    root: Path,
    manifest: Any,
    *,
    expected_archive_sha256: str | None = None,
) -> bool:
    """Return whether an install receipt authenticates the current payload stats."""
    descriptor = -1
    try:
        resource_root = Path(root).resolve(strict=True)
        if not resource_root.is_dir() or resource_root.is_symlink():
            return False
        metadata_path = resource_root / "DB_METADATA.json"
        descriptor = os.open(
            metadata_path,
            os.O_RDONLY
            | getattr(os, "O_CLOEXEC", 0)
            | getattr(os, "O_NOFOLLOW", 0),
        )
        opened = os.fstat(descriptor)
        current = metadata_path.lstat()
        if (
            not stat.S_ISREG(opened.st_mode)
            or opened.st_nlink != 1
            or opened.st_size <= 0
            or opened.st_size > INSTALL_METADATA_MAX_BYTES
            or not stat.S_ISREG(current.st_mode)
            or current.st_nlink != 1
            or (current.st_dev, current.st_ino) != (opened.st_dev, opened.st_ino)
        ):
            return False
        with os.fdopen(descriptor, "r", encoding="utf-8") as handle:
            descriptor = -1
            payload = json.load(handle)
    except (OSError, TypeError, ValueError, UnicodeDecodeError, json.JSONDecodeError):
        return False
    finally:
        if descriptor >= 0:
            os.close(descriptor)

    expected_paths = tuple(item.path for item in manifest.files)
    if (
        not isinstance(payload, dict)
        or payload.get("component") != "virosync_core"
        or payload.get("manifest_sha256") != manifest.manifest_sha256
        or payload.get("version") != manifest.version
        or payload.get("required_files") != list(expected_paths)
        or (
            expected_archive_sha256 is not None
            and payload.get("archive_sha256") != expected_archive_sha256
        )
    ):
        return False
    verified_files = payload.get("verified_files")
    if not isinstance(verified_files, dict) or set(verified_files) != set(expected_paths):
        return False

    allowed_files = {*expected_paths, RESOURCE_MANIFEST_NAME, "DB_METADATA.json"}
    allowed_directories = {
        parent.as_posix()
        for relative in allowed_files
        for parent in PurePosixPath(relative).parents
        if parent.as_posix() != "."
    }
    try:
        root_metadata = resource_root.lstat()
        if (
            not stat.S_ISDIR(root_metadata.st_mode)
            or root_metadata.st_mode & 0o222
        ):
            return False
        actual_files: set[str] = set()
        actual_directories: set[str] = set()
        for candidate in resource_root.rglob("*"):
            relative = candidate.relative_to(resource_root).as_posix()
            metadata = candidate.lstat()
            if metadata.st_mode & 0o222:
                return False
            if stat.S_ISREG(metadata.st_mode):
                actual_files.add(relative)
            elif stat.S_ISDIR(metadata.st_mode):
                actual_directories.add(relative)
            else:
                return False
    except OSError:
        return False
    if actual_files != allowed_files or actual_directories != allowed_directories:
        return False

    expected_fields = {"device", "inode", "size", "mtime_ns", "ctime_ns"}
    for item in manifest.files:
        try:
            candidate = resource_root / item.path
            metadata = candidate.lstat()
        except OSError:
            return False
        if (
            not stat.S_ISREG(metadata.st_mode)
            or candidate.is_symlink()
            or metadata.st_nlink != 1
            or metadata.st_size != item.size
        ):
            return False
        receipt = verified_files.get(item.path)
        if not isinstance(receipt, dict) or set(receipt) != expected_fields:
            return False
        observed = {
            "device": metadata.st_dev,
            "inode": metadata.st_ino,
            "size": metadata.st_size,
            "mtime_ns": metadata.st_mtime_ns,
            "ctime_ns": metadata.st_ctime_ns,
        }
        if receipt != observed:
            return False
    return True


def _install_receipt_is_current(
    root: Path,
    source: ResourceSource,
    *,
    manifest: Any | None = None,
) -> bool:
    if manifest is None:
        manifest = _load_core_install_manifest(root, source)
    return verified_install_receipt(
        root,
        manifest,
        expected_archive_sha256=source.archive_sha256,
    )


def _load_core_install_manifest(root: Path, source: ResourceSource) -> Any:
    manifest = load_resource_manifest(
        root,
        expected_version=source.version,
        expected_manifest_sha256=source.manifest_sha256,
    )
    if manifest.bundle_kind == "source":
        raise ResourceInstallError(
            "A schema-v2 source/repair bundle cannot be installed as "
            "ViroSync core runtime resources"
        )
    return manifest


def _candidate_path(target: Path, source: ResourceSource) -> Path:
    safe_version = re.sub(r"[^A-Za-z0-9._-]+", "-", source.version).strip("-")
    if not safe_version:
        raise ResourceInstallError(f"Invalid resource version: {source.version!r}")
    return target.parent / f"{target.name}-{safe_version}-{source.manifest_sha256[:16]}"


def _active_candidate(target: Path, candidate: Path) -> Path | None:
    if not target.is_symlink() or candidate.is_symlink() or not candidate.is_dir():
        return None
    link = os.readlink(target)
    if Path(link).is_absolute() or Path(link).name != link or link != candidate.name:
        return None
    try:
        if target.resolve(strict=True) != candidate.resolve(strict=True):
            return None
    except OSError:
        return None
    return candidate


def active_installed_candidate(
    target: Path,
    source: ResourceSource,
) -> Path | None:
    """Return the exact relative sibling candidate selected for *source*."""
    target = Path(target)
    return _active_candidate(target, _candidate_path(target, source))


def install_core_resources(
    target: Path,
    source: ResourceSource,
    *,
    copy_archive: ArchiveCopier,
    verify_tree: TreeVerifier,
    required_files: list[str],
    full: bool,
    reuse_existing: bool = False,
    reject_invalid_existing: bool = False,
    semantic_runner: Callable[..., subprocess.CompletedProcess] | None = None,
    fault_injector: FaultInjector | None = None,
    progress_callback: Callable[[float, str], None] | None = None,
) -> Path:
    """Authenticate, stage, validate, finalize, and atomically activate a bundle."""
    target = Path(target)

    def report(progress: float, stage: str) -> None:
        if progress_callback is not None:
            progress_callback(progress, stage)

    source = ResourceSource(
        version=source.version,
        source=source.source,
        filename=Path(source.filename).name,
        archive_sha256=normalize_sha256(source.archive_sha256, "archive_sha256"),
        manifest_sha256=normalize_sha256(source.manifest_sha256, "manifest_sha256"),
    )
    if not source.filename or source.filename in {".", ".."}:
        raise ResourceInstallError("Resource archive filename is invalid")

    target.parent.mkdir(parents=True, exist_ok=True)
    with sibling_install_lock(target):
        recover_pending_install(target)
        candidate = _candidate_path(target, source)
        verification_kwargs: dict[str, Any] = {
            "expected_version": source.version,
            "expected_manifest_sha256": source.manifest_sha256,
            "verify_hashes": True,
            "full": full,
        }
        if semantic_runner is not None:
            verification_kwargs["command_runner"] = semantic_runner

        if reuse_existing:
            active = _active_candidate(target, candidate)
            if active is not None:
                manifest = _load_core_install_manifest(active, source)
                required_files = [item.path for item in manifest.files]
                receipt_is_current = _install_receipt_is_current(
                    active,
                    source,
                    manifest=manifest,
                )
                if receipt_is_current:
                    report(100, "core resources already installed")
                    return target
                report(60, "validating existing resources")
                verify_tree(active, **verification_kwargs)
                _make_tree_immutable(active, finalize_root=False)
                _write_install_metadata(active, source, required_files)
                _chmod_if_needed(active, 0o555)
                _fsync_tree(active)
                _fsync_directory(target.parent)
                report(100, "core resources ready")
                return target
            if (
                reject_invalid_existing
                and target.exists()
                and not target.is_symlink()
            ):
                verify_tree(
                    target,
                    expected_version=source.version,
                    expected_manifest_sha256=None,
                    verify_hashes=False,
                    full=False,
                )

        if not full:
            raise ResourceInstallError(
                "Full semantic validation is required before resource activation; "
                "use the fast verify command only for an existing active release"
            )

        stage_root = Path(
            tempfile.mkdtemp(prefix=f".{target.name}.stage-", dir=target.parent)
        )
        payload_root = stage_root / "payload"
        archive_path = stage_root / source.filename
        promoted = False
        try:
            report(5, "downloading core resources")
            copy_archive(source.source, archive_path)
            _invoke_fault(fault_injector, "after_download")
            if not archive_path.is_file() or archive_path.stat().st_size == 0:
                raise ResourceInstallError("Downloaded core-resource archive is empty")
            report(35, "verifying archive checksum")
            actual_archive_sha = sha256_file(archive_path)
            if actual_archive_sha != source.archive_sha256:
                raise ResourceInstallError(
                    "Core-resource archive SHA-256 mismatch: "
                    f"expected {source.archive_sha256}, found {actual_archive_sha}"
                )
            _invoke_fault(fault_injector, "after_archive_verify")
            report(45, "extracting archive")
            safe_extract_archive(archive_path, payload_root)
            manifest = _load_core_install_manifest(payload_root, source)
            required_files = [item.path for item in manifest.files]
            _validate_extracted_inventory(payload_root, required_files)
            _invoke_fault(fault_injector, "after_extract")
            report(70, "validating extracted resources")
            verify_tree(payload_root, **verification_kwargs)
            _invoke_fault(fault_injector, "after_stage_validate")

            if candidate.exists() or candidate.is_symlink():
                if candidate.is_symlink() or not candidate.is_dir():
                    raise ResourceInstallError(
                        f"Immutable candidate path has an invalid type: {candidate}"
                )
                verification_kwargs["verify_hashes"] = True
                verify_tree(candidate, **verification_kwargs)
                receipt_is_current = _install_receipt_is_current(
                    candidate,
                    source,
                    manifest=manifest,
                )
                _make_tree_immutable(candidate, finalize_root=False)
                if not receipt_is_current:
                    _write_install_metadata(candidate, source, required_files)
                _chmod_if_needed(candidate, 0o555)
                _fsync_tree(candidate)
                shutil.rmtree(payload_root)
            else:
                os.rename(payload_root, candidate)
                promoted = True
                receipt_is_current = _install_receipt_is_current(
                    candidate,
                    source,
                    manifest=manifest,
                )
                _make_tree_immutable(candidate, finalize_root=False)
                if not receipt_is_current:
                    _write_install_metadata(candidate, source, required_files)
                _chmod_if_needed(candidate, 0o555)
                _fsync_tree(candidate)
            _fsync_directory(target.parent)
            _invoke_fault(fault_injector, "after_candidate_promote")
            report(90, "activating resources")
            _activate_candidate(
                target,
                candidate,
                fault_injector=fault_injector,
            )
            report(100, "core resources ready")
            return target
        except BaseException:
            if promoted and not target.is_symlink():
                # A fully validated candidate may remain for forensics/retry. It is
                # never visible at the stable path unless activation succeeds.
                logger.info("Retaining validated inactive resource candidate: %s", candidate)
            raise
        finally:
            if stage_root.exists():
                _make_tree_removable(stage_root)
                shutil.rmtree(stage_root, ignore_errors=True)
