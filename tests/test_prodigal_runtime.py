"""Exercise native installation boundaries without downloads, compilers or gene calling."""

from __future__ import annotations

import io
import json
import os
import subprocess
import sys
import tempfile
import zipfile
from collections.abc import Iterator
from pathlib import Path
from types import SimpleNamespace

import pytest

from virosync.utils import prodigal_runtime as runtime


@pytest.fixture
def packaged_recipe(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[bytes]:
    """Provide a small packaged recipe and source archive through real file readers."""
    archive = io.BytesIO()
    with zipfile.ZipFile(archive, "w") as source:
        source.writestr("upstream/main.c", "fixture source\n")
    resources = tmp_path / "package/data/prodigal_gv"
    resources.mkdir(parents=True)
    (resources / "capacity.patch").write_bytes(b"fixture patch\n")
    (resources / "pixi.toml").write_bytes(b"fixture manifest\n")
    (resources / "pixi.lock").write_bytes(b"fixture lock\n")
    (resources / "LICENSE").write_bytes(b"fixture license\n")
    (resources / "recipe.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "archive_url": "https://example.invalid/source.zip",
                "archive_sha256": runtime._digest(archive.getvalue()),
                "source_directory": "upstream",
                "patch_sha256": runtime._digest(b"fixture patch\n"),
                "compiler": "fixture-cc",
                "cflags": ["-O3"],
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(runtime, "files", lambda _package: tmp_path / "package")
    monkeypatch.setattr(runtime.urllib.request, "urlopen", lambda *_args, **_kwargs: io.BytesIO(archive.getvalue()))
    # The real prefix limit also applies when pytest's base directory is deeply nested.
    with tempfile.TemporaryDirectory(prefix="vs-native-test-", dir="/tmp") as directory:
        monkeypatch.setenv("VIROSYNC_PRODIGAL_RUNTIME", str(Path(directory) / "owned runtime"))
        yield archive.getvalue()


@pytest.fixture
def native_commands(monkeypatch: pytest.MonkeyPatch) -> list[list[str]]:
    """Replace only expensive external commands, retaining installation and receipt logic."""
    commands: list[list[str]] = []

    def invoke(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        commands.append(command)
        if command[:2] == ["pixi", "install"]:
            build = Path(command[-1]).parent
            environment = build / ".pixi/envs/default"
            (environment / "bin").mkdir(parents=True)
            (environment / "bin/fixture-cc").write_bytes(b"compiler fixture")
            (environment / "lib").mkdir()
            (environment / "lib/libz.so.1").write_bytes(b"zlib fixture")
        elif "make" in command:
            source = Path(command[command.index("-C") + 1])
            executable = source / "prodigal-gv"
            executable.write_bytes(b"corrected executable fixture")
            executable.chmod(0o755)
        elif command[-1] == "-v":
            assert "LD_LIBRARY_PATH" not in kwargs["env"]
            return subprocess.CompletedProcess(command, 0, "", "Prodigal V2.11.0-gv\n")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(runtime.subprocess, "run", invoke)
    return commands


def test_setup_keeps_prefix_and_reuses_verified_binary(
    packaged_recipe: bytes, native_commands: list[list[str]], monkeypatch: pytest.MonkeyPatch
) -> None:
    """One setup feeds the same resolver without build activation or a second build."""
    monkeypatch.setenv("LD_LIBRARY_PATH", "/transient/build/library/path")
    executable = runtime.setup_prodigal_runtime()
    install = executable.parents[1]
    prefix = install / "build/.pixi/envs/default"
    assert prefix.is_dir()
    assert len(native_commands) == 4
    assert native_commands[1][5:9] == ["patch", "--batch", "--forward", "--fuzz=0"]
    compile_command = native_commands[2]
    assert f"CC='{prefix / 'bin/fixture-cc'}'" in compile_command
    assert any(str(prefix / "lib") in arg and "-rpath" in arg for arg in compile_command)
    receipt = json.loads((install / "receipt.json").read_text(encoding="utf-8"))
    assert receipt["runtime_directory"] == str(install)
    assert receipt["binary"] == runtime._file_identity(executable)
    assert receipt["normal_runtime_version"] == "Prodigal V2.11.0-gv"
    assert runtime.resolve_prodigal_executable() == executable
    assert runtime.setup_prodigal_runtime() == executable
    assert len(native_commands) == 4


@pytest.mark.parametrize("change", ["binary", "build_lock", "environment", "detached", "library", "location"])
def test_changed_installation_cannot_be_selected(
    packaged_recipe: bytes, native_commands: list[list[str]], change: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Reject damaged or relocated installations rather than using an old PATH caller."""
    executable = runtime.setup_prodigal_runtime()
    install = executable.parents[1]
    if change == "binary":
        executable.write_bytes(b"different binary")
    elif change == "build_lock":
        (install / "build/pixi.lock").write_bytes(b"different lock")
    elif change == "environment":
        (install / "build/.pixi/envs/default").rename(install / "lost-environment")
    elif change == "detached":
        environment = install / "build/.pixi/envs/default"
        detached = install / "detached-environment"
        environment.rename(detached)
        environment.symlink_to(detached, target_is_directory=True)
    elif change == "library":
        (install / "build/.pixi/envs/default/lib/libz.so.1").write_bytes(b"changed zlib")
    else:
        moved = install.parent.with_name("moved-runtime")
        install.parent.rename(moved)
        monkeypatch.setenv("VIROSYNC_PRODIGAL_RUNTIME", str(moved))
    with pytest.raises(RuntimeError):
        runtime.resolve_prodigal_executable()
    assert len(native_commands) == 4


def test_detached_environment_stops_before_patch_or_compile(
    packaged_recipe: bytes, native_commands: list[list[str]], monkeypatch: pytest.MonkeyPatch
) -> None:
    """A Pixi detached prefix cannot enter a build with a different retention contract."""
    invoke = runtime.subprocess.run

    def install_detached(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        completed = invoke(command, **kwargs)
        environment = Path(command[-1]).parent / ".pixi/envs/default"
        detached = environment.parent / "detached"
        environment.rename(detached)
        environment.symlink_to(detached, target_is_directory=True)
        return completed

    monkeypatch.setattr(runtime.subprocess, "run", install_detached)
    with pytest.raises(RuntimeError, match="retained environment"):
        runtime.setup_prodigal_runtime()
    assert len(native_commands) == 1
    assert native_commands[0][:2] == ["pixi", "install"]


def test_missing_corrected_runtime_never_uses_path(packaged_recipe: bytes, monkeypatch: pytest.MonkeyPatch) -> None:
    """A vulnerable executable on PATH cannot satisfy the installation requirement."""
    monkeypatch.setattr(runtime.shutil, "which", lambda _name: pytest.fail("unexpected PATH lookup"))
    with pytest.raises(RuntimeError, match="not installed.*python -m virosync"):
        runtime.resolve_prodigal_executable()


def test_archive_mismatch_retains_incomplete_evidence(
    packaged_recipe: bytes, native_commands: list[list[str]], monkeypatch: pytest.MonkeyPatch
) -> None:
    """A bad download never reaches compilation, publication or automatic replacement."""
    monkeypatch.setattr(runtime.urllib.request, "urlopen", lambda *_args, **_kwargs: io.BytesIO(b"wrong archive"))
    with pytest.raises(RuntimeError, match="archive hash mismatch"):
        runtime.setup_prodigal_runtime()
    assert native_commands == []
    with pytest.raises(RuntimeError, match="Incomplete native build retained"):
        runtime.setup_prodigal_runtime()


def test_packaged_patch_must_match_recipe(packaged_recipe: bytes) -> None:
    """Reject modified package inputs before creating a native installation."""
    runtime.files("virosync").joinpath("data/prodigal_gv/capacity.patch").write_bytes(b"changed patch")
    with pytest.raises(RuntimeError, match="capacity patch is inconsistent"):
        runtime.setup_prodigal_runtime()


def test_relative_runtime_root_is_rejected(packaged_recipe: bytes, monkeypatch: pytest.MonkeyPatch) -> None:
    """Do not select a different installation when the working directory changes."""
    monkeypatch.setenv("VIROSYNC_PRODIGAL_RUNTIME", "relative/runtime")
    with pytest.raises(ValueError, match="absolute directory"):
        runtime.resolve_prodigal_executable()


def compiler_prefix_root(prefix_bytes: int, basename: str) -> Path:
    """Construct a runtime root whose complete compiler prefix has the requested byte length."""
    recipe = runtime._load_recipe()
    root = Path(os.environ["VIROSYNC_PRODIGAL_RUNTIME"]).parent / basename
    suffix = Path(recipe.sha256) / "build/.pixi/envs/default"
    padding = prefix_bytes - len(os.fsencode(root / suffix))
    assert padding >= 0
    root = root.with_name(root.name + "x" * padding)
    assert len(os.fsencode(root / suffix)) == prefix_bytes
    return root


def test_compiler_prefix_accepts_255_bytes(
    packaged_recipe: bytes, native_commands: list[list[str]], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Accept the locked compiler's exact prefix limit for setup and resolution."""
    root = compiler_prefix_root(255, "runtime")
    monkeypatch.setenv("VIROSYNC_PRODIGAL_RUNTIME", str(root))
    executable = runtime.setup_prodigal_runtime()
    assert runtime.resolve_prodigal_executable() == executable
    assert len(native_commands) == 4


@pytest.mark.parametrize("basename", ["runtime", "éé"])
@pytest.mark.parametrize("entrypoint", ["setup_prodigal_runtime", "resolve_prodigal_executable"])
def test_overlong_compiler_prefix_precedes_side_effects(
    packaged_recipe: bytes,
    native_commands: list[list[str]],
    monkeypatch: pytest.MonkeyPatch,
    basename: str,
    entrypoint: str,
) -> None:
    """Reject 256 filesystem bytes before creating files or invoking external commands."""
    root = compiler_prefix_root(256, basename)
    monkeypatch.setenv("VIROSYNC_PRODIGAL_RUNTIME", str(root))
    monkeypatch.setattr(
        runtime.urllib.request,
        "urlopen",
        lambda *_args, **_kwargs: pytest.fail("overlong prefix attempted a download"),
    )
    with pytest.raises(RuntimeError, match="256 bytes.*at most 255.*VIROSYNC_PRODIGAL_RUNTIME"):
        getattr(runtime, entrypoint)()
    assert not root.exists()
    assert native_commands == []


def test_version_probe_uses_resolver(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The auxiliary provenance probe must report the same corrected executable."""
    from virosync.utils import provenance

    executable = tmp_path / "corrected/prodigal-gv"
    observed: list[list[str]] = []
    monkeypatch.setitem(sys.modules, "torch", SimpleNamespace(__version__="fixture"))
    monkeypatch.setitem(sys.modules, "torch_geometric", SimpleNamespace(__version__="fixture"))
    monkeypatch.setitem(sys.modules, "pyhmmer", SimpleNamespace(__version__="fixture"))
    monkeypatch.setattr(provenance, "resolve_prodigal_executable", lambda: executable)

    def probe(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        observed.append(command)
        return subprocess.CompletedProcess(command, 0, "fixture", "")

    monkeypatch.setattr(provenance.subprocess, "run", probe)
    assert provenance.capture_tool_versions()["prodigal-gv"] == "fixture"
    assert [str(executable), "-v"] in observed
