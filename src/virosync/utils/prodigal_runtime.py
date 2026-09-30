"""Build and select ViroSync's source-pinned Prodigal-GV runtime on Linux x86-64."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import logging
import os
import platform
import shlex
import shutil
import subprocess
import urllib.request
import zipfile
from collections.abc import Sequence
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path

from virosync.utils.atomic_write import atomic_write_context

LOGGER = logging.getLogger(__name__)
RESOURCE_NAMES = ("recipe.json", "capacity.patch", "pixi.toml", "pixi.lock", "LICENSE")
SETUP_COMMAND = "python -m virosync.utils.prodigal_runtime setup"
# The locked ld_impl_linux-64 binary has a 255-byte relocation placeholder.
MAX_COMPILER_PREFIX_BYTES = 255


@dataclass(frozen=True)
class NativeRecipe:
    """Packaged source and build inputs for one immutable runtime installation."""

    archive_url: str
    archive_sha256: str
    source_directory: str
    compiler: str
    cflags: tuple[str, ...]
    resources: dict[str, bytes]
    sha256: str


def resolve_prodigal_executable() -> Path:
    """Return the verified corrected caller, without searching PATH or installing it."""
    recipe = _load_recipe()
    return _verify_installation(_runtime_directory(recipe), recipe)


def setup_prodigal_runtime() -> Path:
    """Build once in a stable owned prefix, retaining failed setup evidence."""
    recipe = _load_recipe()
    runtime = _runtime_directory(recipe)
    root = runtime.parent
    root.mkdir(parents=True, exist_ok=True, mode=0o700)
    if root.stat().st_uid != os.getuid():
        raise RuntimeError(f"Native setup requires a runtime directory owned by the current user: {root}")
    with (root / ".setup.lock").open("a", encoding="utf-8") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if (runtime / "receipt.json").exists():
            return _verify_installation(runtime, recipe)
        if runtime.exists():
            raise RuntimeError(
                f"Incomplete native build retained at {runtime}; move it aside before rerunning {SETUP_COMMAND}"
            )
        runtime.mkdir(mode=0o700)
        LOGGER.info("Building corrected Prodigal-GV in %s", runtime)
        _build_runtime(runtime, recipe)
        return _verify_installation(runtime, recipe)


def _load_recipe() -> NativeRecipe:
    if platform.system() != "Linux" or platform.machine() not in {"x86_64", "AMD64"}:
        raise RuntimeError("The corrected Prodigal-GV runtime requires Linux x86-64 and Pixi")
    directory = files("virosync").joinpath("data/prodigal_gv")
    resources = {name: directory.joinpath(name).read_bytes() for name in RESOURCE_NAMES}
    document = json.loads(resources["recipe.json"])
    if document["schema_version"] != 1 or _digest(resources["capacity.patch"]) != document["patch_sha256"]:
        raise RuntimeError("Packaged Prodigal-GV recipe or capacity patch is inconsistent")
    identity = {name: _digest(content) for name, content in resources.items()}
    return NativeRecipe(
        archive_url=document["archive_url"],
        archive_sha256=document["archive_sha256"],
        source_directory=document["source_directory"],
        compiler=document["compiler"],
        cflags=tuple(document["cflags"]),
        resources=resources,
        sha256=_digest(json.dumps(identity, sort_keys=True).encode()),
    )


def _runtime_root() -> Path:
    selected = os.environ.get("VIROSYNC_PRODIGAL_RUNTIME")
    root = Path(selected).expanduser() if selected else Path.home() / ".local/share/virosync/prodigal-gv"
    if not root.is_absolute():
        raise ValueError("VIROSYNC_PRODIGAL_RUNTIME must name an absolute directory")
    return root.resolve()


def _runtime_directory(recipe: NativeRecipe) -> Path:
    """Check the locked compiler's binary-prefix limit before setup creates files."""
    runtime = _runtime_root() / recipe.sha256
    environment = (runtime / "build/.pixi/envs/default").resolve()
    prefix_bytes = len(os.fsencode(environment))
    if prefix_bytes > MAX_COMPILER_PREFIX_BYTES:
        raise RuntimeError(
            f"Native compiler environment prefix is {prefix_bytes} bytes; "
            f"the locked compiler allows at most {MAX_COMPILER_PREFIX_BYTES}: {environment}. "
            "Set VIROSYNC_PRODIGAL_RUNTIME to a shorter absolute path in an owned permanent "
            "directory before setup, and keep that value for subsequent runs."
        )
    return runtime


def _verify_installation(runtime: Path, recipe: NativeRecipe) -> Path:
    receipt_path = runtime / "receipt.json"
    if not receipt_path.is_file():
        raise RuntimeError(f"Corrected Prodigal-GV is not installed at {runtime}; run {SETUP_COMMAND}")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("schema_version") != 1 or receipt.get("recipe_sha256") != recipe.sha256:
        raise RuntimeError(f"Native runtime receipt does not match the packaged recipe: {runtime}")
    if receipt.get("runtime_directory") != str(runtime):
        raise RuntimeError(f"Native runtime was moved; its compiler/library prefix must remain at {runtime}")
    for name in RESOURCE_NAMES:
        path = runtime / "build" / name
        if path.is_symlink() or not path.is_file() or path.read_bytes() != recipe.resources[name]:
            raise RuntimeError(f"Native runtime build inputs changed: {path}")
    executable = runtime / "bin/prodigal-gv"
    if executable.is_symlink() or not executable.is_file() or not os.access(executable, os.X_OK):
        raise RuntimeError(f"Corrected Prodigal-GV executable is missing or invalid: {executable}")
    if _file_identity(executable) != receipt["binary"]:
        raise RuntimeError(f"Corrected Prodigal-GV executable changed: {executable}")
    environment = runtime / "build/.pixi/envs/default"
    _require_stable_environment(environment)
    zlib = environment / "lib/libz.so.1"
    if not zlib.is_file() or {"path": str(zlib.resolve()), **_file_identity(zlib)} != receipt["zlib"]:
        raise RuntimeError(f"Native runtime zlib changed or is missing: {zlib}")
    return executable


def _require_stable_environment(environment: Path) -> None:
    if not environment.is_dir() or environment.resolve() != environment:
        raise RuntimeError(f"Native runtime requires its retained environment at {environment}")


def _build_runtime(runtime: Path, recipe: NativeRecipe) -> None:
    build = runtime / "build"
    build.mkdir()
    for name, content in recipe.resources.items():
        (build / name).write_bytes(content)
    archive = build / "source.zip"
    with urllib.request.urlopen(recipe.archive_url, timeout=60) as response, archive.open("xb") as output:
        shutil.copyfileobj(response, output)
    if _file_identity(archive)["sha256"] != recipe.archive_sha256:
        raise RuntimeError(f"Prodigal-GV source archive hash mismatch: {archive}")
    with zipfile.ZipFile(archive) as compressed:
        compressed.extractall(build)
    source = build / recipe.source_directory
    environment = build / ".pixi/envs/default"
    compiler = environment / "bin" / recipe.compiler
    manifest = build / "pixi.toml"
    commands = [
        ["pixi", "install", "--locked", "--manifest-path", str(manifest)],
        [
            "pixi",
            "run",
            "--locked",
            "--manifest-path",
            str(manifest),
            "patch",
            "--batch",
            "--forward",
            "--fuzz=0",
            "--directory",
            str(source),
            "-p1",
            "-i",
            str(build / "capacity.patch"),
        ],
        _make_command(source, environment, recipe),
    ]
    with (runtime / "build.log").open("w", encoding="utf-8") as log:
        for command in commands:
            log.write(shlex.join(command) + "\n")
            log.flush()
            subprocess.run(command, cwd=source, check=True, stdout=log, stderr=subprocess.STDOUT)
            _require_stable_environment(environment)
    destination = runtime / "bin"
    destination.mkdir()
    executable = destination / "prodigal-gv"
    shutil.copy2(source / "prodigal-gv", executable)
    version = _normal_runtime_version(executable)
    receipt = {
        "schema_version": 1,
        "recipe_sha256": recipe.sha256,
        "runtime_directory": str(runtime),
        "source_archive_sha256": recipe.archive_sha256,
        "compiler": {"path": str(compiler), **_file_identity(compiler)},
        "zlib": {
            "path": str((environment / "lib/libz.so.1").resolve()),
            **_file_identity(environment / "lib/libz.so.1"),
        },
        "command": commands[-1],
        "binary": _file_identity(executable),
        "normal_runtime_version": version,
    }
    with atomic_write_context(runtime / "receipt.json") as handle:
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _make_command(source: Path, environment: Path, recipe: NativeRecipe) -> list[str]:
    compiler = environment / "bin" / recipe.compiler
    include = f"-I{environment / 'include'}"
    libraries = [f"-L{environment / 'lib'}", f"-Wl,-rpath,{environment / 'lib'}"]
    return [
        "pixi",
        "run",
        "--locked",
        "--manifest-path",
        str(environment.parents[2] / "pixi.toml"),
        "make",
        "-C",
        str(source),
        "-j1",
        f"CC={shlex.quote(str(compiler))}",
        f"CFLAGS={shlex.join([*recipe.cflags, include])}",
        f"LDFLAGS={shlex.join(libraries)}",
    ]


def _normal_runtime_version(executable: Path) -> str:
    environment = dict(os.environ)
    environment.pop("LD_LIBRARY_PATH", None)
    completed = subprocess.run(
        [str(executable), "-v"], check=True, capture_output=True, text=True, env=environment, timeout=30
    )
    return (completed.stdout + completed.stderr).strip()


def _file_identity(path: Path) -> dict[str, int | str]:
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    return {"size": path.stat().st_size, "sha256": digest}


def _digest(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def main(argv: Sequence[str] | None = None) -> int:
    """Install the packaged native runtime explicitly; never run a genome analysis."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["setup"], help="Build or verify the corrected native runtime.")
    parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    try:
        executable = setup_prodigal_runtime()
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        LOGGER.error("Native runtime setup failed: %s", error)
        return 1
    LOGGER.info("Corrected Prodigal-GV ready: %s", executable)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
