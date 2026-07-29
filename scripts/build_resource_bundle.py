#!/usr/bin/env python3
"""Build an authenticated, deterministic ViroSync core-resource archive.

The source ``resources/virosync`` tree is read-only input. Generated metadata
and optional regenerated HMM/DIAMOND indices are staged elsewhere and supplied
to the archive as overrides, so a build never rewrites the installed resources.

Example
-------
    pixi run python scripts/build_resource_bundle.py \
        --version v1.0.6 --threads 16 \
        --skip-hmmpress --skip-marker-dmnd
"""

from __future__ import annotations

import argparse
from collections.abc import Callable, Mapping
from dataclasses import dataclass
import gzip
import io
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile

# Keep the documented ``pixi run python scripts/...`` invocation independent of
# whether the project has been installed in editable mode.
_SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

from virosync.utils.resource_manifest import (
    CORE_RESOURCE_FILES,
    RUNTIME_RESOURCE_FILES,
    RESOURCE_MANIFEST_NAME,
    SOURCE_RESOURCE_FILES,
    DiamondSequenceCounter,
    ResourceManifestError,
    ResourcePayload,
    build_resource_manifest,
    build_split_resource_manifests,
    diamond_sequence_count,
    sha256_file,
)

HMM_INDEX_SUFFIXES = (".h3f", ".h3i", ".h3m", ".h3p")
TAR_FORMAT = tarfile.GNU_FORMAT
_VERSION_RE = re.compile(r"^v[0-9]+\.[0-9]+\.[0-9]+$")


@dataclass(frozen=True)
class BundleBuildResult:
    """Identities of a completed bundle build."""

    output: Path
    version: str
    archive_sha256: str
    manifest_sha256: str


@dataclass(frozen=True)
class SplitBundleBuildResult:
    """Identities of completed runtime and source bundle builds."""

    runtime: BundleBuildResult
    source: BundleBuildResult


def _run(
    command: list[str],
    command_runner: Callable[..., subprocess.CompletedProcess[str]],
) -> None:
    print(f"  $ {' '.join(command)}", flush=True)
    command_runner(command, check=True)


def _require_tool(name: str, tool_finder: Callable[[str], str | None]) -> str:
    path = tool_finder(name)
    if path is None:
        raise ResourceManifestError(f"required tool {name!r} was not found on PATH")
    return path


def _count_prefixed_lines(path: Path, prefix: bytes) -> int:
    with path.open("rb") as handle:
        return sum(1 for line in handle if line.startswith(prefix))


def _readme_bytes(
    resources_dir: Path,
    version: str,
    *,
    split_runtime: bool = False,
) -> bytes:
    hmm_count = _count_prefixed_lines(resources_dir / "models/combined.hmm", b"NAME")
    marker_count = _count_prefixed_lines(resources_dir / "marker/marker.faa", b">")
    if split_runtime:
        content = (
            f"ViroSync core runtime resource bundle {version}\n"
            f"{'=' * 48}\n\n"
            f"HMM profiles (models/combined.hmm): {hmm_count}\n"
            f"TIER-1 marker proteins (source artifact): {marker_count}\n\n"
            "Runtime payload:\n"
            "  models/combined.hmm                         Phase 1 marker HMM library\n"
            "  models/model_annotations_with_interpro.tsv  one row per HMM profile\n"
            "  models/og_marker_name_map.tsv               VS<->OG marker name map\n"
            "  marker/marker.dmnd                           TIER-1 validation database\n"
            "  genomes/combined_proteome.dmnd               TIER-2 taxonomy database\n"
            "  taxonomy/labels.tsv                          genome-to-lineage labels\n\n"
            "Source/repair artifact:\n"
            "  models/combined.hmm.h3f/.h3i/.h3m/.h3p      HMMER indices\n"
            "  marker/marker.faa                            marker protein source\n\n"
            f"Payload integrity and semantic counts are recorded in {RESOURCE_MANIFEST_NAME}.\n"
        )
        return content.encode("utf-8")
    content = (
        f"ViroSync core resource bundle {version}\n"
        f"{'=' * 40}\n\n"
        f"HMM profiles (models/combined.hmm): {hmm_count}\n"
        f"TIER-1 marker proteins (marker/marker.faa): {marker_count}\n\n"
        "Layout:\n"
        "  models/combined.hmm(+.h3f/.h3i/.h3m/.h3p)  Phase 1 marker HMM library\n"
        "  models/model_annotations_with_interpro.tsv  one row per HMM profile\n"
        "  models/og_marker_name_map.tsv               VS<->OG marker name map\n"
        "  marker/marker.faa, marker/marker.dmnd        TIER-1 validation database\n"
        "  genomes/combined_proteome.dmnd               TIER-2 taxonomy database\n"
        "  taxonomy/labels.tsv                          genome-to-lineage labels\n\n"
        f"Payload integrity and semantic counts are recorded in {RESOURCE_MANIFEST_NAME}.\n"
    )
    return content.encode("utf-8")


def _derived_overrides(
    resources_dir: Path,
    staging_dir: Path,
    *,
    threads: int,
    skip_hmmpress: bool,
    skip_marker_dmnd: bool,
    command_runner: Callable[..., subprocess.CompletedProcess[str]],
    tool_finder: Callable[[str], str | None],
) -> dict[str, ResourcePayload]:
    overrides: dict[str, ResourcePayload] = {}
    combined_hmm = resources_dir / "models/combined.hmm"
    marker_faa = resources_dir / "marker/marker.faa"
    if not combined_hmm.is_file():
        raise ResourceManifestError(f"missing resource payload: {combined_hmm}")
    if not marker_faa.is_file():
        raise ResourceManifestError(f"missing resource payload: {marker_faa}")

    if skip_hmmpress:
        missing = [
            f"models/combined.hmm{suffix}"
            for suffix in HMM_INDEX_SUFFIXES
            if not Path(f"{combined_hmm}{suffix}").is_file()
        ]
        if missing:
            raise ResourceManifestError(
                "--skip-hmmpress requires all four existing indices; "
                f"missing={missing}"
            )
    else:
        hmmpress = _require_tool("hmmpress", tool_finder)
        staged_hmm = staging_dir / "models" / "combined.hmm"
        staged_hmm.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(combined_hmm, staged_hmm)
        _run([hmmpress, "-f", str(staged_hmm)], command_runner)
        for suffix in HMM_INDEX_SUFFIXES:
            relative = f"models/combined.hmm{suffix}"
            staged_index = Path(f"{staged_hmm}{suffix}")
            if not staged_index.is_file() or staged_index.stat().st_size <= 0:
                raise ResourceManifestError(
                    f"hmmpress did not create a non-empty {relative}"
                )
            overrides[relative] = staged_index

    marker_dmnd = resources_dir / "marker/marker.dmnd"
    if skip_marker_dmnd:
        if not marker_dmnd.is_file() or marker_dmnd.stat().st_size <= 0:
            raise ResourceManifestError(
                "--skip-marker-dmnd requires an existing non-empty marker/marker.dmnd"
            )
    else:
        diamond = _require_tool("diamond", tool_finder)
        staged_stem = staging_dir / "marker" / "marker"
        staged_stem.parent.mkdir(parents=True, exist_ok=True)
        _run(
            [
                diamond,
                "makedb",
                "--in",
                str(marker_faa),
                "--db",
                str(staged_stem),
                "--threads",
                str(threads),
            ],
            command_runner,
        )
        staged_dmnd = staged_stem.with_suffix(".dmnd")
        if not staged_dmnd.is_file() or staged_dmnd.stat().st_size <= 0:
            raise ResourceManifestError(
                "diamond makedb did not create a non-empty marker/marker.dmnd"
            )
        overrides["marker/marker.dmnd"] = staged_dmnd

    return overrides


def _payload(
    resources_dir: Path,
    relative: str,
    overrides: Mapping[str, ResourcePayload],
) -> ResourcePayload:
    return overrides.get(relative, resources_dir / relative)


def _payload_size(payload: ResourcePayload) -> int:
    return len(payload) if isinstance(payload, bytes) else payload.stat().st_size


def _add_regular_member(
    archive: tarfile.TarFile,
    arcname: str,
    payload: ResourcePayload,
) -> None:
    info = tarfile.TarInfo(arcname)
    info.size = _payload_size(payload)
    info.mode = 0o644
    info.mtime = 0
    info.uid = 0
    info.gid = 0
    info.uname = ""
    info.gname = ""
    if isinstance(payload, bytes):
        archive.addfile(info, io.BytesIO(payload))
    else:
        with payload.open("rb") as handle:
            archive.addfile(info, handle)


def create_deterministic_archive(
    resources_dir: Path,
    output: Path,
    *,
    overrides: Mapping[str, ResourcePayload],
    manifest_bytes: bytes,
    payload_files: tuple[str, ...] = CORE_RESOURCE_FILES,
) -> None:
    """Write a canonical gzip archive and atomically publish it."""

    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.",
        suffix=".tmp",
        dir=output.parent,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with temporary.open("wb") as raw:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                compresslevel=9,
                fileobj=raw,
                mtime=0,
            ) as compressed:
                with tarfile.open(
                    fileobj=compressed,
                    mode="w",
                    # combined_proteome.dmnd exceeds USTAR's 8 GiB size field.
                    # GNU tar uses deterministic base-256 encoding for it.
                    format=TAR_FORMAT,
                ) as archive:
                    for relative in payload_files:
                        _add_regular_member(
                            archive,
                            f"virosync/{relative}",
                            _payload(resources_dir, relative, overrides),
                        )
                    _add_regular_member(
                        archive,
                        f"virosync/{RESOURCE_MANIFEST_NAME}",
                        manifest_bytes,
                    )
            raw.flush()
            os.fsync(raw.fileno())
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def build_resource_bundle(
    resources_dir: Path,
    output: Path,
    version: str,
    *,
    threads: int = 8,
    skip_hmmpress: bool = False,
    skip_marker_dmnd: bool = False,
    command_runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
    tool_finder: Callable[[str], str | None] = shutil.which,
    diamond_sequence_counter: DiamondSequenceCounter | None = None,
) -> BundleBuildResult:
    """Build a complete bundle without modifying ``resources_dir``."""

    if _VERSION_RE.fullmatch(version) is None:
        raise ResourceManifestError("version must have the form vMAJOR.MINOR.PATCH")
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise ResourceManifestError("threads must be a positive integer")

    logical_resources_dir = Path(resources_dir).expanduser()
    if logical_resources_dir.name != "virosync":
        raise ResourceManifestError(
            "resources directory must be addressed through a path named 'virosync': "
            f"{logical_resources_dir}"
        )
    resources_dir = logical_resources_dir.resolve(strict=True)
    if not resources_dir.is_dir():
        raise ResourceManifestError(
            f"resources directory must be an existing directory: {resources_dir}"
        )
    output = Path(output).expanduser().resolve(strict=False)
    if output == resources_dir or resources_dir in output.parents:
        raise ResourceManifestError(
            "bundle output must be outside the resources directory"
        )
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(
        prefix=".virosync-bundle-stage-",
        dir=output.parent,
    ) as temporary_dir:
        overrides = _derived_overrides(
            resources_dir,
            Path(temporary_dir),
            threads=threads,
            skip_hmmpress=skip_hmmpress,
            skip_marker_dmnd=skip_marker_dmnd,
            command_runner=command_runner,
            tool_finder=tool_finder,
        )
        overrides["DB_VERSION"] = f"{version}\n".encode("utf-8")
        overrides["DATABASE_README.txt"] = _readme_bytes(resources_dir, version)

        if diamond_sequence_counter is None:

            def _count(payload: ResourcePayload) -> int:
                if not isinstance(payload, Path):
                    raise ResourceManifestError(
                        "DIAMOND database override must be a file path"
                    )
                return diamond_sequence_count(payload, command_runner=command_runner)

            diamond_sequence_counter = _count

        manifest, manifest_bytes = build_resource_manifest(
            resources_dir,
            version,
            overrides=overrides,
            diamond_sequence_counter=diamond_sequence_counter,
        )
        create_deterministic_archive(
            resources_dir,
            output,
            overrides=overrides,
            manifest_bytes=manifest_bytes,
        )

    return BundleBuildResult(
        output=output,
        version=version,
        archive_sha256=sha256_file(output),
        manifest_sha256=manifest.manifest_sha256,
    )


def build_split_resource_bundles(
    resources_dir: Path,
    runtime_output: Path,
    source_output: Path,
    version: str,
    *,
    threads: int = 8,
    skip_hmmpress: bool = False,
    skip_marker_dmnd: bool = False,
    command_runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
    tool_finder: Callable[[str], str | None] = shutil.which,
    diamond_sequence_counter: DiamondSequenceCounter | None = None,
) -> SplitBundleBuildResult:
    """Build bound runtime and source bundles without modifying ``resources_dir``."""

    if _VERSION_RE.fullmatch(version) is None:
        raise ResourceManifestError("version must have the form vMAJOR.MINOR.PATCH")
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise ResourceManifestError("threads must be a positive integer")

    logical_resources_dir = Path(resources_dir).expanduser()
    if logical_resources_dir.name != "virosync":
        raise ResourceManifestError(
            f"resources directory must be addressed through a path named 'virosync': {logical_resources_dir}"
        )
    resources_dir = logical_resources_dir.resolve(strict=True)
    if not resources_dir.is_dir():
        raise ResourceManifestError(
            f"resources directory must be an existing directory: {resources_dir}"
        )

    runtime_output = Path(runtime_output).expanduser().resolve(strict=False)
    source_output = Path(source_output).expanduser().resolve(strict=False)
    if runtime_output == source_output:
        raise ResourceManifestError("runtime and source bundle outputs must differ")
    for output in (runtime_output, source_output):
        if output == resources_dir or resources_dir in output.parents:
            raise ResourceManifestError(
                "bundle output must be outside the resources directory"
            )
        if output.exists() and output.is_dir():
            raise ResourceManifestError(
                f"bundle output must not be a directory: {output}"
            )
        output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(
        prefix=".virosync-bundle-stage-",
        dir=runtime_output.parent,
    ) as temporary_dir:
        overrides = _derived_overrides(
            resources_dir,
            Path(temporary_dir),
            threads=threads,
            skip_hmmpress=skip_hmmpress,
            skip_marker_dmnd=skip_marker_dmnd,
            command_runner=command_runner,
            tool_finder=tool_finder,
        )
        overrides["DB_VERSION"] = f"{version}\n".encode("utf-8")
        overrides["DATABASE_README.txt"] = _readme_bytes(
            resources_dir,
            version,
            split_runtime=True,
        )

        if diamond_sequence_counter is None:

            def _count(payload: ResourcePayload) -> int:
                if not isinstance(payload, Path):
                    raise ResourceManifestError(
                        "DIAMOND database override must be a file path"
                    )
                return diamond_sequence_count(payload, command_runner=command_runner)

            diamond_sequence_counter = _count

        (
            (runtime_manifest, runtime_manifest_bytes),
            (source_manifest, source_manifest_bytes),
        ) = build_split_resource_manifests(
            resources_dir,
            version,
            overrides=overrides,
            diamond_sequence_counter=diamond_sequence_counter,
        )

        with (
            tempfile.TemporaryDirectory(
                prefix=f".{runtime_output.name}.",
                dir=runtime_output.parent,
            ) as runtime_prepared_dir,
            tempfile.TemporaryDirectory(
                prefix=f".{source_output.name}.",
                dir=source_output.parent,
            ) as source_prepared_dir,
        ):
            runtime_temporary = Path(runtime_prepared_dir) / "bundle.tar.gz"
            source_temporary = Path(source_prepared_dir) / "bundle.tar.gz"
            create_deterministic_archive(
                resources_dir,
                runtime_temporary,
                overrides=overrides,
                manifest_bytes=runtime_manifest_bytes,
                payload_files=RUNTIME_RESOURCE_FILES,
            )
            create_deterministic_archive(
                resources_dir,
                source_temporary,
                overrides=overrides,
                manifest_bytes=source_manifest_bytes,
                payload_files=SOURCE_RESOURCE_FILES,
            )
            runtime_archive_sha256 = sha256_file(runtime_temporary)
            source_archive_sha256 = sha256_file(source_temporary)
            os.replace(runtime_temporary, runtime_output)
            os.replace(source_temporary, source_output)

    return SplitBundleBuildResult(
        runtime=BundleBuildResult(
            output=runtime_output,
            version=version,
            archive_sha256=runtime_archive_sha256,
            manifest_sha256=runtime_manifest.manifest_sha256,
        ),
        source=BundleBuildResult(
            output=source_output,
            version=version,
            archive_sha256=source_archive_sha256,
            manifest_sha256=source_manifest.manifest_sha256,
        ),
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--resources-dir",
        type=Path,
        default=Path("resources/virosync"),
        help="read-only resources/virosync tree to package",
    )
    parser.add_argument("--version", default="v1.0.6", help="bundle version")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="legacy or runtime tarball output",
    )
    parser.add_argument(
        "--split",
        action="store_true",
        help="build separate schema-v2 runtime and source bundles",
    )
    parser.add_argument(
        "--source-output",
        type=Path,
        default=None,
        help="source tarball output for --split",
    )
    parser.add_argument(
        "--threads", type=int, default=8, help="threads for diamond makedb"
    )
    parser.add_argument(
        "--skip-hmmpress", action="store_true", help="reuse existing .h3* indices"
    )
    parser.add_argument(
        "--skip-marker-dmnd",
        action="store_true",
        help="reuse the existing marker.dmnd",
    )
    parser.add_argument(
        "--no-proteome",
        action="store_true",
        help="retained compatibility flag; incomplete release bundles are rejected",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    if args.no_proteome:
        parser.error(
            "--no-proteome cannot produce a schema-v1 release bundle; "
            "genomes/combined_proteome.dmnd is required"
        )
    if args.source_output is not None and not args.split:
        parser.error("--source-output requires --split")
    try:
        if args.split:
            runtime_output = args.output or Path(
                f"resources_{args.version.replace('.', '_')}_runtime.tar.gz"
            )
            source_output = args.source_output or Path(
                f"resources_{args.version.replace('.', '_')}_source.tar.gz"
            )
            split_result = build_split_resource_bundles(
                args.resources_dir,
                runtime_output,
                source_output,
                args.version,
                threads=args.threads,
                skip_hmmpress=args.skip_hmmpress,
                skip_marker_dmnd=args.skip_marker_dmnd,
            )
        else:
            output = args.output or Path(
                f"resources_{args.version.replace('.', '_')}.tar.gz"
            )
            result = build_resource_bundle(
                args.resources_dir,
                output,
                args.version,
                threads=args.threads,
                skip_hmmpress=args.skip_hmmpress,
                skip_marker_dmnd=args.skip_marker_dmnd,
            )
    except (OSError, ResourceManifestError, subprocess.SubprocessError) as exc:
        parser.exit(2, f"ERROR: {exc}\n")

    if args.split:
        for label, bundle in (
            ("Runtime bundle", split_result.runtime),
            ("Source bundle", split_result.source),
        ):
            size_gb = bundle.output.stat().st_size / 1e9
            print(f"{label}: {bundle.output} ({size_gb:.2f} GB)")
            print(f"archive sha256: {bundle.archive_sha256}")
            print(f"manifest sha256: {bundle.manifest_sha256}")
    else:
        size_gb = result.output.stat().st_size / 1e9
        print(f"Bundle: {result.output} ({size_gb:.2f} GB)")
        print(f"archive sha256: {result.archive_sha256}")
        print(f"manifest sha256: {result.manifest_sha256}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
