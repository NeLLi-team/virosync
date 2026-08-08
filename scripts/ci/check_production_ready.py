#!/usr/bin/env python3
"""Validate ViroSync public-release surfaces from the locked Pixi environment."""

from __future__ import annotations

import os
import re
import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from virosync.config.application_config import ApplicationConfig
from virosync.config.pipeline_config import ConfigError
from virosync.utils.database_manager import ViroSyncDatabaseManager
from virosync.utils.resource_manifest import (
    ResourceManifestError,
    load_resource_manifest,
)
SOFTWARE_VERSION = "1.0.0"
DATABASE_VERSION = "v1.0.7"
RESOURCE_ARCHIVE = "resources_v1_0_7_runtime.tar.gz"
RESOURCE_URL = f"https://dl.newlineages.com/virosync/{RESOURCE_ARCHIVE}"
RESOURCE_ARCHIVE_SHA256 = (
    "57daed0b39bf2bc4c4f84ec3b612c6034a3d26ea38e7ec5fba4f4469da36e9a2"
)
RESOURCE_MANIFEST_SHA256 = (
    "f3aeed77045f4728207c6997f5986ed155056e2b4b2a297574d57686982a18b3"
)
RELEASE_MANIFEST_PATH = Path(
    "release-manifests/resources_v1_0_7/RESOURCE_MANIFEST.json"
)
RESOURCE_IDENTITY = (
    RESOURCE_URL,
    DATABASE_VERSION,
    RESOURCE_ARCHIVE_SHA256,
    RESOURCE_MANIFEST_SHA256,
)

FORBIDDEN_TRACKED_PATTERNS = (
    "resources/",
    "results/",
    ".memd/",
    "tasks/",
    "docs/syn2_real_benchmark_results_writeup.md",
    "docs/syn2_real_benchmark_results_writeup.qmd",
    "config/orchestration_vs791_benchmark.yaml",
    "tests/test_repeat_spike_benchmark.py",
    "benchmarking/",
)

EXPECTED_IGNORES = (
    "/benchmarking/",
    "/docs/ms/",
    "resources/",
    "results/",
    "tasks/",
    ".memd/",
    "resources_v*.tar.gz",
    "config/orchestration_*_benchmark.yaml",
)


def read_text(path: str) -> str:
    return (ROOT / path).read_text()


def read_toml(path: str) -> dict:
    with (ROOT / path).open("rb") as handle:
        return tomllib.load(handle)


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def require_contains(path: str, needle: str, failures: list[str]) -> None:
    text = read_text(path)
    require(
        needle in text, f"{path} does not contain expected text: {needle}", failures
    )


def require_regex(path: str, pattern: str, message: str, failures: list[str]) -> None:
    text = read_text(path)
    require(
        re.search(pattern, text, re.MULTILINE | re.DOTALL) is not None,
        message,
        failures,
    )


def check_software_version(failures: list[str]) -> None:
    pyproject = read_toml("pyproject.toml")
    pixi = read_toml("pixi.toml")

    require(
        pyproject.get("project", {}).get("version") == SOFTWARE_VERSION,
        f"pyproject.toml project.version must be {SOFTWARE_VERSION}",
        failures,
    )
    require(
        pixi.get("workspace", {}).get("version") == SOFTWARE_VERSION,
        f"pixi.toml workspace.version must be {SOFTWARE_VERSION}",
        failures,
    )
    require_regex(
        "src/virosync/__init__.py",
        rf'^__version__ = "{re.escape(SOFTWARE_VERSION)}"$',
        f"src/virosync/__init__.py __version__ must be {SOFTWARE_VERSION}",
        failures,
    )
    require_contains(
        "README.md",
        f"https://img.shields.io/badge/version-{SOFTWARE_VERSION}-blue",
        failures,
    )

    changelog = read_text("CHANGELOG.md")
    first_release = re.search(
        r"^## \[([^\]]+)\] - \d{4}-\d{2}-\d{2}$", changelog, re.MULTILINE
    )
    require(
        first_release is not None, "CHANGELOG.md is missing a release header", failures
    )
    if first_release is not None:
        require(
            first_release.group(1) == SOFTWARE_VERSION,
            f"CHANGELOG.md latest release must be {SOFTWARE_VERSION}",
            failures,
        )


def resource_identity(config: ApplicationConfig) -> tuple[object, ...]:
    """Return the typed core-resource identity from one shipped config."""

    orchestration = config.orchestration
    return (
        orchestration.core_resources_url,
        orchestration.core_resources_version,
        orchestration.core_resources_sha256,
        orchestration.core_resources_manifest_sha256,
    )


def check_resource_version(failures: list[str]) -> None:
    config_paths = (
        "config/orchestration.yaml",
        "config/orchestration_archaeal.yaml",
    )
    parsed_configs: dict[str, ApplicationConfig] = {}
    for relative in config_paths:
        try:
            parsed_configs[relative] = ApplicationConfig.from_yaml(ROOT / relative)
        except (ConfigError, OSError) as exc:
            failures.append(f"{relative} does not parse as a typed config: {exc}")

    for relative, config in parsed_configs.items():
        require(
            resource_identity(config) == RESOURCE_IDENTITY,
            f"{relative} core-resource identity differs from the release pin",
            failures,
        )
    if set(parsed_configs) == set(config_paths):
        require(
            resource_identity(parsed_configs[config_paths[0]])
            == resource_identity(parsed_configs[config_paths[1]]),
            "shipped configs disagree on the core-resource identity",
            failures,
        )

    require(
        bool(ViroSyncDatabaseManager.DATABASE_SOURCES),
        "DATABASE_SOURCES has no authenticated source",
        failures,
    )
    if ViroSyncDatabaseManager.DATABASE_SOURCES:
        source = ViroSyncDatabaseManager.DATABASE_SOURCES[0]
        manager_identity = (
            source.get("source"),
            source.get("version"),
            source.get("archive_sha256"),
            source.get("manifest_sha256"),
        )
        require(
            manager_identity == RESOURCE_IDENTITY,
            "DATABASE_SOURCES first source differs from the shipped config identity",
            failures,
        )
        require(
            source.get("filename") == RESOURCE_ARCHIVE,
            f"DATABASE_SOURCES first filename must be {RESOURCE_ARCHIVE}",
            failures,
        )
        require(
            ViroSyncDatabaseManager.DATABASE_VERSION == DATABASE_VERSION,
            f"database manager version must be {DATABASE_VERSION}",
            failures,
        )

    try:
        manifest = load_resource_manifest(
            ROOT / RELEASE_MANIFEST_PATH,
            expected_version=DATABASE_VERSION,
            expected_manifest_sha256=RESOURCE_MANIFEST_SHA256,
        )
    except (OSError, ResourceManifestError) as exc:
        failures.append(f"tracked release manifest is invalid: {exc}")
        manifest = None
    if manifest is not None:
        require(
            manifest.schema_version == 2,
            "tracked release manifest schema must be 2",
            failures,
        )
        require(
            manifest.bundle_kind == "runtime",
            "tracked release manifest must describe a runtime bundle",
            failures,
        )
        require(
            len(manifest.files) == 9,
            "tracked release manifest must authenticate 9 payload files",
            failures,
        )
        require(
            manifest.semantic_counts.get("hmm_models") == 1053,
            "tracked release manifest HMM count must be 1053",
            failures,
        )
        require(
            manifest.semantic_counts.get("marker_proteins") == 3864300,
            "tracked release manifest marker count must be 3864300",
            failures,
        )
        require(
            manifest.semantic_counts.get("pfam_models") == 937,
            "tracked release manifest Pfam count must be 937",
            failures,
        )
        require(
            manifest.semantic_counts.get("taxonomy_labels") == 190317,
            "tracked release manifest taxonomy count must be 190317",
            failures,
        )

    manager = read_text("src/virosync/utils/database_manager.py")
    require(
        "resources_v1_0_3.tar.gz" not in manager
        and "resources_v1_0_2.tar.gz" not in manager
        and "resources_v1_0_1.tar.gz" not in manager
        and "resources_v1_0_0.tar.gz" not in manager,
        "DATABASE_SOURCES must not retain legacy v1.0.0-v1.0.3 fallback archives",
        failures,
    )
    require(
        "combined_ga.hmm" not in manager,
        "database_manager.py must not accept legacy models/combined_ga.hmm resources",
        failures,
    )

    for path in ("README.md", "docs/METHODS.md", "docs/RESOURCE_BUNDLE.md"):
        require_contains(path, DATABASE_VERSION, failures)
        require_contains(path, RESOURCE_ARCHIVE, failures)

    # The smoke workflow runs on the self-hosted runner with resources
    # provisioned out-of-band, so it no longer caches by key; it must still pin
    # the expected DB version (asserted at runtime against the provisioned bundle).
    require_regex(
        ".github/workflows/example-smoke.yml",
        rf"^\s*EXPECTED_DB_VERSION:\s*{re.escape(DATABASE_VERSION)}\s*$",
        "example-smoke EXPECTED_DB_VERSION differs from the release pin",
        failures,
    )


def check_resource_documentation(failures: list[str]) -> None:
    bundle_doc = read_text("docs/RESOURCE_BUNDLE.md")
    expected_values = (
        "HMM models: `1053`",
        "Pfam models: `937`",
        "marker annotation table lines: `1054`",
        "OG marker map lines: `809`",
        "taxonomy label lines: `190317`",
        "models/model_annotations_with_interpro.tsv",
        "models/og_marker_name_map.tsv",
    )
    for value in expected_values:
        require(
            value in bundle_doc, f"docs/RESOURCE_BUNDLE.md missing {value}", failures
        )


def is_forbidden_tracked_path(path: str) -> bool:
    """Return whether *path* is outside the curated release allowlist."""

    # The benchmark harness, its results, and the manuscript live in the
    # virosync-bench repository; nothing under these prefixes belongs here.
    if path.startswith("benchmarking/"):
        return True
    if path.startswith("docs/ms/"):
        return True

    return any(
        path == pattern.rstrip("/") or path.startswith(pattern)
        for pattern in FORBIDDEN_TRACKED_PATTERNS
    ) or re.fullmatch(r"resources_v.*\.tar\.gz", path) is not None


def check_artifact_exclusion(failures: list[str]) -> None:
    gitignore = read_text(".gitignore")
    for pattern in EXPECTED_IGNORES:
        require(pattern in gitignore, f".gitignore is missing {pattern}", failures)

    result = subprocess.run(
        ["git", "ls-files", "-s", "-z"],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        failures.append(f"git ls-files failed: {result.stderr.strip()}")
        return

    tracked_entries: list[tuple[str, str]] = []
    for entry in result.stdout.split("\0"):
        if not entry:
            continue
        metadata, separator, path = entry.partition("\t")
        if not separator:
            failures.append(f"git ls-files emitted an invalid entry: {entry!r}")
            continue
        mode = metadata.split(maxsplit=1)[0]
        tracked_entries.append((mode, path))
    forbidden = [path for _, path in tracked_entries if is_forbidden_tracked_path(path)]
    require(
        not forbidden,
        "Forbidden generated/private paths are tracked: " + ", ".join(forbidden),
        failures,
    )
    unsafe_modes = [
        f"{path} ({mode})"
        for mode, path in tracked_entries
        if mode not in {"100644", "100755"}
    ]
    require(
        not unsafe_modes,
        "Benchmark/manuscript source contains a symlink or submodule: "
        + ", ".join(unsafe_modes),
        failures,
    )


def check_github_tag(failures: list[str]) -> None:
    ref_type = os.environ.get("GITHUB_REF_TYPE")
    ref_name = os.environ.get("GITHUB_REF_NAME")
    if ref_type != "tag":
        return

    expected_tag = f"v{SOFTWARE_VERSION}"
    require(
        ref_name == expected_tag,
        f"GitHub tag {ref_name!r} must match package version tag {expected_tag!r}",
        failures,
    )


def main() -> int:
    failures: list[str] = []
    check_software_version(failures)
    check_resource_version(failures)
    check_resource_documentation(failures)
    check_artifact_exclusion(failures)
    check_github_tag(failures)

    if failures:
        print("Production readiness checks failed:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1

    print(
        f"Production readiness checks passed for ViroSync {SOFTWARE_VERSION} / DB {DATABASE_VERSION}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
