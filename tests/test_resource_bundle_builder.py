from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tarfile

import pytest

from virosync.utils.resource_manifest import (
    CORE_RESOURCE_FILES,
    RUNTIME_RESOURCE_FILES,
    RESOURCE_MANIFEST_NAME,
    SOURCE_RESOURCE_FILES,
    ResourceManifestError,
)

_SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts/build_resource_bundle.py"
_SCRIPT_SPEC = importlib.util.spec_from_file_location(
    "build_resource_bundle", _SCRIPT_PATH
)
assert _SCRIPT_SPEC is not None and _SCRIPT_SPEC.loader is not None
_SCRIPT_MODULE = importlib.util.module_from_spec(_SCRIPT_SPEC)
sys.modules[_SCRIPT_SPEC.name] = _SCRIPT_MODULE
_SCRIPT_SPEC.loader.exec_module(_SCRIPT_MODULE)
build_resource_bundle = _SCRIPT_MODULE.build_resource_bundle
build_split_resource_bundles = _SCRIPT_MODULE.build_split_resource_bundles
TAR_FORMAT = _SCRIPT_MODULE.TAR_FORMAT


def _resource_tree(root: Path) -> Path:
    payloads = {
        "DB_VERSION": b"stale-source-version\n",
        "DATABASE_README.txt": b"stale source readme\n",
        "models/combined.hmm": b"HMMER3/f\nNAME  VS000001\n//\n",
        "models/combined.hmm.h3f": b"h3f\n",
        "models/combined.hmm.h3i": b"h3i\n",
        "models/combined.hmm.h3m": b"h3m\n",
        "models/combined.hmm.h3p": b"h3p\n",
        "models/model_annotations_with_interpro.tsv": (
            b"model\tannotation\nVS000001\tsynthetic\n"
        ),
        "models/og_marker_name_map.tsv": b"model\tmarker\nVS000001\tOG1\n",
        "marker/marker.faa": b">marker-one\nMPEP\n",
        "marker/marker.dmnd": b"marker diamond\n",
        "genomes/combined_proteome.dmnd": b"proteome diamond\n",
        "taxonomy/labels.tsv": b"genome-one\tNCLDV\n",
        "DB_METADATA.json": b'{"install_local": true}\n',
        "unrelated.txt": b"must not be packaged\n",
    }
    for relative, content in payloads.items():
        destination = root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    return root


def _snapshot(root: Path) -> dict[str, tuple[int, int, str]]:
    snapshot = {}
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix()
        metadata = path.stat()
        snapshot[relative] = (
            metadata.st_size,
            metadata.st_mtime_ns,
            hashlib.sha256(path.read_bytes()).hexdigest(),
        )
    return snapshot


def _diamond_count(payload: Path | bytes) -> int:
    assert isinstance(payload, Path)
    return 1 if payload.name == "marker.dmnd" else 2


def test_builder_is_source_immutable_and_archive_is_deterministic(
    tmp_path: Path,
) -> None:
    resources = _resource_tree(tmp_path / "source" / "virosync")
    before = _snapshot(resources)
    output_one = tmp_path / "one.tar.gz"
    output_two = tmp_path / "two.tar.gz"

    first = build_resource_bundle(
        resources,
        output_one,
        "v1.0.6",
        skip_hmmpress=True,
        skip_marker_dmnd=True,
        diamond_sequence_counter=_diamond_count,
    )
    second = build_resource_bundle(
        resources,
        output_two,
        "v1.0.6",
        skip_hmmpress=True,
        skip_marker_dmnd=True,
        diamond_sequence_counter=_diamond_count,
    )

    assert _snapshot(resources) == before
    assert output_one.read_bytes() == output_two.read_bytes()
    assert output_one.read_bytes()[4:8] == b"\0\0\0\0"
    assert first.archive_sha256 == second.archive_sha256
    assert first.manifest_sha256 == second.manifest_sha256

    with tarfile.open(output_one, "r:gz") as archive:
        members = archive.getmembers()
        expected_names = [f"virosync/{path}" for path in CORE_RESOURCE_FILES]
        expected_names.append(f"virosync/{RESOURCE_MANIFEST_NAME}")
        assert [member.name for member in members] == expected_names
        assert all(member.isfile() for member in members)
        assert all(member.mtime == 0 for member in members)
        assert all(member.uid == 0 and member.gid == 0 for member in members)
        assert all(member.mode == 0o644 for member in members)
        assert archive.extractfile("virosync/DB_VERSION").read() == b"v1.0.6\n"
        readme = archive.extractfile("virosync/DATABASE_README.txt").read()
        assert b"ViroSync core resource bundle v1.0.6" in readme
        manifest = json.load(archive.extractfile(f"virosync/{RESOURCE_MANIFEST_NAME}"))

    assert [entry["path"] for entry in manifest["files"]] == list(CORE_RESOURCE_FILES)
    assert RESOURCE_MANIFEST_NAME not in {entry["path"] for entry in manifest["files"]}
    assert "DB_METADATA.json" not in {entry["path"] for entry in manifest["files"]}
    assert manifest["semantic_counts"]["marker_diamond_sequences"] == 1
    assert manifest["semantic_counts"]["proteome_diamond_sequences"] == 2
    assert (resources / "DB_VERSION").read_bytes() == b"stale-source-version\n"


def test_requested_regeneration_fails_cleanly_when_tool_is_absent(
    tmp_path: Path,
) -> None:
    resources = _resource_tree(tmp_path / "source" / "virosync")
    before = _snapshot(resources)

    with pytest.raises(ResourceManifestError, match="hmmpress.*not found"):
        build_resource_bundle(
            resources,
            tmp_path / "bundle.tar.gz",
            "v1.0.6",
            skip_hmmpress=False,
            skip_marker_dmnd=True,
            tool_finder=lambda _name: None,
            diamond_sequence_counter=_diamond_count,
        )

    assert _snapshot(resources) == before
    assert not (tmp_path / "bundle.tar.gz").exists()


def test_output_inside_source_tree_is_rejected_before_writing(tmp_path: Path) -> None:
    resources = _resource_tree(tmp_path / "source" / "virosync")
    output = resources / "bundle.tar.gz"

    with pytest.raises(ResourceManifestError, match="outside the resources directory"):
        build_resource_bundle(
            resources,
            output,
            "v1.0.6",
            skip_hmmpress=True,
            skip_marker_dmnd=True,
            diamond_sequence_counter=_diamond_count,
        )

    assert not output.exists()


def test_builder_accepts_the_stable_virosync_symlink(tmp_path: Path) -> None:
    versioned = _resource_tree(
        tmp_path / "resources" / "virosync-v1.0.6-0123456789abcdef"
    )
    stable = versioned.parent / "virosync"
    stable.symlink_to(versioned.name)
    before = _snapshot(versioned)
    output = tmp_path / "bundle.tar.gz"

    result = build_resource_bundle(
        stable,
        output,
        "v1.0.6",
        skip_hmmpress=True,
        skip_marker_dmnd=True,
        diamond_sequence_counter=_diamond_count,
    )

    assert result.output == output
    assert output.is_file()
    assert _snapshot(versioned) == before


def test_builder_rejects_a_stale_reused_marker_database(tmp_path: Path) -> None:
    resources = _resource_tree(tmp_path / "source" / "virosync")

    def mismatched_diamond_count(payload: Path | bytes) -> int:
        assert isinstance(payload, Path)
        return 99 if payload.name == "marker.dmnd" else 2

    with pytest.raises(
        ResourceManifestError,
        match="marker protein and marker DIAMOND sequence counts differ",
    ):
        build_resource_bundle(
            resources,
            tmp_path / "bundle.tar.gz",
            "v1.0.6",
            skip_hmmpress=True,
            skip_marker_dmnd=True,
            diamond_sequence_counter=mismatched_diamond_count,
        )

    assert not (tmp_path / "bundle.tar.gz").exists()


def test_skip_flags_require_complete_existing_derived_files(tmp_path: Path) -> None:
    resources = _resource_tree(tmp_path / "source" / "virosync")
    (resources / "models/combined.hmm.h3p").unlink()

    with pytest.raises(ResourceManifestError, match="--skip-hmmpress.*missing"):
        build_resource_bundle(
            resources,
            tmp_path / "bundle.tar.gz",
            "v1.0.6",
            skip_hmmpress=True,
            skip_marker_dmnd=True,
            diamond_sequence_counter=_diamond_count,
        )


def test_archive_format_supports_real_proteome_database_size() -> None:
    member = tarfile.TarInfo("virosync/genomes/combined_proteome.dmnd")
    member.size = 11_269_683_342

    with pytest.raises(ValueError, match="overflow"):
        member.tobuf(format=tarfile.USTAR_FORMAT)
    encoded = member.tobuf(format=TAR_FORMAT)
    assert len(encoded) == tarfile.BLOCKSIZE


def test_split_builder_is_deterministic_exact_and_bound(tmp_path: Path) -> None:
    resources = _resource_tree(tmp_path / "source" / "virosync")
    before = _snapshot(resources)
    first_runtime = tmp_path / "first-runtime.tar.gz"
    first_source = tmp_path / "first-source.tar.gz"
    second_runtime = tmp_path / "second-runtime.tar.gz"
    second_source = tmp_path / "second-source.tar.gz"

    first = build_split_resource_bundles(
        resources,
        first_runtime,
        first_source,
        "v1.0.6",
        skip_hmmpress=True,
        skip_marker_dmnd=True,
        diamond_sequence_counter=_diamond_count,
    )
    second = build_split_resource_bundles(
        resources,
        second_runtime,
        second_source,
        "v1.0.6",
        skip_hmmpress=True,
        skip_marker_dmnd=True,
        diamond_sequence_counter=_diamond_count,
    )

    assert _snapshot(resources) == before
    assert first_runtime.read_bytes() == second_runtime.read_bytes()
    assert first_source.read_bytes() == second_source.read_bytes()
    assert first.runtime.archive_sha256 == second.runtime.archive_sha256
    assert first.source.archive_sha256 == second.source.archive_sha256

    manifests = {}
    for kind, archive_path, expected_payloads in (
        ("runtime", first_runtime, RUNTIME_RESOURCE_FILES),
        ("source", first_source, SOURCE_RESOURCE_FILES),
    ):
        with tarfile.open(archive_path, "r:gz") as archive:
            expected_names = [f"virosync/{relative}" for relative in expected_payloads]
            expected_names.append(f"virosync/{RESOURCE_MANIFEST_NAME}")
            assert [member.name for member in archive.getmembers()] == expected_names
            manifests[kind] = json.load(
                archive.extractfile(f"virosync/{RESOURCE_MANIFEST_NAME}")
            )

    runtime_manifest = manifests["runtime"]
    source_manifest = manifests["source"]
    with tarfile.open(first_runtime, "r:gz") as archive:
        runtime_readme = archive.extractfile(
            "virosync/DATABASE_README.txt"
        ).read()
    assert b"Source/repair artifact:" in runtime_readme
    assert runtime_manifest["schema_version"] == 2
    assert runtime_manifest["bundle_kind"] == "runtime"
    assert "runtime_manifest_sha256" not in runtime_manifest
    assert source_manifest["schema_version"] == 2
    assert source_manifest["bundle_kind"] == "source"
    assert source_manifest["runtime_manifest_sha256"] == first.runtime.manifest_sha256
    assert runtime_manifest["semantic_counts"] == source_manifest["semantic_counts"]
    assert runtime_manifest["semantic_counts"]["marker_proteins"] == 1
    assert runtime_manifest["semantic_counts"]["marker_diamond_sequences"] == 1


def test_split_prep_failure_publishes_neither_artifact(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    resources = _resource_tree(tmp_path / "source" / "virosync")
    runtime_output = tmp_path / "runtime.tar.gz"
    source_output = tmp_path / "source.tar.gz"
    runtime_output.write_bytes(b"existing runtime\n")
    source_output.write_bytes(b"existing source\n")
    real_create = _SCRIPT_MODULE.create_deterministic_archive
    calls = 0

    def fail_second_archive(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise ResourceManifestError("synthetic source prep failure")
        return real_create(*args, **kwargs)

    monkeypatch.setattr(
        _SCRIPT_MODULE,
        "create_deterministic_archive",
        fail_second_archive,
    )

    with pytest.raises(ResourceManifestError, match="source prep failure"):
        build_split_resource_bundles(
            resources,
            runtime_output,
            source_output,
            "v1.0.6",
            skip_hmmpress=True,
            skip_marker_dmnd=True,
            diamond_sequence_counter=_diamond_count,
        )

    assert runtime_output.read_bytes() == b"existing runtime\n"
    assert source_output.read_bytes() == b"existing source\n"
