from __future__ import annotations

import json
from pathlib import Path
import subprocess

import pytest

import virosync.utils.resource_manifest as resource_manifest
from virosync.utils.resource_manifest import (
    CORE_RESOURCE_FILES,
    RUNTIME_RESOURCE_FILES,
    RESOURCE_MANIFEST_NAME,
    SEMANTIC_COUNT_KEYS,
    SOURCE_RESOURCE_FILES,
    ResourceManifestError,
    build_resource_manifest,
    build_split_resource_manifests,
    load_resource_manifest,
    validate_resource_tree,
)


def _payloads(version: str = "v1.0.6") -> dict[str, bytes]:
    return {
        "DB_VERSION": f"{version}\n".encode(),
        "DATABASE_README.txt": b"synthetic resources\n",
        "models/combined.hmm": (b"HMMER3/f\nNAME  VS000001\n//\nNAME  VS000002\n//\n"),
        "models/combined.hmm.h3f": b"h3f\n",
        "models/combined.hmm.h3i": b"h3i\n",
        "models/combined.hmm.h3m": b"h3m\n",
        "models/combined.hmm.h3p": b"h3p\n",
        "models/model_annotations_with_interpro.tsv": (
            b"model\tannotation\nVS000001\tone\nVS000002\ttwo\n"
        ),
        "models/og_marker_name_map.tsv": (
            b"model\tmarker\nVS000001\tOG1\nVS000002\tOG2\n"
        ),
        "marker/marker.faa": b">one\nMPEP\n>two\nMPEP\n>three\nMPEP\n",
        "marker/marker.dmnd": b"synthetic marker diamond\n",
        "genomes/combined_proteome.dmnd": b"synthetic proteome diamond\n",
        "taxonomy/labels.tsv": (
            b"genome\tlineage\ngenome-one\tNCLDV\ngenome-two\tPLV\n"
        ),
    }


def _write_tree(root: Path, version: str = "v1.0.6") -> tuple[Path, str]:
    for relative, content in _payloads(version).items():
        destination = root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)

    def diamond_count(payload: Path | bytes) -> int:
        assert isinstance(payload, Path)
        return 3 if payload.name == "marker.dmnd" else 5

    manifest, content = build_resource_manifest(
        root,
        version,
        diamond_sequence_counter=diamond_count,
    )
    (root / RESOURCE_MANIFEST_NAME).write_bytes(content)
    return root, manifest.manifest_sha256


def test_manifest_schema_and_fast_validation_use_no_child_process(
    tmp_path: Path,
) -> None:
    root, manifest_sha256 = _write_tree(tmp_path / "virosync")
    assert (
        manifest_sha256
        == "3c3976abfea9dc7e75bc8491a7c125a6519f79856de923cacc06da4262c47c2b"
    )

    def reject_runner(*_args, **_kwargs):
        raise AssertionError("fast validation must not invoke a child process")

    result = validate_resource_tree(
        root,
        expected_version="v1.0.6",
        expected_manifest_sha256=manifest_sha256,
        command_runner=reject_runner,
    )
    manifest = load_resource_manifest(root)

    assert result.version == "v1.0.6"
    assert result.manifest_sha256 == manifest_sha256
    assert result.files_verified == 13
    assert result.full is False
    assert tuple(item.path for item in manifest.files) == CORE_RESOURCE_FILES
    assert tuple(manifest.semantic_counts) == SEMANTIC_COUNT_KEYS
    assert manifest.semantic_counts == {
        "hmm_models": 2,
        "hmm_index_files": 4,
        "model_annotations": 2,
        "og_marker_mappings": 2,
        "marker_proteins": 3,
        "marker_diamond_sequences": 3,
        "proteome_diamond_sequences": 5,
        "taxonomy_labels": 2,
    }


def test_metadata_only_validation_does_not_hash_or_scan_payloads(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root, manifest_sha256 = _write_tree(tmp_path / "virosync")

    def reject(*_args, **_kwargs):
        raise AssertionError("metadata-only validation read payload contents")

    monkeypatch.setattr(resource_manifest, "sha256_file", reject)
    monkeypatch.setattr(resource_manifest, "compute_semantic_counts", reject)

    result = validate_resource_tree(
        root,
        expected_manifest_sha256=manifest_sha256,
        verify_hashes=False,
        full=False,
    )

    assert result.manifest_sha256 == manifest_sha256


@pytest.mark.parametrize("kind", ["file", "directory", "symlink"])
def test_unexpected_resource_tree_paths_are_rejected(
    tmp_path: Path,
    kind: str,
) -> None:
    root, _ = _write_tree(tmp_path / "virosync")
    unexpected = root / "unexpected"
    if kind == "file":
        unexpected.write_text("extra\n")
    elif kind == "directory":
        unexpected.mkdir()
    else:
        unexpected.symlink_to("DB_VERSION")

    with pytest.raises(ResourceManifestError, match="unexpected"):
        validate_resource_tree(root, verify_hashes=False)


def test_hmm_and_annotation_identifiers_must_match(tmp_path: Path) -> None:
    root, _ = _write_tree(tmp_path / "virosync")
    annotation_path = root / "models/model_annotations_with_interpro.tsv"
    original = annotation_path.read_text()
    replacement = original.replace("VS000002", "VS999999")
    annotation_path.write_text(replacement)
    manifest_path = root / RESOURCE_MANIFEST_NAME
    document = json.loads(manifest_path.read_text())
    entry = next(
        item
        for item in document["files"]
        if item["path"] == "models/model_annotations_with_interpro.tsv"
    )
    entry["sha256"] = resource_manifest.sha256_file(annotation_path)
    manifest_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")

    with pytest.raises(ResourceManifestError, match="identifiers disagree"):
        validate_resource_tree(root)


def test_full_validation_checks_both_diamond_database_counts(tmp_path: Path) -> None:
    root, manifest_sha256 = _write_tree(tmp_path / "virosync")
    calls: list[list[str]] = []

    def fake_run(command, **kwargs):
        calls.append(list(command))
        count = 3 if Path(command[-1]).name == "marker.dmnd" else 5
        return subprocess.CompletedProcess(
            command,
            0,
            stdout=f"Database type  Diamond database\nSequences  {count}\n",
            stderr="",
        )

    result = validate_resource_tree(
        root,
        expected_manifest_sha256=manifest_sha256,
        full=True,
        command_runner=fake_run,
    )

    assert result.full is True
    assert [command[:3] for command in calls] == [
        ["diamond", "dbinfo", "--db"],
        ["diamond", "dbinfo", "--db"],
    ]


def test_marker_fasta_and_diamond_counts_must_match(tmp_path: Path) -> None:
    root, _ = _write_tree(tmp_path / "virosync")
    manifest_path = root / RESOURCE_MANIFEST_NAME
    document = json.loads(manifest_path.read_text(encoding="utf-8"))
    document["semantic_counts"]["marker_diamond_sequences"] = 2
    manifest_path.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ResourceManifestError, match="marker.*count"):
        load_resource_manifest(root)


def test_same_size_corruption_is_rejected(tmp_path: Path) -> None:
    root, _ = _write_tree(tmp_path / "virosync")
    target = root / "models/og_marker_name_map.tsv"
    original = target.read_bytes()
    target.write_bytes(bytes([original[0] ^ 1]) + original[1:])

    with pytest.raises(ResourceManifestError, match="SHA-256 mismatch"):
        validate_resource_tree(root)


def test_text_semantic_count_mismatch_is_rejected(tmp_path: Path) -> None:
    root, _ = _write_tree(tmp_path / "virosync")
    manifest_path = root / RESOURCE_MANIFEST_NAME
    document = json.loads(manifest_path.read_text())
    document["semantic_counts"]["hmm_models"] = 99
    manifest_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")

    with pytest.raises(
        ResourceManifestError, match="semantic count mismatch for hmm_models"
    ):
        validate_resource_tree(root)


def test_payload_symlink_is_rejected_even_when_content_matches(tmp_path: Path) -> None:
    root, _ = _write_tree(tmp_path / "virosync")
    target = root / "marker/marker.faa"
    real = root / "marker/real-marker.faa"
    target.rename(real)
    target.symlink_to(real.name)

    with pytest.raises(ResourceManifestError, match="symlink"):
        validate_resource_tree(root)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda doc: doc["files"].pop(), "payload set"),
        (
            lambda doc: doc["semantic_counts"].pop("taxonomy_labels"),
            "semantic_counts keys",
        ),
        (lambda doc: doc.update(schema_version=2), "schema_version"),
        (lambda doc: doc.update(version="v9.9.9"), "version disagreement"),
    ],
)
def test_malformed_or_incomplete_manifest_is_rejected(
    tmp_path: Path,
    mutation,
    message: str,
) -> None:
    root, _ = _write_tree(tmp_path / "virosync")
    manifest_path = root / RESOURCE_MANIFEST_NAME
    document = json.loads(manifest_path.read_text())
    mutation(document)
    manifest_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")

    with pytest.raises(ResourceManifestError, match=message):
        load_resource_manifest(root)


def test_expected_version_and_manifest_digest_are_enforced(tmp_path: Path) -> None:
    root, _ = _write_tree(tmp_path / "virosync")

    with pytest.raises(ResourceManifestError, match="resource version mismatch"):
        load_resource_manifest(root, expected_version="v1.0.7")
    with pytest.raises(ResourceManifestError, match="manifest SHA-256 mismatch"):
        load_resource_manifest(root, expected_manifest_sha256="0" * 64)


def test_legacy_six_count_manifest_is_accepted_but_unknown_counts_are_not(
    tmp_path: Path,
) -> None:
    root, _ = _write_tree(tmp_path / "virosync")
    manifest_path = root / RESOURCE_MANIFEST_NAME
    document = json.loads(manifest_path.read_text())
    document["semantic_counts"].pop("hmm_index_files")
    document["semantic_counts"].pop("og_marker_mappings")
    manifest_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")

    result = validate_resource_tree(root)
    assert set(result.semantic_counts) == {
        "hmm_models",
        "model_annotations",
        "marker_proteins",
        "marker_diamond_sequences",
        "proteome_diamond_sequences",
        "taxonomy_labels",
    }

    document["semantic_counts"]["unknown_count"] = 1
    manifest_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    with pytest.raises(ResourceManifestError, match="unexpected=.*unknown_count"):
        load_resource_manifest(root)


def test_noncanonical_payload_role_is_rejected(tmp_path: Path) -> None:
    root, _ = _write_tree(tmp_path / "virosync")
    manifest_path = root / RESOURCE_MANIFEST_NAME
    document = json.loads(manifest_path.read_text())
    document["files"][0]["role"] = "noncanonical"
    manifest_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")

    with pytest.raises(ResourceManifestError, match="role mismatch for DB_VERSION"):
        load_resource_manifest(root)


def test_schema_v2_manifests_are_strict_bound_views_of_union(
    tmp_path: Path,
) -> None:
    root = tmp_path / "full" / "virosync"
    for relative, content in _payloads().items():
        destination = root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)

    def diamond_count(payload: Path | bytes) -> int:
        assert isinstance(payload, Path)
        return 3 if payload.name == "marker.dmnd" else 5

    (runtime, runtime_bytes), (source, source_bytes) = build_split_resource_manifests(
        root,
        "v1.0.6",
        diamond_sequence_counter=diamond_count,
    )

    assert runtime.schema_version == source.schema_version == 2
    assert runtime.bundle_kind == "runtime"
    assert source.bundle_kind == "source"
    assert tuple(item.path for item in runtime.files) == RUNTIME_RESOURCE_FILES
    assert tuple(item.path for item in source.files) == SOURCE_RESOURCE_FILES
    assert set(RUNTIME_RESOURCE_FILES).isdisjoint(SOURCE_RESOURCE_FILES)
    assert set(RUNTIME_RESOURCE_FILES) | set(SOURCE_RESOURCE_FILES) == set(
        CORE_RESOURCE_FILES
    )
    assert runtime.semantic_counts == source.semantic_counts
    assert tuple(runtime.semantic_counts) == SEMANTIC_COUNT_KEYS
    assert source.runtime_manifest_sha256 == runtime.manifest_sha256

    runtime_path = tmp_path / "runtime-manifest.json"
    source_path = tmp_path / "source-manifest.json"
    runtime_path.write_bytes(runtime_bytes)
    source_path.write_bytes(source_bytes)
    assert load_resource_manifest(runtime_path) == runtime
    assert (
        load_resource_manifest(
            source_path,
            expected_runtime_manifest_sha256=runtime.manifest_sha256,
        )
        == source
    )
    with pytest.raises(ResourceManifestError, match="runtime manifest SHA-256"):
        load_resource_manifest(
            source_path,
            expected_runtime_manifest_sha256="0" * 64,
        )

    for kind, payload_files, manifest_bytes in (
        ("runtime", RUNTIME_RESOURCE_FILES, runtime_bytes),
        ("source", SOURCE_RESOURCE_FILES, source_bytes),
    ):
        extracted = tmp_path / kind / "virosync"
        for relative in payload_files:
            destination = extracted / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes((root / relative).read_bytes())
        (extracted / RESOURCE_MANIFEST_NAME).write_bytes(manifest_bytes)
        result = validate_resource_tree(
            extracted,
            expected_runtime_manifest_sha256=(
                runtime.manifest_sha256 if kind == "source" else None
            ),
        )
        assert result.files_verified == len(payload_files)


@pytest.mark.parametrize(
    ("manifest_kind", "mutation", "message"),
    [
        (
            "runtime",
            lambda doc: doc["files"].pop(),
            "payload set",
        ),
        (
            "source",
            lambda doc: doc["semantic_counts"].pop("hmm_index_files"),
            "semantic_counts keys",
        ),
        (
            "source",
            lambda doc: doc.update(runtime_manifest_sha256="not-a-digest"),
            "runtime_manifest_sha256",
        ),
        (
            "runtime",
            lambda doc: doc.update(runtime_manifest_sha256="0" * 64),
            "extra=.*runtime_manifest_sha256",
        ),
        (
            "runtime",
            lambda doc: doc["files"][5].update(role="diamond_database"),
            "role mismatch",
        ),
        (
            "runtime",
            lambda doc: doc.update(bundle_kind=[]),
            "bundle_kind",
        ),
    ],
)
def test_schema_v2_rejects_incomplete_or_unexpected_contract_fields(
    tmp_path: Path,
    manifest_kind: str,
    mutation,
    message: str,
) -> None:
    root = tmp_path / "source" / "virosync"
    for relative, content in _payloads().items():
        destination = root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)

    def diamond_count(payload: Path | bytes) -> int:
        assert isinstance(payload, Path)
        return 3 if payload.name == "marker.dmnd" else 5

    (_, runtime_bytes), (_, source_bytes) = build_split_resource_manifests(
        root,
        "v1.0.6",
        diamond_sequence_counter=diamond_count,
    )
    content = runtime_bytes if manifest_kind == "runtime" else source_bytes
    document = json.loads(content)
    mutation(document)
    manifest_path = tmp_path / "mutated-manifest.json"
    manifest_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")

    with pytest.raises(ResourceManifestError, match=message):
        load_resource_manifest(manifest_path)
