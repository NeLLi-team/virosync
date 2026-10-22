from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

import virosync.utils.resource_manifest as resource_manifest
from virosync.utils.resource_manifest import (
    CORE_RESOURCE_FILES,
    GENE_FUNCTION_RESOURCE_FILES,
    LEGACY_SEMANTIC_COUNT_KEYS,
    RESOURCE_MANIFEST_NAME,
    RUNTIME_RESOURCE_FILES,
    SEMANTIC_COUNT_KEYS,
    SOURCE_RESOURCE_FILES,
    ResourceManifestError,
    build_resource_manifest,
    build_split_resource_manifests,
    load_resource_manifest,
    runtime_resource_files,
    semantic_count_keys,
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
        "models/model_annotations_with_interpro.tsv": (b"model\tannotation\nVS000001\tone\nVS000002\ttwo\n"),
        "models/og_marker_name_map.tsv": (b"model\tmarker\nVS000001\tOG1\nVS000002\tOG2\n"),
        "marker/marker.faa": b">one\nMPEP\n>two\nMPEP\n>three\nMPEP\n",
        "marker/marker.dmnd": b"synthetic marker diamond\n",
        "genomes/combined_proteome.dmnd": b"synthetic proteome diamond\n",
        "taxonomy/labels.tsv": (b"genome\tlineage\ngenome-one\tNCLDV\ngenome-two\tPLV\n"),
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


def _add_pfam_hmm(root: Path) -> None:
    path = root / "models/pfam_virosync_screening.hmm"
    path.write_bytes(b"HMMER3/f\nNAME  PfamOne\nGA    10.0 10.0;\n//\n")


def test_manifest_schema_and_fast_validation_use_no_child_process(
    tmp_path: Path,
) -> None:
    root, manifest_sha256 = _write_tree(tmp_path / "virosync")
    assert manifest_sha256 == "3c3976abfea9dc7e75bc8491a7c125a6519f79856de923cacc06da4262c47c2b"

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
    assert tuple(manifest.semantic_counts) == LEGACY_SEMANTIC_COUNT_KEYS
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
    entry = next(item for item in document["files"] if item["path"] == "models/model_annotations_with_interpro.tsv")
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

    with pytest.raises(ResourceManifestError, match="semantic count mismatch for hmm_models"):
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


def test_legacy_six_count_manifest_is_accepted_but_new_counts_are_not(
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

    for key in ("pfam_models", "unknown_count"):
        document["semantic_counts"][key] = 1
        manifest_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
        with pytest.raises(ResourceManifestError, match=f"unexpected=.*{key}"):
            load_resource_manifest(root)
        document["semantic_counts"].pop(key)


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
    _add_pfam_hmm(root)

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
    assert set(RUNTIME_RESOURCE_FILES) | set(SOURCE_RESOURCE_FILES) == {
        *CORE_RESOURCE_FILES,
        "models/pfam_virosync_screening.hmm",
    }
    assert runtime.semantic_counts == source.semantic_counts
    assert tuple(runtime.semantic_counts) == SEMANTIC_COUNT_KEYS
    assert runtime.semantic_counts["pfam_models"] == 1
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
            expected_runtime_manifest_sha256=(runtime.manifest_sha256 if kind == "source" else None),
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
    _add_pfam_hmm(root)

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


GENE_FUNCTION_LABELS_HEADER = b"profile\tsource\tlabel\tmodule\tmax_evalue\tdescription\n"


def _add_gene_function_files(root: Path, profile_names: tuple[str, str] = ("PfamOne", "VOG00001")) -> None:
    first, second = profile_names
    (root / "models/gene_function.hmm").write_bytes(
        b"HMMER3/f\nNAME  PfamOne\n//\nNAME  VOG00001\n//\n",
    )
    (root / "models/gene_function_labels.tsv").write_bytes(
        GENE_FUNCTION_LABELS_HEADER
        + f"{first}\tPfam 38.2\tMTase\tother\t1e-05\tone\n".encode()
        + f"{second}\tVOGDB 228\tMCP\tcapsid\t1e-09\ttwo\n".encode()
    )


def _diamond_count(payload: Path | bytes) -> int:
    assert isinstance(payload, Path)
    return 3 if payload.name == "marker.dmnd" else 5


def _runtime_tree(tmp_path: Path, version: str) -> tuple[Path, resource_manifest.ResourceManifest, bytes]:
    """Build a full tree with the gene-function files; return its extracted runtime view."""
    full = tmp_path / "full" / "virosync"
    for relative, content in _payloads(version).items():
        destination = full / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    _add_pfam_hmm(full)
    _add_gene_function_files(full)
    (runtime, runtime_bytes), _ = build_split_resource_manifests(
        full,
        version,
        diamond_sequence_counter=_diamond_count,
    )
    extracted = tmp_path / "runtime" / "virosync"
    for item in runtime.files:
        destination = extracted / item.path
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes((full / item.path).read_bytes())
    (extracted / RESOURCE_MANIFEST_NAME).write_bytes(runtime_bytes)
    return extracted, runtime, runtime_bytes


@pytest.mark.parametrize("version", ["v1.0.7", "v1.1.0", "v1.1.99"])
def test_runtime_resource_files_before_v1_2_0_keep_the_nine_file_list(version: str) -> None:
    assert runtime_resource_files(version) == RUNTIME_RESOURCE_FILES
    assert len(RUNTIME_RESOURCE_FILES) == 9


@pytest.mark.parametrize("version", ["v1.2.0", "v1.10.0", "v2.0.0"])
def test_runtime_resource_files_from_v1_2_0_add_the_gene_function_files(version: str) -> None:
    assert runtime_resource_files(version) == (
        *RUNTIME_RESOURCE_FILES,
        "models/gene_function.hmm",
        "models/gene_function_labels.tsv",
    )
    assert GENE_FUNCTION_RESOURCE_FILES == ("models/gene_function.hmm", "models/gene_function_labels.tsv")


@pytest.mark.parametrize("version", ["v1.0.7", "v1.1.0"])
def test_semantic_count_keys_before_v1_2_0_keep_the_nine_keys(version: str) -> None:
    assert semantic_count_keys(version) == SEMANTIC_COUNT_KEYS


def test_semantic_count_keys_from_v1_2_0_add_the_gene_function_counts() -> None:
    assert semantic_count_keys("v1.2.0") == (*SEMANTIC_COUNT_KEYS, "gene_function_models", "gene_function_labels")


def test_v1_2_0_runtime_manifest_lists_the_gene_function_files_with_their_own_roles(tmp_path: Path) -> None:
    extracted, runtime, _ = _runtime_tree(tmp_path, "v1.2.0")

    roles = {item.path: item.role for item in runtime.files}
    assert tuple(roles) == runtime_resource_files("v1.2.0")
    assert roles["models/gene_function.hmm"] == "gene_function_models"
    assert roles["models/gene_function_labels.tsv"] == "gene_function_labels"
    assert runtime.semantic_counts["gene_function_models"] == 2
    assert runtime.semantic_counts["gene_function_labels"] == 2
    assert load_resource_manifest(extracted) == runtime
    assert validate_resource_tree(extracted).files_verified == 11


def test_v1_2_0_source_manifest_carries_the_same_counts_as_the_runtime_manifest(tmp_path: Path) -> None:
    full = tmp_path / "full" / "virosync"
    for relative, content in _payloads("v1.2.0").items():
        destination = full / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    _add_pfam_hmm(full)
    _add_gene_function_files(full)

    (runtime, _), (source, _) = build_split_resource_manifests(
        full,
        "v1.2.0",
        diamond_sequence_counter=_diamond_count,
    )

    assert tuple(item.path for item in source.files) == SOURCE_RESOURCE_FILES
    assert source.semantic_counts == runtime.semantic_counts


def test_runtime_manifest_before_v1_2_0_ignores_gene_function_files_in_the_tree(tmp_path: Path) -> None:
    _, runtime, _ = _runtime_tree(tmp_path, "v1.1.0")

    assert tuple(item.path for item in runtime.files) == RUNTIME_RESOURCE_FILES
    assert tuple(runtime.semantic_counts) == SEMANTIC_COUNT_KEYS


def test_building_a_v1_2_0_manifest_requires_the_gene_function_files(tmp_path: Path) -> None:
    full = tmp_path / "full" / "virosync"
    for relative, content in _payloads("v1.2.0").items():
        destination = full / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    _add_pfam_hmm(full)

    with pytest.raises(ResourceManifestError, match="gene_function_models"):
        build_split_resource_manifests(full, "v1.2.0", diamond_sequence_counter=_diamond_count)


def test_building_a_v1_2_0_manifest_rejects_unequal_gene_function_counts(tmp_path: Path) -> None:
    full = tmp_path / "full" / "virosync"
    for relative, content in _payloads("v1.2.0").items():
        destination = full / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    _add_pfam_hmm(full)
    _add_gene_function_files(full)
    (full / "models/gene_function_labels.tsv").write_bytes(
        GENE_FUNCTION_LABELS_HEADER + b"PfamOne\tPfam 38.2\tMTase\tother\t1e-05\tone\n"
    )

    with pytest.raises(ResourceManifestError, match="gene-function model and label counts differ"):
        build_split_resource_manifests(full, "v1.2.0", diamond_sequence_counter=_diamond_count)


def test_building_a_v1_2_0_manifest_rejects_labels_for_other_profiles(tmp_path: Path) -> None:
    full = tmp_path / "full" / "virosync"
    for relative, content in _payloads("v1.2.0").items():
        destination = full / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    _add_pfam_hmm(full)
    _add_gene_function_files(full, ("PfamOne", "VOG99999"))

    with pytest.raises(ResourceManifestError, match="gene-function profile identifiers disagree"):
        build_split_resource_manifests(full, "v1.2.0", diamond_sequence_counter=_diamond_count)


def test_validation_rejects_labels_for_other_profiles(tmp_path: Path) -> None:
    extracted, _, _ = _runtime_tree(tmp_path, "v1.2.0")
    labels = extracted / "models/gene_function_labels.tsv"
    labels.write_bytes(labels.read_bytes().replace(b"VOG00001", b"VOG99999"))
    manifest_path = extracted / RESOURCE_MANIFEST_NAME
    document = json.loads(manifest_path.read_text())
    entry = next(item for item in document["files"] if item["path"] == "models/gene_function_labels.tsv")
    entry["sha256"] = resource_manifest.sha256_file(labels)
    manifest_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")

    with pytest.raises(ResourceManifestError, match="gene-function profile identifiers disagree"):
        validate_resource_tree(extracted)


def test_validation_rejects_a_gene_function_count_that_differs_from_the_manifest(tmp_path: Path) -> None:
    extracted, _, _ = _runtime_tree(tmp_path, "v1.2.0")
    manifest_path = extracted / RESOURCE_MANIFEST_NAME
    document = json.loads(manifest_path.read_text())
    document["semantic_counts"]["gene_function_models"] = 3
    document["semantic_counts"]["gene_function_labels"] = 3
    manifest_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")

    with pytest.raises(ResourceManifestError, match="semantic count mismatch for gene_function_models"):
        validate_resource_tree(extracted)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda doc: doc["files"].pop(), "payload set"),
        (lambda doc: doc["files"].pop(9), "payload set"),
        (lambda doc: doc["semantic_counts"].pop("gene_function_labels"), "semantic_counts keys"),
        (lambda doc: doc["files"][9].update(role="hmm_models"), "role mismatch"),
        (lambda doc: doc["semantic_counts"].update(gene_function_labels=3), "gene-function model and label counts"),
    ],
)
def test_v1_2_0_runtime_manifest_rejects_a_broken_gene_function_contract(
    tmp_path: Path,
    mutation,
    message: str,
) -> None:
    extracted, _, runtime_bytes = _runtime_tree(tmp_path, "v1.2.0")
    document = json.loads(runtime_bytes)
    mutation(document)
    manifest_path = tmp_path / "mutated-manifest.json"
    manifest_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")

    with pytest.raises(ResourceManifestError, match=message):
        load_resource_manifest(manifest_path)


def test_a_v1_1_0_manifest_cannot_carry_the_gene_function_files(tmp_path: Path) -> None:
    _, _, runtime_bytes = _runtime_tree(tmp_path, "v1.2.0")
    document = json.loads(runtime_bytes)
    document.update(version="v1.1.0", resource_version="v1.1.0")
    manifest_path = tmp_path / "mutated-manifest.json"
    manifest_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")

    with pytest.raises(ResourceManifestError, match="payload set.*unexpected=.*gene_function"):
        load_resource_manifest(manifest_path)
