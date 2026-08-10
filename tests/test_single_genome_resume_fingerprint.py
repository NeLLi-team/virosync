"""Fingerprint-strict schema-v3 resume and legacy compatibility coverage.

Covers the four acceptance areas of the v1.0.5 resume hardening:
  A1  a missing config_fingerprint means STALE unless the legacy opt-in is set
  A2  every completion path (incl. early exits) records the fingerprint
  A3  run identity separates scalar config, resources, environment, and runtime
      controls, and its allowlists cover every FIELD_SPEC
  A4  v1/v2 manifests stay readable only through the internal compatibility opt-in
"""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import pytest

import virosync.output_contract as output_contract
from virosync.ablation import ABLATION_CONTRACT_SHA256
from virosync.config import MaskingBackend, MaskingConfig
from virosync.config.pipeline_config import FIELD_SPECS
from virosync.pipeline.phase0.masking import mask_genome_pipeline
from virosync.orchestration._flows.single_genome import (
    orchestrator as orchestrator_module,
)
from virosync.orchestration._flows.single_genome import run_state as run_state_module
from virosync.orchestration._flows.single_genome.manifest import (
    _FINGERPRINT_CONFIG_FIELDS,
    _FINGERPRINT_ENVIRONMENT_FIELDS,
    _FINGERPRINT_INPUT_FIELDS,
    _FINGERPRINT_RESOURCE_FIELDS,
    _FINGERPRINT_RUNTIME_ONLY_FIELDS,
    _compute_config_fingerprint,
    _write_completion_manifest,
    _write_empty_run_log,
)
from virosync.orchestration._flows.single_genome.orchestrator import (
    _enabled_executable_identities,
    _enabled_model_identities,
    _enabled_resource_identities,
    _masking_request_identity,
)
from virosync.orchestration._flows.single_genome.resume import (
    _completed_run_artifacts,
    _manifest_is_stale,
    _valid_completion_manifest,
)
from virosync.orchestration._flows.single_genome.run_state import (
    build_code_identity,
    build_environment_identity,
    build_input_identity,
    canonical_sha256,
    compute_run_fingerprint,
)
from test_single_genome_resume import (
    _publish_schema3_success as _publish_closed_schema3_success,
    _start_schema3_run,
)


def _seed_outputs(output_dir: Path) -> None:
    """Write the minimal valid run.log + prediction tables (no manifest yet)."""
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "virosync_predictions.tsv").write_text("eve_id\n")
    (output_dir / "virosync_predictions_detailed.tsv").write_text("eve_id\n")
    (output_dir / "run.log").write_text(
        "# ViroSync Run Log: demo\n\n## Results Summary\nGEVEs detected: 0\n"
    )
    _seed_masking_status(output_dir)


def _seed_masking_status(output_dir: Path) -> None:
    input_fasta = output_dir / "input.fna"
    input_fasta.parent.mkdir(parents=True, exist_ok=True)
    input_fasta.write_text(">demo\nACGT\n")
    mask_genome_pipeline(
        input_fasta,
        output_dir / "phase0" / "masking",
        config=MaskingConfig(),
    )


def _run_identity(
    *,
    seed: str = "demo",
    genome_id: str = "demo",
    input_path: Path | None = None,
    output_dir: Path | None = None,
    config_sha256: str | None = None,
    environment: dict | None = None,
    resources: list[dict] | None = None,
    coordinate_schema_version: int | None = None,
    coordinate_convention: str | None = None,
    output_schema_version: int | None = None,
    summary_schema_version: int = 3,
) -> dict:
    digest = canonical_sha256({"seed": seed})
    if output_dir is None:
        output_dir = Path("/tmp") / "virosync-schema3-fixtures" / seed
    if input_path is None:
        input_path = output_dir / "input.fna"
    if environment is None:
        environment_payload = {
            "lock_sha256": digest,
            "runtime_sha256": canonical_sha256(
                {"seed": seed, "kind": "runtime"}
            ),
            "requested_device": "cpu",
            "effective_device": "cpu",
        }
        environment = {
            **environment_payload,
            "sha256": canonical_sha256(environment_payload),
        }
    return {
        "genome_id": genome_id,
        "input_path": str(input_path.absolute()),
        "output_dir": str(output_dir.absolute()),
        "input": {"size": 4, "sha256": digest},
        "config": {
            "sha256": config_sha256 or digest,
            "ablation_id": "A0",
            "ablation_contract_sha256": ABLATION_CONTRACT_SHA256,
        },
        "code": {"version": "test", "source_sha256": digest},
        "environment": environment,
        "coordinate_schema_version": (
            output_contract.COORDINATE_SCHEMA_VERSION
            if coordinate_schema_version is None
            else coordinate_schema_version
        ),
        "coordinate_convention": (
            output_contract.COORDINATE_CONVENTION
            if coordinate_convention is None
            else coordinate_convention
        ),
        "output_schema_version": (
            output_contract.OUTPUT_SCHEMA_VERSION
            if output_schema_version is None
            else output_schema_version
        ),
        "summary_schema_version": summary_schema_version,
        "requested_masking": _masking_request_identity(MaskingConfig()),
        "resources": resources or [],
    }


def _publish_schema3_success(output_dir: Path, *, seed: str = "demo") -> str:
    run_fingerprint = _start_schema3_run(output_dir, identity_seed=seed)
    (output_dir / "virosync_predictions.tsv").write_text("eve_id\n")
    _publish_closed_schema3_success(output_dir, run_fingerprint)
    return run_fingerprint


# --- A1: missing fingerprint => stale unless legacy opt-in --------------------

def test_resume_missing_fingerprint_is_stale_under_expected(tmp_path: Path) -> None:
    _seed_outputs(tmp_path)
    _write_completion_manifest(tmp_path, genome_id="demo", status="success")  # no fp
    assert _completed_run_artifacts(tmp_path, expected_fingerprint="abc123") is None


def test_resume_missing_fingerprint_accepted_with_legacy_optin(tmp_path: Path) -> None:
    _seed_outputs(tmp_path)
    _write_completion_manifest(tmp_path, genome_id="demo", status="success")  # no fp
    assert (
        _completed_run_artifacts(
            tmp_path,
            expected_fingerprint="abc123",
            allow_missing_fingerprint=True,
        )
        is None
    )
    assert (
        _completed_run_artifacts(
            tmp_path,
            expected_fingerprint="abc123",
            allow_missing_fingerprint=True,
            allow_legacy_schema=True,
        )
        is not None
    )


def test_resume_matching_fingerprint_accepted(tmp_path: Path) -> None:
    run_fingerprint = _publish_schema3_success(tmp_path)

    assert (
        _completed_run_artifacts(
            tmp_path,
            expected_fingerprint=run_fingerprint,
        )
        is not None
    )


def test_resume_differing_fingerprint_is_stale(tmp_path: Path) -> None:
    run_fingerprint = _publish_schema3_success(tmp_path)
    different = canonical_sha256({"different": run_fingerprint})

    assert _completed_run_artifacts(tmp_path, expected_fingerprint=different) is None
    assert (
        _completed_run_artifacts(
            tmp_path,
            expected_fingerprint=different,
            allow_missing_fingerprint=True,
            allow_legacy_schema=True,
        )
        is None
    )


def test_schema2_without_expected_fingerprint_is_default_stale(tmp_path: Path) -> None:
    _seed_outputs(tmp_path)
    _write_completion_manifest(tmp_path, genome_id="demo", status="success")  # no fp

    assert _completed_run_artifacts(tmp_path) is None
    assert (
        _completed_run_artifacts(
            tmp_path,
            allow_missing_fingerprint=True,
            allow_legacy_schema=True,
        )
        is not None
    )


# --- A2: every completion path writes the fingerprint ------------------------

def test_write_empty_run_log_records_fingerprint(tmp_path: Path) -> None:
    (tmp_path / "virosync_predictions_detailed.tsv").write_text("eve_id\n")
    _seed_masking_status(tmp_path)
    _write_empty_run_log(
        output_dir=tmp_path,
        genome_id="demo",
        reason="no HMM hits",
        elapsed_sec=0.0,
        fingerprint="fp-early-exit",
    )
    payload = json.loads((tmp_path / "virosync_run_complete.json").read_text())
    assert payload["config_fingerprint"] == "fp-early-exit"
    assert payload["schema_version"] == 2
    assert payload["coordinate_schema_version"] == output_contract.COORDINATE_SCHEMA_VERSION
    assert payload["output_schema_version"] == output_contract.OUTPUT_SCHEMA_VERSION
    assert payload["coordinate_convention"] == output_contract.COORDINATE_CONVENTION
    assert not (tmp_path / "virosync_run_state.json").exists()


# --- A3: fingerprint sensitivity ---------------------------------------------

_BASE = {
    "ablation_id": "A0",
    "ablation_contract_sha256": ABLATION_CONTRACT_SHA256,
    "assembly_mode": "default",
    "extension_kb": 5,
    "marker_floor_priority_only": 0.55,
    "boundary_host_trim_enabled": False,
    "boundary_diamond_random_seed": 42,
    "search_backend": "diamond",
    "threads": 8,
    "device": "cuda",
    "resume": True,
    "extended_output": True,
}


def _fp(overrides: dict | None = None) -> str:
    cfg = dict(_BASE)
    if overrides:
        cfg.update(overrides)
    return _compute_config_fingerprint(cfg)


@pytest.mark.parametrize(
    ("field", "new_value"),
    [
        ("coordinate_schema_version", output_contract.COORDINATE_SCHEMA_VERSION + 1),
        ("output_schema_version", output_contract.OUTPUT_SCHEMA_VERSION + 1),
        ("summary_schema_version", 4),
    ],
)
def test_output_contract_identity_changes_fingerprint(
    field: str,
    new_value: int | str,
) -> None:
    baseline_identity = _run_identity()
    changed_identity = dict(baseline_identity)
    changed_identity[field] = new_value

    assert compute_run_fingerprint(changed_identity) != compute_run_fingerprint(
        baseline_identity
    )


def test_noncanonical_coordinate_convention_is_rejected() -> None:
    identity = _run_identity(coordinate_convention="1-based, closed [start, end]")

    with pytest.raises(ValueError, match="coordinate_convention"):
        compute_run_fingerprint(identity)


@pytest.mark.parametrize(
    "field,new",
    [
        ("assembly_mode", "permissive"),
        ("marker_floor_priority_only", 0.6),
        ("boundary_host_trim_enabled", True),
        ("boundary_diamond_random_seed", 7),
        # The backend is hashed even though Diamond is the only valid value.
        ("search_backend", "other"),
        ("extension_kb", 8),
        ("extended_output", False),  # output-artifact schema -> fingerprinted
        ("masking", MaskingConfig(backend=MaskingBackend.TRF)),
    ],
)
def test_output_determining_knob_changes_fingerprint(field, new) -> None:
    assert _fp() != _fp({field: new})


def test_ablation_id_changes_config_and_full_run_fingerprints() -> None:
    baseline_config = _fp()
    ablated_config = _fp({"ablation_id": "A6"})
    baseline_identity = _run_identity(config_sha256=baseline_config)
    ablated_identity = _run_identity(config_sha256=ablated_config)
    ablated_identity["config"]["ablation_id"] = "A6"

    assert baseline_config != ablated_config
    assert compute_run_fingerprint(baseline_identity) != compute_run_fingerprint(
        ablated_identity
    )


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("ablation_id", "A7", "ablation_id must be one of A0-A6"),
        (
            "ablation_contract_sha256",
            "0" * 64,
            "ablation contract does not match",
        ),
    ],
)
def test_run_identity_rejects_invalid_ablation_binding(
    field: str,
    value: str,
    message: str,
) -> None:
    identity = _run_identity()
    identity["config"][field] = value

    with pytest.raises(ValueError, match=message):
        compute_run_fingerprint(identity)


def test_config_fingerprint_is_full_scalar_only_sha256() -> None:
    fingerprint = _fp()

    assert len(fingerprint) == 64
    assert set(fingerprint) <= set("0123456789abcdef")


@pytest.mark.parametrize(
    "field,new",
    [("threads", 32), ("device", "cpu"), ("resume", False), ("max_concurrent_genomes", 2)],
)
def test_non_scalar_knob_does_not_change_config_fingerprint(field, new) -> None:
    assert _fp() == _fp({field: new})


def _run_fingerprint_with_resources(flat_config: dict) -> str:
    resources = [asdict(item) for item in _enabled_resource_identities(flat_config)]
    return compute_run_fingerprint(_run_identity(resources=resources))


def test_environment_device_change_changes_full_run_fingerprint(tmp_path: Path) -> None:
    lock = tmp_path / "pixi.lock"
    lock.write_text("version = 6\n")
    cpu = asdict(
        build_environment_identity(
            lock,
            requested_device="cpu",
            effective_device="cpu",
        )
    )
    cuda = asdict(
        build_environment_identity(
            lock,
            requested_device="cuda",
            effective_device="cuda",
        )
    )

    assert compute_run_fingerprint(
        _run_identity(environment=cpu)
    ) != compute_run_fingerprint(_run_identity(environment=cuda))


def test_same_path_same_size_input_mutation_changes_full_run_fingerprint(
    tmp_path: Path,
) -> None:
    genome = tmp_path / "same.fna"
    genome.write_bytes(b">x\nAAAA\n")
    before = asdict(build_input_identity(genome))
    genome.write_bytes(b">x\nAAAT\n")
    after = asdict(build_input_identity(genome))

    before_identity = _run_identity()
    before_identity["input"] = before
    after_identity = _run_identity()
    after_identity["input"] = after
    assert compute_run_fingerprint(before_identity) != compute_run_fingerprint(
        after_identity
    )


def test_same_size_installed_source_mutation_changes_full_run_fingerprint(
    tmp_path: Path,
) -> None:
    source_root = tmp_path / "virosync"
    source_root.mkdir()
    module = source_root / "module.py"
    module.write_text("VALUE = 1\n")
    before = asdict(build_code_identity(source_root, version="test"))
    module.write_text("VALUE = 2\n")
    after = asdict(build_code_identity(source_root, version="test"))

    before_identity = _run_identity()
    before_identity["code"] = before
    after_identity = _run_identity()
    after_identity["code"] = after
    assert compute_run_fingerprint(before_identity) != compute_run_fingerprint(
        after_identity
    )


def test_code_identity_reads_source_on_each_call(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_root = tmp_path / "virosync"
    source_root.mkdir()
    (source_root / "module.py").write_text("VALUE = 1\n")
    original_hash = run_state_module._hash_relative_file
    calls = 0

    def tracked_hash(root: Path, relative_path: str) -> tuple[int, str]:
        nonlocal calls
        calls += 1
        return original_hash(root, relative_path)

    monkeypatch.setattr(run_state_module, "_hash_relative_file", tracked_hash)
    build_code_identity(source_root, version="test")
    build_code_identity(source_root, version="test")

    assert calls == 2


def test_manifestless_resource_identity_reads_members_on_each_call(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    resource_root = tmp_path / "optional-resource"
    resource_root.mkdir()
    (resource_root / "payload.tsv").write_text("id\tvalue\n1\tA\n")
    original_hash = run_state_module._hash_relative_file
    calls = 0

    def tracked_hash(root: Path, relative_path: str) -> tuple[int, str]:
        nonlocal calls
        calls += 1
        return original_hash(root, relative_path)

    monkeypatch.setattr(run_state_module, "_hash_relative_file", tracked_hash)
    for _ in range(2):
        run_state_module.build_resource_identity(
            "optional",
            "test",
            resource_root,
            kind="optional",
        )

    assert calls == 2


def test_executable_identity_reads_binary_on_each_call(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    binary = tmp_path / "tool"
    binary.write_bytes(b"TOOL0001")
    original_builder = orchestrator_module.build_resource_set_identity
    calls = 0

    def tracked_builder(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original_builder(*args, **kwargs)

    monkeypatch.setattr(
        orchestrator_module,
        "build_resource_set_identity",
        tracked_builder,
    )
    orchestrator_module._executable_path_identity("tool", binary)
    orchestrator_module._executable_path_identity("tool", binary)

    assert calls == 2


def test_db_version_bump_changes_full_run_fingerprint_at_identical_size(
    tmp_path: Path,
) -> None:
    bundle = tmp_path / "virosync"
    (bundle / "marker").mkdir(parents=True)
    db = bundle / "marker" / "marker.dmnd"
    db.write_bytes(b"x" * 100)
    (bundle / "DB_VERSION").write_text("v1.0.5\n")
    fp_old = _run_fingerprint_with_resources({"marker_db": db})
    (bundle / "DB_VERSION").write_text("v1.0.6\n")  # same byte size, new version
    fp_new = _run_fingerprint_with_resources({"marker_db": db})
    assert fp_old != fp_new


def test_same_name_same_size_database_byte_changes_full_run_fingerprint(
    tmp_path: Path,
) -> None:
    bundle = tmp_path / "virosync"
    (bundle / "marker").mkdir(parents=True)
    db = bundle / "marker" / "marker.dmnd"
    db.write_bytes(b"A" * 100)
    (bundle / "DB_VERSION").write_text("v1.0.6\n")
    before = _run_fingerprint_with_resources({"marker_db": db})
    db.write_bytes(b"B" + b"A" * 99)
    after = _run_fingerprint_with_resources({"marker_db": db})

    assert before != after


def test_prebuilt_marker_db_excludes_unused_marker_source_from_fingerprint(
    tmp_path: Path,
) -> None:
    marker_db = tmp_path / "marker.dmnd"
    marker_db.write_bytes(b"prebuilt")
    baseline = _run_fingerprint_with_resources({"marker_db": marker_db})
    with_unused_source = _run_fingerprint_with_resources(
        {
            "marker_db": marker_db,
            "marker_faa_db": tmp_path / "missing-unused-marker.faa",
        }
    )

    assert with_unused_source == baseline


def test_selected_marker_source_remains_in_fingerprint(tmp_path: Path) -> None:
    marker_faa = tmp_path / "marker.faa"
    marker_faa.write_bytes(b">marker\nAAAA\n")
    before = _run_fingerprint_with_resources({"marker_faa_db": marker_faa})
    marker_faa.write_bytes(b">marker\nAAAT\n")
    after = _run_fingerprint_with_resources({"marker_faa_db": marker_faa})

    assert after != before


def test_phylo_db_identity_gated_on_enable_phylogenetic(tmp_path: Path) -> None:
    db = tmp_path / "gvclass.dmnd"
    db.write_bytes(b"y" * 50)
    off_with_db = _run_fingerprint_with_resources(
        {"gvclass_db": db, "enable_phylogenetic": False, "run_gvclass": False}
    )
    off_no_db = _run_fingerprint_with_resources(
        {"gvclass_db": None, "enable_phylogenetic": False, "run_gvclass": False}
    )
    assert off_with_db == off_no_db  # path ignored while phylo is off
    on_with_db = _run_fingerprint_with_resources(
        {"gvclass_db": db, "enable_phylogenetic": True, "run_gvclass": False}
    )
    assert on_with_db != off_with_db


def test_tmvec_database_dir_gated_on_use_tmvec_database(tmp_path: Path) -> None:
    dbdir = tmp_path / "tmvec_db"
    dbdir.mkdir()
    (dbdir / "tmvec.dmnd").write_bytes(b"tmvec")
    # off: tmvec dir change is invisible
    off_a = _run_fingerprint_with_resources(
        {"tmvec_database_dir": dbdir, "use_tmvec_database": False}
    )
    off_b = _run_fingerprint_with_resources(
        {"tmvec_database_dir": None, "use_tmvec_database": False}
    )
    assert off_a == off_b
    on_db = _run_fingerprint_with_resources(
        {"tmvec_database_dir": dbdir, "use_tmvec_database": True}
    )
    on_none = _run_fingerprint_with_resources(
        {"tmvec_database_dir": None, "use_tmvec_database": True}
    )
    assert on_db != on_none


def test_gvclass_db_gated_by_run_gvclass_not_only_phylogenetic(tmp_path: Path) -> None:
    db = tmp_path / "gvclass.dmnd"
    db.write_bytes(b"z" * 40)
    base = {"enable_phylogenetic": False, "run_gvclass": True}
    with_db = _run_fingerprint_with_resources({**base, "gvclass_db": db})
    without_db = _run_fingerprint_with_resources({**base, "gvclass_db": None})
    assert with_db != without_db


def test_viral_structure_db_gated_by_use_boltz(tmp_path: Path) -> None:
    db = tmp_path / "struct.db"
    db.write_bytes(b"s" * 30)
    base = {"use_boltz": True}
    with_db = _run_fingerprint_with_resources({**base, "viral_structure_db": db})
    without_db = _run_fingerprint_with_resources(
        {**base, "viral_structure_db": None}
    )
    assert with_db != without_db


def _run_fingerprint_for_identity_items(items) -> str:
    return compute_run_fingerprint(
        _run_identity(resources=[asdict(item) for item in items])
    )


def _identity_digest(items, name: str) -> str | None:
    return next(
        (item.manifest_sha256 for item in items if item.name == name),
        None,
    )


def test_gvclass_path_executable_is_feature_gated_and_content_bound(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    binary = tmp_path / "gvclass"
    binary.write_bytes(b"GVCLASS1")
    monkeypatch.setattr(
        orchestrator_module.shutil,
        "which",
        lambda name: str(binary) if name == "gvclass" else None,
    )
    masking = MaskingConfig()

    off_before = _enabled_executable_identities(
        {"enable_phylogenetic": False, "skip_structural": True},
        masking,
    )
    on_before = _enabled_executable_identities(
        {"enable_phylogenetic": True, "skip_structural": True},
        masking,
    )
    binary.write_bytes(b"GVCLASS2")
    off_after = _enabled_executable_identities(
        {"enable_phylogenetic": False, "skip_structural": True},
        masking,
    )
    on_after = _enabled_executable_identities(
        {"enable_phylogenetic": True, "skip_structural": True},
        masking,
    )

    assert _identity_digest(off_before, "executable:gvclass") is None
    assert _run_fingerprint_for_identity_items(off_before) == (
        _run_fingerprint_for_identity_items(off_after)
    )
    assert _identity_digest(on_before, "executable:gvclass") != (
        _identity_digest(on_after, "executable:gvclass")
    )


def test_gvclass_batch_executable_is_run_gvclass_gated(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    binary_dir = tmp_path / "gvclass-install"
    binary_dir.mkdir()
    binary = binary_dir / "gvclass"
    binary.write_bytes(b"BATCH001")
    monkeypatch.setattr(orchestrator_module.shutil, "which", lambda name: None)
    masking = MaskingConfig()

    off = _enabled_executable_identities(
        {
            "run_gvclass": False,
            "gvclass_path": binary_dir,
            "skip_structural": True,
        },
        masking,
    )
    on_before = _enabled_executable_identities(
        {
            "run_gvclass": True,
            "gvclass_path": binary_dir,
            "skip_structural": True,
        },
        masking,
    )
    binary.write_bytes(b"BATCH002")
    on_after = _enabled_executable_identities(
        {
            "run_gvclass": True,
            "gvclass_path": binary_dir,
            "skip_structural": True,
        },
        masking,
    )

    assert _identity_digest(off, "executable:gvclass-batch") is None
    assert _identity_digest(on_before, "executable:gvclass-batch") != (
        _identity_digest(on_after, "executable:gvclass-batch")
    )


def test_boltz_executable_is_feature_gated_and_content_bound(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from virosync.utils import executables

    binary = tmp_path / "boltz"
    binary.write_bytes(b"BOLTZ001")
    monkeypatch.setattr(orchestrator_module.shutil, "which", lambda name: None)
    monkeypatch.setattr(
        executables,
        "resolve_boltz_executable",
        lambda: binary,
    )
    masking = MaskingConfig()

    off_before = _enabled_executable_identities(
        {"use_boltz": False, "skip_structural": True},
        masking,
    )
    on_before = _enabled_executable_identities(
        {"use_boltz": True, "skip_structural": True},
        masking,
    )
    binary.write_bytes(b"BOLTZ002")
    off_after = _enabled_executable_identities(
        {"use_boltz": False, "skip_structural": True},
        masking,
    )
    on_after = _enabled_executable_identities(
        {"use_boltz": True, "skip_structural": True},
        masking,
    )

    assert _identity_digest(off_before, "executable:boltz") is None
    assert _run_fingerprint_for_identity_items(off_before) == (
        _run_fingerprint_for_identity_items(off_after)
    )
    assert _identity_digest(on_before, "executable:boltz") != (
        _identity_digest(on_after, "executable:boltz")
    )


def test_skani_executable_is_always_content_bound(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    binary = tmp_path / "skani"
    binary.write_bytes(b"SKANI001")
    monkeypatch.setattr(
        orchestrator_module.shutil,
        "which",
        lambda name: str(binary) if name == "skani" else None,
    )
    before = _enabled_executable_identities(
        {"skip_structural": True},
        MaskingConfig(),
    )
    binary.write_bytes(b"SKANI002")
    after = _enabled_executable_identities(
        {"skip_structural": True},
        MaskingConfig(),
    )

    assert _identity_digest(before, "executable:skani") != (
        _identity_digest(after, "executable:skani")
    )
    assert _run_fingerprint_for_identity_items(before) != (
        _run_fingerprint_for_identity_items(after)
    )


def test_tmvec_model_revisions_are_feature_gated_and_fingerprinted(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from virosync.pipeline.phase3 import tmvec_predictor

    assert _enabled_model_identities({"use_tmvec_database": False}) == []
    before = _enabled_model_identities({"use_tmvec_database": True})
    monkeypatch.setattr(
        tmvec_predictor,
        "PROTTRANS_MODEL_REVISION",
        "revision-after-test-mutation",
    )
    after = _enabled_model_identities({"use_tmvec_database": True})

    assert _identity_digest(
        before,
        f"model:{tmvec_predictor.PROTTRANS_MODEL_ID}",
    ) != _identity_digest(
        after,
        f"model:{tmvec_predictor.PROTTRANS_MODEL_ID}",
    )
    assert _run_fingerprint_for_identity_items(before) != (
        _run_fingerprint_for_identity_items(after)
    )


@pytest.mark.parametrize(
    "surface",
    ["config", "source", "core resource", "optional resource", "masking"],
)
def test_closed_run_identity_invalidates_each_output_surface(
    surface: str,
) -> None:
    resources = [
        {
            "name": "core",
            "kind": "core",
            "version": "v1.0.6",
            "manifest_sha256": "a" * 64,
        },
        {
            "name": "optional",
            "kind": "optional",
            "version": "fixture",
            "manifest_sha256": "b" * 64,
        },
    ]
    baseline = _run_identity(resources=resources)
    changed = json.loads(json.dumps(baseline))
    if surface == "config":
        changed["config"]["sha256"] = "c" * 64
    elif surface == "source":
        changed["code"]["source_sha256"] = "c" * 64
    elif surface == "core resource":
        changed["resources"][0]["manifest_sha256"] = "c" * 64
    elif surface == "optional resource":
        changed["resources"][1]["manifest_sha256"] = "c" * 64
    else:
        changed["requested_masking"]["backend"] = "trf"

    assert compute_run_fingerprint(baseline) != compute_run_fingerprint(
        changed
    )


def test_fingerprint_inputs_are_real_impl_parameters() -> None:
    """Every fingerprinted flat field must be an actual _single_genome_flow_impl param,
    else locals() filtering yields None and the field is silently dropped from the hash.
    """
    import inspect

    from virosync.orchestration._flows.single_genome.orchestrator import (
        _single_genome_flow_impl,
    )

    params = set(inspect.signature(_single_genome_flow_impl).parameters)
    missing = _FINGERPRINT_INPUT_FIELDS - params
    assert not missing, f"fingerprint reads non-impl params (silently None): {sorted(missing)}"


# --- A3 drift guard: every FIELD_SPEC is classified --------------------------

def test_fingerprint_allowlist_covers_all_field_specs() -> None:
    all_flats = {spec.flat for spec in FIELD_SPECS}
    categories = {
        "config": _FINGERPRINT_CONFIG_FIELDS,
        "resources": _FINGERPRINT_RESOURCE_FIELDS,
        "environment": _FINGERPRINT_ENVIRONMENT_FIELDS,
        "runtime": _FINGERPRINT_RUNTIME_ONLY_FIELDS,
    }
    classified = frozenset().union(*categories.values())
    overlaps = {
        field: sorted(name for name, values in categories.items() if field in values)
        for field in classified
        if sum(field in values for values in categories.values()) > 1
    }
    assert not overlaps, f"fields belong to multiple identity categories: {overlaps}"
    assert _FINGERPRINT_INPUT_FIELDS == (
        _FINGERPRINT_CONFIG_FIELDS
        | _FINGERPRINT_RESOURCE_FIELDS
        | _FINGERPRINT_ENVIRONMENT_FIELDS
    )
    unclassified = all_flats - classified
    assert not unclassified, (
        "unclassified config knobs (triage in manifest.py): "
        f"{sorted(unclassified)}"
    )
    stale = classified - all_flats
    assert not stale, f"classified fields no longer in FIELD_SPECS: {sorted(stale)}"


# --- A4: schema-version backward read + legacy manifests ---------------------

def test_v1_manifest_without_fingerprint_honours_legacy_optin(tmp_path: Path) -> None:
    _seed_outputs(tmp_path)
    (tmp_path / "virosync_run_complete.json").write_text(
        json.dumps({"schema_version": 1, "genome_id": "demo", "status": "success"})
    )
    assert _completed_run_artifacts(tmp_path, expected_fingerprint="x") is None
    assert (
        _completed_run_artifacts(
            tmp_path,
            expected_fingerprint="x",
            allow_missing_fingerprint=True,
        )
        is None
    )
    assert (
        _completed_run_artifacts(
            tmp_path,
            expected_fingerprint="x",
            allow_missing_fingerprint=True,
            allow_legacy_schema=True,
        )
        is not None
    )


def test_unknown_schema_version_rejected(tmp_path: Path) -> None:
    manifest = tmp_path / "virosync_run_complete.json"
    manifest.write_text(
        json.dumps({"schema_version": 99, "genome_id": "demo", "status": "success"})
    )
    assert not _valid_completion_manifest(manifest, expected_fingerprint=None)


# --- CRITICAL: stale COMPLETED run must wipe phase-level caches, not just the -----
# --- top-level short-circuit. _manifest_is_stale drives the genome-dir discard. ----

def test_manifest_is_stale_detects_completed_run_with_differing_fingerprint(
    tmp_path: Path,
) -> None:
    run_fingerprint = _publish_schema3_success(tmp_path)

    assert _manifest_is_stale(
        tmp_path,
        expected_fingerprint=canonical_sha256({"different": run_fingerprint}),
    ) is True


def test_manifest_is_stale_false_for_matching_fingerprint(tmp_path: Path) -> None:
    run_fingerprint = _publish_schema3_success(tmp_path)

    assert (
        _manifest_is_stale(
            tmp_path,
            expected_fingerprint=run_fingerprint,
        )
        is False
    )


def test_schema2_manifest_is_default_stale_even_when_fingerprint_matches(
    tmp_path: Path,
) -> None:
    _seed_masking_status(tmp_path)
    _write_completion_manifest(
        tmp_path,
        genome_id="demo",
        status="success",
        fingerprint="legacy-fingerprint",
    )

    assert _manifest_is_stale(
        tmp_path,
        expected_fingerprint="legacy-fingerprint",
    ) is True


def test_manifest_is_stale_false_for_interrupted_run_without_manifest(tmp_path: Path) -> None:
    # An interrupted run leaves intermediates but no completion manifest -> resumable,
    # NOT a stale completion (so we must not wipe it).
    assert _manifest_is_stale(tmp_path, expected_fingerprint="x") is False


def test_manifest_is_stale_false_for_corrupt_manifest(tmp_path: Path) -> None:
    (tmp_path / "virosync_run_complete.json").write_text("{not json")
    assert _manifest_is_stale(tmp_path, expected_fingerprint="x") is False
