"""Single-genome pipeline orchestrator.

Owns the public ``single_genome_flow`` entry point, its implementation
``_single_genome_flow_impl``, the Phase 0 subflow, and compatibility exports.
Threads per-genome state across the phase subflows imported from the sibling
phase modules.
"""

import inspect
import json
import logging
import os
import shutil
import subprocess
from dataclasses import asdict
from functools import wraps
from pathlib import Path
from typing import Callable, Optional

import virosync
from virosync.ablation import (
    ABLATION_CONTRACT_SHA256,
    AblationCounters,
    AblationEvents,
    AblationID,
    InterventionCounts,
    validate_ablation_events_bytes,
)
from virosync.config import MaskingBackend, MaskingConfig, PipelineConfig
from virosync.output_contract import (
    COORDINATE_CONVENTION,
    COORDINATE_SCHEMA_VERSION,
    OUTPUT_SCHEMA_VERSION,
    EFFECTIVE_EVE_CLASS_COUNT_KEYS,
)
from virosync.orchestration.runtime import call_task, get_orchestration_logger
from virosync.orchestration.tasks import (
    mask_genome_task,
    generate_proteome_task,
    generate_outputs_task,
    create_summary_artifact_task,
)
from virosync.pipeline.phase0.masking import (
    MaskingResult,
    load_masking_result,
)
from virosync.orchestration._flows.utils import (
    _detect_explicit_overrides,
    _filter_kwargs_to_signature,
    _merge_config_with_kwargs,
)
from virosync.utils.atomic_write import atomic_write_context
from virosync.utils.path_safety import require_strict_child, validate_path_component
from virosync.utils.resource_installer import verified_install_receipt
from virosync.validation.tsv_invariants import (
    TSVInvariantError,
    enforce_tsv_invariants,
)

from .loaders import (
    _count_fasta_records,
)
from .manifest import (
    _FINGERPRINT_INPUT_FIELDS,
    _FINGERPRINT_RESOURCE_FIELDS,
    _FINGERPRINT_RESOURCE_GATED,
    _clear_success_markers,
    _compute_config_fingerprint,
    _empty_prediction_summary,
    _summarize_prediction_outputs,
    _write_completion_manifest,
)
from .reports import _generate_required_reports
from .resume import (
    _completed_run_artifacts,
)
from .run_state import (
    ConfigIdentity,
    EnvironmentIdentity,
    PHASE_MARKER_FILENAMES,
    RUN_STATE_FILENAME,
    ResourceIdentity,
    build_artifact_identity,
    build_code_identity,
    build_environment_identity,
    build_input_identity,
    build_phase_record,
    build_resource_identity,
    build_resource_set_identity,
    canonical_sha256,
    compute_run_fingerprint,
    invalidate_from_phase,
    load_run_state,
    marker_sha256,
    plan_resume,
    publish_phase_completion,
    publish_run_failed,
    publish_run_started,
    publish_run_success,
    runtime_environment_sha256,
    sibling_run_lock,
)
from .phase_state import PHASE2_STATE_SCHEMA
from .phase1_state import PHASE1_STATE_SCHEMA
from .phase2_resume_state import PHASE2_RESUME_STATE_SCHEMA
from .phase1 import _run_phase1_subflow
from .phase2 import _run_phase2_subflow
from .phase3 import _run_phase3_subflow

logger = logging.getLogger(__name__)

_RUN_SUMMARY_SCHEMA_VERSION = 3
_VERIFIED_CORE_INVENTORIES: dict[Path, tuple[object, ...]] = {}
_EXECUTABLE_IDENTITY_CACHE: dict[
    tuple[str, Path],
    tuple[tuple[int, int, int, int, int], object],
] = {}

_FINAL_ROOT_NAMES = frozenset(
    {
        "ablation_events.json",
        "run.log",
        "virosync_run_complete.json",
        "virosync_predictions.tsv",
        "virosync_predictions.bed",
        "virosync_predictions.gff3",
        "virosync_predictions_detailed.tsv",
        "virosync_summary.json",
        "virosync_tsv_invariant_report.tsv",
        "gvclass_results.tsv",
        "host_signature_model.png",
    }
)


def _is_owned_root_final(path: Path) -> bool:
    return (
        path.name in _FINAL_ROOT_NAMES
        or path.name.endswith("_eves.fna")
    )


def _masking_request_identity(masking: MaskingConfig) -> dict[str, object]:
    """Return the canonical pre-run masking request."""

    library = masking.repeatmasker_library
    library_sha256 = None
    if library is not None and Path(library).is_file():
        library_sha256 = build_input_identity(Path(library)).sha256
    return {
        "backend": masking.backend.value,
        "failure_policy": masking.failure_policy.value,
        "fallback_backend": (
            masking.fallback_backend.value
            if masking.fallback_backend is not None
            else None
        ),
        "repeatmasker_species": masking.repeatmasker_species,
        "repeatmasker_library": str(library) if library is not None else None,
        "repeatmasker_library_sha256": library_sha256,
    }


def _project_lock_path() -> Path | None:
    """Find the Pixi lock for a source checkout, if one is installed."""

    candidates: list[Path] = []
    pixi_project_root = os.environ.get("PIXI_PROJECT_ROOT")
    if pixi_project_root:
        candidates.append(Path(pixi_project_root))
    candidates.append(Path(virosync.__file__).resolve().parent)
    seen: set[Path] = set()
    for start in candidates:
        for directory in (start, *start.parents):
            if directory in seen:
                continue
            seen.add(directory)
            lock = directory / "pixi.lock"
            if lock.is_file():
                return lock
    return None


def _environment_identity(device: str) -> EnvironmentIdentity:
    """Bind Pixi when available, with an installed-distribution fallback."""

    import platform

    effective_device = f"cpu:{platform.machine() or 'unknown'}"
    if device.lower() == "cuda":
        try:
            import torch

            if torch.cuda.is_available():
                index = torch.cuda.current_device()
                properties = torch.cuda.get_device_properties(index)
                major, minor = torch.cuda.get_device_capability(index)
                driver_version = "unknown"
                nvidia_smi = shutil.which("nvidia-smi")
                if nvidia_smi is not None:
                    try:
                        probe = subprocess.run(
                            [
                                nvidia_smi,
                                "--query-gpu=driver_version",
                                "--format=csv,noheader",
                            ],
                            check=True,
                            capture_output=True,
                            text=True,
                            timeout=5,
                        )
                    except (OSError, subprocess.SubprocessError):
                        pass
                    else:
                        versions = sorted(
                            {
                                line.strip()
                                for line in probe.stdout.splitlines()
                                if line.strip()
                            }
                        )
                        if versions:
                            driver_version = "+".join(versions)
                cudnn_version = torch.backends.cudnn.version() or "unavailable"
                effective_device = (
                    f"cuda:{properties.name}:sm_{major}{minor}:"
                    f"runtime_{torch.version.cuda or 'unknown'}:"
                    f"driver_{driver_version}:cudnn_{cudnn_version}"
                )
            else:
                effective_device = "cuda:unavailable"
        except Exception:
            effective_device = "cuda:unavailable"

    lock_path = _project_lock_path()
    if lock_path is not None:
        return build_environment_identity(
            lock_path,
            requested_device=device,
            effective_device=effective_device,
        )

    from importlib.metadata import distributions

    installed = sorted(
        (
            distribution.metadata.get("Name", "unknown"),
            distribution.version,
        )
        for distribution in distributions()
    )
    lock_sha256 = canonical_sha256(
        {"schema_version": 1, "installed_distributions": installed}
    )
    payload = {
        "lock_sha256": lock_sha256,
        "runtime_sha256": runtime_environment_sha256(),
        "requested_device": device,
        "effective_device": effective_device,
    }
    return EnvironmentIdentity(**payload, sha256=canonical_sha256(payload))


def _nearest_resource_version(path: Path) -> str:
    start = path if path.is_dir() else path.parent
    for directory in (start, *start.parents):
        version_file = directory / "DB_VERSION"
        if not version_file.is_file():
            continue
        try:
            value = version_file.read_text(encoding="utf-8").strip()
        except (OSError, UnicodeError):
            continue
        if value:
            return value
    return "unversioned"


def _authenticated_core_manifest(path: Path):
    """Return the R7 manifest only when *path* is one authenticated payload."""

    from virosync.utils.resource_manifest import (
        RESOURCE_MANIFEST_NAME,
        ResourceManifestError,
        load_resource_manifest,
    )

    resolved = path.resolve(strict=True)
    start = resolved if resolved.is_dir() else resolved.parent
    for directory in (start, *start.parents):
        manifest_path = directory / RESOURCE_MANIFEST_NAME
        if not manifest_path.is_file():
            continue
        try:
            manifest = load_resource_manifest(manifest_path)
        except ResourceManifestError as exc:
            raise ValueError(f"invalid resource manifest above {path}: {exc}") from exc
        if resolved.is_file():
            try:
                relative = resolved.relative_to(directory.resolve(strict=True)).as_posix()
            except ValueError:
                continue
            if relative in {item.path for item in manifest.files}:
                _verify_core_manifest_payload(directory, manifest_path, manifest)
                return manifest_path, manifest
        return None
    return None


def _core_inventory_signature(root: Path, manifest) -> tuple[object, ...]:
    signature: list[object] = [manifest.manifest_sha256]
    for item in manifest.files:
        candidate = root / item.path
        metadata = candidate.lstat()
        if candidate.is_symlink() or not candidate.is_file():
            raise ValueError(f"core resource payload is not regular: {item.path}")
        signature.append(
            (
                item.path,
                metadata.st_dev,
                metadata.st_ino,
                metadata.st_size,
                metadata.st_mtime_ns,
                metadata.st_ctime_ns,
            )
        )
    return tuple(signature)


def _verify_core_manifest_payload(root: Path, manifest_path: Path, manifest) -> None:
    signature = _core_inventory_signature(root, manifest)
    cache_key = manifest_path.resolve(strict=True)
    if _VERIFIED_CORE_INVENTORIES.get(cache_key) == signature:
        return
    if not verified_install_receipt(root, manifest):
        from virosync.utils.resource_manifest import validate_resource_tree

        validate_resource_tree(
            root,
            expected_version=manifest.version,
            expected_manifest_sha256=manifest.manifest_sha256,
            verify_hashes=True,
            full=False,
        )
    _VERIFIED_CORE_INVENTORIES[cache_key] = signature


def _enabled_resource_identities(flat_config: dict) -> list:
    """Build identities for every enabled resource-valued configuration field."""

    gated = dict(_FINGERPRINT_RESOURCE_GATED)
    marker_source_fields = {"faa_dir", "marker_faa_db", "marker_faa_dir"}
    prebuilt_marker_selected = bool(flat_config.get("marker_db"))
    identities = []
    seen: set[tuple[str, str, str, str]] = set()
    for field in sorted(_FINGERPRINT_RESOURCE_FIELDS):
        if prebuilt_marker_selected and field in marker_source_fields:
            continue
        gates = gated.get(field)
        if gates is not None and not any(bool(flat_config.get(gate)) for gate in gates):
            continue
        value = flat_config.get(field)
        if not value:
            continue
        selected = Path(value).resolve(strict=True)
        authenticated = _authenticated_core_manifest(selected)
        if authenticated is not None:
            manifest_path, manifest = authenticated
            identity = build_resource_identity(
                "core",
                manifest.version,
                manifest_path,
                kind="core",
            )
        elif field == "hmm_database" and selected.is_file():
            dependency_paths = {selected.name: selected}
            for suffix in (".h3f", ".h3i", ".h3m", ".h3p"):
                sidecar = selected.with_name(f"{selected.name}{suffix}")
                if sidecar.is_file():
                    dependency_paths[sidecar.name] = sidecar
            for sibling_name in (
                "model_annotations_with_interpro.tsv",
                "og_marker_name_map.tsv",
            ):
                sibling = selected.parent / sibling_name
                if sibling.is_file():
                    dependency_paths[sibling.name] = sibling
            identity = build_resource_set_identity(
                field,
                _nearest_resource_version(selected),
                dependency_paths,
            )
        elif field == "viral_structure_db" and selected.is_file():
            dependency_paths = {
                sibling.name: sibling
                for sibling in selected.parent.glob(f"{selected.name}*")
                if sibling.is_file() and not sibling.is_symlink()
            }
            identity = build_resource_set_identity(
                field,
                _nearest_resource_version(selected),
                dependency_paths,
            )
        else:
            identity = build_resource_identity(
                field,
                _nearest_resource_version(selected),
                selected,
                kind="optional",
            )
        key = (
            identity.name,
            identity.kind,
            identity.version,
            identity.manifest_sha256,
        )
        if key not in seen:
            identities.append(identity)
            seen.add(key)
    return identities


def _executable_path_identity(name: str, executable: str | Path):
    path = Path(executable).resolve(strict=True)
    metadata = path.stat()
    signature = (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
    )
    cache_key = (name, path)
    cached = _EXECUTABLE_IDENTITY_CACHE.get(cache_key)
    if cached is not None and cached[0] == signature:
        return cached[1]
    identity = build_resource_set_identity(
        f"executable:{name}",
        "resolved-binary",
        {path.name: path},
    )
    _EXECUTABLE_IDENTITY_CACHE[cache_key] = (signature, identity)
    return identity


def _executable_resource_identity(name: str):
    executable = shutil.which(name)
    if executable is None:
        return None
    return _executable_path_identity(name, executable)


def _enabled_executable_identities(
    flat_config: dict,
    masking: MaskingConfig,
) -> list:
    names = {"diamond", "prodigal-gv", "skani"}
    if masking.backend in {MaskingBackend.TRF, MaskingBackend.TRF_REPEATMASKER}:
        names.add("trf")
    if masking.backend in {
        MaskingBackend.REPEATMASKER,
        MaskingBackend.TRF_REPEATMASKER,
    }:
        names.add("RepeatMasker")
    if not bool(flat_config.get("skip_structural")):
        names.add("foldseek")
    if bool(flat_config.get("enable_phylogenetic")):
        names.add("gvclass")
    identities = [
        identity
        for name in sorted(names)
        if (identity := _executable_resource_identity(name)) is not None
    ]
    if bool(flat_config.get("run_gvclass")) and flat_config.get("gvclass_path"):
        gvclass = Path(flat_config["gvclass_path"]) / "gvclass"
        if gvclass.is_file() and not gvclass.is_symlink():
            identities.append(
                _executable_path_identity("gvclass-batch", gvclass)
            )
    if bool(flat_config.get("use_boltz")):
        from virosync.utils.executables import resolve_boltz_executable

        boltz = resolve_boltz_executable()
        if boltz is not None:
            identities.append(_executable_path_identity("boltz", boltz))
    if masking.backend in {
        MaskingBackend.REPEATMASKER,
        MaskingBackend.TRF_REPEATMASKER,
    }:
        repeatmasker = shutil.which("RepeatMasker")
        if repeatmasker is not None:
            executable = Path(repeatmasker).resolve(strict=True)
            for library in (
                executable.parent / "Libraries",
                executable.parent.parent / "share" / "RepeatMasker" / "Libraries",
            ):
                if library.is_dir() and not library.is_symlink():
                    identities.append(
                        build_resource_identity(
                            "repeatmasker-library",
                            _nearest_resource_version(library),
                            library,
                            kind="optional",
                        )
                    )
                    break
    return identities


def _enabled_model_identities(flat_config: dict) -> list[ResourceIdentity]:
    if not bool(flat_config.get("use_tmvec_database")):
        return []
    from virosync.pipeline.phase3.tmvec_predictor import (
        PROTTRANS_MODEL_ID,
        PROTTRANS_MODEL_REVISION,
        PROTTRANS_PROXY_MODEL_ID,
        PROTTRANS_PROXY_MODEL_REVISION,
    )

    return [
        ResourceIdentity(
            name=f"model:{model_id}",
            kind="optional",
            version=revision,
            manifest_sha256=canonical_sha256(
                {
                    "schema_version": 1,
                    "model_id": model_id,
                    "revision": revision,
                }
            ),
        )
        for model_id, revision in (
            (PROTTRANS_MODEL_ID, PROTTRANS_MODEL_REVISION),
            (PROTTRANS_PROXY_MODEL_ID, PROTTRANS_PROXY_MODEL_REVISION),
        )
    ]


def _build_run_identity(
    *,
    genome_path: Path,
    genome_id: str,
    output_dir: Path,
    flat_config: dict,
    masking: MaskingConfig,
    device: str,
) -> tuple[dict[str, object], str]:
    """Build the complete immutable identity and its full run fingerprint."""

    ablation_id = AblationID(flat_config["ablation_id"])
    if flat_config["ablation_contract_sha256"] != ABLATION_CONTRACT_SHA256:
        raise ValueError("ablation contract SHA-256 differs from this ViroSync build")
    source_root = Path(virosync.__file__).resolve().parent
    identity = {
        "genome_id": genome_id,
        "input_path": str(genome_path.resolve(strict=True)),
        "output_dir": str(output_dir.resolve(strict=False)),
        "input": asdict(build_input_identity(genome_path)),
        "config": asdict(
            ConfigIdentity(
                sha256=_compute_config_fingerprint(flat_config),
                ablation_id=ablation_id.value,
                ablation_contract_sha256=ABLATION_CONTRACT_SHA256,
            )
        ),
        "code": asdict(
            build_code_identity(
                source_root,
                version=getattr(virosync, "__version__", "unknown"),
            )
        ),
        "environment": asdict(_environment_identity(device)),
        "resources": [
            asdict(item)
            for item in (
                *_enabled_resource_identities(flat_config),
                *_enabled_executable_identities(flat_config, masking),
                *_enabled_model_identities(flat_config),
            )
        ],
        "coordinate_schema_version": COORDINATE_SCHEMA_VERSION,
        "coordinate_convention": COORDINATE_CONVENTION,
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "summary_schema_version": _RUN_SUMMARY_SCHEMA_VERSION,
        "requested_masking": _masking_request_identity(masking),
    }
    return identity, compute_run_fingerprint(identity)


def _artifact_schema(path: Path, output_dir: Path) -> str:
    relative = path.relative_to(output_dir).as_posix()
    if path.name == "ablation_events.json":
        return "virosync.ablation_events/v1"
    if path.name == "virosync_predictions.tsv":
        return "canonical-predictions-v4"
    if path.name == "virosync_predictions_detailed.tsv":
        return "detailed-predictions-v4"
    if path.name == "virosync_predictions.bed":
        return "canonical-predictions-bed-v1"
    if path.name == "virosync_predictions.gff3":
        return "canonical-predictions-gff3-v1"
    if path.name == "eve_ani_edges.tsv":
        return "eve-ani-edges-v1"
    if path.name == "virosync_summary.json":
        return "virosync-summary-v3"
    if path.name == "masking_status.json":
        return "masking-status-v1"
    if path.name == "refined_state.json":
        return PHASE2_STATE_SCHEMA
    if relative == "phase1/resume_state.json":
        return PHASE1_STATE_SCHEMA
    if relative == "phase2/resume_state.json":
        return PHASE2_RESUME_STATE_SCHEMA
    if path.name == "virosync_tsv_invariant_report.tsv":
        return "tsv-invariant-report-v1"
    if path.name == "virosync_run_complete.json":
        return "completion-manifest-v2"
    if relative == "notebooks/jupyter/eve_analysis.ipynb":
        return "eve-analysis-notebook-v1"
    if path.name == "run.log":
        return "run-log-v1"
    suffix = path.suffix.lower().lstrip(".") or "binary"
    return f"virosync-{suffix}-v1:{relative}"


def _artifact_identities(
    output_dir: Path,
    paths: list[Path],
) -> tuple:
    identities = []
    seen: set[str] = set()
    for path in sorted(paths, key=lambda item: item.as_posix()):
        try:
            path.lstat()
        except FileNotFoundError:
            continue
        if path.is_symlink():
            raise ValueError(f"run artifact must not be a symlink: {path}")
        if not path.is_file():
            continue
        relative = path.relative_to(output_dir).as_posix()
        if relative in seen:
            continue
        identities.append(
            build_artifact_identity(
                relative,
                root=output_dir,
                schema=_artifact_schema(path, output_dir),
            )
        )
        seen.add(relative)
    return tuple(identities)


def _phase_artifacts(output_dir: Path, phase: int) -> tuple:
    relative_paths = {
        0: (
            "phase0/ablation_events.json",
            "phase0/masking/masking_status.json",
            "phase0/proteome.fasta",
            "phase0/genes.gff",
        ),
        1: (
            "phase1/ablation_events.json",
            "phase1/resume_state.json",
        ),
        2: (
            "phase2/ablation_events.json",
            "phase2/refined_boundaries.bed",
            "phase2/refined_state.json",
            "phase2/resume_state.json",
            "phase2/superset_diamond/full_proteome.tsv",
        ),
        3: ("phase3/ablation_events.json",),
    }[phase]
    paths = [output_dir / relative for relative in relative_paths]
    if phase == 0:
        status_path = output_dir / "phase0" / "masking" / "masking_status.json"
        if status_path.is_file():
            masking_result = load_masking_result(status_path)
            masked_path = Path(masking_result.output_path)
            try:
                relative_masked_path = masked_path.resolve(strict=True).relative_to(
                    output_dir.resolve(strict=True)
                )
            except (OSError, ValueError):
                pass
            else:
                paths.append(output_dir / relative_masked_path)
    return _artifact_identities(
        output_dir,
        [path for path in paths if path.is_file()],
    )


def _write_ablation_events(
    output_dir: Path,
    *,
    ablation_id: AblationID,
    counters: AblationCounters | None = None,
    phase: int | None = None,
) -> Path:
    """Atomically write one canonical cumulative ablation event document."""

    if counters is None:
        counters = AblationCounters.for_ablation(ablation_id)
    events = AblationEvents(ablation_id=ablation_id, counters=counters)
    parent = output_dir if phase is None else output_dir / f"phase{phase}"
    path = parent / "ablation_events.json"
    with atomic_write_context(path, "wb") as handle:
        handle.write(events.to_bytes())
    return path


def _merge_ablation_counts(
    *,
    ablation_id: AblationID,
    current: AblationCounters,
    additional: InterventionCounts | None,
) -> AblationCounters:
    """Add one phase-local count group to the selected cumulative arm."""

    if additional is None:
        return current
    if not isinstance(additional, InterventionCounts):
        raise TypeError("phase ablation_counts must be InterventionCounts")
    return AblationCounters.for_ablation(
        ablation_id,
        opportunities=current.total_opportunities + additional.opportunities,
        interventions=current.total_interventions + additional.interventions,
        changed=current.total_changed + additional.changed,
    )


def _validate_ablation_events_file(
    path: Path,
    *,
    expected_ablation_id: AblationID,
) -> AblationEvents:
    """Validate canonical event bytes and their selected benchmark arm."""

    events = validate_ablation_events_bytes(path.read_bytes())
    if events.ablation_id is not expected_ablation_id:
        raise ValueError(
            "ablation event ID differs from the authenticated run identity"
        )
    return events


def _final_artifacts(
    output_dir: Path,
    *,
    ablation_id: AblationID,
) -> tuple:
    _validate_ablation_events_file(
        output_dir / "ablation_events.json",
        expected_ablation_id=ablation_id,
    )
    paths: list[Path] = []
    for candidate in output_dir.iterdir():
        if candidate.is_symlink():
            if _is_owned_root_final(candidate) or candidate.name == "notebooks":
                raise ValueError(f"final output must not be a symlink: {candidate}")
            continue
        if candidate.is_file() and _is_owned_root_final(candidate):
            paths.append(candidate)
    synthesis = output_dir / "phase3_synthesis"
    paths.extend(
        candidate
        for name in (
            "virosync_predictions.tsv",
            "virosync_predictions_detailed.tsv",
            "virosync_predictions.bed",
            "virosync_predictions.gff3",
            "virosync_summary.json",
            "eve_ani_edges.tsv",
        )
        if (candidate := synthesis / name).is_file()
    )
    notebook = output_dir / "notebooks" / "jupyter" / "eve_analysis.ipynb"
    if notebook.is_file():
        paths.append(notebook)
    return _artifact_identities(output_dir, paths)


def _result_identity(
    summary: dict,
    *,
    terminal_phase: int | None,
    benchmark_eligible: bool,
    promoted_low_rows: int,
) -> dict:
    low_rows = int(summary.get("low_tier", 0) or 0)
    if (
        type(promoted_low_rows) is not int
        or promoted_low_rows < 0
        or promoted_low_rows > low_rows
    ):
        raise ValueError(
            "promoted_low_rows must count a subset of canonical LOW rows"
        )
    return {
        "terminal_phase": terminal_phase,
        "canonical_rows": int(summary.get("accepted", 0) or 0),
        "detailed_rows": int(summary.get("predictions", 0) or 0),
        "accepted_bp": int(summary.get("accepted_bp", 0) or 0),
        "promoted_low_rows": promoted_low_rows,
        "class_counts": {
            eve_class: int(summary.get(count_key, 0) or 0)
            for eve_class, count_key in EFFECTIVE_EVE_CLASS_COUNT_KEYS.items()
        },
        "tier_counts": {
            "HIGH": int(summary.get("high_tier", 0) or 0),
            "MEDIUM": int(summary.get("medium_tier", 0) or 0),
            "LOW": int(summary.get("low_tier", 0) or 0),
        },
        "benchmark_eligible": benchmark_eligible,
    }


def _run_benchmark_eligible(output_dir: Path) -> bool:
    status_path = output_dir / "phase0" / "masking" / "masking_status.json"
    return bool(load_masking_result(status_path).benchmark_eligible)


def _authenticated_output_files(output_dir: Path) -> dict[str, object]:
    """Return the canonical public mapping from a published schema-v3 success."""

    output_dir = Path(output_dir)
    state = load_run_state(output_dir)
    if state.status != "success":
        raise ValueError("output files require a published schema-v3 success state")
    return {
        "run_state": str(output_dir / RUN_STATE_FILENAME),
        "artifacts": {
            artifact.relative_path: str(output_dir / artifact.relative_path)
            for artifact in sorted(
                state.artifacts,
                key=lambda item: item.relative_path,
            )
        },
    }


def _publish_phase_state(
    *,
    output_dir: Path,
    phase: int,
    run_fingerprint: str,
    artifacts: tuple,
    outcome: str = "complete",
    masking_result: MaskingResult | None = None,
    requested_masking: dict[str, object] | None = None,
) -> Path:
    dependency = (
        run_fingerprint
        if phase == 0
        else marker_sha256(output_dir, phase - 1)
    )
    kwargs = {}
    if phase == 0:
        if masking_result is None or requested_masking is None:
            raise ValueError("Phase 0 state requires requested and actual masking")
        status_path = Path(masking_result.status_path)
        actual = json.loads(status_path.read_text(encoding="utf-8"))
        actual["status_sha256"] = masking_result.status_sha256
        kwargs = {
            "requested_masking": requested_masking,
            "actual_masking": actual,
        }
    record = build_phase_record(
        phase=phase,
        run_fingerprint=run_fingerprint,
        dependency_sha256=dependency,
        artifacts=artifacts,
        outcome=outcome,
        **kwargs,
    )
    return publish_phase_completion(output_dir, record)


def _load_verified_phase0(
    *,
    genome_path: Path,
    output_dir: Path,
    masking: MaskingConfig,
) -> dict:
    status_path = output_dir / "phase0" / "masking" / "masking_status.json"
    result = load_masking_result(
        status_path,
        expected_config=masking,
        expected_input=genome_path,
    )
    proteome_path = output_dir / "phase0" / "proteome.fasta"
    genes_path = output_dir / "phase0" / "genes.gff"
    if not proteome_path.is_file() or not genes_path.is_file():
        raise ValueError("verified Phase 0 marker lacks reloadable outputs")
    return {
        "masked_path": result.output_path,
        "repeat_regions": list(result.repeat_regions),
        "proteome_path": proteome_path,
        "n_genes": _count_fasta_records(proteome_path),
        "elapsed": 0.0,
        "masking_result": result,
    }


def _publish_terminal_success(
    *,
    output_dir: Path,
    phase: int,
    run_fingerprint: str,
    result: dict,
    ablation_id: AblationID,
    summary: dict | None = None,
    outcome: str = "terminal_zero",
) -> dict:
    current_events = _validate_ablation_events_file(
        output_dir / "ablation_events.json",
        expected_ablation_id=ablation_id,
    )
    _write_ablation_events(
        output_dir,
        ablation_id=ablation_id,
        counters=current_events.counters,
        phase=phase,
    )
    if summary is None:
        summary = _empty_prediction_summary()
    benchmark_eligible = _run_benchmark_eligible(output_dir)
    final_artifacts = _final_artifacts(
        output_dir,
        ablation_id=ablation_id,
    )
    phase_artifacts = tuple(
        sorted(
            (*_phase_artifacts(output_dir, phase), *final_artifacts),
            key=lambda artifact: artifact.relative_path,
        )
    )
    _publish_phase_state(
        output_dir=output_dir,
        phase=phase,
        run_fingerprint=run_fingerprint,
        artifacts=phase_artifacts,
        outcome=outcome,
    )
    state = publish_run_success(
        output_dir,
        run_fingerprint=run_fingerprint,
        artifacts=final_artifacts,
        result=_result_identity(
            summary,
            terminal_phase=phase,
            benchmark_eligible=benchmark_eligible,
            promoted_low_rows=0,
        ),
    )
    return {
        **result,
        "benchmark_eligible": benchmark_eligible,
        "legacy_resume": False,
        "output_files": _authenticated_output_files(output_dir),
        "run_state": asdict(state),
    }


def _write_combined_eve_fasta(
    *,
    output_dir: Path,
    genome_id: str,
    genome_path: Path,
    results: list | tuple,
) -> Path | None:
    """Write the root combined EVE FASTA shared by terminal output paths."""

    if not genome_path.is_file() or not results:
        return None

    from Bio import SeqIO
    from virosync.pipeline.phase3.output_generator import OutputGenerator

    genome_seqs = {
        str(record.id): str(record.seq)
        for record in SeqIO.parse(genome_path, "fasta")
    }
    eve_fasta = output_dir / f"{genome_id}_eves.fna"
    generator = OutputGenerator(output_dir, genome_fasta=genome_path)
    generator.genome_sequences = genome_seqs
    generator.write_combined_eve_fasta(list(results), eve_fasta)
    if not eve_fasta.is_file():
        raise RuntimeError(
            "combined EVE FASTA writer produced no file for nonempty results"
        )
    return eve_fasta


def _publish_a1_seed_surface(
    *,
    output_dir: Path,
    genome_id: str,
    run_fingerprint: str,
    merged_seeds: list,
    masked_path: Path,
    proteome_path: Path,
    taxonomy_labels_file: Path | None,
    seed_marker_allowlist: list[str] | None,
    extended_output: bool,
    export_all_eve_sequences: bool,
    logger,
    current_counters: AblationCounters,
) -> dict:
    """Publish A1's nonzero Phase-1 seed surface without Phase 2 or 3."""

    from virosync.pipeline.phase3.phase1_surface import (
        build_phase1_seed_surface,
    )

    surface = build_phase1_seed_surface(merged_seeds)
    counters = _merge_ablation_counts(
        ablation_id=AblationID.A1,
        current=current_counters,
        additional=surface.intervention_counts,
    )
    _write_ablation_events(
        output_dir,
        ablation_id=AblationID.A1,
        counters=counters,
    )
    output_files = call_task(
        generate_outputs_task,
        verification_results=list(surface.detailed_results),
        canonical_results=list(surface.canonical_results),
        output_dir=output_dir,
        genome_path=masked_path,
        proteome_path=proteome_path,
        accepted_only=False,
        extended_output=extended_output,
        seed_marker_allowlist=seed_marker_allowlist,
        export_all_eve_sequences=export_all_eve_sequences,
        promoted_low_results=[],
    )
    canonical_path = Path(output_files["predictions_tsv"])
    detailed_path = Path(output_files["predictions_detailed_tsv"])
    invariant_path = output_dir / "virosync_tsv_invariant_report.tsv"
    enforce_tsv_invariants(
        detailed_tsv=detailed_path,
        report_out=invariant_path,
    )
    output_files["tsv_invariant_report"] = invariant_path
    eve_fasta = _write_combined_eve_fasta(
        output_dir=output_dir,
        genome_id=genome_id,
        genome_path=masked_path,
        results=surface.canonical_results,
    )
    if eve_fasta is not None:
        output_files["eve_fasta"] = eve_fasta
    output_files.update(
        _generate_required_reports(
            output_dir=output_dir,
            genome_id=genome_id,
            taxonomy_labels_file=taxonomy_labels_file,
            logger=logger,
        )
    )
    persisted_summary = _summarize_prediction_outputs(
        canonical_path,
        detailed_path,
    )
    run_log_path = output_dir / "run.log"
    with atomic_write_context(run_log_path, "w") as handle:
        handle.write(f"# ViroSync Run Log: {genome_id}\n\n")
        handle.write("## Results Summary\n")
        handle.write("Prediction stage: Phase-1 seed surface (A1)\n")
        handle.write(
            f"Canonical EVEs: {persisted_summary['accepted']}\n"
        )
        handle.write(
            f"Detailed candidates: {persisted_summary['predictions']}\n"
        )
        handle.write("Confidence kind: not_scored\n")
    output_files["run_log"] = run_log_path
    completion = _write_completion_manifest(
        output_dir=output_dir,
        genome_id=genome_id,
        status="success",
        output_files=output_files,
        fingerprint=run_fingerprint,
    )
    output_files["completion_manifest"] = completion
    return _publish_terminal_success(
        output_dir=output_dir,
        phase=1,
        run_fingerprint=run_fingerprint,
        result={
            "genome_id": genome_id,
            "success": True,
            "output_files": output_files,
            **persisted_summary,
        },
        ablation_id=AblationID.A1,
        summary=persisted_summary,
        outcome="terminal_ablation",
    )


def _mark_failed_attempt(output_dir: Path, exc: BaseException) -> None:
    try:
        state = load_run_state(output_dir)
    except (FileNotFoundError, OSError, TypeError, ValueError):
        return
    if state.status != "running":
        return
    publish_run_failed(
        output_dir,
        run_fingerprint=state.run_fingerprint,
        error_type=type(exc).__name__,
        message=str(exc),
    )


def _serialized_run(function):
    """Hold one sibling lock and persist uncaught/returned failures."""

    @wraps(function)
    def _wrapped(*args, **kwargs):
        bound = inspect.signature(function).bind(*args, **kwargs)
        output_dir = _validate_clean_run_target(
            Path(bound.arguments["output_dir"]),
            bound.arguments["genome_id"],
        )
        with sibling_run_lock(output_dir):
            try:
                result = function(*args, **kwargs)
            except BaseException as exc:
                _mark_failed_attempt(output_dir, exc)
                raise
            if isinstance(result, dict) and result.get("success") is False:
                _mark_failed_attempt(
                    output_dir,
                    RuntimeError(str(result.get("error") or "pipeline failed")),
                )
            return result

    return _wrapped


def _validate_clean_run_target(output_dir: Path, genome_id: str) -> Path:
    """Return the lexical per-genome output path after resolved-path validation."""
    validate_path_component(genome_id, "genome ID")
    output_dir = Path(output_dir)
    if ".." in output_dir.parts:
        raise ValueError(f"refusing output path with a parent segment: {output_dir}")
    if output_dir.is_symlink():
        raise ValueError(f"refusing to use symlink output directory: {output_dir}")
    if output_dir.name != genome_id:
        raise ValueError(
            "output directory final component must exactly match genome ID: "
            f"{output_dir.name!r} != {genome_id!r}"
        )
    require_strict_child(output_dir.parent, output_dir)
    return output_dir


def _remove_output_dir(output_dir: Path, genome_id: str) -> None:
    """Remove a pre-existing per-genome output after an immediate safety check."""
    output_dir = _validate_clean_run_target(output_dir, genome_id)
    shutil.rmtree(require_strict_child(output_dir.parent, output_dir))


def _revalidate_completed_run(
    output_dir: Path,
    completed_artifacts: dict[str, Path],
):
    """Recheck cached detailed output before returning a resumed success."""
    output_dir = Path(output_dir)
    invariant_report_path = output_dir / "virosync_tsv_invariant_report.tsv"
    gene_taxonomy = (
        output_dir / "phase3_synthesis" / "gene_taxonomy" / "gene_taxonomy_all.tsv"
    )
    try:
        report = enforce_tsv_invariants(
            detailed_tsv=completed_artifacts["predictions_detailed"],
            report_out=invariant_report_path,
            gene_taxonomy_all_tsv=gene_taxonomy if gene_taxonomy.exists() else None,
        )
    except Exception:
        _clear_success_markers(output_dir)
        raise
    completed_artifacts["tsv_invariant_report"] = invariant_report_path
    return report


def _pin_cuda_device(device: str, logger) -> None:
    """Pin torch to requested GPU in worker process if one is configured."""
    if (device or "").lower() != "cuda":
        return
    try:
        import os
        import torch
    except Exception:
        return
    if not torch.cuda.is_available():
        return
    visible = (os.environ.get("CUDA_VISIBLE_DEVICES") or "").strip()
    requested = (os.environ.get("VIROSYNC_GPU") or "").strip()
    target_idx: Optional[int] = None
    if visible:
        first = visible.split(",", 1)[0].strip()
        if torch.cuda.device_count() == 1:
            target_idx = 0
        elif first.isdigit() and int(first) < torch.cuda.device_count():
            target_idx = int(first)
        elif requested.isdigit() and int(requested) < torch.cuda.device_count():
            # CUDA UUID selectors can leave multiple devices visible on some setups;
            # fall back to the explicit requested index when available.
            target_idx = int(requested)
    elif requested.isdigit() and int(requested) < torch.cuda.device_count():
        target_idx = int(requested)
    if target_idx is None:
        return
    try:
        torch.cuda.set_device(target_idx)
        logger.info("Pinned CUDA device: cuda:%d", target_idx)
    except Exception as exc:
        logger.warning("Failed to pin CUDA device: %s", exc)


def _run_phase0_subflow(
    genome_path: Path,
    output_dir: Path,
    genome_id: str,
    threads: int,
    skip_masking: Optional[bool],
    resume: bool,
    logger,
    masking: Optional[MaskingConfig] = None,
) -> dict:
    """
    Phase 0: Preprocessing (repeat masking + gene prediction).

    Args:
        genome_path: Path to input genome FASTA
        output_dir: Base output directory
        genome_id: Genome identifier for logging
        threads: Number of threads for parallel processing
        skip_masking: Deprecated compatibility override
        logger: Logger instance

    Returns:
        dict with keys:
            - masked_path: Path to masked genome FASTA
            - repeat_regions: List of RepeatRegion objects
            - proteome_path: Path to protein FASTA
            - n_genes: Number of predicted genes
    """
    import time

    phase0_start = time.time()
    logger.info("-" * 60)
    logger.info("Phase 0: Preprocessing (masking + proteome)")

    phase0_dir = output_dir / "phase0"
    cached_proteome = phase0_dir / "proteome.fasta"
    cached_genes = phase0_dir / "genes.gff"
    if resume and cached_proteome.exists() and cached_genes.exists():
        logger.info(
            "Phase 0 partial-cache reuse is disabled until schema-v3 run state; "
            "recomputing masking and proteome"
        )

    masking = masking or MaskingConfig()
    if skip_masking is not None:
        if skip_masking:
            masking = masking.with_backend(MaskingBackend.OFF)
        elif masking.backend is MaskingBackend.OFF:
            masking = masking.with_backend(MaskingBackend.TRF_REPEATMASKER)

    masking_task_result = call_task(
        mask_genome_task,
        genome_path=genome_path,
        output_dir=output_dir / "phase0",
        threads=threads,
        masking=masking,
    )
    status_path = phase0_dir / "masking" / "masking_status.json"
    if isinstance(masking_task_result, MaskingResult):
        try:
            if (
                masking_task_result.status_path is None
                or masking_task_result.status_path.resolve() != status_path.resolve()
            ):
                raise ValueError(
                    "returned status path is not phase0/masking/masking_status.json"
                )
            masking_result = load_masking_result(
                status_path,
                repeat_regions=masking_task_result.repeat_regions,
                expected_config=masking,
                expected_input=genome_path,
            )
            if (
                masking_result.to_status_payload()
                != masking_task_result.to_status_payload()
                or masking_task_result.status_sha256
                != masking_result.status_sha256
            ):
                raise ValueError("returned result disagrees with persisted status")
        except (OSError, ValueError, KeyError) as exc:
            raise ValueError(f"masking status mismatch: {exc}") from exc
    else:
        try:
            masked_path, repeat_regions = masking_task_result
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "masking task returned neither MaskingResult nor a legacy path tuple"
            ) from exc
        try:
            masking_result = load_masking_result(
                status_path,
                repeat_regions=tuple(repeat_regions),
                expected_config=masking,
                expected_input=genome_path,
            )
        except (OSError, ValueError, KeyError) as exc:
            raise ValueError(f"masking status mismatch: {exc}") from exc
        if masking_result.output_path.resolve() != Path(masked_path).resolve():
            raise ValueError(
                "masking status output path mismatch with the sequence passed to Prodigal"
            )

    masked_path = masking_result.output_path
    repeat_regions = list(masking_result.repeat_regions)

    proteome_result = call_task(
        generate_proteome_task,
        genome_path=masked_path,
        output_dir=output_dir / "phase0",
        threads=threads,
    )
    proteome_path, n_genes = proteome_result

    logger.info(f"Generated {n_genes} genes (prodigal-gv)")

    phase0_elapsed = time.time() - phase0_start
    logger.info(f"Phase 0 complete: {phase0_elapsed:.1f}s")

    return {
        "masked_path": masked_path,
        "repeat_regions": repeat_regions,
        "proteome_path": proteome_path,
        "n_genes": n_genes,
        "elapsed": phase0_elapsed,
        "masking_result": masking_result,
    }


@_serialized_run
def _single_genome_flow_impl(
    genome_path: Path,
    output_dir: Path,
    genome_id: str,
    ablation_id: str,
    ablation_contract_sha256: str,
    # Configuration object (recommended - use instead of individual kwargs)
    # Database paths (can override config)
    hmm_database: Optional[Path],
    hmm_allowlist: Optional[Path],
    seed_marker_allowlist: Optional[list[str]],
    marker_faa_db: Optional[Path],
    marker_db: Optional[Path],
    gene_taxonomy_faa_db: Optional[Path],
    marker_faa_dir: Optional[Path],
    faa_dir: Optional[Path],
    gvclass_db: Optional[Path],
    diamond_db: Optional[Path],
    enable_phylogenetic: bool,
    # Taxonomy lookup for host signature comparison
    taxonomy_labels_file: Optional[Path],
    # Host taxonomy configuration
    host_prefixes: Optional[list[str]],
    host_label: str,
    high_pident_host_threshold: float,
    # Parameters
    threads: int,
    max_threads: Optional[int],
    device: str,
    masking: MaskingConfig,
    skip_structural: bool,
    use_boltz: bool,
    boltz_mcp_only: bool,
    boltz_use_msa_server: bool,
    boltz_min_seq_len: int,
    boltz_max_seq_len: int,
    boltz_no_kernels: bool,
    use_tmvec_database: bool,
    tmvec_require_gpu: bool,
    tmvec_databases: Optional[list[str]],
    tmvec_database_dir: Optional[Path],
    tmvec_min_score: float,
    viral_structure_db: Optional[Path],
    assembly_mode: str,
    high_tier_threshold: float,
    low_tier_threshold: float,
    use_crf_in_final_score: bool,
    priority_marker_list: Optional[list[str]],
    marker_floor_priority_only: float,
    marker_floor_priority_plus_family: float,
    marker_floor_priority_multi_family: float,
    marker_family_bonus_per_family: float,
    marker_multi_family_bonus: float,
    hmm_chunk_size: Optional[int],
    gene_taxonomy_threads: Optional[int],
    interproscan_enabled: bool,
    interproscan_dir: Optional[Path],
    interproscan_keywords: Optional[list[str]],
    interproscan_threads: Optional[int],
    interproscan_applications: Optional[list[str]],
    extended_output: bool,
    export_all_eve_sequences: bool,
    # GVClass batch classification
    run_gvclass: bool,
    gvclass_path: Optional[Path],
    # HMM-gated workflow options
    rebuild_db: bool,
    initial_window_bp: int,
    initial_window_genes: int,
    min_markers_initial: int,
    extension_kb: int,
    merge_distance: int,
    host_taxonomy_deviation_enabled: bool,
    host_taxonomy_deviation_allow_seeds: bool,
    host_taxonomy_deviation_min_token_len: int,
    host_taxonomy_deviation_min_tokens: int,
    host_taxonomy_deviation_overlap_threshold: float,
    host_taxonomy_deviation_max_pident: float,
    host_taxonomy_deviation_max_hits: int,
    host_taxonomy_deviation_window_bp: int,
    host_taxonomy_deviation_window_count: int,
    host_taxonomy_deviation_window_seed: int,
    host_taxonomy_deviation_window_min_markers: int,
    host_taxonomy_deviation_seed_window_bp: int,
    host_taxonomy_deviation_seed_min_markers: int,
    marker_validation_top_k: int,
    novel_marker_min_score: float,
    novel_marker_min_coverage: float,
    novel_marker_require_cluster: bool,
    boundary_host_trim_enabled: bool,
    boundary_host_trim_window_bp: int,
    boundary_host_trim_step_bp: int,
    boundary_host_trim_max_host_fraction: float,
    boundary_host_trim_min_viral_fraction: float,
    boundary_host_trim_score_threshold: float,
    host_signature_evidence_threshold: float,
    boundary_host_trim_buffer_kb: int,
    boundary_host_signature_min_token_len: int,
    boundary_host_trim_min_overlap_score: float,
    boundary_taxonomy_ml_enabled: bool,
    boundary_taxonomy_ml_model: str,
    boundary_taxonomy_ml_threshold: float,
    boundary_taxonomy_ml_neighbor_window: int,
    taxonomy_weight_mode: str,
    # Batched Diamond configuration (Phase 2)
    # Note: Batched Diamond is always enabled (use_batched_boundary_diamond removed Jan 2026)
    boundary_diamond_flank_genes: int,
    boundary_diamond_control_sample_size: int,
    boundary_diamond_control_min_distance: int,
    boundary_diamond_top_k: int,
    boundary_diamond_chunk_size: int,
    boundary_diamond_random_seed: int,
    boundary_diamond_superset_prototype_enabled: bool,
    resume: bool,
    # Search backend (diamond -- sole backend)
    search_backend: str = "diamond",
    progress_callback: Optional[Callable[[float, str, bool], None]] = None,
) -> dict:
    """
    Implementation of single genome flow (internal).

    Called by the public wrapper only.
    All parameters passed explicitly (no defaults, no config merging).
    """
    import time
    from datetime import datetime, timezone

    selected_ablation = AblationID(ablation_id)
    ablation_id = selected_ablation.value
    if ablation_contract_sha256 != ABLATION_CONTRACT_SHA256:
        raise ValueError("ablation contract SHA-256 differs from this ViroSync build")
    if host_prefixes is None:
        host_prefixes = ["EUK__", "MITO__", "PLASTID__"]
    host_label = (host_label or "EUK").upper()

    logger = get_orchestration_logger(__name__)

    def report_progress(percent: float, stage: str) -> None:
        if progress_callback is not None:
            progress_callback(percent, stage, False)

    report_progress(0, "starting")
    _pin_cuda_device(device=device, logger=logger)
    genome_start_time = time.time()
    logger.info("=" * 60)
    logger.info(f"GENOME: {genome_id}")
    logger.info("=" * 60)
    logger.info(f"Input: {genome_path}")
    if (gvclass_db or diamond_db) and not enable_phylogenetic:
        logger.info("Phylogenetic validation disabled; GVClass/Diamond inputs ignored.")

    # Ensure paths are Path objects
    # Resolve an input symlink once so identity hashing and every phase consume
    # the same pinned regular file even if the link is retargeted mid-run.
    genome_path = Path(genome_path).resolve(strict=True)
    output_dir = _validate_clean_run_target(Path(output_dir), genome_id)
    # A clean run (resume=False, e.g. --clean-run) must start from a pristine
    # output directory: file-existence-based resume checks would otherwise let
    # stale artifacts from an aborted prior run be mistaken for fresh results.
    if not resume and output_dir.exists():
        logger.info("Clean run: removing existing output directory %s", output_dir)
        _remove_output_dir(output_dir, genome_id)
    output_dir.mkdir(parents=True, exist_ok=True)

    if max_threads:
        threads = min(threads, max_threads)

    background = None
    gene_data = None

    flat_identity_config = {
        key: value
        for key, value in locals().items()
        if key in _FINGERPRINT_INPUT_FIELDS
    }
    run_identity, run_fingerprint = _build_run_identity(
        genome_path=genome_path,
        genome_id=genome_id,
        output_dir=output_dir,
        flat_config=flat_identity_config,
        masking=masking,
        device=device,
    )
    # Retain the compatibility field name in schema-v2 reports, but bind it to
    # the complete schema-v3 identity rather than the old 16-character config hash.
    config_fingerprint = run_fingerprint

    completed_artifacts = None
    if resume:
        completed_artifacts = _completed_run_artifacts(
            output_dir,
            expected_fingerprint=run_fingerprint,
            expected_input=genome_path,
        )
    phase3_predictions = output_dir / "phase3_synthesis" / "virosync_predictions.tsv"
    refined_bed = output_dir / "phase2" / "refined_boundaries.bed"
    validated_hits_tsv = output_dir / "phase1" / "marker_validation" / "validated_marker_hits.tsv"

    if resume and completed_artifacts is not None:
        resume_counts = _empty_prediction_summary()
        try:
            resume_counts = _summarize_prediction_outputs(
                completed_artifacts["phase3_predictions"],
                completed_artifacts["predictions_detailed"],
            )
        except Exception:
            _clear_success_markers(output_dir)
            logger.exception(
                "Failed to summarize completed predictions from %s",
                completed_artifacts["phase3_predictions"],
            )
            raise
        state = load_run_state(output_dir)
        benchmark_eligible = _run_benchmark_eligible(output_dir)
        expected_result = _result_identity(
            resume_counts,
            terminal_phase=state.result.get("terminal_phase") if state.result else None,
            benchmark_eligible=benchmark_eligible,
            promoted_low_rows=(
                int(state.result.get("promoted_low_rows", resume_counts["low_tier"]))
                if state.result
                else resume_counts["low_tier"]
            ),
        )
        observed_result = dict(state.result or {})
        observed_result.setdefault(
            "promoted_low_rows",
            resume_counts["low_tier"],
        )
        observed_result = {
            key: observed_result.get(key) for key in expected_result
        }
        if observed_result != expected_result:
            logger.warning(
                "Resume: persisted result counts disagree with authenticated "
                "outputs; invalidating the run."
            )
            invalidate_from_phase(output_dir, from_phase=0)
            completed_artifacts = None
        if completed_artifacts is not None:
            logger.info("Resume: completed schema-3 outputs validated")
            report_progress(100, "complete")
            return {
                "genome_id": genome_id,
                "success": True,
                "benchmark_eligible": benchmark_eligible,
                "legacy_resume": False,
                "output_files": _authenticated_output_files(output_dir),
                "elapsed_sec": 0.0,
                **resume_counts,
            }

    resume_plan = (
        plan_resume(output_dir, expected_run_fingerprint=run_fingerprint)
        if resume
        else None
    )
    reusable_phases = set(resume_plan.reusable_phases if resume_plan else ())
    restart_phase = resume_plan.restart_phase if resume_plan else 0
    recoverable_completion = resume_plan is not None and (
        resume_plan.terminal_phase is not None
        or restart_phase == len(PHASE_MARKER_FILENAMES)
    )
    if not recoverable_completion and restart_phase < len(PHASE_MARKER_FILENAMES):
        invalidate_from_phase(output_dir, from_phase=restart_phase)
    ablation_counters = AblationCounters.for_ablation(selected_ablation)
    if reusable_phases:
        previous_events = _validate_ablation_events_file(
            output_dir / f"phase{max(reusable_phases)}" / "ablation_events.json",
            expected_ablation_id=selected_ablation,
        )
        ablation_counters = previous_events.counters
    publish_run_started(
        output_dir,
        run_fingerprint=run_fingerprint,
        identities=run_identity,
        preserve_success_artifacts=recoverable_completion,
    )
    if not recoverable_completion:
        _write_ablation_events(
            output_dir,
            ablation_id=selected_ablation,
            counters=ablation_counters,
        )

    # A crash can leave the final/terminal phase marker durable while the last
    # success-state replace never happened. Revalidate and promote without rerunning.
    if recoverable_completion:
        final_artifacts = _final_artifacts(
            output_dir,
            ablation_id=selected_ablation,
        )
        final_paths = {
            artifact.relative_path: output_dir / artifact.relative_path
            for artifact in final_artifacts
        }

        def _recovered_path(*relative_paths: str) -> Path:
            for relative_path in relative_paths:
                if relative_path in final_paths:
                    return final_paths[relative_path]
            raise ValueError(
                f"recovered completion is missing one of {relative_paths!r}"
            )

        recovered_canonical = _recovered_path(
            "phase3_synthesis/virosync_predictions.tsv",
            "virosync_predictions.tsv",
        )
        recovered_detailed = _recovered_path(
            "virosync_predictions_detailed.tsv",
            "phase3_synthesis/virosync_predictions_detailed.tsv",
        )
        recovered_summary = _summarize_prediction_outputs(
            recovered_canonical,
            recovered_detailed,
        )
        recovered_summary_document = json.loads(
            _recovered_path(
                "phase3_synthesis/virosync_summary.json",
                "virosync_summary.json",
            ).read_text(encoding="utf-8")
        )
        recovered_statistics = recovered_summary_document.get("statistics")
        if not isinstance(recovered_statistics, dict):
            raise ValueError("recovered ViroSync summary has no statistics object")
        recovered_promoted_low = recovered_statistics.get(
            "promoted_low_confidence"
        )
        if type(recovered_promoted_low) is not int:
            raise ValueError(
                "recovered ViroSync summary has no promoted LOW count"
            )
        benchmark_eligible = _run_benchmark_eligible(output_dir)
        publish_run_success(
            output_dir,
            run_fingerprint=run_fingerprint,
            artifacts=final_artifacts,
            result=_result_identity(
                recovered_summary,
                terminal_phase=resume_plan.terminal_phase,
                benchmark_eligible=benchmark_eligible,
                promoted_low_rows=recovered_promoted_low,
            ),
        )
        report_progress(100, "complete")
        return {
            "genome_id": genome_id,
            "success": True,
            "benchmark_eligible": benchmark_eligible,
            "legacy_resume": False,
            "output_files": _authenticated_output_files(output_dir),
            "elapsed_sec": 0.0,
            **recovered_summary,
        }

    # === PHASE 0: Preprocessing ===
    if 0 in reusable_phases:
        phase0_result = _load_verified_phase0(
            genome_path=genome_path,
            output_dir=output_dir,
            masking=masking,
        )
        logger.info("Resume: reused authenticated Phase 0")
    else:
        phase0_result = _run_phase0_subflow(
            genome_path=genome_path,
            output_dir=output_dir,
            genome_id=genome_id,
            threads=threads,
            skip_masking=None,
            resume=False,
            logger=logger,
            masking=masking,
        )

    # Extract Phase 0 outputs
    masked_path = phase0_result["masked_path"]
    repeat_regions = phase0_result["repeat_regions"]
    proteome_path = phase0_result["proteome_path"]
    n_genes = phase0_result["n_genes"]
    phase0_elapsed = phase0_result["elapsed"]
    masking_result = phase0_result["masking_result"]
    if 0 not in reusable_phases:
        _write_ablation_events(
            output_dir,
            ablation_id=selected_ablation,
            counters=ablation_counters,
            phase=0,
        )
        _publish_phase_state(
            output_dir=output_dir,
            phase=0,
            run_fingerprint=run_fingerprint,
            artifacts=_phase_artifacts(output_dir, 0),
            masking_result=masking_result,
            requested_masking=run_identity["requested_masking"],
        )
    report_progress(20, "phase 0 complete")

    # === PHASE 1: Seeding ===
    phase1_result = _run_phase1_subflow(
        # Core inputs from Phase 0
        masked_path=masked_path,
        proteome_path=proteome_path,
        repeat_regions=repeat_regions,
        # Core identifiers
        output_dir=output_dir,
        genome_id=genome_id,
        # HMM & Database parameters
        hmm_database=hmm_database,
        hmm_allowlist=hmm_allowlist,
        hmm_chunk_size=hmm_chunk_size,
        marker_faa_db=marker_faa_db,
        marker_faa_dir=marker_faa_dir,
        marker_db=marker_db,
        faa_dir=faa_dir,
        gene_taxonomy_faa_db=gene_taxonomy_faa_db,
        # Taxonomy parameters
        taxonomy_labels_file=taxonomy_labels_file,
        host_prefixes=host_prefixes,
        host_label=host_label,
        taxonomy_weight_mode=taxonomy_weight_mode,
        # Host taxonomy deviation parameters
        host_taxonomy_deviation_enabled=host_taxonomy_deviation_enabled,
        host_taxonomy_deviation_allow_seeds=host_taxonomy_deviation_allow_seeds,
        host_taxonomy_deviation_min_token_len=host_taxonomy_deviation_min_token_len,
        host_taxonomy_deviation_min_tokens=host_taxonomy_deviation_min_tokens,
        host_taxonomy_deviation_overlap_threshold=host_taxonomy_deviation_overlap_threshold,
        host_taxonomy_deviation_max_pident=host_taxonomy_deviation_max_pident,
        host_taxonomy_deviation_max_hits=host_taxonomy_deviation_max_hits,
        host_taxonomy_deviation_window_bp=host_taxonomy_deviation_window_bp,
        host_taxonomy_deviation_window_count=host_taxonomy_deviation_window_count,
        host_taxonomy_deviation_window_seed=host_taxonomy_deviation_window_seed,
        host_taxonomy_deviation_window_min_markers=host_taxonomy_deviation_window_min_markers,
        host_taxonomy_deviation_seed_window_bp=host_taxonomy_deviation_seed_window_bp,
        host_taxonomy_deviation_seed_min_markers=host_taxonomy_deviation_seed_min_markers,
        marker_validation_top_k=marker_validation_top_k,
        novel_marker_min_score=novel_marker_min_score,
        novel_marker_min_coverage=novel_marker_min_coverage,
        novel_marker_require_cluster=novel_marker_require_cluster,
        # Region assembly parameters
        initial_window_bp=initial_window_bp,
        initial_window_genes=initial_window_genes,
        min_markers_initial=min_markers_initial,
        extension_kb=extension_kb,
        merge_distance=merge_distance,
        boundary_host_signature_min_token_len=boundary_host_signature_min_token_len,
        # Workflow configuration
        rebuild_db=rebuild_db,
        assembly_mode=assembly_mode,
        extended_output=extended_output,
        resume=1 in reusable_phases,
        # Threading
        threads=threads,
        # Search backend
        search_backend=search_backend,
        # Logger
        logger=logger,
        # Resume fingerprint for early-exit completion manifests
        config_fingerprint=config_fingerprint,
        resume_authorized=1 in reusable_phases,
        ablation_id=selected_ablation,
    )
    report_progress(45, "phase 1 complete")

    # Handle Phase 1 early exit (error or no predictions)
    if "success" in phase1_result:
        if phase1_result.get("success") is True:
            ablation_counters = _merge_ablation_counts(
                ablation_id=selected_ablation,
                current=ablation_counters,
                additional=phase1_result.get("ablation_counts"),
            )
            _write_ablation_events(
                output_dir,
                ablation_id=selected_ablation,
                counters=ablation_counters,
            )
            return _publish_terminal_success(
                output_dir=output_dir,
                phase=1,
                run_fingerprint=run_fingerprint,
                result=phase1_result,
                ablation_id=selected_ablation,
            )
        return phase1_result

    # Extract Phase 1 outputs
    merged_seeds = phase1_result["merged_seeds"]
    validated_markers = phase1_result["validated_markers"]
    host_signature_model = phase1_result["host_signature_model"]
    host_signatures = phase1_result["host_signatures"]
    background = phase1_result["background"]
    gene_data = phase1_result["gene_data"]
    host_deviation_summary = phase1_result["host_deviation_summary"]
    phase1_elapsed = phase1_result["elapsed"]
    host_signature_model_payload = host_signature_model.to_dict() if host_signature_model else None
    if 1 not in reusable_phases:
        ablation_counters = _merge_ablation_counts(
            ablation_id=selected_ablation,
            current=ablation_counters,
            additional=phase1_result.get("ablation_counts"),
        )
        _write_ablation_events(
            output_dir,
            ablation_id=selected_ablation,
            counters=ablation_counters,
        )
    if selected_ablation is AblationID.A1:
        return _publish_a1_seed_surface(
            output_dir=output_dir,
            genome_id=genome_id,
            run_fingerprint=run_fingerprint,
            merged_seeds=merged_seeds,
            masked_path=masked_path,
            proteome_path=proteome_path,
            taxonomy_labels_file=taxonomy_labels_file,
            seed_marker_allowlist=seed_marker_allowlist,
            extended_output=extended_output,
            export_all_eve_sequences=export_all_eve_sequences,
            logger=logger,
            current_counters=ablation_counters,
        )
    if 1 not in reusable_phases:
        _write_ablation_events(
            output_dir,
            ablation_id=selected_ablation,
            counters=ablation_counters,
            phase=1,
        )
        _publish_phase_state(
            output_dir=output_dir,
            phase=1,
            run_fingerprint=run_fingerprint,
            artifacts=_phase_artifacts(output_dir, 1),
        )

    # === PHASE 2: Boundary Refinement ===
    from virosync.orchestration.resource_monitor import ResourceMonitor

    with ResourceMonitor(
        task_name="boundary_refinement",
        genome_id=genome_id,
        phase="phase2",
        output_dir=Path(output_dir) / "phase2",
        threads=threads,
    ):
        phase2_result = _run_phase2_subflow(
            # Core inputs from Phase 0
            masked_path=masked_path,
            proteome_path=proteome_path,
            # Core inputs from Phase 1
            merged_seeds=merged_seeds,
            validated_markers=validated_markers,
            host_signature_model=host_signature_model,
            # Core identifiers
            output_dir=output_dir,
            genome_id=genome_id,
            # Resume configuration
            resume=2 in reusable_phases,
            refined_bed=refined_bed,
            # Database parameters
            gene_taxonomy_faa_db=gene_taxonomy_faa_db,
            marker_db=marker_db,
            taxonomy_labels_file=taxonomy_labels_file,
            # Host configuration
            host_prefixes=host_prefixes,
            host_label=host_label,
            high_pident_host_threshold=high_pident_host_threshold,
            # Phase 2a: Host-signature trimming parameters
            boundary_host_trim_enabled=boundary_host_trim_enabled,
            boundary_host_trim_window_bp=boundary_host_trim_window_bp,
            boundary_host_trim_step_bp=boundary_host_trim_step_bp,
            boundary_host_trim_max_host_fraction=boundary_host_trim_max_host_fraction,
            boundary_host_trim_min_viral_fraction=boundary_host_trim_min_viral_fraction,
            boundary_host_trim_score_threshold=boundary_host_trim_score_threshold,
            boundary_host_trim_buffer_kb=boundary_host_trim_buffer_kb,
            boundary_host_trim_min_overlap_score=boundary_host_trim_min_overlap_score,
            boundary_host_signature_min_token_len=boundary_host_signature_min_token_len,
            taxonomy_weight_mode=taxonomy_weight_mode,
            boundary_taxonomy_ml_enabled=boundary_taxonomy_ml_enabled,
            boundary_taxonomy_ml_model=boundary_taxonomy_ml_model,
            boundary_taxonomy_ml_threshold=boundary_taxonomy_ml_threshold,
            boundary_taxonomy_ml_neighbor_window=boundary_taxonomy_ml_neighbor_window,
            # Phase 2b: Batched Diamond parameters
            boundary_diamond_flank_genes=boundary_diamond_flank_genes,
            boundary_diamond_control_sample_size=boundary_diamond_control_sample_size,
            boundary_diamond_control_min_distance=boundary_diamond_control_min_distance,
            boundary_diamond_top_k=boundary_diamond_top_k,
            boundary_diamond_chunk_size=boundary_diamond_chunk_size,
            boundary_diamond_random_seed=boundary_diamond_random_seed,
            boundary_diamond_superset_prototype_enabled=(
                boundary_diamond_superset_prototype_enabled
            ),
            # Threading
            threads=threads,
            gene_taxonomy_threads=gene_taxonomy_threads,
            # Output configuration
            extended_output=extended_output,
            # Search backend
            search_backend=search_backend,
            # Timing reference
            genome_start_time=genome_start_time,
            # Logger
            logger=logger,
            # Resume fingerprint for early-exit completion manifests
            config_fingerprint=config_fingerprint,
            resume_authorized=2 in reusable_phases,
            ablation_id=selected_ablation,
        )
    report_progress(75, "phase 2 complete")

    # Handle Phase 2 early exit (no seeds or no boundaries)
    if "success" in phase2_result:
        if phase2_result.get("success") is True:
            ablation_counters = _merge_ablation_counts(
                ablation_id=selected_ablation,
                current=ablation_counters,
                additional=phase2_result.get("ablation_counts"),
            )
            _write_ablation_events(
                output_dir,
                ablation_id=selected_ablation,
                counters=ablation_counters,
            )
            return _publish_terminal_success(
                output_dir=output_dir,
                phase=2,
                run_fingerprint=run_fingerprint,
                result=phase2_result,
                ablation_id=selected_ablation,
            )
        return phase2_result

    # Extract Phase 2 outputs
    refined_boundaries = phase2_result["refined_boundaries"]
    boundary_taxonomy_map = phase2_result["boundary_taxonomy_map"]
    boundary_control_stats = phase2_result["boundary_control_stats"]
    boundary_diamond_query = phase2_result["boundary_diamond_query"]
    proteome_index = phase2_result["proteome_index"]
    goto_phase3 = phase2_result["goto_phase3"]
    boundaries_bed = phase2_result["boundaries_bed"]
    phase2_elapsed = phase2_result["elapsed"]
    if 2 not in reusable_phases:
        ablation_counters = _merge_ablation_counts(
            ablation_id=selected_ablation,
            current=ablation_counters,
            additional=phase2_result.get("ablation_counts"),
        )
        _write_ablation_events(
            output_dir,
            ablation_id=selected_ablation,
            counters=ablation_counters,
        )
        _write_ablation_events(
            output_dir,
            ablation_id=selected_ablation,
            counters=ablation_counters,
            phase=2,
        )
        _publish_phase_state(
            output_dir=output_dir,
            phase=2,
            run_fingerprint=run_fingerprint,
            artifacts=_phase_artifacts(output_dir, 2),
            outcome=phase2_result.get("phase_outcome", "complete"),
        )

    # === PHASE 3: Evidence Synthesis ===
    phase3_result = _run_phase3_subflow(
        # Core inputs from Phase 0
        masked_path=masked_path,
        proteome_path=proteome_path,
        # Core inputs from Phase 1
        validated_markers=validated_markers,
        host_signatures=host_signatures,
        host_signature_model=host_signature_model,
        host_signature_model_payload=host_signature_model_payload,
        # Core inputs from Phase 2
        refined_boundaries=refined_boundaries,
        boundary_taxonomy_map=boundary_taxonomy_map,
        boundary_control_stats=boundary_control_stats,
        boundary_diamond_query=boundary_diamond_query,
        proteome_index=proteome_index,
        boundaries_bed=boundaries_bed,
        merged_seeds=merged_seeds,
        # Core identifiers
        output_dir=output_dir,
        genome_id=genome_id,
        # Resume configuration
        resume=False,
        validated_hits_tsv=validated_hits_tsv,
        # Database parameters
        gene_taxonomy_faa_db=gene_taxonomy_faa_db,
        marker_db=marker_db,
        marker_faa_db=marker_faa_db,
        marker_faa_dir=marker_faa_dir,
        faa_dir=faa_dir,
        diamond_db=diamond_db,
        gvclass_db=gvclass_db,
        hmm_database=hmm_database,
        viral_structure_db=viral_structure_db,
        taxonomy_labels_file=taxonomy_labels_file,
        # Host configuration
        host_prefixes=host_prefixes,
        host_label=host_label,
        high_pident_host_threshold=high_pident_host_threshold,
        boundary_host_trim_score_threshold=boundary_host_trim_score_threshold,
        host_signature_evidence_threshold=host_signature_evidence_threshold,
        boundary_diamond_flank_genes=boundary_diamond_flank_genes,
        # Verification parameters
        high_tier_threshold=high_tier_threshold,
        low_tier_threshold=low_tier_threshold,
        use_crf_in_final_score=use_crf_in_final_score,
        priority_marker_list=priority_marker_list,
        marker_floor_priority_only=marker_floor_priority_only,
        marker_floor_priority_plus_family=marker_floor_priority_plus_family,
        marker_floor_priority_multi_family=marker_floor_priority_multi_family,
        marker_family_bonus_per_family=marker_family_bonus_per_family,
        marker_multi_family_bonus=marker_multi_family_bonus,
        enable_phylogenetic=enable_phylogenetic,
        # Structural analysis parameters
        skip_structural=skip_structural,
        use_boltz=use_boltz,
        boltz_mcp_only=boltz_mcp_only,
        boltz_use_msa_server=boltz_use_msa_server,
        boltz_min_seq_len=boltz_min_seq_len,
        boltz_max_seq_len=boltz_max_seq_len,
        boltz_no_kernels=boltz_no_kernels,
        use_tmvec_database=use_tmvec_database,
        tmvec_require_gpu=tmvec_require_gpu,
        tmvec_databases=tmvec_databases,
        tmvec_database_dir=tmvec_database_dir,
        tmvec_min_score=tmvec_min_score,
        device=device,
        # InterProScan parameters
        interproscan_enabled=interproscan_enabled,
        interproscan_dir=interproscan_dir,
        interproscan_keywords=interproscan_keywords,
        interproscan_threads=interproscan_threads,
        interproscan_applications=interproscan_applications,
        # Database rebuild (for Phase 3 fallback)
        rebuild_db=rebuild_db,
        # Threading
        threads=threads,
        gene_taxonomy_threads=gene_taxonomy_threads,
        # Logger
        logger=logger,
        resume_authorized=False,
        ablation_id=selected_ablation,
        assembly_mode=assembly_mode,
    )
    report_progress(90, "phase 3 complete")

    # Release GPU memory after Phase 3 (TMVec models ~45 GiB).
    try:
        from virosync.utils.gpu import release_gpu_memory
        release_gpu_memory()
    except Exception:
        pass  # Non-critical; best-effort cleanup

    # Extract Phase 3 outputs
    verification_results = phase3_result["verification_results"]
    accepted_results = phase3_result["accepted_results"]
    promoted_low_results = phase3_result.get("promoted_low_results", [])
    classification_stats = phase3_result["classification_stats"]
    phase3_elapsed = phase3_result["elapsed"]
    precomputed_tmvec = phase3_result["precomputed_tmvec"]
    ablation_counters = _merge_ablation_counts(
        ablation_id=selected_ablation,
        current=ablation_counters,
        additional=phase3_result.get("ablation_counts"),
    )

    # Write evidence_graph.json with coherence analyses from all candidates
    from virosync.pipeline.phase3.evidence_graph import write_evidence_graph_json

    coherence_analyses = [
        r.coherence_analysis for r in verification_results if r.coherence_analysis
    ]
    if coherence_analyses:
        evidence_graph_path = output_dir / "phase3_synthesis" / "evidence_graph.json"
        evidence_graph_path.parent.mkdir(parents=True, exist_ok=True)
        write_evidence_graph_json(coherence_analyses, evidence_graph_path, genome_id)

    # Classify MCP proteins as DJR/SJR using multi-signal approach
    marker_validation_dir = output_dir / "phase1" / "marker_validation"
    marker_hits_path = marker_validation_dir / "validated_marker_hits.tsv"
    sequences_path = marker_validation_dir / "hmm_hit_porfs.faa"
    jelly_roll_output = output_dir / "phase3_synthesis" / "virosync_jelly_roll_proteins.tsv"

    if jelly_roll_output.exists():
        logger.info("Jelly roll classification already available: %s", jelly_roll_output)
    elif marker_hits_path.exists() and sequences_path.exists():
        from virosync.orchestration.tasks import classify_jelly_roll_task

        # Build paths for optional signal sources
        interproscan_batch_path = output_dir / "phase3" / "interproscan" / "interproscan_batch.tsv"
        tmvec_results_path = output_dir / "phase3_synthesis" / "virosync_tmvec_proteins.tsv"
        foldseek_results_path = output_dir / "structural_analysis" / "foldseek_pdb_results.tsv"

        logger.info("Classifying MCP proteins as DJR/SJR")
        call_task(
            classify_jelly_roll_task,
            marker_hits_path=marker_hits_path,
            sequences_path=sequences_path,
            output_path=jelly_roll_output,
            interproscan_path=interproscan_batch_path if interproscan_batch_path.exists() else None,
            tmvec_results_path=tmvec_results_path if tmvec_results_path.exists() else None,
            foldseek_results_path=foldseek_results_path if foldseek_results_path.exists() else None,
        )
    else:
        logger.info("Skipping jelly roll classification (marker validation files not found)")

    # Generate ALL outputs to phase3_synthesis/ (intermediate files)
    candidate_output_dir = output_dir / "phase3_synthesis"
    output_files_all = call_task(
        generate_outputs_task,
        verification_results=verification_results,
        output_dir=candidate_output_dir,
        genome_path=masked_path,
        proteome_path=proteome_path,
        accepted_only=False,
        extended_output=extended_output,
        seed_marker_allowlist=seed_marker_allowlist,
        export_all_eve_sequences=export_all_eve_sequences,
        canonical_results=accepted_results,
        promoted_low_results=promoted_low_results,
    )

    # === FINAL OUTPUTS TO ROOT ===
    final_outputs = {}
    invariant_report = None
    invariant_report_path = None

    # 1. Detailed predictions TSV (copy from phase3_synthesis to root)
    src = candidate_output_dir / "virosync_predictions_detailed.tsv"
    dst = output_dir / "virosync_predictions_detailed.tsv"
    detailed_generated = src.is_file()
    if detailed_generated:
        shutil.copy(src, dst)
        final_outputs["predictions_detailed"] = str(dst)

    gene_taxonomy_all = output_files_all.get("gene_taxonomy_all")
    gene_taxonomy_all_path = (
        Path(gene_taxonomy_all) if isinstance(gene_taxonomy_all, (str, Path)) else None
    )
    if gene_taxonomy_all_path is not None and not gene_taxonomy_all_path.exists():
        gene_taxonomy_all_path = None

    invariant_report_path = output_dir / "virosync_tsv_invariant_report.tsv"
    detailed_for_check = dst if detailed_generated else src
    try:
        invariant_report = enforce_tsv_invariants(
            detailed_tsv=detailed_for_check,
            report_out=invariant_report_path,
            gene_taxonomy_all_tsv=gene_taxonomy_all_path,
        )
    except TSVInvariantError as exc:
        _clear_success_markers(output_dir)
        logger.error("%s", exc)
        raise
    final_outputs["tsv_invariant_report"] = str(invariant_report_path)

    if invariant_report.warning_count:
        preview = "; ".join(
            f"{issue.eve_id}:{issue.check}"
            for issue in invariant_report.warning_issues[:5]
        )
        logger.warning(
            "Detailed TSV invariant check passed with warnings: "
            "rows=%d warnings=%d (%s)",
            invariant_report.rows_checked,
            invariant_report.warning_count,
            preview if preview else "no preview",
        )
    else:
        logger.info(
            "Detailed TSV invariant check passed: rows=%d issues=%d",
            invariant_report.rows_checked,
            invariant_report.issue_count,
        )

    in_memory_accepted = int(phase3_result["accepted"])
    if in_memory_accepted != len(accepted_results):
        raise RuntimeError(
            "in-memory Phase 3 accepted count disagrees with accepted results: "
            f"count={in_memory_accepted} results={len(accepted_results)}"
        )
    persisted_summary = _summarize_prediction_outputs(
        candidate_output_dir / "virosync_predictions.tsv",
        candidate_output_dir / "virosync_predictions_detailed.tsv",
        expected_accepted=len(accepted_results),
        expected_candidates=len(verification_results),
    )
    for eve_class, count_key in EFFECTIVE_EVE_CLASS_COUNT_KEYS.items():
        expected_count = int(classification_stats.get(eve_class, 0) or 0)
        if persisted_summary[count_key] != expected_count:
            raise RuntimeError(
                "persisted effective-class count disagrees with in-memory Phase 3: "
                f"class={eve_class} expected={expected_count} "
                f"persisted={persisted_summary[count_key]}"
            )

    accepted = persisted_summary["accepted"]
    tier_counts = {
        "HIGH": persisted_summary["high_tier"],
        "MEDIUM": persisted_summary["medium_tier"],
        "LOW": persisted_summary["low_tier"],
    }
    candidate_tier_counts = {
        "HIGH": persisted_summary["candidate_high_tier"],
        "MEDIUM": persisted_summary["candidate_medium_tier"],
        "LOW": persisted_summary["candidate_low_tier"],
    }
    quality_gate_dropped = persisted_summary["quality_gate_dropped"]

    # 2. Combined EVE FASTA: {genome_id}_eves.fna
    eve_fasta = _write_combined_eve_fasta(
        output_dir=output_dir,
        genome_id=genome_id,
        genome_path=masked_path,
        results=accepted_results,
    )
    if eve_fasta is not None:
        final_outputs["eve_fasta"] = str(eve_fasta)

    # 3. GVClass results (if enabled)
    gvclass_results_path = None
    if run_gvclass and gvclass_path:
        from virosync.pipeline.phase3.gvclass_runner import (
            load_gvclass_id_map,
            parse_gvclass_results,
            run_gvclass_batch,
            write_gvclass_results_tsv,
        )
        gvclass_input = candidate_output_dir / "gvclass_input" / "nucleotide"
        if gvclass_input.exists() and list(gvclass_input.glob("*.fna")):
            logger.info(f"Running GVClass on {len(list(gvclass_input.glob('*.fna')))} EVE sequences")
            gvclass_out = candidate_output_dir / "gvclass_output"
            summary = run_gvclass_batch(
                gvclass_input, gvclass_out, Path(gvclass_path),
                threads=threads, gvclass_db=gvclass_db
            )
            if summary:
                id_map = load_gvclass_id_map(
                    candidate_output_dir / "gvclass_input" / "manifest.tsv"
                )
                results_dict = parse_gvclass_results(summary, id_map=id_map)
                gvclass_tsv = output_dir / "gvclass_results.tsv"
                write_gvclass_results_tsv(results_dict, gvclass_tsv)
                final_outputs["gvclass_results"] = str(gvclass_tsv)
                gvclass_results_path = gvclass_tsv
                logger.info(f"GVClass results written: {gvclass_tsv}")
        else:
            logger.warning("GVClass skipped: no EVE sequences found in gvclass_input/nucleotide/")

    output_files = {
        "candidates": output_files_all,
        "final": final_outputs,
        "phase2_boundaries": str(boundaries_bed),
        "evidence_graph": str(output_dir / "phase3_synthesis" / "evidence_graph.json"),
        "masking_status": str(masking_result.status_path),
    }

    output_files.update(
        _generate_required_reports(
            output_dir=output_dir,
            genome_id=genome_id,
            taxonomy_labels_file=taxonomy_labels_file,
            logger=logger,
        )
    )

    # Create summary artifact placeholder
    call_task(
        create_summary_artifact_task,
        genome_id=genome_id,
        n_seeds=len(merged_seeds),
        n_boundaries=len(refined_boundaries),
        n_verified=len(verification_results),
        tier_counts=tier_counts,
        output_files=output_files,
    )

    # Calculate final timing
    total_elapsed = time.time() - genome_start_time

    # Write per-genome run log for batch tracking. This file is the resume
    # completion signal: its presence is checked by batch orchestration to
    # decide whether to skip a genome. Written atomically so that a SIGKILL
    # mid-write cannot leave a truncated "complete" marker.
    run_log_path = output_dir / "run.log"
    with atomic_write_context(run_log_path, "w") as f:
        f.write(f"# ViroSync Run Log: {genome_id}\n")
        f.write(f"# Generated: {datetime.now(timezone.utc).isoformat()}\n")
        f.write(f"# Input: {genome_path}\n")
        f.write(f"# Output: {output_dir}\n")
        f.write(f"# Total time: {total_elapsed:.1f}s\n")
        f.write("#" + "=" * 59 + "\n\n")
        f.write("## Timing\n")
        f.write(f"Phase 0 (masking + proteome): {phase0_elapsed:.1f}s\n")
        f.write(f"Phase 1 (seeding): {phase1_elapsed:.1f}s\n")
        f.write(f"Phase 2 (boundary refinement): {phase2_elapsed:.1f}s\n")
        f.write(f"Phase 3 (verification): {phase3_elapsed:.1f}s\n")
        f.write(f"Total: {total_elapsed:.1f}s\n\n")
        f.write("## Results Summary\n")
        f.write(f"Seeds identified: {len(merged_seeds)}\n")
        f.write(f"Boundaries refined: {len(refined_boundaries)}\n")
        f.write(f"Phase-3 candidates: {len(verification_results)}\n")
        f.write(f"Canonical EVEs: {accepted}\n")
        f.write(f"Quality-gate dropped: {quality_gate_dropped}\n")
        f.write(
            "  Candidate tiers: "
            f"HIGH={candidate_tier_counts['HIGH']}, "
            f"MEDIUM={candidate_tier_counts['MEDIUM']}, "
            f"LOW={candidate_tier_counts['LOW']}\n"
        )
        f.write(
            "  Canonical tiers: "
            f"HIGH={tier_counts['HIGH']}, "
            f"MEDIUM={tier_counts['MEDIUM']}, "
            f"LOW={tier_counts['LOW']} "
            f"(normal-gate-promoted={len(promoted_low_results)})\n\n"
        )
        if invariant_report is not None:
            f.write("## Detailed TSV Invariant Check\n")
            f.write(f"Status: {invariant_report.status}\n")
            f.write(f"Rows checked: {invariant_report.rows_checked}\n")
            f.write(f"Issues: {invariant_report.issue_count}\n")
            f.write(f"Errors: {invariant_report.error_count}\n")
            f.write(f"Warnings: {invariant_report.warning_count}\n")
            if invariant_report_path is not None:
                f.write(f"Report: {invariant_report_path}\n")
            if invariant_report.issue_count:
                for issue in invariant_report.issues[:10]:
                    f.write(
                        f"  {issue.eve_id}\t{issue.check}\t{issue.message}\n"
                    )
            f.write("\n")
        if host_deviation_summary:
            f.write("## Phase 1 Host-Taxonomy Deviation\n")
            f.write("Enabled: true\n")
            f.write(f"Deviation markers: {host_deviation_summary.get('markers_total', 0)}\n")
            f.write(f"Seedable deviation markers: {host_deviation_summary.get('markers_seedable', 0)}\n")
            baseline = host_deviation_summary.get("baseline") or {}
            if baseline:
                baseline_source = baseline.get("baseline_source", "unknown")
                f.write(f"Baseline source: {baseline_source}\n")
                f.write(
                    f"Baseline markers: {baseline.get('baseline_marker_hits', 0)} "
                    f"tokens: {baseline.get('baseline_token_count', 0)}\n"
                )
                if "background_window_used" in baseline:
                    f.write(
                        "Background windows: "
                        f"{baseline.get('background_window_used', 0)}/"
                        f"{baseline.get('background_window_count', 0)} "
                        f"size={baseline.get('background_window_bp', 0)}bp\n"
                    )
            f.write(f"Report: {host_deviation_summary.get('report_path', '')}\n\n")
        f.write("## Output Files\n")
        for key, path in output_files.items():
            if isinstance(path, dict):
                for subkey, subpath in path.items():
                    f.write(f"  {key}.{subkey}: {subpath}\n")
            else:
                f.write(f"  {key}: {path}\n")
        f.write("\n## Detected EVEs by Tier\n")
        for tier in ["HIGH", "MEDIUM", "LOW"]:
            tier_results = [vr for vr in accepted_results if vr.confidence_tier == tier]
            if tier_results:
                f.write(f"\n### {tier} Confidence ({len(tier_results)})\n")
                for vr in tier_results:
                    region_bp = vr.end - vr.start
                    genes = getattr(vr, 'gene_count', 0)
                    f.write(f"  {vr.eve_id}: {vr.start}-{vr.end} ({region_bp:,} bp, {genes} genes)\n")
                    f.write(f"    confidence: {vr.final_confidence:.3f}\n")
    logger.info(f"Run log written: {run_log_path}")

    # Capture tool/DB/input provenance for reproducibility (best-effort).
    try:
        from virosync.utils.provenance import write_provenance
        write_provenance(
            output_dir,
            {
                "hmm_database": str(hmm_database) if hmm_database else None,
                "marker_db": str(marker_db) if marker_db else None,
                "gene_taxonomy_faa_db": str(gene_taxonomy_faa_db) if gene_taxonomy_faa_db else None,
                "taxonomy_labels": str(taxonomy_labels_file) if taxonomy_labels_file else None,
                "use_tmvec_database": use_tmvec_database,
                "use_interproscan": interproscan_enabled,
                "masking": masking,
                "masking_status_path": str(masking_result.status_path),
                "masking_status_sha256": masking_result.status_sha256,
            },
            input_genome=genome_path,
        )
    except Exception as exc:  # provenance is non-critical; never fail the run for it
        logger.warning("Provenance capture failed (non-fatal): %s", exc)

    completion_manifest = _write_completion_manifest(
        output_dir=output_dir,
        genome_id=genome_id,
        status="success",
        output_files=output_files,
        fingerprint=config_fingerprint,
    )
    output_files["completion_manifest"] = str(completion_manifest)
    logger.info("Completion manifest written: %s", completion_manifest)

    _write_ablation_events(
        output_dir,
        ablation_id=selected_ablation,
        counters=ablation_counters,
        phase=3,
    )
    _write_ablation_events(
        output_dir,
        ablation_id=selected_ablation,
        counters=ablation_counters,
    )
    final_artifacts = _final_artifacts(
        output_dir,
        ablation_id=selected_ablation,
    )
    phase3_artifacts_by_path = {
        artifact.relative_path: artifact
        for artifact in (*_phase_artifacts(output_dir, 3), *final_artifacts)
    }
    _publish_phase_state(
        output_dir=output_dir,
        phase=3,
        run_fingerprint=run_fingerprint,
        artifacts=tuple(
            phase3_artifacts_by_path[path]
            for path in sorted(phase3_artifacts_by_path)
        ),
    )
    publish_run_success(
        output_dir,
        run_fingerprint=run_fingerprint,
        artifacts=final_artifacts,
        result=_result_identity(
            persisted_summary,
            terminal_phase=None,
            benchmark_eligible=masking_result.benchmark_eligible,
            promoted_low_rows=len(promoted_low_results),
        ),
    )
    output_files = _authenticated_output_files(output_dir)
    logger.info("Schema-v3 run state written: %s", output_files["run_state"])

    # Final summary
    logger.info("=" * 60)
    logger.info(f"GENOME COMPLETE: {genome_id}")
    logger.info(f"  Total time: {total_elapsed:.1f}s")
    logger.info(
        f"  Result: {accepted} canonical EVEs "
        f"({len(verification_results)} candidates; "
        f"HIGH={tier_counts['HIGH']}, MEDIUM={tier_counts['MEDIUM']}, "
        f"LOW={tier_counts['LOW']}, "
        f"normal-gate-promoted={len(promoted_low_results)})"
    )
    logger.info("=" * 60)

    report_progress(100, "complete")
    return {
        "genome_id": genome_id,
        "success": True,
        "benchmark_eligible": masking_result.benchmark_eligible,
        "legacy_resume": False,
        "output_files": output_files,
        "elapsed_sec": total_elapsed,
        **persisted_summary,
        "tsv_invariant_issues": invariant_report.issue_count if invariant_report else 0,
    }


def single_genome_flow(
    genome_path: Path,
    output_dir: Path,
    genome_id: str,
    # Configuration object (recommended - use instead of individual kwargs)
    config: Optional[PipelineConfig] = None,
    ablation_id: str = "A0",
    # Database paths (can override config)
    hmm_database: Optional[Path] = None,
    hmm_allowlist: Optional[Path] = None,
    seed_marker_allowlist: Optional[list[str]] = None,
    marker_faa_db: Optional[Path] = None,
    marker_db: Optional[Path] = None,
    gene_taxonomy_faa_db: Optional[Path] = None,
    marker_faa_dir: Optional[Path] = None,
    faa_dir: Optional[Path] = None,
    gvclass_db: Optional[Path] = None,
    diamond_db: Optional[Path] = None,
    enable_phylogenetic: bool = False,
    # Taxonomy lookup for host signature comparison
    taxonomy_labels_file: Optional[Path] = None,
    # Host taxonomy configuration
    host_prefixes: Optional[list[str]] = None,
    host_label: str = "EUK",
    high_pident_host_threshold: float = 70.0,
    # Parameters
    threads: int = 8,
    max_threads: Optional[int] = None,
    device: str = "cpu",
    search_backend: str = "diamond",
    masking: Optional[MaskingConfig] = None,
    skip_masking: Optional[bool] = None,
    skip_structural: bool = True,
    use_boltz: bool = False,
    boltz_mcp_only: bool = True,
    boltz_use_msa_server: bool = False,
    boltz_min_seq_len: int = 100,
    boltz_max_seq_len: int = 1000,
    boltz_no_kernels: bool = True,
    use_tmvec_database: bool = False,
    tmvec_require_gpu: bool = False,
    tmvec_databases: Optional[list[str]] = None,
    tmvec_database_dir: Optional[Path] = None,
    tmvec_min_score: float = 0.5,
    viral_structure_db: Optional[Path] = None,
    assembly_mode: str = "default",
    high_tier_threshold: float = 0.7,
    low_tier_threshold: float = 0.2,
    use_crf_in_final_score: bool = False,
    priority_marker_list: Optional[list[str]] = None,
    marker_floor_priority_only: float = 0.55,
    marker_floor_priority_plus_family: float = 0.70,
    marker_floor_priority_multi_family: float = 0.80,
    marker_family_bonus_per_family: float = 0.06,
    marker_multi_family_bonus: float = 0.08,
    hmm_chunk_size: Optional[int] = None,
    gene_taxonomy_threads: Optional[int] = None,
    interproscan_enabled: bool = False,
    interproscan_dir: Optional[Path] = None,
    interproscan_keywords: Optional[list[str]] = None,
    interproscan_threads: Optional[int] = None,
    interproscan_applications: Optional[list[str]] = None,
    extended_output: bool = True,
    export_all_eve_sequences: bool = True,
    # GVClass batch classification
    run_gvclass: bool = False,
    gvclass_path: Optional[Path] = None,
    # HMM-gated workflow options
    rebuild_db: bool = False,
    initial_window_bp: int = 10000,
    initial_window_genes: int = 5,
    min_markers_initial: int = 1,
    extension_kb: int = 5,
    merge_distance: int = 1000,
    host_taxonomy_deviation_enabled: bool = False,
    host_taxonomy_deviation_allow_seeds: bool = False,
    host_taxonomy_deviation_min_token_len: int = 6,
    host_taxonomy_deviation_min_tokens: int = 3,
    host_taxonomy_deviation_overlap_threshold: float = 0.3,
    host_taxonomy_deviation_max_pident: float = 70.0,
    host_taxonomy_deviation_max_hits: int = 5,
    host_taxonomy_deviation_window_bp: int = 5000,
    host_taxonomy_deviation_window_count: int = 25,
    host_taxonomy_deviation_window_seed: int = 13,
    host_taxonomy_deviation_window_min_markers: int = 1,
    host_taxonomy_deviation_seed_window_bp: int = 10000,
    host_taxonomy_deviation_seed_min_markers: int = 3,
    marker_validation_top_k: int = 10,
    novel_marker_min_score: float = 30.0,
    novel_marker_min_coverage: float = 0.5,
    novel_marker_require_cluster: bool = True,
    boundary_host_trim_enabled: bool = True,
    boundary_host_trim_window_bp: int = 5000,
    boundary_host_trim_step_bp: int = 1000,
    boundary_host_trim_max_host_fraction: float = 0.3,
    boundary_host_trim_min_viral_fraction: float = 0.05,
    boundary_host_trim_score_threshold: float = 0.3,
    host_signature_evidence_threshold: float = 0.3,
    boundary_host_trim_buffer_kb: int = 5,
    boundary_host_signature_min_token_len: int = 3,
    boundary_host_trim_min_overlap_score: float = 0.40,
    boundary_taxonomy_ml_enabled: bool = False,
    boundary_taxonomy_ml_model: str = "logreg",
    boundary_taxonomy_ml_threshold: float = 0.5,
    boundary_taxonomy_ml_neighbor_window: int = 3,
    taxonomy_weight_mode: str = "rank",
    # Batched Diamond configuration (Phase 2)
    # Note: Batched Diamond is always enabled (use_batched_boundary_diamond removed Jan 2026)
    boundary_diamond_flank_genes: int = 10,  # Seeds already extended by ±5 genes
    boundary_diamond_control_sample_size: int = 100,
    boundary_diamond_control_min_distance: int = 30,
    boundary_diamond_top_k: int = 10,
    boundary_diamond_chunk_size: int = 10000,
    boundary_diamond_random_seed: int = 42,
    boundary_diamond_superset_prototype_enabled: bool = False,
    resume: bool = True,
    progress_callback: Optional[Callable[[float, str, bool], None]] = None,
) -> dict:
    """
    Process a single genome through all ViroSync phases.

    This function orchestrates:
    - Phase 0: Preprocessing (masking + proteome generation)
    - Phase 1: HMM scan -> marker validation -> region assembly
    - Phase 2: Boundary refinement (parallel per seed)
    - Phase 3: Evidence synthesis (parallel per boundary)

    Args:
        genome_path: Path to input genome FASTA
        output_dir: Output directory for results
        genome_id: Unique identifier for this genome
        config: Optional PipelineConfig (recommended - overrides individual params)
        ... (additional parameters - see implementation)

    Returns:
        Dictionary with results summary
    """
    # Fast path: no config, pass params directly
    if config is None:
        effective_masking = masking or MaskingConfig()
        if skip_masking is not None:
            effective_masking = (
                effective_masking.with_backend(MaskingBackend.OFF)
                if skip_masking
                else effective_masking.with_backend(
                    MaskingBackend.TRF_REPEATMASKER
                )
                if effective_masking.backend is MaskingBackend.OFF
                else effective_masking
            )
        return _single_genome_flow_impl(
            genome_path=genome_path,
            output_dir=output_dir,
            genome_id=genome_id,
            ablation_id=ablation_id,
            ablation_contract_sha256=ABLATION_CONTRACT_SHA256,
            hmm_database=hmm_database,
            hmm_allowlist=hmm_allowlist,
            seed_marker_allowlist=seed_marker_allowlist,
            marker_faa_db=marker_faa_db,
            marker_db=marker_db,
            gene_taxonomy_faa_db=gene_taxonomy_faa_db,
            marker_faa_dir=marker_faa_dir,
            faa_dir=faa_dir,
            gvclass_db=gvclass_db,
            diamond_db=diamond_db,
            enable_phylogenetic=enable_phylogenetic,
            taxonomy_labels_file=taxonomy_labels_file,
            host_prefixes=host_prefixes,
            host_label=host_label,
            high_pident_host_threshold=high_pident_host_threshold,
            threads=threads,
            max_threads=max_threads,
            device=device,
            search_backend=search_backend,
            masking=effective_masking,
            skip_structural=skip_structural,
            use_boltz=use_boltz,
            boltz_mcp_only=boltz_mcp_only,
            boltz_use_msa_server=boltz_use_msa_server,
            boltz_min_seq_len=boltz_min_seq_len,
            boltz_max_seq_len=boltz_max_seq_len,
            boltz_no_kernels=boltz_no_kernels,
            use_tmvec_database=use_tmvec_database,
            tmvec_require_gpu=tmvec_require_gpu,
            tmvec_databases=tmvec_databases,
            tmvec_database_dir=tmvec_database_dir,
            tmvec_min_score=tmvec_min_score,
            viral_structure_db=viral_structure_db,
            assembly_mode=assembly_mode,
            high_tier_threshold=high_tier_threshold,
            low_tier_threshold=low_tier_threshold,
            use_crf_in_final_score=use_crf_in_final_score,
            priority_marker_list=priority_marker_list,
            marker_floor_priority_only=marker_floor_priority_only,
            marker_floor_priority_plus_family=marker_floor_priority_plus_family,
            marker_floor_priority_multi_family=marker_floor_priority_multi_family,
            marker_family_bonus_per_family=marker_family_bonus_per_family,
            marker_multi_family_bonus=marker_multi_family_bonus,
            hmm_chunk_size=hmm_chunk_size,
            gene_taxonomy_threads=gene_taxonomy_threads,
            interproscan_enabled=interproscan_enabled,
            interproscan_dir=interproscan_dir,
            interproscan_keywords=interproscan_keywords,
            interproscan_threads=interproscan_threads,
            interproscan_applications=interproscan_applications,
            extended_output=extended_output,
            export_all_eve_sequences=export_all_eve_sequences,
            run_gvclass=run_gvclass,
            gvclass_path=gvclass_path,
            rebuild_db=rebuild_db,
            initial_window_bp=initial_window_bp,
            initial_window_genes=initial_window_genes,
            min_markers_initial=min_markers_initial,
            extension_kb=extension_kb,
            merge_distance=merge_distance,
            host_taxonomy_deviation_enabled=host_taxonomy_deviation_enabled,
            host_taxonomy_deviation_allow_seeds=host_taxonomy_deviation_allow_seeds,
            host_taxonomy_deviation_min_token_len=host_taxonomy_deviation_min_token_len,
            host_taxonomy_deviation_min_tokens=host_taxonomy_deviation_min_tokens,
            host_taxonomy_deviation_overlap_threshold=host_taxonomy_deviation_overlap_threshold,
            host_taxonomy_deviation_max_pident=host_taxonomy_deviation_max_pident,
            host_taxonomy_deviation_max_hits=host_taxonomy_deviation_max_hits,
            host_taxonomy_deviation_window_bp=host_taxonomy_deviation_window_bp,
            host_taxonomy_deviation_window_count=host_taxonomy_deviation_window_count,
            host_taxonomy_deviation_window_seed=host_taxonomy_deviation_window_seed,
            host_taxonomy_deviation_window_min_markers=host_taxonomy_deviation_window_min_markers,
            host_taxonomy_deviation_seed_window_bp=host_taxonomy_deviation_seed_window_bp,
            host_taxonomy_deviation_seed_min_markers=host_taxonomy_deviation_seed_min_markers,
            marker_validation_top_k=marker_validation_top_k,
            novel_marker_min_score=novel_marker_min_score,
            novel_marker_min_coverage=novel_marker_min_coverage,
            novel_marker_require_cluster=novel_marker_require_cluster,
            boundary_host_trim_enabled=boundary_host_trim_enabled,
            boundary_host_trim_window_bp=boundary_host_trim_window_bp,
            boundary_host_trim_step_bp=boundary_host_trim_step_bp,
            boundary_host_trim_max_host_fraction=boundary_host_trim_max_host_fraction,
            boundary_host_trim_min_viral_fraction=boundary_host_trim_min_viral_fraction,
            boundary_host_trim_score_threshold=boundary_host_trim_score_threshold,
            host_signature_evidence_threshold=host_signature_evidence_threshold,
            boundary_host_trim_buffer_kb=boundary_host_trim_buffer_kb,
            boundary_host_signature_min_token_len=boundary_host_signature_min_token_len,
            boundary_host_trim_min_overlap_score=boundary_host_trim_min_overlap_score,
            boundary_taxonomy_ml_enabled=boundary_taxonomy_ml_enabled,
            boundary_taxonomy_ml_model=boundary_taxonomy_ml_model,
            boundary_taxonomy_ml_threshold=boundary_taxonomy_ml_threshold,
            boundary_taxonomy_ml_neighbor_window=boundary_taxonomy_ml_neighbor_window,
            taxonomy_weight_mode=taxonomy_weight_mode,
            boundary_diamond_flank_genes=boundary_diamond_flank_genes,
            boundary_diamond_control_sample_size=boundary_diamond_control_sample_size,
            boundary_diamond_control_min_distance=boundary_diamond_control_min_distance,
            boundary_diamond_top_k=boundary_diamond_top_k,
            boundary_diamond_chunk_size=boundary_diamond_chunk_size,
            boundary_diamond_random_seed=boundary_diamond_random_seed,
            boundary_diamond_superset_prototype_enabled=(
                boundary_diamond_superset_prototype_enabled
            ),
            resume=resume,
            progress_callback=progress_callback,
        )

    # Detect explicit overrides (values different from signature defaults)
    explicit_overrides = _detect_explicit_overrides(
        signature=inspect.signature(single_genome_flow),
        passed_kwargs=locals().copy(),
        exclude_keys={
            'config',
            'genome_path',
            'output_dir',
            'genome_id',
            'progress_callback',
        },
    )
    if "skip_masking" in explicit_overrides:
        if "masking" in explicit_overrides:
            config = config.with_overrides(
                masking=explicit_overrides.pop("masking")
            )
        config = config.with_overrides(
            skip_masking=explicit_overrides.pop("skip_masking")
        )

    # Merge config with explicit overrides
    merged = _merge_config_with_kwargs(
        config=config,
        explicit_kwargs=explicit_overrides,
        exclude_keys={
            'config',
            'genome_path',
            'output_dir',
            'genome_id',
            'progress_callback',
        },
    )

    # Filter merged kwargs to match impl signature and fill missing params with defaults
    filtered = _filter_kwargs_to_signature(
        kwargs=merged,
        target_signature=inspect.signature(_single_genome_flow_impl),
        defaults_signature=inspect.signature(single_genome_flow),
    )
    filtered.pop("progress_callback", None)

    # Call implementation with filtered values
    return _single_genome_flow_impl(
        genome_path=genome_path,
        output_dir=output_dir,
        genome_id=genome_id,
        progress_callback=progress_callback,
        **filtered,
    )


run_single_genome_task = single_genome_flow
