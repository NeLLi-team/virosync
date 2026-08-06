"""ViroSync orchestration CLI commands."""

import click
import yaml
import copy
import hashlib
import importlib.util
import logging
import os
import sys
import shutil
import subprocess
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional
from click.core import ParameterSource

import virosync
from virosync.orchestration.python_runner import (
    BatchProgress,
    _batch_result_status,
    _preflight_genome_runs,
    run_batch_python,
)
from virosync.config import (
    ApplicationConfig,
    ConfigError,
    FeatureResolution,
    PipelineConfig,
)
from virosync.utils.database_manager import ViroSyncDatabaseManager
from virosync.utils.executables import resolve_boltz_executable

# Supported genome file extensions
GENOME_EXTENSIONS = {".fna", ".fasta", ".fa"}
GVCLASS_PATH_ENV_VAR = "VIROSYNC_GVCLASS_PATH"
_TMVEC_CHECKPOINT_MD5 = {
    "tm_vec_swiss_model_large_params.json": "fbb1f2288be74ad6c5ac1c05a19f876d",
    "tm_vec_swiss_model_large.ckpt": "69d8ef7a3286b8f6077fb89310ac19dd",
}


def _collect_genome_paths(input_path: Path) -> list[Path]:
    """
    Collect genome paths from input.

    Input can be:
    - A single genome file (.fna, .fasta, .fa)
    - A directory containing genome files
    - A text file with genome paths (one per line)

    Returns:
        List of Path objects to genome files
    """
    input_path = Path(input_path)

    if input_path.is_file():
        # Check if it's a genome file or a list file
        suffix = "".join(input_path.suffixes).lower()
        if any(suffix.endswith(ext) for ext in GENOME_EXTENSIONS):
            # Single genome file
            return [input_path]
        else:
            # Assume it's a list file with genome paths
            genome_paths = [
                Path(line.strip())
                for line in input_path.read_text().splitlines()
                if line.strip() and not line.startswith("#")
            ]
            return genome_paths
    elif input_path.is_dir():
        # Directory: find all genome files
        genome_paths = []
        for ext in GENOME_EXTENSIONS:
            genome_paths.extend(input_path.glob(f"*{ext}"))
        return sorted(genome_paths)
    else:
        raise click.ClickException(f"Input path does not exist: {input_path}")


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def orchestrate():
    """Orchestration commands for ViroSync pipeline execution."""
    pass


def _load_config(path: Optional[Path]) -> ApplicationConfig:
    """Decode configuration without resolving, checking, or installing resources."""
    try:
        if path is None:
            return ApplicationConfig.from_dict({"schema_version": 1})
        return ApplicationConfig.from_yaml(Path(path))
    except ConfigError as exc:
        raise click.ClickException(str(exc)) from exc


def _gpu_uuid_from_index(index: int) -> Optional[str]:
    try:
        proc = subprocess.run(
            ["nvidia-smi", "--query-gpu=index,uuid", "--format=csv,noheader"],
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return None
    for line in proc.stdout.splitlines():
        parts = [p.strip() for p in line.split(",", 1)]
        if len(parts) != 2:
            continue
        if parts[0].isdigit() and int(parts[0]) == index:
            return parts[1]
    return None


def _apply_gpu_id_env(
    gpu_id: Optional[int],
    configured_gpu_id: Optional[int],
) -> Optional[str]:
    resolved = gpu_id if gpu_id is not None else configured_gpu_id
    if resolved is None:
        return None
    requested = str(resolved).strip()
    selector = requested
    if requested.isdigit():
        uuid = _gpu_uuid_from_index(int(requested))
        if uuid:
            selector = uuid
    os.environ["CUDA_VISIBLE_DEVICES"] = selector
    os.environ["VIROSYNC_GPU"] = requested
    return f"gpu-id={requested}, CUDA_VISIBLE_DEVICES={selector}"


def _effective_config_payload(
    application_config: ApplicationConfig,
    optional_features: dict[str, FeatureResolution],
) -> dict:
    """Build the sole effective configuration artifact representation."""
    return application_config.effective_payload(optional_features)


def _resolve_extension_kb(cli_value: Optional[int], phase1_config: dict) -> int:
    """Resolve extension parameter from config, supporting both extension_bp and extension_kb.

    Priority: CLI value > extension_bp (converted) > extension_kb > default (5)
    """
    if cli_value is not None:
        return cli_value
    # Use is-not-None checks (not just key presence): a config with an explicit
    # ``extension_bp: null`` / ``extension_kb: null`` must fall through to the default,
    # not return None (which crashes the downstream region-extension multiply).
    extension_bp = phase1_config.get("extension_bp")
    if extension_bp is not None:
        return extension_bp // 1000  # convert bp -> kb
    extension_kb = phase1_config.get("extension_kb")
    if extension_kb is not None:
        return extension_kb
    return 5  # Phase1Config.extension_kb default


def _cap_threads_per_worker(
    threads_per_worker: int,
    max_threads: Optional[int],
    max_concurrent_genomes: Optional[int],
) -> tuple[int, Optional[str]]:
    """Cap threads-per-genome so concurrent genomes don't oversubscribe max_threads.

    The genome pool runs up to ``max_concurrent_genomes`` genomes concurrently (the
    python_runner ThreadPoolExecutor), each using ``threads_per_worker`` threads, so the
    budget is ``max_threads // max_concurrent_genomes`` -- NOT ``// workers`` (the prior
    bug oversubscribed whenever max_concurrent_genomes exceeded workers). Returns the
    (possibly reduced) value plus a warning message when it was reduced.
    """
    if not (max_threads and max_concurrent_genomes and max_concurrent_genomes > 0):
        return threads_per_worker, None
    worker_budget = max(1, max_threads // max_concurrent_genomes)
    if threads_per_worker > worker_budget:
        warning = (
            f"Reducing threads-per-genome {threads_per_worker} -> {worker_budget} to keep "
            f"max_threads={max_threads} across {max_concurrent_genomes} concurrent genomes"
        )
        return worker_budget, warning
    return threads_per_worker, None


def _prompt_optional_archive_choice(
    name: str,
    default_target: Path,
    default_source: Optional[str],
) -> tuple[bool, Path, Optional[str]]:
    """Interactive prompt for optional archive-backed resources."""
    should_setup = click.confirm(
        f"Set up {name} now?",
        default=bool(default_source),
    )
    if not should_setup:
        return False, default_target, None

    target_input = click.prompt(f"{name} install location", default=str(default_target))
    target_path = ViroSyncDatabaseManager.normalize_path(target_input)

    if default_source:
        source_input = click.prompt(f"{name} archive URL/path", default=default_source)
    else:
        source_input = click.prompt(
            f"{name} archive URL/path (leave empty to use existing files only)",
            default="",
            show_default=False,
        )
    source = source_input.strip() or None
    return True, target_path, source


def _warn_optional_resource(message: str) -> None:
    click.echo(click.style(f"Warning: {message}", fg="yellow"), err=True)


def _command_output_flags(local_verbose: bool = False) -> tuple[bool, bool]:
    """Resolve root and command-local verbosity and configure logging."""
    root_obj = click.get_current_context().find_root().obj or {}
    verbose = bool(local_verbose or root_obj.get("verbose"))
    quiet = bool(root_obj.get("quiet")) and not verbose
    logging.getLogger().setLevel(
        logging.DEBUG if verbose else logging.ERROR
    )
    return verbose, quiet


def _database_version(config: PipelineConfig) -> str:
    """Return the installed database version, never the requested version."""
    for value in (
        config.databases.hmm_database,
        config.databases.marker_db,
        config.databases.gene_taxonomy_faa_db,
    ):
        if value is None:
            continue
        path = Path(value)
        for candidate in (path.parent, path.parent.parent):
            if (candidate / "DB_VERSION").is_file():
                return ViroSyncDatabaseManager.get_database_version(candidate)
    return "unknown"


def _print_banner(database_version: str) -> None:
    from virosync.cli.main import print_banner

    print_banner(database_version)


def _tmvec_runtime_issues(config: PipelineConfig) -> list[str]:
    """Return local-only TMVec readiness failures without downloading assets."""
    issues = []
    missing_databases = ViroSyncDatabaseManager.missing_tmvec_files(
        tmvec_root=config.phase3.tmvec_database_dir,
        databases=config.phase3.tmvec_databases,
    )
    if missing_databases:
        root = (
            config.phase3.tmvec_database_dir
            or ViroSyncDatabaseManager.default_tmvec_path()
        )
        issues.append(
            f"TMVec databases not available under {root}: "
            + "; ".join(missing_databases)
        )

    dependency_modules = {
        "torch": "torch",
        "transformers": "transformers",
        "pytorch-lightning": "pytorch_lightning",
        "sentencepiece": "sentencepiece",
    }
    missing_dependencies = [
        label
        for label, module in dependency_modules.items()
        if importlib.util.find_spec(module) is None
    ]
    if missing_dependencies:
        issues.append(
            "missing dependencies: " + ", ".join(sorted(missing_dependencies))
        )

    cache_dir = Path.home() / ".cache" / "virosync" / "tmvec"
    checkpoint_files = (
        cache_dir / "tm_vec_swiss_model_large_params.json",
        cache_dir / "tm_vec_swiss_model_large.ckpt",
    )
    def checkpoint_ready(path: Path) -> bool:
        try:
            return path.is_file() and path.stat().st_size > 0
        except OSError:
            return False

    missing_checkpoints = [
        path.name for path in checkpoint_files if not checkpoint_ready(path)
    ]
    if missing_checkpoints:
        issues.append(
            f"missing trained checkpoint files under {cache_dir}: "
            + ", ".join(missing_checkpoints)
        )

    if not missing_checkpoints:
        for path in checkpoint_files:
            digest = hashlib.md5(usedforsecurity=False)
            try:
                with path.open("rb") as handle:
                    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                        digest.update(chunk)
            except OSError as exc:
                issues.append(f"cannot read trained checkpoint {path}: {exc}")
                continue
            if digest.hexdigest() != _TMVEC_CHECKPOINT_MD5[path.name]:
                issues.append(f"trained checkpoint checksum mismatch: {path}")

    if config.compute.device.value != "cuda":
        if config.phase3.tmvec_require_gpu:
            issues.append("compute.device must be 'cuda' when tmvec_require_gpu is true")
        return issues

    if "torch" not in missing_dependencies:
        try:
            import torch

            if not torch.cuda.is_available() or torch.cuda.device_count() < 1:
                issues.append("selected CUDA device is not available")
            else:
                torch.cuda.get_device_properties(0)
        except Exception as exc:
            issues.append(f"selected CUDA device validation failed: {exc}")
    return issues


def _resolve_optional_features(
    config: PipelineConfig,
) -> tuple[PipelineConfig, dict[str, FeatureResolution]]:
    """Resolve optional features once, recording requested and effective state."""
    resolved = copy.deepcopy(config)
    states: dict[str, FeatureResolution] = {}

    boltz_requested = resolved.phase3.use_boltz
    boltz_issues = []
    if boltz_requested:
        boltz_db = resolved.phase3.viral_structure_db
        if boltz_db is None:
            boltz_issues.append("viral_structure_db is not set")
        else:
            boltz_path = Path(boltz_db)
            if boltz_path.is_dir():
                boltz_issues.append("viral_structure_db must be a FoldSeek prefix")
            elif (
                not boltz_path.exists()
                and not boltz_path.with_suffix(".dbtype").exists()
            ):
                boltz_issues.append(f"FoldSeek database not found: {boltz_path}")
        if resolve_boltz_executable() is None:
            boltz_issues.append("Boltz executable not found")
        if shutil.which("foldseek") is None:
            boltz_issues.append("FoldSeek executable not found")
        if not resolved.phase3.boltz_use_msa_server:
            boltz_issues.append("boltz_use_msa_server is false")
    if boltz_issues:
        reason = "; ".join(boltz_issues)
        _warn_optional_resource(f"Boltz disabled: {reason}.")
        resolved.phase3.use_boltz = False
    else:
        reason = None
    states["boltz"] = FeatureResolution(
        requested=boltz_requested,
        required=False,
        enabled=resolved.phase3.use_boltz,
        reason_code="unavailable" if boltz_issues else None,
        details=tuple(boltz_issues),
    )

    tmvec_requested = resolved.phase3.use_tmvec_database
    tmvec_issues = _tmvec_runtime_issues(resolved) if tmvec_requested else []
    tmvec_reason = "; ".join(tmvec_issues) or None
    if tmvec_issues and resolved.phase3.tmvec_require_gpu:
        raise click.ClickException("TMVec strict requirements not met: " + tmvec_reason)
    if tmvec_issues:
        _warn_optional_resource(f"TMVec disabled: {tmvec_reason}.")
        resolved.phase3.use_tmvec_database = False
    states["tmvec"] = FeatureResolution(
        requested=tmvec_requested,
        required=resolved.phase3.tmvec_require_gpu,
        enabled=resolved.phase3.use_tmvec_database,
        reason_code="requirements_unmet" if tmvec_issues else None,
        details=tuple(tmvec_issues),
    )

    interpro_requested = resolved.phase3.interproscan_enabled
    interpro_reason = None
    if interpro_requested and not ViroSyncDatabaseManager.interproscan_available(
        resolved.phase3.interproscan_dir
    ):
        interpro_reason = (
            f"InterProScan not available at {resolved.phase3.interproscan_dir}"
        )
        _warn_optional_resource(f"InterProScan disabled: {interpro_reason}.")
        resolved.phase3.interproscan_enabled = False
    states["interproscan"] = FeatureResolution(
        requested=interpro_requested,
        required=False,
        enabled=resolved.phase3.interproscan_enabled,
        reason_code="unavailable" if interpro_reason else None,
        details=(interpro_reason,) if interpro_reason else (),
    )
    return resolved, states


def _coerce_application_config(
    value: ApplicationConfig | dict,
) -> ApplicationConfig:
    """Accept the public type plus the pre-v1 test/helper mapping shape."""
    if isinstance(value, ApplicationConfig):
        return value
    if not isinstance(value, dict):
        raise ConfigError("Application configuration must be a mapping")
    if "schema_version" in value or "orchestration" in value:
        data = dict(value)
        data.setdefault("schema_version", 1)
        return ApplicationConfig.from_dict(data)
    pipeline_sections = {
        key: item
        for key, item in value.items()
        if key
        in {
            "ablation",
            "databases",
            "compute",
            "host",
            "phase1",
            "phase2",
            "phase3",
            "execution",
        }
    }
    legacy_orchestration = {
        key: item for key, item in value.items() if key not in pipeline_sections
    }
    return ApplicationConfig.from_dict(
        {
            "schema_version": 1,
            "orchestration": legacy_orchestration,
            **pipeline_sections,
        }
    )


def _build_pipeline_config(
    yaml_config: ApplicationConfig | dict,
    clean_run: bool,
    hmm_db: Optional[Path] = None,
    hmm_allowlist: Optional[Path] = None,
    marker_faa_db: Optional[Path] = None,
    marker_faa_dir: Optional[Path] = None,
    marker_db: Optional[Path] = None,
    faa_dir: Optional[Path] = None,
    gvclass_db: Optional[Path] = None,
    gvclass_path: Optional[Path] = None,
    diamond_db: Optional[Path] = None,
    threads: Optional[int] = None,
    max_threads: Optional[int] = None,
    device: Optional[str] = None,
    search_backend: Optional[str] = None,
    max_concurrent_genomes: Optional[int] = None,
    assembly_mode: Optional[str] = None,
    hmm_chunk_size: Optional[int] = None,
    rebuild_db: Optional[bool] = None,
    phase1_initial_window_bp: Optional[int] = None,
    phase1_initial_window_genes: Optional[int] = None,
    phase1_min_markers_initial: Optional[int] = None,
    phase1_extension_kb: Optional[int] = None,
    phase1_merge_distance: Optional[int] = None,
    frameshift_screening_enabled: Optional[bool] = None,
    enable_phylogenetic: Optional[bool] = None,
    skip_masking: Optional[bool] = None,
    skip_structural: Optional[bool] = None,
    boltz: Optional[bool] = None,
    tmvec: Optional[bool] = None,
    tmvec_gpu: Optional[bool] = None,
    interproscan: Optional[bool] = None,
    high_tier_threshold: Optional[float] = None,
    low_tier_threshold: Optional[float] = None,
    use_taxonomy_ml: Optional[bool] = None,
    taxonomy_ml_model: Optional[str] = None,
) -> PipelineConfig:
    """Apply explicit CLI overrides exactly once to a decoded pipeline."""
    application = _coerce_application_config(yaml_config)
    if tmvec is False and tmvec_gpu is True:
        raise click.UsageError("--no-tmvec cannot be combined with --tmvec-gpu")
    if tmvec_gpu is True and device not in {None, "cuda"}:
        raise click.UsageError("--device cpu cannot be combined with --tmvec-gpu")
    overrides = {
        "hmm_database": hmm_db,
        "hmm_allowlist": hmm_allowlist,
        "marker_faa_db": marker_faa_db,
        "marker_faa_dir": marker_faa_dir,
        "marker_db": marker_db,
        "faa_dir": faa_dir,
        "gvclass_db": gvclass_db,
        "gvclass_path": gvclass_path,
        "diamond_db": diamond_db,
        "threads": threads,
        "max_threads": max_threads,
        "device": device,
        "search_backend": search_backend,
        "assembly_mode": assembly_mode,
        "hmm_chunk_size": hmm_chunk_size,
        "initial_window_bp": phase1_initial_window_bp,
        "initial_window_genes": phase1_initial_window_genes,
        "min_markers_initial": phase1_min_markers_initial,
        "extension_kb": phase1_extension_kb,
        "merge_distance": phase1_merge_distance,
        "frameshift_screening_enabled": frameshift_screening_enabled,
        "skip_masking": skip_masking,
        "skip_structural": skip_structural,
        "high_tier_threshold": high_tier_threshold,
        "low_tier_threshold": low_tier_threshold,
        "boundary_taxonomy_ml_enabled": use_taxonomy_ml,
        "boundary_taxonomy_ml_model": taxonomy_ml_model,
    }
    if rebuild_db is not None:
        overrides["rebuild_db"] = rebuild_db
    if enable_phylogenetic is not None:
        overrides["enable_phylogenetic"] = enable_phylogenetic
    if gvclass_path is not None:
        overrides["run_gvclass"] = True
    if boltz is not None:
        overrides["use_boltz"] = boltz
        if boltz and skip_structural is None:
            overrides["skip_structural"] = False
    if tmvec is not None:
        overrides["use_tmvec_database"] = tmvec
    if tmvec_gpu is True:
        overrides["use_tmvec_database"] = True
        overrides["tmvec_require_gpu"] = True
        overrides["device"] = "cuda"
    elif tmvec_gpu is False:
        overrides["tmvec_require_gpu"] = False
    if interproscan is not None:
        overrides["interproscan_enabled"] = interproscan
    if clean_run:
        overrides["resume"] = False

    # Process-level concurrency is intentionally not accepted as a pipeline
    # override. Keeping the argument avoids breaking direct helper callers.
    del max_concurrent_genomes
    try:
        return application.pipeline.with_overrides(**overrides)
    except ConfigError as exc:
        raise click.ClickException(str(exc)) from exc


def _resolve_pipeline_resources(
    config: PipelineConfig,
    orchestration_config,
    config_path: Optional[Path],
) -> PipelineConfig:
    """Resolve/install core resources at the sole side-effecting boundary."""
    payload = config.to_flow_kwargs()
    payload.update(
        {
            "database_root": orchestration_config.database_root,
            "core_resources_url": orchestration_config.core_resources_url,
            "core_resources_version": orchestration_config.core_resources_version,
            "core_resources_sha256": orchestration_config.core_resources_sha256,
            "core_resources_manifest_sha256": (
                orchestration_config.core_resources_manifest_sha256
            ),
        }
    )
    resolved = ViroSyncDatabaseManager.resolve_config_paths(payload, config_path)
    overrides = {
        "hmm_database": resolved.get("hmm_database") or resolved.get("hmm_db"),
        "marker_faa_db": resolved.get("marker_faa_db"),
        "marker_db": resolved.get("marker_db"),
        "gene_taxonomy_faa_db": resolved.get("gene_taxonomy_faa_db"),
        "taxonomy_labels_file": resolved.get("taxonomy_labels_file"),
        "tmvec_database_dir": resolved.get("tmvec_database_dir"),
        "interproscan_dir": resolved.get("interproscan_dir"),
    }
    updated = config.with_overrides(**overrides)
    if resolved.get("marker_faa_db") is None:
        updated = replace(
            updated,
            databases=replace(updated.databases, marker_faa_db=None),
        )
    return updated


def _validate_runtime_config(config: PipelineConfig) -> None:
    """Fail with all resource/config errors after explicit resolution."""
    errors = config.validate()
    errors.extend(config.validate_database_paths(check_files=True))
    if config.phase1.frameshift_screening_enabled:
        missing_bath_tools = [
            name
            for name in ("bathconvert", "bathsearch")
            if shutil.which(name) is None
        ]
        if missing_bath_tools:
            errors.append(
                "phase1.frameshift_screening_enabled requires commands on PATH: "
                + ", ".join(missing_bath_tools)
                + " (see docs/FRAMESHIFT_SCREENING.md)"
            )
    if config.phase3.run_gvclass:
        if config.phase3.gvclass_path is None:
            errors.append("phase3.run_gvclass requires phase3.gvclass_path")
        elif not config.phase3.gvclass_path.exists():
            errors.append(f"GVClass path does not exist: {config.phase3.gvclass_path}")
    if errors:
        raise click.ClickException(
            "Invalid runtime configuration: " + "; ".join(errors)
        )


@orchestrate.command("setup")
@click.option(
    "--config",
    "config_path",
    type=click.Path(path_type=Path),
    default=Path("config/orchestration.yaml"),
    show_default=True,
    help="Config file to read/write setup paths",
)
@click.option(
    "--db-root",
    type=click.Path(path_type=Path),
    default=None,
    help="Root directory for ViroSync core resources",
)
@click.option(
    "--core-resource",
    type=str,
    default=None,
    help="Core resources archive path or URL",
)
@click.option(
    "--core-version",
    type=str,
    default=None,
    help="Optional core resources version string override",
)
@click.option(
    "--core-resource-sha256",
    type=str,
    default=None,
    help="Pinned SHA-256 of the complete core archive",
)
@click.option(
    "--core-manifest-sha256",
    type=str,
    default=None,
    help="Pinned SHA-256 of RESOURCE_MANIFEST.json inside the archive",
)
@click.option(
    "--tmvec-url",
    type=str,
    default=None,
    help="Optional TMVec archive URL/path",
)
@click.option(
    "--tmvec-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="TMVec target directory",
)
@click.option(
    "--interproscan-url",
    type=str,
    default=None,
    help="Optional InterProScan archive URL/path",
)
@click.option(
    "--interproscan-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="InterProScan target directory",
)
@click.option(
    "--boltz-db-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Optional FoldSeek viral structure database path for Boltz",
)
@click.option(
    "--interactive-optional/--no-interactive-optional",
    default=False,
    show_default=True,
    help="Prompt to set up optional TMVec/InterProScan and Boltz DB path",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Reinstall resources even if files already exist",
)
@click.option(
    "--write-config/--no-write-config",
    default=True,
    show_default=True,
    help="Write resolved resource paths back to config file",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="Show setup details and diagnostic logs",
)
def setup(
    config_path: Path,
    db_root: Optional[Path],
    core_resource: Optional[str],
    core_version: Optional[str],
    core_resource_sha256: Optional[str],
    core_manifest_sha256: Optional[str],
    tmvec_url: Optional[str],
    tmvec_dir: Optional[Path],
    interproscan_url: Optional[str],
    interproscan_dir: Optional[Path],
    boltz_db_dir: Optional[Path],
    interactive_optional: bool,
    force: bool,
    write_config: bool,
    verbose: bool,
):
    """Install ViroSync resources and optional TMVec/InterProScan assets."""
    verbose, quiet = _command_output_flags(verbose)
    application_config = (
        _load_config(config_path)
        if config_path.exists()
        else ApplicationConfig.from_dict({"schema_version": 1})
    )
    config_data = application_config.to_dict()
    orchestration_cfg = dict(config_data.get("orchestration", {}))
    phase3_cfg = dict(config_data.get("phase3", {}))

    env_db_root = os.environ.get("VIROSYNC_DB_ROOT")
    root = Path(
        db_root
        or env_db_root
        or orchestration_cfg.get("database_root")
        or ViroSyncDatabaseManager.default_database_path()
    )
    root = ViroSyncDatabaseManager.normalize_path(root)
    source = (
        core_resource
        or orchestration_cfg.get("core_resources_url")
        or ViroSyncDatabaseManager.DATABASE_SOURCES[0]["source"]
    )
    source_record = ViroSyncDatabaseManager._record_for_source(source)
    use_config_identity = (
        core_resource is None
        or core_resource == orchestration_cfg.get("core_resources_url")
    )
    selected_version = (
        core_version
        or (
            orchestration_cfg.get("core_resources_version")
            if use_config_identity
            else None
        )
        or (source_record or {}).get("version")
    )
    selected_archive_sha256 = (
        core_resource_sha256
        or (
            orchestration_cfg.get("core_resources_sha256")
            if use_config_identity
            else None
        )
        or (source_record or {}).get("archive_sha256")
    )
    selected_manifest_sha256 = (
        core_manifest_sha256
        or (
            orchestration_cfg.get("core_resources_manifest_sha256")
            if use_config_identity
            else None
        )
        or (source_record or {}).get("manifest_sha256")
    )
    if not quiet:
        installed_version = ViroSyncDatabaseManager.get_database_version(root)
        _print_banner(
            installed_version
            if installed_version != "unknown"
            else selected_version or "not installed"
        )
    progress = BatchProgress(1, unit="resource set") if not quiet else None

    if db_root is None and not env_db_root and not orchestration_cfg.get("database_root"):
        # Skip prompt if databases already exist at the default location
        if not ViroSyncDatabaseManager._check_missing_files(root):
            if verbose:
                click.echo(f"Using existing database at: {root}")
        else:
            click.echo("ViroSync core resources:")
            click.echo("  Download size: ~6.8 GB (compressed)")
            click.echo("  Installed size: ~15 GB")
            click.echo(f"  Default location: {root}")
            click.echo()
            user_path = click.prompt("Install location", default=str(root))
            root = ViroSyncDatabaseManager.normalize_path(user_path)
            if not click.confirm(f"Download and install to {root}?", default=True):
                raise SystemExit("Setup cancelled.")

    if progress is not None:
        progress.update("resources", 0, "starting")

    if verbose:
        click.echo(f"Installing ViroSync core resources into: {root}")
        click.echo(f"Core source: {source}")

    try:
        installed_root = ViroSyncDatabaseManager.setup_database(
            database_path=str(root),
            source=source,
            version=selected_version,
            archive_sha256=selected_archive_sha256,
            manifest_sha256=selected_manifest_sha256,
            force=force,
            full=True,
            progress_callback=(
                (
                    lambda percent, stage: progress.update(
                        "resources",
                        percent * 0.85,
                        stage,
                    )
                )
                if progress is not None
                else None
            ),
        )
    except Exception as exc:
        if progress is not None:
            progress.update("resources", 100, "failed", True)
            progress.finish(False)
        click.echo(click.style(f"Core setup failed: {exc}", fg="red"), err=True)
        raise SystemExit(1)
    if progress is not None:
        progress.update("resources", 85, "checking optional resources")

    defaults = ViroSyncDatabaseManager.default_paths(installed_root)
    interactive_optional = interactive_optional and sys.stdin.isatty()

    tmvec_target = ViroSyncDatabaseManager.normalize_path(
        tmvec_dir
        or phase3_cfg.get("tmvec_database_dir")
        or defaults["tmvec_database_dir"]
    )
    tmvec_source = tmvec_url or orchestration_cfg.get("tmvec_resources_url")
    tmvec_requested = bool(
        tmvec_url is not None
        or tmvec_dir is not None
        or tmvec_source
        or phase3_cfg.get("use_tmvec_database")
    )
    tmvec_databases = (
        phase3_cfg.get("tmvec_databases")
        or orchestration_cfg.get("tmvec_databases")
        or ["bfvd"]
    )
    if isinstance(tmvec_databases, str):
        tmvec_databases = [tmvec_databases]
    tmvec_databases = [str(db_name) for db_name in tmvec_databases]

    interpro_target = ViroSyncDatabaseManager.normalize_path(
        interproscan_dir
        or phase3_cfg.get("interproscan_dir")
        or ViroSyncDatabaseManager.default_interproscan_path(installed_root)
    )
    interpro_source = interproscan_url or orchestration_cfg.get("interproscan_resources_url")
    interpro_requested = bool(
        interproscan_url is not None
        or interproscan_dir is not None
        or interpro_source
        or phase3_cfg.get("interproscan_enabled")
    )

    boltz_db_path = None
    if boltz_db_dir or phase3_cfg.get("viral_structure_db"):
        boltz_db_path = ViroSyncDatabaseManager.normalize_path(
            boltz_db_dir or phase3_cfg.get("viral_structure_db")
        )

    if interactive_optional:
        if progress is not None and progress.is_tty:
            click.echo()
        if tmvec_url is None and tmvec_dir is None:
            tmvec_requested, tmvec_target, tmvec_source = _prompt_optional_archive_choice(
                name="TMVec",
                default_target=tmvec_target,
                default_source=tmvec_source,
            )
        if interproscan_url is None and interproscan_dir is None:
            interpro_requested, interpro_target, interpro_source = _prompt_optional_archive_choice(
                name="InterProScan",
                default_target=interpro_target,
                default_source=interpro_source,
            )
        if boltz_db_dir is None:
            enable_boltz = click.confirm(
                "Configure Boltz/FoldSeek viral structure DB path now?",
                default=bool(boltz_db_path),
            )
            if enable_boltz:
                boltz_default = boltz_db_path or (installed_root / "viral_structure_db")
                boltz_input = click.prompt(
                    "Boltz viral structure DB path",
                    default=str(boltz_default),
                )
                boltz_db_path = ViroSyncDatabaseManager.normalize_path(boltz_input)

    tmvec_required: list[str] = []
    for db_name in tmvec_databases:
        tmvec_required.extend(
            ViroSyncDatabaseManager.TMVEC_REQUIRED_FILES.get(db_name, [])
        )
    if not tmvec_required:
        tmvec_required = ViroSyncDatabaseManager.TMVEC_REQUIRED_FILES["bfvd"]

    if tmvec_requested and tmvec_source:
        tmvec_ok = ViroSyncDatabaseManager.setup_optional_archive(
            name="tmvec",
            target_path=tmvec_target,
            source=tmvec_source,
            required_files=tmvec_required,
            force=force,
            progress_callback=(
                (
                    lambda percent, stage: progress.update(
                        "resources",
                        85 + percent * 0.07,
                        f"TMVec {stage}",
                    )
                )
                if progress is not None
                else None
            ),
        )
    else:
        tmvec_ok = not ViroSyncDatabaseManager.missing_tmvec_files(
            tmvec_root=tmvec_target,
            databases=tmvec_databases,
        )

    if tmvec_ok and verbose:
        click.echo(
            click.style(
                f"TMVec ready ({','.join(tmvec_databases)}): {tmvec_target}",
                fg="green",
            )
        )
    elif not tmvec_ok and (tmvec_requested or verbose):
        if tmvec_requested:
            click.echo(
                click.style(
                    "TMVec unavailable; runtime will proceed with TMVec disabled.",
                    fg="yellow",
                ),
                err=True,
            )
        else:
            click.echo(
                click.style(
                    "TMVec setup skipped; runtime will proceed with TMVec disabled.",
                    fg="yellow",
                ),
                err=True,
            )
    if progress is not None:
        progress.update("resources", 92, "TMVec check complete")

    if interpro_requested:
        interpro_ok = ViroSyncDatabaseManager.setup_optional_archive(
            name="interproscan",
            target_path=interpro_target,
            source=interpro_source,
            required_files=ViroSyncDatabaseManager.INTERPROSCAN_REQUIRED_FILES,
            force=force,
            progress_callback=(
                (
                    lambda percent, stage: progress.update(
                        "resources",
                        92 + percent * 0.05,
                        f"InterProScan {stage}",
                    )
                )
                if progress is not None
                else None
            ),
        )
    else:
        interpro_ok = ViroSyncDatabaseManager.interproscan_available(interpro_target)

    if interpro_ok and verbose:
        click.echo(click.style(f"InterProScan ready: {interpro_target}", fg="green"))
    elif not interpro_ok and (interpro_requested or verbose):
        if interpro_requested:
            click.echo(
                click.style(
                    "InterProScan unavailable; runtime will proceed with InterProScan disabled.",
                    fg="yellow",
                ),
                err=True,
            )
        else:
            click.echo(
                click.style(
                    "InterProScan setup skipped; runtime will proceed with InterProScan disabled.",
                    fg="yellow",
                ),
                err=True,
            )
    if progress is not None:
        progress.update("resources", 97, "InterProScan check complete")

    if write_config:
        cfg = config_data
        orch = cfg.setdefault("orchestration", {})
        databases = cfg.setdefault("databases", {})
        p3 = cfg.setdefault("phase3", {})

        orch["database_root"] = str(installed_root)
        orch["core_resources_url"] = source
        orch["core_resources_version"] = selected_version
        orch["core_resources_sha256"] = selected_archive_sha256
        orch["core_resources_manifest_sha256"] = selected_manifest_sha256
        databases["hmm_database"] = str(defaults["hmm_db"])
        databases["marker_faa_db"] = (
            str(defaults["marker_faa_db"])
            if defaults["marker_faa_db"] is not None
            else None
        )
        databases["marker_db"] = str(defaults["marker_db"])
        databases["gene_taxonomy_faa_db"] = str(defaults["gene_taxonomy_faa_db"])
        databases["taxonomy_labels_file"] = str(defaults["taxonomy_labels_file"])

        p3["tmvec_database_dir"] = str(tmvec_target)
        p3["tmvec_databases"] = tmvec_databases
        p3["use_tmvec_database"] = bool(p3.get("use_tmvec_database") and tmvec_ok)

        p3["interproscan_dir"] = str(interpro_target)
        p3["interproscan_enabled"] = bool(
            p3.get("interproscan_enabled") and interpro_ok
        )

        if boltz_db_path is not None:
            p3["viral_structure_db"] = str(boltz_db_path)

        try:
            ApplicationConfig.from_dict(cfg).to_yaml(config_path)
        except ConfigError as exc:
            raise click.ClickException(str(exc)) from exc
        if verbose:
            click.echo(f"Updated config: {config_path}")

    if progress is not None:
        progress.update("resources", 100, "complete")
        progress.finish(True)
    if not quiet:
        click.echo(click.style("Setup complete.", fg="green"))


@orchestrate.group("resources")
def resources_command():
    """Verify installed core resources."""


@resources_command.command("verify")
@click.option(
    "--config",
    "config_path",
    type=click.Path(path_type=Path),
    default=Path("config/orchestration.yaml"),
    show_default=True,
    help="Configuration containing the expected resource identity",
)
@click.option(
    "--db-root",
    type=click.Path(path_type=Path),
    default=None,
    help="Stable core-resource path to verify",
)
@click.option(
    "--full",
    is_flag=True,
    default=False,
    help="Hash every payload and run DIAMOND dbinfo checks",
)
def verify_resources(config_path: Path, db_root: Optional[Path], full: bool) -> None:
    """Verify pinned metadata; use --full for payload hashes and semantic probes."""
    application_config = (
        _load_config(config_path)
        if config_path.exists()
        else ApplicationConfig.from_dict({"schema_version": 1})
    )
    orchestration_cfg = application_config.orchestration
    source_record = ViroSyncDatabaseManager._record_for_source(
        orchestration_cfg.core_resources_url
        or ViroSyncDatabaseManager.DATABASE_SOURCES[0]["source"]
    )
    root = ViroSyncDatabaseManager.normalize_path(
        db_root
        or os.environ.get("VIROSYNC_DB_ROOT")
        or orchestration_cfg.database_root
        or ViroSyncDatabaseManager.default_database_path()
    )
    try:
        result = ViroSyncDatabaseManager.verify_database(
            root,
            expected_version=(
                orchestration_cfg.core_resources_version
                or (source_record or {}).get("version")
            ),
            manifest_sha256=(
                orchestration_cfg.core_resources_manifest_sha256
                or (source_record or {}).get("manifest_sha256")
            ),
            full=full,
        )
    except Exception as exc:
        raise click.ClickException(f"Core resource verification failed: {exc}") from exc
    mode = "full" if full else "fast"
    click.echo(f"Core resources verified ({mode}): {root}")
    click.echo(f"Version: {result.version}")
    click.echo(f"Manifest SHA-256: {result.manifest_sha256}")
    click.echo(f"Authenticated payloads: {result.files_verified}")


@orchestrate.command("run")
@click.option(
    "--input", "-i", "input_path",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input: genome file (.fna/.fasta), directory with genomes, or list file",
)
@click.option(
    "--output", "-o", "output_dir",
    required=True,
    type=click.Path(path_type=Path),
    help="Base output directory (subdirs created per genome)",
)
@click.option(
    "--config",
    "config_path",
    type=click.Path(exists=True, path_type=Path),
    help="YAML config file for orchestration defaults",
)
@click.option(
    "--clean-run",
    is_flag=True,
    default=False,
    help="Ignore existing outputs and start from scratch",
)
@click.option(
    "--workers", "-w",
    default=None,
    type=click.IntRange(min=1),
    help="Number of parallel genome slots (default: config value or 4)",
)
@click.option(
    "--threads-per-worker",
    default=None,
    type=click.IntRange(min=1),
    help="Threads per genome slot (default: config value or 8)",
)
@click.option(
    "--max-concurrent-genomes",
    default=None,
    type=click.IntRange(min=1),
    help="Maximum genomes processed in parallel (default: equals workers)",
)
@click.option(
    "--hmm-db",
    type=click.Path(exists=True, path_type=Path),
    help="Path to HMM database",
)
@click.option(
    "--hmm-allowlist",
    type=click.Path(exists=True, path_type=Path),
    help="Path to HMM allowlist file",
)
@click.option(
    "--marker-faa-db",
    type=click.Path(exists=True, path_type=Path),
    help="Marker FAA database file (e.g., marker.faa). Overrides --marker-faa-dir.",
)
@click.option(
    "--marker-faa-dir",
    type=click.Path(exists=True, path_type=Path),
    help="Directory with marker FAA files",
)
@click.option(
    "--marker-db",
    type=click.Path(exists=True, path_type=Path),
    help="Marker Diamond database file (e.g., marker.dmnd). Used for Phase 1 marker validation.",
)
@click.option(
    "--faa-dir",
    type=click.Path(exists=True, path_type=Path),
    help="Directory with FAA files used to build marker-validation inputs when needed",
)
@click.option(
    "--gvclass-db",
    type=click.Path(exists=True, path_type=Path),
    help="Path to GVClass database",
)
@click.option(
    "--gvclass",
    "gvclass_path",
    type=click.Path(exists=True, path_type=Path),
    help=f"Path to GVClass installation directory (defaults to ${GVCLASS_PATH_ENV_VAR} when set)",
)
@click.option(
    "--diamond-db",
    type=click.Path(exists=True, path_type=Path),
    help="Path to Diamond database",
)
@click.option(
    "--enable-phylogenetic/--disable-phylogenetic",
    default=None,
    help="Enable/disable GVClass/Diamond phylogenetic validation",
)
@click.option(
    "--assembly-mode",
    default=None,
    type=click.Choice(["default", "fragmented", "relaxed", "strict"]),
    help="Assembly mode for HHG seeding (default: config value or default)",
)
@click.option(
    "--high-tier-threshold",
    default=None,
    type=click.FloatRange(min=0.0, max=1.0),
    help="Confidence threshold for HIGH tier (default: config value or 0.7)",
)
@click.option(
    "--low-tier-threshold",
    default=None,
    type=click.FloatRange(min=0.0, max=1.0),
    help="Confidence threshold for LOW tier (default: config value or 0.2)",
)
@click.option(
    "--hmm-chunk-size",
    default=None,
    type=click.IntRange(min=1),
    help="Chunk size (pORFs) for HMM search to reduce memory use",
)
@click.option(
    "--rebuild-db/--no-rebuild-db",
    default=None,
    help="Enable/disable rebuilding marker.faa/combined.faa from inputs",
)
@click.option(
    "--phase1-initial-window-bp",
    default=None,
    type=click.IntRange(min=1),
    help="Phase 1: initial marker clustering window (bp, default: config value or 10000)",
)
@click.option(
    "--phase1-initial-window-genes",
    default=None,
    type=click.IntRange(min=1),
    help="Phase 1: initial marker clustering window (genes, default: config value or 5)",
)
@click.option(
    "--phase1-min-markers-initial",
    default=None,
    type=click.IntRange(min=1),
    help="Phase 1: minimum markers to form initial cluster (default: config value or 2)",
)
@click.option(
    "--phase1-extension-kb",
    default=None,
    type=click.IntRange(min=0),
    help="Phase 1: extension distance from outermost markers (kb, default: config value or 5)",
)
@click.option(
    "--phase1-merge-distance",
    default=None,
    type=click.IntRange(min=0),
    help="Phase 1: max gap to merge overlapping regions (bp, default: config value or 1000)",
)
@click.option(
    "--frameshift-screening/--no-frameshift-screening",
    "frameshift_screening_enabled",
    default=None,
    help="Enable/disable frameshift-sensitive marker rescue (requires BATH)",
)
@click.option(
    "--device",
    default=None,
    type=click.Choice(["cuda", "cpu"]),
    help="Device for GPU tasks (default: config value or cpu)",
)
@click.option(
    "--search-backend",
    default=None,
    type=click.Choice(["diamond"]),
    help="Sequence search backend (Diamond only)",
)
@click.option(
    "--gpu-id",
    default=None,
    type=click.IntRange(min=0),
    help="Select GPU by ID (sets CUDA_VISIBLE_DEVICES/VIROSYNC_GPU)",
)
@click.option(
    "--skip-masking/--no-skip-masking",
    default=None,
    help="Skip repeat masking step (default: use config value)",
)
@click.option(
    "--skip-structural/--no-skip-structural",
    default=None,
    help="Skip Boltz/FoldSeek structural homology (default: use config value)",
)
@click.option(
    "--boltz/--no-boltz",
    default=None,
    help="Enable/disable Boltz + FoldSeek structural homology (default: use config value)",
)
@click.option(
    "--tmvec/--no-tmvec",
    default=None,
    help="Enable/disable TMVec structural similarity search",
)
@click.option(
    "--tmvec-gpu/--no-tmvec-gpu",
    default=None,
    help="Require/do not require GPU for TMVec",
)
@click.option(
    "--interproscan/--no-interproscan",
    default=None,
    help="Enable/disable InterProScan annotation",
)
@click.option(
    "--use-taxonomy-ml/--no-taxonomy-ml",
    "use_taxonomy_ml",
    default=None,
    help="Enable taxonomy boundary ML refinement in Phase 2",
)
@click.option(
    "--taxonomy-ml-model",
    default=None,
    type=click.Choice(["logreg", "gbdt", "xgboost"]),
    help="Taxonomy boundary ML model choice (logreg|gbdt|xgboost)",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="Show configuration and diagnostic logs",
)
def run(
    input_path: Path,
    output_dir: Path,
    workers: Optional[int],
    threads_per_worker: Optional[int],
    max_concurrent_genomes: Optional[int],
    hmm_db: Optional[Path],
    hmm_allowlist: Optional[Path],
    config_path: Optional[Path],
    clean_run: bool,
    marker_faa_db: Optional[Path],
    marker_faa_dir: Optional[Path],
    marker_db: Optional[Path],
    faa_dir: Optional[Path],
    gvclass_db: Optional[Path],
    gvclass_path: Optional[Path],
    diamond_db: Optional[Path],
    enable_phylogenetic: Optional[bool],
    assembly_mode: Optional[str],
    high_tier_threshold: Optional[float],
    low_tier_threshold: Optional[float],
    hmm_chunk_size: Optional[int],
    rebuild_db: Optional[bool],
    phase1_initial_window_bp: Optional[int],
    phase1_initial_window_genes: Optional[int],
    phase1_min_markers_initial: Optional[int],
    phase1_extension_kb: Optional[int],
    phase1_merge_distance: Optional[int],
    frameshift_screening_enabled: Optional[bool],
    device: Optional[str],
    search_backend: Optional[str],
    gpu_id: Optional[int],
    skip_masking: Optional[bool],
    skip_structural: Optional[bool],
    boltz: Optional[bool],
    tmvec: Optional[bool],
    tmvec_gpu: Optional[bool],
    interproscan: Optional[bool],
    use_taxonomy_ml: Optional[bool],
    taxonomy_ml_model: Optional[str],
    verbose: bool,
    ):
    """Run ViroSync pipeline with standard-library Python parallelization.

    Processes one or more genomes in parallel using Python worker threads.
    Each genome uses its configured tool thread budget internally. Output is
    organized into per-genome subdirectories.

    Input can be:
    - A single genome file (.fna, .fasta, .fa)
    - A directory containing genome files
    - A text file with genome paths (one per line, # for comments)

    Examples:
        # Single genome
        virosync run -i genome.fna -o results/

        # Directory of genomes
        virosync run -i genomes/ -o results/ -w 8

        # List file
        virosync run -i genomes.txt -o results/ -w 8
    """
    verbose, quiet = _command_output_flags(verbose)
    ctx = click.get_current_context()

    def explicit(name: str, value):
        if ctx.get_parameter_source(name) == ParameterSource.COMMANDLINE:
            return value
        return None

    application_config = _load_config(config_path)
    worker_override = explicit("workers", workers)
    concurrency_override = explicit("max_concurrent_genomes", max_concurrent_genomes)
    if (
        worker_override is not None
        and concurrency_override is not None
        and worker_override != concurrency_override
    ):
        raise click.UsageError(
            "--workers and --max-concurrent-genomes must match when both are set"
        )
    effective_concurrency = (
        concurrency_override
        if concurrency_override is not None
        else (
            worker_override
            if worker_override is not None
            else application_config.orchestration.max_concurrent_genomes
        )
    )

    pipeline_config = _build_pipeline_config(
        yaml_config=application_config,
        clean_run=bool(explicit("clean_run", clean_run)),
        hmm_db=hmm_db,
        hmm_allowlist=hmm_allowlist,
        marker_faa_db=marker_faa_db,
        marker_faa_dir=marker_faa_dir,
        marker_db=marker_db,
        faa_dir=faa_dir,
        gvclass_db=gvclass_db,
        gvclass_path=gvclass_path,
        diamond_db=diamond_db,
        threads=explicit("threads_per_worker", threads_per_worker),
        device=explicit("device", device),
        search_backend=explicit("search_backend", search_backend),
        assembly_mode=explicit("assembly_mode", assembly_mode),
        hmm_chunk_size=hmm_chunk_size,
        rebuild_db=explicit("rebuild_db", rebuild_db),
        phase1_initial_window_bp=explicit(
            "phase1_initial_window_bp", phase1_initial_window_bp
        ),
        phase1_initial_window_genes=explicit(
            "phase1_initial_window_genes", phase1_initial_window_genes
        ),
        phase1_min_markers_initial=explicit(
            "phase1_min_markers_initial", phase1_min_markers_initial
        ),
        phase1_extension_kb=explicit("phase1_extension_kb", phase1_extension_kb),
        phase1_merge_distance=explicit("phase1_merge_distance", phase1_merge_distance),
        frameshift_screening_enabled=explicit(
            "frameshift_screening_enabled", frameshift_screening_enabled
        ),
        enable_phylogenetic=explicit("enable_phylogenetic", enable_phylogenetic),
        skip_masking=explicit("skip_masking", skip_masking),
        skip_structural=explicit("skip_structural", skip_structural),
        boltz=explicit("boltz", boltz),
        tmvec=explicit("tmvec", tmvec),
        tmvec_gpu=explicit("tmvec_gpu", tmvec_gpu),
        interproscan=explicit("interproscan", interproscan),
        high_tier_threshold=explicit("high_tier_threshold", high_tier_threshold),
        low_tier_threshold=explicit("low_tier_threshold", low_tier_threshold),
        use_taxonomy_ml=explicit("use_taxonomy_ml", use_taxonomy_ml),
        taxonomy_ml_model=explicit("taxonomy_ml_model", taxonomy_ml_model),
    )

    effective_threads, cap_warning = _cap_threads_per_worker(
        pipeline_config.compute.threads,
        pipeline_config.compute.max_threads,
        effective_concurrency,
    )
    if cap_warning:
        click.echo(click.style(cap_warning, fg="yellow"), err=True)
        pipeline_config = pipeline_config.with_overrides(threads=effective_threads)
    effective_orchestration = replace(
        application_config.orchestration,
        max_concurrent_genomes=effective_concurrency,
        gpu_id=(
            explicit("gpu_id", gpu_id)
            if explicit("gpu_id", gpu_id) is not None
            else application_config.orchestration.gpu_id
        ),
    )
    application_config = replace(
        application_config,
        orchestration=effective_orchestration,
        pipeline=pipeline_config,
    )

    # R1 safety preflight is complete before resource resolution or output writes.
    genome_paths = _collect_genome_paths(input_path)
    if not genome_paths:
        raise click.ClickException("No genome files found in input")
    try:
        _preflight_genome_runs(genome_paths, output_dir)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    if verbose:
        click.echo(f"Run start (UTC): {timestamp}")
        click.echo(f"Pipeline: ViroSync {virosync.__version__}")
        if config_path:
            click.echo(f"Config file: {config_path}")

    selected_gpu = _apply_gpu_id_env(
        explicit("gpu_id", gpu_id),
        application_config.orchestration.gpu_id,
    )
    if selected_gpu is not None and verbose:
        click.echo(f"GPU selection: {selected_gpu}")

    try:
        pipeline_config = _resolve_pipeline_resources(
            pipeline_config,
            application_config.orchestration,
            config_path,
        )
        pipeline_config, optional_features = _resolve_optional_features(pipeline_config)
        _validate_runtime_config(pipeline_config)
    except ConfigError as exc:
        raise click.ClickException(str(exc)) from exc
    if config_path and verbose:
        click.echo("Config loaded and validated.")
    application_config = replace(
        application_config,
        pipeline=pipeline_config,
    )
    effective_payload = _effective_config_payload(
        application_config,
        optional_features,
    )

    if not quiet:
        _print_banner(
            _database_version(pipeline_config)
        )
    if verbose:
        click.echo(
            f"Processing {len(genome_paths)} genomes with "
            f"{effective_concurrency} concurrent genomes"
        )
        click.echo(f"Threads per genome: {pipeline_config.compute.threads}")
        click.echo(f"Output directory: {output_dir}")
        click.echo("Effective config (CLI overrides applied):")
        click.echo(yaml.safe_dump(effective_payload, sort_keys=False).strip())

    progress = (
        BatchProgress(len(genome_paths))
        if not verbose and not quiet
        else None
    )

    results = run_batch_python(
        genome_paths=genome_paths,
        output_base_dir=output_dir,
        config=pipeline_config,
        max_concurrent_genomes=effective_concurrency,
        retries=application_config.orchestration.retries,
        retry_delay_seconds=(application_config.orchestration.retry_delay_seconds),
        effective_config=effective_payload,
        progress=progress,
    )

    # Summary
    successful = sum(1 for r in results if r.get("success", False))
    failed_results = [r for r in results if not r.get("success", False)]
    ineligible_results = [
        result
        for result in results
        if _batch_result_status(result) == "success_with_warnings"
    ]
    total_accepted = sum(r.get("accepted", 0) for r in results)
    total_candidates = sum(r.get("predictions", 0) for r in results)
    total_time = sum(r.get("elapsed_sec", 0) for r in results)

    if not quiet:
        click.echo("")
        click.echo("=" * 50)
    if failed_results:
        heading = "Batch Processing Failed"
    elif ineligible_results:
        heading = "Batch Processing Completed with Warnings"
    else:
        heading = "Batch Processing Complete"
    if not quiet:
        click.echo(click.style(heading, bold=True))
        click.echo(f"Successful: {successful}/{len(results)} genomes")
        click.echo(
            f"Benchmark eligible: {successful - len(ineligible_results)}/"
            f"{successful} successful genomes"
        )
        click.echo(
            click.style(
                f"Total EVEs: {total_accepted} canonical "
                f"({total_candidates} candidates)",
                fg="green",
            )
        )
        click.echo(f"Total time: {total_time:.0f}s")

    # Show per-genome summary
    if not quiet:
        click.echo("")
        click.echo("Per-genome results:")
        for result in results:
            genome_id = result.get("genome_id", "unknown")
            if result.get("success", False):
                accepted = result.get("accepted", 0)
                candidates = result.get("predictions", 0)
                elapsed = result.get("elapsed_sec", 0)
                warning = ""
                if _batch_result_status(result) == "success_with_warnings":
                    warning = (
                        " [SUCCESS WITH WARNINGS: benchmark_eligible=false, "
                        "legacy_resume="
                        f"{str(result.get('legacy_resume') is True).lower()}]"
                    )
                click.echo(
                    f"  {genome_id}: {accepted} canonical EVEs "
                    f"({candidates} candidates, {elapsed:.0f}s){warning}"
                )
            else:
                error = result.get("error", "Unknown error")
                click.echo(
                    click.style(
                        f"  {genome_id}: FAILED - {error}",
                        fg="red",
                    )
                )

    summary_path = output_dir / "batch_summary.tsv"
    report_path = output_dir / "batch_report.md"
    if not quiet:
        click.echo("")
        click.echo(f"Batch summary: {summary_path}")
        click.echo(f"Batch report: {report_path}")

    if failed_results or ineligible_results:
        problems = []
        if failed_results:
            problems.append(f"{len(failed_results)}/{len(results)} genomes failed")
        if ineligible_results:
            problems.append(
                f"{len(ineligible_results)}/{successful} successful genomes are "
                "benchmark-ineligible"
            )
        raise click.ClickException(
            "; ".join(problems)
            + f". Summary: {summary_path}; report: {report_path}"
        )


@orchestrate.command("info")
def info():
    """Show orchestration system information."""
    click.echo("ViroSync Orchestration System")
    click.echo("=" * 40)
    click.echo("Backend: plain Python ThreadPoolExecutor")
    click.echo(f"ViroSync version: {virosync.__version__}")
    click.echo("")
    click.echo("Available commands:")
    click.echo("  setup           - Install core resources and optional TMVec/InterProScan")
    click.echo("  run             - Process genome(s) with local Python parallelism")
    click.echo("")
    click.echo("Input options:")
    click.echo("  --input/-i      - Single genome file, directory, or list file")
