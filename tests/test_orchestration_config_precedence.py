from __future__ import annotations

import inspect
import json
from pathlib import Path

import click
import numpy as np
import pytest
from click.testing import CliRunner

import virosync.orchestration.cli as orchestration_cli
from virosync.config import (
    ApplicationConfig,
    ConfigError,
    FeatureResolution,
    PipelineConfig,
)
from virosync.orchestration.cli import (
    _build_pipeline_config,
    _resolve_optional_features,
    _validate_runtime_config,
    orchestrate,
)


def _application() -> ApplicationConfig:
    return ApplicationConfig.from_dict(
        {
            "schema_version": 1,
            "orchestration": {
                "max_concurrent_genomes": 6,
                "retries": 2,
                "retry_delay_seconds": 3,
            },
            "compute": {"threads": 16, "max_threads": 48, "device": "cpu"},
            "phase1": {
                "assembly_mode": "relaxed",
                "initial_window_bp": 12000,
                "extension_kb": 9,
                "frameshift_screening_enabled": True,
            },
            "phase3": {
                "high_tier_threshold": 0.8,
                "low_tier_threshold": 0.1,
                "export_all_eve_sequences": True,
            },
        }
    )


def _disabled_feature_states() -> dict[str, FeatureResolution]:
    return {
        name: FeatureResolution(False, False, False)
        for name in ("boltz", "tmvec", "interproscan")
    }


def _stub_ready_tmvec_preflight_assets(monkeypatch) -> None:
    monkeypatch.setattr(
        orchestration_cli.ViroSyncDatabaseManager,
        "load_tmvec_manifest",
        lambda *args, **kwargs: {"model": {"family": "tmvec2"}},
    )
    monkeypatch.setattr(
        orchestration_cli.importlib.util,
        "find_spec",
        lambda name: object(),
    )


def test_build_pipeline_config_preserves_yaml_when_cli_values_absent(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        orchestration_cli.ViroSyncDatabaseManager,
        "resolve_config_paths",
        lambda *args, **kwargs: pytest.fail("pure merge reached resource resolver"),
    )

    config = _build_pipeline_config(
        yaml_config=_application(),
        clean_run=False,
    )

    assert config.compute.threads == 16
    assert config.phase1.assembly_mode.value == "relaxed"
    assert config.phase1.initial_window_bp == 12000
    assert config.phase1.extension_kb == 9
    assert config.phase1.frameshift_screening_enabled is True
    assert config.phase3.high_tier_threshold == 0.8
    assert config.phase3.export_all_eve_sequences is True
    assert config.execution.resume is True


def test_legacy_resume_is_absent_from_cli_yaml_and_flow_surfaces() -> None:
    help_result = CliRunner().invoke(orchestrate, ["run", "--help"])

    assert help_result.exit_code == 0, help_result.output
    assert "--resume-allow-legacy" not in help_result.output
    with pytest.raises(ConfigError, match="resume_allow_legacy_fingerprint"):
        ApplicationConfig.from_dict(
            {
                "schema_version": 1,
                "execution": {"resume_allow_legacy_fingerprint": True},
            }
        )

    from virosync.orchestration._flows.single_genome.orchestrator import (
        _single_genome_flow_impl,
        single_genome_flow,
    )

    assert "resume_allow_legacy_fingerprint" not in inspect.signature(
        single_genome_flow
    ).parameters
    assert "resume_allow_legacy_fingerprint" not in inspect.signature(
        _single_genome_flow_impl
    ).parameters


def test_build_pipeline_config_applies_explicit_cli_values_once() -> None:
    config = _build_pipeline_config(
        yaml_config=_application(),
        clean_run=True,
        threads=8,
        device="cuda",
        assembly_mode="strict",
        phase1_initial_window_bp=10000,
        phase1_extension_kb=0,
        frameshift_screening_enabled=False,
        high_tier_threshold=0.7,
        low_tier_threshold=0.2,
        tmvec=False,
    )

    assert config.compute.threads == 8
    assert config.compute.device.value == "cuda"
    assert config.phase1.assembly_mode.value == "strict"
    assert config.phase1.initial_window_bp == 10000
    assert config.phase1.extension_kb == 0
    assert config.phase1.frameshift_screening_enabled is False
    assert config.phase3.high_tier_threshold == 0.7
    assert config.phase3.low_tier_threshold == 0.2
    assert config.phase3.use_tmvec_database is False
    assert config.execution.resume is False
    assert not hasattr(config.compute, "max_concurrent_genomes")


def test_invalid_cross_field_cli_override_fails_before_resource_resolution(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        orchestration_cli.ViroSyncDatabaseManager,
        "resolve_config_paths",
        lambda *args, **kwargs: pytest.fail("resource resolver called"),
    )

    with pytest.raises(click.ClickException, match="high_tier_threshold"):
        _build_pipeline_config(
            yaml_config=_application(),
            clean_run=False,
            high_tier_threshold=0.05,
        )


def test_paired_boolean_overrides_work_in_both_directions() -> None:
    enabled = ApplicationConfig.from_dict(
        {
            "schema_version": 1,
            "compute": {"device": "cuda"},
            "phase1": {
                "rebuild_db": True,
                "frameshift_screening_enabled": True,
            },
            "phase3": {
                "enable_phylogenetic": True,
                "use_tmvec_database": True,
                "tmvec_require_gpu": True,
                "interproscan_enabled": True,
            },
        }
    )

    disabled = _build_pipeline_config(
        yaml_config=enabled,
        clean_run=False,
        rebuild_db=False,
        frameshift_screening_enabled=False,
        enable_phylogenetic=False,
        tmvec=False,
        tmvec_gpu=False,
        interproscan=False,
    )
    reenabled = _build_pipeline_config(
        yaml_config=ApplicationConfig.from_dict({"schema_version": 1}),
        clean_run=False,
        rebuild_db=True,
        frameshift_screening_enabled=True,
        enable_phylogenetic=True,
        tmvec=True,
        tmvec_gpu=True,
        interproscan=True,
    )

    assert disabled.phase1.rebuild_db is False
    assert disabled.phase1.frameshift_screening_enabled is False
    assert disabled.phase3.enable_phylogenetic is False
    assert disabled.phase3.use_tmvec_database is False
    assert disabled.phase3.tmvec_require_gpu is False
    assert disabled.phase3.interproscan_enabled is False
    assert reenabled.phase1.rebuild_db is True
    assert reenabled.phase1.frameshift_screening_enabled is True
    assert reenabled.phase3.enable_phylogenetic is True
    assert reenabled.phase3.use_tmvec_database is True
    assert reenabled.phase3.tmvec_require_gpu is True
    assert reenabled.phase3.interproscan_enabled is True


def test_no_tmvec_conflicts_with_required_tmvec_gpu() -> None:
    with pytest.raises(click.UsageError, match="cannot be combined"):
        _build_pipeline_config(
            yaml_config=_application(),
            clean_run=False,
            tmvec=False,
            tmvec_gpu=True,
        )


def test_explicit_cpu_conflicts_with_required_tmvec_gpu() -> None:
    with pytest.raises(click.UsageError, match="--device cpu"):
        _build_pipeline_config(
            yaml_config=_application(),
            clean_run=False,
            device="cpu",
            tmvec_gpu=True,
        )


@pytest.mark.parametrize(
    ("flags", "message"),
    [
        (["--no-tmvec", "--tmvec-gpu"], "--no-tmvec"),
        (["--device", "cpu", "--tmvec-gpu"], "--device cpu"),
    ],
)
def test_contradictory_tmvec_flags_fail_at_real_cli_boundary(
    tmp_path: Path,
    monkeypatch,
    flags: list[str],
    message: str,
) -> None:
    input_fasta = tmp_path / "input.fna"
    input_fasta.write_text(">scaffold_1\nACGT\n")
    monkeypatch.setattr(orchestration_cli, "_load_config", lambda path: _application())
    monkeypatch.setattr(
        orchestration_cli,
        "_resolve_pipeline_resources",
        lambda *args: pytest.fail("resource resolver called"),
    )

    result = CliRunner().invoke(
        orchestrate,
        ["run", "-i", str(input_fasta), "-o", str(tmp_path / "results"), *flags],
    )

    assert result.exit_code == 2
    assert message in result.output
    assert not (tmp_path / "results").exists()


def test_retired_tier_model_options_are_absent_from_cli_help() -> None:
    result = CliRunner().invoke(orchestrate, ["run", "--help"])

    assert result.exit_code == 0
    assert "--tier1-model" not in result.output
    assert "--tier2-model" not in result.output


def test_tmvec_request_fails_instead_of_disabling(monkeypatch) -> None:
    config = PipelineConfig.from_dict({"phase3": {"use_tmvec_database": True}})
    monkeypatch.setattr(
        orchestration_cli,
        "_tmvec_runtime_issues",
        lambda config: ["TMVec databases not available under /missing"],
    )

    with pytest.raises(click.ClickException, match="TMVec2 requirements not met"):
        _resolve_optional_features(config)


def test_tmvec_strict_mode_fails_instead_of_disabling(monkeypatch) -> None:
    config = PipelineConfig.from_dict(
        {
            "compute": {"device": "cuda"},
            "phase3": {
                "use_tmvec_database": True,
                "tmvec_require_gpu": True,
            },
        }
    )
    monkeypatch.setattr(
        orchestration_cli,
        "_tmvec_runtime_issues",
        lambda config: ["CUDA device is not available"],
    )

    with pytest.raises(click.ClickException, match="TMVec2 requirements not met"):
        _resolve_optional_features(config)


def test_cuda_tmvec_request_fails_when_cuda_is_unavailable(
    monkeypatch,
) -> None:
    import torch

    _stub_ready_tmvec_preflight_assets(monkeypatch)
    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    config = PipelineConfig.from_dict(
        {
            "compute": {"device": "cuda"},
            "phase3": {"use_tmvec_database": True},
        }
    )

    with pytest.raises(click.ClickException, match="selected CUDA device is not available"):
        _resolve_optional_features(config)


def test_plain_cpu_tmvec_manifest_checksum_failure_is_fatal(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        orchestration_cli,
        "_tmvec_runtime_issues",
        lambda config: ["TMVec file checksum mismatch"],
    )
    config = PipelineConfig.from_dict(
        {
            "compute": {"device": "cpu"},
            "phase3": {"use_tmvec_database": True},
        }
    )

    with pytest.raises(click.ClickException, match="checksum mismatch"):
        _resolve_optional_features(config)


def _stub_runtime_boundaries(monkeypatch, application, received) -> None:
    monkeypatch.setattr(orchestration_cli, "_load_config", lambda path: application)
    monkeypatch.setattr(
        orchestration_cli,
        "_resolve_pipeline_resources",
        lambda config, orchestration, path: config,
    )
    monkeypatch.setattr(
        orchestration_cli,
        "_resolve_optional_features",
        lambda config: (
            config,
            _disabled_feature_states(),
        ),
    )
    monkeypatch.setattr(
        orchestration_cli, "_validate_runtime_config", lambda config: None
    )

    def fake_runner(**kwargs):
        received.update(kwargs)
        return [
            {
                "success": True,
                "benchmark_eligible": True,
                "legacy_resume": False,
                "accepted": 0,
                "predictions": 0,
                "elapsed_sec": 0.0,
                "genome_id": "input",
            }
        ]

    monkeypatch.setattr(orchestration_cli, "run_batch_python", fake_runner)


@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("--max-concurrent-genomes", 7),
        ("--workers", 1),
    ],
)
def test_run_resolves_process_concurrency_once(
    tmp_path: Path,
    monkeypatch,
    option: str,
    value: int,
) -> None:
    received: dict[str, object] = {}
    input_fasta = tmp_path / "input.fna"
    input_fasta.write_text(">scaffold_1\nACGT\n")
    _stub_runtime_boundaries(monkeypatch, _application(), received)

    result = CliRunner().invoke(
        orchestrate,
        [
            "run",
            "-i",
            str(input_fasta),
            "-o",
            str(tmp_path / "results"),
            option,
            str(value),
        ],
    )

    assert result.exit_code == 0, result.output
    assert received["max_concurrent_genomes"] == value
    assert received["retries"] == 2
    assert received["retry_delay_seconds"] == 3
    effective = received["effective_config"]
    assert effective["orchestration"]["max_concurrent_genomes"] == value
    assert effective["phase3"]["export_all_eve_sequences"] is True


def test_scalar_cli_overrides_reach_the_runner_config(
    tmp_path: Path,
    monkeypatch,
) -> None:
    received: dict[str, object] = {}
    input_fasta = tmp_path / "input.fna"
    input_fasta.write_text(">scaffold_1\nACGT\n")
    _stub_runtime_boundaries(monkeypatch, _application(), received)

    result = CliRunner().invoke(
        orchestrate,
        [
            "run",
            "-i",
            str(input_fasta),
            "-o",
            str(tmp_path / "results"),
            "--hmm-chunk-size",
            "123",
            "--phase1-initial-window-genes",
            "7",
            "--phase1-min-markers-initial",
            "3",
            "--phase1-merge-distance",
            "42",
            "--search-backend",
            "diamond",
            "--use-taxonomy-ml",
            "--taxonomy-ml-model",
            "gbdt",
        ],
    )

    assert result.exit_code == 0, result.output
    config = received["config"]
    assert config.phase1.hmm_chunk_size == 123
    assert config.phase1.initial_window_genes == 7
    assert config.phase1.min_markers_initial == 3
    assert config.phase1.merge_distance == 42
    assert config.compute.search_backend.value == "diamond"
    assert config.phase2.taxonomy_ml_enabled is True
    assert config.phase2.taxonomy_ml_model == "gbdt"


@pytest.mark.parametrize(
    ("flags", "use_boltz", "skip_structural"),
    [
        (["--boltz"], True, False),
        (["--boltz", "--skip-structural"], True, True),
        (["--no-boltz", "--no-skip-structural"], False, False),
    ],
)
def test_boltz_and_structural_skip_cli_interaction(
    tmp_path: Path,
    monkeypatch,
    flags: list[str],
    use_boltz: bool,
    skip_structural: bool,
) -> None:
    received: dict[str, object] = {}
    input_fasta = tmp_path / "input.fna"
    input_fasta.write_text(">scaffold_1\nACGT\n")
    _stub_runtime_boundaries(monkeypatch, _application(), received)

    result = CliRunner().invoke(
        orchestrate,
        [
            "run",
            "-i",
            str(input_fasta),
            "-o",
            str(tmp_path / "results"),
            *flags,
        ],
    )

    assert result.exit_code == 0, result.output
    config = received["config"]
    assert config.phase3.use_boltz is use_boltz
    assert config.phase3.skip_structural is skip_structural


@pytest.mark.parametrize("root_flag", ["--quiet", "--verbose"])
def test_root_output_flags_reach_run_command(
    tmp_path: Path,
    monkeypatch,
    root_flag: str,
) -> None:
    received: dict[str, object] = {}
    input_fasta = tmp_path / "input.fna"
    input_fasta.write_text(">scaffold_1\nACGT\n")
    _stub_runtime_boundaries(monkeypatch, _application(), received)

    from virosync.cli.main import cli as root_cli

    result = CliRunner().invoke(
        root_cli,
        [
            root_flag,
            "run",
            "-i",
            str(input_fasta),
            "-o",
            str(tmp_path / "results"),
        ],
    )

    assert result.exit_code == 0, result.output
    if root_flag == "--quiet":
        assert "Batch Processing Complete" not in result.output
        assert "Effective config" not in result.output
    else:
        assert "Effective config" in result.output
        assert "Batch Processing Complete" in result.output


def test_run_rejects_conflicting_concurrency_aliases_before_resources(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = tmp_path / "input.fna"
    input_fasta.write_text(">scaffold_1\nACGT\n")
    monkeypatch.setattr(orchestration_cli, "_load_config", lambda path: _application())
    monkeypatch.setattr(
        orchestration_cli,
        "_resolve_pipeline_resources",
        lambda *args: pytest.fail("resource resolver called"),
    )

    result = CliRunner().invoke(
        orchestrate,
        [
            "run",
            "-i",
            str(input_fasta),
            "-o",
            str(tmp_path / "results"),
            "--workers",
            "2",
            "--max-concurrent-genomes",
            "3",
        ],
    )

    assert result.exit_code == 2
    assert "must match" in result.output
    assert not (tmp_path / "results").exists()


@pytest.mark.parametrize("enable", [False, True])
def test_paired_boolean_flags_cross_the_real_cli_boundary(
    tmp_path: Path,
    monkeypatch,
    enable: bool,
) -> None:
    received: dict[str, object] = {}
    input_fasta = tmp_path / "input.fna"
    input_fasta.write_text(">scaffold_1\nACGT\n")
    application = ApplicationConfig.from_dict(
        {
            "schema_version": 1,
            "compute": {"device": "cuda" if not enable else "cpu"},
            "phase1": {
                "rebuild_db": not enable,
                "frameshift_screening_enabled": not enable,
            },
            "phase3": {
                "enable_phylogenetic": not enable,
                "use_tmvec_database": not enable,
                "tmvec_require_gpu": not enable,
                "interproscan_enabled": not enable,
            },
        }
    )
    _stub_runtime_boundaries(monkeypatch, application, received)
    if enable:
        flags = [
            "--rebuild-db",
            "--frameshift-screening",
            "--enable-phylogenetic",
            "--tmvec",
            "--tmvec-gpu",
            "--interproscan",
        ]
    else:
        flags = [
            "--no-rebuild-db",
            "--no-frameshift-screening",
            "--disable-phylogenetic",
            "--no-tmvec",
            "--no-tmvec-gpu",
            "--no-interproscan",
        ]

    result = CliRunner().invoke(
        orchestrate,
        [
            "run",
            "-i",
            str(input_fasta),
            "-o",
            str(tmp_path / "results"),
            *flags,
        ],
    )

    assert result.exit_code == 0, result.output
    config = received["config"]
    assert config.phase1.rebuild_db is enable
    assert config.phase1.frameshift_screening_enabled is enable
    assert config.phase3.enable_phylogenetic is enable
    assert config.phase3.use_tmvec_database is enable
    assert config.phase3.tmvec_require_gpu is enable
    assert config.phase3.interproscan_enabled is enable
    assert config.compute.device.value == "cuda"


@pytest.mark.parametrize(
    ("option", "filename", "field", "is_dir"),
    [
        ("--hmm-db", "markers.hmm", "hmm_database", False),
        ("--hmm-allowlist", "allowlist.txt", "hmm_allowlist", False),
        ("--marker-faa-db", "marker.faa", "marker_faa_db", False),
        ("--marker-faa-dir", "marker-parts", "marker_faa_dir", True),
        ("--marker-db", "marker.dmnd", "marker_db", False),
        ("--faa-dir", "proteins", "faa_dir", True),
        ("--gvclass-db", "gvclass.dmnd", "gvclass_db", False),
        ("--diamond-db", "phylogenetic.dmnd", "diamond_db", False),
    ],
)
def test_relative_database_path_override_uses_cli_working_directory(
    tmp_path: Path,
    monkeypatch,
    option: str,
    filename: str,
    field: str,
    is_dir: bool,
) -> None:
    received: dict[str, object] = {}
    input_fasta = tmp_path / "input.fna"
    marker_path = tmp_path / filename
    input_fasta.write_text(">scaffold_1\nACGT\n")
    if is_dir:
        marker_path.mkdir()
    else:
        marker_path.write_text(">marker\nAAAA\n")
    _stub_runtime_boundaries(monkeypatch, _application(), received)
    monkeypatch.chdir(tmp_path)

    result = CliRunner().invoke(
        orchestrate,
        [
            "run",
            "-i",
            "input.fna",
            "-o",
            "results",
            option,
            filename,
        ],
    )

    assert result.exit_code == 0, result.output
    config = received["config"]
    assert getattr(config.databases, field) == marker_path.resolve()


def test_path_preflight_precedes_resource_resolution(
    tmp_path: Path,
    monkeypatch,
) -> None:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "same.fna").write_text(">a\nACGT\n")
    (inputs / "same.fa").write_text(">b\nACGT\n")
    monkeypatch.setattr(orchestration_cli, "_load_config", lambda path: _application())
    monkeypatch.setattr(
        orchestration_cli,
        "_resolve_pipeline_resources",
        lambda *args: pytest.fail("resource resolver called"),
    )

    result = CliRunner().invoke(
        orchestrate,
        ["run", "-i", str(inputs), "-o", str(tmp_path / "results")],
    )

    assert result.exit_code == 1
    assert "duplicate genome ID" in result.output
    assert not (tmp_path / "results").exists()


def test_config_success_message_and_runner_wait_for_runtime_validation(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_fasta = tmp_path / "input.fna"
    input_fasta.write_text(">scaffold_1\nACGT\n")
    config_path = tmp_path / "config.yaml"
    config_path.write_text("schema_version: 1\n")
    output_dir = tmp_path / "results"
    monkeypatch.setattr(orchestration_cli, "_load_config", lambda path: _application())
    monkeypatch.setattr(
        orchestration_cli,
        "_resolve_pipeline_resources",
        lambda config, orchestration, path: config,
    )
    monkeypatch.setattr(
        orchestration_cli,
        "_resolve_optional_features",
        lambda config: (config, {}),
    )

    def fail_validation(config):
        raise click.ClickException("injected runtime validation failure")

    monkeypatch.setattr(
        orchestration_cli,
        "_validate_runtime_config",
        fail_validation,
    )
    monkeypatch.setattr(
        orchestration_cli,
        "run_batch_python",
        lambda **kwargs: pytest.fail("runner called before validation"),
    )

    result = CliRunner().invoke(
        orchestrate,
        [
            "run",
            "-i",
            str(input_fasta),
            "-o",
            str(output_dir),
            "--config",
            str(config_path),
        ],
    )

    assert result.exit_code == 1
    assert "injected runtime validation failure" in result.output
    assert "Config loaded" not in result.output
    assert not output_dir.exists()


@pytest.mark.parametrize("invalid_kind", ["empty", "directory"])
def test_strict_tmvec_invalid_assets_fail_before_config_loaded(
    tmp_path: Path,
    monkeypatch,
    invalid_kind: str,
) -> None:
    input_fasta = tmp_path / "input.fna"
    input_fasta.write_text(">scaffold_1\nACGT\n")
    config_path = tmp_path / "config.yaml"
    config_path.write_text("schema_version: 1\n")
    tmvec_root = tmp_path / "tmvec"
    tmvec_root.mkdir()
    invalid = tmvec_root / "bfvd_embeddings.npy"
    if invalid_kind == "directory":
        invalid.mkdir()
    else:
        invalid.write_bytes(b"")
    np.save(
        tmvec_root / "bfvd_annotations.npy",
        np.array([{"id": "valid"}], dtype=object),
        allow_pickle=True,
    )
    application = ApplicationConfig.from_dict(
        {
            "schema_version": 1,
            "compute": {"device": "cuda"},
            "phase3": {
                "use_tmvec_database": True,
                "tmvec_require_gpu": True,
                "tmvec_database_dir": str(tmvec_root),
                "tmvec_databases": ["bfvd"],
            },
        }
    )
    monkeypatch.setattr(orchestration_cli, "_load_config", lambda path: application)
    monkeypatch.setattr(
        orchestration_cli,
        "_resolve_pipeline_resources",
        lambda config, orchestration, path: config,
    )
    monkeypatch.setattr(orchestration_cli.importlib.util, "find_spec", lambda name: None)

    result = CliRunner().invoke(
        orchestrate,
        [
            "run",
            "-i",
            str(input_fasta),
            "-o",
            str(tmp_path / "results"),
            "--config",
            str(config_path),
        ],
    )

    assert result.exit_code == 1
    assert "TMVec2 requirements not met" in result.output
    assert "TMVEC_MANIFEST.json" in result.output
    assert "Config loaded" not in result.output
    assert not (tmp_path / "results").exists()


@pytest.mark.parametrize(
    ("database_updates", "phase3", "message"),
    [
        ({"hmm_allowlist": "missing-allowlist.txt"}, {}, "HMM allowlist"),
        ({"marker_db": "missing-marker.dmnd"}, {}, "Marker Diamond database"),
        (
            {"gvclass_db": "missing-gvclass-db"},
            {"enable_phylogenetic": True},
            "GVClass database",
        ),
        (
            {"diamond_db": "missing-phylo.dmnd"},
            {"enable_phylogenetic": True},
            "Phylogenetic Diamond database",
        ),
    ],
)
def test_selected_configured_resource_paths_fail_before_run(
    tmp_path: Path,
    database_updates: dict[str, str],
    phase3: dict,
    message: str,
) -> None:
    hmm = tmp_path / "markers.hmm"
    marker = tmp_path / "marker.dmnd"
    gvclass = tmp_path / "gvclass-db"
    diamond = tmp_path / "phylo.dmnd"
    for path in (hmm, marker, gvclass, diamond):
        path.write_bytes(b"db")
    databases = {
        "hmm_database": str(hmm),
        "marker_db": str(marker),
        "gvclass_db": str(gvclass),
        "diamond_db": str(diamond),
        **database_updates,
    }
    config = PipelineConfig.from_dict(
        {
            "databases": databases,
            "phase1": {"frameshift_screening_enabled": False},
            "phase3": phase3,
        }
    )

    with pytest.raises(click.ClickException, match=message):
        _validate_runtime_config(config)


def test_unused_marker_build_inputs_are_not_broadly_existence_checked(
    tmp_path: Path,
) -> None:
    hmm = tmp_path / "markers.hmm"
    marker = tmp_path / "marker.dmnd"
    hmm.write_bytes(b"hmm")
    marker.write_bytes(b"db")
    config = PipelineConfig.from_dict(
        {
            "databases": {
                "hmm_database": str(hmm),
                "marker_db": str(marker),
                "marker_faa_db": str(tmp_path / "unused-missing.faa"),
                "faa_dir": str(tmp_path / "unused-missing-dir"),
            },
            "phase1": {"frameshift_screening_enabled": False},
        }
    )

    _validate_runtime_config(config)


def test_runtime_config_rejects_graphviz_render_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    hmm = tmp_path / "markers.hmm"
    marker = tmp_path / "marker.dmnd"
    hmm.write_bytes(b"hmm")
    marker.write_bytes(b"db")
    config = PipelineConfig.from_dict(
        {
            "databases": {
                "hmm_database": str(hmm),
                "marker_db": str(marker),
            },
            "phase1": {"frameshift_screening_enabled": False},
        }
    )
    monkeypatch.setattr(
        orchestration_cli,
        "graphviz_runtime_error",
        lambda: "Graphviz cannot render the required sfdp PNG report",
    )

    with pytest.raises(click.ClickException, match="Graphviz cannot render"):
        _validate_runtime_config(config)


def test_selected_marker_build_inputs_are_existence_checked(tmp_path: Path) -> None:
    hmm = tmp_path / "markers.hmm"
    faa_dir = tmp_path / "faa"
    hmm.write_bytes(b"hmm")
    faa_dir.mkdir()
    config = PipelineConfig.from_dict(
        {
            "databases": {
                "hmm_database": str(hmm),
                "marker_faa_db": str(tmp_path / "missing-marker.faa"),
                "faa_dir": str(faa_dir),
            },
            "phase1": {"frameshift_screening_enabled": False},
        }
    )

    with pytest.raises(click.ClickException, match="Marker FAA database"):
        _validate_runtime_config(config)


def test_rebuild_requires_and_validates_marker_build_inputs(tmp_path: Path) -> None:
    hmm = tmp_path / "markers.hmm"
    marker_faa = tmp_path / "marker.faa"
    faa_dir = tmp_path / "faa"
    hmm.write_bytes(b"hmm")
    marker_faa.write_text(">marker\nAAAA\n")
    faa_dir.mkdir()

    valid = PipelineConfig.from_dict(
        {
            "databases": {
                "hmm_database": str(hmm),
                "marker_db": str(tmp_path / "unused-missing-marker.dmnd"),
                "marker_faa_db": str(marker_faa),
                "faa_dir": str(faa_dir),
            },
            "phase1": {
                "rebuild_db": True,
                "frameshift_screening_enabled": False,
            },
        }
    )
    _validate_runtime_config(valid)

    missing_sources = PipelineConfig.from_dict(
        {
            "databases": {
                "hmm_database": str(hmm),
                "marker_db": str(tmp_path / "installed-marker.dmnd"),
            },
            "phase1": {
                "rebuild_db": True,
                "frameshift_screening_enabled": False,
            },
        }
    )
    with pytest.raises(click.ClickException, match=r"rebuild_db requires.*faa_dir"):
        _validate_runtime_config(missing_sources)


def test_enabled_frameshift_screening_requires_bath_tools_before_run(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    hmm = tmp_path / "markers.hmm"
    marker = tmp_path / "marker.dmnd"
    hmm.write_bytes(b"hmm")
    marker.write_bytes(b"db")
    config = PipelineConfig.from_dict(
        {
            "databases": {
                "hmm_database": str(hmm),
                "marker_db": str(marker),
            },
            "phase1": {"frameshift_screening_enabled": True},
        }
    )
    monkeypatch.setattr(
        orchestration_cli.shutil,
        "which",
        lambda name: "/mock/bin/bathconvert" if name == "bathconvert" else None,
    )

    with pytest.raises(
        click.ClickException,
        match="requires commands on PATH: bathsearch",
    ):
        _validate_runtime_config(config)


@pytest.mark.parametrize(
    ("config_path", "golden_path", "host_label", "host_prefixes"),
    [
        (
            Path("config/orchestration.yaml"),
            Path("tests/data/effective_config_default.json"),
            "EUK",
            ["EUK__", "MITO__", "PLASTID__"],
        ),
        (
            Path("config/orchestration_archaeal.yaml"),
            Path("tests/data/effective_config_archaeal.json"),
            "ARC",
            ["ARC__"],
        ),
    ],
)
def test_shipped_yaml_reaches_real_cli_runner_boundary(
    tmp_path: Path,
    monkeypatch,
    config_path: Path,
    golden_path: Path,
    host_label: str,
    host_prefixes: list[str],
) -> None:
    received: dict[str, object] = {}
    input_fasta = tmp_path / "input.fna"
    input_fasta.write_text(">scaffold_1\nACGT\n")
    monkeypatch.setattr(
        orchestration_cli,
        "_resolve_pipeline_resources",
        lambda config, orchestration, path: config,
    )
    monkeypatch.setattr(
        orchestration_cli,
        "_resolve_optional_features",
        lambda config: (
            config,
            _disabled_feature_states(),
        ),
    )
    monkeypatch.setattr(
        orchestration_cli, "_validate_runtime_config", lambda config: None
    )

    def fake_runner(**kwargs):
        received.update(kwargs)
        return [
            {
                "success": True,
                "benchmark_eligible": True,
                "legacy_resume": False,
                "accepted": 0,
                "predictions": 0,
                "elapsed_sec": 0.0,
                "genome_id": "input",
            }
        ]

    monkeypatch.setattr(orchestration_cli, "run_batch_python", fake_runner)

    result = CliRunner().invoke(
        orchestrate,
        [
            "run",
            "-i",
            str(input_fasta),
            "-o",
            str(tmp_path / "results"),
            "--config",
            str(config_path.resolve()),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Config loaded and validated." not in result.output
    assert "Software version" in result.output
    assert "Database version" in result.output
    assert isinstance(
        received["progress"],
        orchestration_cli.BatchProgress,
    )
    pipeline = received["config"]
    assert pipeline.host.label == host_label
    assert pipeline.host.prefixes == host_prefixes
    assert pipeline.phase3.export_all_eve_sequences is True
    effective = received["effective_config"]
    assert effective["host"]["label"] == host_label
    assert effective["phase3"]["export_all_eve_sequences"] is True
    assert len(effective["effective_config_sha256"]) == 64
    assert effective == json.loads(golden_path.read_text())

    received.clear()
    verbose_result = CliRunner().invoke(
        orchestrate,
        [
            "run",
            "-i",
            str(input_fasta),
            "-o",
            str(tmp_path / "verbose-results"),
            "--config",
            str(config_path.resolve()),
            "--verbose",
        ],
    )

    assert verbose_result.exit_code == 0, verbose_result.output
    assert "Config loaded and validated." in verbose_result.output
    assert "Effective config (CLI overrides applied):" in verbose_result.output
    assert received["progress"] is None


# The pipeline used to have two scientific configurations: `orchestrate run` has no
# `--config` default (unlike `orchestrate setup`), so omitting the flag silently
# selected the dataclass defaults instead of the shipped, benchmarked YAML. Four
# result-affecting settings disagreed -- min_markers_initial, host_trim_enabled,
# skip_structural, export_all_eve_sequences -- so a run without `--config` produced
# different regions with no warning. The defaults were aligned to the YAML; this
# test fails if they drift apart again.
RESULT_AFFECTING_DEFAULTS = (
    ("phase1", "frameshift_screening_enabled"),
    ("phase1", "min_markers_initial"),
    ("phase2", "host_trim_enabled"),
    ("phase3", "skip_structural"),
    ("phase3", "export_all_eve_sequences"),
)


def test_result_affecting_defaults_match_the_shipped_config():
    shipped = ApplicationConfig.from_yaml(Path("config/orchestration.yaml")).pipeline
    bare = ApplicationConfig.from_dict({"schema_version": 1}).pipeline

    mismatched = {
        f"{block}.{field}": (
            getattr(getattr(bare, block), field),
            getattr(getattr(shipped, block), field),
        )
        for block, field in RESULT_AFFECTING_DEFAULTS
        if getattr(getattr(bare, block), field) != getattr(getattr(shipped, block), field)
    }
    assert not mismatched, (
        "running without --config would produce different results than the shipped "
        f"config/orchestration.yaml; (no-config, yaml) = {mismatched}"
    )


# single_genome_flow duplicates these settings as keyword defaults, so it carried a
# third copy of the scientific configuration: aligning the dataclasses alone left the
# exported Python API on the old values. Flow argument names differ from the config
# field names, hence the explicit mapping.
FLOW_ARG_FOR_FIELD = {
    ("phase1", "frameshift_screening_enabled"): "frameshift_screening_enabled",
    ("phase1", "min_markers_initial"): "min_markers_initial",
    ("phase2", "host_trim_enabled"): "boundary_host_trim_enabled",
    ("phase3", "skip_structural"): "skip_structural",
    ("phase3", "export_all_eve_sequences"): "export_all_eve_sequences",
}


def test_flow_signature_defaults_match_the_config_defaults():
    from virosync.orchestration._flows.single_genome.orchestrator import single_genome_flow

    defaults = ApplicationConfig.from_dict({"schema_version": 1}).pipeline
    signature = inspect.signature(single_genome_flow).parameters

    mismatched = {
        arg: (signature[arg].default, getattr(getattr(defaults, block), field))
        for (block, field), arg in FLOW_ARG_FOR_FIELD.items()
        if signature[arg].default != getattr(getattr(defaults, block), field)
    }
    assert not mismatched, (
        "single_genome_flow defaults disagree with the pipeline config defaults, so "
        f"calling the flow directly would produce different results; (flow, config) = {mismatched}"
    )
