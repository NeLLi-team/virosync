"""Strict v1 configuration contract tests."""

import hashlib
import inspect
import json
from collections import Counter
from dataclasses import fields
from pathlib import Path

import pytest
import yaml

from virosync.ablation import ABLATION_CONTRACT_SHA256, AblationID
from virosync.config import ApplicationConfig, ConfigError, PipelineConfig
from virosync.config.pipeline_config import FIELD_SPECS, _RETIRED_PIPELINE_KEYS
from virosync.orchestration._flows.single_genome import (
    orchestrator as orchestrator_module,
)
from virosync.orchestration._flows.single_genome.orchestrator import (
    _single_genome_flow_impl,
    single_genome_flow,
)


_RETIRED_FLAT_KEYS = {
    "tier1_model_path",
    "tier2_model_path",
    "composition_window_bp",
    "composition_step_bp",
    "composition_classifier_threshold",
    "composition_training_step_bp",
    "composition_expansion_step_bp",
    "composition_max_gap_without_marker_kb",
    "composition_max_region_bp",
    "boundary_refinement_mode",
    "boundary_window_bp",
    "boundary_step_bp",
    "boundary_classifier_threshold",
    "boundary_ensemble_strategy",
    "boundary_refinement_margin_kb",
    "boundary_ani_threshold",
    "boundary_min_cluster_bp",
}


def test_pipeline_default_round_trip_covers_every_flow_field() -> None:
    config = PipelineConfig.from_dict({})
    kwargs = config.to_flow_kwargs()

    assert set(kwargs) == {spec.flat for spec in FIELD_SPECS if spec.emit}
    assert kwargs["device"] == "cpu"
    assert kwargs["marker_validation_top_k"] == 10
    assert kwargs["frameshift_screening_enabled"] is False
    enabled = PipelineConfig.from_dict(
        {"phase1": {"frameshift_screening_enabled": True}}
    )
    assert enabled.phase1.frameshift_screening_enabled is True
    assert "max_concurrent_genomes" not in kwargs
    assert "use_structural_homology" not in kwargs


def test_every_emitted_pipeline_key_has_exactly_one_runner_parameter() -> None:
    emitted = {spec.flat for spec in FIELD_SPECS if spec.emit}
    runner_parameters = set(inspect.signature(_single_genome_flow_impl).parameters) - {
        "genome_path",
        "output_dir",
        "genome_id",
        "progress_callback",
    }

    assert emitted == runner_parameters


def test_progress_callback_is_forwarded_once_with_config(
    tmp_path: Path,
    monkeypatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_impl(**kwargs):
        captured.update(kwargs)
        return {"success": True}

    fake_impl.__signature__ = inspect.signature(_single_genome_flow_impl)
    monkeypatch.setattr(
        orchestrator_module,
        "_single_genome_flow_impl",
        fake_impl,
    )

    def callback(percent, stage, failed):
        return None

    result = single_genome_flow(
        genome_path=tmp_path / "input.fna",
        output_dir=tmp_path / "output",
        genome_id="input",
        config=PipelineConfig(),
        progress_callback=callback,
    )

    assert result == {"success": True}
    assert captured["progress_callback"] is callback


def test_every_pipeline_dataclass_field_has_exactly_one_field_spec() -> None:
    config = PipelineConfig()
    declared = {
        (section.name, item.name)
        for section in fields(config)
        for item in fields(getattr(config, section.name))
    }
    mapped = Counter((spec.section, spec.field) for spec in FIELD_SPECS)

    assert set(mapped) == declared
    assert all(count == 1 for count in mapped.values())


def test_retired_noop_controls_are_absent_from_every_public_schema_surface() -> None:
    config = PipelineConfig()
    declared = {
        f"{section.name}.{item.name}"
        for section in fields(config)
        for item in fields(getattr(config, section.name))
    }
    mapped = {spec.flat for spec in FIELD_SPECS}
    impl_parameters = set(inspect.signature(_single_genome_flow_impl).parameters)
    wrapper_parameters = set(inspect.signature(single_genome_flow).parameters)

    assert len(_RETIRED_PIPELINE_KEYS) == 17
    assert not (_RETIRED_PIPELINE_KEYS & declared)
    assert not (_RETIRED_FLAT_KEYS & mapped)
    assert not (_RETIRED_FLAT_KEYS & impl_parameters)
    assert not (_RETIRED_FLAT_KEYS & wrapper_parameters)

    root = Path(__file__).parents[1]
    for relative in ("config/orchestration.yaml", "config/orchestration_archaeal.yaml"):
        raw = yaml.safe_load((root / relative).read_text())
        present = {
            dotted
            for dotted in _RETIRED_PIPELINE_KEYS
            if dotted.split(".", 1)[1] in raw.get(dotted.split(".", 1)[0], {})
        }
        assert not present, f"{relative} contains retired controls: {sorted(present)}"


def test_pipeline_parser_accepts_only_canonical_sections() -> None:
    config = PipelineConfig.from_dict(
        {
            "compute": {"threads": 6, "device": "cpu"},
            "phase2": {"diamond_flank_genes": 0},
            "phase3": {"export_all_eve_sequences": True},
        }
    )

    assert config.compute.threads == 6
    assert config.phase2.diamond_flank_genes == 0
    assert config.phase3.export_all_eve_sequences is True


def test_ablation_configuration_is_closed_and_digest_bound() -> None:
    config = PipelineConfig.from_dict(
        {
            "ablation": {
                "id": "A6",
                "contract_sha256": ABLATION_CONTRACT_SHA256,
            }
        }
    )

    assert config.ablation.id is AblationID.A6
    assert config.to_flow_kwargs()["ablation_id"] == "A6"
    assert (
        config.to_flow_kwargs()["ablation_contract_sha256"]
        == ABLATION_CONTRACT_SHA256
    )

    with pytest.raises(ConfigError, match="ablation.id.*must be one of"):
        PipelineConfig.from_dict({"ablation": {"id": "A7"}})
    with pytest.raises(ConfigError, match="contract_sha256.*does not match"):
        PipelineConfig.from_dict(
            {"ablation": {"contract_sha256": "0" * 64}}
        )
    with pytest.raises(ConfigError, match="A4 cannot use.*taxonomy_ml_enabled"):
        PipelineConfig.from_dict(
            {
                "ablation": {"id": "A4"},
                "phase2": {"taxonomy_ml_enabled": True},
            }
        )


def test_retired_flat_parser_is_targeted_error() -> None:
    with pytest.raises(ConfigError, match="Flat pipeline configuration"):
        PipelineConfig._from_flat_dict({"threads": 4})


@pytest.mark.parametrize(
    ("data", "message"),
    [
        ({"compute": {"threads": True}}, "must be an integer"),
        ({"phase3": {"use_boltz": 1}}, "must be a boolean"),
        ({"phase3": {"tmvec_databases": "bfvd"}}, "must be a list"),
        ({"compute": {"device": "gpu"}}, "must be one of"),
        ({"phase2": {"taxonomy_ml_threshold": 1.1}}, "between 0.0 and 1.0"),
        ({"phase2": {"diamond_control_min_distance": -1}}, "must be >= 0"),
        ({"phase2": {"diamond_flank_genes": -1}}, "must be >= 0"),
    ],
)
def test_strict_types_and_ranges(data: dict, message: str) -> None:
    with pytest.raises(ConfigError, match=message):
        PipelineConfig.from_dict(data)


def test_unknown_key_reports_dotted_path_and_suggestion() -> None:
    with pytest.raises(
        ConfigError,
        match=r"phase3\.export_all_eve_sequence.*export_all_eve_sequences",
    ):
        PipelineConfig.from_dict({"phase3": {"export_all_eve_sequence": True}})


@pytest.mark.parametrize("dotted", sorted(_RETIRED_PIPELINE_KEYS))
def test_retired_noop_controls_have_targeted_dotted_errors(dotted: str) -> None:
    section, key = dotted.split(".", 1)

    with pytest.raises(
        ConfigError,
        match=rf"Unknown configuration key '{dotted}'.*never affected runtime",
    ):
        PipelineConfig.from_dict({section: {key: None}})


@pytest.mark.parametrize("alias", ["tier1_model_path", "tier2_model_path"])
def test_retired_tier_model_aliases_have_targeted_errors(alias: str) -> None:
    with pytest.raises(
        ConfigError,
        match=rf"orchestration\.{alias}.*retired alias",
    ):
        ApplicationConfig.from_dict({"orchestration": {alias: "legacy.joblib"}})


def test_unknown_orchestration_key_suggests_accepted_legacy_alias() -> None:
    with pytest.raises(ConfigError, match="threads_per_worker"):
        ApplicationConfig.from_dict({"orchestration": {"threads_per_workr": 4}})


@pytest.mark.parametrize(
    ("orchestration", "message"),
    [
        (
            {"core_resources_url": "https://example.invalid/custom.tar.gz"},
            "requires core_resources_version",
        ),
        (
            {
                "core_resources_url": "https://example.invalid/custom.tar.gz",
                "core_resources_version": "v1.0.6",
                "core_resources_sha256": "0" * 64,
            },
            "must be configured together",
        ),
        (
            {
                "core_resources_url": "https://example.invalid/custom.tar.gz",
                "core_resources_version": "1.0.6",
                "core_resources_sha256": "0" * 64,
                "core_resources_manifest_sha256": "1" * 64,
            },
            "form vMAJOR.MINOR.PATCH",
        ),
        (
            {
                "core_resources_url": "https://example.invalid/custom.tar.gz",
                "core_resources_version": "v1.0.6",
                "core_resources_sha256": "A" * 64,
                "core_resources_manifest_sha256": "1" * 64,
            },
            "lowercase 64-character SHA-256",
        ),
    ],
)
def test_core_resource_identity_is_strict_and_all_or_none(
    orchestration: dict,
    message: str,
) -> None:
    with pytest.raises(ConfigError, match=message):
        ApplicationConfig.from_dict({"orchestration": orchestration})


def test_application_normalizes_documented_aliases_once() -> None:
    config = ApplicationConfig.from_dict(
        {
            "schema_version": 1,
            "orchestration": {
                "threads_per_worker": 12,
                "workers": 3,
                "use_structural_homology": True,
            },
            "phase1": {
                "hhg_seeding": {
                    "marker_top_k": 7,
                    "novel_marker_min_score": 42.0,
                    "novel_marker_min_coverage": 0.7,
                    "novel_marker_require_cluster": False,
                }
            },
            "phase2": {
                "boundary_diamond": {
                    "flank_genes": 0,
                    "control_min_distance": 0,
                    "top_k": 4,
                }
            },
            "phase3": {"high_pident_euk_threshold": 80.0},
        }
    )

    assert config.pipeline.compute.threads == 12
    assert config.orchestration.max_concurrent_genomes == 3
    assert config.pipeline.phase3.use_boltz is True
    assert config.pipeline.phase1.marker_validation_top_k == 7
    assert config.pipeline.phase1.novel_marker_min_score == 42.0
    assert config.pipeline.phase1.novel_marker_min_coverage == 0.7
    assert config.pipeline.phase1.novel_marker_require_cluster is False
    assert config.pipeline.phase2.diamond_flank_genes == 0
    assert config.pipeline.phase2.diamond_control_min_distance == 0
    assert config.pipeline.phase2.diamond_top_k == 4
    assert config.pipeline.host.high_pident_threshold == 80.0


def test_unversioned_legacy_mapping_defaults_to_v1() -> None:
    config = ApplicationConfig.from_dict(
        {"orchestration": {"threads_per_worker": 5, "workers": 2}}
    )

    assert config.schema_version == 1
    assert config.pipeline.compute.threads == 5
    assert config.orchestration.max_concurrent_genomes == 2


def test_explicit_unsupported_schema_version_is_rejected() -> None:
    with pytest.raises(ConfigError, match="Unsupported schema_version"):
        ApplicationConfig.from_dict({"schema_version": 2})


def test_tmvec_gpu_requirement_cannot_be_enabled_without_tmvec() -> None:
    with pytest.raises(ConfigError, match="tmvec_require_gpu.*use_tmvec_database"):
        PipelineConfig.from_dict(
            {
                "phase3": {
                    "use_tmvec_database": False,
                    "tmvec_require_gpu": True,
                }
            }
        )


def test_alias_and_canonical_collision_is_error() -> None:
    with pytest.raises(ConfigError, match="conflicts with canonical key"):
        ApplicationConfig.from_dict(
            {
                "schema_version": 1,
                "compute": {"threads": 8},
                "orchestration": {"threads_per_worker": 4},
            }
        )


@pytest.mark.parametrize(
    ("section", "message"),
    [
        (
            {"phase1": {"hhg_seeding": {"validate_hmm_only": True}}},
            "validate_hmm_only",
        ),
        (
            {"phase2": {"boundary_diamond": {"threads": 8}}},
            "compute.gene_taxonomy_threads",
        ),
        (
            {"phase2": {"boundary_diamond": {"host_prefix": "ARC__"}}},
            "derived from host.label",
        ),
        (
            {"phase2": {"host_trim": {"use_control_baseline": True}}},
            "not implemented",
        ),
        (
            {"phase1": {"composition_clusterer": "hdbscan"}},
            "composition_clusterer",
        ),
    ],
)
def test_removed_keys_are_targeted_errors(section: dict, message: str) -> None:
    with pytest.raises(ConfigError, match=message):
        ApplicationConfig.from_dict({"schema_version": 1, **section})


def test_duplicate_yaml_key_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "duplicate.yaml"
    path.write_text("schema_version: 1\ncompute:\n  threads: 4\n  threads: 8\n")

    with pytest.raises(ConfigError, match="Duplicate YAML key 'threads'"):
        ApplicationConfig.from_yaml(path)


def test_paths_are_config_relative_and_effective_payload_keeps_nulls(
    tmp_path: Path,
) -> None:
    config_dir = tmp_path / "cfg"
    config_dir.mkdir()
    path = config_dir / "config.yaml"
    path.write_text("schema_version: 1\ndatabases:\n  hmm_database: data/markers.hmm\n")

    config = ApplicationConfig.from_yaml(path)
    payload = config.effective_payload()

    assert config.pipeline.databases.hmm_database == (
        config_dir / "data" / "markers.hmm"
    )
    assert payload["databases"]["hmm_allowlist"] is None
    assert payload["optional_features"] == {}
    digest = payload.pop("effective_config_sha256")
    canonical = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    assert digest == hashlib.sha256(canonical).hexdigest()


def test_shipped_configs_are_canonical_and_portable() -> None:
    default = ApplicationConfig.from_yaml(Path("config/orchestration.yaml"))
    archaeal = ApplicationConfig.from_yaml(Path("config/orchestration_archaeal.yaml"))

    assert default.pipeline.phase3.export_all_eve_sequences is True
    assert archaeal.pipeline.host.label == "ARC"
    assert archaeal.pipeline.host.prefixes == ["ARC__"]
    assert archaeal.pipeline.phase3.export_all_eve_sequences is True
    assert archaeal.pipeline.phase3.use_tmvec_database is False
    assert archaeal.pipeline.phase3.interproscan_enabled is False
    assert "/media/" not in str(archaeal.to_dict())


@pytest.mark.parametrize(
    "path",
    [
        Path("config/orchestration.yaml"),
        Path("config/orchestration_archaeal.yaml"),
    ],
)
def test_every_shipped_yaml_key_survives_strict_decoding(path: Path) -> None:
    declared = yaml.safe_load(path.read_text())
    effective = ApplicationConfig.from_yaml(path).to_dict()

    def assert_consumed(source, decoded, dotted="") -> None:
        if isinstance(source, dict):
            for key, value in source.items():
                child = f"{dotted}.{key}" if dotted else key
                assert key in decoded, f"dropped configuration key: {child}"
                assert_consumed(value, decoded[key], child)
        else:
            assert decoded == source, f"changed configuration value: {dotted}"

    assert_consumed(declared, effective)
