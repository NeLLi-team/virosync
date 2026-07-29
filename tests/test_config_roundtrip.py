"""Config round-trip tests.

Guards against YAML keys that are silently dropped by the parser (i.e. present
in config but never read, so edits have no effect). This is exactly the failure
mode that left ``host_trim_min_overlap_score`` stuck at its default.
"""

from __future__ import annotations

from pathlib import Path

import yaml

from virosync.config import ApplicationConfig, PipelineConfig

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO_ROOT / "config" / "orchestration.yaml"


def test_default_config_loads_and_round_trips() -> None:
    """The shipped config must load and survive a config -> flow-kwargs round-trip."""
    config = PipelineConfig.from_yaml(DEFAULT_CONFIG)
    flow_kwargs = config.to_flow_kwargs()
    assert isinstance(flow_kwargs, dict) and flow_kwargs


def test_host_trim_min_overlap_score_is_not_silently_dropped(tmp_path: Path) -> None:
    """A non-default phase2.host_trim_min_overlap_score must reach the config + flow kwargs."""
    data = yaml.safe_load(DEFAULT_CONFIG.read_text())
    data.setdefault("phase2", {})["host_trim_min_overlap_score"] = 0.77

    cfg = ApplicationConfig.from_dict(data).pipeline
    assert cfg.phase2.host_trim_min_overlap_score == 0.77

    flow_kwargs = cfg.to_flow_kwargs()
    assert flow_kwargs["boundary_host_trim_min_overlap_score"] == 0.77


def test_taxonomy_ml_enabled_round_trips(tmp_path: Path) -> None:
    """taxonomy_ml_enabled must be readable from YAML (was undocumented/silent-default)."""
    data = yaml.safe_load(DEFAULT_CONFIG.read_text())
    data.setdefault("phase2", {})["taxonomy_ml_enabled"] = True

    cfg = ApplicationConfig.from_dict(data).pipeline
    assert cfg.phase2.taxonomy_ml_enabled is True
    assert cfg.to_flow_kwargs()["boundary_taxonomy_ml_enabled"] is True


def test_pipeline_to_yaml_writes_readable_application_schema(tmp_path: Path) -> None:
    path = tmp_path / "pipeline.yaml"
    original = PipelineConfig.from_dict(
        {
            "compute": {"threads": 3, "device": "cpu"},
            "host": {"label": "ARC", "prefixes": ["ARC__"]},
            "phase3": {"export_all_eve_sequences": True},
        }
    )

    original.to_yaml(path)
    loaded = PipelineConfig.from_yaml(path)

    assert loaded.compute.threads == 3
    assert loaded.host.label == "ARC"
    assert loaded.phase3.export_all_eve_sequences is True


def test_default_config_has_no_dead_phase1_hhg_block() -> None:
    """The retired phase1.hhg_seeding block must stay removed (it was never parsed)."""
    data = yaml.safe_load(DEFAULT_CONFIG.read_text())
    assert "hhg_seeding" not in data.get("phase1", {})
    assert "host_trim" not in data.get("phase2", {})
