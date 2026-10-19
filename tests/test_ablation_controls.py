"""Tests for the two ablation controls.

Terminal-repeat refinement and Pfam arbitration always ran, so neither could be
measured against a run without it. Both now carry a switch that defaults to on,
which leaves every shipped configuration behaving as before.
"""

from __future__ import annotations

from virosync.config.pipeline_config import PipelineConfig


def test_both_ablation_controls_default_to_on() -> None:
    config = PipelineConfig()

    assert config.phase2.terminal_repeat_enabled is True
    assert config.phase1.pfam_arbitration_enabled is True


def test_the_shipped_configuration_leaves_both_steps_on() -> None:
    config = PipelineConfig.from_yaml("config/orchestration.yaml")

    assert config.phase2.terminal_repeat_enabled is True
    assert config.phase1.pfam_arbitration_enabled is True


def test_each_control_can_be_turned_off_on_its_own() -> None:
    config = PipelineConfig()
    config.phase2.terminal_repeat_enabled = False

    assert config.phase2.terminal_repeat_enabled is False
    assert config.phase1.pfam_arbitration_enabled is True
