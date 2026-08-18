"""Regression tests for GVClass flag/path resolution in orchestration CLI."""

from pathlib import Path

import click
import pytest

from virosync.orchestration.cli import (
    _build_pipeline_config,
    _validate_runtime_config,
    run,
)


def test_gvclass_enabled_in_yaml_requires_path() -> None:
    config = _build_pipeline_config(
        yaml_config={"phase3": {"run_gvclass": True}},
        clean_run=True,
    )

    with pytest.raises(click.ClickException, match="run_gvclass requires"):
        _validate_runtime_config(config)


def test_gvclass_cli_path_enables_run(tmp_path: Path) -> None:
    gvclass_root = tmp_path / "gvclass"
    gvclass_root.mkdir()

    cfg = _build_pipeline_config(
        yaml_config={},
        clean_run=True,
        gvclass_path=gvclass_root,
    )

    assert cfg.phase3.run_gvclass is True
    assert cfg.phase3.gvclass_path == gvclass_root


def test_gvclass_yaml_path_enables_run(tmp_path: Path) -> None:
    gvclass_root = tmp_path / "gvclass"
    gvclass_root.mkdir()
    yaml_config = {
        "phase3": {
            "run_gvclass": True,
            "gvclass_path": str(gvclass_root),
        }
    }

    cfg = _build_pipeline_config(yaml_config=yaml_config, clean_run=True)

    assert cfg.phase3.run_gvclass is True
    assert cfg.phase3.gvclass_path == gvclass_root


def test_gvclass_option_reads_legacy_environment_path() -> None:
    option = next(param for param in run.params if param.name == "gvclass_path")
    assert option.envvar == "VIROSYNC_GVCLASS_PATH"


def test_gvclass_validation_requires_executable(tmp_path: Path) -> None:
    gvclass_root = tmp_path / "gvclass"
    gvclass_root.mkdir()
    config = _build_pipeline_config(
        yaml_config={},
        clean_run=True,
        gvclass_path=gvclass_root,
    )

    with pytest.raises(click.ClickException, match="missing or not executable"):
        _validate_runtime_config(config)

    executable = gvclass_root / "gvclass"
    executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    executable.chmod(0o755)
    with pytest.raises(click.ClickException) as exc_info:
        _validate_runtime_config(config)
    assert "GVClass executable" not in str(exc_info.value)
