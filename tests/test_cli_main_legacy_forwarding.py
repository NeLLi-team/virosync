"""Tests for bare -i/-o CLI invocation forwarding."""

import logging

from click.testing import CliRunner

from virosync.cli.main import _configure_logging, _inject_run, _is_bare_run, cli


def test_bare_run_detected_for_bare_input_output() -> None:
    args = [
        "-i", "example/test-1.fna",
        "-o", "results/example_manual",
        "--config", "config/orchestration.yaml",
    ]
    assert _is_bare_run(args)


def test_bare_run_detected_for_long_options_with_equals() -> None:
    assert _is_bare_run(["--input=example/test-1.fna", "--output=results/test"])


def test_bare_run_not_detected_when_command_present() -> None:
    args = [
        "run",
        "-i", "example/test-1.fna",
        "-o", "results/example_manual",
    ]
    assert not _is_bare_run(args)


def test_bare_run_not_detected_for_help() -> None:
    assert not _is_bare_run(["--help"])


def test_inject_run_after_global_flags() -> None:
    args = [
        "-v",
        "-i", "example/test-1.fna",
        "-o", "results/example_manual",
    ]
    assert _inject_run(args) == [
        "-v",
        "run",
        "-i",
        "example/test-1.fna",
        "-o",
        "results/example_manual",
    ]


def test_orchestrate_help_lists_run_command() -> None:
    result = CliRunner().invoke(cli, ["orchestrate", "--help"])

    assert result.exit_code == 0, result.output
    assert "run" in result.output


def test_top_level_run_without_required_options_fails() -> None:
    result = CliRunner().invoke(cli, ["run"])

    assert result.exit_code == 2
    assert "Missing option '--input'" in result.output


def test_verbose_logging_does_not_enable_third_party_debug_payloads() -> None:
    _configure_logging(True)
    try:
        assert logging.getLogger("virosync.report.generate").isEnabledFor(
            logging.DEBUG
        )
        assert not logging.getLogger("papermill").isEnabledFor(logging.DEBUG)
    finally:
        _configure_logging(False)
