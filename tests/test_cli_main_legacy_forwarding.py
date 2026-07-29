"""Tests for bare -i/-o CLI invocation forwarding."""

from virosync.cli.main import _inject_run, _is_bare_run


def test_bare_run_detected_for_bare_input_output() -> None:
    args = [
        "-i", "example/test-1.fna",
        "-o", "results/example_manual",
        "--config", "config/orchestration.yaml",
    ]
    assert _is_bare_run(args)


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
