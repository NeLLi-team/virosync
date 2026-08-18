"""Keep the public Click surface in the command-line reference."""

from __future__ import annotations

from pathlib import Path
import re

import click
import pytest

from virosync.cli.main import cli
from virosync.orchestration.cli import orchestrate


REFERENCE_PATH = Path(__file__).parents[1] / "docs" / "reference" / "cli.md"


def _reference_section(marker: str) -> str:
    text = REFERENCE_PATH.read_text(encoding="utf-8")
    start_token = f"<!-- cli-reference:{marker} -->"
    start = text.index(start_token) + len(start_token)
    end = text.find("<!-- cli-reference:", start)
    return text[start:] if end == -1 else text[start:end]


def _option_is_present(option: str, section: str) -> bool:
    pattern = rf"(?<![A-Za-z0-9_-]){re.escape(option)}(?![A-Za-z0-9_-])"
    return re.search(pattern, section) is not None


resources = orchestrate.commands["resources"]
PUBLIC_COMMAND_SECTIONS = (
    (cli, "virosync"),
    (cli.commands["info"], "virosync-info"),
    (orchestrate.commands["run"], "virosync-run"),
    (orchestrate, "virosync-orchestrate"),
    (orchestrate.commands["setup"], "virosync-orchestrate-setup"),
    (resources, "virosync-orchestrate-resources"),
    (resources.commands["verify"], "virosync-orchestrate-resources-verify"),
    (orchestrate.commands["info"], "virosync-orchestrate-info"),
)


@pytest.mark.parametrize(("command", "marker"), PUBLIC_COMMAND_SECTIONS)
def test_each_public_click_option_is_in_its_reference_section(
    command: click.Command,
    marker: str,
) -> None:
    section = _reference_section(marker)
    missing = []
    for parameter in command.params:
        if not isinstance(parameter, click.Option) or parameter.hidden:
            continue
        for option in (*parameter.opts, *parameter.secondary_opts):
            if not _option_is_present(option, section):
                missing.append(option)

    assert not missing, f"{command.name} options absent from {marker}: {missing}"


def test_documented_command_tree_matches_click_groups() -> None:
    assert set(cli.commands) == {"info", "orchestrate", "run"}
    assert set(orchestrate.commands) == {"info", "resources", "run", "setup"}
    assert set(resources.commands) == {"verify"}
