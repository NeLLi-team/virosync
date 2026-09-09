#!/usr/bin/env python3
"""Run ViroSync commands with deferred imports for optional runtimes.

Usage:
    virosync run -i genome.fasta -o results/ --config config/orchestration.yaml
    virosync run -i genomes/ -o results/ -w 4 --threads-per-worker 8
"""

import logging
import os
import sys
from pathlib import Path
from typing import Any

import click

from virosync import __version__

logging.basicConfig(
    level=logging.ERROR,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("virosync")

_GLOBAL_CLI_FLAGS = {"-v", "--verbose", "-q", "--quiet"}


def _configure_logging(verbose: bool) -> None:
    """Show ViroSync diagnostics without third-party debug payloads."""
    logging.getLogger().setLevel(logging.ERROR)
    logging.getLogger("virosync").setLevel(logging.DEBUG if verbose else logging.ERROR)


def print_banner(database_version: str = "not resolved") -> None:
    """Print ViroSync banner."""
    banner = f"""
╔══════════════════════════════════════════════════════════════════════╗
║                                                                      ║
║   ██╗   ██╗██╗██████╗  ██████╗ ███████╗██╗   ██╗███╗   ██╗ ██████╗   ║
║   ██║   ██║██║██╔══██╗██╔═══██╗██╔════╝╚██╗ ██╔╝████╗  ██║██╔════╝   ║
║   ██║   ██║██║██████╔╝██║   ██║███████╗ ╚████╔╝ ██╔██╗ ██║██║        ║
║   ╚██╗ ██╔╝██║██╔══██╗██║   ██║╚════██║  ╚██╔╝  ██║╚██╗██║██║        ║
║    ╚████╔╝ ██║██║  ██║╚██████╔╝███████║   ██║   ██║ ╚████║╚██████╗   ║
║     ╚═══╝  ╚═╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝   ╚═╝   ╚═╝  ╚═══╝ ╚═════╝   ║
║                                                                      ║
║   Giant Endogenous Viral Element Detection Framework                 ║
║   Software version {__version__:<50}║
║   Database version {database_version:<50}║
╚══════════════════════════════════════════════════════════════════════╝
"""
    click.echo(banner)


CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


@click.group(context_settings=CONTEXT_SETTINGS, invoke_without_command=True)
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose output")
@click.option("--quiet", "-q", is_flag=True, help="Suppress non-essential output")
@click.version_option(version=__version__, prog_name="virosync")
@click.pass_context
def cli(ctx: click.Context, verbose: bool, quiet: bool) -> None:
    """ViroSync - Detect Giant Endogenous Viral Elements in genomes."""
    ctx.ensure_object(dict)

    if verbose:
        _configure_logging(True)
        ctx.obj["verbose"] = True
    elif quiet:
        _configure_logging(False)
        ctx.obj["quiet"] = True
    else:
        _configure_logging(False)
        ctx.obj["verbose"] = False
        ctx.obj["quiet"] = False

    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@cli.command()
@click.option(
    "--config",
    "config_path",
    type=click.Path(path_type=Path),
    default=Path("config/orchestration.yaml"),
    show_default=True,
    help="Configuration used to resolve the database root",
)
def info(config_path: Path) -> None:
    """Show ViroSync configuration and system info."""
    import psutil
    import torch
    import yaml

    from virosync.utils.database_manager import ViroSyncDatabaseManager

    configured_root = None
    if config_path.is_file():
        try:
            payload = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
            configured_root = (payload.get("orchestration") or {}).get("database_root")
        except (OSError, TypeError, yaml.YAMLError):
            configured_root = None
    database_root_value = (
        os.environ.get("VIROSYNC_DB_ROOT") or configured_root or ViroSyncDatabaseManager.default_database_path()
    )
    database_root = Path(database_root_value)
    if configured_root and not database_root.is_absolute():
        database_root = config_path.parent / database_root
    database_root = ViroSyncDatabaseManager.normalize_path(database_root)
    print_banner(ViroSyncDatabaseManager.get_database_version(database_root))

    click.echo("System Information:")
    click.echo(f"  Python: {sys.version.split()[0]}")
    click.echo(f"  PyTorch: {torch.__version__}")
    click.echo(f"  CUDA available: {torch.cuda.is_available()}")

    if torch.cuda.is_available():
        click.echo(f"  CUDA version: {torch.version.cuda}")
        click.echo(f"  GPU: {torch.cuda.get_device_name(0)}")
        gpu_mem = torch.cuda.get_device_properties(0).total_memory / (1024**3)
        click.echo(f"  GPU memory: {gpu_mem:.1f} GB")

    click.echo(f"\n  CPU cores: {psutil.cpu_count()}")
    click.echo(f"  RAM: {psutil.virtual_memory().total / (1024**3):.1f} GB")


class _LazyRunCommand(click.Command):
    """Lazy proxy that exposes ``orchestrate run`` as top-level ``virosync run``."""

    def __init__(self) -> None:
        super().__init__(
            name="run",
            help="Run ViroSync on one or more genomes.",
            context_settings=CONTEXT_SETTINGS,
        )
        self._real: click.Command | None = None

    def _load(self) -> click.Command:
        if self._real is None:
            from virosync.orchestration.cli import orchestrate

            self._real = orchestrate.commands["run"]
        return self._real

    def get_params(self, ctx: click.Context) -> list[click.Parameter]:
        return self._load().get_params(ctx)

    def parse_args(self, ctx: click.Context, args: list[str]) -> list[str]:
        return self._load().parse_args(ctx, args)

    def invoke(self, ctx: click.Context) -> object:
        return self._load().invoke(ctx)

    def format_help(self, ctx: click.Context, formatter: click.HelpFormatter) -> None:
        self._load().format_help(ctx, formatter)


cli.add_command(_LazyRunCommand())


class _LazyOrchestrateGroup(click.Group):
    """Lazy-loading wrapper so top-level help skips orchestration imports."""

    # Static metadata so help renders without importing orchestration module.
    _STATIC_COMMANDS = {
        "info": "Show orchestration system information.",
        "resources": "Verify installed core resources.",
        "run": "Run ViroSync on one or more genomes.",
        "setup": "Install ViroSync resources and optional databases.",
    }

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self._real: click.Group | None = None

    def _load(self) -> click.Group:
        if self._real is None:
            from virosync.orchestration.cli import orchestrate

            self._real = orchestrate
        return self._real

    def parse_args(self, ctx: click.Context, args: list[str]) -> list[str]:
        if not args:
            click.echo(ctx.get_help())
            ctx.exit(0)
        return super().parse_args(ctx, args)

    def list_commands(self, ctx: click.Context) -> list[str]:
        if self._real is None:
            return sorted(self._STATIC_COMMANDS)
        return self._real.list_commands(ctx)

    def get_command(self, ctx: click.Context, cmd_name: str) -> click.Command | None:
        return self._load().get_command(ctx, cmd_name)

    def format_commands(self, ctx: click.Context, formatter: click.HelpFormatter) -> None:
        if self._real is not None:
            super().format_commands(ctx, formatter)
            return
        with formatter.section("Commands"):
            formatter.write_dl(sorted(self._STATIC_COMMANDS.items()))


cli.add_command(
    _LazyOrchestrateGroup(
        name="orchestrate",
        help="Orchestration commands for ViroSync pipeline execution.",
        context_settings=CONTEXT_SETTINGS,
    )
)


def _first_non_global_token(argv: list[str]) -> str | None:
    """Return first argv token after top-level global flags."""
    idx = 0
    while idx < len(argv) and argv[idx] in _GLOBAL_CLI_FLAGS:
        idx += 1
    return argv[idx] if idx < len(argv) else None


def _is_bare_run(argv: list[str]) -> bool:
    """Detect bare ``-i/-o`` invocation and treat it as ``run``."""
    first_token = _first_non_global_token(argv)
    if first_token is None or not first_token.startswith("-"):
        return False

    has_input = any(token in {"-i", "--input"} or token.startswith("--input=") for token in argv)
    has_output = any(token in {"-o", "--output"} or token.startswith("--output=") for token in argv)
    return has_input and has_output


def _inject_run(argv: list[str]) -> list[str]:
    """Insert ``run`` after global flags."""
    idx = 0
    while idx < len(argv) and argv[idx] in _GLOBAL_CLI_FLAGS:
        idx += 1
    return [*argv[:idx], "run", *argv[idx:]]


def main() -> int | None:
    """Forward bare input/output options to ``run`` and invoke the CLI."""
    argv = sys.argv[1:]
    if _is_bare_run(argv):
        sys.argv = [sys.argv[0], *_inject_run(argv)]
    return cli(prog_name="virosync", obj={})


if __name__ == "__main__":
    raise SystemExit(main())
