#!/usr/bin/env python3
"""
ViroSync CLI - Main Entry Point.

A next-generation framework for detecting Giant Endogenous Viral Elements
(GEVEs) in eukaryotic genomes.

Usage:
    virosync run -i genome.fasta -o results/ --config config/orchestration.yaml
    virosync run -i genomes/ -o results/ -w 4 --threads-per-worker 8
"""

import logging
import sys
from typing import Optional

import click
from virosync import __version__

# Configure root logger
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("virosync")

_GLOBAL_CLI_FLAGS = {"-v", "--verbose", "-q", "--quiet"}


def print_banner():
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
║   Version {__version__:<59}║
╚══════════════════════════════════════════════════════════════════════╝
"""
    click.echo(banner)


CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


@click.group(context_settings=CONTEXT_SETTINGS, invoke_without_command=True)
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose output")
@click.option("--quiet", "-q", is_flag=True, help="Suppress non-essential output")
@click.version_option(version=__version__, prog_name="virosync")
@click.pass_context
def cli(ctx, verbose: bool, quiet: bool):
    """ViroSync - Detect Giant Endogenous Viral Elements in genomes."""
    ctx.ensure_object(dict)

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        ctx.obj["verbose"] = True
    elif quiet:
        logging.getLogger().setLevel(logging.WARNING)
        ctx.obj["quiet"] = True
    else:
        ctx.obj["verbose"] = False
        ctx.obj["quiet"] = False

    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@cli.command()
def info():
    """Show ViroSync configuration and system info."""
    import torch
    import psutil

    print_banner()

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


# ---------------------------------------------------------------------------
# Lazy-loading wrappers: defer orchestration imports until actually needed
# ---------------------------------------------------------------------------

class _LazyRunCommand(click.Command):
    """Lazy proxy that exposes ``orchestrate run`` as top-level ``virosync run``."""

    def __init__(self):
        super().__init__(
            name="run",
            help="Run ViroSync pipeline with Python parallelization.",
            context_settings=CONTEXT_SETTINGS,
        )
        self._real = None

    def _load(self):
        if self._real is None:
            from virosync.orchestration.cli import orchestrate
            self._real = orchestrate.commands["run"]

    # -- Click overrides that delegate to the real command --

    def get_params(self, ctx):
        self._load()
        return self._real.get_params(ctx)

    def parse_args(self, ctx, args):
        if not args:
            self._load()
            click.echo(ctx.get_help())
            ctx.exit(0)
        self._load()
        return self._real.parse_args(ctx, args)

    def invoke(self, ctx):
        self._load()
        return self._real.invoke(ctx)

    def format_help(self, ctx, formatter):
        self._load()
        return self._real.format_help(ctx, formatter)


cli.add_command(_LazyRunCommand())


class _LazyOrchestrateGroup(click.Group):
    """Lazy-loading wrapper so top-level help skips orchestration imports."""

    # Static metadata so help renders without importing orchestration module.
    _STATIC_COMMANDS = {
        "info": "Show orchestration system information.",
        "resources": "Verify installed core resources.",
        "setup": "Install ViroSync resources and optional databases.",
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._real = None

    def _load(self):
        if self._real is None:
            from virosync.orchestration.cli import orchestrate
            self._real = orchestrate

    def parse_args(self, ctx, args):
        if not args:
            click.echo(ctx.get_help())
            ctx.exit(0)
        return super().parse_args(ctx, args)

    def list_commands(self, ctx):
        if self._real is None:
            return sorted(self._STATIC_COMMANDS)
        return self._real.list_commands(ctx)

    def get_command(self, ctx, cmd_name):
        self._load()
        return self._real.get_command(ctx, cmd_name)

    def format_commands(self, ctx, formatter):
        if self._real is None:
            commands = [
                (name, desc) for name, desc in sorted(self._STATIC_COMMANDS.items())
            ]
            if commands:
                with formatter.section("Commands"):
                    formatter.write_dl(commands)
        else:
            super().format_commands(ctx, formatter)


cli.add_command(
    _LazyOrchestrateGroup(
        name="orchestrate",
        help="Orchestration commands for ViroSync pipeline execution.",
        context_settings=CONTEXT_SETTINGS,
    )
)


# ---------------------------------------------------------------------------
# Bare ``virosync -i … -o …`` shortcut  →  ``virosync run -i … -o …``
# ---------------------------------------------------------------------------

def _first_non_global_token(argv: list[str]) -> Optional[str]:
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

    has_input = any(token in {"-i", "--input"} for token in argv)
    has_output = any(token in {"-o", "--output"} for token in argv)
    return has_input and has_output


def _inject_run(argv: list[str]) -> list[str]:
    """Insert ``run`` after global flags."""
    idx = 0
    while idx < len(argv) and argv[idx] in _GLOBAL_CLI_FLAGS:
        idx += 1
    return [*argv[:idx], "run", *argv[idx:]]


def main():
    """Main entry point."""
    argv = sys.argv[1:]
    if _is_bare_run(argv):
        sys.argv = [sys.argv[0], *_inject_run(argv)]
    return cli(prog_name="virosync", obj={})


if __name__ == "__main__":
    sys.exit(main())
