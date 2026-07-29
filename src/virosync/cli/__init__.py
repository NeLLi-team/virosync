"""
ViroSync CLI Module.

Provides command-line interface for running the ViroSync pipeline.
"""


def __getattr__(name: str):
    if name in {"cli", "main"}:
        from virosync.cli import main as cli_main
        return getattr(cli_main, name)
    raise AttributeError(f"module 'virosync.cli' has no attribute {name!r}")
