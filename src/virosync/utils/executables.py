"""Utility helpers for locating external executables."""

import shutil
from pathlib import Path
from typing import Optional


def resolve_boltz_executable() -> Optional[str]:
    """Return the path to the boltz executable, or None if not found.

    Checks:
      1. ``boltz`` on PATH (via ``shutil.which``)
      2. A local ``scripts/boltz`` relative to the package root
    """
    on_path = shutil.which("boltz")
    if on_path is not None:
        return on_path

    local_script = Path(__file__).resolve().parents[3] / "scripts" / "boltz"
    if local_script.is_file():
        return str(local_script)

    return None
