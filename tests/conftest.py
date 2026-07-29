"""Test configuration for local package imports.

Pytest in CI and local pixi runs should be able to import from ``src/``
without requiring an editable install first.
"""

from __future__ import annotations

import sys
from pathlib import Path


SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


import pytest


@pytest.fixture(autouse=True)
def _block_database_auto_download(monkeypatch):
    """Globally neutralize ViroSyncDatabaseManager auto-download in unit tests.

    ``ViroSyncDatabaseManager.resolve_config_paths()`` will, on a
    missing-database check, call ``setup_database()`` to download the
    ~6.8 GB ViroSync resource tarball. Local dev boxes have the DB cached
    so the call is a no-op, but on a fresh GitHub Actions runner the
    download exhausts the ~14 GB free disk and fails the entire test
    session. Unit tests for CLI/YAML precedence, GVClass flag handling,
    etc. have no business pulling down the production DB. Patch the
    dangerous methods to safe no-ops by default; tests that genuinely
    need the real implementation can override in their own monkeypatch
    or fixture.
    """
    try:
        from virosync.utils.database_manager import ViroSyncDatabaseManager
    except ImportError:
        return

    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "resolve_config_paths",
        staticmethod(lambda config, config_path=None: config),
    )
    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "setup_database",
        classmethod(lambda cls, *args, **kwargs: Path("/nonexistent/skip-download")),
    )
