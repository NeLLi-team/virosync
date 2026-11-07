"""A finished DIAMOND call returns at once instead of at the next watchdog poll."""

from __future__ import annotations

import os
import stat
import time
from pathlib import Path

import pytest

from virosync.pipeline import search_backend


def _install_fake_diamond(directory: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    script = directory / "diamond"
    script.write_text('#!/bin/sh\nsleep 0.3\nwhile [ "$1" != "--out" ]; do shift; done\nshift\n: > "$1"\n')
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", f"{directory}{os.pathsep}{os.environ['PATH']}")


def test_run_sequence_search_returns_when_diamond_exits_not_at_poll_interval(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _install_fake_diamond(tmp_path, monkeypatch)
    monkeypatch.setattr(search_backend, "DIAMOND_WATCHDOG_POLL_INTERVAL", 10)
    query = tmp_path / "query.faa"
    query.write_text(">q1\nMKVLAAGIVGLLLAQ\n")
    output = tmp_path / "out" / "hits.tsv"

    started = time.monotonic()
    search_backend.run_sequence_search(query, tmp_path / "db.dmnd", output, threads=1)
    elapsed = time.monotonic() - started

    assert output.exists()
    assert elapsed < 3
