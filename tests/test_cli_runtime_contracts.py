"""Phase B: runtime-contract bug fixes.

B1  _resolve_extension_kb returns the real default (5) instead of None
B2  threads-per-genome is capped by max_threads // max_concurrent_genomes
B3  resource download uses bounded retries + timeouts (curl -f fails on 4xx/5xx)
"""

from __future__ import annotations

import subprocess

import pytest

from virosync.orchestration.cli import _cap_threads_per_worker, _resolve_extension_kb
from virosync.utils import database_manager as dm


# --- B1: extension_kb None-default -------------------------------------------

@pytest.mark.parametrize(
    "cli_value,phase1_config,expected",
    [
        (None, {}, 5),  # the bug: returned None -> TypeError on later multiply
        (None, {"extension_bp": 8000}, 8),
        (None, {"extension_kb": 12}, 12),
        (None, {"extension_bp": 3000, "extension_kb": 99}, 3),  # extension_bp wins
        (7, {}, 7),  # explicit CLI value wins
        (7, {"extension_kb": 12}, 7),
        (None, {"extension_kb": None}, 5),  # explicit null -> default, not None
        (None, {"extension_bp": None}, 5),  # explicit null -> default, not TypeError
        (None, {"extension_bp": None, "extension_kb": 12}, 12),  # null bp falls through
    ],
)
def test_resolve_extension_kb(cli_value, phase1_config, expected) -> None:
    assert _resolve_extension_kb(cli_value, phase1_config) == expected


# --- B2: thread cap by max_concurrent_genomes --------------------------------

def test_cap_threads_reduces_when_oversubscribed() -> None:
    # 48 threads across 6 concurrent genomes -> 8 each; 16 requested -> reduced + warned
    capped, warning = _cap_threads_per_worker(16, max_threads=48, max_concurrent_genomes=6)
    assert capped == 8
    assert warning is not None and "8" in warning


def test_cap_threads_keeps_fitting_value() -> None:
    # recommended single-genome config: budget 48 >= 32 requested -> unchanged
    capped, warning = _cap_threads_per_worker(32, max_threads=48, max_concurrent_genomes=1)
    assert capped == 32
    assert warning is None


def test_cap_threads_noop_without_max_threads() -> None:
    capped, warning = _cap_threads_per_worker(16, max_threads=None, max_concurrent_genomes=6)
    assert capped == 16
    assert warning is None


def test_cap_threads_floor_is_one() -> None:
    capped, _ = _cap_threads_per_worker(16, max_threads=4, max_concurrent_genomes=100)
    assert capped == 1


# --- B3: hardened download ---------------------------------------------------

def test_download_archive_uses_hardened_flags(tmp_path, monkeypatch) -> None:
    captured: list[list[str]] = []

    def fake_run(command, **kwargs):
        captured.append(list(command))
        return subprocess.CompletedProcess(args=command, returncode=1, stdout="", stderr="boom")

    monkeypatch.setattr(dm.subprocess, "run", fake_run)
    monkeypatch.setattr(
        dm.ViroSyncDatabaseManager, "_certifi_ca_bundle", classmethod(lambda cls: None)
    )

    with pytest.raises(RuntimeError):
        dm.ViroSyncDatabaseManager._copy_or_download_archive(
            "https://example.com/resources.tar.gz", tmp_path / "out.tar.gz"
        )

    wget = next(c for c in captured if c and c[0] == "wget")
    curl = next(c for c in captured if c and c[0] == "curl")
    assert "--tries=3" in wget and "--timeout=120" in wget
    assert "-fL" in curl  # -f => fail on HTTP error instead of saving an error page
    assert curl[curl.index("--retry") + 1] == "3"
    assert curl[curl.index("--retry-delay") + 1] == "5"
    assert curl[curl.index("--connect-timeout") + 1] == "30"
