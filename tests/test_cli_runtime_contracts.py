"""Phase B: runtime-contract bug fixes.

B2  threads-per-genome is capped by max_threads // max_concurrent_genomes
B3  resource download uses bounded retries + timeouts (curl -f fails on 4xx/5xx)
"""

from __future__ import annotations

import io
from pathlib import Path
import subprocess

import pytest

from virosync.orchestration import cli as orchestration_cli
from virosync.orchestration.cli import (
    _apply_gpu_id_env,
    _cap_threads_per_worker,
    _collect_genome_paths,
)
from virosync.utils import database_manager as dm
from virosync.utils.resource_installer import (
    _download_error_detail,
    copy_or_download_archive,
)


# --- B2: thread cap by max_concurrent_genomes --------------------------------


def test_directory_input_uses_supported_files_at_top_level_only(tmp_path: Path) -> None:
    top = tmp_path / "inputs"
    nested = top / "nested"
    nested.mkdir(parents=True)
    expected = [top / "a.fa", top / "b.fasta", top / "c.fna"]
    for path in [*expected, nested / "hidden.fna", top / "ignored.txt"]:
        path.write_text(">seq\nACGT\n", encoding="utf-8")

    assert _collect_genome_paths(top) == expected


def test_list_input_keeps_paths_relative_to_working_directory(tmp_path: Path) -> None:
    list_path = tmp_path / "genomes.txt"
    list_path.write_text(
        "# comment\nrelative/input.fna\n\n/absolute/input.fna\n",
        encoding="utf-8",
    )

    assert _collect_genome_paths(list_path) == [
        Path("relative/input.fna"),
        Path("/absolute/input.fna"),
    ]


def test_gpu_id_sets_uuid_selector_and_requested_index(monkeypatch) -> None:
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    monkeypatch.delenv("VIROSYNC_GPU", raising=False)
    monkeypatch.setattr(
        orchestration_cli,
        "_gpu_uuid_from_index",
        lambda _index: "GPU-test-uuid",
    )

    detail = _apply_gpu_id_env(2, None)

    assert detail == "gpu-id=2, CUDA_VISIBLE_DEVICES=GPU-test-uuid"
    assert orchestration_cli.os.environ["CUDA_VISIBLE_DEVICES"] == "GPU-test-uuid"
    assert orchestration_cli.os.environ["VIROSYNC_GPU"] == "2"


def test_gpu_id_uses_numeric_selector_when_uuid_lookup_fails(monkeypatch) -> None:
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    monkeypatch.delenv("VIROSYNC_GPU", raising=False)
    monkeypatch.setattr(
        orchestration_cli,
        "_gpu_uuid_from_index",
        lambda _index: None,
    )

    _apply_gpu_id_env(None, 3)

    assert orchestration_cli.os.environ["CUDA_VISIBLE_DEVICES"] == "3"
    assert orchestration_cli.os.environ["VIROSYNC_GPU"] == "3"

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
    captured_kwargs: list[dict] = []

    def fake_run(command, **kwargs):
        captured.append(list(command))
        captured_kwargs.append(kwargs)
        return subprocess.CompletedProcess(args=command, returncode=1, stdout="", stderr="boom")

    monkeypatch.setattr(dm.subprocess, "run", fake_run)
    monkeypatch.setattr(
        dm.ViroSyncDatabaseManager, "_certifi_ca_bundle", classmethod(lambda cls: None)
    )

    with pytest.raises(RuntimeError) as exc_info:
        dm.ViroSyncDatabaseManager._copy_or_download_archive(
            "https://example.com/resources.tar.gz", tmp_path / "out.tar.gz"
        )

    assert "boom" in str(exc_info.value)
    wget = next(c for c in captured if c and c[0] == "wget")
    curl = next(c for c in captured if c and c[0] == "curl")
    assert "--tries=3" in wget and "--timeout=120" in wget
    assert "--quiet" not in wget
    assert "--show-progress" not in wget
    assert "--progress=bar:force:noscroll" in wget
    assert "-fL" in curl  # -f => fail on HTTP error instead of saving an error page
    assert curl[curl.index("--retry") + 1] == "3"
    assert curl[curl.index("--retry-delay") + 1] == "5"
    assert curl[curl.index("--connect-timeout") + 1] == "30"
    assert "--progress-bar" in curl and "--show-error" in curl
    assert all(kwargs["stderr"] is subprocess.PIPE for kwargs in captured_kwargs)


def test_download_archive_captures_tool_bar_and_forwards_percentages(
    tmp_path,
) -> None:
    destination = tmp_path / "resources.tar.gz"
    observed: list[float] = []

    class FakeProcess:
        def __init__(self):
            self.stderr = io.BytesIO(b"resources  10%\rresources  55%\rresources 100%\n")

        def wait(self):
            return 0

    def fake_popen(command, **_kwargs):
        destination.write_bytes(b"archive")
        return FakeProcess()

    copy_or_download_archive(
        "https://example.com/resources.tar.gz",
        destination,
        popen_factory=fake_popen,
        progress_callback=observed.append,
    )

    assert observed == [10.0, 55.0, 100.0]


def test_download_archive_progress_does_not_regress_on_curl_fallback(
    tmp_path,
) -> None:
    destination = tmp_path / "resources.tar.gz"
    observed: list[float] = []

    class FakeProcess:
        def __init__(self, stderr: bytes, returncode: int):
            self.stderr = io.BytesIO(stderr)
            self._returncode = returncode

        def wait(self):
            return self._returncode

    def fake_popen(command, **_kwargs):
        if command[0] == "wget":
            return FakeProcess(b"resources 10%\rresources 40%\nfailed\n", 1)
        destination.write_bytes(b"archive")
        return FakeProcess(b"resources 5%\rresources 30%\rresources 100%\n", 0)

    copy_or_download_archive(
        "https://example.com/resources.tar.gz",
        destination,
        popen_factory=fake_popen,
        progress_callback=observed.append,
    )

    assert observed == [10.0, 40.0, 100.0]


def test_streamed_download_terminates_child_when_progress_callback_fails(
    tmp_path,
) -> None:
    destination = tmp_path / "resources.tar.gz"

    class FakeProcess:
        def __init__(self):
            self.stderr = io.BytesIO(b"resources 10%\n")
            self.terminated = False
            self.wait_calls: list[int | None] = []

        def terminate(self):
            self.terminated = True

        def wait(self, timeout=None):
            self.wait_calls.append(timeout)
            return 0

        def kill(self):
            raise AssertionError("terminate should be sufficient")

    process = FakeProcess()

    def fail_progress(_percent: float) -> None:
        raise BrokenPipeError("closed output")

    with pytest.raises(BrokenPipeError):
        copy_or_download_archive(
            "https://example.com/resources.tar.gz",
            destination,
            popen_factory=lambda *_args, **_kwargs: process,
            progress_callback=fail_progress,
        )

    assert process.terminated
    assert process.wait_calls == [5]
    assert process.stderr.closed


def test_download_error_detail_redacts_urls_hosts_and_ips() -> None:
    detail = _download_error_detail(
        "Resolving secret.example.org (secret.example.org)... 192.0.2.1\n"
        "curl: (6) Could not resolve host: https://secret.example.org/file?token=abc\n"
    )

    assert "secret.example.org" not in detail
    assert "192.0.2.1" not in detail
    assert "token=abc" not in detail
    assert "<url>" in detail
