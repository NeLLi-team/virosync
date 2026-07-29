from __future__ import annotations

import subprocess
from pathlib import Path

from virosync.utils.database_manager import ViroSyncDatabaseManager
from virosync.utils.ssl_env import clear_stale_ssl_env_vars


def test_clear_stale_ssl_env_vars_removes_missing_paths(tmp_path: Path) -> None:
    existing_cert = tmp_path / "cert.pem"
    existing_cert.write_text("dummy")
    env = {
        "SSL_CERT_FILE": str(tmp_path / "missing.pem"),
        "CURL_CA_BUNDLE": str(existing_cert),
        "REQUESTS_CA_BUNDLE": str(tmp_path / "missing-requests.pem"),
    }

    removed = clear_stale_ssl_env_vars(env)

    assert removed == {
        "SSL_CERT_FILE": str(tmp_path / "missing.pem"),
        "REQUESTS_CA_BUNDLE": str(tmp_path / "missing-requests.pem"),
    }
    assert "SSL_CERT_FILE" not in env
    assert "REQUESTS_CA_BUNDLE" not in env
    assert env["CURL_CA_BUNDLE"] == str(existing_cert)


def test_download_uses_certifi_bundle_and_sanitized_env(
    tmp_path: Path, monkeypatch
) -> None:
    archive_path = tmp_path / "archive.tar.gz"
    ca_bundle = tmp_path / "cacert.pem"
    ca_bundle.write_text("ca")
    monkeypatch.setenv("SSL_CERT_FILE", str(tmp_path / "stale-cert.pem"))

    captured_calls: list[tuple[list[str], dict[str, str]]] = []

    def fake_run(command, env=None, stdout=None, stderr=None):
        assert stdout == subprocess.DEVNULL
        captured_calls.append((command, env))
        output_flag = "-O" if "-O" in command else "-o"
        output_path = Path(command[command.index(output_flag) + 1])
        output_path.write_bytes(b"archive")
        return subprocess.CompletedProcess(args=command, returncode=0)

    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "_certifi_ca_bundle",
        staticmethod(lambda: str(ca_bundle)),
    )
    monkeypatch.setattr("virosync.utils.database_manager.subprocess.run", fake_run)

    ViroSyncDatabaseManager._copy_or_download_archive(
        "https://example.com/archive.tar.gz",
        archive_path,
    )

    assert archive_path.exists()
    assert len(captured_calls) == 1
    command, env = captured_calls[0]
    assert "--ca-certificate" in command
    assert str(ca_bundle) in command
    assert "SSL_CERT_FILE" not in env
