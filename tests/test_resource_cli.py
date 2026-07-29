from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner
import yaml

from virosync.config import ApplicationConfig
from virosync.orchestration.cli import _resolve_pipeline_resources, orchestrate
from virosync.utils.database_manager import ViroSyncDatabaseManager
from virosync.utils.resource_manifest import ResourceValidationResult


MANIFEST_SHA256 = (
    "7c845e29ff44b141b946864291b61eb6eefc0c695b901ad6d7351f62988f226f"
)
ARCHIVE_SHA256 = "1e513c922fd45f9e46ab672558c136713990082d51ef5875c4d705797c5a035a"
_REAL_RESOLVE_CONFIG_PATHS = ViroSyncDatabaseManager.resolve_config_paths.__func__


def test_run_resource_resolution_honors_environment_database_root(
    tmp_path: Path,
    monkeypatch,
) -> None:
    database_root = tmp_path / "isolated-resources"
    calls: list[Path] = []

    def fake_setup(cls, **kwargs):
        calls.append(Path(kwargs["database_path"]))
        return database_root

    monkeypatch.setenv("VIROSYNC_DB_ROOT", str(database_root))
    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "setup_database",
        classmethod(fake_setup),
    )

    resolved = _REAL_RESOLVE_CONFIG_PATHS(
        ViroSyncDatabaseManager,
        {
            "database_root": None,
            "core_resources_url": "fixture.tar.gz",
            "core_resources_version": "v1.0.6",
            "core_resources_sha256": ARCHIVE_SHA256,
            "core_resources_manifest_sha256": MANIFEST_SHA256,
            "hmm_database": None,
            "marker_db": None,
            "gene_taxonomy_faa_db": None,
        }
    )

    assert calls == [database_root]
    assert resolved["database_root"] == str(database_root)
    assert resolved["hmm_db"] == str(database_root / "models/combined.hmm")
    assert resolved["marker_faa_db"] is None


def test_resource_verify_cli_is_fast_by_default(
    tmp_path: Path,
    monkeypatch,
) -> None:
    calls: list[dict] = []

    def fake_verify(cls, database_path, **kwargs):
        calls.append({"database_path": Path(database_path), **kwargs})
        return ResourceValidationResult(
            version="v1.0.6",
            manifest_sha256=MANIFEST_SHA256,
            semantic_counts={},
            files_verified=13,
            full=kwargs["full"],
        )

    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "verify_database",
        classmethod(fake_verify),
    )
    result = CliRunner().invoke(
        orchestrate,
        [
            "resources",
            "verify",
            "--config",
            "config/orchestration.yaml",
            "--db-root",
            str(tmp_path / "virosync"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert calls == [
        {
            "database_path": tmp_path / "virosync",
            "expected_version": "v1.0.6",
            "manifest_sha256": MANIFEST_SHA256,
            "full": False,
        }
    ]
    assert "Core resources verified (fast)" in result.output
    assert "Authenticated payloads: 13" in result.output


def test_resource_verify_cli_full_flag_reaches_verifier(
    tmp_path: Path,
    monkeypatch,
) -> None:
    full_values: list[bool] = []

    def fake_verify(cls, database_path, **kwargs):
        full_values.append(kwargs["full"])
        return ResourceValidationResult(
            version="v1.0.6",
            manifest_sha256=MANIFEST_SHA256,
            semantic_counts={},
            files_verified=13,
            full=True,
        )

    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "verify_database",
        classmethod(fake_verify),
    )
    result = CliRunner().invoke(
        orchestrate,
        [
            "resources",
            "verify",
            "--config",
            "config/orchestration.yaml",
            "--db-root",
            str(tmp_path / "virosync"),
            "--full",
        ],
    )

    assert result.exit_code == 0, result.output
    assert full_values == [True]
    assert "Core resources verified (full)" in result.output


def test_setup_cli_exposes_digest_controls_and_has_no_fast_activation() -> None:
    result = CliRunner().invoke(orchestrate, ["setup", "--help"])

    assert result.exit_code == 0, result.output
    assert "--core-resource-sha256" in result.output
    assert "--core-manifest-sha256" in result.output
    assert "--fast" not in result.output


def test_custom_setup_source_does_not_inherit_shipped_digest_pins(
    tmp_path: Path,
    monkeypatch,
) -> None:
    captured: list[dict] = []

    def fake_setup(cls, **kwargs):
        captured.append(kwargs)
        return Path(kwargs["database_path"])

    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "setup_database",
        classmethod(fake_setup),
    )
    result = CliRunner().invoke(
        orchestrate,
        [
            "setup",
            "--config",
            "config/orchestration.yaml",
            "--db-root",
            str(tmp_path / "resources" / "virosync"),
            "--core-resource",
            str(tmp_path / "custom.tar.gz"),
            "--no-write-config",
        ],
    )

    assert result.exit_code == 0, result.output
    assert len(captured) == 1
    assert captured[0]["source"] == str(tmp_path / "custom.tar.gz")
    assert captured[0]["version"] is None
    assert captured[0]["archive_sha256"] is None
    assert captured[0]["manifest_sha256"] is None
    assert captured[0]["full"] is True


def test_setup_writes_null_when_runtime_bundle_has_no_marker_source(
    tmp_path: Path,
    monkeypatch,
) -> None:
    database_root = tmp_path / "resources" / "virosync"
    database_root.mkdir(parents=True)
    config_path = tmp_path / "orchestration.yaml"

    def fake_setup(cls, **_kwargs):
        return database_root

    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "setup_database",
        classmethod(fake_setup),
    )
    result = CliRunner().invoke(
        orchestrate,
        [
            "setup",
            "--config",
            str(config_path),
            "--db-root",
            str(database_root),
            "--core-resource",
            str(tmp_path / "runtime.tar.gz"),
            "--core-version",
            "v1.0.6",
            "--core-resource-sha256",
            ARCHIVE_SHA256,
            "--core-manifest-sha256",
            MANIFEST_SHA256,
        ],
    )

    assert result.exit_code == 0, result.output
    payload = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    assert payload["databases"]["marker_faa_db"] is None
    assert payload["databases"]["marker_db"] == str(
        database_root / "marker" / "marker.dmnd"
    )


def test_pipeline_resolution_clears_absent_auto_marker_source(
    tmp_path: Path,
    monkeypatch,
) -> None:
    application = ApplicationConfig.from_dict(
        {
            "databases": {
                "marker_faa_db": str(tmp_path / "old-marker.faa"),
            }
        }
    )

    def fake_resolve(cls, _payload, _config_path):
        return {
            "hmm_database": str(tmp_path / "combined.hmm"),
            "marker_faa_db": None,
            "marker_db": str(tmp_path / "marker.dmnd"),
            "gene_taxonomy_faa_db": str(tmp_path / "combined-proteome.dmnd"),
        }

    monkeypatch.setattr(
        ViroSyncDatabaseManager,
        "resolve_config_paths",
        classmethod(fake_resolve),
    )
    resolved = _resolve_pipeline_resources(
        application.pipeline,
        application.orchestration,
        None,
    )

    assert resolved.databases.marker_faa_db is None
