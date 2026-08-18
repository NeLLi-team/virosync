"""Checks for the optional runtime preflight command."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np

import scripts.check_structural_runtime as structural_runtime
from scripts.check_structural_runtime import main


def test_interproscan_preflight_rejects_non_executable_script(
    tmp_path: Path,
    capsys,
) -> None:
    interproscan_dir = tmp_path / "interproscan"
    interproscan_dir.mkdir()
    script = interproscan_dir / "interproscan.sh"
    script.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    script.chmod(0o644)
    config = tmp_path / "config.yaml"
    config.write_text(
        "phase3:\n"
        f"  interproscan_dir: {interproscan_dir}\n",
        encoding="utf-8",
    )

    exit_code = main(["--config", str(config), "--require-interproscan"])

    assert exit_code == 1
    assert "missing or not executable" in capsys.readouterr().out


def test_require_all_optional_runs_each_preflight_group(
    tmp_path: Path,
    capsys,
    monkeypatch,
) -> None:
    config = tmp_path / "config.yaml"
    config.write_text("phase3: {}\n", encoding="utf-8")
    monkeypatch.setattr(
        structural_runtime,
        "_check_module",
        lambda name: (False, f"missing {name}"),
    )
    monkeypatch.setattr(
        structural_runtime,
        "_check_command",
        lambda name: (False, f"missing {name}"),
    )
    monkeypatch.setattr(
        structural_runtime,
        "_check_boltz_command",
        lambda: (False, "missing boltz"),
    )

    exit_code = main(["--config", str(config), "--require-all-optional"])
    output = capsys.readouterr().out

    assert exit_code == 1
    assert "python module lobster" in output
    assert "command foldseek" in output
    assert "phase3.interproscan_dir is not configured" in output
    assert "[SKIP] TMVec" not in output
    assert "[SKIP] Boltz" not in output
    assert "[SKIP] InterProScan" not in output


def test_tmvec_preflight_rejects_unknown_device(
    tmp_path: Path,
    capsys,
    monkeypatch,
) -> None:
    config = tmp_path / "config.yaml"
    config.write_text(
        "compute:\n"
        "  device: gpu\n"
        "phase3:\n"
        "  use_tmvec_database: true\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        structural_runtime,
        "_check_module",
        lambda name: (True, f"found {name}"),
    )

    exit_code = main(["--config", str(config), "--require-tmvec"])
    output = capsys.readouterr().out

    assert exit_code == 1
    assert "compute.device must be cpu or cuda" in output


def test_tmvec_preflight_uses_production_search_and_upstream_reference(
    tmp_path: Path,
    monkeypatch,
) -> None:
    cpu_reference = tmp_path / "cpu-reference.npy"
    cuda_reference = tmp_path / "cuda-reference.npy"
    np.save(cpu_reference, np.ones(512, dtype=np.float32))
    np.save(cuda_reference, np.ones(512, dtype=np.float32))
    manifest = {
        "smoke_query": {
            "id": "query-1",
            "sequence": "M" * 60,
            "expected_target_id": "BFVD-1",
            "expected_score": 0.9,
            "score_tolerance": 0.01,
            "reference_embeddings": {
                "cpu": {
                    "path": cpu_reference.name,
                    "atol": 1e-5,
                    "rtol": 1e-5,
                },
                "cuda": {
                    "path": cuda_reference.name,
                    "atol": 1e-5,
                    "rtol": 1e-5,
                },
            },
        }
    }
    monkeypatch.setattr(
        structural_runtime.ViroSyncDatabaseManager,
        "load_tmvec_manifest",
        lambda *args, **kwargs: manifest,
    )
    calls = []

    class Predictor:
        def embed_batch(self, sequences):
            calls.append(("embed_batch", sequences))
            return np.ones((1, 512), dtype=np.float32)

        def release(self):
            calls.append(("release",))

    class Searcher:
        def __init__(self, **kwargs):
            calls.append(("init", kwargs))
            self.predictor = Predictor()

        def search_batch(self, proteins, databases):
            calls.append(("search_batch", proteins, databases))
            return {
                "query-1": {
                    "bfvd": SimpleNamespace(target_id="BFVD-1", tm_score=0.9)
                }
            }

    monkeypatch.setattr(
        "virosync.pipeline.phase3.tmvec_database.TMVecDatabaseSearch",
        Searcher,
    )

    parity_ok, parity_message, query_ok, query_message = (
        structural_runtime._check_tmvec_real_query(
            tmp_path,
            device="cpu",
            require_gpu=False,
        )
    )

    assert parity_ok is True
    assert query_ok is True
    assert "matches upstream" in parity_message
    assert query_message == "BFVD top hit BFVD-1 at score 0.900000"
    assert ("search_batch", [("query-1", "M" * 60)], ["bfvd"]) in calls
    assert calls[-1] == ("release",)
