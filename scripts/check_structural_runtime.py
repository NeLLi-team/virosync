#!/usr/bin/env python3
"""Preflight checks for optional structural runtime (TMVec + Boltz).

This is intentionally operator-facing: it verifies the exact Python modules,
executables, and configured resource paths used by the optional structural
layers before a run is started with `--tmvec` and/or `--boltz`.
"""

from __future__ import annotations

import argparse
import importlib
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

from virosync.utils.database_manager import ViroSyncDatabaseManager
from virosync.utils.executables import resolve_boltz_executable


def _load_phase3(config_path: Path) -> dict:
    data = yaml.safe_load(config_path.read_text()) or {}
    if "phase3" in data and isinstance(data["phase3"], dict):
        return dict(data["phase3"])
    return {}


def _load_compute(config_path: Path) -> dict:
    data = yaml.safe_load(config_path.read_text()) or {}
    if "compute" in data and isinstance(data["compute"], dict):
        return dict(data["compute"])
    return {}


def _check_module(module_name: str) -> tuple[bool, str]:
    try:
        mod = importlib.import_module(module_name)
    except Exception as exc:
        return False, str(exc)

    version = getattr(mod, "__version__", None)
    if version is None:
        return True, "installed"
    return True, str(version)


def _check_command(name: str) -> tuple[bool, str]:
    path = shutil.which(name)
    if not path:
        return False, "not found on PATH"
    return True, path


def _check_boltz_command() -> tuple[bool, str]:
    boltz_executable = resolve_boltz_executable()
    if boltz_executable is None:
        return False, "not found on PATH and scripts/boltz missing"

    try:
        proc = subprocess.run(
            [boltz_executable, "--help"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=20,
            check=False,
        )
    except Exception as exc:
        return False, f"failed to execute: {exc}"

    if proc.returncode != 0:
        return False, f"boltz --help failed (exit={proc.returncode})"
    return True, f"ok ({boltz_executable})"


def _check_tmvec_real_query(
    tmvec_dir: str | Path,
    *,
    device: str,
    require_gpu: bool,
) -> tuple[bool, str, bool, str]:
    searcher = None
    predictor = None
    try:
        from virosync.pipeline.phase3.tmvec_database import TMVecDatabaseSearch

        root = Path(tmvec_dir)
        manifest = ViroSyncDatabaseManager.load_tmvec_manifest(
            root,
            verify_hashes=True,
            databases=["bfvd"],
        )
        smoke = manifest["smoke_query"]
        reference_contract = smoke["reference_embeddings"][device]
        searcher = TMVecDatabaseSearch(
            database_root=root,
            databases=["bfvd"],
            device=device,
            require_gpu=require_gpu,
            fail_on_unavailable=True,
        )
        predictor = searcher.predictor
        runtime_embedding = np.asarray(
            predictor.embed_batch([smoke["sequence"]]),
            dtype=np.float32,
        )
        reference_embedding = np.asarray(
            np.load(root / reference_contract["path"], allow_pickle=False),
            dtype=np.float32,
        ).reshape(1, -1)
        if runtime_embedding.shape != reference_embedding.shape:
            return (
                False,
                f"runtime shape {runtime_embedding.shape} != reference shape "
                f"{reference_embedding.shape}",
                False,
                "not run",
            )
        if not np.allclose(
            runtime_embedding,
            reference_embedding,
            atol=float(reference_contract["atol"]),
            rtol=float(reference_contract["rtol"]),
        ):
            max_error = float(np.max(np.abs(runtime_embedding - reference_embedding)))
            return (
                False,
                f"runtime embedding differs from the upstream reference "
                f"(max_abs_error={max_error:.6g})",
                False,
                "not run",
            )

        results = searcher.search_batch(
            [(smoke["id"], smoke["sequence"])],
            databases=["bfvd"],
        )
        hit = results.get(smoke["id"], {}).get("bfvd")
        if hit is None:
            return (
                True,
                "runtime embedding matches upstream reference",
                False,
                "no BFVD hit",
            )
        score_error = abs(float(hit.tm_score) - float(smoke["expected_score"]))
        if hit.target_id != smoke["expected_target_id"]:
            return (
                True,
                "runtime embedding matches upstream reference",
                False,
                f"top hit {hit.target_id!r} != {smoke['expected_target_id']!r}",
            )
        if score_error > float(smoke["score_tolerance"]):
            return (
                True,
                "runtime embedding matches upstream reference",
                False,
                f"score {hit.tm_score:.6f} differs from expected "
                f"{float(smoke['expected_score']):.6f} by {score_error:.6f}",
            )
    except Exception as exc:
        return False, str(exc), False, "not run"
    finally:
        if predictor is not None:
            predictor.release()

    return (
        True,
        "runtime embedding matches upstream reference",
        True,
        f"BFVD top hit {hit.target_id} at score {hit.tm_score:.6f}",
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/orchestration.yaml"),
        help="Path to orchestration YAML config",
    )
    parser.add_argument(
        "--require-tmvec",
        action="store_true",
        help="Fail unless TMVec2 and the configured BFVD resources pass a real query.",
    )
    parser.add_argument(
        "--require-boltz",
        action="store_true",
        help="Fail unless the Boltz wrapper/runtime, Foldseek, and configured Foldseek DB are available.",
    )
    parser.add_argument(
        "--require-interproscan",
        action="store_true",
        help="Fail unless the configured InterProScan installation is available.",
    )
    parser.add_argument(
        "--require-all-optional",
        action="store_true",
        help="Equivalent to --require-tmvec --require-boltz --require-interproscan.",
    )
    args = parser.parse_args(argv)

    config_path = args.config.resolve()
    if not config_path.exists():
        print(f"[FAIL] config not found: {config_path}")
        return 1

    phase3 = _load_phase3(config_path)
    compute = _load_compute(config_path)
    tmvec_dir = phase3.get("tmvec_database_dir") or str(
        ViroSyncDatabaseManager.default_tmvec_path()
    )
    tmvec_databases = phase3.get("tmvec_databases") or ["bfvd"]
    tmvec_require_gpu = bool(phase3.get("tmvec_require_gpu"))
    tmvec_device = str(compute.get("device") or ("cuda" if tmvec_require_gpu else "cpu"))
    interpro_dir = phase3.get("interproscan_dir")
    viral_structure_db = phase3.get("viral_structure_db")
    boltz_use_msa_server = bool(phase3.get("boltz_use_msa_server"))

    failures = 0
    require_tmvec = args.require_all_optional or args.require_tmvec
    require_boltz = args.require_all_optional or args.require_boltz
    require_interproscan = args.require_all_optional or args.require_interproscan
    check_tmvec = require_tmvec or bool(phase3.get("use_tmvec_database"))
    check_boltz = require_boltz or bool(phase3.get("use_boltz"))
    check_interproscan = require_interproscan or bool(phase3.get("interproscan_enabled"))

    print(f"[INFO] config: {config_path}")
    print(f"[INFO] phase3.tmvec_database_dir: {tmvec_dir}")
    print(f"[INFO] phase3.tmvec_databases: {tmvec_databases}")
    print(f"[INFO] TMVec2 device: {tmvec_device}")
    print(f"[INFO] phase3.interproscan_dir: {interpro_dir}")
    print(f"[INFO] phase3.viral_structure_db: {viral_structure_db}")
    print(f"[INFO] phase3.boltz_use_msa_server: {boltz_use_msa_server}")

    if not any([check_tmvec, check_boltz, check_interproscan]):
        print("[INFO] no optional structural/domain layers are enabled or required")

    tmvec_ready_for_model_smoke = True

    if check_tmvec:
        if tmvec_device not in {"cpu", "cuda"}:
            print(
                "[FAIL] compute.device must be cpu or cuda for TMVec2; "
                f"found {tmvec_device}"
            )
            failures += 1
            tmvec_ready_for_model_smoke = False
        required_modules = [
            "torch",
            "transformers",
            "lightning",
            "lobster",
        ]
        for name in required_modules:
            ok, msg = _check_module(name)
            state = "OK" if ok else "FAIL"
            print(f"[{state}] python module {name}: {msg}")
            if not ok:
                failures += 1
                tmvec_ready_for_model_smoke = False

        if tmvec_require_gpu and tmvec_device != "cuda":
            print("[FAIL] TMVec2 GPU mode requires compute.device: cuda")
            failures += 1
            tmvec_ready_for_model_smoke = False
        if tmvec_device == "cuda":
            try:
                import torch

                cuda_ok = bool(torch.cuda.is_available())
                if cuda_ok:
                    gpu_count = torch.cuda.device_count()
                    print(f"[OK] CUDA runtime: available ({gpu_count} device(s))")
                else:
                    print("[FAIL] CUDA runtime: torch.cuda.is_available() == False")
                    failures += 1
                    tmvec_ready_for_model_smoke = False
            except Exception as exc:
                print(f"[FAIL] CUDA runtime: unable to query torch ({exc})")
                failures += 1
                tmvec_ready_for_model_smoke = False
        elif tmvec_device == "cpu":
            print("[OK] TMVec2 device: CPU")
    else:
        print("[SKIP] TMVec runtime checks (not enabled; pass --require-tmvec to force)")

    if check_boltz:
        if not boltz_use_msa_server:
            print(
                "[FAIL] phase3.boltz_use_msa_server is false; ViroSync's "
                "ViroSync writes sequence-only Boltz inputs and needs MSA-server mode."
            )
            failures += 1

        ok, msg = _check_command("foldseek")
        state = "OK" if ok else "FAIL"
        print(f"[{state}] command foldseek: {msg}")
        if not ok:
            failures += 1

        ok, msg = _check_boltz_command()
        state = "OK" if ok else "FAIL"
        print(f"[{state}] command boltz: {msg}")
        if not ok:
            failures += 1
    else:
        print("[SKIP] Boltz/Foldseek runtime checks (not enabled; pass --require-boltz to force)")

    if check_interproscan:
        if interpro_dir:
            interpro_path = Path(interpro_dir)
            interpro_exec = interpro_path / "interproscan.sh"
            if ViroSyncDatabaseManager.interproscan_available(interpro_path):
                print(f"[OK] interproscan.sh: {interpro_exec}")
            else:
                print(
                    "[FAIL] interproscan.sh missing or not executable: "
                    f"{interpro_exec}"
                )
                failures += 1
        else:
            print("[FAIL] phase3.interproscan_dir is not configured")
            failures += 1
    else:
        print("[SKIP] InterProScan checks (not enabled; pass --require-interproscan to force)")

    if check_tmvec:
        if tmvec_dir:
            missing_tmvec = ViroSyncDatabaseManager.missing_tmvec_files(
                tmvec_root=Path(tmvec_dir),
                databases=[str(db) for db in tmvec_databases],
            )
            if missing_tmvec:
                print("[FAIL] TMVec database files missing:")
                for entry in missing_tmvec:
                    print(f"  - {entry}")
                failures += 1
                tmvec_ready_for_model_smoke = False
            else:
                print("[OK] TMVec2 manifest and BFVD files: complete")
        else:
            print("[FAIL] phase3.tmvec_database_dir is not configured")
            failures += 1
            tmvec_ready_for_model_smoke = False

        if tmvec_ready_for_model_smoke:
            parity_ok, parity_msg, query_ok, query_msg = _check_tmvec_real_query(
                tmvec_dir,
                device=tmvec_device,
                require_gpu=tmvec_require_gpu,
            )
            parity_state = "OK" if parity_ok else "FAIL"
            print(f"[{parity_state}] TMVec2 upstream parity: {parity_msg}")
            query_state = "OK" if query_ok else "FAIL"
            print(f"[{query_state}] TMVec2 BFVD smoke query: {query_msg}")
            if not parity_ok or not query_ok:
                failures += 1

    if check_boltz:
        if viral_structure_db:
            prefix = Path(viral_structure_db)
            if prefix.exists() or prefix.with_suffix(".dbtype").exists():
                print(f"[OK] Boltz/FoldSeek DB prefix: {prefix}")
            else:
                print(f"[FAIL] Boltz/FoldSeek DB prefix missing: {prefix} or {prefix.with_suffix('.dbtype')}")
                failures += 1
        else:
            print("[FAIL] phase3.viral_structure_db is not configured")
            failures += 1

    if failures:
        print(f"[FAIL] structural runtime preflight failed ({failures} issue(s))")
        return 1

    print("[OK] structural runtime preflight passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
