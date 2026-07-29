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

import yaml

from virosync.utils.database_manager import ViroSyncDatabaseManager
from virosync.utils.executables import resolve_boltz_executable


def _load_phase3(config_path: Path) -> dict:
    data = yaml.safe_load(config_path.read_text()) or {}
    if "phase3" in data and isinstance(data["phase3"], dict):
        return dict(data["phase3"])
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


def _check_tmvec_model_smoke() -> tuple[bool, str]:
    try:
        from virosync.pipeline.phase3.tmvec_predictor import get_tmvec_predictor

        predictor = get_tmvec_predictor(device="cuda", require_gpu=True)
        try:
            embeddings = predictor.embed_batch(["M" * 60], max_residues=100, max_batch=1)
        finally:
            if hasattr(predictor, "release"):
                predictor.release()
    except Exception as exc:
        return False, str(exc)

    if embeddings.shape != (1, 512):
        return False, f"unexpected embedding shape: {embeddings.shape}"
    return True, "single-sequence embedding produced"


def _check_tmvec_database_integrity(
    tmvec_dir: str | Path,
    tmvec_databases: list[str],
) -> tuple[bool, str]:
    try:
        from virosync.pipeline.phase3.tmvec_database import TMVecDatabaseSearch

        searcher = TMVecDatabaseSearch(
            database_root=Path(tmvec_dir),
            databases=tmvec_databases,
            require_gpu=True,
        )
        for name in tmvec_databases:
            paths = searcher._db_paths.get(name)
            if not paths or name != "bfvd":
                continue
            emb_path = paths.get("embeddings")
            if emb_path is None or not emb_path.exists():
                continue
            legacy_reason = searcher._legacy_untrained_bfvd_reason(emb_path)
            if legacy_reason:
                return False, legacy_reason
    except Exception as exc:
        return False, str(exc)

    return True, "configured databases passed compatibility checks"


def main() -> int:
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
        help="Fail unless the TMVec runtime, CUDA, and configured TMVec databases are available.",
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
    args = parser.parse_args()

    config_path = args.config.resolve()
    if not config_path.exists():
        print(f"[FAIL] config not found: {config_path}")
        return 1

    phase3 = _load_phase3(config_path)
    tmvec_dir = phase3.get("tmvec_database_dir")
    tmvec_databases = phase3.get("tmvec_databases") or ["bfvd"]
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
    print(f"[INFO] phase3.interproscan_dir: {interpro_dir}")
    print(f"[INFO] phase3.viral_structure_db: {viral_structure_db}")
    print(f"[INFO] phase3.boltz_use_msa_server: {boltz_use_msa_server}")

    if not any([check_tmvec, check_boltz, check_interproscan]):
        print("[INFO] no optional structural/domain layers are enabled or required")

    tmvec_ready_for_model_smoke = True

    if check_tmvec:
        required_modules = [
            "torch",
            "tm_vec",
            "transformers",
            "pytorch_lightning",
        ]
        for name in required_modules:
            ok, msg = _check_module(name)
            state = "OK" if ok else "FAIL"
            print(f"[{state}] python module {name}: {msg}")
            if not ok:
                failures += 1
                tmvec_ready_for_model_smoke = False

        try:
            import torch

            cuda_ok = bool(torch.cuda.is_available())
            if cuda_ok:
                gpu_count = torch.cuda.device_count()
                print(f"[OK] CUDA runtime: available ({gpu_count} device(s))")
            else:
                print("[FAIL] CUDA runtime: torch.cuda.is_available() == False (TMVec will be disabled)")
                failures += 1
                tmvec_ready_for_model_smoke = False
        except Exception as exc:
            print(f"[FAIL] CUDA runtime: unable to query torch ({exc})")
            failures += 1
            tmvec_ready_for_model_smoke = False
    else:
        print("[SKIP] TMVec runtime checks (not enabled; pass --require-tmvec to force)")

    if check_boltz:
        if not boltz_use_msa_server:
            print(
                "[FAIL] phase3.boltz_use_msa_server is false; ViroSync's "
                "Boltz integration does not currently write local MSA inputs."
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
            if interpro_exec.exists() and interpro_exec.is_file():
                print(f"[OK] interproscan.sh: {interpro_exec}")
            else:
                print(f"[FAIL] interproscan.sh missing: {interpro_exec}")
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
                print("[OK] TMVec database files: complete")
                ok, msg = _check_tmvec_database_integrity(
                    tmvec_dir=tmvec_dir,
                    tmvec_databases=[str(db) for db in tmvec_databases],
                )
                state = "OK" if ok else "FAIL"
                print(f"[{state}] TMVec database compatibility: {msg}")
                if not ok:
                    failures += 1
                    tmvec_ready_for_model_smoke = False
        else:
            print("[FAIL] phase3.tmvec_database_dir is not configured")
            failures += 1
            tmvec_ready_for_model_smoke = False

        if tmvec_ready_for_model_smoke:
            ok, msg = _check_tmvec_model_smoke()
            state = "OK" if ok else "FAIL"
            print(f"[{state}] TMVec model smoke: {msg}")
            if not ok:
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
