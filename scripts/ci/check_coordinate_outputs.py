#!/usr/bin/env python3
"""Reject output roots that do not declare the current coordinate contract."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Iterable, Sequence
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from virosync.output_contract import coordinate_contract_metadata


METADATA_FILENAMES = (
    "virosync_run_state.json",
    "virosync_run_complete.json",
    "virosync_summary.json",
)
COORDINATE_OUTPUT_FILENAMES = (
    "virosync_predictions.tsv",
    "virosync_predictions_detailed.tsv",
    "virosync_predictions.bed",
    "virosync_predictions.gff3",
    "refined_boundaries.bed",
    "hhg_seeds.bed",
    "marker_seed_regions.bed",
    "marker_seed_regions.tsv",
    "validated_marker_hits.tsv",
)
CLEAN_RUN_GUIDANCE = (
    "Regenerate ViroSync outputs from raw inputs with --clean-run before "
    "R4 validation."
)
_MISSING = object()


def _metadata_paths(root: Path) -> list[Path]:
    if root.is_file():
        return [root] if root.name in METADATA_FILENAMES else []
    if not root.is_dir():
        return []
    return sorted(
        path
        for filename in METADATA_FILENAMES
        for path in root.rglob(filename)
        if path.is_file()
    )


def _coordinate_output_paths(root: Path) -> list[Path]:
    if not root.is_dir():
        return []
    return sorted(
        path
        for filename in COORDINATE_OUTPUT_FILENAMES
        for path in root.rglob(filename)
        if path.is_file()
    )


def _owning_run_dir(path: Path, run_dirs: set[Path]) -> Path | None:
    owners = [run_dir for run_dir in run_dirs if run_dir in path.parents]
    return max(owners, key=lambda owner: len(owner.parts), default=None)


def _metadata_errors(path: Path) -> list[str]:
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        return [f"{path}: invalid JSON ({exc}). {CLEAN_RUN_GUIDANCE}"]
    if not isinstance(payload, dict):
        return [
            f"{path}: metadata root must be a JSON object. "
            f"{CLEAN_RUN_GUIDANCE}"
        ]

    errors = []
    contract_payload = payload
    if path.name == "virosync_run_state.json":
        if payload.get("schema_version") != 3 or payload.get("status") != "success":
            errors.append(
                f"{path}: authoritative run state is not schema-v3 success. "
                f"{CLEAN_RUN_GUIDANCE}"
            )
        identities = payload.get("identities")
        contract_payload = identities if isinstance(identities, dict) else {}
    for field, expected in coordinate_contract_metadata().items():
        actual = contract_payload.get(field, _MISSING)
        if type(actual) is type(expected) and actual == expected:
            continue
        rendered = "<missing>" if actual is _MISSING else repr(actual)
        errors.append(
            f"{path}: {field}={rendered}; expected {expected!r}. "
            f"{CLEAN_RUN_GUIDANCE}"
        )
    return errors


def check_coordinate_output_roots(roots: Iterable[Path]) -> list[str]:
    """Return completion/summary contract errors below ``roots``."""
    errors = []
    for raw_root in roots:
        root = Path(raw_root)
        metadata_paths = _metadata_paths(root)
        coordinate_paths = _coordinate_output_paths(root)
        if not metadata_paths:
            errors.append(
                f"{root}: no completion or summary metadata found. "
                f"{CLEAN_RUN_GUIDANCE}"
            )
        for path in metadata_paths:
            errors.extend(_metadata_errors(path))

        manifest_paths = [
            path
            for path in metadata_paths
            if path.name in {"virosync_run_state.json", "virosync_run_complete.json"}
        ]
        summary_paths = [
            path
            for path in metadata_paths
            if path.name == "virosync_summary.json"
        ]
        run_dirs = {path.parent for path in manifest_paths}
        owned_summaries: dict[Path, list[Path]] = {
            run_dir: [] for run_dir in run_dirs
        }
        owned_coordinates: dict[Path, list[Path]] = {
            run_dir: [] for run_dir in run_dirs
        }

        for path in summary_paths:
            owner = _owning_run_dir(path, run_dirs)
            if owner is None:
                errors.append(
                    f"{path}: no enclosing ViroSync completion metadata "
                    "for this summary. "
                    f"{CLEAN_RUN_GUIDANCE}"
                )
            else:
                owned_summaries[owner].append(path)

        for path in coordinate_paths:
            owner = _owning_run_dir(path, run_dirs)
            if owner is None:
                errors.append(
                    f"{path}: no enclosing ViroSync completion metadata "
                    "for this coordinate output. "
                    f"{CLEAN_RUN_GUIDANCE}"
                )
            else:
                owned_coordinates[owner].append(path)

        for run_dir, paths in owned_coordinates.items():
            if paths and not owned_summaries[run_dir]:
                errors.append(
                    f"{run_dir}: coordinate outputs exist but "
                    "virosync_summary.json is missing. "
                    f"{CLEAN_RUN_GUIDANCE}"
                )
    return errors


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "roots",
        nargs="+",
        type=Path,
        help="Result or benchmark roots to inspect recursively",
    )
    args = parser.parse_args(argv)

    errors = check_coordinate_output_roots(args.roots)
    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1
    print(
        "Coordinate output metadata contract verified for "
        f"{len(args.roots)} root(s)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
