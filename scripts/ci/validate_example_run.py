#!/usr/bin/env python3
"""Validate a completed example run and snapshot authenticated outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path, PurePath
import re
import sys
from typing import Sequence


ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from virosync.orchestration._flows.single_genome.run_state import (
    RUN_STATE_FILENAME,
    RUN_STATE_SCHEMA_VERSION,
    load_run_state,
    plan_resume,
)


SNAPSHOT_SCHEMA_VERSION = 4
PREDICTION_TABLES = {
    "canonical_eve_ids": Path("phase3_synthesis/virosync_predictions.tsv"),
    "detailed_eve_ids": Path("phase3_synthesis/virosync_predictions_detailed.tsv"),
}
CLASS_FIELDS = {
    "ncldv": "NCLDV",
    "mirus": "MIRUS",
    "ppv": "PPV",
    "cress": "CRESS",
    "phage": "PHAGE",
    "viral_unknown": "VIRAL_UNKNOWN",
    "unknown": "UNKNOWN",
}
TIER_FIELDS = {
    "high_tier": "HIGH",
    "medium_tier": "MEDIUM",
    "low_tier": "LOW",
}
SNAPSHOT_SUMMARY_FIELDS = (
    "genome_id",
    "status",
    "benchmark_eligible",
    "legacy_resume",
    "predictions",
    "accepted",
    "high_tier",
    "medium_tier",
    "low_tier",
    "ncldv",
    "mirus",
    "ppv",
    "cress",
    "phage",
    "viral_unknown",
    "unknown",
    "total_bp",
    "genes",
    "hallmarks",
)


class ExampleValidationError(ValueError):
    """The example output does not satisfy the release-smoke contract."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _nonnegative_int(row: dict[str, str], field: str) -> int:
    try:
        value = int(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        raise ExampleValidationError(
            f"batch summary field {field!r} must be an integer"
        ) from exc
    if value < 0:
        raise ExampleValidationError(
            f"batch summary field {field!r} must be nonnegative"
        )
    return value


def _elapsed_seconds(row: dict[str, str]) -> float:
    try:
        value = float(row["elapsed_sec"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ExampleValidationError(
            "batch summary elapsed_sec must be numeric"
        ) from exc
    if value < 0:
        raise ExampleValidationError(
            "batch summary elapsed_sec must be nonnegative"
        )
    return value


def _load_batch_rows(output_root: Path) -> list[dict[str, str]]:
    summary = output_root / "batch_summary.tsv"
    try:
        with summary.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except OSError as exc:
        raise ExampleValidationError(f"cannot read {summary.name}: {exc}") from exc
    if not rows:
        raise ExampleValidationError("batch_summary.tsv has no result rows")
    genome_ids = [row.get("genome_id", "") for row in rows]
    if not all(genome_ids) or len(set(genome_ids)) != len(genome_ids):
        raise ExampleValidationError(
            "batch_summary.tsv genome IDs must be nonempty and unique"
        )
    return rows


def _load_eve_ids(run_dir: Path, relative_path: Path) -> list[str]:
    path = run_dir / relative_path
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if "eve_id" not in (reader.fieldnames or ()):
                raise ExampleValidationError(f"{relative_path} lacks an eve_id column")
            eve_ids = [row["eve_id"] for row in reader]
    except OSError as exc:
        raise ExampleValidationError(f"cannot read {relative_path}: {exc}") from exc
    if not all(eve_ids) or len(set(eve_ids)) != len(eve_ids):
        raise ExampleValidationError(f"{relative_path} EVE IDs must be nonempty and unique")
    return sorted(eve_ids)


def _validated_genome_dir(output_root: Path, genome_id: str) -> Path:
    component = PurePath(genome_id)
    if (
        component.is_absolute()
        or len(component.parts) != 1
        or component.parts[0] in {"", ".", ".."}
    ):
        raise ExampleValidationError(
            f"batch summary contains unsafe genome ID {genome_id!r}"
        )
    run_dir = output_root / genome_id
    if not run_dir.is_dir():
        raise ExampleValidationError(f"missing result directory for {genome_id}")
    return run_dir


def _require_batch_contract(
    row: dict[str, str],
    result: dict[str, object],
    *,
    require_resume: bool,
) -> None:
    genome_id = row["genome_id"]
    if row.get("status") != "success":
        raise ExampleValidationError(f"{genome_id}: batch status is not success")
    if row.get("benchmark_eligible") != "true":
        raise ExampleValidationError(f"{genome_id}: run is not benchmark eligible")
    if row.get("legacy_resume") != "false":
        raise ExampleValidationError(f"{genome_id}: legacy resume was reported")
    if row.get("error", ""):
        raise ExampleValidationError(f"{genome_id}: batch summary reports an error")
    if result.get("benchmark_eligible") is not True:
        raise ExampleValidationError(
            f"{genome_id}: authoritative state is not benchmark eligible"
        )

    expected_counts = {
        "predictions": result.get("detailed_rows"),
        "accepted": result.get("canonical_rows"),
        "total_bp": result.get("accepted_bp"),
    }
    tier_counts = result.get("tier_counts")
    class_counts = result.get("class_counts")
    if not isinstance(tier_counts, dict) or not isinstance(class_counts, dict):
        raise ExampleValidationError(
            f"{genome_id}: authoritative state lacks tier/class counts"
        )
    expected_counts.update(
        {field: tier_counts.get(tier) for field, tier in TIER_FIELDS.items()}
    )
    expected_counts.update(
        {field: class_counts.get(label) for field, label in CLASS_FIELDS.items()}
    )
    for field, expected in expected_counts.items():
        actual = _nonnegative_int(row, field)
        if type(expected) is not int or actual != expected:
            raise ExampleValidationError(
                f"{genome_id}: batch {field}={actual} differs from state {expected!r}"
            )
    for field in ("genes", "hallmarks"):
        _nonnegative_int(row, field)
    elapsed = _elapsed_seconds(row)
    if require_resume and elapsed != 0:
        raise ExampleValidationError(
            f"{genome_id}: unchanged resume elapsed_sec must be 0, got {elapsed}"
        )


def build_snapshot(
    output_root: Path,
    *,
    require_resume: bool = False,
) -> dict[str, object]:
    """Validate *output_root* and return a deterministic, path-free snapshot."""

    output_root = Path(output_root)
    if not output_root.is_dir():
        raise ExampleValidationError(
            f"example output root is not a directory: {output_root}"
        )
    runs: dict[str, object] = {}
    for row in _load_batch_rows(output_root):
        genome_id = row["genome_id"]
        run_dir = _validated_genome_dir(output_root, genome_id)
        try:
            state = load_run_state(run_dir)
        except (OSError, TypeError, ValueError) as exc:
            raise ExampleValidationError(
                f"{genome_id}: invalid schema-v3 run state: {exc}"
            ) from exc
        if state.schema_version != RUN_STATE_SCHEMA_VERSION or state.status != "success":
            raise ExampleValidationError(
                f"{genome_id}: authoritative run state is not schema-v3 success"
            )
        if re.fullmatch(r"[0-9a-f]{64}", state.run_fingerprint) is None:
            raise ExampleValidationError(
                f"{genome_id}: run fingerprint is not a SHA-256 digest"
            )
        plan = plan_resume(
            run_dir,
            expected_run_fingerprint=state.run_fingerprint,
        )
        if not plan.completed:
            raise ExampleValidationError(
                f"{genome_id}: completion artifacts are stale: {plan.reason}"
            )
        if not isinstance(state.result, dict):
            raise ExampleValidationError(
                f"{genome_id}: success state has no result payload"
            )
        _require_batch_contract(row, state.result, require_resume=require_resume)

        artifacts = [
            {
                "relative_path": artifact.relative_path,
                "size": artifact.size,
                "sha256": artifact.sha256,
                "schema": artifact.schema,
                "row_count": artifact.row_count,
            }
            for artifact in sorted(
                state.artifacts,
                key=lambda item: item.relative_path,
            )
        ]
        runs[genome_id] = {
            "run_fingerprint": state.run_fingerprint,
            "attempt": state.attempt,
            "result": state.result,
            "artifacts": artifacts,
            **{
                field: _load_eve_ids(run_dir, relative_path)
                for field, relative_path in PREDICTION_TABLES.items()
            },
            "run_state_sha256": _sha256(run_dir / RUN_STATE_FILENAME),
            "batch": {field: row[field] for field in SNAPSHOT_SUMMARY_FIELDS},
        }
    return {
        "schema_version": SNAPSHOT_SCHEMA_VERSION,
        "runs": runs,
    }


def _load_snapshot(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ExampleValidationError(f"cannot read snapshot {path}: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != SNAPSHOT_SCHEMA_VERSION:
        raise ExampleValidationError(
            f"snapshot {path} is not schema v{SNAPSHOT_SCHEMA_VERSION}"
        )
    return payload


def _require_expected_totals(
    snapshot: dict[str, object],
    *,
    predictions: int | None,
    accepted: int | None,
) -> None:
    runs = snapshot["runs"]
    if not isinstance(runs, dict):
        raise ExampleValidationError("snapshot runs payload must be an object")
    totals = {"predictions": 0, "accepted": 0}
    for run in runs.values():
        if not isinstance(run, dict) or not isinstance(run.get("batch"), dict):
            raise ExampleValidationError("snapshot run lacks batch counts")
        batch = run["batch"]
        for field in totals:
            totals[field] += int(batch[field])
    for field, expected in (("predictions", predictions), ("accepted", accepted)):
        if expected is not None and totals[field] != expected:
            raise ExampleValidationError(
                f"example {field} total {totals[field]} differs from expected {expected}"
            )


def _require_expected_eve_ids(
    snapshot: dict[str, object],
    *,
    canonical: list[str],
    detailed: list[str],
) -> None:
    if not canonical and not detailed:
        return
    runs = snapshot["runs"]
    if not isinstance(runs, dict) or len(runs) != 1:
        raise ExampleValidationError("exact EVE ID checks require an output root with one genome")
    run = next(iter(runs.values()))
    if not isinstance(run, dict):
        raise ExampleValidationError("snapshot run payload must be an object")
    observed = {field: run.get(field) for field in PREDICTION_TABLES}
    if not all(isinstance(ids, list) for ids in observed.values()):
        raise ExampleValidationError("snapshot run lacks EVE ID lists")
    for label, expected, field in (
        ("canonical", canonical, "canonical_eve_ids"),
        ("detailed", detailed, "detailed_eve_ids"),
    ):
        if expected and sorted(expected) != observed[field]:
            raise ExampleValidationError(
                f"example {label} EVE IDs {observed[field]!r} differ from expected {sorted(expected)!r}"
            )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_root", type=Path)
    snapshot = parser.add_mutually_exclusive_group()
    snapshot.add_argument("--write-snapshot", type=Path)
    snapshot.add_argument("--compare-snapshot", type=Path)
    parser.add_argument("--require-resume", action="store_true")
    parser.add_argument("--expect-predictions", type=int)
    parser.add_argument("--expect-accepted", type=int)
    parser.add_argument("--expect-canonical-eve-id", action="append", default=[])
    parser.add_argument("--expect-detailed-eve-id", action="append", default=[])
    args = parser.parse_args(argv)
    if args.require_resume and args.compare_snapshot is None:
        parser.error("--require-resume requires --compare-snapshot")
    for name in ("expect_predictions", "expect_accepted"):
        if getattr(args, name) is not None and getattr(args, name) < 0:
            parser.error(f"--{name.replace('_', '-')} must be nonnegative")

    try:
        current = build_snapshot(
            args.output_root,
            require_resume=args.require_resume,
        )
        _require_expected_totals(
            current,
            predictions=args.expect_predictions,
            accepted=args.expect_accepted,
        )
        _require_expected_eve_ids(
            current,
            canonical=args.expect_canonical_eve_id,
            detailed=args.expect_detailed_eve_id,
        )
        if args.write_snapshot is not None:
            args.write_snapshot.parent.mkdir(parents=True, exist_ok=True)
            args.write_snapshot.write_text(
                json.dumps(current, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        if args.compare_snapshot is not None:
            expected = _load_snapshot(args.compare_snapshot)
            if current != expected:
                raise ExampleValidationError(
                    "authenticated run fingerprint, counts, or artifact identities changed"
                )
    except ExampleValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    artifact_count = sum(
        len(run["artifacts"]) for run in current["runs"].values()  # type: ignore[union-attr,index]
    )
    print(
        f"Validated {len(current['runs'])} schema-v3 run(s) and "
        f"{artifact_count} authenticated artifact(s)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
