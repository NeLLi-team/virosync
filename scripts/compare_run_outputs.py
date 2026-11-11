#!/usr/bin/env python3
"""Compare the scientific outputs of two ViroSync genome output folders.

The script compares the prediction tables, BED and GFF3 exports, the summary,
evidence profiles, gene-function tables, the Phase 2 taxonomy table and boundary
BED, the host signature model and the Phase 1 checkpoint. Tables must match byte
for byte. JSON files must match as documents: key order is free, list order is
not, and the key ``generated_at`` is ignored. The script prints each difference
and exits with status 1 when any exists, 2 when an input is unusable.

Example
-------
    pixi run python scripts/compare_run_outputs.py results_fresh/GENOME results_new/GENOME
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Iterator, Sequence
from itertools import islice
from pathlib import Path

# Each tuple lists the places where one output can lie. The first one that exists counts.
SCIENTIFIC_OUTPUTS: tuple[tuple[str, ...], ...] = (
    ("phase3_synthesis/virosync_predictions.tsv", "virosync_predictions.tsv"),
    ("virosync_predictions_detailed.tsv", "phase3_synthesis/virosync_predictions_detailed.tsv"),
    ("virosync_repeat_candidates.tsv", "phase3_synthesis/virosync_repeat_candidates.tsv"),
    ("phase3_synthesis/virosync_predictions.bed", "virosync_predictions.bed"),
    ("phase3_synthesis/virosync_predictions.gff3", "virosync_predictions.gff3"),
    ("phase3_synthesis/virosync_summary.json", "virosync_summary.json"),
    ("phase3_synthesis/evidence_profiles.json",),
    ("phase3_synthesis/gene_functions.tsv",),
    ("phase3_synthesis/gene_function_length_exclusions.tsv",),
    ("phase3_synthesis/eve_ani_edges.tsv",),
    ("phase2/boundary_diamond/taxonomy_map.tsv",),
    ("phase2/refined_boundaries.bed",),
    ("phase1/marker_validation/host_signature_model.json",),
    ("phase1/resume_state.json",),
)
# The only fields ignored in JSON files: they record when the file was written.
VOLATILE_JSON_KEYS = frozenset({"generated_at"})
MAX_JSON_DIFFERENCES = 10
MAX_LINE_CHARS = 200


def compare_runs(first: Path, second: Path) -> list[str]:
    """Return one message per difference between the scientific outputs of two folders."""
    differences: list[str] = []
    for alternatives in SCIENTIFIC_OUTPUTS:
        left, right = _locate(first, alternatives), _locate(second, alternatives)
        if left is None and right is None:
            continue
        if left is None or right is None:
            name = left or right
            differences.append(f"{name}: only in {'the first' if right is None else 'the second'} folder")
            continue
        if left != right:
            differences.append(f"{alternatives[0]}: location differs: {left} != {right}")
        differences.extend(_file_differences(left, first / left, second / right))
    return differences


def count_outputs(folder: Path) -> int:
    """Return how many scientific outputs the folder holds."""
    return sum(_locate(folder, alternatives) is not None for alternatives in SCIENTIFIC_OUTPUTS)


def _locate(folder: Path, alternatives: Sequence[str]) -> str | None:
    return next((relative for relative in alternatives if (folder / relative).is_file()), None)


def _file_differences(name: str, left: Path, right: Path) -> list[str]:
    if left.suffix == ".json":
        found = _json_differences(json.loads(left.read_text("utf-8")), json.loads(right.read_text("utf-8")))
        return [f"{name}: {message}" for message in islice(found, MAX_JSON_DIFFERENCES)]
    if left.read_bytes() == right.read_bytes():
        return []
    return [f"{name}: {_first_line_difference(left, right)}"]


def _first_line_difference(left: Path, right: Path) -> str:
    left_lines = left.read_text("utf-8", errors="replace").splitlines()
    right_lines = right.read_text("utf-8", errors="replace").splitlines()
    for number, (a, b) in enumerate(zip(left_lines, right_lines, strict=False), start=1):
        if a != b:
            return f"line {number} differs: {a[:MAX_LINE_CHARS]} != {b[:MAX_LINE_CHARS]}"
    return f"line count differs: {len(left_lines)} != {len(right_lines)}"


def _json_differences(left: object, right: object, path: str = "") -> Iterator[str]:
    if isinstance(left, dict) and isinstance(right, dict):
        for key in sorted((set(left) | set(right)) - VOLATILE_JSON_KEYS):
            child = f"{path}.{key}" if path else key
            if key not in left or key not in right:
                yield f"{child} only in {'the first' if key in left else 'the second'} folder"
            else:
                yield from _json_differences(left[key], right[key], child)
    elif isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            yield f"{path} has {len(left)} items, then {len(right)}"
        for index, (a, b) in enumerate(zip(left, right, strict=False)):
            yield from _json_differences(a, b, f"{path}[{index}]")
    elif left != right or type(left) is not type(right):
        yield f"{path or 'value'} differs: {left!r} != {right!r}"


def main(argv: Sequence[str] | None = None) -> int:
    """Compare two genome output folders and return the exit status."""
    parser = argparse.ArgumentParser(description="Compare the scientific outputs of two ViroSync genome folders.")
    parser.add_argument("first", type=Path, help="first genome output folder")
    parser.add_argument("second", type=Path, help="second genome output folder")
    args = parser.parse_args(argv)
    for folder in (args.first, args.second):
        if not folder.is_dir():
            print(f"error: not a directory: {folder}", file=sys.stderr)
            return 2
    if not (count_outputs(args.first) or count_outputs(args.second)):
        print("error: neither folder holds a scientific output", file=sys.stderr)
        return 2
    differences = compare_runs(args.first, args.second)
    for message in differences:
        print(message)
    if differences:
        print(f"{len(differences)} difference(s)")
        return 1
    print(f"identical: {max(count_outputs(args.first), count_outputs(args.second))} outputs compared")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
