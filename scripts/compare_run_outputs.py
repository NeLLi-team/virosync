#!/usr/bin/env python3
"""Compare the scientific outputs of two ViroSync genome output folders.

The script reads ``virosync_run_state.json`` in each folder. A run state other than
``success`` is a difference. A successful run must hold every output that its state
requires: the six exports and ``evidence_profiles.json`` for each run; the Phase 1
checkpoint, the host signature model and the three Phase 2 files for a run that did
not end early; and ``<genome_id>_eves.fna``, ``gene_taxonomy/gene_taxonomy_all.tsv``
and ``gvclass_input/manifest.tsv`` for a run with accepted EVEs. A missing required
output is a difference even when both folders lack it.

The script then compares the prediction tables, BED and GFF3 exports, the summary,
the evidence profiles and graph, the gene-function tables, the combined EVE FASTA,
the Phase 2 taxonomy table, boundary BED and checkpoint, the host signature model
and the Phase 1 checkpoint. It compares the directories ``gene_taxonomy``,
``gvclass_input`` and ``eve_sequences_all`` by member set and content. Tables and
sequence files must match byte for byte. JSON files must match as documents: key
order is free, list order is not, and the key ``generated_at`` is ignored. The
writers of these outputs hold no other field that records a time.

The script prints each difference and exits with status 1 when any exists, and with
status 2 when an input is unusable.

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
from typing import NamedTuple

# Each entry lists the places where one output can lie. The first one that exists counts.
# ``{genome_id}`` stands for the genome ID in the run state.
SCIENTIFIC_OUTPUTS: dict[str, tuple[str, ...]] = {
    "canonical predictions": ("phase3_synthesis/virosync_predictions.tsv", "virosync_predictions.tsv"),
    "detailed predictions": ("virosync_predictions_detailed.tsv", "phase3_synthesis/virosync_predictions_detailed.tsv"),
    "repeat candidates": ("virosync_repeat_candidates.tsv", "phase3_synthesis/virosync_repeat_candidates.tsv"),
    "BED": ("phase3_synthesis/virosync_predictions.bed", "virosync_predictions.bed"),
    "GFF3": ("phase3_synthesis/virosync_predictions.gff3", "virosync_predictions.gff3"),
    "summary": ("phase3_synthesis/virosync_summary.json", "virosync_summary.json"),
    "evidence profiles": ("phase3_synthesis/evidence_profiles.json", "evidence_profiles.json"),
    "evidence graph": ("phase3_synthesis/evidence_graph.json",),
    "gene functions": ("phase3_synthesis/gene_functions.tsv",),
    "gene function exclusions": ("phase3_synthesis/gene_function_length_exclusions.tsv",),
    "ANI edges": ("phase3_synthesis/eve_ani_edges.tsv",),
    "EVE sequences": ("{genome_id}_eves.fna",),
    "phase 2 taxonomy table": ("phase2/boundary_diamond/taxonomy_map.tsv",),
    "phase 2 boundaries": ("phase2/refined_boundaries.bed",),
    "phase 2 checkpoint": ("phase2/refined_state.json",),
    "host signature model": ("phase1/marker_validation/host_signature_model.json",),
    "phase 1 checkpoint": ("phase1/resume_state.json",),
}
SCIENTIFIC_DIRECTORIES: tuple[str, ...] = (
    "phase3_synthesis/gene_taxonomy",
    "phase3_synthesis/gvclass_input",
    "phase3_synthesis/eve_sequences_all",
)
_REQUIRED_EVERY_RUN = (
    "canonical predictions",
    "detailed predictions",
    "repeat candidates",
    "BED",
    "GFF3",
    "summary",
    "evidence profiles",
)
_REQUIRED_FULL_RUN = (
    "phase 1 checkpoint",
    "host signature model",
    "phase 2 taxonomy table",
    "phase 2 boundaries",
    "phase 2 checkpoint",
)
_REQUIRED_WITH_EVES = ("EVE sequences",)
_REQUIRED_WITH_EVES_FILES = (
    "phase3_synthesis/gene_taxonomy/gene_taxonomy_all.tsv",
    "phase3_synthesis/gvclass_input/manifest.tsv",
)
# The only field ignored in JSON files: it records when the file was written.
VOLATILE_JSON_KEYS = frozenset({"generated_at"})
RUN_STATE_FILENAME = "virosync_run_state.json"
MAX_DIFFERENCES_PER_FILE = 10
MAX_LINE_CHARS = 200


class RunSummary(NamedTuple):
    """The parts of a run state that decide which outputs a run must hold."""

    status: str
    genome_id: str
    terminal_phase: int | None
    accepted: int


def read_run_summary(folder: Path) -> RunSummary:
    """Read the completion state of a genome output folder."""
    path = folder / RUN_STATE_FILENAME
    try:
        payload = json.loads(path.read_text("utf-8"))
        result = payload.get("result") or {}
        return RunSummary(
            status=str(payload["status"]),
            genome_id=str(payload["identities"]["genome_id"]),
            terminal_phase=result.get("terminal_phase"),
            accepted=int(result.get("canonical_rows", 0)),
        )
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        raise ValueError(f"no readable {RUN_STATE_FILENAME} in {folder}: {exc}") from exc


def compare_runs(first: Path, second: Path) -> list[str]:
    """Return one message per difference between the scientific outputs of two folders."""
    summaries = {first: read_run_summary(first), second: read_run_summary(second)}
    differences: list[str] = []
    for folder, summary in summaries.items():
        differences.extend(_completion_problems(folder, summary))
    if summaries[first].terminal_phase != summaries[second].terminal_phase:
        differences.append(f"the first run ended {_ended(summaries[first])} and the second {_ended(summaries[second])}")
    for name, alternatives in SCIENTIFIC_OUTPUTS.items():
        left = _locate(first, alternatives, summaries[first].genome_id)
        right = _locate(second, alternatives, summaries[second].genome_id)
        if left is None and right is None:
            continue
        if left is None or right is None:
            side = "the first" if right is None else "the second"
            differences.append(f"{left or right}: only in {side} folder")
            continue
        if left != right:
            differences.append(f"{name}: location differs: {left} != {right}")
        differences.extend(_file_differences(left, first / left, second / right))
    for directory in SCIENTIFIC_DIRECTORIES:
        differences.extend(_directory_differences(directory, first / directory, second / directory))
    return differences


def count_outputs(folder: Path) -> int:
    """Return how many scientific outputs and directories the folder holds."""
    genome_id = read_run_summary(folder).genome_id
    files = sum(_locate(folder, alternatives, genome_id) is not None for alternatives in SCIENTIFIC_OUTPUTS.values())
    return files + sum((folder / directory).is_dir() for directory in SCIENTIFIC_DIRECTORIES)


def _ended(summary: RunSummary) -> str:
    return "after Phase 3" if summary.terminal_phase is None else f"in Phase {summary.terminal_phase}"


def _completion_problems(folder: Path, summary: RunSummary) -> list[str]:
    if summary.status != "success":
        return [f"{folder}: run state is {summary.status!r}, not 'success'"]
    required = list(_REQUIRED_EVERY_RUN)
    if summary.terminal_phase is None:
        required.extend(_REQUIRED_FULL_RUN)
    if summary.accepted:
        required.extend(_REQUIRED_WITH_EVES)
    problems = [
        f"{folder}: required output is missing: {name} ({SCIENTIFIC_OUTPUTS[name][0].format(genome_id=summary.genome_id)})"
        for name in required
        if _locate(folder, SCIENTIFIC_OUTPUTS[name], summary.genome_id) is None
    ]
    if summary.accepted:
        problems.extend(
            f"{folder}: required output is missing: {relative}"
            for relative in _REQUIRED_WITH_EVES_FILES
            if not (folder / relative).is_file()
        )
    return problems


def _locate(folder: Path, alternatives: Sequence[str], genome_id: str) -> str | None:
    candidates = (alternative.format(genome_id=genome_id) for alternative in alternatives)
    return next((relative for relative in candidates if (folder / relative).is_file()), None)


def _file_differences(name: str, left: Path, right: Path) -> list[str]:
    if left.suffix == ".json":
        found = _json_differences(json.loads(left.read_text("utf-8")), json.loads(right.read_text("utf-8")))
        return [f"{name}: {message}" for message in islice(found, MAX_DIFFERENCES_PER_FILE)]
    if left.read_bytes() == right.read_bytes():
        return []
    return [f"{name}: {_first_line_difference(left, right)}"]


def _directory_differences(directory: str, left: Path, right: Path) -> list[str]:
    if not left.is_dir() and not right.is_dir():
        return []
    if not left.is_dir() or not right.is_dir():
        return [f"{directory}: only in {'the first' if left.is_dir() else 'the second'} folder"]
    left_members, right_members = _members(left), _members(right)
    differences = [
        f"{directory}/{member}: only in {'the first' if member in left_members else 'the second'} folder"
        for member in sorted(left_members.keys() ^ right_members.keys())
    ]
    for member in sorted(left_members.keys() & right_members.keys()):
        if left_members[member].read_bytes() != right_members[member].read_bytes():
            differences.append(
                f"{directory}/{member}: {_first_line_difference(left_members[member], right_members[member])}"
            )
    return differences


def _members(directory: Path) -> dict[str, Path]:
    return {path.relative_to(directory).as_posix(): path for path in sorted(directory.rglob("*")) if path.is_file()}


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
    try:
        differences = compare_runs(args.first, args.second)
        compared = max(count_outputs(args.first), count_outputs(args.second))
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    for message in differences:
        print(message)
    if differences:
        print(f"{len(differences)} difference(s)")
        return 1
    print(f"identical: {compared} outputs compared")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
