"""Start a run from the Phase 0 and Phase 1 results of an earlier run.

A source change after Phase 1 changes the run fingerprint, so a normal resume
repeats the Phase 1 searches. This module validates an earlier run of the same
genome, copies its ``phase0/`` and ``phase1/`` files into a fresh output folder
and rebuilds the host signature model with the current code. The orchestrator
then publishes both phases under the new fingerprint and runs Phase 2 and
Phase 3 in full. The earlier folder is only read.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
from collections.abc import Collection, Iterable, Mapping
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import cast

from Bio import SeqIO

from virosync.ablation import InterventionCounts
from virosync.config import MaskingBackend, PipelineConfig
from virosync.pipeline.host_signatures import TaxonomyLabelLookup, set_taxonomy_lookup
from virosync.pipeline.phase0.masking import load_masking_result, write_masking_status
from virosync.pipeline.phase1.frameshift_screening import is_rescued_protein_id
from virosync.pipeline.phase1.hhg_seeding import HMMHit
from virosync.pipeline.phase1.marker_validation import (
    NovelMarkerCriteria,
    ValidatedMarkerHit,
    filter_validated_markers,
)
from virosync.utils.atomic_write import atomic_write

from .phase1 import Phase1Result, Phase1Terminal, _build_host_model, _finish_without_seeds
from .phase1_state import Phase1ResumeState, load_phase1_state, write_phase1_state
from .run_state import (
    RUN_STATE_FILENAME,
    _validate_run_identity,
    load_run_state,
    marker_sha256,
    validate_phase_prefix,
)

PROVENANCE_SCHEMA = "virosync.phase_reuse/v1"
PROVENANCE_FILENAME = "reuse_provenance.json"

_MASKING_STATUS = "phase0/masking/masking_status.json"
_MODEL_JSON = "phase1/marker_validation/host_signature_model.json"
_MODEL_DEBUG = "phase1/marker_validation/host_signature_model_debug.tsv"
_CHECKPOINT = "phase1/resume_state.json"
_HOST_FASTA = "phase1/marker_validation/host_hmm_hit_porfs.faa"
_HOST_TSV = "phase1/marker_validation/host_diamond_top10.tsv"
_RESCUE_HOST_FASTA = "phase1/frameshift_screening/validation/host_hmm_hit_porfs.faa"
# Required ordinary-marker search files, authenticated by new Phase 1 markers. The protein search
# writes all four whenever it finds a hit. Without them, jelly-roll classification, the viral
# top-ten flags or the gene-function evidence drop out without an error.
_SEARCH_FILES = (
    "phase1/marker_validation/validated_marker_hits.tsv",
    "phase1/marker_validation/diamond_top10_taxonomy.tsv",
    "phase1/marker_validation/hmm_hit_porfs.faa",
    "phase1/hmm/hmm_hits_validated.tsv",
)
PHASE1_CONSUMED_FILES = (
    *_SEARCH_FILES,
    "phase1/novelty/diamond/diamond_combined.tsv",
    "phase1/hhg/hmm_hits_validated.tsv",
    _MODEL_JSON,
    _MODEL_DEBUG,
    _HOST_FASTA,
    _HOST_TSV,
    _RESCUE_HOST_FASTA,
    "phase1/frameshift_screening/validation/host_diamond_top10.tsv",
    "phase1/frameshift_screening/confirmed_frameshift_markers.tsv",
    "phase1/frameshift_screening/confirmed_frameshift_proteins.faa",
    "phase1/frameshift_screening/frameshift_hits.tsv",
    "phase1/frameshift_screening/frameshift_events.tsv",
    "phase1/pfam_arbitration.tsv",
)
# Phase 1 files that the host signature model changes. They are rebuilt, not copied.
_PHASE1_REBUILT = (_MODEL_JSON, _MODEL_DEBUG, _CHECKPOINT)
# An earlier import record is never copied: the new record names its source run instead.
_PHASE0_NOT_COPIED = frozenset({_MASKING_STATUS, f"phase0/{PROVENANCE_FILENAME}"})
_PHASE1_NOT_COPIED = frozenset({*_PHASE1_REBUILT, f"phase1/{PROVENANCE_FILENAME}"})
# Run identity components that may differ: the paths of the run and the code.
_UNCOMPARED_COMPONENTS = frozenset({"input_path", "output_dir", "code"})
# The host model reads neither the model name nor the score of a host-search query.
_HOST_QUERY_MODEL = "host_search_query"
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_COPY_CHUNK_BYTES = 1024 * 1024


class Phase1ReuseError(ValueError):
    """The earlier run, or its folder, cannot serve as the source of this run."""


@dataclass(frozen=True, slots=True)
class Phase1Reuse:
    """Request to start from the Phase 0 and Phase 1 results of an earlier run.

    ``old_base`` is the ``-o`` folder of the earlier run, which holds one subfolder
    per genome. ``code_sha256`` is the source SHA-256 that the earlier run recorded.
    """

    old_base: Path
    code_sha256: str

    def __post_init__(self) -> None:
        if _SHA256.fullmatch(self.code_sha256) is None:
            raise Phase1ReuseError("the code SHA-256 must be 64 lowercase hexadecimal characters")

    def old_dir(self, genome_id: str) -> Path:
        """Return the folder of the earlier run of one genome."""
        return self.old_base / genome_id

    def require_disjoint(self, output_base: Path, genome_id: str | None = None) -> None:
        """Refuse an output folder that equals, contains or lies inside the earlier folder.

        Checks the base folders and, with ``genome_id``, the two genome folders.
        """
        _require_disjoint(self.old_base, output_base)
        if genome_id is not None:
            _require_disjoint(self.old_dir(genome_id), Path(output_base) / genome_id)


@dataclass(frozen=True, slots=True)
class ReusedRun:
    """An earlier run whose Phase 0 and Phase 1 passed validation.

    ``markers`` holds the SHA-256 of its Phase 0 and Phase 1 markers.
    ``terminal_reason`` is the early-exit reason of a run that ended in Phase 1.
    ``checkpoint`` is the Phase 1 checkpoint of a run that did not end there.
    """

    root: Path
    run_fingerprint: str
    code_sha256: str
    markers: Mapping[int, str]
    terminal_reason: str | None = None
    checkpoint: Phase1ResumeState | None = None
    artifacts: Mapping[str, str] = field(default_factory=dict)


def _require_disjoint(old: Path, new: Path) -> None:
    old_path, new_path = Path(old).resolve(), Path(new).resolve()
    if old_path == new_path:
        raise Phase1ReuseError(f"the new output folder {new} equals the earlier folder {old}")
    if old_path in new_path.parents:
        raise Phase1ReuseError(f"the new output folder {new} lies inside the earlier folder {old}")
    if new_path in old_path.parents:
        raise Phase1ReuseError(f"the new output folder {new} contains the earlier folder {old}")


def open_reused_run(request: Phase1Reuse, *, genome_id: str, identity: Mapping[str, object]) -> ReusedRun:
    """Validate the earlier run of one genome against the identity of the new run.

    Reads the earlier folder only. Raises ``Phase1ReuseError`` with the reason when
    the run state, the Phase 0 and Phase 1 records, the recorded artifact hashes,
    the code SHA-256, the ablation arm or any identity component other than the
    code and the two paths does not match, when its masking fell back, or when it
    lacks a Phase 1 file that Phase 2 or Phase 3 read.
    """
    root = request.old_dir(genome_id)
    for path in (request.old_base, root, root / "phase0", root / "phase1"):
        if path.is_symlink():
            raise Phase1ReuseError(f"the earlier run contains a symlink root: {path}")
    try:
        state = load_run_state(root)
    except (OSError, TypeError, ValueError) as exc:
        hint = ""
        if (request.old_base / RUN_STATE_FILENAME).exists():
            hint = "; pass the base folder of the earlier run (its -o value), not a genome folder"
        raise Phase1ReuseError(f"no earlier run at {root}: {exc}{hint}") from exc
    if state.status == "running":
        raise Phase1ReuseError("the earlier run is still running")
    code = cast(Mapping[str, str], state.identities["code"])["source_sha256"]
    if code != request.code_sha256:
        raise Phase1ReuseError(f"the earlier run records code SHA-256 {code}, not {request.code_sha256}")
    differing = _differing_components(state.identities, _validate_run_identity(identity))
    if differing:
        raise Phase1ReuseError("the run identity differs from the earlier run in: " + ", ".join(differing))
    ablation = cast(Mapping[str, str], state.identities["config"])["ablation_id"]
    if ablation != "A0":
        raise Phase1ReuseError(f"only A0 runs can be reused; the earlier run used {ablation}")
    records, reason = validate_phase_prefix(root, expected_run_fingerprint=state.run_fingerprint, through_phase=1)
    if len(records) < 2:
        raise Phase1ReuseError(
            f"the earlier run at {root} has no valid Phase 0 and Phase 1: {reason or 'Phase 0 ends the run'}. "
            "The validators open the paths that the earlier run recorded, so its folder and its input genome "
            f"{state.identities['input_path']} must still lie at the recorded paths"
        )
    artifacts = {}
    for record in records:
        for artifact in record.artifacts:
            previous = artifacts.setdefault(artifact.relative_path, artifact.sha256)
            if previous != artifact.sha256:
                raise Phase1ReuseError(f"conflicting authenticated hashes: {artifact.relative_path}")
    for relative in PHASE1_CONSUMED_FILES:
        if (root / relative).exists():
            _require_authenticated(root, relative, artifacts)
    _require_authenticated(root, _MASKING_STATUS, artifacts)
    _require_masking_without_fallback(root)
    terminal = records[1].outcome == "terminal_zero"
    checkpoint = None if terminal else _load_complete_checkpoint(root, artifacts)
    if terminal:
        _require_authenticated(root, "virosync_run_complete.json", artifacts)
    return ReusedRun(
        root=root,
        run_fingerprint=state.run_fingerprint,
        code_sha256=code,
        markers={phase: marker_sha256(root, phase) for phase in (0, 1)},
        terminal_reason=_terminal_reason(root) if terminal else None,
        checkpoint=checkpoint,
        artifacts=artifacts,
    )


def _require_masking_without_fallback(root: Path) -> None:
    """Refuse a run whose masking fell back: a new run could mask with the requested backend."""
    status = json.loads((root / _MASKING_STATUS).read_text(encoding="utf-8"))
    if status["status"] == "fallback":
        raise Phase1ReuseError(
            f"the earlier run's masking fell back from {status['requested_backend']} to "
            f"{status['effective_backend']} ({status['fallback_reason']}); a new run could mask with "
            "the requested backend and give another proteome"
        )


def _load_complete_checkpoint(root: Path, artifacts: Mapping[str, str]) -> Phase1ResumeState:
    """Load the Phase 1 checkpoint and require the search files that a run with a protein hit wrote."""
    _require_authenticated(root, _CHECKPOINT, artifacts)
    checkpoint = load_phase1_state(root / _CHECKPOINT)
    _require_authenticated(root, _CHECKPOINT, artifacts)
    if not checkpoint.validated_markers:
        raise Phase1ReuseError("empty nonterminal Phase 1 checkpoint cannot be reused")
    if all(is_rescued_protein_id(marker.query_porf) for marker in checkpoint.validated_markers):
        return checkpoint
    missing = [relative for relative in _SEARCH_FILES if not (root / relative).is_file()]
    if missing:
        raise Phase1ReuseError(
            f"the earlier run at {root} lacks Phase 1 files that Phase 2 or Phase 3 read: " + ", ".join(missing)
        )
    for relative in _SEARCH_FILES:
        _require_authenticated(root, relative, artifacts)
    return checkpoint


def _require_authenticated(root: Path, relative: str, artifacts: Mapping[str, str]) -> None:
    expected = artifacts.get(relative)
    if expected is None:
        raise Phase1ReuseError(f"no authenticated old phase artifact for {relative}; run a fresh Phase 1")
    path = root / relative
    if path.is_symlink() or not path.is_file():
        raise Phase1ReuseError(f"authenticated input is not a regular file: {relative}")
    with path.open("rb") as handle:
        actual = hashlib.file_digest(handle, "sha256").hexdigest()
    if actual != expected:
        raise Phase1ReuseError(f"authenticated input changed: {relative}")


def require_matching_phase0_origin(source: ReusedRun, output_dir: Path) -> None:
    """Reject a reusable Phase 0 from a fresh run or a different import origin."""
    path = output_dir / "phase0" / PROVENANCE_FILENAME
    try:
        provenance = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise Phase1ReuseError("reusable Phase 0 has no authenticated import origin") from exc
    if (
        provenance.get("old_run_folder") != str(source.root.resolve())
        or provenance.get("old_phase_marker_sha256") != source.markers[0]
        or provenance.get("old_run_fingerprint") != source.run_fingerprint
    ):
        raise Phase1ReuseError("reusable Phase 0 import origin differs")


def _differing_components(earlier: Mapping[str, object], current: Mapping[str, object]) -> list[str]:
    differing: list[str] = []
    for key in sorted(set(earlier) | set(current)):
        before, after = earlier.get(key), current.get(key)
        if key in _UNCOMPARED_COMPONENTS or before == after:
            continue
        if key == "resources":
            changed = _resource_keys(before) ^ _resource_keys(after)
            differing.append("resources (" + ", ".join(sorted({name for name, *_ in changed})) + ")")
        elif isinstance(before, dict) and isinstance(after, dict):
            differing.extend(
                f"{key}.{field}" for field in sorted(set(before) | set(after)) if before.get(field) != after.get(field)
            )
        else:
            differing.append(key)
    return differing


def _resource_keys(resources: object) -> set[tuple[str, str, str]]:
    items = cast(Iterable[Mapping[str, str]], resources or ())
    return {(item["name"], item["version"], item["manifest_sha256"]) for item in items}


def _terminal_reason(root: Path) -> str:
    manifest = root / "virosync_run_complete.json"
    try:
        reason = json.loads(manifest.read_text(encoding="utf-8"))["reason"]
    except (OSError, KeyError, ValueError) as exc:
        raise Phase1ReuseError(f"cannot read the early-exit reason of the earlier run from {manifest}: {exc}") from exc
    if not isinstance(reason, str):
        raise Phase1ReuseError(f"the early-exit reason in {manifest} is not a string")
    return reason


def import_phase0(source: ReusedRun, *, output_dir: Path, genome_path: Path) -> None:
    """Copy ``phase0/`` and point its masking status at the new run.

    The status records the absolute path of the masked genome, which lies in the
    run folder, and a fingerprint of that record.
    """
    copied = _copy_tree(source, output_dir, "phase0", skip=_PHASE0_NOT_COPIED)
    _rebase_masking_status(source, output_dir, genome_path)
    _write_provenance(source, output_dir, 0, copied, _rewritten(source, output_dir, [_MASKING_STATUS]))


def _rebase_masking_status(source: ReusedRun, output_dir: Path, genome_path: Path) -> None:
    _require_authenticated(source.root, _MASKING_STATUS, source.artifacts)
    status = source.root / _MASKING_STATUS
    earlier = load_masking_result(status)
    region_count = json.loads(status.read_text(encoding="utf-8"))["repeat_region_count"]
    if earlier.effective_backend is MaskingBackend.OFF:
        output_path = Path(genome_path)
    else:
        try:
            relative = earlier.output_path.resolve().relative_to(source.root.resolve())
        except ValueError as exc:
            raise Phase1ReuseError(f"the masked genome {earlier.output_path} is outside the earlier run") from exc
        output_path = Path(output_dir).resolve() / relative
    _require_authenticated(source.root, _MASKING_STATUS, source.artifacts)
    write_masking_status(
        replace(earlier, output_path=output_path),
        Path(output_dir) / "phase0" / "masking",
        repeat_region_count=region_count,
    )


def import_phase1(
    source: ReusedRun,
    *,
    output_dir: Path,
    genome_id: str,
    config: PipelineConfig,
    masked_path: Path,
    proteome_path: Path,
    config_fingerprint: str,
    logger,
) -> Phase1Result | Phase1Terminal:
    """Copy ``phase1/``, rebuild the host signature model and return the Phase 1 outcome.

    A run that ended in Phase 1 without seeds writes its empty outputs again, with
    the early-exit reason of the earlier run. Every other run loads its lossless
    checkpoint and replaces the model, the signature set and their files.
    """
    copied = _copy_tree(source, output_dir, "phase1", skip=_PHASE1_NOT_COPIED)
    model_inputs = {
        "output_dir": output_dir,
        "config": config,
        "masked_path": masked_path,
        "proteome_path": proteome_path,
        "logger": logger,
    }
    if source.terminal_reason is not None:
        if (source.root / _MODEL_JSON).is_file():
            _rebuild_host_model(**model_inputs)
        _write_provenance(source, output_dir, 1, copied, _rewritten(source, output_dir, _PHASE1_REBUILT))
        logger.info("Phase 1 imported from %s (%d files copied)", source.root, len(copied))
        return _finish_without_seeds(
            reason=source.terminal_reason,
            masked_path=masked_path,
            proteome_path=proteome_path,
            output_dir=output_dir,
            genome_id=genome_id,
            config=config,
            logger=logger,
            config_fingerprint=config_fingerprint,
            elapsed_sec=0.0,
            ablation_counts=InterventionCounts(),
        )
    if source.checkpoint is None:
        raise Phase1ReuseError("nonterminal import lacks its validated checkpoint")
    _require_authenticated(source.root, _CHECKPOINT, source.artifacts)
    state = source.checkpoint
    confirmed = {
        m.query_porf.split("|aa", 1)[0] for m in state.validated_markers if is_rescued_protein_id(m.query_porf)
    }
    check_rescue_host_search(output_dir / _RESCUE_HOST_FASTA, confirmed)
    host_signatures, model = _rebuild_host_model(**model_inputs)
    write_phase1_state(
        output_dir / _CHECKPOINT,
        validated_markers=state.validated_markers,
        hallmark_evidence=state.hallmark_evidence,
        hallmark_witnesses=state.hallmark_witnesses,
        merged_seeds=state.merged_seeds,
        host_signature_model=model,
        host_signatures=host_signatures,
        host_deviation_summary=state.host_deviation_summary,
    )
    _write_provenance(source, output_dir, 1, copied, _rewritten(source, output_dir, _PHASE1_REBUILT))
    logger.info(
        "Phase 1 imported from %s (%d files copied, %d markers, %d seeds)",
        source.root,
        len(copied),
        len(state.validated_markers),
        len(state.merged_seeds),
    )
    return Phase1Result(
        merged_seeds=state.merged_seeds,
        validated_markers=state.validated_markers,
        hallmark_evidence=state.hallmark_evidence,
        hallmark_witnesses=state.hallmark_witnesses,
        host_signature_model=model,
        host_signatures=host_signatures,
        host_deviation_summary=state.host_deviation_summary,
        elapsed=0.0,
    )


def _rebuild_host_model(*, output_dir: Path, config: PipelineConfig, masked_path: Path, proteome_path: Path, logger):
    labels = config.databases.taxonomy_labels_file
    lookup = TaxonomyLabelLookup.load(Path(labels)) if labels and Path(labels).exists() else None
    set_taxonomy_lookup(lookup)
    phase1 = config.phase1
    records = rebuild_host_records(
        host_query_fasta=output_dir / _HOST_FASTA,
        host_diamond_tsv=output_dir / _HOST_TSV,
        proteome_path=proteome_path,
        genome_path=masked_path,
        taxonomy_lookup=lookup,
        taxonomy_weight_mode=config.phase2.taxonomy_weight_mode,
        novel_criteria=NovelMarkerCriteria(
            min_hmm_score=phase1.novel_marker_min_score,
            min_hmm_coverage=phase1.novel_marker_min_coverage,
            require_cluster=phase1.novel_marker_require_cluster,
            initial_window_bp=phase1.initial_window_bp,
        ),
    )
    return _build_host_model(
        records,
        host_prefixes=config.host.prefixes,
        min_token_length=config.phase2.host_signature_min_token_len,
        weight_mode=config.phase2.taxonomy_weight_mode,
        output_dir=output_dir,
        logger=logger,
    )


def rebuild_host_records(
    *,
    host_query_fasta: Path,
    host_diamond_tsv: Path,
    proteome_path: Path,
    genome_path: Path,
    taxonomy_lookup: dict | None,
    taxonomy_weight_mode: str,
    novel_criteria: NovelMarkerCriteria,
) -> list[ValidatedMarkerHit]:
    """Rebuild the host records of the Phase 1 host search from its saved files.

    The host signature model reads, for each host-searched protein, its scaffold and
    its ten DIAMOND matches. ``host_hmm_hit_porfs.faa`` lists the proteins once, in
    the order of their first selected hit, which is the order of the search records.
    ``host_diamond_top10.tsv`` holds the matches. The scaffold comes from the
    proteome. One placeholder hit per protein therefore yields records with the same
    model-relevant content as the search records, built by the same function.
    """
    for path in (host_query_fasta, host_diamond_tsv):
        if not Path(path).is_file():
            raise Phase1ReuseError(f"cannot rebuild the host records: {Path(path).name} is missing")
    hits = [
        HMMHit(
            query_name=record.id,
            target_name=_HOST_QUERY_MODEL,
            score=0.0,
            evalue=1.0,
            domain_score=0.0,
            query_start=1,
            query_end=len(record.seq),
        )
        for record in SeqIO.parse(host_query_fasta, "fasta")
    ]
    return filter_validated_markers(
        hmm_hits=hits,
        diamond_output_tsv=Path(host_diamond_tsv),
        proteome_fasta=Path(proteome_path),
        genome_fasta=Path(genome_path),
        output_dir=None,
        novel_criteria=novel_criteria,
        taxonomy_lookup=taxonomy_lookup,
        taxonomy_weight_mode=taxonomy_weight_mode,
        max_seqs=10,
    )


def check_rescue_host_search(rescue_host_faa: Path, confirmed_ids: Collection[str]) -> None:
    """Stop when a confirmed frameshift-rescue protein is in the rescue host search.

    The search excludes proteins with viral evidence, so a confirmed protein adds no
    host record. A record that did would need hits that this mode does not save.
    """
    if not confirmed_ids:
        return
    if not Path(rescue_host_faa).is_file():
        raise Phase1ReuseError(f"cannot rebuild the host records: {Path(rescue_host_faa).name} is missing")
    present = sorted({record.id for record in SeqIO.parse(rescue_host_faa, "fasta")} & set(confirmed_ids))
    if present:
        raise Phase1ReuseError(
            "cannot rebuild the host records exactly: confirmed rescue proteins are in the rescue host search: "
            + ", ".join(present)
        )


def _copy_tree(source: ReusedRun, output_dir: Path, top: str, *, skip: frozenset[str]) -> list[dict[str, object]]:
    """Copy only authenticated regular files below ``top`` and check their saved hashes."""
    source_root = source.root
    if (source_root / top).is_symlink():
        raise Phase1ReuseError(f"the earlier run contains a symlink root: {top}")
    copied: list[dict[str, object]] = []
    pending = [Path(top)]
    while pending:
        directory = pending.pop()
        (Path(output_dir) / directory).mkdir(parents=True, exist_ok=True)
        with os.scandir(Path(source_root) / directory) as entries:
            for entry in sorted(entries, key=lambda item: item.name):
                relative = directory / entry.name
                if entry.is_symlink():
                    raise Phase1ReuseError(f"the earlier run contains a symlink: {relative}")
                if entry.is_dir():
                    pending.append(relative)
                elif not entry.is_file():
                    raise Phase1ReuseError(f"the earlier run contains a special file: {relative}")
                elif relative.as_posix() not in skip and relative.as_posix() in source.artifacts:
                    copied.append(
                        _copy_file(
                            Path(entry.path),
                            Path(output_dir) / relative,
                            relative,
                            source.artifacts[relative.as_posix()],
                        )
                    )
    return sorted(copied, key=lambda item: str(item["path"]))


def _copy_file(source: Path, destination: Path, relative: Path, expected_sha256: str) -> dict[str, object]:
    digest = hashlib.sha256()
    size = 0
    with source.open("rb") as reader, destination.open("xb") as writer:
        for chunk in iter(lambda: reader.read(_COPY_CHUNK_BYTES), b""):
            digest.update(chunk)
            writer.write(chunk)
            size += len(chunk)
    if digest.hexdigest() != expected_sha256:
        raise Phase1ReuseError(f"authenticated input changed during copy: {relative}")
    return {"path": relative.as_posix(), "size": size, "sha256": digest.hexdigest()}


def _rewritten(source: ReusedRun, output_dir: Path, relatives: Iterable[str]) -> list[dict[str, object]]:
    entries: list[dict[str, object]] = []
    for relative in relatives:
        earlier, current = _sha256_or_none(source.root / relative), _sha256_or_none(Path(output_dir) / relative)
        if earlier is not None or current is not None:
            entries.append({"path": relative, "old_sha256": earlier, "sha256": current})
    return entries


def _sha256_or_none(path: Path) -> str | None:
    if not path.is_file():
        return None
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _write_provenance(
    source: ReusedRun,
    output_dir: Path,
    phase: int,
    copied: list[dict[str, object]],
    rewritten: list[dict[str, object]],
) -> None:
    document = {
        "schema": PROVENANCE_SCHEMA,
        "phase": phase,
        "old_run_folder": str(source.root.resolve()),
        "old_run_fingerprint": source.run_fingerprint,
        "old_code_sha256": source.code_sha256,
        "old_phase_marker_sha256": source.markers[phase],
        "copied": copied,
        "authentication": "old phase artifact hashes; unauthenticated files were not copied",
        "parent_reuse_provenance": {
            relative: sha for relative, sha in source.artifacts.items() if relative.endswith("/" + PROVENANCE_FILENAME)
        },
        "rewritten": rewritten,
    }
    atomic_write(
        Path(output_dir) / f"phase{phase}" / PROVENANCE_FILENAME,
        json.dumps(document, indent=2, sort_keys=True) + "\n",
    )
