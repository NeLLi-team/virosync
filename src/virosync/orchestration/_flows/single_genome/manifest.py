"""Run-log, completion-manifest, and config-fingerprint helpers.

Used by the phase subflows (early-exit paths) and the orchestrator to write the
resume markers and summarize prediction tables. Depends on ``loaders._safe_int``
but never on the phase modules or orchestrator.
"""

import csv
import hashlib
import json
from dataclasses import fields, is_dataclass
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Optional

import virosync.output_contract as output_contract
from virosync.utils.atomic_write import atomic_write_context
from virosync.validation.tsv_invariants import enforce_tsv_invariants

from .loaders import _safe_int

# Manifest schema: bumped 1 -> 2 for fingerprint-strict resume (Phase A).
# Schema-v1 compatibility readers remain available to explicit migration tools
# and regression tests, but normal orchestration never enables their reuse.
_MANIFEST_SCHEMA_VERSION = 2


def _clear_success_markers(output_dir: Path) -> None:
    """Remove stale completion signals without touching diagnostic artifacts."""
    output_dir = Path(output_dir)
    for name in ("run.log", "virosync_run_complete.json"):
        (output_dir / name).unlink(missing_ok=True)


def _empty_prediction_summary() -> dict[str, int]:
    """Return the stable numeric result surface for a zero-prediction run."""
    return {
        "predictions": 0,
        "accepted": 0,
        "high_tier": 0,
        "medium_tier": 0,
        "low_tier": 0,
        "candidate_high_tier": 0,
        "candidate_medium_tier": 0,
        "candidate_low_tier": 0,
        "accepted_bp": 0,
        "total_genes": 0,
        "total_hallmarks": 0,
        **output_contract.empty_effective_eve_class_counts(),
        "quality_gate_dropped": 0,
    }


def _write_empty_run_log(
    output_dir: Path,
    genome_id: str,
    reason: str,
    elapsed_sec: float,
    input_path: Optional[Path] = None,
    output_files: Optional[dict] = None,
    fingerprint: Optional[str] = None,
) -> Path:
    """Write a minimal run.log so resume treats a zero-result run as complete.

    Phase 1 early exits (no HMM hits, no validated markers, no seeds) and the
    Phase 2 "no refined boundaries" exit used to return without writing
    ``run.log``, so any subsequent resume would redo the whole pipeline for
    the same genome. Emitting a compact completion marker here lets batch
    orchestration skip these genomes on restart.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    _clear_success_markers(output_dir)

    invariant_report_path = output_dir / "virosync_tsv_invariant_report.tsv"
    invariant_report = enforce_tsv_invariants(
        detailed_tsv=output_dir / "virosync_predictions_detailed.tsv",
        report_out=invariant_report_path,
    )
    if output_files is None:
        output_files = {}
    output_files["tsv_invariant_report"] = str(invariant_report_path)

    run_log_path = output_dir / "run.log"
    with atomic_write_context(run_log_path, "w") as f:
        f.write(f"# ViroSync Run Log: {genome_id}\n")
        f.write(f"# Generated: {datetime.now(timezone.utc).isoformat()}\n")
        if input_path is not None:
            f.write(f"# Input: {input_path}\n")
        f.write(f"# Output: {output_dir}\n")
        f.write(f"# Total time: {elapsed_sec:.1f}s\n")
        f.write("#" + "=" * 59 + "\n\n")
        f.write("## Results Summary\n")
        f.write("Seeds identified: 0\n")
        f.write("Boundaries refined: 0\n")
        f.write("EVEs detected: 0\n")
        f.write(f"Early-exit reason: {reason}\n")
        f.write("\n## Detailed TSV Invariant Check\n")
        f.write(f"Status: {invariant_report.status}\n")
        f.write(f"Rows checked: {invariant_report.rows_checked}\n")
        f.write(f"Issues: {invariant_report.issue_count}\n")
        f.write(f"Errors: {invariant_report.error_count}\n")
        f.write(f"Warnings: {invariant_report.warning_count}\n")
        f.write(f"Report: {invariant_report_path}\n")
    _write_completion_manifest(
        output_dir=output_dir,
        genome_id=genome_id,
        status="success",
        reason=reason,
        output_files=output_files,
        fingerprint=fingerprint,
    )
    return run_log_path


# --- Schema-v3 run-identity field partition ---------------------------------
#
# Schema v3 keeps scalar configuration, resources, execution environment, and
# runtime-only controls as separate identities.  The partition is hard-coded so
# a new PipelineConfig field fails the drift guard until it is classified.
#
# Names are the FLAT kwargs of ``_single_genome_flow_impl`` -- i.e. the *effective*
# config it actually runs (already resolved from YAML + CLI). Fingerprinting those
# directly avoids the ``PipelineConfig.with_overrides`` round-trip blind spot for
# knobs that are not overridable via flat kwargs (boundary_diamond_*), which
# would otherwise read back as defaults.

# Resource paths are excluded from the scalar config digest.  ``run_state`` binds
# authenticated core manifests and deterministic content manifests instead.
_FINGERPRINT_RESOURCE_FIELDS = (
    "hmm_database",
    "hmm_allowlist",
    "marker_db",
    "gene_taxonomy_faa_db",
    "taxonomy_labels_file",
    # marker source build-inputs: when rebuild_db builds marker_db from these, the
    # built DB identity is computed BEFORE the build, so the sources must be hashed too.
    "marker_faa_db",
    "marker_faa_dir",
    "faa_dir",
)
_FINGERPRINT_RESOURCE_GATED = (
    ("gvclass_db", ("enable_phylogenetic", "run_gvclass")),
    ("diamond_db", ("enable_phylogenetic",)),
    ("viral_structure_db", ("use_boltz",)),
    ("tmvec_database_dir", ("use_tmvec_database",)),
    ("interproscan_dir", ("interproscan_enabled",)),
    ("gvclass_path", ("run_gvclass",)),
)
# Scalar output-determining fields -> hashed by canonical value (enums already arrive
# as their ``.value`` strings; lists are sorted; None -> "none").
_FINGERPRINT_SCALAR_FIELDS = (
    "ablation_id",
    "ablation_contract_sha256",
    "seed_marker_allowlist",
    # compute
    "search_backend",
    # phase1
    "rebuild_db",
    "assembly_mode",
    "hmm_chunk_size",
    "initial_window_bp",
    "initial_window_genes",
    "min_markers_initial",
    "extension_kb",
    "merge_distance",
    "host_taxonomy_deviation_enabled",
    "host_taxonomy_deviation_allow_seeds",
    "host_taxonomy_deviation_min_token_len",
    "host_taxonomy_deviation_min_tokens",
    "host_taxonomy_deviation_overlap_threshold",
    "host_taxonomy_deviation_max_pident",
    "host_taxonomy_deviation_max_hits",
    "host_taxonomy_deviation_window_bp",
    "host_taxonomy_deviation_window_count",
    "host_taxonomy_deviation_window_seed",
    "host_taxonomy_deviation_window_min_markers",
    "host_taxonomy_deviation_seed_window_bp",
    "host_taxonomy_deviation_seed_min_markers",
    "marker_validation_top_k",
    "novel_marker_min_score",
    "novel_marker_min_coverage",
    "novel_marker_require_cluster",
    # phase2
    "taxonomy_weight_mode",
    "boundary_taxonomy_ml_enabled",
    "boundary_taxonomy_ml_model",
    "boundary_taxonomy_ml_threshold",
    "boundary_taxonomy_ml_neighbor_window",
    "boundary_host_trim_enabled",
    "boundary_host_trim_window_bp",
    "boundary_host_trim_step_bp",
    "boundary_host_trim_max_host_fraction",
    "boundary_host_trim_min_viral_fraction",
    "boundary_host_trim_score_threshold",
    "boundary_host_trim_buffer_kb",
    "boundary_host_trim_min_overlap_score",
    "boundary_host_signature_min_token_len",
    "boundary_diamond_flank_genes",
    "boundary_diamond_control_sample_size",
    "boundary_diamond_control_min_distance",
    "boundary_diamond_top_k",
    "boundary_diamond_chunk_size",
    "boundary_diamond_random_seed",
    "boundary_diamond_superset_prototype_enabled",
    # phase3
    "high_tier_threshold",
    "low_tier_threshold",
    "use_crf_in_final_score",
    "priority_marker_list",
    "marker_floor_priority_only",
    "marker_floor_priority_plus_family",
    "marker_floor_priority_multi_family",
    "marker_family_bonus_per_family",
    "marker_multi_family_bonus",
    "enable_phylogenetic",
    "skip_structural",
    "use_boltz",
    "boltz_mcp_only",
    "boltz_use_msa_server",
    "boltz_min_seq_len",
    "boltz_max_seq_len",
    "boltz_no_kernels",
    "use_tmvec_database",
    "tmvec_require_gpu",
    "tmvec_databases",
    "tmvec_min_score",
    "interproscan_enabled",
    "interproscan_applications",
    "interproscan_keywords",
    "export_all_eve_sequences",
    "extended_output",
    "run_gvclass",
    "host_signature_evidence_threshold",
    # host
    "host_prefixes",
    "host_label",
    "high_pident_host_threshold",
    # execution
    "masking",
)
_FINGERPRINT_CONFIG_FIELDS = frozenset(_FINGERPRINT_SCALAR_FIELDS)
_FINGERPRINT_RESOURCE_FIELDS = frozenset(
    _FINGERPRINT_RESOURCE_FIELDS
    + tuple(field for field, _gate in _FINGERPRINT_RESOURCE_GATED)
)
_FINGERPRINT_ENVIRONMENT_FIELDS = frozenset({"device"})

# Flat fields deliberately EXCLUDED from the fingerprint: runtime/IO-only knobs that
# never change which EVEs are accepted. Kept explicit so the drift-guard test can prove
# every PipelineConfig FIELD_SPEC is classified (fingerprinted XOR excluded) -- a newly
# added output-determining knob then fails that test until it is triaged here.
_FINGERPRINT_RUNTIME_ONLY_FIELDS = frozenset(
    {
        # compute / parallelism (never output-determining)
        "threads",
        "max_threads",
        "gene_taxonomy_threads",
        "interproscan_threads",
        # resume control (does not change accepted EVEs)
        "resume",
    }
)

# Flat locals needed to build every schema-v3 identity.  This compatibility
# alias remains the orchestrator's explicit locals filter.
_FINGERPRINT_INPUT_FIELDS = frozenset(
    _FINGERPRINT_CONFIG_FIELDS
    | _FINGERPRINT_RESOURCE_FIELDS
    | _FINGERPRINT_ENVIRONMENT_FIELDS
)


def _canonical_config_value(value):
    """Return a typed, unambiguous JSON value for config fingerprinting."""

    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, Enum):
        return _canonical_config_value(value.value)
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value):
        payload = {
            item.name: getattr(value, item.name)
            for item in fields(value)
        }
        library = payload.get("repeatmasker_library")
        if library is not None:
            library_path = Path(library)
            payload["repeatmasker_library"] = {
                "path": str(library_path),
                "sha256": (
                    hashlib.sha256(library_path.read_bytes()).hexdigest()
                    if library_path.is_file()
                    else "missing"
                ),
            }
        return _canonical_config_value(payload)
    if isinstance(value, dict):
        return {
            str(key): _canonical_config_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_canonical_config_value(item) for item in value]
    raise TypeError(f"unsupported config fingerprint value: {type(value)!r}")


def _compute_config_fingerprint(flat_config: dict) -> str:
    """Return the full SHA-256 of effective scalar output configuration.

    Resources and device/environment are deliberately absent: schema-v3 run
    identity binds those categories independently. Runtime controls such as
    thread count and resume mode never affect this digest.
    """
    document = {
        "schema_version": 1,
        "fields": {
            field: _canonical_config_value(flat_config.get(field))
            for field in sorted(_FINGERPRINT_SCALAR_FIELDS)
        },
    }
    encoded = json.dumps(
        document,
        allow_nan=False,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _masking_status_identity(output_dir: Path) -> dict | None:
    """Return a verified status-file identity for completion/provenance binding."""
    from virosync.pipeline.phase0.masking import load_masking_result

    output_dir = Path(output_dir)
    status_path = output_dir / "phase0" / "masking" / "masking_status.json"
    if not status_path.is_file():
        return None
    result = load_masking_result(status_path)
    payload = json.loads(status_path.read_text(encoding="utf-8"))
    return {
        "path": "phase0/masking/masking_status.json",
        "sha256": result.status_sha256,
        "result_fingerprint": payload["result_fingerprint"],
        "benchmark_eligible": bool(payload["benchmark_eligible"]),
        "status": payload["status"],
    }


def _effective_masking_fingerprint(
    requested_fingerprint: str,
    status_sha256: str,
) -> str:
    """Bind the pre-run requested identity to the post-mask status file."""
    return hashlib.sha256(
        f"{requested_fingerprint}|{status_sha256}".encode()
    ).hexdigest()


def _write_completion_manifest(
    output_dir: Path,
    genome_id: str,
    status: str,
    *,
    reason: str | None = None,
    output_files: dict | None = None,
    fingerprint: str | None = None,
) -> Path:
    """Atomically write the final resume completion marker."""
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "virosync_run_complete.json"
    payload = {
        "schema_version": _MANIFEST_SCHEMA_VERSION,
        "genome_id": genome_id,
        "status": status,
        "reason": reason or "",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        **output_contract.coordinate_contract_metadata(),
    }
    if output_files is not None:
        payload["output_files"] = _json_safe(output_files)
    if fingerprint is not None:
        payload["config_fingerprint"] = fingerprint
    masking_status = None
    if status == "success":
        try:
            masking_status = _masking_status_identity(output_dir)
        except (KeyError, OSError, TypeError, ValueError) as exc:
            _clear_success_markers(output_dir)
            raise ValueError(
                "cannot write success completion manifest without a valid "
                "phase0/masking/masking_status.json"
            ) from exc
        if masking_status is None:
            _clear_success_markers(output_dir)
            raise ValueError(
                "cannot write success completion manifest without a valid "
                "phase0/masking/masking_status.json"
            )
    if masking_status is not None:
        payload["masking_status"] = masking_status
        if fingerprint is not None:
            payload["effective_masking_fingerprint"] = (
                _effective_masking_fingerprint(
                    fingerprint,
                    masking_status["sha256"],
                )
            )
    with atomic_write_context(manifest_path, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest_path


def _json_safe(value):
    """Convert nested values used in run manifests to JSON-safe scalars."""
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def _summarize_predictions_tsv(
    predictions_tsv: Path,
    *,
    canonical: bool = True,
) -> dict[str, int]:
    """Summarize prediction rows.

    Canonical ``virosync_predictions.tsv`` rows have already been selected by
    the active acceptance policy, including any canonical LOW-confidence rows.
    Detailed TSV rows are all candidates and should be summarized with
    ``canonical=False``.
    """
    stats = _empty_prediction_summary()
    with predictions_tsv.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or [])
        has_persisted_class = "effective_eve_class" in fieldnames
        gene_count_field = (
            "gene_taxonomy_total"
            if "gene_taxonomy_total" in fieldnames
            else "total_proteins"
            if "total_proteins" in fieldnames
            else None
        )
        hallmark_field = (
            "hallmark_total"
            if "hallmark_total" in fieldnames
            else "hallmark_count"
            if "hallmark_count" in fieldnames
            else None
        )
        for row in reader:
            stats["predictions"] += 1
            tier = (row.get("confidence_tier") or "").upper()
            if tier == "HIGH":
                stats["high_tier"] += 1
            elif tier == "MEDIUM":
                stats["medium_tier"] += 1
            elif tier == "LOW":
                stats["low_tier"] += 1
            if canonical:
                stats["accepted"] += 1
                stats["accepted_bp"] += _safe_int(row.get("length"))
                if gene_count_field is not None:
                    stats["total_genes"] += _safe_int(row.get(gene_count_field))
                if hallmark_field is not None:
                    stats["total_hallmarks"] += _safe_int(row.get(hallmark_field))
                if has_persisted_class:
                    effective_class = output_contract.normalize_effective_eve_class(
                        row.get("effective_eve_class")
                    )
                else:
                    # The resolver answers in the GATE vocabulary, which still
                    # carries MIXED, so a legacy tree without a persisted class
                    # needs the published fold before it can index the partition.
                    effective_class = output_contract.normalize_effective_eve_class(
                        output_contract.resolve_effective_eve_class(
                            confidence_tier=tier,
                            region_classification=row.get("region_classification"),
                            classification=row.get("classification"),
                            likely_family=row.get("likely_family"),
                        )
                    )
                class_key = output_contract.EFFECTIVE_EVE_CLASS_COUNT_KEYS[
                    effective_class
                ]
                stats[class_key] += 1
    if canonical and output_contract.effective_eve_class_count_total(stats) != stats["accepted"]:
        raise ValueError(
            "exclusive effective-class counts do not sum to accepted predictions"
        )
    return stats


def _summarize_prediction_outputs(
    canonical_tsv: Path,
    detailed_tsv: Path,
    *,
    expected_accepted: int | None = None,
    expected_candidates: int | None = None,
) -> dict[str, int]:
    """Build the public result counts from persisted canonical and detailed TSVs."""
    summary = _summarize_predictions_tsv(canonical_tsv)
    candidates = _summarize_predictions_tsv(detailed_tsv, canonical=False)

    accepted = summary["accepted"]
    candidate_count = candidates["predictions"]
    if candidate_count < accepted:
        raise ValueError(
            "persisted candidate count is smaller than accepted count: "
            f"accepted={accepted} candidates={candidate_count}"
        )
    if expected_accepted is not None and accepted != expected_accepted:
        raise ValueError(
            "persisted accepted count disagrees with the in-memory gate: "
            f"expected={expected_accepted} persisted={accepted}"
        )
    if expected_candidates is not None and candidate_count != expected_candidates:
        raise ValueError(
            "persisted candidate count disagrees with in-memory Phase 3: "
            f"expected={expected_candidates} persisted={candidate_count}"
        )

    summary["predictions"] = candidate_count
    summary["candidate_high_tier"] = candidates["high_tier"]
    summary["candidate_medium_tier"] = candidates["medium_tier"]
    summary["candidate_low_tier"] = candidates["low_tier"]
    summary["quality_gate_dropped"] = candidate_count - accepted
    return summary
