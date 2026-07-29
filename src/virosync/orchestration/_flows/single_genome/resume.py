"""Resume short-circuit detection and artifact validation helpers.

Self-contained: validates run logs, completion manifests, and prediction TSVs
to decide whether a prior run can be reused. No sibling-module dependencies.
"""

import csv
import hashlib
import json
from pathlib import Path
from typing import Optional

from .run_state import (
    RUN_STATE_FILENAME,
    load_run_state,
    plan_resume,
)


def _completed_run_artifacts(
    output_dir: Path,
    expected_fingerprint: str | None = None,
    *,
    allow_missing_fingerprint: bool = False,
    allow_legacy_schema: bool = False,
    expected_input: Path | None = None,
) -> dict[str, Path] | None:
    """Return the exact artifact set for a safe full-resume short-circuit.

    Schema-v3 run state is authoritative. Schema-v1/v2 completion manifests are
    considered only under the explicit legacy opt-in; they are never promoted to
    schema-v3 state or benchmark eligibility.
    """
    output_dir = Path(output_dir)
    state_path = output_dir / RUN_STATE_FILENAME
    if state_path.is_file():
        try:
            state = load_run_state(output_dir)
            fingerprint = expected_fingerprint or state.run_fingerprint
            plan = plan_resume(
                output_dir,
                expected_run_fingerprint=fingerprint,
            )
        except (OSError, TypeError, ValueError):
            return None
        if not plan.completed or state.status != "success":
            return None

        paths = {
            artifact.relative_path: output_dir / artifact.relative_path
            for artifact in state.artifacts
        }

        def _select(*relative_paths: str) -> Path | None:
            for relative in relative_paths:
                candidate = paths.get(relative)
                if candidate is not None:
                    return candidate
            return None

        canonical = _select(
            "phase3_synthesis/virosync_predictions.tsv",
            "virosync_predictions.tsv",
        )
        detailed = _select(
            "virosync_predictions_detailed.tsv",
            "phase3_synthesis/virosync_predictions_detailed.tsv",
        )
        if canonical is None or detailed is None:
            return None
        artifacts = {
            "phase3_predictions": canonical,
            "predictions_detailed": detailed,
            "run_state": state_path,
        }
        run_log = paths.get("run.log")
        if run_log is not None:
            artifacts["run_log"] = run_log
        compatibility = paths.get("virosync_run_complete.json")
        if compatibility is not None:
            artifacts["completion_manifest"] = compatibility
        masking_status = output_dir / "phase0" / "masking" / "masking_status.json"
        if masking_status.is_file():
            artifacts["masking_status"] = masking_status
        return artifacts

    if not allow_legacy_schema:
        return None

    run_log = output_dir / "run.log"
    manifest = output_dir / "virosync_run_complete.json"
    phase3_predictions = output_dir / "phase3_synthesis" / "virosync_predictions.tsv"
    root_predictions = output_dir / "virosync_predictions.tsv"
    root_detailed = output_dir / "virosync_predictions_detailed.tsv"
    phase3_detailed = output_dir / "phase3_synthesis" / "virosync_predictions_detailed.tsv"

    if not _valid_resume_run_log(run_log):
        return None
    if not _valid_completion_manifest(
        manifest,
        # Explicit legacy reuse cannot prove the schema-v3 fingerprint. The
        # compatibility manifest's own masking/config binding is still checked.
        expected_fingerprint=None,
        allow_missing_fingerprint=allow_missing_fingerprint,
        expected_input=expected_input,
    ):
        return None

    predictions = _first_valid_tsv(
        [phase3_predictions, root_predictions],
        required_fields={"eve_id"},
    )
    detailed = _first_valid_tsv(
        [root_detailed, phase3_detailed],
        required_fields={"eve_id"},
    )
    if predictions is None or detailed is None:
        return None

    artifacts = {
        "phase3_predictions": predictions,
        "predictions_detailed": detailed,
        "run_log": run_log,
        "completion_manifest": manifest,
    }
    payload = _load_completion_payload(manifest) or {}
    if payload.get("masking_status"):
        artifacts["masking_status"] = (
            output_dir / "phase0" / "masking" / "masking_status.json"
        )
    return artifacts


def _first_valid_tsv(
    paths: list[Path],
    *,
    required_fields: set[str],
) -> Path | None:
    for path in paths:
        if _valid_tsv_header(path, required_fields=required_fields):
            return path
    return None


def _valid_tsv_header(path: Path, *, required_fields: set[str]) -> bool:
    if not path.exists() or not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fieldnames = set(reader.fieldnames or [])
    except OSError:
        return False
    return required_fields.issubset(fieldnames)


def _valid_resume_run_log(path: Path) -> bool:
    if not path.exists() or not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return False
    return text.startswith("# ViroSync Run Log:") and "## Results Summary" in text


# Manifest schema versions this compatibility reader understands. Production
# orchestration never enables schema-v1/v2 reuse; the legacy branch remains for
# explicit migration tooling and regression tests only.
_READABLE_SCHEMA_VERSIONS = frozenset({1, 2})


def _load_completion_payload(path: Path) -> dict | None:
    """Parse a STRUCTURALLY valid completion manifest, else None.

    Structural validity = readable JSON with a known schema_version, a success status,
    and a genome_id. Deliberately ignores the config fingerprint (see _fingerprint_ok)
    so callers can separate "is this a real ViroSync completion" from "does it match
    the current config".
    """
    if not path.exists() or not path.is_file() or path.stat().st_size == 0:
        return None
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if (
        payload.get("schema_version") in _READABLE_SCHEMA_VERSIONS
        and payload.get("status") == "success"
        and bool(payload.get("genome_id"))
    ):
        return payload
    return None


def _fingerprint_ok(
    payload: dict,
    expected_fingerprint: str | None,
    allow_missing_fingerprint: bool,
) -> bool:
    """Fingerprint-strict resume gate (Phase A):
    - no expected fingerprint    -> caller does not enforce it (accept)
    - manifest lacks fingerprint -> accept only under the legacy opt-in
    - both present               -> accept iff equal
    """
    manifest_fp = payload.get("config_fingerprint")
    if expected_fingerprint is None:
        return True
    if manifest_fp is None:
        return allow_missing_fingerprint
    return manifest_fp == expected_fingerprint


def _masking_identity_ok(
    payload: dict,
    output_dir: Path,
    *,
    allow_missing_identity: bool,
    expected_input: Path | None = None,
) -> bool:
    """Verify a recorded whole-file status and its semantic result identity."""
    identity = payload.get("masking_status")
    if identity is None:
        return allow_missing_identity
    if not isinstance(identity, dict):
        return False
    if identity.get("path") != "phase0/masking/masking_status.json":
        return False
    status_path = output_dir / "phase0" / "masking" / "masking_status.json"
    try:
        from virosync.pipeline.phase0.masking import load_masking_result

        result = load_masking_result(status_path, expected_input=expected_input)
        status_payload = json.loads(status_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return False
    if identity.get("sha256") != result.status_sha256:
        return False
    if identity.get("result_fingerprint") != status_payload.get("result_fingerprint"):
        return False
    if identity.get("benchmark_eligible") != status_payload.get("benchmark_eligible"):
        return False
    if identity.get("status") != status_payload.get("status"):
        return False
    requested_fingerprint = payload.get("config_fingerprint")
    effective_fingerprint = payload.get("effective_masking_fingerprint")
    if requested_fingerprint is not None:
        expected_effective = hashlib.sha256(
            f"{requested_fingerprint}|{result.status_sha256}".encode()
        ).hexdigest()
        if effective_fingerprint != expected_effective:
            return False
    return True


def _valid_completion_manifest(
    path: Path,
    expected_fingerprint: str | None = None,
    *,
    allow_missing_fingerprint: bool = False,
    expected_input: Path | None = None,
) -> bool:
    payload = _load_completion_payload(path)
    if payload is None:
        return False
    return _fingerprint_ok(
        payload,
        expected_fingerprint,
        allow_missing_fingerprint,
    ) and _masking_identity_ok(
        payload,
        path.parent,
        allow_missing_identity=allow_missing_fingerprint,
        expected_input=expected_input,
    )


def _manifest_is_stale(
    output_dir: Path,
    *,
    expected_fingerprint: str | None,
    allow_missing_fingerprint: bool = False,
    expected_input: Path | None = None,
) -> bool:
    """True iff a prior run COMPLETED here under a different effective config.

    A structurally valid completion manifest that FAILS the fingerprint gate means the
    per-genome output (and its phase-level intermediate caches) was produced under a
    different config/database and must NOT be reused -- the file-existence resume in
    Phase 0-3 would otherwise fold stale intermediates into a "fresh" result. Returns
    False when no/absent/corrupt manifest is present (an interrupted run with no
    completion marker is still legitimately resumable).
    """
    output_dir = Path(output_dir)
    state_path = output_dir / RUN_STATE_FILENAME
    if state_path.is_file():
        try:
            state = load_run_state(output_dir)
            fingerprint = expected_fingerprint or state.run_fingerprint
            return not plan_resume(
                output_dir,
                expected_run_fingerprint=fingerprint,
            ).completed
        except (OSError, TypeError, ValueError):
            return True

    payload = _load_completion_payload(output_dir / "virosync_run_complete.json")
    if payload is None:
        return False
    if not allow_missing_fingerprint:
        return True
    return not (
        _fingerprint_ok(payload, expected_fingerprint, allow_missing_fingerprint)
        and _masking_identity_ok(
            payload,
            output_dir,
            allow_missing_identity=allow_missing_fingerprint,
            expected_input=expected_input,
        )
    )


def _require_phase2b_gene_taxonomy_db(
    gene_taxonomy_faa_db: Optional[Path],
    *,
    has_seeds: bool,
) -> Optional[Path]:
    """Return the required Phase 2b taxonomy DB, or fail before boundary work."""
    if not has_seeds:
        return None
    if gene_taxonomy_faa_db is None:
        raise ValueError(
            "Phase 2b requires gene_taxonomy_faa_db. "
            "The marker validation database is not a valid fallback for "
            "boundary taxonomy refinement."
        )
    phase2b_db = Path(gene_taxonomy_faa_db)
    if not phase2b_db.exists():
        raise FileNotFoundError(
            f"Phase 2b gene taxonomy database not found: {phase2b_db}"
        )
    return phase2b_db
