"""Plain-Python batch orchestration for ViroSync."""

from __future__ import annotations

import csv
import io
import json
import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from virosync.config import PipelineConfig
from virosync.orchestration._flows.single_genome import single_genome_flow
from virosync.output_contract import effective_eve_class_count_total
from virosync.utils.atomic_write import atomic_write
from virosync.utils.path_safety import require_strict_child, validate_path_component
from virosync.validation.tsv_invariants import TSVInvariantError

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class GenomeRunSpec:
    """Preflighted input-to-output mapping for one genome run."""

    input_path: Path
    genome_id: str
    output_dir: Path


def _exclusive_class_count_violation(result: dict) -> str | None:
    """Return why a successful result's class partition is not exhaustive."""
    if not result.get("success", False):
        return None
    accepted = int(result.get("accepted", 0) or 0)
    classified = effective_eve_class_count_total(result)
    if classified == accepted:
        return None
    return (
        "exclusive effective-class counts do not sum to accepted predictions "
        f"for {result.get('genome_id', '?')}: "
        f"accepted={accepted} classified={classified}"
    )


def _enforce_exclusive_class_counts(results: list[dict]) -> list[dict]:
    """Demote results whose class partition is not exhaustive to explicit failures.

    One corrupt genome must never suppress the whole batch's summary, and it must
    never be summarized as a credible success either, so it is reported as a
    failed row carrying the invariant violation.
    """
    enforced: list[dict] = []
    for result in results:
        violation = _exclusive_class_count_violation(result)
        if violation is None:
            enforced.append(result)
            continue
        logger.error("%s", violation)
        enforced.append(
            _failure_result(result.get("genome_id", "?"), ValueError(violation))
        )
    return enforced


def _preflight_genome_runs(
    genome_paths: Iterable[Path],
    output_base_dir: Path,
) -> list[GenomeRunSpec]:
    """Validate the complete batch mapping before any worker or output is created."""
    output_base_dir = Path(output_base_dir)
    specs: list[GenomeRunSpec] = []
    invalid: list[str] = []
    sources_by_id: dict[str, list[Path]] = {}

    for raw_path in genome_paths:
        input_path = Path(raw_path)
        genome_id = input_path.stem
        try:
            validate_path_component(genome_id, "genome ID")
            output_dir = output_base_dir / genome_id
            if ".." in output_dir.parts:
                raise ValueError(f"output path contains a parent segment: {output_dir}")
            require_strict_child(output_base_dir, output_dir)
            if output_dir.is_symlink():
                raise ValueError(f"output path is a symlink: {output_dir}")
        except (TypeError, ValueError) as exc:
            invalid.append(f"- {input_path}: {exc}")
            continue

        specs.append(
            GenomeRunSpec(
                input_path=input_path,
                genome_id=genome_id,
                output_dir=output_dir,
            )
        )
        sources_by_id.setdefault(genome_id, []).append(input_path)

    conflicts = {
        genome_id: source_paths
        for genome_id, source_paths in sources_by_id.items()
        if len(source_paths) > 1
    }
    if invalid or conflicts:
        lines = ["Unsafe or ambiguous genome input mapping:"]
        lines.extend(invalid)
        for genome_id, source_paths in sorted(conflicts.items()):
            lines.append(f"- duplicate genome ID {genome_id!r}:")
            lines.extend(f"  - {source_path}" for source_path in source_paths)
        raise ValueError("\n".join(lines))

    return specs


def _single_genome_callable():
    return getattr(single_genome_flow, "fn", single_genome_flow)


def _failure_result(genome_id: str, exc: BaseException) -> dict:
    return {
        "genome_id": genome_id,
        "success": False,
        "benchmark_eligible": False,
        "legacy_resume": False,
        "error": str(exc),
        "predictions": 0,
        "accepted": 0,
        "elapsed_sec": 0,
    }


def _normalize_worker_result(genome_id: str, result: object) -> dict:
    """Return one well-formed result record for every requested input."""
    if not isinstance(result, dict):
        return _failure_result(
            genome_id,
            TypeError(
                "single-genome flow returned "
                f"{type(result).__name__}, expected a result mapping"
            ),
        )
    normalized = dict(result)
    normalized.setdefault("genome_id", genome_id)
    normalized.setdefault("benchmark_eligible", False)
    normalized.setdefault("legacy_resume", False)
    return normalized


def _is_benchmark_eligible(result: dict) -> bool:
    """Return whether a successful result is admissible to the benchmark."""

    return bool(
        result.get("success", False)
        and result.get("benchmark_eligible") is True
        and result.get("legacy_resume") is False
    )


def _batch_result_status(result: dict) -> str:
    """Return the explicit batch status for one normalized worker result."""

    if not result.get("success", False):
        return "failed"
    if _is_benchmark_eligible(result):
        return "success"
    return "success_with_warnings"


def _run_one_genome(
    *,
    spec: GenomeRunSpec,
    config: PipelineConfig,
    retries: int = 1,
    retry_delay_seconds: int = 60,
) -> dict:
    run_single = _single_genome_callable()
    last_exc: BaseException | None = None

    for attempt in range(retries + 1):
        try:
            if attempt:
                logger.info(
                    "%s: retry %d/%d after %.0fs",
                    spec.genome_id,
                    attempt,
                    retries,
                    retry_delay_seconds,
                )
            attempt_config = (
                config if attempt == 0 else config.with_overrides(resume=True)
            )
            return _normalize_worker_result(
                spec.genome_id,
                run_single(
                    genome_path=spec.input_path,
                    output_dir=spec.output_dir,
                    genome_id=spec.genome_id,
                    config=attempt_config,
                ),
            )
        except TSVInvariantError as exc:
            logger.error(
                "%s: deterministic invariant failure: %s",
                spec.genome_id,
                exc,
            )
            return _failure_result(spec.genome_id, exc)
        except Exception as exc:
            last_exc = exc
            logger.exception(
                "%s: run attempt %d/%d failed",
                spec.genome_id,
                attempt + 1,
                retries + 1,
            )
            if attempt < retries:
                time.sleep(retry_delay_seconds)

    assert last_exc is not None
    return _failure_result(spec.genome_id, last_exc)


def _write_batch_summary(output_base_dir: Path, results: list[dict]) -> Path:
    results = _enforce_exclusive_class_counts(results)
    summary_path = output_base_dir / "batch_summary.tsv"
    fields = [
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
        "vp",
        "plv",
        "mirus",
        "mixed",
        "total_bp",
        "genes",
        "hallmarks",
        "elapsed_sec",
        "error",
        "ppv",
        "unknown",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    for result in results:
        writer.writerow(
            {
                "genome_id": result.get("genome_id", "?"),
                "status": _batch_result_status(result),
                "benchmark_eligible": (
                    "true" if _is_benchmark_eligible(result) else "false"
                ),
                "legacy_resume": (
                    "true" if result.get("legacy_resume") is True else "false"
                ),
                "predictions": result.get("predictions", 0),
                "accepted": result.get("accepted", 0),
                "high_tier": result.get("high_tier", 0),
                "medium_tier": result.get("medium_tier", 0),
                "low_tier": result.get("low_tier", 0),
                "ncldv": result.get("ncldv_count", 0),
                "vp": result.get("vp_count", 0),
                "plv": result.get("plv_count", 0),
                "mirus": result.get("mirus_count", 0),
                "mixed": result.get("mixed_count", 0),
                "total_bp": result.get("accepted_bp", 0),
                "genes": result.get("total_genes", 0),
                "hallmarks": result.get("total_hallmarks", 0),
                "elapsed_sec": f"{result.get('elapsed_sec', 0):.0f}",
                "error": result.get("error", "") or "",
                "ppv": result.get("ppv_count", 0),
                "unknown": result.get("unknown_count", 0),
            }
        )
    atomic_write(summary_path, buffer.getvalue())
    return summary_path


def _write_batch_report(output_base_dir: Path, results: list[dict]) -> Path:
    results = _enforce_exclusive_class_counts(results)
    total_genomes = len(results)
    successful = sum(1 for r in results if r.get("success", False))
    benchmark_eligible = sum(
        1 for result in results if _batch_result_status(result) == "success"
    )
    success_with_warnings = sum(
        1
        for result in results
        if _batch_result_status(result) == "success_with_warnings"
    )
    legacy_resumes = sum(
        1
        for result in results
        if result.get("success", False) and result.get("legacy_resume") is True
    )
    total_accepted = sum(r.get("accepted", 0) for r in results)
    total_predictions = sum(r.get("predictions", 0) for r in results)
    total_high = sum(r.get("high_tier", 0) for r in results)
    total_medium = sum(r.get("medium_tier", 0) for r in results)
    total_low = sum(r.get("low_tier", 0) for r in results)
    total_time = sum(r.get("elapsed_sec", 0) for r in results)
    total_ncldv = sum(r.get("ncldv_count", 0) for r in results)
    total_vp = sum(r.get("vp_count", 0) for r in results)
    total_plv = sum(r.get("plv_count", 0) for r in results)
    total_mirus = sum(r.get("mirus_count", 0) for r in results)
    total_mixed = sum(r.get("mixed_count", 0) for r in results)
    total_ppv = sum(r.get("ppv_count", 0) for r in results)
    total_unknown = sum(r.get("unknown_count", 0) for r in results)
    total_classified = (
        total_ncldv
        + total_vp
        + total_plv
        + total_mirus
        + total_mixed
        + total_ppv
        + total_unknown
    )
    if total_classified != total_accepted:
        raise ValueError(
            "exclusive effective-class counts do not sum to accepted predictions: "
            f"accepted={total_accepted} classified={total_classified}"
        )
    total_bp = sum(r.get("accepted_bp", 0) for r in results)
    total_genes = sum(r.get("total_genes", 0) for r in results)
    total_hallmarks = sum(r.get("total_hallmarks", 0) for r in results)
    failed_results = [r for r in results if not r.get("success", False)]

    report_path = output_base_dir / "batch_report.md"
    with report_path.open("w") as handle:
        timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
        handle.write("# ViroSync Batch Report\n\n")
        handle.write(f"**Generated:** {timestamp}\n\n")
        handle.write("## Summary\n\n")
        handle.write(
            f"- **Genomes processed:** {total_genomes} "
            f"({successful} successful, {total_genomes - successful} failed)\n"
        )
        handle.write(
            f"- **Benchmark eligibility:** {benchmark_eligible} eligible, "
            f"{success_with_warnings} success with warnings\n"
        )
        handle.write(f"- **Legacy resumes:** {legacy_resumes}\n")
        handle.write(
            f"- **Total GEVEs:** {total_accepted} canonical "
            f"({total_predictions} candidates)\n"
        )
        handle.write(
            f"- **Canonical tiers:** HIGH={total_high}, MEDIUM={total_medium}, "
            f"LOW={total_low}\n"
        )
        handle.write(f"- **Total processing time:** {total_time:.0f}s\n\n")
        handle.write("### GEVE Classification\n\n")
        handle.write("| Category | Count | Description |\n")
        handle.write("|----------|-------|-------------|\n")
        handle.write(f"| NCLDV | {total_ncldv} | Nucleocytoviricota (giant viruses) |\n")
        handle.write(f"| VP | {total_vp} | Virophages |\n")
        handle.write(f"| PLV | {total_plv} | Polinton-like viruses |\n")
        handle.write(f"| MIRUS | {total_mirus} | Mirusviricota |\n")
        handle.write(f"| MIXED | {total_mixed} | Multiple viral lineages |\n")
        handle.write(f"| PPV | {total_ppv} | Preplasmiviricota |\n")
        handle.write(f"| UNKNOWN | {total_unknown} | Unrecognized effective class |\n")
        handle.write(f"| **Total** | **{total_classified}** | |\n\n")
        handle.write("### Region Statistics (Canonical GEVEs)\n\n")
        handle.write(f"- **Total bp:** {total_bp:,}\n")
        handle.write(f"- **Total genes:** {total_genes:,}\n")
        handle.write(f"- **Total hallmark markers:** {total_hallmarks:,}\n")
        if total_accepted > 0:
            handle.write(f"- **Average region size:** {total_bp // total_accepted:,} bp\n")
            handle.write(f"- **Average genes per region:** {total_genes / total_accepted:.1f}\n")
            handle.write(
                f"- **Average hallmarks per region:** "
                f"{total_hallmarks / total_accepted:.1f}\n"
            )
        handle.write("\n")
        handle.write("## Per-Genome Results\n\n")
        handle.write(
            "| Genome | Status | Benchmark eligible | Legacy resume | HIGH | MED | LOW | NCLDV | VP | PLV | MIRUS | Mixed | PPV | Unknown | bp | Genes | Time |\n"
        )
        handle.write(
            "|--------|--------|--------------------|---------------|------|-----|-----|-------|----|-----|-------|-------|-----|---------|----|-------|------|\n"
        )
        for result in sorted(results, key=lambda r: (-r.get("predictions", 0), r.get("genome_id", ""))):
            gid = result.get("genome_id", "?")
            if result.get("success", False):
                handle.write(
                    f"| {gid} | {_batch_result_status(result)} | "
                    f"{'yes' if _is_benchmark_eligible(result) else 'no'} | "
                    f"{'yes' if result.get('legacy_resume') is True else 'no'} | "
                    f"{result.get('high_tier', 0)} | "
                    f"{result.get('medium_tier', 0)} | {result.get('low_tier', 0)} | "
                    f"{result.get('ncldv_count', 0)} | {result.get('vp_count', 0)} | "
                    f"{result.get('plv_count', 0)} | {result.get('mirus_count', 0)} | "
                    f"{result.get('mixed_count', 0)} | {result.get('ppv_count', 0)} | "
                    f"{result.get('unknown_count', 0)} | {result.get('accepted_bp', 0):,} | "
                    f"{result.get('total_genes', 0)} | {result.get('elapsed_sec', 0):.0f}s |\n"
                )
            else:
                handle.write(
                    f"| {gid} | failed | no | no | - | - | - | - | - | - | - | - | - | - | - | - | FAILED |\n"
                )
        if failed_results:
            handle.write("\n## Failed Genomes\n\n")
            for result in failed_results:
                handle.write(
                    f"- **{result.get('genome_id', '?')}**: "
                    f"{result.get('error', 'Unknown error')}\n"
                )
    return report_path


def run_batch_python(
    *,
    genome_paths: Iterable[Path],
    output_base_dir: Path,
    config: PipelineConfig,
    max_concurrent_genomes: int,
    retries: int = 1,
    retry_delay_seconds: int = 60,
    effective_config: dict | None = None,
) -> list[dict]:
    """Run genomes concurrently with standard-library Python primitives."""
    genome_paths = [Path(path) for path in genome_paths]
    output_base_dir = Path(output_base_dir)
    run_specs = _preflight_genome_runs(genome_paths, output_base_dir)
    output_base_dir.mkdir(parents=True, exist_ok=True)
    if effective_config is not None:
        effective_path = output_base_dir / "effective_config.json"
        atomic_write(
            effective_path,
            json.dumps(effective_config, indent=2, sort_keys=True) + "\n",
        )

    if max_concurrent_genomes < 1:
        raise ValueError("max_concurrent_genomes must be >= 1")

    logger.info("=" * 60)
    logger.info("ViroSync Python batch processing: %d genomes", len(run_specs))
    logger.info("Max concurrent genomes: %d", max_concurrent_genomes)
    logger.info("=" * 60)

    results_by_index: dict[int, dict] = {}
    with ThreadPoolExecutor(max_workers=max_concurrent_genomes) as executor:
        futures = {
            executor.submit(
                _run_one_genome,
                spec=spec,
                config=config,
                retries=retries,
                retry_delay_seconds=retry_delay_seconds,
            ): (idx, spec.genome_id)
            for idx, spec in enumerate(run_specs)
        }
        completed = 0
        for future in as_completed(futures):
            idx, genome_id = futures[future]
            completed += 1
            try:
                result = _normalize_worker_result(genome_id, future.result())
            except Exception as exc:
                logger.exception(
                    "[%d/%d] %s: failed after runner exception",
                    completed,
                    len(run_specs),
                    genome_id,
                )
                result = _failure_result(genome_id, exc)
            results_by_index[idx] = result
            status = "ok" if result.get("success", False) else "failed"
            logger.info(
                "[%d/%d] %s: %s accepted=%s predictions=%s elapsed=%.0fs",
                completed,
                len(run_specs),
                genome_id,
                status,
                result.get("accepted", 0),
                result.get("predictions", 0),
                result.get("elapsed_sec", 0),
            )

    results = _enforce_exclusive_class_counts(
        [results_by_index[idx] for idx in range(len(run_specs))]
    )
    summary_path = _write_batch_summary(output_base_dir, results)
    report_path = _write_batch_report(output_base_dir, results)
    logger.info("Batch summary written: %s", summary_path)
    logger.info("Batch report written: %s", report_path)
    return results
