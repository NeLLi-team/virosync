"""
Invariant checks for ``virosync_predictions_detailed.tsv`` outputs.

This module is used in two ways:
1) Automatically from the orchestration flow after each run.
2) Manually as a small regression-check script via ``python -m``.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

EMPTY_VALUES = {"", ".", "NA", "None", "null", "NULL"}
BOOL_TRUE_VALUES = {"1", "true", "True", "yes", "YES"}


@dataclass(frozen=True)
class InvariantIssue:
    eve_id: str
    check: str
    severity: str
    message: str


class TSVInvariantError(RuntimeError):
    """Raised after fatal detailed-TSV diagnostics have been written."""

    def __init__(self, report: "InvariantReport", report_path: Path):
        self.report = report
        self.report_path = Path(report_path)
        preview = ", ".join(
            f"{issue.eve_id}:{issue.check}" for issue in report.fatal_issues[:5]
        )
        detail = f" ({preview})" if preview else ""
        super().__init__(
            "Detailed TSV invariant check failed: "
            f"errors={report.error_count}, warnings={report.warning_count}; "
            f"report={self.report_path}{detail}"
        )


@dataclass
class InvariantReport:
    rows_checked: int
    issues: list[InvariantIssue]

    @property
    def issue_count(self) -> int:
        return len(self.issues)

    @property
    def warning_issues(self) -> list[InvariantIssue]:
        return [
            issue
            for issue in self.issues
            if issue.severity.strip().lower() == "warning"
        ]

    @property
    def fatal_issues(self) -> list[InvariantIssue]:
        # Fail closed: a misspelled or future severity must not silently turn a
        # scientific inconsistency into a successful run.
        return [
            issue
            for issue in self.issues
            if issue.severity.strip().lower() != "warning"
        ]

    @property
    def warning_count(self) -> int:
        return len(self.warning_issues)

    @property
    def error_count(self) -> int:
        return len(self.fatal_issues)

    @property
    def passed(self) -> bool:
        return self.error_count == 0

    @property
    def status(self) -> str:
        if not self.passed:
            return "FAIL"
        if self.warning_count:
            return "PASS_WITH_WARNINGS"
        return "PASS"


def _is_empty(value: object) -> bool:
    if value is None:
        return True
    return str(value).strip() in EMPTY_VALUES


def _parse_int(value: object) -> Optional[int]:
    if _is_empty(value):
        return 0
    try:
        return int(str(value))
    except (TypeError, ValueError):
        return None


def _parse_float(value: object) -> Optional[float]:
    if _is_empty(value):
        return 0.0
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return None


def _parse_taxonomy_counts(raw: str) -> tuple[dict[str, int], Optional[str]]:
    counts: dict[str, int] = {}
    if _is_empty(raw):
        return counts, None

    for token in str(raw).split(";"):
        token = token.strip()
        if not token:
            continue
        if ":" not in token:
            return {}, f"Invalid token '{token}' in taxonomy_best_hits"
        key, value = token.split(":", 1)
        key = key.strip()
        value = value.strip()
        try:
            counts[key] = int(value)
        except ValueError:
            return {}, f"Invalid count '{value}' for taxonomy key '{key}'"
    return counts, None


def _has_marker_token(marker_blob: str, pattern: re.Pattern[str]) -> bool:
    if _is_empty(marker_blob):
        return False
    return bool(pattern.search(marker_blob))


def _count_name_tokens(raw: object, sep: str = ",") -> int:
    if _is_empty(raw):
        return 0
    return len([token.strip() for token in str(raw).split(sep) if token.strip()])


def _parse_counted_name_blob(raw: object, sep: str = ",") -> tuple[dict[str, int], Optional[str]]:
    if _is_empty(raw):
        return {}, None

    counts: dict[str, int] = {}
    for token in str(raw).split(sep):
        token = token.strip()
        if not token:
            continue
        if ":" not in token:
            return {}, f"Invalid token '{token}' in counted name list"
        name, count_raw = token.rsplit(":", 1)
        name = name.strip()
        if not name:
            return {}, f"Missing model name in token '{token}'"
        try:
            count = int(count_raw)
        except ValueError:
            return {}, f"Invalid count '{count_raw}' for model '{name}'"
        counts[name] = count
    return counts, None


def _load_gene_taxonomy_totals(
    gene_taxonomy_all_tsv: Optional[Path],
) -> dict[str, dict[str, int]]:
    """
    Load per-EVE interior-gene taxonomy summaries from gene_taxonomy_all.tsv.
    """
    if gene_taxonomy_all_tsv is None or not gene_taxonomy_all_tsv.exists():
        return {}

    by_eve: dict[str, dict[str, int]] = {}
    with gene_taxonomy_all_tsv.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            eve_id = (row.get("eve_id") or "").strip()
            if not eve_id:
                continue

            is_flanking = str(row.get("is_flanking", "")).strip() in BOOL_TRUE_VALUES
            if is_flanking:
                continue

            state = by_eve.setdefault(
                eve_id,
                {
                    "total_interior": 0,
                    "ncldv_top10": 0,
                    "mirus_top10": 0,
                },
            )

            state["total_interior"] += 1
            top10_raw = str(row.get("top10_origins", "")).strip()
            if top10_raw and top10_raw not in EMPTY_VALUES:
                top10_tokens = {token.strip() for token in top10_raw.split(",") if token.strip()}
                if "NCLDV" in top10_tokens:
                    state["ncldv_top10"] += 1
                if "MIRUS" in top10_tokens:
                    state["mirus_top10"] += 1
            else:
                # Backward-compatible fallback for older exports without top10_origins.
                if str(row.get("has_ncldv_top10", "")).strip() in BOOL_TRUE_VALUES:
                    state["ncldv_top10"] += 1
                if str(row.get("has_mirus_top10", "")).strip() in BOOL_TRUE_VALUES:
                    state["mirus_top10"] += 1

    return by_eve


def run_tsv_invariant_checks(
    detailed_tsv: Path,
    gene_taxonomy_all_tsv: Optional[Path] = None,
) -> InvariantReport:
    """
    Check core TSV invariants that previously regressed.
    """
    issues: list[InvariantIssue] = []
    rows_checked = 0

    og_pattern = re.compile(r"\bOG\d+\b", re.IGNORECASE)
    gvogm_pattern = re.compile(r"\bGVOGM\d+\b", re.IGNORECASE)

    gene_tax_summaries = _load_gene_taxonomy_totals(gene_taxonomy_all_tsv)

    with detailed_tsv.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or [])
        required_support_fields = {
            "total_proteins",
            "ncldv_top10_proteins",
            "mirus_top10_proteins",
            "ppv_top10_proteins",
            "cress_top10_proteins",
        }
        missing_support_fields = sorted(required_support_fields - fieldnames)
        if "eve_id" not in fieldnames:
            issues.append(
                InvariantIssue(
                    eve_id=".",
                    check="detailed_tsv_schema",
                    severity="error",
                    message="Detailed TSV header is missing required column: eve_id",
                )
            )
        for row in reader:
            rows_checked += 1
            eve_id = (row.get("eve_id") or f"row_{rows_checked}").strip()

            def add_issue(check: str, message: str, severity: str = "error") -> None:
                issues.append(
                    InvariantIssue(
                        eve_id=eve_id,
                        check=check,
                        severity=severity,
                        message=message,
                    )
                )

            total_proteins = _parse_int(row.get("total_proteins"))
            og_count = _parse_int(row.get("og_count"))
            gvogm_count = _parse_int(row.get("gvogm_count"))
            og_unvalidated_count = _parse_int(row.get("og_unvalidated_count"))
            gvogm_unvalidated_count = _parse_int(row.get("gvogm_unvalidated_count"))
            ncldv_top_count = _parse_int(row.get("ncldv_top10_proteins"))
            mirus_top_count = _parse_int(row.get("mirus_top10_proteins"))
            ppv_top_count = _parse_int(row.get("ppv_top10_proteins"))
            cress_top_count = _parse_int(row.get("cress_top10_proteins"))
            host_sig_count = _parse_int(row.get("host_signature_gene_count"))
            host_sig_fraction = _parse_float(row.get("host_signature_fraction"))

            for field_name, value in (
                ("total_proteins", total_proteins),
                ("og_count", og_count),
                ("gvogm_count", gvogm_count),
                ("og_unvalidated_count", og_unvalidated_count),
                ("gvogm_unvalidated_count", gvogm_unvalidated_count),
                ("ncldv_top10_proteins", ncldv_top_count),
                ("mirus_top10_proteins", mirus_top_count),
                ("ppv_top10_proteins", ppv_top_count),
                ("cress_top10_proteins", cress_top_count),
                ("host_signature_gene_count", host_sig_count),
                ("host_signature_fraction", host_sig_fraction),
            ):
                if value is None:
                    add_issue("parse_error", f"Field '{field_name}' is not numeric")

            if total_proteins is None:
                total_proteins = 0
            if og_count is None:
                og_count = 0
            if gvogm_count is None:
                gvogm_count = 0
            if og_unvalidated_count is None:
                og_unvalidated_count = 0
            if gvogm_unvalidated_count is None:
                gvogm_unvalidated_count = 0
            if ncldv_top_count is None:
                ncldv_top_count = 0
            if mirus_top_count is None:
                mirus_top_count = 0
            if ppv_top_count is None:
                ppv_top_count = 0
            if cress_top_count is None:
                cress_top_count = 0
            if host_sig_count is None:
                host_sig_count = 0
            if host_sig_fraction is None:
                host_sig_fraction = 0.0

            taxonomy_counts, taxonomy_error = _parse_taxonomy_counts(
                str(row.get("taxonomy_best_hits", "")).strip()
            )
            if taxonomy_error:
                add_issue("taxonomy_best_hits_format", taxonomy_error)
            elif taxonomy_counts:
                observed_total = sum(taxonomy_counts.values())
                if total_proteins != observed_total:
                    add_issue(
                        "taxonomy_total_mismatch",
                        f"total_proteins={total_proteins}, taxonomy_sum={observed_total}",
                    )

            if ncldv_top_count > total_proteins:
                add_issue(
                    "ncldv_top_count_out_of_range",
                    f"ncldv_top10_proteins={ncldv_top_count} > total_proteins={total_proteins}",
                )
            if mirus_top_count > total_proteins:
                add_issue(
                    "mirus_top_count_out_of_range",
                    f"mirus_top10_proteins={mirus_top_count} > total_proteins={total_proteins}",
                )
            if ppv_top_count > total_proteins:
                add_issue(
                    "ppv_top_count_out_of_range",
                    f"ppv_top10_proteins={ppv_top_count} > total_proteins={total_proteins}",
                )
            if cress_top_count > total_proteins:
                add_issue(
                    "cress_top_count_out_of_range",
                    f"cress_top10_proteins={cress_top_count} > total_proteins={total_proteins}",
                )

            if host_sig_count > total_proteins:
                add_issue(
                    "host_signature_count_out_of_range",
                    f"host_signature_gene_count={host_sig_count} > total_proteins={total_proteins}",
                )
            if host_sig_count > 0 and host_sig_fraction <= 0.0:
                add_issue(
                    "host_signature_fraction_zero",
                    "host_signature_gene_count > 0 but host_signature_fraction <= 0",
                )
            if host_sig_count == 0 and host_sig_fraction > 0.0:
                add_issue(
                    "host_signature_fraction_nonzero",
                    "host_signature_gene_count == 0 but host_signature_fraction > 0",
                )
            if host_sig_fraction > 1.0:
                add_issue(
                    "host_signature_fraction_out_of_range",
                    f"host_signature_fraction={host_sig_fraction} > 1.0",
                )

            og_names = row.get("og_names", "")
            gvogm_names = row.get("gvogm_names", "")
            og_unvalidated_names = row.get("og_unvalidated_names", "")
            gvogm_unvalidated_names = row.get("gvogm_unvalidated_names", "")
            og_unvalidated_name_count = _count_name_tokens(og_unvalidated_names)
            gvogm_unvalidated_name_count = _count_name_tokens(gvogm_unvalidated_names)
            og_name_counts, og_name_error = _parse_counted_name_blob(og_names)
            gvogm_name_counts, gvogm_name_error = _parse_counted_name_blob(gvogm_names)

            if og_name_error:
                add_issue(
                    "og_names_format",
                    og_name_error,
                )
            if gvogm_name_error:
                add_issue(
                    "gvogm_names_format",
                    gvogm_name_error,
                )
            if og_count > 0 and _is_empty(og_names):
                add_issue("og_names_missing", "og_count > 0 but og_names is empty")
            if gvogm_count > 0 and _is_empty(gvogm_names):
                add_issue("gvogm_names_missing", "gvogm_count > 0 but gvogm_names is empty")
            if og_count == 0 and not _is_empty(og_names):
                add_issue("og_names_nonzero", "og_count == 0 but og_names is populated")
            if gvogm_count == 0 and not _is_empty(gvogm_names):
                add_issue(
                    "gvogm_names_nonzero",
                    "gvogm_count == 0 but gvogm_names is populated",
                )
            if og_count > sum(og_name_counts.values()):
                add_issue(
                    "og_count_names_mismatch",
                    f"og_count={og_count}, og_names_total={sum(og_name_counts.values())}",
                )
            if gvogm_count > sum(gvogm_name_counts.values()):
                add_issue(
                    "gvogm_count_names_mismatch",
                    f"gvogm_count={gvogm_count}, gvogm_names_total={sum(gvogm_name_counts.values())}",
                )
            if og_unvalidated_count > 0 and _is_empty(og_unvalidated_names):
                add_issue(
                    "og_unvalidated_names_missing",
                    "og_unvalidated_count > 0 but og_unvalidated_names is empty",
                )
            if gvogm_unvalidated_count > 0 and _is_empty(gvogm_unvalidated_names):
                add_issue(
                    "gvogm_unvalidated_names_missing",
                    "gvogm_unvalidated_count > 0 but gvogm_unvalidated_names is empty",
                )
            if og_unvalidated_count == 0 and not _is_empty(og_unvalidated_names):
                add_issue(
                    "og_unvalidated_names_nonzero",
                    "og_unvalidated_count == 0 but og_unvalidated_names is populated",
                )
            if gvogm_unvalidated_count == 0 and not _is_empty(gvogm_unvalidated_names):
                add_issue(
                    "gvogm_unvalidated_names_nonzero",
                    "gvogm_unvalidated_count == 0 but gvogm_unvalidated_names is populated",
                )
            if og_unvalidated_count != og_unvalidated_name_count:
                add_issue(
                    "og_unvalidated_count_names_mismatch",
                    f"og_unvalidated_count={og_unvalidated_count}, og_unvalidated_names_count={og_unvalidated_name_count}",
                )
            if gvogm_unvalidated_count != gvogm_unvalidated_name_count:
                add_issue(
                    "gvogm_unvalidated_count_names_mismatch",
                    f"gvogm_unvalidated_count={gvogm_unvalidated_count}, gvogm_unvalidated_names_count={gvogm_unvalidated_name_count}",
                )

            marker_blob = "|".join(
                [
                    str(row.get("seed_marker_names", "")),
                    str(row.get("other_marker_names", "")),
                    str(row.get("seed_marker_patterns", "")),
                    str(row.get("other_marker_patterns", "")),
                ]
            )
            if (
                _has_marker_token(marker_blob, og_pattern)
                and (og_count + og_unvalidated_count) == 0
            ):
                add_issue(
                    "og_marker_evidence_zero_count",
                    "Marker columns contain OG identifiers but both og_count and og_unvalidated_count are 0",
                )
            if (
                _has_marker_token(marker_blob, gvogm_pattern)
                and (gvogm_count + gvogm_unvalidated_count) == 0
            ):
                add_issue(
                    "gvogm_marker_evidence_zero_count",
                    "Marker columns contain GVOGm identifiers but both gvogm_count and gvogm_unvalidated_count are 0",
                )

            seed_marker_names = row.get("seed_marker_names", "")
            seed_marker_patterns = row.get("seed_marker_patterns", "")
            if not _is_empty(seed_marker_patterns) and _is_empty(seed_marker_names):
                add_issue(
                    "seed_marker_names_missing",
                    "seed_marker_patterns is populated but seed_marker_names is empty",
                )
            if not _is_empty(seed_marker_names) and _is_empty(seed_marker_patterns):
                add_issue(
                    "seed_marker_patterns_missing",
                    "seed_marker_names is populated but seed_marker_patterns is empty",
                )

            gene_tax = gene_tax_summaries.get(eve_id)
            if gene_tax:
                interior_total = gene_tax["total_interior"]
                interior_ncldv = gene_tax["ncldv_top10"]
                interior_mirus = gene_tax["mirus_top10"]

                if total_proteins != interior_total:
                    add_issue(
                        "gene_taxonomy_total_mismatch",
                        f"total_proteins={total_proteins}, gene_taxonomy_interior={interior_total}",
                    )
                if ncldv_top_count != interior_ncldv:
                    add_issue(
                        "gene_taxonomy_ncldv_mismatch",
                        f"ncldv_top10_proteins={ncldv_top_count}, interior_ncldv_top10={interior_ncldv}",
                    )
                if mirus_top_count != interior_mirus:
                    add_issue(
                        "gene_taxonomy_mirus_mismatch",
                        f"mirus_top10_proteins={mirus_top_count}, interior_mirus_top10={interior_mirus}",
                    )

        if rows_checked and missing_support_fields:
            issues.append(
                InvariantIssue(
                    eve_id=".",
                    check="detailed_tsv_schema",
                    severity="error",
                    message=(
                        "Non-empty detailed TSV header is missing required columns: "
                        + ", ".join(missing_support_fields)
                    ),
                )
            )

    return InvariantReport(rows_checked=rows_checked, issues=issues)


def write_tsv_invariant_report(report: InvariantReport, output_path: Path) -> Path:
    """
    Write invariant-check results as a TSV report.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "status",
                "rows_checked",
                "issue_count",
                "error_count",
                "warning_count",
            ]
        )
        writer.writerow(
            [
                report.status,
                report.rows_checked,
                report.issue_count,
                report.error_count,
                report.warning_count,
            ]
        )
        writer.writerow([])
        writer.writerow(["eve_id", "check", "severity", "message"])
        for issue in report.issues:
            writer.writerow([issue.eve_id, issue.check, issue.severity, issue.message])
    return output_path


def enforce_tsv_invariants(
    detailed_tsv: Path,
    report_out: Path,
    gene_taxonomy_all_tsv: Optional[Path] = None,
) -> InvariantReport:
    """Write diagnostics and raise only for fatal detailed-TSV issues.

    Missing, unreadable, and malformed inputs become reportable scientific
    failures rather than bypassing the completion gate. Warning-only reports
    remain successful and are labeled distinctly.
    """
    detailed_tsv = Path(detailed_tsv)
    report_out = Path(report_out)

    if not detailed_tsv.is_file():
        report = InvariantReport(
            rows_checked=0,
            issues=[
                InvariantIssue(
                    eve_id=".",
                    check="detailed_tsv_missing",
                    severity="error",
                    message=f"Detailed TSV not found: {detailed_tsv}",
                )
            ],
        )
    else:
        try:
            report = run_tsv_invariant_checks(
                detailed_tsv=detailed_tsv,
                gene_taxonomy_all_tsv=gene_taxonomy_all_tsv,
            )
        except (OSError, UnicodeError, csv.Error) as exc:
            report = InvariantReport(
                rows_checked=0,
                issues=[
                    InvariantIssue(
                        eve_id=".",
                        check="detailed_tsv_unreadable",
                        severity="error",
                        message=f"Could not read detailed TSV {detailed_tsv}: {exc}",
                    )
                ],
            )

    write_tsv_invariant_report(report, report_out)
    if not report.passed:
        raise TSVInvariantError(report, report_out)
    return report


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check invariants in virosync detailed TSV outputs."
    )
    parser.add_argument(
        "--detailed-tsv",
        required=True,
        type=Path,
        help="Path to virosync_predictions_detailed.tsv",
    )
    parser.add_argument(
        "--gene-taxonomy-all-tsv",
        default=None,
        type=Path,
        help="Optional path to gene_taxonomy/gene_taxonomy_all.tsv",
    )
    parser.add_argument(
        "--report-out",
        default=None,
        type=Path,
        help="Output report TSV path (default: sibling virosync_tsv_invariant_report.tsv)",
    )
    parser.add_argument(
        "--fail-on-issues",
        action="store_true",
        help="Exit non-zero when error-severity invariant issues are found",
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = _parse_args(argv)
    detailed_tsv = args.detailed_tsv
    if not detailed_tsv.exists():
        raise SystemExit(f"Detailed TSV not found: {detailed_tsv}")

    report_out = args.report_out
    if report_out is None:
        report_out = detailed_tsv.parent / "virosync_tsv_invariant_report.tsv"

    gene_taxonomy_all = args.gene_taxonomy_all_tsv
    if gene_taxonomy_all is not None and not gene_taxonomy_all.exists():
        gene_taxonomy_all = None

    report = run_tsv_invariant_checks(
        detailed_tsv=detailed_tsv,
        gene_taxonomy_all_tsv=gene_taxonomy_all,
    )
    write_tsv_invariant_report(report, report_out)

    print(
        f"TSV invariant check: {report.status} "
        f"(rows={report.rows_checked}, issues={report.issue_count}, "
        f"errors={report.error_count}, warnings={report.warning_count})"
    )
    print(f"Report: {report_out}")

    if args.fail_on_issues and report.error_count > 0:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
