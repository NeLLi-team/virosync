"""Optional per-genome report generation (Jupyter notebook)."""

from pathlib import Path


def _generate_required_reports(
    output_dir: Path,
    genome_id: str,
    taxonomy_labels_file: Path | None,
    logger,
    host_credit_threshold: float,
    generate_plots: bool,
) -> dict[str, str]:
    """Render the per-genome notebook only when plotting is enabled."""
    if not generate_plots:
        logger.info("Plotting disabled; scientific exports and audits remain enabled")
        return {}

    from virosync.report.generate import generate_eve_report

    report_paths = generate_eve_report(
        output_dir=output_dir,
        genome_id=genome_id,
        host_credit_threshold=host_credit_threshold,
        tax_labels_path=(Path(taxonomy_labels_file) if taxonomy_labels_file else None),
    )
    logger.info("EVE analysis notebook: %s", report_paths.jupyter)
    return {
        "eve_analysis_notebook": str(report_paths.jupyter),
    }
