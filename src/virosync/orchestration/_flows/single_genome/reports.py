"""Required per-genome report generation (Jupyter notebook)."""

from pathlib import Path
from typing import Optional


def _generate_required_reports(
    output_dir: Path,
    genome_id: str,
    taxonomy_labels_file: Optional[Path],
    logger,
) -> dict[str, str]:
    """Generate the required per-genome Jupyter EVE analysis notebook."""
    from virosync.report.generate import generate_eve_report

    report_paths = generate_eve_report(
        output_dir=output_dir,
        genome_id=genome_id,
        tax_labels_path=(
            Path(taxonomy_labels_file) if taxonomy_labels_file else None
        ),
    )
    logger.info("EVE analysis notebook: %s", report_paths.jupyter)
    return {
        "eve_analysis_notebook": str(report_paths.jupyter),
    }
