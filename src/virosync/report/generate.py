"""Generate the per-genome EVE analysis Jupyter notebook.

The notebook source is maintained as a single jupytext (percent-format) Python
file, ``eve_analysis.py`` -- the one source of truth. At runtime it is converted
to a notebook, parameterized for the genome, and executed with papermill so the
shipped ``.ipynb`` carries rendered figures and outputs.
"""

import logging
import tempfile
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

# Single source of truth: jupytext percent-format notebook source.
_SOURCE = Path(__file__).parent / "eve_analysis.py"
_KERNEL = "python3"
# The kernel starts in the output directory, where a relative PYTHONPATH finds nothing.
_PACKAGE_ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class EveReportPaths:
    """Path to the generated per-genome report notebook."""

    jupyter: Path


def generate_eve_report(
    output_dir: Path,
    genome_id: str,
    tax_labels_path: Path | None = None,
    host_credit_threshold: float = 0.35,
) -> EveReportPaths:
    """Render the per-genome EVE analysis notebook.

    The jupytext source is converted to a notebook, parameterized for this
    genome, and executed with papermill. When plotting is enabled, the executed
    notebook is required, so failures are raised rather than swallowed.
    """
    output_dir = Path(output_dir).resolve()
    if not _SOURCE.exists():
        raise FileNotFoundError(f"Notebook source not found: {_SOURCE}")

    try:
        import jupytext
        import papermill as pm
    except ImportError as exc:
        raise RuntimeError("jupytext and papermill are required to generate the EVE report") from exc

    jupyter_dir = output_dir / "notebooks" / "jupyter"
    jupyter_dir.mkdir(parents=True, exist_ok=True)
    output_path = jupyter_dir / "eve_analysis.ipynb"

    parameters = {
        "RESULTS_DIR": str(output_dir),
        "GENOME_ID": genome_id,
        "HOST_CREDIT_THRESHOLD": host_credit_threshold,
        "VIROSYNC_SRC": str(_PACKAGE_ROOT),
    }
    if tax_labels_path:
        parameters["TAX_LABELS_PATH"] = str(Path(tax_labels_path).resolve())

    # jupytext source (.py) -> notebook -> parameterized execution (papermill).
    notebook = jupytext.read(_SOURCE)
    with tempfile.TemporaryDirectory() as _tmp:
        source_ipynb = Path(_tmp) / "eve_analysis_source.ipynb"
        jupytext.write(notebook, source_ipynb)
        pm.execute_notebook(
            str(source_ipynb),
            str(output_path),
            parameters=parameters,
            cwd=str(output_dir),
            kernel_name=_KERNEL,
            progress_bar=logger.isEnabledFor(logging.DEBUG),
        )

    if not output_path.exists():
        raise RuntimeError(f"Jupyter report was not written: {output_path}")
    logger.info("EVE analysis notebook written: %s", output_path)

    return EveReportPaths(jupyter=output_path)
