"""Single-genome processing package.

The implementation is split across submodules:

- ``orchestrator`` -- the ``single_genome_flow`` entry point and implementation
- ``phase1`` / ``phase2`` / ``phase3`` -- per-phase subflows
- ``loaders`` / ``manifest`` / ``resume`` / ``reports`` -- helper functions

Import from ``virosync.orchestration.flows`` for the public interface. This
package re-exports the historical module-level names so the import path
``virosync.orchestration._flows.single_genome`` stays stable for callers/tests.
"""

from .orchestrator import (
    _pin_cuda_device,
    _run_phase0_subflow,
    _single_genome_flow_impl,
    run_single_genome_task,
    single_genome_flow,
)
from .phase1 import _run_phase1_subflow
from .phase2 import _run_phase2_subflow
from .phase3 import _run_phase3_subflow
from .loaders import (
    _build_merged_seeds_from_regions,
    _count_fasta_records,
    _load_interproscan_summary,
    _load_tmvec_cache,
    _safe_int,
    _serialize_tmvec_cache,
)
from .manifest import (
    _compute_config_fingerprint,
    _json_safe,
    _summarize_predictions_tsv,
    _write_completion_manifest,
    _write_empty_run_log,
)
from .reports import _generate_required_reports
from .resume import (
    _completed_run_artifacts,
    _first_valid_tsv,
    _require_phase2b_gene_taxonomy_db,
    _valid_completion_manifest,
    _valid_resume_run_log,
    _valid_tsv_header,
)

__all__ = [
    "single_genome_flow",
    "run_single_genome_task",
    "_single_genome_flow_impl",
    "_run_phase0_subflow",
    "_run_phase1_subflow",
    "_run_phase2_subflow",
    "_run_phase3_subflow",
    "_completed_run_artifacts",
    "_require_phase2b_gene_taxonomy_db",
    "_summarize_predictions_tsv",
    "_write_completion_manifest",
]
