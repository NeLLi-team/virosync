"""
ViroSync orchestration entry points.

This module is the canonical import location for orchestration functions.
Internal implementations are in _flows/ package.

Example:
    from virosync.orchestration.flows import single_genome_flow, batch_genome_flow

    # Process a single genome
    result = single_genome_flow(
        genome_path=Path("genome.fasta"),
        output_dir=Path("output/my_genome"),
        genome_id="my_genome",
        config=PipelineConfig(),
    )

    # Process multiple genomes in parallel
    results = batch_genome_flow(
        genome_paths=[Path("g1.fasta"), Path("g2.fasta")],
        output_base_dir=Path("results/"),
        config=PipelineConfig(),
    )
"""

# Re-export from internal package
from virosync.orchestration._flows.single_genome import single_genome_flow
from virosync.orchestration._flows.batch_genome import batch_genome_flow
from virosync.orchestration._flows.utils import (
    _merge_config_with_kwargs,
    _matches_allowlist,
    build_marker_faa,
    ensure_combined_faa,
)

__all__ = [
    # Main flows
    "single_genome_flow",
    "batch_genome_flow",
    # Utilities (for backwards compatibility)
    "_merge_config_with_kwargs",
    "_matches_allowlist",
    "build_marker_faa",
    "ensure_combined_faa",
]
