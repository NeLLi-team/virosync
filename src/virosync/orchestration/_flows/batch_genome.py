"""Batch genome processing compatibility wrapper."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from virosync.config import PipelineConfig
from virosync.orchestration.python_runner import run_batch_python


def batch_genome_flow(
    genome_paths: list[Path],
    output_base_dir: Path,
    config: Optional[PipelineConfig] = None,
    threads_per_genome: Optional[int] = None,
    max_concurrent_genomes: Optional[int] = None,
    **kwargs,
) -> list[dict]:
    """Process multiple genomes with the plain-Python batch runner.

    The historical public name is retained for compatibility with callers that
    imported ``batch_genome_flow`` directly.
    """
    pipeline_config = config or PipelineConfig()
    overrides = {key: value for key, value in kwargs.items() if value is not None}
    if threads_per_genome is not None:
        overrides["threads"] = threads_per_genome
    if overrides:
        pipeline_config = pipeline_config.with_overrides(**overrides)

    concurrency = max_concurrent_genomes or 4
    if concurrency < 1:
        raise ValueError("max_concurrent_genomes must be >= 1")

    return run_batch_python(
        genome_paths=genome_paths,
        output_base_dir=output_base_dir,
        config=pipeline_config,
        max_concurrent_genomes=concurrency,
    )
