"""ViroSync orchestration module.

This module provides plain-Python orchestration for the ViroSync pipeline.

Key components:
- tasks: Granular Python task functions for each pipeline phase
- flows: Pipeline functions for single and batch genome processing
- cli: Command-line interface for orchestrated runs
- utils: Helper functions for data wiring
"""

from virosync.utils.ssl_env import clear_stale_ssl_env_vars

# Clear stale SSL overrides first so relocated environments still work.
clear_stale_ssl_env_vars()

from virosync.orchestration.flows import (
    single_genome_flow,
)
from virosync.orchestration.tasks import (
    generate_outputs_task,
    generate_proteome_task,
    hhg_seeding_task,
    mask_genome_task,
    verify_eve_candidates_batched_task,
    verify_eve_task,
)

__all__ = [
    # Tasks
    "mask_genome_task",
    "generate_proteome_task",
    "hhg_seeding_task",
    "verify_eve_task",
    "verify_eve_candidates_batched_task",
    "generate_outputs_task",
    # Flows
    "single_genome_flow",
]
