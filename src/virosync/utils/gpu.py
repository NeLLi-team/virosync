"""GPU memory management and device selection utilities.

Provides functions to select the best available GPU and release GPU memory
between pipeline phases to prevent CUDA OOM errors when worker threads process
multiple genomes sequentially.
"""

from __future__ import annotations

import gc
import logging

logger = logging.getLogger(__name__)




def release_gpu_memory() -> None:
    """Release GPU memory held by cached models and tensors.

    Clears the global TMVec predictor (the largest GPU consumer), runs Python
    garbage collection, and empties the CUDA memory cache. Call this between
    pipeline phases to prevent OOM when a worker leaves Phase 3.
    """
    try:
        from virosync.pipeline.phase3.tmvec_predictor import release_tmvec_predictor
        release_tmvec_predictor()
    except Exception as exc:
        logger.debug("TMVec predictor release skipped: %s", exc)

    gc.collect()

    try:
        import torch
        if torch.cuda.is_available():
            before = torch.cuda.memory_allocated()
            torch.cuda.empty_cache()
            after = torch.cuda.memory_allocated()
            freed_mib = (before - after) / (1024 ** 2)
            if freed_mib > 1:
                logger.info("Released %.0f MiB GPU memory", freed_mib)
    except ImportError:
        pass

