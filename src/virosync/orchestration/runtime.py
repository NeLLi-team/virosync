"""Runtime helpers for plain-Python orchestration."""

from __future__ import annotations

import logging
from typing import Any, Callable


def get_orchestration_logger(name: str = "virosync.orchestration") -> logging.Logger:
    """Return a standard Python logger for orchestration code."""
    return logging.getLogger(name)


def call_task(task_or_func: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
    """Call a task-like object or plain function without orchestration side effects."""
    raw_func = getattr(task_or_func, "fn", task_or_func)
    return raw_func(*args, **kwargs)
