"""
Lightweight per-task resource monitoring for Python orchestration runs.
"""

from __future__ import annotations

import json
import logging
import socket
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import psutil

from virosync.utils.path_safety import require_strict_child, safe_filename_component


@dataclass
class ResourceMetrics:
    task_name: str
    genome_id: str
    phase: str
    task_id: str
    wall_time_sec: float
    cpu_user_sec: float
    cpu_system_sec: float
    max_rss_bytes: int
    max_cpu_percent: float
    threads: int
    host: str
    pid: int
    timestamp: str


class ResourceMonitor:
    """
    Track wall time, CPU time, peak RSS, and CPU utilization for a task.
    """

    def __init__(
        self,
        task_name: str,
        genome_id: str,
        phase: str,
        output_dir: Path,
        threads: int,
        task_id: Optional[str] = None,
        sample_interval: float = 0.5,
    ) -> None:
        self.task_name = task_name
        self.genome_id = genome_id
        self.phase = phase
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.task_id = task_id
        self.sample_interval = sample_interval
        self._proc = psutil.Process()
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._max_rss = 0
        self._max_cpu_percent = 0.0
        self._start_time = 0.0
        self._start_cpu_user = 0.0
        self._start_cpu_system = 0.0

    def __enter__(self) -> "ResourceMonitor":
        logger = logging.getLogger(__name__)
        logger.info(
            "Resource monitor start: task=%s phase=%s genome=%s threads=%s",
            self.task_name,
            self.phase,
            self.genome_id,
            self.threads,
        )
        self._start_time = time.time()
        cpu_times = self._proc.cpu_times()
        self._start_cpu_user = cpu_times.user
        self._start_cpu_system = cpu_times.system
        self._prime_cpu_percent()
        self._thread = threading.Thread(target=self._sample_loop, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self._stop.set()
        if self._thread:
            self._thread.join(timeout=2)
        self._write_metrics()

    def _prime_cpu_percent(self) -> None:
        try:
            self._proc.cpu_percent(None)
            for child in self._proc.children(recursive=True):
                child.cpu_percent(None)
        except psutil.Error:
            return

    def _sample_loop(self) -> None:
        while not self._stop.is_set():
            try:
                rss = self._rss_total()
                self._max_rss = max(self._max_rss, rss)
                cpu_percent = self._cpu_percent_total()
                self._max_cpu_percent = max(self._max_cpu_percent, cpu_percent)
            except psutil.Error:
                pass
            time.sleep(self.sample_interval)

    def _rss_total(self) -> int:
        total = 0
        processes = [self._proc] + self._proc.children(recursive=True)
        for proc in processes:
            try:
                total += proc.memory_info().rss
            except psutil.Error:
                continue
        return total

    def _cpu_percent_total(self) -> float:
        total = 0.0
        processes = [self._proc] + self._proc.children(recursive=True)
        for proc in processes:
            try:
                total += proc.cpu_percent(None)
            except psutil.Error:
                continue
        return total

    def _write_metrics(self) -> None:
        end_time = time.time()
        cpu_times = self._proc.cpu_times()
        metrics = ResourceMetrics(
            task_name=self.task_name,
            genome_id=self.genome_id,
            phase=self.phase,
            task_id=self.task_id or self.genome_id,
            wall_time_sec=end_time - self._start_time,
            cpu_user_sec=cpu_times.user - self._start_cpu_user,
            cpu_system_sec=cpu_times.system - self._start_cpu_system,
            max_rss_bytes=self._max_rss,
            max_cpu_percent=self._max_cpu_percent,
            threads=self.threads,
            host=socket.gethostname(),
            pid=self._proc.pid,
            timestamp=time.strftime("%Y-%m-%dT%H:%M:%S"),
        )

        task_tag = safe_filename_component(self.task_id or self.genome_id)
        metrics_dir = self.output_dir / "resource_metrics"
        out_path = metrics_dir / f"{self.phase}_{self.task_name}_{task_tag}_{int(time.time())}.json"
        require_strict_child(metrics_dir, out_path)
        metrics_dir.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as f:
            json.dump(metrics.__dict__, f, indent=2)
        logger = logging.getLogger(__name__)
        logger.info(
            "Resource monitor end: task=%s phase=%s wall=%.1fs rss=%.2fGB cpu=%.1f%% out=%s",
            metrics.task_name,
            metrics.phase,
            metrics.wall_time_sec,
            metrics.max_rss_bytes / (1024 ** 3),
            metrics.max_cpu_percent,
            out_path,
        )
