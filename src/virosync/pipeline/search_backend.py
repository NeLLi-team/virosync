"""
Unified sequence search backend for ViroSync.

Provides a single entry point for running protein sequence searches
against Diamond databases. All Diamond call sites delegate to
run_sequence_search(); Diamond is the sole supported backend (since v1.0.5).
"""

import logging
import os
import signal
import subprocess
import tempfile
import time
from pathlib import Path

logger = logging.getLogger(__name__)

# Default output columns matching Diamond outfmt 6
DEFAULT_COLUMNS = ["qseqid", "sseqid", "evalue", "bitscore", "pident", "qcovhsp"]

DEFAULT_SEARCH_TIMEOUT = 3600


def _env_int(name: str, default: int) -> int:
    value = os.environ.get(name)
    if value is None or value == "":
        return default
    try:
        parsed = int(value)
    except ValueError:
        logger.warning("Invalid integer for %s=%r; using %d", name, value, default)
        return default
    return parsed


def run_sequence_search(
    query_fasta: Path,
    db_path: Path,
    output_tsv: Path,
    threads: int,
    backend: str = "diamond",
    evalue: float = 1e-5,
    max_target_seqs: int = 10,
    output_columns: list[str] | None = None,
    extra_flags: list[str] | None = None,
    timeout: int = DEFAULT_SEARCH_TIMEOUT,
) -> Path:
    """
    Run a protein sequence search using the specified backend.

    Args:
        query_fasta: Path to query FASTA file.
        db_path: Path to the Diamond database file (.dmnd).
        output_tsv: Path for the output TSV file.
        threads: Number of CPU threads.
        backend: must be "diamond" (the only supported backend).
        evalue: E-value threshold.
        max_target_seqs: Maximum hits per query.
        output_columns: Column names in Diamond outfmt 6 convention.
                        Default: qseqid sseqid evalue bitscore pident qcovhsp.
        extra_flags: Additional CLI flags passed verbatim to the tool.
        timeout: Timeout in seconds.

    Returns:
        Path to the output TSV file.
    """
    if output_columns is None:
        output_columns = list(DEFAULT_COLUMNS)

    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    if backend != "diamond":
        raise ValueError(
            f"Unsupported search backend: {backend!r} (only 'diamond' is supported)"
        )
    return _run_diamond(
        query_fasta=query_fasta,
        db_path=db_path,
        output_tsv=output_tsv,
        threads=threads,
        evalue=evalue,
        max_target_seqs=max_target_seqs,
        output_columns=output_columns,
        extra_flags=extra_flags,
        timeout=timeout,
    )


# Diamond 2.1.21 occasionally hits a finalization-phase futex deadlock when
# invoked many times in sequence from a long-running parent (Apr/May 2026).
# The same input can stop making CPU/IO progress in the orchestrator while
# finishing quickly in a fresh process. The defenses below contain the blast
# radius:
#   * per-call ``tempfile.TemporaryDirectory`` passed via Diamond's ``--tmpdir``
#     so each invocation gets a clean workspace and inherits no filesystem
#     state from prior runs
#   * process-group isolation so a killed parent cannot leave Diamond workers
#     running in the background
#   * wall-clock and no-progress watchdogs, followed by one reduced-thread
#     retry, which clears the futex contention seen in benchmark reruns
# 0 = auto: derive a reduced-thread retry from the original thread count
# (max(8, threads//2), never exceeding the original). Set a positive value to
# force a specific retry thread count.
DIAMOND_RETRY_THREADS = _env_int("VIROSYNC_DIAMOND_RETRY_THREADS", 0)
DIAMOND_DEFAULT_TIMEOUT = _env_int("VIROSYNC_DIAMOND_TIMEOUT_SEC", 900)
DIAMOND_NO_PROGRESS_TIMEOUT = _env_int(
    "VIROSYNC_DIAMOND_NO_PROGRESS_TIMEOUT_SEC", 60
)
DIAMOND_NO_PROGRESS_MIN_RUNTIME = _env_int(
    "VIROSYNC_DIAMOND_NO_PROGRESS_MIN_RUNTIME_SEC", 120
)
DIAMOND_WATCHDOG_POLL_INTERVAL = _env_int(
    "VIROSYNC_DIAMOND_WATCHDOG_POLL_INTERVAL_SEC", 10
)


class DiamondNoProgressTimeout(subprocess.TimeoutExpired):
    """Raised when Diamond remains alive but CPU and IO counters stop moving."""


def _read_process_activity(pid: int) -> tuple[int, int, int, int, int, int, int] | None:
    """Return cumulative CPU and IO counters for watchdog progress checks."""
    try:
        stat_text = Path(f"/proc/{pid}/stat").read_text()
        after_comm = stat_text.rsplit(") ", 1)[1].split()
        cpu_ticks = int(after_comm[11]) + int(after_comm[12])

        io_values: dict[str, int] = {}
        for line in Path(f"/proc/{pid}/io").read_text().splitlines():
            key, value = line.split(":", 1)
            io_values[key] = int(value.strip())
    except (OSError, IndexError, ValueError):
        return None

    return (
        cpu_ticks,
        io_values.get("rchar", 0),
        io_values.get("wchar", 0),
        io_values.get("syscr", 0),
        io_values.get("syscw", 0),
        io_values.get("read_bytes", 0),
        io_values.get("write_bytes", 0),
    )


def _terminate_process_group(proc: subprocess.Popen) -> None:
    """Terminate a start_new_session=True process and its children."""
    if proc.poll() is not None:
        return
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        proc.wait(timeout=10)
        return
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        return
    proc.wait(timeout=10)


def _run_diamond(
    query_fasta: Path,
    db_path: Path,
    output_tsv: Path,
    threads: int,
    evalue: float,
    max_target_seqs: int,
    output_columns: list[str],
    extra_flags: list[str] | None,
    timeout: int,
) -> Path:
    """Run Diamond blastp search.

    Wraps Diamond with isolation (per-call tempdir, new session, scrubbed
    env) and an automatic single retry at reduced thread count if the call
    times out. Caller-facing exceptions are unchanged: ``RuntimeError`` on
    final timeout, ``CalledProcessError`` on non-zero exit.
    """
    db_prefix = db_path
    if db_prefix.suffix == ".dmnd":
        db_prefix = db_prefix.with_suffix("")

    base_cmd = [
        "diamond", "blastp",
        "--query", str(query_fasta),
        "--db", str(db_prefix),
        "--out", str(output_tsv),
        "--evalue", str(evalue),
        "--max-target-seqs", str(max_target_seqs),
        "--outfmt", "6", *output_columns,
    ]
    if extra_flags:
        base_cmd.extend(extra_flags)

    # Clear environment variables that multiprocessing workers may set which
    # restrict Diamond's internal thread pool and cause deadlocks.
    env = os.environ.copy()
    env.pop("OMP_NUM_THREADS", None)
    env.pop("OPENBLAS_NUM_THREADS", None)
    env.pop("MKL_NUM_THREADS", None)

    def _invoke(t: int) -> None:
        # Isolated temp directory per call. DIAMOND writes intermediate index chunks
        # here and we want them gone (and uncontended) the moment the call
        # completes, regardless of where the caller's TMPDIR points.
        with tempfile.TemporaryDirectory(prefix="diamond_call_") as tdir:
            cmd = base_cmd + ["--threads", str(t), "--tmpdir", tdir]
            logger.debug("Running Diamond: %s", " ".join(cmd))
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                start_new_session=True,
                env=env,
            )
            start = time.monotonic()
            effective_timeout = (
                DIAMOND_DEFAULT_TIMEOUT
                if timeout == DEFAULT_SEARCH_TIMEOUT
                else timeout
            )
            last_activity = _read_process_activity(proc.pid)
            last_activity_at = start
            poll_interval = max(1, DIAMOND_WATCHDOG_POLL_INTERVAL)

            while True:
                returncode = proc.poll()
                if returncode is not None:
                    if returncode != 0:
                        raise subprocess.CalledProcessError(returncode, cmd)
                    return

                now = time.monotonic()
                elapsed = now - start
                if effective_timeout > 0 and elapsed >= effective_timeout:
                    _terminate_process_group(proc)
                    raise subprocess.TimeoutExpired(cmd, effective_timeout)

                current_activity = _read_process_activity(proc.pid)
                if current_activity is None:
                    # If /proc is unavailable, fall back to wall-clock timeout.
                    last_activity_at = now
                elif last_activity is None or current_activity != last_activity:
                    last_activity = current_activity
                    last_activity_at = now
                elif (
                    DIAMOND_NO_PROGRESS_TIMEOUT > 0
                    and elapsed >= DIAMOND_NO_PROGRESS_MIN_RUNTIME
                    and now - last_activity_at >= DIAMOND_NO_PROGRESS_TIMEOUT
                ):
                    _terminate_process_group(proc)
                    raise DiamondNoProgressTimeout(
                        cmd,
                        DIAMOND_NO_PROGRESS_TIMEOUT,
                    )

                sleep_for = poll_interval
                if effective_timeout > 0:
                    sleep_for = min(sleep_for, max(0.1, effective_timeout - elapsed))
                time.sleep(sleep_for)

    start_time = time.time()
    try:
        _invoke(threads)
        elapsed = time.time() - start_time
        logger.debug("Diamond completed in %.2f seconds", elapsed)
    except subprocess.TimeoutExpired as e:
        elapsed = time.time() - start_time
        reason = (
            "made no CPU/IO progress"
            if isinstance(e, DiamondNoProgressTimeout)
            else "timed out"
        )
        retry_threads = (
            DIAMOND_RETRY_THREADS
            if DIAMOND_RETRY_THREADS > 0
            else max(1, min(threads, max(8, threads // 2)))
        )
        logger.warning(
            "Diamond blastp %s after %.1fs at threads=%d; retrying once "
            "with threads=%d",
            reason,
            elapsed, threads, retry_threads,
        )
        if output_tsv.exists():
            output_tsv.unlink()
        retry_start = time.time()
        try:
            _invoke(retry_threads)
            logger.info(
                "Diamond blastp retry succeeded in %.1fs (initial attempt "
                "stalled; reduced-thread retry cleared it)",
                time.time() - retry_start,
            )
        except subprocess.TimeoutExpired as retry_error:
            total = time.time() - start_time
            retry_reason = (
                "NO-PROGRESS"
                if isinstance(retry_error, DiamondNoProgressTimeout)
                else "TIMEOUT"
            )
            logger.error(
                "Diamond blastp %s on both attempts (total %.1fs)",
                retry_reason, total,
            )
            raise RuntimeError(
                f"Diamond blastp timeout after {total:.2f}s including retry. "
                "Check worker environment and system resources."
            )
        except subprocess.CalledProcessError as e:
            logger.error(
                "Diamond blastp retry failed with exit code %d", e.returncode,
            )
            raise
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        logger.error(
            "Diamond blastp failed with exit code %d after %.2fs",
            e.returncode, elapsed,
        )
        raise
    return output_tsv
