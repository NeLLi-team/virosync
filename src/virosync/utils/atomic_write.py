"""
Atomic file write utilities to prevent corruption on crash.

Usage:
    from virosync.utils.atomic_write import atomic_write

    # Write text file atomically
    atomic_write(output_path, content_string)

    # Write with context manager
    with atomic_write_context(output_path) as f:
        f.write("content")
"""

import logging
import tempfile
from pathlib import Path
from typing import Union

logger = logging.getLogger(__name__)


def atomic_write(path: Union[Path, str], content: str, encoding: str = "utf-8") -> None:
    """
    Write file atomically to prevent corruption on crash.

    Writes to a temporary file in the same directory, then atomically
    renames it to the target path. This ensures the file is either
    fully written or not written at all.

    Args:
        path: Target file path
        content: Content to write
        encoding: Character encoding (default: utf-8)
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Create temp file in same directory (important for atomic rename)
    tmp = tempfile.NamedTemporaryFile(
        mode="w",
        encoding=encoding,
        dir=path.parent,
        delete=False,
        prefix=f".tmp_{path.name}_",
        suffix=".tmp",
    )
    tmp_path = Path(tmp.name)

    try:
        # Write and close inside the try so a failed write (e.g. ENOSPC) or a
        # failed flush at close() still cleans up the temp file
        with tmp:
            tmp.write(content)

        # Atomic rename (POSIX guarantees atomicity)
        tmp_path.rename(path)
    except Exception:
        # Clean up temp file on failure
        tmp_path.unlink(missing_ok=True)
        raise


class atomic_write_context:
    """
    Context manager for atomic file writes.

    Usage:
        with atomic_write_context(output_path) as f:
            f.write("line 1\n")
            f.write("line 2\n")
    """

    def __init__(self, path: Union[Path, str], mode: str = "w", encoding: str = "utf-8"):
        self.path = Path(path)
        self.mode = mode
        self.encoding = encoding
        self.tmp_file = None
        self.tmp_path = None

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)

        # Create temp file
        self.tmp_file = tempfile.NamedTemporaryFile(
            mode=self.mode,
            encoding=self.encoding if "b" not in self.mode else None,
            dir=self.path.parent,
            delete=False,
            prefix=f".tmp_{self.path.name}_",
            suffix=".tmp",
        )
        self.tmp_path = Path(self.tmp_file.name)
        return self.tmp_file

    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            self.tmp_file.close()

            if exc_type is None:
                # Success - atomic rename
                self.tmp_path.rename(self.path)
            else:
                # Failure - clean up temp file
                self.tmp_path.unlink(missing_ok=True)
        except Exception as e:
            logger.error("Failed to complete atomic write for %s: %s", self.path, e)
            self.tmp_path.unlink(missing_ok=True)
            raise

        return False  # Don't suppress exceptions
