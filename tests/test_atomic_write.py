from __future__ import annotations

import errno
import tempfile
from pathlib import Path

import pytest

from virosync.utils import atomic_write as atomic_write_module
from virosync.utils.atomic_write import atomic_write


def test_atomic_write_failure_leaves_no_temp_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A failed write (e.g. ENOSPC mid-batch) must not orphan a ``.tmp_*`` file."""
    target = tmp_path / "phase2_complete.json"
    real_named_temporary_file = tempfile.NamedTemporaryFile

    def failing_named_temporary_file(*args, **kwargs):
        handle = real_named_temporary_file(*args, **kwargs)

        def fail_write(_content: str) -> int:
            raise OSError(errno.ENOSPC, "No space left on device")

        handle.write = fail_write
        return handle

    monkeypatch.setattr(
        atomic_write_module.tempfile,
        "NamedTemporaryFile",
        failing_named_temporary_file,
    )

    with pytest.raises(OSError, match="No space left on device"):
        atomic_write(target, '{"status": "complete"}')

    assert not target.exists()
    assert not list(tmp_path.glob(".tmp_*"))
    assert list(tmp_path.iterdir()) == []
