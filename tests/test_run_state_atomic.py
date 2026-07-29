from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from virosync.orchestration._flows.single_genome import run_state


def _payload(value: str) -> dict[str, object]:
    return {"schema_version": 3, "value": value}


@pytest.mark.parametrize("preexisting", [False, True])
def test_atomic_json_replace_failure_leaves_old_valid_or_absent(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    preexisting: bool,
) -> None:
    target = tmp_path / "state.json"
    if preexisting:
        run_state.atomic_write_json(target, _payload("old"))

    def fail_replace(source, destination):
        raise OSError("injected before replace")

    monkeypatch.setattr(os, "replace", fail_replace)
    with pytest.raises(OSError, match="before replace"):
        run_state.atomic_write_json(target, _payload("new"))

    if preexisting:
        assert json.loads(target.read_text()) == _payload("old")
    else:
        assert not target.exists()
    assert not list(tmp_path.glob(".state.json.*.tmp"))


def test_atomic_json_directory_fsync_failure_leaves_new_valid_state(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "state.json"
    run_state.atomic_write_json(target, _payload("old"))
    real_fsync = os.fsync
    calls = 0

    def fail_directory_fsync(descriptor: int) -> None:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected after replace")
        real_fsync(descriptor)

    monkeypatch.setattr(os, "fsync", fail_directory_fsync)
    with pytest.raises(OSError, match="after replace"):
        run_state.atomic_write_json(target, _payload("new"))

    assert json.loads(target.read_text()) == _payload("new")
    assert not list(tmp_path.glob(".state.json.*.tmp"))


def test_atomic_json_serialization_failure_never_touches_existing_state(
    tmp_path: Path,
) -> None:
    target = tmp_path / "state.json"
    run_state.atomic_write_json(target, _payload("old"))

    with pytest.raises(ValueError, match="NaN"):
        run_state.atomic_write_json(target, {"invalid": float("nan")})

    assert json.loads(target.read_text()) == _payload("old")


def test_phase_artifact_fsyncs_file_and_directory_chain(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "result"
    state_path = root / "phase2" / "nested" / "resume_state.json"
    state_path.parent.mkdir(parents=True)
    state_path.write_text('{"schema_version": 1}\n', encoding="utf-8")
    artifact = run_state.build_artifact_identity(
        state_path,
        root=root,
        schema="phase2-resume-state-v1",
    )
    synced: list[Path] = []
    real_fsync = os.fsync

    def _record_fsync(descriptor: int) -> None:
        synced.append(Path(os.readlink(f"/proc/self/fd/{descriptor}")))
        real_fsync(descriptor)

    monkeypatch.setattr(os, "fsync", _record_fsync)

    run_state._fsync_artifact(root, artifact)

    assert synced == [
        state_path,
        state_path.parent,
        state_path.parent.parent,
        root,
    ]
