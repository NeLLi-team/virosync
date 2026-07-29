from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from virosync.pipeline.phase3 import structural_homology
from virosync.pipeline.phase3.structural_homology import (
    BoltzFoldSeekAnalyzer,
    FoldSeekSearcher,
)
from virosync.utils.executables import resolve_boltz_executable


def test_resolve_boltz_executable_finds_repo_wrapper(monkeypatch) -> None:
    monkeypatch.setenv("PATH", "")

    resolved = resolve_boltz_executable()

    assert resolved is not None
    assert Path(resolved).name == "boltz"
    assert Path(resolved).parent.name == "scripts"


def test_boltz_analyzer_requires_msa_server_for_sequence_only_inputs(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        structural_homology,
        "resolve_boltz_executable",
        lambda: "/bin/boltz",
    )
    monkeypatch.setattr(
        structural_homology.shutil,
        "which",
        lambda name: f"/bin/{name}",
    )

    assert not BoltzFoldSeekAnalyzer(
        viral_db_path=tmp_path / "bfvd",
        use_msa_server=False,
    ).available()
    assert BoltzFoldSeekAnalyzer(
        viral_db_path=tmp_path / "bfvd",
        use_msa_server=True,
    ).available()


def test_foldseek_search_requests_alignment_backtrace(
    tmp_path: Path,
    monkeypatch,
) -> None:
    calls: list[list[str]] = []

    def _fake_run(cmd, **kwargs):
        calls.append([str(item) for item in cmd])
        if len(cmd) > 1 and cmd[1] == "convertalis":
            Path(cmd[5]).write_text("")
        return SimpleNamespace(stdout="foldseek test", stderr="", returncode=0)

    monkeypatch.setattr(structural_homology.subprocess, "run", _fake_run)
    query = tmp_path / "query.cif"
    query.write_text("data_query\n")

    searcher = FoldSeekSearcher(database_path=tmp_path / "bfvd", threads=1)
    assert searcher.search(query) == []

    search_calls = [cmd for cmd in calls if len(cmd) > 1 and cmd[1] == "search"]
    assert search_calls
    assert "-a" in search_calls[0]
