"""Reuse the Phase 2 gene-taxonomy search of an earlier run instead of repeating it."""

from __future__ import annotations

import json
import logging
import os
import time
from dataclasses import replace
from pathlib import Path

import pytest
from test_phase1_reuse import (
    OLD_CODE,
    _flow,
    _invoke_run,
    _make_worlds,
    _provenance,
    _sha256_of,
    _tree,
    _write_snapshot,
)
from test_phase2_superset_prototype import _write_proteome
from test_run_state_orchestration import _MockPipeline

from virosync.config import PipelineConfig
from virosync.orchestration._flows.single_genome import orchestrator, phase1_reuse, phase2
from virosync.orchestration._flows.single_genome.phase1_reuse import Phase1Reuse, Phase1ReuseError
from virosync.pipeline.host_signatures import HostSignatureModel
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
from virosync.pipeline.phase2 import boundary_diamond as bd
from virosync.pipeline.phase2.boundary_diamond import SavedBoundarySearch, SavedHitTable, SavedSearchError

QUERY = "phase2/boundary_diamond/genome_boundary_query.faa"
HITS = "phase2/boundary_diamond/boundary_diamond.tsv"
CHUNKS = ("phase2/boundary_diamond/diamond_chunk_1.tsv", "phase2/boundary_diamond/diamond_chunk_2.tsv")
SAVED_FILES = (QUERY, HITS)
SEARCH_REUSE = "phase2/boundary_diamond/search_reuse.json"
HIT_LINES = (
    "p1\tNCLDV__virus|x\t1e-30\t120.0\t60.0\t95.0\n",
    "p2\tEUK__host|y\t1e-20\t100.0\t50.0\t90.0\n",
)
HIT_TABLE = "".join(HIT_LINES)
Worlds = tuple[_MockPipeline, _MockPipeline]


@pytest.fixture
def worlds(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Worlds:
    """An earlier finished run that saved its Phase 2 search, and a second pipeline with another code identity."""
    old, new = _make_worlds(tmp_path, monkeypatch)
    for relative, content in ((QUERY, ">p1\nMMM\n>p2\nMMM\n"), (HITS, HIT_TABLE)):
        path = old.output_dir / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    database = _database(old)
    database.parent.mkdir(parents=True)
    database.write_bytes(b"database")
    _write_provenance(old, {"path": str(database), "exists": True, "size_bytes": 8})
    return old, new


def _database(old: _MockPipeline) -> Path:
    """The gene-taxonomy database that the earlier run searched and the new run uses."""
    return old.output_dir.parent.parent / "bundle/genomes/combined_proteome.dmnd"


def _write_provenance(old: _MockPipeline, record: dict | None) -> None:
    """Write the ``provenance.json`` of the earlier run with ``record`` as its gene-taxonomy database."""
    databases = {} if record is None else {"gene_taxonomy_db": record}
    (old.output_dir / "provenance.json").write_text(json.dumps({"databases": databases}), encoding="utf-8")


def _run(old: _MockPipeline, new: _MockPipeline, request: Phase1Reuse, **kwargs) -> dict:
    return _flow(new, request, gene_taxonomy_faa_db=_database(old), **kwargs)


def _request(
    old: _MockPipeline, *, manifest: tuple[Path, str] | None = None, phase2_search: bool = True
) -> Phase1Reuse:
    request = Phase1Reuse(old_base=old.output_dir.parent, code_sha256=OLD_CODE, reuse_phase2_search=phase2_search)
    if manifest is None:
        return request
    return replace(request, snapshot_manifest=manifest[0], snapshot_sha256=manifest[1])


def _listing(tmp_path: Path, old: _MockPipeline, *, hashes: dict[str, str] | None = None) -> tuple[Path, str]:
    """Write a snapshot of the saved Phase 2 files; ``hashes`` replaces the listed hash of a file."""
    hashes = hashes or {}
    folder = old.output_dir / "phase2/boundary_diamond"
    relatives = sorted(path.relative_to(old.output_dir).as_posix() for path in folder.iterdir())
    lines = [f"{hashes.get(rel, _sha256_of(old.output_dir / rel))}  demo/{rel}" for rel in relatives]
    return _write_snapshot(tmp_path, lines)


def _use_chunks(old: _MockPipeline) -> None:
    """Make the earlier search a two-chunk search: one table per query protein."""
    (old.output_dir / HITS).unlink()
    for relative, line in zip(CHUNKS, HIT_LINES, strict=True):
        (old.output_dir / relative).write_text(line, encoding="utf-8")


def _open(
    old: _MockPipeline,
    new: _MockPipeline,
    request: Phase1Reuse,
    *,
    superset: bool = False,
    chunk_size: int = 100,
    database: Path | None = None,
):
    config = PipelineConfig()
    config.databases.gene_taxonomy_faa_db = database or _database(old)
    config.phase2.diamond_superset_prototype_enabled = superset
    config.phase2.diamond_chunk_size = chunk_size
    return phase1_reuse.open_reused_run(request, genome_id="demo", identity=new.identity, config=config)


def _saved_search(old: _MockPipeline, tables: tuple[str, ...] = (HITS,)) -> SavedBoundarySearch:
    return SavedBoundarySearch(
        query_sha256=_sha256_of(old.output_dir / QUERY),
        hit_tables=tuple(SavedHitTable(old.output_dir / rel, _sha256_of(old.output_dir / rel)) for rel in tables),
    )


def _spy_phase2(new: _MockPipeline, monkeypatch: pytest.MonkeyPatch) -> list[dict]:
    """Record the arguments of the Phase 2 call and let the mock run it."""
    seen: list[dict] = []

    def spy(**kwargs):
        seen.append(kwargs)
        return new._phase2(**kwargs)

    monkeypatch.setattr(orchestrator, "_run_phase2_subflow", spy)
    return seen


# --- command line ---------------------------------------------------------------------------------------------


def test_cli_refuses_the_flag_without_the_earlier_folder(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    result, received = _invoke_run(tmp_path, monkeypatch, "--reuse-phase2-search")

    assert result.exit_code != 0
    assert "--reuse-phase2-search: needs --reuse-phase1-from" in result.output
    assert not received


def test_cli_carries_the_flag_in_the_request(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    old_base = tmp_path / "campaign"
    old_base.mkdir()

    result, received = _invoke_run(
        tmp_path,
        monkeypatch,
        *("--reuse-phase1-from", str(old_base), "--reuse-phase1-code-sha256", OLD_CODE),
        "--reuse-phase2-search",
    )

    assert result.exit_code == 0, result.output
    assert received["phase1_reuse"] == Phase1Reuse(old_base=old_base, code_sha256=OLD_CODE, reuse_phase2_search=True)


def test_a_request_does_not_reuse_the_search_by_default(tmp_path: Path) -> None:
    assert Phase1Reuse(old_base=tmp_path, code_sha256=OLD_CODE).reuse_phase2_search is False


# --- authentication at import ---------------------------------------------------------------------------------


def test_a_snapshot_authenticates_the_saved_search(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds

    source = _open(old, new, _request(old, manifest=_listing(tmp_path, old)))

    assert source.phase2_search == _saved_search(old)


@pytest.mark.parametrize("unlisted", SAVED_FILES)
def test_a_saved_file_without_a_hash_is_refused_and_named(worlds: Worlds, tmp_path: Path, unlisted: str) -> None:
    old, new = worlds
    kept = next(rel for rel in SAVED_FILES if rel != unlisted)
    manifest = _write_snapshot(tmp_path, [f"{_sha256_of(old.output_dir / kept)}  demo/{kept}"])

    with pytest.raises(Phase1ReuseError, match=unlisted):
        _open(old, new, _request(old, manifest=manifest))


def test_the_saved_search_needs_a_hash_even_without_a_snapshot(worlds: Worlds) -> None:
    old, new = worlds

    with pytest.raises(Phase1ReuseError, match=QUERY):
        _open(old, new, _request(old))


@pytest.mark.parametrize("changed", SAVED_FILES)
def test_a_listed_hash_that_differs_from_the_file_is_refused(worlds: Worlds, tmp_path: Path, changed: str) -> None:
    old, new = worlds
    manifest = _listing(tmp_path, old, hashes={changed: "0" * 64})

    with pytest.raises(Phase1ReuseError, match=f"authenticated input changed: {changed}"):
        _open(old, new, _request(old, manifest=manifest))


@pytest.mark.parametrize("linked", SAVED_FILES)
def test_a_symlinked_saved_file_is_refused(worlds: Worlds, tmp_path: Path, linked: str) -> None:
    old, new = worlds
    manifest = _listing(tmp_path, old)
    path = old.output_dir / linked
    target = tmp_path / "elsewhere"
    target.write_bytes(path.read_bytes())
    path.unlink()
    path.symlink_to(target)

    with pytest.raises(Phase1ReuseError, match="not a regular file"):
        _open(old, new, _request(old, manifest=manifest))


def test_a_symlinked_search_folder_is_refused(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    manifest = _listing(tmp_path, old)
    folder = old.output_dir / "phase2/boundary_diamond"
    moved = tmp_path / "moved"
    folder.rename(moved)
    folder.symlink_to(moved)

    with pytest.raises(Phase1ReuseError, match="symlink"):
        _open(old, new, _request(old, manifest=manifest))


@pytest.mark.parametrize("missing", SAVED_FILES)
def test_an_earlier_run_with_only_one_saved_file_is_refused(worlds: Worlds, tmp_path: Path, missing: str) -> None:
    old, new = worlds
    manifest = _listing(tmp_path, old)
    (old.output_dir / missing).unlink()

    with pytest.raises(Phase1ReuseError, match=missing):
        _open(old, new, _request(old, manifest=manifest))


def test_a_snapshot_authenticates_a_chunked_saved_search(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    _use_chunks(old)

    source = _open(old, new, _request(old, manifest=_listing(tmp_path, old)), chunk_size=1)

    assert source.phase2_search == _saved_search(old, CHUNKS)


def test_a_missing_chunk_table_is_refused(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    _use_chunks(old)
    manifest = _listing(tmp_path, old)
    (old.output_dir / CHUNKS[1]).unlink()

    with pytest.raises(Phase1ReuseError, match=f"lacks {CHUNKS[1]}"):
        _open(old, new, _request(old, manifest=manifest), chunk_size=1)


@pytest.mark.parametrize("extra", ["diamond_chunk_3.tsv", "boundary_diamond.tsv"], ids=["beyond-the-last", "single"])
def test_a_hit_table_outside_the_expected_chunk_set_is_refused(worlds: Worlds, tmp_path: Path, extra: str) -> None:
    old, new = worlds
    _use_chunks(old)
    (old.output_dir / "phase2/boundary_diamond" / extra).write_text(HIT_TABLE, encoding="utf-8")
    manifest = _listing(tmp_path, old)

    with pytest.raises(Phase1ReuseError, match=f"outside the expected set.*{extra}"):
        _open(old, new, _request(old, manifest=manifest), chunk_size=1)


def test_a_chunk_table_beside_a_single_table_is_refused(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    (old.output_dir / CHUNKS[0]).write_text(HIT_TABLE, encoding="utf-8")
    manifest = _listing(tmp_path, old)

    with pytest.raises(Phase1ReuseError, match=f"outside the expected set.*{CHUNKS[0]}"):
        _open(old, new, _request(old, manifest=manifest))


def test_a_chunk_table_without_a_hash_is_refused_and_named(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    _use_chunks(old)
    lines = [f"{_sha256_of(old.output_dir / rel)}  demo/{rel}" for rel in (QUERY, CHUNKS[0])]

    with pytest.raises(Phase1ReuseError, match=CHUNKS[1]):
        _open(old, new, _request(old, manifest=_write_snapshot(tmp_path, lines)), chunk_size=1)


def test_a_chunk_table_listing_authenticates_nothing_without_the_flag(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    _use_chunks(old)

    source = _open(old, new, _request(old, manifest=_listing(tmp_path, old), phase2_search=False), chunk_size=1)

    assert source.snapshot_authenticated == {}


def test_an_earlier_run_without_a_saved_search_is_imported_and_records_that(worlds: Worlds) -> None:
    old, new = worlds
    for relative in SAVED_FILES:
        (old.output_dir / relative).unlink()

    source = _open(old, new, _request(old))

    assert source.phase2_search is None
    assert source.phase2_search_requested is True


def test_the_full_proteome_variant_is_refused(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds

    with pytest.raises(Phase1ReuseError, match="diamond_superset_prototype_enabled"):
        _open(old, new, _request(old, manifest=_listing(tmp_path, old)), superset=True)


def test_the_full_proteome_variant_is_accepted_without_the_flag(worlds: Worlds) -> None:
    old, new = worlds

    source = _open(old, new, _request(old, phase2_search=False), superset=True)

    assert source.phase2_search is None


# --- search database --------------------------------------------------------------------------------------------


def test_the_same_database_in_another_bundle_is_accepted(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    moved = tmp_path / "other-bundle/genomes/combined_proteome.dmnd"
    moved.parent.mkdir(parents=True)
    moved.write_bytes(_database(old).read_bytes())

    source = _open(old, new, _request(old, manifest=_listing(tmp_path, old)), database=moved)

    assert source.phase2_search == _saved_search(old)


def test_a_relative_database_path_naming_the_recorded_file_is_accepted(
    worlds: Worlds, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    old, new = worlds
    monkeypatch.chdir(_database(old).parent)

    source = _open(old, new, _request(old, manifest=_listing(tmp_path, old)), database=Path("combined_proteome.dmnd"))

    assert source.phase2_search == _saved_search(old)


def test_a_database_with_another_name_is_refused(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    other = tmp_path / "bundle/marker/marker.dmnd"
    other.parent.mkdir(parents=True)
    other.write_bytes(_database(old).read_bytes())

    with pytest.raises(
        Phase1ReuseError, match=r"database differs.*marker/marker\.dmnd.*genomes/combined_proteome\.dmnd"
    ):
        _open(old, new, _request(old, manifest=_listing(tmp_path, old)), database=other)


def test_a_database_with_another_size_is_refused(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    _database(old).write_bytes(b"a larger database")

    with pytest.raises(Phase1ReuseError, match=r"database differs.*17 bytes.*8 bytes"):
        _open(old, new, _request(old, manifest=_listing(tmp_path, old)))


def test_a_new_run_without_a_database_is_refused(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    config = PipelineConfig()

    with pytest.raises(Phase1ReuseError, match="no gene-taxonomy database"):
        phase1_reuse.open_reused_run(
            _request(old, manifest=_listing(tmp_path, old)), genome_id="demo", identity=new.identity, config=config
        )


def test_a_missing_provenance_file_is_refused(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    manifest = _listing(tmp_path, old)
    (old.output_dir / "provenance.json").unlink()

    with pytest.raises(Phase1ReuseError, match="no regular provenance.json"):
        _open(old, new, _request(old, manifest=manifest))


def test_a_symlinked_provenance_file_is_refused(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    manifest = _listing(tmp_path, old)
    path = old.output_dir / "provenance.json"
    target = tmp_path / "elsewhere.json"
    target.write_bytes(path.read_bytes())
    path.unlink()
    path.symlink_to(target)

    with pytest.raises(Phase1ReuseError, match="no regular provenance.json"):
        _open(old, new, _request(old, manifest=manifest))


def test_a_provenance_file_that_is_not_json_is_refused(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    manifest = _listing(tmp_path, old)
    (old.output_dir / "provenance.json").write_text("{", encoding="utf-8")

    with pytest.raises(Phase1ReuseError, match="cannot read"):
        _open(old, new, _request(old, manifest=manifest))


@pytest.mark.parametrize("record", [None, {"path": "x/genomes/combined_proteome.dmnd"}, {"size_bytes": 8}])
def test_a_provenance_file_without_the_database_record_is_refused(
    worlds: Worlds, tmp_path: Path, record: dict | None
) -> None:
    old, new = worlds
    manifest = _listing(tmp_path, old)
    _write_provenance(old, record)

    with pytest.raises(Phase1ReuseError, match="lacks the record databases.gene_taxonomy_db"):
        _open(old, new, _request(old, manifest=manifest))


def test_the_database_is_not_compared_without_the_flag(worlds: Worlds) -> None:
    old, new = worlds
    (old.output_dir / "provenance.json").unlink()

    assert _open(old, new, _request(old, phase2_search=False)).phase2_search is None


# --- Phase 2 marker of the earlier run -----------------------------------------------------------------------


def _write_phase2_marker(old: _MockPipeline, listed: dict[str, str]) -> None:
    artifacts = [{"relative_path": rel, "sha256": sha} for rel, sha in listed.items()]
    (old.output_dir / "phase2.complete.json").write_text(json.dumps({"artifacts": artifacts}), encoding="utf-8")


def test_a_phase2_marker_hash_that_differs_from_the_listed_hash_is_refused(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    manifest = _listing(tmp_path, old)
    _write_phase2_marker(old, {HITS: "0" * 64})

    with pytest.raises(Phase1ReuseError, match=f"snapshot conflicts with the recorded hash: {HITS}"):
        _open(old, new, _request(old, manifest=manifest))


def test_a_phase2_marker_hash_that_equals_the_listed_hash_is_accepted(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    manifest = _listing(tmp_path, old)
    _write_phase2_marker(old, {rel: _sha256_of(old.output_dir / rel) for rel in SAVED_FILES})

    assert _open(old, new, _request(old, manifest=manifest)).phase2_search == _saved_search(old)


def test_a_phase2_marker_that_lists_a_file_twice_is_refused_when_any_entry_differs(
    worlds: Worlds, tmp_path: Path
) -> None:
    old, new = worlds
    manifest = _listing(tmp_path, old)
    right = _sha256_of(old.output_dir / HITS)
    artifacts = [{"relative_path": HITS, "sha256": sha} for sha in ("0" * 64, right)]
    (old.output_dir / "phase2.complete.json").write_text(json.dumps({"artifacts": artifacts}), encoding="utf-8")

    with pytest.raises(Phase1ReuseError, match=f"snapshot conflicts with the recorded hash: {HITS}"):
        _open(old, new, _request(old, manifest=manifest))


def test_a_phase2_marker_never_authenticates_a_file(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    manifest = _write_snapshot(tmp_path, [f"{_sha256_of(old.output_dir / QUERY)}  demo/{QUERY}"])
    _write_phase2_marker(old, {HITS: _sha256_of(old.output_dir / HITS)})

    with pytest.raises(Phase1ReuseError, match=f"no authenticated hash for {HITS}"):
        _open(old, new, _request(old, manifest=manifest))


@pytest.mark.parametrize(
    "content",
    [
        "{",
        "[]",
        '{"artifacts": {}}',
        '{"artifacts": [{"relative_path": "x"}]}',
        '{"artifacts": ["x"]}',
        '{"artifacts": [{"relative_path": 1, "sha256": "x"}]}',
        '{"artifacts": [{"relative_path": "x", "sha256": 1}]}',
    ],
)
def test_a_malformed_phase2_marker_is_refused(worlds: Worlds, tmp_path: Path, content: str) -> None:
    old, new = worlds
    manifest = _listing(tmp_path, old)
    (old.output_dir / "phase2.complete.json").write_text(content, encoding="utf-8")

    with pytest.raises(Phase1ReuseError, match="Phase 2 marker"):
        _open(old, new, _request(old, manifest=manifest))


def test_a_symlinked_phase2_marker_is_refused(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    manifest = _listing(tmp_path, old)
    target = tmp_path / "marker.json"
    target.write_text('{"artifacts": []}', encoding="utf-8")
    (old.output_dir / "phase2.complete.json").unlink()
    (old.output_dir / "phase2.complete.json").symlink_to(target)

    with pytest.raises(Phase1ReuseError, match="Phase 2 marker"):
        _open(old, new, _request(old, manifest=manifest))


# --- set boundaries and genome ids ---------------------------------------------------------------------------


def test_a_query_with_exactly_the_chunk_size_imports_the_single_table(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds

    source = _open(old, new, _request(old, manifest=_listing(tmp_path, old)), chunk_size=2)

    assert source.phase2_search == _saved_search(old)


def test_a_query_with_twice_the_chunk_size_imports_two_chunk_tables_and_refuses_a_third(
    worlds: Worlds, tmp_path: Path
) -> None:
    old, new = worlds
    _use_chunks(old)
    (old.output_dir / QUERY).write_text(">p1\nM\n>p2\nM\n>p3\nM\n>p4\nM\n", encoding="utf-8")
    manifest = _listing(tmp_path, old)

    assert _open(old, new, _request(old, manifest=manifest), chunk_size=2).phase2_search == _saved_search(old, CHUNKS)

    third = old.output_dir / "phase2/boundary_diamond/diamond_chunk_3.tsv"
    third.write_text(HIT_TABLE, encoding="utf-8")
    with pytest.raises(Phase1ReuseError, match=f"outside the expected set.*{third.name}"):
        _open(old, new, _request(old, manifest=_listing(tmp_path, old)), chunk_size=2)


def test_a_snapshot_line_of_a_genome_whose_id_is_a_prefix_authenticates_nothing(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    lines = [f"{_sha256_of(old.output_dir / rel)}  dem/{rel}" for rel in SAVED_FILES]

    with pytest.raises(Phase1ReuseError, match="no authenticated hash"):
        _open(old, new, _request(old, manifest=_write_snapshot(tmp_path, lines)))


# --- snapshot compatibility -----------------------------------------------------------------------------------


def test_a_listing_with_phase2_lines_authenticates_nothing_without_the_flag(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds

    source = _open(old, new, _request(old, manifest=_listing(tmp_path, old), phase2_search=False))

    assert source.phase2_search is None
    assert source.snapshot_authenticated == {}
    assert not any(relative.startswith("phase2/") for relative in source.artifacts)


def test_a_listing_with_phase2_lines_keeps_the_phase1_provenance_without_the_flag(
    worlds: Worlds, tmp_path: Path
) -> None:
    old, new = worlds

    _run(old, new, _request(old, manifest=_listing(tmp_path, old), phase2_search=False))

    assert "snapshot_authenticated" not in _provenance(new, 1)


# --- the import in the new run --------------------------------------------------------------------------------


def test_the_flow_hands_the_saved_search_to_phase2(worlds: Worlds, tmp_path: Path, monkeypatch) -> None:
    old, new = worlds
    seen = _spy_phase2(new, monkeypatch)

    result = _run(old, new, _request(old, manifest=_listing(tmp_path, old)))

    assert result["success"] is True
    assert seen[0]["saved_boundary_search"] == _saved_search(old)
    assert seen[0]["require_saved_boundary_search"] is True


def test_the_flow_without_the_flag_hands_phase2_no_saved_search(worlds: Worlds, tmp_path: Path, monkeypatch) -> None:
    old, new = worlds
    seen = _spy_phase2(new, monkeypatch)

    _run(old, new, _request(old, manifest=_listing(tmp_path, old), phase2_search=False))

    assert seen[0]["saved_boundary_search"] is None
    assert seen[0]["require_saved_boundary_search"] is False
    assert not (new.output_dir / SEARCH_REUSE).exists()


def test_the_flow_without_saved_files_requires_a_saved_search_in_phase2(worlds: Worlds, monkeypatch) -> None:
    old, new = worlds
    for relative in SAVED_FILES:
        (old.output_dir / relative).unlink()
    seen = _spy_phase2(new, monkeypatch)

    result = _run(old, new, _request(old))

    assert result["success"] is True
    assert seen[0]["saved_boundary_search"] is None
    assert seen[0]["require_saved_boundary_search"] is True
    assert not (new.output_dir / SEARCH_REUSE).exists()


def test_the_import_record_names_the_files_and_the_snapshot(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    manifest = _listing(tmp_path, old)

    _run(old, new, _request(old, manifest=manifest))

    assert json.loads((new.output_dir / SEARCH_REUSE).read_text()) == {
        "schema": "virosync.phase2_search_reuse/v1",
        "old_run_folder": str(old.output_dir.resolve()),
        "files": [
            {"path": QUERY, "sha256": _sha256_of(old.output_dir / QUERY), "authenticated_by": "snapshot"},
            {"path": HITS, "sha256": _sha256_of(old.output_dir / HITS), "authenticated_by": "snapshot"},
        ],
        "snapshot_manifest_sha256": manifest[1],
    }


def test_the_import_record_lists_every_chunk_table(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    _use_chunks(old)

    _run(old, new, _request(old, manifest=_listing(tmp_path, old)), boundary_diamond_chunk_size=1)

    record = json.loads((new.output_dir / SEARCH_REUSE).read_text())
    assert [item["path"] for item in record["files"]] == [QUERY, *CHUNKS]


def _imported_source(old: _MockPipeline) -> phase1_reuse.ReusedRun:
    """An opened earlier run with its saved search."""
    return phase1_reuse.ReusedRun(
        root=old.output_dir,
        run_fingerprint="1" * 64,
        code_sha256=OLD_CODE,
        markers={},
        snapshot_sha256="2" * 64,
        phase2_search_requested=True,
        phase2_search=_saved_search(old),
    )


def test_the_import_is_logged(worlds: Worlds, tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    old, _ = worlds
    caplog.set_level(logging.INFO)

    phase1_reuse.record_phase2_search_import(_imported_source(old), tmp_path / "new", logging.getLogger(__name__))

    assert [record.getMessage() for record in caplog.records] == [
        f"Phase 2 search imported from {old.output_dir.resolve()}"
    ]


@pytest.mark.parametrize("phase", [0, 1])
def test_the_provenance_records_say_that_the_import_was_requested(worlds: Worlds, tmp_path: Path, phase: int) -> None:
    old, new = worlds

    _run(old, new, _request(old, manifest=_listing(tmp_path, old)))

    assert _provenance(new, phase)["phase2_search_reuse_requested"] is True


@pytest.mark.parametrize("phase", [0, 1])
def test_the_provenance_records_have_no_request_field_without_the_flag(worlds: Worlds, phase: int) -> None:
    old, new = worlds

    _run(old, new, _request(old, phase2_search=False))

    assert "phase2_search_reuse_requested" not in _provenance(new, phase)


def test_the_import_leaves_the_earlier_folder_unchanged(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    manifest = _listing(tmp_path, old)
    before = _tree(old.output_dir.parent)

    _run(old, new, _request(old, manifest=manifest))

    assert _tree(old.output_dir.parent) == before


def test_a_resumed_run_imports_again(worlds: Worlds, tmp_path: Path) -> None:
    old, new = worlds
    request = _request(old, manifest=_listing(tmp_path, old))
    _run(old, new, request)

    assert _run(old, new, request, resume=True)["success"] is True
    assert (new.output_dir / SEARCH_REUSE).is_file()


# --- the search in Phase 2 ------------------------------------------------------------------------------------


def _bd_inputs(tmp_path: Path) -> dict:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">p1 # 1 # 90 # + # ID=contig_1;\nMPEPTIDE\n>p2 # 101 # 190 # + # ID=contig_2;\nMPEPTIDE\n", encoding="utf-8"
    )
    return {
        "query": bd.GenomeDiamondQuery(eve_porf_ids={"seed": ["p1", "p2"]}, all_porf_ids=["p1", "p2"]),
        "proteome_fasta": proteome,
        "diamond_db": tmp_path / "db.dmnd",
        "proteome_index": {"contig": [bd.pORF("p1", "contig", 0, 90), bd.pORF("p2", "contig", 100, 190)]},
    }


def _run_search(tmp_path: Path, output: str, *, chunk_size: int = 10, **extra):
    return bd.run_batched_diamond(
        output_dir=tmp_path / output,
        config=bd.BoundaryDiamondConfig(chunk_size=chunk_size, threads=1),
        **_bd_inputs(tmp_path),
        **extra,
    )


def _fake_search(table: str):
    def search(**kwargs) -> None:
        Path(kwargs["output"]).write_text(table, encoding="utf-8")

    return search


def _forbid_search(monkeypatch: pytest.MonkeyPatch) -> None:
    def fail(**_kwargs):
        raise AssertionError("a DIAMOND search started")

    monkeypatch.setattr(bd, "run_diamond_blastp", fail)
    monkeypatch.setattr(bd, "run_diamond_chunked", fail)


def _earlier_search(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[dict, SavedBoundarySearch]:
    """Run the search once and describe its saved files."""
    monkeypatch.setattr(bd, "run_diamond_blastp", _fake_search(HIT_TABLE))
    taxonomy = _run_search(tmp_path, "earlier")
    saved = SavedBoundarySearch(
        query_sha256=_sha256_of(tmp_path / "earlier/genome_boundary_query.faa"),
        hit_tables=(
            SavedHitTable(
                tmp_path / "earlier/boundary_diamond.tsv", _sha256_of(tmp_path / "earlier/boundary_diamond.tsv")
            ),
        ),
    )
    return taxonomy, saved


def test_a_saved_search_gives_the_classification_of_the_search_without_a_search(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    searched, saved = _earlier_search(tmp_path, monkeypatch)
    _forbid_search(monkeypatch)

    imported = _run_search(tmp_path, "later", saved_search=saved)

    assert imported == searched
    assert (tmp_path / "later/boundary_diamond.tsv").read_text() == HIT_TABLE


def test_a_query_that_differs_from_the_saved_query_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _, saved = _earlier_search(tmp_path, monkeypatch)
    _forbid_search(monkeypatch)

    with pytest.raises(SavedSearchError, match="query proteins differ from the earlier run"):
        _run_search(tmp_path, "later", saved_search=replace(saved, query_sha256="0" * 64))

    assert not (tmp_path / "later/boundary_diamond.tsv").exists()


def _earlier_chunked_search(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[dict, SavedBoundarySearch]:
    """Run a two-chunk search once, with one hit per chunk, and describe its saved files."""
    tables = iter(HIT_LINES)

    def chunk_search(**kwargs) -> None:
        Path(kwargs["output"]).write_text(next(tables), encoding="utf-8")

    monkeypatch.setattr(bd, "run_diamond_blastp", chunk_search)
    taxonomy = _run_search(tmp_path, "earlier", chunk_size=1)
    folder = tmp_path / "earlier"
    saved = SavedBoundarySearch(
        query_sha256=_sha256_of(folder / "genome_boundary_query.faa"),
        hit_tables=tuple(
            SavedHitTable(folder / name, _sha256_of(folder / name))
            for name in ("diamond_chunk_1.tsv", "diamond_chunk_2.tsv")
        ),
    )
    return taxonomy, saved


def test_a_saved_chunked_search_gives_the_merged_classification_without_a_search(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    searched, saved = _earlier_chunked_search(tmp_path, monkeypatch)
    _forbid_search(monkeypatch)
    monkeypatch.setattr(bd, "split_fasta", lambda *_args, **_kwargs: pytest.fail("the query was split"))

    imported = _run_search(tmp_path, "later", chunk_size=1, saved_search=saved)

    assert imported == searched
    assert [(tmp_path / "later" / f"diamond_chunk_{n}.tsv").read_text() for n in (1, 2)] == list(HIT_LINES)
    assert not (tmp_path / "later/boundary_diamond.tsv").exists()


def test_saved_tables_that_do_not_fit_the_query_and_chunk_size_are_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, saved = _earlier_search(tmp_path, monkeypatch)
    _forbid_search(monkeypatch)

    with pytest.raises(SavedSearchError, match="do not match diamond_chunk_1.tsv, diamond_chunk_2.tsv"):
        _run_search(tmp_path, "later", chunk_size=1, saved_search=saved)

    assert not list((tmp_path / "later").glob("*.tsv"))


def test_saved_chunk_tables_that_do_not_fit_a_single_search_are_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, saved = _earlier_chunked_search(tmp_path, monkeypatch)
    _forbid_search(monkeypatch)

    with pytest.raises(SavedSearchError, match="do not match boundary_diamond.tsv"):
        _run_search(tmp_path, "later", chunk_size=10, saved_search=saved)


def test_a_hit_table_that_differs_from_its_authenticated_hash_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, saved = _earlier_search(tmp_path, monkeypatch)
    tampered = replace(saved, hit_tables=(replace(saved.hit_tables[0], sha256="0" * 64),))
    _forbid_search(monkeypatch)

    with pytest.raises(SavedSearchError, match="not the authenticated"):
        _run_search(tmp_path, "later", saved_search=tampered)

    assert not (tmp_path / "later/boundary_diamond.tsv").exists()


def test_a_saved_table_replaced_by_a_symlink_is_refused_at_the_copy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, saved = _earlier_search(tmp_path, monkeypatch)
    table = saved.hit_tables[0].path
    identical = tmp_path / "identical.tsv"
    identical.write_bytes(table.read_bytes())
    table.unlink()
    table.symlink_to(identical)
    _forbid_search(monkeypatch)

    with pytest.raises(SavedSearchError, match="saved Phase 2 hit table"):
        _run_search(tmp_path, "later", saved_search=saved)

    assert not (tmp_path / "later/boundary_diamond.tsv").exists()


def _saved_in_genome_layout(tmp_path: Path, saved: SavedBoundarySearch) -> SavedBoundarySearch:
    """Copy the saved tables to ``genome/phase2/boundary_diamond/`` and point the search there."""
    folder = tmp_path / "genome/phase2/boundary_diamond"
    folder.mkdir(parents=True)
    for table in saved.hit_tables:
        (folder / table.path.name).write_bytes(table.path.read_bytes())
    return replace(saved, hit_tables=tuple(replace(table, path=folder / table.path.name) for table in saved.hit_tables))


@pytest.mark.parametrize("replaced", ["genome/phase2/boundary_diamond", "genome/phase2"])
def test_a_saved_table_whose_parent_folder_was_replaced_by_a_symlink_is_refused_at_the_copy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, replaced: str
) -> None:
    _, saved = _earlier_search(tmp_path, monkeypatch)
    saved = _saved_in_genome_layout(tmp_path, saved)
    folder = tmp_path / replaced
    identical = tmp_path / "identical-folder"
    folder.rename(identical)
    folder.symlink_to(identical)
    _forbid_search(monkeypatch)

    with pytest.raises(SavedSearchError, match="cannot open the saved Phase 2 hit table"):
        _run_search(tmp_path, "later", saved_search=saved)

    assert not (tmp_path / "later/boundary_diamond.tsv").exists()


def test_a_fifo_in_place_of_a_saved_table_is_refused_without_blocking(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, saved = _earlier_search(tmp_path, monkeypatch)
    table = saved.hit_tables[0].path
    table.unlink()
    os.mkfifo(table)
    _forbid_search(monkeypatch)

    with pytest.raises(SavedSearchError, match="is not a regular file"):
        _run_search(tmp_path, "later", saved_search=saved)


def test_a_saved_table_that_is_not_a_regular_file_is_refused_at_the_copy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, saved = _earlier_search(tmp_path, monkeypatch)
    table = saved.hit_tables[0].path
    table.unlink()
    table.mkdir()
    _forbid_search(monkeypatch)

    with pytest.raises(SavedSearchError, match="saved Phase 2 hit table"):
        _run_search(tmp_path, "later", saved_search=saved)


def test_without_a_saved_search_the_search_runs_and_writes_its_table(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    searches: list[Path] = []

    def search(**kwargs) -> None:
        searches.append(Path(kwargs["output"]))
        _fake_search(HIT_TABLE)(**kwargs)

    monkeypatch.setattr(bd, "run_diamond_blastp", search)

    _run_search(tmp_path, "plain")

    assert searches == [tmp_path / "plain/boundary_diamond.tsv"]
    assert (tmp_path / "plain/boundary_diamond.tsv").read_text() == HIT_TABLE


def _phase2_result(output_dir: Path, tmp_path: Path, **kwargs) -> phase2.Phase2Result | phase2.Phase2Terminal:
    masked = tmp_path / "masked.fna"
    masked.write_text(">contig\n" + "ACGT" * 200 + "\n")
    proteome = tmp_path / "proteome.faa"
    _write_proteome(proteome)
    database = tmp_path / "combined.dmnd"
    database.write_bytes(b"database")
    config = PipelineConfig().with_overrides(
        resume=False,
        gene_taxonomy_faa_db=database,
        host_prefixes=["EUK__"],
        boundary_host_trim_score_threshold=0.5,
        boundary_taxonomy_ml_neighbor_window=1,
        threads=1,
        gene_taxonomy_threads=1,
        extended_output=False,
    )
    config.phase2.diamond_flank_genes = 1
    config.phase2.diamond_control_sample_size = 1
    config.phase2.diamond_control_min_distance = 1
    config.phase2.diamond_chunk_size = 100
    seed = MergedSeed(
        scaffold="contig",
        start=100,
        end=350,
        seed_id="seed-1",
        sources=["hhg", "marker_validation"],
        confidence="high",
        confirmed_seeds=(ConfirmedSeed("seed-1", "contig", 100, 350),),
    )
    return phase2._run_phase2_subflow(
        masked_path=masked,
        proteome_path=proteome,
        merged_seeds=[seed],
        validated_markers=[],
        host_signature_model=HostSignatureModel(token_weights={"host": 1.0}, lineage_weights={"euk|host": 1.0}),
        output_dir=output_dir,
        genome_id="genome",
        config=config,
        genome_start_time=time.time(),
        logger=logging.getLogger(__name__),
        metadata_families={},
        **kwargs,
    )


def _host_hit_for_every_query(**kwargs) -> None:
    ids = [line[1:].strip() for line in Path(kwargs["query"]).read_text().splitlines() if line.startswith(">")]
    Path(kwargs["output"]).write_text("".join(f"{i}\tEUK__host|p\t1e-20\t100.0\t80.0\t90.0\n" for i in ids))


def test_phase2_with_a_saved_search_matches_phase2_with_the_search(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(bd, "run_diamond_blastp", _host_hit_for_every_query)
    searched = _phase2_result(tmp_path / "earlier", tmp_path)
    folder = tmp_path / "earlier/phase2/boundary_diamond"
    saved = SavedBoundarySearch(
        query_sha256=_sha256_of(folder / "genome_boundary_query.faa"),
        hit_tables=(SavedHitTable(folder / "boundary_diamond.tsv", _sha256_of(folder / "boundary_diamond.tsv")),),
    )
    _forbid_search(monkeypatch)

    imported = _phase2_result(
        tmp_path / "later", tmp_path, saved_boundary_search=saved, require_saved_boundary_search=True
    )

    assert imported.boundary_taxonomy_map == searched.boundary_taxonomy_map
    assert imported.boundary_control_stats == searched.boundary_control_stats


def test_phase2_that_requires_a_saved_search_stops_when_it_has_none(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _forbid_search(monkeypatch)

    with pytest.raises(Phase1ReuseError, match="saved no Phase 2 search"):
        _phase2_result(tmp_path / "later", tmp_path, require_saved_boundary_search=True)
