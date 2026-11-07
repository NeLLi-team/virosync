"""Keep native marker search depth and provenance aligned with configuration."""

import inspect
import json
from pathlib import Path

import pytest

from virosync.config import PipelineConfig
from virosync.orchestration._flows.single_genome.orchestrator import single_genome_flow
from virosync.orchestration.tasks import hhg_seeding_task
from virosync.pipeline.phase1 import hhg_seeding
from virosync.utils import provenance


def test_hhg_default_matches_effective_config():
    assert (
        inspect.signature(hhg_seeding_task).parameters["marker_top_k"].default
        == PipelineConfig().phase1.marker_validation_top_k
        == 20
    )
    assert inspect.signature(single_genome_flow).parameters["marker_validation_top_k"].default == 20


@pytest.mark.parametrize("top_k", [7, 20, 35])
def test_legacy_marker_search_forwards_requested_depth(tmp_path, monkeypatch, top_k):
    from virosync.pipeline import search_backend

    target = tmp_path / "marker.faa"
    target.write_text(">ref\nM\n")
    (tmp_path / "tmp").mkdir()
    (tmp_path / "tmp" / "marker.dmnd").touch()
    calls = []
    monkeypatch.setattr(search_backend, "run_genome_filtered_search", lambda **kw: calls.append(kw))
    hhg_seeding._run_diamond_marker_search(
        query_fasta=Path("query.faa"), target_faa=target, output_file=tmp_path / "hits.tsv", threads=1, max_seqs=top_k
    )
    assert calls[0]["max_target_seqs"] == top_k


@pytest.mark.parametrize("config, expected", [({}, 20), ({"marker_validation_top_k": 35}, 35)])
def test_provenance_records_marker_depth(tmp_path, monkeypatch, config, expected):
    monkeypatch.setattr(provenance, "capture_tool_versions", lambda: {})
    monkeypatch.setattr(provenance, "capture_gene_caller_identity", lambda version: {})
    provenance.write_provenance(tmp_path, config)
    assert (
        json.loads((tmp_path / "provenance.json").read_text())["configuration"]["marker_validation_top_k"] == expected
    )


@pytest.mark.parametrize("top_k", [7, 35])
def test_legacy_validation_passes_requested_depth_to_search(tmp_path, monkeypatch, top_k):
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">s_1 # 1 # 3 # 1 # ID=1_1\nM\n")
    (tmp_path / "ATPase.faa").write_text(">PLV__ref|p\nM\n")
    calls = []
    monkeypatch.setattr(hhg_seeding, "_run_diamond_marker_search", lambda **kw: calls.append(kw))
    hhg_seeding.validate_hmm_hits_with_markers(
        [hhg_seeding.HMMHit("s_1", "ATPase", 100, 1e-20, 100, 1, 1)],
        proteome,
        tmp_path,
        top_k=top_k,
    )
    assert calls[0]["max_seqs"] == top_k
