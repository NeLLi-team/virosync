from pathlib import Path

from virosync.orchestration import tasks
from virosync.pipeline.phase1.marker_validation import NovelMarkerCriteria


def test_marker_validation_task_forwards_top_k_and_novel_criteria(
    tmp_path: Path,
    monkeypatch,
) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">p1\nM\n")
    marker_db = tmp_path / "marker.dmnd"
    marker_db.write_bytes(b"")
    seen: dict[str, object] = {}

    monkeypatch.setattr(
        "virosync.pipeline.phase1.extract_hmm_hit_sequences",
        lambda **kwargs: 1,
    )

    def fake_diamond(**kwargs):
        seen["diamond_max_seqs"] = kwargs["max_seqs"]
        kwargs["output_tsv"].write_text("")

    def fake_filter(**kwargs):
        seen["filter_max_seqs"] = kwargs["max_seqs"]
        seen["novel_criteria"] = kwargs["novel_criteria"]
        return []

    monkeypatch.setattr(
        "virosync.pipeline.phase1.run_diamond_on_hmm_hits",
        fake_diamond,
    )
    monkeypatch.setattr(
        "virosync.pipeline.phase1.filter_validated_markers",
        fake_filter,
    )
    criteria = NovelMarkerCriteria(
        min_hmm_score=41.0,
        min_hmm_coverage=0.65,
        require_cluster=False,
    )
    run_task = getattr(tasks.marker_validation_task, "fn", tasks.marker_validation_task)

    result = run_task(
        hmm_hits=[object()],
        proteome_path=proteome,
        marker_db=marker_db,
        output_dir=tmp_path / "out",
        threads=1,
        max_seqs=7,
        novel_criteria=criteria,
    )

    assert result == []
    assert seen["diamond_max_seqs"] == 7
    assert seen["filter_max_seqs"] == 7
    assert seen["novel_criteria"] is criteria
