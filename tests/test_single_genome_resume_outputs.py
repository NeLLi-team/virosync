from __future__ import annotations

from virosync.orchestration._flows.single_genome import (
    _completed_run_artifacts,
    _summarize_predictions_tsv,
)
from test_single_genome_resume import (
    _publish_schema3_success,
    _start_schema3_run,
)


def test_completed_run_artifacts_requires_final_root_outputs(tmp_path) -> None:
    output_dir = tmp_path / "genome"
    run_fingerprint = _start_schema3_run(output_dir)
    phase3_dir = output_dir / "phase3_synthesis"
    phase3_dir.mkdir(parents=True)
    (phase3_dir / "virosync_predictions.tsv").write_text("eve_id\n")

    assert _completed_run_artifacts(output_dir) is None

    (output_dir / "virosync_predictions_detailed.tsv").write_text("eve_id\n")
    assert _completed_run_artifacts(output_dir) is None

    (output_dir / "run.log").write_text(
        "# ViroSync Run Log: demo\n\n## Results Summary\nCanonical GEVEs: 0\n"
    )
    assert _completed_run_artifacts(output_dir) is None

    _publish_schema3_success(output_dir, run_fingerprint)
    artifacts = _completed_run_artifacts(output_dir)

    assert artifacts is not None
    assert set(artifacts) == {
        "phase3_predictions",
        "predictions_detailed",
        "run_log",
        "completion_manifest",
        "masking_status",
        "run_state",
    }


def test_completed_run_artifacts_rejects_mutated_prediction_headers(tmp_path) -> None:
    output_dir = tmp_path / "genome"
    run_fingerprint = _start_schema3_run(output_dir)
    phase3_dir = output_dir / "phase3_synthesis"
    phase3_dir.mkdir(parents=True)
    predictions = phase3_dir / "virosync_predictions.tsv"
    predictions.write_text("eve_id\n")
    (output_dir / "virosync_predictions_detailed.tsv").write_text("eve_id\n")
    (output_dir / "run.log").write_text(
        "# ViroSync Run Log: demo\n\n## Results Summary\nCanonical GEVEs: 0\n"
    )
    _publish_schema3_success(output_dir, run_fingerprint)
    predictions.write_text("not_eve_id\n")

    assert _completed_run_artifacts(output_dir) is None


def test_summarize_predictions_tsv_counts_canonical_rows_as_accepted(tmp_path) -> None:
    predictions_tsv = tmp_path / "virosync_predictions.tsv"
    predictions_tsv.write_text(
        "\t".join(
            [
                "eve_id",
                "length",
                "confidence_tier",
                "hallmark_total",
                "total_proteins",
                "classification",
            ]
        )
        + "\n"
        + "\t".join(["EVE_1", "100", "HIGH", "3", "5", "NCLDV"])
        + "\n"
        + "\t".join(["EVE_2", "50", "LOW", "1", "2", "VP"])
        + "\n"
        + "\t".join(["EVE_3", "80", "MEDIUM", "2", "4", "MIRUS"])
        + "\n"
    )

    stats = _summarize_predictions_tsv(predictions_tsv)

    assert stats["predictions"] == 3
    assert stats["accepted"] == 3
    assert stats["high_tier"] == 1
    assert stats["medium_tier"] == 1
    assert stats["low_tier"] == 1
    assert stats["accepted_bp"] == 230
    assert stats["total_genes"] == 11
    assert stats["total_hallmarks"] == 6
    assert stats["ncldv_count"] == 1
    # Legacy "VP" rows roll up under the unified Preplasmiviricota class.
    assert stats["ppv_count"] == 1
    assert stats["vp_count"] == 0
    assert stats["mirus_count"] == 1


def test_summarize_predictions_tsv_counts_detailed_candidates_without_accepting(tmp_path) -> None:
    detailed_tsv = tmp_path / "virosync_predictions_detailed.tsv"
    detailed_tsv.write_text(
        "eve_id\tlength\tconfidence_tier\thallmark_total\ttotal_proteins\tclassification\n"
        "EVE_1\t100\tHIGH\t3\t5\tNCLDV\n"
        "EVE_2\t50\tLOW\t1\t2\tVP\n"
    )

    stats = _summarize_predictions_tsv(detailed_tsv, canonical=False)

    assert stats["predictions"] == 2
    assert stats["accepted"] == 0
    assert stats["high_tier"] == 1
    assert stats["low_tier"] == 1
