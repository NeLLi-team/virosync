"""Tests for the plain-Python orchestration backend."""

import csv
import io
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from virosync.orchestration import python_runner
from virosync.orchestration.runtime import call_task
from virosync.validation.tsv_invariants import (
    InvariantIssue,
    InvariantReport,
    TSVInvariantError,
)


class _TaskLike:
    def __init__(self):
        self.called = False
        self.fn_called = False

    def __call__(self, value: int) -> int:
        self.called = True
        return value + 100

    def fn(self, value: int) -> int:
        self.fn_called = True
        return value + 1


def test_call_task_uses_raw_function_for_task_like_objects() -> None:
    task = _TaskLike()

    assert call_task(task, 4) == 5
    assert task.fn_called is True
    assert task.called is False


def test_run_batch_python_writes_summary_and_preserves_input_order(
    tmp_path: Path,
    monkeypatch,
) -> None:
    inputs = []
    for name in ["b.fna", "a.fna"]:
        path = tmp_path / name
        path.write_text(">scaffold\nACGT\n")
        inputs.append(path)

    def _fake_single(
        *,
        genome_path,
        output_dir,
        genome_id,
        config,
    ):
        output_dir.mkdir(parents=True, exist_ok=True)
        return {
            "genome_id": genome_id,
            "success": True,
            "benchmark_eligible": True,
            "legacy_resume": False,
            "predictions": 2,
            "accepted": 1,
            "high_tier": 1,
            "medium_tier": 0,
            "low_tier": 0,
            "ncldv_count": 1,
            "mirus_count": 0,
            "ppv_count": 0,
            "cress_count": 0,
            "phage_count": 0,
            "viral_unknown_count": 0,
            "unknown_count": 0,
            "accepted_bp": 100,
            "total_genes": 5,
            "total_hallmarks": 2,
            "elapsed_sec": 3.2,
        }

    monkeypatch.setattr(python_runner, "_single_genome_callable", lambda: _fake_single)

    results = python_runner.run_batch_python(
        genome_paths=inputs,
        output_base_dir=tmp_path / "out",
        config=SimpleNamespace(),
        max_concurrent_genomes=2,
        retry_delay_seconds=0,
        effective_config={"schema_version": 1, "effective_config_sha256": "abc"},
    )

    assert [result["genome_id"] for result in results] == ["b", "a"]
    summary = (tmp_path / "out" / "batch_summary.tsv").read_text().splitlines()
    assert summary[0].startswith(
        "genome_id\tstatus\tbenchmark_eligible\tlegacy_resume\t"
        "predictions\taccepted"
    )
    assert summary[0].endswith(
        "\tncldv\tmirus\tppv\tcress\tphage\tviral_unknown\tunknown\t"
        "total_bp\tgenes\thallmarks\telapsed_sec\terror"
    )
    assert summary[1].startswith("b\tsuccess\ttrue\tfalse\t2\t1\t1")
    assert summary[2].startswith("a\tsuccess\ttrue\tfalse\t2\t1\t1")
    assert (tmp_path / "out" / "batch_report.md").exists()
    assert json.loads((tmp_path / "out" / "effective_config.json").read_text()) == {
        "effective_config_sha256": "abc",
        "schema_version": 1,
    }


def test_preflight_preserves_ordinary_genome_output_name(tmp_path: Path) -> None:
    input_path = tmp_path / "sample.fna"

    specs = python_runner._preflight_genome_runs([input_path], tmp_path / "out")

    assert specs == [
        python_runner.GenomeRunSpec(
            input_path=input_path,
            genome_id="sample",
            output_dir=(tmp_path / "out" / "sample").resolve(),
        )
    ]


def test_preflight_preserves_relative_output_path(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    input_path = Path("sample.fna")

    specs = python_runner._preflight_genome_runs([input_path], Path("results/run"))

    assert specs[0].output_dir == Path("results/run/sample")
    assert specs[0].output_dir.is_absolute() is False


def test_output_parent_segment_fails_before_executor_or_output_creation(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    input_path = Path("sample.fna")
    input_path.write_text(">scaffold\nACGT\n")
    sentinel = Path("sentinel.txt")
    sentinel.write_bytes(b"parent sentinel\n")
    executor_created = False

    def _unexpected_executor(*args, **kwargs):
        nonlocal executor_created
        executor_created = True
        raise AssertionError("executor must not be created for an invalid batch")

    monkeypatch.setattr(python_runner, "ThreadPoolExecutor", _unexpected_executor)

    with pytest.raises(ValueError, match="parent segment"):
        python_runner.run_batch_python(
            genome_paths=[input_path],
            output_base_dir=Path("results") / ".." / "out",
            config=SimpleNamespace(),
            max_concurrent_genomes=1,
            effective_config={"schema_version": 1},
        )

    assert executor_created is False
    assert Path("out").exists() is False
    assert sentinel.read_bytes() == b"parent sentinel\n"


@pytest.mark.parametrize("filename", ["...fna", "..fna"])
def test_unsafe_stem_fails_before_executor_or_output_creation(
    tmp_path: Path,
    monkeypatch,
    filename: str,
) -> None:
    input_path = tmp_path / filename
    input_path.write_text(">scaffold\nACGT\n")
    sentinel = tmp_path / "sentinel.txt"
    sentinel.write_bytes(b"parent sentinel\n")
    output_root = tmp_path / "out"
    executor_created = False

    def _unexpected_executor(*args, **kwargs):
        nonlocal executor_created
        executor_created = True
        raise AssertionError("executor must not be created for an invalid batch")

    monkeypatch.setattr(python_runner, "ThreadPoolExecutor", _unexpected_executor)

    with pytest.raises(ValueError, match="Unsafe or ambiguous"):
        python_runner.run_batch_python(
            genome_paths=[input_path],
            output_base_dir=output_root,
            config=SimpleNamespace(),
            max_concurrent_genomes=1,
        )

    assert executor_created is False
    assert output_root.exists() is False
    assert sentinel.read_bytes() == b"parent sentinel\n"


@pytest.mark.parametrize(
    "case",
    ["same_path", "different_extensions", "different_directories"],
)
def test_duplicate_genome_ids_list_every_source_and_create_no_output(
    tmp_path: Path,
    case: str,
) -> None:
    if case == "same_path":
        path = tmp_path / "sample.fna"
        paths = [path, path]
    elif case == "different_extensions":
        paths = [tmp_path / "sample.fa", tmp_path / "sample.fna"]
    else:
        paths = [tmp_path / "left" / "sample.fna", tmp_path / "right" / "sample.fna"]
    for path in set(paths):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(">scaffold\nACGT\n")
    output_root = tmp_path / "out"

    with pytest.raises(ValueError, match="duplicate genome ID") as exc_info:
        python_runner.run_batch_python(
            genome_paths=paths,
            output_base_dir=output_root,
            config=SimpleNamespace(),
            max_concurrent_genomes=2,
        )

    message = str(exc_info.value)
    for path in paths:
        assert str(path) in message
    assert output_root.exists() is False


def test_preexisting_output_symlink_is_rejected_before_worker(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "sample.fna"
    input_path.write_text(">scaffold\nACGT\n")
    output_root = tmp_path / "out"
    output_root.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    sentinel = outside / "sentinel.txt"
    sentinel.write_bytes(b"outside sentinel\n")
    (output_root / "sample").symlink_to(outside, target_is_directory=True)

    with pytest.raises(ValueError, match="Unsafe or ambiguous"):
        python_runner.run_batch_python(
            genome_paths=[input_path],
            output_base_dir=output_root,
            config=SimpleNamespace(),
            max_concurrent_genomes=1,
        )

    assert (output_root / "sample").is_symlink()
    assert sentinel.read_bytes() == b"outside sentinel\n"


def test_deterministic_invariant_error_is_not_retried(
    tmp_path: Path,
    monkeypatch,
) -> None:
    calls = 0
    invariant_error = TSVInvariantError(
        InvariantReport(
            rows_checked=1,
            issues=[InvariantIssue("EVE_1", "broken", "error", "injected")],
        ),
        tmp_path / "invariant.tsv",
    )

    def _fake_single(**kwargs):
        nonlocal calls
        calls += 1
        raise invariant_error

    def _unexpected_sleep(seconds):
        raise AssertionError("deterministic invariant errors must not sleep or retry")

    monkeypatch.setattr(python_runner, "_single_genome_callable", lambda: _fake_single)
    monkeypatch.setattr(python_runner.time, "sleep", _unexpected_sleep)
    spec = python_runner.GenomeRunSpec(
        input_path=tmp_path / "sample.fna",
        genome_id="sample",
        output_dir=tmp_path / "out" / "sample",
    )

    result = python_runner._run_one_genome(
        spec=spec,
        config=SimpleNamespace(),
        retries=3,
        retry_delay_seconds=60,
    )

    assert calls == 1
    assert result["success"] is False
    assert "Detailed TSV invariant check failed" in result["error"]


def test_clean_first_attempt_switches_to_schema_resume_on_retry(
    tmp_path: Path,
    monkeypatch,
) -> None:
    resume_values: list[bool] = []
    output_dir = tmp_path / "out" / "sample"

    def _fake_single(*, config, **kwargs):
        resume_values.append(config.execution.resume)
        if len(resume_values) == 1:
            output_dir.mkdir(parents=True)
            (output_dir / "phase1.complete.json").write_text("checkpoint\n")
            raise RuntimeError("transient failure after Phase 1")
        assert (output_dir / "phase1.complete.json").read_text() == "checkpoint\n"
        return {
            "genome_id": "sample",
            "success": True,
            "benchmark_eligible": True,
            "legacy_resume": False,
        }

    monkeypatch.setattr(python_runner, "_single_genome_callable", lambda: _fake_single)
    monkeypatch.setattr(python_runner.time, "sleep", lambda _seconds: None)
    config = python_runner.PipelineConfig()
    config.execution.resume = False
    spec = python_runner.GenomeRunSpec(
        input_path=tmp_path / "sample.fna",
        genome_id="sample",
        output_dir=output_dir,
    )

    result = python_runner._run_one_genome(
        spec=spec,
        config=config,
        retries=1,
        retry_delay_seconds=0,
    )

    assert result["success"] is True
    assert resume_values == [False, True]
    assert config.execution.resume is False


def test_non_mapping_worker_result_still_writes_one_failed_row(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_path = tmp_path / "sample.fna"
    input_path.write_text(">scaffold\nACGT\n")
    monkeypatch.setattr(
        python_runner,
        "_single_genome_callable",
        lambda: (lambda **kwargs: None),
    )

    results = python_runner.run_batch_python(
        genome_paths=[input_path],
        output_base_dir=tmp_path / "out",
        config=SimpleNamespace(),
        max_concurrent_genomes=1,
    )

    assert len(results) == 1
    assert results[0]["genome_id"] == "sample"
    assert results[0]["success"] is False
    with (tmp_path / "out" / "batch_summary.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["genome_id"] == "sample"
    assert rows[0]["status"] == "failed"
    assert "returned NoneType" in rows[0]["error"]


def test_batch_summary_quotes_multiline_tabbed_errors(tmp_path: Path) -> None:
    error = "first\tfield\nsecond line"

    summary_path = python_runner._write_batch_summary(
        tmp_path,
        [
            {
                "genome_id": "sample",
                "success": False,
                "error": error,
                "elapsed_sec": 0.0,
            }
        ],
    )

    with summary_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["error"] == error


def test_batch_outputs_fold_legacy_subtypes_and_include_cress(tmp_path: Path) -> None:
    legacy_result = {
        "genome_id": "synthetic",
        "success": True,
        "benchmark_eligible": True,
        "legacy_resume": False,
        "predictions": 8,
        "accepted": 8,
        "high_tier": 8,
        "medium_tier": 0,
        "low_tier": 0,
        "ncldv_count": 1,
        "vp_count": 1,
        "plv_count": 1,
        "mirus_count": 1,
        "mixed_count": 1,
        "ppv_count": 1,
        "cress_count": 1,
        "unknown_count": 1,
        "accepted_bp": 21007,
        "total_genes": 7,
        "total_hallmarks": 7,
        "elapsed_sec": 1,
    }
    result = python_runner._normalize_worker_result("synthetic", legacy_result)

    summary_path = python_runner._write_batch_summary(tmp_path, [result])
    report_path = python_runner._write_batch_report(tmp_path, [result])

    with summary_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    assert reader.fieldnames is not None
    assert "vp" not in reader.fieldnames
    assert "plv" not in reader.fieldnames
    assert "mixed" not in reader.fieldnames
    assert rows[0]["ppv"] == "3"
    assert rows[0]["viral_unknown"] == "1"
    assert rows[0]["cress"] == rows[0]["unknown"] == "1"

    report = report_path.read_text()
    assert "| PPV | 3 |" in report
    assert "| CRESS | 1 |" in report
    assert "| VIRAL_UNKNOWN | 1 |" in report
    assert "| UNKNOWN | 1 |" in report
    assert "| VP |" not in report
    assert "| PLV |" not in report
    assert "| MIXED |" not in report
    assert "| **Total** | **8** |" in report
    assert (
        "| synthetic | success | yes | no | 8 | 0 | 0 | 1 | 1 | 3 | 1 | 0 | 1 | 1 |"
        in report
    )


def test_batch_progress_is_monotonic_and_finishes_failed_queries() -> None:
    stream = io.StringIO()
    progress = python_runner.BatchProgress(
        2,
        stream=stream,
        is_tty=False,
        unit="genomes",
    )

    progress.update("first", 20, "phase 0")
    progress.update("first", 10, "stale update")
    progress.update("second", 40, "phase 1")
    progress.update("first", 50, "failed", failed=True)
    progress.update("second", 100, "complete")
    progress.finish(False)

    lines = [line for line in stream.getvalue().splitlines() if line]
    assert [line.split("] ", 1)[1].split("%", 1)[0].strip() for line in lines] == [
        "10",
        "30",
        "70",
        "100",
        "100",
    ]
    assert "2/2 genomes" in lines[-1]
    assert lines[-1].endswith("finished with failures")


def test_batch_outputs_fail_offsetting_per_genome_class_mismatches(
    tmp_path: Path,
) -> None:
    results = [
        {
            "genome_id": "over",
            "success": True,
            "accepted": 1,
            "ncldv_count": 2,
        },
        {
            "genome_id": "under",
            "success": True,
            "accepted": 1,
        },
        {
            "genome_id": "clean",
            "success": True,
            "benchmark_eligible": True,
            "legacy_resume": False,
            "predictions": 1,
            "accepted": 1,
            "high_tier": 1,
            "ncldv_count": 1,
            "accepted_bp": 400,
        },
    ]

    summary_path = python_runner._write_batch_summary(tmp_path, results)
    report_path = python_runner._write_batch_report(tmp_path, results)

    with summary_path.open(newline="") as handle:
        rows = {
            row["genome_id"]: row for row in csv.DictReader(handle, delimiter="\t")
        }
    assert set(rows) == {"over", "under", "clean"}
    for genome_id in ("over", "under"):
        assert rows[genome_id]["status"] == "failed"
        assert "do not sum to accepted predictions" in rows[genome_id]["error"]
        assert rows[genome_id]["accepted"] == "0"
        assert rows[genome_id]["ncldv"] == "0"
    assert rows["clean"]["status"] == "success"
    assert rows["clean"]["accepted"] == rows["clean"]["ncldv"] == "1"

    report = report_path.read_text()
    assert "| over | failed |" in report
    assert "| under | failed |" in report
    assert "for over: accepted=1 classified=2" in report
    assert "| **Total** | **1** |" in report
    assert "(1 successful, 2 failed)" in report
