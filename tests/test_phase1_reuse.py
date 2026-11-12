"""Start a run from the Phase 0 and Phase 1 results of an earlier run."""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import asdict
from pathlib import Path

import pytest
from click.testing import CliRunner
from test_phase1_resume_state import _phase1_kwargs
from test_run_state_orchestration import _MockPipeline, _NullResourceMonitor

from virosync.config import MaskingBackend, MaskingConfig, MaskingFailurePolicy, PipelineConfig
from virosync.orchestration import cli as orchestration_cli
from virosync.orchestration import python_runner
from virosync.orchestration._flows.single_genome import orchestrator, phase1_reuse
from virosync.orchestration._flows.single_genome import phase1 as phase1_module
from virosync.orchestration._flows.single_genome.phase1 import Phase1Result, _run_phase1_subflow
from virosync.orchestration._flows.single_genome.phase1_reuse import (
    Phase1Reuse,
    Phase1ReuseError,
    ReusedRun,
    check_rescue_host_search,
    import_phase0,
    import_phase1,
    rebuild_host_records,
)
from virosync.orchestration._flows.single_genome.phase1_state import load_phase1_state
from virosync.orchestration._flows.single_genome.run_state import (
    build_input_identity,
    canonical_sha256,
    compute_run_fingerprint,
    load_run_state,
    plan_resume,
)
from virosync.orchestration.cli import orchestrate
from virosync.pipeline import host_signatures
from virosync.pipeline.host_signatures import HostSignatureModel, TaxonomyLabelLookup
from virosync.pipeline.phase0 import masking as masking_module
from virosync.pipeline.phase0.masking import (
    MaskingBackendError,
    MaskingResult,
    load_masking_result,
    write_masking_status,
)
from virosync.pipeline.phase1.hhg_seeding import HMMHit
from virosync.pipeline.phase1.marker_validation import (
    NovelMarkerCriteria,
    collect_host_signatures,
    extract_hmm_hit_sequences,
    filter_validated_markers,
)

OLD_CODE = "3" * 64
NEW_CODE = "4" * 64


def _tree(root: Path) -> dict[str, str | None]:
    """Map every path below ``root`` to the SHA-256 of its bytes (None for a directory)."""
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        for path in sorted(root.rglob("*"))
    }


def _make_worlds(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    search_files: bool = True,
    terminal_phase: int | None = None,
    masking: MaskingConfig | None = None,
) -> tuple[_MockPipeline, _MockPipeline]:
    """Finish an earlier run, then build a second pipeline whose code identity differs."""
    (tmp_path / "old").mkdir()
    (tmp_path / "new").mkdir()
    old = _MockPipeline(tmp_path / "old", monkeypatch)
    old.search_files = search_files
    old.terminal_phase = terminal_phase
    if masking is not None:
        old.masking = masking
        old.identity["requested_masking"] = orchestrator._masking_request_identity(masking)
        old.fingerprint = compute_run_fingerprint(old.identity)
    assert old.run()["success"] is True
    new = _MockPipeline(tmp_path / "new", monkeypatch)
    new.identity["requested_masking"] = old.identity["requested_masking"]
    new.identity["code"]["source_sha256"] = NEW_CODE
    new.fingerprint = compute_run_fingerprint(new.identity)
    return old, new


@pytest.fixture
def worlds(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[_MockPipeline, _MockPipeline]:
    """An earlier finished run and a second pipeline whose code identity differs."""
    return _make_worlds(tmp_path, monkeypatch)


@pytest.fixture
def search_worlds(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[_MockPipeline, _MockPipeline]:
    """The same, with an ordinary marker in the old checkpoint and the Phase 1 search files on disk."""
    return _make_worlds(tmp_path, monkeypatch, search_files=True)


def _reuse(old: _MockPipeline, code_sha256: str = OLD_CODE) -> Phase1Reuse:
    return Phase1Reuse(old_base=old.output_dir.parent, code_sha256=code_sha256)


def _flow(new: _MockPipeline, reuse: Phase1Reuse, *, output_dir: Path | None = None, **kwargs) -> dict:
    return orchestrator.single_genome_flow(
        genome_path=new.genome,
        output_dir=output_dir or new.output_dir,
        genome_id="demo",
        device="cpu",
        phase1_reuse=reuse,
        **kwargs,
    )


def _sentinel_model() -> HostSignatureModel:
    return HostSignatureModel(
        token_weights={"sentinel": 3.0},
        token_counts={"sentinel": 1},
        token_bits={"sentinel": [3.0]},
        max_weight=3.0,
        min_token_length=3,
        host_prefixes=["EUK__", "MITO__", "PLASTID__"],
        weight_mode="rank",
    )


def test_reuse_runs_phase2_and_phase3_and_loads_neither_phase0_nor_phase1(
    worlds: tuple[_MockPipeline, _MockPipeline],
) -> None:
    old, new = worlds

    result = _flow(new, _reuse(old))

    assert result["success"] is True
    assert new.calls == {"phase0": 0, "phase1": 0, "phase2": 1, "phase3": 1}
    assert new.resume_flags["phase2"] == [False]
    assert new.resume_flags["phase3"] == [False]
    assert new.events == ["running", "phase0", "phase1", "phase2", "phase3", "success"]
    assert plan_resume(new.output_dir, expected_run_fingerprint=new.fingerprint).completed
    state = load_run_state(new.output_dir)
    assert state.run_fingerprint == new.fingerprint
    assert state.run_fingerprint != old.fingerprint
    assert state.identities["code"]["source_sha256"] == NEW_CODE


def test_reuse_leaves_the_old_folder_byte_identical(worlds: tuple[_MockPipeline, _MockPipeline]) -> None:
    old, new = worlds
    before = _tree(old.output_dir.parent)

    _flow(new, _reuse(old))

    assert _tree(old.output_dir.parent) == before


def test_reuse_copies_the_phase_files_and_records_their_provenance(
    worlds: tuple[_MockPipeline, _MockPipeline],
) -> None:
    old, new = worlds
    _flow(new, _reuse(old))

    copied = "phase1/marker_validation/host_diamond_top10.tsv"
    assert (new.output_dir / copied).read_bytes() == (old.output_dir / copied).read_bytes()
    provenance = json.loads((new.output_dir / "phase1/reuse_provenance.json").read_text())
    assert provenance["old_run_fingerprint"] == old.fingerprint
    assert provenance["old_code_sha256"] == OLD_CODE
    assert provenance["old_run_folder"] == str(old.output_dir.resolve())
    entry = next(item for item in provenance["copied"] if item["path"] == copied)
    assert entry["sha256"] == hashlib.sha256((new.output_dir / copied).read_bytes()).hexdigest()
    rewritten = {item["path"]: item for item in provenance["rewritten"]}
    checkpoint = "phase1/resume_state.json"
    assert rewritten[checkpoint]["old_sha256"] == hashlib.sha256((old.output_dir / checkpoint).read_bytes()).hexdigest()
    assert rewritten[checkpoint]["sha256"] == hashlib.sha256((new.output_dir / checkpoint).read_bytes()).hexdigest()
    phase0 = json.loads((new.output_dir / "phase0/reuse_provenance.json").read_text())
    assert phase0["old_run_fingerprint"] == old.fingerprint
    for phase in (0, 1):
        marker = json.loads((new.output_dir / f"phase{phase}.complete.json").read_text())
        assert f"phase{phase}/reuse_provenance.json" in {item["relative_path"] for item in marker["artifacts"]}


def test_host_model_is_rebuilt_by_the_current_builder_and_the_old_model_is_not_used(
    worlds: tuple[_MockPipeline, _MockPipeline], monkeypatch: pytest.MonkeyPatch
) -> None:
    old, new = worlds
    sentinel = _sentinel_model()
    builder_calls: list[tuple[list, dict]] = []

    def builder(records, **kwargs):
        builder_calls.append((list(records), kwargs))
        return sentinel

    monkeypatch.setattr(host_signatures, "build_host_signature_model", builder)
    seen: dict[str, object] = {}

    def phase2(**kwargs):
        seen["phase2"] = kwargs["host_signature_model"]
        return new._phase2(**kwargs)

    def phase3(**kwargs):
        seen["phase3"] = kwargs["host_signature_model_payload"]
        return new._phase3(**kwargs)

    monkeypatch.setattr(orchestrator, "_run_phase2_subflow", phase2)
    monkeypatch.setattr(orchestrator, "_run_phase3_subflow", phase3)

    _flow(new, _reuse(old))

    assert len(builder_calls) == 1
    assert builder_calls[0][1] == {
        "min_token_length": 3,
        "host_prefixes": {"EUK__", "MITO__", "PLASTID__"},
        "weight_mode": "rank",
    }
    assert load_phase1_state(old.output_dir / "phase1/resume_state.json").host_signature_model.token_weights == {}
    rebuilt = load_phase1_state(new.output_dir / "phase1/resume_state.json").host_signature_model
    assert rebuilt.token_weights == {"sentinel": 3.0}
    model_json = new.output_dir / "phase1/marker_validation/host_signature_model.json"
    assert json.loads(model_json.read_text()) == sentinel.to_dict()
    assert seen["phase2"] is sentinel
    assert seen["phase3"] == sentinel.to_dict()


def test_an_interrupted_reused_run_resumes_inside_the_new_folder(
    worlds: tuple[_MockPipeline, _MockPipeline],
) -> None:
    old, new = worlds
    new.fail_phase3_once = True
    with pytest.raises(RuntimeError, match="injected Phase 3 failure"):
        _flow(new, _reuse(old))

    result = _flow(new, _reuse(old))

    assert result["success"] is True
    # The retry loads the checkpoints of the new folder; it neither imports nor searches again.
    assert new.calls == {"phase0": 0, "phase1": 1, "phase2": 2, "phase3": 2}
    assert new.resume_flags["phase1"] == [True]
    assert new.resume_flags["phase2"] == [False, True]
    assert plan_resume(new.output_dir, expected_run_fingerprint=new.fingerprint).completed


def test_a_clean_run_imports_again(worlds: tuple[_MockPipeline, _MockPipeline]) -> None:
    old, new = worlds
    _flow(new, _reuse(old))

    result = _flow(new, _reuse(old), resume=False)

    assert result["success"] is True
    assert new.calls == {"phase0": 0, "phase1": 0, "phase2": 2, "phase3": 2}


def test_an_import_that_fails_half_way_is_repeated_cleanly(
    worlds: tuple[_MockPipeline, _MockPipeline], monkeypatch: pytest.MonkeyPatch
) -> None:
    old, new = worlds
    original = phase1_reuse._rebuild_host_model

    def fail_once(**kwargs):
        monkeypatch.setattr(phase1_reuse, "_rebuild_host_model", original)
        raise RuntimeError("injected import failure")

    monkeypatch.setattr(phase1_reuse, "_rebuild_host_model", fail_once)
    with pytest.raises(RuntimeError, match="injected import failure"):
        _flow(new, _reuse(old))
    assert load_run_state(new.output_dir).status == "failed"
    assert (new.output_dir / "phase0.complete.json").is_file()
    assert not (new.output_dir / "phase1.complete.json").exists()

    assert _flow(new, _reuse(old))["success"] is True
    assert plan_resume(new.output_dir, expected_run_fingerprint=new.fingerprint).completed


SEARCH_FILES = (
    "phase1/marker_validation/validated_marker_hits.tsv",
    "phase1/marker_validation/diamond_top10_taxonomy.tsv",
    "phase1/marker_validation/hmm_hit_porfs.faa",
    "phase1/hmm/hmm_hits_validated.tsv",
)


def test_an_old_run_with_its_search_files_is_reused_and_the_files_are_copied(
    search_worlds: tuple[_MockPipeline, _MockPipeline],
) -> None:
    old, new = search_worlds

    assert _flow(new, _reuse(old))["success"] is True

    provenance = json.loads((new.output_dir / "phase1/reuse_provenance.json").read_text())
    copied = {item["path"] for item in provenance["copied"]}
    assert set(SEARCH_FILES) <= copied
    for relative in SEARCH_FILES:
        assert (new.output_dir / relative).read_bytes() == (old.output_dir / relative).read_bytes()


@pytest.mark.parametrize("missing", SEARCH_FILES)
def test_an_old_run_without_a_file_that_phase3_reads_is_refused(
    search_worlds: tuple[_MockPipeline, _MockPipeline], missing: str
) -> None:
    old, new = search_worlds
    (old.output_dir / missing).unlink()

    with pytest.raises(Phase1ReuseError, match="valid Phase 0 and Phase 1"):
        _flow(new, _reuse(old))

    assert not new.output_dir.exists()


def test_an_old_run_that_ended_in_phase1_needs_no_search_files(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    old, new = _make_worlds(tmp_path, monkeypatch, search_files=False, terminal_phase=1)
    monkeypatch.setattr(phase1_module, "call_task", new._call_task)
    monkeypatch.setattr(phase1_module, "_generate_required_reports", new._generate_required_reports)

    assert _flow(new, _reuse(old))["success"] is True


def test_an_old_run_whose_masking_fell_back_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config = MaskingConfig(
        backend=MaskingBackend.TRF,
        failure_policy=MaskingFailurePolicy.FALLBACK,
        fallback_backend=MaskingBackend.OFF,
    )

    def failing_backend(**kwargs):
        raise MaskingBackendError(MaskingBackend.TRF, "injected failure")

    monkeypatch.setattr(masking_module, "_run_backend_attempt", failing_backend)
    old, new = _make_worlds(tmp_path, monkeypatch, masking=config)
    status = json.loads((old.output_dir / "phase0/masking/masking_status.json").read_text())
    assert status["status"] == "fallback"

    with pytest.raises(Phase1ReuseError, match=r"masking fell back.*trf.*off"):
        _flow(new, _reuse(old), masking=config)

    assert not new.output_dir.exists()


def test_a_moved_input_is_refused_with_its_recorded_path(worlds: tuple[_MockPipeline, _MockPipeline]) -> None:
    old, new = worlds
    recorded = str(old.genome.resolve())
    old.genome.rename(old.genome.with_name("moved.fna"))

    with pytest.raises(Phase1ReuseError, match="recorded") as refusal:
        _flow(new, _reuse(old))

    assert recorded in str(refusal.value)


def test_the_option_help_and_the_docs_state_the_code_precondition() -> None:
    result = CliRunner().invoke(orchestrate, ["run", "--help"])
    help_text = " ".join(result.output.split())
    assert "identical to the current source, except for the host signature model" in help_text
    assert "cannot check" in help_text
    root = Path(__file__).parents[1]
    for document in ("docs/METHODS.md", "docs/reference/cli.md"):
        assert "except for the host signature model" in " ".join((root / document).read_text().split())


def test_old_run_that_ended_in_phase1_gives_the_same_terminal_result(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    (tmp_path / "old").mkdir()
    (tmp_path / "new").mkdir()
    old = _MockPipeline(tmp_path / "old", monkeypatch)
    old.terminal_phase = 1
    assert old.run()["success"] is True
    new = _MockPipeline(tmp_path / "new", monkeypatch)
    new.identity["code"]["source_sha256"] = NEW_CODE
    new.fingerprint = compute_run_fingerprint(new.identity)
    monkeypatch.setattr(phase1_module, "call_task", new._call_task)
    monkeypatch.setattr(phase1_module, "_generate_required_reports", new._generate_required_reports)
    stale_model = old.output_dir / "phase1/marker_validation/host_signature_model.json"
    stale_model.write_text('{"stale": true}\n')
    # Publish this test's synthetic source model as part of its authenticated fixture.
    marker_path = old.output_dir / "phase1.complete.json"
    marker = json.loads(marker_path.read_text())
    artifacts = {item["relative_path"]: item for item in marker["artifacts"]}
    artifacts.update({item.relative_path: asdict(item) for item in orchestrator._phase_artifacts(old.output_dir, 1)})
    marker["artifacts"] = [artifacts[key] for key in sorted(artifacts)]
    marker_path.write_text(json.dumps(marker))
    sentinel = _sentinel_model()
    monkeypatch.setattr(host_signatures, "build_host_signature_model", lambda records, **kwargs: sentinel)
    before = _tree(old.output_dir.parent)

    result = _flow(new, _reuse(old))

    assert result["success"] is True
    assert new.calls == {"phase0": 0, "phase1": 0, "phase2": 0, "phase3": 0}
    new_state = load_run_state(new.output_dir)
    old_state = load_run_state(old.output_dir)
    assert new_state.status == "success"
    assert new_state.result == old_state.result
    assert new_state.result["terminal_phase"] == 1
    assert plan_resume(new.output_dir, expected_run_fingerprint=new.fingerprint).terminal_phase == 1
    for name in ("virosync_predictions.tsv", "virosync_predictions_detailed.tsv", "virosync_predictions.bed"):
        assert (new.output_dir / name).read_bytes() == (old.output_dir / name).read_bytes()
    manifest = json.loads((new.output_dir / "virosync_run_complete.json").read_text())
    assert manifest["reason"] == json.loads((old.output_dir / "virosync_run_complete.json").read_text())["reason"]
    assert manifest["config_fingerprint"] == new.fingerprint
    assert json.loads((new.output_dir / "phase1/marker_validation/host_signature_model.json").read_text()) == (
        sentinel.to_dict()
    )
    assert _tree(old.output_dir.parent) == before


def test_old_terminal_run_without_a_host_model_gets_none(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "old").mkdir()
    (tmp_path / "new").mkdir()
    old = _MockPipeline(tmp_path / "old", monkeypatch)
    old.terminal_phase = 1
    old.run()
    new = _MockPipeline(tmp_path / "new", monkeypatch)
    new.identity["code"]["source_sha256"] = NEW_CODE
    new.fingerprint = compute_run_fingerprint(new.identity)
    monkeypatch.setattr(phase1_module, "call_task", new._call_task)
    monkeypatch.setattr(phase1_module, "_generate_required_reports", new._generate_required_reports)

    assert _flow(new, _reuse(old))["success"] is True

    assert not (new.output_dir / "phase1/marker_validation/host_signature_model.json").exists()


def test_wrong_code_sha256_is_refused_before_anything_is_written(
    worlds: tuple[_MockPipeline, _MockPipeline],
) -> None:
    old, new = worlds

    with pytest.raises(Phase1ReuseError, match=r"code SHA-256.*" + OLD_CODE):
        _flow(new, _reuse(old, "5" * 64))

    assert not new.output_dir.exists()
    assert new.calls == {"phase0": 0, "phase1": 0, "phase2": 0, "phase3": 0}


def _change_input(new: _MockPipeline) -> None:
    new.genome.write_text(">scaffold\nACGTA\n")
    new.identity["input"] = asdict(build_input_identity(new.genome))


def _change_config(new: _MockPipeline) -> None:
    new.identity["config"]["sha256"] = "9" * 64


def _change_resources(new: _MockPipeline) -> None:
    new.identity["resources"] = [{"name": "hmm", "kind": "core", "version": "1", "manifest_sha256": "a" * 64}]


def _change_environment(new: _MockPipeline) -> None:
    environment = new.identity["environment"]
    environment["runtime_sha256"] = "6" * 64
    environment["sha256"] = canonical_sha256({key: value for key, value in environment.items() if key != "sha256"})


def _change_masking(new: _MockPipeline) -> None:
    new.identity["requested_masking"] = {**new.identity["requested_masking"], "backend": "trf"}


@pytest.mark.parametrize(
    ("component", "change"),
    [
        ("input", _change_input),
        ("config", _change_config),
        ("resources", _change_resources),
        ("environment", _change_environment),
        ("requested_masking", _change_masking),
    ],
)
def test_a_differing_identity_component_is_refused_and_named(
    worlds: tuple[_MockPipeline, _MockPipeline], component: str, change
) -> None:
    old, new = worlds
    change(new)
    new.fingerprint = compute_run_fingerprint(new.identity)

    with pytest.raises(Phase1ReuseError, match=rf"identity differs.*{component}"):
        _flow(new, _reuse(old))

    assert not new.output_dir.exists()


def test_only_the_production_arm_can_be_reused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "old").mkdir()
    (tmp_path / "new").mkdir()
    old = _MockPipeline(tmp_path / "old", monkeypatch)
    old.identity["config"]["ablation_id"] = "A2"
    old.fingerprint = compute_run_fingerprint(old.identity)
    orchestrator.single_genome_flow(
        genome_path=old.genome, output_dir=old.output_dir, genome_id="demo", device="cpu", ablation_id="A2"
    )
    new = _MockPipeline(tmp_path / "new", monkeypatch)
    new.identity["config"]["ablation_id"] = "A2"
    new.identity["code"]["source_sha256"] = NEW_CODE
    new.fingerprint = compute_run_fingerprint(new.identity)

    with pytest.raises(Phase1ReuseError, match="A0"):
        _flow(new, _reuse(old), ablation_id="A2")


@pytest.mark.parametrize("artifact", ["phase0/proteome.fasta", "phase1/resume_state.json"])
def test_a_tampered_old_artifact_is_refused(worlds: tuple[_MockPipeline, _MockPipeline], artifact: str) -> None:
    old, new = worlds
    path = old.output_dir / artifact
    path.write_bytes(path.read_bytes() + b" ")

    with pytest.raises(Phase1ReuseError, match="valid Phase 0 and Phase 1"):
        _flow(new, _reuse(old))

    assert not new.output_dir.exists()


def test_a_missing_old_run_is_refused(worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path) -> None:
    old, new = worlds
    empty = tmp_path / "empty"
    empty.mkdir()

    with pytest.raises(Phase1ReuseError, match="no earlier run"):
        _flow(new, Phase1Reuse(old_base=empty, code_sha256=OLD_CODE))


def test_the_genome_folder_is_not_accepted_for_the_base_folder(
    worlds: tuple[_MockPipeline, _MockPipeline],
) -> None:
    old, new = worlds

    with pytest.raises(Phase1ReuseError, match="base folder"):
        _flow(new, Phase1Reuse(old_base=old.output_dir, code_sha256=OLD_CODE))


@pytest.mark.parametrize("layout", ["same", "inside", "contains"])
@pytest.mark.parametrize("resume", [True, False], ids=["resume", "clean-run"])
def test_overlapping_folders_are_refused_and_clean_run_cannot_delete_the_old_run(
    worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path, layout: str, resume: bool
) -> None:
    old, new = worlds
    old_base = old.output_dir.parent
    output_dir = {
        "same": old.output_dir,
        "inside": old_base / "nested" / "demo",
        "contains": old_base.parent / "demo",
    }[layout]
    before = _tree(tmp_path / "old")

    with pytest.raises(Phase1ReuseError, match="inside|contain|equal"):
        _flow(new, _reuse(old), output_dir=output_dir, resume=resume)

    assert _tree(tmp_path / "old") == before


def test_old_folder_for_several_genomes_is_the_base_folder(tmp_path: Path) -> None:
    reuse = Phase1Reuse(old_base=tmp_path / "campaign", code_sha256=OLD_CODE)

    assert reuse.old_dir("g1") == tmp_path / "campaign" / "g1"
    assert reuse.old_dir("g2") == tmp_path / "campaign" / "g2"


def test_the_batch_preflight_refuses_overlapping_bases_before_any_write(tmp_path: Path) -> None:
    old_base = tmp_path / "campaign"
    old_base.mkdir()
    (old_base / "g1.txt").write_text("keep")
    reuse = Phase1Reuse(old_base=old_base, code_sha256=OLD_CODE)
    before = _tree(tmp_path)

    with pytest.raises(ValueError, match="equal"):
        python_runner.run_batch_python(
            genome_paths=[tmp_path / "g1.fna"],
            output_base_dir=old_base,
            config=PipelineConfig(),
            max_concurrent_genomes=1,
            effective_config={"threads": 1},
            phase1_reuse=reuse,
        )

    assert _tree(tmp_path) == before


def test_the_runner_hands_the_request_to_every_flow_call_and_does_not_retry_a_refusal(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reuse = Phase1Reuse(old_base=tmp_path / "campaign", code_sha256=OLD_CODE)
    calls: list[dict] = []

    def refuse(**kwargs):
        calls.append(kwargs)
        raise Phase1ReuseError("not reusable")

    monkeypatch.setattr(python_runner, "single_genome_flow", refuse)
    spec = python_runner.GenomeRunSpec(
        input_path=tmp_path / "g1.fna", genome_id="g1", output_dir=tmp_path / "results" / "g1"
    )

    result = python_runner._run_one_genome(
        spec=spec, config=PipelineConfig(), retries=2, retry_delay_seconds=0, phase1_reuse=reuse
    )

    assert result["success"] is False
    assert "not reusable" in result["error"]
    assert len(calls) == 1
    assert calls[0]["phase1_reuse"] is reuse


def _invoke_run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *options: str, output: Path | None = None):
    received: dict[str, object] = {}
    genome = tmp_path / "input.fna"
    genome.write_text(">scaffold_1\nACGT\n", encoding="utf-8")
    monkeypatch.setattr(orchestration_cli, "_resolve_pipeline_resources", lambda config, orchestration, path: config)
    monkeypatch.setattr(orchestration_cli, "_validate_runtime_config", lambda config: None)
    monkeypatch.setattr(orchestration_cli, "_preflight_native_tools", lambda: None)
    monkeypatch.setattr(
        orchestration_cli,
        "_resolve_optional_features",
        lambda config: (config, {}),
    )
    monkeypatch.setattr(orchestration_cli.ApplicationConfig, "effective_payload", lambda self, features: {})

    def fake_runner(**kwargs: object) -> list[dict]:
        received.update(kwargs)
        return [
            {
                "success": True,
                "benchmark_eligible": True,
                "legacy_resume": False,
                "accepted": 0,
                "predictions": 0,
                "elapsed_sec": 0.0,
                "genome_id": "input",
            }
        ]

    monkeypatch.setattr(orchestration_cli, "run_batch_python", fake_runner)
    result = CliRunner().invoke(
        orchestrate, ["run", "-i", str(genome), "-o", str(output or tmp_path / "results"), *options]
    )
    return result, received


def test_cli_passes_the_base_folder_and_the_code_sha256(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    old_base = tmp_path / "campaign"
    old_base.mkdir()

    result, received = _invoke_run(
        tmp_path,
        monkeypatch,
        "--reuse-phase1-from",
        str(old_base),
        "--reuse-phase1-code-sha256",
        OLD_CODE,
    )

    assert result.exit_code == 0, result.output
    assert received["phase1_reuse"] == Phase1Reuse(old_base=old_base, code_sha256=OLD_CODE)


def test_cli_without_the_options_passes_no_request(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    result, received = _invoke_run(tmp_path, monkeypatch)

    assert result.exit_code == 0, result.output
    assert received["phase1_reuse"] is None


@pytest.mark.parametrize(
    "options",
    [
        ("--reuse-phase1-code-sha256", OLD_CODE),
        ("--reuse-phase1-from", "{old}"),
    ],
    ids=["sha-only", "folder-only"],
)
def test_cli_needs_both_reuse_options(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, options: tuple) -> None:
    old_base = tmp_path / "campaign"
    old_base.mkdir()

    result, received = _invoke_run(tmp_path, monkeypatch, *(item.format(old=old_base) for item in options))

    assert result.exit_code != 0
    assert "required together" in result.output
    assert not received


def test_cli_rejects_a_malformed_code_sha256(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    old_base = tmp_path / "campaign"
    old_base.mkdir()

    result, received = _invoke_run(
        tmp_path, monkeypatch, "--reuse-phase1-from", str(old_base), "--reuse-phase1-code-sha256", "ABC"
    )

    assert result.exit_code != 0
    assert "--reuse-phase1-code-sha256: the code SHA-256" in result.output
    assert not received


def test_cli_refuses_an_output_folder_that_equals_the_old_folder(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    old_base = tmp_path / "campaign"
    old_base.mkdir()
    (old_base / "keep.txt").write_text("keep")
    before = _tree(old_base)

    result, received = _invoke_run(
        tmp_path,
        monkeypatch,
        "--reuse-phase1-from",
        str(old_base),
        "--reuse-phase1-code-sha256",
        OLD_CODE,
        output=old_base,
    )

    assert result.exit_code != 0
    assert "equal" in result.output
    assert not received
    assert _tree(old_base) == before


# --- Host records -------------------------------------------------------------------


_PROTEOME = """\
>s1_1 # 11 # 100 # 1 # ID=1_1;partial=00
MKTAYIAKQR
>s1_2 # 201 # 300 # -1 # ID=1_2;partial=00
MKVLAAGIVG
>s2_1 # 1 # 90 # 1 # ID=2_1;partial=00
MQRSTLLAAG
>s2_2 # 5 # 50 # 1 # ID=2_2;partial=00
MMMMMKKKKK
>s2_3 # 60 # 120 # 1 # ID=2_3;partial=00
MGGGGAAAAA
"""
_GENOME = ">s1\n" + "ACGT" * 100 + "\n>s2\n" + "ACGT" * 50 + "\n"
# (query, target, evalue, bits, pident, qcov); rows of different queries are interleaved.
_HOST_ROWS = [
    ("s1_1", "EUK__A|p5", 1e-30, 300.0, 40.0, 90.0),
    ("s2_1", "EUK__A|p1", 1e-20, 200.0, 35.0, 80.0),
    ("s1_2", "EUK__B|p9", 1e-18, 180.0, 33.0, 80.0),
    ("s2_3", "NCLDV__X|p1", 1e-15, 200.0, 45.0, 80.0),
    ("s2_1", "BAC__C|p3", 1e-12, 150.0, 30.0, 70.0),
    ("s1_1", "EUK__A|p6", 1e-25, 250.0, 38.0, 90.0),
    ("s1_2", "EUK__A|p8", 1e-10, 120.0, 30.0, 70.0),
    ("s2_3", "EUK__A|p2", 1e-8, 100.0, 28.0, 60.0),
    ("s2_1", "EUK__B|p2", 1e-6, 100.0, 27.0, 60.0),
    ("s1_1", "EUK__B|p7", 1e-5, 90.0, 26.0, 50.0),
]
# Hit order matters: s2_1 comes first and s1_2 has two segments around s1_1.
_HITS = [
    HMMHit("s2_1", "GVOGm0004", 70.0, 1e-12, 70.0, 1, 10),
    HMMHit("s1_2", "GVOGm0001", 50.0, 1e-8, 50.0, 1, 5),
    HMMHit("s1_1", "GVOGm0002", 60.0, 1e-9, 60.0, 1, 10),
    HMMHit("s1_2", "GVOGm0010", 55.0, 1e-9, 55.0, 6, 10),
    HMMHit("s2_2", "GVOGm0005", 40.0, 1e-5, 40.0, 1, 10),
    HMMHit("s2_3", "GVOGm0006", 45.0, 1e-6, 45.0, 1, 10),
]


@pytest.fixture
def host_search(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    lookup = TaxonomyLabelLookup(
        {
            "EUK__A": "EUK|Amoebozoa|Entamoeba",
            "EUK__B": "EUK|Alveolata|Ciliophora",
            "BAC__C": "BAC|Proteobacteria|Escherichia",
            "NCLDV__X": "NCLDV|Mimiviridae|Mimivirus",
        }
    )
    monkeypatch.setattr(host_signatures, "_TAXONOMY_LOOKUP", lookup)
    paths = {
        "proteome": tmp_path / "proteome.fasta",
        "genome": tmp_path / "genome.fna",
        "tsv": tmp_path / "host_diamond_top10.tsv",
        "faa": tmp_path / "host_hmm_hit_porfs.faa",
    }
    paths["proteome"].write_text(_PROTEOME)
    paths["genome"].write_text(_GENOME)
    paths["tsv"].write_text("".join("\t".join(map(str, row)) + "\n" for row in _HOST_ROWS))
    extract_hmm_hit_sequences(_HITS, paths["proteome"], paths["faa"])
    return paths


@pytest.mark.parametrize("weight_mode", ["rank", "bitscore"])
def test_rebuilt_host_records_give_the_model_of_the_search_records(
    host_search: dict[str, Path], weight_mode: str
) -> None:
    lookup = host_signatures.get_taxonomy_lookup()
    criteria = NovelMarkerCriteria()
    searched = filter_validated_markers(
        hmm_hits=_HITS,
        diamond_output_tsv=host_search["tsv"],
        proteome_fasta=host_search["proteome"],
        genome_fasta=host_search["genome"],
        output_dir=None,
        taxonomy_lookup=lookup,
        taxonomy_weight_mode=weight_mode,
        novel_criteria=criteria,
        max_seqs=10,
    )

    rebuilt = rebuild_host_records(
        host_query_fasta=host_search["faa"],
        host_diamond_tsv=host_search["tsv"],
        proteome_path=host_search["proteome"],
        genome_path=host_search["genome"],
        taxonomy_lookup=lookup,
        taxonomy_weight_mode=weight_mode,
        novel_criteria=criteria,
    )

    assert [record.query_porf.split("|aa")[0] for record in rebuilt] == ["s2_1", "s1_2", "s1_1", "s2_2", "s2_3"]
    prefixes = {"EUK__"}
    model_args = {"min_token_length": 3, "host_prefixes": prefixes, "weight_mode": weight_mode}
    expected = host_signatures.build_host_signature_model(searched, **model_args)
    actual = host_signatures.build_host_signature_model(rebuilt, **model_args)
    assert expected.token_weights  # the fixture reaches the model
    assert json.dumps(actual.to_dict(), indent=2) == json.dumps(expected.to_dict(), indent=2)
    assert collect_host_signatures(rebuilt, host_prefixes=prefixes) == collect_host_signatures(
        searched, host_prefixes=prefixes
    )


def test_a_missing_host_search_file_stops_the_rebuild(host_search: dict[str, Path]) -> None:
    host_search["tsv"].unlink()

    with pytest.raises(Phase1ReuseError, match="host_diamond_top10.tsv"):
        rebuild_host_records(
            host_query_fasta=host_search["faa"],
            host_diamond_tsv=host_search["tsv"],
            proteome_path=host_search["proteome"],
            genome_path=host_search["genome"],
            taxonomy_lookup=None,
            taxonomy_weight_mode="rank",
            novel_criteria=NovelMarkerCriteria(),
        )


def test_a_confirmed_rescue_protein_in_the_rescue_host_search_stops_the_import(tmp_path: Path) -> None:
    rescue_host_faa = tmp_path / "host_hmm_hit_porfs.faa"
    rescue_host_faa.write_text(">GVOGm0001_VSR0123456789abcdef\nMKV\n>other\nMKV\n")

    check_rescue_host_search(rescue_host_faa, confirmed_ids={"unrelated_VSR0000000000000000"})
    with pytest.raises(Phase1ReuseError, match="GVOGm0001_VSR0123456789abcdef"):
        check_rescue_host_search(rescue_host_faa, confirmed_ids={"GVOGm0001_VSR0123456789abcdef"})


# --- Equivalence with a fresh Phase 1 search ------------------------------------------

# Marker-search matches; s1_1 has viral support, so the host search excludes it.
_MAIN_ROWS = [
    ("s1_1", "NCLDV__X|p1", 1e-40, 400.0, 45.0, 90.0),
    ("s1_1", "EUK__A|p5", 1e-30, 300.0, 40.0, 90.0),
    ("s1_2", "EUK__B|p9", 1e-18, 180.0, 33.0, 80.0),
    ("s2_1", "EUK__A|p1", 1e-20, 200.0, 35.0, 80.0),
    ("s2_3", "EUK__A|p2", 1e-8, 100.0, 28.0, 60.0),
]


def _fasta_ids(path: Path) -> list[str]:
    return [line[1:].split()[0] for line in path.read_text().splitlines() if line.startswith(">")]


def _tsv(rows: list[tuple], ids: list[str]) -> str:
    return "".join("\t".join(map(str, row)) + "\n" for row in rows if row[0] in ids)


def _fresh_phase1(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, assemble_regions: bool
) -> tuple[dict, Phase1Result]:
    """Run the real marker task with canned DIAMOND output; return the Phase 1 keyword arguments and result."""
    from virosync.orchestration import resource_monitor
    from virosync.orchestration.tasks import marker_validation_task, region_assembly_task
    from virosync.pipeline import phase1 as phase1_package
    from virosync.pipeline.phase1 import marker_validation as marker_validation_module

    labels = tmp_path / "labels.tsv"
    labels.write_text(
        "EUK__A\tEUK|Amoebozoa|Entamoeba\nEUK__B\tEUK|Alveolata|Ciliophora\n"
        "BAC__C\tBAC|Proteobacteria|Escherichia\nNCLDV__X\tNCLDV|Mimiviridae|Mimivirus\n"
    )
    kwargs = _phase1_kwargs(tmp_path, tmp_path / "fresh")
    config = kwargs["config"]
    config.execution.resume = False
    config.phase1.pfam_arbitration_enabled = False
    config.phase1.frameshift_screening_enabled = False
    config.databases.hmm_database = tmp_path / "markers.hmm"
    config.databases.marker_db = tmp_path / "markers.dmnd"
    config.databases.taxonomy_labels_file = labels
    config.databases.hmm_database.write_text("HMM\n")
    config.databases.marker_db.write_text("DB\n")
    kwargs["resume_authorized"] = False
    kwargs["masked_path"].write_text(_GENOME)
    kwargs["proteome_path"].write_text(_PROTEOME)

    def marker_search(*, hmm_hit_fasta, output_tsv, **_ignored):
        output_tsv.write_text(_tsv(_MAIN_ROWS, _fasta_ids(hmm_hit_fasta)))
        return output_tsv

    def host_search(hmm_hit_fasta, _database, output_tsv, **_ignored):
        output_tsv.write_text(_tsv(_HOST_ROWS, _fasta_ids(hmm_hit_fasta)))
        return output_tsv

    def run_task(task, **kw):
        if task is phase1_module.hhg_seeding_task:
            return [], list(_HITS)
        if task is marker_validation_task:
            return marker_validation_task(**kw)
        if task is region_assembly_task:
            return region_assembly_task(**kw) if assemble_regions else []
        raise AssertionError(task)

    monkeypatch.setattr(resource_monitor, "ResourceMonitor", _NullResourceMonitor)
    monkeypatch.setattr(phase1_package, "run_diamond_on_hmm_hits", marker_search)
    monkeypatch.setattr(marker_validation_module, "run_host_diamond_on_hmm_hits", host_search)
    monkeypatch.setattr(phase1_module, "call_task", run_task)
    fresh = _run_phase1_subflow(**kwargs, metadata_families={})
    assert isinstance(fresh, Phase1Result)
    return kwargs, fresh


def _import_into(replay_dir: Path, kwargs: dict) -> Phase1Result:
    source = ReusedRun(
        root=kwargs["output_dir"],
        run_fingerprint="1" * 64,
        code_sha256=OLD_CODE,
        markers={0: "a" * 64, 1: "b" * 64},
        checkpoint=load_phase1_state(kwargs["output_dir"] / "phase1/resume_state.json"),
        artifacts={name: digest for name, digest in _tree(kwargs["output_dir"]).items() if digest is not None},
    )
    replay = import_phase1(
        source,
        output_dir=replay_dir,
        genome_id=kwargs["genome_id"],
        config=kwargs["config"],
        masked_path=kwargs["masked_path"],
        proteome_path=kwargs["proteome_path"],
        config_fingerprint=kwargs["config_fingerprint"],
        logger=logging.getLogger("test-phase1-reuse"),
    )
    assert isinstance(replay, Phase1Result)
    return replay


def test_import_gives_the_phase1_files_of_a_fresh_search(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The real marker task writes the host search files; the import must rebuild the same model from them."""
    kwargs, fresh = _fresh_phase1(tmp_path, monkeypatch, assemble_regions=False)
    fresh_dir, replay_dir = kwargs["output_dir"], tmp_path / "replay"
    host_proteins = _fasta_ids(fresh_dir / "phase1/marker_validation/host_hmm_hit_porfs.faa")
    assert host_proteins == ["s2_1", "s1_2", "s2_2", "s2_3"]
    assert fresh.host_signature_model.token_weights  # the fixture reaches the model

    replay = _import_into(replay_dir, kwargs)

    assert replay.host_signature_model == fresh.host_signature_model
    assert replay.host_signatures == fresh.host_signatures
    assert replay.validated_markers == fresh.validated_markers
    for relative in (
        "phase1/marker_validation/host_signature_model.json",
        "phase1/marker_validation/host_signature_model_debug.tsv",
        "phase1/marker_validation/validated_marker_hits.tsv",
        "phase1/marker_validation/host_diamond_top10.tsv",
        "phase1/resume_state.json",
    ):
        assert (replay_dir / relative).read_bytes() == (fresh_dir / relative).read_bytes(), relative


def test_an_assembled_region_of_the_imported_checkpoint_runs_through_phase2_like_the_fresh_one(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from virosync.orchestration._flows.single_genome import phase2
    from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy, GenomeDiamondQuery, SeedGeneMapping, pORF

    kwargs, fresh = _fresh_phase1(tmp_path, monkeypatch, assemble_regions=True)
    replay = _import_into(tmp_path / "replay", kwargs)
    assert fresh.merged_seeds, "the fixture must assemble one region"
    assert replay.merged_seeds == fresh.merged_seeds
    seed = fresh.merged_seeds[0]
    genes = [
        pORF(id=f"{seed.scaffold}_{index}", scaffold=seed.scaffold, start=100 * index, end=100 * index + 100)
        for index in range(4)
    ]
    taxonomy = {
        gene.id: GeneTaxonomy(
            porf_id=gene.id,
            scaffold=gene.scaffold,
            start=gene.start,
            end=gene.end,
            has_hit=True,
            top1_prefix="EUK__",
            top1_target="EUK__A|p",
            top10_prefixes=["EUK__"],
            top10_pidents=[80.0],
            top10_targets=["EUK__A|p"],
            top10_bits=[100.0],
        )
        for gene in genes
    }
    query = GenomeDiamondQuery(
        seed_gene_mappings={
            seed.seed_id: SeedGeneMapping(
                seed_id=seed.seed_id,
                scaffold=seed.scaffold,
                seed_start=seed.start,
                seed_end=seed.end,
                flank_start_bp=0,
                flank_end_bp=400,
            )
        }
    )
    monkeypatch.setattr(phase2, "build_proteome_index", lambda _path: {seed.scaffold: genes})
    monkeypatch.setattr(phase2, "_require_phase2b_gene_taxonomy_db", lambda *_args, **_kwargs: tmp_path / "db")
    monkeypatch.setattr(
        phase2,
        "_run_boundary_diamond",
        lambda **_kwargs: phase2._BoundaryEvidence(taxonomy_map=taxonomy, control_stats=None, query=query),
    )
    config = PipelineConfig().with_overrides(resume=False, threads=1, gene_taxonomy_threads=1)

    def run_phase2(result: Phase1Result, output_dir: Path) -> phase2.Phase2Result:
        outcome = phase2._run_phase2_subflow(
            masked_path=kwargs["masked_path"],
            proteome_path=kwargs["proteome_path"],
            merged_seeds=result.merged_seeds,
            validated_markers=result.validated_markers,
            host_signature_model=result.host_signature_model,
            output_dir=output_dir,
            genome_id="genome-a",
            config=config,
            genome_start_time=0,
            logger=logging.getLogger("test-phase1-reuse"),
            metadata_families={},
        )
        assert isinstance(outcome, phase2.Phase2Result)
        return outcome

    fresh_boundaries = run_phase2(fresh, tmp_path / "phase2-fresh").refined_boundaries
    replay_boundaries = run_phase2(replay, tmp_path / "phase2-replay").refined_boundaries

    assert fresh_boundaries
    assert replay_boundaries == fresh_boundaries
    assert replay_boundaries[0].confirmed_seeds == seed.confirmed_seeds


# --- Phase 0 ------------------------------------------------------------------------


def test_phase0_import_rebases_the_masked_genome_and_keeps_the_old_status(tmp_path: Path) -> None:
    old_root = tmp_path / "old" / "g"
    new_root = tmp_path / "new" / "g"
    genome = tmp_path / "g.fna"
    genome.write_text(">s\nACGTACGT\n")
    masking_dir = old_root / "phase0" / "masking"
    masking_dir.mkdir(parents=True)
    masked = masking_dir / "genome.masked.fna"
    masked.write_text(">s\nACGTNNNN\n")
    result = MaskingResult(
        output_path=masked,
        repeat_regions=(),
        requested_backend=MaskingBackend.TRF,
        effective_backend=MaskingBackend.TRF,
        failure_policy=MaskingFailurePolicy.STRICT,
        status="success",
        legacy_adapter=False,
        backend_versions=(("trf", "4.09.1"),),
        masked_bases=4,
        repeatmasker_species=None,
        repeatmasker_library=None,
        repeatmasker_library_sha256=None,
        configured_fallback_backend=None,
        fallback_backend=None,
        fallback_reason=None,
        input_sha256=hashlib.sha256(genome.read_bytes()).hexdigest(),
        output_sha256=hashlib.sha256(masked.read_bytes()).hexdigest(),
    )
    write_masking_status(result, masking_dir)
    (old_root / "phase0" / "proteome.fasta").write_text(">p\nM\n")
    (old_root / "phase0" / "reuse_provenance.json").write_text('{"earlier": "import"}')
    (old_root / "phase0" / "genes.gff").write_text("s\ttest\tCDS\t1\t3\t.\t+\t0\tID=p\n")
    before = _tree(old_root)
    source = ReusedRun(
        root=old_root,
        run_fingerprint="1" * 64,
        code_sha256=OLD_CODE,
        markers={0: "a" * 64, 1: "b" * 64},
        artifacts={name: digest for name, digest in before.items() if digest is not None},
    )

    import_phase0(source, output_dir=new_root, genome_path=genome)

    status = new_root / "phase0" / "masking" / "masking_status.json"
    loaded = load_masking_result(status, expected_input=genome)
    assert loaded.output_path == (new_root / "phase0" / "masking" / "genome.masked.fna").resolve()
    assert loaded.masked_bases == 4
    assert (new_root / "phase0" / "masking" / "genome.masked.fna").read_bytes() == masked.read_bytes()
    assert (new_root / "phase0" / "proteome.fasta").read_bytes() == b">p\nM\n"
    loaded_phase0 = orchestrator._load_verified_phase0(
        genome_path=genome, output_dir=new_root, masking=MaskingConfig(backend=MaskingBackend.TRF)
    )
    assert loaded_phase0["masked_path"] == loaded.output_path
    assert _tree(old_root) == before
    provenance = json.loads((new_root / "phase0" / "reuse_provenance.json").read_text())
    assert [item["path"] for item in provenance["rewritten"]] == ["phase0/masking/masking_status.json"]
    assert "phase0/masking/masking_status.json" not in {item["path"] for item in provenance["copied"]}
    assert provenance["old_run_fingerprint"] == "1" * 64  # the record of the new import, not the earlier one
    assert "phase0/reuse_provenance.json" not in {item["path"] for item in provenance["copied"]}


def test_phase0_import_refuses_a_symlink_in_the_old_run(tmp_path: Path) -> None:
    old_root = tmp_path / "old" / "g"
    (old_root / "phase0").mkdir(parents=True)
    target = tmp_path / "elsewhere.txt"
    target.write_text("x")
    (old_root / "phase0" / "link").symlink_to(target)
    source = ReusedRun(root=old_root, run_fingerprint="1" * 64, code_sha256=OLD_CODE, markers={})

    with pytest.raises(Phase1ReuseError, match="symlink"):
        import_phase0(source, output_dir=tmp_path / "new" / "g", genome_path=tmp_path / "g.fna")


@pytest.mark.parametrize("relative", SEARCH_FILES)
def test_repair_changed_consumed_file_is_refused(tmp_path, monkeypatch, relative):
    old, new = _make_worlds(tmp_path, monkeypatch, search_files=True)
    path = old.output_dir / relative
    path.write_bytes(path.read_bytes() + b"changed after completion\n")
    with pytest.raises(Phase1ReuseError):
        phase1_reuse.open_reused_run(_reuse(old), genome_id="demo", identity=new.identity)


def test_repair_unbound_consumed_files_are_refused(tmp_path, monkeypatch):
    old, new = _make_worlds(tmp_path, monkeypatch, search_files=True)
    path = old.output_dir / "phase1.complete.json"
    marker = json.loads(path.read_text())
    marker["artifacts"] = [a for a in marker["artifacts"] if a["relative_path"] not in SEARCH_FILES]
    path.write_text(json.dumps(marker))
    with pytest.raises(Phase1ReuseError, match="authenticated"):
        phase1_reuse.open_reused_run(_reuse(old), genome_id="demo", identity=new.identity)


def test_repair_copy_checks_validated_expected_hash(tmp_path, monkeypatch):
    old, new = _make_worlds(tmp_path, monkeypatch, search_files=True)
    source = phase1_reuse.open_reused_run(_reuse(old), genome_id="demo", identity=new.identity)
    path = old.output_dir / "phase0/proteome.fasta"
    path.write_bytes(path.read_bytes() + b"changed between check and copy\n")
    with pytest.raises(Phase1ReuseError, match="changed"):
        import_phase0(source, output_dir=new.output_dir, genome_path=new.genome)


def test_repair_running_old_run_is_refused(tmp_path, monkeypatch):
    old, new = _make_worlds(tmp_path, monkeypatch, search_files=True)
    path = old.output_dir / "virosync_run_state.json"
    state = json.loads(path.read_text())
    state["status"] = "running"
    path.write_text(json.dumps(state))
    with pytest.raises(Phase1ReuseError, match="running"):
        phase1_reuse.open_reused_run(_reuse(old), genome_id="demo", identity=new.identity)


@pytest.mark.parametrize("relative", ["phase0", "phase1"])
def test_repair_symlink_phase_root_is_refused(tmp_path, monkeypatch, relative):
    old, new = _make_worlds(tmp_path, monkeypatch, search_files=True)
    path = old.output_dir / relative
    moved = old.output_dir / (relative + "-moved")
    path.rename(moved)
    path.symlink_to(moved, target_is_directory=True)
    with pytest.raises(Phase1ReuseError, match="symlink"):
        phase1_reuse.open_reused_run(_reuse(old), genome_id="demo", identity=new.identity)


def test_repair_fresh_phase0_cannot_receive_imported_phase1(tmp_path, monkeypatch):
    from virosync.orchestration._flows.single_genome.run_state import invalidate_from_phase

    old, new = _make_worlds(tmp_path, monkeypatch, search_files=True)
    new.search_files = True
    assert new.run()["success"]
    invalidate_from_phase(new.output_dir, from_phase=1)
    with pytest.raises(Phase1ReuseError, match="origin"):
        _flow(new, _reuse(old))


def test_repair_empty_nonterminal_checkpoint_is_refused(tmp_path, monkeypatch):
    old, new = _make_worlds(tmp_path, monkeypatch, search_files=False)
    with pytest.raises(Phase1ReuseError, match="empty"):
        phase1_reuse.open_reused_run(_reuse(old), genome_id="demo", identity=new.identity)


def test_repair_output_lookup_cannot_cross_run_root(tmp_path):
    from virosync.pipeline.phase3.output_generator import OutputGenerator

    output = tmp_path / "new" / "phase3_synthesis"
    output.mkdir(parents=True)
    stray = tmp_path / "phase1/marker_validation/diamond_top10_taxonomy.tsv"
    stray.parent.mkdir(parents=True)
    stray.write_text("unowned search\n")
    generator = OutputGenerator(output_dir=output, run_root=output.parent)
    assert generator._find_diamond_results() is None


@pytest.mark.parametrize(
    "relative",
    ["phase1/novelty/diamond/diamond_combined.tsv", "phase1/hhg/hmm_hits_validated.tsv", "phase1/pfam_arbitration.tsv"],
)
def test_repair_existing_unbound_fallback_is_refused(worlds, relative):
    old, new = worlds
    path = old.output_dir / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("unbound\n")
    with pytest.raises(Phase1ReuseError, match="authenticated"):
        phase1_reuse.open_reused_run(_reuse(old), genome_id="demo", identity=new.identity)


def _drop_search_hashes(old: _MockPipeline) -> None:
    """Make the earlier run look like one made by an older source: its markers lack the search file hashes."""
    path = old.output_dir / "phase1.complete.json"
    marker = json.loads(path.read_text())
    marker["artifacts"] = [a for a in marker["artifacts"] if a["relative_path"] not in SEARCH_FILES]
    path.write_text(json.dumps(marker))


def _sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _snapshot_lines(old: _MockPipeline, *, genome: str = "demo", prefix: str = "", mark: str = "  ") -> list[str]:
    return [f"{_sha256_of(old.output_dir / relative)}{mark}{prefix}{genome}/{relative}" for relative in SEARCH_FILES]


def _write_snapshot(tmp_path: Path, lines: list[str]) -> tuple[Path, str]:
    path = tmp_path / "snapshot.sha256"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path, _sha256_of(path)


def _snapshot_reuse(old: _MockPipeline, manifest: Path, manifest_sha256: str) -> Phase1Reuse:
    return Phase1Reuse(
        old_base=old.output_dir.parent,
        code_sha256=OLD_CODE,
        snapshot_manifest=manifest,
        snapshot_sha256=manifest_sha256,
    )


def _open_with_snapshot(old: _MockPipeline, new: _MockPipeline, manifest: Path, manifest_sha256: str):
    request = _snapshot_reuse(old, manifest, manifest_sha256)
    return phase1_reuse.open_reused_run(request, genome_id="demo", identity=new.identity)


def _provenance(run: _MockPipeline, phase: int) -> dict:
    return json.loads((run.output_dir / f"phase{phase}/reuse_provenance.json").read_text())


def test_an_old_run_without_search_hashes_is_refused_without_a_snapshot(
    search_worlds: tuple[_MockPipeline, _MockPipeline],
) -> None:
    old, new = search_worlds
    _drop_search_hashes(old)

    with pytest.raises(Phase1ReuseError, match="no authenticated old phase artifact"):
        _flow(new, _reuse(old))


def test_a_snapshot_lets_an_old_run_without_search_hashes_be_reused_and_copies_the_files(
    search_worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = search_worlds
    _drop_search_hashes(old)
    manifest, manifest_sha256 = _write_snapshot(tmp_path, _snapshot_lines(old))

    result = _flow(new, _snapshot_reuse(old, manifest, manifest_sha256))

    assert result["success"] is True
    for relative in SEARCH_FILES:
        assert (new.output_dir / relative).read_bytes() == (old.output_dir / relative).read_bytes()


def test_the_phase1_provenance_names_the_manifest_and_the_files_it_authenticated(
    search_worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = search_worlds
    _drop_search_hashes(old)
    manifest, manifest_sha256 = _write_snapshot(tmp_path, _snapshot_lines(old))

    _flow(new, _snapshot_reuse(old, manifest, manifest_sha256))

    provenance = _provenance(new, 1)
    assert provenance["snapshot_manifest_sha256"] == manifest_sha256
    assert provenance["snapshot_authenticated"] == {
        relative: _sha256_of(old.output_dir / relative) for relative in SEARCH_FILES
    }


def test_the_phase1_provenance_says_that_a_snapshot_authenticated_the_files(
    search_worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = search_worlds
    _drop_search_hashes(old)
    manifest, manifest_sha256 = _write_snapshot(tmp_path, _snapshot_lines(old))

    _flow(new, _snapshot_reuse(old, manifest, manifest_sha256))

    assert "snapshot" in _provenance(new, 1)["authentication"]
    assert _provenance(new, 1)["schema"] == phase1_reuse.PROVENANCE_SCHEMA


def test_the_phase0_provenance_has_no_snapshot_fields(
    search_worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = search_worlds
    _drop_search_hashes(old)
    manifest, manifest_sha256 = _write_snapshot(tmp_path, _snapshot_lines(old))

    _flow(new, _snapshot_reuse(old, manifest, manifest_sha256))

    phase0 = _provenance(new, 0)
    assert "snapshot_manifest_sha256" not in phase0
    assert "snapshot_authenticated" not in phase0
    assert "snapshot" not in phase0["authentication"]


@pytest.mark.parametrize(("prefix", "mark"), [("./", "  "), ("", " *")], ids=["dot-slash", "binary-mark"])
def test_a_snapshot_line_may_start_with_dot_slash_or_use_the_binary_mark(
    search_worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path, prefix: str, mark: str
) -> None:
    old, new = search_worlds
    _drop_search_hashes(old)
    manifest, manifest_sha256 = _write_snapshot(tmp_path, _snapshot_lines(old, prefix=prefix, mark=mark))

    source = _open_with_snapshot(old, new, manifest, manifest_sha256)

    assert source.snapshot_authenticated == {
        relative: _sha256_of(old.output_dir / relative) for relative in SEARCH_FILES
    }


def test_a_snapshot_hash_that_differs_from_the_file_is_refused(
    search_worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = search_worlds
    _drop_search_hashes(old)
    lines = _snapshot_lines(old)
    lines[0] = "0" * 64 + lines[0][64:]
    manifest, manifest_sha256 = _write_snapshot(tmp_path, lines)

    with pytest.raises(Phase1ReuseError, match="authenticated input changed"):
        _flow(new, _snapshot_reuse(old, manifest, manifest_sha256))

    assert not new.output_dir.exists()


def test_a_file_that_changes_between_validation_and_copy_is_refused_with_a_snapshot(
    search_worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = search_worlds
    _drop_search_hashes(old)
    manifest, manifest_sha256 = _write_snapshot(tmp_path, _snapshot_lines(old))
    source = _open_with_snapshot(old, new, manifest, manifest_sha256)
    path = old.output_dir / SEARCH_FILES[0]
    path.write_bytes(path.read_bytes() + b"changed between check and copy\n")

    with pytest.raises(Phase1ReuseError, match="changed during copy"):
        import_phase1(
            source,
            output_dir=new.output_dir,
            genome_id="demo",
            config=PipelineConfig(),
            masked_path=new.genome,
            proteome_path=new.genome,
            config_fingerprint=new.fingerprint,
            logger=logging.getLogger("test"),
        )


def test_a_manifest_with_another_sha256_than_the_stated_one_is_refused(
    search_worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = search_worlds
    _drop_search_hashes(old)
    manifest, _ = _write_snapshot(tmp_path, _snapshot_lines(old))

    with pytest.raises(Phase1ReuseError, match="SHA-256"):
        _flow(new, _snapshot_reuse(old, manifest, "0" * 64))

    assert not new.output_dir.exists()


def test_a_snapshot_hash_that_conflicts_with_a_recorded_hash_is_refused(
    search_worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = search_worlds
    lines = _snapshot_lines(old)
    lines[1] = "0" * 64 + lines[1][64:]
    manifest, manifest_sha256 = _write_snapshot(tmp_path, lines)

    with pytest.raises(Phase1ReuseError, match="snapshot conflicts with the recorded hash"):
        _flow(new, _snapshot_reuse(old, manifest, manifest_sha256))


def test_a_snapshot_of_another_genome_does_not_authenticate_this_genome(
    search_worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = search_worlds
    _drop_search_hashes(old)
    manifest, manifest_sha256 = _write_snapshot(tmp_path, _snapshot_lines(old, genome="other"))

    with pytest.raises(Phase1ReuseError, match="no authenticated old phase artifact"):
        _flow(new, _snapshot_reuse(old, manifest, manifest_sha256))


HOST_MODEL = "phase1/marker_validation/host_signature_model.json"


def _add_unrecorded_host_model(old: _MockPipeline) -> Path:
    """Add a consumed Phase 1 file that the phase markers of the earlier run do not record."""
    path = old.output_dir / HOST_MODEL
    path.write_text("{}\n", encoding="utf-8")
    return path


def test_a_consumed_file_without_a_recorded_hash_is_refused_without_a_snapshot_entry(
    worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = worlds
    _add_unrecorded_host_model(old)
    manifest, manifest_sha256 = _write_snapshot(tmp_path, [f"{'a' * 64}  demo/phase1/other.tsv"])

    with pytest.raises(Phase1ReuseError, match=f"no authenticated old phase artifact for {HOST_MODEL}"):
        _flow(new, _snapshot_reuse(old, manifest, manifest_sha256))


def test_a_consumed_file_without_a_recorded_hash_is_accepted_with_a_snapshot_entry(
    worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = worlds
    model = _add_unrecorded_host_model(old)
    manifest, manifest_sha256 = _write_snapshot(tmp_path, [f"{_sha256_of(model)}  demo/{HOST_MODEL}"])

    assert _flow(new, _snapshot_reuse(old, manifest, manifest_sha256))["success"] is True


def test_the_phase1_provenance_lists_a_rebuilt_file_that_the_snapshot_authenticated(
    worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = worlds
    model = _add_unrecorded_host_model(old)
    manifest, manifest_sha256 = _write_snapshot(tmp_path, [f"{_sha256_of(model)}  demo/{HOST_MODEL}"])

    _flow(new, _snapshot_reuse(old, manifest, manifest_sha256))

    provenance = _provenance(new, 1)
    assert provenance["snapshot_authenticated"] == {HOST_MODEL: _sha256_of(model)}
    assert HOST_MODEL not in {item["path"] for item in provenance["copied"]}


def test_a_snapshot_entry_for_a_file_outside_the_consumed_files_authenticates_nothing(
    worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = worlds
    unconsumed = old.output_dir / "phase1/unconsumed.tsv"
    unconsumed.write_text("unconsumed\n")
    manifest, manifest_sha256 = _write_snapshot(tmp_path, [f"{_sha256_of(unconsumed)}  demo/phase1/unconsumed.tsv"])

    source = _open_with_snapshot(old, new, manifest, manifest_sha256)

    assert "phase1/unconsumed.tsv" not in source.artifacts
    assert source.snapshot_authenticated == {}


def test_a_snapshot_entry_for_an_absent_consumed_file_authenticates_nothing(
    worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = worlds
    absent = "phase1/novelty/diamond/diamond_combined.tsv"
    manifest, manifest_sha256 = _write_snapshot(tmp_path, [f"{'0' * 64}  demo/{absent}"])

    source = _open_with_snapshot(old, new, manifest, manifest_sha256)

    assert not (old.output_dir / absent).exists()
    assert source.snapshot_authenticated == {}


def test_a_snapshot_that_matches_the_recorded_hashes_changes_no_provenance_byte(
    search_worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = search_worlds
    manifest, manifest_sha256 = _write_snapshot(tmp_path, _snapshot_lines(old))
    _flow(new, _reuse(old))
    plain = {phase: (new.output_dir / f"phase{phase}/reuse_provenance.json").read_bytes() for phase in (0, 1)}

    _flow(new, _snapshot_reuse(old, manifest, manifest_sha256), resume=False)

    with_snapshot = {phase: (new.output_dir / f"phase{phase}/reuse_provenance.json").read_bytes() for phase in (0, 1)}
    assert with_snapshot == plain


@pytest.mark.parametrize(
    "line",
    [
        "{sha}",
        "{sha} demo/phase1/hmm/hmm_hits_validated.tsv",
        "{upper}  demo/phase1/hmm/hmm_hits_validated.tsv",
        "{short}  demo/phase1/hmm/hmm_hits_validated.tsv",
        "{sha}  ",
        "",
    ],
    ids=["no-path", "one-space", "uppercase-hash", "short-hash", "empty-path", "blank-line"],
)
def test_a_malformed_snapshot_line_is_refused(
    search_worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path, line: str
) -> None:
    old, new = search_worlds
    _drop_search_hashes(old)
    sha = "a" * 64
    manifest, manifest_sha256 = _write_snapshot(
        tmp_path, [*_snapshot_lines(old), line.format(sha=sha, upper=sha.upper(), short=sha[:63])]
    )

    with pytest.raises(Phase1ReuseError, match="malformed snapshot line"):
        _open_with_snapshot(old, new, manifest, manifest_sha256)


def test_a_path_listed_twice_in_a_snapshot_with_different_hashes_is_refused(
    search_worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = search_worlds
    _drop_search_hashes(old)
    lines = _snapshot_lines(old)
    manifest, manifest_sha256 = _write_snapshot(tmp_path, [*lines, "0" * 64 + "  ./" + lines[0][66:]])

    with pytest.raises(Phase1ReuseError, match="twice with different hashes"):
        _open_with_snapshot(old, new, manifest, manifest_sha256)


def test_a_symlinked_snapshot_is_refused(search_worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path) -> None:
    old, new = search_worlds
    _drop_search_hashes(old)
    manifest, manifest_sha256 = _write_snapshot(tmp_path, _snapshot_lines(old))
    link = tmp_path / "link.sha256"
    link.symlink_to(manifest)

    with pytest.raises(Phase1ReuseError, match="regular file"):
        _open_with_snapshot(old, new, link, manifest_sha256)


def test_a_snapshot_that_is_a_folder_is_refused(
    search_worlds: tuple[_MockPipeline, _MockPipeline], tmp_path: Path
) -> None:
    old, new = search_worlds
    _drop_search_hashes(old)

    with pytest.raises(Phase1ReuseError, match="regular file"):
        _open_with_snapshot(old, new, tmp_path, "0" * 64)


@pytest.mark.parametrize(
    "fields",
    [{"snapshot_manifest": Path("m.sha256")}, {"snapshot_sha256": "0" * 64}],
    ids=["manifest-only", "sha-only"],
)
def test_a_request_needs_both_snapshot_fields(tmp_path: Path, fields: dict) -> None:
    with pytest.raises(Phase1ReuseError, match="together"):
        Phase1Reuse(old_base=tmp_path, code_sha256=OLD_CODE, **fields)


def test_a_request_rejects_a_malformed_snapshot_sha256(tmp_path: Path) -> None:
    with pytest.raises(Phase1ReuseError, match="snapshot SHA-256"):
        Phase1Reuse(old_base=tmp_path, code_sha256=OLD_CODE, snapshot_manifest=Path("m.sha256"), snapshot_sha256="ABC")


def _cli_snapshot(tmp_path: Path) -> tuple[Path, Path]:
    old_base = tmp_path / "campaign"
    old_base.mkdir()
    manifest = tmp_path / "snapshot.sha256"
    manifest.write_text(f"{'a' * 64}  g/phase1/hmm/hmm_hits_validated.tsv\n", encoding="utf-8")
    return old_base, manifest


def test_cli_passes_the_snapshot_to_the_request(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    old_base, manifest = _cli_snapshot(tmp_path)

    result, received = _invoke_run(
        tmp_path,
        monkeypatch,
        *("--reuse-phase1-from", str(old_base), "--reuse-phase1-code-sha256", OLD_CODE),
        *("--reuse-phase1-snapshot", str(manifest), "--reuse-phase1-snapshot-sha256", _sha256_of(manifest)),
    )

    assert result.exit_code == 0, result.output
    assert received["phase1_reuse"] == Phase1Reuse(
        old_base=old_base,
        code_sha256=OLD_CODE,
        snapshot_manifest=manifest,
        snapshot_sha256=_sha256_of(manifest),
    )


def test_cli_refuses_a_snapshot_with_another_sha256_before_any_genome_runs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    old_base, manifest = _cli_snapshot(tmp_path)

    result, received = _invoke_run(
        tmp_path,
        monkeypatch,
        *("--reuse-phase1-from", str(old_base), "--reuse-phase1-code-sha256", OLD_CODE),
        *("--reuse-phase1-snapshot", str(manifest), "--reuse-phase1-snapshot-sha256", "0" * 64),
    )

    assert result.exit_code != 0
    assert "does not match --reuse-phase1-snapshot-sha256" in result.output
    assert not received


@pytest.mark.parametrize(
    "options",
    [
        ("--reuse-phase1-snapshot", "{manifest}"),
        ("--reuse-phase1-snapshot-sha256", "0" * 64),
    ],
    ids=["path-only", "sha-only"],
)
def test_cli_needs_both_snapshot_options(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, options: tuple) -> None:
    old_base, manifest = _cli_snapshot(tmp_path)

    result, received = _invoke_run(
        tmp_path,
        monkeypatch,
        *("--reuse-phase1-from", str(old_base), "--reuse-phase1-code-sha256", OLD_CODE),
        *(item.format(manifest=manifest) for item in options),
    )

    assert result.exit_code != 0
    assert "required together" in result.output
    assert not received


def test_cli_needs_the_reuse_folder_for_the_snapshot_options(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _, manifest = _cli_snapshot(tmp_path)

    result, received = _invoke_run(
        tmp_path,
        monkeypatch,
        *("--reuse-phase1-snapshot", str(manifest), "--reuse-phase1-snapshot-sha256", _sha256_of(manifest)),
    )

    assert result.exit_code != 0
    assert "need --reuse-phase1-from" in result.output
    assert not received
