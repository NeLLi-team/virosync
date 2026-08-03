from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from virosync.ablation import (
    ABLATION_CONTRACT_SHA256,
    AblationCounters,
    AblationEvents,
    AblationID,
    validate_ablation_events_bytes,
)
from virosync.config import MaskingBackend, MaskingConfig, MaskingFailurePolicy
from virosync.orchestration._flows.single_genome import orchestrator
from virosync.orchestration._flows.single_genome import (
    run_state as run_state_module,
)
from virosync.orchestration._flows.single_genome.manifest import (
    _empty_prediction_summary,
    _write_empty_run_log,
)
from virosync.orchestration._flows.single_genome.phase_state import (
    load_phase2_state,
    write_phase2_state,
)
from virosync.orchestration._flows.single_genome.phase1_state import (
    load_phase1_state,
    write_phase1_state,
)
from virosync.orchestration._flows.single_genome.phase2_resume_state import (
    write_phase2_resume_state,
)
from virosync.orchestration._flows.single_genome.run_state import (
    PHASE_MARKER_FILENAMES,
    RUN_STATE_FILENAME,
    atomic_write_json,
    build_artifact_identity,
    build_input_identity,
    canonical_sha256,
    compute_run_fingerprint,
    load_run_state,
    marker_sha256,
    plan_resume,
)
from virosync.output_contract import (
    EFFECTIVE_EVE_CLASSES,
    OUTPUT_SCHEMA_VERSION,
    normalize_effective_eve_class,
)
from virosync.pipeline.phase0.masking import (
    MaskingResult,
    load_masking_result,
    mask_genome_pipeline,
    write_masking_status,
)
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.host_signatures import HostSignatureModel
from virosync.pipeline.phase1.seed_merger import MergedSeed


_PREDICTION_HEADER = (
    "eve_id\tscaffold\tstart\tend\tconfidence_tier\tlength\t"
    "gene_taxonomy_total\thallmark_total\tfinal_confidence\t"
    "effective_eve_class\n"
)
class _NullResourceMonitor:
    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        return False


class _FakeOutputGenerator:
    def __init__(self, output_dir, genome_fasta=None):
        self.genome_sequences = {}

    def write_combined_eve_fasta(self, results, path):
        Path(path).write_text(">EVE_1\nACGT\n")


class _MockPipeline:
    """Small deterministic phase implementation for orchestration-state tests."""

    def __init__(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
        self.genome = tmp_path / "demo.fna"
        self.genome.write_text(">scaffold\nACGT\n")
        self.output_dir = tmp_path / "results" / "demo"
        self.calls = {f"phase{phase}": 0 for phase in range(4)}
        self.resume_flags = {"phase1": [], "phase2": [], "phase3": []}
        self.events: list[str] = []
        self.fail_phase3_once = False
        self.terminal_phase: int | None = None
        self.candidate_only = False
        self.effective_eve_class = "NCLDV"
        self.phase3_boundaries: list[RefinedBoundary] = []
        self.require_combined_fasta_before_report = False

        self.identity, self.fingerprint = self._run_identity()
        monkeypatch.setattr(
            orchestrator,
            "_build_run_identity",
            lambda **kwargs: (self.identity, self.fingerprint),
        )
        monkeypatch.setattr(
            orchestrator,
            "_run_phase0_subflow",
            self._phase0,
        )
        monkeypatch.setattr(
            orchestrator,
            "_run_phase1_subflow",
            self._phase1,
        )
        monkeypatch.setattr(
            orchestrator,
            "_run_phase2_subflow",
            self._phase2,
        )
        monkeypatch.setattr(
            orchestrator,
            "_run_phase3_subflow",
            self._phase3,
        )
        monkeypatch.setattr(
            orchestrator,
            "_generate_required_reports",
            self._generate_required_reports,
        )
        monkeypatch.setattr(orchestrator, "call_task", self._call_task)

        from virosync.orchestration import resource_monitor
        from virosync.pipeline.phase3 import output_generator
        from virosync.utils import gpu, provenance

        monkeypatch.setattr(resource_monitor, "ResourceMonitor", _NullResourceMonitor)
        monkeypatch.setattr(output_generator, "OutputGenerator", _FakeOutputGenerator)
        monkeypatch.setattr(gpu, "release_gpu_memory", lambda: None)
        monkeypatch.setattr(provenance, "write_provenance", lambda *args, **kwargs: None)
        monkeypatch.setattr(
            orchestrator,
            "enforce_tsv_invariants",
            self._enforce_tsv_invariants,
        )
        self._install_state_spies(monkeypatch)

    def _run_identity(self) -> tuple[dict[str, object], str]:
        environment_payload = {
            "lock_sha256": "4" * 64,
            "runtime_sha256": "5" * 64,
            "requested_device": "cpu",
            "effective_device": "cpu",
        }
        identity = {
            "genome_id": "demo",
            "input_path": str(self.genome.resolve()),
            "output_dir": str(self.output_dir.resolve()),
            "input": asdict(build_input_identity(self.genome)),
            "config": {
                "sha256": "2" * 64,
                "ablation_id": "A0",
                "ablation_contract_sha256": ABLATION_CONTRACT_SHA256,
            },
            "code": {"version": "test", "source_sha256": "3" * 64},
            "environment": {
                **environment_payload,
                "sha256": canonical_sha256(environment_payload),
            },
            "coordinate_schema_version": 2,
            "coordinate_convention": "0-based, half-open [start, end)",
            "output_schema_version": OUTPUT_SCHEMA_VERSION,
            "summary_schema_version": 3,
            "requested_masking": orchestrator._masking_request_identity(
                MaskingConfig()
            ),
            "resources": [],
        }
        return identity, compute_run_fingerprint(identity)

    def _install_state_spies(self, monkeypatch: pytest.MonkeyPatch) -> None:
        original_started = orchestrator.publish_run_started
        original_phase = orchestrator.publish_phase_completion
        original_success = orchestrator.publish_run_success

        def publish_started(*args, **kwargs):
            state = original_started(*args, **kwargs)
            self.events.append("running")
            return state

        def publish_phase(output_dir, record):
            assert record.artifacts
            assert all(
                (Path(output_dir) / artifact.relative_path).is_file()
                for artifact in record.artifacts
            )
            path = original_phase(output_dir, record)
            self.events.append(f"phase{record.phase}")
            return path

        def publish_success(output_dir, **kwargs):
            assert (Path(output_dir) / PHASE_MARKER_FILENAMES[-1]).is_file() or (
                kwargs["result"]["terminal_phase"] is not None
            )
            assert all(
                (Path(output_dir) / artifact.relative_path).is_file()
                for artifact in kwargs["artifacts"]
            )
            state = original_success(output_dir, **kwargs)
            self.events.append("success")
            return state

        monkeypatch.setattr(orchestrator, "publish_run_started", publish_started)
        monkeypatch.setattr(orchestrator, "publish_phase_completion", publish_phase)
        monkeypatch.setattr(orchestrator, "publish_run_success", publish_success)

    def _phase0(self, **kwargs):
        self.calls["phase0"] += 1
        phase_dir = Path(kwargs["output_dir"]) / "phase0"
        masking_result = mask_genome_pipeline(
            self.genome,
            phase_dir / "masking",
            config=MaskingConfig(),
        )
        proteome = phase_dir / "proteome.fasta"
        genes = phase_dir / "genes.gff"
        proteome.write_text(">p1\nM\n")
        genes.write_text("scaffold\ttest\tCDS\t1\t4\t.\t+\t0\tID=p1\n")
        return {
            "masked_path": masking_result.output_path,
            "repeat_regions": [],
            "proteome_path": proteome,
            "n_genes": 1,
            "elapsed": 0.0,
            "masking_result": load_masking_result(masking_result.status_path),
        }

    def _phase1(self, **kwargs):
        self.calls["phase1"] += 1
        self.resume_flags["phase1"].append(bool(kwargs["resume"]))
        phase_dir = Path(kwargs["output_dir"]) / "phase1"
        phase_dir.mkdir(parents=True, exist_ok=True)
        state_path = phase_dir / "resume_state.json"
        merged_seeds = [
            MergedSeed(
                scaffold="scaffold",
                start=0,
                end=4,
                seed_id="seed-1",
                sources=["integration-test"],
            )
        ]
        host_model = HostSignatureModel()
        if kwargs["resume"]:
            state = load_phase1_state(state_path)
            merged_seeds = state.merged_seeds
            host_model = state.host_signature_model
        else:
            write_phase1_state(
                state_path,
                validated_markers=[],
                merged_seeds=merged_seeds,
                host_signature_model=host_model,
                host_signatures=set(),
                host_deviation_summary=None,
            )
        if self.terminal_phase == 1:
            self._write_predictions(Path(kwargs["output_dir"]), rows=0)
            reports = self._generate_required_reports(
                output_dir=kwargs["output_dir"]
            )
            _write_empty_run_log(
                output_dir=Path(kwargs["output_dir"]),
                genome_id=kwargs["genome_id"],
                reason="no markers",
                elapsed_sec=0.0,
                input_path=self.genome,
                output_files=reports,
                fingerprint=self.fingerprint,
            )
            return {
                "genome_id": kwargs["genome_id"],
                "success": True,
                **_empty_prediction_summary(),
                "output_files": {},
                "elapsed_sec": 0.0,
            }
        return {
            "merged_seeds": merged_seeds,
            "validated_markers": [],
            "host_signature_model": host_model,
            "host_signatures": set(),
            "background": None,
            "gene_data": {},
            "host_deviation_summary": {},
            "elapsed": 0.0,
        }

    def _phase2(self, **kwargs):
        self.calls["phase2"] += 1
        resumed = bool(kwargs["resume"])
        self.resume_flags["phase2"].append(resumed)
        phase_dir = Path(kwargs["output_dir"]) / "phase2"
        phase_dir.mkdir(parents=True, exist_ok=True)
        state_path = phase_dir / "refined_state.json"
        bed_path = phase_dir / "refined_boundaries.bed"
        if resumed:
            boundaries = load_phase2_state(state_path)
        elif self.terminal_phase == 2:
            boundaries = []
            write_phase2_state(state_path, boundaries)
            write_phase2_resume_state(
                phase_dir / "resume_state.json",
                refined_boundaries=boundaries,
                boundary_taxonomy_map={},
                boundary_control_stats=None,
                boundary_diamond_query=None,
            )
            bed_path.write_text("")
            self._write_predictions(Path(kwargs["output_dir"]), rows=0)
            reports = self._generate_required_reports(
                output_dir=kwargs["output_dir"]
            )
            _write_empty_run_log(
                output_dir=Path(kwargs["output_dir"]),
                genome_id=kwargs["genome_id"],
                reason="no boundaries",
                elapsed_sec=0.0,
                input_path=self.genome,
                output_files=reports,
                fingerprint=self.fingerprint,
            )
            return {
                "genome_id": kwargs["genome_id"],
                "success": True,
                **_empty_prediction_summary(),
                "output_files": {},
                "elapsed_sec": 0.0,
            }
        else:
            boundaries = [
                RefinedBoundary(
                    scaffold="scaffold",
                    start=0,
                    end=4,
                    seed_id="lossless-seed",
                    seed_sources=["integration-test"],
                    hallmark_genes=["MCP"],
                    confidence=0.91,
                )
            ]
            write_phase2_state(state_path, boundaries)
            write_phase2_resume_state(
                phase_dir / "resume_state.json",
                refined_boundaries=boundaries,
                boundary_taxonomy_map={},
                boundary_control_stats=None,
                boundary_diamond_query=None,
            )
            bed_path.write_text(
                "scaffold\t0\t4\tEVE_scaffold_0-4\t910\t.\n"
            )
        return {
            "refined_boundaries": boundaries,
            "boundary_taxonomy_map": {},
            "boundary_control_stats": None,
            "boundary_diamond_query": None,
            "proteome_index": {},
            "goto_phase3": resumed,
            "boundaries_bed": bed_path,
            "elapsed": 0.0,
        }

    def _phase3(self, **kwargs):
        self.calls["phase3"] += 1
        self.resume_flags["phase3"].append(bool(kwargs["resume"]))
        self.phase3_boundaries = list(kwargs["refined_boundaries"])
        phase_dir = Path(kwargs["output_dir"]) / "phase3"
        phase_dir.mkdir(parents=True, exist_ok=True)
        (phase_dir / "state.tsv").write_text("state\nready\n")
        if self.fail_phase3_once:
            self.fail_phase3_once = False
            raise RuntimeError("injected Phase 3 failure")
        result = SimpleNamespace(
            coherence_analysis=None,
            confidence_tier="HIGH",
            eve_id="EVE_1",
            start=0,
            end=4,
            gene_count=1,
            final_confidence=0.91,
        )
        accepted_results = [] if self.candidate_only else [result]
        return {
            "verification_results": [result],
            "accepted_results": accepted_results,
            "classification_stats": (
                {}
                if self.candidate_only
                else {normalize_effective_eve_class(self.effective_eve_class): 1}
            ),
            "accepted": len(accepted_results),
            "elapsed": 0.0,
            "precomputed_tmvec": None,
        }

    def _call_task(self, task_fn, *args, **kwargs):
        if task_fn is orchestrator.generate_outputs_task:
            output_dir = Path(kwargs["output_dir"])
            result = kwargs["verification_results"][0]
            self._write_predictions(
                output_dir,
                rows=0 if self.candidate_only else 1,
                detailed_rows=1,
                confidence_tier=result.confidence_tier,
                final_confidence=result.final_confidence,
                effective_eve_class=self.effective_eve_class,
            )
            return {
                "predictions_tsv": output_dir / "virosync_predictions.tsv",
                "predictions_detailed_tsv": (
                    output_dir / "virosync_predictions_detailed.tsv"
                ),
                "predictions_bed": output_dir / "virosync_predictions.bed",
                "predictions_gff": output_dir / "virosync_predictions.gff3",
                "summary_json": output_dir / "virosync_summary.json",
            }
        if task_fn is orchestrator.create_summary_artifact_task:
            return None
        raise AssertionError(f"unexpected task in mocked run: {task_fn}")

    @staticmethod
    def _enforce_tsv_invariants(*, report_out, **kwargs):
        detailed = Path(kwargs["detailed_tsv"])
        rows_checked = max(0, len(detailed.read_text().splitlines()) - 1)
        Path(report_out).write_text(
            "status\trows_checked\tissue_count\terror_count\twarning_count\n"
            f"PASS\t{rows_checked}\t0\t0\t0\n\n"
            "eve_id\tcheck\tseverity\tmessage\n"
        )
        return SimpleNamespace(
            warning_count=0,
            warning_issues=[],
            rows_checked=rows_checked,
            issue_count=0,
            error_count=0,
            status="PASS",
            issues=[],
        )

    @staticmethod
    def _write_predictions(
        directory: Path,
        *,
        rows: int,
        detailed_rows: int | None = None,
        confidence_tier: str = "HIGH",
        final_confidence: float = 0.91,
        effective_eve_class: str = "NCLDV",
    ) -> None:
        directory.mkdir(parents=True, exist_ok=True)
        if detailed_rows is None:
            detailed_rows = rows
        prediction_row = (
            "EVE_1\tscaffold\t0\t4\t"
            f"{confidence_tier}\t4\t1\t1\t{final_confidence}\t{effective_eve_class}\n"
        )
        canonical_content = _PREDICTION_HEADER + (prediction_row * rows)
        detailed_content = _PREDICTION_HEADER + (
            prediction_row * detailed_rows
        )
        export_score = int(min(1000, final_confidence * 1000))
        (directory / "virosync_predictions.tsv").write_text(canonical_content)
        (directory / "virosync_predictions_detailed.tsv").write_text(
            detailed_content
        )
        (directory / "virosync_predictions.bed").write_text(
            f"scaffold\t0\t4\tEVE_1\t{export_score}\t.\n" if rows else ""
        )
        (directory / "virosync_predictions.gff3").write_text(
            "##gff-version 3\n"
            + (
                f"scaffold\tViroSync\tEVE\t1\t4\t{export_score}\t.\t.\t"
                f"ID=EVE_1;Name=EVE_1;confidence={final_confidence}\n"
                if rows
                else ""
            )
        )
        (directory / "virosync_summary.json").write_text(
            json.dumps(
                {
                    "virosync_version": "test",
                    "coordinate_schema_version": 2,
                    "coordinate_convention": "0-based, half-open [start, end)",
                    "output_schema_version": OUTPUT_SCHEMA_VERSION,
                    "statistics": {
                        "canonical_predictions": rows,
                        "total_candidates": detailed_rows,
                        "total_accepted_length_bp": 4 if rows else 0,
                        "high_confidence": (
                            rows if confidence_tier == "HIGH" else 0
                        ),
                        "medium_confidence": (
                            rows if confidence_tier == "MEDIUM" else 0
                        ),
                        "low_confidence": (
                            rows if confidence_tier == "LOW" else 0
                        ),
                        "promoted_low_confidence": 0,
                    },
                    "per_scaffold": {},
                }
            )
            + "\n"
        )

    def _generate_required_reports(self, *, output_dir, **kwargs):
        output_dir = Path(output_dir)
        if self.require_combined_fasta_before_report:
            assert (output_dir / "demo_eves.fna").is_file()
        notebook = Path(output_dir) / "notebooks" / "jupyter" / "eve_analysis.ipynb"
        notebook.parent.mkdir(parents=True, exist_ok=True)
        notebook.write_text(
            '{"cells":[],"metadata":{},"nbformat":4,"nbformat_minor":5}\n'
        )
        (output_dir / "host_signature_model.png").write_bytes(b"synthetic-png\n")
        (output_dir / "gvclass_results.tsv").write_text(
            "eve_id\tclassification\nEVE_1\tNCLDV\n",
            encoding="utf-8",
        )
        return {"eve_analysis_notebook": str(notebook)}

    def run(self) -> dict:
        return orchestrator.single_genome_flow(
            genome_path=self.genome,
            output_dir=self.output_dir,
            genome_id="demo",
            device="cpu",
            resume=True,
        )


@pytest.fixture
def mocked_pipeline(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> _MockPipeline:
    return _MockPipeline(tmp_path, monkeypatch)


def test_a1_publishes_nonzero_phase1_terminal_and_resumes_without_phase2_or_3(
    mocked_pipeline: _MockPipeline,
) -> None:
    mocked_pipeline.require_combined_fasta_before_report = True
    mocked_pipeline.identity["config"]["ablation_id"] = "A1"
    mocked_pipeline.fingerprint = compute_run_fingerprint(
        mocked_pipeline.identity
    )

    result = orchestrator.single_genome_flow(
        genome_path=mocked_pipeline.genome,
        output_dir=mocked_pipeline.output_dir,
        genome_id="demo",
        device="cpu",
        resume=True,
        ablation_id="A1",
    )

    assert result["success"] is True
    assert result["accepted"] == 1
    assert result["predictions"] == 1
    assert mocked_pipeline.calls == {
        "phase0": 1,
        "phase1": 1,
        "phase2": 0,
        "phase3": 0,
    }
    marker = json.loads(
        (mocked_pipeline.output_dir / PHASE_MARKER_FILENAMES[1]).read_text()
    )
    assert marker["outcome"] == "terminal_ablation"
    assert "demo_eves.fna" in {
        artifact["relative_path"] for artifact in marker["artifacts"]
    }
    assert not (
        mocked_pipeline.output_dir / PHASE_MARKER_FILENAMES[2]
    ).exists()
    state = load_run_state(mocked_pipeline.output_dir)
    assert state.result["terminal_phase"] == 1
    assert state.result["canonical_rows"] == 1
    assert state.result["tier_counts"] == {
        "HIGH": 0,
        "MEDIUM": 0,
        "LOW": 1,
    }
    assert state.result["promoted_low_rows"] == 0
    assert "demo_eves.fna" in {
        artifact.relative_path for artifact in state.artifacts
    }
    events = validate_ablation_events_bytes(
        (mocked_pipeline.output_dir / "ablation_events.json").read_bytes()
    )
    assert events.ablation_id.value == "A1"
    assert (
        events.counters.total_opportunities,
        events.counters.total_interventions,
        events.counters.total_changed,
    ) == (1, 1, 0)

    resumed = orchestrator.single_genome_flow(
        genome_path=mocked_pipeline.genome,
        output_dir=mocked_pipeline.output_dir,
        genome_id="demo",
        device="cpu",
        resume=True,
        ablation_id="A1",
    )
    assert resumed["accepted"] == 1
    assert mocked_pipeline.calls == {
        "phase0": 1,
        "phase1": 1,
        "phase2": 0,
        "phase3": 0,
    }

    eve_fasta = mocked_pipeline.output_dir / "demo_eves.fna"
    original = eve_fasta.read_bytes()
    eve_fasta.write_bytes(original.replace(b"ACGT", b"TGCA"))
    assert eve_fasta.stat().st_size == len(original)

    repaired = orchestrator.single_genome_flow(
        genome_path=mocked_pipeline.genome,
        output_dir=mocked_pipeline.output_dir,
        genome_id="demo",
        device="cpu",
        resume=True,
        ablation_id="A1",
    )
    assert repaired["accepted"] == 1
    assert mocked_pipeline.calls == {
        "phase0": 1,
        "phase1": 2,
        "phase2": 0,
        "phase3": 0,
    }
    assert eve_fasta.read_bytes() == original


_REQUIRED_FINAL_CLASSES = (
    ("ablation events", ("ablation_events.json",)),
    (
        "canonical predictions",
        (
            "phase3_synthesis/virosync_predictions.tsv",
            "virosync_predictions.tsv",
        ),
    ),
    (
        "detailed predictions",
        (
            "virosync_predictions_detailed.tsv",
            "phase3_synthesis/virosync_predictions_detailed.tsv",
        ),
    ),
    (
        "BED export",
        (
            "phase3_synthesis/virosync_predictions.bed",
            "virosync_predictions.bed",
        ),
    ),
    (
        "GFF3 export",
        (
            "phase3_synthesis/virosync_predictions.gff3",
            "virosync_predictions.gff3",
        ),
    ),
    (
        "summary",
        (
            "phase3_synthesis/virosync_summary.json",
            "virosync_summary.json",
        ),
    ),
    ("invariant report", ("virosync_tsv_invariant_report.tsv",)),
    ("run log", ("run.log",)),
    ("completion manifest", ("virosync_run_complete.json",)),
    ("executed notebook", ("notebooks/jupyter/eve_analysis.ipynb",)),
)

_REQUIRED_PHASE_ARTIFACTS = (
    (0, "phase0/ablation_events.json"),
    (0, "phase0/masking/masking_status.json"),
    (0, "phase0/proteome.fasta"),
    (0, "phase0/genes.gff"),
    (1, "phase1/ablation_events.json"),
    (1, "phase1/resume_state.json"),
    (2, "phase2/ablation_events.json"),
    (2, "phase2/refined_boundaries.bed"),
    (2, "phase2/refined_state.json"),
    (2, "phase2/resume_state.json"),
    (3, "phase3/ablation_events.json"),
)


def _mutate_recorded_file(path: Path, original: bytes, operation: str) -> None:
    if operation == "delete":
        path.unlink()
        return
    if operation == "truncate":
        path.write_bytes(original[: len(original) // 2])
        return
    if operation != "same-size":
        raise AssertionError(f"unknown mutation: {operation}")
    if not original:
        raise AssertionError(
            f"same-size mutation needs non-empty bytes: {path}"
        )
    changed = bytearray(original)
    changed[0] ^= 1
    path.write_bytes(changed)


def _assert_mutations_invalidate_at_phase(
    *,
    path: Path,
    expected_phase: int,
    mocked_pipeline: _MockPipeline,
) -> None:
    original = path.read_bytes()
    for operation in ("same-size", "delete", "truncate"):
        _mutate_recorded_file(path, original, operation)
        stale = plan_resume(
            mocked_pipeline.output_dir,
            expected_run_fingerprint=mocked_pipeline.fingerprint,
        )
        assert stale.completed is False, operation
        assert stale.reusable_phases == tuple(range(expected_phase)), operation
        assert stale.restart_phase == expected_phase, operation
        path.write_bytes(original)
        restored = plan_resume(
            mocked_pipeline.output_dir,
            expected_run_fingerprint=mocked_pipeline.fingerprint,
        )
        assert restored.completed is True, operation


def _rebind_phase_artifact(
    output_dir: Path,
    *,
    phase: int,
    relative_path: str,
) -> None:
    marker_path = output_dir / PHASE_MARKER_FILENAMES[phase]
    marker = json.loads(marker_path.read_text())
    for index, artifact in enumerate(marker["artifacts"]):
        if artifact["relative_path"] != relative_path:
            continue
        rebound = build_artifact_identity(
            output_dir / relative_path,
            root=output_dir,
            schema=artifact["schema"],
        )
        marker["artifacts"][index] = asdict(rebound)
        atomic_write_json(marker_path, marker)
        return
    raise AssertionError(f"phase marker does not record {relative_path}")


@pytest.mark.parametrize(
    ("mutation", "reason"),
    [
        ("wrong-id", "ablation event ID differs from the run identity"),
        ("wrong-digest", "ablation contract SHA256 mismatch"),
        ("noncanonical", "canonical form"),
        ("negative", "must be nonnegative"),
        ("boolean", "must be an integer"),
    ],
)
def test_semantically_invalid_ablation_fragment_restarts_its_phase(
    mocked_pipeline: _MockPipeline,
    mutation: str,
    reason: str,
) -> None:
    mocked_pipeline.run()
    relative_path = "phase1/ablation_events.json"
    path = mocked_pipeline.output_dir / relative_path
    document = json.loads(path.read_bytes())
    counter = document["counters"]["phase1_seed_surface_export"]
    if mutation == "wrong-id":
        document["ablation_id"] = "A1"
    elif mutation == "wrong-digest":
        document["contract_sha256"] = "0" * 64
    elif mutation == "negative":
        counter["opportunities"] = -1
    elif mutation == "boolean":
        counter["opportunities"] = True
    elif mutation != "noncanonical":
        raise AssertionError(f"unknown mutation: {mutation}")
    content = json.dumps(
        document,
        indent=2 if mutation == "noncanonical" else None,
        separators=None if mutation == "noncanonical" else (",", ":"),
        sort_keys=mutation != "noncanonical",
    ).encode("utf-8")
    path.write_bytes(content)
    _rebind_phase_artifact(
        mocked_pipeline.output_dir,
        phase=1,
        relative_path=relative_path,
    )

    plan = plan_resume(
        mocked_pipeline.output_dir,
        expected_run_fingerprint=mocked_pipeline.fingerprint,
    )
    _, prefix_reason = run_state_module._validated_phase_prefix(
        mocked_pipeline.output_dir,
        mocked_pipeline.fingerprint,
    )

    assert plan.reusable_phases == (0,)
    assert plan.restart_phase == 1
    assert reason in (prefix_reason or "")


def test_decreasing_ablation_fragment_counters_restart_at_first_regression(
    mocked_pipeline: _MockPipeline,
) -> None:
    mocked_pipeline.identity["config"]["ablation_id"] = "A2"
    mocked_pipeline.fingerprint = compute_run_fingerprint(
        mocked_pipeline.identity
    )
    orchestrator.single_genome_flow(
        genome_path=mocked_pipeline.genome,
        output_dir=mocked_pipeline.output_dir,
        genome_id="demo",
        device="cpu",
        resume=True,
        ablation_id="A2",
    )
    phase1_events = AblationEvents(
        ablation_id=AblationID.A2,
        counters=AblationCounters.for_ablation(
            AblationID.A2,
            opportunities=3,
            interventions=2,
            changed=1,
        ),
    )
    relative_path = "phase1/ablation_events.json"
    (mocked_pipeline.output_dir / relative_path).write_bytes(
        phase1_events.to_bytes()
    )
    _rebind_phase_artifact(
        mocked_pipeline.output_dir,
        phase=1,
        relative_path=relative_path,
    )
    phase2_marker = mocked_pipeline.output_dir / PHASE_MARKER_FILENAMES[2]
    phase2 = json.loads(phase2_marker.read_text())
    phase2["dependency_sha256"] = marker_sha256(
        mocked_pipeline.output_dir,
        phase=1,
    )
    atomic_write_json(phase2_marker, phase2)

    plan = plan_resume(
        mocked_pipeline.output_dir,
        expected_run_fingerprint=mocked_pipeline.fingerprint,
    )
    _, prefix_reason = run_state_module._validated_phase_prefix(
        mocked_pipeline.output_dir,
        mocked_pipeline.fingerprint,
    )

    assert plan.reusable_phases == (0, 1)
    assert plan.restart_phase == 2
    assert "ablation counters decreased" in (prefix_reason or "")


def test_ablation_counters_cannot_appear_before_the_owner_phase(
    mocked_pipeline: _MockPipeline,
) -> None:
    mocked_pipeline.identity["config"]["ablation_id"] = "A5"
    mocked_pipeline.fingerprint = compute_run_fingerprint(
        mocked_pipeline.identity
    )
    orchestrator.single_genome_flow(
        genome_path=mocked_pipeline.genome,
        output_dir=mocked_pipeline.output_dir,
        genome_id="demo",
        device="cpu",
        resume=True,
        ablation_id="A5",
    )
    early_events = AblationEvents(
        ablation_id=AblationID.A5,
        counters=AblationCounters.for_ablation(
            AblationID.A5,
            opportunities=1,
            interventions=1,
            changed=0,
        ),
    )
    relative_path = "phase0/ablation_events.json"
    (mocked_pipeline.output_dir / relative_path).write_bytes(
        early_events.to_bytes()
    )
    _rebind_phase_artifact(
        mocked_pipeline.output_dir,
        phase=0,
        relative_path=relative_path,
    )

    records, reason = run_state_module._validated_phase_prefix(
        mocked_pipeline.output_dir,
        mocked_pipeline.fingerprint,
    )

    assert records == ()
    assert "nonzero before their owner phase" in (reason or "")


def test_ablation_counters_cannot_change_after_the_owner_phase(
    mocked_pipeline: _MockPipeline,
) -> None:
    mocked_pipeline.identity["config"]["ablation_id"] = "A2"
    mocked_pipeline.fingerprint = compute_run_fingerprint(
        mocked_pipeline.identity
    )
    orchestrator.single_genome_flow(
        genome_path=mocked_pipeline.genome,
        output_dir=mocked_pipeline.output_dir,
        genome_id="demo",
        device="cpu",
        resume=True,
        ablation_id="A2",
    )
    for phase, counts in (
        (1, (2, 1, 1)),
        (2, (3, 2, 1)),
    ):
        events = AblationEvents(
            ablation_id=AblationID.A2,
            counters=AblationCounters.for_ablation(
                AblationID.A2,
                opportunities=counts[0],
                interventions=counts[1],
                changed=counts[2],
            ),
        )
        relative_path = f"phase{phase}/ablation_events.json"
        (mocked_pipeline.output_dir / relative_path).write_bytes(
            events.to_bytes()
        )
        _rebind_phase_artifact(
            mocked_pipeline.output_dir,
            phase=phase,
            relative_path=relative_path,
        )
        if phase < 2:
            next_marker = mocked_pipeline.output_dir / PHASE_MARKER_FILENAMES[2]
            payload = json.loads(next_marker.read_text())
            payload["dependency_sha256"] = marker_sha256(
                mocked_pipeline.output_dir,
                phase=1,
            )
            atomic_write_json(next_marker, payload)

    records, reason = run_state_module._validated_phase_prefix(
        mocked_pipeline.output_dir,
        mocked_pipeline.fingerprint,
    )

    assert tuple(record.phase for record in records) == (0, 1)
    assert "changed after their owner phase" in (reason or "")


def test_root_ablation_events_must_equal_the_terminal_phase_fragment(
    mocked_pipeline: _MockPipeline,
) -> None:
    mocked_pipeline.identity["config"]["ablation_id"] = "A6"
    mocked_pipeline.fingerprint = compute_run_fingerprint(
        mocked_pipeline.identity
    )
    orchestrator.single_genome_flow(
        genome_path=mocked_pipeline.genome,
        output_dir=mocked_pipeline.output_dir,
        genome_id="demo",
        device="cpu",
        resume=True,
        ablation_id="A6",
    )
    root_events = AblationEvents(
        ablation_id=AblationID.A6,
        counters=AblationCounters.for_ablation(
            AblationID.A6,
            opportunities=1,
            interventions=1,
            changed=1,
        ),
    )
    (mocked_pipeline.output_dir / "ablation_events.json").write_bytes(
        root_events.to_bytes()
    )
    _rebind_phase_artifact(
        mocked_pipeline.output_dir,
        phase=3,
        relative_path="ablation_events.json",
    )

    records, reason = run_state_module._validated_phase_prefix(
        mocked_pipeline.output_dir,
        mocked_pipeline.fingerprint,
    )

    assert tuple(record.phase for record in records) == (0, 1, 2)
    assert "final ablation events differ" in (reason or "")


def test_fresh_run_publishes_artifact_backed_state_in_order(
    mocked_pipeline: _MockPipeline,
) -> None:
    result = mocked_pipeline.run()

    assert result["success"] is True
    assert mocked_pipeline.events == [
        "running",
        "phase0",
        "phase1",
        "phase2",
        "phase3",
        "success",
    ]
    assert load_run_state(mocked_pipeline.output_dir).status == "success"
    assert all(
        (mocked_pipeline.output_dir / marker).is_file()
        for marker in PHASE_MARKER_FILENAMES
    )


@pytest.mark.parametrize(
    ("artifact_class", "path_choices"),
    _REQUIRED_FINAL_CLASSES,
    ids=[item[0].replace(" ", "-") for item in _REQUIRED_FINAL_CLASSES],
)
def test_every_required_final_artifact_class_is_content_authenticated(
    mocked_pipeline: _MockPipeline,
    artifact_class: str,
    path_choices: tuple[str, ...],
) -> None:
    mocked_pipeline.run()
    recorded = {
        artifact.relative_path
        for artifact in load_run_state(mocked_pipeline.output_dir).artifacts
    }
    relative_path = next(
        (choice for choice in path_choices if choice in recorded),
        None,
    )
    assert relative_path is not None, artifact_class

    _assert_mutations_invalidate_at_phase(
        path=mocked_pipeline.output_dir / relative_path,
        expected_phase=3,
        mocked_pipeline=mocked_pipeline,
    )


@pytest.mark.parametrize(
    ("phase", "relative_path"),
    _REQUIRED_PHASE_ARTIFACTS,
    ids=[
        f"phase{phase}-{Path(path).name}"
        for phase, path in _REQUIRED_PHASE_ARTIFACTS
    ],
)
def test_every_required_phase_artifact_is_content_authenticated(
    mocked_pipeline: _MockPipeline,
    phase: int,
    relative_path: str,
) -> None:
    mocked_pipeline.run()

    _assert_mutations_invalidate_at_phase(
        path=mocked_pipeline.output_dir / relative_path,
        expected_phase=phase,
        mocked_pipeline=mocked_pipeline,
    )


def test_every_recorded_artifact_rejects_all_content_mutations(
    mocked_pipeline: _MockPipeline,
) -> None:
    mocked_pipeline.run()
    recorded: dict[str, int] = {}
    for phase, marker_name in enumerate(PHASE_MARKER_FILENAMES):
        record = run_state_module._load_phase_record(
            mocked_pipeline.output_dir / marker_name
        )
        for artifact in record.artifacts:
            recorded.setdefault(artifact.relative_path, phase)
    state = load_run_state(mocked_pipeline.output_dir)
    for artifact in state.artifacts:
        recorded.setdefault(artifact.relative_path, 3)

    assert {
        "demo_eves.fna",
        "gvclass_results.tsv",
        "host_signature_model.png",
        "phase3_synthesis/virosync_predictions.tsv",
        "phase3_synthesis/virosync_predictions_detailed.tsv",
        "virosync_predictions_detailed.tsv",
    }.issubset(recorded)
    for relative_path, expected_phase in sorted(recorded.items()):
        _assert_mutations_invalidate_at_phase(
            path=mocked_pipeline.output_dir / relative_path,
            expected_phase=expected_phase,
            mocked_pipeline=mocked_pipeline,
        )


def test_dynamic_masked_fasta_is_recorded_and_content_authenticated(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "result"
    masking_dir = output_dir / "phase0" / "masking"
    masking_dir.mkdir(parents=True)
    input_fasta = tmp_path / "input.fna"
    input_fasta.write_text(">scaffold\nACGT\n", encoding="utf-8")
    masked_fasta = masking_dir / "masked.fasta"
    masked_fasta.write_text(">scaffold\nACGT\n", encoding="utf-8")

    def _sha256(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    write_masking_status(
        MaskingResult(
            output_path=masked_fasta,
            repeat_regions=(),
            requested_backend=MaskingBackend.TRF,
            effective_backend=MaskingBackend.TRF,
            failure_policy=MaskingFailurePolicy.STRICT,
            status="success",
            legacy_adapter=False,
            backend_versions=(("trf", "4.10.0rc2"),),
            masked_bases=0,
            repeatmasker_species=None,
            repeatmasker_library=None,
            repeatmasker_library_sha256=None,
            configured_fallback_backend=None,
            fallback_backend=None,
            fallback_reason=None,
            input_sha256=_sha256(input_fasta),
            output_sha256=_sha256(masked_fasta),
        ),
        masking_dir,
    )
    artifacts = orchestrator._phase_artifacts(output_dir, 0)
    identity = next(
        artifact
        for artifact in artifacts
        if artifact.relative_path == "phase0/masking/masked.fasta"
    )
    original = masked_fasta.read_bytes()

    for operation in ("same-size", "delete", "truncate"):
        _mutate_recorded_file(masked_fasta, original, operation)
        assert not run_state_module.validate_artifact_identity(
            identity,
            root=output_dir,
        )
        masked_fasta.write_bytes(original)
        assert run_state_module.validate_artifact_identity(
            identity,
            root=output_dir,
        )


def test_phase_artifacts_support_relative_output_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    output_dir = Path("results/demo")
    proteome = output_dir / "phase0/proteome.fasta"
    proteome.parent.mkdir(parents=True)
    proteome.write_text(">gene_1\nMPEPTIDE\n", encoding="utf-8")
    input_fasta = Path("input.fna")
    input_fasta.write_text(">scaffold\nACGT\n", encoding="utf-8")
    masking_dir = output_dir / "phase0/masking"
    masked_fasta = masking_dir / "masked.fasta"
    masked_fasta.parent.mkdir(parents=True)
    masked_fasta.write_text(">scaffold\nACGT\n", encoding="utf-8")
    input_sha256 = hashlib.sha256(input_fasta.read_bytes()).hexdigest()
    masked_sha256 = hashlib.sha256(masked_fasta.read_bytes()).hexdigest()
    write_masking_status(
        MaskingResult(
            output_path=masked_fasta,
            repeat_regions=(),
            requested_backend=MaskingBackend.TRF,
            effective_backend=MaskingBackend.TRF,
            failure_policy=MaskingFailurePolicy.STRICT,
            status="success",
            legacy_adapter=False,
            backend_versions=(("trf", "4.10.0rc2"),),
            masked_bases=0,
            repeatmasker_species=None,
            repeatmasker_library=None,
            repeatmasker_library_sha256=None,
            configured_fallback_backend=None,
            fallback_backend=None,
            fallback_reason=None,
            input_sha256=input_sha256,
            output_sha256=masked_sha256,
        ),
        masking_dir,
    )

    artifacts = orchestrator._phase_artifacts(output_dir, 0)
    identities = {artifact.relative_path: artifact for artifact in artifacts}

    assert identities["phase0/proteome.fasta"].sha256 == hashlib.sha256(
        proteome.read_bytes()
    ).hexdigest()
    assert identities["phase0/masking/masked.fasta"].sha256 == masked_sha256


@pytest.mark.parametrize("duplicate_kind", ["canonical", "summary"])
def test_duplicate_final_outputs_must_agree(
    mocked_pipeline: _MockPipeline,
    duplicate_kind: str,
) -> None:
    mocked_pipeline.run()
    state = load_run_state(mocked_pipeline.output_dir)
    assert state.result is not None
    recorded = {artifact.relative_path for artifact in state.artifacts}

    if duplicate_kind == "canonical":
        source_relative = next(
            path
            for path in (
                "phase3_synthesis/virosync_predictions.tsv",
                "virosync_predictions.tsv",
            )
            if path in recorded
        )
        duplicate_relative = (
            "virosync_predictions.tsv"
            if source_relative.startswith("phase3_synthesis/")
            else "phase3_synthesis/virosync_predictions.tsv"
        )
        content = (
            mocked_pipeline.output_dir / source_relative
        ).read_text().replace("0.91", "0.92")
        schema = "canonical-predictions-v4"
        expected_error = "duplicate canonical prediction tables disagree"
    else:
        source_relative = next(
            path
            for path in (
                "phase3_synthesis/virosync_summary.json",
                "virosync_summary.json",
            )
            if path in recorded
        )
        duplicate_relative = (
            "virosync_summary.json"
            if source_relative.startswith("phase3_synthesis/")
            else "phase3_synthesis/virosync_summary.json"
        )
        payload = json.loads(
            (mocked_pipeline.output_dir / source_relative).read_text()
        )
        payload["per_scaffold"] = {"different": {}}
        content = json.dumps(payload, sort_keys=True) + "\n"
        schema = "virosync-summary-v3"
        expected_error = "duplicate ViroSync summaries disagree"

    duplicate = mocked_pipeline.output_dir / duplicate_relative
    duplicate.parent.mkdir(parents=True, exist_ok=True)
    duplicate.write_text(content)
    duplicate_identity = build_artifact_identity(
        duplicate,
        root=mocked_pipeline.output_dir,
        schema=schema,
    )

    with pytest.raises(ValueError, match=expected_error):
        run_state_module._validate_success_artifacts(
            mocked_pipeline.output_dir,
            (*state.artifacts, duplicate_identity),
            state.result,
            run_fingerprint=state.run_fingerprint,
            identities=state.identities,
        )


def test_exception_retains_prefix_and_retry_reuses_only_verified_phases(
    mocked_pipeline: _MockPipeline,
) -> None:
    mocked_pipeline.fail_phase3_once = True
    with pytest.raises(RuntimeError, match="injected Phase 3 failure"):
        mocked_pipeline.run()

    state = load_run_state(mocked_pipeline.output_dir)
    plan = plan_resume(
        mocked_pipeline.output_dir,
        expected_run_fingerprint=mocked_pipeline.fingerprint,
    )
    assert state.status == "failed"
    assert plan.reusable_phases == (0, 1, 2)
    assert plan.restart_phase == 3
    assert not (mocked_pipeline.output_dir / PHASE_MARKER_FILENAMES[3]).exists()

    result = mocked_pipeline.run()

    assert result["success"] is True
    assert mocked_pipeline.calls == {
        "phase0": 1,
        "phase1": 2,
        "phase2": 2,
        "phase3": 2,
    }
    assert mocked_pipeline.resume_flags["phase1"] == [False, True]
    assert mocked_pipeline.resume_flags["phase2"] == [False, True]
    assert mocked_pipeline.resume_flags["phase3"] == [False, False]
    assert mocked_pipeline.phase3_boundaries[0].seed_id == "lossless-seed"
    assert mocked_pipeline.phase3_boundaries[0].hallmark_genes == ["MCP"]


def test_unmarked_partial_phase_files_are_recomputed(
    mocked_pipeline: _MockPipeline,
) -> None:
    phase0 = mocked_pipeline.output_dir / "phase0"
    phase0.mkdir(parents=True)
    stale = phase0 / "proteome.fasta"
    stale.write_text(">stale\nX\n")

    mocked_pipeline.run()

    assert mocked_pipeline.calls["phase0"] == 1
    assert stale.read_text() == ">p1\nM\n"
    assert load_run_state(mocked_pipeline.output_dir).status == "success"


@pytest.mark.parametrize(
    ("phase", "relative_path", "partial_content"),
    [
        (0, "phase0/genes.gff", "##gff-version 3\n"),
        (1, "phase1/resume_state.json", '{"schema_version": 1}\n'),
        (2, "phase2/refined_boundaries.bed", "# scaffold start end\n"),
        (3, "phase3/state.tsv", "state\n"),
    ],
)
def test_header_only_unmarked_partial_is_never_reused(
    mocked_pipeline: _MockPipeline,
    phase: int,
    relative_path: str,
    partial_content: str,
) -> None:
    partial = mocked_pipeline.output_dir / relative_path
    partial.parent.mkdir(parents=True, exist_ok=True)
    partial.write_text(partial_content)

    mocked_pipeline.run()

    assert mocked_pipeline.calls[f"phase{phase}"] == 1
    assert partial.read_text() != partial_content
    assert plan_resume(
        mocked_pipeline.output_dir,
        expected_run_fingerprint=mocked_pipeline.fingerprint,
    ).completed


def test_state_json_writer_and_reader_enforce_one_size_limit(
    tmp_path: Path,
) -> None:
    limit = run_state_module._MAX_STATE_BYTES
    destination = tmp_path / "state.json"
    original = b'{"sentinel": true}\n'
    destination.write_bytes(original)

    with pytest.raises(ValueError, match="exceeds.*reader limit"):
        atomic_write_json(destination, {"payload": "x" * limit})
    assert destination.read_bytes() == original

    output_dir = tmp_path / "oversize-read"
    output_dir.mkdir()
    oversized = output_dir / RUN_STATE_FILENAME
    with oversized.open("wb") as handle:
        handle.truncate(limit + 1)
    with pytest.raises(ValueError, match="state JSON has invalid size"):
        load_run_state(output_dir)


@pytest.mark.parametrize(
    ("phase", "artifact", "operation"),
    [
        (0, "phase0/proteome.fasta", "mutate"),
        (1, "phase1/resume_state.json", "mutate"),
        (2, "phase2/refined_state.json", "delete"),
    ],
)
def test_stale_phase_artifact_restarts_that_phase_and_downstream(
    mocked_pipeline: _MockPipeline,
    phase: int,
    artifact: str,
    operation: str,
) -> None:
    mocked_pipeline.fail_phase3_once = True
    with pytest.raises(RuntimeError):
        mocked_pipeline.run()

    artifact_path = mocked_pipeline.output_dir / artifact
    if operation == "delete":
        artifact_path.unlink()
    else:
        artifact_path.write_text("mutated\n")

    stale_plan = plan_resume(
        mocked_pipeline.output_dir,
        expected_run_fingerprint=mocked_pipeline.fingerprint,
    )
    assert stale_plan.reusable_phases == tuple(range(phase))
    assert stale_plan.restart_phase == phase

    result = mocked_pipeline.run()

    assert result["success"] is True
    assert mocked_pipeline.calls["phase0"] == (2 if phase == 0 else 1)
    assert mocked_pipeline.resume_flags["phase1"][-1] is (phase > 1)
    assert mocked_pipeline.resume_flags["phase2"][-1] is (phase > 2)
    assert mocked_pipeline.resume_flags["phase3"][-1] is False


@pytest.mark.parametrize("terminal_phase", [1, 2])
def test_terminal_zero_publishes_exact_prefix_and_zero_success(
    mocked_pipeline: _MockPipeline,
    terminal_phase: int,
) -> None:
    mocked_pipeline.terminal_phase = terminal_phase

    result = mocked_pipeline.run()

    state = load_run_state(mocked_pipeline.output_dir)
    assert result["success"] is True
    assert state.status == "success"
    assert state.result is not None
    assert state.result["terminal_phase"] == terminal_phase
    assert state.result["canonical_rows"] == 0
    assert state.result["detailed_rows"] == 0
    assert state.result["accepted_bp"] == 0
    assert not any(state.result["class_counts"].values())
    assert not any(state.result["tier_counts"].values())
    assert {
        path.name
        for path in mocked_pipeline.output_dir.glob("phase*.complete.json")
    } == set(PHASE_MARKER_FILENAMES[: terminal_phase + 1])


def test_terminal_zero_rejects_dangling_downstream_marker_symlink(
    mocked_pipeline: _MockPipeline,
) -> None:
    mocked_pipeline.terminal_phase = 1
    mocked_pipeline.run()
    downstream = mocked_pipeline.output_dir / PHASE_MARKER_FILENAMES[2]
    downstream.symlink_to("missing-phase2-marker.json")

    stale = plan_resume(
        mocked_pipeline.output_dir,
        expected_run_fingerprint=mocked_pipeline.fingerprint,
    )

    assert stale.completed is False
    assert stale.reusable_phases == (0,)
    assert stale.restart_phase == 1


@pytest.mark.parametrize("terminal_phase", [1, 2])
def test_terminal_marker_promotes_after_crash_without_phase_rerun(
    mocked_pipeline: _MockPipeline,
    monkeypatch: pytest.MonkeyPatch,
    terminal_phase: int,
) -> None:
    mocked_pipeline.terminal_phase = terminal_phase
    original_success = orchestrator.publish_run_success

    def crash_before_success(*args, **kwargs):
        raise RuntimeError("injected crash after terminal marker")

    monkeypatch.setattr(orchestrator, "publish_run_success", crash_before_success)
    with pytest.raises(RuntimeError, match="after terminal marker"):
        mocked_pipeline.run()

    assert (
        mocked_pipeline.output_dir / PHASE_MARKER_FILENAMES[terminal_phase]
    ).is_file()
    calls_before = dict(mocked_pipeline.calls)
    monkeypatch.setattr(orchestrator, "publish_run_success", original_success)

    resumed = mocked_pipeline.run()

    assert resumed["success"] is True
    assert mocked_pipeline.calls == calls_before
    state = load_run_state(mocked_pipeline.output_dir)
    assert state.status == "success"
    assert state.result["terminal_phase"] == terminal_phase


def test_phase3_marker_survives_failed_success_publish_and_promotes_without_rerun(
    mocked_pipeline: _MockPipeline,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original_success = orchestrator.publish_run_success

    def crash_before_success(*args, **kwargs):
        raise RuntimeError("injected crash after Phase 3 marker")

    monkeypatch.setattr(orchestrator, "publish_run_success", crash_before_success)
    with pytest.raises(RuntimeError, match="after Phase 3 marker"):
        mocked_pipeline.run()

    assert (mocked_pipeline.output_dir / PHASE_MARKER_FILENAMES[3]).is_file()
    calls_before = dict(mocked_pipeline.calls)
    monkeypatch.setattr(orchestrator, "publish_run_success", original_success)

    result = mocked_pipeline.run()

    assert result["success"] is True
    assert mocked_pipeline.calls == calls_before
    assert load_run_state(mocked_pipeline.output_dir).status == "success"


def test_completed_run_phase3_artifact_mutation_reuses_phase0_to_phase2(
    mocked_pipeline: _MockPipeline,
) -> None:
    mocked_pipeline.run()
    calls_after_fresh = dict(mocked_pipeline.calls)
    (
        mocked_pipeline.output_dir
        / "phase3_synthesis"
        / "virosync_predictions.tsv"
    ).write_text("mutated\n")

    stale = plan_resume(
        mocked_pipeline.output_dir,
        expected_run_fingerprint=mocked_pipeline.fingerprint,
    )
    assert stale.reusable_phases == (0, 1, 2)
    assert stale.restart_phase == 3

    result = mocked_pipeline.run()

    assert result["success"] is True
    assert mocked_pipeline.calls["phase0"] == calls_after_fresh["phase0"]
    assert mocked_pipeline.calls["phase1"] == calls_after_fresh["phase1"] + 1
    assert mocked_pipeline.calls["phase2"] == calls_after_fresh["phase2"] + 1
    assert mocked_pipeline.calls["phase3"] == calls_after_fresh["phase3"] + 1


def test_unowned_phase3_diagnostic_mutation_does_not_invalidate_success(
    mocked_pipeline: _MockPipeline,
) -> None:
    mocked_pipeline.run()
    (mocked_pipeline.output_dir / "phase3" / "state.tsv").write_text(
        "mutated diagnostic\n"
    )

    plan = plan_resume(
        mocked_pipeline.output_dir,
        expected_run_fingerprint=mocked_pipeline.fingerprint,
    )

    assert plan.completed is True
    assert plan.reusable_phases == (0, 1, 2, 3)


def test_fresh_and_completed_resume_return_identical_summary_counts(
    mocked_pipeline: _MockPipeline,
) -> None:
    fresh = mocked_pipeline.run()
    calls_after_fresh = dict(mocked_pipeline.calls)

    resumed = mocked_pipeline.run()

    summary_fields = set(_empty_prediction_summary())
    assert {field: fresh[field] for field in summary_fields} == {
        field: resumed[field] for field in summary_fields
    }
    assert fresh["accepted"] == 1
    assert fresh["predictions"] == 1
    assert fresh["accepted_bp"] == 4
    assert mocked_pipeline.calls == calls_after_fresh
    assert Path(resumed["output_files"]["run_state"]).name == RUN_STATE_FILENAME


def test_completed_resume_normalizes_legacy_missing_promoted_low_count(
    mocked_pipeline: _MockPipeline,
) -> None:
    mocked_pipeline.run()
    calls_after_fresh = dict(mocked_pipeline.calls)
    state_path = mocked_pipeline.output_dir / RUN_STATE_FILENAME
    document = json.loads(state_path.read_text())
    document["result"].pop("promoted_low_rows")
    atomic_write_json(state_path, document)

    resumed = mocked_pipeline.run()

    assert resumed["success"] is True
    assert mocked_pipeline.calls == calls_after_fresh


@pytest.mark.parametrize("legacy_class", ["VP", "PLV", "MIXED"])
def test_legacy_class_tokens_resume_without_recompute(
    mocked_pipeline: _MockPipeline,
    legacy_class: str,
) -> None:
    mocked_pipeline.effective_eve_class = legacy_class

    fresh = mocked_pipeline.run()
    calls_after_fresh = dict(mocked_pipeline.calls)
    resumed = mocked_pipeline.run()

    assert fresh["success"] is True
    assert resumed["success"] is True
    assert mocked_pipeline.calls == calls_after_fresh
    plan = plan_resume(
        mocked_pipeline.output_dir,
        expected_run_fingerprint=mocked_pipeline.fingerprint,
    )
    assert plan.completed is True


def _result_with_promoted_low(value: object = 1) -> dict[str, object]:
    return {
        "terminal_phase": None,
        "canonical_rows": 2,
        "detailed_rows": 2,
        "accepted_bp": 8,
        "class_counts": {
            eve_class: 2 if eve_class == "NCLDV" else 0
            for eve_class in EFFECTIVE_EVE_CLASSES
        },
        "tier_counts": {"HIGH": 0, "MEDIUM": 0, "LOW": 2},
        "promoted_low_rows": value,
        "benchmark_eligible": True,
    }


@pytest.mark.parametrize("value", [-1, True, 3])
def test_result_identity_rejects_invalid_promoted_low_count(value: object) -> None:
    with pytest.raises((TypeError, ValueError), match="promoted_low_rows"):
        run_state_module._validated_result(_result_with_promoted_low(value))


def test_legacy_result_identity_defaults_promoted_count_to_low_rows() -> None:
    result = _result_with_promoted_low()
    result.pop("promoted_low_rows")

    normalized = run_state_module._validated_result(result)

    assert normalized["promoted_low_rows"] == 2


@pytest.mark.parametrize(
    "field",
    ["low_confidence", "promoted_low_confidence"],
)
def test_success_validation_rejects_independent_low_summary_mismatch(
    mocked_pipeline: _MockPipeline,
    field: str,
) -> None:
    mocked_pipeline.run()
    state = load_run_state(mocked_pipeline.output_dir)
    summary_artifact = next(
        artifact
        for artifact in state.artifacts
        if artifact.relative_path.endswith("virosync_summary.json")
    )
    summary_path = mocked_pipeline.output_dir / summary_artifact.relative_path
    summary = json.loads(summary_path.read_text())
    summary["statistics"][field] = 1
    summary_path.write_text(json.dumps(summary) + "\n")
    rebound_artifacts = tuple(
        build_artifact_identity(
            summary_path,
            root=mocked_pipeline.output_dir,
            schema=summary_artifact.schema,
        )
        if artifact.relative_path == summary_artifact.relative_path
        else artifact
        for artifact in state.artifacts
    )

    with pytest.raises(ValueError, match=f"summary {field} disagrees"):
        run_state_module._validate_success_artifacts(
            mocked_pipeline.output_dir,
            rebound_artifacts,
            state.result,
            run_fingerprint=state.run_fingerprint,
            identities=state.identities,
        )


def test_candidates_without_acceptance_publish_and_resume_exactly(
    mocked_pipeline: _MockPipeline,
) -> None:
    mocked_pipeline.candidate_only = True

    fresh = mocked_pipeline.run()
    resumed = mocked_pipeline.run()
    state = load_run_state(mocked_pipeline.output_dir)

    assert fresh["success"] is True
    assert fresh["predictions"] == resumed["predictions"] == 1
    assert fresh["accepted"] == resumed["accepted"] == 0
    assert (
        fresh["quality_gate_dropped"]
        == resumed["quality_gate_dropped"]
        == 1
    )
    assert state.result is not None
    assert state.result["canonical_rows"] == 0
    assert state.result["detailed_rows"] == 1
    assert state.result["accepted_bp"] == 0
    assert not any(state.result["class_counts"].values())
    assert not any(state.result["tier_counts"].values())
    assert plan_resume(
        mocked_pipeline.output_dir,
        expected_run_fingerprint=mocked_pipeline.fingerprint,
    ).completed


@pytest.mark.parametrize("terminal_phase", [None, 1, 2])
def test_output_files_are_identical_for_fresh_resumed_and_recovered_success(
    mocked_pipeline: _MockPipeline,
    terminal_phase: int | None,
) -> None:
    mocked_pipeline.terminal_phase = terminal_phase
    fresh = mocked_pipeline.run()
    resumed = mocked_pipeline.run()

    orchestrator.publish_run_started(
        mocked_pipeline.output_dir,
        run_fingerprint=mocked_pipeline.fingerprint,
        identities=mocked_pipeline.identity,
        preserve_success_artifacts=True,
    )
    recovered = mocked_pipeline.run()

    assert fresh["output_files"] == resumed["output_files"] == recovered["output_files"]
    assert set(fresh["output_files"]) == {"run_state", "artifacts"}
    state = load_run_state(mocked_pipeline.output_dir)
    assert fresh["output_files"]["artifacts"] == {
        artifact.relative_path: str(
            mocked_pipeline.output_dir / artifact.relative_path
        )
        for artifact in sorted(
            state.artifacts,
            key=lambda item: item.relative_path,
        )
    }
