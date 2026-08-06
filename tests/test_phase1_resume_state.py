from __future__ import annotations

import copy
import json
import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

from virosync.orchestration._flows.single_genome import phase1 as phase1_module
from virosync.orchestration._flows.single_genome.phase1 import (
    _run_phase1_subflow,
)
from virosync.orchestration._flows.single_genome.phase1_state import (
    PHASE1_STATE_ARTIFACT_TYPE,
    PHASE1_STATE_FILENAME,
    PHASE1_STATE_SCHEMA_VERSION,
    Phase1StateError,
    load_phase1_state,
    phase1_state_from_document,
    phase1_state_to_document,
    write_phase1_state,
)
from virosync.pipeline.host_signatures import HostSignatureModel
from virosync.pipeline.phase1.frameshift_screening import (
    ANNOTATION_CLASS,
    FrameshiftHit,
    rescued_protein_id,
)
from virosync.pipeline.phase1.hhg_seeding import Anchor, HMMHit
from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit
from virosync.pipeline.phase1.seed_merger import MergedSeed


def _complete_marker() -> ValidatedMarkerHit:
    return ValidatedMarkerHit(
        query_porf="scaffold/alpha_17",
        scaffold="scaffold/alpha",
        start=101,
        end=999,
        strand="-",
        hmm_target="GVOGm0003",
        hmm_score=123.456789012345,
        hmm_evalue=1.23456789012345e-87,
        validation_status="validated",
        top10_prefixes="EUK__,NCLDV__,PLV__",
        best_hit_target="NCLDV__exact|protein-1",
        best_hit_pident=33.333333333333336,
        best_hit_bits=456.7890123456789,
        has_ncldv=1,
        has_mirus=0,
        has_plv=1,
        has_vp=1,
        has_viral=1,
        top10_targets="host-1,virus-1",
        top10_pidents="99.123456789,33.333333333333336",
        top10_bitscores="500.000000001,456.7890123456789",
        top10_evalues="0.0,1.23456789012345e-87",
        taxonomy_substring_counts="Amoebozoa:7.123456789",
        taxonomy_raw_counts="Amoebozoa:8",
    )


def _complete_seed() -> MergedSeed:
    anchor = Anchor(
        porf_id="scaffold/alpha_17",
        scaffold="scaffold/alpha",
        start=101,
        end=999,
        strand="-",
        hallmark_gene="GVOGm0003",
        score=123.456789012345,
        evalue=1.23456789012345e-87,
    )
    second_anchor = Anchor(
        porf_id="scaffold/alpha_18",
        scaffold="scaffold/alpha",
        start=1001,
        end=1550,
        strand="+",
        hallmark_gene="GVOGm0054",
        score=87.6543210987654,
        evalue=9.87654321098765e-41,
    )
    return MergedSeed(
        scaffold="scaffold/alpha",
        start=51,
        end=1601,
        seed_id="seed_0_scaffold/alpha_51",
        sources=["hhg", "marker_validation"],
        hhg_score=123.456789012345,
        novelty_score=0.123456789012345,
        compositional_score=0.234567890123456,
        mean_kfd=0.345678901234567,
        mean_composite=0.456789012345678,
        max_kfd=0.567890123456789,
        max_composite=0.678901234567891,
        gc_deviation=-0.123456789012345,
        cub_deviation=0.012345678901234,
        n_windows=7,
        cluster_ids=[3, 9],
        anchors=[anchor],
        hhg_anchors=[second_anchor],
        priority=0.789012345678901,
        confidence="high",
        score=0.890123456789012,
        predicted_family="MIXED",
        region_classification_ncldv_markers=4,
        region_classification_vp_plv_markers=2,
        region_classification_mirus_markers=1,
        host_trim_original_start=41,
        host_trim_original_end=1611,
        host_trimmed_start=51,
        host_trimmed_end=1601,
        host_trim_reason="host-taxonomy",
        host_trim_common_euk_taxonomy="Eukaryota;Amoebozoa",
    )


def _complete_host_model() -> HostSignatureModel:
    return HostSignatureModel(
        token_weights={"amoebozoa": 7.123456789012345},
        token_counts={"amoebozoa": 8},
        token_bits={"amoebozoa": [99.12345678901234, 88.98765432109876]},
        max_weight=7.123456789012345,
        min_token_length=4,
        host_prefixes=["EUK__", "HOST__"],
        weight_mode="bitscore",
    )


def _complete_deviation_summary() -> dict[str, object]:
    return {
        "enabled": True,
        "markers_total": 5,
        "markers_seedable": 2,
        "baseline": {
            "fractions": [0.123456789012345, 0.987654321098765],
            "nested": {"count": 17, "usable": True, "label": None},
        },
        "report_path": "phase1/marker_validation/host_taxonomy_deviation.tsv",
    }


def _write_complete_state(path: Path) -> None:
    write_phase1_state(
        path,
        validated_markers=[_complete_marker()],
        merged_seeds=[_complete_seed()],
        host_signature_model=_complete_host_model(),
        host_signatures={"Amoebozoa", "Eukaryota"},
        host_deviation_summary=_complete_deviation_summary(),
    )


def _phase1_kwargs(tmp_path: Path, output_dir: Path) -> dict[str, object]:
    return {
        "masked_path": tmp_path / "masked.fna",
        "proteome_path": tmp_path / "proteome.faa",
        "repeat_regions": [],
        "output_dir": output_dir,
        "genome_id": "genome-a",
        "hmm_database": None,
        "hmm_allowlist": None,
        "hmm_chunk_size": None,
        "frameshift_screening_enabled": False,
        "marker_faa_db": None,
        "marker_faa_dir": None,
        "marker_db": None,
        "faa_dir": None,
        "gene_taxonomy_faa_db": None,
        "taxonomy_labels_file": None,
        "host_prefixes": ["EUK__"],
        "host_label": "EUK__",
        "taxonomy_weight_mode": "rank",
        "host_taxonomy_deviation_enabled": False,
        "host_taxonomy_deviation_allow_seeds": False,
        "host_taxonomy_deviation_min_token_len": 3,
        "host_taxonomy_deviation_min_tokens": 2,
        "host_taxonomy_deviation_overlap_threshold": 0.5,
        "host_taxonomy_deviation_max_pident": 80.0,
        "host_taxonomy_deviation_max_hits": 10,
        "host_taxonomy_deviation_window_bp": 10000,
        "host_taxonomy_deviation_window_count": 10,
        "host_taxonomy_deviation_window_seed": 17,
        "host_taxonomy_deviation_window_min_markers": 2,
        "host_taxonomy_deviation_seed_window_bp": 10000,
        "host_taxonomy_deviation_seed_min_markers": 2,
        "marker_validation_top_k": 10,
        "novel_marker_min_score": 30.0,
        "novel_marker_min_coverage": 0.5,
        "novel_marker_require_cluster": True,
        "initial_window_bp": 20000,
        "initial_window_genes": 10,
        "min_markers_initial": 1,
        "extension_kb": 10,
        "merge_distance": 10000,
        "boundary_host_signature_min_token_len": 3,
        "rebuild_db": False,
        "assembly_mode": "default",
        "extended_output": False,
        "resume": True,
        "threads": 1,
        "search_backend": "diamond",
        "logger": logging.getLogger("test-phase1-resume-state"),
        "config_fingerprint": "f" * 64,
        "resume_authorized": True,
    }


def test_phase1_state_file_round_trip_preserves_every_consumed_field(
    tmp_path: Path,
) -> None:
    state_path = tmp_path / "phase1" / PHASE1_STATE_FILENAME

    _write_complete_state(state_path)
    loaded = load_phase1_state(state_path)

    assert loaded.validated_markers == [_complete_marker()]
    assert loaded.merged_seeds == [_complete_seed()]
    assert loaded.host_signature_model == _complete_host_model()
    assert loaded.host_signatures == {"Amoebozoa", "Eukaryota"}
    assert loaded.host_deviation_summary == _complete_deviation_summary()
    payload = json.loads(state_path.read_text())
    assert payload["schema_version"] == PHASE1_STATE_SCHEMA_VERSION
    assert payload["artifact_type"] == PHASE1_STATE_ARTIFACT_TYPE
    assert payload["validated_markers"][0]["hmm_score"] == 123.456789012345
    assert payload["validated_markers"][0]["has_plv"] == 1
    assert payload["validated_markers"][0]["has_vp"] == 1


def test_phase1_state_rejects_schema_drift_and_lossy_values() -> None:
    document = phase1_state_to_document(
        validated_markers=[_complete_marker()],
        merged_seeds=[_complete_seed()],
        host_signature_model=_complete_host_model(),
        host_signatures={"Amoebozoa", "Eukaryota"},
        host_deviation_summary=_complete_deviation_summary(),
    )

    unknown_schema = copy.deepcopy(document)
    unknown_schema["schema_version"] = 2
    with pytest.raises(Phase1StateError, match="unsupported.*schema_version"):
        phase1_state_from_document(unknown_schema)

    missing_marker_field = copy.deepcopy(document)
    del missing_marker_field["validated_markers"][0]["hmm_evalue"]
    with pytest.raises(Phase1StateError, match="missing=.*hmm_evalue"):
        phase1_state_from_document(missing_marker_field)

    extra_seed_field = copy.deepcopy(document)
    extra_seed_field["merged_seeds"][0]["python_type"] = "arbitrary.Class"
    with pytest.raises(Phase1StateError, match="extra=.*python_type"):
        phase1_state_from_document(extra_seed_field)

    rounded_boolean = copy.deepcopy(document)
    rounded_boolean["validated_markers"][0]["has_plv"] = True
    with pytest.raises(Phase1StateError, match="has_plv must be an integer"):
        phase1_state_from_document(rounded_boolean)

    nonfinite_model = copy.deepcopy(document)
    nonfinite_model["host_signature_model"]["token_bits"]["amoebozoa"][0] = float(
        "nan"
    )
    with pytest.raises(Phase1StateError, match="token_bits.*must be finite"):
        phase1_state_from_document(nonfinite_model)

    unsorted_signatures = copy.deepcopy(document)
    unsorted_signatures["host_signatures"] = ["Eukaryota", "Amoebozoa"]
    with pytest.raises(Phase1StateError, match="sorted and unique"):
        phase1_state_from_document(unsorted_signatures)

    open_deviation_schema = copy.deepcopy(document)
    open_deviation_schema["host_deviation_summary"]["unexpected"] = 1
    with pytest.raises(Phase1StateError, match="extra=.*unexpected"):
        phase1_state_from_document(open_deviation_schema)


def test_phase1_state_file_rejects_duplicate_keys_and_nonstandard_numbers(
    tmp_path: Path,
) -> None:
    duplicate_path = tmp_path / "duplicate.json"
    duplicate_path.write_text(
        '{"artifact_type":"virosync.phase1.resume_state",'
        '"schema_version":1,"schema_version":1,'
        '"validated_markers":[],"merged_seeds":[],'
        '"host_signature_model":{},"host_signatures":[],'
        '"host_deviation_summary":null}'
    )
    with pytest.raises(Phase1StateError, match="duplicate JSON key"):
        load_phase1_state(duplicate_path)

    nonfinite_path = tmp_path / "nonfinite.json"
    nonfinite_path.write_text(
        '{"artifact_type":"virosync.phase1.resume_state",'
        '"schema_version":1,"validated_markers":NaN,"merged_seeds":[],'
        '"host_signature_model":{},"host_signatures":[],'
        '"host_deviation_summary":null}'
    )
    with pytest.raises(Phase1StateError, match="non-finite JSON number"):
        load_phase1_state(nonfinite_path)


def test_phase1_state_failed_validation_preserves_existing_file(
    tmp_path: Path,
) -> None:
    state_path = tmp_path / PHASE1_STATE_FILENAME
    state_path.write_text("previous valid state\n")
    marker = _complete_marker()
    marker.hmm_score = float("nan")

    with pytest.raises(Phase1StateError, match="hmm_score must be finite"):
        write_phase1_state(
            state_path,
            validated_markers=[marker],
            merged_seeds=[_complete_seed()],
            host_signature_model=_complete_host_model(),
            host_signatures=set(),
            host_deviation_summary=None,
        )

    assert state_path.read_text() == "previous valid state\n"


def test_phase1_state_replace_failure_preserves_existing_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    state_path = tmp_path / PHASE1_STATE_FILENAME
    state_path.write_text("previous valid state\n")
    original_rename = Path.rename

    def fail_target_replace(source: Path, target: Path) -> Path:
        if Path(target) == state_path:
            raise OSError("injected replace failure")
        return original_rename(source, target)

    monkeypatch.setattr(Path, "rename", fail_target_replace)

    with pytest.raises(OSError, match="injected replace failure"):
        _write_complete_state(state_path)

    assert state_path.read_text() == "previous valid state\n"
    assert not list(tmp_path.glob(f".tmp_{PHASE1_STATE_FILENAME}_*.tmp"))


def test_authenticated_phase1_resume_loads_only_exact_state(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "output"
    state_path = output_dir / "phase1" / PHASE1_STATE_FILENAME
    _write_complete_state(state_path)

    result = _run_phase1_subflow(**_phase1_kwargs(tmp_path, output_dir))

    assert result["validated_markers"] == [_complete_marker()]
    assert result["merged_seeds"] == [_complete_seed()]
    assert result["host_signature_model"] == _complete_host_model()
    assert result["host_signatures"] == {"Amoebozoa", "Eukaryota"}
    assert result["host_deviation_summary"] == _complete_deviation_summary()
    assert not (output_dir / "phase1" / "marker_validation").exists()
    assert not (output_dir / "phase1" / "region_assembly").exists()


def test_enabled_phase1_resume_requires_confirmed_frameshift_faa(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "output"
    state_path = output_dir / "phase1" / PHASE1_STATE_FILENAME
    _write_complete_state(state_path)
    kwargs = _phase1_kwargs(tmp_path, output_dir)
    kwargs["frameshift_screening_enabled"] = True

    with pytest.raises(ValueError, match="confirmed frameshift protein FAA"):
        _run_phase1_subflow(**kwargs)

    confirmed_faa = (
        output_dir
        / "phase1"
        / "frameshift_screening"
        / "confirmed_frameshift_proteins.faa"
    )
    confirmed_faa.parent.mkdir(parents=True)
    confirmed_faa.write_text("")
    with pytest.raises(ValueError, match="confirmed frameshift marker table"):
        _run_phase1_subflow(**kwargs)


@pytest.mark.parametrize(
    ("enabled", "expected_prefix"),
    [
        (False, ["hmm"]),
        (True, ["frameshift", "hmm"]),
    ],
)
def test_frameshift_screening_runs_before_zero_hmm_hit_return_only_when_enabled(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    enabled: bool,
    expected_prefix: list[str],
) -> None:
    output_dir = tmp_path / "output"
    kwargs = _phase1_kwargs(tmp_path, output_dir)
    kwargs["resume"] = False
    kwargs["resume_authorized"] = False
    kwargs["frameshift_screening_enabled"] = enabled
    hmm_database = tmp_path / "markers.hmm"
    marker_db = tmp_path / "markers.dmnd"
    hmm_database.write_text("HMM\n")
    marker_db.write_text("DB\n")
    kwargs["hmm_database"] = hmm_database
    kwargs["marker_db"] = marker_db
    events: list[str] = []

    def fake_call_task(task, **task_kwargs):
        if task is phase1_module.frameshift_screening_task:
            assert task_kwargs == {
                "masked_fasta": kwargs["masked_path"],
                "hmm_database": hmm_database,
                "output_dir": output_dir / "phase1" / "frameshift_screening",
                "threads": kwargs["threads"],
            }
            events.append("frameshift")
            return []
        if task is phase1_module.hhg_seeding_task:
            events.append("hmm")
            return [], []
        if task is phase1_module.generate_outputs_task:
            events.append("outputs")
            return {}
        raise AssertionError(f"unexpected task: {task}")

    monkeypatch.setattr(phase1_module, "call_task", fake_call_task)
    monkeypatch.setattr(phase1_module, "_generate_required_reports", lambda **kwargs: {})
    monkeypatch.setattr(phase1_module, "_write_empty_run_log", lambda **kwargs: None)

    result = _run_phase1_subflow(**kwargs)

    assert result["success"] is True
    assert events[: len(expected_prefix)] == expected_prefix
    assert events[-1] == "outputs"


def test_pfam_contradiction_stops_before_marker_validation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_dir = tmp_path / "output"
    kwargs = _phase1_kwargs(tmp_path, output_dir)
    kwargs["resume"] = False
    kwargs["resume_authorized"] = False
    hmm_database = tmp_path / "models" / "combined.hmm"
    hmm_database.parent.mkdir()
    hmm_database.write_text("HMM\n")
    (hmm_database.parent / "pfam_virosync_screening.hmm").write_text("PFAM\n")
    marker_db = tmp_path / "markers.dmnd"
    marker_db.write_text("DB\n")
    kwargs["hmm_database"] = hmm_database
    kwargs["marker_db"] = marker_db
    hits = [
        HMMHit("protein-1", "VS000806", 100.0, 1e-20, 90.0, 1, 100),
        HMMHit("protein-1", "VS000369", 90.0, 1e-18, 80.0, 1, 100),
    ]
    events = []

    def fake_call_task(task, **task_kwargs):
        if task is phase1_module.hhg_seeding_task:
            events.append("hmm")
            return [], hits
        if task is phase1_module.pfam_arbitration_task:
            events.append("pfam")
            assert task_kwargs["hmm_hits"] == hits
            assert task_kwargs["proteins"] == {"protein-1"}
            return []
        if task is phase1_module.generate_outputs_task:
            events.append("outputs")
            return {}
        if task.__name__ == "marker_validation_task":
            raise AssertionError("contradicted hits reached marker validation")
        raise AssertionError(f"unexpected task: {task}")

    monkeypatch.setattr(phase1_module, "call_task", fake_call_task)
    monkeypatch.setattr(phase1_module, "_generate_required_reports", lambda **kwargs: {})
    monkeypatch.setattr(phase1_module, "_write_empty_run_log", lambda **kwargs: None)

    result = _run_phase1_subflow(**kwargs)

    assert result["success"] is True
    assert events == ["hmm", "pfam", "outputs"]


def test_missing_pfam_resource_skips_arbitration_and_preserves_marker_validation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_dir = tmp_path / "output"
    kwargs = _phase1_kwargs(tmp_path, output_dir)
    kwargs["resume"] = False
    kwargs["resume_authorized"] = False
    hmm_database = tmp_path / "models" / "combined.hmm"
    hmm_database.parent.mkdir()
    hmm_database.write_text("HMM\n")
    marker_db = tmp_path / "markers.dmnd"
    marker_db.write_text("DB\n")
    kwargs["hmm_database"] = hmm_database
    kwargs["marker_db"] = marker_db
    hits = [
        HMMHit("protein-1", "model-a", 100.0, 1e-20, 90.0, 1, 100),
        HMMHit("protein-1", "model-b", 90.0, 1e-18, 80.0, 1, 100),
    ]
    events = []

    def fake_call_task(task, **task_kwargs):
        if task is phase1_module.hhg_seeding_task:
            events.append("hmm")
            return [], hits
        if task is phase1_module.pfam_arbitration_task:
            raise AssertionError("missing Pfam resource reached arbitration")
        if task.__name__ == "marker_validation_task":
            events.append("marker")
            assert task_kwargs["hmm_hits"] == hits
            return []
        if task is phase1_module.generate_outputs_task:
            events.append("outputs")
            return {}
        raise AssertionError(f"unexpected task: {task}")

    monkeypatch.setattr(phase1_module, "call_task", fake_call_task)
    monkeypatch.setattr(phase1_module, "_generate_required_reports", lambda **kwargs: {})
    monkeypatch.setattr(phase1_module, "_write_empty_run_log", lambda **kwargs: None)

    result = _run_phase1_subflow(**kwargs)

    assert result["success"] is True
    assert events == ["hmm", "marker", "outputs"]


def test_confirmed_frameshift_marker_can_seed_without_a_protein_hmm_hit(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_dir = tmp_path / "output"
    kwargs = _phase1_kwargs(tmp_path, output_dir)
    kwargs["resume"] = False
    kwargs["resume_authorized"] = False
    kwargs["frameshift_screening_enabled"] = True
    kwargs["masked_path"].write_text(">contig_1\n" + "A" * 500 + "\n")
    kwargs["proteome_path"].write_text("")
    hmm_database = tmp_path / "markers.hmm"
    marker_db = tmp_path / "markers.dmnd"
    hmm_database.write_text("HMM\n")
    marker_db.write_text("DB\n")
    kwargs["hmm_database"] = hmm_database
    kwargs["marker_db"] = marker_db
    hit = FrameshiftHit(
        annotation_class=ANNOTATION_CLASS,
        hit_id="1",
        target_name="contig_1",
        target_accession="-",
        query_name="VS000001",
        query_accession="-",
        hmm_len=100,
        hmm_from=10,
        hmm_to=90,
        seq_len=500,
        ali_start=100,
        ali_end=300,
        strand="+",
        evalue=1e-20,
        score=90.0,
        bias=0.0,
        pid=40.0,
        shifts=1,
        stops=0,
        description="synthetic rescue",
    )
    protein_id = rescued_protein_id(hit)
    marker = ValidatedMarkerHit(
        query_porf=f"{protein_id}|aa1-80",
        scaffold="contig_1",
        start=100,
        end=300,
        strand="+",
        hmm_target="VS000001",
        hmm_score=90.0,
        hmm_evalue=1e-20,
        validation_status="validated",
        top10_prefixes="NCLDV__",
        best_hit_target="NCLDV__reference",
        best_hit_pident=35.0,
        best_hit_bits=80.0,
        has_ncldv=1,
        has_mirus=0,
        has_plv=0,
        has_vp=0,
        has_viral=1,
    )

    def fake_call_task(task, **task_kwargs):
        if task is phase1_module.frameshift_screening_task:
            frameshift_dir = task_kwargs["output_dir"]
            frameshift_dir.mkdir(parents=True)
            (frameshift_dir / "frameshift_candidates.faa").write_text(
                f">{protein_id} # 101 # 300 # 1 # "
                f"ID=0_{protein_id.rsplit('_', 1)[1]};"
                "annotation=frameshift_rescued_domain\n"
                + "M" * 80
                + "\n"
            )
            return [hit]
        if task is phase1_module.hhg_seeding_task:
            return [], []
        if task.__name__ == "marker_validation_task":
            if task_kwargs["output_dir"].name == "marker_validation":
                return []
            task_kwargs["output_dir"].mkdir(parents=True)
            (task_kwargs["output_dir"] / "diamond_top10.tsv").write_text(
                f"{protein_id}\tNCLDV__reference\t1e-20\t80\t35\t75\n"
            )
            return [marker]
        if task.__name__ == "region_assembly_task":
            assert task_kwargs["validated_markers"] == [marker]
            return [
                SimpleNamespace(
                    scaffold="contig_1",
                    start=50,
                    end=350,
                    length=300,
                    marker_count=1,
                    markers=[marker],
                    predicted_family="NCLDV",
                )
            ]
        raise AssertionError(f"unexpected task: {task}")

    monkeypatch.setattr(phase1_module, "call_task", fake_call_task)

    result = _run_phase1_subflow(**kwargs)

    assert result["validated_markers"] == [marker]
    assert len(result["merged_seeds"]) == 1
    assert result["merged_seeds"][0].anchors[0].porf_id == marker.query_porf
    assert result["merged_seeds"][0].sources == ["frameshift_rescue"]
    confirmed_faa = (
        output_dir
        / "phase1"
        / "frameshift_screening"
        / "confirmed_frameshift_proteins.faa"
    )
    assert protein_id in confirmed_faa.read_text()


def test_fresh_phase1_writes_exact_state_after_classification(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_dir = tmp_path / "output"
    kwargs = _phase1_kwargs(tmp_path, output_dir)
    kwargs["resume"] = False
    kwargs["resume_authorized"] = False
    hmm_database = tmp_path / "markers.hmm"
    marker_db = tmp_path / "markers.dmnd"
    hmm_database.write_text("HMM\n")
    marker_db.write_text("DB\n")
    kwargs["hmm_database"] = hmm_database
    kwargs["marker_db"] = marker_db
    marker = _complete_marker()

    def fake_call_task(task, **task_kwargs):
        if task is phase1_module.hhg_seeding_task:
            return [], [
                HMMHit(
                    query_name="scaffold/alpha_17",
                    target_name="GVOGm0003",
                    score=123.0,
                    evalue=1e-40,
                    domain_score=120.0,
                    query_start=1,
                    query_end=100,
                )
            ]
        if task.__name__ == "marker_validation_task":
            return [marker]
        if task.__name__ == "region_assembly_task":
            return [
                SimpleNamespace(
                    scaffold=marker.scaffold,
                    start=51,
                    end=1601,
                    length=1550,
                    marker_count=1,
                    markers=[marker],
                )
            ]
        raise AssertionError(f"unexpected task: {task}")

    monkeypatch.setattr(phase1_module, "call_task", fake_call_task)

    result = _run_phase1_subflow(**kwargs)
    state_path = output_dir / "phase1" / PHASE1_STATE_FILENAME
    loaded = load_phase1_state(state_path)

    assert state_path.is_file()
    assert loaded.validated_markers == result["validated_markers"]
    assert loaded.merged_seeds == result["merged_seeds"]
    assert loaded.host_signature_model == result["host_signature_model"]
    assert loaded.host_signatures == result["host_signatures"]
    assert loaded.host_deviation_summary == result["host_deviation_summary"]
    assert loaded.validated_markers[0].hmm_score == 123.456789012345
    assert loaded.validated_markers[0].has_plv == 1
    assert loaded.validated_markers[0].has_vp == 1
