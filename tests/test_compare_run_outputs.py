"""Compare the scientific outputs of two genome output folders."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts import compare_run_outputs as compare

PREDICTIONS = "phase3_synthesis/virosync_predictions.tsv"
DETAILED = "virosync_predictions_detailed.tsv"
SUMMARY = "phase3_synthesis/virosync_summary.json"
MODEL = "phase1/marker_validation/host_signature_model.json"
EVIDENCE = "phase3_synthesis/evidence_profiles.json"
REFINED_STATE = "phase2/refined_state.json"
EVES = "demo_eves.fna"
GENE_TAXONOMY = "phase3_synthesis/gene_taxonomy"
GVCLASS = "phase3_synthesis/gvclass_input"
ALL_EVES = "phase3_synthesis/eve_sequences_all"


def _write(root: Path, relative: str, content: str) -> None:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _state(root: Path, *, status: str = "success", terminal_phase: int | None = None, rows: int = 1) -> None:
    result = {"terminal_phase": terminal_phase, "canonical_rows": rows} if status == "success" else None
    _write(
        root,
        "virosync_run_state.json",
        json.dumps({"status": status, "identities": {"genome_id": "demo"}, "result": result}),
    )


def _exports(root: Path, prefix: str, generated_at: str) -> None:
    """The six exports that every completed run holds, in ``prefix`` (a folder or the root)."""
    _write(root, f"{prefix}virosync_predictions.tsv", "eve_id\tstart\nEVE_1\t10\n")
    _write(root, f"{prefix}virosync_predictions.bed", "scaffold\t10\t20\tEVE_1\t900\t.\n")
    _write(root, f"{prefix}virosync_predictions.gff3", "##gff-version 3\n")
    _write(
        root,
        f"{prefix}virosync_summary.json",
        json.dumps({"generated_at": generated_at, "statistics": {"canonical_predictions": 1}}),
    )
    _write(root, "virosync_predictions_detailed.tsv", "eve_id\tstart\nEVE_1\t10\nEVE_2\t50\n")
    _write(root, "virosync_repeat_candidates.tsv", "eve_id\n")


def _folder(root: Path, *, generated_at: str = "2026-11-10T10:00:00", rows: int = 1) -> Path:
    """A completed run that reached Phase 3, with accepted EVEs when ``rows`` is positive."""
    _state(root, rows=rows)
    _exports(root, "phase3_synthesis/", generated_at)
    _write(root, EVIDENCE, json.dumps({"EVE_1": {"final_confidence": 0.9}}))
    _write(root, MODEL, json.dumps({"token_weights": {"a": 1.5, "b": 2.5}, "token_bits": {"a": [1.5, 0.5]}}))
    _write(root, "phase1/resume_state.json", json.dumps({"schema_version": 5}))
    _write(root, REFINED_STATE, json.dumps({"boundaries": [{"start": 10, "end": 20}]}))
    _write(root, "phase2/refined_boundaries.bed", "scaffold\t10\t20\tEVE_1\t900\t.\n")
    _write(root, "phase2/boundary_diamond/taxonomy_map.tsv", "porf\ttaxon\np1\tEUK\n")
    _write(root, "run.log", f"# Generated: {generated_at}\n")
    if rows:
        _write(root, EVES, ">EVE_1\nACGT\n")
        _write(root, f"{GENE_TAXONOMY}/gene_taxonomy_all.tsv", "eve_id\tporf\nEVE_1\tp1\n")
        _write(root, f"{GENE_TAXONOMY}/EVE_1_gene_taxonomy.tsv", "porf\tprefix\np1\tEUK\n")
        _write(root, f"{GVCLASS}/manifest.tsv", "eve_id\tnucleotide_fasta\nEVE_1\tnucleotide/EVE_1.fna\n")
        _write(root, f"{GVCLASS}/nucleotide/EVE_1.fna", ">EVE_1\nACGT\n")
    return root


def _early_exit(root: Path, *, evidence: str = '{"eve": 1}') -> Path:
    """A completed run that ended in Phase 1 without seeds: empty exports at the root."""
    _state(root, terminal_phase=1, rows=0)
    _exports(root, "", "2026-11-10T10:00:00")
    _write(root, "evidence_profiles.json", evidence)
    return root


def test_identical_folders_have_no_difference(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")

    assert compare.compare_runs(first, second) == []


def test_main_prints_and_exits_zero_for_identical_folders(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")

    assert compare.main([str(first), str(second)]) == 0

    assert "identical" in capsys.readouterr().out


def test_a_table_difference_names_the_file_and_the_first_differing_line(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    _write(second, DETAILED, "eve_id\tstart\nEVE_1\t10\nEVE_2\t51\n")

    differences = compare.compare_runs(first, second)

    assert len(differences) == 1
    assert DETAILED in differences[0]
    assert "line 3" in differences[0]
    assert "EVE_2\t50" in differences[0]
    assert "EVE_2\t51" in differences[0]


def test_main_exits_one_and_prints_each_difference(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    _write(second, DETAILED, "changed\n")
    _write(second, PREDICTIONS, "changed\n")

    assert compare.main([str(first), str(second)]) == 1

    output = capsys.readouterr().out
    assert DETAILED in output
    assert PREDICTIONS in output


def test_the_time_stamp_of_a_json_file_is_ignored_and_any_other_key_is_not(tmp_path: Path) -> None:
    first = _folder(tmp_path / "a", generated_at="2026-11-10T10:00:00")
    second = _folder(tmp_path / "b", generated_at="2026-11-11T12:30:00")
    assert compare.compare_runs(first, second) == []

    _write(second, SUMMARY, json.dumps({"generated_at": "x", "statistics": {"canonical_predictions": 2}}))
    differences = compare.compare_runs(first, second)

    assert len(differences) == 1
    assert "statistics.canonical_predictions" in differences[0]


def test_json_key_order_is_ignored_but_list_order_is_not(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    _write(second, MODEL, json.dumps({"token_bits": {"a": [1.5, 0.5]}, "token_weights": {"b": 2.5, "a": 1.5}}))
    assert compare.compare_runs(first, second) == []

    _write(second, MODEL, json.dumps({"token_weights": {"a": 1.5, "b": 2.5}, "token_bits": {"a": [0.5, 1.5]}}))
    differences = compare.compare_runs(first, second)

    assert len(differences) == 2
    assert "token_bits.a[0]" in differences[0]
    assert "token_bits.a[1]" in differences[1]


def test_a_file_in_one_folder_only_is_a_difference(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    (second / DETAILED).unlink()

    differences = compare.compare_runs(first, second)

    assert any("only in" in message and DETAILED in message for message in differences)


def test_the_same_table_in_another_location_is_a_difference(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    (second / PREDICTIONS).rename(second / "virosync_predictions.tsv")

    differences = compare.compare_runs(first, second)

    assert len(differences) == 1
    assert "virosync_predictions.tsv" in differences[0]
    assert "location" in differences[0]


def test_files_outside_the_scientific_set_are_not_compared(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    _write(second, "run.log", "# Generated: another day\n")
    _write(second, "phase1/reuse_provenance.json", "{}")

    assert compare.compare_runs(first, second) == []


@pytest.mark.parametrize(
    ("relative", "content"),
    [
        (EVES, ">EVE_1\nACGA\n"),
        (f"{GENE_TAXONOMY}/gene_taxonomy_all.tsv", "eve_id\tporf\nEVE_1\tp2\n"),
        (f"{GENE_TAXONOMY}/EVE_1_gene_taxonomy.tsv", "porf\tprefix\np1\tBAC\n"),
        (f"{GVCLASS}/nucleotide/EVE_1.fna", ">EVE_1\nTTTT\n"),
        (f"{GVCLASS}/manifest.tsv", "eve_id\tnucleotide_fasta\nEVE_2\t\n"),
        (REFINED_STATE, json.dumps({"boundaries": [{"start": 10, "end": 21}]})),
        (EVIDENCE, json.dumps({"EVE_1": {"final_confidence": 0.8}})),
    ],
)
def test_sequence_taxonomy_and_checkpoint_outputs_are_compared(tmp_path: Path, relative: str, content: str) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    _write(second, relative, content)

    differences = compare.compare_runs(first, second)

    assert differences
    assert all(relative.split("/")[-1] in message or relative in message for message in differences)


@pytest.mark.parametrize("directory", [GENE_TAXONOMY, GVCLASS])
def test_a_directory_member_in_one_folder_only_is_a_difference(tmp_path: Path, directory: str) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    _write(second, f"{directory}/extra_member.tsv", "x\n")

    differences = compare.compare_runs(first, second)

    assert len(differences) == 1
    assert directory in differences[0]
    assert "extra_member.tsv" in differences[0]


def test_the_all_eve_sequence_directory_is_compared_when_present(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    _write(first, f"{ALL_EVES}/nucleotide/EVE_1.fna", ">EVE_1\nACGT\n")
    _write(second, f"{ALL_EVES}/nucleotide/EVE_1.fna", ">EVE_1\nACGG\n")

    differences = compare.compare_runs(first, second)

    assert len(differences) == 1
    assert ALL_EVES in differences[0]


def test_the_root_evidence_profiles_of_a_run_that_ended_early_are_compared(tmp_path: Path) -> None:
    first, second = _early_exit(tmp_path / "a"), _early_exit(tmp_path / "b")
    assert compare.compare_runs(first, second) == []

    _write(second, "evidence_profiles.json", '{"eve": 2}')
    differences = compare.compare_runs(first, second)

    assert len(differences) == 1
    assert "evidence_profiles.json" in differences[0]


def test_a_run_that_ended_early_needs_only_the_exports(tmp_path: Path) -> None:
    first, second = _early_exit(tmp_path / "a"), _early_exit(tmp_path / "b")

    assert compare.compare_runs(first, second) == []


def test_runs_that_ended_in_different_places_are_a_difference(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a", rows=0), _early_exit(tmp_path / "b")

    differences = compare.compare_runs(first, second)

    assert any("ended" in message for message in differences)


def test_two_incomplete_completed_runs_do_not_pass(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    for folder in (first, second):
        (folder / REFINED_STATE).unlink()
        (folder / EVES).unlink()

    differences = compare.compare_runs(first, second)

    assert len(differences) == 4
    assert all("required output is missing" in message for message in differences)
    assert any(REFINED_STATE in message for message in differences)
    assert any(EVES in message for message in differences)


@pytest.mark.parametrize(
    "required",
    [
        PREDICTIONS,
        "phase3_synthesis/virosync_predictions.bed",
        "phase3_synthesis/virosync_predictions.gff3",
        SUMMARY,
        DETAILED,
        "virosync_repeat_candidates.tsv",
        EVIDENCE,
        "phase1/resume_state.json",
        MODEL,
        "phase2/refined_boundaries.bed",
        "phase2/boundary_diamond/taxonomy_map.tsv",
        f"{GENE_TAXONOMY}/gene_taxonomy_all.tsv",
        f"{GVCLASS}/manifest.tsv",
    ],
)
def test_each_required_output_of_a_completed_run_is_enforced(tmp_path: Path, required: str) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    for folder in (first, second):
        (folder / required).unlink()

    differences = compare.compare_runs(first, second)

    assert differences
    assert any("required output is missing" in message for message in differences)


def test_a_run_without_accepted_eves_needs_no_sequence_outputs(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a", rows=0), _folder(tmp_path / "b", rows=0)

    assert compare.compare_runs(first, second) == []


def test_a_run_that_did_not_succeed_is_a_difference(tmp_path: Path) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    _state(second, status="failed")

    differences = compare.compare_runs(first, second)

    assert any("'failed'" in message and "success" in message for message in differences)


def test_a_missing_folder_exits_two(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    first = _folder(tmp_path / "a")

    assert compare.main([str(first), str(tmp_path / "missing")]) == 2

    assert "missing" in capsys.readouterr().err


def test_a_folder_without_a_run_state_exits_two(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    first, second = _folder(tmp_path / "a"), _folder(tmp_path / "b")
    (second / "virosync_run_state.json").unlink()

    assert compare.main([str(first), str(second)]) == 2

    assert "virosync_run_state.json" in capsys.readouterr().err
