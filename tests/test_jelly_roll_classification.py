import csv
import logging
import subprocess
import sys
from pathlib import Path

import pytest

from virosync.orchestration.tasks import (
    _build_jelly_roll_summary_for_boundary,
    classify_jelly_roll_task,
)
from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    EvidenceSynthesizerConfig,
    VerificationResult,
    calculate_eve_confidence,
    load_jelly_roll_data,
)
from virosync.utils.path_safety import safe_filename_component

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "classify_jelly_roll.py"
sys.path.insert(0, str(SCRIPT_PATH.parent))

from classify_jelly_roll import load_sequences  # noqa: E402


def _write_inputs(
    tmp_path: Path,
    hits: list[tuple[str, str, str]],
    sequences: list[tuple[str, str]],
) -> tuple[Path, Path]:
    marker_hits_path = tmp_path / "validated_marker_hits.tsv"
    marker_hits_path.write_text(
        "query_porf\thmm_target\tvalidation_status\n"
        + "".join(f"{porf_id}\t{marker}\t{status}\n" for porf_id, marker, status in hits)
    )
    sequences_path = tmp_path / "hmm_hit_porfs.faa"
    sequences_path.write_text("".join(f">{porf_id}\n{sequence}\n" for porf_id, sequence in sequences))
    return marker_hits_path, sequences_path


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_load_sequences_supports_exact_base_and_filter_modes(tmp_path: Path) -> None:
    _, sequences_path = _write_inputs(
        tmp_path,
        [],
        [
            ("base", "A" * 20),
            ("domain|aa4-12", "D" * 9),
            ("plain", "P" * 15),
            ("unrelated", "U" * 10),
        ],
    )

    requested = {"base|aa4-12", "domain|aa4-12", "plain", "missing|aa1-8"}

    assert set(load_sequences(sequences_path, requested)) == {
        "base",
        "domain|aa4-12",
        "plain",
    }
    assert load_sequences(sequences_path, set()) == {}
    assert set(load_sequences(sequences_path)) == {
        "base",
        "domain|aa4-12",
        "plain",
        "unrelated",
    }


def test_task_handles_sequence_forms_lengths_and_one_result_per_base(tmp_path: Path) -> None:
    marker_hits_path, sequences_path = _write_inputs(
        tmp_path,
        [
            ("full|aa1-200", "PLV_MCP", "validated"),
            ("full|aa250-460", "PLV_MCP", "validated"),
            ("domain|aa101-480", "PLV_MCP", "validated"),
            ("plain", "PLV_MCP", "validated"),
            ("empty|aa10-300", "PLV_MCP", "validated"),
            ("missing|aa1-100", "PLV_MCP", "validated"),
        ],
        [
            ("full", "F" * 520),
            ("domain|aa101-480", "D" * 80),
            ("plain", "P" * 410),
            ("empty|aa10-300", ""),
            ("empty", "E" * 430),
            ("unrelated", "U" * 600),
        ],
    )
    output_path = tmp_path / "out" / "jelly_roll.tsv"

    assert classify_jelly_roll_task(marker_hits_path, sequences_path, output_path) == output_path

    rows = {row["protein_id"]: row for row in _read_rows(output_path)}
    assert set(rows) == {"full|aa1-200", "domain|aa101-480", "plain", "empty|aa10-300"}
    assert rows["full|aa1-200"]["length"] == "520"
    assert rows["full|aa1-200"]["evidence"] == "hmm:2_domains"
    assert rows["domain|aa101-480"]["length"] == "480"
    assert rows["plain"]["length"] == "410"
    assert rows["empty|aa10-300"]["length"] == "430"
    assert {row["mcp_support"] for row in rows.values()} == {"sequence_supported"}
    assert {row["validation_status"] for row in rows.values()} == {"validated"}


def test_task_loader_summary_and_scorer_keep_support_separate_from_fold(tmp_path: Path) -> None:
    cases = [
        ("weak_mirus", "Mirus_MCP", "unvalidated", 92, "UNKNOWN", "candidate", False),
        ("valid_mirus", "Mirus_MCP", "validated_novel", 92, "HK97", "sequence_supported", True),
        ("explicit_hk97", "HK97_capsid", "validated", 92, "HK97", "sequence_supported", True),
        ("mirus_substring", "NotMirus_MCP", "validated", 92, "UNKNOWN", "sequence_supported", True),
        ("valid_generic", "MCP", "validated", 92, "UNKNOWN", "sequence_supported", True),
        ("long_candidate", "MCP", "unvalidated", 500, "DJR", "candidate", False),
        ("valid_djr", "PLV_MCP", "validated", 500, "DJR", "sequence_supported", True),
    ]
    marker_hits_path, sequences_path = _write_inputs(
        tmp_path,
        [(porf_id, marker, status) for porf_id, marker, status, *_rest in cases],
        [(porf_id, "A" * length) for porf_id, _marker, _status, length, *_rest in cases],
    )
    output_path = tmp_path / "support.tsv"

    classify_jelly_roll_task(marker_hits_path, sequences_path, output_path)

    output_rows = {row["protein_id"]: row for row in _read_rows(output_path)}
    loaded = load_jelly_roll_data(output_path)
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())

    for porf_id, _marker, status, _length, fold, support, promotes_mcp in cases:
        row = output_rows[porf_id]
        assert row["type"] == fold
        assert row["mcp_support"] == support
        assert row["validation_status"] == status

        summary = _build_jelly_roll_summary_for_boundary([(porf_id, "A")], loaded)
        assert summary is not None
        assert summary["supported_mcp_count"] == int(promotes_mcp)
        assert summary["mcp_proteins"][0]["mcp_support"] == support
        assert summary["mcp_proteins"][0]["validation_status"] == status
        result = VerificationResult(eve_id=porf_id, scaffold="scaf", start=0, end=1000)
        if "mirus" in porf_id:
            result.region_classification_mirus_markers = 1
        synthesizer._apply_jelly_roll_summary(result, summary)
        score = calculate_eve_confidence(result, crf_confidence=0.0)
        assert result.has_mcp is promotes_mcp
        profile = result.to_dict()
        assert profile["jelly_roll_supported_mcp_count"] == int(promotes_mcp)
        assert profile["jelly_roll_hk97_count"] == int(fold == "HK97")
        if porf_id == "weak_mirus":
            assert result.score_components["priority_floor"] == 0.0
            assert score < 0.7
        if porf_id == "long_candidate":
            assert result.jelly_roll_confidence_bonus == 0.0
        if porf_id == "valid_djr":
            assert result.jelly_roll_confidence_bonus > 0.0

    combined = _build_jelly_roll_summary_for_boundary(
        [("valid_djr", "A"), ("long_candidate", "A")],
        loaded,
    )
    assert combined is not None
    assert combined["djr_count"] == 2
    assert combined["supported_mcp_count"] == 1
    assert combined["total_mcp"] == 2
    supported_confidence = float(output_rows["valid_djr"]["confidence"])
    all_confidences = [float(output_rows[porf_id]["confidence"]) for porf_id in ("valid_djr", "long_candidate")]
    assert combined["avg_confidence"] == pytest.approx(sum(all_confidences) / 2)
    assert combined["confidence_bonus"] == pytest.approx(0.05 * supported_confidence * (1 / combined["total_mcp"]))

    preexisting = VerificationResult(eve_id="preexisting", scaffold="scaf", start=0, end=1000, has_mcp=True)
    weak_summary = _build_jelly_roll_summary_for_boundary([("weak_mirus", "A")], loaded)
    synthesizer._apply_jelly_roll_summary(preexisting, weak_summary)
    assert preexisting.has_mcp is True


def test_multi_domain_uses_supported_status_without_changing_output_id(tmp_path: Path) -> None:
    marker_hits_path, sequences_path = _write_inputs(
        tmp_path,
        [
            ("multi|aa1-100", "Mirus_MCP", "unvalidated"),
            ("multi|aa150-250", "Mirus_MCP", "validated_novel"),
            ("mixed|aa1-100", "Mirus_MCP", "unvalidated"),
            ("mixed|aa150-250", "MCP", "validated"),
        ],
        [("multi", "M" * 300), ("mixed", "M" * 300)],
    )
    output_path = tmp_path / "multi.tsv"

    classify_jelly_roll_task(marker_hits_path, sequences_path, output_path)

    rows = {row["protein_id"]: row for row in _read_rows(output_path)}
    assert rows["multi|aa1-100"]["validation_status"] == "validated_novel"
    assert rows["multi|aa1-100"]["mcp_support"] == "sequence_supported"
    assert rows["multi|aa1-100"]["type"] == "HK97"
    assert rows["mixed|aa1-100"]["validation_status"] == "validated"
    assert rows["mixed|aa1-100"]["mcp_support"] == "sequence_supported"
    assert rows["mixed|aa1-100"]["type"] == "UNKNOWN"


def test_foldseek_hk97_requires_quality_and_exact_query_match(tmp_path: Path) -> None:
    cases = [
        ("qualified", "HK97", "structure_supported"),
        ("weak_tm", "HK97", "candidate"),
        ("weak_evalue", "HK97", "candidate"),
        ("bad_number", "HK97", "candidate"),
        ("negative_evalue", "HK97", "candidate"),
        ("bad_range", "HK97", "candidate"),
        ("bad_field", "HK97", "candidate"),
        ("legacy", "HK97", "candidate"),
        ("arbitrary_name", "UNKNOWN", "candidate"),
        ("generic_multi|aa1-100", "HK97", "structure_supported"),
        ("gene1", "UNKNOWN", "candidate"),
        ("contig_1", "UNKNOWN", "candidate"),
        ("ambiguous", "UNKNOWN", "candidate"),
        ("ambiguous_model", "UNKNOWN", "candidate"),
        ("ambiguous_model_0", "UNKNOWN", "candidate"),
        ("exact_raw", "HK97", "structure_supported"),
        ("domain_raw|aa1-100", "HK97", "structure_supported"),
        ("encoded|aa1-100", "HK97", "structure_supported"),
    ]
    marker_hits_path, sequences_path = _write_inputs(
        tmp_path,
        [(porf_id, "MCP", "unvalidated") for porf_id, *_rest in cases]
        + [("generic_multi|aa150-250", "MCP", "unvalidated")],
        [(porf_id, "A" * 100) for porf_id, *_rest in cases],
    )
    foldseek_path = tmp_path / "foldseek.tsv"
    foldseek_path.write_text(
        "qualified_model\tpdb100#1OHG_A\t1e-10\t100\t0.8\t0.8\t0.9\t0.7\n"
        "weak_tm_model\tpdb100#1OHG_A\t1e-10\t100\t0.8\t0.8\t0.9\t0.49\n"
        "weak_evalue_model\tpdb100#1OHG_A\t0.01\t100\t0.8\t0.8\t0.9\t0.8\n"
        "bad_number_model\tpdb100#1OHG_A\tnan\t100\t0.8\t0.8\t0.9\tinf\n"
        "negative_evalue_model\tpdb100#1OHG_A\t-1e-10\t100\t0.8\t0.8\t0.9\t0.8\n"
        "bad_range_model\tpdb100#1OHG_A\t1e-10\t100\t0.8\t0.8\t0.9\t1.1\n"
        "bad_field_model\tpdb100#1OHG_A\t1e-10\t100\tnot-a-number\t0.8\t0.9\t0.8\n"
        "legacy_model\tpdb100#1OHG_A\n"
        "arbitrary_name_model\tsome_HK97_capsid\t1e-10\t100\t0.8\t0.8\t0.9\t0.8\n"
        "generic_multi_model\tpdb100#1OHG_A\t1e-10\t100\t0.8\t0.8\t0.9\t0.8\n"
        "gene10_model\tpdb100#1OHG_A\t1e-10\t100\t0.8\t0.8\t0.9\t0.8\n"
        "prefix_contig_1_model\tpdb100#1OHG_A\t1e-10\t100\t0.8\t0.8\t0.9\t0.8\n"
        "99_contig_1_model_0\tpdb100#1OHG_A\t1e-10\t100\t0.8\t0.8\t0.9\t0.8\n"
        "ambiguous_model\tpdb100#1OHG_A\t1e-10\t100\t0.8\t0.8\t0.9\t0.8\n"
        "ambiguous_model_0\tpdb100#1OHG_A\t1e-10\t100\t0.8\t0.8\t0.9\t0.8\n"
        "exact_raw\tpdb100#1OHG_A\t1e-10\t100\t0.8\t0.8\t0.9\t0.8\n"
        "domain_raw|aa1-100\tpdb100#1OHG_A\t1e-10\t100\t0.8\t0.8\t0.9\t0.8\n"
        f"{safe_filename_component('encoded|aa1-100')}_model_0\tpdb100#1OHG_A\t1e-10\t100\t0.8\t0.8\t0.9\t0.8\n"
    )
    output_path = tmp_path / "foldseek_output.tsv"

    classify_jelly_roll_task(
        marker_hits_path,
        sequences_path,
        output_path,
        foldseek_results_path=foldseek_path,
    )

    rows = {row["protein_id"]: row for row in _read_rows(output_path)}
    for porf_id, fold, support in cases:
        assert rows[porf_id]["type"] == fold
        assert rows[porf_id]["mcp_support"] == support


def test_old_tsv_defaults_to_candidate_once(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    old_path = tmp_path / "old.tsv"
    old_path.write_text(
        "protein_id\ttype\tconfidence\tlength\tmarker\tevidence\n"
        "old1\tDJR\t0.9\t500\tPLV_MCP\tlength\n"
        "old2\tUNKNOWN\t0.0\t100\tMCP\tinsufficient_evidence\n"
    )

    logger = logging.getLogger(load_jelly_roll_data.__module__)
    prior_disabled = logger.disabled
    prior_propagate = logger.propagate
    prior_global_disable = logging.root.manager.disable
    logger.disabled = False
    logger.propagate = True
    logging.disable(logging.NOTSET)
    try:
        with caplog.at_level(logging.INFO, logger=logger.name):
            loaded = load_jelly_roll_data(old_path)
    finally:
        logger.disabled = prior_disabled
        logger.propagate = prior_propagate
        logging.disable(prior_global_disable)

    assert loaded["old1"][0]["mcp_support"] == "candidate"
    notices = [record for record in caplog.records if "mcp_support column" in record.message]
    assert len(notices) == 1


def test_cli_loads_base_id_and_prefers_exact_nonempty_domain_sequence(tmp_path: Path) -> None:
    marker_hits_path, sequences_path = _write_inputs(
        tmp_path,
        [
            ("base_only|aa101-480", "PLV_MCP", "validated_novel"),
            ("both|aa101-480", "PLV_MCP", "validated"),
            ("polb", "PolB", "validated"),
        ],
        [
            ("base_only", "F" * 520),
            ("both", "B" * 600),
            ("both|aa101-480", "D" * 80),
            ("polb", "P" * 100),
            ("unrelated", "U" * 500),
        ],
    )
    output_path = tmp_path / "cli.tsv"

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--marker-hits",
            str(marker_hits_path),
            "--sequences",
            str(sequences_path),
            "--output",
            str(output_path),
            "--include-sequences",
            "--all-markers",
        ],
        cwd=REPO_ROOT,
        check=True,
    )

    rows = {row["protein_id"]: row for row in _read_rows(output_path)}
    with output_path.open() as handle:
        header = handle.readline().rstrip().split("\t")
    assert header[:6] == ["protein_id", "type", "confidence", "length", "marker", "evidence"]
    assert header[6:] == ["mcp_support", "validation_status"]
    assert set(rows) == {"base_only|aa101-480", "both|aa101-480", "polb"}
    assert rows["base_only|aa101-480"]["length"] == "520"
    assert rows["both|aa101-480"]["length"] == "480"
    assert rows["base_only|aa101-480"]["mcp_support"] == "sequence_supported"
    assert rows["base_only|aa101-480"]["validation_status"] == "validated_novel"
    assert rows["polb"]["type"] == "UNKNOWN"
    assert rows["polb"]["mcp_support"] == "candidate"
    assert rows["polb"]["validation_status"] == "validated"

    with output_path.with_suffix(".with_sequences.tsv").open() as handle:
        header = handle.readline().rstrip().split("\t")
    assert header[:7] == [
        "protein_id",
        "type",
        "confidence",
        "length",
        "marker",
        "evidence",
        "sequence",
    ]
    assert header[7:] == ["mcp_support", "validation_status"]
    sequence_rows = {row["protein_id"]: row for row in _read_rows(output_path.with_suffix(".with_sequences.tsv"))}
    assert sequence_rows["base_only|aa101-480"]["sequence"] == "F" * 520
    assert sequence_rows["both|aa101-480"]["sequence"] == "D" * 80
    assert sequence_rows["both|aa101-480"]["mcp_support"] == "sequence_supported"
    assert sequence_rows["both|aa101-480"]["validation_status"] == "validated"
