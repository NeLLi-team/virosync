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
        "query_porf\thmm_target\tvalidation_status\thmm_score\n"
        + "".join(f"{porf_id}\t{marker}\t{status}\t50.0\n" for porf_id, marker, status in hits)
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

        summary = _build_jelly_roll_summary_for_boundary([(porf_id, "A")], loaded, {})
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
        assert result.has_capsid_evidence is promotes_mcp
        profile = result.to_dict()
        assert profile["jelly_roll_supported_mcp_count"] == int(promotes_mcp)
        assert profile["jelly_roll_hk97_count"] == int(fold == "HK97")
        # Fold annotation alone is not retained catalogue HMM evidence.
        assert result.score_components["hallmark_summary"]["core_count"] == 0
        assert result.score_components["hallmark_floor"] == 0.0
        assert score < 0.2

    combined = _build_jelly_roll_summary_for_boundary(
        [("valid_djr", "A"), ("long_candidate", "A")],
        loaded,
        {},
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
    weak_summary = _build_jelly_roll_summary_for_boundary([("weak_mirus", "A")], loaded, {})
    synthesizer._apply_jelly_roll_summary(preexisting, weak_summary)
    assert preexisting.has_mcp is True


def test_curated_structure_mcp_reaches_jelly_roll_selection(tmp_path: Path) -> None:
    marker_hits_path, sequences_path = _write_inputs(
        tmp_path,
        [
            ("curated_mcp", "EVE_CAPSID_NOVEL_0029", "validated_novel"),
            ("other_capsid", "EVE_CAPSID_NOVEL_0028", "validated_novel"),
            ("substring_mcp", "NotMirus_MCP", "validated_novel"),
        ],
        [
            ("curated_mcp", "M" * 100),
            ("other_capsid", "C" * 100),
            ("substring_mcp", "S" * 100),
        ],
    )
    output_path = tmp_path / "curated.tsv"

    classify_jelly_roll_task(marker_hits_path, sequences_path, output_path)

    rows = _read_rows(output_path)
    assert len(rows) == 1
    assert rows[0]["marker"] == "EVE_CAPSID_NOVEL_0029"
    assert rows[0]["mcp_support"] == "sequence_supported"


def test_batched_scoring_counts_rescue_linked_fragments_once_for_capsid_fold(tmp_path, monkeypatch):
    from virosync.orchestration import tasks, utils
    from virosync.pipeline.phase1.seed_merger import ConfirmedSeed
    from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary, boundary_candidate_id

    genome = tmp_path / "genome.fna"
    genome.write_text(">s\n" + "ACGT" * 500 + "\n")
    boundary = RefinedBoundary("s", 100, 1000, confirmed_seeds=(ConfirmedSeed("seed", "s", 100, 1000),))
    identity = boundary_candidate_id(boundary)
    rescued_id = "s_VSR0123456789abcdef"
    markers = [
        {
            "porf_id": gene,
            "hallmark_gene": "VS000309",
            "hmm_score": score,
            "scoring_locus_id": rescued_id,
            "validation_status": "validated",
            "top10_prefixes": "NCLDV__",
            "top10_pidents": "80",
        }
        for gene, score in [("left", 65), ("right", 60), (rescued_id, 144)]
    ]
    folds = {
        gene: [
            {
                "porf_id": f"{gene}|aa1-20",
                "classification": "DJR",
                "confidence": confidence,
                "mcp_support": "sequence_supported",
                "validation_status": "validated",
            }
        ]
        for gene, confidence in [("left", 0.8), ("right", 0.9)]
    }
    monkeypatch.setattr(utils, "get_genes_for_boundary", lambda **kwargs: [("left", "AA"), ("right", "AA")])
    [result] = tasks.verify_eve_candidates_batched_task(
        boundaries=[boundary],
        genome_path=genome,
        work_dir=tmp_path / "work",
        proteome_path=tmp_path / "unused.faa",
        hallmark_hits_map={identity: markers},
        hallmark_evidence_map={},
        hallmark_witnesses_map={},
        candidate_annotations={},
        jelly_roll_map=folds,
        max_workers=1,
    )
    assert result.jelly_roll_total_mcp == result.jelly_roll_supported_mcp_count == 1
    assert result.jelly_roll_djr_count == 1
    assert result.jelly_roll_confidence_bonus == pytest.approx(0.05 * 0.9)
    assert {record["porf_id"] for record in result.jelly_roll_mcp_proteins} == {"left|aa1-20", "right|aa1-20"}


@pytest.mark.parametrize(
    "other_fold,other_confidence,other_support,linked_bonus",
    [("SJR", 0.95, "sequence_supported", 0.0), ("DJR", 0.8, "candidate", 0.045), ("DJR", 0.95, "candidate", 0.045)],
)
def test_fold_alias_recomputes_counts_without_discarding_raw_annotations(
    other_fold, other_confidence, other_support, linked_bonus
):
    folds = {
        "left": [{"porf_id": "left", "classification": "DJR", "confidence": 0.9, "mcp_support": "sequence_supported"}],
        "right": [
            {
                "porf_id": "right",
                "classification": other_fold,
                "confidence": other_confidence,
                "mcp_support": other_support,
            }
        ],
    }
    sequences = [("left", "AA"), ("right", "AA")]
    separate = _build_jelly_roll_summary_for_boundary(sequences, folds, {})
    linked = _build_jelly_roll_summary_for_boundary(sequences, folds, {"left": "locus", "right": "locus"})
    assert separate["total_mcp"] == 2
    assert separate["confidence_bonus"] == pytest.approx(0.0225)
    assert linked["total_mcp"] == 1
    assert linked["confidence_bonus"] == pytest.approx(linked_bonus)
    assert linked["supported_mcp_count"] == 1
    assert linked["avg_confidence"] == (max(0.9, other_confidence) if other_support != "candidate" else 0.9)
    assert linked["mcp_proteins"] == separate["mcp_proteins"]


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


@pytest.mark.parametrize(
    "status,expected_fold", [("validated", "HK97"), ("validated_novel", "HK97"), ("unvalidated", "UNKNOWN")]
)
def test_mirus_vs_alias_cannot_use_length_or_multiple_domains_as_djr(
    tmp_path: Path, status: str, expected_fold: str
) -> None:
    marker_hits_path, sequences_path = _write_inputs(
        tmp_path,
        [
            ("long", "VS000005", status),
            ("multi|aa1-100", "VS000005", status),
            ("multi|aa150-250", "VS000005", status),
        ],
        [("long", "M" * 500), ("multi", "M" * 300)],
    )
    output_path = tmp_path / "mirus_alias.tsv"

    classify_jelly_roll_task(marker_hits_path, sequences_path, output_path)

    rows = {row["protein_id"]: row for row in _read_rows(output_path)}
    assert rows["long"]["type"] == expected_fold
    assert rows["multi|aa1-100"]["type"] == expected_fold


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


@pytest.mark.parametrize("score", ["9.999", "9.9996", "bad", "", "nan", "inf"])
@pytest.mark.parametrize("model", ["PLV_MCP", "Mirus_MCP"])
def test_sequence_support_requires_full_score_floor(tmp_path, score, model):
    from classify_jelly_roll import classify_proteins, load_marker_hits

    path = tmp_path / "hits.tsv"
    path.write_text(f"query_porf\thmm_target\thmm_score\tvalidation_status\nprotein\t{model}\t{score}\tvalidated\n")
    (result,) = classify_proteins(load_marker_hits(path), {"protein": "A" * 500})
    assert result.mcp_support == "candidate"
    assert result.validation_status == "validated"
    if model == "Mirus_MCP":
        assert result.jelly_roll_type != "HK97"


def test_low_validated_hit_preserves_independent_foldseek_support(tmp_path):
    from classify_jelly_roll import classify_proteins, load_foldseek_results, load_marker_hits

    path = tmp_path / "hits.tsv"
    path.write_text("query_porf\thmm_target\thmm_score\tvalidation_status\nprotein\tMirus_MCP\t9.999\tvalidated\n")
    foldseek = tmp_path / "foldseek.tsv"
    foldseek.write_text("protein_model\tpdb100#1OHG_A\t1e-10\t100\t0.8\t0.8\t0.9\t0.7\n")
    (result,) = classify_proteins(
        load_marker_hits(path), {"protein": "A" * 500}, foldseek_results=load_foldseek_results(foldseek)
    )
    assert result.mcp_support == "structure_supported"
    assert result.jelly_roll_type == "HK97"
    assert result.validation_status == "validated"


@pytest.mark.parametrize("reverse", [False, True])
def test_eligible_alternate_hit_supports_same_base_protein(tmp_path, reverse):
    from classify_jelly_roll import classify_proteins, load_marker_hits

    rows = ["protein|aa1-100\tPLV_MCP\t9.999\tvalidated\n", "protein|aa200-300\tMirus_MCP\t10\tvalidated_novel\n"]
    if reverse:
        rows.reverse()
    path = tmp_path / "hits.tsv"
    path.write_text("query_porf\thmm_target\thmm_score\tvalidation_status\n" + "".join(rows))
    (result,) = classify_proteins(load_marker_hits(path), {"protein": "A" * 500})
    assert result.protein_id == rows[0].split("\t")[0]
    assert result.mcp_support == "sequence_supported"
    assert result.validation_status == rows[0].strip().split("\t")[-1]


@pytest.mark.parametrize("score, support", [(9.9996, "candidate"), (10.0, "sequence_supported")])
def test_validation_writer_preserves_capsid_score_boundary(tmp_path, score, support):
    from classify_jelly_roll import classify_proteins, load_marker_hits

    from virosync.pipeline.phase1.hhg_seeding import HMMHit
    from virosync.pipeline.phase1.marker_validation import filter_validated_markers

    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">scaffold_1 # 101 # 1600 # 1 # ID=1_1\n" + "M" * 500 + "\n")
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("scaffold_1\tNCLDV__virus|p1\t1e-30\t100\t90\t100\n")
    hit = HMMHit("scaffold_1", "PLV_MCP", score, 1e-8, score, 1, 100)
    (marker,) = filter_validated_markers([hit], diamond, proteome, output_dir=tmp_path / "validated")
    assert marker.validation_status == "validated"
    path = tmp_path / "validated/validated_marker_hits.tsv"
    assert float(_read_rows(path)[0]["hmm_score"]) == score
    (result,) = classify_proteins(load_marker_hits(path), {"scaffold_1": "M" * 500})
    assert result.mcp_support == support
    assert result.validation_status == marker.validation_status


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("alternate", ["weak", "blocked", "eligible_hk97"])
def test_duplicate_query_alternatives_preserve_eligible_support(tmp_path, reverse, alternate):
    from classify_jelly_roll import classify_proteins, count_hmm_domains_per_protein, load_marker_hits

    other = {
        "weak": "Mirus_MCP\t9.9996\tvalidated\t",
        "blocked": "Mirus_MCP\t50\tvalidated\tcellular",
        "eligible_hk97": "Mirus_MCP\t10\tvalidated\t",
    }[alternate]
    rows = ["protein|aa1-100\tPLV_MCP\t10\tvalidated\t\n", f"protein|aa1-100\t{other}\n"]
    if reverse:
        rows.reverse()
    path = tmp_path / "hits.tsv"
    path.write_text("query_porf\thmm_target\thmm_score\tvalidation_status\tpfam_override\n" + "".join(rows))
    hits = load_marker_hits(path)
    assert len(hits["protein|aa1-100"]) == 2
    assert count_hmm_domains_per_protein(hits) == {"protein": 1}
    (result,) = classify_proteins(hits, {"protein": "A" * 500})
    assert result.mcp_support == "sequence_supported"
    assert result.validation_status == "validated"
    assert result.protein_id == "protein|aa1-100"
    assert result.marker == rows[0].split("\t")[1]
    if alternate == "eligible_hk97":
        assert result.jelly_roll_type == "HK97"
    else:
        assert result.jelly_roll_type != "HK97"


@pytest.mark.parametrize("structural", [False, True])
def test_blocked_validated_capsid_needs_independent_structural_support(tmp_path, structural):
    from classify_jelly_roll import classify_proteins, load_foldseek_results, load_marker_hits

    path = tmp_path / "hits.tsv"
    path.write_text(
        "query_porf\thmm_target\thmm_score\tvalidation_status\tpfam_override\nprotein\tMirus_MCP\t50\tvalidated\tcellular\n"
    )
    foldseek = tmp_path / "foldseek.tsv"
    foldseek.write_text("protein_model\tpdb100#1OHG_A\t1e-10\t100\t0.8\t0.8\t0.9\t0.7\n")
    (result,) = classify_proteins(
        load_marker_hits(path),
        {"protein": "A" * 500},
        foldseek_results=load_foldseek_results(foldseek) if structural else None,
    )
    assert result.mcp_support == ("structure_supported" if structural else "candidate")
    assert result.validation_status == "validated"
    assert result.jelly_roll_type == ("HK97" if structural else "UNKNOWN")
