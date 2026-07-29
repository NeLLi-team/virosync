from __future__ import annotations

"""Schema regression tests for ViroSync detailed TSV exports."""

import csv
from pathlib import Path

from virosync.pipeline.phase3.evidence_synthesizer import VerificationResult
from virosync.pipeline.phase3.output_generator import OutputGenerator
from virosync.validation.tsv_invariants import run_tsv_invariant_checks


def _write_validated_marker_hits(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "query_porf\tscaffold\tstart\tend\tstrand\thmm_target\thmm_score\tvalidation_status",
                "protA\tscaf1\t110\t150\t+\tGVOGm0003\t50\tvalidated",
                "protA\tscaf1\t110\t150\t+\tOG1000\t40\tvalidated",
                "protB\tscaf1\t160\t200\t+\tGVOGm0003\t60\tvalidated",
                "protC\tscaf1\t210\t250\t+\tGVOGm0007\t70\tvalidated",
                "protD\tscaf1\t260\t300\t+\tOG2000\t80\tvalidated",
                "protE\tscaf1\t310\t350\t+\tOG2000\t90\tvalidated_novel",
                "protF\tscaf1\t360\t390\t+\tGVOGm0007\t20\tunvalidated",
            ]
        )
        + "\n"
    )


def test_write_predictions_detailed_tsv_uses_simplified_schema_and_marker_counts(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "phase3_synthesis"
    output_dir.mkdir()
    _write_validated_marker_hits(
        tmp_path / "phase1" / "marker_validation" / "validated_marker_hits.tsv"
    )

    generator = OutputGenerator(output_dir=output_dir, extended_output=True)
    result = VerificationResult(
        eve_id="EVE_scaf1_100-400",
        scaffold="scaf1",
        start=100,
        end=400,
        confidence_tier="HIGH",
        final_confidence=0.95,
        likely_family="VP",
        gene_count=5,
        gene_taxonomy_total=5,
        gene_taxonomy_records=[
            {
                "top1_target": "EUK__host_a",
                "top1_prefix": "EUK",
                "top10_prefixes": ["NCLDV"],
                "is_flanking": False,
            },
            {
                "top1_target": "MIRUS__hit",
                "top1_prefix": "MIRUS",
                "top10_prefixes": ["MIRUS"],
                "is_flanking": False,
            },
            {
                "top1_target": "VP__hit",
                "top1_prefix": "VP",
                "top10_prefixes": ["VP"],
                "is_flanking": False,
            },
            {
                "top1_target": "PLV__hit",
                "top1_prefix": "PLV",
                "top10_prefixes": ["PLV"],
                "is_flanking": False,
            },
            {
                "top1_target": "",
                "top1_prefix": "UNKNOWN",
                "top10_prefixes": [],
                "is_flanking": False,
            },
        ],
        host_signature_gene_count=1,
        host_signature_fraction=0.2,
        host_signature_weighted_mean=0.4,
    )

    detailed_tsv = generator.write_predictions_detailed_tsv([result])

    with detailed_tsv.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 1
    row = rows[0]
    ordered_fieldnames = list(reader.fieldnames or [])
    fieldnames = set(ordered_fieldnames)

    assert "classification" in fieldnames
    assert "likely_group" in fieldnames  # capscan MCP group column
    assert "hallmark_total" in fieldnames
    assert "hallmark_unique" in fieldnames
    assert "plv_top10_proteins" in fieldnames
    assert "vp_top10_proteins" in fieldnames
    assert ordered_fieldnames[-1] == "effective_eve_class"

    assert "status" not in fieldnames
    assert "likely_family" not in fieldnames
    assert "crf_confidence" not in fieldnames
    assert "crf_posterior" not in fieldnames
    assert "seed_compositional_score" not in fieldnames
    assert "cub_deviation" not in fieldnames
    assert "coherence_score" not in fieldnames
    assert "structural_score" not in fieldnames
    assert "hallmark_count" not in fieldnames
    assert "hallmark_diversity" not in fieldnames
    assert "gvogm_ncldv_mirus_names" not in fieldnames
    assert "og_ncldv_mirus_names" not in fieldnames
    assert "interproscan_numt_hits" not in fieldnames
    assert "interproscan_numt_markers" not in fieldnames
    assert "numt_flag" not in fieldnames
    assert "candidate_common_euk_taxonomy" not in fieldnames
    assert "phase1_host_signature_host_prefixes" not in fieldnames
    assert "phase1_host_signature_top_tokens" not in fieldnames
    assert "phase1_host_signature_token_count" not in fieldnames

    # GVClass unified VP and PLV into Preplasmiviricota; a legacy "VP" label
    # resolves to the PPV lineage in the public class column.
    assert row["classification"] == "VP"
    assert row["effective_eve_class"] == "PPV"
    assert row["likely_group"] == "."  # no capscan group-defining marker -> "."
    assert row["hallmark_total"] == "5"
    assert row["hallmark_unique"] == "4"
    assert row["gvogm_count"] == "3"
    assert row["gvogm_names"] == "GVOGm0003:2,GVOGm0007:1"
    assert row["og_count"] == "3"
    assert row["og_names"] == "OG1000:1,OG2000:2"
    assert row["gvogm_unvalidated_count"] == "1"
    assert row["gvogm_unvalidated_names"] == "GVOGm0007"
    assert row["ncldv_top10_proteins"] == "1"
    assert row["mirus_top10_proteins"] == "1"
    assert row["plv_top10_proteins"] == "1"
    assert row["vp_top10_proteins"] == "1"
    assert row["host_signature_fraction"] == "0.2000"


def test_detailed_tsv_invariants_allow_counted_names_to_exceed_unique_protein_count(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "phase3_synthesis"
    output_dir.mkdir()
    _write_validated_marker_hits(
        tmp_path / "phase1" / "marker_validation" / "validated_marker_hits.tsv"
    )

    generator = OutputGenerator(output_dir=output_dir, extended_output=True)
    result = VerificationResult(
        eve_id="EVE_scaf1_100-400",
        scaffold="scaf1",
        start=100,
        end=400,
        confidence_tier="HIGH",
        final_confidence=0.95,
        likely_family="VP",
        gene_count=5,
        gene_taxonomy_total=5,
        gene_taxonomy_records=[
            {
                "top1_target": "EUK__host_a",
                "top1_prefix": "EUK",
                "top10_prefixes": ["NCLDV"],
                "is_flanking": False,
            },
            {
                "top1_target": "MIRUS__hit",
                "top1_prefix": "MIRUS",
                "top10_prefixes": ["MIRUS"],
                "is_flanking": False,
            },
            {
                "top1_target": "VP__hit",
                "top1_prefix": "VP",
                "top10_prefixes": ["VP"],
                "is_flanking": False,
            },
            {
                "top1_target": "PLV__hit",
                "top1_prefix": "PLV",
                "top10_prefixes": ["PLV"],
                "is_flanking": False,
            },
            {
                "top1_target": "",
                "top1_prefix": "UNKNOWN",
                "top10_prefixes": [],
                "is_flanking": False,
            },
        ],
        host_signature_gene_count=1,
        host_signature_fraction=0.2,
        host_signature_weighted_mean=0.4,
    )

    detailed_tsv = generator.write_predictions_detailed_tsv([result])
    report = run_tsv_invariant_checks(detailed_tsv)

    assert report.passed, report.issues


def test_detailed_tsv_ignores_stale_taxonomy_total_without_records(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "phase3_synthesis"
    output_dir.mkdir()
    proteome = tmp_path / "phase0" / "proteome.fasta"
    proteome.parent.mkdir()
    proteome.write_text(
        "\n".join(
            [
                ">scaf1_1 # 1 # 100 # 1 # ID=scaf1_1;",
                "MAAA",
                ">scaf1_2 # 150 # 250 # 1 # ID=scaf1_2;",
                "MAAA",
                ">scaf1_3 # 301 # 350 # 1 # ID=scaf1_3;",
                "MAAA",
            ]
        )
        + "\n"
    )

    generator = OutputGenerator(
        output_dir=output_dir,
        proteome_fasta=proteome,
        extended_output=True,
    )
    result = VerificationResult(
        eve_id="EVE_scaf1_100-300",
        scaffold="scaf1",
        start=100,
        end=300,
        confidence_tier="LOW",
        final_confidence=0.1,
        likely_family="UNKNOWN",
        gene_count=3,
        gene_taxonomy_total=3,
        gene_taxonomy_records=[],
    )

    detailed_tsv = generator.write_predictions_detailed_tsv([result])

    with detailed_tsv.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))

    assert row["total_proteins"] == "1"
    assert "UNK:1" in row["taxonomy_best_hits"]
    assert "NO_HITS:0" in row["taxonomy_best_hits"]

    report = run_tsv_invariant_checks(detailed_tsv)
    assert report.passed, report.issues
