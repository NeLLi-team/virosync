from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path

import pytest
from Bio import SeqIO

from virosync.ablation import AblationID
from virosync.config.pipeline_config import PipelineConfig
from virosync.orchestration._flows.single_genome.phase3 import _annotate_integration_genes
from virosync.output_contract import (
    COORDINATE_CONVENTION,
    COORDINATE_SCHEMA_VERSION,
    DETAILED_PREDICTION_COLUMNS,
    DETAILED_PREDICTION_EXTENDED_COLUMNS,
    INTEGRATION_EVIDENCE_COLUMNS,
    OUTPUT_SCHEMA_VERSION,
)
from virosync.pipeline.phase2.boundary_refiner import (
    RefinedBoundary,
    merge_adjacent_viral_boundaries,
)
from virosync.pipeline.phase3.acceptance_selection import select_phase3_acceptance
from virosync.pipeline.phase3.evidence_graph import (
    CoherenceAnalysis,
    EvidenceCorrelationGraph,
    EvidenceProfile,
    EvidenceType,
    WindowEvidence,
)
from virosync.pipeline.phase3.evidence_synthesizer import (
    VerificationResult,
    VerificationStatus,
)
from virosync.pipeline.phase3.output_generator import (
    OutputGenerator,
    evaluate_v2_quality_gate,
)

CANONICAL_BASE_FIELDS = (
    "eve_id",
    "scaffold",
    "start",
    "end",
    "length",
    "confidence_tier",
    "final_confidence",
    "region_classification",
    "region_classification_ncldv_markers",
    "region_classification_vp_plv_markers",
    "region_classification_mirus_markers",
    "classification",
    "likely_group",
    "kfd",
    "gc_deviation",
    "hallmark_total",
    "hallmark_unique",
    "hallmark_non_atpase",
    "has_virus_specific",
    "has_structural_support",
    "mcp_gene_ids",
    "predicted_taxonomy",
    "taxonomy_confidence",
    "gene_taxonomy_total",
    "gene_taxonomy_ncldv_top10",
    "gene_taxonomy_mirus_top10",
    "gene_taxonomy_phage_top10",
    "gene_taxonomy_viral_top10",
    "gene_taxonomy_total_with_flanking",
    "gene_taxonomy_flanking_count",
    "gene_taxonomy_viral_interior",
    "gene_taxonomy_viral_flanking",
    "gene_taxonomy_cellular",
    "gene_taxonomy_unknown",
    "gene_taxonomy_has_ncldv_mirus",
    "interproscan_total_hits",
    "interproscan_viral_hits",
    "interproscan_keyword_hits",
    "candidate_start",
    "candidate_end",
    "candidate_length",
    "candidate_reduction_bp",
    "candidate_reduction_reason",
    *INTEGRATION_EVIDENCE_COLUMNS,
)
CANONICAL_EXTENDED_FIELDS = (
    "interproscan_category_hits",
    "interproscan_family_hits",
    "interproscan_category_score",
    "interproscan_score",
    "gene_taxonomy_vp_plv_top10",
    "gene_taxonomy_dominant_family",
    "gene_taxonomy_dominant_fraction",
    "ppv_subtype",
    "host_signature_gene_count",
    "host_signature_fraction",
    "host_signature_weighted_mean",
    "marker_category_hits",
    "marker_family_hits",
    "marker_complement_score",
    "family_consistency_score",
    "vp_completeness",
    "ppv_completeness",
    "ncldv_completeness",
    "mirus_completeness",
    "seed_marker_names",
    "other_marker_names",
    "seed_marker_patterns",
    "other_marker_patterns",
)


def _write_genome_fasta(path: Path) -> None:
    path.write_text(">contig_1\n" + ("A" * 160) + "\n>contig_10\n" + ("C" * 120) + "\n")


def _write_proteome_fasta(path: Path) -> None:
    path.write_text(
        ">contig_1_1 # 10 # 40 # + # ID=1_1;\n"
        "MPEPTIDE\n"
        ">contig_1_2 # 80 # 110 # + # ID=1_2;\n"
        "MKSECOND\n"
        ">contig_10_1 # 12 # 35 # + # ID=10_1;\n"
        "MOTHERSEQ\n"
    )


def _build_result(
    *,
    eve_id: str,
    scaffold: str,
    start: int,
    end: int,
    confidence_tier: str,
    status: VerificationStatus,
    region_classification: str = "NCLDV",
    likely_family: str = "NCLDV",
    taxonomy_class: str = "NCLDV",
    hallmark_count: int = 0,
    has_mcp: bool = False,
) -> VerificationResult:
    return VerificationResult(
        eve_id=eve_id,
        scaffold=scaffold,
        start=start,
        end=end,
        confidence_tier=confidence_tier,
        status=status,
        final_confidence=0.9 if confidence_tier != "LOW" else 0.1,
        region_classification=region_classification,
        likely_family=likely_family,
        taxonomy_class=taxonomy_class,
        hallmark_count=hallmark_count,
        has_mcp=has_mcp,
    )


def test_write_gvclass_export_uses_coordinate_overlap(tmp_path: Path) -> None:
    genome_fasta = tmp_path / "genome.fna"
    proteome_fasta = tmp_path / "proteome.faa"
    _write_genome_fasta(genome_fasta)
    _write_proteome_fasta(proteome_fasta)

    generator = OutputGenerator(
        output_dir=tmp_path,
        genome_fasta=genome_fasta,
        proteome_fasta=proteome_fasta,
    )
    results = [
        _build_result(
            eve_id="EVE_1",
            scaffold="contig_1",
            start=0,
            end=50,
            confidence_tier="HIGH",
            status=VerificationStatus.HIGH_CONFIDENCE,
        ),
        _build_result(
            eve_id="EVE_2",
            scaffold="contig_1",
            start=60,
            end=120,
            confidence_tier="LOW",
            status=VerificationStatus.AMBIGUOUS,
        ),
    ]

    export_dir = generator.write_gvclass_export(results, tmp_path / "gvclass")

    protein_records = list(SeqIO.parse(export_dir / "protein" / "EVE_1.faa", "fasta"))

    assert [record.id for record in protein_records] == ["contig_1_1"]


def test_write_gvclass_export_includes_confirmed_frameshift_proteins(
    tmp_path: Path,
) -> None:
    genome_fasta = tmp_path / "genome.fna"
    proteome_fasta = tmp_path / "proteome.faa"
    confirmed_faa = tmp_path / "phase1" / "frameshift_screening" / "confirmed_frameshift_proteins.faa"
    confirmed_faa.parent.mkdir(parents=True)
    _write_genome_fasta(genome_fasta)
    _write_proteome_fasta(proteome_fasta)
    confirmed_faa.write_text(
        ">contig_1_VSR0123456789abcdef # 20 # 45 # 1 # "
        "ID=0_VSR0123456789abcdef;annotation=frameshift_rescued_domain;"
        "model=VS000001;shifts=1;stops=0;literal_stops=0\n"
        "MRESCUED\n"
    )
    confirmed_faa.with_name("confirmed_frameshift_markers.tsv").write_text(
        "query_porf\tscaffold\tstart\tend\tstrand\thmm_target\thmm_score\t"
        "validation_status\n"
        "contig_1_VSR0123456789abcdef|aa1-8\tcontig_1\t19\t45\t+\t"
        "VS000001\t80\tvalidated\n"
    )
    generator = OutputGenerator(
        output_dir=tmp_path,
        genome_fasta=genome_fasta,
        proteome_fasta=proteome_fasta,
    )
    result = _build_result(
        eve_id="EVE_1",
        scaffold="contig_1",
        start=0,
        end=50,
        confidence_tier="HIGH",
        status=VerificationStatus.HIGH_CONFIDENCE,
    )

    export_dir = generator.write_gvclass_export([result], tmp_path / "gvclass")
    records = list(SeqIO.parse(export_dir / "protein" / "EVE_1.faa", "fasta"))

    assert [record.id for record in records] == [
        "contig_1_1",
        "contig_1_VSR0123456789abcdef",
    ]
    assert "annotation=frameshift_rescued_domain" in records[1].description
    assert generator._marker_names_for_region(
        "contig_1",
        0,
        50,
        status_filter="validated",
    ) == ["VS000001"]


def test_write_gvclass_export_writes_rescue_when_ordinary_proteome_is_empty(
    tmp_path: Path,
) -> None:
    genome_fasta = tmp_path / "genome.fna"
    proteome_fasta = tmp_path / "proteome.faa"
    confirmed_faa = tmp_path / "phase1" / "frameshift_screening" / "confirmed_frameshift_proteins.faa"
    confirmed_faa.parent.mkdir(parents=True)
    _write_genome_fasta(genome_fasta)
    proteome_fasta.write_text("")
    confirmed_faa.write_text(
        ">contig_1_VSR0123456789abcdef # 20 # 45 # 1 # "
        "ID=0_VSR0123456789abcdef;annotation=frameshift_rescued_domain;"
        "model=VS000001;shifts=1;stops=0;literal_stops=0\n"
        "MRESCUED\n"
    )
    generator = OutputGenerator(
        output_dir=tmp_path,
        genome_fasta=genome_fasta,
        proteome_fasta=proteome_fasta,
    )
    result = _build_result(
        eve_id="EVE_1",
        scaffold="contig_1",
        start=0,
        end=50,
        confidence_tier="HIGH",
        status=VerificationStatus.HIGH_CONFIDENCE,
    )

    export_dir = generator.write_gvclass_export([result], tmp_path / "gvclass")
    records = list(SeqIO.parse(export_dir / "protein" / "EVE_1.faa", "fasta"))

    assert [record.id for record in records] == ["contig_1_VSR0123456789abcdef"]


def test_detailed_tsv_does_not_count_rescue_as_an_ordinary_protein(
    tmp_path: Path,
) -> None:
    genome_fasta = tmp_path / "genome.fna"
    proteome_fasta = tmp_path / "proteome.faa"
    confirmed_faa = tmp_path / "phase1" / "frameshift_screening" / "confirmed_frameshift_proteins.faa"
    confirmed_faa.parent.mkdir(parents=True)
    _write_genome_fasta(genome_fasta)
    _write_proteome_fasta(proteome_fasta)
    confirmed_faa.write_text(
        ">contig_1_VSR0123456789abcdef # 20 # 45 # 1 # "
        "ID=0_VSR0123456789abcdef;annotation=frameshift_rescued_domain\n"
        "MRESCUED\n"
    )
    confirmed_faa.with_name("confirmed_frameshift_markers.tsv").write_text(
        "query_porf\tscaffold\tstart\tend\tstrand\thmm_target\thmm_score\t"
        "validation_status\n"
        "contig_1_VSR0123456789abcdef|aa1-8\tcontig_1\t19\t45\t+\t"
        "VS000001\t80\tvalidated\n"
    )
    generator = OutputGenerator(
        output_dir=tmp_path,
        genome_fasta=genome_fasta,
        proteome_fasta=proteome_fasta,
        extended_output=True,
    )
    result = _build_result(
        eve_id="EVE_1",
        scaffold="contig_1",
        start=0,
        end=50,
        confidence_tier="HIGH",
        status=VerificationStatus.HIGH_CONFIDENCE,
    )

    detailed = generator.write_predictions_detailed_tsv([result])

    with detailed.open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["total_proteins"] == "1"
    assert row["hallmark_total"] == "1"


def test_write_eve_sequences_uses_coordinate_overlap(tmp_path: Path) -> None:
    genome_fasta = tmp_path / "genome.fna"
    proteome_fasta = tmp_path / "proteome.faa"
    _write_genome_fasta(genome_fasta)
    _write_proteome_fasta(proteome_fasta)

    generator = OutputGenerator(
        output_dir=tmp_path,
        genome_fasta=genome_fasta,
        proteome_fasta=proteome_fasta,
    )
    results = [
        _build_result(
            eve_id="EVE_1",
            scaffold="contig_1",
            start=0,
            end=50,
            confidence_tier="HIGH",
            status=VerificationStatus.HIGH_CONFIDENCE,
        ),
        _build_result(
            eve_id="EVE_2",
            scaffold="contig_1",
            start=60,
            end=120,
            confidence_tier="LOW",
            status=VerificationStatus.AMBIGUOUS,
        ),
    ]

    export_dir = generator.write_eve_sequences(results, tmp_path / "all_eves")

    first_records = list(SeqIO.parse(export_dir / "protein" / "EVE_1.faa", "fasta"))
    second_records = list(SeqIO.parse(export_dir / "protein" / "EVE_2.faa", "fasta"))

    assert [record.id for record in first_records] == ["contig_1_1"]
    assert [record.id for record in second_records] == ["contig_1_2"]


def test_same_coordinate_candidate_ids_remain_unique_across_exports(
    tmp_path: Path,
) -> None:
    genome_fasta = tmp_path / "genome.fna"
    proteome_fasta = tmp_path / "proteome.faa"
    _write_genome_fasta(genome_fasta)
    _write_proteome_fasta(proteome_fasta)
    generator = OutputGenerator(
        output_dir=tmp_path,
        genome_fasta=genome_fasta,
        proteome_fasta=proteome_fasta,
    )
    candidate_ids = [
        "EVE_contig_1_0-50-cordinary",
        "EVE_contig_1_0-50-crescue",
    ]
    results = [
        _build_result(
            eve_id=candidate_id,
            scaffold="contig_1",
            start=0,
            end=50,
            confidence_tier="HIGH",
            status=VerificationStatus.HIGH_CONFIDENCE,
        )
        for candidate_id in candidate_ids
    ]

    canonical_tsv = generator.write_predictions_tsv(results)
    detailed_tsv = generator.write_predictions_detailed_tsv(results)
    bed_path = generator.write_predictions_bed(results)
    gff_path = generator.write_predictions_gff(results)
    sequence_dir = generator.write_eve_sequences(results, tmp_path / "sequences")

    with canonical_tsv.open(newline="") as handle:
        canonical_ids = [row["eve_id"] for row in csv.DictReader(handle, delimiter="\t")]
    with detailed_tsv.open(newline="") as handle:
        detailed_ids = [row["eve_id"] for row in csv.DictReader(handle, delimiter="\t")]
    bed_ids = [line.split("\t")[3] for line in bed_path.read_text().splitlines()]
    gff_ids = [
        line.split("ID=", 1)[1].split(";", 1)[0]
        for line in gff_path.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    sequence_ids = [path.stem for path in (sequence_dir / "nucleotide").glob("*.fna")]

    assert set(canonical_ids) == set(candidate_ids)
    assert len(canonical_ids) == 2
    assert set(detailed_ids) == set(candidate_ids)
    assert len(detailed_ids) == 2
    assert set(bed_ids) == set(candidate_ids)
    assert len(bed_ids) == 2
    assert set(gff_ids) == set(candidate_ids)
    assert len(gff_ids) == 2
    assert set(sequence_ids) == set(candidate_ids)
    assert len(sequence_ids) == 2


@pytest.mark.parametrize("writer_name", ["write_gvclass_export", "write_eve_sequences"])
def test_per_eve_exports_encode_paths_but_preserve_raw_manifest_ids(
    tmp_path: Path,
    writer_name: str,
) -> None:
    genome_fasta = tmp_path / "genome.fna"
    _write_genome_fasta(genome_fasta)
    generator = OutputGenerator(output_dir=tmp_path, genome_fasta=genome_fasta)
    raw_ids = [
        "EVE_NODE/1",
        "../EVE_NODE",
        "EVE_NODE 1",
        "EVE_λ",
        "EVE_control\n",
    ]
    results = [
        _build_result(
            eve_id=raw_id,
            scaffold="contig_1",
            start=index,
            end=index + 20,
            confidence_tier="HIGH",
            status=VerificationStatus.HIGH_CONFIDENCE,
        )
        for index, raw_id in enumerate(raw_ids)
    ]
    export_dir = tmp_path / writer_name

    getattr(generator, writer_name)(results, export_dir)

    with (export_dir / "manifest.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["eve_id"] for row in rows] == raw_ids
    relative_paths = [row["nucleotide_fasta"] for row in rows]
    assert len(set(relative_paths)) == len(raw_ids)
    assert all(row["protein_fasta"] == "" for row in rows)
    for relative_path in relative_paths:
        assert relative_path
        exported = (export_dir / relative_path).resolve()
        assert export_dir.resolve() in exported.parents
        assert exported.is_file()
        assert "/" not in Path(relative_path).name.removesuffix(".fna")


def test_gene_taxonomy_uses_safe_filename_without_changing_raw_output_key(
    tmp_path: Path,
) -> None:
    generator = OutputGenerator(output_dir=tmp_path)
    raw_id = "EVE_NODE/1 λ"
    result = _build_result(
        eve_id=raw_id,
        scaffold="contig_1",
        start=0,
        end=20,
        confidence_tier="HIGH",
        status=VerificationStatus.HIGH_CONFIDENCE,
    )
    result.gene_taxonomy_records = [{"porf_id": "p1", "scaffold": "contig_1"}]

    output_files = generator.write_gene_taxonomy([result])

    output_path = output_files[f"gene_taxonomy_{raw_id}"]
    assert output_path.is_file()
    assert (tmp_path / "gene_taxonomy").resolve() in output_path.resolve().parents
    assert "/" not in output_path.name


@pytest.mark.parametrize(
    ("record", "expected_start", "expected_evalue", "expected_combined_score"),
    [
        pytest.param(
            {
                "porf_id": "phase2-current",
                "scaffold": "contig_1",
                "start": 0,
                "end": 100,
                "top1_prefix": "EUK",
                "top1_target": "current-hit",
                "top1_pident": 87.5,
                "top1_evalue": 1e-20,
            },
            "0",
            "1e-20",
            "87.5",
            id="phase2-current-coordinates",
        ),
        pytest.param(
            {
                "porf_id": "phase2-alias",
                "scaffold": "contig_1",
                "porf_start": 0,
                "porf_end": 100,
                "top1_prefix": "EUK",
                "top1_target": "alias-hit",
                "top1_pident": 76.5,
                "top1_evalue": 2e-30,
            },
            "0",
            "2e-30",
            "76.5",
            id="phase2-zero-coordinate-alias",
        ),
        pytest.param(
            {
                "porf_id": "legacy",
                "scaffold": "contig_1",
                "start": 0,
                "end": 100,
                "best_hit_origin": "EUK",
                "best_hit_target": "legacy-hit",
                "best_hit_evalue": 3e-40,
            },
            "0",
            "3e-40",
            "3e-40",
            id="legacy-evalue",
        ),
        pytest.param(
            {
                "porf_id": "hybrid-legacy-row",
                "top1_prefix": "EUK",
                "top1_target": "phase2-hit",
                "top1_pident": 66.5,
                "top1_evalue": 4e-50,
                "best_hit_origin": "LEGACY",
                "best_hit_target": "legacy-hit",
                "best_hit_evalue": 9e-9,
            },
            "",
            "9e-09",
            "66.5",
            id="legacy-per-eve-with-phase2-combined-fields",
        ),
    ],
)
def test_gene_taxonomy_preserves_metric_and_zero_coordinate_contract(
    tmp_path: Path,
    record: dict[str, object],
    expected_start: str,
    expected_evalue: str,
    expected_combined_score: str,
) -> None:
    result = _build_result(
        eve_id="EVE_METRIC",
        scaffold="contig_1",
        start=0,
        end=100,
        confidence_tier="HIGH",
        status=VerificationStatus.HIGH_CONFIDENCE,
    )
    result.gene_taxonomy_records = [record]

    output_files = OutputGenerator(output_dir=tmp_path).write_gene_taxonomy([result])

    with output_files["gene_taxonomy_EVE_METRIC"].open(encoding="utf-8", newline="") as handle:
        per_eve = next(csv.DictReader(handle, delimiter="\t"))
    with output_files["gene_taxonomy_all"].open(encoding="utf-8", newline="") as handle:
        combined = next(csv.DictReader(handle, delimiter="\t"))
    assert per_eve["start"] == expected_start
    assert combined["start"] == expected_start
    assert per_eve["best_hit_evalue"] == expected_evalue
    assert combined["best_hit_score"] == expected_combined_score


def test_per_eve_export_rejects_duplicate_raw_ids_before_writing(
    tmp_path: Path,
) -> None:
    generator = OutputGenerator(output_dir=tmp_path)
    results = [
        _build_result(
            eve_id="EVE_duplicate",
            scaffold="contig_1",
            start=index,
            end=index + 20,
            confidence_tier="HIGH",
            status=VerificationStatus.HIGH_CONFIDENCE,
        )
        for index in range(2)
    ]
    export_dir = tmp_path / "all_eves"

    with pytest.raises(ValueError, match="duplicate EVE ID"):
        generator.write_eve_sequences(results, export_dir)

    assert export_dir.exists() is False


def test_generate_all_counts_promoted_low_as_canonical_prediction(tmp_path: Path) -> None:
    generator = OutputGenerator(output_dir=tmp_path)
    result = _build_result(
        eve_id="EVE_LOW_PROMOTED",
        scaffold="contig_1",
        start=0,
        end=6001,
        confidence_tier="LOW",
        status=VerificationStatus.AMBIGUOUS,
        region_classification="UNKNOWN",
        likely_family="NCLDV",
        hallmark_count=2,
    )

    output_files = generator.generate_all([result])

    with Path(output_files["predictions_tsv"]).open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    summary = json.loads(Path(output_files["summary_json"]).read_text())

    assert len(rows) == 1
    assert rows[0]["confidence_tier"] == "LOW"
    assert rows[0]["effective_eve_class"] == "NCLDV"
    assert summary["statistics"]["total_candidates"] == 1
    assert summary["statistics"]["canonical_predictions"] == 1
    assert summary["statistics"]["promoted_low_confidence"] == 1
    assert summary["statistics"]["total_accepted_length_bp"] == 6001


def test_generate_all_drops_high_medium_that_fail_v2_gate(tmp_path: Path) -> None:
    generator = OutputGenerator(output_dir=tmp_path)
    result = _build_result(
        eve_id="EVE_HIGH_DROPPED",
        scaffold="contig_1",
        start=0,
        end=1000,
        confidence_tier="HIGH",
        status=VerificationStatus.HIGH_CONFIDENCE,
        region_classification="NCLDV",
        hallmark_count=0,
        has_mcp=False,
    )

    output_files = generator.generate_all([result])

    with Path(output_files["predictions_tsv"]).open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    summary = json.loads(Path(output_files["summary_json"]).read_text())

    assert rows == []
    assert summary["statistics"]["total_candidates"] == 1
    assert summary["statistics"]["canonical_predictions"] == 0
    assert summary["statistics"]["total_accepted_length_bp"] == 0


def test_a6_preselected_surface_reaches_canonical_output_unchanged(
    tmp_path: Path,
) -> None:
    generator = OutputGenerator(output_dir=tmp_path)
    normal_keep = _build_result(
        eve_id="EVE_NORMAL_KEEP",
        scaffold="contig_1",
        start=0,
        end=6001,
        confidence_tier="HIGH",
        status=VerificationStatus.HIGH_CONFIDENCE,
        region_classification="NCLDV",
        hallmark_count=1,
        has_mcp=True,
    )
    normal_drop = _build_result(
        eve_id="EVE_A6_RETAIN",
        scaffold="contig_1",
        start=7000,
        end=8000,
        confidence_tier="LOW",
        status=VerificationStatus.AMBIGUOUS,
        region_classification="UNKNOWN",
        hallmark_count=0,
        has_mcp=False,
    )
    selection = select_phase3_acceptance(
        [normal_keep, normal_drop],
        AblationID.A6,
    )

    output_files = generator.generate_all(
        list(selection.detailed_results),
        canonical_results=list(selection.canonical_results),
        promoted_low_results=list(selection.promoted_low_results),
    )

    with Path(output_files["predictions_tsv"]).open() as handle:
        canonical_rows = list(csv.DictReader(handle, delimiter="\t"))
    with Path(output_files["predictions_detailed_tsv"]).open() as handle:
        detailed_rows = list(csv.DictReader(handle, delimiter="\t"))
    summary = json.loads(Path(output_files["summary_json"]).read_text())
    bed_rows = Path(output_files["predictions_bed"]).read_text().splitlines()
    gff_rows = [
        line
        for line in Path(output_files["predictions_gff"]).read_text().splitlines()
        if line and not line.startswith("#")
    ]
    evidence = json.loads(Path(output_files["evidence_json"]).read_text())

    assert [row["eve_id"] for row in canonical_rows] == [
        "EVE_NORMAL_KEEP",
        "EVE_A6_RETAIN",
    ]
    assert [row["eve_id"] for row in detailed_rows] == [
        "EVE_NORMAL_KEEP",
        "EVE_A6_RETAIN",
    ]
    assert summary["statistics"]["total_candidates"] == 2
    assert summary["statistics"]["canonical_predictions"] == 2
    assert summary["statistics"]["low_confidence"] == 1
    assert summary["statistics"]["promoted_low_confidence"] == 0
    assert len(bed_rows) == 2
    assert len(gff_rows) == 2
    assert list(evidence) == ["EVE_NORMAL_KEEP", "EVE_A6_RETAIN"]


def test_a0_preselected_output_matches_generator_gate_fallback(
    tmp_path: Path,
) -> None:
    normal_keep = _build_result(
        eve_id="EVE_NORMAL_KEEP",
        scaffold="contig_1",
        start=0,
        end=6001,
        confidence_tier="HIGH",
        status=VerificationStatus.HIGH_CONFIDENCE,
        region_classification="NCLDV",
        hallmark_count=1,
        has_mcp=True,
    )
    promoted_low = _build_result(
        eve_id="EVE_LOW_PROMOTED",
        scaffold="contig_1",
        start=7000,
        end=14001,
        confidence_tier="LOW",
        status=VerificationStatus.AMBIGUOUS,
        region_classification="UNKNOWN",
        likely_family="VP",
        hallmark_count=2,
    )
    promoted_low.hallmark_genes = ["VP_MCP"]
    normal_drop = _build_result(
        eve_id="EVE_NORMAL_DROP",
        scaffold="contig_1",
        start=15000,
        end=16000,
        confidence_tier="HIGH",
        status=VerificationStatus.HIGH_CONFIDENCE,
        region_classification="NCLDV",
    )
    results = [normal_keep, promoted_low, normal_drop]
    selection = select_phase3_acceptance(results, AblationID.A0)

    fallback_files = OutputGenerator(output_dir=tmp_path / "fallback").generate_all(results)
    preselected_files = OutputGenerator(output_dir=tmp_path / "preselected").generate_all(
        list(selection.detailed_results),
        canonical_results=list(selection.canonical_results),
        promoted_low_results=list(selection.promoted_low_results),
    )

    assert [result.eve_id for result in selection.canonical_results] == [
        "EVE_NORMAL_KEEP",
        "EVE_LOW_PROMOTED",
    ]
    for key in ("predictions_tsv", "predictions_bed", "predictions_detailed_tsv"):
        assert Path(preselected_files[key]).read_bytes() == Path(fallback_files[key]).read_bytes()
    fallback_gff = [
        line for line in Path(fallback_files["predictions_gff"]).read_text().splitlines() if not line.startswith("#")
    ]
    preselected_gff = [
        line for line in Path(preselected_files["predictions_gff"]).read_text().splitlines() if not line.startswith("#")
    ]
    assert preselected_gff == fallback_gff
    assert json.loads(Path(preselected_files["evidence_json"]).read_text()) == json.loads(
        Path(fallback_files["evidence_json"]).read_text()
    )
    fallback_summary = json.loads(Path(fallback_files["summary_json"]).read_text())
    preselected_summary = json.loads(Path(preselected_files["summary_json"]).read_text())
    assert preselected_summary["statistics"] == fallback_summary["statistics"]
    assert preselected_summary["per_scaffold"] == fallback_summary["per_scaffold"]


def test_protein_counts_use_full_half_open_hit_overlap(tmp_path: Path) -> None:
    phase3_dir = tmp_path / "phase3"
    phase3_dir.mkdir()
    hits_path = tmp_path / "phase1" / "marker_validation" / "validated_marker_hits.tsv"
    hits_path.parent.mkdir(parents=True)
    hits_path.write_text(
        "query_porf\tscaffold\tstart\tend\tstrand\thmm_target\thmm_score\tvalidation_status\n"
        "left_touch\tctg\t10\t20\t+\tPLV_MCP_1\t100\tvalidated\n"
        "inside\tctg\t20\t25\t+\tVP_MCP_1\t100\tvalidated\n"
        "left_crossing\tctg\t15\t22\t+\tVP_Penton_1\t100\tvalidated\n"
        "right_touch\tctg\t30\t40\t+\tVP_ATPase_1\t100\tvalidated\n"
    )
    generator = OutputGenerator(output_dir=phase3_dir)

    observed = generator._load_protein_counts_by_region([("ctg", 20, 30)])

    assert observed == {
        ("ctg", 20, 30): {
            "VP_MCP": 1,
            "VP_Penton": 1,
        }
    }


def test_malformed_confirmed_table_does_not_discard_ordinary_marker_hits(
    tmp_path: Path,
) -> None:
    phase3_dir = tmp_path / "phase3"
    phase3_dir.mkdir()
    ordinary = tmp_path / "phase1" / "marker_validation" / "validated_marker_hits.tsv"
    ordinary.parent.mkdir(parents=True)
    ordinary.write_text(
        "query_porf\tscaffold\tstart\tend\tstrand\thmm_target\thmm_score\t"
        "validation_status\n"
        "ordinary\tctg\t20\t30\t+\tOG000001\t100\tvalidated\n"
    )
    confirmed = tmp_path / "phase1" / "frameshift_screening" / "confirmed_frameshift_markers.tsv"
    confirmed.parent.mkdir(parents=True)
    confirmed.write_text("bad_header\nvalue\n")
    generator = OutputGenerator(output_dir=phase3_dir)

    hits, protein_models = generator._parse_validated_marker_hits()

    assert hits["ctg"] == [(20, 30, "OG000001", "validated", "ordinary")]
    assert protein_models == {"ordinary": {"OG000001"}}


@pytest.mark.parametrize(
    "writer_name",
    ["write_predictions_tsv", "write_predictions_detailed_tsv"],
)
def test_prediction_tsv_writers_preserve_zero_candidate_start(
    tmp_path: Path,
    writer_name: str,
) -> None:
    result = _build_result(
        eve_id="EVE_ZERO",
        scaffold="ctg",
        start=0,
        end=20,
        confidence_tier="HIGH",
        status=VerificationStatus.HIGH_CONFIDENCE,
    )
    result.candidate_start = 0
    result.candidate_end = 20
    result.candidate_length = 20
    generator = OutputGenerator(output_dir=tmp_path, extended_output=False)

    output_path = getattr(generator, writer_name)([result])

    with output_path.open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["candidate_start"] == "0"
    assert row["candidate_end"] == "20"


@pytest.mark.parametrize(
    "writer_name",
    ["write_predictions_tsv", "write_predictions_detailed_tsv"],
)
@pytest.mark.parametrize("extended_output", [False, True])
def test_prediction_tsv_writers_publish_the_taxonomy_class(
    tmp_path: Path,
    writer_name: str,
    extended_output: bool,
) -> None:
    # The gate resolves this region to NCLDV, but publication reports the
    # taxonomy consensus, so the two must be able to disagree in the column.
    result = _build_result(
        eve_id="EVE_LOW_CONFLICT",
        scaffold="ctg",
        start=0,
        end=6001,
        confidence_tier="LOW",
        status=VerificationStatus.AMBIGUOUS,
        region_classification="PPV",
        likely_family="NCLDV",
        taxonomy_class="MIRUS",
        hallmark_count=2,
    )
    generator = OutputGenerator(output_dir=tmp_path, extended_output=extended_output)

    output_path = getattr(generator, writer_name)([result])

    with output_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    assert reader.fieldnames is not None
    if writer_name == "write_predictions_tsv":
        expected_fields = (
            CANONICAL_BASE_FIELDS
            + (CANONICAL_EXTENDED_FIELDS if extended_output else ("interproscan_score",))
            + ("effective_eve_class",)
        )
    else:
        expected_fields = tuple(
            column
            for column in DETAILED_PREDICTION_COLUMNS
            if extended_output or column not in DETAILED_PREDICTION_EXTENDED_COLUMNS
        )
    assert reader.fieldnames == list(expected_fields)
    assert rows[0]["effective_eve_class"] == "MIRUS"
    assert evaluate_v2_quality_gate(result).effective_class == "NCLDV"
    if writer_name == "write_predictions_tsv":
        assert rows[0]["region_classification"] == "PPV"
        assert rows[0]["classification"] == "NCLDV"
    else:
        assert rows[0]["likely_family"] == "NCLDV"

    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    empty_output = getattr(
        OutputGenerator(output_dir=empty_dir, extended_output=extended_output),
        writer_name,
    )([])
    empty_header = empty_output.read_text().splitlines()[0].split("\t")
    assert empty_header == list(expected_fields)


def test_canonical_tsv_uses_dot_for_unassigned_ppv_subtype(
    tmp_path: Path,
) -> None:
    result = _build_result(
        eve_id="EVE_PPV",
        scaffold="ctg",
        start=0,
        end=6001,
        confidence_tier="HIGH",
        status=VerificationStatus.HIGH_CONFIDENCE,
        region_classification="PPV",
        likely_family="PPV",
        taxonomy_class="PPV",
        hallmark_count=2,
    )
    result.ppv_subtype = ""

    output_path = OutputGenerator(
        output_dir=tmp_path,
        extended_output=True,
    ).write_predictions_tsv([result])

    with output_path.open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["effective_eve_class"] == "PPV"
    assert row["ppv_subtype"] == "."


def test_summary_records_coordinate_contract(tmp_path: Path) -> None:
    summary_path = OutputGenerator(output_dir=tmp_path).write_summary([])

    summary = json.loads(summary_path.read_text())

    assert summary["coordinate_schema_version"] == COORDINATE_SCHEMA_VERSION
    assert summary["output_schema_version"] == OUTPUT_SCHEMA_VERSION
    assert summary["coordinate_convention"] == COORDINATE_CONVENTION


def test_evidence_profile_serialization_does_not_follow_set_iteration_order() -> None:
    """Key and list order in evidence_profiles.json / evidence_graph.json.

    ``WindowEvidence.evidence_types`` is a set and ``Enum.__hash__`` is
    ``hash(self._name_)``, so anything that inherits its iteration order varies
    with PYTHONHASHSEED between runs.
    """
    etypes = [
        EvidenceType.HALLMARK_MCP,
        EvidenceType.HALLMARK_A32,
        EvidenceType.HIGH_KFD,
        EvidenceType.HIGH_NOVELTY,
        EvidenceType.CRF_CORE_VIRAL,
        EvidenceType.VIRAL_STRUCTURE,
        EvidenceType.ANOMALOUS_GC,
    ]
    profile = EvidenceProfile(eve_id="EVE_1", scaffold="scaffold", start=0, end=1000)
    for i in range(4):
        window = WindowEvidence(scaffold="scaffold", start=i * 250, end=(i + 1) * 250)
        for etype in etypes:
            window.add_evidence(etype)
        profile.windows.append(window)
    profile.compute_aggregates()

    analysis = CoherenceAnalysis(
        eve_id="EVE_1",
        profile=profile,
        graph=EvidenceCorrelationGraph(),
    )
    analysis.compute_coherence()
    emitted = analysis.to_dict()

    expected = sorted(e.value for e in etypes)
    assert list(emitted["evidence_profile"]["evidence_counts"]) == expected
    assert list(emitted["evidence_profile"]["evidence_coverage"]) == expected
    assert emitted["windows"][0]["evidence_types"] == expected
    assert emitted["graph_summary"]["evidence_types"] == expected

    # Every type co-occurs in every window, so all edge weights tie exactly.
    # Truncating to the strongest 5 must still keep the same 5 pairs each run.
    connections = emitted["graph_summary"]["strongest_connections"]
    assert [(u, v) for u, v, _ in connections] == list(itertools.combinations(expected, 2))[:5]


def test_merged_seed_sources_are_ordered_like_the_detailed_tsv_field() -> None:
    """``write_predictions_detailed_tsv`` emits ``sorted(r.seed_sources)`` while
    ``VerificationResult.to_dict`` emits the list as built, so the Phase 2 merge
    must build it sorted or the two files disagree for identical data.
    """
    boundaries = [
        RefinedBoundary(
            scaffold="scaffold",
            start=0,
            end=1000,
            original_end=1000,
            seed_sources=["novelty", "hhg", "compositional"],
        ),
        RefinedBoundary(
            scaffold="scaffold",
            start=500,
            end=1500,
            original_end=1500,
            seed_sources=["marker", "compositional"],
        ),
    ]

    merged = merge_adjacent_viral_boundaries(
        boundaries,
        taxonomy_map={},
        proteome_index={},
    )

    assert len(merged) == 1
    assert merged[0].seed_sources == ["compositional", "hhg", "marker", "novelty"]


@pytest.mark.parametrize("extended_output", [False, True])
def test_integration_evidence_preserves_boundaries_and_gene_provenance(tmp_path: Path, extended_output: bool) -> None:
    """Both public tables and JSON retain repeat and enzyme evidence."""
    result = VerificationResult(
        eve_id="EVE_contig_1_0-140",
        scaffold="contig_1",
        start=0,
        end=140,
        tir_present=True,
        tir_status="detected",
        tir_left_start=0,
        tir_left_end=30,
        tir_right_start=110,
        tir_right_end=140,
        tir_identity=0.9666666667,
        tir_alignment_capped=False,
        tir_alignment_length=30,
        tir_boundary_override=True,
        pre_tir_start=20,
        pre_tir_end=120,
        integration_gene_hits=[
            {
                "protein_id": "contig_1_1",
                "mechanism": "tyrosine_recombinase",
                "location": "interior",
                "source": "pfam_hmm",
            },
            {
                "protein_id": "contig_1_1",
                "mechanism": "tyrosine_recombinase",
                "location": "interior",
                "source": "marker_annotation",
            },
            {"protein_id": "contig_1_2", "mechanism": "dde_integrase", "location": "downstream", "source": "pfam_hmm"},
        ],
    )
    generator = OutputGenerator(output_dir=tmp_path, extended_output=extended_output)
    for writer in (generator.write_predictions_tsv, generator.write_predictions_detailed_tsv):
        output_path = writer([result])
        with output_path.open() as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        assert (row["start"], row["end"], row["length"]) == ("0", "140", "140")
        assert (row["tir_present"], row["tir_boundary_override"]) == ("1", "1")
        assert row["tir_alignment_capped"] == "0"
        assert row["tir_alignment_length"] == "30"
        assert row["tir_left_start"] == "0"
        assert (row["pre_tir_start"], row["pre_tir_end"]) == ("20", "120")
        assert row["recombinase_genes"] == "contig_1_1"
        assert json.loads(row["integration_gene_evidence"]) == result.integration_gene_hits
        assert row["integration_hmm_status"] == "not_assessed"
        assert row["integration_hmm_unsearched"] == "[]"
    profiles = json.loads(generator.write_evidence_profiles([result]).read_text())
    assert profiles[result.eve_id]["tir_left_start"] == 0
    assert profiles[result.eve_id]["integration_gene_hits"] == result.integration_gene_hits
    assert profiles[result.eve_id]["integration_hmm_status"] == "not_assessed"
    assert profiles[result.eve_id]["integration_hmm_unsearched"] == []


@pytest.mark.parametrize("independent_annotation", [False, True])
def test_incomplete_hmm_screen_preserves_acceptance_and_all_exports(
    tmp_path: Path, independent_annotation: bool
) -> None:
    """Unsearched proteins remain distinct from optional independent positive evidence."""
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">contig_1_1 # 1 # 300003 # -1 # ID=1_1;partial=00\n" + "M" * 100_001 + "\n", encoding="utf-8")
    original_bytes = proteome.read_bytes()
    result = _build_result(
        eve_id="EVE_contig_1_0-300003",
        scaffold="contig_1",
        start=0,
        end=300_003,
        confidence_tier="HIGH",
        status=VerificationStatus.HIGH_CONFIDENCE,
        hallmark_count=3,
        has_mcp=True,
    )
    assert result.integration_hmm_status == "not_assessed"
    assert result.integration_hmm_unsearched == []
    acceptance_before = select_phase3_acceptance([result], AblationID.A0)
    assert acceptance_before.canonical_results == (result,)
    before = result.to_dict()
    config = PipelineConfig()
    config.compute.threads = 1
    config.databases.hmm_database = ""
    config.phase3.interproscan_enabled = independent_annotation
    interpro = tmp_path / "phase3" / "interproscan" / "interproscan_batch.tsv"
    interpro.parent.mkdir(parents=True)
    interpro.write_text(
        "EVE_contig_1_0-300003|contig_1_1\tmd5\t100001\tPfam\tPF00589\tPhage integrase\t3\t170\t1E-30"
        "\tT\t2026-01-01\tIPR002104\tTyrosine recombinase\n",
        encoding="utf-8",
    )

    _annotate_integration_genes(
        [result], proteome_path=proteome, validated_markers=[], output_dir=tmp_path, config=config
    )

    assert result.integration_hmm_status == "incomplete_sequence_length"
    expected_unsearched = [
        {
            "region_id": result.eve_id,
            "protein_id": "contig_1_1",
            "length_aa": 100_001,
            "scaffold": "contig_1",
            "start": 0,
            "end": 300_003,
            "strand": "-",
            "location": "interior",
            "reason": "sequence_length_limit",
            "limit_aa": 100_000,
        }
    ]
    assert result.integration_hmm_unsearched == expected_unsearched
    assert len(result.integration_gene_hits) == int(independent_annotation)
    assert all(hit["source"] == "interproscan" for hit in result.integration_gene_hits)
    after = result.to_dict()
    changed = {key for key in after if after[key] != before[key]}
    assert changed <= {"integration_gene_hits", "integration_hmm_status", "integration_hmm_unsearched"}
    acceptance_after = select_phase3_acceptance([result], AblationID.A0)
    assert acceptance_after == acceptance_before
    generator = OutputGenerator(output_dir=tmp_path)
    accepted_path = generator.write_predictions_tsv(list(acceptance_after.canonical_results))
    detailed_path = generator.write_predictions_detailed_tsv([result])
    with accepted_path.open(encoding="utf-8") as handle:
        accepted = next(csv.DictReader(handle, delimiter="\t"))
    with detailed_path.open(encoding="utf-8") as handle:
        detailed = next(csv.DictReader(handle, delimiter="\t"))
    evidence = json.loads(generator.write_evidence_profiles([result]).read_text(encoding="utf-8"))[result.eve_id]
    assert (
        accepted["integration_hmm_status"] == detailed["integration_hmm_status"] == evidence["integration_hmm_status"]
    )
    assert json.loads(accepted["integration_hmm_unsearched"]) == expected_unsearched
    assert (
        json.loads(detailed["integration_hmm_unsearched"])
        == evidence["integration_hmm_unsearched"]
        == expected_unsearched
    )
    assert accepted["recombinase_genes"] == ("contig_1_1" if independent_annotation else ".")
    assert proteome.read_bytes() == original_bytes


def test_empty_integration_selection_exports_complete_status(tmp_path: Path) -> None:
    """An EVE with no selected proteins still has an assessed, empty screen."""
    proteome = tmp_path / "empty.faa"
    proteome.write_text(">other_1 # 1 # 90 # 1 # ID=1_1;partial=00\n" + "M" * 30 + "\n", encoding="utf-8")
    result = VerificationResult(eve_id="eve-1", scaffold="ctg", start=0, end=100)
    config = PipelineConfig()
    config.compute.threads = 1
    config.databases.hmm_database = ""

    _annotate_integration_genes(
        [result], proteome_path=proteome, validated_markers=[], output_dir=tmp_path, config=config
    )

    assert result.integration_hmm_status == "complete"
    assert result.integration_hmm_unsearched == []
    path = OutputGenerator(output_dir=tmp_path).write_predictions_detailed_tsv([result])
    with path.open(encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["integration_hmm_status"] == "complete"
    assert row["integration_hmm_unsearched"] == "[]"
