from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path

import pytest
from Bio import SeqIO

from virosync.ablation import AblationID
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
from virosync.pipeline.phase3.output_generator import OutputGenerator
from virosync.output_contract import (
    COORDINATE_CONVENTION,
    COORDINATE_SCHEMA_VERSION,
    DETAILED_PREDICTION_COLUMNS,
    DETAILED_PREDICTION_EXTENDED_COLUMNS,
    OUTPUT_SCHEMA_VERSION,
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
    path.write_text(
        ">contig_1\n"
        + ("A" * 160)
        + "\n>contig_10\n"
        + ("C" * 120)
        + "\n"
    )


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

    fallback_files = OutputGenerator(
        output_dir=tmp_path / "fallback"
    ).generate_all(results)
    preselected_files = OutputGenerator(
        output_dir=tmp_path / "preselected"
    ).generate_all(
        list(selection.detailed_results),
        canonical_results=list(selection.canonical_results),
        promoted_low_results=list(selection.promoted_low_results),
    )

    assert [result.eve_id for result in selection.canonical_results] == [
        "EVE_NORMAL_KEEP",
        "EVE_LOW_PROMOTED",
    ]
    for key in ("predictions_tsv", "predictions_bed", "predictions_detailed_tsv"):
        assert Path(preselected_files[key]).read_bytes() == Path(
            fallback_files[key]
        ).read_bytes()
    fallback_gff = [
        line
        for line in Path(fallback_files["predictions_gff"]).read_text().splitlines()
        if not line.startswith("#")
    ]
    preselected_gff = [
        line
        for line in Path(preselected_files["predictions_gff"]).read_text().splitlines()
        if not line.startswith("#")
    ]
    assert preselected_gff == fallback_gff
    assert json.loads(Path(preselected_files["evidence_json"]).read_text()) == json.loads(
        Path(fallback_files["evidence_json"]).read_text()
    )
    fallback_summary = json.loads(Path(fallback_files["summary_json"]).read_text())
    preselected_summary = json.loads(
        Path(preselected_files["summary_json"]).read_text()
    )
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
def test_prediction_tsv_writers_append_tier_aware_effective_class(
    tmp_path: Path,
    writer_name: str,
    extended_output: bool,
) -> None:
    result = _build_result(
        eve_id="EVE_LOW_CONFLICT",
        scaffold="ctg",
        start=0,
        end=6001,
        confidence_tier="LOW",
        status=VerificationStatus.AMBIGUOUS,
        region_classification="PPV",
        likely_family="NCLDV",
        hallmark_count=2,
    )
    generator = OutputGenerator(output_dir=tmp_path, extended_output=extended_output)

    output_path = getattr(generator, writer_name)([result])

    with output_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    assert reader.fieldnames is not None
    if writer_name == "write_predictions_tsv":
        expected_fields = CANONICAL_BASE_FIELDS + (
            CANONICAL_EXTENDED_FIELDS
            if extended_output
            else ("interproscan_score",)
        ) + ("effective_eve_class",)
    else:
        expected_fields = tuple(
            column
            for column in DETAILED_PREDICTION_COLUMNS
            if extended_output
            or column not in DETAILED_PREDICTION_EXTENDED_COLUMNS
        )
    assert reader.fieldnames == list(expected_fields)
    assert rows[0]["effective_eve_class"] == "NCLDV"
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
    assert [(u, v) for u, v, _ in connections] == list(
        itertools.combinations(expected, 2)
    )[:5]


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

    merged = merge_adjacent_viral_boundaries(boundaries, taxonomy_map={})

    assert len(merged) == 1
    assert merged[0].seed_sources == ["compositional", "hhg", "marker", "novelty"]
