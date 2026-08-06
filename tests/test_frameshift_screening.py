from __future__ import annotations

import subprocess
from dataclasses import replace
from pathlib import Path

import pytest
from Bio import SeqIO

from virosync.config import MaskingConfig
from virosync.orchestration._flows.single_genome import orchestrator
from virosync.pipeline.phase1 import frameshift_screening
from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit


def _hmm_record(name: str) -> str:
    return f"HMMER3/f [synthetic]\nNAME  {name}\nLENG  5\n//\n"


def _event_row(
    hit_id: str,
    ali_from: int,
    ali_to: int,
    shifts: int,
    stops: int,
    description: str = "synthetic hit",
) -> str:
    return (
        f"{hit_id} contig_1 - VS000001 - 100 2 90 1000 "
        f"{ali_from} {ali_to} 1e-20 75.5 0.1 42.0 {shifts} {stops} {description}\n"
    )


def _bath_alignment_report(
    *,
    query: str = "VS000001",
    target: str = "contig_1",
    ali_from: int = 10,
    ali_to: int = 40,
    amino_acids: str = "M * a - G",
) -> str:
    return (
        f"Query:       {query}  [M=100]\n"
        "Annotation for each hit (and alignments):\n"
        f">> {target}  synthetic target\n"
        "  Alignment:\n"
        f"                       {amino_acids}\n"
        f"  {target} {ali_from}    ATG  TGA  GCT  ---  GGT    {ali_to}\n"
        "//\n"
    )


def test_filter_vs_profiles_keeps_only_complete_exact_vs_records(tmp_path: Path) -> None:
    source = tmp_path / "combined.hmm"
    output = tmp_path / "filtered.hmm"
    source.write_text(
        _hmm_record("VS000001")
        + _hmm_record("NCLDV_MCP")
        + _hmm_record("VS123456")
        + _hmm_record("VS123456_extra")
    )

    count = frameshift_screening.filter_vs_profiles(source, output)

    assert count == 2
    assert output.read_text() == _hmm_record("VS000001") + _hmm_record("VS123456")


def test_filter_vs_profiles_rejects_duplicate_and_zero_selection(
    tmp_path: Path,
) -> None:
    duplicate = tmp_path / "duplicate.hmm"
    duplicate.write_text(_hmm_record("VS000001") * 2)
    with pytest.raises(ValueError, match="Duplicate VS profile NAME.*VS000001"):
        frameshift_screening.filter_vs_profiles(
            duplicate,
            tmp_path / "duplicate.filtered.hmm",
        )

    no_vs = tmp_path / "no-vs.hmm"
    no_vs.write_text(_hmm_record("NCLDV_MCP"))
    with pytest.raises(ValueError, match="No VS###### profiles"):
        frameshift_screening.filter_vs_profiles(
            no_vs,
            tmp_path / "no-vs.filtered.hmm",
        )

    incomplete = tmp_path / "incomplete.hmm"
    incomplete.write_text(_hmm_record("VS000001").removesuffix("//\n"))
    with pytest.raises(ValueError, match="Incomplete HMMER record"):
        frameshift_screening.filter_vs_profiles(
            incomplete,
            tmp_path / "incomplete.filtered.hmm",
        )


def test_parser_normalizes_coordinates_strand_and_filters_event_free_rows(
    tmp_path: Path,
) -> None:
    tblout = tmp_path / "bathsearch.tblout"
    normalized = tmp_path / "frameshift_hits.tsv"
    tblout.write_text(
        "# BATH output\n"
        + _event_row("forward", 10, 40, shifts=1, stops=0)
        + _event_row("reverse", 90, 50, shifts=0, stops=2, description="-")
        + _event_row("empty-description", 200, 250, shifts=1, stops=0, description="")
        + _event_row("intact", 100, 150, shifts=0, stops=0)
    )

    hits = frameshift_screening.parse_bathsearch_tblout(tblout)
    frameshift_screening.write_frameshift_hits(hits, normalized)

    assert [hit.hit_id for hit in hits] == [
        "forward",
        "reverse",
        "empty-description",
    ]
    assert (hits[0].ali_start, hits[0].ali_end, hits[0].strand) == (9, 40, "+")
    assert (hits[1].ali_start, hits[1].ali_end, hits[1].strand) == (49, 90, "-")
    assert hits[1].description == "-"
    assert hits[2].description == ""
    assert {hit.annotation_class for hit in hits} == {
        frameshift_screening.ANNOTATION_CLASS
    }
    lines = normalized.read_text().splitlines()
    assert lines[0].split("\t") == [
        "annotation_class",
        "hit_id",
        "target_name",
        "target_accession",
        "query_name",
        "query_accession",
        "hmm_len",
        "hmm_from",
        "hmm_to",
        "seq_len",
        "ali_start",
        "ali_end",
        "strand",
        "evalue",
        "score",
        "bias",
        "pid",
        "shifts",
        "stops",
        "description",
    ]
    assert len(lines) == 4


def test_parser_rejects_malformed_rows(tmp_path: Path) -> None:
    tblout = tmp_path / "malformed.tblout"
    tblout.write_text("hit with too few columns\n")

    with pytest.raises(ValueError, match="expected 18 columns"):
        frameshift_screening.parse_bathsearch_tblout(tblout)

    invalid_numeric = tmp_path / "invalid-numeric.tblout"
    invalid_numeric.write_text(
        _event_row("invalid", 10, 40, shifts=1, stops=0).replace(
            " 100 2 90 1000 ",
            " invalid 2 90 1000 ",
        )
    )
    with pytest.raises(ValueError, match="Invalid numeric value"):
        frameshift_screening.parse_bathsearch_tblout(invalid_numeric)


def test_parser_matches_real_bath_tblout_fixture() -> None:
    fixture = Path(__file__).parent / "data" / "bathsearch_vs000001.tblout"
    header = fixture.read_text().splitlines()[0].removeprefix("#").split()

    assert header[:4] == ["hit", "ID", "target", "name"]
    assert header[-4:] == ["stops", "description", "of", "target"]
    hits = frameshift_screening.parse_bathsearch_tblout(fixture)
    assert len(hits) == 1
    hit = hits[0]
    assert (
        hit.target_name,
        hit.query_name,
        hit.hmm_len,
        hit.hmm_from,
        hit.hmm_to,
        hit.seq_len,
        hit.ali_start,
        hit.ali_end,
        hit.strand,
        hit.evalue,
        hit.score,
        hit.bias,
        hit.pid,
        hit.shifts,
        hit.stops,
    ) == (
        "VS000001_consensus_plus_one_base",
        "VS000001",
        179,
        1,
        178,
        538,
        0,
        533,
        "+",
        1.2e-85,
        272.6,
        0.2,
        99.44,
        2,
        0,
    )


def test_translation_parser_and_candidate_faa_preserve_coordinates_and_provenance(
    tmp_path: Path,
) -> None:
    report = tmp_path / "bathsearch.txt"
    report.write_text(_bath_alignment_report())
    translations = frameshift_screening.parse_bathsearch_translations(report)
    hit = frameshift_screening.parse_bathsearch_tblout(
        _write_text(
            tmp_path / "bathsearch.tblout",
            _event_row("1", 10, 40, shifts=1, stops=1),
        )
    )[0]

    candidate_faa = tmp_path / "frameshift_candidates.faa"
    indexed = frameshift_screening.write_frameshift_candidate_faa(
        [hit],
        translations,
        candidate_faa,
    )

    record = next(SeqIO.parse(candidate_faa, "fasta"))
    assert str(record.seq) == "MXAG"
    assert record.id in indexed
    assert record.id.startswith("contig_1_VSR")
    assert "annotation=frameshift_rescued_domain" in record.description
    assert "literal_stops=1" in record.description
    assert frameshift_screening.is_rescued_protein_id(record.id)
    assert frameshift_screening.is_rescued_protein_id(f"{record.id}|aa1-4")
    assert not frameshift_screening.is_rescued_protein_id("contig_1_VSRnot-generated")
    assert not frameshift_screening.is_rescued_protein_id(
        "contig_1_VSR0123456789ABCDEF"
    )


def _write_text(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def _validated_marker(hit: frameshift_screening.FrameshiftHit) -> ValidatedMarkerHit:
    return ValidatedMarkerHit(
        query_porf=f"{frameshift_screening.rescued_protein_id(hit)}|aa1-40",
        scaffold=hit.target_name,
        start=hit.ali_start,
        end=hit.ali_end,
        strand=hit.strand,
        hmm_target=hit.query_name,
        hmm_score=hit.score,
        hmm_evalue=hit.evalue,
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


def test_confirmed_rescue_selection_requires_coverage_and_deduplicates_loci(
    tmp_path: Path,
) -> None:
    base = frameshift_screening.parse_bathsearch_tblout(
        _write_text(
            tmp_path / "hits.tblout",
            _event_row("1", 10, 100, shifts=1, stops=0),
        )
    )[0]
    best = replace(base, score=90.0)
    overlapping = replace(base, hit_id="2", ali_start=20, ali_end=90, score=70.0)
    low_model_coverage = replace(
        base,
        hit_id="3",
        ali_start=200,
        ali_end=260,
        hmm_from=1,
        hmm_to=30,
        score=85.0,
    )
    low_query_coverage = replace(
        base,
        hit_id="4",
        ali_start=300,
        ali_end=390,
        score=88.0,
    )
    novel_only = replace(
        base,
        hit_id="5",
        ali_start=400,
        ali_end=490,
        score=89.0,
    )
    hits = [best, overlapping, low_model_coverage, low_query_coverage, novel_only]
    markers = [_validated_marker(hit) for hit in hits]
    markers[-1] = replace(markers[-1], validation_status="validated_novel")
    diamond = tmp_path / "diamond_top10.tsv"
    diamond.write_text(
        "".join(
            f"{frameshift_screening.rescued_protein_id(hit)}\tNCLDV__ref\t"
            f"1e-20\t80\t35\t{40 if hit is low_query_coverage else 75}\n"
            for hit in hits
        )
    )

    confirmed = frameshift_screening.select_confirmed_frameshift_markers(
        markers,
        {frameshift_screening.rescued_protein_id(hit): hit for hit in hits},
        diamond,
        validated_prefixes={"NCLDV__"},
        min_pident=25.0,
    )

    assert [marker.hmm_score for marker in confirmed] == [90.0]


def test_confirmed_marker_write_preserves_previous_file_on_failure(
    tmp_path: Path,
) -> None:
    output = tmp_path / "confirmed_frameshift_markers.tsv"
    output.write_text("previous complete artifact\n")

    class BrokenMarker:
        query_porf = "contig_1_VSR0123456789abcdef"

        @property
        def scaffold(self) -> str:
            raise RuntimeError("simulated write failure")

    with pytest.raises(RuntimeError, match="simulated write failure"):
        frameshift_screening.write_confirmed_frameshift_markers(
            [BrokenMarker()],
            output,
        )

    assert output.read_text() == "previous complete artifact\n"
    assert not list(tmp_path.glob(".tmp_confirmed_frameshift_markers.tsv_*"))


def test_runner_uses_frameshift_and_threshold_options_and_keeps_raw_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    hmm_database = tmp_path / "combined.hmm"
    masked_fasta = tmp_path / "masked.fna"
    output_dir = tmp_path / "output"
    hmm_database.write_text(_hmm_record("VS000001") + _hmm_record("NCLDV_MCP"))
    for suffix in ("h3f", "h3i", "h3m", "h3p"):
        Path(f"{hmm_database}.{suffix}").write_text("stale pressed sidecar\n")
    masked_fasta.write_text(">contig_1\nACGT\n")
    commands: list[tuple[list[str], dict]] = []
    converted_input = ""

    monkeypatch.setattr(
        frameshift_screening.shutil,
        "which",
        lambda name: f"/mock/bin/{name}",
    )

    def fake_run(args: list[str], **kwargs):
        nonlocal converted_input
        commands.append((args, kwargs))
        if Path(args[0]).name == "bathconvert":
            assert Path(args[2]) != hmm_database
            assert not any(
                Path(f"{args[2]}.{suffix}").exists()
                for suffix in ("h3f", "h3i", "h3m", "h3p")
            )
            converted_input = Path(args[2]).read_text()
            Path(args[1]).write_text("converted models\n")
        else:
            Path(args[args.index("--tblout") + 1]).write_text(
                _event_row("1", 10, 40, shifts=1, stops=0)
            )
            Path(args[args.index("--fstblout") + 1]).write_text("raw fs output\n")
            Path(args[args.index("-o") + 1]).write_text(_bath_alignment_report())
        return subprocess.CompletedProcess(args, 0, stdout="captured", stderr="")

    monkeypatch.setattr(frameshift_screening.subprocess, "run", fake_run)

    hits = frameshift_screening.run_frameshift_screening(
        masked_fasta=masked_fasta,
        hmm_database=hmm_database,
        output_dir=output_dir,
        threads=6,
    )

    assert [hit.hit_id for hit in hits] == ["1"]
    assert converted_input == _hmm_record("VS000001")
    assert len(commands) == 2
    bathsearch_args, bathsearch_kwargs = commands[1]
    assert bathsearch_args.count("--fs") == 1
    assert bathsearch_args[bathsearch_args.index("--cpu") + 1] == "6"
    assert bathsearch_args[bathsearch_args.index("-E") + 1] == "1e-5"
    assert bathsearch_args[bathsearch_args.index("--incE") + 1] == "1e-5"
    assert bathsearch_kwargs == {
        "capture_output": True,
        "text": True,
        "check": False,
        "timeout": 3600,
    }
    assert (output_dir / "bathsearch.tblout").is_file()
    assert (output_dir / "bathsearch.fstblout").is_file()
    assert (output_dir / "bathsearch.txt").is_file()
    assert (output_dir / "frameshift_hits.tsv").is_file()
    assert (output_dir / "frameshift_candidates.faa").is_file()
    assert not (output_dir / "vs_profiles.bhmm").exists()


def test_runner_fails_clearly_for_missing_tool_and_command_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        frameshift_screening.shutil,
        "which",
        lambda name: "/mock/bin/bathconvert" if name == "bathconvert" else None,
    )
    with pytest.raises(RuntimeError, match="requires 'bathsearch' on PATH"):
        frameshift_screening.run_frameshift_screening(
            tmp_path / "masked.fna",
            tmp_path / "combined.hmm",
            tmp_path / "missing-tool-output",
            threads=1,
        )

    hmm_database = tmp_path / "combined.hmm"
    hmm_database.write_text(_hmm_record("VS000001"))
    monkeypatch.setattr(
        frameshift_screening.shutil,
        "which",
        lambda name: f"/mock/bin/{name}",
    )
    monkeypatch.setattr(
        frameshift_screening.subprocess,
        "run",
        lambda args, **kwargs: subprocess.CompletedProcess(
            args,
            7,
            stdout="",
            stderr="conversion failed",
        ),
    )
    with pytest.raises(RuntimeError, match="bathconvert failed with exit code 7.*conversion failed"):
        frameshift_screening.run_frameshift_screening(
            tmp_path / "masked.fna",
            hmm_database,
            tmp_path / "failed-command-output",
            threads=1,
        )

    def timeout(args, **kwargs):
        raise subprocess.TimeoutExpired(args, kwargs["timeout"])

    monkeypatch.setattr(frameshift_screening.subprocess, "run", timeout)
    with pytest.raises(RuntimeError, match="bathconvert timed out after 3600 seconds"):
        frameshift_screening.run_frameshift_screening(
            tmp_path / "masked.fna",
            hmm_database,
            tmp_path / "timed-out-command-output",
            threads=1,
        )


def test_frameshift_optional_runtime_identities_and_phase_artifact(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    requested: list[str] = []

    def fake_identity(name: str):
        requested.append(name)
        return None

    monkeypatch.setattr(orchestrator, "_executable_resource_identity", fake_identity)
    orchestrator._enabled_executable_identities(
        {"frameshift_screening_enabled": False},
        MaskingConfig(),
    )
    assert "bathconvert" not in requested
    assert "bathsearch" not in requested

    requested.clear()
    orchestrator._enabled_executable_identities(
        {"frameshift_screening_enabled": True},
        MaskingConfig(),
    )
    assert "bathconvert" in requested
    assert "bathsearch" in requested

    diagnostic = tmp_path / "phase1" / "frameshift_screening" / "frameshift_hits.tsv"
    diagnostic.parent.mkdir(parents=True)
    diagnostic.write_text("annotation_class\n")
    confirmed = diagnostic.parent / "confirmed_frameshift_proteins.faa"
    confirmed.write_text(
        ">contig_1_VSR0123456789abcdef # 1 # 30 # 1 # "
        "ID=0_VSR0123456789abcdef;annotation=frameshift_rescued_domain\nMPEPTIDE\n"
    )
    confirmed.with_name("confirmed_frameshift_markers.tsv").write_text(
        "query_porf\tscaffold\tstart\tend\tstrand\thmm_target\thmm_score\t"
        "validation_status\n"
    )
    identities = orchestrator._phase_artifacts(tmp_path, 1)
    assert {identity.relative_path for identity in identities} == {
        "phase1/frameshift_screening/confirmed_frameshift_markers.tsv",
        "phase1/frameshift_screening/confirmed_frameshift_proteins.faa",
        "phase1/frameshift_screening/frameshift_hits.tsv",
    }
    assert {identity.schema for identity in identities} == {
        "frameshift-hits-v1",
        "frameshift-rescued-markers-v1",
        "frameshift-rescued-proteins-v1",
    }
