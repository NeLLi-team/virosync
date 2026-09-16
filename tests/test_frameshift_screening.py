from __future__ import annotations

import csv
from collections import Counter
from dataclasses import replace
from pathlib import Path

import pyhmmer
import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from virosync.config import MaskingConfig
from virosync.orchestration._flows.single_genome import orchestrator
from virosync.pipeline.phase1 import frameshift_screening
from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit

_PROTEIN = "ACDEFGHIKLMNPQRSTVWY" * 4
_CODONS = {
    "A": "GCT",
    "C": "TGT",
    "D": "GAT",
    "E": "GAA",
    "F": "TTT",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "K": "AAA",
    "L": "CTG",
    "M": "ATG",
    "N": "AAT",
    "P": "CCT",
    "Q": "CAA",
    "R": "CGT",
    "S": "TCT",
    "T": "ACT",
    "V": "GTT",
    "W": "TGG",
    "Y": "TAT",
}


def _build_hmm(name: str, protein: str = _PROTEIN) -> pyhmmer.plan7.HMM:
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequence = pyhmmer.easel.TextSequence(
        name=name.encode(),
        sequence=protein,
    ).digitize(alphabet)
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(
        sequence,
        pyhmmer.plan7.Background(alphabet),
    )
    return hmm


def _write_hmms(path: Path, *hmms: pyhmmer.plan7.HMM) -> Path:
    with path.open("wb") as handle:
        for hmm in hmms:
            hmm.write(handle)
    return path


def _coding_dna(protein: str = _PROTEIN) -> str:
    return "".join(_CODONS[amino_acid] for amino_acid in protein)


def _write_fasta(path: Path, records: dict[str, str]) -> Path:
    path.write_text(
        "".join(f">{name}\n{sequence}\n" for name, sequence in records.items()),
        encoding="utf-8",
    )
    return path


def _run_screen(
    tmp_path: Path,
    records: dict[str, str],
) -> tuple[list[frameshift_screening.FrameshiftHit], Path]:
    hmm_path = _write_hmms(tmp_path / "profiles.hmm", _build_hmm("VS000001"))
    fasta_path = _write_fasta(tmp_path / "genomes.fna", records)
    output_dir = tmp_path / "output"
    hits = frameshift_screening.run_frameshift_screening(
        fasta_path,
        hmm_path,
        output_dir,
        threads=2,
    )
    return hits, output_dir


def test_native_screen_recovers_supported_events_and_rejects_controls(
    tmp_path: Path,
) -> None:
    coding_dna = _coding_dna()
    insertion = coding_dna[:120] + "A" + coding_dna[120:]
    codon_deletion = coding_dna[:120] + coding_dna[123:]
    stop = coding_dna[:120] + "TAA" + coding_dna[123:]
    reverse = str(Seq(insertion).reverse_complement())
    padding = "N" * 60
    records = {
        "plus": padding + insertion + padding,
        "minus": padding + reverse + padding,
        "intact": padding + coding_dna + padding,
        "codon_deletion": padding + codon_deletion + padding,
        "stop": padding + stop + padding,
        "multiple": (padding + insertion + "N" * 90 + coding_dna + "N" * 90 + insertion + padding),
    }

    hits, output_dir = _run_screen(tmp_path, records)

    assert Counter(hit.target_name for hit in hits) == {
        "plus": 1,
        "minus": 1,
        "stop": 1,
        "multiple": 2,
    }
    assert all(hit.backend == "native_codon" for hit in hits)
    assert all(hit.pid is None for hit in hits)
    assert all(hit.evalue <= frameshift_screening.PEPTIDE_MAX_EVALUE for hit in hits)
    assert all((hit.hmm_to - hit.hmm_from + 1) / hit.hmm_len >= frameshift_screening.MIN_MODEL_COVERAGE for hit in hits)
    plus_hit = next(hit for hit in hits if hit.target_name == "plus")
    minus_hit = next(hit for hit in hits if hit.target_name == "minus")
    stop_hit = next(hit for hit in hits if hit.target_name == "stop")
    assert (plus_hit.ali_start, plus_hit.ali_end, plus_hit.strand) == (60, 301, "+")
    assert (minus_hit.ali_start, minus_hit.ali_end, minus_hit.strand) == (60, 301, "-")
    assert plus_hit.native_score - plus_hit.no_frameshift_score >= 6.0
    assert (stop_hit.shifts, stop_hit.stops) == (0, 1)
    assert stop_hit.native_score == pytest.approx(stop_hit.no_frameshift_score)

    with (output_dir / "frameshift_events.tsv").open(
        newline="",
        encoding="utf-8",
    ) as handle:
        events = list(csv.DictReader(handle, delimiter="\t"))
    assert len(events) == 5
    plus_event = next(event for event in events if event["target_name"] == "plus")
    minus_event = next(event for event in events if event["target_name"] == "minus")
    assert (plus_event["event_kind"], plus_event["genomic_start"], plus_event["genomic_end"]) == (
        "insertion",
        "177",
        "181",
    )
    assert (minus_event["event_kind"], minus_event["genomic_start"], minus_event["genomic_end"]) == (
        "insertion",
        "180",
        "184",
    )
    candidate_records = list(SeqIO.parse(output_dir / "frameshift_candidates.faa", "fasta"))
    assert len(candidate_records) == 5
    assert all("annotation=frameshift_rescued_domain" in record.description for record in candidate_records)
    assert all(frameshift_screening.is_rescued_protein_id(record.id) for record in candidate_records)


def test_native_screen_writes_header_only_empty_artifacts(tmp_path: Path) -> None:
    padding = "N" * 60
    hits, output_dir = _run_screen(
        tmp_path,
        {"intact": padding + _coding_dna() + padding},
    )

    assert hits == []
    hit_lines = (output_dir / "frameshift_hits.tsv").read_text(encoding="utf-8").splitlines()
    event_lines = (output_dir / "frameshift_events.tsv").read_text(encoding="utf-8").splitlines()
    assert len(hit_lines) == 1
    assert hit_lines[0].endswith("backend\tnative_score\tno_frameshift_score\tevalue_interpretation")
    assert event_lines == [
        "candidate_id\ttarget_name\tquery_name\tstrand\tevent_kind\t"
        "genomic_start\tgenomic_end\tsize\tmodel_position\tprotein_position"
    ]
    assert (output_dir / "frameshift_candidates.faa").read_text(encoding="utf-8") == ""


def test_seed_search_bounds_long_stop_rich_targets_and_finds_boundary_domain() -> None:
    protein = _PROTEIN * 3
    coding_start = 2_900
    coding_end = coding_start + len(protein) * 3
    dna = "TAA" * 4_000
    dna = dna[:coding_start] + _coding_dna(protein) + dna[coding_end:]
    chunks = frameshift_screening._translate_chunks({"long_contig": dna})

    assert max(chunk.chunk_end - chunk.chunk_start for chunk in chunks) <= 3_000
    seeds = frameshift_screening._search_seeds(
        [_build_hmm("VS000001", protein)],
        chunks,
        threads=1,
    )

    assert any(
        seed.target_name == "long_contig" and seed.strand == "+" and seed.start < coding_end and coding_start < seed.end
        for seed in seeds
    )


def test_profile_loader_validates_exact_names_duplicates_and_malformed_input(
    tmp_path: Path,
) -> None:
    selected = _write_hmms(
        tmp_path / "selected.hmm",
        _build_hmm("OTHER"),
        _build_hmm("VS000001"),
    )
    assert [hmm.name for hmm in frameshift_screening._load_vs_profiles(selected)] == [b"VS000001"]

    duplicate = _write_hmms(
        tmp_path / "duplicate.hmm",
        _build_hmm("VS000001"),
        _build_hmm("VS000001"),
    )
    with pytest.raises(ValueError, match="Duplicate VS profile NAME.*VS000001"):
        frameshift_screening._load_vs_profiles(duplicate)

    no_models = _write_hmms(
        tmp_path / "no-vs.hmm",
        _build_hmm("OTHER"),
    )
    with pytest.raises(ValueError, match="No VS###### profiles"):
        frameshift_screening._load_vs_profiles(no_models)

    malformed = tmp_path / "malformed.hmm"
    malformed.write_text("not an HMM\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Malformed HMM database"):
        frameshift_screening._load_vs_profiles(malformed)


def test_event_support_requires_observed_residues_on_both_flanks() -> None:
    event = frameshift_screening.FrameshiftEvent(
        kind="stop",
        nt_start=30,
        nt_end=33,
        model_position=10,
        size=0,
        protein_position=10,
    )
    domain = frameshift_screening._PeptideDomain(
        sequence_score=50.0,
        bias=0.0,
        evalue=1e-20,
        hmm_from=1,
        hmm_to=20,
        target_start=0,
        target_end=20,
    )

    assert frameshift_screening._event_is_internal(
        event,
        domain,
        "A" * 10 + "X" + "A" * 9,
        0,
        60,
    )
    assert not frameshift_screening._event_is_internal(
        event,
        domain,
        "AXXXXXXXXA" + "X" + "A" * 9,
        0,
        60,
    )


def _hit(
    *,
    hit_id: str = "one",
    start: int = 9,
    end: int = 100,
    score: float = 75.5,
    hmm_from: int = 1,
    hmm_to: int = 90,
) -> frameshift_screening.FrameshiftHit:
    return frameshift_screening.FrameshiftHit(
        annotation_class=frameshift_screening.ANNOTATION_CLASS,
        hit_id=hit_id,
        target_name="contig_1",
        target_accession="-",
        query_name="VS000001",
        query_accession="-",
        hmm_len=100,
        hmm_from=hmm_from,
        hmm_to=hmm_to,
        seq_len=1_000,
        ali_start=start,
        ali_end=end,
        strand="+",
        evalue=1e-20,
        score=score,
        bias=0.1,
        pid=0.0,
        shifts=1,
        stops=0,
        description="synthetic native hit",
        native_score=80.0,
        no_frameshift_score=60.0,
    )


def _validated_marker(
    hit: frameshift_screening.FrameshiftHit,
) -> ValidatedMarkerHit:
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


def test_confirmed_rescue_selection_and_writers_preserve_contract(
    tmp_path: Path,
) -> None:
    base = _hit(score=90.0)
    overlapping = replace(
        base,
        hit_id="two",
        ali_start=20,
        ali_end=90,
        score=70.0,
    )
    low_model_coverage = replace(
        base,
        hit_id="three",
        ali_start=200,
        ali_end=260,
        hmm_to=30,
        score=85.0,
    )
    low_query_coverage = replace(
        base,
        hit_id="four",
        ali_start=300,
        ali_end=390,
        score=88.0,
    )
    novel_only = replace(
        base,
        hit_id="five",
        ali_start=400,
        ali_end=490,
        score=89.0,
    )
    hits = [base, overlapping, low_model_coverage, low_query_coverage, novel_only]
    markers = [_validated_marker(hit) for hit in hits]
    markers[-1] = replace(markers[-1], validation_status="validated_novel")
    diamond = tmp_path / "diamond_top10.tsv"
    diamond.write_text(
        "".join(
            f"{frameshift_screening.rescued_protein_id(hit)}\tNCLDV__ref\t"
            f"1e-20\t80\t35\t{40 if hit is low_query_coverage else 75}\n"
            for hit in hits
        ),
        encoding="utf-8",
    )

    confirmed = frameshift_screening.select_confirmed_frameshift_markers(
        markers,
        {frameshift_screening.rescued_protein_id(hit): hit for hit in hits},
        diamond,
        validated_prefixes={"NCLDV__"},
        min_pident=25.0,
    )

    assert [marker.hmm_score for marker in confirmed] == [90.0]
    candidate_faa = tmp_path / "frameshift_candidates.faa"
    frameshift_screening.write_frameshift_candidate_faa(
        hits,
        {frameshift_screening.rescued_protein_id(hit): "MPEPTIDE" for hit in hits},
        candidate_faa,
    )
    confirmed_faa = tmp_path / "confirmed_frameshift_proteins.faa"
    assert (
        frameshift_screening.write_confirmed_frameshift_faa(
            candidate_faa,
            confirmed,
            confirmed_faa,
        )
        == 1
    )
    assert [record.id for record in SeqIO.parse(confirmed_faa, "fasta")] == [
        frameshift_screening.rescued_protein_id(base)
    ]


def test_confirmed_marker_write_preserves_previous_file_on_failure(
    tmp_path: Path,
) -> None:
    output = tmp_path / "confirmed_frameshift_markers.tsv"
    output.write_text("previous complete artifact\n", encoding="utf-8")

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

    assert output.read_text(encoding="utf-8") == "previous complete artifact\n"
    assert not list(tmp_path.glob(".tmp_confirmed_frameshift_markers.tsv_*"))


def test_frameshift_runtime_has_no_bath_identity_and_tracks_native_artifacts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    requested: list[str] = []

    def fake_identity(name: str):
        requested.append(name)
        return None

    monkeypatch.setattr(orchestrator, "_executable_resource_identity", fake_identity)
    orchestrator._enabled_executable_identities(
        {"frameshift_screening_enabled": True},
        MaskingConfig(),
    )
    assert "bathconvert" not in requested
    assert "bathsearch" not in requested

    diagnostic = tmp_path / "phase1" / "frameshift_screening" / "frameshift_hits.tsv"
    diagnostic.parent.mkdir(parents=True)
    diagnostic.write_text("annotation_class\n", encoding="utf-8")
    diagnostic.with_name("frameshift_events.tsv").write_text(
        "candidate_id\n",
        encoding="utf-8",
    )
    confirmed = diagnostic.parent / "confirmed_frameshift_proteins.faa"
    confirmed.write_text(
        ">contig_1_VSR0123456789abcdef # 1 # 30 # 1 # "
        "ID=0_VSR0123456789abcdef;annotation=frameshift_rescued_domain\nMPEPTIDE\n",
        encoding="utf-8",
    )
    confirmed.with_name("confirmed_frameshift_markers.tsv").write_text(
        "query_porf\tscaffold\tstart\tend\tstrand\thmm_target\thmm_score\tvalidation_status\n",
        encoding="utf-8",
    )

    identities = orchestrator._phase_artifacts(tmp_path, 1)

    assert {identity.relative_path for identity in identities} == {
        "phase1/frameshift_screening/confirmed_frameshift_markers.tsv",
        "phase1/frameshift_screening/confirmed_frameshift_proteins.faa",
        "phase1/frameshift_screening/frameshift_events.tsv",
        "phase1/frameshift_screening/frameshift_hits.tsv",
    }
    schemas = {identity.relative_path: identity.schema for identity in identities}
    assert schemas["phase1/frameshift_screening/frameshift_hits.tsv"] == "frameshift-hits-v2"
    assert schemas["phase1/frameshift_screening/frameshift_events.tsv"] == "frameshift-events-v1"
