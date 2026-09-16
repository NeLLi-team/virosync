"""Find and validate frameshift-disrupted VS marker domains in genomic DNA."""

from __future__ import annotations

import csv
import hashlib
import logging
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, fields, replace
from pathlib import Path
from typing import Protocol, TypeVar

import pyhmmer
from Bio import SeqIO
from Bio.Seq import Seq

from virosync.pipeline.phase1.native_frameshift import (
    CodonAlignment,
    FrameshiftEvent,
    align_frameshift_profile,
)
from virosync.utils.atomic_write import atomic_write_context

logger = logging.getLogger(__name__)


ANNOTATION_CLASS = "frameshift_rescue_candidate"
VS_NAME_PATTERN = re.compile(r"VS[0-9]{6}")
RESCUED_PROTEIN_ID_PATTERN = re.compile(r"_VSR[0-9a-f]{16}$")
MIN_MODEL_COVERAGE = 0.5
MIN_DIAMOND_QUERY_COVERAGE = 50.0

SEED_CHUNK_NT = 3_000
SEED_OVERLAP_NT = 1_500
SEED_MAX_EVALUE = 10.0
SEED_MIN_DOMAIN_BITS = 8.0
SEED_FILTER_F1 = 0.1
SEED_FILTER_F2 = 0.1
SEED_FILTER_F3 = 0.1
ALIGNMENT_WINDOW_MAX_NT = 12_000
ALIGNMENT_EXTRA_FLANK_NT = 90
EVENT_FLANK_AA = 8
MIN_FRAMESHIFT_GAIN_BITS = 6.0
PEPTIDE_MAX_EVALUE = 1e-5
EVALUE_INTERPRETATION = "conditional peptide support; not calibrated native DNA significance"


@dataclass(frozen=True, slots=True)
class FrameshiftHit:
    """One validated event-bearing marker hit in genomic coordinates."""

    annotation_class: str
    hit_id: str
    target_name: str
    target_accession: str
    query_name: str
    query_accession: str
    hmm_len: int
    hmm_from: int
    hmm_to: int
    seq_len: int
    ali_start: int
    ali_end: int
    strand: str
    evalue: float
    score: float
    bias: float
    pid: float | None
    shifts: int
    stops: int
    description: str
    backend: str = "native_codon"
    native_score: float = 0.0
    no_frameshift_score: float = 0.0
    evalue_interpretation: str = EVALUE_INTERPRETATION


@dataclass(frozen=True, slots=True)
class FrameshiftEventRecord:
    """One retained disruption mapped to genomic coordinates."""

    candidate_id: str
    target_name: str
    query_name: str
    strand: str
    event_kind: str
    genomic_start: int
    genomic_end: int
    size: int
    model_position: int
    protein_position: int


@dataclass(frozen=True, slots=True)
class _TranslatedChunk:
    """Coordinate map for one translated frame of a bounded DNA chunk."""

    name: str
    contig_name: str
    chunk_start: int
    chunk_end: int
    strand: str
    frame: int
    sequence: str


@dataclass(frozen=True, slots=True)
class _Seed:
    """One permissive protein-profile seed in genomic coordinates."""

    query_name: str
    target_name: str
    strand: str
    start: int
    end: int


@dataclass(frozen=True, slots=True)
class _CandidateWindow:
    """One bounded native-alignment window and the seeds it contains."""

    query_name: str
    target_name: str
    strand: str
    start: int
    end: int
    seed_spans: tuple[tuple[int, int], ...]


@dataclass(frozen=True, slots=True)
class _PeptideDomain:
    """Best strictly supported domain on a native reconstructed peptide."""

    sequence_score: float
    bias: float
    evalue: float
    hmm_from: int
    hmm_to: int
    target_start: int
    target_end: int


@dataclass(frozen=True, slots=True)
class _Candidate:
    """A validated output hit with its peptide and retained events."""

    hit: FrameshiftHit
    sequence: str
    events: tuple[FrameshiftEventRecord, ...]


class _Marker(Protocol):
    """Fields used to confirm a validated rescue marker."""

    query_porf: str
    scaffold: str
    start: int
    end: int
    strand: str
    hmm_score: float
    validation_status: str


_MarkerT = TypeVar("_MarkerT", bound=_Marker)


def rescued_protein_id(hit: FrameshiftHit) -> str:
    """Return a stable ID whose Prodigal-style suffix preserves the scaffold."""
    token = hashlib.sha256(
        (f"{hit.query_name}\0{hit.target_name}\0{hit.ali_start}\0{hit.ali_end}\0{hit.strand}").encode()
    ).hexdigest()[:16]
    return f"{hit.target_name}_VSR{token}"


def is_rescued_protein_id(value: str) -> bool:
    """Return whether a protein ID has ViroSync's generated rescue suffix."""
    return bool(isinstance(value, str) and RESCUED_PROTEIN_ID_PATTERN.search(value.split("|aa", 1)[0]))


def write_frameshift_hits(hits: Iterable[FrameshiftHit], output_path: Path) -> None:
    """Write the normalized native-screening diagnostic table."""
    field_names = [field.name for field in fields(FrameshiftHit)]
    with atomic_write_context(output_path, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(field_names)
        for hit in hits:
            writer.writerow([getattr(hit, name) for name in field_names])


def write_frameshift_events(
    events: Iterable[FrameshiftEventRecord],
    output_path: Path,
) -> None:
    """Write retained native events with zero-based half-open coordinates."""
    field_names = [field.name for field in fields(FrameshiftEventRecord)]
    with atomic_write_context(output_path, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(field_names)
        for event in events:
            writer.writerow([getattr(event, name) for name in field_names])


def write_frameshift_candidate_faa(
    hits: Iterable[FrameshiftHit],
    sequences_by_protein_id: Mapping[str, str],
    output_path: Path,
) -> dict[str, FrameshiftHit]:
    """Write validated reconstructed domains as pseudo-protein candidates."""
    by_protein_id: dict[str, FrameshiftHit] = {}
    with atomic_write_context(output_path, "w") as handle:
        for hit in hits:
            protein_id = rescued_protein_id(hit)
            sequence = sequences_by_protein_id.get(protein_id)
            if sequence is None:
                raise ValueError(f"No native translation found for event hit {hit.hit_id}")
            if protein_id in by_protein_id:
                raise ValueError(f"Duplicate rescued protein ID: {protein_id}")
            suffix = protein_id.rsplit("_", 1)[1]
            strand = "1" if hit.strand == "+" else "-1"
            handle.write(
                f">{protein_id} # {hit.ali_start + 1} # {hit.ali_end} # {strand} # "
                f"ID=0_{suffix};annotation=frameshift_rescued_domain;"
                f"model={hit.query_name};shifts={hit.shifts};stops={hit.stops};"
                f"literal_stops={hit.stops}\n"
            )
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset : offset + 60] + "\n")
            by_protein_id[protein_id] = hit
    return by_protein_id


def diamond_query_coverages(
    path: Path,
    *,
    min_pident: float,
    validated_prefixes: set[str],
) -> dict[str, float]:
    """Return maximum qualifying viral DIAMOND query coverage per query."""
    coverages: dict[str, float] = {}
    if not Path(path).is_file():
        return coverages
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 6:
                continue
            query, target, _evalue, _bits, pident, qcov = columns[:6]
            prefix = target.split("__", 1)[0] + "__" if "__" in target else ""
            try:
                pident_value = float(pident)
                qcov_value = float(qcov)
            except ValueError as error:
                raise ValueError(f"Invalid DIAMOND numeric value at {path}:{line_number}") from error
            if prefix in validated_prefixes and pident_value >= min_pident:
                coverages[query] = max(coverages.get(query, 0.0), qcov_value)
    return coverages


def select_confirmed_frameshift_markers(
    markers: Iterable[_MarkerT],
    hits_by_protein_id: dict[str, FrameshiftHit],
    diamond_output: Path,
    *,
    validated_prefixes: set[str],
    min_pident: float,
    min_model_coverage: float = MIN_MODEL_COVERAGE,
    min_diamond_query_coverage: float = MIN_DIAMOND_QUERY_COVERAGE,
) -> list[_MarkerT]:
    """Retain DIAMOND-confirmed, sufficiently covered, non-overlapping loci."""
    qcov_by_query = diamond_query_coverages(
        diamond_output,
        min_pident=min_pident,
        validated_prefixes=validated_prefixes,
    )
    eligible: list[_MarkerT] = []
    for marker in markers:
        base_id = str(marker.query_porf).split("|aa", 1)[0]
        hit = hits_by_protein_id.get(base_id)
        if hit is None or marker.validation_status != "validated":
            continue
        model_coverage = (hit.hmm_to - hit.hmm_from + 1) / hit.hmm_len
        if model_coverage < min_model_coverage:
            continue
        if qcov_by_query.get(base_id, 0.0) < min_diamond_query_coverage:
            continue
        eligible.append(marker)

    sorted_markers = sorted(
        eligible,
        key=lambda item: (
            item.scaffold,
            item.strand,
            item.start,
            item.end,
            -item.hmm_score,
        ),
    )
    confirmed: list[_MarkerT] = []
    cluster: list[_MarkerT] = []
    cluster_end = -1

    def finish_cluster() -> None:
        if cluster:
            confirmed.append(max(cluster, key=lambda item: item.hmm_score))

    for marker in sorted_markers:
        same_locus_group = (
            cluster
            and cluster[0].scaffold == marker.scaffold
            and cluster[0].strand == marker.strand
            and marker.start < cluster_end
        )
        if not same_locus_group:
            finish_cluster()
            cluster = [marker]
            cluster_end = marker.end
            continue
        cluster.append(marker)
        cluster_end = max(cluster_end, marker.end)
    finish_cluster()
    return sorted(
        confirmed,
        key=lambda item: (item.scaffold, item.start, item.end),
    )


def write_confirmed_frameshift_faa(
    candidate_faa: Path,
    confirmed_markers: Iterable[_Marker],
    output_path: Path,
) -> int:
    """Write candidate records retained as confirmed rescue seed markers."""
    confirmed_ids = {str(marker.query_porf).split("|aa", 1)[0] for marker in confirmed_markers}
    records = [record for record in SeqIO.parse(candidate_faa, "fasta") if record.id in confirmed_ids]
    with atomic_write_context(output_path, "w") as handle:
        SeqIO.write(records, handle, "fasta")
    return len(records)


def write_confirmed_frameshift_markers(
    markers: Iterable[object],
    output_path: Path,
) -> None:
    """Write the confirmed rescue markers consumed by EVE reporting."""
    columns = (
        "query_porf",
        "scaffold",
        "start",
        "end",
        "strand",
        "hmm_target",
        "hmm_score",
        "validation_status",
    )
    with atomic_write_context(output_path, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for marker in markers:
            writer.writerow([getattr(marker, column) for column in columns])


def run_frameshift_screening(
    masked_fasta: Path,
    hmm_database: Path,
    output_dir: Path,
    threads: int,
) -> list[FrameshiftHit]:
    """Run native codon-aware screening and write event-bearing candidates."""
    if threads < 1:
        raise ValueError(f"threads must be positive, got {threads}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    hmms = _load_vs_profiles(hmm_database)
    contigs = _load_contigs(masked_fasta)
    translated_chunks = _translate_chunks(contigs)
    seeds = _search_seeds(hmms, translated_chunks, threads)
    windows = _merge_seed_windows(seeds, contigs, hmms)
    candidates = _screen_windows(windows, contigs, hmms)

    hits = [candidate.hit for candidate in candidates]
    sequences = {rescued_protein_id(candidate.hit): candidate.sequence for candidate in candidates}
    events = [event for candidate in candidates for event in candidate.events]
    write_frameshift_hits(hits, output_dir / "frameshift_hits.tsv")
    write_frameshift_events(events, output_dir / "frameshift_events.tsv")
    write_frameshift_candidate_faa(
        hits,
        sequences,
        output_dir / "frameshift_candidates.faa",
    )
    logger.info(
        "Native frameshift screening retained %d candidates and %d events",
        len(hits),
        len(events),
    )
    return hits


def _load_vs_profiles(hmm_path: Path) -> list[pyhmmer.plan7.HMM]:
    """Load exact VS marker profiles and reject ambiguous profile databases."""
    selected: list[pyhmmer.plan7.HMM] = []
    selected_names: set[str] = set()
    try:
        with pyhmmer.plan7.HMMFile(hmm_path) as handle:
            for hmm in handle:
                name = _decode(hmm.name)
                if not VS_NAME_PATTERN.fullmatch(name):
                    continue
                if name in selected_names:
                    raise ValueError(f"Duplicate VS profile NAME in {hmm_path}: {name}")
                selected.append(hmm)
                selected_names.add(name)
    except (EOFError, OSError, ValueError) as error:
        if isinstance(error, ValueError) and str(error).startswith("Duplicate VS profile NAME"):
            raise
        raise ValueError(f"Malformed HMM database {hmm_path}: {error}") from error
    if not selected:
        raise ValueError(f"No VS###### profiles found in {hmm_path}")
    return selected


def _load_contigs(fasta_path: Path) -> dict[str, str]:
    """Load unique named contigs as uppercase DNA."""
    contigs: dict[str, str] = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        if record.id in contigs:
            raise ValueError(f"Duplicate FASTA record ID in {fasta_path}: {record.id}")
        contigs[record.id] = str(record.seq).upper().replace("U", "T")
    return contigs


def _translate_chunks(
    contigs: Mapping[str, str],
) -> list[_TranslatedChunk]:
    """Translate bounded overlapping DNA chunks in all six reading frames."""
    chunks: list[_TranslatedChunk] = []
    chunk_number = 0
    for contig_name, dna in contigs.items():
        for chunk_start, chunk_end in _chunk_spans(
            len(dna),
            SEED_CHUNK_NT,
            SEED_OVERLAP_NT,
        ):
            forward = dna[chunk_start:chunk_end]
            for strand, oriented in (("+", forward), ("-", _reverse_complement(forward))):
                for frame in range(3):
                    protein = _translate_frame(oriented, frame)
                    if protein == "":
                        continue
                    chunks.append(
                        _TranslatedChunk(
                            name=f"frame_{chunk_number}",
                            contig_name=contig_name,
                            chunk_start=chunk_start,
                            chunk_end=chunk_end,
                            strand=strand,
                            frame=frame,
                            sequence=protein,
                        )
                    )
                    chunk_number += 1
    return chunks


def _chunk_spans(
    sequence_length: int,
    chunk_size: int,
    overlap: int,
) -> Iterable[tuple[int, int]]:
    """Yield bounded spans that cover a sequence without gaps."""
    start = 0
    while start < sequence_length:
        end = min(start + chunk_size, sequence_length)
        yield start, end
        if end == sequence_length:
            return
        start = end - overlap


def _translate_frame(dna: str, frame: int) -> str:
    """Translate one frame, representing stops as unknown amino acids."""
    usable_length = ((len(dna) - frame) // 3) * 3
    if usable_length <= 0:
        return ""
    sequence = str(Seq(dna[frame : frame + usable_length]).translate())
    return sequence.replace("*", "X")


def _search_seeds(
    hmms: list[pyhmmer.plan7.HMM],
    chunks: list[_TranslatedChunk],
    threads: int,
) -> list[_Seed]:
    """Find permissive profile seeds in translated chunks with PyHMMER."""
    if not chunks:
        return []
    alphabet = pyhmmer.easel.Alphabet.amino()
    metadata = {chunk.name: chunk for chunk in chunks}
    sequences = [
        pyhmmer.easel.TextSequence(
            name=chunk.name.encode(),
            sequence=chunk.sequence,
        ).digitize(alphabet)
        for chunk in chunks
    ]
    results = pyhmmer.hmmsearch(
        hmms,
        sequences,
        cpus=threads,
        parallel="targets",
        E=SEED_MAX_EVALUE,
        domE=SEED_MAX_EVALUE,
        incE=SEED_MAX_EVALUE,
        incdomE=SEED_MAX_EVALUE,
        F1=SEED_FILTER_F1,
        F2=SEED_FILTER_F2,
        F3=SEED_FILTER_F3,
    )
    seeds: set[_Seed] = set()
    for hmm, top_hits in zip(hmms, results, strict=True):
        query_name = _decode(hmm.name)
        for hit in top_hits:
            chunk = metadata[_decode(hit.name)]
            for domain in hit.domains:
                if domain.i_evalue > SEED_MAX_EVALUE or domain.score < SEED_MIN_DOMAIN_BITS:
                    continue
                alignment = domain.alignment
                target_start = alignment.target_from - 1
                target_end = alignment.target_to
                start, end = _translated_to_genomic(
                    chunk,
                    target_start,
                    target_end,
                )
                seeds.add(
                    _Seed(
                        query_name=query_name,
                        target_name=chunk.contig_name,
                        strand=chunk.strand,
                        start=start,
                        end=end,
                    )
                )
    return sorted(
        seeds,
        key=lambda seed: (
            seed.query_name,
            seed.target_name,
            seed.strand,
            seed.start,
            seed.end,
        ),
    )


def _translated_to_genomic(
    chunk: _TranslatedChunk,
    amino_start: int,
    amino_end: int,
) -> tuple[int, int]:
    """Map a translated target interval to genomic coordinates."""
    oriented_start = chunk.frame + amino_start * 3
    oriented_end = chunk.frame + amino_end * 3
    if chunk.strand == "+":
        return (
            chunk.chunk_start + oriented_start,
            chunk.chunk_start + oriented_end,
        )
    return (
        chunk.chunk_end - oriented_end,
        chunk.chunk_end - oriented_start,
    )


def _merge_seed_windows(
    seeds: Iterable[_Seed],
    contigs: Mapping[str, str],
    hmms: Iterable[pyhmmer.plan7.HMM],
) -> list[_CandidateWindow]:
    """Merge overlapping seed windows while keeping native DP bounded."""
    hmms_by_name = {_decode(hmm.name): hmm for hmm in hmms}
    proposed: list[_CandidateWindow] = []
    for seed in seeds:
        hmm = hmms_by_name[seed.query_name]
        flank = hmm.M * 3 + ALIGNMENT_EXTRA_FLANK_NT
        contig_length = len(contigs[seed.target_name])
        proposed.append(
            _CandidateWindow(
                query_name=seed.query_name,
                target_name=seed.target_name,
                strand=seed.strand,
                start=max(0, seed.start - flank),
                end=min(contig_length, seed.end + flank),
                seed_spans=((seed.start, seed.end),),
            )
        )

    proposed.sort(
        key=lambda window: (
            window.query_name,
            window.target_name,
            window.strand,
            window.start,
            window.end,
        )
    )
    merged: list[_CandidateWindow] = []
    for window in proposed:
        if not merged or not _can_merge_windows(
            merged[-1],
            window,
            hmms_by_name[window.query_name].M,
        ):
            merged.append(window)
            continue
        previous = merged[-1]
        merged[-1] = replace(
            previous,
            end=max(previous.end, window.end),
            seed_spans=tuple(sorted(set(previous.seed_spans + window.seed_spans))),
        )
    return merged


def _can_merge_windows(
    left: _CandidateWindow,
    right: _CandidateWindow,
    model_length: int,
) -> bool:
    """Return whether two windows share a group and fit the DP size bound."""
    same_group = (
        left.query_name == right.query_name and left.target_name == right.target_name and left.strand == right.strand
    )
    max_length = max(
        ALIGNMENT_WINDOW_MAX_NT,
        model_length * 9 + 2 * ALIGNMENT_EXTRA_FLANK_NT,
    )
    return same_group and right.start <= left.end and max(left.end, right.end) - left.start <= max_length


def _screen_windows(
    windows: Iterable[_CandidateWindow],
    contigs: Mapping[str, str],
    hmms: Iterable[pyhmmer.plan7.HMM],
) -> list[_Candidate]:
    """Run native alignment, peptide validation, event cropping, and deduplication."""
    hmms_by_name = {_decode(hmm.name): hmm for hmm in hmms}
    candidates: list[_Candidate] = []
    for window in windows:
        hmm = hmms_by_name[window.query_name]
        genomic_dna = contigs[window.target_name][window.start : window.end]
        oriented_dna = genomic_dna if window.strand == "+" else _reverse_complement(genomic_dna)
        seed_spans = tuple(_orient_seed_span(window, span) for span in window.seed_spans)
        for alignment in _enumerate_alignments(hmm, oriented_dna, seed_spans):
            candidate = _validate_candidate(
                window,
                len(contigs[window.target_name]),
                hmm,
                oriented_dna,
                alignment,
            )
            if candidate is not None:
                candidates.append(candidate)

    best_by_locus: dict[tuple[str, str, str, int, int], _Candidate] = {}
    for candidate in candidates:
        hit = candidate.hit
        key = (
            hit.query_name,
            hit.target_name,
            hit.strand,
            hit.ali_start,
            hit.ali_end,
        )
        previous = best_by_locus.get(key)
        if previous is None or hit.native_score > previous.hit.native_score:
            best_by_locus[key] = candidate
    return sorted(
        best_by_locus.values(),
        key=lambda candidate: (
            candidate.hit.target_name,
            candidate.hit.ali_start,
            candidate.hit.ali_end,
            candidate.hit.query_name,
            candidate.hit.strand,
        ),
    )


def _orient_seed_span(
    window: _CandidateWindow,
    span: tuple[int, int],
) -> tuple[int, int]:
    """Map a genomic seed span into an oriented candidate window."""
    start, end = span
    if window.strand == "+":
        return start - window.start, end - window.start
    return window.end - end, window.end - start


def _enumerate_alignments(
    hmm: pyhmmer.plan7.HMM,
    dna: str,
    seed_spans: tuple[tuple[int, int], ...],
) -> Iterable[CodonAlignment]:
    """Yield non-overlapping local alignments until every seed is resolved."""
    pending = [(0, len(dna), seed_spans)]
    while pending:
        fragment_start, fragment_end, fragment_seeds = pending.pop()
        if not fragment_seeds or fragment_end <= fragment_start:
            continue
        alignment = align_frameshift_profile(
            hmm,
            dna[fragment_start:fragment_end],
        )
        if alignment is None or alignment.nt_end <= alignment.nt_start:
            continue
        absolute = _offset_alignment(alignment, fragment_start)
        yield absolute

        left_seeds: list[tuple[int, int]] = []
        right_seeds: list[tuple[int, int]] = []
        for seed_start, seed_end in fragment_seeds:
            seed_center = (seed_start + seed_end) // 2
            if seed_center < absolute.nt_start:
                left_seeds.append((seed_start, seed_end))
            elif seed_center >= absolute.nt_end:
                right_seeds.append((seed_start, seed_end))
        if left_seeds:
            pending.append((fragment_start, absolute.nt_start, tuple(left_seeds)))
        if right_seeds:
            pending.append((absolute.nt_end, fragment_end, tuple(right_seeds)))


def _offset_alignment(alignment: CodonAlignment, offset: int) -> CodonAlignment:
    """Shift an alignment from a fragment into its parent oriented window."""
    return replace(
        alignment,
        nt_start=alignment.nt_start + offset,
        nt_end=alignment.nt_end + offset,
        events=tuple(
            replace(
                event,
                nt_start=event.nt_start + offset,
                nt_end=event.nt_end + offset,
            )
            for event in alignment.events
        ),
        codon_spans=tuple((start + offset, end + offset) for start, end in alignment.codon_spans),
    )


def _validate_candidate(
    window: _CandidateWindow,
    contig_length: int,
    hmm: pyhmmer.plan7.HMM,
    oriented_dna: str,
    alignment: CodonAlignment,
) -> _Candidate | None:
    """Revalidate, crop, and retain only supported internal events."""
    if len(alignment.sequence) != len(alignment.codon_spans):
        raise ValueError(f"Native alignment sequence/span mismatch for {_decode(hmm.name)}")
    domain = _best_peptide_domain(hmm, alignment.sequence)
    if domain is None:
        return None

    selected_spans = alignment.codon_spans[domain.target_start : domain.target_end]
    if not selected_spans:
        return None
    crop_start = selected_spans[0][0]
    crop_end = selected_spans[-1][1]
    retained_events = tuple(
        event
        for event in alignment.events
        if _event_is_internal(
            event,
            domain,
            alignment.sequence,
            crop_start,
            crop_end,
        )
    )
    if not retained_events:
        return None
    scores = _cropped_native_scores(
        hmm,
        oriented_dna,
        alignment,
        crop_start,
        crop_end,
        retained_events,
    )
    if scores is None:
        return None
    native_score, baseline_score = scores
    has_frameshift = any(event.kind != "stop" for event in retained_events)
    if has_frameshift and native_score - baseline_score < MIN_FRAMESHIFT_GAIN_BITS:
        return None

    genomic_start, genomic_end = _oriented_to_genomic(
        window,
        crop_start,
        crop_end,
    )
    query_name = _decode(hmm.name)
    hit = FrameshiftHit(
        annotation_class=ANNOTATION_CLASS,
        hit_id=_native_hit_id(
            query_name,
            window.target_name,
            genomic_start,
            genomic_end,
            window.strand,
        ),
        target_name=window.target_name,
        target_accession="-",
        query_name=query_name,
        query_accession=_decode(hmm.accession) or "-",
        hmm_len=hmm.M,
        hmm_from=domain.hmm_from,
        hmm_to=domain.hmm_to,
        seq_len=contig_length,
        ali_start=genomic_start,
        ali_end=genomic_end,
        strand=window.strand,
        evalue=domain.evalue,
        score=domain.sequence_score,
        bias=domain.bias,
        pid=None,
        shifts=sum(event.kind != "stop" for event in retained_events),
        stops=sum(event.kind == "stop" for event in retained_events),
        description="native codon alignment with conditional peptide HMM support",
        native_score=native_score,
        no_frameshift_score=baseline_score,
    )
    candidate_id = rescued_protein_id(hit)
    events = tuple(
        _event_record(
            candidate_id,
            window,
            query_name,
            event,
            domain.target_start,
        )
        for event in retained_events
    )
    sequence = alignment.sequence[domain.target_start : domain.target_end]
    return _Candidate(hit=hit, sequence=sequence, events=events)


def _cropped_native_scores(
    hmm: pyhmmer.plan7.HMM,
    dna: str,
    alignment: CodonAlignment,
    crop_start: int,
    crop_end: int,
    retained_events: tuple[FrameshiftEvent, ...],
) -> tuple[float, float] | None:
    """Return native and ordinary-codon scores for the reported peptide span."""
    if (crop_start, crop_end) == (alignment.nt_start, alignment.nt_end):
        return alignment.score, alignment.baseline_score
    rescored = align_frameshift_profile(hmm, dna[crop_start:crop_end])
    if rescored is None:
        return None
    scored_start = crop_start + rescored.nt_start
    scored_end = crop_start + rescored.nt_end
    if any(event.nt_start < scored_start or event.nt_end > scored_end for event in retained_events):
        return None
    rescored_events = {
        (
            crop_start + event.nt_start,
            crop_start + event.nt_end,
            event.kind,
            event.size,
        )
        for event in rescored.events
    }
    if any((event.nt_start, event.nt_end, event.kind, event.size) not in rescored_events for event in retained_events):
        return None
    return rescored.score, rescored.baseline_score


def _best_peptide_domain(
    hmm: pyhmmer.plan7.HMM,
    sequence: str,
) -> _PeptideDomain | None:
    """Return the strongest sufficiently covered strict peptide HMM domain."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    digital = pyhmmer.easel.TextSequence(
        name=b"native_candidate",
        sequence=sequence,
    ).digitize(alphabet)
    top_hits = next(
        pyhmmer.hmmsearch(
            [hmm],
            [digital],
            cpus=1,
            E=PEPTIDE_MAX_EVALUE,
            domE=PEPTIDE_MAX_EVALUE,
            incE=PEPTIDE_MAX_EVALUE,
            incdomE=PEPTIDE_MAX_EVALUE,
        )
    )
    if not top_hits:
        return None
    hit = top_hits[0]
    if hit.evalue > PEPTIDE_MAX_EVALUE:
        return None
    supported = [
        domain
        for domain in hit.domains
        if domain.i_evalue <= PEPTIDE_MAX_EVALUE
        and (domain.alignment.hmm_to - domain.alignment.hmm_from + 1) / hmm.M >= MIN_MODEL_COVERAGE
    ]
    if not supported:
        return None
    domain = max(supported, key=lambda item: (item.score, -item.i_evalue))
    return _PeptideDomain(
        sequence_score=hit.score,
        bias=hit.bias,
        evalue=hit.evalue,
        hmm_from=domain.alignment.hmm_from,
        hmm_to=domain.alignment.hmm_to,
        target_start=domain.alignment.target_from - 1,
        target_end=domain.alignment.target_to,
    )


def _event_is_internal(
    event: FrameshiftEvent,
    domain: _PeptideDomain,
    sequence: str,
    crop_start: int,
    crop_end: int,
) -> bool:
    """Return whether an event has peptide support and context on both sides."""
    protein_position = event.protein_position
    left_flank = sequence[domain.target_start : protein_position]
    right_flank = sequence[protein_position + 1 : domain.target_end]
    return (
        sum(amino_acid != "X" for amino_acid in left_flank) >= EVENT_FLANK_AA
        and sum(amino_acid != "X" for amino_acid in right_flank) >= EVENT_FLANK_AA
        and crop_start <= event.nt_start
        and event.nt_end <= crop_end
    )


def _event_record(
    candidate_id: str,
    window: _CandidateWindow,
    query_name: str,
    event: FrameshiftEvent,
    peptide_start: int,
) -> FrameshiftEventRecord:
    """Map one retained native event into its output record."""
    genomic_start, genomic_end = _oriented_to_genomic(
        window,
        event.nt_start,
        event.nt_end,
    )
    return FrameshiftEventRecord(
        candidate_id=candidate_id,
        target_name=window.target_name,
        query_name=query_name,
        strand=window.strand,
        event_kind=event.kind,
        genomic_start=genomic_start,
        genomic_end=genomic_end,
        size=event.size,
        model_position=event.model_position,
        protein_position=event.protein_position - peptide_start,
    )


def _oriented_to_genomic(
    window: _CandidateWindow,
    oriented_start: int,
    oriented_end: int,
) -> tuple[int, int]:
    """Map an oriented half-open window span to genomic coordinates."""
    if window.strand == "+":
        return window.start + oriented_start, window.start + oriented_end
    return window.end - oriented_end, window.end - oriented_start


def _native_hit_id(
    query_name: str,
    target_name: str,
    start: int,
    end: int,
    strand: str,
) -> str:
    """Return a stable compact identifier for one native hit."""
    token = hashlib.sha256(f"{query_name}\0{target_name}\0{start}\0{end}\0{strand}".encode()).hexdigest()[:16]
    return f"NFS{token}"


def _reverse_complement(dna: str) -> str:
    """Return uppercase reverse-complemented DNA."""
    return str(Seq(dna).reverse_complement()).upper()


def _decode(value: bytes | str | None) -> str:
    """Decode an optional HMMER text field."""
    if value is None:
        return ""
    return value.decode() if isinstance(value, bytes) else value
