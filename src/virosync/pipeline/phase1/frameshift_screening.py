"""Frameshift-sensitive screening and rescue of VS marker domains."""

from __future__ import annotations

import csv
import hashlib
import logging
import re
import shutil
import subprocess
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass, fields
from pathlib import Path

from virosync.utils.atomic_write import atomic_write_context


logger = logging.getLogger(__name__)


ANNOTATION_CLASS = "frameshift_rescue_candidate"
COMMAND_TIMEOUT_SECONDS = 3600
VS_NAME_PATTERN = re.compile(r"VS[0-9]{6}")
RESCUED_PROTEIN_ID_PATTERN = re.compile(r"_VSR[0-9a-f]{16}$")
MIN_MODEL_COVERAGE = 0.5
MIN_DIAMOND_QUERY_COVERAGE = 50.0


@dataclass(frozen=True)
class FrameshiftHit:
    """One event-bearing BATH hit in normalized nucleotide coordinates."""

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
    pid: float
    shifts: int
    stops: int
    description: str


@dataclass(frozen=True)
class BathTranslation:
    """One model-conditioned amino-acid alignment from BATH's text report."""

    query_name: str
    hit_id: str
    target_name: str
    ali_from: int
    ali_to: int
    sequence: str
    literal_stops: int


def rescued_protein_id(hit: FrameshiftHit) -> str:
    """Return a stable ID whose Prodigal-style suffix preserves the scaffold."""

    token = hashlib.sha256(
        f"{hit.query_name}\0{hit.target_name}\0{hit.ali_start}\0"
        f"{hit.ali_end}\0{hit.strand}".encode("utf-8")
    ).hexdigest()[:16]
    return f"{hit.target_name}_VSR{token}"


def is_rescued_protein_id(value: str) -> bool:
    """Return whether a protein ID has ViroSync's generated rescue suffix."""

    return bool(
        isinstance(value, str)
        and RESCUED_PROTEIN_ID_PATTERN.search(value.split("|aa", 1)[0])
    )


def filter_vs_profiles(hmm_path: Path, output_path: Path) -> int:
    """Stream complete HMMER records and write exact ``VS######`` profiles."""

    hmm_path = Path(hmm_path)
    output_path = Path(output_path)
    selected_names: set[str] = set()
    record: list[str] = []

    with hmm_path.open(encoding="utf-8") as source, output_path.open(
        "x", encoding="utf-8"
    ) as output:
        for line in source:
            record.append(line)
            if line.strip() != "//":
                continue

            names = [
                parts[1]
                for record_line in record
                if len(parts := record_line.split()) >= 2 and parts[0] == "NAME"
            ]
            if len(names) == 1 and VS_NAME_PATTERN.fullmatch(names[0]):
                name = names[0]
                if name in selected_names:
                    raise ValueError(f"Duplicate VS profile NAME in {hmm_path}: {name}")
                selected_names.add(name)
                output.writelines(record)
            record = []

        if any(line.strip() for line in record):
            raise ValueError(f"Incomplete HMMER record in {hmm_path}")
        if not selected_names:
            raise ValueError(f"No VS###### profiles found in {hmm_path}")

    return len(selected_names)


def parse_bathsearch_tblout(path: Path) -> list[FrameshiftHit]:
    """Parse event-bearing BATH hits and normalize alignment coordinates."""

    hits: list[FrameshiftHit] = []
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split(maxsplit=17)
            if len(columns) == 17:
                columns.append("")
            if len(columns) != 18:
                raise ValueError(
                    f"Invalid bathsearch tblout row at {path}:{line_number}: "
                    f"expected 18 columns, found {len(columns)}"
                )
            (
                hit_id,
                target_name,
                target_accession,
                query_name,
                query_accession,
                hmm_len,
                hmm_from,
                hmm_to,
                seq_len,
                ali_from,
                ali_to,
                evalue,
                score,
                bias,
                pid,
                shifts,
                stops,
                description,
            ) = columns
            try:
                ali_from_int = int(ali_from)
                ali_to_int = int(ali_to)
                shifts_int = int(shifts)
                stops_int = int(stops)
                numeric = {
                    "hmm_len": int(hmm_len),
                    "hmm_from": int(hmm_from),
                    "hmm_to": int(hmm_to),
                    "seq_len": int(seq_len),
                    "evalue": float(evalue),
                    "score": float(score),
                    "bias": float(bias),
                    "pid": float(pid),
                }
            except ValueError as error:
                raise ValueError(
                    f"Invalid numeric value in bathsearch tblout at {path}:{line_number}"
                ) from error
            if shifts_int <= 0 and stops_int <= 0:
                continue
            hits.append(
                FrameshiftHit(
                    annotation_class=ANNOTATION_CLASS,
                    hit_id=hit_id,
                    target_name=target_name,
                    target_accession=target_accession,
                    query_name=query_name,
                    query_accession=query_accession,
                    ali_start=min(ali_from_int, ali_to_int) - 1,
                    ali_end=max(ali_from_int, ali_to_int),
                    strand="+" if ali_from_int <= ali_to_int else "-",
                    shifts=shifts_int,
                    stops=stops_int,
                    description=description,
                    **numeric,
                )
            )
    return hits


def write_frameshift_hits(hits: list[FrameshiftHit], output_path: Path) -> None:
    """Write the stable normalized diagnostic TSV."""

    output_path = Path(output_path)
    with atomic_write_context(output_path, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        field_names = [field.name for field in fields(FrameshiftHit)]
        writer.writerow(field_names)
        for hit in hits:
            writer.writerow([getattr(hit, name) for name in field_names])


def parse_bathsearch_translations(path: Path) -> dict[tuple[str, str], BathTranslation]:
    """Parse BATH's model-conditioned target amino-acid alignment rows.

    BATH has no machine-readable protein export. Its text report prints the
    amino-acid row immediately before each nucleotide target row. Hit IDs are
    ordinals within each query and match ``--tblout`` ordering.
    """

    translations: dict[tuple[str, str], BathTranslation] = {}
    query_name = ""
    target_name = ""
    hit_ordinal = 0
    amino_acids: list[str] = []
    ali_from: int | None = None
    ali_to: int | None = None
    previous_line = ""

    def finish_hit() -> None:
        nonlocal target_name, amino_acids, ali_from, ali_to
        if not target_name:
            return
        if not amino_acids or ali_from is None or ali_to is None:
            raise ValueError(
                f"Missing translated alignment for {query_name} hit {hit_ordinal} "
                f"({target_name}) in {path}"
            )
        raw_sequence = "".join(amino_acids).replace("-", "").upper()
        if not raw_sequence:
            raise ValueError(
                f"Empty translated alignment for {query_name} hit {hit_ordinal} in {path}"
            )
        invalid = sorted(set(raw_sequence) - set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*"))
        if invalid:
            raise ValueError(
                f"Invalid translated residues for {query_name} hit {hit_ordinal} "
                f"in {path}: {''.join(invalid)}"
            )
        key = (query_name, str(hit_ordinal))
        if key in translations:
            raise ValueError(f"Duplicate BATH translation key {key} in {path}")
        translations[key] = BathTranslation(
            query_name=query_name,
            hit_id=str(hit_ordinal),
            target_name=target_name,
            ali_from=ali_from,
            ali_to=ali_to,
            sequence=raw_sequence.replace("*", "X"),
            literal_stops=raw_sequence.count("*"),
        )
        target_name = ""
        amino_acids = []
        ali_from = None
        ali_to = None

    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("Query:"):
                finish_hit()
                parts = line.split()
                if len(parts) < 2:
                    raise ValueError(f"Invalid BATH query line at {path}:{line_number}")
                query_name = parts[1]
                hit_ordinal = 0
                previous_line = line
                continue
            if line.startswith(">> "):
                finish_hit()
                if not query_name:
                    raise ValueError(f"BATH hit precedes query at {path}:{line_number}")
                parts = line[3:].split()
                if not parts:
                    raise ValueError(f"Invalid BATH hit line at {path}:{line_number}")
                target_name = parts[0]
                hit_ordinal += 1
                previous_line = line
                continue

            if target_name:
                tokens = line.split()
                is_target_row = (
                    len(tokens) >= 4
                    and tokens[0] == target_name
                    and tokens[1].lstrip("-").isdigit()
                    and tokens[-1].lstrip("-").isdigit()
                    and all(
                        re.fullmatch(r"[ACGTUNacgtun-]+", token)
                        for token in tokens[2:-1]
                    )
                )
                if is_target_row:
                    aa_tokens = previous_line.split()
                    if not aa_tokens or not all(
                        re.fullmatch(r"[A-Za-z*-]", token) for token in aa_tokens
                    ):
                        raise ValueError(
                            f"Invalid BATH amino-acid row before {path}:{line_number}"
                        )
                    amino_acids.extend(aa_tokens)
                    if ali_from is None:
                        ali_from = int(tokens[1])
                    ali_to = int(tokens[-1])
            previous_line = line
    finish_hit()
    if not translations:
        raise ValueError(f"No translated BATH alignments found in {path}")
    return translations


def write_frameshift_candidate_faa(
    hits: Iterable[FrameshiftHit],
    translations: dict[tuple[str, str], BathTranslation],
    output_path: Path,
) -> dict[str, FrameshiftHit]:
    """Write event-bearing BATH aligned domains as pseudo-protein candidates."""

    by_protein_id: dict[str, FrameshiftHit] = {}
    with atomic_write_context(output_path, "w") as handle:
        for hit in hits:
            key = (hit.query_name, hit.hit_id)
            translation = translations.get(key)
            if translation is None:
                raise ValueError(f"No BATH translation found for event hit {key}")
            raw_from = hit.ali_start + 1 if hit.strand == "+" else hit.ali_end
            raw_to = hit.ali_end if hit.strand == "+" else hit.ali_start + 1
            if (
                translation.target_name != hit.target_name
                or translation.ali_from != raw_from
                or translation.ali_to != raw_to
            ):
                raise ValueError(
                    f"BATH text/tblout mismatch for {key}: "
                    f"text={translation.target_name}:{translation.ali_from}-{translation.ali_to}, "
                    f"table={hit.target_name}:{raw_from}-{raw_to}"
                )
            protein_id = rescued_protein_id(hit)
            if protein_id in by_protein_id:
                raise ValueError(f"Duplicate rescued protein ID: {protein_id}")
            suffix = protein_id.rsplit("_", 1)[1]
            strand = "1" if hit.strand == "+" else "-1"
            handle.write(
                f">{protein_id} # {hit.ali_start + 1} # {hit.ali_end} # {strand} # "
                f"ID=0_{suffix};annotation=frameshift_rescued_domain;"
                f"model={hit.query_name};shifts={hit.shifts};stops={hit.stops};"
                f"literal_stops={translation.literal_stops}\n"
            )
            sequence = translation.sequence
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
                raise ValueError(
                    f"Invalid DIAMOND numeric value at {path}:{line_number}"
                ) from error
            if prefix in validated_prefixes and pident_value >= min_pident:
                coverages[query] = max(coverages.get(query, 0.0), qcov_value)
    return coverages


def select_confirmed_frameshift_markers(
    markers: Iterable[object],
    hits_by_protein_id: dict[str, FrameshiftHit],
    diamond_output: Path,
    *,
    validated_prefixes: set[str],
    min_pident: float,
    min_model_coverage: float = MIN_MODEL_COVERAGE,
    min_diamond_query_coverage: float = MIN_DIAMOND_QUERY_COVERAGE,
) -> list[object]:
    """Retain DIAMOND-confirmed, sufficiently covered, non-overlapping loci."""

    qcov_by_query = diamond_query_coverages(
        diamond_output,
        min_pident=min_pident,
        validated_prefixes=validated_prefixes,
    )
    eligible: list[object] = []
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
    confirmed: list[object] = []
    cluster: list[object] = []
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
    return sorted(confirmed, key=lambda item: (item.scaffold, item.start, item.end))


def write_confirmed_frameshift_faa(
    candidate_faa: Path,
    confirmed_markers: Iterable[object],
    output_path: Path,
) -> int:
    """Write candidate records retained as confirmed rescue seed markers."""

    from Bio import SeqIO

    confirmed_ids = {
        str(marker.query_porf).split("|aa", 1)[0] for marker in confirmed_markers
    }
    records = [
        record
        for record in SeqIO.parse(candidate_faa, "fasta")
        if record.id in confirmed_ids
    ]
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


def _required_tool(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        raise RuntimeError(f"Frameshift screening requires '{name}' on PATH")
    return executable


def _run_command(args: list[str], label: str) -> subprocess.CompletedProcess[str]:
    try:
        result = subprocess.run(
            args,
            capture_output=True,
            text=True,
            check=False,
            timeout=COMMAND_TIMEOUT_SECONDS,
        )
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(
            f"{label} timed out after {COMMAND_TIMEOUT_SECONDS} seconds"
        ) from error
    except OSError as error:
        raise RuntimeError(f"Failed to start {label}: {error}") from error
    if result.returncode != 0:
        detail = (result.stderr or result.stdout).strip()
        suffix = f": {detail}" if detail else ""
        raise RuntimeError(f"{label} failed with exit code {result.returncode}{suffix}")
    return result


def run_frameshift_screening(
    masked_fasta: Path,
    hmm_database: Path,
    output_dir: Path,
    threads: int,
) -> list[FrameshiftHit]:
    """Run the VS-profile BATH screen and write event-bearing candidate domains."""

    bathconvert = _required_tool("bathconvert")
    bathsearch = _required_tool("bathsearch")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tblout_path = output_dir / "bathsearch.tblout"
    fstblout_path = output_dir / "bathsearch.fstblout"
    main_output_path = output_dir / "bathsearch.txt"
    hits_path = output_dir / "frameshift_hits.tsv"
    candidates_path = output_dir / "frameshift_candidates.faa"
    with tempfile.TemporaryDirectory(prefix="virosync_frameshift_") as temp_dir:
        filtered_hmm = Path(temp_dir) / "vs_profiles.hmm"
        models_path = Path(temp_dir) / "vs_profiles.bhmm"
        selected_profile_count = filter_vs_profiles(hmm_database, filtered_hmm)
        logger.info("Selected %d VS marker profiles for BATH", selected_profile_count)
        _run_command(
            [bathconvert, str(models_path), str(filtered_hmm)],
            "bathconvert",
        )
        _run_command(
            [
                bathsearch,
                "--fs",
                "--cpu",
                str(threads),
                "-E",
                "1e-5",
                "--incE",
                "1e-5",
                "--tblout",
                str(tblout_path),
                "--fstblout",
                str(fstblout_path),
                "-o",
                str(main_output_path),
                str(models_path),
                str(masked_fasta),
            ],
            "bathsearch",
        )
    hits = parse_bathsearch_tblout(tblout_path)
    write_frameshift_hits(hits, hits_path)
    if hits:
        translations = parse_bathsearch_translations(main_output_path)
        write_frameshift_candidate_faa(hits, translations, candidates_path)
    else:
        candidates_path.write_text("", encoding="utf-8")
    return hits
