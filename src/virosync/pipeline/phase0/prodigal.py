"""
Genome-wide gene calling with Prodigal-GV.

Supports parallel gene calling by splitting scaffolds into chunks
and running multiple prodigal-gv CLI processes.
"""

from concurrent.futures import ThreadPoolExecutor
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import json
import logging
import re
import shutil
import signal
import subprocess
import tempfile
from typing import Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


_LONG_SCAFFOLD_BP = 10_000_000
_TILE_CORE_BP = 1_000_000
# Cleanup-abort recovery relies on this disposable flank. Do not shrink it
# without re-running the real long-scaffold regression gate.
_TILE_OVERLAP_BP = 25_000
_TILE_ID_PREFIX = "__virosync_tile_"
_SEQUENCE_DATA_RE = re.compile(
    r'^# Sequence Data: seqnum=\d+;seqlen=(\d+);seqhdr="(.*)"\s*$'
)
_KNOWN_CLEANUP_ERRORS = (
    "double free or corruption",
    "free(): invalid pointer",
    "munmap_chunk(): invalid pointer",
)


@dataclass
class GenePrediction:
    gene_id: str
    scaffold: str
    start: int
    end: int
    strand: str
    protein: str


@dataclass(frozen=True)
class _ProdigalValidation:
    gff_count: int
    faa_count: int
    discarded_coordinates: tuple[tuple[str, int, int, str], ...] = ()


def parse_prodigal_header(header: str, record_id: str) -> Optional[tuple[str, int, int, str]]:
    """Parse and normalize coordinates from a Prodigal-GV protein header.

    Prodigal header format example:
      >contig_1 # 143 # 456 # - # ID=1_1;...

    Prodigal reports 1-based inclusive coordinates. Returned coordinates use
    ViroSync's 0-based half-open ``[start, end)`` convention and strands are
    normalized to ``+`` or ``-``. Structurally unrelated headers return
    ``None``; malformed Prodigal-like metadata raises ``ValueError``.
    """
    parts = header.split(" # ")
    if len(parts) == 1:
        return None
    if len(parts) < 4:
        raise ValueError(
            f"Malformed Prodigal header for {record_id!r}: expected start, end, and strand"
        )
    try:
        start = int(parts[1])
        end = int(parts[2])
    except ValueError as exc:
        raise ValueError(
            f"Invalid Prodigal coordinate for {record_id!r}: "
            f"start={parts[1]!r}, end={parts[2]!r}"
        ) from exc

    if start < 1:
        raise ValueError(f"Invalid Prodigal start for {record_id!r}: {start}; expected >= 1")
    if end < start:
        raise ValueError(
            f"Invalid Prodigal end for {record_id!r}: {end}; expected >= start ({start})"
        )

    strand_token = parts[3].strip()
    strand_map = {"1": "+", "+": "+", "-1": "-", "-": "-"}
    try:
        strand = strand_map[strand_token]
    except KeyError as exc:
        raise ValueError(
            f"Invalid Prodigal strand for {record_id!r}: {strand_token!r}"
        ) from exc

    scaffold_fields = parts[0].split()
    if not scaffold_fields:
        raise ValueError(f"Invalid Prodigal scaffold for {record_id!r}: empty scaffold")
    scaffold_token = scaffold_fields[0]
    scaffold = scaffold_token
    gene_index = None
    if "ID=" in header:
        try:
            id_part = header.split("ID=", 1)[1].split(";", 1)[0]
            gene_index = id_part.split("_")[-1]
        except IndexError:
            gene_index = None
    if gene_index and scaffold_token.endswith(f"_{gene_index}"):
        scaffold = scaffold_token[: -(len(gene_index) + 1)]
    return scaffold, start - 1, end, strand


def _rebase_interval(start: int, end: int, offset: int) -> tuple[int, int]:
    """Move a normalized half-open interval into its parent scaffold."""
    return start + offset, end + offset


def _owns_midpoint(
    start: int,
    end: int,
    core_start: int,
    core_end: int,
) -> bool:
    """Return whether a half-open interval's midpoint belongs to a core."""
    return core_start <= (start + end) // 2 < core_end


def _format_prodigal_header(
    gene_id: str,
    start: int,
    end: int,
    strand: str,
    attributes: list[str],
) -> str:
    """Serialize a normalized interval using Prodigal's header coordinates."""
    strand_token = "1" if strand == "+" else "-1"
    return (
        f">{gene_id} # {start + 1} # {end} # "
        f"{strand_token} # {';'.join(attributes)}\n"
    )


def _contains_internal_tile(fasta_path: str) -> bool:
    """Return whether a FASTA contains a ViroSync tile record."""
    prefix = f">{_TILE_ID_PREFIX}"
    with open(fasta_path) as handle:
        return any(line.startswith(prefix) for line in handle if line.startswith(">"))


def _validate_tiled_prodigal_output(
    input_fasta: Path,
    proteins_faa: Path,
    genes_gff: Path,
    tile_cores: Optional[dict[str, tuple[int, int]]] = None,
    allow_unowned_gff_only: bool = False,
) -> _ProdigalValidation:
    """Require complete input metadata and matching FAA/GFF CDS coordinates."""
    if not proteins_faa.exists():
        raise RuntimeError(f"missing protein FASTA: {proteins_faa}")
    if not genes_gff.exists():
        raise RuntimeError(f"missing GFF: {genes_gff}")
    has_final_newline = True
    if proteins_faa.stat().st_size:
        with proteins_faa.open("rb") as handle:
            handle.seek(-1, 2)
            has_final_newline = handle.read(1) == b"\n"

    expected_sequences: Counter[tuple[str, int]] = Counter()
    input_ids: list[str] = []
    input_count = 0
    for record in SeqIO.parse(input_fasta, "fasta"):
        expected_sequences[(record.id, len(record.seq))] += 1
        input_ids.append(record.id)
        input_count += 1
    if len(expected_sequences) != input_count:
        raise RuntimeError(f"duplicate input FASTA IDs in {input_fasta}")

    observed_sequences: Counter[tuple[str, int]] = Counter()
    gff_coordinates: Counter[tuple[str, int, int, str]] = Counter()
    gff_order: list[tuple[str, int, int, str]] = []
    with genes_gff.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("# Sequence Data:"):
                match = _SEQUENCE_DATA_RE.match(line)
                if not match:
                    raise RuntimeError(
                        f"unparseable Sequence Data at {genes_gff}:{line_number}"
                    )
                length, header = match.groups()
                record_id = header.split(maxsplit=1)[0]
                observed_sequences[(record_id, int(length))] += 1
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            try:
                coordinate = (
                    fields[0],
                    int(fields[3]) - 1,
                    int(fields[4]),
                    fields[6],
                )
                gff_coordinates[coordinate] += 1
                gff_order.append(coordinate)
            except ValueError as exc:
                raise RuntimeError(
                    f"invalid CDS coordinates at {genes_gff}:{line_number}"
                ) from exc

    if observed_sequences != expected_sequences:
        raise RuntimeError(
            "GFF Sequence Data does not match the input FASTA: "
            f"expected {expected_sequences}, observed {observed_sequences}"
        )

    faa_coordinates: Counter[tuple[str, int, int, str]] = Counter()
    faa_order: list[tuple[str, int, int, str]] = []
    invalid_proteins: list[tuple[tuple[str, int, int, str], str]] = []
    malformed_final_header = False
    faa_records = list(SeqIO.parse(proteins_faa, "fasta"))
    for index, record in enumerate(faa_records):
        try:
            parsed = parse_prodigal_header(record.description, record.id)
        except ValueError as exc:
            if allow_unowned_gff_only and index == len(faa_records) - 1:
                malformed_final_header = True
                continue
            raise RuntimeError(
                f"malformed Prodigal protein header: {record.id}"
            ) from exc
        if parsed is None:
            raise RuntimeError(f"unparseable Prodigal protein header: {record.id}")
        scaffold, start, end, strand = parsed
        coordinate = (scaffold, start, end, strand)
        faa_coordinates[coordinate] += 1
        faa_order.append(coordinate)
        nucleotide_length = end - start
        protein = str(record.seq)
        if nucleotide_length % 3 or len(protein) != nucleotide_length // 3:
            invalid_proteins.append(
                (
                    coordinate,
                    f"protein length does not match CDS span for {record.id}: "
                    f"protein={len(protein)}, CDS={nucleotide_length} bp",
                )
            )
        if "*" in protein[:-1]:
            invalid_proteins.append(
                (coordinate, f"internal stop codon in Prodigal protein: {record.id}")
            )

    if not allow_unowned_gff_only:
        if not has_final_newline:
            raise RuntimeError(f"protein FASTA lacks a final newline: {proteins_faa}")
        if invalid_proteins:
            raise RuntimeError(invalid_proteins[0][1])
        if gff_coordinates != faa_coordinates:
            raise RuntimeError(
                "GFF and protein FASTA CDS coordinates differ: "
                f"GFF={sum(gff_coordinates.values())}, "
                f"FAA={sum(faa_coordinates.values())}"
            )
        return _ProdigalValidation(
            gff_count=sum(gff_coordinates.values()),
            faa_count=sum(faa_coordinates.values()),
        )

    if input_count != 1:
        raise RuntimeError("cleanup-abort recovery requires one isolated input record")
    tile_cores = tile_cores or {}
    record_id = input_ids[0]
    core = tile_cores.get(record_id)
    if core is None:
        raise RuntimeError(
            "cleanup-abort output for an untiled record cannot be accepted safely: "
            f"{record_id}"
        )
    faa_only = faa_coordinates - gff_coordinates
    if faa_only:
        raise RuntimeError(
            "protein FASTA contains coordinates absent from the GFF: "
            f"{list(faa_only)[:3]}"
        )

    affected = set(gff_coordinates - faa_coordinates)
    affected.update(coordinate for coordinate, _ in invalid_proteins)
    if not has_final_newline and not malformed_final_header:
        if not faa_order:
            raise RuntimeError("nonempty cleanup-abort output contains no proteins")
        affected.add(faa_order[-1])
    if malformed_final_header and not affected:
        raise RuntimeError(
            "malformed final protein header has no matching GFF suffix"
        )
    if not affected:
        return _ProdigalValidation(
            gff_count=sum(gff_coordinates.values()),
            faa_count=sum(faa_coordinates.values()),
        )

    if any(count != 1 for count in gff_coordinates.values()):
        raise RuntimeError("cleanup-abort recovery does not allow duplicate GFF coordinates")
    affected_indices = [
        index for index, coordinate in enumerate(gff_order) if coordinate in affected
    ]
    if not affected_indices or set(gff_order[affected_indices[0] :]) != affected:
        raise RuntimeError("cleanup-abort output loss is not a contiguous GFF suffix")

    core_start, core_end = core
    for scaffold, start, end, strand in affected:
        if scaffold != record_id:
            raise RuntimeError(
                f"cleanup-abort coordinate maps to unexpected record: {scaffold}"
            )
        if _owns_midpoint(start, end, core_start, core_end):
            raise RuntimeError(
                "cleanup-abort output is incomplete inside the owned core: "
                f"{scaffold}:{start}-{end}({strand})"
            )

    discarded = tuple(gff_order[affected_indices[0] :])
    return _ProdigalValidation(
        gff_count=sum(gff_coordinates.values()),
        faa_count=sum(faa_coordinates.values()),
        discarded_coordinates=discarded,
    )


def _retain_failed_prodigal_attempt(
    failure_dir: Path,
    paths: list[Path],
) -> Path:
    """Copy diagnostic artifacts outside the temporary worker directory."""
    failure_dir.mkdir(parents=True, exist_ok=True)
    for path in paths:
        if path.exists():
            shutil.copy2(path, failure_dir / path.name)
    return failure_dir


def _known_cleanup_failure(returncode: int, stderr_path: Path) -> bool:
    """Recognize the observed Prodigal-GV allocator cleanup abort."""
    if returncode != -signal.SIGABRT or not stderr_path.exists():
        return False
    stderr = stderr_path.read_text(errors="replace").lower()
    return any(signature in stderr for signature in _KNOWN_CLEANUP_ERRORS)


def _write_cleanup_abort_audit(
    audit_dir: Path,
    record_id: str,
    returncode: int,
    stderr_path: Path,
    validation: _ProdigalValidation,
) -> Path:
    """Record every accepted partial-output cleanup abort."""
    audit_dir.mkdir(parents=True, exist_ok=True)
    stderr = stderr_path.read_text(errors="replace").strip()
    payload = {
        "record_id": record_id,
        "returncode": returncode,
        "stderr": stderr,
        "gff_count": validation.gff_count,
        "faa_count": validation.faa_count,
        "discarded_coordinates": [
            {
                "scaffold": scaffold,
                "start_0based": start,
                "end_exclusive": end,
                "strand": strand,
            }
            for scaffold, start, end, strand in validation.discarded_coordinates
        ],
    }
    audit_path = audit_dir / f"{record_id}.json"
    temporary_path = audit_path.with_suffix(".json.tmp")
    temporary_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary_path.replace(audit_path)
    return audit_path


def _remove_discarded_proteins(
    proteins_faa: Path,
    discarded_coordinates: tuple[tuple[str, int, int, str], ...],
) -> None:
    """Remove an aborted output suffix before concatenating retry outputs."""
    if not discarded_coordinates:
        return
    discarded = set(discarded_coordinates)
    retained = []
    records = list(SeqIO.parse(proteins_faa, "fasta"))
    for index, record in enumerate(records):
        try:
            parsed = parse_prodigal_header(record.description, record.id)
        except ValueError:
            if index == len(records) - 1:
                continue
            raise
        if parsed is None:
            raise RuntimeError(f"unparseable Prodigal protein header: {record.id}")
        scaffold, start, end, strand = parsed
        if (scaffold, start, end, strand) not in discarded:
            retained.append(record)
    temporary_path = proteins_faa.with_suffix(".clean.faa")
    SeqIO.write(retained, temporary_path, "fasta")
    temporary_path.replace(proteins_faa)


def _run_prodigal_on_chunk(
    chunk_fasta: str,
    chunk_out: str,
    has_tiles: Optional[bool] = None,
    tile_cores: Optional[dict[str, tuple[int, int]]] = None,
) -> str:
    """Run prodigal-gv CLI on a single chunk FASTA. Returns protein FASTA path."""
    if has_tiles is None:
        has_tiles = _contains_internal_tile(chunk_fasta)
    tile_cores = tile_cores or {}
    cmd = [
        "prodigal-gv",
        "-i", chunk_fasta,
        "-a", chunk_out,
        "-o", str(Path(chunk_out).with_suffix(".gff")) if has_tiles else "/dev/null",
        "-f", "gff",
        "-p", "meta",
        "-q",
    ]
    if has_tiles:
        stderr_path = Path(chunk_out).with_suffix(".stderr")
        with stderr_path.open("w") as stderr:
            completed = subprocess.run(
                cmd,
                check=False,
                stdout=subprocess.DEVNULL,
                stderr=stderr,
            )
        if completed.returncode == 0:
            try:
                _validate_tiled_prodigal_output(
                    Path(chunk_fasta),
                    Path(chunk_out),
                    Path(cmd[cmd.index("-o") + 1]),
                )
            except RuntimeError:
                failure_dir = _retain_failed_prodigal_attempt(
                    Path(chunk_out).parent.parent
                    / "prodigal_failures"
                    / Path(chunk_out).stem,
                    [
                        Path(chunk_fasta),
                        Path(chunk_out),
                        Path(cmd[cmd.index("-o") + 1]),
                        stderr_path,
                    ],
                )
                logger.error("Retained invalid Prodigal output in %s", failure_dir)
                raise
            return chunk_out

        logger.warning(
            "prodigal-gv exited nonzero for tiled worker %s; retrying "
            "its records separately",
            chunk_fasta,
        )
        with tempfile.TemporaryDirectory(
            dir=Path(chunk_out).parent,
        ) as retry_dir_name:
            retry_dir = Path(retry_dir_name)
            retry_outputs: list[Path] = []
            for index, record in enumerate(SeqIO.parse(chunk_fasta, "fasta")):
                retry_input = retry_dir / f"record_{index}.fasta"
                retry_output = retry_dir / f"record_{index}.faa"
                retry_gff = retry_dir / f"record_{index}.gff"
                retry_stderr = retry_dir / f"record_{index}.stderr"
                SeqIO.write([record], retry_input, "fasta")
                retry_cmd = [
                    "prodigal-gv",
                    "-i", str(retry_input),
                    "-a", str(retry_output),
                    "-o", str(retry_gff),
                    "-f", "gff",
                    "-p", "meta",
                    "-q",
                ]
                with retry_stderr.open("w") as stderr:
                    retry_completed = subprocess.run(
                        retry_cmd,
                        check=False,
                        stdout=subprocess.DEVNULL,
                        stderr=stderr,
                    )
                try:
                    known_cleanup = (
                        retry_completed.returncode != 0
                        and _known_cleanup_failure(
                            retry_completed.returncode,
                            retry_stderr,
                        )
                    )
                    if retry_completed.returncode != 0 and not known_cleanup:
                        raise RuntimeError(
                            "unrecognized nonzero Prodigal-GV exit: "
                            f"{retry_completed.returncode}"
                        )
                    validation = _validate_tiled_prodigal_output(
                        retry_input,
                        retry_output,
                        retry_gff,
                        tile_cores=tile_cores,
                        allow_unowned_gff_only=known_cleanup,
                    )
                    if known_cleanup:
                        if retry_output.stat().st_size == 0:
                            raise RuntimeError(
                                "nonzero Prodigal-GV exit produced no proteins"
                            )
                        audit_path = _write_cleanup_abort_audit(
                            Path(chunk_out).parent.parent
                            / "accepted_cleanup_aborts",
                            record.id,
                            retry_completed.returncode,
                            retry_stderr,
                            validation,
                        )
                        _remove_discarded_proteins(
                            retry_output,
                            validation.discarded_coordinates,
                        )
                        logger.warning(
                            "Accepted owned-core-complete output after known "
                            "Prodigal-GV cleanup abort for %s; audit: %s",
                            record.id,
                            audit_path,
                        )
                except (RuntimeError, ValueError):
                    failure_dir = _retain_failed_prodigal_attempt(
                        Path(chunk_out).parent.parent
                        / "prodigal_failures"
                        / Path(chunk_out).stem
                        / f"record_{index}",
                        [retry_input, retry_output, retry_gff, retry_stderr],
                    )
                    logger.error("Retained failed Prodigal attempt in %s", failure_dir)
                    raise
                retry_outputs.append(retry_output)

            with open(chunk_out, "w") as merged:
                for retry_output in retry_outputs:
                    with retry_output.open() as handle:
                        shutil.copyfileobj(handle, merged)
        return chunk_out

    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        if Path(chunk_out).exists() and Path(chunk_out).stat().st_size > 0:
            logger.warning(
                "prodigal-gv exited nonzero after writing %s; accepting the "
                "existing untiled output for compatibility",
                chunk_out,
            )
        else:
            raise
    return chunk_out


def run_prodigal_genome(
    genome_fasta: Path,
    output_dir: Path,
    mode: str = "meta",
    threads: int = 1,
) -> tuple[Path, list[GenePrediction]]:
    """
    Run prodigal-gv on a genome and parse gene predictions.

    Long scaffolds are tiled even with one requested thread. With multiple
    threads, work records are balanced across parallel prodigal-gv processes.

    Returns:
        (proteins_faa_path, list of GenePrediction)
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    proteins_faa = output_dir / "proteome.fasta"
    gff_path = output_dir / "genes.gff"

    if shutil.which("prodigal-gv") is None:
        raise RuntimeError("prodigal-gv not found in PATH.")

    use_parallel = threads > 1
    if not use_parallel:
        with genome_fasta.open() as handle:
            use_parallel = any(
                len(record.seq) > _LONG_SCAFFOLD_BP
                for record in SeqIO.parse(handle, "fasta")
            )
    if use_parallel:
        return _run_prodigal_parallel(
            genome_fasta, output_dir, proteins_faa, gff_path, max(1, threads)
        )
    return _run_prodigal_single(
        genome_fasta, output_dir, proteins_faa, gff_path, mode
    )


def _run_prodigal_parallel(
    genome_fasta: Path,
    output_dir: Path,
    proteins_faa: Path,
    gff_path: Path,
    threads: int,
) -> tuple[Path, list[GenePrediction]]:
    """Split genome into chunks and run prodigal-gv in parallel."""
    # Read all scaffolds and split long records into bounded work units. The
    # merged proteome is restored to genome-global coordinates before Phase 1.
    records = list(SeqIO.parse(genome_fasta, "fasta"))
    reserved_ids = [
        record.id for record in records if record.id.startswith(_TILE_ID_PREFIX)
    ]
    if reserved_ids:
        raise RuntimeError(
            "input scaffold ID uses ViroSync's reserved tile prefix: "
            f"{reserved_ids[0]}"
        )
    scaffold_order = {record.id: index for index, record in enumerate(records)}
    work_records: list[SeqRecord] = []
    tile_sources: dict[str, tuple[str, int, int, int]] = {}
    tiled_scaffolds = 0

    for record_index, record in enumerate(records):
        if len(record.seq) <= _LONG_SCAFFOLD_BP:
            work_records.append(record)
            continue

        tiled_scaffolds += 1
        tile_count = (len(record.seq) + _TILE_CORE_BP - 1) // _TILE_CORE_BP
        for tile_index in range(tile_count):
            # Even partitions avoid a tiny, context-poor terminal tile.
            core_start = len(record.seq) * tile_index // tile_count
            core_end = len(record.seq) * (tile_index + 1) // tile_count
            tile_start = max(0, core_start - _TILE_OVERLAP_BP)
            tile_end = min(len(record.seq), core_end + _TILE_OVERLAP_BP)
            tile_id = (
                f"{_TILE_ID_PREFIX}{record_index:06d}_{tile_index:06d}"
            )
            work_records.append(
                SeqRecord(
                    record.seq[tile_start:tile_end],
                    id=tile_id,
                    description="",
                )
            )
            tile_sources[tile_id] = (
                record.id,
                tile_start,
                core_start,
                core_end,
            )

    n_chunks = min(threads, len(work_records))

    logger.info(
        "prodigal-gv parallel: %d scaffolds (%d tiled) split into %d chunks "
        "(%d work records, %d threads)",
        len(records), tiled_scaffolds, n_chunks, len(work_records), threads,
    )

    with tempfile.TemporaryDirectory(dir=output_dir) as tmpdir:
        tmpdir = Path(tmpdir)
        chunk_inputs = []
        chunk_outputs = []
        chunk_has_tiles = []
        chunk_tile_cores = []

        # Distribute scaffolds round-robin by size (largest first)
        records_sorted = sorted(work_records, key=lambda r: len(r.seq), reverse=True)
        chunks: list[list] = [[] for _ in range(n_chunks)]
        chunk_sizes = [0] * n_chunks
        for record in records_sorted:
            # Assign to smallest chunk (greedy balancing)
            idx = chunk_sizes.index(min(chunk_sizes))
            chunks[idx].append(record)
            chunk_sizes[idx] += len(record.seq)

        for i, chunk_records in enumerate(chunks):
            if not chunk_records:
                continue
            chunk_fasta = tmpdir / f"chunk_{i}.fasta"
            chunk_out = tmpdir / f"chunk_{i}.faa"
            SeqIO.write(chunk_records, chunk_fasta, "fasta")
            chunk_inputs.append(str(chunk_fasta))
            chunk_outputs.append(str(chunk_out))
            cores = {
                record.id: (
                    tile_sources[record.id][2] - tile_sources[record.id][1],
                    tile_sources[record.id][3] - tile_sources[record.id][1],
                )
                for record in chunk_records
                if record.id in tile_sources
            }
            # Strict FAA/GFF validation is genome-scoped, so a short scaffold's
            # behavior cannot change with thread-dependent chunk packing.
            chunk_has_tiles.append(bool(tile_sources))
            chunk_tile_cores.append(cores)

        # Run prodigal-gv in parallel
        with ThreadPoolExecutor(max_workers=n_chunks) as executor:
            list(
                executor.map(
                    _run_prodigal_on_chunk,
                    chunk_inputs,
                    chunk_outputs,
                    chunk_has_tiles,
                    chunk_tile_cores,
                )
            )

        genes: list[GenePrediction] = []
        if not tile_sources:
            # Preserve the existing output path exactly for ordinary genomes.
            with open(proteins_faa, "w") as out_f:
                for chunk_out in chunk_outputs:
                    if not Path(chunk_out).exists():
                        continue
                    with open(chunk_out) as in_f:
                        out_f.write(in_f.read())

            for record in SeqIO.parse(proteins_faa, "fasta"):
                parsed = parse_prodigal_header(record.description, record.id)
                if not parsed:
                    continue
                scaffold, start, end, strand = parsed
                genes.append(
                    GenePrediction(
                        gene_id=record.id,
                        scaffold=scaffold,
                        start=start,
                        end=end,
                        strand=strand,
                        protein=str(record.seq),
                    )
                )
        else:
            predictions: dict[
                str, list[tuple[int, int, str, str, str]]
            ] = {}
            for chunk_out in chunk_outputs:
                if not Path(chunk_out).exists():
                    continue
                for record in SeqIO.parse(chunk_out, "fasta"):
                    parsed = parse_prodigal_header(record.description, record.id)
                    if not parsed:
                        continue
                    scaffold, start, end, strand = parsed
                    metadata_parts = record.description.split(" # ", 4)
                    metadata = metadata_parts[4] if len(metadata_parts) == 5 else ""

                    tile_source = tile_sources.get(scaffold)
                    if tile_source:
                        scaffold, tile_start, core_start, core_end = tile_source
                        start, end = _rebase_interval(start, end, tile_start)
                        if not _owns_midpoint(start, end, core_start, core_end):
                            continue

                    predictions.setdefault(scaffold, []).append(
                        (start, end, strand, str(record.seq), metadata)
                    )

            unknown_scaffolds = sorted(set(predictions) - set(scaffold_order))
            if unknown_scaffolds:
                raise RuntimeError(
                    "prodigal-gv output could not be mapped to input scaffolds: "
                    + ", ".join(unknown_scaffolds[:5])
                )

            with open(proteins_faa, "w") as out_f:
                for scaffold in sorted(
                    predictions,
                    key=lambda item: scaffold_order[item],
                ):
                    scaffold_predictions = sorted(
                        predictions[scaffold],
                        key=lambda item: (item[0], item[1], item[2], item[3]),
                    )
                    for gene_index, (
                        start,
                        end,
                        strand,
                        protein,
                        metadata,
                    ) in enumerate(scaffold_predictions, start=1):
                        gene_id = f"{scaffold}_{gene_index}"
                        attributes = metadata.split(";") if metadata else []
                        for index, attribute in enumerate(attributes):
                            if attribute.startswith("ID="):
                                attributes[index] = f"ID=1_{gene_index}"
                                break
                        else:
                            attributes.insert(0, f"ID=1_{gene_index}")
                        out_f.write(
                            _format_prodigal_header(
                                gene_id,
                                start,
                                end,
                                strand,
                                attributes,
                            )
                        )
                        for offset in range(0, len(protein), 60):
                            out_f.write(protein[offset : offset + 60] + "\n")
                        genes.append(
                            GenePrediction(
                                gene_id=gene_id,
                                scaffold=scaffold,
                                start=start,
                                end=end,
                                strand=strand,
                                protein=protein,
                            )
                        )

    # Write an empty GFF placeholder (parallel mode doesn't merge GFFs)
    gff_path.touch()

    logger.info("prodigal-gv parallel: predicted %d genes", len(genes))
    return proteins_faa, genes


def _run_prodigal_single(
    genome_fasta: Path,
    output_dir: Path,
    proteins_faa: Path,
    gff_path: Path,
    mode: str,
) -> tuple[Path, list[GenePrediction]]:
    """Run prodigal-gv CLI on entire genome (single process)."""
    # genes.gff is an upstream-native Prodigal artifact and therefore retains
    # standard 1-based GFF coordinates. Only header-derived in-memory records
    # are normalized to ViroSync's 0-based half-open convention.
    cmd = [
        "prodigal-gv",
        "-i", str(genome_fasta),
        "-a", str(proteins_faa),
        "-o", str(gff_path),
        "-f", "gff",
        "-p", mode,
        "-q",
    ]
    logger.info("Running prodigal-gv: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        if proteins_faa.exists() and proteins_faa.stat().st_size > 0:
            logger.warning(
                "prodigal-gv exited with error %s, but output files exist. "
                "Treating as success (known cleanup issue).", e.returncode
            )
        else:
            raise

    genes: list[GenePrediction] = []
    for record in SeqIO.parse(proteins_faa, "fasta"):
        parsed = parse_prodigal_header(record.description, record.id)
        if not parsed:
            continue
        scaffold, start, end, strand = parsed
        genes.append(
            GenePrediction(
                gene_id=record.id,
                scaffold=scaffold,
                start=start,
                end=end,
                strand=strand,
                protein=str(record.seq),
            )
        )

    logger.info("Parsed %d prodigal genes", len(genes))
    return proteins_faa, genes


def load_gene_predictions(proteins_faa: Path) -> dict[str, list[GenePrediction]]:
    """
    Load gene predictions from a prodigal-gv proteins FASTA.
    """
    by_scaffold: dict[str, list[GenePrediction]] = {}
    for record in SeqIO.parse(proteins_faa, "fasta"):
        parsed = parse_prodigal_header(record.description, record.id)
        if not parsed:
            continue
        scaffold, start, end, strand = parsed
        by_scaffold.setdefault(scaffold, []).append(
            GenePrediction(
                gene_id=record.id,
                scaffold=scaffold,
                start=start,
                end=end,
                strand=strand,
                protein=str(record.seq),
            )
        )
    for scaffold in by_scaffold:
        by_scaffold[scaffold].sort(key=lambda g: g.start)
    return by_scaffold
