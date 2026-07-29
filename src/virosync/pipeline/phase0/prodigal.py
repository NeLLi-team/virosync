"""
Genome-wide gene calling with Prodigal-GV.

Supports parallel gene calling by splitting scaffolds into chunks
and running multiple prodigal-gv CLI processes.
"""

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
import logging
import shutil
import subprocess
import tempfile
from typing import Optional

from Bio import SeqIO

logger = logging.getLogger(__name__)


@dataclass
class GenePrediction:
    gene_id: str
    scaffold: str
    start: int
    end: int
    strand: str
    protein: str


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


def _run_prodigal_on_chunk(chunk_fasta: str, chunk_out: str) -> str:
    """Run prodigal-gv CLI on a single chunk FASTA. Returns protein FASTA path."""
    cmd = [
        "prodigal-gv",
        "-i", chunk_fasta,
        "-a", chunk_out,
        "-o", "/dev/null",
        "-f", "gff",
        "-p", "meta",
        "-q",
    ]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        if Path(chunk_out).exists() and Path(chunk_out).stat().st_size > 0:
            pass  # known cleanup issue
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

    When threads > 1, splits scaffolds into chunks and runs parallel
    prodigal-gv CLI processes for ~linear speedup.

    Returns:
        (proteins_faa_path, list of GenePrediction)
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    proteins_faa = output_dir / "proteome.fasta"
    gff_path = output_dir / "genes.gff"

    if shutil.which("prodigal-gv") is None:
        raise RuntimeError("prodigal-gv not found in PATH.")

    if threads > 1:
        return _run_prodigal_parallel(
            genome_fasta, output_dir, proteins_faa, gff_path, threads
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
    # Read all scaffolds and split into roughly equal chunks
    records = list(SeqIO.parse(genome_fasta, "fasta"))
    n_chunks = min(threads, len(records))

    logger.info(
        "prodigal-gv parallel: %d scaffolds split into %d chunks (%d threads)",
        len(records), n_chunks, threads,
    )

    with tempfile.TemporaryDirectory(dir=output_dir) as tmpdir:
        tmpdir = Path(tmpdir)
        chunk_inputs = []
        chunk_outputs = []

        # Distribute scaffolds round-robin by size (largest first)
        records_sorted = sorted(records, key=lambda r: len(r.seq), reverse=True)
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

        # Run prodigal-gv in parallel
        with ThreadPoolExecutor(max_workers=n_chunks) as executor:
            list(executor.map(_run_prodigal_on_chunk, chunk_inputs, chunk_outputs))

        # Concatenate results and parse
        genes: list[GenePrediction] = []
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
