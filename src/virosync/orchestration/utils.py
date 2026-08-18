"""
ViroSync Orchestration Utilities.

Helper functions for data wiring between orchestration task functions.
Shared helpers keep task inputs on disk instead of serializing large objects.
"""

from bisect import bisect_left
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from threading import Lock
from typing import Any, Callable, Optional

from Bio import SeqIO
from virosync.pipeline.phase0.prodigal import parse_prodigal_header
from virosync.orchestration.resource_monitor import ResourceMonitor


_ProteomeRecord = tuple[str, str, int, int, str]
_ScaffoldIndex = tuple[
    tuple[int, ...],
    tuple[tuple[int, _ProteomeRecord], ...],
    int,
]
_ProteomeData = tuple[
    tuple[_ProteomeRecord, ...],
    dict[str, _ScaffoldIndex],
]


def _parse_proteome_records(
    proteome_path: str, _size: int, _mtime: int
) -> tuple[_ProteomeRecord, ...]:
    """Parse a proteome FASTA once into (id, scaffold, start, end, seq) tuples.

    The loader caches these records by (path, size, mtime), so per-boundary
    lookups parse the proteome once. Coordinates and IDs match the prior inline
    parse.
    """
    records: list[tuple[str, str, int, int, str]] = []
    for record in SeqIO.parse(proteome_path, "fasta"):
        parsed = parse_prodigal_header(record.description, record.id)
        if parsed:
            scaffold, start, end, _strand = parsed
        else:
            # Parse gene header to extract coordinates
            header_parts = record.id.split("_")
            try:
                if header_parts[0].lower() == "porf":
                    scaffold = header_parts[1]
                    start = int(header_parts[2])
                    end = int(header_parts[3])
                else:
                    scaffold = "_".join(header_parts[:-4])
                    start = int(header_parts[-4])
                    end = int(header_parts[-3])
            except (ValueError, IndexError):
                scaffold = record.id
                start = 0
                end = len(record.seq) * 3
        records.append((record.id, scaffold, start, end, str(record.seq)))
    return tuple(records)


def _build_proteome_index(
    records: tuple[_ProteomeRecord, ...],
) -> dict[str, _ScaffoldIndex]:
    by_scaffold: dict[str, list[tuple[int, _ProteomeRecord]]] = defaultdict(list)
    for file_index, record in enumerate(records):
        by_scaffold[record[1]].append((file_index, record))

    index = {}
    for scaffold, file_ordered in by_scaffold.items():
        coordinate_ordered = tuple(
            sorted(file_ordered, key=lambda item: (item[1][2], item[0]))
        )
        index[scaffold] = (
            tuple(record[2] for _file_index, record in coordinate_ordered),
            coordinate_ordered,
            max(
                (
                    max(0, record[3] - record[2])
                    for _file_index, record in file_ordered
                ),
                default=0,
            ),
        )
    return index


@lru_cache(maxsize=4)
def _load_proteome_data(
    proteome_path: str, size: int, mtime: int
) -> _ProteomeData:
    records = _parse_proteome_records(proteome_path, size, mtime)
    return records, _build_proteome_index(records)


_PROTEOME_CACHE_LOCK = Lock()


def _cached_proteome_data(
    proteome_path: Path,
) -> _ProteomeData:
    p = Path(proteome_path)
    try:
        st = p.stat()
    except OSError:
        records = _parse_proteome_records(str(p), 0, 0)
        return records, _build_proteome_index(records)

    # functools.lru_cache can execute the wrapped function more than once when
    # concurrent misses race. Hold this lock around the cache-wrapper call so a
    # 32-thread Phase 3 start parses one large proteome exactly once.
    with _PROTEOME_CACHE_LOCK:
        return _load_proteome_data(str(p), st.st_size, int(st.st_mtime))


def get_overlapping_genes(
    proteome_path: Path,
    boundary_scaffold: Optional[str] = None,
    boundary_start: Optional[int] = None,
    boundary_end: Optional[int] = None,
) -> dict[str, list[tuple[str, str]]]:
    """
    Get gene sequences grouped by scaffold.

    If boundary coordinates are provided, filters to genes overlapping that region.

    Args:
        proteome_path: Path to proteome FASTA file
        boundary_scaffold: If provided, filter to this scaffold
        boundary_start: If provided with boundary_end, filter to overlapping region
        boundary_end: If provided with boundary_start, filter to overlapping region

    Returns:
        Dictionary mapping scaffold ID to list of (gene_id, sequence) tuples
    """
    genes_by_scaffold: dict[str, list[tuple[str, str]]] = defaultdict(list)
    records, scaffold_index = _cached_proteome_data(proteome_path)

    if boundary_scaffold:
        indexed = scaffold_index.get(boundary_scaffold)
        if indexed is None:
            return {}
        starts, coordinate_ordered, max_gene_length = indexed
        if boundary_start is not None and boundary_end is not None:
            # Every overlapping record starts after this exact lower bound.
            left = bisect_left(starts, boundary_start - max_gene_length)
            right = bisect_left(starts, boundary_end)
            selected = [
                item
                for item in coordinate_ordered[left:right]
                if item[1][3] > boundary_start
            ]
            selected.sort(key=lambda item: item[0])
            records = tuple(record for _file_index, record in selected)

    for record_id, scaffold, start, end, seq in records:
        # Filter by scaffold if specified
        if boundary_scaffold and scaffold != boundary_scaffold:
            continue

        # Filter by overlap if boundary specified
        if boundary_start is not None and boundary_end is not None:
            if end <= boundary_start or start >= boundary_end:
                continue  # No overlap

        genes_by_scaffold[scaffold].append((record_id, seq))

    return dict(genes_by_scaffold)


def get_genes_for_boundary(
    proteome_path: Path,
    scaffold: str,
    start: int,
    end: int,
    max_porfs: int = 10000,
) -> list[tuple[str, str]]:
    """
    Get gene sequences overlapping a specific boundary.

    Args:
        proteome_path: Path to proteome FASTA file
        scaffold: Scaffold/contig ID
        start: Boundary start coordinate
        end: Boundary end coordinate
        max_porfs: Maximum number of genes to return (default: 10000).
                   TMVec batch processing can handle thousands efficiently.
                   Note: Boltz structural analysis has a separate limit
                   (structural_max_porfs) applied in evidence_synthesizer.py.

    Returns:
        List of (gene_id, sequence) tuples
    """
    genes = get_overlapping_genes(
        proteome_path,
        boundary_scaffold=scaffold,
        boundary_start=start,
        boundary_end=end,
    )

    scaffold_genes = genes.get(scaffold, [])

    # Return all genes up to max_porfs
    # TMVec batch processing handles thousands of proteins efficiently
    # Boltz has a separate limit (structural_max_porfs=50) applied later
    return scaffold_genes[:max_porfs]


def run_with_monitor(
    task_name: str,
    func: Callable[..., Any],
    monitor_output_dir: Path,
    threads: int,
    phase: str = "phase3",
    task_id: Optional[str] = None,
    genome_path: Optional[Path] = None,
    genome_fasta: Optional[Path] = None,
    **kwargs: Any,
) -> Any:
    """
    Run a function under ResourceMonitor with best-effort genome context.
    """
    genome_source = genome_path or genome_fasta
    genome_id = Path(genome_source).stem if genome_source else "unknown_genome"
    monitor_output_dir = Path(monitor_output_dir)
    monitor_output_dir.mkdir(parents=True, exist_ok=True)

    if genome_fasta is not None and "genome_fasta" not in kwargs:
        kwargs["genome_fasta"] = genome_fasta
    if genome_path is not None and "genome_path" not in kwargs:
        kwargs["genome_path"] = genome_path
    # Pass threads to the function if it accepts it
    if "threads" not in kwargs:
        kwargs["threads"] = threads

    with ResourceMonitor(
        task_name=task_name,
        genome_id=genome_id,
        phase=phase,
        output_dir=monitor_output_dir,
        threads=threads,
        task_id=task_id,
    ):
        return func(**kwargs)
