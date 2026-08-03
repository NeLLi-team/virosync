"""
Gene-level taxonomy via Diamond for final EVE confidence scoring.

After candidate EVE regions are identified, this module runs Diamond on ALL genes
within each region to provide per-gene taxonomy. This enables:

1. Boosting confidence if NCLDV/MIRUS genes found in region
2. Penalizing confidence if high-identity EUK genes found
3. Detecting chimeric regions (mixed viral + host genes)

This is Step 9 from PIPELINE_HMM_GATED_PLAN.md.
"""

import logging
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from Bio import SeqIO

from virosync.output_contract import canonical_family
from virosync.pipeline.phase0.prodigal import parse_prodigal_header
from virosync.pipeline.taxonomy_utils import resolve_org_id
from virosync.utils.path_safety import (
    require_strict_child,
    safe_filename_component,
    safe_filename_components,
)

logger = logging.getLogger(__name__)

VIRAL_PREFIXES = {"NCLDV", "MIRUS", "VP", "PLV", "PPV", "CRESS", "GVMAG", "PHAGE"}
MIN_VIRAL_HIT_PIDENT = 25.0

# Reference namespaces that name a viral lineage. GVMAG and PHAGE are namespaces
# rather than published classes, so callers that need a class fold them via
# viral_hit_categories().
PUBLISHED_VIRAL_PREFIXES = frozenset(
    {"NCLDV", "MIRUS", "PPV", "CRESS", "GVMAG", "PHAGE"}
)
_PREPLASMIVIRICOTA_TOKEN = "preplasmiviricota"


@dataclass
class GeneDiamondHit:
    """Represents a single Diamond hit for a gene."""

    porf_id: str
    target: str
    evalue: float
    bits: float
    pident: float
    qcov: float
    rank: int = 1  # Rank in top-10 hits


@dataclass
class GeneTaxonomy:
    """Per-gene taxonomy summary."""

    porf_id: str
    porf_start: int
    porf_end: int
    top1_target: str = "."
    top1_prefix: str = "UNKNOWN"
    top1_pident: float = 0.0
    top1_evalue: float = 1.0
    top10_prefixes: list[str] = None
    top10_raw_prefixes: list[str] = None
    top10_targets: list[str] = None
    top10_bitscores: list[float] = None
    top10_pidents: list[float] = None
    top10_evalues: list[float] = None
    has_ncldv_mirus: bool = False  # Strict NCLDV/MIRUS only
    has_vp_plv: bool = False
    has_viral: bool = False
    is_high_pident_euk: bool = False

    def __post_init__(self):
        if self.top10_prefixes is None:
            self.top10_prefixes = []
        if self.top10_raw_prefixes is None:
            self.top10_raw_prefixes = []
        if self.top10_targets is None:
            self.top10_targets = []
        if self.top10_bitscores is None:
            self.top10_bitscores = []
        if self.top10_pidents is None:
            self.top10_pidents = []
        if self.top10_evalues is None:
            self.top10_evalues = []


def extract_raw_prefix(target_id: str) -> str:
    """Return the raw taxonomy prefix for a target id (e.g., VP, PLV)."""
    if "__" not in target_id:
        return "UNKNOWN"
    return target_id.split("__", 1)[0]


def extract_prefix(target_id: str) -> str:
    """
    Extract taxonomic prefix from target ID.

    Prefixes:
    - EUK__: Eukaryote
    - BAC__: Bacteria
    - ARC__: Archaea
    - MITO__: Mitochondria (maps to EUK)
    - PLASTID__: Plastids (maps to EUK)
    - NCLDV__: Giant viruses (Nucleocytoviricota)
    - MIRUS__: Mirusviricota
    - GVMAG__: Giant virus MAGs
    - PHAGE__: Bacteriophages
    - VP__: Virophages
    - PLV__: Polinton-like viruses
    - CRESS__: CRESS/ssDNA viruses

    Returns:
        Prefix string (e.g., "NCLDV", "EUK", "PLV", "VP", "CRESS", "UNKNOWN")
    """
    if target_id.startswith("EUK__"):
        return "EUK"
    elif target_id.startswith("BAC__"):
        return "BAC"
    elif target_id.startswith("ARC__"):
        return "ARC"
    elif target_id.startswith("MITO__"):
        return "EUK"  # Mitochondria are eukaryotic
    elif target_id.startswith("PLASTID__"):
        return "EUK"  # Plastids are eukaryotic
    elif target_id.startswith("NCLDV__"):
        return "NCLDV"
    elif target_id.startswith("MIRUS__"):
        return "MIRUS"
    elif target_id.startswith("GVMAG__"):
        return "GVMAG"
    elif target_id.startswith("PHAGE__"):
        return "PHAGE"
    elif target_id.startswith("VP__"):
        return "VP"
    elif target_id.startswith("PLV__"):
        return "PLV"
    elif target_id.startswith("PPV__"):
        return "PPV"
    elif target_id.startswith("CRESS__"):
        return "CRESS"
    else:
        return "UNKNOWN"


def summarize_dominant_family(
    taxonomies: list[GeneTaxonomy],
) -> tuple[str, float]:
    """
    Return the dominant viral family for a region and the fraction of genes supporting it.

    Families are counted by top-10 prefix membership. A region whose genes carry no
    viral prefix has no dominant family, so it reports ("UNKNOWN", 0.0) rather than the
    first family in iteration order. Callers export this label and feed it to family
    agreement scoring, so a placeholder family would read as real evidence.
    """
    if not taxonomies:
        return "UNKNOWN", 0.0
    family_counts = {
        family: 0 for family in ("NCLDV", "MIRUS", "PPV", "CRESS")
    }
    for taxonomy in taxonomies:
        supported = {
            canonical_family(prefix.rstrip("_"))
            for prefix, pident in zip(
                taxonomy.top10_prefixes or [],
                taxonomy.top10_pidents or [],
            )
            if pident >= MIN_VIRAL_HIT_PIDENT
        }
        for family in family_counts:
            if family in supported:
                family_counts[family] += 1
    dominant_family = max(family_counts, key=family_counts.get)
    dominant_count = family_counts[dominant_family]
    if dominant_count == 0:
        return "UNKNOWN", 0.0
    return dominant_family, dominant_count / len(taxonomies)


def _split_top10_field(value: object) -> list[str]:
    """Split a top-10 field that may be a comma-joined string or a list."""
    if isinstance(value, str):
        return value.split(",")
    if isinstance(value, (list, tuple)):
        return list(value)
    return []


def qualified_viral_hits(record: dict) -> list[tuple[str, str]]:
    """Return ``(prefix, target)`` for identity-qualified viral top-10 hits.

    One place canonicalizes a top-10 hit list into viral reference namespaces:
    VP and PLV fold onto PPV, and a hit below ``MIN_VIRAL_HIT_PIDENT`` does not
    qualify. ``record`` needs ``top10_prefixes`` and ``top10_pidents``;
    ``top10_targets`` is optional and yields ``""`` when absent.
    """
    prefixes = [
        canonical_family(str(prefix).rstrip("_"))
        for prefix in _split_top10_field(record.get("top10_prefixes"))
        if str(prefix).strip()
    ]
    pidents: list[float] = []
    for value in _split_top10_field(record.get("top10_pidents")):
        try:
            pidents.append(float(value))
        except (TypeError, ValueError):
            pidents.append(0.0)
    targets = [
        str(target).strip()
        for target in _split_top10_field(record.get("top10_targets"))
    ]
    hits: list[tuple[str, str]] = []
    for index, (prefix, pident) in enumerate(zip(prefixes, pidents)):
        if prefix not in PUBLISHED_VIRAL_PREFIXES:
            continue
        if pident < MIN_VIRAL_HIT_PIDENT:
            continue
        hits.append((prefix, targets[index] if index < len(targets) else ""))
    return hits


def viral_hit_categories(
    record: dict,
    taxonomy_lookup: Optional[dict] = None,
) -> list[str]:
    """Return the published viral class of each qualified top-10 hit.

    GVMAG is the giant-virus MAG namespace, so it reports NCLDV: without the
    fold a marker with 8 GVMAG and 2 NCLDV hits would read as two lineages. A
    PHAGE target reports PPV when its resolved lineage is Preplasmiviricota (the
    legacy ``PHAGE__VARDNA__`` namespace holds virophage-like genomes) and PHAGE
    otherwise. ``taxonomy_lookup`` of ``None`` skips lineage resolution.
    """
    categories: list[str] = []
    for prefix, target in qualified_viral_hits(record):
        if prefix == "GVMAG":
            categories.append("NCLDV")
        elif prefix == "PHAGE":
            categories.append(
                "PPV"
                if _is_preplasmiviricota(target, taxonomy_lookup)
                else "PHAGE"
            )
        else:
            categories.append(prefix)
    return categories


def _is_preplasmiviricota(target: str, taxonomy_lookup: Optional[dict]) -> bool:
    if not target or not taxonomy_lookup:
        return False
    org_id = resolve_org_id(target, taxonomy_lookup)
    lineage = str(taxonomy_lookup.get(org_id) or "").lower()
    return _PREPLASMIVIRICOTA_TOKEN in lineage


def _resolve_diamond_db_prefix(target_path: Path) -> Path:
    """
    Resolve a Diamond database prefix from a provided path.

    Accepts:
    - A .dmnd path
    - A prefix (without .dmnd)
    - A FASTA path if a matching .dmnd exists alongside it
    """
    if target_path.suffix == ".dmnd":
        prefix = target_path.with_suffix("")
    elif target_path.suffix in (".faa", ".fasta", ".fa"):
        dmnd_path = target_path.with_suffix(".dmnd")
        if not dmnd_path.exists():
            raise FileNotFoundError(f"Diamond DB not found for FASTA: {dmnd_path}")
        prefix = dmnd_path.with_suffix("")
    else:
        dmnd_path = Path(str(target_path) + ".dmnd")
        prefix = target_path
        if not dmnd_path.exists():
            raise FileNotFoundError(f"Diamond DB not found: {dmnd_path}")

    dmnd_path = Path(str(prefix) + ".dmnd")
    if not dmnd_path.exists():
        raise FileNotFoundError(f"Diamond DB not found: {dmnd_path}")
    return prefix


def run_diamond_blastp(
    query_fasta: Path,
    target_db: Path,
    output_file: Path,
    threads: int = 8,
    evalue: float = 1e-5,
    max_seqs: int = 10,
    chunk_size: int = 10000,
    fast: bool = False,
    min_threads_per_chunk: int = 4,
    search_backend: str = "diamond",
) -> Path:
    """
    Run Diamond blastp with top-10 hits.

    For large query sets (>chunk_size), queries are chunked and processed in parallel.
    Parallelism is determined by: max_parallel = threads // min_threads_per_chunk
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)
    db_prefix = _resolve_diamond_db_prefix(target_db)

    # Count queries to decide if chunking is needed
    query_count = sum(1 for line in open(query_fasta) if line.startswith(">"))
    if output_file.exists():
        output_file.unlink()

    if query_count <= chunk_size:
        _run_single_diamond_blastp(
            query_fasta=query_fasta,
            db_prefix=db_prefix,
            output_file=output_file,
            threads=threads,
            evalue=evalue,
            max_seqs=max_seqs,
            fast=fast,
            search_backend=search_backend,
        )
        return output_file

    # Calculate parallel execution parameters
    max_parallel = max(1, threads // min_threads_per_chunk)
    threads_per_chunk = max(1, threads // max_parallel)

    queries = list(SeqIO.parse(query_fasta, "fasta"))
    total_chunks = (len(queries) + chunk_size - 1) // chunk_size

    logger.info(
        "Chunking %d queries into %d batches of %d, running %d in parallel with %d threads each",
        query_count, total_chunks, chunk_size, max_parallel, threads_per_chunk,
    )

    with tempfile.TemporaryDirectory(dir=output_file.parent) as tmp_dir:
        tmp_path = Path(tmp_dir)

        # Prepare all chunk FASTAs first
        chunk_files = []
        for i in range(0, len(queries), chunk_size):
            chunk = queries[i : i + chunk_size]
            chunk_num = i // chunk_size + 1
            chunk_fasta = tmp_path / f"chunk_{chunk_num}.faa"
            chunk_output = tmp_path / f"chunk_{chunk_num}.tsv"

            with chunk_fasta.open("w") as handle:
                for record in chunk:
                    handle.write(f">{record.id}\n{record.seq}\n")

            chunk_files.append((chunk_num, chunk_fasta, chunk_output, len(chunk)))

        def process_chunk(args):
            chunk_num, chunk_fasta, chunk_output, chunk_len = args
            logger.info("Processing chunk %d/%d (%d queries)", chunk_num, total_chunks, chunk_len)
            _run_single_diamond_blastp(
                query_fasta=chunk_fasta,
                db_prefix=db_prefix,
                output_file=chunk_output,
                threads=threads_per_chunk,
                evalue=evalue,
                max_seqs=max_seqs,
                fast=fast,
                search_backend=search_backend,
            )
            return chunk_num, chunk_output

        # Process chunks in parallel
        completed_chunks = []
        with ThreadPoolExecutor(max_workers=max_parallel) as executor:
            futures = {executor.submit(process_chunk, cf): cf[0] for cf in chunk_files}
            for future in as_completed(futures):
                chunk_num, chunk_output = future.result()
                completed_chunks.append((chunk_num, chunk_output))
                logger.info("Completed chunk %d/%d", chunk_num, total_chunks)

        # Merge results in order
        completed_chunks.sort(key=lambda x: x[0])
        for chunk_num, chunk_output in completed_chunks:
            if chunk_output.exists():
                with output_file.open("a") as out_handle, chunk_output.open() as in_handle:
                    out_handle.write(in_handle.read())

    return output_file


def _run_single_diamond_blastp(
    query_fasta: Path,
    db_prefix: Path,
    output_file: Path,
    threads: int,
    evalue: float,
    max_seqs: int,
    fast: bool,
    search_backend: str = "diamond",
) -> None:
    from virosync.pipeline.search_backend import run_sequence_search

    extra_flags = ["--fast"] if fast else None  # Diamond is the only backend

    run_sequence_search(
        query_fasta=query_fasta,
        db_path=db_prefix,
        output_tsv=output_file,
        threads=threads,
        backend=search_backend,
        evalue=evalue,
        max_target_seqs=max_seqs,
        extra_flags=extra_flags,
    )


def parse_diamond_top10(diamond_output: Path) -> dict[str, list[GeneDiamondHit]]:
    """
    Parse Diamond output and extract top-10 hits per query.

    Args:
        diamond_output: Path to Diamond output TSV

    Returns:
        Dict mapping porf_id -> list of GeneDiamondHit (sorted by bits, top-10)
    """
    hits_by_porf: dict[str, list[GeneDiamondHit]] = {}

    if not diamond_output.exists():
        logger.warning(f"Diamond output not found: {diamond_output}")
        return hits_by_porf

    with open(diamond_output) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue

            query, target, evalue, bits, pident, qcov = parts[:6]
            try:
                hit = GeneDiamondHit(
                    porf_id=query,
                    target=target,
                    evalue=float(evalue),
                    bits=float(bits),
                    pident=float(pident),
                    qcov=float(qcov),
                )
                hits_by_porf.setdefault(query, []).append(hit)
            except ValueError:
                continue

    # Sort by bits score and keep top-10
    for porf_id, hits in hits_by_porf.items():
        hits.sort(key=lambda h: h.bits, reverse=True)
        hits_by_porf[porf_id] = hits[:10]
        # Assign ranks
        for i, hit in enumerate(hits_by_porf[porf_id], start=1):
            hit.rank = i

    return hits_by_porf


def extract_porf_info(porf_id: str, description: Optional[str] = None) -> tuple[int, int]:
    """
    Extract start/end coordinates from pORF ID.

    Expected format: scaffold_start_end_strand_frame

    Returns:
        Tuple of (start, end)
    """
    if description:
        parsed = parse_prodigal_header(description, porf_id)
        if parsed:
            _scaffold, start, end, _strand = parsed
            return start, end
    return 0, 0


def classify_gene_taxonomy(
    porf_id: str,
    hits: list[GeneDiamondHit],
    high_pident_euk_threshold: float = 70.0,
    start: Optional[int] = None,
    end: Optional[int] = None,
) -> GeneTaxonomy:
    """
    Classify gene taxonomy based on top-10 Diamond hits.

    Args:
        porf_id: pORF identifier
        hits: List of top-10 hits (sorted by bits)
        high_pident_euk_threshold: % identity threshold for flagging high-confidence EUK genes

    Returns:
        GeneTaxonomy object with classification
    """
    if start is None or end is None:
        start, end = extract_porf_info(porf_id)

    if not hits:
        return GeneTaxonomy(
            porf_id=porf_id,
            porf_start=start,
            porf_end=end,
        )

    # Extract prefixes from top-10
    top10_prefixes = [extract_prefix(hit.target) for hit in hits]
    top10_raw_prefixes = [extract_raw_prefix(hit.target) for hit in hits]
    top10_targets = [hit.target for hit in hits]
    top10_bitscores = [hit.bits for hit in hits]
    top10_pidents = [hit.pident for hit in hits]
    top10_evalues = [hit.evalue for hit in hits]
    top1_prefix = top10_prefixes[0] if top10_prefixes else "UNKNOWN"
    top1_target = hits[0].target if hits else "."
    top1_pident = hits[0].pident if hits else 0.0
    top1_evalue = hits[0].evalue if hits else 1.0

    # Check for viral families in top-10
    has_vp = "VP" in top10_prefixes
    has_plv = "PLV" in top10_prefixes
    has_vp_plv = has_vp or has_plv or ("PPV" in top10_prefixes)
    has_ncldv = "NCLDV" in top10_prefixes
    has_mirus = "MIRUS" in top10_prefixes
    has_viral = any(
        prefix in VIRAL_PREFIXES
        and idx < len(top10_pidents)
        and top10_pidents[idx] >= MIN_VIRAL_HIT_PIDENT
        for idx, prefix in enumerate(top10_prefixes)
    )
    has_ncldv_mirus = has_ncldv or has_mirus

    # Check for high-identity EUK gene
    is_high_pident_euk = (
        top1_prefix == "EUK" and top1_pident >= high_pident_euk_threshold
    )

    return GeneTaxonomy(
        porf_id=porf_id,
        porf_start=start,
        porf_end=end,
        top1_target=top1_target,
        top1_prefix=top1_prefix,
        top1_pident=top1_pident,
        top1_evalue=top1_evalue,
        top10_prefixes=top10_prefixes,
        top10_raw_prefixes=top10_raw_prefixes,
        top10_targets=top10_targets,
        top10_bitscores=top10_bitscores,
        top10_pidents=top10_pidents,
        top10_evalues=top10_evalues,
        has_ncldv_mirus=has_ncldv_mirus,
        has_vp_plv=has_vp_plv,
        has_viral=has_viral,
        is_high_pident_euk=is_high_pident_euk,
    )


def run_gene_taxonomy_diamond(
    eve_id: str,
    scaffold: str,
    start: int,
    end: int,
    proteome_fasta: Path,
    combined_faa_db: Path,
    output_dir: Path,
    threads: int = 4,
    high_pident_euk_threshold: float = 70.0,
    search_backend: str = "diamond",
) -> tuple[list[GeneTaxonomy], dict]:
    """
    Run Diamond on all genes within an EVE region for taxonomy assignment.

    This is the main function implementing Step 9 from PIPELINE_HMM_GATED_PLAN.md.

    Args:
        eve_id: EVE identifier (e.g., "EVE_scaffold_12345")
        scaffold: Scaffold name
        start: Region start (bp)
        end: Region end (bp)
        proteome_fasta: Path to proteome FASTA (all pORFs)
        combined_faa_db: Path to combined Diamond database (cellular + gvmags + phage)
        output_dir: Output directory for per-EVE files
        threads: Number of threads for Diamond
        high_pident_euk_threshold: % identity threshold for flagging EUK genes

    Returns:
        Tuple of (list of GeneTaxonomy objects, summary dict)
    """
    filename_component = safe_filename_component(eve_id)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Extract all pORFs within region boundaries
    logger.info(f"{eve_id}: Extracting genes from {scaffold}:{start}-{end}")

    region_porfs = []
    region_coords: dict[str, tuple[int, int]] = {}
    tmp_faa = tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False)
    tmp_faa_path = Path(tmp_faa.name)
    try:
        with tmp_faa:
            for record in SeqIO.parse(proteome_fasta, "fasta"):
                parsed = parse_prodigal_header(record.description, record.id)
                if not parsed:
                    continue
                gene_scaffold, gene_start, gene_end, _strand = parsed
                if gene_scaffold != scaffold:
                    continue
                if gene_start < end and gene_end > start:
                    tmp_faa.write(f">{record.id}\n{record.seq}\n")
                    region_porfs.append(record.id)
                    region_coords[record.id] = (gene_start, gene_end)

        if not region_porfs:
            logger.warning(f"{eve_id}: No pORFs found in region")
            return [], {
                "total": 0,
                "ncldv_mirus": 0,
                "vp_plv": 0,
                "viral_top10": 0,
                "high_pident_euk": 0,
                "has_ncldv_mirus": False,
                "has_vp_plv": False,
                "dominant_family": "UNKNOWN",
                "dominant_fraction": 0.0,
            }

        logger.info(f"{eve_id}: Extracted {len(region_porfs)} pORFs")

        # Step 2: Run Diamond against combined database
        diamond_output = output_dir / f"{filename_component}_diamond.tsv"
        require_strict_child(output_dir, diamond_output)
        try:
            logger.info(
                "%s: Running Diamond against %s (%d pORFs, threads=%d)",
                eve_id,
                combined_faa_db,
                len(region_porfs),
                threads,
            )
            run_diamond_blastp(
                query_fasta=tmp_faa_path,
                target_db=combined_faa_db,
                output_file=diamond_output,
                threads=threads,
                fast=False,
                search_backend=search_backend,
            )
        except subprocess.CalledProcessError as e:
            # Do NOT convert a Diamond failure into a zero-viral summary: a fabricated
            # all-host result is indistinguishable from a genuine host region and would
            # silently produce false negatives. Fail loud so the genome is marked failed.
            logger.error("%s: Diamond search failed; aborting gene taxonomy: %s", eve_id, e)
            raise
    finally:
        tmp_faa_path.unlink(missing_ok=True)

    # Step 3: Parse top-10 hits
    hits_by_porf = parse_diamond_top10(diamond_output)

    # Step 4: Classify each gene
    gene_taxonomies = []
    for porf_id in region_porfs:
        hits = hits_by_porf.get(porf_id, [])
        start_end = region_coords.get(porf_id, (0, 0))
        taxonomy = classify_gene_taxonomy(
            porf_id,
            hits,
            high_pident_euk_threshold,
            start=start_end[0],
            end=start_end[1],
        )
        gene_taxonomies.append(taxonomy)

    # Step 5: Write per-EVE output file
    output_tsv = output_dir / f"{filename_component}.tsv"
    require_strict_child(output_dir, output_tsv)
    with open(output_tsv, "w") as f:
        # Write header
        f.write(
            "porf_id\tporf_start\tporf_end\ttop1_target\ttop1_prefix\t"
            "top1_pident\ttop10_prefixes\ttop10_raw_prefixes\thas_ncldv_mirus\thas_vp_plv\t"
            "has_viral\tis_high_pident_euk\n"
        )
        # Write data
        for tax in gene_taxonomies:
            f.write(
                f"{tax.porf_id}\t{tax.porf_start}\t{tax.porf_end}\t"
                f"{tax.top1_target}\t{tax.top1_prefix}\t{tax.top1_pident:.2f}\t"
                f"{','.join(tax.top10_prefixes)}\t{','.join(tax.top10_raw_prefixes)}\t"
                f"{int(tax.has_ncldv_mirus)}\t{int(tax.has_vp_plv)}\t"
                f"{int(tax.has_viral)}\t{int(tax.is_high_pident_euk)}\n"
            )

    # Step 6: Generate summary statistics
    ncldv_mirus = sum(1 for t in gene_taxonomies if t.has_ncldv_mirus)
    vp_plv = sum(1 for t in gene_taxonomies if t.has_vp_plv)
    viral_top10 = sum(1 for t in gene_taxonomies if t.has_viral)
    dominant_family, dominant_fraction = summarize_dominant_family(gene_taxonomies)

    summary = {
        "total": len(gene_taxonomies),
        "ncldv_mirus": ncldv_mirus,
        "vp_plv": vp_plv,
        "viral_top10": viral_top10,
        "high_pident_euk": sum(1 for t in gene_taxonomies if t.is_high_pident_euk),
        "has_ncldv_mirus": any(t.has_ncldv_mirus for t in gene_taxonomies),
        "has_vp_plv": vp_plv > 0,
        "dominant_family": dominant_family,
        "dominant_fraction": dominant_fraction,
    }

    logger.info(
        f"{eve_id}: Gene taxonomy complete - "
        f"{summary['total']} genes, "
        f"{summary['ncldv_mirus']} NCLDV/MIRUS, "
        f"{summary['vp_plv']} VP/PLV, "
        f"{summary['high_pident_euk']} high-identity EUK"
    )

    return gene_taxonomies, summary


def materialize_gene_taxonomy_batch_from_cached_hits(
    regions: list[dict],
    proteome_fasta: Path,
    diamond_hits: dict[str, list],
    output_dir: Path,
    high_pident_euk_threshold: float = 70.0,
) -> dict[str, tuple[list[GeneTaxonomy], dict]]:
    """Rebuild Phase 2a region records from raw-pORF superset hits.

    Overlapping regions intentionally receive separate taxonomy records for the
    same pORF, matching the legacy region-prefixed query behavior without a
    second DIAMOND search.
    """

    filename_components = safe_filename_components(
        (region["eve_id"] for region in regions),
        label="EVE ID",
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    regions_by_scaffold: dict[str, list[dict]] = {}
    for region in regions:
        regions_by_scaffold.setdefault(region["scaffold"], []).append(region)
    for scaffold_regions in regions_by_scaffold.values():
        scaffold_regions.sort(key=lambda region: region["start"])

    results_with_hits: dict[str, list[GeneTaxonomy]] = {
        region["eve_id"]: [] for region in regions
    }
    results_without_hits: dict[str, list[GeneTaxonomy]] = {
        region["eve_id"]: [] for region in regions
    }
    for record in SeqIO.parse(proteome_fasta, "fasta"):
        parsed = parse_prodigal_header(record.description, record.id)
        if not parsed:
            continue
        scaffold, start, end, _strand = parsed
        for region in regions_by_scaffold.get(scaffold, []):
            if start < region["end"] and end > region["start"]:
                hits = diamond_hits.get(record.id, [])[:10]
                taxonomy = classify_gene_taxonomy(
                    record.id,
                    hits,
                    high_pident_euk_threshold,
                    start=start,
                    end=end,
                )
                destination = results_with_hits if hits else results_without_hits
                destination[region["eve_id"]].append(taxonomy)

    # Match the legacy two-pass layout: hit-bearing records first, then no-hit
    # genes. The legacy hit group follows DIAMOND output order, while this
    # cached path follows proteome order; downstream logic keys records by ID.
    results = {
        eve_id: results_with_hits[eve_id] + results_without_hits[eve_id]
        for eve_id in results_with_hits
    }

    summaries: dict[str, dict] = {}
    for eve_id, taxonomies in results.items():
        dominant_family, dominant_fraction = summarize_dominant_family(taxonomies)

        ncldv_mirus = sum(t.has_ncldv_mirus for t in taxonomies)
        vp_plv = sum(t.has_vp_plv for t in taxonomies)
        summaries[eve_id] = {
            "total": len(taxonomies),
            "ncldv_mirus": ncldv_mirus,
            "vp_plv": vp_plv,
            "viral_top10": sum(t.has_viral for t in taxonomies),
            "high_pident_euk": sum(t.is_high_pident_euk for t in taxonomies),
            "has_ncldv_mirus": ncldv_mirus > 0,
            "has_vp_plv": vp_plv > 0,
            "dominant_family": dominant_family,
            "dominant_fraction": dominant_fraction,
        }

        try:
            output_tsv = output_dir / f"{filename_components[eve_id]}.tsv"
            require_strict_child(output_dir, output_tsv)
            with output_tsv.open("w") as handle:
                handle.write(
                    "porf_id\tstart\tend\ttop1_prefix\ttop1_target\t"
                    "top1_pident\ttop10_prefixes\ttop10_raw_prefixes\t"
                    "has_ncldv_mirus\thas_vp_plv\thas_viral\t"
                    "is_high_pident_euk\n"
                )
                for taxonomy in taxonomies:
                    handle.write(
                        f"{taxonomy.porf_id}\t{taxonomy.porf_start}\t"
                        f"{taxonomy.porf_end}\t{taxonomy.top1_prefix}\t"
                        f"{taxonomy.top1_target}\t{taxonomy.top1_pident:.1f}\t"
                        f"{','.join(taxonomy.top10_prefixes or [])}\t"
                        f"{','.join(taxonomy.top10_raw_prefixes or [])}\t"
                        f"{int(taxonomy.has_ncldv_mirus)}\t"
                        f"{int(taxonomy.has_vp_plv)}\t"
                        f"{int(taxonomy.has_viral)}\t"
                        f"{int(taxonomy.is_high_pident_euk)}\n"
                    )
        except Exception as exc:
            logger.warning(
                "Gene taxonomy superset: failed to write TSV for %s: %s",
                eve_id,
                exc,
            )

    logger.info(
        "Gene taxonomy superset: materialized %d records across %d regions",
        sum(len(taxonomies) for taxonomies in results.values()),
        len(results),
    )
    return {
        eve_id: (results[eve_id], summaries[eve_id])
        for eve_id in results
    }


def run_gene_taxonomy_diamond_batch(
    regions: list[dict],
    proteome_fasta: Path,
    combined_faa_db: Path,
    output_dir: Path,
    threads: int = 4,
    high_pident_euk_threshold: float = 70.0,
    search_backend: str = "diamond",
) -> dict[str, tuple[list[GeneTaxonomy], dict]]:
    """
    Run Diamond gene taxonomy for all candidate regions in one batch.

    Uses genome-wide prodigal-gv genes from Phase 0 (no re-running prodigal).

    Args:
        regions: List of dicts with keys: eve_id, scaffold, start, end
        proteome_fasta: Path to prodigal-gv proteins FASTA (genome-wide genes)
        combined_faa_db: Path to combined Diamond database
        output_dir: Output directory for batch outputs
        threads: Threads for Diamond
        high_pident_euk_threshold: % identity threshold for EUK

    Returns:
        Dict mapping eve_id -> (gene_taxonomies, summary)
    """
    filename_components = safe_filename_components(
        (region["eve_id"] for region in regions),
        label="EVE ID",
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    regions_by_scaffold: dict[str, list[dict]] = {}
    for region in regions:
        regions_by_scaffold.setdefault(region["scaffold"], []).append(region)

    for scaffold in regions_by_scaffold:
        regions_by_scaffold[scaffold].sort(key=lambda r: r["start"])

    region_genes: dict[str, list[str]] = {r["eve_id"]: [] for r in regions}
    query_map: dict[str, tuple[str, str]] = {}
    query_coords: dict[str, tuple[int, int]] = {}
    total_queries = 0

    tmp_faa = tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False)
    tmp_faa_path = Path(tmp_faa.name)
    try:
        with tmp_faa:
            for record in SeqIO.parse(proteome_fasta, "fasta"):
                parsed = parse_prodigal_header(record.description, record.id)
                if not parsed:
                    continue
                scaffold, start, end, _strand = parsed
                for region in regions_by_scaffold.get(scaffold, []):
                    if start < region["end"] and end > region["start"]:
                        query_id = f"{region['eve_id']}|{record.id}"
                        tmp_faa.write(f">{query_id}\n{record.seq}\n")
                        query_map[query_id] = (region["eve_id"], record.id)
                        query_coords[query_id] = (start, end)
                        region_genes.setdefault(region["eve_id"], []).append(record.id)
                        total_queries += 1

        total_regions = len(regions)
        total_region_genes = sum(len(v) for v in region_genes.values())
        if total_queries == 0:
            logger.warning("Gene taxonomy batch: no prodigal genes found in any region")
            return {
                region["eve_id"]: ([], {
                    "total": 0,
                "ncldv_mirus": 0,
                "vp_plv": 0,
                "viral_top10": 0,
                "high_pident_euk": 0,
                "has_ncldv_mirus": False,
                "has_vp_plv": False,
                "dominant_family": "UNKNOWN",
                "dominant_fraction": 0.0,
            })
                for region in regions
            }

        logger.info(
            "Gene taxonomy batch: regions=%d total_region_genes=%d queries=%d",
            total_regions,
            total_region_genes,
            total_queries,
        )
        logger.info(
            "Gene taxonomy batch: running Diamond vs %s (threads=%d)",
            combined_faa_db,
            threads,
        )

        diamond_output = output_dir / "batch_diamond.tsv"
        try:
            run_diamond_blastp(
                query_fasta=tmp_faa_path,
                target_db=combined_faa_db,
                output_file=diamond_output,
                threads=threads,
                fast=False,
                search_backend=search_backend,
            )
        except subprocess.CalledProcessError as e:
            # Fail loud rather than fabricating zero-viral summaries for every region:
            # a swallowed Diamond failure would silently over-trim/false-negative EVEs.
            logger.error("Gene taxonomy batch: Diamond search failed; aborting: %s", e)
            raise
    finally:
        tmp_faa_path.unlink(missing_ok=True)

    hits_by_query = parse_diamond_top10(diamond_output)
    logger.info(
        "Gene taxonomy batch: parsed %d queries with hits from %s",
        len(hits_by_query),
        diamond_output,
    )

    results: dict[str, list[GeneTaxonomy]] = {r["eve_id"]: [] for r in regions}
    seen_genes: dict[str, set[str]] = {r["eve_id"]: set() for r in regions}
    matched_queries = 0
    unmatched_queries = 0
    for query_id, hits in hits_by_query.items():
        mapping = query_map.get(query_id)
        if not mapping:
            unmatched_queries += 1
            if unmatched_queries <= 3:
                logger.warning("Gene taxonomy batch: unmatched query_id=%s", query_id)
            continue
        matched_queries += 1
        eve_id, porf_id = mapping
        start, end = query_coords.get(query_id, (0, 0))
        if not hits:
            taxonomy = GeneTaxonomy(
                porf_id=porf_id,
                porf_start=start,
                porf_end=end,
            )
        else:
            top10_prefixes = [extract_prefix(hit.target) for hit in hits]
            top10_raw_prefixes = [extract_raw_prefix(hit.target) for hit in hits]
            top10_targets = [hit.target for hit in hits]
            top10_bitscores = [hit.bits for hit in hits]
            top10_pidents = [hit.pident for hit in hits]
            top10_evalues = [hit.evalue for hit in hits]
            top1_prefix = top10_prefixes[0] if top10_prefixes else "UNKNOWN"
            top1_target = hits[0].target
            top1_pident = hits[0].pident
            has_vp = "VP" in top10_prefixes
            has_plv = "PLV" in top10_prefixes
            has_vp_plv = has_vp or has_plv or ("PPV" in top10_prefixes)
            has_viral = any(
                prefix in VIRAL_PREFIXES
                and idx < len(top10_pidents)
                and top10_pidents[idx] >= MIN_VIRAL_HIT_PIDENT
                for idx, prefix in enumerate(top10_prefixes)
            )
            has_ncldv_mirus = "NCLDV" in top10_prefixes or "MIRUS" in top10_prefixes
            is_high_pident_euk = (
                top1_prefix == "EUK" and top1_pident >= high_pident_euk_threshold
            )
            taxonomy = GeneTaxonomy(
                porf_id=porf_id,
                porf_start=start,
                porf_end=end,
                top1_target=top1_target,
                top1_prefix=top1_prefix,
                top1_pident=top1_pident,
                top1_evalue=hits[0].evalue,
                top10_prefixes=top10_prefixes,
                top10_raw_prefixes=top10_raw_prefixes,
                top10_targets=top10_targets,
                top10_bitscores=top10_bitscores,
                top10_pidents=top10_pidents,
                top10_evalues=top10_evalues,
                has_ncldv_mirus=has_ncldv_mirus,
                has_vp_plv=has_vp_plv,
                has_viral=has_viral,
                is_high_pident_euk=is_high_pident_euk,
            )
        results[eve_id].append(taxonomy)
        seen_genes[eve_id].add(porf_id)

    logger.info(
        "Gene taxonomy batch: matched=%d unmatched=%d (query_map has %d entries)",
        matched_queries,
        unmatched_queries,
        len(query_map),
    )

    # Add genes with no Diamond hits
    for eve_id, genes in region_genes.items():
        missing = [g for g in genes if g not in seen_genes.get(eve_id, set())]
        for porf_id in missing:
            start, end = 0, 0
            for query_id, mapping in query_map.items():
                mapped_eve, mapped_id = mapping
                if mapped_eve == eve_id and mapped_id == porf_id:
                    start, end = query_coords.get(query_id, (0, 0))
                    break
            results[eve_id].append(
                GeneTaxonomy(
                    porf_id=porf_id,
                    porf_start=start,
                    porf_end=end,
                )
            )

    # Log aggregation stats
    total_records = sum(len(r) for r in results.values())
    logger.info(
        "Gene taxonomy batch: aggregated %d records across %d regions",
        total_records,
        len(results),
    )

    summaries: dict[str, dict] = {}
    for eve_id, taxonomies in results.items():
        ncldv_mirus = sum(1 for t in taxonomies if t.has_ncldv_mirus)
        vp_plv = sum(1 for t in taxonomies if t.has_vp_plv)
        viral_top10 = sum(1 for t in taxonomies if t.has_viral)
        dominant_family, dominant_fraction = summarize_dominant_family(taxonomies)
        high_pident_euk = sum(1 for t in taxonomies if t.is_high_pident_euk)
        summaries[eve_id] = {
            "total": len(taxonomies),
            "ncldv_mirus": ncldv_mirus,
            "vp_plv": vp_plv,
            "viral_top10": viral_top10,
            "high_pident_euk": high_pident_euk,
            "has_ncldv_mirus": ncldv_mirus > 0,
            "has_vp_plv": vp_plv > 0,
            "dominant_family": dominant_family,
            "dominant_fraction": dominant_fraction,
        }

        # Write per-EVE TSV (non-fatal if fails)
        try:
            output_tsv = output_dir / f"{filename_components[eve_id]}.tsv"
            require_strict_child(output_dir, output_tsv)
            with open(output_tsv, "w") as f:
                f.write(
                    "porf_id\tstart\tend\ttop1_prefix\ttop1_target\t"
                    "top1_pident\ttop10_prefixes\ttop10_raw_prefixes\thas_ncldv_mirus\thas_vp_plv\t"
                    "has_viral\tis_high_pident_euk\n"
                )
                for t in taxonomies:
                    f.write(
                        f"{t.porf_id}\t{t.porf_start}\t{t.porf_end}\t{t.top1_prefix}\t"
                        f"{t.top1_target}\t{t.top1_pident:.1f}\t"
                        f"{','.join(t.top10_prefixes or [])}\t"
                        f"{','.join(t.top10_raw_prefixes or [])}\t{int(t.has_ncldv_mirus)}\t"
                        f"{int(t.has_vp_plv)}\t{int(t.has_viral)}\t"
                        f"{int(t.is_high_pident_euk)}\n"
                    )
        except Exception as e:
            logger.warning("Gene taxonomy batch: failed to write TSV for %s: %s", eve_id, e)

    logger.info("Gene taxonomy batch: completed with %d per-region TSVs", len(results))

    return {eve_id: (results[eve_id], summaries[eve_id]) for eve_id in results}
