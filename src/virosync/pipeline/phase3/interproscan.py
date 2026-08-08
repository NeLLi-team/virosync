"""
InterProScan annotation for EVE candidate regions.

Runs InterProScan on prodigal-gv proteins overlapping candidate regions,
then summarizes viral/hallmark-related annotations to boost confidence.
"""

from __future__ import annotations

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable, Optional

from Bio import SeqIO

from virosync.pipeline.phase0.prodigal import parse_prodigal_header

logger = logging.getLogger(__name__)


DEFAULT_VIRAL_KEYWORDS = [
    "viral",
    "virion",
    "virophage",
    "capsid",
    "major capsid",
    "mcp",
    "portal",
    "terminase",
    "packaging",
    "ncldv",
    "mirus",
    "pox",
    "polymerase",
    "helicase",
    "transcription factor",
    "vltf",
]

NUMT_KEYWORDS = [
    "mitochondrial",
    "mitochondrion",
    "cytochrome",
    "cox",
    "cytochrome c oxidase",
    "nadh dehydrogenase",
    "complex i",
    "nadh-ubiquinone",
    "atp synthase",
    "atpase f0",
    "atpase f1",
    "electron transport",
    "oxidative phosphorylation",
]

INTERPRO_CATEGORY_KEYWORDS = {
    "capsid": ["capsid", "mcp", "major capsid", "minor capsid", "coat"],
    "portal": ["portal"],
    "terminase": ["terminase"],
    "packaging_atpase": ["packaging", "a32", "atpase"],
    "penton": ["penton"],
    "triplex": ["triplex"],
    "protease": ["protease", "peptidase"],
    "polymerase": ["polymerase", "polb", "dna polymerase"],
    "helicase": ["helicase", "primase", "d5"],
    "transcription": ["vltf", "transcription factor", "rna polymerase", "rnap"],
    "numt": NUMT_KEYWORDS,
}

INTERPRO_FAMILY_KEYWORDS = {
    # Virophage and Polinton-like keywords both denote Preplasmiviricota.
    "PPV": ["virophage", "polinton", "polinton-like", "plv"],
    "MIRUS": ["mirus"],
    "NCLDV": ["ncldv", "pox", "nucleocytoviricota", "mimivirus"],
}




def _keyword_match(text: str, keywords: Iterable[str]) -> bool:
    text = text.lower()
    return any(k.lower() in text for k in keywords)


def _categorize_interpro(text: str) -> set[str]:
    text = text.lower()
    categories = set()
    for category, keywords in INTERPRO_CATEGORY_KEYWORDS.items():
        if any(k in text for k in keywords):
            categories.add(category)
    return categories


def _infer_family_from_interpro(text: str) -> set[str]:
    text = text.lower()
    families = set()
    for family, keywords in INTERPRO_FAMILY_KEYWORDS.items():
        if any(k in text for k in keywords):
            families.add(family)
    return families


def run_interproscan_batch(
    regions: list[dict],
    proteome_fasta: Path,
    interproscan_dir: Path,
    output_dir: Path,
    threads: int = 4,
    keywords: Optional[list[str]] = None,
    applications: Optional[list[str]] = None,
) -> dict[str, dict]:
    """
    Run InterProScan on all candidate regions in one batch.

    Returns:
        Mapping of eve_id -> summary dict with counts and keyword hits.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    keywords = keywords or DEFAULT_VIRAL_KEYWORDS

    interproscan_exec = Path(interproscan_dir) / "interproscan.sh"
    if not interproscan_exec.exists():
        raise FileNotFoundError(f"interproscan.sh not found at {interproscan_exec}")

    regions_by_scaffold: dict[str, list[dict]] = {}
    for region in regions:
        regions_by_scaffold.setdefault(region["scaffold"], []).append(region)
    for scaffold in regions_by_scaffold:
        regions_by_scaffold[scaffold].sort(key=lambda r: r["start"])

    region_genes: dict[str, list[str]] = {r["eve_id"]: [] for r in regions}
    query_map: dict[str, str] = {}
    total_queries = 0

    stripped = 0
    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False) as tmp_faa:
        tmp_faa_path = Path(tmp_faa.name)
        for record in SeqIO.parse(proteome_fasta, "fasta"):
            parsed = parse_prodigal_header(record.description, record.id)
            if not parsed:
                continue
            scaffold, start, end, _strand = parsed
            for region in regions_by_scaffold.get(scaffold, []):
                if start < region["end"] and end > region["start"]:
                    query_id = f"{region['eve_id']}|{record.id}"
                    seq = str(record.seq).upper()
                    cleaned = "".join(c for c in seq if c in "ACDEFGHIKLMNPQRSTVWY")
                    if cleaned != seq:
                        stripped += 1
                    tmp_faa.write(f">{query_id}\n{cleaned}\n")
                    query_map[query_id] = region["eve_id"]
                    region_genes.setdefault(region["eve_id"], []).append(record.id)
                    total_queries += 1

    if total_queries == 0:
        logger.warning("InterProScan batch: no prodigal genes found in any region")
        tmp_faa_path.unlink(missing_ok=True)
        return {
            region["eve_id"]: {
                "total_hits": 0,
                "viral_hits": 0,
                "keyword_hits": [],
            }
            for region in regions
        }

    output_tsv = output_dir / "interproscan_batch.tsv"
    cmd = [
        str(interproscan_exec),
        "-i", str(tmp_faa_path),
        "-f", "tsv",
        "-o", str(output_tsv),
        "--cpu", str(threads),
    ]

    if applications:
        cmd.extend(["-appl", ",".join(applications)])

    if stripped:
        logger.info("InterProScan batch: stripped invalid residues from %d queries", stripped)

    logger.info("InterProScan batch: running %s", " ".join(cmd))
    try:
        subprocess.run(
            cmd,
            check=True,
            text=True,
            timeout=3600,  # 1 hour max
            capture_output=True,
        )
    except subprocess.TimeoutExpired:
        logger.error("InterProScan timed out after 1 hour")
        tmp_faa_path.unlink(missing_ok=True)
        raise
    except subprocess.CalledProcessError as e:
        logger.error("InterProScan failed: %s", e.stderr)
        tmp_faa_path.unlink(missing_ok=True)
        raise
    finally:
        tmp_faa_path.unlink(missing_ok=True)

    summaries: dict[str, dict] = {
        region["eve_id"]: {
            "total_hits": 0,
            "viral_hits": 0,
            "keyword_hits": [],
            "category_hits": [],
            "family_hits": [],
            "category_score": 0.0,
            "numt_hits": 0,
            "numt_markers": [],
        }
        for region in regions
    }
    keyword_hits: dict[str, set[str]] = {r["eve_id"]: set() for r in regions}
    category_hits: dict[str, set[str]] = {r["eve_id"]: set() for r in regions}
    family_hits: dict[str, set[str]] = {r["eve_id"]: set() for r in regions}
    numt_hits: dict[str, set[str]] = {r["eve_id"]: set() for r in regions}

    with output_tsv.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 13:
                continue
            query_id = parts[0]
            signature = parts[4]
            signature_desc = parts[5]
            interpro_id = parts[11] if len(parts) > 11 else ""
            interpro_desc = parts[12] if len(parts) > 12 else ""

            eve_id = query_map.get(query_id)
            if not eve_id:
                continue
            summaries[eve_id]["total_hits"] += 1
            text = " ".join([signature, signature_desc, interpro_id, interpro_desc]).strip()
            if _keyword_match(text, keywords):
                summaries[eve_id]["viral_hits"] += 1
                keyword_hits[eve_id].add(signature_desc or signature or interpro_id)
            categories = _categorize_interpro(text)
            if categories:
                category_hits[eve_id].update(categories)
            families = _infer_family_from_interpro(text)
            if families:
                family_hits[eve_id].update(families)
            # Check for NUMT (mitochondrial) markers
            if _keyword_match(text, NUMT_KEYWORDS):
                summaries[eve_id]["numt_hits"] += 1
                numt_hits[eve_id].add(signature_desc or signature or interpro_id)

    for eve_id, hits in keyword_hits.items():
        summaries[eve_id]["keyword_hits"] = sorted(hits)
    for eve_id, hits in category_hits.items():
        summaries[eve_id]["category_hits"] = sorted(hits)
        if hits:
            summaries[eve_id]["category_score"] = min(len(hits) / 4.0, 1.0)
    for eve_id, hits in family_hits.items():
        summaries[eve_id]["family_hits"] = sorted(hits)
    for eve_id, hits in numt_hits.items():
        summaries[eve_id]["numt_markers"] = sorted(hits)

    summary_path = output_dir / "interproscan_summary.tsv"
    with summary_path.open("w") as handle:
        handle.write(
            "eve_id\ttotal_hits\tviral_hits\tkeyword_hits\tcategory_hits\tfamily_hits\tcategory_score\tnumt_hits\tnumt_markers\n"
        )
        for eve_id, summary in summaries.items():
            handle.write(
                f"{eve_id}\t{summary['total_hits']}\t{summary['viral_hits']}\t"
                f"{'|'.join(summary['keyword_hits']) if summary['keyword_hits'] else '.'}\t"
                f"{'|'.join(summary['category_hits']) if summary['category_hits'] else '.'}\t"
                f"{'|'.join(summary['family_hits']) if summary['family_hits'] else '.'}\t"
                f"{summary['category_score']:.4f}\t"
                f"{summary['numt_hits']}\t"
                f"{'|'.join(summary['numt_markers']) if summary['numt_markers'] else '.'}\n"
            )

    logger.info("InterProScan batch: wrote %s", summary_path)
    return summaries
