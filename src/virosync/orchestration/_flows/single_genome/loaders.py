"""On-disk artifact loaders and parsers for the single-genome flow.

Leaf module: pure parse/load helpers used by the phase subflows (resume paths)
and the orchestrator. Depends only on leaf pipeline types, never on the phase
modules or orchestrator, so it introduces no import cycle.
"""

import csv
import json
import logging
from pathlib import Path

from virosync.utils.atomic_write import atomic_write_context
from virosync.pipeline.phase1.hhg_seeding import Anchor
from virosync.pipeline.phase3.mcp_detection import is_mcp_gene
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy, ControlStats

# Module-level logger for parse helpers.
logger = logging.getLogger(__name__)


def _count_fasta_records(fasta_path: Path) -> int:
    """Count FASTA entries without loading sequences into memory."""
    count = 0
    with fasta_path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def _load_refined_boundaries_from_bed(boundaries_bed: Path) -> list[RefinedBoundary]:
    """Load refined boundaries from BED."""
    boundaries: list[RefinedBoundary] = []
    with boundaries_bed.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            scaffold = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3] if len(parts) > 3 else ""
            score = int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else 0
            confidence = score / 1000.0
            seed_id = name if name.startswith("seed_") else ""
            boundaries.append(
                RefinedBoundary(
                    scaffold=scaffold,
                    start=start,
                    end=end,
                    seed_id=seed_id,
                    original_start=start,
                    original_end=end,
                    confidence=confidence,
                    posterior_probability=confidence,
                )
            )
    return boundaries


def _load_taxonomy_map_from_tsv(taxonomy_tsv: Path) -> dict[str, GeneTaxonomy]:
    """Load Phase 2b taxonomy TSV written by write_taxonomy_map()."""
    taxonomy_map: dict[str, GeneTaxonomy] = {}
    with taxonomy_tsv.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if "top10_targets" not in fieldnames:
            logger.warning(
                "Taxonomy TSV %s uses old format (missing top10 hit columns). "
                "Host signature scoring will be degraded. "
                "Delete phase2/ and re-run to regenerate.",
                taxonomy_tsv,
            )
        for row in reader:
            porf_id = (row.get("porf_id") or "").strip()
            if not porf_id:
                continue
            top10_prefixes = [p for p in (row.get("top10_prefixes") or "").split(",") if p]
            top10_targets = [t for t in (row.get("top10_targets") or "").split(",") if t]
            top10_bits_raw = (row.get("top10_bits") or "").split(",")
            top10_bits = [float(b) for b in top10_bits_raw if b]
            top10_pidents_raw = (row.get("top10_pidents") or "").split(",")
            top10_pidents = [float(p) for p in top10_pidents_raw if p]
            top10_evalues_raw = (row.get("top10_evalues") or "").split(",")
            top10_evalues = [float(e) for e in top10_evalues_raw if e]
            taxonomy_map[porf_id] = GeneTaxonomy(
                porf_id=porf_id,
                scaffold=(row.get("scaffold") or "").strip(),
                start=int(row.get("start") or 0),
                end=int(row.get("end") or 0),
                top1_target=(row.get("top1_target") or "").strip(),
                top1_prefix=(row.get("top1_prefix") or "UNKNOWN").strip() or "UNKNOWN",
                top1_pident=float(row.get("top1_pident") or 0.0),
                top1_evalue=float(row.get("top1_evalue") or 1.0),
                top10_prefixes=top10_prefixes,
                top10_targets=top10_targets,
                top10_bits=top10_bits,
                top10_pidents=top10_pidents,
                top10_evalues=top10_evalues,
                has_ncldv_mirus=bool(int(row.get("has_ncldv_mirus") or 0)),
                has_vp_plv=bool(int(row.get("has_vp_plv") or 0)),
                has_viral=bool(int(row.get("has_viral") or 0)),
                has_hit=bool(int(row.get("has_hit") or 0)),
            )
    return taxonomy_map


def _load_control_stats_from_json(control_stats_json: Path) -> ControlStats:
    """Load Phase 2b control stats JSON written by write_control_stats()."""
    with control_stats_json.open() as handle:
        payload = json.load(handle)
    return ControlStats(
        n_genes=int(payload.get("n_genes", 0)),
        n_no_hits=int(payload.get("n_no_hits", 0)),
        no_hit_frequency=float(payload.get("no_hit_frequency", 0.0)),
        host_frequency=float(payload.get("host_frequency", 0.0)),
        mean_pident=float(payload.get("mean_pident", 0.0)),
        dominant_organism=str(payload.get("dominant_organism", "unknown")),
        host_prefix=str(payload.get("host_prefix", "EUK__")),
    )


def _load_marker_seed_regions(seed_regions_tsv: Path) -> list[dict[str, int | str]]:
    """Load marker seed regions TSV written by region_assembly."""
    regions: list[dict[str, int | str]] = []
    with seed_regions_tsv.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            scaffold = (row.get("scaffold") or "").strip()
            if not scaffold:
                continue
            regions.append(
                {
                    "region_id": (row.get("region_id") or "").strip(),
                    "scaffold": scaffold,
                    "start": int(row.get("start") or 0),
                    "end": int(row.get("end") or 0),
                }
            )
    regions.sort(key=lambda r: (str(r["scaffold"]), int(r["start"]), int(r["end"])))
    return regions


def _build_merged_seeds_from_regions(
    regions: list[dict[str, int | str]],
    validated_markers: list,
) -> list:
    """Rebuild MergedSeed list from saved Phase 1 regions and validated markers."""
    from virosync.pipeline.phase1.seed_merger import MergedSeed

    markers_by_scaffold: dict[str, list] = {}
    for marker in validated_markers:
        markers_by_scaffold.setdefault(marker.scaffold, []).append(marker)

    merged_seeds = []
    for idx, region in enumerate(regions):
        scaffold = str(region["scaffold"])
        start = int(region["start"])
        end = int(region["end"])
        anchors: list[Anchor] = []
        for marker in markers_by_scaffold.get(scaffold, []):
            if marker.start >= end or marker.end <= start:
                continue
            porf_id = getattr(marker, "porf_id", None) or getattr(marker, "query_porf", "")
            anchors.append(
                Anchor(
                    porf_id=porf_id,
                    scaffold=marker.scaffold,
                    start=marker.start,
                    end=marker.end,
                    strand=getattr(marker, "strand", "+"),
                    hallmark_gene=marker.hmm_target,
                    score=marker.hmm_score,
                    evalue=getattr(marker, "hmm_evalue", 0.0),
                )
            )

        seed = MergedSeed(
            scaffold=scaffold,
            start=start,
            end=end,
            seed_id=f"seed_{idx}_{scaffold}_{start}",
            sources=["hhg", "marker_validation"],
            confidence="high"
            if any(is_mcp_gene(a.hallmark_gene) for a in anchors)
            else "medium",
            hhg_score=len(anchors) * 10.0,
            novelty_score=0.0,
            compositional_score=0.0,
            mean_kfd=0.0,
            max_kfd=0.0,
            mean_composite=0.0,
            max_composite=0.0,
            gc_deviation=0.0,
            n_windows=0,
            cluster_ids=[],
            anchors=anchors,
            hhg_anchors=anchors,
        )
        merged_seeds.append(seed)
    return merged_seeds


def _load_interproscan_summary(summary_path: Path) -> dict[str, dict]:
    """Load InterProScan summary TSV into eve_id -> summary dict."""
    summaries: dict[str, dict] = {}
    with summary_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            eve_id = (row.get("eve_id") or "").strip()
            if not eve_id:
                continue
            summaries[eve_id] = {
                "total_hits": int(row.get("total_hits") or 0),
                "viral_hits": int(row.get("viral_hits") or 0),
                "keyword_hits": []
                if (row.get("keyword_hits") or ".") == "."
                else [k for k in (row.get("keyword_hits") or "").split("|") if k],
                "category_hits": []
                if (row.get("category_hits") or ".") == "."
                else [k for k in (row.get("category_hits") or "").split("|") if k],
                "family_hits": []
                if (row.get("family_hits") or ".") == "."
                else [k for k in (row.get("family_hits") or "").split("|") if k],
                "category_score": float(row.get("category_score") or 0.0),
                "numt_hits": int(row.get("numt_hits") or 0),
                "numt_markers": []
                if (row.get("numt_markers") or ".") == "."
                else [k for k in (row.get("numt_markers") or "").split("|") if k],
            }
    return summaries


def _serialize_tmvec_cache(precomputed_tmvec: dict[str, dict], output_path: Path) -> None:
    """Persist precomputed TMVec hits for resume."""
    payload: dict[str, dict[str, dict | None]] = {}
    for porf_id, db_hits in precomputed_tmvec.items():
        db_payload: dict[str, dict | None] = {}
        for db_name, hit in db_hits.items():
            if hit is None:
                db_payload[db_name] = None
                continue
            db_payload[db_name] = {
                "target_id": hit.target_id,
                "tm_score": hit.tm_score,
                "database": hit.database,
                "protein_name": hit.protein_name,
                "organism": hit.organism,
                "lineage": hit.lineage,
                "keywords": hit.keywords,
            }
        payload[porf_id] = db_payload
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with atomic_write_context(output_path, "w") as handle:
        json.dump(payload, handle, indent=2)


def _load_tmvec_cache(cache_path: Path) -> dict[str, dict]:
    """Load persisted TMVec cache from JSON."""
    from virosync.pipeline.phase3.tmvec_database import TMVecHit

    with cache_path.open() as handle:
        payload = json.load(handle)
    precomputed_tmvec: dict[str, dict] = {}
    for porf_id, db_hits in payload.items():
        hit_map: dict[str, object | None] = {}
        for db_name, hit in (db_hits or {}).items():
            if not hit:
                hit_map[db_name] = None
                continue
            hit_map[db_name] = TMVecHit(
                target_id=str(hit.get("target_id", "")),
                tm_score=float(hit.get("tm_score", 0.0)),
                database=str(hit.get("database", db_name)),
                protein_name=hit.get("protein_name"),
                organism=hit.get("organism"),
                lineage=hit.get("lineage"),
                keywords=hit.get("keywords"),
            )
        precomputed_tmvec[str(porf_id)] = hit_map
    return precomputed_tmvec


def _safe_int(value: object, default: int = 0) -> int:
    """Parse integer from scalar-like value with fallback."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return default
