"""Competitive original-top-ten evidence for the coreless LOW acceptance check."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import TYPE_CHECKING

from virosync.output_contract import canonical_family
from virosync.pipeline.phase1.viral_markers import base_marker_gene_id
from virosync.pipeline.phase2.boundary_diamond import MIN_VIRAL_HIT_PIDENT

if TYPE_CHECKING:
    from virosync.pipeline.phase3.evidence_synthesizer import VerificationResult

VIRAL_CLASSES = frozenset({"NCLDV", "MIRUS", "PPV", "CRESS", "GVMAG", "EUK-PEVE"})
CELLULAR_CLASSES = frozenset({"EUK", "BAC", "ARC"})


@dataclass(frozen=True, slots=True)
class CompetitiveViralSupport:
    """Supporting native locus IDs (or base pORF IDs) and assessment counts."""

    supporting_gene_ids: tuple[str, ...]
    assessed_gene_count: int
    unassessed_gene_count: int

    @property
    def unavailable(self) -> bool:
        """Distinguish missing evidence from a searched gene with no retained hits."""
        return self.unassessed_gene_count > 0 or self.assessed_gene_count == 0


def assess_competitive_viral_support(result: VerificationResult) -> CompetitiveViralSupport:
    """Compare retained viral and cellular scores inside the final half-open span.

    Exact physical aliases, shared base pORF IDs and native scoring-locus links
    count once. Pool their competitors before comparison, so a strong cellular
    fragment cannot be hidden by another fragment. Viral hits require 25%
    identity; deeper hits and phage-only evidence never supply positive support.
    """
    records = []
    unassessed = 0
    for record in result.gene_taxonomy_records:
        start, end = record.get("start"), record.get("end")
        if (
            not isinstance(start, int)
            or not isinstance(end, int)
            or not record.get("porf_id")
            or not record.get("scaffold")
        ):
            unassessed += 1
            continue
        if record.get("scaffold") == result.scaffold and result.start <= start < end <= result.end:
            records.append(record)
    loci = result.scoring_locus_ids
    supporting = []
    assessed = 0
    for members in _gene_groups(records, loci):
        scores = [_competitive_scores(record) for record in members]
        if any(score is None for score in scores):
            unassessed += 1
            continue
        assessed += 1
        viral = [score[0] for score in scores if score is not None and score[0] is not None]
        cellular = [score[1] for score in scores if score is not None and score[1] is not None]
        if viral and (not cellular or max(viral) > max(cellular)):
            base_ids = {base_marker_gene_id(r["porf_id"]) for r in members}
            gene_id = min({loci[gene] for gene in base_ids if gene in loci} or base_ids)
            supporting.append(gene_id)
    return CompetitiveViralSupport(tuple(sorted(set(supporting))), assessed, unassessed)


def _gene_groups(records: list[dict], loci: dict[str, str]) -> list[list[dict]]:
    """Pool transitive gene, coordinate and scoring-locus aliases in linear space."""
    parents: dict[tuple, tuple] = {}

    def root(key: tuple) -> tuple:
        parents.setdefault(key, key)
        while parents[key] != key:
            parents[key] = parents[parents[key]]
            key = parents[key]
        return key

    for record in records:
        gene = base_marker_gene_id(record["porf_id"])
        keys = [("gene", gene), ("position", record["start"], record["end"])]
        if gene in loci:
            keys.append(("locus", loci[gene]))
        group = root(keys[0])
        for key in keys[1:]:
            parents[root(key)] = group
    groups: dict[tuple, list[dict]] = {}
    for record in records:
        group = root(("gene", base_marker_gene_id(record["porf_id"])))
        groups.setdefault(group, []).append(record)
    return list(groups.values())


def _competitive_scores(record: dict) -> tuple[float | None, float | None] | None:
    """Return best viral/cellular scores, or unavailable for malformed retained hits."""
    prefixes = record.get("top10_prefixes")
    bits = record.get("top10_bitscores")
    targets = record.get("top10_targets")
    pidents = record.get("top10_pidents")
    if not all(isinstance(values, (list, tuple)) for values in (prefixes, bits, targets, pidents)):
        return None
    if not len(prefixes) == len(bits) == len(targets) == len(pidents) or len(prefixes) > 10:
        return None
    if not prefixes:
        return (None, None) if record.get("has_hit") is False else None
    viral, cellular = [], []
    for prefix, raw_score, target, raw_identity in zip(prefixes, bits, targets, pidents, strict=True):
        if isinstance(raw_score, bool):
            return None
        try:
            score = float(raw_score)
            identity = float(raw_identity)
        except (TypeError, ValueError):
            return None
        if not isfinite(score) or score < 0 or not isfinite(identity) or not 0 <= identity <= 100:
            return None
        if not isinstance(prefix, str) or not prefix.strip():
            return None
        namespace = canonical_family(prefix.rstrip("_"))
        # Some readers fold organelle namespaces into a cellular prefix.
        if str(target).startswith(("MITO__", "PLASTID__")):
            continue
        if namespace in VIRAL_CLASSES and identity >= MIN_VIRAL_HIT_PIDENT:
            viral.append(score)
        elif namespace in CELLULAR_CLASSES:
            cellular.append(score)
    return (max(viral) if viral else None, max(cellular) if cellular else None)
