"""Bounded eligibility and ownership checks for rejected adjacent seed pairs."""

from __future__ import annotations

import json
from dataclasses import replace

from virosync.output_contract import resolve_effective_eve_class
from virosync.pipeline.hallmark_evidence import assess_hallmarks
from virosync.pipeline.host_signatures import HostSignatureModel, host_lineage_share
from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit
from virosync.pipeline.phase2.boundary_diamond import (
    GeneTaxonomy,
    GenomeDiamondQuery,
    has_identity_qualified_viral_hit,
    pORF,
)
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary, boundary_candidate_id
from virosync.pipeline.phase2.seed_extension import _merge_pair
from virosync.pipeline.phase2.terminal_repeats import find_terminal_inverted_repeats
from virosync.pipeline.phase3.evidence_synthesizer import VerificationResult


def eligible_joint_pairs(
    boundaries: list[RefinedBoundary],
    results: list[VerificationResult],
    query: GenomeDiamondQuery | None,
    proteome_index: dict[str, list[pORF]],
    taxonomy: dict[str, GeneTaxonomy],
    host_model: HostSignatureModel | None,
    host_threshold: float,
    *,
    accepted_ids: frozenset[str],
) -> list[tuple[RefinedBoundary, RefinedBoundary]]:
    """Select adjacent rejected pairs in overlapping searched contexts."""
    if query is None or host_model is None or not host_model.token_weights:
        return []
    by_id = {result.eve_id: result for result in results}
    ordered = sorted(boundaries, key=lambda b: (b.scaffold, b.start, b.end, b.seed_id))
    pairs = []
    searched = set(query.all_porf_ids)
    for left, right in zip(ordered, ordered[1:]):
        children = (left, right)
        if left.scaffold != right.scaffold or left.end > right.start:
            continue
        scored = [by_id[boundary_candidate_id(child)] for child in children]
        classes = [effective_class(result) for result in scored]
        if classes[0] not in {"NCLDV", "MIRUS"} or classes[0] != classes[1]:
            continue
        if any(result.eve_id in accepted_ids for result in scored):
            continue
        if any(result.confidence_tier != "LOW" for result in scored):
            continue
        if any(not child.confirmed_seeds or child.joint_seed_provenance for child in children):
            continue
        if any(_crosses_stop(child, left.start, right.end) for child in children):
            continue
        if ("frameshift_rescue" in left.seed_sources) != ("frameshift_rescue" in right.seed_sources):
            continue
        if _has_third_overlap(ordered, children):
            continue
        if all(_has_multi_marker_witness(result, classes[0]) for result in scored):
            continue
        mappings = [
            [query.seed_gene_mappings.get(core.seed_id) for core in child.confirmed_seeds] for child in children
        ]
        if any(mapping is None for group in mappings for mapping in group):
            continue
        envelopes = [
            (min(m.flank_start_bp for m in group if m is not None), max(m.flank_end_bp for m in group if m is not None))
            for group in mappings
        ]
        if max(span[0] for span in envelopes) >= min(span[1] for span in envelopes):
            continue
        genes = [g for g in proteome_index.get(left.scaffold, []) if g.start < right.end and g.end > left.start]
        if any(g.id not in searched or g.id not in taxonomy for g in genes):
            continue
        if _gap_has_host_barrier(children, genes, taxonomy, host_model, host_threshold):
            continue
        pairs.append(children)
    return pairs


def effective_class(result: VerificationResult) -> str:
    """Resolve the ordinary final gate class without borrowing a neighbor's label."""
    return resolve_effective_eve_class(
        confidence_tier=result.confidence_tier,
        region_classification=result.region_classification,
        classification=getattr(result, "classification", ""),
        likely_family=result.likely_family,
    )


def joint_boundary(children: tuple[RefinedBoundary, RefinedBoundary]) -> RefinedBoundary:
    """Preserve exact cores and record the rejected geometry before recomputation."""
    provenance = {
        "rule": "rejected_adjacent_r2_v1",
        "children": [
            {
                "candidate_id": boundary_candidate_id(child),
                "start": child.start,
                "end": child.end,
                "seed_ids": [core.seed_id for core in child.confirmed_seeds],
                "left_stop_reason": child.left_stop_reason,
                "right_stop_reason": child.right_stop_reason,
                "left_stop": child.left_stop,
                "right_stop": child.right_stop,
            }
            for child in children
        ],
    }
    return replace(
        _merge_pair(*children), candidate_id="", joint_seed_provenance=json.dumps(provenance, sort_keys=True)
    )


def joint_result_supported(result: VerificationResult, children: tuple[RefinedBoundary, RefinedBoundary]) -> bool:
    """Require unchanged R2 support contributed by both original seed sets."""
    support = result.supplementary_viral_evidence
    if support is None or not support.eligible or support.group != effective_class(result):
        return False
    records = {record["porf_id"]: record for record in result.gene_taxonomy_records}
    return all(
        any(
            core.start <= records[gene]["start"] < records[gene]["end"] <= core.end
            for core in child.confirmed_seeds
            for gene in support.counted_gene_ids
        )
        for child in children
    )


def has_clear_internal_boundary_scan(
    children: tuple[RefinedBoundary, RefinedBoundary],
    sequence: str | None,
    markers: list[ValidatedMarkerHit],
    extension_bp: int,
) -> bool:
    """Require complete native scans without a child-bounding internal repeat pair."""
    if sequence is None:
        return False
    left, right = children
    gap_start = max(core.end for core in left.confirmed_seeds)
    gap_end = min(core.start for core in right.confirmed_seeds)
    if len(sequence) < right.end or not set(sequence[gap_start:gap_end].upper()) <= set("ACGT"):
        return False
    scan_start, scan_end = max(0, left.start - extension_bp), min(len(sequence), right.end + extension_bp)
    spans = []
    for child in children:
        retained = [
            m
            for m in markers
            if m.scaffold == child.scaffold
            and m.validation_status in {"validated", "validated_novel"}
            and child.start <= m.start < m.end <= child.end
        ]
        if not retained:
            return False
        spans.append((min(m.start for m in retained), max(m.end for m in retained)))
    try:
        for index, (start, end) in enumerate(spans):
            scan = find_terminal_inverted_repeats(
                sequence,
                scan_start=scan_start,
                scan_end=scan_end,
                marker_start=start,
                marker_end=end,
            )
            if scan.pairing_capped:
                return False
            other_start, other_end = spans[1 - index]
            if any(
                (repeat.right_end <= other_start and left.start < repeat.right_end < right.end)
                or (repeat.left_start >= other_end and left.start < repeat.left_start < right.end)
                for repeat in scan.candidates
            ):
                return False
    except (OSError, RuntimeError, ValueError):
        return False
    return True


def _crosses_stop(child: RefinedBoundary, start: int, end: int) -> bool:
    return (child.left_stop is not None and start < child.left_stop) or (
        child.right_stop is not None and end > child.right_stop
    )


def _has_third_overlap(boundaries: list[RefinedBoundary], children: tuple[RefinedBoundary, RefinedBoundary]) -> bool:
    left, right = children
    return any(
        other is not left
        and other is not right
        and other.scaffold == left.scaffold
        and other.start < right.end
        and other.end > left.start
        for other in boundaries
    )


def _has_multi_marker_witness(result: VerificationResult, family: str) -> bool:
    panel = {"NCLDV": "Nucleocytoviricota", "MIRUS": "Mirusviricota"}[family]
    return any(
        w.panel == panel and w.rule in {"capsid_partner_5kb", "three_core_15kb"}
        for w in assess_hallmarks(result).witnesses
    )


def _gap_has_host_barrier(
    children: tuple[RefinedBoundary, RefinedBoundary],
    genes: list[pORF],
    taxonomy: dict[str, GeneTaxonomy],
    model: HostSignatureModel,
    threshold: float,
) -> bool:
    left, right = children
    run = 0
    for gene in sorted(genes, key=lambda g: (g.start, g.end)):
        if gene.end <= left.end or gene.start >= right.start:
            continue
        record = taxonomy[gene.id]
        host = host_lineage_share(record, model) >= threshold
        viral = has_identity_qualified_viral_hit(record.top10_prefixes, record.top10_pidents)
        run = run + 1 if host and not viral else 0
        if run == 3:
            return True
    return False
