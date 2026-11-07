"""Bounded native assessment of adjacent accepted or rejected seed pairs."""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass, replace
from typing import Literal

from virosync.output_contract import resolve_effective_eve_class
from virosync.pipeline.eve_type_rules import composition, may_merge
from virosync.pipeline.hallmark_evidence import HallmarkEvidence, assess_hallmarks
from virosync.pipeline.host_signatures import HostSignatureModel, host_lineage_share, host_stop_available
from virosync.pipeline.phase1.marker_roles import is_ordinary_marker_evidence
from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit
from virosync.pipeline.phase2.boundary_diamond import (
    GeneTaxonomy,
    GenomeDiamondQuery,
    pORF,
)
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary, boundary_candidate_id
from virosync.pipeline.phase2.repeat_evidence import EndAssessment, assess_repeat_evidence
from virosync.pipeline.phase2.seed_extension import _merge_pair, is_extension_host_gene
from virosync.pipeline.phase2.terminal_repeats import _is_low_complexity, find_terminal_inverted_repeats
from virosync.pipeline.phase3.evidence_synthesizer import VerificationResult

_STRICT_CORE_FAMILIES = frozenset({"PolB", "A32", "RNAPL", "TFIIB", "TopoII", "VLTF3"})
_LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class JointRepeatAssessment:
    """Record repeat uncertainty without claiming a complete negative scan."""

    status: Literal["completed", "accepted_child_ambiguity"]
    new_gap: tuple[int, int]
    inter_core_ambiguous_bp: int
    gap_sides: tuple[EndAssessment, EndAssessment]


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
    viral_marker_genes: frozenset[str] = frozenset(),
) -> list[tuple[RefinedBoundary, RefinedBoundary]]:
    """Select pairs from original boundaries without consuming earlier unions."""
    if query is None or not host_stop_available(host_model):
        return []
    by_id = {result.eve_id: result for result in results}
    ordered = sorted(boundaries, key=lambda b: (b.scaffold, b.start, b.end, b.seed_id))
    pairs = []
    for left, right in zip(ordered, ordered[1:]):
        children = (left, right)
        if left.scaffold != right.scaffold or left.end > right.start:
            continue
        scored = [by_id[boundary_candidate_id(child)] for child in children]
        classes = [effective_class(result) for result in scored]
        if not set(classes) <= {"NCLDV", "MIRUS"} and not classes[0] == classes[1] == "PPV":
            continue
        compositions = [
            composition((item.gene, item.model) for item, _ in _qualified_core_items(result)) for result in scored
        ]
        if not may_merge(*compositions):
            continue
        accepted = [result.eve_id in accepted_ids for result in scored]
        if any(accepted) != all(accepted):
            continue
        if all(accepted):
            families = [qualified_core_families(result) for result in scored]
            if right.start - left.end > 5000 or not any(families):
                continue
        elif any(result.confidence_tier != "LOW" for result in scored):
            continue
        if any(not child.confirmed_seeds or child.joint_seed_provenance for child in children):
            continue
        if any(_crosses_stop(child, left.start, right.end) for child in children):
            continue
        if ("frameshift_rescue" in left.seed_sources) != ("frameshift_rescue" in right.seed_sources):
            continue
        if _has_third_overlap(ordered, children):
            continue
        if all(_has_multi_marker_witness(result, cls) for result, cls in zip(scored, classes, strict=True)):
            continue
        if not _hull_fully_searched(children, query, proteome_index, taxonomy):
            continue
        mappings = [
            [query.seed_gene_mappings.get(core.seed_id) for core in child.confirmed_seeds] for child in children
        ]
        envelopes = [
            (min(m.flank_start_bp for m in group if m is not None), max(m.flank_end_bp for m in group if m is not None))
            for group in mappings
        ]
        if not all(accepted) and max(span[0] for span in envelopes) >= min(span[1] for span in envelopes):
            continue
        genes = [g for g in proteome_index.get(left.scaffold, []) if g.start < right.end and g.end > left.start]
        if _gap_has_host_barrier(
            children,
            genes,
            taxonomy,
            host_model,
            host_threshold,
            viral_marker_genes,
        ):
            continue
        pairs.append(children)
    return pairs


def _hull_fully_searched(
    children: tuple[RefinedBoundary, RefinedBoundary],
    query: GenomeDiamondQuery | None,
    proteome_index: dict[str, list[pORF]],
    taxonomy: dict[str, GeneTaxonomy],
) -> bool:
    """Prove gene-level coverage before treating hull padding as searched DNA."""
    left, right = children
    if query is None or left.scaffold != right.scaffold or any(not child.confirmed_seeds for child in children):
        return False
    if any(query.seed_gene_mappings.get(core.seed_id) is None for child in children for core in child.confirmed_seeds):
        return False
    genes = [gene for gene in proteome_index.get(left.scaffold, []) if gene.start < right.end and gene.end > left.start]
    searched = set(query.all_porf_ids)
    return bool(genes) and all(gene.id in searched and gene.id in taxonomy for gene in genes)


def effective_class(result: VerificationResult) -> str:
    """Resolve the ordinary final gate class without borrowing a neighbor's label."""
    return resolve_effective_eve_class(
        confidence_tier=result.confidence_tier,
        region_classification=result.region_classification,
        classification=getattr(result, "classification", ""),
        likely_family=result.likely_family,
    )


def joint_boundary(
    children: tuple[RefinedBoundary, RefinedBoundary],
    *,
    rule: str = "rejected_adjacent_r2_v1",
    repeat_assessment: JointRepeatAssessment | None = None,
) -> RefinedBoundary:
    """Preserve exact cores and record the child geometry before recomputation."""
    provenance = {
        "rule": rule,
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
    if repeat_assessment is not None:
        provenance["repeat_assessment"] = asdict(repeat_assessment)
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


def qualified_core_families(result: VerificationResult) -> frozenset[str]:
    """Resolve qualifying core identities from selected-panel catalogue roles.

    MIRUS RNAP roles pool two subunits. A model's NCLDV catalogue link preserves
    their identity without interpreting every MIRUS RNAP hit as the large subunit.
    """
    families = set()
    for item, roles in _qualified_core_items(result):
        for role in roles:
            family = role["marker"]
            if role["panel"] == "Mirusviricota" and family == "RNAP":
                subunits = {
                    link["marker"]
                    for link in item.model_roles
                    if link["panel"] == "Nucleocytoviricota"
                    and link["marker"] in {"RNAPL", "RNAPS"}
                    and not link["excluded_from_scoring"]
                }
                families.update(subunits or {family})
            else:
                families.add(family)
    return frozenset(families)


def qualified_panel_families(result: VerificationResult) -> frozenset[str]:
    """Retain original selected-panel names alongside resolved core identities."""
    return frozenset(role["marker"] for _, roles in _qualified_core_items(result) for role in roles)


def accepted_joint_supported(result: VerificationResult, children: tuple[RefinedBoundary, RefinedBoundary]) -> bool:
    """Require native qualifying evidence wholly inside each original interval."""
    intervals = [(item.start, item.end) for item, _ in _qualified_core_items(result)]
    support = result.supplementary_viral_evidence
    if support is not None and support.eligible and support.group == effective_class(result):
        counted = set(support.counted_gene_ids)
        intervals.extend(
            (record["start"], record["end"]) for record in result.gene_taxonomy_records if record["porf_id"] in counted
        )
    return all(any(child.start <= start < end <= child.end for start, end in intervals) for child in children)


def _qualified_core_items(result: VerificationResult) -> list[tuple[HallmarkEvidence, list[dict]]]:
    summary = assess_hallmarks(result).summary
    qualified = set(summary.qualifying_evidence_ids)
    retained = []
    for item in result.hallmark_evidence:
        if item.evidence_id not in qualified:
            continue
        roles = [
            role
            for role in item.model_roles
            if role["panel"] == summary.panel and role["role"] == "core_marker" and not role["excluded_from_scoring"]
        ]
        if roles:
            retained.append((item, roles))
    return retained


def assess_joint_repeat_scan(
    children: tuple[RefinedBoundary, RefinedBoundary],
    sequence: str | None,
    markers: list[ValidatedMarkerHit],
    extension_bp: int,
    *,
    query: GenomeDiamondQuery | None,
    proteome_index: dict[str, list[pORF]],
    taxonomy: dict[str, GeneTaxonomy],
    accepted_ids: frozenset[str] = frozenset(),
) -> JointRepeatAssessment | None:
    """Withhold separating repeats; record permitted uncertainty in accepted children."""
    if sequence is None or query is None or not _hull_fully_searched(children, query, proteome_index, taxonomy):
        return None
    left, right = children
    if len(sequence) < right.end:
        return None
    accepted_children = all(boundary_candidate_id(child) in accepted_ids for child in children)
    gap_start = max(core.end for core in left.confirmed_seeds)
    gap_end = min(core.start for core in right.confirmed_seeds)
    inter_core_ambiguous_bp = sum(base not in "ACGT" for base in sequence[gap_start:gap_end].upper())
    if inter_core_ambiguous_bp and not accepted_children:
        return None
    if not set(sequence[left.end : right.start].upper()) <= set("ACGT"):
        return None
    scan_start, scan_end = max(0, left.start - extension_bp), min(len(sequence), right.end + extension_bp)
    spans = []
    for child in children:
        retained = [
            m
            for m in markers
            if m.scaffold == child.scaffold
            and m.validation_status in {"validated", "validated_novel"}
            and is_ordinary_marker_evidence(m)
            and child.start <= m.start < m.end <= child.end
        ]
        if not retained:
            return None
        spans.append((min(m.start for m in retained), max(m.end for m in retained)))
    gap_sides = []
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
                _LOGGER.info(
                    "withheld_unresolved_repeat type=TIR reason=pairing_capped children=%s",
                    [boundary_candidate_id(c) for c in children],
                )
                return None
            other_start, other_end = spans[1 - index]
            for repeat in scan.candidates:
                if (repeat.right_end <= other_start and left.start < repeat.right_end < right.end) or (
                    repeat.left_start >= other_end and left.start < repeat.left_start < right.end
                ):
                    _LOGGER.info(
                        "withheld_unresolved_repeat type=TIR children=%s arms=%s",
                        [boundary_candidate_id(c) for c in children],
                        repeat,
                    )
                    return None
        coverage = tuple(
            (mapping.flank_start_bp, mapping.flank_end_bp)
            for mapping in query.seed_gene_mappings.values()
            if mapping.scaffold == left.scaffold
        ) + ((left.start, right.end),)
        for index, child in enumerate(children):
            direct = assess_repeat_evidence(
                sequence,
                scaffold=child.scaffold,
                input_id=boundary_candidate_id(child),
                start=child.start,
                end=child.end,
                parent_start=child.start,
                parent_end=child.end,
                coverage_intervals=coverage,
                extension_bp=extension_bp,
            )
            gap_side = direct.right if index == 0 else direct.left
            if gap_side.status != "completed" and not (
                accepted_children and _ambiguity_confined_to_children(gap_side, sequence, children)
            ):
                _LOGGER.info(
                    "withheld_unresolved_repeat type=DTR children=%s gap_side=%s",
                    [boundary_candidate_id(c) for c in children],
                    gap_side,
                )
                return None
            gap_sides.append(gap_side)
            start, end = spans[index]
            other_start, other_end = spans[1 - index]
            for repeat in direct.candidates:
                encloses = repeat.left_end <= start and end <= repeat.right_start
                separates = (index == 0 and repeat.right_end <= other_start) or (
                    index == 1 and repeat.left_start >= other_end
                )
                if encloses and separates:
                    _LOGGER.info(
                        "withheld_unresolved_repeat type=DTR children=%s arms=%s",
                        [boundary_candidate_id(c) for c in children],
                        repeat,
                    )
                    return None
    except (OSError, RuntimeError, ValueError):
        return None
    uncertain = inter_core_ambiguous_bp > 0 or any(side.status != "completed" for side in gap_sides)
    return JointRepeatAssessment(
        "accepted_child_ambiguity" if uncertain else "completed",
        (left.end, right.start),
        inter_core_ambiguous_bp,
        (gap_sides[0], gap_sides[1]),
    )


def _ambiguity_confined_to_children(
    side: EndAssessment,
    sequence: str,
    children: tuple[RefinedBoundary, RefinedBoundary],
) -> bool:
    """Allow only ambiguous child sequence, never another incomplete-scan reason."""
    if (
        side.status != "incomplete"
        or set(side.reasons) - {"midpoint_partition"} != {"ambiguous_bases"}
        or side.window_start is None
        or side.window_end is None
    ):
        return False
    window = sequence[side.window_start : side.window_end].upper()
    # Native end assessment skips complexity when Ns are present. The exception
    # must not turn a low-complexity window into apparently usable evidence.
    if _is_low_complexity("".join(base for base in window if base in "ACGT")):
        return False
    return all(
        base in "ACGT" or any(child.start <= position < child.end for child in children)
        for position, base in enumerate(window, side.window_start)
    )


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
    panels = {
        "NCLDV": {"Nucleocytoviricota"},
        "MIRUS": {"Mirusviricota"},
        "PPV": {"Polinton_like_viruses", "Maverick_Polintons", "Virophages"},
    }[family]
    return any(
        w.panel in panels and w.rule in {"capsid_partner_5kb", "three_core_15kb"}
        for w in assess_hallmarks(result).witnesses
    )


def _gap_has_host_barrier(
    children: tuple[RefinedBoundary, RefinedBoundary],
    genes: list[pORF],
    taxonomy: dict[str, GeneTaxonomy],
    model: HostSignatureModel,
    threshold: float,
    viral_marker_genes: frozenset[str],
) -> bool:
    left, right = children
    run = 0
    for gene in sorted(genes, key=lambda g: (g.start, g.end)):
        if gene.end <= left.end or gene.start >= right.start:
            continue
        record = taxonomy[gene.id]
        host = is_extension_host_gene(
            record,
            host_share_ok=host_lineage_share(record, model) >= threshold,
            marker=gene.id in viral_marker_genes,
        )
        run = run + 1 if host else 0
        if run == 3:
            return True
    return False
