"""Retained profile evidence for graduated hallmark scores and spatial seeds."""

from __future__ import annotations

import hashlib
import json
from bisect import bisect_left, bisect_right
from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from functools import lru_cache
from itertools import accumulate, combinations
from pathlib import Path
from typing import TYPE_CHECKING

from virosync.pipeline.phase1.hhg_seeding import HMMHit
from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed

if TYPE_CHECKING:
    from virosync.pipeline.phase3.evidence_synthesizer import VerificationResult

CATALOGUE_PATH = Path(__file__).resolve().parents[1] / "data" / "hallmark_panels.json"
STRONG_EVALUE = 1e-5
WEAK_EVALUE = 0.01
VIRAL_PREFIXES = ("NCLDV__", "MIRUS__", "PLV__", "VP__", "PPV__", "CRESS__", "EUK-pEVE")
EXCLUDED_FAMILIES = {"SF1/SF3_fusion", "PrimPol/SF3"}


@dataclass(frozen=True, slots=True)
class HallmarkEvidence:
    """One retained model record with unchanged validation and measured diagnostics."""

    query_porf: str
    gene: str
    scaffold: str
    start: int
    end: int
    model: str
    evalue: float
    score: float
    domain_score: float | None
    query_start: int | None
    query_end: int | None
    validation_status: str
    viral_top10_max_pident: float
    source: str

    @property
    def evidence_id(self) -> str:
        """Return a deterministic identity for this retained model and envelope."""
        key = (self.query_porf, self.model, self.scaffold, self.start, self.end)
        return hashlib.sha256(json.dumps(key).encode()).hexdigest()[:24]

    @property
    def model_roles(self) -> list[dict]:
        """Return the canonical panel roles without qualifying or validating the hit."""
        fields = ("panel", "marker", "module", "role", "capsid", "structural_anchor", "excluded_from_scoring")
        return [{name: link[name] for name in fields} for link in _catalogue().get(("phase1_marker", self.model), ())]


@dataclass(frozen=True, slots=True)
class HallmarkWitness:
    """A spatial rule and its exact retained model members."""

    scaffold: str
    panel: str
    rule: str
    tier: str
    start: int
    end: int
    evidence_ids: tuple[str, ...]

    @property
    def witness_id(self) -> str:
        """Return a stable source identity independent of input order."""
        key = (self.scaffold, self.panel, self.rule, self.evidence_ids)
        return "hallmark_" + hashlib.sha256(json.dumps(key).encode()).hexdigest()[:24]


@dataclass(frozen=True, slots=True)
class HallmarkSupport:
    """A contained functional or integration profile supplying context only."""

    gene: str
    model: str
    evalue: float
    stage: str


@dataclass(frozen=True, slots=True)
class HallmarkSummary:
    """Selected compatible panel and its graduated score components."""

    panel: str = "none"
    core_count: int = 0
    strong_count: int = 0
    weighted_core_count: float = 0.0
    compatible_count: int = 0
    capsid_fraction: float = 0.0
    marker_score: float = 0.0
    coherence_bonus: float = 0.0
    floor_tier: str = ""
    qualifying_evidence_ids: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class HallmarkAssessment:
    """One final-coordinate assessment shared by scoring, gates and export."""

    scaffold: str
    start: int
    end: int
    summary: HallmarkSummary
    witnesses: tuple[HallmarkWitness, ...]


def assess_hallmarks(result: VerificationResult) -> HallmarkAssessment:
    """Assess only complete genes in the candidate's final retained interval."""
    retained = [
        h
        for h in result.hallmark_evidence
        if h.scaffold == result.scaffold and result.start <= h.start and h.end <= result.end
    ]
    return HallmarkAssessment(
        result.scaffold,
        result.start,
        result.end,
        summarize_hallmarks(retained, result.hallmark_support),
        tuple(contained_witnesses(retained, result.hallmark_witnesses, result.scaffold, result.start, result.end)),
    )


def collect_hallmark_evidence(
    markers: Sequence[ValidatedMarkerHit],
    raw_hits: Sequence[HMMHit],
    *,
    frameshift_ids: set[str] | None = None,
) -> list[HallmarkEvidence]:
    """Join exact retained identities; never restore discarded model alternatives."""
    raw = {(f"{h.query_name}|aa{h.query_start}-{h.query_end}", h.target_name): h for h in raw_hits}
    catalogue = _catalogue()
    evidence = []
    for marker in markers:
        if ("phase1_marker", marker.hmm_target) not in catalogue:
            continue
        gene = marker.query_porf.split("|aa", 1)[0]
        measured = raw.get((marker.query_porf, marker.hmm_target))
        frameshift = gene in (frameshift_ids or set())
        if measured is None and not frameshift:
            raise ValueError(f"Retained hallmark lacks exact HMM source: {marker.query_porf} {marker.hmm_target}")
        identities = zip(marker.top10_prefixes.split(","), marker.top10_pidents.split(","), strict=False)
        viral_identity = max(
            (float(value) for prefix, value in identities if value and prefix.startswith(VIRAL_PREFIXES)), default=0.0
        )
        evidence.append(
            HallmarkEvidence(
                query_porf=marker.query_porf,
                gene=gene,
                scaffold=marker.scaffold,
                start=marker.start,
                end=marker.end,
                model=marker.hmm_target,
                evalue=marker.hmm_evalue,
                score=marker.hmm_score,
                domain_score=None if frameshift else measured.domain_score,
                query_start=None if frameshift else measured.query_start,
                query_end=None if frameshift else measured.query_end,
                validation_status=marker.validation_status,
                viral_top10_max_pident=viral_identity,
                source="frameshift" if frameshift else "retained_hmm",
            )
        )
    return sorted(evidence, key=lambda h: (h.scaffold, h.start, h.end, h.query_porf, h.model))


def find_spatial_witnesses(evidence: Sequence[HallmarkEvidence]) -> list[HallmarkWitness]:
    """Find independent full-ORF witnesses, including contained three-core windows."""
    grouped = defaultdict(list)
    for hit, link in _linked(evidence):
        if link["role"] == "core_marker":
            grouped[hit.scaffold, link["panel"]].append((hit, link))
    witnesses = {}

    def add(panel: str, selected: list[tuple[HallmarkEvidence, dict]], rule: str, tier: str) -> None:
        members = tuple(
            h.evidence_id
            for h, _ in sorted(selected, key=lambda pair: (pair[0].gene, pair[1]["marker"], pair[0].model))
        )
        witness = HallmarkWitness(
            selected[0][0].scaffold,
            panel,
            rule,
            tier,
            min(h.start for h, _ in selected),
            max(h.end for h, _ in selected),
            members,
        )
        witnesses[witness.witness_id] = witness

    for (_, panel), group in sorted(grouped.items()):
        ordered = sorted(range(len(group)), key=lambda i: group[i][0].start)
        starts = [group[i][0].start for i in ordered]
        prefix_ends = list(accumulate((group[i][0].end for i in ordered), max))
        families = {link["marker"] for _, link in group}
        if len(families) >= 3:
            for left in sorted(set(starts)):
                lower, upper = bisect_left(starts, left), bisect_right(starts, left + 15000)
                # Preserve input order for equal-score/model matching ties.
                local = [group[i] for i in sorted(ordered[lower:upper]) if group[i][0].end <= left + 15000]
                if len({link["marker"] for _, link in local}) < 3:
                    continue
                for right in sorted({h.end for h, _ in local}):
                    selected = _assignment([(h, link) for h, link in local if h.end <= right])
                    if len(selected) >= 3:
                        add(panel, selected, "three_core_15kb", "MEDIUM")
                        break
        for capsid, link in group:
            if not link["capsid"]:
                continue
            add(panel, [(capsid, link)], "capsid", "MEDIUM")
            if len(families) < 2:
                continue
            # Prefix maxima retain long genes starting well before the capsid.
            lower = bisect_left(prefix_ends, capsid.start - 5000)
            upper = bisect_right(starts, capsid.end + 5000)
            for index in ordered[lower:upper]:
                partner, other = group[index]
                if partner.gene == capsid.gene or other["marker"] == link["marker"]:
                    continue
                gap = max(0, partner.start - capsid.end, capsid.start - partner.end)
                if gap <= 5000:
                    add(panel, [(capsid, link), (partner, other)], "capsid_partner_5kb", "HIGH")
    return sorted(witnesses.values(), key=lambda w: (w.scaffold, w.start, w.end, w.panel, w.rule, w.evidence_ids))


def add_rescue_seeds(seeds: list[MergedSeed], witnesses: Sequence[HallmarkWitness]) -> list[MergedSeed]:
    """Append exact witness cores without changing or renumbering ordinary seeds."""
    output = list(seeds)
    candidates = [
        w
        for w in witnesses
        if not any(seed.scaffold == w.scaffold and seed.start <= w.start and seed.end >= w.end for seed in seeds)
    ]
    for witness in candidates:
        if any(
            other.scaffold == witness.scaffold
            and other.start <= witness.start
            and other.end >= witness.end
            and ((other.start, other.end) != (witness.start, witness.end) or other.witness_id < witness.witness_id)
            for other in candidates
        ):
            continue
        identity = witness.witness_id
        output.append(
            MergedSeed(
                witness.scaffold,
                witness.start,
                witness.end,
                seed_id=identity,
                confirmed_seeds=(ConfirmedSeed(identity, witness.scaffold, witness.start, witness.end),),
                sources=["hallmark_rescue"],
                confidence=witness.tier.lower(),
            )
        )
    return output


def contained_witnesses(
    evidence: Sequence[HallmarkEvidence], witnesses: Sequence[HallmarkWitness], scaffold: str, start: int, end: int
) -> list[HallmarkWitness]:
    """Retain original witnesses whose complete valid members remain inside."""
    retained = [h for h in evidence if h.scaffold == scaffold and start <= h.start and h.end <= end]
    by_id = {h.evidence_id: h for h in retained}
    contained = []
    for witness in witnesses:
        if witness.scaffold != scaffold or not set(witness.evidence_ids).issubset(by_id):
            continue
        selected = [by_id[key] for key in witness.evidence_ids]
        if witness in find_spatial_witnesses(selected):
            contained.append(witness)
    return contained


def has_spatial_rescue(result: VerificationResult) -> bool:
    """Return whether a scored candidate retains an authenticated spatial floor."""
    if result.confidence_tier not in {"MEDIUM", "HIGH"}:
        return False
    assessment = getattr(result, "hallmark_assessment", None)
    if not assessment or not assessment.witnesses:
        return False
    if (assessment.scaffold, assessment.start, assessment.end) != (result.scaffold, result.start, result.end):
        return False
    if any(w.tier == "HIGH" for w in assessment.witnesses) and result.confidence_tier != "HIGH":
        return False
    floor = result.score_components.get("hallmark_floor")
    return floor is not None and result.final_confidence >= floor


def summarize_hallmarks(
    evidence: Sequence[HallmarkEvidence], support: Sequence[HallmarkSupport] = ()
) -> HallmarkSummary:
    """Select one compatible panel, then apply the approved graduated weights."""
    core, context = defaultdict(list), defaultdict(list)
    for hit, link in _linked(evidence):
        if (
            hit.validation_status != "validated_novel"
            and hit.viral_top10_max_pident < 25
            and not link["structural_anchor"]
        ):
            continue
        context[link["panel"]].append((hit, link))
        if link["role"] == "core_marker":
            core[link["panel"]].append((hit, link))
    for hit in support:
        if hit.evalue > STRONG_EVALUE:
            continue
        for link in _catalogue().get((hit.stage, hit.model), ()):
            if not link["excluded_from_scoring"] and link["marker"] not in EXCLUDED_FAMILIES:
                context[link["panel"]].append((hit, link))
    choices = []
    for panel, hits in core.items():
        total = len(_assignment(hits))
        strong = len(_assignment([(h, link) for h, link in hits if h.evalue <= STRONG_EVALUE]))
        any_capsid = any(link["capsid"] for _, link in hits)
        strong_capsid = any(link["capsid"] and h.evalue <= STRONG_EVALUE for h, link in hits)
        compatible = len(_assignment(context[panel]))
        if total >= 3:
            effective, weak_weight = float(total), 1.0
        elif total == 2:
            pair_strong = max(
                (
                    int(a.evalue <= STRONG_EVALUE) + int(b.evalue <= STRONG_EVALUE)
                    for (a, left), (b, right) in combinations(hits, 2)
                    if a.gene != b.gene and left["marker"] != right["marker"]
                ),
                default=0,
            )
            effective, weak_weight = pair_strong + (2 - pair_strong) * 0.5, 0.5
        else:
            effective, weak_weight = float(strong) + (total - strong) * 0.25, 0.25
        capsid_fraction = 1.0 if strong_capsid else weak_weight if any_capsid else 0.0
        # The prototype adds a fractional capsid contribution after the capped core term.
        marker_score = min(0.2 * effective + (0.4 if capsid_fraction == 1.0 else 0.0), 1.0)
        if 0 < capsid_fraction < 1:
            marker_score += 0.4 * capsid_fraction
        choices.append(
            HallmarkSummary(
                panel,
                total,
                strong,
                effective,
                compatible,
                capsid_fraction,
                marker_score,
                0.05 * min(max(compatible - 1, 0), 3) * effective / total,
                "HIGH" if any_capsid and total >= 2 else "MEDIUM" if any_capsid else "",
                tuple(sorted({h.evidence_id for h, _ in hits})),
            )
        )
    return max(
        choices,
        key=lambda c: (c.core_count, c.strong_count, bool(c.capsid_fraction), c.compatible_count, c.panel),
        default=HallmarkSummary(),
    )


@lru_cache(maxsize=1)
def _catalogue() -> dict[tuple[str, str], list[dict]]:
    indexed = defaultdict(list)
    for link in json.loads(CATALOGUE_PATH.read_text())["links"]:
        indexed[link["stage"], link["model"]].append(link)
    return dict(indexed)


def _linked(evidence: Sequence[HallmarkEvidence]) -> Iterable[tuple[HallmarkEvidence, dict]]:
    for hit in evidence:
        if hit.evalue > WEAK_EVALUE:
            continue
        for link in hit.model_roles:
            if not link["excluded_from_scoring"] and link["marker"] not in EXCLUDED_FAMILIES:
                yield hit, link


def _assignment(
    hits: Sequence[tuple[HallmarkEvidence | HallmarkSupport, dict]],
) -> list[tuple[HallmarkEvidence | HallmarkSupport, dict]]:
    edges = defaultdict(dict)
    for hit, link in sorted(hits, key=lambda pair: (pair[0].evalue, pair[0].model)):
        edges[link["marker"]].setdefault(hit.gene, (hit, link))
    matched = {}

    def assign(family: str, seen: set[str]) -> bool:
        for gene in sorted(edges[family]):
            if gene in seen:
                continue
            seen.add(gene)
            if gene not in matched or assign(matched[gene], seen):
                matched[gene] = family
                return True
        return False

    for family in sorted(edges):
        assign(family, set())
    return [edges[family][gene] for gene, family in sorted(matched.items())]
