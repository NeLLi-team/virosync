"""Detect integration-associated genes in EVE candidates and their flanks."""

from __future__ import annotations

import csv
import json
import re
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import pyhmmer
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMMFile

from virosync.pipeline.phase0.prodigal import parse_prodigal_header
from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit
from virosync.pipeline.phase1.viral_markers import base_marker_gene_id

__all__ = [
    "INTEGRATION_PROFILE_MANIFEST_PATH",
    "INTEGRATION_PROFILE_PATH",
    "IntegrationGeneHit",
    "IntegrationRegion",
    "IntegrationScanResult",
    "IntegrationUnsearchedProtein",
    "collect_annotated_integration_genes",
    "scan_integration_genes",
]

INTEGRATION_PROFILE_PATH = Path(__file__).parents[2] / "data" / "integration_profiles.hmm"
INTEGRATION_PROFILE_MANIFEST_PATH = Path(__file__).parents[2] / "data" / "integration_profiles.json"
MAX_INTEGRATION_PROTEIN_LENGTH = 100_000

IntegrationLocation = Literal["interior", "upstream", "downstream", "boundary_overlap"]

_ANNOTATION_PATTERNS = (
    (re.compile(r"\btyrosine[\s-]+recombinase\b", re.IGNORECASE), "tyrosine_recombinase"),
    (re.compile(r"\bserine[\s-]+recombinase\b", re.IGNORECASE), "serine_recombinase"),
    (re.compile(r"\bintegrase\b", re.IGNORECASE), "integrase"),
    (re.compile(r"\brecombinase\b", re.IGNORECASE), "recombinase"),
)


@dataclass(frozen=True, slots=True)
class IntegrationRegion:
    """One final EVE interval and the wider interval screened around it.

    Coordinates use ViroSync's 0-based half-open convention.
    """

    region_id: str
    scaffold: str
    scan_start: int
    scan_end: int
    interior_start: int
    interior_end: int


@dataclass(frozen=True, slots=True)
class IntegrationGeneHit:
    """One integration-associated gene annotation tied to an EVE region.

    Genomic coordinates are 0-based half-open. Domain coordinates are the
    1-based inclusive amino-acid envelope coordinates reported by HMMER.
    """

    region_id: str
    protein_id: str
    scaffold: str
    start: int
    end: int
    strand: str
    location: IntegrationLocation
    mechanism: str
    source: str
    profile: str
    accession: str
    sequence_score: float | None = None
    domain_score: float | None = None
    evalue: float | None = None
    domain_start: int | None = None
    domain_end: int | None = None
    validation_status: str = ""
    annotation: str = ""


@dataclass(frozen=True, slots=True)
class IntegrationUnsearchedProtein:
    """One original protein excluded from HMM screening in an EVE context.

    Genomic coordinates are 0-based half-open. This is missing assessment,
    not evidence that the protein lacks an integration domain.
    """

    region_id: str
    protein_id: str
    length_aa: int
    scaffold: str
    start: int
    end: int
    strand: str
    location: IntegrationLocation
    reason: Literal["sequence_length_limit"] = "sequence_length_limit"
    limit_aa: int = MAX_INTEGRATION_PROTEIN_LENGTH


@dataclass(frozen=True, slots=True)
class IntegrationScanResult:
    """Separate positive HMM evidence from proteins the engine cannot assess."""

    hits: list[IntegrationGeneHit]
    unsearched: list[IntegrationUnsearchedProtein]


@dataclass(frozen=True, slots=True)
class _GeneContext:
    """Candidate gene coordinates in one search region."""

    region: IntegrationRegion
    protein_id: str
    scaffold: str
    start: int
    end: int
    strand: str


@dataclass(frozen=True, slots=True)
class _ProfileMetadata:
    """Manifest metadata for one bundled integration profile."""

    accession: str
    name: str
    role: str
    description: str


def scan_integration_genes(
    proteome_path: Path,
    regions: Sequence[IntegrationRegion],
    *,
    threads: int,
    hmm_path: Path = INTEGRATION_PROFILE_PATH,
) -> IntegrationScanResult:
    """Search candidate and flanking proteins with curated profile GA cutoffs.

    Selected proteins of at most 100,000 residues are searched in one batch.
    Longer proteins are recorded per EVE context without modifying sequences.
    A hit must pass both sequence and best-domain gathering thresholds.
    Empty selections have neither hits nor unsearched records.
    """
    if not regions:
        return IntegrationScanResult([], [])

    contexts_by_protein, sequences = _load_candidate_genes(proteome_path, regions)
    supported_sequences: list[pyhmmer.easel.DigitalSequence] = []
    unsearched: list[IntegrationUnsearchedProtein] = []
    for sequence in sequences:
        if len(sequence) <= MAX_INTEGRATION_PROTEIN_LENGTH:
            supported_sequences.append(sequence)
            continue
        for context in contexts_by_protein[_decode(sequence.name)]:
            unsearched.append(
                IntegrationUnsearchedProtein(
                    region_id=context.region.region_id,
                    protein_id=context.protein_id,
                    length_aa=len(sequence),
                    scaffold=context.scaffold,
                    start=context.start,
                    end=context.end,
                    strand=context.strand,
                    location=_gene_location(context.start, context.end, context.region),
                )
            )
    if not supported_sequences:
        return IntegrationScanResult([], unsearched)

    with HMMFile(hmm_path) as handle:
        profiles = list(handle)
    metadata_by_name, _metadata_by_accession = _load_profile_metadata()
    _validate_profiles(profiles, metadata_by_name)

    hits: list[IntegrationGeneHit] = []
    # Count original targets, not EVE contexts. Excluded comparisons are not
    # performed; this conservative normalization retains their multiplicity.
    search_space = {"Z": len(contexts_by_protein)} if unsearched else {}
    searches = pyhmmer.hmmsearch(
        profiles,
        supported_sequences,
        cpus=threads,
        parallel="targets",
        E=float("inf"),
        domE=float("inf"),
        incE=float("inf"),
        incdomE=float("inf"),
        **search_space,
    )
    for profile, top_hits in zip(profiles, searches, strict=True):
        profile_name = _decode(profile.name)
        metadata = metadata_by_name[profile_name]
        accession = _decode(profile.accession)
        sequence_cutoff, domain_cutoff = _gathering_cutoffs(profile)
        for top_hit in top_hits:
            best_domain = max(top_hit.domains, key=lambda domain: domain.score)
            if top_hit.score < sequence_cutoff or best_domain.score < domain_cutoff:
                continue
            protein_id = _decode(top_hit.name)
            for context in contexts_by_protein.get(protein_id, ()):
                hits.append(
                    _make_hit(
                        context,
                        mechanism=metadata.role,
                        source="pfam_hmm",
                        profile=profile_name,
                        accession=accession,
                        sequence_score=float(top_hit.score),
                        domain_score=float(best_domain.score),
                        evalue=float(top_hit.evalue),
                        domain_start=int(best_domain.env_from),
                        domain_end=int(best_domain.env_to),
                        annotation=metadata.description,
                    )
                )
    return IntegrationScanResult(_sorted_unique_hits(hits), unsearched)


def collect_annotated_integration_genes(
    proteome_path: Path,
    validated_markers: Sequence[ValidatedMarkerHit],
    regions: Sequence[IntegrationRegion],
    *,
    model_annotations_path: Path | None,
    interproscan_path: Path | None = None,
) -> list[IntegrationGeneHit]:
    """Collect integration annotations already present in marker metadata.

    Model annotations use only explicit Pfam signatures and the description or
    majority-annotation text. The incidental ``pfam_top_domains`` count column
    is deliberately ignored. Optional raw InterProScan rows add independent
    source records, including reported domain coordinates. Pfam's score field
    is retained as an E-value; other member databases do not share that score
    meaning and are kept as annotation evidence without a numeric score.
    """
    if not regions:
        return []

    metadata_by_name, metadata_by_accession = _load_profile_metadata()
    hits: list[IntegrationGeneHit] = []
    if model_annotations_path is not None and model_annotations_path.exists():
        annotations = _load_model_annotations(model_annotations_path)
        hits.extend(
            _marker_annotation_hits(
                validated_markers,
                regions,
                annotations,
                metadata_by_name,
                metadata_by_accession,
            )
        )
    if interproscan_path is not None and interproscan_path.exists():
        contexts_by_protein, _sequences = _load_candidate_genes(proteome_path, regions)
        hits.extend(
            _interproscan_hits(
                interproscan_path,
                regions,
                contexts_by_protein,
                metadata_by_name,
                metadata_by_accession,
            )
        )
    return _sorted_unique_hits(hits)


def _load_candidate_genes(
    proteome_path: Path,
    regions: Sequence[IntegrationRegion],
) -> tuple[dict[str, list[_GeneContext]], list[pyhmmer.easel.DigitalSequence]]:
    regions_by_scaffold: dict[str, list[IntegrationRegion]] = defaultdict(list)
    for region in regions:
        regions_by_scaffold[region.scaffold].append(region)

    contexts_by_protein: dict[str, list[_GeneContext]] = defaultdict(list)
    sequences: list[pyhmmer.easel.DigitalSequence] = []
    alphabet = pyhmmer.easel.Alphabet.amino()
    with SequenceFile(proteome_path, digital=True, alphabet=alphabet) as handle:
        for sequence in handle:
            protein_id = _decode(sequence.name)
            description = _decode(sequence.description)
            parsed = parse_prodigal_header(f"{protein_id} {description}".rstrip(), protein_id)
            if parsed is None:
                continue
            scaffold, start, end, strand = parsed
            matching_regions = _regions_for_gene(regions_by_scaffold.get(scaffold, ()), start, end)
            if not matching_regions:
                continue
            sequences.append(sequence)
            for region in matching_regions:
                contexts_by_protein[protein_id].append(
                    _GeneContext(
                        region=region,
                        protein_id=protein_id,
                        scaffold=scaffold,
                        start=start,
                        end=end,
                        strand=strand,
                    )
                )
    return dict(contexts_by_protein), sequences


def _regions_for_gene(
    regions: Sequence[IntegrationRegion],
    start: int,
    end: int,
) -> list[IntegrationRegion]:
    return [region for region in regions if start < region.scan_end and end > region.scan_start]


def _load_profile_metadata() -> tuple[dict[str, _ProfileMetadata], dict[str, _ProfileMetadata]]:
    document = json.loads(INTEGRATION_PROFILE_MANIFEST_PATH.read_text(encoding="utf-8"))
    metadata = [
        _ProfileMetadata(
            accession=profile["accession"],
            name=profile["name"],
            role=profile["role"],
            description=profile["description"],
        )
        for profile in document["profiles"]
    ]
    return (
        {profile.name: profile for profile in metadata},
        {profile.accession.split(".", 1)[0]: profile for profile in metadata},
    )


def _validate_profiles(
    profiles: Sequence[pyhmmer.plan7.HMM],
    metadata_by_name: dict[str, _ProfileMetadata],
) -> None:
    for profile in profiles:
        name = _decode(profile.name)
        if name not in metadata_by_name:
            raise ValueError(f"unsupported integration profile {name!r}")
        if not profile.cutoffs.gathering_available():
            raise ValueError(f"integration profile {name!r} has no gathering cutoff")


def _gathering_cutoffs(profile: pyhmmer.plan7.HMM) -> tuple[float, float]:
    cutoffs = profile.cutoffs.gathering
    if cutoffs is None:
        raise ValueError(f"integration profile {_decode(profile.name)!r} has no gathering cutoff")
    return cutoffs


def _load_model_annotations(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"model_name", "description", "majority_annotation", "pfam_signature"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"model annotation table is missing columns: {sorted(missing)}")
        return {row["model_name"].strip(): row for row in reader if row["model_name"].strip()}


def _marker_annotation_hits(
    markers: Sequence[ValidatedMarkerHit],
    regions: Sequence[IntegrationRegion],
    annotations: dict[str, dict[str, str]],
    metadata_by_name: dict[str, _ProfileMetadata],
    metadata_by_accession: dict[str, _ProfileMetadata],
) -> list[IntegrationGeneHit]:
    regions_by_scaffold: dict[str, list[IntegrationRegion]] = defaultdict(list)
    for region in regions:
        regions_by_scaffold[region.scaffold].append(region)

    hits: list[IntegrationGeneHit] = []
    for marker in markers:
        model = marker.hmm_target
        row = annotations.get(model)
        if row is None:
            continue
        classification = _classify_annotation(
            pfam_signature=row["pfam_signature"],
            annotation=" ".join((row["majority_annotation"], row["description"])),
            metadata_by_name=metadata_by_name,
            metadata_by_accession=metadata_by_accession,
        )
        if classification is None:
            continue
        mechanism, accession, annotation = classification
        scaffold = marker.scaffold
        start = marker.start
        end = marker.end
        protein_id = base_marker_gene_id(marker.query_porf)
        for region in _regions_for_gene(regions_by_scaffold.get(scaffold, ()), start, end):
            context = _GeneContext(
                region=region,
                protein_id=protein_id,
                scaffold=scaffold,
                start=start,
                end=end,
                strand=marker.strand,
            )
            hits.append(
                _make_hit(
                    context,
                    mechanism=mechanism,
                    source="marker_annotation",
                    profile=model,
                    accession=accession,
                    sequence_score=marker.hmm_score,
                    evalue=marker.hmm_evalue,
                    validation_status=marker.validation_status,
                    annotation=annotation,
                )
            )
    return hits


def _interproscan_hits(
    path: Path,
    regions: Sequence[IntegrationRegion],
    contexts_by_protein: dict[str, list[_GeneContext]],
    metadata_by_name: dict[str, _ProfileMetadata],
    metadata_by_accession: dict[str, _ProfileMetadata],
) -> list[IntegrationGeneHit]:
    region_prefixes = sorted(((f"{region.region_id}|", region.region_id) for region in regions), reverse=True)
    hits: list[IntegrationGeneHit] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13:
                continue
            query_id = fields[0]
            split_query = _split_interproscan_query(query_id, region_prefixes)
            if split_query is None:
                continue
            protein_id, query_region_id = split_query
            classification = _classify_annotation(
                pfam_signature=fields[4],
                annotation=" ".join((fields[5], fields[11], fields[12])),
                metadata_by_name=metadata_by_name,
                metadata_by_accession=metadata_by_accession,
            )
            if classification is None:
                continue
            mechanism, accession, annotation = classification
            evalue = _optional_float(fields[8]) if fields[3].casefold() == "pfam" else None
            for context in contexts_by_protein.get(protein_id, ()):
                if context.region.region_id != query_region_id:
                    continue
                hits.append(
                    _make_hit(
                        context,
                        mechanism=mechanism,
                        source="interproscan",
                        profile=fields[4] or fields[3],
                        accession=accession,
                        evalue=evalue,
                        domain_start=_optional_int(fields[6]),
                        domain_end=_optional_int(fields[7]),
                        annotation=annotation,
                    )
                )
    return hits


def _split_interproscan_query(
    query_id: str,
    prefixes: Sequence[tuple[str, str]],
) -> tuple[str, str] | None:
    for prefix, region_id in prefixes:
        if query_id.startswith(prefix):
            return query_id[len(prefix) :], region_id
    return None


def _classify_annotation(
    *,
    pfam_signature: str,
    annotation: str,
    metadata_by_name: dict[str, _ProfileMetadata],
    metadata_by_accession: dict[str, _ProfileMetadata],
) -> tuple[str, str, str] | None:
    normalized_annotation = annotation.strip()
    if re.search(r"\byqaj\b", normalized_annotation, re.IGNORECASE):
        return "repair_recombinase", "", normalized_annotation

    for token in pfam_signature.split(";"):
        normalized = token.strip()
        accession = normalized.split(".", 1)[0]
        if normalized.casefold() == "yqaj":
            return "repair_recombinase", "PF09588", normalized_annotation or "YqaJ recombinase"
        if accession == "PF09588":
            return "repair_recombinase", normalized, normalized_annotation or "YqaJ recombinase"
        metadata = metadata_by_accession.get(accession) or metadata_by_name.get(normalized)
        if metadata is not None:
            return metadata.role, metadata.accession, normalized_annotation or metadata.description

    for pattern, mechanism in _ANNOTATION_PATTERNS:
        if pattern.search(normalized_annotation):
            return mechanism, "", normalized_annotation
    return None


def _make_hit(
    context: _GeneContext,
    *,
    mechanism: str,
    source: str,
    profile: str,
    accession: str,
    sequence_score: float | None = None,
    domain_score: float | None = None,
    evalue: float | None = None,
    domain_start: int | None = None,
    domain_end: int | None = None,
    validation_status: str = "",
    annotation: str = "",
) -> IntegrationGeneHit:
    return IntegrationGeneHit(
        region_id=context.region.region_id,
        protein_id=context.protein_id,
        scaffold=context.scaffold,
        start=context.start,
        end=context.end,
        strand=context.strand,
        location=_gene_location(context.start, context.end, context.region),
        mechanism=mechanism,
        source=source,
        profile=profile,
        accession=accession,
        sequence_score=sequence_score,
        domain_score=domain_score,
        evalue=evalue,
        domain_start=domain_start,
        domain_end=domain_end,
        validation_status=validation_status,
        annotation=annotation,
    )


def _gene_location(start: int, end: int, region: IntegrationRegion) -> IntegrationLocation:
    if end <= region.interior_start:
        return "upstream"
    if start >= region.interior_end:
        return "downstream"
    if start >= region.interior_start and end <= region.interior_end:
        return "interior"
    return "boundary_overlap"


def _sorted_unique_hits(hits: Sequence[IntegrationGeneHit]) -> list[IntegrationGeneHit]:
    best_by_key: dict[tuple[str, str, str, str], IntegrationGeneHit] = {}
    for hit in hits:
        key = (hit.region_id, hit.protein_id, hit.source, hit.profile)
        previous = best_by_key.get(key)
        if previous is None or _hit_score(hit) > _hit_score(previous):
            best_by_key[key] = hit
    return sorted(
        best_by_key.values(),
        key=lambda hit: (
            hit.region_id,
            hit.scaffold,
            hit.start,
            hit.end,
            hit.protein_id,
            hit.source,
            hit.profile,
        ),
    )


def _hit_score(hit: IntegrationGeneHit) -> tuple[float, float, float]:
    return (
        hit.domain_score if hit.domain_score is not None else float("-inf"),
        hit.sequence_score if hit.sequence_score is not None else float("-inf"),
        -hit.evalue if hit.evalue is not None else float("-inf"),
    )


def _optional_float(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


def _optional_int(value: str) -> int | None:
    try:
        return int(value)
    except ValueError:
        return None


def _decode(value: bytes | None) -> str:
    return value.decode() if value is not None else ""
