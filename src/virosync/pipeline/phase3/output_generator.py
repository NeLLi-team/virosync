"""Output Generator for EVE Predictions.

Generates standardized outputs compatible with downstream tools:
- GVClass-compatible FASTA and TSV files
- BED/GFF3 annotation files
- Evidence profiles (JSON/HDF5)
- Summary reports
"""

import csv
import json
import logging
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import cast
from urllib.parse import quote

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from virosync import __version__
from virosync.output_contract import (
    CONCRETE_EVE_CLASSES,
    DETAILED_PREDICTION_COLUMNS,
    DETAILED_PREDICTION_EXTENDED_COLUMNS,
    DETAILED_TAXONOMY_PARTITION,
    INTEGRATION_EVIDENCE_COLUMNS,
    canonical_family,
    coordinate_contract_metadata,
    normalize_effective_eve_class,
    resolve_effective_eve_class,
)
from virosync.utils.atomic_write import atomic_write_context
from virosync.utils.path_safety import require_strict_child, safe_filename_components

from .evidence_synthesizer import VerificationResult
from .gene_taxonomy import qualified_viral_hits

logger = logging.getLogger(__name__)


# Concrete viral families used for eve_class resolution/precedence. NOTE: the
# gate ALSO accepts "MIXED" (multi-family regions) as a first-class category via
# its own scoring branch; MIXED is deliberately kept out of this set so a
# concrete family always wins label resolution, but it is NOT disqualified.
_V2_EVE_CLASSES = CONCRETE_EVE_CLASSES

_GENE_TAXONOMY_COLUMNS = (
    "porf_id",
    "scaffold",
    "start",
    "end",
    "best_hit_origin",
    "best_hit_target",
    "best_hit_evalue",
    "top10_origins",
    "has_viral_neighbor",
    "has_ncldv_top10",
    "has_mirus_top10",
    "has_vp_plv_top10",
)
_COMBINED_GENE_TAXONOMY_COLUMNS = (
    "eve_id",
    "contig",
    "porf_id",
    "start",
    "end",
    "best_hit_origin",
    "best_hit_target",
    "best_hit_score",
    "top10_origins",
    "has_viral_neighbor",
    "has_ncldv_top10",
    "has_mirus_top10",
    "has_vp_plv_top10",
    "is_high_pident_euk",
    "is_flanking",
    "flank_position",
)


@dataclass(frozen=True, slots=True)
class _MarkerDisplayGroups:
    """Marker labels in legacy and protein-pattern groupings."""

    seed_legacy: tuple[str, ...]
    other_legacy: tuple[str, ...]
    seed_patterns: tuple[str, ...]
    other_patterns: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class _DetailedOutputSources:
    """Read-only source indexes used to construct detailed prediction rows."""

    porfs_by_scaffold: dict[str, list[tuple[int, int, str]]]
    rescued_protein_ids: set[str]
    diamond_flags: dict[str, tuple[bool, bool]]
    hmm_targets: dict[str, set[str]]
    marker_hits_by_scaffold: dict[str, list[tuple[int, int, str, str, str, float]]]
    genome_gc: float


@dataclass(frozen=True, slots=True)
class _GeneTaxonomyOutputRecord:
    """One normalized taxonomy record for both exported TSV schemas."""

    porf_id: str
    per_eve_scaffold: str
    combined_scaffold: str
    per_eve_start: str
    per_eve_end: str
    combined_start: str
    combined_end: str
    best_hit_origin: str
    best_hit_target: str
    best_hit_evalue: str
    best_hit_score: str
    top10_origins: str
    per_eve_has_viral_neighbor: bool
    combined_has_viral_neighbor: bool
    has_ncldv_top10: bool
    has_mirus_top10: bool
    per_eve_has_vp_plv_top10: bool
    combined_has_vp_plv_top10: bool
    is_high_pident_euk: bool
    is_flanking: bool
    flank_position: str
    legacy_per_eve_fields: tuple[str, ...] | None = None

    def per_eve_row(self) -> list[str]:
        """Return fields in the per-EVE taxonomy schema."""
        if self.legacy_per_eve_fields is not None:
            return list(self.legacy_per_eve_fields)
        return [
            self.porf_id,
            self.per_eve_scaffold,
            self.per_eve_start,
            self.per_eve_end,
            self.best_hit_origin,
            self.best_hit_target,
            self.best_hit_evalue,
            self.top10_origins,
            _tsv_flag(self.per_eve_has_viral_neighbor),
            _tsv_flag(self.has_ncldv_top10),
            _tsv_flag(self.has_mirus_top10),
            _tsv_flag(self.per_eve_has_vp_plv_top10),
        ]

    def combined_row(self, eve_id: str) -> list[str]:
        """Return fields in the aggregate taxonomy schema."""
        return [
            eve_id,
            self.combined_scaffold,
            self.porf_id,
            self.combined_start,
            self.combined_end,
            self.best_hit_origin,
            self.best_hit_target,
            self.best_hit_score,
            self.top10_origins,
            _tsv_flag(self.combined_has_viral_neighbor),
            _tsv_flag(self.has_ncldv_top10),
            _tsv_flag(self.has_mirus_top10),
            _tsv_flag(self.combined_has_vp_plv_top10),
            _tsv_flag(self.is_high_pident_euk),
            _tsv_flag(self.is_flanking),
            self.flank_position,
        ]


def _tsv_flag(value: object) -> str:
    """Serialize a truth value as the TSV contract's integer flag."""
    return "1" if value else "0"


def _integration_evidence_row(result: VerificationResult) -> dict[str, str]:
    """Format the same integration evidence for accepted and detailed tables."""
    recombinases = sorted(
        {str(hit["protein_id"]) for hit in result.integration_gene_hits if "recombinase" in str(hit["mechanism"])}
    )
    row = {
        "tir_present": _tsv_flag(result.tir_present),
        "tir_status": result.tir_status,
        "tir_candidate_count": str(result.tir_candidate_count),
        "tir_identity": f"{result.tir_identity:.4f}" if result.tir_present else ".",
        "tir_alignment_capped": _tsv_flag(result.tir_alignment_capped),
        "tir_alignment_length": str(result.tir_alignment_length) if result.tir_present else ".",
        "tir_boundary_override": _tsv_flag(result.tir_boundary_override),
        "tsd_sequence": result.tsd_sequence or ".",
        "recombinase_genes": "|".join(recombinases) or ".",
        "integration_gene_evidence": (
            json.dumps(result.integration_gene_hits, separators=(",", ":")) if result.integration_gene_hits else "."
        ),
    }
    for name in (
        "tir_scan_start",
        "tir_scan_end",
        "tir_left_start",
        "tir_left_end",
        "tir_right_start",
        "tir_right_end",
        "pre_tir_start",
        "pre_tir_end",
    ):
        value = getattr(result, name)
        row[name] = str(value) if value is not None else "."
    return row


def _per_eve_coordinate(record: Mapping[str, object], alias: str, canonical: str) -> object:
    """Apply the per-EVE alias fallback while retaining a numeric zero."""
    alias_value = record.get(alias)
    if alias_value == 0 and alias_value is not False:
        return alias_value
    return alias_value or record.get(canonical, "")


def _normalize_gene_taxonomy_record(record: Mapping[str, object]) -> _GeneTaxonomyOutputRecord:
    """Interpret a Phase 2 or legacy Phase 3 taxonomy record once."""
    if "top1_prefix" not in record:
        scaffold = cast(str, record.get("scaffold", "."))
        best_hit_evalue = str(record.get("best_hit_evalue", ""))
        return _GeneTaxonomyOutputRecord(
            porf_id=cast(str, record.get("porf_id", ".")),
            per_eve_scaffold=scaffold,
            combined_scaffold=scaffold,
            per_eve_start=str(record.get("start", "")),
            per_eve_end=str(record.get("end", "")),
            combined_start=str(record.get("start", "")),
            combined_end=str(record.get("end", "")),
            best_hit_origin=cast(str, record.get("best_hit_origin", ".")),
            best_hit_target=cast(str, record.get("best_hit_target", ".")),
            best_hit_evalue=best_hit_evalue,
            best_hit_score=best_hit_evalue,
            top10_origins=cast(str, record.get("top10_origins", "")),
            per_eve_has_viral_neighbor=bool(record.get("has_viral_neighbor")),
            combined_has_viral_neighbor=bool(record.get("has_viral_neighbor")),
            has_ncldv_top10=bool(record.get("has_ncldv_top10")),
            has_mirus_top10=bool(record.get("has_mirus_top10")),
            per_eve_has_vp_plv_top10=bool(record.get("has_vp_plv_top10")),
            combined_has_vp_plv_top10=bool(record.get("has_vp_plv_top10")),
            is_high_pident_euk=False,
            is_flanking=False,
            flank_position=".",
        )

    porf_id = cast(str, record.get("porf_id", "."))
    top10 = record.get("top10_prefixes", [])
    top10_origins = ",".join(top10) if isinstance(top10, list) else str(top10)
    top10_prefixes = top10_origins.split(",") if top10_origins else []
    has_mirus = "MIRUS" in top10_prefixes
    has_vp_plv = any(prefix in {"VP", "PLV", "PPV"} for prefix in top10_prefixes)
    flank_position = record.get("flank_position", "")
    legacy_per_eve_fields = None
    if "start" not in record and "porf_start" not in record:
        legacy_per_eve_fields = (
            porf_id,
            cast(str, record.get("scaffold", ".")),
            str(record.get("start", "")),
            str(record.get("end", "")),
            cast(str, record.get("best_hit_origin", ".")),
            cast(str, record.get("best_hit_target", ".")),
            str(record.get("best_hit_evalue", "")),
            cast(str, record.get("top10_origins", "")),
            _tsv_flag(record.get("has_viral_neighbor")),
            _tsv_flag(record.get("has_ncldv_top10")),
            _tsv_flag(record.get("has_mirus_top10")),
            _tsv_flag(record.get("has_vp_plv_top10")),
        )
    return _GeneTaxonomyOutputRecord(
        porf_id=porf_id,
        per_eve_scaffold=cast(
            str,
            record.get("scaffold") or (porf_id.split("|", 1)[0] if "|" in porf_id else porf_id),
        ),
        combined_scaffold=cast(str, record.get("scaffold", ".")),
        per_eve_start=str(_per_eve_coordinate(record, "porf_start", "start")),
        per_eve_end=str(_per_eve_coordinate(record, "porf_end", "end")),
        combined_start=str(record.get("porf_start", record.get("start", ""))),
        combined_end=str(record.get("porf_end", record.get("end", ""))),
        best_hit_origin=cast(str, record.get("top1_prefix", ".")),
        best_hit_target=cast(str, record.get("top1_target", ".")),
        best_hit_evalue=str(record.get("top1_evalue", "")),
        best_hit_score=str(record.get("top1_pident", "")),
        top10_origins=top10_origins,
        per_eve_has_viral_neighbor=bool(record.get("has_ncldv_mirus")),
        combined_has_viral_neighbor=bool(record.get("has_viral")),
        has_ncldv_top10=bool(record.get("has_ncldv_mirus")),
        has_mirus_top10=has_mirus,
        per_eve_has_vp_plv_top10=bool(record.get("has_vp_plv") or has_vp_plv) if top10_origins else False,
        combined_has_vp_plv_top10=has_vp_plv,
        is_high_pident_euk=bool(record.get("is_high_pident_euk")),
        is_flanking=bool(record.get("is_flanking", False)),
        flank_position=cast(str, flank_position if flank_position else "."),
        legacy_per_eve_fields=legacy_per_eve_fields,
    )


def _taxonomy_prefixes(value: object) -> list[str]:
    if isinstance(value, str):
        raw_prefixes = value.split(",")
    elif isinstance(value, (list, tuple)):
        raw_prefixes = value
    else:
        raw_prefixes = []
    return [canonical_family(str(prefix).rstrip("_")) for prefix in raw_prefixes if str(prefix).strip()]


def _qualified_viral_prefixes(record: dict) -> list[str]:
    """Reference namespaces of a record's identity-qualified viral top-10 hits.

    The taxonomy_best_hits partition reports GVMAG and PHAGE as their own
    buckets, so it uses the raw namespaces rather than the published classes
    that ``viral_hit_categories`` folds them onto.
    """
    return [prefix for prefix, _target in qualified_viral_hits(record)]


def _taxonomy_partition_bucket(record: dict) -> str:
    top1_target = str(record.get("top1_target") or "")
    if not top1_target or top1_target in {".", "NA", "None"}:
        return "NO_HITS"

    qualified_viral = _qualified_viral_prefixes(record)
    if qualified_viral:
        return qualified_viral[0]

    raw_top1 = top1_target.split("__", 1)[0].upper() if "__" in top1_target else "UNKNOWN"
    if raw_top1 in {"MITO", "PLASTID"}:
        return raw_top1
    top1 = canonical_family(record.get("top1_prefix") or "UNKNOWN")
    if top1 in {"EUK", "BAC", "ARC"}:
        return top1
    return "UNK"


def _gff3_escape(value: object) -> str:
    """Percent-encode a raw biological identifier for a GFF3 field."""
    return quote(
        str(value),
        safe=("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._:-|"),
    )


@dataclass(frozen=True)
class QualityGateDecision:
    """Decision made by the canonical v2 output quality gate."""

    kept: bool
    effective_class: str
    reason: str
    promoted_low: bool = False


def _resolve_eve_class(result: VerificationResult) -> str:
    """Resolve the HIGH/MEDIUM ``eve_class`` for the v2 quality gate.

    Returns the effective EVE class for v2 filter evaluation, preferring a
    concrete family from ``region_classification`` (authoritative when it is a
    known class), then ``classification`` / ``likely_family``. A ``MIXED``
    region (multiple viral families seeded together) is surfaced as
    ``"MIXED"`` so the gate can score it under the normal viral rule instead of
    dropping it; anything else becomes ``"UNKNOWN"`` and is rejected.
    """
    return resolve_effective_eve_class(
        region_classification=getattr(result, "region_classification", ""),
        classification=getattr(result, "classification", ""),
        likely_family=getattr(result, "likely_family", ""),
    )


def _published_eve_class(result: VerificationResult) -> str:
    """Return the PUBLISHED class for one result.

    Publication reports the taxonomy consensus over the validated markers' own
    top-10 hits, not the gate's acceptance label: the gate's vocabulary still
    carries MIXED and drives acceptance only.
    """
    return normalize_effective_eve_class(getattr(result, "taxonomy_class", "UNKNOWN"))


_ATPASE_HALLMARK_NAMES = {"plv_pc_054", "gvogm0760"}


def _is_atpase_marker(name: str) -> bool:
    """True if a hallmark marker is a packaging ATPase (A32 / FtsK-HerA family).

    These cross-hit ubiquitous cellular P-loop NTPases, so an ATPase-only PLV/VP
    region is unreliable and must not, on its own, support an accepted call.
    """
    n = (name or "").lower()
    return n in _ATPASE_HALLMARK_NAMES or "atpase" in n


def evaluate_v2_quality_gate(result: VerificationResult) -> QualityGateDecision:
    """Evaluate one result against the canonical v2 quality gate."""
    tier = (getattr(result, "confidence_tier", "") or "").upper()
    length = max(0, int(getattr(result, "end", 0)) - int(getattr(result, "start", 0)))
    hallmark = int(getattr(result, "hallmark_count", 0) or 0)
    has_mcp = bool(getattr(result, "has_mcp", False))
    # Non-ATPase hallmark count: a PLV/VP region supported only by the broad packaging
    # ATPase (which cross-hits cellular NTPases) is unreliable. Require MCP or >=1
    # non-ATPase hallmark for acceptance; ATPase-only stays discovery-only (gated out).
    hallmark_genes = getattr(result, "hallmark_genes", []) or []
    non_atpase_hallmark = sum(1 for g in hallmark_genes if not _is_atpase_marker(g))
    eve_class = resolve_effective_eve_class(
        confidence_tier=tier,
        region_classification=getattr(result, "region_classification", ""),
        classification=getattr(result, "classification", ""),
        likely_family=getattr(result, "likely_family", ""),
    )

    if tier in ("HIGH", "MEDIUM"):
        if eve_class == "MIXED":
            # MIXED = multiple viral families seeded in one region, the expected
            # signature of NCLDV-adjacent capscan PLVs (Aquintoviricetes "Near-"
            # groups, "NCV-like" groups) that carry NCLDV-family hallmark hits
            # alongside their own capsid, so no single family wins the classifier
            # tie-break. MIXED is a first-class accepted category, scored under
            # the same rule as PLV/VP/PPV: an MCP is strong evidence but NOT
            # required -- >=2 hallmarks with >=1 non-ATPase also qualifies. A
            # high-scoring MIXED region is kept, not disqualified.
            kept = length > 2000 and (has_mcp or (hallmark >= 2 and non_atpase_hallmark >= 1))
            reason = "mixed_high_medium_pass" if kept else "mixed_high_medium_gate"
            return QualityGateDecision(kept, eve_class, reason)
        if eve_class not in _V2_EVE_CLASSES:
            return QualityGateDecision(False, eve_class, "unsupported_class")
        if eve_class == "CRESS":
            identity_marker_support = "CRESS" in (getattr(result, "marker_family_hits", []) or [])
            kept = identity_marker_support and length > 0 and hallmark >= 1
            reason = "cress_identity_high_medium_pass" if kept else "cress_identity_required"
            return QualityGateDecision(kept, eve_class, reason)
        if eve_class == "PPV":
            kept = length > 2000 and (has_mcp or (hallmark >= 2 and non_atpase_hallmark >= 1))
            reason = "small_dna_high_medium_pass" if kept else "small_dna_high_medium_gate"
            return QualityGateDecision(kept, eve_class, reason)
        if eve_class in ("NCLDV", "MIRUS"):
            kept = length > 5000 or has_mcp
            reason = "ncldv_mirus_high_medium_pass" if kept else "ncldv_mirus_high_medium_gate"
            return QualityGateDecision(kept, eve_class, reason)

    if tier == "LOW":
        family = eve_class
        if family in ("NCLDV", "MIRUS"):
            kept = length > 5000 and hallmark >= 2
            reason = "ncldv_mirus_low_promoted" if kept else "ncldv_mirus_low_gate"
            return QualityGateDecision(kept, family, reason, promoted_low=kept)
        # LOW must never be LOOSER than the region's own HIGH/MEDIUM rule, or
        # raising a region's confidence would remove it from the output. These
        # two branches used to omit the `hallmark >= 2` conjunct, so a region
        # with one non-ATPase hallmark and no MCP was published at LOW and
        # dropped at MEDIUM. Unlike NCLDV/MIRUS above, the bare-MCP shortcut is
        # kept here: an MCP is the primary diagnostic for Preplasmiviricota, and
        # dropping it would reject 58% of published LOW PPV calls, which is a
        # sensitivity change rather than a correctness fix.
        if family == "CRESS":
            identity_marker_support = "CRESS" in (getattr(result, "marker_family_hits", []) or [])
            kept = identity_marker_support and length > 0 and hallmark >= 1
            reason = "cress_identity_low_promoted" if kept else "cress_identity_required"
            return QualityGateDecision(
                kept,
                family,
                reason,
                promoted_low=kept,
            )
        if family == "PPV":
            kept = length > 2000 and (has_mcp or (hallmark >= 2 and non_atpase_hallmark >= 1))
            reason = "small_dna_low_promoted" if kept else "small_dna_low_gate"
            return QualityGateDecision(kept, family, reason, promoted_low=kept)
        if family == "MIXED":
            kept = length > 2000 and (has_mcp or (hallmark >= 2 and non_atpase_hallmark >= 1))
            reason = "mixed_low_promoted" if kept else "mixed_low_gate"
            return QualityGateDecision(kept, family, reason, promoted_low=kept)
        return QualityGateDecision(False, family, "low_unsupported_family")

    return QualityGateDecision(False, eve_class, "unsupported_tier")


class OutputGenerator:
    """Generates all output files for verified EVE predictions."""

    def __init__(
        self,
        output_dir: Path,
        genome_fasta: Path | None = None,
        proteome_fasta: Path | None = None,
        extended_output: bool = True,
        seed_marker_allowlist: list[str] | None = None,
        export_all_eve_sequences: bool = False,
    ):
        """Initialize output generator.

        Args:
            output_dir: Base output directory
            genome_fasta: Path to genome FASTA (for extracting sequences)
            proteome_fasta: Path to proteome FASTA (for extracting proteins)
        """
        self.output_dir = Path(output_dir)
        self.genome_fasta = genome_fasta
        self.proteome_fasta = proteome_fasta
        self.extended_output = extended_output
        self.seed_marker_allowlist = seed_marker_allowlist or []
        self.export_all_eve_sequences = export_all_eve_sequences

        # Load sequences if provided
        self._genome_sequences = None
        self._proteome_sequences = None
        self._marker_hits = None
        self._protein_counts_cache = None
        self._protein_to_models_cache = None  # NEW: Cache for protein-to-models mapping
        self._porfs_by_scaffold = None
        self._rescued_protein_sequences = None
        self._rescued_protein_descriptions = None

    @staticmethod
    def _eve_filename_components(
        results: list[VerificationResult],
    ) -> dict[str, str]:
        """Map raw EVE IDs to distinct filesystem-safe filename components."""
        return safe_filename_components(
            (result.eve_id for result in results),
            label="EVE ID",
        )

    @property
    def genome_sequences(self) -> dict[str, str]:
        """Lazy load genome sequences."""
        if self._genome_sequences is None and self.genome_fasta:
            self._genome_sequences = {}
            for record in SeqIO.parse(self.genome_fasta, "fasta"):
                self._genome_sequences[record.id] = str(record.seq)
        return self._genome_sequences or {}

    @property
    def proteome_sequences(self) -> dict[str, str]:
        """Lazy load proteome sequences."""
        if self._proteome_sequences is None and self.proteome_fasta:
            self._proteome_sequences = {}
            for record in SeqIO.parse(self.proteome_fasta, "fasta"):
                self._proteome_sequences[record.id] = str(record.seq)
        return self._proteome_sequences or {}

    def _find_file(self, relative_paths: list[Path]) -> Path | None:
        for rel in relative_paths:
            candidate = (self.output_dir / rel).resolve()
            if candidate.exists():
                return candidate
        return None

    def _find_diamond_results(self) -> Path | None:
        return self._find_file(
            [
                Path("phase1/marker_validation/diamond_top10_taxonomy.tsv"),
                Path("../phase1/marker_validation/diamond_top10_taxonomy.tsv"),
                Path("../../phase1/marker_validation/diamond_top10_taxonomy.tsv"),
                Path("phase1/novelty/diamond/diamond_combined.tsv"),
                Path("../phase1/novelty/diamond/diamond_combined.tsv"),
                Path("../../phase1/novelty/diamond/diamond_combined.tsv"),
            ]
        )

    @staticmethod
    def _base_porf_id(porf_id: str) -> str:
        """Normalize pORF IDs by stripping optional domain suffixes."""
        return porf_id.split("|aa", 1)[0] if "|aa" in porf_id else porf_id

    @staticmethod
    def _is_seed_marker_name(marker_name: str) -> bool:
        """Heuristic fallback for seed marker classification.

        Used when no explicit seed-marker allowlist is provided.
        """
        key = marker_name.lower()
        if key in {"og1352", "og484"}:
            return True
        return key.startswith(("gvogm", "gamadvirusmcp", "plv_", "vp_", "mirus_"))

    @staticmethod
    def _normalize_model_name(marker_name: str) -> str:
        """Normalize marker IDs to stable display names."""
        marker = (marker_name or "").strip()
        upper = marker.upper()
        if upper.startswith("GVOGM"):
            return f"GVOGm{upper[5:]}"
        if upper.startswith("OG"):
            return f"OG{upper[2:]}"
        return marker

    def _ensure_protein_to_models_cache(
        self,
        eve_regions: list[tuple[str, int, int]],
    ) -> dict[tuple[str, int, int], dict[str, set[str]]]:
        """Populate region-level protein-to-model mapping when missing."""
        if self._protein_to_models_cache is None:
            hits_by_scaffold, protein_to_models_global = self._parse_validated_marker_hits()
            self._protein_to_models_cache = self._build_region_protein_mapping(
                eve_regions, hits_by_scaffold, protein_to_models_global
            )
        return self._protein_to_models_cache

    def _protein_model_summary(
        self,
        protein_to_models: dict[str, set[str]],
        prefix: str | None = None,
    ) -> tuple[int, Counter]:
        """Summarize model support as unique-protein totals plus per-model counts."""
        total = 0
        counts: Counter[str] = Counter()
        normalized_prefix = (prefix or "").upper()

        for models in protein_to_models.values():
            matching = {
                self._normalize_model_name(model)
                for model in models
                if not normalized_prefix or model.upper().startswith(normalized_prefix)
            }
            if not matching:
                continue
            total += 1
            for model in matching:
                counts[model] += 1

        return total, counts

    @staticmethod
    def _format_counted_names(model_counts: Counter[str]) -> str:
        """Format per-model protein counts as MODEL:n tokens."""
        if not model_counts:
            return "."
        items = sorted(model_counts.items(), key=lambda item: item[0])
        return ",".join(f"{model}:{count}" for model, count in items)

    def _group_marker_names_for_display(
        self,
        marker_names: list[str],
        scaffold: str,
        start: int,
        end: int,
        use_protein_patterns: bool = False,
    ) -> list[str]:
        """Group markers by protein hit patterns (if enabled and available).

        Args:
            marker_names: Marker names to group
            scaffold, start, end: EVE region
            use_protein_patterns: If True, use protein-pattern grouping

        Returns:
            List of formatted marker strings
        """
        if not marker_names:
            return []

        # Legacy grouping (always available as fallback)
        if not use_protein_patterns or self._protein_to_models_cache is None:
            return self._group_marker_names_legacy(marker_names, scaffold, start, end)

        region_key = (scaffold, start, end)
        protein_to_models = self._protein_to_models_cache.get(region_key, {})

        if not protein_to_models:
            return self._group_marker_names_legacy(marker_names, scaffold, start, end)

        # Build reverse mapping: frozenset(models) -> protein_count
        pattern_counts = Counter()

        # Normalize case for comparison
        marker_set = {m.upper() for m in marker_names}
        mapped_markers = set()

        for protein, models in protein_to_models.items():
            relevant_models = models & marker_set
            if relevant_models:
                pattern_counts[frozenset(relevant_models)] += 1
                mapped_markers |= relevant_models

        # Check for unmapped markers
        unmapped = marker_set - mapped_markers

        # Format mapped patterns
        result = []
        sorted_patterns = sorted(pattern_counts.items(), key=lambda x: (-x[1], -len(x[0]), sorted(x[0])[0]))

        for model_pattern, protein_count in sorted_patterns:
            sorted_models = sorted(model_pattern)
            model_str = "/".join(sorted_models)
            result.append(f"{model_str}:{protein_count}")

        # Add unmapped markers using legacy grouping (safety net)
        if unmapped:
            logger.debug(f"Region {scaffold}:{start}-{end} has {len(unmapped)} unmapped markers")
            # sorted(), not list(): set iteration order varies between runs under
            # string hash randomization and made this column non-deterministic.
            legacy_unmapped = self._group_marker_names_legacy(sorted(unmapped), scaffold, start, end)
            result.extend(legacy_unmapped)

        return result

    def _group_marker_names_legacy(
        self,
        marker_names: list[str],
        scaffold: str,
        start: int,
        end: int,
    ) -> list[str]:
        """Legacy functional grouping (fallback).

        Preserves exact existing behavior for backward compatibility.
        Groups PLV_MCP_1-10 as "PLV_MCP(n proteins)" where n is unique protein count.
        Falls back to hit counts if protein data unavailable.
        """
        if not marker_names:
            return []

        # Group markers by functional category
        plv_mcp_hits = []
        vp_mcp_hits = []
        vp_atpase_hits = []
        vp_penton_hits = []
        other_markers = []

        for marker in marker_names:
            key = marker.lower()
            if key.startswith("plv_mcp"):
                plv_mcp_hits.append(marker)
            elif key.startswith("vp_mcp"):
                vp_mcp_hits.append(marker)
            elif key.startswith("vp_atpase"):
                vp_atpase_hits.append(marker)
            elif key.startswith("vp_penton"):
                vp_penton_hits.append(marker)
            else:
                other_markers.append(marker)

        # Build deduplicated list
        result = []

        # Try to get protein counts from cached Phase 1 data
        region_key = (scaffold, start, end)
        protein_counts = {}
        if self._protein_counts_cache is not None and region_key in self._protein_counts_cache:
            protein_counts = self._protein_counts_cache[region_key]

        # Add grouped markers with protein counts (or hit counts as fallback)
        if plv_mcp_hits:
            if "PLV_MCP" in protein_counts:
                count = protein_counts["PLV_MCP"]
                if count == 1:
                    result.append("PLV_MCP(1 protein)")
                else:
                    result.append(f"PLV_MCP({count} proteins)")
            else:
                count = len(plv_mcp_hits)
                result.append(f"PLV_MCP({count}x)")

        if vp_mcp_hits:
            if "VP_MCP" in protein_counts:
                count = protein_counts["VP_MCP"]
                if count == 1:
                    result.append("VP_MCP(1 protein)")
                else:
                    result.append(f"VP_MCP({count} proteins)")
            else:
                count = len(vp_mcp_hits)
                result.append(f"VP_MCP({count}x)")

        if vp_atpase_hits:
            if "VP_ATPase" in protein_counts:
                count = protein_counts["VP_ATPase"]
                if count == 1:
                    result.append("VP_ATPase(1 protein)")
                else:
                    result.append(f"VP_ATPase({count} proteins)")
            else:
                count = len(vp_atpase_hits)
                result.append(f"VP_ATPase({count}x)")

        if vp_penton_hits:
            if "VP_Penton" in protein_counts:
                count = protein_counts["VP_Penton"]
                if count == 1:
                    result.append("VP_Penton(1 protein)")
                else:
                    result.append(f"VP_Penton({count} proteins)")
            else:
                count = len(vp_penton_hits)
                result.append(f"VP_Penton({count}x)")

        # Add other markers (deduplicated with counts)
        other_counts = Counter(other_markers)
        for marker in sorted(set(other_markers), key=other_markers.index):
            count = other_counts[marker]
            if count > 1:
                result.append(f"{marker}({count}x)")
            else:
                result.append(marker)

        return result

    def _marker_display_groups(self, result: VerificationResult) -> _MarkerDisplayGroups:
        """Build both marker label groupings for one prediction row."""
        marker_names = self._marker_names_for_region(result.scaffold, result.start, result.end)
        if self.seed_marker_allowlist:
            allowlist = {marker.lower() for marker in self.seed_marker_allowlist}
            seed_markers = [marker for marker in marker_names if marker.lower() in allowlist]
        else:
            seed_markers = [marker for marker in marker_names if self._is_seed_marker_name(marker)]
        other_markers = [marker for marker in marker_names if marker not in seed_markers]

        return _MarkerDisplayGroups(
            seed_legacy=tuple(
                self._group_marker_names_for_display(
                    seed_markers,
                    result.scaffold,
                    result.start,
                    result.end,
                    use_protein_patterns=False,
                )
            ),
            other_legacy=tuple(
                self._group_marker_names_for_display(
                    other_markers,
                    result.scaffold,
                    result.start,
                    result.end,
                    use_protein_patterns=False,
                )
            ),
            seed_patterns=tuple(
                self._group_marker_names_for_display(
                    seed_markers,
                    result.scaffold,
                    result.start,
                    result.end,
                    use_protein_patterns=True,
                )
            ),
            other_patterns=tuple(
                self._group_marker_names_for_display(
                    other_markers,
                    result.scaffold,
                    result.start,
                    result.end,
                    use_protein_patterns=True,
                )
            ),
        )

    def _load_protein_counts_by_region(
        self, eve_regions: list[tuple[str, int, int]]
    ) -> dict[tuple[str, int, int], dict[str, int]]:
        """Count unique proteins per EVE region per marker group from Phase 1 data.

        Args:
            eve_regions: List of (scaffold, start, end) tuples for EVE regions

        Returns:
            Dict mapping (scaffold, start, end) -> {marker_group: protein_count}
        """
        if self._protein_counts_cache is not None:
            return self._protein_counts_cache

        validated_hits_paths = self._validated_marker_hit_paths()
        if not validated_hits_paths:
            logger.warning(
                "Phase 1 marker-hit tables not found - marker display will show hit counts instead of protein counts"
            )
            self._protein_counts_cache = {}
            return self._protein_counts_cache

        # First collect all hits by interval
        hits_by_interval = []  # [(scaffold, start, end, base_protein, marker_group)]

        for validated_hits_path in validated_hits_paths:
            try:
                with open(validated_hits_path) as f:
                    # Skip header
                    next(f)
                    for line in f:
                        fields = line.strip().split("\t")
                        if len(fields) < 6:
                            continue

                        query_porf = fields[0]
                        scaffold = fields[1]
                        start = int(fields[2])
                        end = int(fields[3])
                        marker = fields[5]

                        if "|aa" in query_porf:
                            base_protein = query_porf.rsplit("|aa", 1)[0]
                        else:
                            base_protein = query_porf

                        marker_lower = marker.lower()
                        if marker_lower.startswith("plv_mcp"):
                            marker_group = "PLV_MCP"
                        elif marker_lower.startswith("vp_mcp"):
                            marker_group = "VP_MCP"
                        elif marker_lower.startswith("vp_atpase"):
                            marker_group = "VP_ATPase"
                        elif marker_lower.startswith("vp_penton"):
                            marker_group = "VP_Penton"
                        else:
                            continue

                        hits_by_interval.append((scaffold, start, end, base_protein, marker_group))
            except (OSError, StopIteration, ValueError) as error:
                logger.warning(
                    "Error loading protein counts from %s: %s",
                    validated_hits_path,
                    error,
                )

        # Now assign hits to EVE regions
        protein_counts = {}

        for eve_scaffold, eve_start, eve_end in eve_regions:
            region_key = (eve_scaffold, eve_start, eve_end)
            protein_counts[region_key] = {}

            for hit_scaffold, hit_start, hit_end, base_protein, marker_group in hits_by_interval:
                if hit_scaffold == eve_scaffold and hit_start < eve_end and hit_end > eve_start:
                    if marker_group not in protein_counts[region_key]:
                        protein_counts[region_key][marker_group] = set()
                    protein_counts[region_key][marker_group].add(base_protein)

        # Convert sets to counts
        for region_key in protein_counts:
            for marker_group in protein_counts[region_key]:
                protein_counts[region_key][marker_group] = len(protein_counts[region_key][marker_group])

        self._protein_counts_cache = protein_counts
        logger.info(f"Loaded protein counts for {len(protein_counts)} EVE regions from Phase 1 data")
        return protein_counts

    def _parse_validated_marker_hits(
        self,
    ) -> tuple[dict[str, list[tuple]], dict[str, set[str]]]:
        """Parse validated_marker_hits.tsv once for both:
        1. Hits by scaffold for overlap checking
        2. Protein-to-models mapping

        Returns:
            (hits_by_scaffold, base_protein_to_models)

        hits_by_scaffold: {scaffold: [(start, end, marker, validation_status, base_protein)]}
        base_protein_to_models: {base_protein_id: {model1, model2, ...}}
        """
        import csv
        from collections import defaultdict

        validated_hits_paths = self._validated_marker_hit_paths()
        if not validated_hits_paths:
            logger.warning("validated_marker_hits.tsv not found - protein-pattern grouping unavailable")
            return {}, {}

        hits_by_scaffold = defaultdict(list)
        protein_to_models = defaultdict(set)

        for validated_hits_path in validated_hits_paths:
            file_hits = defaultdict(list)
            file_models = defaultdict(set)
            try:
                with open(validated_hits_path) as f:
                    reader = csv.DictReader(f, delimiter="\t")
                    for row in reader:
                        query_porf = row["query_porf"]
                        scaffold = row["scaffold"]
                        start = int(row["start"])
                        end = int(row["end"])
                        marker = row["hmm_target"]
                        status = row.get("validation_status", "")

                        # Only include validated markers
                        if status not in ("validated", "validated_novel"):
                            continue

                        # Extract base protein ID
                        base_protein = query_porf.rsplit("|aa", 1)[0] if "|aa" in query_porf else query_porf

                        # Index by scaffold for fast lookup
                        file_hits[scaffold].append((start, end, marker, status, base_protein))

                        # Build protein-to-models mapping
                        file_models[base_protein].add(self._normalize_model_name(marker))
            except (KeyError, OSError, TypeError, ValueError) as error:
                logger.error(
                    "Error parsing marker hits from %s: %s",
                    validated_hits_path,
                    error,
                )
                continue
            for scaffold, hits in file_hits.items():
                hits_by_scaffold[scaffold].extend(hits)
            for protein, models in file_models.items():
                protein_to_models[protein].update(models)

        # Sort hits by start position for efficient overlap checking
        for scaffold in hits_by_scaffold:
            hits_by_scaffold[scaffold].sort(key=lambda x: x[0])

        logger.info(
            "Parsed %d validated marker hits from %d files",
            sum(len(v) for v in hits_by_scaffold.values()),
            len(validated_hits_paths),
        )
        return dict(hits_by_scaffold), dict(protein_to_models)

    def _build_region_protein_mapping(
        self,
        eve_regions: list[tuple[str, int, int]],
        hits_by_scaffold: dict[str, list[tuple]],
        protein_to_models: dict[str, set[str]],
    ) -> dict[tuple[str, int, int], dict[str, set[str]]]:
        """Assign proteins to EVE regions using proper overlap logic.

        Args:
            eve_regions: [(scaffold, start, end), ...]
            hits_by_scaffold: From _parse_validated_marker_hits
            protein_to_models: From _parse_validated_marker_hits

        Returns:
            {(scaffold, start, end): {base_protein: {models}}}
        """
        from collections import defaultdict

        region_mapping = {}

        for eve_scaffold, eve_start, eve_end in eve_regions:
            region_key = (eve_scaffold, eve_start, eve_end)
            proteins_in_region = defaultdict(set)

            # Only check hits on the same scaffold
            scaffold_hits = hits_by_scaffold.get(eve_scaffold, [])

            for hit_start, hit_end, marker, status, base_protein in scaffold_hits:
                # Proper overlap check (from existing code)
                if hit_start < eve_end and hit_end > eve_start:
                    proteins_in_region[base_protein] |= protein_to_models.get(base_protein, set())

            region_mapping[region_key] = dict(proteins_in_region)

        logger.info(f"Built protein-to-models mapping for {len(region_mapping)} EVE regions")
        return region_mapping

    def _find_hmm_hits(self) -> Path | None:
        return self._find_file(
            [
                Path("phase1/marker_validation/validated_marker_hits.tsv"),
                Path("../phase1/marker_validation/validated_marker_hits.tsv"),
                Path("../../phase1/marker_validation/validated_marker_hits.tsv"),
                Path("phase1/hhg/hmm_hits_validated.tsv"),
                Path("../phase1/hhg/hmm_hits_validated.tsv"),
                Path("../../phase1/hhg/hmm_hits_validated.tsv"),
                Path("hmm_hits_validated.tsv"),
            ]
        )

    def _find_validated_marker_hits(self) -> Path | None:
        return self._find_file(
            [
                Path("phase1/marker_validation/validated_marker_hits.tsv"),
                Path("../phase1/marker_validation/validated_marker_hits.tsv"),
                Path("../../phase1/marker_validation/validated_marker_hits.tsv"),
            ]
        )

    def _find_confirmed_frameshift_marker_hits(self) -> Path | None:
        return self._find_file(
            [
                Path("phase1/frameshift_screening/confirmed_frameshift_markers.tsv"),
                Path("../phase1/frameshift_screening/confirmed_frameshift_markers.tsv"),
                Path("../../phase1/frameshift_screening/confirmed_frameshift_markers.tsv"),
            ]
        )

    def _validated_marker_hit_paths(self) -> list[Path]:
        return [
            path
            for path in (
                self._find_validated_marker_hits(),
                self._find_confirmed_frameshift_marker_hits(),
            )
            if path is not None
        ]

    def _find_confirmed_frameshift_proteins(self) -> Path | None:
        return self._find_file(
            [
                Path("phase1/frameshift_screening/confirmed_frameshift_proteins.faa"),
                Path("../phase1/frameshift_screening/confirmed_frameshift_proteins.faa"),
                Path("../../phase1/frameshift_screening/confirmed_frameshift_proteins.faa"),
            ]
        )

    def _load_marker_hits(self) -> dict[str, list[tuple[int, int, str, str, str, float]]]:
        if self._marker_hits is not None:
            return self._marker_hits

        hits_by_scaffold: dict[str, list[tuple[int, int, str, str, str, float]]] = {}
        hits_paths = self._validated_marker_hit_paths()
        if not hits_paths:
            self._marker_hits = hits_by_scaffold
            return hits_by_scaffold

        for hits_path in hits_paths:
            with open(hits_path) as handle:
                header = handle.readline().rstrip("\n").split("\t")
                idx = {name: i for i, name in enumerate(header)}
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 7:
                        continue
                    scaffold = parts[idx.get("scaffold", 1)]
                    try:
                        start = int(parts[idx.get("start", 2)])
                        end = int(parts[idx.get("end", 3)])
                    except ValueError:
                        continue
                    target = parts[idx.get("hmm_target", 5)]
                    status = parts[idx.get("validation_status", 7)]
                    porf_id = parts[idx.get("query_porf", 0)]
                    try:
                        hmm_score = float(parts[idx.get("hmm_score", 6)])
                    except (ValueError, IndexError):
                        hmm_score = 0.0
                    hits_by_scaffold.setdefault(scaffold, []).append((start, end, target, status, porf_id, hmm_score))

        for scaffold in hits_by_scaffold:
            hits_by_scaffold[scaffold].sort(key=lambda x: x[0])
        self._marker_hits = hits_by_scaffold
        return hits_by_scaffold

    def _marker_names_for_region(
        self,
        scaffold: str,
        start: int,
        end: int,
        status_filter: str | None = None,
    ) -> list[str]:
        """Return deduplicated marker names for a region.

        When multiple HMM profiles hit the same protein, only the
        best-scoring model (highest hmm_score) is returned.
        """
        hits_by_scaffold = self._load_marker_hits()
        hits = hits_by_scaffold.get(scaffold, [])
        # Collect hits overlapping the region, then deduplicate by base gene
        by_gene: dict[str, tuple[str, float]] = {}  # base_gene -> (target, score)
        for h_start, h_end, target, status, porf_id, hmm_score in hits:
            if h_start < end and h_end > start:
                if status_filter == "validated" and status not in ("validated", "validated_novel"):
                    continue
                if status_filter == "unvalidated" and status in ("validated", "validated_novel"):
                    continue
                base_gene = porf_id.split("|aa")[0] if porf_id else target
                if base_gene not in by_gene or hmm_score > by_gene[base_gene][1]:
                    by_gene[base_gene] = (target, hmm_score)
        return [target for target, _ in by_gene.values()]

    def _load_porfs_by_scaffold(self) -> dict[str, list[tuple[int, int, str]]]:
        from virosync.pipeline.phase0.prodigal import parse_prodigal_header

        if self._porfs_by_scaffold is not None:
            return self._porfs_by_scaffold

        porfs: dict[str, list[tuple[int, int, str]]] = {}
        rescued_sequences: dict[str, str] = {}
        rescued_descriptions: dict[str, str] = {}
        sources: list[tuple[Path, bool]] = []
        if self.proteome_fasta:
            sources.append((Path(self.proteome_fasta), False))
        rescued_fasta = self._find_confirmed_frameshift_proteins()
        if rescued_fasta:
            sources.append((rescued_fasta, True))
        for source, is_rescued in sources:
            for record in SeqIO.parse(source, "fasta"):
                parsed = parse_prodigal_header(record.description, record.id)
                if not parsed:
                    continue
                scaffold, start, end, _strand = parsed
                porfs.setdefault(scaffold, []).append((start, end, record.id))
                if is_rescued:
                    rescued_sequences[record.id] = str(record.seq)
                    rescued_descriptions[record.id] = record.description[len(record.id) :].strip()
        for scaffold in porfs:
            porfs[scaffold].sort()
        self._porfs_by_scaffold = porfs
        self._rescued_protein_sequences = rescued_sequences
        self._rescued_protein_descriptions = rescued_descriptions
        return porfs

    def _protein_records_for_region(
        self,
        scaffold: str,
        start: int,
        end: int,
    ) -> list[SeqRecord]:
        """Return protein records whose coordinates overlap a region."""
        eve_proteins = []
        porfs_by_scaffold = self._load_porfs_by_scaffold()
        for porf_start, porf_end, porf_id in porfs_by_scaffold.get(scaffold, []):
            if porf_start >= end or porf_end <= start:
                continue
            seq = self.proteome_sequences.get(porf_id)
            description = ""
            if not seq and self._rescued_protein_sequences:
                seq = self._rescued_protein_sequences.get(porf_id)
                description = (self._rescued_protein_descriptions or {}).get(
                    porf_id,
                    "frameshift_rescued_domain",
                )
            if not seq:
                continue
            eve_proteins.append(
                SeqRecord(
                    Seq(seq),
                    id=porf_id,
                    description=description,
                )
            )
        return eve_proteins

    def _load_diamond_top10_flags(self) -> dict[str, tuple[bool, bool]]:
        """Return query -> (has_ncldv, has_mirus) based on top-10 Diamond hits."""
        diamond_path = self._find_diamond_results()
        if not diamond_path or not diamond_path.exists():
            return {}

        if diamond_path.name == "diamond_top10_taxonomy.tsv":
            flags: dict[str, tuple[bool, bool]] = {}
            with open(diamond_path) as handle:
                header = handle.readline()
                if not header:
                    return flags
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 7:
                        continue
                    query = parts[0]
                    base_query = self._base_porf_id(query)
                    has_ncldv = parts[4].strip() == "1"
                    has_mirus = parts[5].strip() == "1"
                    value = (has_ncldv, has_mirus)
                    flags[query] = value
                    flags[base_query] = value
            return flags

        top10: dict[str, list[tuple[float, str]]] = {}
        with open(diamond_path) as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                query, target, _, bits = parts[0], parts[1], parts[2], parts[3]
                try:
                    bits_val = float(bits)
                except ValueError:
                    bits_val = 0.0
                hits = top10.setdefault(query, [])
                hits.append((bits_val, target))
                if len(hits) > 10:
                    hits.sort(key=lambda x: x[0], reverse=True)
                    del hits[10:]

        flags: dict[str, tuple[bool, bool]] = {}
        for query, hits in top10.items():
            hits.sort(key=lambda x: x[0], reverse=True)
            top_targets = [t for _, t in hits[:10]]
            value = (
                any(t.startswith("NCLDV__") for t in top_targets),
                any(t.startswith("MIRUS__") for t in top_targets),
            )
            flags[query] = value
            flags[self._base_porf_id(query)] = value
        return flags

    def _load_hmm_targets(self) -> dict[str, set[str]]:
        """Return query -> set of HMM targets from validated hits."""
        hmm_path = self._find_hmm_hits()
        if not hmm_path or not hmm_path.exists():
            return {}
        targets: dict[str, set[str]] = {}
        with open(hmm_path) as handle:
            header = handle.readline()
            if not header:
                return targets
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                if hmm_path.name == "validated_marker_hits.tsv":
                    query = parts[0]
                    target = parts[5] if len(parts) > 5 else ""
                    validation_status = parts[7] if len(parts) > 7 else ""
                    if validation_status not in ("validated", "validated_novel"):
                        continue
                else:
                    query, target = parts[0], parts[1]
                targets.setdefault(query, set()).add(target)
                targets.setdefault(self._base_porf_id(query), set()).add(target)
        return targets

    def _genome_gc(self) -> float:
        if not self.genome_sequences:
            return 0.0
        seq = "".join(self.genome_sequences.values()).upper()
        if not seq:
            return 0.0
        gc = seq.count("G") + seq.count("C")
        return (gc / len(seq)) * 100.0

    def _region_gc(self, scaffold: str, start: int, end: int) -> float:
        seq = self.genome_sequences.get(scaffold, "")
        if not seq:
            return 0.0
        region = seq[start:end].upper()
        if not region:
            return 0.0
        gc = region.count("G") + region.count("C")
        return (gc / len(region)) * 100.0

    def generate_all(
        self,
        results: list[VerificationResult],
        accepted_only: bool = False,
        apply_v2_gate: bool = True,
        canonical_results: list[VerificationResult] | None = None,
        promoted_low_results: list[VerificationResult] | None = None,
    ) -> dict[str, Path]:
        """Generate all output files.

        Args:
            results: List of VerificationResult objects
            accepted_only: Only include accepted predictions (default: False, include all)
            apply_v2_gate: When True, apply the v2 class/length/marker
                quality gate (the canonical acceptance gate) to the
                canonical ``accepted`` output artifacts. Default True; set
                False to preserve the old
                "emit everything after the LOW filter" behavior.
            canonical_results: Explicit preselected canonical surface. When
                provided, it must contain the same result objects as
                ``results`` and no gate or legacy LOW filter is applied.
            promoted_low_results: Required with ``canonical_results``. The
                identity-preserving subset that the normal gate promoted from
                LOW confidence.

        Returns:
            Dictionary mapping output type to file path
        """
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Keep full list - all candidates are reported with confidence tiers
        all_results = list(results)

        # Build caches ONCE for all results (not subsets) to avoid incomplete caches
        all_regions = [(r.scaffold, r.start, r.end) for r in all_results]

        # Parse validated marker hits once and build region-specific mapping
        hits_by_scaffold, protein_to_models_global = self._parse_validated_marker_hits()
        self._protein_to_models_cache = self._build_region_protein_mapping(
            all_regions, hits_by_scaffold, protein_to_models_global
        )

        # Load protein counts (existing functionality)
        self._load_protein_counts_by_region(all_regions)

        promoted_low_count = 0
        if canonical_results is not None:
            if accepted_only:
                raise ValueError("accepted_only cannot accompany preselected canonical results")
            if promoted_low_results is None:
                raise ValueError("promoted_low_results is required with canonical_results")
            remaining = Counter(id(result) for result in all_results)
            for result in canonical_results:
                identity = id(result)
                if remaining[identity] <= 0:
                    raise ValueError("canonical_results must be an identity-preserving subset")
                remaining[identity] -= 1
            canonical_remaining = Counter(id(result) for result in canonical_results)
            for result in promoted_low_results:
                identity = id(result)
                if canonical_remaining[identity] <= 0:
                    raise ValueError("promoted_low_results must be an identity-preserving subset of canonical_results")
                if (result.confidence_tier or "").upper() != "LOW":
                    raise ValueError("promoted_low_results must contain only LOW results")
                canonical_remaining[identity] -= 1
            promoted_low_count = len(promoted_low_results)
            results = list(canonical_results)
        else:
            if promoted_low_results is not None:
                raise ValueError("promoted_low_results requires preselected canonical_results")
            if accepted_only and not apply_v2_gate:
                results = [r for r in results if r.is_accepted]
            if apply_v2_gate:
                # The v2 gate supersedes the legacy LOW prefilter (which was
                # stricter and would drop promotable PLV/VP or NCLDV/MIRUS LOW
                # calls before the v2 logic could promote them). The v2 gate is
                # the single, canonical acceptance gate.
                decisions = [evaluate_v2_quality_gate(r) for r in results]
                results = [r for r, decision in zip(results, decisions) if decision.kept]
                dropped = sum(1 for decision in decisions if not decision.kept)
                promoted_low_count = sum(1 for decision in decisions if decision.promoted_low)
                if dropped > 0:
                    logger.info(
                        "v2 quality gate dropped %d predictions (class/length/marker rules)",
                        dropped,
                    )
                if promoted_low_count > 0:
                    logger.info(
                        "v2 quality gate promoted %d LOW-confidence predictions",
                        promoted_low_count,
                    )
            else:
                # Legacy behavior: keep LOW only if MCP or ≥3 hallmarks.
                n_before = len(results)
                results = [r for r in results if r.confidence_tier != "LOW" or r.has_mcp or r.hallmark_count >= 3]
                n_filtered = n_before - len(results)
                if n_filtered > 0:
                    logger.info(
                        "Filtered %d LOW-confidence predictions without MCP or >=3 hallmarks",
                        n_filtered,
                    )

        if not results:
            logger.warning("No results to output; writing empty output files")
            output_files = {}
            output_files["predictions_tsv"] = self.write_predictions_tsv([])
            output_files["predictions_bed"] = self.write_predictions_bed([])
            output_files["predictions_gff"] = self.write_predictions_gff([])
            output_files["predictions_detailed_tsv"] = self.write_predictions_detailed_tsv(all_results)
            output_files["interproscan_summary_tsv"] = self.write_interproscan_summary([])
            output_files["evidence_json"] = self.write_evidence_profiles([])
            output_files["tmvec_proteins_tsv"] = self.write_tmvec_proteins_tsv(all_results)
            output_files["summary_json"] = self.write_summary(
                [],
                total_candidates=len(all_results),
                promoted_low_confidence=promoted_low_count,
            )
            return output_files

        output_files = {}

        # Generate each output type
        output_files["predictions_tsv"] = self.write_predictions_tsv(results)
        output_files["predictions_bed"] = self.write_predictions_bed(results)
        output_files["predictions_gff"] = self.write_predictions_gff(results)
        output_files["predictions_detailed_tsv"] = self.write_predictions_detailed_tsv(all_results)
        output_files["interproscan_summary_tsv"] = self.write_interproscan_summary(results)
        output_files["evidence_json"] = self.write_evidence_profiles(results)
        output_files["tmvec_proteins_tsv"] = self.write_tmvec_proteins_tsv(all_results)
        output_files["summary_json"] = self.write_summary(
            results,
            total_candidates=len(all_results),
            promoted_low_confidence=promoted_low_count,
        )
        output_files.update(self.write_gene_taxonomy(all_results))  # Write gene taxonomy for ALL EVEs

        # GVClass-compatible outputs
        if self.genome_sequences:
            gvclass_dir = self.output_dir / "gvclass_input"
            output_files["gvclass_dir"] = self.write_gvclass_export(results, gvclass_dir)

        if self.export_all_eve_sequences and self.genome_sequences:
            all_dir = self.output_dir / "eve_sequences_all"
            output_files["eve_sequences_all"] = self.write_eve_sequences(all_results, all_dir)

        logger.info(f"Generated {len(output_files)} output files in {self.output_dir}")

        return output_files

    def write_interproscan_summary(self, results: list[VerificationResult]) -> Path:
        """Write InterProScan annotation summary per region."""
        output_path = self.output_dir / "interproscan_summary.tsv"
        with atomic_write_context(output_path, "w") as f:
            if self.extended_output:
                f.write(
                    "eve_id\tinterproscan_total_hits\tinterproscan_viral_hits\t"
                    "interproscan_keywords\tinterproscan_categories\tinterproscan_families\t"
                    "interproscan_category_score\tinterproscan_score\n"
                )
            else:
                f.write(
                    "eve_id\tinterproscan_total_hits\tinterproscan_viral_hits\t"
                    "interproscan_keywords\tinterproscan_score\n"
                )
            for r in results:
                if self.extended_output:
                    f.write(
                        f"{r.eve_id}\t{r.interproscan_total_hits}\t{r.interproscan_viral_hits}\t"
                        f"{'|'.join(r.interproscan_keyword_hits) if r.interproscan_keyword_hits else '.'}\t"
                        f"{'|'.join(r.interproscan_category_hits) if r.interproscan_category_hits else '.'}\t"
                        f"{'|'.join(r.interproscan_family_hits) if r.interproscan_family_hits else '.'}\t"
                        f"{r.interproscan_category_score:.4f}\t{r.interproscan_score:.4f}\n"
                    )
                else:
                    f.write(
                        f"{r.eve_id}\t{r.interproscan_total_hits}\t{r.interproscan_viral_hits}\t"
                        f"{'|'.join(r.interproscan_keyword_hits) if r.interproscan_keyword_hits else '.'}\t"
                        f"{r.interproscan_score:.4f}\n"
                    )
        logger.info("Wrote InterProScan summary to %s", output_path)
        return output_path

    def _prediction_row(self, result: VerificationResult) -> list[str]:
        """Build one row for the canonical prediction table."""
        region_key = (result.scaffold, result.start, result.end)
        protein_to_models = dict(self._protein_to_models_cache.get(region_key, {}))
        hallmark_total, hallmark_model_counts = self._protein_model_summary(protein_to_models)
        hallmark_unique = len(hallmark_model_counts)

        row = [
            result.eve_id,
            result.scaffold,
            str(result.start),
            str(result.end),
            str(result.length),
            result.confidence_tier or "UNKNOWN",
            f"{result.final_confidence:.4f}",
            result.region_classification or ".",
            str(result.region_classification_ncldv_markers),
            str(result.region_classification_vp_plv_markers),
            str(result.region_classification_mirus_markers),
            normalize_effective_eve_class(result.likely_family),
            getattr(result, "likely_group", "") or ".",
            f"{result.kfd:.4f}",
            f"{result.gc_deviation:.4f}",
            str(hallmark_total if hallmark_total else result.hallmark_count),
            str(hallmark_unique if hallmark_unique else result.hallmark_diversity),
            str(sum(1 for gene in (result.hallmark_genes or []) if not _is_atpase_marker(gene))),
            _tsv_flag(result.has_virus_specific_marker),
            _tsv_flag(result.has_structural_support),
            "|".join(result.mcp_gene_ids) if result.mcp_gene_ids else ".",
            result.predicted_taxonomy or ".",
            f"{result.taxonomy_confidence:.4f}" if result.taxonomy_confidence else ".",
            str(result.gene_taxonomy_total),
            str(result.gene_taxonomy_ncldv_top10),
            str(result.gene_taxonomy_mirus_top10),
            str(result.gene_taxonomy_phage_top10),
            str(result.gene_taxonomy_viral_top10),
            str(result.gene_taxonomy_total_with_flanking),
            str(result.gene_taxonomy_flanking_count),
            str(result.gene_taxonomy_viral_interior),
            str(result.gene_taxonomy_viral_flanking),
            str(result.gene_taxonomy_cellular),
            str(result.gene_taxonomy_unknown),
            _tsv_flag(result.gene_taxonomy_has_ncldv_mirus),
            str(result.interproscan_total_hits),
            str(result.interproscan_viral_hits),
            "|".join(result.interproscan_keyword_hits) if result.interproscan_keyword_hits else ".",
            str(result.candidate_start) if result.candidate_start is not None else ".",
            str(result.candidate_end) if result.candidate_end is not None else ".",
            str(result.candidate_length or 0),
            str(result.candidate_reduction_bp or 0),
            result.candidate_reduction_reason or ".",
        ]
        integration = _integration_evidence_row(result)
        row.extend(integration[column] for column in INTEGRATION_EVIDENCE_COLUMNS)
        if not self.extended_output:
            row.extend([f"{result.interproscan_score:.4f}", _published_eve_class(result)])
            return row

        marker_groups = self._marker_display_groups(result)
        row.extend(
            [
                "|".join(result.interproscan_category_hits) if result.interproscan_category_hits else ".",
                "|".join(result.interproscan_family_hits) if result.interproscan_family_hits else ".",
                f"{result.interproscan_category_score:.4f}",
                f"{result.interproscan_score:.4f}",
                str(result.gene_taxonomy_vp_plv_top10),
                result.gene_taxonomy_dominant_family or "UNKNOWN",
                f"{result.gene_taxonomy_dominant_fraction:.4f}",
                (
                    result.ppv_subtype
                    if (_published_eve_class(result) == "PPV" and result.ppv_subtype in {"VP", "PLV"})
                    else "."
                ),
                str(result.host_signature_gene_count),
                f"{result.host_signature_fraction:.4f}",
                f"{getattr(result, 'host_signature_weighted_mean', 0.0):.4f}",
                "|".join(result.marker_category_hits) if result.marker_category_hits else ".",
                "|".join(result.marker_family_hits) if result.marker_family_hits else ".",
                f"{result.marker_complement_score:.4f}",
                f"{result.family_consistency_score:.4f}",
                result.vp_completeness,
                result.ppv_completeness,
                result.ncldv_completeness,
                result.mirus_completeness,
                "|".join(marker_groups.seed_legacy) if marker_groups.seed_legacy else ".",
                "|".join(marker_groups.other_legacy) if marker_groups.other_legacy else ".",
                "|".join(marker_groups.seed_patterns) if marker_groups.seed_patterns else ".",
                "|".join(marker_groups.other_patterns) if marker_groups.other_patterns else ".",
                _published_eve_class(result),
            ]
        )
        return row

    def write_predictions_tsv(self, results: list[VerificationResult]) -> Path:
        """Write predictions to TSV format.

        This is the main summary file with all metrics.
        """
        output_path = self.output_dir / "virosync_predictions.tsv"
        results = sorted(
            results,
            key=lambda r: (r.final_confidence, r.eve_id),
            reverse=True,
        )

        # Load protein counts from Phase 1 for all EVE regions
        eve_regions = [(r.scaffold, r.start, r.end) for r in results]
        self._load_protein_counts_by_region(eve_regions)
        self._ensure_protein_to_models_cache(eve_regions)

        columns = [
            "eve_id",
            "scaffold",
            "start",
            "end",
            "length",
            "confidence_tier",
            "final_confidence",
            "region_classification",
            "region_classification_ncldv_markers",
            "region_classification_vp_plv_markers",
            "region_classification_mirus_markers",
            "classification",
            "likely_group",
            "kfd",
            "gc_deviation",
            "hallmark_total",
            "hallmark_unique",
            "hallmark_non_atpase",
            "has_virus_specific",
            "has_structural_support",
            "mcp_gene_ids",
            "predicted_taxonomy",
            "taxonomy_confidence",
            "gene_taxonomy_total",
            "gene_taxonomy_ncldv_top10",
            "gene_taxonomy_mirus_top10",
            "gene_taxonomy_phage_top10",
            "gene_taxonomy_viral_top10",
            "gene_taxonomy_total_with_flanking",
            "gene_taxonomy_flanking_count",
            "gene_taxonomy_viral_interior",
            "gene_taxonomy_viral_flanking",
            "gene_taxonomy_cellular",
            "gene_taxonomy_unknown",
            "gene_taxonomy_has_ncldv_mirus",
            "interproscan_total_hits",
            "interproscan_viral_hits",
            "interproscan_keyword_hits",
            "candidate_start",
            "candidate_end",
            "candidate_length",
            "candidate_reduction_bp",
            "candidate_reduction_reason",
        ]
        columns.extend(INTEGRATION_EVIDENCE_COLUMNS)
        if self.extended_output:
            columns.extend(
                [
                    "interproscan_category_hits",
                    "interproscan_family_hits",
                    "interproscan_category_score",
                    "interproscan_score",
                    "gene_taxonomy_vp_plv_top10",
                    "gene_taxonomy_dominant_family",
                    "gene_taxonomy_dominant_fraction",
                    "ppv_subtype",
                    "host_signature_gene_count",
                    "host_signature_fraction",
                    "host_signature_weighted_mean",
                    "marker_category_hits",
                    "marker_family_hits",
                    "marker_complement_score",
                    "family_consistency_score",
                    "vp_completeness",
                    "ppv_completeness",
                    "ncldv_completeness",
                    "mirus_completeness",
                    "seed_marker_names",
                    "other_marker_names",
                    "seed_marker_patterns",  # NEW: protein-pattern grouping
                    "other_marker_patterns",  # NEW: protein-pattern grouping
                ]
            )
        else:
            columns.append("interproscan_score")
        columns.append("effective_eve_class")

        with atomic_write_context(output_path, "w") as handle:
            handle.write("\t".join(columns) + "\n")

            for result in results:
                handle.write("\t".join(self._prediction_row(result)) + "\n")

        logger.info(f"Wrote {len(results)} predictions to {output_path}")
        return output_path

    def _detailed_prediction_row(
        self,
        result: VerificationResult,
        sources: _DetailedOutputSources,
        columns: list[str],
    ) -> list[str]:
        """Build one row for the detailed prediction table."""
        porfs = [
            porf_id
            for porf_start, porf_end, porf_id in sources.porfs_by_scaffold.get(result.scaffold, [])
            if porf_id not in sources.rescued_protein_ids and porf_start < result.end and porf_end > result.start
        ]

        region_key = (result.scaffold, result.start, result.end)
        protein_to_models = dict(self._protein_to_models_cache.get(region_key, {}))
        if not protein_to_models and sources.hmm_targets and porfs:
            for porf_id in porfs:
                base_porf_id = self._base_porf_id(porf_id)
                models = {
                    self._normalize_model_name(model)
                    for model in (
                        sources.hmm_targets.get(base_porf_id, set()) | sources.hmm_targets.get(porf_id, set())
                    )
                }
                if models:
                    protein_to_models[base_porf_id] = models

        hallmark_total, hallmark_model_counts = self._protein_model_summary(protein_to_models)
        hallmark_unique = len(hallmark_model_counts)
        gvogm_count, gvogm_model_counts = self._protein_model_summary(protein_to_models, prefix="GVOGM")
        og_count, og_model_counts = self._protein_model_summary(protein_to_models, prefix="OG")

        unvalidated_gvogm_names: set[str] = set()
        unvalidated_og_names: set[str] = set()
        for hit_start, hit_end, target, status, _porf_id, _score in sources.marker_hits_by_scaffold.get(
            result.scaffold, []
        ):
            if hit_start >= result.end or hit_end <= result.start:
                continue
            if status in ("validated", "validated_novel"):
                continue
            target_upper = target.upper()
            if target_upper.startswith("GVOGM"):
                unvalidated_gvogm_names.add(target)
            elif target_upper.startswith("OG"):
                unvalidated_og_names.add(target)

        taxonomy_counts = {prefix: 0 for prefix in DETAILED_TAXONOMY_PARTITION}
        top10_support = {prefix: 0 for prefix in ("NCLDV", "MIRUS", "PPV", "CRESS")}

        all_gene_tax_records = getattr(result, "gene_taxonomy_records", []) or []
        gene_tax_records: list[dict] = []
        for record in all_gene_tax_records:
            is_flanking = (
                record.get("is_flanking", False) if isinstance(record, dict) else getattr(record, "is_flanking", False)
            )
            if not is_flanking:
                gene_tax_records.append(record if isinstance(record, dict) else getattr(record, "__dict__", {}))
        summary_gene_tax_total = getattr(result, "gene_taxonomy_total", 0) or 0
        for record in gene_tax_records:
            raw_viral = set(_taxonomy_prefixes(record.get("top10_prefixes")))
            for family in top10_support:
                if family in raw_viral:
                    top10_support[family] += 1
            taxonomy_counts[_taxonomy_partition_bucket(record)] += 1

        if not gene_tax_records:
            for porf_id in porfs:
                base_porf_id = self._base_porf_id(porf_id)
                flags = sources.diamond_flags.get(base_porf_id) or sources.diamond_flags.get(porf_id)
                if not flags:
                    continue
                has_ncldv, has_mirus = flags
                if has_ncldv:
                    top10_support["NCLDV"] += 1
                if has_mirus:
                    top10_support["MIRUS"] += 1
            taxonomy_counts["UNK"] = len(porfs)

        total_proteins = max(summary_gene_tax_total, len(gene_tax_records)) if gene_tax_records else len(porfs)
        observed_total = sum(taxonomy_counts.values())
        if observed_total < total_proteins:
            taxonomy_counts["UNK"] += total_proteins - observed_total

        region_gc = self._region_gc(result.scaffold, result.start, result.end)
        taxonomy_summary = ";".join(f"{prefix}:{taxonomy_counts[prefix]}" for prefix in DETAILED_TAXONOMY_PARTITION)
        marker_groups = (
            self._marker_display_groups(result) if self.extended_output else _MarkerDisplayGroups((), (), (), ())
        )
        published_class = _published_eve_class(result)
        ppv_subtype = result.ppv_subtype if published_class == "PPV" and result.ppv_subtype in {"VP", "PLV"} else "."
        row_values = {
            "eve_id": result.eve_id,
            "scaffold": result.scaffold,
            "start": str(result.start),
            "end": str(result.end),
            "length": str(result.length),
            "confidence_tier": result.confidence_tier or "UNKNOWN",
            "final_confidence": f"{result.final_confidence:.4f}",
            "effective_eve_class": published_class,
            "likely_family": normalize_effective_eve_class(result.likely_family),
            "ppv_subtype": ppv_subtype,
            "likely_group": getattr(result, "likely_group", "") or ".",
            "candidate_start": str(result.candidate_start) if result.candidate_start is not None else ".",
            "candidate_end": str(result.candidate_end) if result.candidate_end is not None else ".",
            "candidate_length": str(result.candidate_length or 0),
            "candidate_reduction_bp": str(result.candidate_reduction_bp or 0),
            "candidate_reduction_reason": result.candidate_reduction_reason or ".",
            "seed_sources": "|".join(sorted(result.seed_sources)) if result.seed_sources else ".",
            "canonical_selection_outcome": result.canonical_selection_outcome or ".",
            **_integration_evidence_row(result),
            "hallmark_total": str(hallmark_total if hallmark_total else result.hallmark_count),
            "hallmark_unique": str(hallmark_unique if hallmark_unique else result.hallmark_diversity),
            "mcp_gene_ids": "|".join(result.mcp_gene_ids) if result.mcp_gene_ids else ".",
            "tier1_bypassed_marker_count": str(len(result.tier1_bypassed_marker_ids)),
            "tier1_bypassed_marker_ids": (
                "|".join(result.tier1_bypassed_marker_ids) if result.tier1_bypassed_marker_ids else "."
            ),
            "tier1_bypassed_marker_models": (
                "|".join(result.tier1_bypassed_marker_models) if result.tier1_bypassed_marker_models else "."
            ),
            "gvogm_count": str(gvogm_count),
            "gvogm_names": self._format_counted_names(gvogm_model_counts),
            "og_count": str(og_count),
            "og_names": self._format_counted_names(og_model_counts),
            "gvogm_unvalidated_count": str(len(unvalidated_gvogm_names)),
            "gvogm_unvalidated_names": (",".join(sorted(unvalidated_gvogm_names)) if unvalidated_gvogm_names else "."),
            "og_unvalidated_count": str(len(unvalidated_og_names)),
            "og_unvalidated_names": ",".join(sorted(unvalidated_og_names)) if unvalidated_og_names else ".",
            "marker_complement_score": f"{result.marker_complement_score:.4f}",
            "family_consistency_score": f"{result.family_consistency_score:.4f}",
            "seed_marker_names": "|".join(marker_groups.seed_legacy) if marker_groups.seed_legacy else ".",
            "other_marker_names": "|".join(marker_groups.other_legacy) if marker_groups.other_legacy else ".",
            "seed_marker_patterns": "|".join(marker_groups.seed_patterns) if marker_groups.seed_patterns else ".",
            "other_marker_patterns": "|".join(marker_groups.other_patterns) if marker_groups.other_patterns else ".",
            "total_proteins": str(total_proteins),
            "ncldv_top10_proteins": str(top10_support["NCLDV"]),
            "mirus_top10_proteins": str(top10_support["MIRUS"]),
            "ppv_top10_proteins": str(top10_support["PPV"]),
            "cress_top10_proteins": str(top10_support["CRESS"]),
            "taxonomy_best_hits": taxonomy_summary,
            "kfd": f"{result.kfd:.4f}",
            "gc_deviation": f"{result.gc_deviation:.4f}",
            "region_gc_percent": f"{region_gc:.3f}",
            "genome_gc_percent": f"{sources.genome_gc:.3f}",
            "gc_delta": f"{(region_gc - sources.genome_gc):.3f}",
            "host_signature_gene_count": str(result.host_signature_gene_count),
            "host_signature_fraction": f"{result.host_signature_fraction:.4f}",
            "host_signature_weighted_mean": f"{getattr(result, 'host_signature_weighted_mean', 0.0):.4f}",
            "interproscan_total_hits": str(result.interproscan_total_hits),
            "interproscan_viral_hits": str(result.interproscan_viral_hits),
            "interproscan_keyword_hits": (
                "|".join(result.interproscan_keyword_hits) if result.interproscan_keyword_hits else "."
            ),
            "interproscan_category_score": f"{result.interproscan_category_score:.4f}",
            "interproscan_score": f"{result.interproscan_score:.4f}",
            "vp_completeness": result.vp_completeness,
            "ppv_completeness": result.ppv_completeness,
            "ncldv_completeness": result.ncldv_completeness,
            "mirus_completeness": result.mirus_completeness,
            "ani_cluster_id": str(result.cluster_id) if result.cluster_id >= 0 else ".",
            "ani_cluster_size": str(result.cluster_size),
            "ani_max_percent": f"{result.max_cluster_ani:.4f}" if result.max_cluster_ani > 0 else ".",
            "taxonomy_class_before_ani": (
                normalize_effective_eve_class(result.taxonomy_class_before_ani)
                if result.taxonomy_class_before_ani
                else "."
            ),
            "taxonomy_class_propagated_from": result.taxonomy_class_propagated_from or ".",
        }
        return [row_values[column] for column in columns]

    def write_predictions_detailed_tsv(self, results: list[VerificationResult]) -> Path:
        """Write detailed predictions with GVOG/OG counts, protein counts, and GC stats."""
        output_path = self.output_dir / "virosync_predictions_detailed.tsv"
        results = sorted(
            results,
            key=lambda result: (result.final_confidence, result.eve_id),
            reverse=True,
        )
        sources = _DetailedOutputSources(
            porfs_by_scaffold=self._load_porfs_by_scaffold(),
            rescued_protein_ids=set(self._rescued_protein_sequences or {}),
            diamond_flags=self._load_diamond_top10_flags(),
            hmm_targets=self._load_hmm_targets(),
            marker_hits_by_scaffold=self._load_marker_hits(),
            genome_gc=self._genome_gc(),
        )

        eve_regions = [(result.scaffold, result.start, result.end) for result in results]
        self._load_protein_counts_by_region(eve_regions)
        self._ensure_protein_to_models_cache(eve_regions)

        columns = [
            column
            for column in DETAILED_PREDICTION_COLUMNS
            if self.extended_output or column not in DETAILED_PREDICTION_EXTENDED_COLUMNS
        ]

        with atomic_write_context(output_path, "w") as handle:
            handle.write("\t".join(columns) + "\n")
            for result in results:
                row = self._detailed_prediction_row(result, sources, columns)
                handle.write("\t".join(row) + "\n")

        logger.info(f"Wrote {len(results)} detailed predictions to {output_path}")
        return output_path

    def write_predictions_bed(self, results: list[VerificationResult]) -> Path:
        """Write predictions to BED6 format."""
        output_path = self.output_dir / "virosync_predictions.bed"

        with atomic_write_context(output_path, "w") as f:
            for r in results:
                persisted_confidence = float(f"{r.final_confidence:.4f}")
                score = int(min(1000, persisted_confidence * 1000))
                strand = "."
                f.write(f"{r.scaffold}\t{r.start}\t{r.end}\t{r.eve_id}\t{score}\t{strand}\n")

        logger.info(f"Wrote {len(results)} predictions to {output_path}")
        return output_path

    def write_predictions_gff(self, results: list[VerificationResult]) -> Path:
        """Write predictions to GFF3 format with attributes."""
        output_path = self.output_dir / "virosync_predictions.gff3"

        with atomic_write_context(output_path, "w") as f:
            f.write("##gff-version 3\n")
            # Version, not a wall-clock timestamp: a generation time in a data file
            # makes two runs of the same input differ, and virosync_summary.json
            # already records generated_at for provenance.
            f.write(f"# ViroSync predictions, virosync {__version__}\n")

            for r in results:
                persisted_confidence = float(f"{r.final_confidence:.4f}")
                score = int(min(1000, persisted_confidence * 1000))
                published_class = _published_eve_class(r)

                # Build attributes
                attrs = [
                    f"ID={_gff3_escape(r.eve_id)}",
                    f"Name={_gff3_escape(r.eve_id)}",
                    f"confidence={r.final_confidence:.4f}",
                    f"status={_gff3_escape(r.status.value)}",
                    f"hallmark_diversity={r.hallmark_diversity}",
                    f"effective_eve_class={_gff3_escape(published_class)}",
                ]
                if published_class == "PPV" and r.ppv_subtype in {"VP", "PLV"}:
                    attrs.append(f"ppv_subtype={_gff3_escape(r.ppv_subtype)}")

                if r.region_classification:
                    attrs.append(f"region_classification={_gff3_escape(r.region_classification)}")
                if r.has_virus_specific_marker:
                    attrs.append("has_virus_specific=true")
                if r.has_structural_support:
                    attrs.append("has_structural_support=true")
                if r.predicted_taxonomy:
                    attrs.append(f"taxonomy={_gff3_escape(r.predicted_taxonomy)}")
                if r.gene_taxonomy_total:
                    attrs.append(f"gene_taxonomy_total={r.gene_taxonomy_total}")
                    attrs.append(f"gene_taxonomy_viral_top10={r.gene_taxonomy_viral_top10}")
                    attrs.append(f"gene_taxonomy_total_with_flanking={r.gene_taxonomy_total_with_flanking}")
                    attrs.append(f"gene_taxonomy_flanking_count={r.gene_taxonomy_flanking_count}")
                    attrs.append(f"gene_taxonomy_viral_interior={r.gene_taxonomy_viral_interior}")
                    attrs.append(f"gene_taxonomy_viral_flanking={r.gene_taxonomy_viral_flanking}")

                attr_str = ";".join(attrs)

                # GFF columns: seqid source type start end score strand phase attributes
                f.write(
                    f"{_gff3_escape(r.scaffold)}\tViroSync\tEVE\t{r.start + 1}\t{r.end}\t{score}\t.\t.\t{attr_str}\n"
                )

        logger.info(f"Wrote {len(results)} predictions to {output_path}")
        return output_path

    def write_gene_taxonomy(self, results: list[VerificationResult]) -> dict[str, Path]:
        """Write per-candidate and aggregate gene taxonomy tables."""
        filename_components = self._eve_filename_components(results)
        output_dir = self.output_dir / "gene_taxonomy"
        output_dir.mkdir(parents=True, exist_ok=True)

        taxonomy_tables = [
            (
                result,
                [_normalize_gene_taxonomy_record(record) for record in result.gene_taxonomy_records],
            )
            for result in results
            if result.gene_taxonomy_records
        ]

        output_files: dict[str, Path] = {}
        for result, records in taxonomy_tables:
            output_path = output_dir / f"{filename_components[result.eve_id]}_gene_taxonomy.tsv"
            require_strict_child(output_dir, output_path)
            output_files[f"gene_taxonomy_{result.eve_id}"] = output_path
            with atomic_write_context(output_path, "w", encoding="utf-8") as handle:
                handle.write("\t".join(_GENE_TAXONOMY_COLUMNS) + "\n")
                for record in records:
                    handle.write("\t".join(record.per_eve_row()) + "\n")

        if output_files:
            logger.info("Wrote %d gene taxonomy tables to %s", len(output_files), output_dir)

        all_records_path = output_dir / "gene_taxonomy_all.tsv"
        all_records_count = 0
        with atomic_write_context(all_records_path, "w", encoding="utf-8") as handle:
            handle.write("\t".join(_COMBINED_GENE_TAXONOMY_COLUMNS) + "\n")
            for result, records in taxonomy_tables:
                for record in records:
                    handle.write("\t".join(record.combined_row(result.eve_id)) + "\n")
                    all_records_count += 1

        if all_records_count > 0:
            output_files["gene_taxonomy_all"] = all_records_path
            logger.info("Wrote combined gene taxonomy: %d records to %s", all_records_count, all_records_path)

        return output_files

    def write_evidence_profiles(self, results: list[VerificationResult]) -> Path:
        """Write detailed evidence profiles to JSON."""
        output_path = self.output_dir / "evidence_profiles.json"

        profiles = {}
        for r in results:
            profile = r.to_dict()

            # Add coherence details if available
            if r.coherence_analysis:
                profile["coherence_details"] = {
                    "interpretation": r.coherence_analysis.interpretation,
                    "confidence_level": r.coherence_analysis.confidence_level,
                }
                if r.coherence_analysis.profile:
                    profile["evidence_coverage"] = {
                        k.value: v for k, v in r.coherence_analysis.profile.evidence_coverage.items()
                    }

            # Add structural details
            if r.structural_results:
                profile["structural_hits"] = [
                    {
                        "porf_id": sr.porf_id,
                        "supports_viral": sr.supports_viral_origin,
                        "score": sr.structural_evidence_score,
                        "prediction_plddt": sr.prediction.mean_plddt if sr.prediction else None,
                    }
                    for sr in r.structural_results[:10]  # Limit
                ]

            profiles[r.eve_id] = profile

        with atomic_write_context(output_path, "w") as f:
            json.dump(profiles, f, indent=2, default=str)

        logger.info(f"Wrote evidence profiles to {output_path}")
        return output_path

    def write_tmvec_proteins_tsv(self, results: list[VerificationResult]) -> Path:
        """Write per-protein TMVec2 BFVD hits."""
        output_path = self.output_dir / "virosync_tmvec_proteins.tsv"
        header = [
            "eve_id",
            "porf_id",
            "length",
            "tmvec_bfvd_score",
            "tmvec_bfvd_hit",
            "tmvec_bfvd_annotation",
            "tmvec_bfvd_organism",
            "tmvec_bfvd_lineage",
            "tmvec_bfvd_keywords",
        ]
        with atomic_write_context(output_path, "w") as f:
            f.write("\t".join(header) + "\n")
            for result in results:
                for record in result.tmvec_all_proteins:
                    row = [str(record.get(col, "")) for col in header]
                    f.write("\t".join(row) + "\n")
        logger.info("Wrote TMVec per-protein hits to %s", output_path)
        return output_path

    def write_summary(
        self,
        results: list[VerificationResult],
        *,
        total_candidates: int | None = None,
        promoted_low_confidence: int = 0,
    ) -> Path:
        """Write summary statistics for the selected canonical predictions."""
        output_path = self.output_dir / "virosync_summary.json"

        total = len(results)
        candidate_total = total if total_candidates is None else total_candidates
        high_conf = sum(1 for r in results if (getattr(r, "confidence_tier", "") or "").upper() == "HIGH")
        medium_conf = sum(1 for r in results if (getattr(r, "confidence_tier", "") or "").upper() == "MEDIUM")
        low_conf = sum(1 for r in results if (getattr(r, "confidence_tier", "") or "").upper() == "LOW")
        if (
            type(promoted_low_confidence) is not int
            or promoted_low_confidence < 0
            or promoted_low_confidence > low_conf
        ):
            raise ValueError("promoted_low_confidence must be a nonnegative count of canonical LOW results")

        total_length = sum(r.length for r in results)
        with_virus_specific = sum(1 for r in results if r.has_virus_specific_marker)
        with_structural = sum(1 for r in results if r.has_structural_support)

        confidences = [r.final_confidence for r in results]

        summary = {
            "generated_at": datetime.now().isoformat(),
            "virosync_version": __version__,
            **coordinate_contract_metadata(),
            "statistics": {
                "total_candidates": candidate_total,
                "canonical_predictions": total,
                "high_medium_confidence": high_conf + medium_conf,
                "high_confidence": high_conf,
                "medium_confidence": medium_conf,
                "low_confidence": low_conf,
                "promoted_low_confidence": promoted_low_confidence,
                "total_accepted_length_bp": total_length,
                "with_virus_specific_markers": with_virus_specific,
                "with_structural_support": with_structural,
                "mean_confidence": float(np.mean(confidences)) if confidences else 0.0,
                "median_confidence": float(np.median(confidences)) if confidences else 0.0,
            },
            "per_scaffold": self._summarize_per_scaffold(results),
        }

        with atomic_write_context(output_path, "w") as f:
            json.dump(summary, f, indent=2)

        logger.info(f"Wrote summary to {output_path}")
        return output_path

    def _summarize_per_scaffold(self, results: list[VerificationResult]) -> dict:
        """Summarize canonical v2-gated results per scaffold."""
        per_scaffold = {}

        for r in results:
            if r.scaffold not in per_scaffold:
                per_scaffold[r.scaffold] = {
                    "count": 0,
                    "canonical_predictions": 0,
                    "high_medium_confidence": 0,
                    "total_length": 0,
                }

            per_scaffold[r.scaffold]["count"] += 1
            per_scaffold[r.scaffold]["canonical_predictions"] += 1
            per_scaffold[r.scaffold]["total_length"] += r.length
            if (getattr(r, "confidence_tier", "") or "").upper() in {"HIGH", "MEDIUM"}:
                per_scaffold[r.scaffold]["high_medium_confidence"] += 1

        return per_scaffold

    def write_gvclass_export(
        self,
        results: list[VerificationResult],
        output_dir: Path,
    ) -> Path:
        """Write GVClass-compatible export.

        Creates per-element FASTA files that can be directly input to GVClass.
        """
        # ``results`` is already the canonical v2-gated accepted set.
        accepted = list(results)
        filename_components = self._eve_filename_components(accepted)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Write nucleotide sequences
        nuc_dir = output_dir / "nucleotide"
        nuc_dir.mkdir(exist_ok=True)

        for r in accepted:
            if r.scaffold in self.genome_sequences:
                seq = self.genome_sequences[r.scaffold][r.start : r.end]
                record = SeqRecord(
                    Seq(seq),
                    id=r.eve_id,
                    description=f"scaffold={r.scaffold} start={r.start} end={r.end} confidence={r.final_confidence:.4f}",
                )
                output_path = nuc_dir / f"{filename_components[r.eve_id]}.fna"
                require_strict_child(nuc_dir, output_path)
                SeqIO.write([record], output_path, "fasta")

        # Write ordinary proteins and any confirmed frameshift-rescued domains.
        if self.proteome_sequences or self._find_confirmed_frameshift_proteins():
            prot_dir = output_dir / "protein"
            prot_dir.mkdir(exist_ok=True)

            for r in accepted:
                eve_proteins = self._protein_records_for_region(
                    r.scaffold,
                    r.start,
                    r.end,
                )
                if eve_proteins:
                    output_path = prot_dir / f"{filename_components[r.eve_id]}.faa"
                    require_strict_child(prot_dir, output_path)
                    SeqIO.write(eve_proteins, output_path, "fasta")

        # Write manifest
        manifest_path = output_dir / "manifest.tsv"
        with open(manifest_path, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")
            writer.writerow(["eve_id", "nucleotide_fasta", "protein_fasta", "confidence"])
            for r in accepted:
                component = filename_components[r.eve_id]
                nuc_path = (Path("nucleotide") / f"{component}.fna").as_posix()
                prot_path = (Path("protein") / f"{component}.faa").as_posix()
                if not (output_dir / nuc_path).exists():
                    nuc_path = ""
                if not (output_dir / prot_path).exists():
                    prot_path = ""
                writer.writerow([r.eve_id, nuc_path, prot_path, f"{r.final_confidence:.4f}"])

        logger.info(f"Wrote GVClass export for {len(accepted)} EVEs to {output_dir}")
        return output_dir

    def write_eve_sequences(
        self,
        results: list[VerificationResult],
        output_dir: Path,
    ) -> Path:
        """Write per-EVE nucleotide/protein FASTA files for all results.

        This export is intended for downstream manual inspection and
        includes high/medium confidence and low confidence regions.
        """
        filename_components = self._eve_filename_components(results)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Write nucleotide sequences
        nuc_dir = output_dir / "nucleotide"
        nuc_dir.mkdir(exist_ok=True)

        for r in results:
            if r.scaffold in self.genome_sequences:
                seq = self.genome_sequences[r.scaffold][r.start : r.end]
                record = SeqRecord(
                    Seq(seq),
                    id=r.eve_id,
                    description=(
                        f"scaffold={r.scaffold} start={r.start} end={r.end} "
                        f"status={r.status.value} confidence={r.final_confidence:.4f}"
                    ),
                )
                output_path = nuc_dir / f"{filename_components[r.eve_id]}.fna"
                require_strict_child(nuc_dir, output_path)
                SeqIO.write([record], output_path, "fasta")

        # Write ordinary proteins and any confirmed frameshift-rescued domains.
        if self.proteome_sequences or self._find_confirmed_frameshift_proteins():
            prot_dir = output_dir / "protein"
            prot_dir.mkdir(exist_ok=True)

            for r in results:
                eve_proteins = self._protein_records_for_region(
                    r.scaffold,
                    r.start,
                    r.end,
                )
                if eve_proteins:
                    output_path = prot_dir / f"{filename_components[r.eve_id]}.faa"
                    require_strict_child(prot_dir, output_path)
                    SeqIO.write(eve_proteins, output_path, "fasta")

        # Write manifest
        manifest_path = output_dir / "manifest.tsv"
        with open(manifest_path, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")
            writer.writerow(["eve_id", "nucleotide_fasta", "protein_fasta", "status", "confidence"])
            for r in results:
                component = filename_components[r.eve_id]
                nuc_path = (Path("nucleotide") / f"{component}.fna").as_posix()
                prot_path = (Path("protein") / f"{component}.faa").as_posix()
                if not (output_dir / nuc_path).exists():
                    nuc_path = ""
                if not (output_dir / prot_path).exists():
                    prot_path = ""
                writer.writerow(
                    [
                        r.eve_id,
                        nuc_path,
                        prot_path,
                        r.status.value,
                        f"{r.final_confidence:.4f}",
                    ]
                )

        logger.info(f"Wrote all-EVE sequence export for {len(results)} EVEs to {output_dir}")
        return output_dir

    @genome_sequences.setter
    def genome_sequences(self, value: dict[str, str]):
        """Allow setting genome sequences directly."""
        self._genome_sequences = value

    def write_combined_eve_fasta(
        self,
        results: list[VerificationResult],
        output_path: Path,
    ) -> Path:
        """Write all EVE sequences to a single multi-FASTA file.

        Args:
            results: List of VerificationResult objects
            output_path: Output FASTA path

        Returns:
            Path to written FASTA file
        """
        records = []
        for r in results:
            if r.scaffold in self.genome_sequences:
                seq = self.genome_sequences[r.scaffold][r.start : r.end]
                record = SeqRecord(
                    Seq(seq),
                    id=r.eve_id,
                    description=(
                        f"scaffold={r.scaffold} start={r.start} end={r.end} "
                        f"tier={r.confidence_tier} confidence={r.final_confidence:.4f}"
                    ),
                )
                records.append(record)

        if records:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            SeqIO.write(records, output_path, "fasta")
            logger.info(f"Wrote {len(records)} EVE sequences to {output_path}")
        else:
            logger.warning(f"No EVE sequences to write to {output_path}")

        return output_path
