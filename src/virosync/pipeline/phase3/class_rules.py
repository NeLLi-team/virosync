"""Name accepted calls by their capsid profiles.

ViroSync takes the class of a region from the top-ten hits of its genes. A capsid
gene whose viral hits are only EUK-pEVE references, or whose hits span two capsid
classes, leaves the region unclassified even when the capsid profile that found it
belongs to one lineage. Two rules use the capsid profiles instead:

- A call whose capsid profiles name only Mirusviricota is MIRUS, whatever its other
  markers or gene-taxonomy class say. Mirusviricota combine an HK97-fold capsid with
  NCLDV-like informational genes, so NCLDV markers beside a Mirus capsid are expected.
- An unclassified call takes the class its capsid profiles name, when they name
  exactly one class and its marker families name no other.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

from virosync.output_contract import normalize_effective_eve_class
from virosync.pipeline.phase1.marker_roles import is_ordinary_marker_evidence

from .evidence_synthesizer import VerificationResult
from .mcp_detection import capsid_profile_class

__all__ = ["apply_capsid_class_rules", "capsid_class"]

UNCLASSIFIED_CLASSES = ("VIRAL_UNKNOWN", "UNKNOWN")


def apply_capsid_class_rules(results: Iterable[VerificationResult], markers: Sequence[object]) -> int:
    """Set ``taxonomy_class`` of the results the capsid rules name; return how many changed.

    ``markers`` are the Phase 1 marker hits of the genome. A call reads the best
    model of each gene that overlaps it.
    """
    markers_by_scaffold: dict[str, list[object]] = {}
    for marker in markers:
        markers_by_scaffold.setdefault(marker.scaffold, []).append(marker)
    changed = 0
    for result in results:
        names = _region_marker_names(result, markers_by_scaffold.get(result.scaffold, []))
        named = _ruled_class(
            normalize_effective_eve_class(result.taxonomy_class),
            capsid_class(names, result.marker_family_hits),
        )
        if named:
            result.taxonomy_class = named
            changed += 1
    return changed


def capsid_class(marker_names: Iterable[str], marker_families: Iterable[str]) -> str:
    """Return the class the capsid profiles name, or ``""``.

    The capsids must name one class. The marker families must name no other,
    except beside a Mirus capsid.
    """
    classes = {capsid_profile_class(name) for name in marker_names} - {""}
    if len(classes) != 1:
        return ""
    (named,) = classes
    if named != "MIRUS" and set(marker_families) - {named}:
        return ""
    return named


def _ruled_class(current: str, by_capsid: str) -> str:
    if by_capsid == "MIRUS" and current != "MIRUS":
        return "MIRUS"
    if current in UNCLASSIFIED_CLASSES and by_capsid:
        return by_capsid
    return ""


def _region_marker_names(result: VerificationResult, markers: Sequence[object]) -> list[str]:
    """Return the best model of each gene overlapping the call: ordinary evidence first, then score."""
    best_by_gene: dict[str, tuple[tuple[bool, float], str]] = {}
    for marker in markers:
        if not (marker.start < result.end and marker.end > result.start):
            continue
        gene = marker.query_porf.split("|aa")[0]
        rank = (is_ordinary_marker_evidence(marker), float(marker.hmm_score))
        if gene not in best_by_gene or rank > best_by_gene[gene][0]:
            best_by_gene[gene] = (rank, marker.hmm_target)
    return [model for _, model in best_by_gene.values()]
