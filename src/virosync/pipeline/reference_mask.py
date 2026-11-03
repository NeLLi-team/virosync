"""Reference mask: the genomes whose database proteins a run must not use as evidence.

The reference databases hold the proteomes of genomes that ViroSync later screens,
and of benchmark hosts. Their own endogenous viral elements carry host or
``EUK-pEVE`` labels, so a hit to the genome under test either hides or inflates
an element. A run therefore ignores every hit to its own assembly and, when
asked, to the genomes of named taxa. The databases stay complete: the search
wrapper removes the masked rows from DIAMOND's results.

A mask is a set of genome keys, the part of a subject ID between ``LABEL__`` and
``|``. A key masks the genome under every label (``EUK__``, ``EUK-pEVE__``,
``MITO__``, ``PLASTID__``, viral labels).
"""

import logging
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

logger = logging.getLogger(__name__)

ASSEMBLY_ACCESSION = re.compile(r"(GC[AF])_(\d{9})[._](\d+)")
EUKPROT_MARKER_KEY = re.compile(r"(EP[0-9]{5})_.+")


@dataclass(frozen=True, slots=True)
class ReferenceMask:
    """Genome keys to ignore in reference searches, and how they were chosen."""

    genome_keys: frozenset[str] = frozenset()
    taxa: tuple[str, ...] = ()
    own_assembly_keys: tuple[str, ...] = ()

    def to_record(self) -> dict[str, object]:
        """Return the JSON-safe summary that provenance.json stores."""
        return {
            "n_genome_keys": len(self.genome_keys),
            "taxa": list(self.taxa),
            "own_assembly_keys": list(self.own_assembly_keys),
        }


def genome_key(subject_id: str) -> str:
    """Return the genome key of a subject ID or organism key.

    ``EUK__GCA-052324655-1|XP_1`` and ``EUK__GCA-052324655-1`` give
    ``GCA-052324655-1``. An ID without a ``LABEL__`` head is its own key up to ``|``.
    """
    head = subject_id.partition("|")[0]
    _label, separator, key = head.partition("__")
    return key if separator else head


def is_masked(subject_id: str, genome_keys: frozenset[str]) -> bool:
    """Whether a subject ID belongs to a masked genome.

    Legacy EukProt markers can carry an organism name after their five-digit key
    (``EUK__EP00151_Rhizophagus_irregularis``) where the labels use the bare key
    (``EUK__EP00151``). Other canonical keys match exactly; an underscore suffix
    can identify a different assembly.
    """
    key = genome_key(subject_id)
    if key in genome_keys:
        return True
    legacy = EUKPROT_MARKER_KEY.fullmatch(key)
    return legacy is not None and legacy.group(1) in genome_keys


def build_reference_mask(
    genome_id: str,
    *,
    mask_own_assembly: bool,
    taxa: Sequence[str],
    labels: Mapping[str, str] | None,
    self_keys: Sequence[str] = (),
) -> ReferenceMask:
    """Build the mask of one genome run.

    Args:
        genome_id: Run identifier; a GenBank or RefSeq assembly accession inside
            it names the own assembly.
        mask_own_assembly: Mask the own assembly when ``genome_id`` holds an accession.
        taxa: Taxon names; each masks every genome whose lineage holds that name
            at any rank, compared ignoring case and surrounding blanks.
        labels: Taxonomy labels, ``LABEL__genomekey`` to ``ROOT|...|species``.
            Only needed when ``taxa`` is not empty.
        self_keys: Reference keys of the genome's own proteomes, bare or with a
            ``LABEL__`` head. They replace the accession read from ``genome_id``.

    Raises:
        ValueError: If ``taxa`` is not empty and ``labels`` is ``None``.
    """
    if not mask_own_assembly:
        own_assembly_keys: tuple[str, ...] = ()
    elif self_keys:
        own_assembly_keys = tuple(dict.fromkeys(twin for key in self_keys for twin in _explicit_assembly_keys(key)))
    else:
        own_assembly_keys = _own_assembly_keys(genome_id)
    return ReferenceMask(
        genome_keys=frozenset(own_assembly_keys) | _taxon_genome_keys(taxa, labels),
        taxa=tuple(taxa),
        own_assembly_keys=own_assembly_keys,
    )


def _explicit_assembly_keys(subject_id: str) -> tuple[str, ...]:
    """Keep exact keys, adding only the proven NCBI accession spelling twin."""
    key = genome_key(subject_id)
    match = re.fullmatch(r"(GC[AF])-([0-9]{9})-([0-9]+)", key) or re.fullmatch(r"(GC[AF])_([0-9]{9})\.([0-9]+)", key)
    if match is None:
        return (key,)
    prefix, number, version = match.groups()
    return f"{prefix}-{number}-{version}", f"{prefix}_{number}.{version}"


def _own_assembly_keys(genome_id: str) -> tuple[str, ...]:
    """Return the keys of the assembly accession in ``genome_id``, if it has one.

    The label form (``GCA-052324655-1``) keys nuclear proteomes. The accession
    itself (``GCF_000001215.4``) keys the organelle rows of the same assembly.
    """
    match = ASSEMBLY_ACCESSION.search(genome_id)
    if match is None:
        return ()
    prefix, number, version = match.groups()
    return f"{prefix}-{number}-{version}", f"{prefix}_{number}.{version}"


def _taxon_genome_keys(taxa: Sequence[str], labels: Mapping[str, str] | None) -> frozenset[str]:
    if not taxa:
        return frozenset()
    if labels is None:
        raise ValueError("reference_mask_taxa needs the taxonomy labels, but none were loaded")
    wanted = {_normalise(taxon) for taxon in taxa}
    found: set[str] = set()
    keys: set[str] = set()
    for organism_key, lineage in labels.items():
        matched = wanted & {_normalise(rank) for rank in lineage.split("|")}
        if matched:
            found |= matched
            keys.add(genome_key(organism_key))
    for taxon in taxa:
        if _normalise(taxon) not in found:
            logger.warning("Reference mask taxon %r matches no genome in the taxonomy labels", taxon)
    return frozenset(keys)


def _normalise(name: str) -> str:
    return name.strip().casefold()
