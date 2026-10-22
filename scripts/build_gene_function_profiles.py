#!/usr/bin/env python3
"""Build the gene-function profile set of a ViroSync v1.2.0 runtime bundle.

The set holds the Pfam families of the gene vocabulary (Bellas and Sommaruga 2026,
Fig. 4) and the VOGDB groups whose consensus description names a gene of that
vocabulary, in any virus lineage. The script writes ``gene_function.hmm`` (HMMER3
text) and ``gene_function_labels.tsv`` (one row per profile) into ``--output-dir``.
It streams both HMM files, so the 4 GB VOGDB library never sits in memory.

Example
-------
    pixi run python scripts/build_gene_function_profiles.py \
        --pfam Pfam-A.hmm --pfam-release 38.2 \
        --vog vogdb.hmm --vog-annotations vog.annotations.tsv \
        --vog-release latest-2025-09-30 --output-dir resources/gene_function
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from collections import Counter
from collections.abc import Collection, Sequence
from dataclasses import dataclass
from pathlib import Path

import pyhmmer

# Keep the documented ``pixi run python scripts/...`` invocation independent of
# whether the project has been installed in editable mode.
_SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

from virosync.pipeline.phase3.gene_functions import (  # noqa: E402  # needs the src path above
    DOMAIN_LABELS,
    LABEL_MODULES,
    MAX_PFAM_EVALUE,
    MAX_VOG_EVALUE,
    Profile,
    vog_label,
    write_gene_function_profiles,
)
from virosync.utils.atomic_write import atomic_write_context  # noqa: E402  # needs the src path above

HMM_NAME = "gene_function.hmm"
LABELS_NAME = "gene_function_labels.tsv"
# VOGDB writes every description that comes from RefSeq with this prefix.
REFSEQ_PREFIX = "REFSEQ "
VOG_ANNOTATION_HEADER = ["#GroupName", "ProteinCount", "SpeciesCount", "FunctionalCategory"]

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class BuildRequest:
    """Inputs and output folder of one profile-set build."""

    pfam: Path
    vog: Path
    vog_annotations: Path
    pfam_release: str
    vog_release: str
    output_dir: Path


def build_gene_function_profiles(request: BuildRequest) -> list[Profile]:
    """Write the HMM library and the labels table, and return the selected profiles.

    The library lists the Pfam profiles first, then the VOG profiles, each in the
    order of its source file; the labels table has the same order.

    Raises:
        ValueError: If a Pfam family of the vocabulary or a selected VOG has no model
            in its library, or the VOG annotation table has another layout.
    """
    vog_profiles = select_vog_profiles(request.vog_annotations, request.vog_release)
    pfam_hmms = _stream_selected(request.pfam, DOMAIN_LABELS)
    _require_all(DOMAIN_LABELS, pfam_hmms, request.pfam)
    vog_hmms = _stream_selected(request.vog, vog_profiles)
    _require_all(vog_profiles, vog_hmms, request.vog)

    profiles = [
        *(_pfam_profile(hmm, request.pfam_release) for hmm in pfam_hmms.values()),
        *(vog_profiles[name] for name in vog_hmms),
    ]
    _write_library([*pfam_hmms.values(), *vog_hmms.values()], request.output_dir / HMM_NAME)
    write_gene_function_profiles(profiles, request.output_dir / LABELS_NAME)
    return profiles


def select_vog_profiles(annotations: Path, release: str) -> dict[str, Profile]:
    """Return the VOG profiles whose consensus description gets a vocabulary label.

    Raises:
        ValueError: If the annotation table does not have the VOGDB header.
    """
    selected: dict[str, Profile] = {}
    with annotations.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = next(reader, [])
        if header[:4] != VOG_ANNOTATION_HEADER:
            raise ValueError(f"{annotations.name} is not a VOGDB annotation table; its header is {header}")
        for vog_id, *_, description in reader:
            description = description.removeprefix(REFSEQ_PREFIX).strip()
            label = vog_label(description)
            if label:
                selected[vog_id] = Profile(
                    vog_id,
                    f"VOGDB {release}",
                    label,
                    LABEL_MODULES[label],
                    MAX_VOG_EVALUE,
                    description,
                )
    return selected


def main(argv: Sequence[str] | None = None) -> int:
    """Parse arguments, build the profile set, print its counts, and return the exit code."""
    request = _parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    try:
        profiles = build_gene_function_profiles(request)
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print("\n".join(_summary_lines(profiles)))
    return 0


def _parse_args(argv: Sequence[str] | None) -> BuildRequest:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pfam", type=Path, required=True, help="Pfam-A.hmm, the HMMER3 text library")
    parser.add_argument("--vog", type=Path, required=True, help="VOGDB HMMER3 text library with all groups")
    parser.add_argument("--vog-annotations", type=Path, required=True, help="VOGDB vog.annotations.tsv")
    parser.add_argument("--pfam-release", required=True, help="Pfam release, for the source column, such as 38.2")
    parser.add_argument("--vog-release", required=True, help="VOGDB release, for the source column")
    parser.add_argument("--output-dir", type=Path, required=True, help="folder for the two output files")
    args = parser.parse_args(argv)
    return BuildRequest(
        pfam=args.pfam,
        vog=args.vog,
        vog_annotations=args.vog_annotations,
        pfam_release=args.pfam_release,
        vog_release=args.vog_release,
        output_dir=args.output_dir,
    )


def _stream_selected(library: Path, wanted: Collection[str]) -> dict[str, pyhmmer.plan7.HMM]:
    """Stream *library* and keep the profiles named in *wanted*, in file order."""
    selected: dict[str, pyhmmer.plan7.HMM] = {}
    with pyhmmer.plan7.HMMFile(library) as handle:
        for hmm in handle:
            name = hmm.name.decode()
            if name not in wanted:
                continue
            if name in selected:
                raise ValueError(f"{library.name} holds profile {name!r} twice")
            selected[name] = hmm
    return selected


def _require_all(wanted: Collection[str], found: Collection[str], library: Path) -> None:
    missing = sorted(set(wanted) - set(found))
    if missing:
        raise ValueError(f"{library.name} has no profile for {missing}")


def _pfam_profile(hmm: pyhmmer.plan7.HMM, release: str) -> Profile:
    name = hmm.name.decode()
    label = DOMAIN_LABELS[name]
    description = hmm.description.decode() if hmm.description else ""
    return Profile(name, f"Pfam {release}", label, LABEL_MODULES[label], MAX_PFAM_EVALUE, description)


def _write_library(hmms: Sequence[pyhmmer.plan7.HMM], path: Path) -> None:
    with atomic_write_context(path, mode="wb") as handle:
        for hmm in hmms:
            hmm.write(handle, binary=False)
    logger.info("wrote %d profiles to %s", len(hmms), path)


def _summary_lines(profiles: Sequence[Profile]) -> list[str]:
    by_database = Counter(profile.database for profile in profiles)
    labels_by_module: dict[str, Counter[str]] = {}
    for profile in profiles:
        labels_by_module.setdefault(profile.module, Counter())[profile.label] += 1
    lines = [f"Pfam profiles: {by_database['Pfam']}", f"VOG profiles: {by_database['VOG']}"]
    lines.extend(
        f"{module}: " + ", ".join(f"{label} {count}" for label, count in sorted(counts.items()))
        for module, counts in labels_by_module.items()
    )
    return lines


if __name__ == "__main__":
    raise SystemExit(main())
