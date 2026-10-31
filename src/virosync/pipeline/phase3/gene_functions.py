"""Gene-function labels for the genes inside accepted EVEs.

A label names a gene with the vocabulary of Bellas and Sommaruga (2026, bioRxiv,
doi:10.64898/2026.06.19.733378, Fig. 4). Catalogue-mapped interior labels can supply compatible context to the hallmark
score; they never create core groups, seeds or boundaries. The evidence for one gene comes in this
order: its best curated ViroSync model (or a renamed OG-backed model whose top Pfam
domain holds in most sampled members), then its strongest Pfam hit, then its
strongest VOG hit. The rules are ported from the real-synteny figure of the benchmark
repository, together with the VOG look-alikes that must stay unlabelled.
"""

from __future__ import annotations

import csv
import io
import logging
from collections import defaultdict
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple, TypeAlias

import pyhmmer

from virosync.utils.atomic_write import atomic_write

logger = logging.getLogger(__name__)

__all__ = [
    "DESCRIBED_MODEL_PREFIXES",
    "DESCRIPTION_RULES",
    "DOMAIN_LABELS",
    "GENE_FUNCTION_COLUMNS",
    "GENE_FUNCTION_LENGTH_EXCLUSION_COLUMNS",
    "LABEL_MODULES",
    "MAX_MODEL_EVALUE",
    "MAX_PFAM_EVALUE",
    "MAX_PROFILE_PROTEIN_LENGTH",
    "MAX_VOG_EVALUE",
    "MODEL_LABELS",
    "PROFILE_COLUMNS",
    "VOG_DESCRIPTION_RULES",
    "GeneFunction",
    "GeneFunctionInputs",
    "GeneFunctionLengthExclusion",
    "GeneFunctionResources",
    "GeneFunctionResult",
    "ModelAnnotation",
    "ModelHit",
    "Profile",
    "ProfileHit",
    "find_gene_function_resources",
    "label_core_genes",
    "label_gene",
    "read_gene_function_profiles",
    "read_gene_functions",
    "read_model_annotations",
    "vog_label",
    "write_gene_function_profiles",
    "write_gene_function_length_exclusions",
    "write_gene_functions",
]

KeywordRules: TypeAlias = tuple[tuple[tuple[str, ...], str], ...]

# Cut-offs of Bellas and Sommaruga (2026): Pfam 1e-5, VOG 1e-9. ViroSync model hits
# take the Pfam cut-off.
MAX_MODEL_EVALUE = 1e-5
MAX_PFAM_EVALUE = 1e-5
MAX_VOG_EVALUE = 1e-9
MAX_PROFILE_PROTEIN_LENGTH = 100_000

# A renamed OG-backed model names a gene only when its top Pfam domain occurs in at
# least this fraction of the sampled members of the model.
OG_MIN_DOMAIN_FRACTION = 0.5
OG_BACKED_SOURCES = ("gvclass_existing_og_to_vs", "OG_backed_relaxed_low_overlap_")

PFAM_DATABASE = "Pfam"
VOG_DATABASE = "VOG"
# The first word of a profile's source string names its database.
_SOURCE_DATABASES = {"Pfam": PFAM_DATABASE, "VOGDB": VOG_DATABASE}

PROFILE_COLUMNS = ("profile", "source", "label", "module", "max_evalue", "description")
GENE_FUNCTION_COLUMNS = ("eve_id", "porf_id", "label", "module", "evidence", "evalue")
GENE_FUNCTION_LENGTH_EXCLUSION_COLUMNS = ("eve_id", "porf_id", "status", "reason", "length_aa", "limit_aa")
_MODEL_HIT_COLUMNS = ("query", "target", "score", "evalue")
_MODEL_ANNOTATION_COLUMNS = (
    "model_name",
    "source",
    "description",
    "majority_annotation",
    "pfam_top_domains",
    "pfam_n_sampled",
)

# Curated models whose name fixes their function; their descriptions do not.
MODEL_LABELS = {
    "GVOGm0003": "MCP",
    "gamadvirusMCP": "MCP",
    "GVOGm0760": "ATPase",
    "ATPase": "ATPase",
    "GVOGm0054": "PolB",
    "GVOGm0023": "RNAPL",
    "GVOGm0890": "VLTF3",
    "VLTF3": "VLTF3",
    "VLTF2": "VLTF2",
    "HUH": "HUH",
    "HUHlong": "HUH",
}
# Curated model families that state their function in their description. Penton
# comes first, because penton models are annotated as minor-capsid-like.
DESCRIBED_MODEL_PREFIXES = ("Mirus_", "VP_", "PLV_", "plv_", "EVE_CAPSID")
DESCRIPTION_RULES: KeywordRules = (
    (("penton", "minor capsid"), "mCP"),
    (("major capsid", "coat protein", "capsid protein"), "MCP"),
    (("packaging atpase", "ftsk-hera"), "ATPase"),
    (("protease",), "Protease"),
    (("portal",), "Portal"),
    (("terminase",), "Terminase"),
    (("triplex",), "Triplex"),
)
# Pfam families of the vocabulary, for genes without a curated model.
DOMAIN_LABELS = {
    "Capsid_NCLDV": "MCP",
    "MCP_V20_C": "MCP",
    "DJR-like_dom": "MCP",
    "Herpes_MCP": "MCP",
    "V18-19": "mCP",
    "Pox_A32": "ATPase",
    "Peptidase_C5": "Protease",
    "DNA_pol_B": "PolB",
    "DNA_pol_B_2": "pPolB",
    "DNA_pol_A": "PolA",
    "RNA_helicase": "SF3",
    "D5_N": "SF3",
    "Viral_helicase1": "SF1",
    "UvrD_C": "SF1",
    "PIF1": "SF1",
    "Prim-Pol": "AEP",
    "PrimPol_C": "AEP",
    "Phage_integrase": "YR",
    "rve": "RVE",
    "rve_2": "RVE",
    "rve_3": "RVE",
    "Integrase_H2C2": "RVE",
    "Pox_VLTF3": "VLTF3",
    "Replitron_HUH": "HUH",
    "Rol_Rep_N": "HUH",
    "Gemini_AL1": "HUH",
    "DNA_methylase": "MTase",
    "N6_N4_Mtase": "MTase",
    "MethyltransfD12": "MTase",
    "N6_Mtase": "MTase",
    "Tet_JBP": "Tet",
    "SLATT_1": "SLATT",
    "SLATT_2": "SLATT",
    "SLATT_3": "SLATT",
    "SLATT_4": "SLATT",
    "Patatin": "Lipase",
    "Lipase_3": "Lipase",
    "TFIIB": "TFIIB",
}
# VOG consensus descriptions of the vocabulary, matched in order. The phrases
# avoid look-alikes found in VOGDB: "DNA polymerase sliding clamp" (PCNA),
# "DNA-packaging protein FI" (a lambda chaperone), "helicase-primase primase
# subunit" (herpesvirus), "TFIIB: glycosyltransferase" and RecT or UvsX
# "recombinase", which does not integrate a genome. An integrase that names no
# type is "Int", as in Bellas and Sommaruga (2026, Fig. 2).
VOG_DESCRIPTION_RULES: KeywordRules = (
    (("minor capsid", "penton"), "mCP"),
    (("major capsid",), "MCP"),
    (("packaging atpase", "opg160"), "ATPase"),
    (("portal protein",), "Portal"),
    (("terminase",), "Terminase"),
    (("dna-directed dna polymerase", "dna polymerase type-b"), "PolB"),
    (("d5-like helicase",), "SF3"),
    (
        ("tyrosine-type recombinase", "recombinase, tyrosine", "tyrosine recombinase"),
        "YR",
    ),
    (("integrase", "site-specific recombinase", "site specific recombinase"), "Int"),
    (("late gene transcription factor 2", "vltf-2"), "VLTF2"),
    (("late gene transcription factor 3", "vltf-3"), "VLTF3"),
    (("initiation factor tfiib", "initiation factor iib"), "TFIIB"),
    (
        ("dna methyltransferase", "dna methylase", "type ii methyltransferase"),
        "MTase",
    ),
)

_MODULE_LABELS = {
    "capsid": ("MCP", "mCP", "ATPase", "Protease", "Portal", "Terminase", "Triplex"),
    "replication": ("PolB", "pPolB", "PolA", "SF1", "SF3", "AEP"),
    "integration": ("YR", "RVE", "Int"),
    "mriya": ("VLTF2", "VLTF3", "HUH"),
    "other": ("MTase", "Tet", "SLATT", "Lipase"),
    "hallmark": ("RNAPL", "TFIIB"),
}
LABEL_MODULES = {label: module for module, labels in _MODULE_LABELS.items() for label in labels}


@dataclass(frozen=True, slots=True)
class Profile:
    """One profile of the gene-function profile set and the label it gives a gene."""

    profile: str
    source: str
    label: str
    module: str
    max_evalue: float
    description: str

    @property
    def database(self) -> str:
        """Database of the profile, "Pfam" or "VOG", from the first word of its source."""
        return _database_of(self.source)


@dataclass(frozen=True, slots=True)
class GeneFunction:
    """The function label of one gene inside an accepted EVE and the evidence for it."""

    eve_id: str
    porf_id: str
    label: str
    module: str
    evidence: str
    evalue: float


@dataclass(frozen=True, slots=True)
class GeneFunctionLengthExclusion:
    """A core protein whose processed length prevents a profile search in one EVE."""

    eve_id: str
    porf_id: str
    length_aa: int


@dataclass(frozen=True, slots=True)
class GeneFunctionResult:
    """Keep positive function evidence separate from missing profile assessment."""

    labels: list[GeneFunction]
    length_exclusions: list[GeneFunctionLengthExclusion]


@dataclass(frozen=True, slots=True)
class ModelHit:
    """One Phase 1 HMM hit of a gene."""

    model: str
    score: float
    evalue: float


@dataclass(frozen=True, slots=True)
class ProfileHit:
    """One hit of a gene against the profile set.

    For a Pfam profile, ``evalue`` is the i-Evalue of the best domain. For a VOG
    profile, it is the full-sequence E-value.
    """

    profile: str
    evalue: float


@dataclass(frozen=True, slots=True)
class ModelAnnotation:
    """The fields of a model-annotation row that decide a gene label."""

    source: str
    description: str
    majority_annotation: str
    pfam_top_domains: str
    pfam_n_sampled: int


@dataclass(frozen=True, slots=True)
class GeneFunctionInputs:
    """Where a run finds the evidence that labels its genes.

    ``hmm_database`` is the run's ``combined.hmm``. The gene-function profile set and
    the model-annotation table sit next to it, and ``model_hits`` is the Phase 1
    validated HMM hit table of the run.
    """

    hmm_database: Path
    model_hits: Path
    threads: int = 1


@dataclass(frozen=True, slots=True)
class GeneFunctionResources:
    """The three bundle files that label a gene."""

    profiles_hmm: Path
    profiles_table: Path
    model_annotations: Path


class _Choice(NamedTuple):
    """A label with the evidence and E-value behind it."""

    label: str
    evidence: str
    evalue: float


def label_gene(
    eve_id: str,
    porf_id: str,
    *,
    model_hits: Sequence[ModelHit],
    profile_hits: Sequence[ProfileHit],
    model_annotations: Mapping[str, ModelAnnotation],
    profiles: Mapping[str, Profile],
) -> GeneFunction | None:
    """Pick the function label of one gene, or return None when the evidence is too weak.

    The first source that gives a label wins: the gene's best-scoring Phase 1 HMM hit
    at E <= 1e-5 (a curated model, or a trusted OG-backed model), then its Pfam hits at
    or below each profile's cut-off, then its VOG hits at or below each profile's
    cut-off. Within Pfam and within VOG, the strongest hit goes first.

    Args:
        eve_id: Identifier of the EVE that holds the gene.
        porf_id: Identifier of the gene.
        model_hits: All Phase 1 HMM hits of the gene.
        profile_hits: All hits of the gene against the gene-function profile set.
        model_annotations: The model-annotation table, by model name.
        profiles: The gene-function profile set, by profile name.
    """
    choice = _model_choice(_best_model_hit(model_hits), model_annotations)
    if choice is None:
        choice = _profile_choice(profile_hits, profiles, PFAM_DATABASE)
    if choice is None:
        choice = _profile_choice(profile_hits, profiles, VOG_DATABASE)
    if choice is None:
        return None
    return GeneFunction(
        eve_id,
        porf_id,
        choice.label,
        LABEL_MODULES[choice.label],
        choice.evidence,
        choice.evalue,
    )


def vog_label(description: str) -> str:
    """Return the vocabulary label of a VOG consensus description, or an empty string."""
    return _keyword_label(description, VOG_DESCRIPTION_RULES)


def read_gene_function_profiles(path: Path) -> dict[str, Profile]:
    """Read ``gene_function_labels.tsv`` into profiles by name, in file order.

    Raises:
        ValueError: If the columns differ, a label is outside the vocabulary, a module
            disagrees with its label, a source names neither Pfam nor VOGDB, or a
            profile repeats.
    """
    profiles: dict[str, Profile] = {}
    for line_number, row in _read_rows(path, PROFILE_COLUMNS):
        label = row["label"]
        if label not in LABEL_MODULES:
            raise ValueError(f"{path.name} line {line_number}: label {label!r} is not in the gene vocabulary")
        module = row["module"]
        if module != LABEL_MODULES[label]:
            raise ValueError(
                f"{path.name} line {line_number}: module {module!r} does not match "
                f"label {label!r} (module {LABEL_MODULES[label]!r})"
            )
        profile = Profile(
            row["profile"],
            row["source"],
            label,
            module,
            float(row["max_evalue"]),
            row["description"],
        )
        _database_of(profile.source)
        if profile.profile in profiles:
            raise ValueError(f"{path.name} line {line_number}: table repeats profile {profile.profile!r}")
        profiles[profile.profile] = profile
    return profiles


def write_gene_function_profiles(profiles: Iterable[Profile], path: Path) -> None:
    """Write profiles to ``gene_function_labels.tsv`` through a temporary file."""
    rows = ((p.profile, p.source, p.label, p.module, f"{p.max_evalue:g}", p.description) for p in profiles)
    _write_rows(path, PROFILE_COLUMNS, rows)


def read_gene_functions(path: Path) -> list[GeneFunction]:
    """Read ``gene_functions.tsv`` into gene-function rows.

    Raises:
        ValueError: If the columns differ from the documented ones.
    """
    return [
        GeneFunction(
            row["eve_id"],
            row["porf_id"],
            row["label"],
            row["module"],
            row["evidence"],
            float(row["evalue"]),
        )
        for _, row in _read_rows(path, GENE_FUNCTION_COLUMNS)
    ]


def write_gene_functions(rows: Iterable[GeneFunction], path: Path) -> None:
    """Write gene-function rows to ``gene_functions.tsv`` through a temporary file."""
    cells = ((r.eve_id, r.porf_id, r.label, r.module, r.evidence, f"{r.evalue:.3g}") for r in rows)
    _write_rows(path, GENE_FUNCTION_COLUMNS, cells)


def write_gene_function_length_exclusions(rows: Iterable[GeneFunctionLengthExclusion], path: Path) -> None:
    """Write length exclusions, including a header for an empty selection, atomically."""
    cells = (
        (
            row.eve_id,
            row.porf_id,
            "unsearched",
            "sequence_length_limit",
            str(row.length_aa),
            str(MAX_PROFILE_PROTEIN_LENGTH),
        )
        for row in rows
    )
    _write_rows(path, GENE_FUNCTION_LENGTH_EXCLUSION_COLUMNS, cells)


def read_model_annotations(path: Path) -> dict[str, ModelAnnotation]:
    """Read the fields of ``model_annotations_with_interpro.tsv`` that decide a label.

    A missing, blank or malformed ``pfam_n_sampled`` reads as 0, which leaves an OG-backed
    model untrusted, so one bad row never stops the labelling of other genes.

    Raises:
        ValueError: If one of the needed columns is missing.
    """
    annotations: dict[str, ModelAnnotation] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        missing = [name for name in _MODEL_ANNOTATION_COLUMNS if name not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path.name} lacks columns {missing}")
        for row in reader:
            n_sampled = (row["pfam_n_sampled"] or "").strip()
            annotations[row["model_name"]] = ModelAnnotation(
                row["source"],
                row["description"],
                row["majority_annotation"],
                row["pfam_top_domains"],
                int(n_sampled) if n_sampled.isdecimal() else 0,
            )
    return annotations


def find_gene_function_resources(hmm_database: Path) -> GeneFunctionResources | None:
    """Return the gene-function files next to *hmm_database*, or None for a bundle before v1.2.0.

    Raises:
        FileNotFoundError: If exactly one of the two gene-function files exists. A bundle
            of v1.2.0 or later carries both, so a lone file marks a corrupt bundle.
    """
    models_dir = Path(hmm_database).parent
    resources = GeneFunctionResources(
        models_dir / "gene_function.hmm",
        models_dir / "gene_function_labels.tsv",
        models_dir / "model_annotations_with_interpro.tsv",
    )
    present = [path for path in (resources.profiles_hmm, resources.profiles_table) if path.is_file()]
    if not present:
        return None
    if len(present) == 1:
        missing = resources.profiles_table if present[0] == resources.profiles_hmm else resources.profiles_hmm
        raise FileNotFoundError(f"gene-function files are incomplete; {missing} is missing next to {present[0].name}")
    return resources


def label_core_genes(
    genes: Sequence[tuple[str, str]],
    resources: GeneFunctionResources,
    *,
    model_hits: Path,
    proteomes: Sequence[Path],
    threads: int,
) -> GeneFunctionResult:
    """Return function labels and length exclusions for ``(eve_id, porf_id)`` pairs.

    The rows keep the order of *genes*. The profile search reports hits down to E 1e-5
    and each profile's own cut-off then applies. E-values follow the hmmscan convention:
    the database size is the number of profiles in the library, so a gene's E-value and
    label do not depend on how many other genes the run searches. A gene without a
    sequence in *proteomes* stays unlabelled and is logged. A missing *model_hits* table
    means no gene has a Phase 1 model hit. Proteins above the profile search's length
    limit retain any curated model label and are recorded separately as unsearched.
    """
    profiles = read_gene_function_profiles(resources.profiles_table)
    annotations = read_model_annotations(resources.model_annotations)
    wanted = {porf_id for _, porf_id in genes}
    model_hits_by_gene = _read_model_hits(model_hits, wanted)
    proteins = _read_proteins(proteomes, wanted)
    supported = pyhmmer.easel.DigitalSequenceBlock(proteins.alphabet)
    excluded_lengths: dict[str, int] = {}
    # Resolve file precedence and remove terminal stops before measuring length.
    for protein in proteins:
        if len(protein) > MAX_PROFILE_PROTEIN_LENGTH:
            excluded_lengths[protein.name.decode()] = len(protein)
        else:
            supported.append(protein)
    if excluded_lengths:
        logger.warning(
            "Gene-function profile search excluded %d proteins longer than %d amino acids",
            len(excluded_lengths),
            MAX_PROFILE_PROTEIN_LENGTH,
        )
    profile_hits_by_gene = _search_profiles(resources.profiles_hmm, supported, profiles, threads)
    labelled = (
        label_gene(
            eve_id,
            porf_id,
            model_hits=model_hits_by_gene.get(porf_id, []),
            profile_hits=profile_hits_by_gene.get(porf_id, []),
            model_annotations=annotations,
            profiles=profiles,
        )
        for eve_id, porf_id in genes
    )
    return GeneFunctionResult(
        labels=[row for row in labelled if row is not None],
        length_exclusions=[
            GeneFunctionLengthExclusion(eve_id, porf_id, excluded_lengths[porf_id])
            for eve_id, porf_id in genes
            if porf_id in excluded_lengths
        ],
    )


def _database_of(source: str) -> str:
    words = source.split()
    database = _SOURCE_DATABASES.get(words[0]) if words else None
    if database is None:
        raise ValueError(f"profile source {source!r} must start with one of {sorted(_SOURCE_DATABASES)}")
    return database


def _keyword_label(text: str, rules: KeywordRules) -> str:
    lowered = text.lower()
    return next(
        (label for keywords, label in rules if any(keyword in lowered for keyword in keywords)),
        "",
    )


def _best_model_hit(model_hits: Sequence[ModelHit]) -> ModelHit | None:
    eligible = [hit for hit in model_hits if hit.evalue <= MAX_MODEL_EVALUE]
    if not eligible:
        return None
    return max(eligible, key=lambda hit: hit.score)


def _model_choice(
    hit: ModelHit | None,
    annotations: Mapping[str, ModelAnnotation],
) -> _Choice | None:
    if hit is None:
        return None
    if hit.model in MODEL_LABELS:
        return _Choice(MODEL_LABELS[hit.model], hit.model, hit.evalue)
    annotation = annotations.get(hit.model)
    if annotation is None:
        return None
    if hit.model.startswith(DESCRIBED_MODEL_PREFIXES):
        label = _keyword_label(f"{annotation.description} {annotation.majority_annotation}", DESCRIPTION_RULES)
        return _Choice(label, hit.model, hit.evalue) if label else None
    if annotation.source.startswith(OG_BACKED_SOURCES):
        return _trusted_og_choice(hit, annotation)
    return None


def _trusted_og_choice(hit: ModelHit, annotation: ModelAnnotation) -> _Choice | None:
    """Label a renamed OG-backed model only by a vocabulary domain most members share.

    Metadata that cannot be trusted (no sampled members, a malformed count, or more
    members with the domain than were sampled) gives no label.
    """
    n_sampled = annotation.pfam_n_sampled
    if n_sampled <= 0:
        return None
    top_domain = _top_domain(annotation.pfam_top_domains)
    if top_domain is None:
        return None
    domain, count = top_domain
    label = DOMAIN_LABELS.get(domain)
    if label is None or count > n_sampled or count / n_sampled < OG_MIN_DOMAIN_FRACTION:
        return None
    evidence = f"{hit.model} (OG; {domain} in {count}/{n_sampled} members)"
    return _Choice(label, evidence, hit.evalue)


def _top_domain(pfam_top_domains: str) -> tuple[str, int] | None:
    """Return the domain with the highest member count from ``Name:count;Name:count``.

    A listing with any entry that has no name or no non-negative integer count is
    malformed and gives None.
    """
    entries = [entry.rpartition(":") for entry in pfam_top_domains.split(";") if entry != ""]
    if not entries or not all(name and count.isdecimal() for name, _, count in entries):
        return None
    return max(((name, int(count)) for name, _, count in entries), key=lambda item: item[1])


def _profile_choice(
    hits: Sequence[ProfileHit],
    profiles: Mapping[str, Profile],
    database: str,
) -> _Choice | None:
    counted = [hit for hit in hits if _is_counted(hit, profiles, database)]
    if not counted:
        return None
    strongest = min(counted, key=lambda hit: (hit.evalue, hit.profile))
    return _Choice(
        profiles[strongest.profile].label,
        f"{database}:{strongest.profile}",
        strongest.evalue,
    )


def _is_counted(hit: ProfileHit, profiles: Mapping[str, Profile], database: str) -> bool:
    profile = profiles.get(hit.profile)
    return profile is not None and profile.database == database and hit.evalue <= profile.max_evalue


def _read_rows(path: Path, columns: Sequence[str]) -> Iterable[tuple[int, dict[str, str]]]:
    """Yield ``(line number, row)`` for a TSV whose header is exactly *columns*."""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        if reader.fieldnames != list(columns):
            raise ValueError(f"{path.name} columns {reader.fieldnames} differ from {list(columns)}")
        for row in reader:
            yield reader.line_num, row


def _write_rows(path: Path, columns: Sequence[str], rows: Iterable[Sequence[str]]) -> None:
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter="\t", quoting=csv.QUOTE_NONE, lineterminator="\n")
    writer.writerow(columns)
    writer.writerows(rows)
    atomic_write(path, buffer.getvalue())


def _read_model_hits(path: Path, wanted: Collection[str]) -> dict[str, list[ModelHit]]:
    """Group the Phase 1 validated HMM hits of the *wanted* genes by gene."""
    hits: dict[str, list[ModelHit]] = defaultdict(list)
    if not path.is_file():
        return hits
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        missing = [name for name in _MODEL_HIT_COLUMNS if name not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path.name} lacks columns {missing}")
        for row in reader:
            if row["query"] in wanted:
                hits[row["query"]].append(ModelHit(row["target"], float(row["score"]), float(row["evalue"])))
    return hits


def _read_proteins(proteomes: Sequence[Path], wanted: Collection[str]) -> pyhmmer.easel.DigitalSequenceBlock:
    """Read the *wanted* proteins from FASTA files of Prodigal proteins, without stop symbols.

    The first file that holds a protein supplies it.
    """
    alphabet = pyhmmer.easel.Alphabet.amino()
    block = pyhmmer.easel.DigitalSequenceBlock(alphabet)
    found: set[str] = set()
    for proteome in proteomes:
        # Phase 1 writes an empty confirmed-rescue file when it confirms no rescue.
        # Any other file that pyhmmer cannot read still raises.
        if Path(proteome).stat().st_size == 0:
            continue
        with pyhmmer.easel.SequenceFile(proteome, format="fasta") as handle:
            for record in handle:
                name = record.name.decode()
                sequence = record.sequence.rstrip("*")
                if name in wanted and name not in found and sequence:
                    block.append(pyhmmer.easel.TextSequence(name=record.name, sequence=sequence).digitize(alphabet))
                    found.add(name)
    missing = set(wanted) - found
    if missing:
        logger.warning(
            "%d core genes have no sequence in the proteome files and stay unlabelled: %s",
            len(missing),
            sorted(missing)[:5],
        )
    return block


def _search_profiles(
    profiles_hmm: Path,
    proteins: pyhmmer.easel.DigitalSequenceBlock,
    profiles: Mapping[str, Profile],
    threads: int,
) -> dict[str, list[ProfileHit]]:
    """Search supported proteins against the profile set and group the hits by gene.

    The effective database size ``Z`` (and ``domZ``) is the number of profiles, as in
    ``hmmscan``, instead of HMMER's ``hmmsearch`` default of the number of searched genes.
    """
    hits: dict[str, list[ProfileHit]] = defaultdict(list)
    if len(proteins) == 0 or not profiles:
        return hits
    with pyhmmer.plan7.HMMFile(profiles_hmm) as library:
        database_size = len(profiles)
        for top_hits in pyhmmer.hmmsearch(
            library, proteins, cpus=threads, E=MAX_PFAM_EVALUE, Z=database_size, domZ=database_size
        ):
            profile = profiles[top_hits.query.name.decode()]
            for hit in top_hits.included:
                hits[hit.name.decode()].append(ProfileHit(profile.profile, _profile_evalue(profile, hit)))
    return hits


def _profile_evalue(profile: Profile, hit: pyhmmer.plan7.Hit) -> float:
    """Return the i-Evalue of the best domain for a Pfam profile, else the sequence E-value."""
    if profile.database == PFAM_DATABASE:
        return min(domain.i_evalue for domain in hit.domains)
    return hit.evalue
