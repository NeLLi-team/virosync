"""Resolve proteins that match more than one ViroSync HMM with Pfam domains."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import pyhmmer
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMMFile

from virosync.pipeline.phase1.hhg_seeding import HMMHit
from virosync.utils.atomic_write import atomic_write_context


CRESS_REP_SCOPE = "CRESS_REP"
CRESS_REP_DISCRIMINATING_DOMAIN = "Gemini_AL1"


@dataclass(frozen=True)
class ModelPfamAnnotation:
    """Pfam signature and source scope for one ViroSync model."""

    signature: frozenset[str]
    source_scope: str


@dataclass(frozen=True)
class PfamArbitrationRecord:
    """One decision for a protein with multiple candidate models."""

    protein: str
    candidates: tuple[tuple[str, float], ...]
    original_model: str
    observed_domains: tuple[str, ...]
    compatible_models: tuple[str, ...]
    final_model: str | None
    outcome: str


def _hit_rank(hit: HMMHit) -> tuple[float, str]:
    """Return a deterministic best-first rank for an HMM hit."""

    return (-hit.score, hit.target_name)


def _best_candidate_hits(hits: list[HMMHit]) -> dict[str, dict[str, HMMHit]]:
    best_by_protein: dict[str, dict[str, HMMHit]] = {}
    for hit in hits:
        best_by_protein.setdefault(hit.query_name, {})[hit.target_name] = hit
    return best_by_protein


def ambiguous_proteins(hits: list[HMMHit]) -> set[str]:
    """Return proteins hit by at least two distinct ViroSync models."""

    return {protein for protein, candidates in _best_candidate_hits(hits).items() if len(candidates) >= 2}


def load_model_pfam_annotations(path: Path) -> dict[str, ModelPfamAnnotation]:
    """Load the enriched model annotation columns needed for arbitration."""

    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"model_name", "pfam_signature", "source_scope"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"model annotation table is missing Pfam columns: {sorted(missing)}")
        annotations = {}
        for row in reader:
            model = row["model_name"].strip()
            if not model:
                continue
            annotations[model] = ModelPfamAnnotation(
                signature=frozenset(domain for domain in row["pfam_signature"].split(";") if domain),
                source_scope=row["source_scope"].strip(),
            )
    return annotations


def scan_pfam_domains(
    proteome_path: Path,
    pfam_hmm_path: Path,
    proteins: set[str],
    threads: int,
) -> dict[str, set[str]]:
    """Scan selected proteins with Pfam gathering thresholds."""

    if not proteins:
        return {}
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences = []
    with SequenceFile(proteome_path, digital=True, alphabet=alphabet) as handle:
        for sequence in handle:
            name = sequence.name.decode()
            if name in proteins:
                sequences.append(sequence)
    with HMMFile(pfam_hmm_path) as handle:
        hmms = list(handle)
    domains_by_protein = {protein: set() for protein in proteins}
    for hmm, top_hits in zip(
        hmms,
        pyhmmer.hmmsearch(
            hmms,
            sequences,
            cpus=threads,
            parallel="targets",
            E=float("inf"),
        ),
    ):
        domain_name = hmm.name.decode()
        sequence_cutoff, domain_cutoff = hmm.cutoffs.gathering
        for hit in top_hits:
            if hit.score >= sequence_cutoff and any(
                domain.score >= domain_cutoff for domain in hit.domains
            ):
                domains_by_protein[hit.name.decode()].add(domain_name)
    return domains_by_protein


def arbitrate_hits(
    hits: list[HMMHit],
    domains_by_protein: dict[str, set[str]],
    annotations: dict[str, ModelPfamAnnotation],
) -> tuple[list[HMMHit], list[PfamArbitrationRecord]]:
    """Apply the handoff's Pfam arbitration rules to ambiguous proteins."""

    best_by_protein = _best_candidate_hits(hits)
    ambiguous = {protein for protein, candidates in best_by_protein.items() if len(candidates) >= 2}
    records = []
    selected_hits: dict[str, HMMHit] = {}

    for protein in sorted(ambiguous):
        candidates = best_by_protein[protein]
        original_hit = min(candidates.values(), key=_hit_rank)
        observed = domains_by_protein.get(protein, set())
        compatible = []
        for model in sorted(candidates):
            annotation = annotations.get(
                model,
                ModelPfamAnnotation(frozenset(), ""),
            )
            required_domains = (
                {CRESS_REP_DISCRIMINATING_DOMAIN}
                if annotation.source_scope == CRESS_REP_SCOPE
                else annotation.signature
            )
            if observed & required_domains:
                compatible.append(model)

        if not observed:
            final_model = original_hit.target_name
            outcome = "unresolved_no_domain"
        elif len(compatible) == 1:
            final_model = compatible[0]
            outcome = "confirmed" if final_model == original_hit.target_name else "reassigned"
        elif len(compatible) > 1:
            final_model = original_hit.target_name
            outcome = "unresolved_shared_domain"
        else:
            final_model = None
            outcome = "contradicted"

        if final_model is not None:
            selected_hits[protein] = candidates[final_model]
        records.append(
            PfamArbitrationRecord(
                protein=protein,
                candidates=tuple((model, candidates[model].score) for model in sorted(candidates)),
                original_model=original_hit.target_name,
                observed_domains=tuple(sorted(observed)),
                compatible_models=tuple(compatible),
                final_model=final_model,
                outcome=outcome,
            )
        )

    retained = [hit for hit in hits if hit.query_name not in ambiguous or selected_hits.get(hit.query_name) is hit]
    return retained, records


def write_pfam_arbitration(
    records: list[PfamArbitrationRecord],
    output_path: Path,
) -> None:
    """Write the per-protein Pfam decision audit."""

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with atomic_write_context(output_path, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "protein",
                "candidates",
                "original_model",
                "observed_domains",
                "compatible_models",
                "final_model",
                "outcome",
            ]
        )
        for record in records:
            writer.writerow(
                [
                    record.protein,
                    ";".join(f"{model}:{score:.3f}" for model, score in record.candidates),
                    record.original_model,
                    ";".join(record.observed_domains),
                    ";".join(record.compatible_models),
                    record.final_model or "",
                    record.outcome,
                ]
            )


def run_pfam_arbitration(
    hits: list[HMMHit],
    proteins: set[str],
    proteome_path: Path,
    pfam_hmm_path: Path,
    model_annotations_path: Path,
    output_path: Path,
    threads: int,
) -> list[HMMHit]:
    """Scan ambiguous proteins, arbitrate their candidates, and write an audit."""

    annotations = load_model_pfam_annotations(model_annotations_path)
    domains = scan_pfam_domains(
        proteome_path,
        pfam_hmm_path,
        proteins,
        threads,
    )
    retained, records = arbitrate_hits(hits, domains, annotations)
    write_pfam_arbitration(records, output_path)
    return retained
