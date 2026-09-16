"""Align nucleotide sequence to a protein profile while permitting frameshifts."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, replace
from itertools import combinations, product
from typing import Literal

import numpy as np
from Bio.Data import CodonTable
from pyhmmer import plan7

__all__ = ["CodonAlignment", "FrameshiftEvent", "align_frameshift_profile"]


_MATCH = 1
_INSERT = 2
_DELETE = 3
_STOP = -2
_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
_DNA_BASES = "ACGT"
_FRAME_STEPS = (3, 1, 2, 4, 5)
_FRAME_PENALTY_BITS = {1: np.log2(0.01), 2: np.log2(0.005)}
_STOP_PENALTY_BITS = np.log2(0.01)

_STANDARD_CODE = CodonTable.unambiguous_dna_by_name["Standard"]
_CODON_TO_AA = dict(_STANDARD_CODE.forward_table)
_STOP_CODONS = frozenset(_STANDARD_CODE.stop_codons)


@dataclass(frozen=True, slots=True)
class FrameshiftEvent:
    """One internal disruption in a codon/profile alignment."""

    kind: Literal["insertion", "deletion", "stop"]
    nt_start: int
    nt_end: int
    model_position: int
    size: int
    protein_position: int


@dataclass(frozen=True, slots=True)
class CodonAlignment:
    """Best local profile alignment in zero-based, half-open coordinates."""

    nt_start: int
    nt_end: int
    model_start: int
    model_end: int
    sequence: str
    score: float
    baseline_score: float
    events: tuple[FrameshiftEvent, ...]
    codon_spans: tuple[tuple[int, int], ...]


def align_frameshift_profile(hmm: plan7.HMM, dna: str) -> CodonAlignment | None:
    """Return the best positive local codon/profile alignment.

    Scores are uncalibrated bits. ``baseline_score`` is the best ordinary-codon
    score within the nucleotide span selected by the unrestricted alignment.
    The input sequence is expected to be uppercase and already oriented.
    """
    if hmm.M == 0 or dna == "":
        return None

    profile_scores = _profile_scores(hmm)
    emission_tables = _profile_emission_tables(profile_scores[0])
    score, endpoint, traces = _run_viterbi(
        dna,
        _FRAME_STEPS,
        profile_scores,
        emission_tables,
        keep_trace=True,
    )
    if endpoint is None or traces is None:
        return None

    alignment = _trace_alignment(dna, score, endpoint, traces)
    selected_dna = dna[alignment.nt_start : alignment.nt_end]
    baseline_score, _, _ = _run_viterbi(
        selected_dna,
        (3,),
        profile_scores,
        emission_tables,
        keep_trace=False,
    )
    return replace(alignment, baseline_score=baseline_score)


def _run_viterbi(
    dna: str,
    steps: tuple[int, ...],
    profile_scores: tuple[np.ndarray, np.ndarray, np.ndarray],
    emission_tables: dict[int, np.ndarray],
    *,
    keep_trace: bool,
) -> tuple[
    float,
    tuple[int, int, int] | None,
    tuple[np.ndarray, np.ndarray, np.ndarray] | None,
]:
    _, insert_scores, transition_scores = profile_scores
    ordinary_codons = _ordinary_codon_indices(dna)
    segment_indices = {step: _segment_indices(dna, step) for step in steps}
    nucleotide_count = len(dna)
    model_length = len(insert_scores) - 1
    negative_infinity = -np.inf

    match_previous = np.full(nucleotide_count + 1, negative_infinity)
    insert_previous = np.full(nucleotide_count + 1, negative_infinity)
    delete_previous = np.full(nucleotide_count + 1, negative_infinity)

    match_trace = None
    insert_trace = None
    delete_trace = None
    if keep_trace:
        shape = (model_length + 1, nucleotide_count + 1)
        match_trace = np.zeros(shape, dtype=np.uint8)
        insert_trace = np.zeros(shape, dtype=np.uint8)
        delete_trace = np.zeros(shape, dtype=np.uint8)

    best_score = 0.0
    endpoint: tuple[int, int, int] | None = None
    endpoint_match_code = 0
    for model_node in range(1, model_length + 1):
        transition_row = transition_scores[model_node - 1]
        delete_current, delete_codes = _delete_row(
            match_previous,
            delete_previous,
            transition_row,
        )
        match_current, match_codes, ordinary_scores, ordinary_codes = _match_row(
            match_previous,
            insert_previous,
            delete_previous,
            transition_row,
            ordinary_codons,
            segment_indices,
            {step: emission_tables[step][model_node] for step in steps},
            steps,
        )
        insert_current, insert_codes = _insert_row(
            match_current,
            insert_scores[model_node],
            transition_scores[model_node],
            ordinary_codons,
        )

        match_best, match_end = _best_match_endpoint(
            ordinary_scores,
            ordinary_codons,
        )
        if match_best > best_score:
            best_score = match_best
            endpoint = (_MATCH, model_node, match_end)
            endpoint_match_code = int(ordinary_codes[match_end])
        insert_end = int(np.argmax(insert_current))
        insert_best = float(insert_current[insert_end])
        if insert_best > best_score:
            best_score = insert_best
            endpoint = (_INSERT, model_node, insert_end)
            endpoint_match_code = 0

        if keep_trace:
            assert match_trace is not None
            assert insert_trace is not None
            assert delete_trace is not None
            match_trace[model_node] = match_codes
            insert_trace[model_node] = insert_codes
            delete_trace[model_node] = delete_codes

        match_previous = match_current
        insert_previous = insert_current
        delete_previous = delete_current

    traces = None
    if keep_trace:
        assert match_trace is not None
        assert insert_trace is not None
        assert delete_trace is not None
        if endpoint is not None and endpoint[0] == _MATCH:
            match_trace[endpoint[1], endpoint[2]] = endpoint_match_code
        traces = (match_trace, insert_trace, delete_trace)
    return best_score, endpoint, traces


def _profile_scores(hmm: plan7.HMM) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    background = np.asarray(plan7.Background(hmm.alphabet).residue_frequencies)
    with np.errstate(divide="ignore"):
        match_scores = np.log2(
            np.asarray(hmm.match_emissions) / background  # type: ignore[attr-defined]  # PyHMMER stubs omit its runtime matrices.
        )
        insert_scores = np.log2(
            np.asarray(hmm.insert_emissions) / background  # type: ignore[attr-defined]  # PyHMMER stubs omit its runtime matrices.
        )
        transition_scores = np.log2(
            np.asarray(hmm.transition_probabilities)  # type: ignore[attr-defined]  # PyHMMER stubs omit its runtime matrices.
        )
    return match_scores, insert_scores, transition_scores


def _profile_emission_tables(match_scores: np.ndarray) -> dict[int, np.ndarray]:
    tables: dict[int, np.ndarray] = {}
    for step, masks in _COMPATIBLE_MASK_TABLES.items():
        table = np.full((len(match_scores), len(masks) + 1), -np.inf)
        for amino_acid in range(len(_AMINO_ACIDS)):
            compatible = np.flatnonzero(masks[:, amino_acid])
            table[:, compatible] = np.maximum(
                table[:, compatible],
                match_scores[:, amino_acid, np.newaxis],
            )
        if step == 3:
            stop_indices = [_SEGMENT_INDICES[3][codon] for codon in _STOP_CODONS]
            table[:, stop_indices] += _STOP_PENALTY_BITS
            table[:, -1] = 0.0
        else:
            table[:, :-1] += _FRAME_PENALTY_BITS[abs(step - 3)]
        tables[step] = table
    return tables


def _delete_row(
    match_previous: np.ndarray,
    delete_previous: np.ndarray,
    transition_row: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    from_match = match_previous + transition_row[plan7.Transitions.MD]
    from_delete = delete_previous + transition_row[plan7.Transitions.DD]
    use_match = from_match >= from_delete
    scores = np.where(use_match, from_match, from_delete)
    codes = np.where(use_match, _MATCH, _DELETE).astype(np.uint8)
    return scores, codes


def _match_row(
    match_previous: np.ndarray,
    insert_previous: np.ndarray,
    delete_previous: np.ndarray,
    transition_row: np.ndarray,
    ordinary_codons: np.ndarray,
    segment_indices: dict[int, np.ndarray],
    emission_rows: dict[int, np.ndarray],
    steps: tuple[int, ...],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    nucleotide_count = len(match_previous) - 1
    scores = np.full(nucleotide_count + 1, -np.inf)
    codes = np.zeros(nucleotide_count + 1, dtype=np.uint8)
    ordinary_scores = np.full(nucleotide_count + 1, -np.inf)
    ordinary_codes = np.zeros(nucleotide_count + 1, dtype=np.uint8)

    for step in steps:
        if nucleotide_count < step:
            continue
        ends = np.arange(step, nucleotide_count + 1)
        starts = ends - step
        local_entry = np.full(len(ends), -np.inf)
        if step == 3:
            local_entry[ordinary_codons != _STOP] = 0.0
        sources = np.vstack(
            (
                match_previous[starts] + transition_row[plan7.Transitions.MM],
                insert_previous[starts] + transition_row[plan7.Transitions.IM],
                delete_previous[starts] + transition_row[plan7.Transitions.DM],
                local_entry,
            )
        )
        source_indices = np.argmax(sources, axis=0)
        source_codes = np.array((_MATCH, _INSERT, _DELETE, 0), dtype=np.uint8)
        candidates = np.max(sources, axis=0) + emission_rows[step][segment_indices[step]]
        candidate_codes = (step << 2) | source_codes[source_indices]
        if step == 3:
            ordinary_scores[ends] = candidates
            ordinary_codes[ends] = candidate_codes
        improved = candidates > scores[ends]
        scores[ends[improved]] = candidates[improved]
        codes[ends[improved]] = candidate_codes[improved]
    return scores, codes, ordinary_scores, ordinary_codes


def _insert_row(
    match_current: np.ndarray,
    emission_row: np.ndarray,
    transition_row: np.ndarray,
    ordinary_codons: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    scores = np.full_like(match_current, -np.inf)
    codes = np.zeros(len(match_current), dtype=np.uint8)
    emissions = np.zeros(len(ordinary_codons))
    sense = ordinary_codons >= 0
    emissions[sense] = emission_row[ordinary_codons[sense]]
    emissions[ordinary_codons == _STOP] = -np.inf

    open_transition = transition_row[plan7.Transitions.MI]
    extend_transition = transition_row[plan7.Transitions.II]
    for phase in range(3):
        ends = np.arange(3 + phase, len(match_current), 3)
        if len(ends) == 0:
            continue
        phase_emissions = emissions[ends - 3]
        _fill_insertion_phase(
            scores,
            codes,
            ends,
            match_current[ends - 3] + open_transition,
            phase_emissions,
            extend_transition,
        )
    return scores, codes


def _fill_insertion_phase(
    scores: np.ndarray,
    codes: np.ndarray,
    ends: np.ndarray,
    opening_scores: np.ndarray,
    emissions: np.ndarray,
    extend_transition: float,
) -> None:
    finite = np.isfinite(emissions)
    if not np.isfinite(extend_transition):
        scores[ends[finite]] = opening_scores[finite] + emissions[finite]
        codes[ends[finite]] = _MATCH
        return

    positions = np.arange(len(ends))
    segment_starts = finite & np.concatenate(([True], ~finite[:-1]))
    segment_start_indices = np.maximum.accumulate(np.where(segment_starts, positions, 0))
    segment_indices = np.cumsum(segment_starts) - 1

    cumulative_before = np.zeros(len(ends))
    if len(ends) > 1:
        extensions = np.where(finite[:-1], emissions[:-1] + extend_transition, 0.0)
        cumulative_before[1:] = np.cumsum(extensions)
    cumulative_before -= cumulative_before[segment_start_indices]
    transformed_openings = opening_scores - cumulative_before

    finite_openings = finite & np.isfinite(transformed_openings)
    if not np.any(finite_openings):
        codes[ends[finite]] = _MATCH
        return
    separation = float(np.ptp(transformed_openings[finite_openings])) + 1.0
    offsets = segment_indices * separation
    adjusted_openings = np.where(
        finite,
        transformed_openings + offsets,
        -np.inf,
    )
    running_adjusted = np.maximum.accumulate(adjusted_openings)
    last_opening = np.maximum.accumulate(np.where(finite_openings, positions, -1))
    has_opening = finite & (last_opening >= segment_start_indices)
    running_openings = running_adjusted - offsets
    scores[ends[has_opening]] = cumulative_before[has_opening] + emissions[has_opening] + running_openings[has_opening]

    prior_adjusted = np.concatenate(([-np.inf], running_adjusted[:-1]))
    prior_adjusted[segment_starts] = -np.inf
    opens_here = adjusted_openings >= prior_adjusted
    codes[ends[finite]] = np.where(
        opens_here[finite] | ~has_opening[finite],
        _MATCH,
        _INSERT,
    )


def _best_match_endpoint(
    scores: np.ndarray,
    ordinary_codons: np.ndarray,
) -> tuple[float, int]:
    ends = np.arange(len(scores))
    eligible = np.isfinite(scores)
    eligible[eligible] = ordinary_codons[ends[eligible] - 3] != _STOP
    eligible_scores = np.where(eligible, scores, -np.inf)
    end = int(np.argmax(eligible_scores))
    return float(eligible_scores[end]), end


def _trace_alignment(
    dna: str,
    score: float,
    endpoint: tuple[int, int, int],
    traces: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> CodonAlignment:
    match_trace, insert_trace, delete_trace = traces
    state, model_node, nucleotide_end = endpoint
    records: list[
        tuple[
            str,
            tuple[int, int],
            tuple[Literal["insertion", "deletion", "stop"], int, int] | None,
        ]
    ] = []
    model_positions: list[int] = []

    while state != 0:
        if state == _MATCH:
            code = int(match_trace[model_node, nucleotide_end])
            step = code >> 2
            previous_state = code & 3
            nucleotide_start = nucleotide_end - step
            model_position = model_node - 1
            amino_acid, event = _matched_residue(
                dna[nucleotide_start:nucleotide_end],
                model_position,
            )
            records.append(
                (
                    amino_acid,
                    (nucleotide_start, nucleotide_end),
                    event,
                )
            )
            model_positions.append(model_position)
            nucleotide_end = nucleotide_start
            model_node -= 1
            state = previous_state
            continue
        if state == _INSERT:
            previous_state = int(insert_trace[model_node, nucleotide_end])
            nucleotide_start = nucleotide_end - 3
            amino_acid = _translate_codon(dna[nucleotide_start:nucleotide_end])
            records.append((amino_acid, (nucleotide_start, nucleotide_end), None))
            nucleotide_end = nucleotide_start
            state = previous_state
            continue

        state = int(delete_trace[model_node, nucleotide_end])
        model_node -= 1

    records.reverse()
    spans = tuple(record[1] for record in records)
    events = tuple(
        FrameshiftEvent(
            kind=event[0],
            nt_start=span[0],
            nt_end=span[1],
            model_position=event[2],
            size=event[1],
            protein_position=protein_position,
        )
        for protein_position, (_, span, event) in enumerate(records)
        if event is not None
    )
    return CodonAlignment(
        nt_start=spans[0][0],
        nt_end=spans[-1][1],
        model_start=min(model_positions),
        model_end=max(model_positions) + 1,
        sequence="".join(record[0] for record in records),
        score=float(score),
        baseline_score=0.0,
        events=events,
        codon_spans=spans,
    )


def _matched_residue(
    segment: str,
    model_position: int,
) -> tuple[
    str,
    tuple[Literal["insertion", "deletion", "stop"], int, int] | None,
]:
    if len(segment) != 3:
        kind: Literal["insertion", "deletion"]
        kind = "insertion" if len(segment) > 3 else "deletion"
        return "X", (kind, len(segment) - 3, model_position)
    if segment in _STOP_CODONS:
        return "X", ("stop", 0, model_position)
    return _translate_codon(segment), None


def _translate_codon(codon: str) -> str:
    return _CODON_TO_AA.get(codon, "X")


def _ordinary_codon_indices(dna: str) -> np.ndarray:
    codon_count = max(len(dna) - 2, 0)
    indices = np.full(codon_count, -1, dtype=np.int8)
    amino_acid_indices = {amino_acid: index for index, amino_acid in enumerate(_AMINO_ACIDS)}
    for start in range(codon_count):
        codon = dna[start : start + 3]
        if codon in _STOP_CODONS:
            indices[start] = _STOP
            continue
        amino_acid = _CODON_TO_AA.get(codon)
        if amino_acid is not None:
            indices[start] = amino_acid_indices[amino_acid]
    return indices


def _segment_indices(dna: str, step: int) -> np.ndarray:
    segment_count = max(len(dna) - step + 1, 0)
    ambiguous_index = len(_SEGMENT_INDICES[step])
    return np.fromiter(
        (
            _SEGMENT_INDICES[step].get(
                dna[start : start + step],
                ambiguous_index,
            )
            for start in range(segment_count)
        ),
        dtype=np.int16,
        count=segment_count,
    )


def _build_compatible_aa_masks() -> dict[int, dict[str, int]]:
    amino_acid_indices = {amino_acid: index for index, amino_acid in enumerate(_AMINO_ACIDS)}
    masks_by_step: dict[int, dict[str, int]] = {step: {} for step in _FRAME_STEPS}
    codons = tuple("".join(bases) for bases in product(_DNA_BASES, repeat=3))

    for step in (1, 2):
        for bases in product(_DNA_BASES, repeat=step):
            segment = "".join(bases)
            compatible = (
                codon
                for codon in codons
                if any(
                    "".join(codon[index] for index in positions) == segment
                    for positions in combinations(range(3), step)
                )
            )
            masks_by_step[step][segment] = _sense_mask(compatible, amino_acid_indices)

    for codon, amino_acid in _CODON_TO_AA.items():
        masks_by_step[3][codon] = 1 << amino_acid_indices[amino_acid]
    for stop_codon in _STOP_CODONS:
        neighbors = (
            stop_codon[:position] + base + stop_codon[position + 1 :]
            for position in range(3)
            for base in _DNA_BASES
            if base != stop_codon[position]
        )
        masks_by_step[3][stop_codon] = _sense_mask(neighbors, amino_acid_indices)

    for step in (4, 5):
        for bases in product(_DNA_BASES, repeat=step):
            segment = "".join(bases)
            codon_subsequences = (
                "".join(segment[index] for index in positions) for positions in combinations(range(step), 3)
            )
            masks_by_step[step][segment] = _sense_mask(
                codon_subsequences,
                amino_acid_indices,
            )
    return masks_by_step


def _sense_mask(
    codons: Iterable[str],
    amino_acid_indices: dict[str, int],
) -> int:
    mask = 0
    for codon in codons:
        amino_acid = _CODON_TO_AA.get(codon)
        if amino_acid is not None:
            mask |= 1 << amino_acid_indices[amino_acid]
    return mask


_COMPATIBLE_AA_MASKS = _build_compatible_aa_masks()
_SEGMENT_INDICES = {
    step: {segment: index for index, segment in enumerate(sorted(masks))}
    for step, masks in _COMPATIBLE_AA_MASKS.items()
}
_AMINO_ACID_BITS = np.left_shift(
    np.uint32(1),
    np.arange(len(_AMINO_ACIDS), dtype=np.uint32),
)
_COMPATIBLE_MASK_TABLES = {
    step: np.asarray(
        [[bool(masks[segment] & bit) for bit in _AMINO_ACID_BITS] for segment in sorted(masks)],
        dtype=bool,
    )
    for step, masks in _COMPATIBLE_AA_MASKS.items()
}
