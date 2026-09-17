from __future__ import annotations

import pytest
from pyhmmer import easel, plan7

from virosync.pipeline.phase1.native_frameshift import (
    CodonAlignment,
    align_frameshift_profile,
)

PROTEIN = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"
CODONS = {
    "A": "GCT",
    "C": "TGC",
    "D": "GAT",
    "E": "GAA",
    "F": "TTT",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "K": "AAA",
    "L": "CTG",
    "M": "ATG",
    "N": "AAT",
    "P": "CCT",
    "Q": "CAA",
    "R": "CGT",
    "S": "TCT",
    "T": "ACT",
    "V": "GTT",
    "W": "TGG",
    "Y": "TAT",
}
DNA = "".join(CODONS[amino_acid] for amino_acid in PROTEIN)


@pytest.fixture(scope="module")
def profile() -> plan7.HMM:
    alphabet = easel.Alphabet.amino()
    sequence = easel.TextSequence(name=b"diverse", sequence=PROTEIN).digitize(alphabet)
    hmm, _, _ = plan7.Builder(alphabet).build(sequence, plan7.Background(alphabet))
    return hmm


def test_intact_sequence_has_complete_ordinary_alignment(profile: plan7.HMM) -> None:
    alignment = _required_alignment(profile, DNA)

    assert alignment.sequence == PROTEIN
    assert (alignment.nt_start, alignment.nt_end) == (0, len(DNA))
    assert (alignment.model_start, alignment.model_end) == (0, len(PROTEIN))
    assert alignment.codon_spans[0] == (0, 3)
    assert alignment.codon_spans[-1] == (len(DNA) - 3, len(DNA))
    assert len(alignment.codon_spans) == len(PROTEIN)
    assert alignment.events == ()
    assert alignment.score == pytest.approx(alignment.baseline_score)


def test_sense_substitutions_remain_ordinary_codons(profile: plan7.HMM) -> None:
    substituted = _replace_codon(_replace_codon(DNA, 8, "GAA"), 31, "GTT")
    alignment = _required_alignment(profile, substituted)

    assert alignment.sequence == "ACDEFGHIELMNPQRSTVWYACDEFGHIKLMVPQRSTVWY"
    assert alignment.events == ()
    assert alignment.score == pytest.approx(alignment.baseline_score)


@pytest.mark.parametrize(
    ("mutated_dna", "expected_sequence"),
    [
        pytest.param(
            DNA[:60] + "GCT" + DNA[60:],
            "ACDEFGHIKLMNPQRSTVWYAACDEFGHIKLMNPQRSTVWY",
            id="three-base-insertion",
        ),
        pytest.param(
            DNA[:60] + "GCTTGC" + DNA[60:],
            "ACDEFGHIKLMNPQRSTVWYACACDEFGHIKLMNPQRSTVWY",
            id="six-base-insertion",
        ),
        pytest.param(
            DNA[:60] + DNA[63:],
            "ACDEFGHIKLMNPQRSTVWYCDEFGHIKLMNPQRSTVWY",
            id="three-base-deletion",
        ),
        pytest.param(
            DNA[:60] + DNA[66:],
            "ACDEFGHIKLMNPQRSTVWYDEFGHIKLMNPQRSTVWY",
            id="six-base-deletion",
        ),
    ],
)
def test_in_frame_indels_use_profile_gaps_without_frame_events(
    profile: plan7.HMM,
    mutated_dna: str,
    expected_sequence: str,
) -> None:
    alignment = _required_alignment(profile, mutated_dna)

    assert alignment.sequence == expected_sequence
    assert alignment.events == ()
    assert alignment.score == pytest.approx(alignment.baseline_score)


@pytest.mark.parametrize(
    ("mutated_dna", "kind", "size", "event_span", "protein_position", "expected_sequence"),
    [
        pytest.param(
            DNA[:60] + "A" + DNA[60:],
            "insertion",
            1,
            (57, 61),
            19,
            "ACDEFGHIKLMNPQRSTVWXACDEFGHIKLMNPQRSTVWY",
            id="one-base-insertion-between-codons",
        ),
        pytest.param(
            DNA[:61] + "A" + DNA[61:],
            "insertion",
            1,
            (60, 64),
            20,
            "ACDEFGHIKLMNPQRSTVWYXCDEFGHIKLMNPQRSTVWY",
            id="one-base-insertion-within-codon",
        ),
        pytest.param(
            DNA[:60] + "AC" + DNA[60:],
            "insertion",
            2,
            (54, 59),
            18,
            "ACDEFGHIKLMNPQRSTVXYACDEFGHIKLMNPQRSTVWY",
            id="two-base-insertion-between-codons",
        ),
        pytest.param(
            DNA[:61] + "AC" + DNA[61:],
            "insertion",
            2,
            (60, 65),
            20,
            "ACDEFGHIKLMNPQRSTVWYXCDEFGHIKLMNPQRSTVWY",
            id="two-base-insertion-within-codon",
        ),
        pytest.param(
            DNA[:60] + DNA[61:],
            "deletion",
            -1,
            (60, 62),
            20,
            "ACDEFGHIKLMNPQRSTVWYXCDEFGHIKLMNPQRSTVWY",
            id="one-base-deletion-between-codons",
        ),
        pytest.param(
            DNA[:61] + DNA[62:],
            "deletion",
            -1,
            (60, 62),
            20,
            "ACDEFGHIKLMNPQRSTVWYXCDEFGHIKLMNPQRSTVWY",
            id="one-base-deletion-within-codon",
        ),
        pytest.param(
            DNA[:60] + DNA[62:],
            "deletion",
            -2,
            (60, 61),
            20,
            "ACDEFGHIKLMNPQRSTVWYXCDEFGHIKLMNPQRSTVWY",
            id="two-base-deletion-between-codons",
        ),
        pytest.param(
            DNA[:61] + DNA[63:],
            "deletion",
            -2,
            (60, 61),
            20,
            "ACDEFGHIKLMNPQRSTVWYXCDEFGHIKLMNPQRSTVWY",
            id="two-base-deletion-within-codon",
        ),
    ],
)
def test_frameshift_indels_export_one_x_with_exact_span(
    profile: plan7.HMM,
    mutated_dna: str,
    kind: str,
    size: int,
    event_span: tuple[int, int],
    protein_position: int,
    expected_sequence: str,
) -> None:
    alignment = _required_alignment(profile, mutated_dna)
    (event,) = alignment.events

    assert alignment.sequence == expected_sequence
    assert event.kind == kind
    assert event.size == size
    assert event.protein_position == protein_position
    assert (event.nt_start, event.nt_end) == event_span
    assert alignment.codon_spans[protein_position] == event_span
    assert event.nt_end - event.nt_start == size + 3
    assert alignment.score >= alignment.baseline_score


def test_internal_stop_exports_x_and_stop_event(profile: plan7.HMM) -> None:
    stopped = _replace_codon(DNA, 20, "TAA")
    alignment = _required_alignment(profile, stopped)
    (event,) = alignment.events

    assert alignment.sequence == "ACDEFGHIKLMNPQRSTVWYXCDEFGHIKLMNPQRSTVWY"
    assert (event.kind, event.size, event.model_position, event.protein_position) == (
        "stop",
        0,
        20,
        20,
    )
    assert (event.nt_start, event.nt_end) == (60, 63)
    assert alignment.score == pytest.approx(alignment.baseline_score)


def test_two_compensating_shifts_are_both_reported(profile: plan7.HMM) -> None:
    first_shift = DNA[:36] + "A" + DNA[36:]
    compensated = first_shift[:85] + first_shift[86:]
    alignment = _required_alignment(profile, compensated)

    assert alignment.sequence == "ACDEFGHIKLMXPQRSTVWYACDEFGHIXLMNPQRSTVWY"
    assert tuple(event.kind for event in alignment.events) == ("insertion", "deletion")
    assert tuple(event.size for event in alignment.events) == (1, -1)
    assert tuple(event.protein_position for event in alignment.events) == (11, 28)
    assert alignment.score > alignment.baseline_score


def test_ambiguous_codons_are_neutral_x_without_events(profile: plan7.HMM) -> None:
    ambiguous = DNA[:54] + "NNNNNN" + DNA[60:]
    alignment = _required_alignment(profile, ambiguous)

    assert alignment.sequence == "ACDEFGHIKLMNPQRSTVXXACDEFGHIKLMNPQRSTVWY"
    assert alignment.codon_spans[18:20] == ((54, 57), (57, 60))
    assert alignment.events == ()
    assert alignment.score == pytest.approx(alignment.baseline_score)


def test_terminal_partial_codons_are_clipped_without_events(profile: plan7.HMM) -> None:
    alignment = _required_alignment(profile, DNA[1:-1])

    assert (alignment.nt_start, alignment.nt_end) == (2, len(DNA) - 4)
    assert (alignment.model_start, alignment.model_end) == (1, len(PROTEIN) - 1)
    assert alignment.sequence == PROTEIN[1:-1]
    assert alignment.events == ()


def test_empty_and_subcodon_inputs_have_no_alignment(profile: plan7.HMM) -> None:
    assert align_frameshift_profile(profile, "") is None
    assert align_frameshift_profile(profile, DNA[:2]) is None


def test_baseline_is_conditioned_on_selected_hit_span(profile: plan7.HMM) -> None:
    damaged = DNA[:60] + "A" + DNA[60:]
    isolated = _required_alignment(profile, damaged)
    partial_paralog = DNA[:66]
    padded = _required_alignment(profile, partial_paralog + "NNNNNN" + damaged)

    assert padded.nt_start == len(partial_paralog) + 6
    assert padded.sequence == isolated.sequence
    assert padded.baseline_score == pytest.approx(isolated.baseline_score)
    assert padded.score == pytest.approx(isolated.score)
    assert padded.score >= padded.baseline_score


def _required_alignment(hmm: plan7.HMM, dna: str) -> CodonAlignment:
    alignment = align_frameshift_profile(hmm, dna)
    assert alignment is not None
    return alignment


def _replace_codon(dna: str, position: int, codon: str) -> str:
    start = position * 3
    return dna[:start] + codon + dna[start + 3 :]
