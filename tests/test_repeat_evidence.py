"""Exercise bounded direct-repeat evidence independently of boundary decisions."""

from __future__ import annotations

import json
import random
from dataclasses import replace

import pytest

from virosync.pipeline.phase2 import repeat_evidence
from virosync.pipeline.phase2.repeat_evidence import assess_repeat_evidence, evidence_from_dict, evidence_to_dict
from virosync.pipeline.phase2.terminal_repeats import (
    TerminalRepeat,
    find_target_site_duplication,
    find_target_site_duplication_at,
)

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _dna(length: int, seed: int = 29) -> str:
    rng = random.Random(seed)
    return "".join(rng.choices("ACGT", k=length))


def _planted(*, arm_length: int = 60, substitutions: tuple[int, ...] = (), inverted: bool = False) -> str:
    arm = _dna(arm_length)
    partner = list(arm)
    for index in substitutions:
        partner[index] = {"A": "C", "C": "G", "G": "T", "T": "A"}[partner[index]]
    paired_arm = "".join(partner)
    if inverted:
        paired_arm = paired_arm.translate(_COMPLEMENT)[::-1]
    return "C" * 100 + arm + "A" * 200 + paired_arm + "G" * 100


def _assess(sequence: str | None, **overrides: object) -> repeat_evidence.RepeatEvidence:
    parameters = {
        "scaffold": "host",
        "input_id": "seed-1",
        "start": 100,
        "end": 420,
        "parent_start": 100,
        "parent_end": 420,
        "coverage_intervals": ((0, len(sequence)),) if sequence is not None else (),
        "extension_bp": 90,
    }
    parameters.update(overrides)
    return assess_repeat_evidence(sequence, **parameters)


@pytest.mark.parametrize(("substitutions", "identity"), [((), 1.0), ((10, 20, 30, 40, 50, 55), 0.9)])
def test_direct_pair_retains_both_arms_without_boundary_authority(
    substitutions: tuple[int, ...], identity: float
) -> None:
    evidence = _assess(_planted(substitutions=substitutions))

    assert len(evidence.candidates) == 1
    pair = evidence.candidates[0]
    assert (pair.left_start, pair.left_end, pair.right_start, pair.right_end) == (100, 160, 360, 420)
    assert pair.identity == pytest.approx(identity)
    assert pair.interpretation == "unresolved"
    assert (evidence.assessed_start, evidence.assessed_end) == (100, 420)
    assert evidence.tsd.status == "not_assessed_no_anchor"


def test_inverted_pair_is_not_reported_as_direct() -> None:
    evidence = _assess(_planted(inverted=True))

    assert evidence.candidates == ()
    assert evidence.tsd.status == "not_assessed_no_anchor"


@pytest.mark.parametrize(("length", "count"), [(49, 0), (50, 1), (615, 1)])
def test_direct_repeat_length_threshold_and_long_arm(length: int, count: int) -> None:
    evidence = _assess(_planted(arm_length=length), end=300 + 2 * length, extension_bp=length + 30)

    assert len(evidence.candidates) == count
    assert not evidence.candidates or evidence.candidates[0].alignment_length == length


def test_outward_pair_and_input_cut_jitter_preserve_geometry() -> None:
    sequence = _planted()
    first = _assess(sequence, start=125, end=395)
    jittered = _assess(sequence, start=130, end=390)

    assert first.candidates[0].left_start == jittered.candidates[0].left_start == 100
    assert first.candidates[0].right_end == jittered.candidates[0].right_end == 420
    assert first.assessed_start == 125
    assert jittered.assessed_start == 130


def test_reverse_complement_preserves_uncapped_pair_geometry() -> None:
    sequence = _planted()
    reverse = sequence.translate(_COMPLEMENT)[::-1]

    forward = _assess(sequence)
    backward = _assess(reverse)

    assert forward.candidates[0].alignment_length == backward.candidates[0].alignment_length == 60
    assert backward.candidates[0].left_start == len(sequence) - forward.candidates[0].right_end
    assert backward.candidates[0].right_end == len(sequence) - forward.candidates[0].left_start


def test_coordinate_zero_and_contig_clipping_retain_positive() -> None:
    sequence = _planted()[100:]

    evidence = _assess(sequence, start=0, end=320, parent_start=0, parent_end=320)

    assert evidence.candidates[0].left_start == 0
    assert evidence.left.status == "incomplete"
    assert "contig_clipped" in evidence.left.reasons


def test_ambiguous_window_and_positive_candidate_coexist() -> None:
    sequence = "N" * 50 + _planted()[50:]

    evidence = _assess(sequence)

    assert len(evidence.candidates) == 1
    assert evidence.left.ambiguous_bp == 40
    assert evidence.right.ambiguous_bp == 0
    assert "ambiguous_bases" in evidence.left.reasons


def test_ambiguous_arm_is_not_bridged_as_mismatch() -> None:
    sequence = _planted()
    sequence = sequence[:130] + "N" + sequence[131:]

    evidence = _assess(sequence)

    assert evidence.candidates == ()
    assert "ambiguous_bases" in evidence.left.reasons


def test_low_complexity_host_sequence_yields_no_candidate() -> None:
    evidence = _assess("ACGT" * 130)

    assert evidence.candidates == ()
    assert "low_complexity" in evidence.left.reasons
    assert "low_complexity" in evidence.right.reasons


def test_repeat_negative_and_missing_input_are_distinct() -> None:
    negative = _assess(_dna(520))
    missing = _assess(None)

    assert negative.candidates == missing.candidates == ()
    assert negative.left.status == "completed"
    assert missing.left.status == "not_assessed"
    assert missing.left.reasons == ("missing_sequence",)


def test_taxonomy_gaps_are_not_filled_by_enclosing_span() -> None:
    evidence = _assess(_planted(), coverage_intervals=((0, 165), (350, 520)))

    assert len(evidence.candidates) == 1
    assert evidence.left.window_end == 165
    assert evidence.right.window_start == 350
    assert "taxonomy_clipped" in evidence.left.reasons
    assert "taxonomy_clipped" in evidence.right.reasons


def test_missing_endpoint_coverage_is_not_assessed() -> None:
    evidence = _assess(_planted(), coverage_intervals=((0, 80), (200, 520)))

    assert evidence.candidates == ()
    assert evidence.left.status == "not_assessed"
    assert evidence.left.reasons == ("missing_coverage",)
    assert evidence.right.status == "incomplete"
    assert "partner_unavailable" in evidence.right.reasons


def test_coverage_traversal_order_does_not_change_ids() -> None:
    sequence = _planted()
    forward = _assess(sequence, coverage_intervals=((0, 260), (260, 520)))
    reversed_coverage = _assess(sequence, coverage_intervals=((260, 520), (0, 260)))

    assert forward == reversed_coverage
    assert _assess(sequence, input_id="seed-2").evidence_id != forward.evidence_id


def test_overlapping_windows_do_not_make_self_match_candidates() -> None:
    evidence = _assess(_dna(520), start=200, end=300, extension_bp=180)

    assert evidence.candidates == ()
    assert evidence.filtered_overlapping == 0
    assert evidence.left.window_end == evidence.right.window_start == 250
    assert evidence.left.reasons == evidence.right.reasons == ("midpoint_partition",)
    assert evidence.left.status == evidence.right.status == "completed"


def test_shared_arm_alternatives_are_retained() -> None:
    arm = _dna(60)
    sequence = "C" * 100 + arm + "A" * 200 + arm + "G" * 40 + arm + "T" * 100

    evidence = _assess(sequence, extension_bp=170)

    assert {(pair.left_start, pair.right_start) for pair in evidence.candidates} == {(100, 360), (100, 460)}
    assert len({pair.candidate_id for pair in evidence.candidates}) == 2


@pytest.mark.parametrize(
    ("limit", "reason"), [("MAX_SEED_PAIRINGS", "seed_pairing_limit"), ("MAX_DIAGONALS", "diagonal_limit")]
)
def test_exhausted_search_limit_is_explicit(monkeypatch: pytest.MonkeyPatch, limit: str, reason: str) -> None:
    monkeypatch.setattr(repeat_evidence, limit, 0)

    evidence = _assess(_planted())

    assert evidence.candidates == ()
    assert reason in evidence.left.reasons
    assert reason in evidence.right.reasons


def test_candidate_limit_keeps_positive_and_marks_incomplete(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(repeat_evidence, "MAX_CANDIDATES", 1)
    arm = _dna(60)
    sequence = "C" * 100 + arm + "A" * 200 + arm + "G" * 40 + arm + "T" * 100

    evidence = _assess(sequence, extension_bp=170)

    assert len(evidence.candidates) == 1
    assert "candidate_limit" in evidence.left.reasons


def test_window_limit_is_recorded_without_enormous_search() -> None:
    evidence = _assess(_planted(), extension_bp=10_001)

    assert "window_limit" in evidence.left.reasons
    assert evidence.max_window_radius_bp == 10_000


def _legacy_tsd(sequence: str, start: int, end: int) -> str:
    """Independent frozen reference for the pre-change short-flank calculation."""
    sequence = sequence.upper()
    available = min(start, len(sequence) - end, 9)
    for length in range(available, 2, -1):
        left, right = sequence[start - length : start], sequence[end : end + length]
        if left == right and set(left) <= set("ACGT") and len(set(left)) > 1:
            return left
    return ""


@pytest.mark.parametrize("duplication", ["ACG", "ACGTA", "ACGTAGTCA", "AAAA", "ACN", "acgt"])
def test_coordinate_core_preserves_legacy_anchored_tsd(duplication: str) -> None:
    sequence = "C" * 30 + duplication + _dna(60) + "G" * 100 + _dna(60, 42) + duplication + "T" * 30
    start, end = 30 + len(duplication), 250 + len(duplication)
    tir = TerminalRepeat(start, start + 60, end - 60, end, 1.0, 60, False)

    assert find_target_site_duplication(sequence, tir) == _legacy_tsd(sequence, start, end)
    assert find_target_site_duplication_at(sequence, start, end) == _legacy_tsd(sequence, start, end)


def test_shifted_tsd_diagnostic_uses_same_translation_offsets() -> None:
    sequence = "C" * 40 + "ACG" + _dna(60) + "G" * 100 + _dna(60, 42) + "ACG" + "T" * 40
    tir = TerminalRepeat(43, 103, 203, 263, 1.0, 60, False)

    evidence = _assess(sequence, start=43, end=263, tir=tir)
    expected = sum(bool(_legacy_tsd(sequence, 43 + offset, 263 + offset)) for offset in (*range(-20, 0), *range(1, 21)))

    assert evidence.tsd.sequence == "ACG"
    assert evidence.tsd.status == "assessed_match"
    assert evidence.tsd.null_assessed == 40
    assert evidence.tsd.null_matches == expected
    assert evidence.tsd.null_geometry == "translate-both-termini;offsets=-20..-1,+1..+20"


def test_tsd_clipped_long_flanks_can_keep_short_match() -> None:
    sequence = "ACG" + _dna(60) + "G" * 100 + _dna(60, 42) + "ACG"
    tir = TerminalRepeat(3, 63, 163, 223, 1.0, 60, False)

    evidence = _assess(sequence, start=3, end=223, tir=tir)

    assert evidence.tsd.sequence == "ACG"
    assert evidence.tsd.status == "incomplete"
    assert evidence.tsd.left_reasons == ("contig_clipped",)
    assert evidence.tsd.right_reasons == ("contig_clipped",)


def test_evidence_json_roundtrip_preserves_candidates_and_limits() -> None:
    evidence = _assess(_planted())

    assert evidence_from_dict(json.loads(json.dumps(evidence_to_dict(evidence)))) == evidence


@pytest.mark.parametrize(
    ("field", "value"),
    [("evidence_id", "wrong"), ("min_identity", 0.85), ("assessed_start", True), ("method", "future-method")],
)
def test_evidence_codec_rejects_corrupt_contract(field: str, value: object) -> None:
    payload = evidence_to_dict(_assess(_planted()))
    payload[field] = value

    with pytest.raises(ValueError):
        evidence_from_dict(payload)


def test_changed_endpoints_keep_parent_pair_without_display_authority() -> None:
    evidence = _assess(_planted())
    inherited = replace(
        evidence,
        left=replace(evidence.left, status="not_assessed", reasons=("endpoint_changed",)),
        display_candidate_id="",
        tsd=repeat_evidence.TsdAssessment("not_assessed_endpoint_changed"),
    )

    decoded = evidence_from_dict(evidence_to_dict(inherited))

    assert decoded == inherited
    assert decoded.candidates == evidence.candidates
    assert decoded.display_candidate_id == ""


@pytest.mark.parametrize("radius", [5_000, 10_000])
def test_realistic_random_windows_complete_without_spending_cap_on_ineligible_diagonals(radius: int) -> None:
    # Diagnostic seed 0 exceeded the old 1024 cap at both window sizes. Match
    # that panel's two consecutive random windows, separated by an unsearched gap.
    width = 2 * radius
    paired_windows = _dna(2 * width, 0)
    sequence = paired_windows[:width] + "N" * width + paired_windows[width:]

    evidence = _assess(sequence, start=radius, end=5 * radius, extension_bp=radius)

    assert evidence.candidates == ()
    assert evidence.left.status == evidence.right.status == "completed"
    assert evidence.left.reasons == evidence.right.reasons == ()


def test_midpoint_partition_excludes_short_eve_upstream_repeat_confounder() -> None:
    arm = _dna(60)
    sequence = _dna(100, 71) + arm + _dna(20, 72) + arm + _dna(500, 73)

    evidence = _assess(sequence, start=300, end=350, extension_bp=300)

    assert evidence.candidates == ()
    assert evidence.left.window_end == evidence.right.window_start == 325
    assert "midpoint_partition" in evidence.left.reasons
    assert "midpoint_partition" in evidence.right.reasons


def test_midpoint_partition_keeps_repeat_outside_both_input_cuts() -> None:
    evidence = _assess(_planted(), start=200, end=300, extension_bp=180)

    assert len(evidence.candidates) == 1
    assert (evidence.candidates[0].left_start, evidence.candidates[0].right_end) == (100, 420)
    assert evidence.left.window_end == evidence.right.window_start == 250
    assert evidence.left.status == evidence.right.status == "completed"
    assert "midpoint_partition" in evidence.left.reasons


def test_pairing_budget_prefers_rare_seeds_over_alphabetical_order(monkeypatch: pytest.MonkeyPatch) -> None:
    index = {"AAAACCC": (0, 1, 2, 3), "TTTCGAC": (10,), "TTTGACG": (30,)}
    monkeypatch.setattr(repeat_evidence, "_kmer_index", lambda sequence: index)
    monkeypatch.setattr(repeat_evidence, "MAX_SEED_PAIRINGS", 2)

    diagonals, limits = repeat_evidence._seed_diagonals("A" * 100, "T" * 100)

    assert diagonals == {0: (10, 30)}
    assert limits == {"seed_pairing_limit"}


def test_unchainable_diagonals_do_not_exhaust_diagonal_budget(monkeypatch: pytest.MonkeyPatch) -> None:
    left_sequence, right_sequence = "A" * 1000, "C" * 1000
    indexes = {
        left_sequence: {"ACGCTGA": (0,), "GATCGTC": (600,)},
        right_sequence: {"ACGCTGA": (10,), "GATCGTC": (610,)},
    }
    monkeypatch.setattr(repeat_evidence, "_kmer_index", indexes.__getitem__)
    monkeypatch.setattr(repeat_evidence, "MAX_DIAGONALS", 0)

    diagonals, limits = repeat_evidence._seed_diagonals(left_sequence, right_sequence)

    assert diagonals == {}
    assert limits == set()


def test_recorded_parameters_are_the_reused_scoring_parameters() -> None:
    evidence = _assess(_planted())

    assert evidence.min_arm_bp == repeat_evidence.MIN_TIR_ARM_BP == 50
    assert evidence.min_identity == repeat_evidence.MIN_TIR_IDENTITY == 0.9
    assert evidence.seed_bp == repeat_evidence._KMER_BP == 7
    assert evidence.max_seed_chain_bp == 493
    assert evidence.min_seed_entropy_bits == 1.2


@pytest.mark.parametrize(
    ("section", "field", "value"),
    [
        ("tsd", "null_geometry", "independent-random-sites"),
        ("tsd", "sequence", "ACG"),
        ("tsd", "anchor_start", 100),
        ("tsd", "null_assessed", 1),
        ("tsd", "left_reasons", ["invented"]),
        ("left", "reasons", ["invented"]),
        ("left", "status", "incomplete"),
        ("left", "requested_start", 11),
    ],
)
def test_codec_rejects_closed_assessment_contract_corruption(section: str, field: str, value: object) -> None:
    payload = evidence_to_dict(_assess(_planted()))
    payload[section][field] = value

    with pytest.raises(ValueError):
        evidence_from_dict(payload)


def test_codec_rejects_partial_tsd_without_recorded_limitation() -> None:
    evidence = _assess(_planted())
    invalid = replace(evidence, tsd=repeat_evidence.TsdAssessment("incomplete", "retained_tir", 100, 420))

    with pytest.raises(ValueError, match="requires a limitation"):
        evidence_from_dict(evidence_to_dict(invalid))


def test_codec_rejects_excess_candidates_before_parameter_validation() -> None:
    evidence = replace(_assess(_planted()), max_candidates=0)

    with pytest.raises(ValueError, match="retention limit"):
        evidence_from_dict(evidence_to_dict(evidence))


def test_codec_preserves_alternative_order_without_changing_display_choice() -> None:
    arm = _dna(60)
    sequence = "C" * 100 + arm + "A" * 200 + arm + "G" * 40 + arm + "T" * 100
    evidence = _assess(sequence, extension_bp=170)
    reordered = replace(evidence, candidates=tuple(reversed(evidence.candidates)))

    assert evidence_from_dict(evidence_to_dict(reordered)) == reordered
    assert reordered.display_candidate_id == evidence.display_candidate_id


def test_endpoint_invalidation_has_one_model_home() -> None:
    evidence = _assess(_planted())

    changed = repeat_evidence.for_changed_endpoints(evidence, start=80, end=420)

    assert changed.assessed_start == 100
    assert changed.left.status == "not_assessed"
    assert changed.right == evidence.right
    assert changed.tsd.status == "not_assessed_endpoint_changed"
    assert changed.display_candidate_id == ""
    assert evidence_from_dict(evidence_to_dict(changed)) == changed


def test_codec_rejects_completed_search_with_a_shortened_method_window() -> None:
    evidence = _assess(_dna(520), start=200, end=300, extension_bp=180)
    shortened = replace(evidence, left=replace(evidence.left, window_start=21))

    with pytest.raises(ValueError, match="whole method window"):
        evidence_from_dict(evidence_to_dict(shortened))


@pytest.mark.parametrize(("start", "end", "reasons"), [(200, 300, ()), (100, 420, ("midpoint_partition",))])
def test_codec_rejects_midpoint_reason_inconsistent_with_geometry(
    start: int, end: int, reasons: tuple[str, ...]
) -> None:
    evidence = _assess(_dna(520), start=start, end=end, extension_bp=90)
    inconsistent = replace(evidence, left=replace(evidence.left, reasons=reasons))

    with pytest.raises(ValueError, match="midpoint partition reason"):
        evidence_from_dict(evidence_to_dict(inconsistent))


def test_codec_requires_reason_for_ambiguous_bases() -> None:
    evidence = _assess(_dna(520))
    unexplained = replace(evidence, left=replace(evidence.left, ambiguous_bp=1))

    with pytest.raises(ValueError, match="lacks its assessment reason"):
        evidence_from_dict(evidence_to_dict(unexplained))


def test_codec_retains_missing_input_without_fabricating_partition_geometry() -> None:
    evidence = _assess(None, start=200, end=300, extension_bp=180)

    assert evidence_from_dict(evidence_to_dict(evidence)) == evidence
    assert evidence.left.window_start is None
    assert evidence.left.reasons == ("missing_sequence",)


def test_codec_retains_changed_endpoint_with_original_partition_geometry() -> None:
    evidence = _assess(_dna(520), start=200, end=300, extension_bp=180)
    changed = repeat_evidence.for_changed_endpoints(evidence, start=190, end=300)

    assert evidence_from_dict(evidence_to_dict(changed)) == changed
    assert changed.left.reasons == ("endpoint_changed", "midpoint_partition")
