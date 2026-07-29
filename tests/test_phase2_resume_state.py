from __future__ import annotations

import copy
import json
from dataclasses import fields
from pathlib import Path

import numpy as np
import pytest

from virosync.features.compositional import WindowFeatures
from virosync.orchestration._flows.single_genome.phase_state import (
    PHASE2_STATE_ARTIFACT_TYPE,
    PHASE2_STATE_SCHEMA_VERSION,
    Phase2StateError,
    load_phase2_state,
    phase2_state_from_document,
    phase2_state_to_document,
    write_phase2_state,
)
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary


def _complete_boundary() -> RefinedBoundary:
    return RefinedBoundary(
        scaffold="scaffold/alpha",
        start=101,
        end=999,
        seed_id="seed_7_scaffold_alpha_90",
        original_start=90,
        original_end=1010,
        candidate_start=80,
        candidate_end=1020,
        host_trim_reason="host-taxonomy",
        host_trim_common_euk_taxonomy="Eukaryota;Viridiplantae",
        seed_sources=["hhg", "novelty", "compositional"],
        seed_confidence="high",
        seed_hhg_score=0.91,
        seed_novelty_score=0.82,
        seed_compositional_score=0.73,
        seed_has_mcp=True,
        predicted_family="NCLDV",
        region_classification_ncldv_markers=4,
        region_classification_vp_plv_markers=2,
        region_classification_mirus_markers=1,
        confidence=0.88,
        posterior_probability=0.93,
        core_viral_start=150,
        core_viral_end=850,
        flank_5_start=101,
        flank_5_end=150,
        flank_3_start=850,
        flank_3_end=999,
        state_sequence=[4, 5],
        state_posteriors=np.array(
            [
                [0.01, 0.02, 0.03, 0.04, 0.40, 0.50],
                [0.02, 0.03, 0.04, 0.01, 0.20, 0.70],
            ],
            dtype=np.float32,
        ),
        hallmark_genes=["MCP", "A32"],
        max_kfd=0.31,
        gc_deviation=0.14,
        cub_deviation=0.07,
        mean_novelty=0.81,
        window_features=[
            WindowFeatures(
                scaffold="scaffold/alpha",
                start=101,
                end=550,
                kfd=0.31,
                cub=0.07,
                gc_content=0.42,
                gc_deviation=0.14,
                porf_density=1.7,
                porf_count=3,
            ),
            WindowFeatures(
                scaffold="scaffold/alpha",
                start=550,
                end=999,
                kfd=0.29,
                cub=0.06,
                gc_content=0.43,
                gc_deviation=0.13,
                porf_density=2.1,
                porf_count=4,
            ),
        ],
    )


def _assert_boundaries_equal(
    expected: RefinedBoundary,
    observed: RefinedBoundary,
) -> None:
    for item in fields(RefinedBoundary):
        if item.name in {"state_posteriors", "window_features"}:
            continue
        assert getattr(observed, item.name) == getattr(expected, item.name)
    assert observed.window_features == expected.window_features
    if expected.state_posteriors is None:
        assert observed.state_posteriors is None
    else:
        assert observed.state_posteriors.dtype == expected.state_posteriors.dtype
        assert observed.state_posteriors.shape == expected.state_posteriors.shape
        np.testing.assert_array_equal(
            observed.state_posteriors,
            expected.state_posteriors,
        )


def test_phase2_state_file_round_trip_preserves_every_boundary_field(
    tmp_path: Path,
) -> None:
    original = _complete_boundary()
    state_path = tmp_path / "phase2" / "refined_state.json"

    write_phase2_state(state_path, [original])
    loaded = load_phase2_state(state_path)

    assert len(loaded) == 1
    _assert_boundaries_equal(original, loaded[0])
    payload = json.loads(state_path.read_text())
    assert payload["schema_version"] == PHASE2_STATE_SCHEMA_VERSION
    assert payload["artifact_type"] == PHASE2_STATE_ARTIFACT_TYPE
    assert payload["boundaries"][0]["state_posteriors"]["dtype"] == "float32"


def test_phase2_state_round_trip_preserves_none_and_empty_posterior_shape() -> None:
    without_posteriors = RefinedBoundary(scaffold="a", start=0, end=10)
    with_empty_posteriors = RefinedBoundary(
        scaffold="b",
        start=10,
        end=20,
        state_posteriors=np.empty((0, 6), dtype=np.float64),
    )

    loaded = phase2_state_from_document(
        phase2_state_to_document([without_posteriors, with_empty_posteriors])
    )

    _assert_boundaries_equal(without_posteriors, loaded[0])
    _assert_boundaries_equal(with_empty_posteriors, loaded[1])


def test_phase2_state_rejects_unknown_schema_and_field_drift() -> None:
    document = phase2_state_to_document([_complete_boundary()])

    unknown_schema = copy.deepcopy(document)
    # Derived, not literal, so a future schema bump cannot turn the "unknown"
    # version into the current one and silently stop testing the rejection.
    unknown_schema["schema_version"] = PHASE2_STATE_SCHEMA_VERSION + 1
    with pytest.raises(Phase2StateError, match="unsupported.*schema_version"):
        phase2_state_from_document(unknown_schema)

    missing_field = copy.deepcopy(document)
    del missing_field["boundaries"][0]["seed_id"]
    with pytest.raises(Phase2StateError, match="missing=.*seed_id"):
        phase2_state_from_document(missing_field)

    extra_field = copy.deepcopy(document)
    extra_field["boundaries"][0]["python_type"] = "arbitrary.Class"
    with pytest.raises(Phase2StateError, match="extra=.*python_type"):
        phase2_state_from_document(extra_field)


@pytest.mark.parametrize(
    ("mutate", "match"),
    [
        (
            lambda boundary: boundary.__setitem__("start", True),
            "start must be an integer",
        ),
        (
            lambda boundary: boundary.__setitem__("confidence", float("nan")),
            "confidence must be finite",
        ),
        (
            lambda boundary: boundary.__setitem__("seed_sources", "hhg"),
            "seed_sources must be a list",
        ),
        (
            lambda boundary: boundary["state_posteriors"].__setitem__(
                "dtype", "object"
            ),
            "dtype must be one of",
        ),
        (
            lambda boundary: boundary["state_posteriors"].__setitem__(
                "shape", [3, 6]
            ),
            "data length.*does not match shape",
        ),
        (
            lambda boundary: boundary["state_posteriors"].update(
                {"shape": [10**100, 0], "data": []}
            ),
            "cannot be reshaped",
        ),
        (
            lambda boundary: boundary["window_features"][0].__setitem__(
                "gc_content", float("inf")
            ),
            "gc_content must be finite",
        ),
    ],
)
def test_phase2_state_rejects_malformed_or_nonfinite_values(
    mutate,
    match: str,
) -> None:
    document = phase2_state_to_document([_complete_boundary()])
    mutate(document["boundaries"][0])

    with pytest.raises(Phase2StateError, match=match):
        phase2_state_from_document(document)


def test_phase2_state_file_rejects_duplicate_keys_and_nonstandard_numbers(
    tmp_path: Path,
) -> None:
    duplicate_path = tmp_path / "duplicate.json"
    duplicate_path.write_text(
        '{"artifact_type":"virosync.phase2.refined_boundaries",'
        '"schema_version":1,"schema_version":1,"boundaries":[]}'
    )
    with pytest.raises(Phase2StateError, match="duplicate JSON key"):
        load_phase2_state(duplicate_path)

    nonfinite_path = tmp_path / "nonfinite.json"
    nonfinite_path.write_text(
        '{"artifact_type":"virosync.phase2.refined_boundaries",'
        '"schema_version":1,"boundaries":NaN}'
    )
    with pytest.raises(Phase2StateError, match="non-finite JSON number"):
        load_phase2_state(nonfinite_path)


def test_phase2_state_never_interprets_executable_type_metadata(
    tmp_path: Path,
) -> None:
    sentinel = tmp_path / "must_not_exist"
    payload = phase2_state_to_document([_complete_boundary()])
    payload["boundaries"][0]["py/object/apply"] = {
        "callable": "pathlib.Path.touch",
        "args": [str(sentinel)],
    }
    state_path = tmp_path / "malicious.json"
    state_path.write_text(json.dumps(payload))

    with pytest.raises(Phase2StateError, match="extra=.*py/object/apply"):
        load_phase2_state(state_path)
    assert not sentinel.exists()


def test_phase2_state_failed_validation_preserves_existing_file(
    tmp_path: Path,
) -> None:
    state_path = tmp_path / "refined_state.json"
    state_path.write_text("previous valid state\n")
    boundary = _complete_boundary()
    boundary.posterior_probability = float("nan")

    with pytest.raises(Phase2StateError, match="posterior_probability must be finite"):
        write_phase2_state(state_path, [boundary])

    assert state_path.read_text() == "previous valid state\n"
