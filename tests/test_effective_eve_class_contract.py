from __future__ import annotations

import ast
import inspect
import textwrap
from types import SimpleNamespace

from virosync.orchestration._flows.single_genome import phase1, phase2
from virosync.orchestration._flows.single_genome.manifest import (
    _empty_prediction_summary,
)
from virosync.output_contract import (
    EFFECTIVE_EVE_CLASSES,
    EFFECTIVE_EVE_CLASS_COUNT_KEYS,
    OUTPUT_SCHEMA_VERSION,
    effective_eve_class_count_total,
    normalize_effective_eve_class,
    resolve_effective_eve_class,
)
from virosync.pipeline.phase3.output_generator import evaluate_v2_quality_gate


def test_effective_class_partition_and_output_schema_are_versioned() -> None:
    assert EFFECTIVE_EVE_CLASSES == (
        "NCLDV",
        "VP",
        "PLV",
        "MIRUS",
        "MIXED",
        "PPV",
        "UNKNOWN",
    )
    assert tuple(EFFECTIVE_EVE_CLASS_COUNT_KEYS) == EFFECTIVE_EVE_CLASSES
    assert OUTPUT_SCHEMA_VERSION == 3


def test_persisted_effective_class_normalization_is_exhaustive() -> None:
    assert normalize_effective_eve_class(" ppv ") == "PPV"
    assert normalize_effective_eve_class("mixed") == "MIXED"
    assert normalize_effective_eve_class("") == "UNKNOWN"
    assert normalize_effective_eve_class("future-lineage") == "UNKNOWN"
    assert normalize_effective_eve_class(None) == "UNKNOWN"


def test_tier_aware_resolver_preserves_low_gate_precedence() -> None:
    labels = {
        "region_classification": "PPV",
        "classification": "NCLDV",
        "likely_family": "NCLDV",
    }
    assert resolve_effective_eve_class(confidence_tier="HIGH", **labels) == "PPV"
    assert resolve_effective_eve_class(confidence_tier="LOW", **labels) == "NCLDV"

    decision = evaluate_v2_quality_gate(
        SimpleNamespace(
            confidence_tier="LOW",
            start=0,
            end=6001,
            hallmark_count=2,
            hallmark_genes=["marker_a", "marker_b"],
            has_mcp=False,
            **labels,
        )
    )
    assert decision.kept is True
    assert decision.effective_class == "NCLDV"


def test_empty_prediction_summary_has_full_exclusive_surface() -> None:
    summary = _empty_prediction_summary()

    assert all(summary[key] == 0 for key in EFFECTIVE_EVE_CLASS_COUNT_KEYS.values())
    assert effective_eve_class_count_total(summary) == summary["accepted"] == 0
    assert summary["high_tier"] == summary["candidate_high_tier"] == 0
    assert summary["medium_tier"] == summary["candidate_medium_tier"] == 0
    assert summary["low_tier"] == summary["candidate_low_tier"] == 0


def _successful_return_dicts(function) -> list[ast.Dict]:
    tree = ast.parse(textwrap.dedent(inspect.getsource(function)))
    successful = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Dict):
            continue
        for key, value in zip(node.keys, node.values):
            if (
                isinstance(key, ast.Constant)
                and key.value == "success"
                and isinstance(value, ast.Constant)
                and value.value is True
            ):
                successful.append(node)
                break
    return successful


def _expands_empty_summary(node: ast.Dict) -> bool:
    return any(
        key is None
        and isinstance(value, ast.Call)
        and isinstance(value.func, ast.Name)
        and value.func.id == "_empty_prediction_summary"
        for key, value in zip(node.keys, node.values)
    )


def _has_literal_key(node: ast.Dict, expected: str) -> bool:
    return any(
        isinstance(key, ast.Constant) and key.value == expected for key in node.keys
    )


def test_all_successful_phase1_phase2_zero_returns_use_full_summary() -> None:
    phase1_returns = _successful_return_dicts(phase1._run_phase1_subflow)
    phase2_returns = _successful_return_dicts(phase2._run_phase2_subflow)

    assert len(phase1_returns) == 2
    assert len(phase2_returns) == 2
    assert all(_expands_empty_summary(node) for node in phase1_returns + phase2_returns)
    assert all(_has_literal_key(node, "elapsed_sec") for node in phase1_returns + phase2_returns)
