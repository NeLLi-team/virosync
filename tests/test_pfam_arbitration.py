from __future__ import annotations

import pytest

from virosync.pipeline.phase1.hhg_seeding import HMMHit
from virosync.pipeline.phase1.pfam_arbitration import (
    ModelPfamAnnotation,
    arbitrate_hits,
)


def _hit(model: str, score: float) -> HMMHit:
    return HMMHit(
        query_name="protein-1",
        target_name=model,
        score=score,
        evalue=1e-20,
        domain_score=score,
        query_start=1,
        query_end=100,
    )


def _annotation(
    *domains: str,
    source_scope: str = "",
) -> ModelPfamAnnotation:
    return ModelPfamAnnotation(frozenset(domains), source_scope)


@pytest.mark.parametrize(
    ("observed", "annotations", "expected_model", "expected_outcome"),
    [
        (
            {"Domain_A"},
            {"model-a": _annotation("Domain_A"), "model-b": _annotation("Domain_B")},
            "model-a",
            "confirmed",
        ),
        (
            {"Domain_B"},
            {"model-a": _annotation("Domain_A"), "model-b": _annotation("Domain_B")},
            "model-b",
            "reassigned",
        ),
        (
            set(),
            {"model-a": _annotation("Domain_A"), "model-b": _annotation("Domain_B")},
            "model-a",
            "unresolved_no_domain",
        ),
        (
            {"Shared"},
            {"model-a": _annotation("Shared"), "model-b": _annotation("Shared")},
            "model-a",
            "unresolved_shared_domain",
        ),
        (
            {"Wrong"},
            {"model-a": _annotation("Domain_A"), "model-b": _annotation("Domain_B")},
            None,
            "contradicted",
        ),
    ],
)
def test_arbitration_outcomes(
    observed: set[str],
    annotations: dict[str, ModelPfamAnnotation],
    expected_model: str | None,
    expected_outcome: str,
) -> None:
    retained, records = arbitrate_hits(
        [_hit("model-a", 100.0), _hit("model-b", 90.0)],
        {"protein-1": observed},
        annotations,
    )

    assert [hit.target_name for hit in retained] == ([expected_model] if expected_model is not None else [])
    assert records[0].final_model == expected_model
    assert records[0].outcome == expected_outcome


def test_cress_rep_requires_catalytic_domain_in_worked_pox_d5_case() -> None:
    annotations = {
        "VS000806": _annotation("AAA", "RNA_helicase", source_scope="CRESS_REP"),
        "VS000369": _annotation("Pox_A32"),
    }

    retained, records = arbitrate_hits(
        [_hit("VS000806", 100.0), _hit("VS000369", 90.0)],
        {"protein-1": {"Pox_D5", "DUF5906", "VapE-like_dom"}},
        annotations,
    )

    assert retained == []
    assert records[0].outcome == "contradicted"


def test_cress_rep_accepts_gemini_al1_but_not_its_generic_signature() -> None:
    annotations = {
        "VS000806": _annotation("AAA", "RNA_helicase", source_scope="CRESS_REP"),
        "other": _annotation("Other"),
    }
    hits = [_hit("VS000806", 100.0), _hit("other", 90.0)]

    catalytic, _ = arbitrate_hits(
        hits,
        {"protein-1": {"Gemini_AL1"}},
        annotations,
    )
    generic, generic_records = arbitrate_hits(
        hits,
        {"protein-1": {"AAA"}},
        annotations,
    )

    assert [hit.target_name for hit in catalytic] == ["VS000806"]
    assert generic == []
    assert generic_records[0].outcome == "contradicted"
