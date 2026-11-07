from __future__ import annotations

from dataclasses import replace

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


def test_reported_penton_survives_stronger_non_capsid_domain_assignment():
    capsid = _hit("EVE_CAPSID_NOVEL_0002", -5)
    noncapsid = _hit("model-a", 200)
    retained, records = arbitrate_hits(
        [noncapsid, capsid],
        {"protein-1": {"HostDomain"}},
        {"model-a": _annotation("HostDomain"), "EVE_CAPSID_NOVEL_0002": _annotation("CapsidDomain")},
    )
    assert retained == [noncapsid, capsid]
    assert records[0].final_model == "model-a"


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


@pytest.mark.parametrize(
    ("observed", "expected_models", "outcome"),
    [
        ({"HostDomain"}, {"GVOGm0003", "GVOGm0054"}, "reassigned"),
        ({"UnrelatedDomain"}, {"GVOGm0003"}, "contradicted"),
    ],
)
def test_overridden_capsid_keeps_only_reported_floor_and_selected_non_capsid(
    tmp_path, observed, expected_models, outcome
):
    from virosync.ablation import AblationID
    from virosync.pipeline.hallmark_evidence import collect_hallmark_evidence, find_spatial_witnesses
    from virosync.pipeline.phase1.marker_roles import MarkerRole, decide_marker_hit_role
    from virosync.pipeline.phase1.marker_validation import filter_validated_markers, load_validated_marker_hits

    capsid = _hit("GVOGm0003", 200)
    noncapsid = _hit("GVOGm0054", 100)
    retained, records = arbitrate_hits(
        [capsid, noncapsid],
        {"protein-1": observed},
        {"GVOGm0003": _annotation("CapsidDomain"), "GVOGm0054": _annotation("HostDomain")},
    )
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">protein-1 # 101 # 400 # 1 # ID=1_1\n" + "M" * 100 + "\n")
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("protein-1\tNCLDV__virus|p1\t1e-30\t100\t90\t100\n")
    markers = filter_validated_markers(retained, diamond, proteome, output_dir=tmp_path / "validation")
    assert records[0].outcome == outcome
    assert {marker.hmm_target for marker in markers} == expected_models
    capsid_marker = next(marker for marker in markers if marker.hmm_target == "GVOGm0003")
    assert capsid_marker.validation_status == "validated"
    assert decide_marker_hit_role(capsid_marker, ablation_id=AblationID.A0).role is MarkerRole.REPORTED_CAPSID
    assert capsid_marker.pfam_override == records[0].outcome
    loaded = load_validated_marker_hits(tmp_path / "validation/validated_marker_hits.tsv")
    assert next(marker for marker in loaded if marker.hmm_target == "GVOGm0003").pfam_override == records[0].outcome
    evidence = collect_hallmark_evidence(markers, retained)
    partner = replace(
        evidence[0],
        query_porf="partner|aa1-100",
        gene="partner",
        model="GVOGm0760",
        start=500,
        end=800,
        pfam_override="",
    )
    witnesses = find_spatial_witnesses([*evidence, partner])
    assert {(w.rule, w.tier) for w in witnesses} == {("reported_capsid", "MEDIUM")}
