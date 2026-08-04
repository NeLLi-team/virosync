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
    normalize_effective_eve_class_counts,
    resolve_effective_eve_class,
)
from virosync.pipeline.phase3.output_generator import evaluate_v2_quality_gate
from virosync.pipeline.phase3.evidence_synthesizer import (
    VerificationResult,
    consensus_taxonomy_class,
    infer_likely_family,
    marker_taxonomy_category,
)


def test_effective_class_partition_and_output_schema_are_versioned() -> None:
    assert EFFECTIVE_EVE_CLASSES == (
        "NCLDV",
        "MIRUS",
        "PPV",
        "CRESS",
        "PHAGE",
        "VIRAL_UNKNOWN",
        "UNKNOWN",
    )
    assert tuple(EFFECTIVE_EVE_CLASS_COUNT_KEYS) == EFFECTIVE_EVE_CLASSES
    assert "MIXED" not in EFFECTIVE_EVE_CLASSES
    assert OUTPUT_SCHEMA_VERSION == 5


def test_persisted_effective_class_normalization_is_exhaustive() -> None:
    assert normalize_effective_eve_class(" ppv ") == "PPV"
    assert normalize_effective_eve_class("vp") == "PPV"
    assert normalize_effective_eve_class("plv") == "PPV"
    assert normalize_effective_eve_class("cress") == "CRESS"
    assert normalize_effective_eve_class("phage") == "PHAGE"
    assert normalize_effective_eve_class("viral_unknown") == "VIRAL_UNKNOWN"
    # MIXED is a retired class, kept as a read alias for schema-4 result trees.
    assert normalize_effective_eve_class("mixed") == "VIRAL_UNKNOWN"
    assert normalize_effective_eve_class("") == "UNKNOWN"
    assert normalize_effective_eve_class("future-lineage") == "UNKNOWN"
    assert normalize_effective_eve_class(None) == "UNKNOWN"


def test_gate_resolver_keeps_mixed_out_of_the_published_partition() -> None:
    # The gate scores MIXED under its own accepting branch, so the resolver must
    # keep answering MIXED even though the published partition retired it.
    labels = {
        "region_classification": "MIXED",
        "classification": "MIXED",
        "likely_family": "MIXED",
    }
    for tier in ("HIGH", "MEDIUM", "LOW"):
        assert resolve_effective_eve_class(confidence_tier=tier, **labels) == "MIXED"
    assert normalize_effective_eve_class("MIXED") == "VIRAL_UNKNOWN"


def test_legacy_class_counts_fold_vp_plv_and_mixed() -> None:
    counts = normalize_effective_eve_class_counts(
        {
            "NCLDV": 1,
            "VP": 2,
            "PLV": 3,
            "MIRUS": 4,
            "MIXED": 5,
            "PPV": 6,
            "UNKNOWN": 7,
        }
    )

    assert counts == {
        "NCLDV": 1,
        "MIRUS": 4,
        "PPV": 11,
        "CRESS": 0,
        "PHAGE": 0,
        "VIRAL_UNKNOWN": 5,
        "UNKNOWN": 7,
    }


def test_schema4_class_counts_fold_mixed_into_viral_unknown() -> None:
    counts = normalize_effective_eve_class_counts(
        {
            "NCLDV": 1,
            "MIRUS": 2,
            "PPV": 3,
            "CRESS": 4,
            "MIXED": 5,
            "UNKNOWN": 6,
        }
    )

    assert counts == {
        "NCLDV": 1,
        "MIRUS": 2,
        "PPV": 3,
        "CRESS": 4,
        "PHAGE": 0,
        "VIRAL_UNKNOWN": 5,
        "UNKNOWN": 6,
    }


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


def test_cress_without_identity_marker_support_is_rejected() -> None:
    decision = evaluate_v2_quality_gate(
        SimpleNamespace(
            confidence_tier="HIGH",
            start=0,
            end=3000,
            hallmark_count=2,
            hallmark_genes=["marker_a", "marker_b"],
            marker_family_hits=[],
            has_mcp=False,
            region_classification="UNKNOWN",
            classification="",
            likely_family="CRESS",
        )
    )

    assert decision.kept is False
    assert decision.effective_class == "CRESS"
    assert decision.reason == "cress_identity_required"


def test_cress_is_reachable_from_identity_qualified_gene_taxonomy() -> None:
    result = VerificationResult(
        eve_id="EVE_CRESS",
        scaffold="ctg",
        start=0,
        end=3000,
        region_classification="UNKNOWN",
        marker_dominant_family="UNKNOWN",
        gene_taxonomy_dominant_family="CRESS",
    )

    assert infer_likely_family(result) == "CRESS"


def test_empty_prediction_summary_has_full_exclusive_surface() -> None:
    summary = _empty_prediction_summary()

    assert all(summary[key] == 0 for key in EFFECTIVE_EVE_CLASS_COUNT_KEYS.values())
    assert effective_eve_class_count_total(summary) == summary["accepted"] == 0
    assert summary["high_tier"] == summary["candidate_high_tier"] == 0
    assert summary["medium_tier"] == summary["candidate_medium_tier"] == 0
    assert summary["low_tier"] == summary["candidate_low_tier"] == 0


def _marker(
    prefixes: str,
    pidents: str,
    *,
    gene: str = "gvogm0100",
    targets: str = "",
    validation_status: str = "validated",
    porf_id: str = "",
) -> dict:
    """One boundary hallmark_hits entry, shaped as build_boundary_evidence packs it."""
    return {
        "hallmark_gene": gene,
        "porf_id": porf_id,
        "top10_prefixes": prefixes,
        "top10_pidents": pidents,
        "top10_targets": targets,
        "validation_status": validation_status,
    }


def _gene(porf_id, prefixes, pidents, targets="", is_flanking=False):
    """One gene_taxonomy_records entry, as build_gene_taxonomy_record packs it."""
    return {
        "porf_id": porf_id,
        "top10_prefixes": prefixes,
        "top10_pidents": pidents,
        "top10_targets": targets,
        "is_flanking": is_flanking,
    }


def test_marker_votes_with_its_own_top10_taxonomy() -> None:
    assert marker_taxonomy_category(_marker("NCLDV,NCLDV", "80,70")) == "NCLDV"
    assert marker_taxonomy_category(_marker("NCLDV,MIRUS", "80,70")) == "VIRAL_UNKNOWN"
    # No qualified viral hit is no vote, not a vote for UNKNOWN.
    assert marker_taxonomy_category(_marker("EUK,BAC", "90,90")) == ""
    assert marker_taxonomy_category(_marker("NCLDV", "24.9")) == ""

    assert consensus_taxonomy_class([_marker("NCLDV,NCLDV", "80,70")]) == "NCLDV"
    assert consensus_taxonomy_class([_marker("NCLDV,MIRUS", "80,70")]) == "VIRAL_UNKNOWN"


def test_markers_without_a_vote_leave_the_denominator() -> None:
    # One NCLDV vote beside one voteless marker is unanimous. Counting the
    # voteless marker would make it 1/2 and demote the EVE to VIRAL_UNKNOWN.
    assert (
        consensus_taxonomy_class(
            [_marker("NCLDV", "80"), _marker("EUK", "90")]
        )
        == "NCLDV"
    )
    # No vote at all: validated markers still say "viral", nothing says which.
    assert consensus_taxonomy_class([_marker("EUK", "90")]) == "VIRAL_UNKNOWN"
    assert (
        consensus_taxonomy_class(
            [_marker("EUK", "90", validation_status="unvalidated")]
        )
        == "UNKNOWN"
    )
    assert consensus_taxonomy_class([]) == "UNKNOWN"


def test_lineage_needs_strictly_more_than_half_the_vote_weight() -> None:
    ncldv = _marker("NCLDV", "80")
    mirus = _marker("MIRUS", "80")

    # Two markers at 2 each. An exact tie is not a call.
    assert consensus_taxonomy_class([ncldv, mirus]) == "VIRAL_UNKNOWN"
    # 4 of 6 weight.
    assert consensus_taxonomy_class([ncldv, ncldv, mirus]) == "NCLDV"
    assert consensus_taxonomy_class([ncldv] * 7 + [mirus] * 3) == "NCLDV"


def test_a_marker_gene_confirmed_by_the_all_gene_search_weighs_three() -> None:
    ncldv = _marker("NCLDV", "80", porf_id="p1")
    mirus = _marker("MIRUS", "80", porf_id="p2")
    # Without the second search: 2 vs 2, a tie.
    assert consensus_taxonomy_class([ncldv, mirus]) == "VIRAL_UNKNOWN"
    # The all-gene search agreeing with the NCLDV marker takes it to 3 of 5.
    confirming = [_gene("p1", "NCLDV", "80")]
    assert (
        consensus_taxonomy_class(
            [ncldv, mirus], gene_taxonomy_records=confirming
        )
        == "NCLDV"
    )
    # Disagreeing instead: NCLDV 2, MIRUS 2, PPV 1. Nothing clears half.
    conflicting = [_gene("p1", "PPV", "80")]
    assert (
        consensus_taxonomy_class(
            [ncldv, mirus], gene_taxonomy_records=conflicting
        )
        == "VIRAL_UNKNOWN"
    )


def test_genes_without_a_marker_vote_once() -> None:
    # No marker at all: the all-gene search alone decides.
    genes = [_gene("g1", "PPV", "80"), _gene("g2", "PPV", "80"), _gene("g3", "NCLDV", "80")]
    assert consensus_taxonomy_class([], gene_taxonomy_records=genes) == "PPV"
    # Flanking genes are outside the element and never vote.
    flanking = [_gene("g9", "NCLDV", "80", is_flanking=True)] * 5
    assert (
        consensus_taxonomy_class([], gene_taxonomy_records=genes + flanking) == "PPV"
    )


def test_a_single_mcp_vote_overrides_every_other_gene() -> None:
    ncldv = _marker("NCLDV", "80")
    mirus_mcp = _marker("MIRUS", "80", gene="mirus_mcp_2")

    # However many other markers disagree, the capsid decides.
    assert consensus_taxonomy_class([ncldv, ncldv, mirus_mcp]) == "MIRUS"
    assert consensus_taxonomy_class([ncldv] * 20 + [mirus_mcp]) == "MIRUS"

    # An MCP whose own top-10 spans two lineages overrides with VIRAL_UNKNOWN.
    mixed_mcp = _marker("NCLDV,MIRUS", "80,80", gene="plv_mcp_1")
    assert consensus_taxonomy_class([ncldv, ncldv, mixed_mcp]) == "VIRAL_UNKNOWN"

    # An MCP with no qualified viral hit casts no vote and overrides nothing.
    silent_mcp = _marker("EUK", "90", gene="plv_mcp_1")
    assert consensus_taxonomy_class([ncldv, ncldv, silent_mcp]) == "NCLDV"


def test_disagreeing_mcp_markers_are_settled_by_the_weighted_vote() -> None:
    """The 5x weight is the MCP tiebreak, not the primary mechanism."""
    ncldv_mcp = _marker("NCLDV", "80", gene="plv_mcp_1", porf_id="m1")
    mirus_mcp = _marker("MIRUS", "80", gene="mirus_mcp_2", porf_id="m2")

    # Two MCPs alone split 5 and 5: a tie is not a call.
    assert consensus_taxonomy_class([ncldv_mcp, mirus_mcp]) == "VIRAL_UNKNOWN"

    # A third gene backing one of them breaks it: NCLDV 6, MIRUS 5.
    assert (
        consensus_taxonomy_class(
            [ncldv_mcp, mirus_mcp],
            gene_taxonomy_records=[_gene("g1", "NCLDV", "80")],
        )
        == "NCLDV"
    )
    # Two MCPs agreeing is not a tie at all; it is a single override.
    assert (
        consensus_taxonomy_class([ncldv_mcp, _marker("NCLDV", "80", gene="vp_mcp_3")])
        == "NCLDV"
    )


def test_gvmag_and_preplasmiviricota_phage_namespaces_fold_to_their_lineage() -> None:
    assert (
        consensus_taxonomy_class([_marker("GVMAG", "80", targets="GVMAG__x|p1")])
        == "NCLDV"
    )
    # A marker mixing the GVMAG namespace with NCLDV proper is still NCLDV.
    assert (
        consensus_taxonomy_class(
            [_marker("GVMAG,NCLDV", "80,80", targets="GVMAG__x|p1,NCLDV__y|p1")]
        )
        == "NCLDV"
    )

    legacy_phage = _marker(
        "PHAGE",
        "80",
        targets="PHAGE__VARDNA__Sputnik|p1",
    )
    lookup = {"PHAGE__Sputnik": "Viruses;Preplasmiviricota;Maveriviricetes"}
    assert consensus_taxonomy_class([legacy_phage], taxonomy_lookup=lookup) == "PPV"
    assert (
        consensus_taxonomy_class(
            [legacy_phage],
            taxonomy_lookup={"PHAGE__Sputnik": "Viruses;Uroviricota;Caudoviricetes"},
        )
        == "PHAGE"
    )
    # Without a lineage lookup the namespace stays PHAGE rather than guessing.
    assert consensus_taxonomy_class([legacy_phage]) == "PHAGE"


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


def test_a_marker_gene_is_never_counted_twice() -> None:
    """The same gene appears in both evidence lists; it is one voter."""
    from virosync.pipeline.phase3.evidence_synthesizer import taxonomy_class_votes

    marker = _marker("NCLDV", "80", porf_id="p1")
    # Confirmed by its own all-gene record: weight 3, not 2 plus a separate 1.
    assert taxonomy_class_votes(
        [marker], gene_taxonomy_records=[_gene("p1", "NCLDV", "80")]
    ) == {"NCLDV": 3}
    # The HMM domain suffix must not break the join.
    assert taxonomy_class_votes(
        [_marker("NCLDV", "80", porf_id="p1|aa1")],
        gene_taxonomy_records=[_gene("p1", "NCLDV", "80")],
    ) == {"NCLDV": 3}
    # A marker with no viral hit of its own votes once, through its gene record.
    assert taxonomy_class_votes(
        [_marker("EUK", "90", porf_id="p1")],
        gene_taxonomy_records=[_gene("p1", "PPV", "80")],
    ) == {"PPV": 1}


def test_an_mcp_weighs_exactly_five_even_when_its_own_gene_search_disagrees() -> None:
    """Two MCPs must meet at 5 each, or one MCP's second search breaks the tie."""
    from virosync.pipeline.phase3.evidence_synthesizer import taxonomy_class_votes

    ncldv_mcp = _marker("NCLDV", "80", gene="plv_mcp_1", porf_id="p1")
    ppv_mcp = _marker("PPV", "80", gene="vp_mcp_3", porf_id="p2")
    # p1's all-gene record contradicts its capsid call and must add nothing.
    records = [_gene("p1", "PPV", "80"), _gene("p2", "PPV", "80")]

    assert taxonomy_class_votes(
        [ncldv_mcp, ppv_mcp], gene_taxonomy_records=records
    ) == {"NCLDV": 5, "PPV": 5}
    assert (
        consensus_taxonomy_class([ncldv_mcp, ppv_mcp], gene_taxonomy_records=records)
        == "VIRAL_UNKNOWN"
    )


def test_an_mcp_profile_is_not_hidden_by_a_higher_scoring_profile() -> None:
    """One protein, two profiles. The capsid must not lose to a bigger score.

    Whether a gene carries an MCP is a property of every profile that hit it.
    Ranking by score alone would cost the gene its weight of 5, its override,
    and its eligibility to donate a class during ANI propagation.
    """
    from virosync.pipeline.phase3.evidence_synthesizer import taxonomy_class_votes

    def hit(gene, score):
        return dict(_marker("PPV", "80", gene=gene, porf_id="p1"), hmm_score=score)

    # Same protein: the MCP wins despite the lower HMM score.
    assert taxonomy_class_votes([hit("vp_mcp_3", 80.0), hit("gvogm0100", 90.0)]) == {
        "PPV": 5
    }
    # Order must not matter.
    assert taxonomy_class_votes([hit("gvogm0100", 90.0), hit("vp_mcp_3", 80.0)]) == {
        "PPV": 5
    }
    # Distinct proteins still vote separately: 5 for the MCP, 2 for the marker.
    other = dict(_marker("PPV", "80", gene="gvogm0100", porf_id="p2"), hmm_score=90.0)
    assert taxonomy_class_votes([hit("vp_mcp_3", 80.0), other]) == {"PPV": 7}
