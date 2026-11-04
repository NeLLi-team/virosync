"""Coreless LOW acceptance uses final-interval, species-masked gene evidence."""

from types import SimpleNamespace
from unittest.mock import patch

import pytest

from virosync.pipeline.phase3.evidence_synthesizer import VerificationResult
from virosync.pipeline.phase3.output_generator import evaluate_v2_quality_gate


def gene(name, *, start=100, end=300, prefixes=None, bits=None, **changes):
    """Make an already searched, retained top-ten record."""
    prefixes = ["EUK-pEVE", "EUK"] if prefixes is None else prefixes
    record = dict(
        porf_id=name,
        scaffold="s",
        start=start,
        end=end,
        top10_prefixes=prefixes,
        top10_bitscores=[120.0, 100.0] if bits is None else bits,
        top10_targets=[prefix + "__reference" for prefix in prefixes],
        top10_pidents=[30.0] * len(prefixes),
        has_hit=bool(prefixes),
        is_flanking=False,
    )
    record.update(changes)
    return record


def candidate(*records, **changes):
    """Make an otherwise accepted LOW NCLDV candidate without a core."""
    values = dict(
        eve_id="eve",
        scaffold="s",
        start=0,
        end=15000,
        confidence_tier="LOW",
        region_classification="NCLDV",
        likely_family="NCLDV",
        hallmark_count=2,
        hallmark_genes=["shared_a", "shared_b"],
        gene_taxonomy_records=list(records),
    )
    values.update(changes)
    return VerificationResult(**values)


def test_short_coreless_low_is_removed_even_with_two_stronger_viral_genes():
    result = candidate(gene("a"), gene("b", start=500, end=700), end=14999)
    decision = evaluate_v2_quality_gate(result)
    assert not decision.kept
    assert decision.reason == "low_short_without_core"
    assert (result.start, result.end, result.confidence_tier) == (0, 14999, "LOW")


def test_exact_15000_bp_with_two_supporting_genes_keeps_original_gate_reason():
    decision = evaluate_v2_quality_gate(candidate(gene("a"), gene("b", start=500, end=700)))
    assert decision.kept and decision.promoted_low
    assert decision.reason == "ncldv_mirus_low_promoted"


@pytest.mark.parametrize("bits", [[100.0, 100.0], [99.0, 100.0]])
def test_ties_and_stronger_cellular_hits_do_not_support(bits):
    decision = evaluate_v2_quality_gate(candidate(gene("a"), gene("b", start=500, end=700, bits=bits)))
    assert not decision.kept
    assert decision.reason == "low_insufficient_competitive_viral_genes"


@pytest.mark.parametrize("prefix", ["EUK-pEVE", "CRESS", "GVMAG", "PLV", "VP", "MIRUS", "NCLDV"])
def test_non_phage_viral_classes_and_peve_count(prefix):
    result = candidate(gene("a"), gene("b", start=500, end=700, prefixes=[prefix], bits=[20.0]))
    assert evaluate_v2_quality_gate(result).kept


@pytest.mark.parametrize("prefix", ["NCLDV", "EUK-pEVE"])
def test_low_rule_requires_25_percent_identity_for_viral_support(prefix):
    result = candidate(gene("a"), gene("b", start=500, end=700, prefixes=[prefix], bits=[120.0], top10_pidents=[24.9]))
    assert not evaluate_v2_quality_gate(result).kept
    result.gene_taxonomy_records[1]["top10_pidents"] = [25.0]
    assert evaluate_v2_quality_gate(result).kept


def test_phage_only_does_not_supply_low_support():
    result = candidate(gene("a"), gene("b", start=500, end=700, prefixes=["PHAGE"], bits=[120.0]))
    assert not evaluate_v2_quality_gate(result).kept


def test_low_support_evidence_profile_is_auditable():
    result = candidate(gene("a"), gene("b", start=500, end=700))
    assert result.to_dict()["competitive_viral_support"] == dict(
        supporting_gene_ids=("a", "b"), assessed_gene_count=2, unassessed_gene_count=0
    )


@pytest.mark.parametrize("prefix", ["EUK", "BAC", "ARC"])
def test_all_cellular_classes_compete(prefix):
    result = candidate(gene("a"), gene("b", start=500, end=700, prefixes=["NCLDV", prefix], bits=[99.0, 100.0]))
    assert not evaluate_v2_quality_gate(result).kept


@pytest.mark.parametrize("prefix", ["UNKNOWN", "MITO", "PLASTID"])
def test_organelle_and_unknown_hits_are_neither_competitors_nor_support(prefix):
    result = candidate(gene("a"), gene("b", start=500, end=700, prefixes=[prefix, "MIRUS"], bits=[300.0, 20.0]))
    assert evaluate_v2_quality_gate(result).kept


@pytest.mark.parametrize(
    "record",
    [
        gene("outside", start=14900, end=15100),
        gene("wrong_scaffold", scaffold="other"),
        gene("unknown", prefixes=["UNKNOWN"], bits=[200.0]),
        gene("no_hit", prefixes=[], bits=[], viral_hit_within_depth=True),
    ],
)
def test_partial_other_scaffold_unknown_and_deep_only_do_not_add_support(record):
    decision = evaluate_v2_quality_gate(candidate(gene("a"), record))
    assert not decision.kept
    assert decision.reason == "low_insufficient_competitive_viral_genes"


def test_final_coordinates_determine_interior_instead_of_stale_flanking_flag():
    result = candidate(gene("a"), gene("b", start=500, end=700, is_flanking=True))
    assert evaluate_v2_quality_gate(result).kept


def test_multiple_models_on_one_gene_are_one_supporting_gene():
    decision = evaluate_v2_quality_gate(candidate(gene("a|aa1-40"), gene("a|aa50-90")))
    assert not decision.kept


def test_linked_fragments_are_one_biological_gene():
    result = candidate(gene("a"), gene("b", start=500, end=700), scoring_locus_ids={"a": "locus", "b": "locus"})
    assert not evaluate_v2_quality_gate(result).kept


def test_physical_aliases_do_not_supply_two_genes():
    assert not evaluate_v2_quality_gate(candidate(gene("a"), gene("alias"))).kept


def test_strong_cellular_fragment_prevents_cherry_picking_viral_fragment():
    result = candidate(
        gene("a"),
        gene("b", start=500, end=700, bits=[100.0, 150.0]),
        gene("c", start=800, end=900),
        scoring_locus_ids={"a": "locus", "b": "locus"},
    )
    assert not evaluate_v2_quality_gate(result).kept


@pytest.mark.parametrize("reverse", [False, True])
def test_base_gene_fragments_pool_competitors_at_different_coordinates(reverse):
    records = [
        gene("a|aa1-40"),
        gene("a|aa50-90", start=500, end=700, bits=[100.0, 150.0]),
        gene("b", start=800, end=900),
    ]
    result = candidate(*(reversed(records) if reverse else records))
    assert evaluate_v2_quality_gate(result).reason == "low_insufficient_competitive_viral_genes"
    assert result.to_dict()["competitive_viral_support"]["assessed_gene_count"] == 2


def test_fragment_physical_alias_and_scoring_locus_links_pool_transitively():
    result = candidate(
        gene("a|aa1-40"),
        gene("a|aa50-90", start=500, end=700),
        gene("alias", start=500, end=700),
        gene("linked", start=800, end=900, bits=[100.0, 150.0]),
        gene("other", start=1000, end=1100),
        scoring_locus_ids={"alias": "locus", "linked": "locus"},
    )
    assert evaluate_v2_quality_gate(result).reason == "low_insufficient_competitive_viral_genes"
    assert result.to_dict()["competitive_viral_support"]["assessed_gene_count"] == 2


@pytest.mark.parametrize("bits", [[120.0], [float("nan"), 100.0], [120.0, "bad"], [120.0, None]])
def test_malformed_hit_arrays_are_unavailable_not_zero_competitors(bits):
    decision = evaluate_v2_quality_gate(candidate(gene("a"), gene("b", start=500, end=700, bits=bits)))
    assert not decision.kept
    assert decision.reason == "low_gene_evidence_unavailable"


def test_missing_taxonomy_differs_from_observed_no_hit():
    assert evaluate_v2_quality_gate(candidate()).reason == "low_gene_evidence_unavailable"
    observed = candidate(gene("a", prefixes=[], bits=[]))
    assert evaluate_v2_quality_gate(observed).reason == "low_insufficient_competitive_viral_genes"


def test_missing_scaffold_is_unavailable_instead_of_silently_outside():
    result = candidate(gene("a"), gene("b", start=500, end=700, scaffold=None))
    assert evaluate_v2_quality_gate(result).reason == "low_gene_evidence_unavailable"
    assert result.to_dict()["competitive_viral_support"]["unassessed_gene_count"] == 1


def test_only_scored_qualifying_core_exempts_low():
    result = candidate(end=6000)
    result.hallmark_assessment = SimpleNamespace(
        summary=SimpleNamespace(qualifying_evidence_ids=("core",), panel="NCLDV")
    )
    result.hallmark_evidence = [
        SimpleNamespace(
            evidence_id="core",
            model_roles=[dict(panel="NCLDV", role="core_marker", excluded_from_scoring=False, marker="A32")],
        )
    ]
    assert evaluate_v2_quality_gate(result).kept
    result.hallmark_evidence[0].model_roles[0]["excluded_from_scoring"] = True
    assert evaluate_v2_quality_gate(result).reason == "low_short_without_core"


@pytest.mark.parametrize("tier", ["MEDIUM", "HIGH"])
def test_extra_low_rule_does_not_touch_higher_tiers(tier):
    assert evaluate_v2_quality_gate(candidate(end=6000, confidence_tier=tier)).kept


def test_spatial_rescue_does_not_bypass_coreless_low_requirement():
    with patch("virosync.pipeline.phase3.output_generator.has_spatial_rescue", return_value=True):
        assert evaluate_v2_quality_gate(candidate(end=6000)).reason == "low_short_without_core"


def test_cress_route_does_not_bypass_coreless_low_requirement():
    result = candidate(end=6000, region_classification="CRESS", likely_family="CRESS", marker_family_hits=["CRESS"])
    assert evaluate_v2_quality_gate(result).reason == "low_short_without_core"
