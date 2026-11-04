"""Accepted pairs retain native class, tier, core ownership and repeat barriers."""

from pathlib import Path

import pytest
from test_joint_seed_rescue import child, rejected
from test_joint_seed_rescue import pair_inputs as pair_inputs
from test_native_hallmark_evidence import hit

from virosync.ablation import AblationID
from virosync.pipeline.host_signatures import HostSignatureModel
from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy, pORF
from virosync.pipeline.phase2.boundary_refiner import boundary_candidate_id
from virosync.pipeline.phase3 import joint_seed_rescue as joint


@pytest.fixture
def accepted_inputs(pair_inputs: dict) -> dict:
    left, right = pair_inputs["results"]
    left.hallmark_evidence = [hit("polb", "GVOGm0054", 200, 400)]
    right.hallmark_evidence = [hit("mcp", "GVOGm0003", 1600, 1800)]
    pair_inputs["accepted_ids"] = frozenset(result.eve_id for result in pair_inputs["results"])
    left.confidence_tier = right.confidence_tier = "MEDIUM"
    return pair_inputs


@pytest.mark.parametrize("tier", ["LOW", "MEDIUM", "HIGH"])
def test_accepted_complementary_pair_can_cross_nonoverlapping_envelopes(accepted_inputs: dict, tier: str) -> None:
    """Independent neighbors without separating evidence can join by design."""
    accepted_inputs["results"][0].confidence_tier = tier
    accepted_inputs["results"][1].confidence_tier = tier
    accepted_inputs["query"].seed_gene_mappings["left"].flank_end_bp = 1000

    assert joint.eligible_joint_pairs(**accepted_inputs) == [tuple(accepted_inputs["boundaries"])]


@pytest.mark.parametrize("gap,allowed", [(5000, True), (5001, False)])
def test_accepted_gap_limit(accepted_inputs: dict, gap: int, allowed: bool) -> None:
    right = child("right", 1000 + gap, 2000 + gap)
    scored = rejected(right)
    scored.confidence_tier = "MEDIUM"
    scored.hallmark_evidence = [hit("mcp", "GVOGm0003", right.start + 100, right.start + 300)]
    accepted_inputs["boundaries"][1] = right
    accepted_inputs["results"][1] = scored
    accepted_inputs["accepted_ids"] = frozenset(result.eve_id for result in accepted_inputs["results"])

    assert bool(joint.eligible_joint_pairs(**accepted_inputs)) is allowed


@pytest.mark.parametrize(
    "family,models,allowed",
    [
        ("NCLDV", ("GVOGm0760", "GVOGm0760"), False),
        ("NCLDV", ("GVOGm0003", "GVOGm0003"), True),
        ("MIRUS", ("GVOGm0023", "VS000307"), False),
        ("MIRUS", ("GVOGm0023", "GVOGm0022"), True),
    ],
)
def test_qualified_catalogue_families_distinguish_rnap_subunits(
    accepted_inputs: dict, family: str, models: tuple[str, str], allowed: bool
) -> None:
    for index, result in enumerate(accepted_inputs["results"]):
        result.likely_family = result.region_classification = family
        result.hallmark_evidence = [hit(f"g{index}", models[index], result.start + 100, result.start + 300)]
    assert bool(joint.eligible_joint_pairs(**accepted_inputs)) is allowed


def test_no_qualifying_core_and_mixed_acceptance_withhold(accepted_inputs: dict) -> None:
    accepted_inputs["results"][0].hallmark_evidence = []
    accepted_inputs["results"][1].hallmark_evidence = []
    assert joint.eligible_joint_pairs(**accepted_inputs) == []
    accepted_inputs["accepted_ids"] = frozenset({accepted_inputs["results"][0].eve_id})
    assert joint.eligible_joint_pairs(**accepted_inputs) == []


def test_hull_gene_outside_gap_must_have_query_and_taxonomy(accepted_inputs: dict) -> None:
    gene = pORF("interior", "s", 300, 400)
    accepted_inputs["proteome_index"] = {"s": [gene]}
    accepted_inputs["taxonomy"] = {gene.id: GeneTaxonomy(gene.id, "s", 300, 400)}
    assert joint.eligible_joint_pairs(**accepted_inputs) == []
    accepted_inputs["query"].all_porf_ids = [gene.id]
    assert len(joint.eligible_joint_pairs(**accepted_inputs)) == 1


@pytest.mark.parametrize("nested", [False, True])
def test_third_boundary_blocks_accepted_pair(accepted_inputs: dict, nested: bool) -> None:
    third = child("third", 200 if nested else 1100, 500 if nested else 1400)
    accepted_inputs["boundaries"].append(third)
    accepted_inputs["results"].append(rejected(third))
    assert joint.eligible_joint_pairs(**accepted_inputs) == []


@pytest.mark.parametrize("failure", [None, "tier", "class", "second_child_class", "missing_child", "gate"])
def test_accepted_union_replaces_only_after_native_checks(failure: str | None, tmp_path: Path) -> None:
    import copy
    import json

    from virosync.orchestration._flows.single_genome.phase3 import _replace_accepted_joints
    from virosync.pipeline.phase3.evidence_synthesizer import calculate_eve_confidence

    children = (child("a", 0, 6000), child("b", 10000, 16000))
    originals = [rejected(boundary) for boundary in children]
    for original in originals:
        original.confidence_tier = "MEDIUM"
    originals[0].hallmark_evidence = [hit("a", "GVOGm0054", 1000, 1300)]
    originals[1].hallmark_evidence = [hit("b", "GVOGm0003", 11000, 11300)]
    boundary = joint.joint_boundary(children, rule="accepted_adjacent_v1")
    result = rejected(boundary)
    result.hallmark_evidence = [hit("a", "GVOGm0054", 1000, 1300), hit("b", "GVOGm0003", 11000, 11300)]
    result.final_confidence = calculate_eve_confidence(result, 0)
    assert result.final_confidence >= 0.7
    result.confidence_tier = "HIGH"
    if failure == "tier":
        originals[0].confidence_tier = "HIGH"
        result.confidence_tier = "MEDIUM"
    elif failure == "class":
        result.region_classification = "MIRUS"
    elif failure == "second_child_class":
        originals[1].region_classification = "MIRUS"
    elif failure == "missing_child":
        result.hallmark_evidence = [hit("a", "GVOGm0054", 1000, 1300), hit("gap", "GVOGm0003", 7000, 7300)]
    elif failure == "gate":
        result.confidence_tier = "LOW"
        originals[0].confidence_tier = originals[1].confidence_tier = "LOW"
    before = copy.deepcopy(originals)

    selected = _replace_accepted_joints(
        originals,
        [result],
        [(boundary, children)],
        ablation_id=AblationID.A0,
        accepted_ids=frozenset(item.eve_id for item in originals),
    )

    assert originals == before
    if failure:
        assert selected == originals
    else:
        assert selected == [result]
        assert result.confirmed_seeds == children[0].confirmed_seeds + children[1].confirmed_seeds
        provenance = json.loads(result.joint_seed_provenance)
        assert provenance["rule"] == "accepted_adjacent_v1"
        assert [item["candidate_id"] for item in provenance["children"]] == [
            boundary_candidate_id(child) for child in children
        ]
        assert provenance["children"][0]["qualified_families"] == ["PolB"]
        assert provenance["children"][1]["qualified_families"] == ["MCP"]
        from virosync.pipeline.phase3.output_generator import OutputGenerator

        writer = OutputGenerator(output_dir=tmp_path)
        exported = json.loads(writer.write_evidence_profiles(selected).read_text())[result.eve_id]
        assert exported["joint_seed_provenance"] == result.joint_seed_provenance
        assert len(exported["confirmed_seeds"]) == 2


@pytest.fixture
def repeat_inputs() -> dict:
    import random

    from test_terminal_repeats import _Marker

    from virosync.pipeline.phase2.boundary_diamond import GenomeDiamondQuery, SeedGeneMapping

    rng = random.Random(927)
    sequence = "".join(rng.choices("ACGT", k=5000))
    return {
        "children": (child("a", 1000, 2000), child("b", 2600, 3500)),
        "sequence": sequence,
        "markers": [_Marker("s", 1300, 1500, "a"), _Marker("s", 2800, 3000, "b")],
        "extension_bp": 500,
        "query": GenomeDiamondQuery(
            all_porf_ids=["a", "b"],
            seed_gene_mappings={
                "a": SeedGeneMapping("a", "s", 1000, 2000, flank_start_bp=0, flank_end_bp=5000),
                "b": SeedGeneMapping("b", "s", 2600, 3500, flank_start_bp=0, flank_end_bp=5000),
            },
        ),
        "proteome_index": {"s": [pORF("a", "s", 1300, 1500), pORF("b", "s", 2800, 3000)]},
        "taxonomy": {"a": GeneTaxonomy("a", "s", 1300, 1500), "b": GeneTaxonomy("b", "s", 2800, 3000)},
    }


@pytest.mark.parametrize("arm_start,gap_start,allowed", [(960, 2200, False), (1350, 1800, True)])
def test_fresh_direct_repeat_requires_complete_marker_span(
    repeat_inputs: dict, arm_start: int, gap_start: int, allowed: bool
) -> None:
    sequence = repeat_inputs["sequence"]
    arm = sequence[arm_start : arm_start + 60]
    repeat_inputs["sequence"] = sequence[:gap_start] + arm + sequence[gap_start + 60 :]

    assert (joint.assess_joint_repeat_scan(**repeat_inputs) is not None) is allowed


@pytest.mark.parametrize("position,allowed", [(700, True), (1800, False)])
def test_only_gap_side_ambiguous_window_withholds(repeat_inputs: dict, position: int, allowed: bool) -> None:
    sequence = repeat_inputs["sequence"]
    repeat_inputs["sequence"] = sequence[:position] + "N" * 70 + sequence[position + 70 :]

    assert (joint.assess_joint_repeat_scan(**repeat_inputs) is not None) is allowed


@pytest.mark.parametrize("limit", ["MAX_CANDIDATES", "MAX_SEED_PAIRINGS"])
def test_direct_repeat_candidate_limit_withholds(
    repeat_inputs: dict, monkeypatch: pytest.MonkeyPatch, limit: str
) -> None:
    from virosync.pipeline.phase2 import repeat_evidence

    sequence = repeat_inputs["sequence"]
    repeat_inputs["sequence"] = sequence[:2200] + sequence[960:1020] + sequence[2260:]
    monkeypatch.setattr(repeat_evidence, limit, 0)

    assert not joint.assess_joint_repeat_scan(**repeat_inputs)


def test_far_side_contig_clipping_does_not_withhold(repeat_inputs: dict) -> None:
    repeat_inputs["children"] = (child("a", 100, 2000), child("b", 2600, 3500))

    assert joint.assess_joint_repeat_scan(**repeat_inputs)


def test_missing_repeat_query_coverage_withholds(repeat_inputs: dict) -> None:
    repeat_inputs["query"] = None

    assert not joint.assess_joint_repeat_scan(**repeat_inputs)


@pytest.mark.parametrize("accepted,position,expected", [(True, 1800, True), (False, 1800, False), (True, 2200, False)])
def test_prepared_joint_records_only_accepted_child_ambiguity(
    repeat_inputs: dict,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    accepted: bool,
    position: int,
    expected: bool,
) -> None:
    """The real repeat scans permit child uncertainty but never an unresolved bridge."""
    import json
    from dataclasses import replace

    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome import phase3

    sequence = repeat_inputs["sequence"]
    sequence = sequence[:position] + "N" * 70 + sequence[position + 70 :]
    raw = tmp_path / "raw.fa"
    raw.write_text(f">s\n{sequence}\n", encoding="utf-8")
    children = repeat_inputs["children"]
    children[0].confirmed_seeds = (replace(children[0].confirmed_seeds[0], start=1300, end=1500),)
    children[1].confirmed_seeds = (replace(children[1].confirmed_seeds[0], start=2800, end=3000),)
    accepted_ids = frozenset(boundary_candidate_id(c) for c in children) if accepted else frozenset()
    config = PipelineConfig()
    config.phase1.extension_kb = 1
    monkeypatch.setattr(phase3, "refine_boundaries_with_terminal_repeats", lambda boundaries, **kwargs: boundaries)
    monkeypatch.setattr(phase3, "refresh_boundary_classification", lambda *args, **kwargs: None)
    monkeypatch.setattr(phase3, "_recalculate_boundary_composition", lambda *args, **kwargs: None)
    monkeypatch.setattr(phase3, "_annotate_boundary_repeats", lambda *args, **kwargs: None)

    prepared = phase3._prepare_joint_boundaries(
        [children],
        raw_genome_path=raw,
        masked_path=raw,
        validated_markers=repeat_inputs["markers"],
        boundary_diamond_query=repeat_inputs["query"],
        boundary_taxonomy_map=repeat_inputs["taxonomy"],
        proteome_index=repeat_inputs["proteome_index"],
        config=config,
        accepted_ids=accepted_ids,
    )

    assert bool(prepared) is expected
    if expected:
        boundary, original_children = prepared[0]
        assert original_children == children
        assert boundary.confirmed_seeds == children[0].confirmed_seeds + children[1].confirmed_seeds
        assessment = json.loads(boundary.joint_seed_provenance)["repeat_assessment"]
        assert assessment["status"] == "accepted_child_ambiguity"
        assert assessment["new_gap"] == [2000, 2600]
        assert assessment["inter_core_ambiguous_bp"] == 70
        assert assessment["gap_sides"][0]["reasons"] == ["ambiguous_bases", "midpoint_partition"]
        assert assessment["gap_sides"][0]["ambiguous_bp"] == 70
        assert assessment["gap_sides"][1]["status"] == "incomplete"
        assert assessment["gap_sides"][1]["ambiguous_bp"] == 70


@pytest.mark.parametrize("arm_start,copy_start", [(960, 2200), (2300, 3700)])
def test_accepted_child_ambiguity_keeps_both_separating_repeat_vetoes(
    repeat_inputs: dict, arm_start: int, copy_start: int
) -> None:
    """A positive pair on either child vetoes even when the left scan is incomplete."""
    sequence = repeat_inputs["sequence"]
    sequence = sequence[:copy_start] + sequence[arm_start : arm_start + 60] + sequence[copy_start + 60 :]
    repeat_inputs["sequence"] = sequence[:1800] + "N" * 70 + sequence[1870:]
    accepted_ids = frozenset(boundary_candidate_id(c) for c in repeat_inputs["children"])

    assert joint.assess_joint_repeat_scan(**repeat_inputs, accepted_ids=accepted_ids) is None


@pytest.mark.parametrize("reason", ["candidate_limit", "seed_pairing_limit", "taxonomy_clipped", "window_limit"])
def test_accepted_ambiguity_keeps_other_incomplete_reasons(
    repeat_inputs: dict, monkeypatch: pytest.MonkeyPatch, reason: str
) -> None:
    """An ambiguity exception must not erase a concurrent native limit or coverage failure."""
    from dataclasses import replace

    from virosync.pipeline.phase2.repeat_evidence import RepeatEvidence

    sequence = repeat_inputs["sequence"]
    repeat_inputs["sequence"] = sequence[:1800] + "N" * 70 + sequence[1870:]
    accepted_ids = frozenset(boundary_candidate_id(c) for c in repeat_inputs["children"])
    native_assess = joint.assess_repeat_evidence

    def limited_assessment(*args: object, **kwargs: object) -> RepeatEvidence:
        evidence = native_assess(*args, **kwargs)
        return replace(evidence, right=replace(evidence.right, reasons=(*evidence.right.reasons, reason)))

    monkeypatch.setattr(joint, "assess_repeat_evidence", limited_assessment)

    assert joint.assess_joint_repeat_scan(**repeat_inputs, accepted_ids=accepted_ids) is None


def test_accepted_ambiguity_outside_children_still_withholds(repeat_inputs: dict) -> None:
    """An endpoint window may extend past the sibling; that exterior DNA is not exempt."""
    sequence = repeat_inputs["sequence"]
    repeat_inputs["sequence"] = sequence[:3550] + "N" * 70 + sequence[3620:]
    repeat_inputs["extension_bp"] = 1600
    accepted_ids = frozenset(boundary_candidate_id(c) for c in repeat_inputs["children"])

    assert joint.assess_joint_repeat_scan(**repeat_inputs, accepted_ids=accepted_ids) is None


def test_accepted_ambiguity_cannot_hide_low_complexity(repeat_inputs: dict) -> None:
    """An N-containing window still needs informative resolved bases."""
    sequence = repeat_inputs["sequence"]
    repeat_inputs["sequence"] = sequence[:1500] + "A" * 300 + "N" * 70 + "A" * 630 + sequence[2500:]
    accepted_ids = frozenset(boundary_candidate_id(c) for c in repeat_inputs["children"])

    assert joint.assess_joint_repeat_scan(**repeat_inputs, accepted_ids=accepted_ids) is None


def test_unavailable_cellular_lineage_withholds(accepted_inputs: dict) -> None:
    accepted_inputs["host_model"] = HostSignatureModel(
        token_weights={"plastid": 1.0}, lineage_weights={"plastid|plastid|plastid": 1.0}
    )

    assert joint.eligible_joint_pairs(**accepted_inputs) == []


@pytest.fixture
def padded_repeat_inputs(repeat_inputs: dict) -> dict:
    from virosync.pipeline.phase2.boundary_diamond import SeedGeneMapping

    genes = [pORF("a", "s", 1300, 1500), pORF("b", "s", 2800, 3000)]
    repeat_inputs["proteome_index"] = {"s": genes}
    repeat_inputs["taxonomy"] = {gene.id: GeneTaxonomy(gene.id, "s", gene.start, gene.end) for gene in genes}
    query = repeat_inputs["query"]
    query.all_porf_ids = [gene.id for gene in genes]
    query.seed_gene_mappings = {
        "a": SeedGeneMapping("a", "s", 1000, 2000, flank_start_bp=1100, flank_end_bp=2500),
        "b": SeedGeneMapping("b", "s", 2600, 3500, flank_start_bp=2100, flank_end_bp=4000),
    }
    return repeat_inputs


def test_gene_free_padding_uses_proven_hull_without_expanding_search(
    padded_repeat_inputs: dict, monkeypatch: pytest.MonkeyPatch
) -> None:
    native = joint.assess_repeat_evidence
    scans = []

    def record(*args: object, **kwargs: object) -> object:
        scan = native(*args, **kwargs)
        scans.append((kwargs, scan))
        return scan

    monkeypatch.setattr(joint, "assess_repeat_evidence", record)
    assert joint.assess_joint_repeat_scan(**padded_repeat_inputs)
    arguments, scan = scans[0]
    assert (1000, 3500) in arguments["coverage_intervals"]
    assert (arguments["start"], arguments["end"]) == (1000, 2000)
    assert scan.left.window_start == 1000
    assert "taxonomy_clipped" in scan.left.reasons
    assert scan.right.status == "completed"
    assert len(scans) == 2


def test_child_bounding_direct_repeat_in_padding_still_withholds(padded_repeat_inputs: dict) -> None:
    sequence = padded_repeat_inputs["sequence"]
    padded_repeat_inputs["sequence"] = sequence[:2200] + sequence[1020:1080] + sequence[2260:]

    assert not joint.assess_joint_repeat_scan(**padded_repeat_inputs)


@pytest.mark.parametrize("missing", ["query", "taxonomy", "mapping", "all_genes"])
def test_repeat_guard_rechecks_hull_completeness(padded_repeat_inputs: dict, missing: str) -> None:
    if missing == "query":
        padded_repeat_inputs["query"].all_porf_ids.remove("a")
    elif missing == "taxonomy":
        del padded_repeat_inputs["taxonomy"]["a"]
    elif missing == "mapping":
        del padded_repeat_inputs["query"].seed_gene_mappings["a"]
    else:
        padded_repeat_inputs["proteome_index"] = {}

    assert not joint.assess_joint_repeat_scan(**padded_repeat_inputs)


def test_hull_proof_does_not_fill_gap_side_coverage_beyond_hull(
    padded_repeat_inputs: dict, monkeypatch: pytest.MonkeyPatch
) -> None:
    native = joint.assess_repeat_evidence
    scans = []
    padded_repeat_inputs["extension_bp"] = 2000
    padded_repeat_inputs["query"].seed_gene_mappings["b"].flank_end_bp = 3500

    def record(*args: object, **kwargs: object) -> object:
        scan = native(*args, **kwargs)
        scans.append(scan)
        return scan

    monkeypatch.setattr(joint, "assess_repeat_evidence", record)
    assert not joint.assess_joint_repeat_scan(**padded_repeat_inputs)
    assert "taxonomy_clipped" in scans[0].right.reasons
    assert scans[0].right.window_end == 3500


@pytest.mark.parametrize("searched", [True, False])
def test_accepted_envelope_hole_requires_queried_gap_gene(accepted_inputs: dict, searched: bool) -> None:
    gap = pORF("hole", "s", 1200, 1300)
    accepted_inputs["query"].seed_gene_mappings["left"].flank_end_bp = 1000
    accepted_inputs["query"].seed_gene_mappings["right"].flank_start_bp = 1500
    accepted_inputs["query"].all_porf_ids = ["hole"] if searched else []
    accepted_inputs["proteome_index"] = {"s": [gap]}
    accepted_inputs["taxonomy"] = {"hole": GeneTaxonomy("hole", "s", 1200, 1300)}

    assert bool(joint.eligible_joint_pairs(**accepted_inputs)) is searched
