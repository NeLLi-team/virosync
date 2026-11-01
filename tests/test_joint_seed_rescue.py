"""Joint rescue preserves rejected seed ownership and biological barriers."""

from dataclasses import replace

import pytest

from virosync.ablation import AblationID
from virosync.pipeline.host_signatures import HostSignatureModel
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed
from virosync.pipeline.phase2.boundary_diamond import GenomeDiamondQuery, SeedGeneMapping, pORF
from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary, boundary_candidate_id
from virosync.pipeline.phase3.evidence_synthesizer import VerificationResult
from virosync.pipeline.phase3.joint_seed_rescue import eligible_joint_pairs


def child(name: str, start: int, end: int) -> RefinedBoundary:
    return RefinedBoundary(
        "s",
        start,
        end,
        seed_id=name,
        original_start=start,
        original_end=end,
        confirmed_seeds=(ConfirmedSeed(name, "s", start, end),),
        left_stop_reason="support_gap",
        right_stop_reason="support_gap",
    )


def rejected(boundary: RefinedBoundary) -> VerificationResult:
    return VerificationResult(
        boundary_candidate_id(boundary),
        "s",
        boundary.start,
        boundary.end,
        confirmed_seeds=boundary.confirmed_seeds,
        confidence_tier="LOW",
        region_classification="NCLDV",
        likely_family="NCLDV",
    )


@pytest.fixture
def pair_inputs() -> dict:
    boundaries = [child("left", 100, 1000), child("right", 1500, 2500)]
    query = GenomeDiamondQuery(
        seed_gene_mappings={
            "left": SeedGeneMapping("left", "s", 100, 1000, flank_start_bp=0, flank_end_bp=1800),
            "right": SeedGeneMapping("right", "s", 1500, 2500, flank_start_bp=800, flank_end_bp=3000),
        }
    )
    return dict(
        boundaries=boundaries,
        results=[rejected(b) for b in boundaries],
        query=query,
        proteome_index={},
        taxonomy={},
        host_model=HostSignatureModel(token_weights={"host": 1.0}),
        host_threshold=0.35,
        accepted_ids=frozenset(),
    )


def test_adjacent_rejected_children_share_context_without_distance_cutoff(pair_inputs: dict) -> None:
    pairs = eligible_joint_pairs(**pair_inputs)
    assert len(pairs) == 1
    assert pairs[0] == tuple(pair_inputs["boundaries"])


def test_accepted_child_is_never_borrowed(pair_inputs: dict) -> None:
    pair_inputs["results"][0].confidence_tier = "HIGH"
    assert eligible_joint_pairs(**pair_inputs) == []


@pytest.mark.parametrize("family", ["MIRUS", "UNKNOWN"])
def test_incompatible_or_unknown_class_abstains(pair_inputs: dict, family: str) -> None:
    pair_inputs["results"][1].region_classification = family
    pair_inputs["results"][1].likely_family = family
    assert eligible_joint_pairs(**pair_inputs) == []


def test_recorded_host_stop_is_preserved(pair_inputs: dict) -> None:
    pair_inputs["boundaries"][0].right_stop = 1200
    assert eligible_joint_pairs(**pair_inputs) == []


def test_nonoverlapping_search_envelopes_abstain(pair_inputs: dict) -> None:
    pair_inputs["query"].seed_gene_mappings["left"].flank_end_bp = 800
    assert eligible_joint_pairs(**pair_inputs) == []


def test_unsearched_gap_orf_abstains(pair_inputs: dict) -> None:
    pair_inputs["proteome_index"] = {"s": [pORF("gap", "s", 1100, 1300, "+")]}
    assert eligible_joint_pairs(**pair_inputs) == []


def test_two_neighbors_are_assessed_before_competition(pair_inputs: dict) -> None:
    third = child("third", 2700, 3500)
    pair_inputs["boundaries"].append(third)
    pair_inputs["results"].append(rejected(third))
    pair_inputs["query"].seed_gene_mappings["third"] = SeedGeneMapping(
        "third",
        "s",
        2700,
        3500,
        flank_start_bp=2000,
        flank_end_bp=4000,
    )
    assert len(eligible_joint_pairs(**pair_inputs)) == 2


def test_intervening_seed_blocks_pair(pair_inputs: dict) -> None:
    third = child("third", 1200, 1400)
    pair_inputs["boundaries"].append(third)
    pair_inputs["results"].append(rejected(third))
    assert eligible_joint_pairs(**pair_inputs) == []


def test_missing_query_and_active_acceptance_abstain(pair_inputs: dict) -> None:
    assert eligible_joint_pairs(**{**pair_inputs, "query": None}) == []
    identity = pair_inputs["results"][0].eve_id
    assert eligible_joint_pairs(**{**pair_inputs, "accepted_ids": frozenset({identity})}) == []


@pytest.mark.parametrize("count,viral_last,blocked", [(2, False, False), (3, False, True), (3, True, False)])
def test_gap_uses_three_host_rule(pair_inputs: dict, monkeypatch, count: int, viral_last: bool, blocked: bool) -> None:
    from types import SimpleNamespace

    from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy
    from virosync.pipeline.phase3 import joint_seed_rescue

    genes = [pORF(f"g{i}", "s", 1050 + i * 100, 1100 + i * 100) for i in range(count)]
    taxonomy = {g.id: GeneTaxonomy(g.id, "s", g.start, g.end) for g in genes}
    if viral_last:
        taxonomy[genes[-1].id].top10_prefixes = ["NCLDV__"]
        taxonomy[genes[-1].id].top10_pidents = [30]
    monkeypatch.setattr(joint_seed_rescue, "host_lineage_share", lambda gene, model: 0.5)
    pair_inputs.update(
        proteome_index={"s": genes}, taxonomy=taxonomy, host_model=SimpleNamespace(token_weights={"host": 1})
    )
    pair_inputs["query"].all_porf_ids = list(taxonomy)
    assert bool(eligible_joint_pairs(**pair_inputs)) is not blocked


def test_native_joint_rescoring_preserves_cores_and_output_resume(tmp_path) -> None:
    import json

    from test_eve_confidence_scoring import _supplementary_gene

    from virosync.orchestration._flows.single_genome.phase3 import _replace_accepted_joints
    from virosync.orchestration._flows.single_genome.phase_state import (
        phase2_state_from_document,
        phase2_state_to_document,
    )
    from virosync.pipeline.phase3.evidence_synthesizer import EvidenceSynthesizer, EvidenceSynthesizerConfig
    from virosync.pipeline.phase3.joint_seed_rescue import joint_boundary, joint_result_supported
    from virosync.pipeline.phase3.output_generator import OutputGenerator

    children = (child("left", 0, 2000), child("right", 3000, 7000))
    children = tuple(
        replace(c, scaffold="scaf", confirmed_seeds=(ConfirmedSeed(c.seed_id, "scaf", c.start, c.end),))
        for c in children
    )
    boundary = joint_boundary(children)
    boundary.predicted_family = "NCLDV"
    restored = phase2_state_from_document(phase2_state_to_document([boundary]))[0]
    synthesizer = EvidenceSynthesizer(
        config=EvidenceSynthesizerConfig(use_boltz=False, use_phylogenetic_validation=False, use_tmvec_database=False),
        work_dir=tmp_path,
    )
    records = [
        _supplementary_gene("a", 100, 200),
        _supplementary_gene("b", 3100, 3200),
        _supplementary_gene("c", 3500, 3600),
    ]
    result = synthesizer.verify_eve(restored, [], gene_taxonomy_records=records, gene_taxonomy_summary={"total": 3})
    assert result.supplementary_viral_evidence.eligible
    assert result.supplementary_viral_evidence.start == 0
    assert result.supplementary_viral_evidence.end == 7000
    assert joint_result_supported(result, children)
    assert result.confirmed_seeds == restored.confirmed_seeds
    assert result.joint_seed_provenance == boundary.joint_seed_provenance
    provenance = json.loads(result.joint_seed_provenance)
    assert [c["candidate_id"] for c in provenance["children"]] == [boundary_candidate_id(c) for c in children]
    assert provenance["children"][0]["right_stop_reason"] == "support_gap"
    originals = [rejected(c) for c in children]
    assert _replace_accepted_joints(originals, [result], [(restored, children)], ablation_id=AblationID.A0) == [result]
    successful = json.loads(result.joint_seed_provenance)
    assert successful["children"][0]["final_confidence"] == originals[0].final_confidence
    assert successful["children"][0]["confidence_tier"] == "LOW"
    assert successful["children"][0]["effective_class"] == "NCLDV"
    path = OutputGenerator(output_dir=tmp_path).write_evidence_profiles([result])
    exported = json.loads(path.read_text())[result.eve_id]
    assert exported["joint_seed_provenance"] == result.joint_seed_provenance
    assert len(exported["confirmed_seeds"]) == 2
    import csv

    detailed = OutputGenerator(output_dir=tmp_path).write_predictions_detailed_tsv([result])
    with detailed.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["joint_seed_provenance"] == result.joint_seed_provenance
    result.supplementary_viral_evidence = replace(result.supplementary_viral_evidence, eligible=False)
    assert _replace_accepted_joints(originals, [result], [(restored, children)], ablation_id=AblationID.A0) == originals


def test_one_core_cannot_donate_all_selected_witnesses() -> None:
    from test_eve_confidence_scoring import _supplementary_candidate

    from virosync.pipeline.phase3.evidence_synthesizer import calculate_eve_confidence
    from virosync.pipeline.phase3.joint_seed_rescue import joint_result_supported

    result = _supplementary_candidate()
    calculate_eve_confidence(result, 0)
    assert result.supplementary_viral_evidence.eligible
    assert not joint_result_supported(result, (child("a", 0, 1000), child("b", 2000, 6000)))


@pytest.mark.parametrize("models,blocked", [(("GVOGm0003", "GVOGm0760"), True), (("GVOGm0760", "GVOGm0760"), False)])
def test_two_native_multi_marker_groups_stay_separate(
    pair_inputs: dict, models: tuple[str, str], blocked: bool
) -> None:
    from test_native_hallmark_evidence import hit

    from virosync.pipeline.hallmark_evidence import find_spatial_witnesses

    left, right = pair_inputs["results"]
    left.hallmark_evidence = [hit("a", models[0], 200, 400), hit("b", models[1], 500, 700)]
    right.hallmark_evidence = [hit("c", models[0], 1600, 1800), hit("d", models[1], 1900, 2100)]
    left.hallmark_witnesses = find_spatial_witnesses(left.hallmark_evidence)
    right.hallmark_witnesses = find_spatial_witnesses(right.hallmark_evidence)
    assert bool(eligible_joint_pairs(**pair_inputs)) is not blocked


def test_native_internal_repeat_veto_and_outer_pair_allowed() -> None:
    from test_terminal_repeats import _adjacent_repeat_sequence, _Marker, _planted_repeat

    from virosync.pipeline.phase3.joint_seed_rescue import has_clear_internal_boundary_scan

    sequence, _ = _adjacent_repeat_sequence()
    children = (child("a", 150, 400), child("b", 570, 820))
    markers = [_Marker("s", 220, 270, "a"), _Marker("s", 640, 690, "b")]
    assert not has_clear_internal_boundary_scan(children, sequence, markers, 5000)
    outer_sequence, _ = _planted_repeat()
    close_children = (child("a", 170, 240), child("b", 260, 350))
    close_markers = [_Marker("s", 190, 210, "a"), _Marker("s", 280, 310, "b")]
    assert has_clear_internal_boundary_scan(close_children, outer_sequence, close_markers, 5000)
    assert not has_clear_internal_boundary_scan(close_children, None, close_markers, 5000)
    assert not has_clear_internal_boundary_scan(close_children, outer_sequence, [], 5000)


def test_pairing_cap_is_unknown_and_blocks_rescue(monkeypatch) -> None:
    from test_terminal_repeats import _Marker, _planted_repeat

    from virosync.pipeline.phase2 import terminal_repeats
    from virosync.pipeline.phase3.joint_seed_rescue import has_clear_internal_boundary_scan

    sequence, _ = _planted_repeat()
    monkeypatch.setattr(terminal_repeats, "_MAX_SEED_PAIRINGS", 0)
    scan = terminal_repeats.find_terminal_inverted_repeats(
        sequence, scan_start=0, scan_end=len(sequence), marker_start=190, marker_end=210
    )
    assert scan.pairing_capped
    assert not has_clear_internal_boundary_scan(
        (child("a", 170, 240), child("b", 260, 350)),
        sequence,
        [_Marker("s", 190, 210, "a"), _Marker("s", 280, 310, "b")],
        5000,
    )


def test_native_joint_finalization_recomputes_class_composition_and_id(tmp_path, monkeypatch) -> None:
    import os
    import sys
    from pathlib import Path

    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome.phase3 import _prepare_joint_boundaries
    from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit

    monkeypatch.setenv("PATH", str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"])
    genome = tmp_path / "genome.fna"
    genome.write_text(">s\n" + "A" * 8000 + "\n")
    children = (child("a", 0, 2000), child("b", 3000, 7000))
    children[0].gc_deviation = 999
    children[0].candidate_id = "stale_child_identity"
    markers = [
        ValidatedMarkerHit(
            "a", "s", 500, 800, "+", "GVOGm0760", 100, 1e-10, "validated", "NCLDV__", "NCLDV__a", 35, 100, 1, 0, 0, 0, 1
        ),
        ValidatedMarkerHit(
            "b",
            "s",
            4000,
            4300,
            "+",
            "VS000151",
            100,
            1e-10,
            "validated",
            "NCLDV__",
            "NCLDV__b",
            35,
            100,
            1,
            0,
            0,
            0,
            1,
        ),
    ]
    query = GenomeDiamondQuery(
        seed_gene_mappings={
            "a": SeedGeneMapping("a", "s", 0, 2000, flank_start_bp=0, flank_end_bp=5000),
            "b": SeedGeneMapping("b", "s", 3000, 7000, flank_start_bp=1000, flank_end_bp=8000),
        }
    )
    prepared = _prepare_joint_boundaries(
        [children],
        raw_genome_path=genome,
        masked_path=genome,
        validated_markers=markers,
        boundary_diamond_query=query,
        boundary_taxonomy_map={},
        proteome_index={},
        config=PipelineConfig(),
    )
    assert len(prepared) == 1
    boundary, retained_children = prepared[0]
    assert retained_children == children
    assert (boundary.start, boundary.end) == (0, 7000)
    assert boundary.candidate_id == "EVE_s_0-7000"
    assert boundary.predicted_family == "NCLDV"
    assert boundary.gc_deviation != 999
    assert boundary.repeat_evidence is not None
    assert boundary.confirmed_seeds == children[0].confirmed_seeds + children[1].confirmed_seeds

    # Execute the same annotation and batched scoring helper called by Phase 3.
    import logging

    from virosync.orchestration._flows.single_genome.phase3 import _verify_native_boundaries
    from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy

    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">a # 501 # 800 # 1 # ID=1\n" + "M" * 100 + "\n"
        ">b # 4001 # 4300 # 1 # ID=2\n" + "M" * 100 + "\n"
        ">c # 4501 # 4800 # 1 # ID=3\n" + "M" * 100 + "\n"
    )
    genes = [pORF("a", "s", 500, 800), pORF("b", "s", 4000, 4300), pORF("c", "s", 4500, 4800)]
    taxonomy = {
        g.id: GeneTaxonomy(
            g.id,
            "s",
            g.start,
            g.end,
            top10_targets=["NCLDV__ref"],
            top10_prefixes=["NCLDV__"],
            top10_pidents=[30],
            top10_evalues=[1e-8],
            top10_bits=[120],
            top10_qcovhsp=[80],
        )
        for g in genes
    }
    config = PipelineConfig()
    config.compute.threads = config.compute.max_threads = 1
    scored = _verify_native_boundaries(
        [boundary],
        boundary_taxonomy_map=taxonomy,
        validated_markers=markers,
        proteome_index={"s": genes},
        boundary_diamond_query=query,
        interproscan_map=None,
        use_precomputed_taxonomy=True,
        proteome_path=proteome,
        config=config,
        logger=logging.getLogger(__name__),
        hallmark_evidence=[],
        hallmark_witnesses=[],
        output_dir=tmp_path,
        jelly_roll_map=None,
        masked_path=genome,
        host_signatures=set(),
        host_signature_model_payload=None,
        precomputed_tmvec=None,
        tmvec_device="cpu",
    )
    assert len(scored) == 1
    assert scored[0].supplementary_viral_evidence.eligible
    assert scored[0].gene_count == 3
    assert scored[0].confirmed_seeds == boundary.confirmed_seeds
    assert scored[0].joint_seed_provenance == boundary.joint_seed_provenance

    from virosync.orchestration._flows.single_genome import phase3

    def unavailable_repeat_tool(*args: object, **kwargs: object) -> None:
        raise RuntimeError("GenomeTools unavailable")

    monkeypatch.setattr(phase3, "refine_boundaries_with_terminal_repeats", unavailable_repeat_tool)
    assert (
        _prepare_joint_boundaries(
            [children],
            raw_genome_path=genome,
            masked_path=genome,
            validated_markers=markers,
            boundary_diamond_query=query,
            boundary_taxonomy_map={},
            proteome_index={},
            config=config,
        )
        == []
    )


def test_internal_repeat_endpoint_in_seed_padding_still_separates() -> None:
    from test_terminal_repeats import _adjacent_repeat_sequence, _Marker

    from virosync.pipeline.phase3.joint_seed_rescue import has_clear_internal_boundary_scan

    sequence, _ = _adjacent_repeat_sequence()
    # First repeat ends at 420, inside the right seed's padded prefix.
    children = (child("a", 100, 350), child("b", 400, 900))
    assert not has_clear_internal_boundary_scan(
        children,
        sequence,
        [_Marker("s", 220, 270, "a"), _Marker("s", 640, 690, "b")],
        5000,
    )
    assert children[1].confirmed_seeds == (ConfirmedSeed("b", "s", 400, 900),)


@pytest.mark.parametrize("both_pass", [False, True])
def test_competition_counts_only_fully_supported_native_results(both_pass: bool) -> None:
    from test_eve_confidence_scoring import _supplementary_gene

    from virosync.orchestration._flows.single_genome.phase3 import _replace_accepted_joints
    from virosync.pipeline.phase3.evidence_synthesizer import calculate_eve_confidence
    from virosync.pipeline.phase3.joint_seed_rescue import joint_boundary

    children = (child("a", 0, 2000), child("b", 3000, 7000), child("c", 8000, 10000))
    children = tuple(
        replace(c, scaffold="scaf", confirmed_seeds=(ConfirmedSeed(c.seed_id, "scaf", c.start, c.end),))
        for c in children
    )
    pairs = [children[:2], children[1:]]
    joints = [(joint_boundary(pair), pair) for pair in pairs]
    first = VerificationResult(
        boundary_candidate_id(joints[0][0]),
        "scaf",
        0,
        7000,
        confirmed_seeds=joints[0][0].confirmed_seeds,
        likely_family="NCLDV",
        region_classification="NCLDV",
        confidence_tier="LOW",
        gene_taxonomy_records=[
            _supplementary_gene("a", 100, 200),
            _supplementary_gene("b", 3100, 3200),
            _supplementary_gene("c", 3500, 3600),
        ],
    )
    second = VerificationResult(
        boundary_candidate_id(joints[1][0]),
        "scaf",
        3000,
        10000,
        confirmed_seeds=joints[1][0].confirmed_seeds,
        likely_family="NCLDV",
        region_classification="NCLDV",
        confidence_tier="LOW",
        gene_taxonomy_records=[
            _supplementary_gene("b", 3100, 3200),
            _supplementary_gene("c", 3500, 3600),
            _supplementary_gene("d", 8500, 8600, qcov=80 if both_pass else 10),
        ],
    )
    calculate_eve_confidence(first, 0)
    calculate_eve_confidence(second, 0)
    originals = [rejected(c) for c in children]
    selected = _replace_accepted_joints(originals, [first, second], joints, ablation_id=AblationID.A0)
    assert selected == originals if both_pass else selected == [originals[2], first]


def test_failed_required_repeat_scan_abstains(monkeypatch) -> None:
    from test_terminal_repeats import _Marker

    from virosync.pipeline.phase3 import joint_seed_rescue

    def fail_scan(*args: object, **kwargs: object) -> None:
        raise RuntimeError("required scan failed")

    monkeypatch.setattr(joint_seed_rescue, "find_terminal_inverted_repeats", fail_scan)
    assert not joint_seed_rescue.has_clear_internal_boundary_scan(
        (child("a", 100, 1000), child("b", 2000, 3000)),
        "A" * 4000,
        [_Marker("s", 400, 500, "a"), _Marker("s", 2400, 2500, "b")],
        5000,
    )


def test_supplementary_support_does_not_bypass_final_length_gate() -> None:
    from test_eve_confidence_scoring import _supplementary_candidate, _supplementary_gene

    from virosync.orchestration._flows.single_genome.phase3 import _replace_accepted_joints
    from virosync.pipeline.phase3.evidence_synthesizer import calculate_eve_confidence
    from virosync.pipeline.phase3.joint_seed_rescue import joint_boundary

    children = (child("a", 0, 1000), child("b", 2000, 5000))
    joint = joint_boundary(children)
    result = _supplementary_candidate(
        end=5000,
        records=[
            _supplementary_gene("a", 100, 200),
            _supplementary_gene("b", 2100, 2200),
            _supplementary_gene("c", 2500, 2600),
        ],
    )
    result.eve_id = boundary_candidate_id(joint)
    result.confirmed_seeds = joint.confirmed_seeds
    calculate_eve_confidence(result, 0)
    assert result.supplementary_viral_evidence.eligible
    originals = [rejected(c) for c in children]
    assert _replace_accepted_joints(originals, [result], [(joint, children)], ablation_id=AblationID.A0) == originals


@pytest.mark.parametrize("model", [None, {}])
def test_unavailable_host_model_abstains(pair_inputs: dict, model: dict | None) -> None:
    from virosync.pipeline.host_signatures import HostSignatureModel

    pair_inputs["host_model"] = None if model is None else HostSignatureModel()
    assert eligible_joint_pairs(**pair_inputs) == []


@pytest.mark.parametrize("gap_base", ["N", "R"])
def test_ambiguous_intercore_dna_abstains(gap_base: str) -> None:
    from test_terminal_repeats import _Marker

    from virosync.pipeline.phase3.joint_seed_rescue import has_clear_internal_boundary_scan

    sequence = "A" * 1400 + gap_base + "A" * 1599
    assert not has_clear_internal_boundary_scan(
        (child("a", 100, 1000), child("b", 2000, 2900)),
        sequence,
        [_Marker("s", 400, 500, "a"), _Marker("s", 2400, 2500, "b")],
        5000,
    )


def test_accepted_overlap_outside_the_pair_blocks_rescue(pair_inputs: dict) -> None:
    third = child("accepted", 0, 500)
    accepted = rejected(third)
    accepted.confidence_tier = "HIGH"
    pair_inputs["boundaries"].append(third)
    pair_inputs["results"].append(accepted)
    assert eligible_joint_pairs(**pair_inputs) == []


def test_joint_output_cannot_be_paired_recursively(pair_inputs: dict) -> None:
    pair_inputs["boundaries"][0].joint_seed_provenance = '{"rule":"rejected_adjacent_r2_v1"}'
    assert eligible_joint_pairs(**pair_inputs) == []


def test_active_selection_rejection_preserves_children_despite_raw_gate_pass() -> None:
    from test_eve_confidence_scoring import _supplementary_candidate

    from virosync.ablation import AblationID
    from virosync.orchestration._flows.single_genome.phase3 import _replace_accepted_joints
    from virosync.pipeline.phase3.acceptance_selection import select_phase3_acceptance
    from virosync.pipeline.phase3.evidence_synthesizer import calculate_eve_confidence
    from virosync.pipeline.phase3.joint_seed_rescue import joint_boundary, joint_result_supported
    from virosync.pipeline.phase3.output_generator import evaluate_v2_quality_gate

    children = (child("a", 0, 250), child("b", 450, 6000))
    children = tuple(
        replace(
            c,
            scaffold="scaf",
            confirmed_seeds=(ConfirmedSeed(c.seed_id, "scaf", c.start, c.end),),
            seed_sources=["frameshift_rescue"],
        )
        for c in children
    )
    boundary = joint_boundary(children)
    result = _supplementary_candidate()
    result.eve_id = boundary_candidate_id(boundary)
    result.confirmed_seeds = boundary.confirmed_seeds
    result.seed_sources = ["frameshift_rescue"]
    calculate_eve_confidence(result, 0)
    originals = [rejected(c) for c in children]

    assert joint_result_supported(result, children)
    assert evaluate_v2_quality_gate(result).kept
    assert not select_phase3_acceptance([result], AblationID.A0).canonical_results
    assert result.canonical_selection_outcome == "rescue_marker_excluded"
    assert _replace_accepted_joints(originals, [result], [(boundary, children)], ablation_id=AblationID.A0) == originals
