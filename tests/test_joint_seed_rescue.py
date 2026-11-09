"""Joint rescue preserves rejected seed ownership and biological barriers."""

from dataclasses import replace

import pytest

from virosync.ablation import AblationID
from virosync.pipeline.host_signatures import HostSignatureModel
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed
from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy, GenomeDiamondQuery, SeedGeneMapping, pORF
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


def repeat_context(length: int) -> dict:
    """Provide explicit searched coverage for planted sequence controls."""
    gene = pORF("coverage", "s", 200, 210)
    return {
        "query": GenomeDiamondQuery(
            all_porf_ids=[gene.id],
            seed_gene_mappings={
                name: SeedGeneMapping(name, "s", 0, length, flank_start_bp=0, flank_end_bp=length)
                for name in ("a", "b")
            },
        ),
        "proteome_index": {"s": [gene]},
        "taxonomy": {gene.id: GeneTaxonomy(gene.id, "s", gene.start, gene.end)},
    }


@pytest.fixture
def pair_inputs(monkeypatch) -> dict:
    from virosync.pipeline import host_signatures

    monkeypatch.setattr(
        host_signatures, "_TAXONOMY_LOOKUP", host_signatures.TaxonomyLabelLookup({"EUK__host": "EUK|host"})
    )
    boundaries = [child("left", 100, 1000), child("right", 1500, 2500)]
    query = GenomeDiamondQuery(
        all_porf_ids=["left_gene", "right_gene"],
        seed_gene_mappings={
            "left": SeedGeneMapping("left", "s", 100, 1000, flank_start_bp=0, flank_end_bp=1800),
            "right": SeedGeneMapping("right", "s", 1500, 2500, flank_start_bp=800, flank_end_bp=3000),
        },
    )
    return dict(
        boundaries=boundaries,
        results=[rejected(b) for b in boundaries],
        query=query,
        proteome_index={"s": [pORF("left_gene", "s", 200, 400), pORF("right_gene", "s", 1600, 1800)]},
        taxonomy={
            "left_gene": GeneTaxonomy("left_gene", "s", 200, 400),
            "right_gene": GeneTaxonomy("right_gene", "s", 1600, 1800),
        },
        host_model=HostSignatureModel(token_weights={"host": 1.0}, lineage_weights={"euk|host": 1.0}),
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


@pytest.mark.parametrize("family,allowed", [("MIRUS", True), ("UNKNOWN", False)])
def test_mixed_ncldv_mirus_allowed_but_unknown_abstains(pair_inputs: dict, family: str, allowed: bool) -> None:
    pair_inputs["results"][1].region_classification = family
    pair_inputs["results"][1].likely_family = family
    assert bool(eligible_joint_pairs(**pair_inputs)) is allowed


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
        proteome_index={"s": genes},
        taxonomy=taxonomy,
        host_model=SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0}),
    )
    pair_inputs["query"].all_porf_ids = list(taxonomy)
    assert bool(eligible_joint_pairs(**pair_inputs)) is not blocked


def test_a_signature_without_a_lineage_abstains(pair_inputs: dict) -> None:
    pair_inputs["host_model"] = HostSignatureModel(
        token_weights={"plastid": 1.0},
        lineage_weights={"plastid|plastid|plastid": 1.0},
    )

    assert eligible_joint_pairs(**pair_inputs) == []


@pytest.mark.parametrize("weight_mode", ["rank", "bitscore"])
@pytest.mark.parametrize(
    "target,count,viral_last,blocked",
    [
        ("EUK__HOST", 2, False, False),
        ("EUK__HOST", 3, False, True),
        ("EUK__HOST", 3, True, False),
        ("EUK__OTHER", 3, False, True),
        ("PLASTID__HOST", 3, False, False),
    ],
)
def test_gap_barrier_uses_the_real_lineage_share(
    pair_inputs: dict, weight_mode: str, target: str, count: int, viral_last: bool, blocked: bool
) -> None:
    from virosync.pipeline import host_signatures
    from virosync.pipeline.host_signatures import TaxonomyLabelLookup, build_host_signature_model
    from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy

    labels = {
        "EUK__HOST": "EUK|Alpha|Sharedclass|Sharedorder|NA|NA|Host one",
        "EUK__OTHER": "EUK|Beta|Sharedclass|Sharedorder|NA|NA|Other one",
        "PLASTID__HOST": "PLASTID|PLASTID|PLASTID|PLASTID|PLASTID|PLASTID|PLASTID",
    }
    marker = {
        "hmm_target": "COG0085",
        "top10_prefixes": "EUK__",
        "top10_targets": "EUK__HOST|p1",
        "top10_bitscores": "100",
        "top10_pidents": "90",
        "top10_evalues": "1e-40",
    }
    host_signatures.set_taxonomy_lookup(TaxonomyLabelLookup(labels))
    model = build_host_signature_model([marker], host_prefixes={"EUK__"}, weight_mode=weight_mode)
    genes = [pORF(f"g{i}", "s", 1050 + i * 100, 1100 + i * 100) for i in range(count)]
    taxonomy = {
        gene.id: GeneTaxonomy(
            gene.id,
            "s",
            gene.start,
            gene.end,
            top10_prefixes=[target.split("__")[0] + "__"],
            top10_targets=[f"{target}|p"],
            top10_pidents=[90.0],
            top10_bits=[100.0],
        )
        for gene in genes
    }
    if viral_last:
        taxonomy[genes[-1].id].top10_prefixes = ["NCLDV__"]
        taxonomy[genes[-1].id].top10_pidents = [30.0]
    pair_inputs.update(proteome_index={"s": genes}, taxonomy=taxonomy, host_model=model)
    pair_inputs["query"].all_porf_ids = list(taxonomy)

    pairs = eligible_joint_pairs(**pair_inputs)

    host_signatures.set_taxonomy_lookup(None)
    assert bool(pairs) is not blocked


@pytest.mark.parametrize("hull_end", [6000, 16000])
def test_native_joint_rescoring_preserves_cores_and_output_resume(tmp_path, hull_end: int) -> None:
    import json

    from test_eve_confidence_scoring import _supplementary_gene

    from virosync.orchestration._flows.single_genome.phase3 import _replace_accepted_joints
    from virosync.orchestration._flows.single_genome.phase_state import (
        phase2_state_from_document,
        phase2_state_to_document,
    )
    from virosync.pipeline.phase3.evidence_synthesizer import EvidenceSynthesizer, EvidenceSynthesizerConfig
    from virosync.pipeline.phase3.joint_seed_rescue import joint_boundary, joint_result_supported
    from virosync.pipeline.phase3.output_generator import (
        OutputGenerator,
        _evaluate_class_quality_gate,
        _qualifying_hallmark_genes,
        evaluate_v2_quality_gate,
    )

    children = (child("left", 0, 2000), child("right", 3000, hull_end))
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
    assert result.supplementary_viral_evidence.end == hull_end
    assert joint_result_supported(result, children)
    assert result.confirmed_seeds == restored.confirmed_seeds
    assert result.joint_seed_provenance == boundary.joint_seed_provenance
    provenance = json.loads(result.joint_seed_provenance)
    assert [c["candidate_id"] for c in provenance["children"]] == [boundary_candidate_id(c) for c in children]
    assert provenance["children"][0]["right_stop_reason"] == "support_gap"
    originals = [rejected(c) for c in children]
    assert result.confidence_tier == "LOW"
    assert _qualifying_hallmark_genes(result) == "."
    assert _evaluate_class_quality_gate(result).kept
    assert evaluate_v2_quality_gate(result).kept
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
    original_locus_ids = result.scoring_locus_ids
    result.scoring_locus_ids = {record["porf_id"]: "one_locus" for record in result.gene_taxonomy_records}
    assert result.supplementary_viral_evidence.eligible
    assert _evaluate_class_quality_gate(result).kept
    assert evaluate_v2_quality_gate(result).reason == "low_insufficient_competitive_viral_genes"
    selected = _replace_accepted_joints(originals, [result], [(restored, children)], ablation_id=AblationID.A0)
    assert all(actual is original for actual, original in zip(selected, originals, strict=True))
    assert tuple(seed for original in selected for seed in original.confirmed_seeds) == boundary.confirmed_seeds
    result.scoring_locus_ids = original_locus_ids
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


@pytest.mark.parametrize(
    "models,status,blocked",
    [
        (("GVOGm0003", "GVOGm0760"), "validated", True),
        (("GVOGm0760", "GVOGm0760"), "validated", False),
        (("GVOGm0003", "GVOGm0760"), "unvalidated", False),
    ],
)
def test_two_native_multi_marker_groups_stay_separate(
    pair_inputs: dict, models: tuple[str, str], status: str, blocked: bool
) -> None:
    from test_native_hallmark_evidence import hit

    from virosync.pipeline.hallmark_evidence import find_spatial_witnesses

    left, right = pair_inputs["results"]
    left.hallmark_evidence = [
        hit("a", models[0], 200, 400, status=status),
        hit("b", models[1], 500, 700, status=status),
    ]
    right.hallmark_evidence = [
        hit("c", models[0], 1600, 1800, status=status),
        hit("d", models[1], 1900, 2100, status=status),
    ]
    left.hallmark_witnesses = find_spatial_witnesses(left.hallmark_evidence)
    right.hallmark_witnesses = find_spatial_witnesses(right.hallmark_evidence)
    assert bool(eligible_joint_pairs(**pair_inputs)) is not blocked


def test_native_internal_repeat_veto_and_outer_pair_allowed() -> None:
    from test_terminal_repeats import _adjacent_repeat_sequence, _Marker, _planted_repeat

    from virosync.pipeline.phase3.joint_seed_rescue import assess_joint_repeat_scan

    sequence, _ = _adjacent_repeat_sequence()
    children = (child("a", 150, 400), child("b", 570, 820))
    markers = [_Marker("s", 220, 270, "a"), _Marker("s", 640, 690, "b")]
    assert not assess_joint_repeat_scan(children, sequence, markers, 5000, **repeat_context(len(sequence)))
    outer_sequence, _ = _planted_repeat()
    close_children = (child("a", 170, 240), child("b", 260, 350))
    close_markers = [_Marker("s", 190, 210, "a"), _Marker("s", 280, 310, "b")]
    assert assess_joint_repeat_scan(
        close_children, outer_sequence, close_markers, 100, **repeat_context(len(outer_sequence))
    )
    assert not assess_joint_repeat_scan(
        close_children, None, close_markers, 5000, **repeat_context(len(outer_sequence))
    )
    assert not assess_joint_repeat_scan(close_children, outer_sequence, [], 5000, **repeat_context(len(outer_sequence)))


def test_pairing_cap_is_unknown_and_blocks_rescue(monkeypatch) -> None:
    from test_terminal_repeats import _Marker, _planted_repeat

    from virosync.pipeline.phase2 import terminal_repeats
    from virosync.pipeline.phase3.joint_seed_rescue import assess_joint_repeat_scan

    sequence, _ = _planted_repeat()
    monkeypatch.setattr(terminal_repeats, "_MAX_SEED_PAIRINGS", 0)
    scan = terminal_repeats.find_terminal_inverted_repeats(
        sequence, scan_start=0, scan_end=len(sequence), marker_start=190, marker_end=210
    )
    assert scan.pairing_capped
    assert not assess_joint_repeat_scan(
        (child("a", 170, 240), child("b", 260, 350)),
        sequence,
        [_Marker("s", 190, 210, "a"), _Marker("s", 280, 310, "b")],
        5000,
        **repeat_context(len(sequence)),
    )


def test_native_joint_finalization_recomputes_class_composition_and_id(tmp_path, monkeypatch) -> None:
    import os
    import random
    import sys
    from pathlib import Path

    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome.phase3 import _prepare_joint_boundaries
    from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit

    monkeypatch.setenv("PATH", str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"])
    genome = tmp_path / "genome.fna"
    genome.write_text(">s\n" + "".join(random.Random(21).choices("ACGT", k=8000)) + "\n")
    config = PipelineConfig()
    config.phase1.extension_kb = 1
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
        all_porf_ids=["a", "b"],
        seed_gene_mappings={
            "a": SeedGeneMapping("a", "s", 0, 2000, flank_start_bp=0, flank_end_bp=5000),
            "b": SeedGeneMapping("b", "s", 3000, 7000, flank_start_bp=1000, flank_end_bp=8000),
        },
    )
    prepared = _prepare_joint_boundaries(
        [children],
        raw_genome_path=genome,
        masked_path=genome,
        validated_markers=markers,
        boundary_diamond_query=query,
        boundary_taxonomy_map={"a": GeneTaxonomy("a", "s", 500, 800), "b": GeneTaxonomy("b", "s", 4000, 4300)},
        proteome_index={"s": [pORF("a", "s", 500, 800), pORF("b", "s", 4000, 4300)]},
        config=config,
        metadata_families={},
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
        marker_annotation_index={},
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
            boundary_taxonomy_map=taxonomy,
            proteome_index={"s": genes},
            config=config,
            metadata_families={},
        )
        == []
    )


def test_internal_repeat_endpoint_in_seed_padding_still_separates() -> None:
    from test_terminal_repeats import _adjacent_repeat_sequence, _Marker

    from virosync.pipeline.phase3.joint_seed_rescue import assess_joint_repeat_scan

    sequence, _ = _adjacent_repeat_sequence()
    # First repeat ends at 420, inside the right seed's padded prefix.
    children = (child("a", 100, 350), child("b", 400, 900))
    assert not assess_joint_repeat_scan(
        children,
        sequence,
        [_Marker("s", 220, 270, "a"), _Marker("s", 640, 690, "b")],
        5000,
        **repeat_context(len(sequence)),
    )
    assert children[1].confirmed_seeds == (ConfirmedSeed("b", "s", 400, 900),)


@pytest.mark.parametrize("both_pass", [False, True])
@pytest.mark.parametrize("accepted", [False, True])
def test_competition_counts_only_fully_supported_native_results(both_pass: bool, accepted: bool) -> None:
    from test_eve_confidence_scoring import _supplementary_gene

    from virosync.orchestration._flows.single_genome.phase3 import _replace_accepted_joints
    from virosync.pipeline.phase3.evidence_synthesizer import calculate_eve_confidence
    from virosync.pipeline.phase3.joint_seed_rescue import joint_boundary

    children = (child("a", 0, 2000), child("b", 3000, 16000), child("c", 17000, 19000))
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
        16000,
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
        19000,
        confirmed_seeds=joints[1][0].confirmed_seeds,
        likely_family="NCLDV",
        region_classification="NCLDV",
        confidence_tier="LOW",
        gene_taxonomy_records=[
            _supplementary_gene("b", 3100, 3200),
            _supplementary_gene("c", 3500, 3600),
            _supplementary_gene("d", 17500, 17600, qcov=80 if both_pass else 10),
        ],
    )
    calculate_eve_confidence(first, 0)
    calculate_eve_confidence(second, 0)
    originals = [rejected(c) for c in children]
    selected = _replace_accepted_joints(
        originals,
        [first, second],
        joints,
        ablation_id=AblationID.A0,
        accepted_ids=frozenset(result.eve_id for result in originals) if accepted else frozenset(),
    )
    assert selected == originals if both_pass else selected == [originals[2], first]


def test_failed_required_repeat_scan_abstains(monkeypatch) -> None:
    from test_terminal_repeats import _Marker

    from virosync.pipeline.phase3 import joint_seed_rescue

    def fail_scan(*args: object, **kwargs: object) -> None:
        raise RuntimeError("required scan failed")

    monkeypatch.setattr(joint_seed_rescue, "find_terminal_inverted_repeats", fail_scan)
    assert not joint_seed_rescue.assess_joint_repeat_scan(
        (child("a", 100, 1000), child("b", 2000, 3000)),
        "A" * 4000,
        [_Marker("s", 400, 500, "a"), _Marker("s", 2400, 2500, "b")],
        5000,
        **repeat_context(4000),
    )


def test_competitive_supplementary_support_can_retain_short_joint() -> None:
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
    retained = _replace_accepted_joints(originals, [result], [(joint, children)], ablation_id=AblationID.A0)
    assert retained == [result]
    assert result.confirmed_seeds == joint.confirmed_seeds
    assert (result.start, result.end) == (0, 5000)


@pytest.mark.parametrize("model", [None, {}])
def test_unavailable_host_model_abstains(pair_inputs: dict, model: dict | None) -> None:
    from virosync.pipeline.host_signatures import HostSignatureModel

    pair_inputs["host_model"] = None if model is None else HostSignatureModel()
    assert eligible_joint_pairs(**pair_inputs) == []


@pytest.mark.parametrize("gap_base", ["N", "R"])
def test_ambiguous_intercore_dna_abstains(gap_base: str) -> None:
    from test_terminal_repeats import _Marker

    from virosync.pipeline.phase3.joint_seed_rescue import assess_joint_repeat_scan

    sequence = "A" * 1400 + gap_base + "A" * 1599
    assert not assess_joint_repeat_scan(
        (child("a", 100, 1000), child("b", 2000, 2900)),
        sequence,
        [_Marker("s", 400, 500, "a"), _Marker("s", 2400, 2500, "b")],
        5000,
        **repeat_context(len(sequence)),
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

    children = (child("a", 0, 250), child("b", 450, 15000))
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
    result = _supplementary_candidate(end=15000)
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


def test_joint_repeat_refinement_keeps_fixed_hull_with_outward_pair(tmp_path) -> None:
    import random

    from test_terminal_repeats import _planted_repeat

    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome.phase3 import _prepare_joint_boundaries
    from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit

    sequence, _ = _planted_repeat()
    rng = random.Random(381)
    sequence = "".join(rng.choices("ACGT", k=6000)) + sequence + "".join(rng.choices("ACGT", k=6000))
    genome = tmp_path / "joint.fna"
    genome.write_text(f">s\n{sequence}\n", encoding="utf-8")
    children = (child("a", 6150, 6200), child("b", 6300, 6350))
    query = GenomeDiamondQuery(
        all_porf_ids=["a", "b"],
        seed_gene_mappings={
            name: SeedGeneMapping(name, "s", start, end, flank_start_bp=0, flank_end_bp=len(sequence))
            for name, start, end in (("a", 6150, 6200), ("b", 6300, 6350))
        },
    )
    markers = [
        ValidatedMarkerHit(
            name,
            "s",
            start,
            end,
            "+",
            "GVOGm0760",
            100,
            1e-10,
            "validated",
            "NCLDV__",
            f"NCLDV__{name}",
            35,
            100,
            1,
            0,
            0,
            0,
            1,
        )
        for name, start, end in (("a", 6160, 6180), ("b", 6310, 6330))
    ]
    prepared = _prepare_joint_boundaries(
        [children],
        raw_genome_path=genome,
        masked_path=genome,
        validated_markers=markers,
        boundary_diamond_query=query,
        boundary_taxonomy_map={"a": GeneTaxonomy("a", "s", 6160, 6180), "b": GeneTaxonomy("b", "s", 6310, 6330)},
        proteome_index={"s": [pORF("a", "s", 6160, 6180), pORF("b", "s", 6310, 6330)]},
        config=PipelineConfig(),
        metadata_families={},
    )
    assert len(prepared) == 1
    joint = prepared[0][0]
    assert (joint.start, joint.end, joint.tir_status, joint.tir_boundary_override) == (6150, 6350, "ambiguous", False)
    assert joint.confirmed_seeds == children[0].confirmed_seeds + children[1].confirmed_seeds


@pytest.mark.parametrize(
    ("deep_viral", "viral_marker", "blocked"),
    [(False, False, True), (True, False, False), (False, True, False)],
    ids=["no_veto", "viral_hit_within_depth", "validated_viral_marker"],
)
def test_gap_barrier_yields_to_a_deep_viral_hit_or_a_validated_viral_marker(
    pair_inputs: dict, monkeypatch, deep_viral: bool, viral_marker: bool, blocked: bool
) -> None:
    from types import SimpleNamespace

    from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy
    from virosync.pipeline.phase3 import joint_seed_rescue

    genes = [pORF(f"g{i}", "s", 1050 + i * 100, 1100 + i * 100) for i in range(3)]
    taxonomy = {g.id: GeneTaxonomy(g.id, "s", g.start, g.end) for g in genes}
    taxonomy["g2"].viral_hit_within_depth = deep_viral
    monkeypatch.setattr(joint_seed_rescue, "host_lineage_share", lambda gene, model: 0.5)
    pair_inputs.update(
        proteome_index={"s": genes},
        taxonomy=taxonomy,
        host_model=SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0}),
        viral_marker_genes=frozenset({"g2"}) if viral_marker else frozenset(),
    )
    pair_inputs["query"].all_porf_ids = list(taxonomy)

    assert bool(eligible_joint_pairs(**pair_inputs)) is not blocked


def test_phase3_flow_gives_the_validated_viral_markers_to_the_joint_rescue_barrier(tmp_path, monkeypatch) -> None:
    import logging
    from types import SimpleNamespace

    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome import phase3

    class StopAfterEligibility(Exception):
        pass

    seen: dict[str, object] = {}

    def capture(*_args, **kwargs):
        seen.update(kwargs)
        raise StopAfterEligibility

    def marker(porf: str, status: str, target: str) -> SimpleNamespace:
        return SimpleNamespace(query_porf=porf, validation_status=status, hmm_target=target, hmm_score=50)

    monkeypatch.setattr(phase3, "eligible_joint_pairs", capture)
    monkeypatch.setattr(phase3, "_run_interproscan", lambda **_kwargs: {})
    monkeypatch.setattr(phase3, "_precompute_tmvec", lambda **_kwargs: (None, None))
    monkeypatch.setattr(phase3, "_load_or_classify_jelly_roll", lambda **_kwargs: {})
    monkeypatch.setattr(phase3, "_verify_native_boundaries", lambda *_args, **_kwargs: [])
    markers = [marker("g1|aa1", "validated", "GVOGm0003"), marker("g2", "unvalidated", "GVOGm0003")]

    with pytest.raises(StopAfterEligibility):
        phase3._run_phase3_subflow(
            masked_path=tmp_path / "masked.fna",
            proteome_path=tmp_path / "proteome.faa",
            validated_markers=markers,
            hallmark_evidence=[],
            hallmark_witnesses=[],
            host_signatures=set(),
            host_signature_model_payload=None,
            refined_boundaries=[],
            boundary_taxonomy_map={},
            boundary_diamond_query=None,
            proteome_index={},
            merged_seeds=[],
            output_dir=tmp_path,
            config=PipelineConfig(),
            validated_hits_tsv=tmp_path / "validated_hits.tsv",
            logger=logging.getLogger(__name__),
            raw_genome_path=tmp_path / "genome.fna",
            metadata_families={},
            marker_annotation_index={},
        )

    assert seen["viral_marker_genes"] == frozenset({"g1"})


@pytest.mark.parametrize("mode,expected", [("relaxed", "NCLDV"), ("default", "UNKNOWN")])
def test_joint_finalization_uses_configured_a2_assembly_floor(tmp_path, monkeypatch, mode, expected):
    import os
    import random
    import sys
    from pathlib import Path

    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome.phase3 import _prepare_joint_boundaries
    from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit

    monkeypatch.setenv("PATH", str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"])
    genome = tmp_path / "genome.fna"
    genome.write_text(">s\n" + "".join(random.Random(21).choices("ACGT", k=8000)) + "\n")
    config = PipelineConfig().with_overrides(ablation_id=AblationID.A2, assembly_mode=mode, extension_kb=1)
    children = (child("a", 0, 2000), child("b", 3000, 7000))
    left_marker = ValidatedMarkerHit(
        "a",
        "s",
        500,
        800,
        "+",
        "unclassified_model",
        100,
        1e-10,
        "validated",
        "NCLDV__",
        "NCLDV__a",
        35,
        100,
        1,
        0,
        0,
        0,
        1,
    )
    right_marker = replace(left_marker, query_porf="b", start=4000, end=4300)
    bypass_marker = replace(
        left_marker,
        query_porf="c",
        start=1000,
        end=1300,
        hmm_target="GVOGm0001",
        hmm_score=80,
        validation_status="supported",
        top10_prefixes="GVMAG__",
    )
    markers = [left_marker, right_marker, bypass_marker]
    query = GenomeDiamondQuery(
        all_porf_ids=["a", "b", "c"],
        seed_gene_mappings={
            "a": SeedGeneMapping("a", "s", 0, 2000, flank_start_bp=0, flank_end_bp=5000),
            "b": SeedGeneMapping("b", "s", 3000, 7000, flank_start_bp=1000, flank_end_bp=8000),
        },
    )
    prepared = _prepare_joint_boundaries(
        [children],
        raw_genome_path=genome,
        masked_path=genome,
        validated_markers=markers,
        boundary_diamond_query=query,
        boundary_taxonomy_map={m.query_porf: GeneTaxonomy(m.query_porf, "s", m.start, m.end) for m in markers},
        proteome_index={"s": [pORF(m.query_porf, "s", m.start, m.end) for m in markers]},
        config=config,
        metadata_families={},
    )
    assert len(prepared) == 1
    boundary, retained_children = prepared[0]
    assert boundary.predicted_family == expected
    assert (boundary.start, boundary.end) == (0, 7000)
    assert boundary.confirmed_seeds == children[0].confirmed_seeds + children[1].confirmed_seeds
    assert retained_children == children
    assert bypass_marker.validation_status == "supported"
    assert not hasattr(bypass_marker, "tier1_bypassed")


@pytest.mark.parametrize("target,blocked", [("EUK__host", True), ("BAC__foreign", False), ("EUK-pEVE__host", False)])
def test_gap_uses_real_token_distribution_and_preserves_children(pair_inputs, monkeypatch, target, blocked):
    """Actual scorer and rescue gate agree without editing confirmed seed spans."""
    from virosync.pipeline import host_signatures

    monkeypatch.setattr(
        host_signatures,
        "_TAXONOMY_LOOKUP",
        {
            "EUK__host": "EUK|Amoebozoa",
            "BAC__foreign": "BAC|Other",
            "EUK-pEVE__host": "EUK|Amoebozoa",
        },
    )
    genes = [pORF(f"g{i}", "s", 1050 + i * 100, 1100 + i * 100) for i in range(3)]
    taxonomy = {
        g.id: GeneTaxonomy(
            g.id,
            "s",
            g.start,
            g.end,
            top10_targets=[target],
            top10_prefixes=[target.split("__")[0] + "__"],
            top10_pidents=[90.0],
        )
        for g in genes
    }
    pair_inputs.update(
        proteome_index={"s": genes},
        taxonomy=taxonomy,
        host_model=HostSignatureModel(token_weights={"eukaryota": 10, "amoebozoa": 10}),
    )
    pair_inputs["query"].all_porf_ids = list(taxonomy)
    before = [boundary.confirmed_seeds for boundary in pair_inputs["boundaries"]]
    assert bool(eligible_joint_pairs(**pair_inputs)) is not blocked
    assert [boundary.confirmed_seeds for boundary in pair_inputs["boundaries"]] == before
