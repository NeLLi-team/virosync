"""Native retained-hallmark scoring and protected spatial rescue contracts."""

import json
from dataclasses import replace

import pytest

from virosync.ablation import AblationID
from virosync.pipeline.hallmark_evidence import (
    HallmarkEvidence,
    add_rescue_seeds,
    collect_hallmark_evidence,
    find_spatial_witnesses,
    summarize_hallmarks,
)
from virosync.pipeline.phase1.hhg_seeding import HMMHit
from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
from virosync.pipeline.phase2.seed_extension import (
    accept_seed_partition,
    extend_confirmed_seed,
    merge_reachable_boundaries,
)
from virosync.pipeline.phase3.acceptance_selection import select_phase3_acceptance
from virosync.pipeline.phase3.eve_ani_clustering import unsupported_eve_ids
from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    EvidenceSynthesizerConfig,
    VerificationResult,
    calculate_eve_confidence,
)
from virosync.pipeline.phase3.output_generator import OutputGenerator


def hit(gene, model, start=100, end=1000, *, evalue=0.001, status="unvalidated", identity=0):
    return HallmarkEvidence(
        query_porf=f"{gene}|aa1-20",
        gene=gene,
        scaffold="s",
        start=start,
        end=end,
        model=model,
        evalue=evalue,
        score=70.123456789,
        domain_score=60.123456789,
        query_start=1,
        query_end=20,
        validation_status=status,
        viral_top10_max_pident=identity,
        source="retained_hmm",
    )


@pytest.mark.parametrize("evalue,weighted", [(1e-5, 1.0), (1.000001e-5, 0.25), (0.01, 0.25), (0.01000001, 0.0)])
def test_sequence_evalue_thresholds_are_inclusive(evalue, weighted):
    assert summarize_hallmarks([hit("g", "GVOGm0760", evalue=evalue, identity=30)]).weighted_core_count == weighted


def test_exact_retained_join_never_restores_discarded_models_and_frameshift_diagnostics_stay_null():
    marker = ValidatedMarkerHit(
        query_porf="g|aa1-20",
        scaffold="s",
        start=100,
        end=1000,
        strand="+",
        hmm_target="GVOGm0003",
        hmm_score=70.123456789,
        hmm_evalue=0.001000000001,
        validation_status="unvalidated",
        top10_prefixes="EUK__",
        top10_pidents="90",
        best_hit_target="EUK__host",
        best_hit_pident=90,
        best_hit_bits=100,
        has_ncldv=0,
        has_mirus=0,
        has_plv=0,
        has_vp=0,
        has_viral=0,
    )
    raw = [
        HMMHit("g", "GVOGm0003", marker.hmm_score, marker.hmm_evalue, 60.123456789, 1, 20),
        HMMHit("g", "GVOGm0760", 60, 1e-8, 59, 1, 20),
    ]
    evidence = collect_hallmark_evidence([marker], raw)
    assert len(evidence) == 1
    assert evidence[0].evalue == 0.001000000001
    assert evidence[0].domain_score == 60.123456789
    assert evidence[0].model == "GVOGm0003"
    assert evidence[0].viral_top10_max_pident == 0
    frameshift = collect_hallmark_evidence([marker], [], frameshift_ids={"g"})
    assert frameshift[0].source == "frameshift"
    assert frameshift[0].domain_score is frameshift[0].query_start is frameshift[0].query_end is None
    assert frameshift[0].evalue == marker.hmm_evalue


def test_graduated_pair_maximizes_strong_valid_assignment():
    hits = [hit("a", "GVOGm0003", evalue=1e-8), hit("b", "GVOGm0003"), hit("a", "GVOGm0760", identity=30)]
    summary = summarize_hallmarks(hits)
    assert summary.core_count == 2
    assert summary.weighted_core_count == 1.0
    assert summary.capsid_fraction == 1.0
    assert summary.marker_score == pytest.approx(0.6)
    assert summary.floor_tier == "HIGH"


def test_spatial_matching_uses_full_gene_span_and_preserves_contained_triple():
    hits = [
        hit("a", "GVOGm0760", 100, 1000),
        hit("b", "GVOGm0054", 5000, 6000),
        hit("c", "GVOGm0013", 14000, 15100),
        hit("d", "GVOGm0461", 14900, 16000),
        hit("copy", "GVOGm0760", 0, 500),
    ]
    witnesses = find_spatial_witnesses(hits)
    assert any(w.start == 100 and w.end == 15100 for w in witnesses)
    assert not find_spatial_witnesses([hits[0], hits[1], replace(hits[2], end=15101)])
    assert not find_spatial_witnesses([hits[0], replace(hits[1], gene="a"), hits[2]])
    assert not find_spatial_witnesses([hits[0], replace(hits[1], model="GVOGm0760"), hits[2]])


def test_capsid_partner_gap_and_nonspatial_floor_are_distinct():
    capsid = hit("a", "GVOGm0003", 100, 1000)
    partner = hit("b", "GVOGm0760", 6000, 7000, identity=30)
    assert any(w.tier == "HIGH" for w in find_spatial_witnesses([capsid, partner]))
    far = replace(partner, start=6001)
    assert all(w.tier == "MEDIUM" for w in find_spatial_witnesses([capsid, far]))
    assert summarize_hallmarks([capsid, far]).floor_tier == "HIGH"


@pytest.mark.parametrize("models", [["GVOGm0760"], ["GVOGm0760", "GVOGm0054", "GVOGm0013"]])
def test_spatial_isolated_groups_skip_impossible_matching(monkeypatch, models):
    from virosync.pipeline import hallmark_evidence

    def impossible_assignment(_hits):
        pytest.fail("A window without three distinct families cannot produce a quorum")

    monkeypatch.setattr(hallmark_evidence, "_assignment", impossible_assignment)
    evidence = [hit(f"g{i}", models[i % len(models)], i * 20000, i * 20000 + 300) for i in range(40)]
    assert find_spatial_witnesses(evidence) == []


def test_spatial_partner_bounds_keep_long_overlaps_and_exact_edges():
    capsid = hit("capsid", "GVOGm0003", 100000, 101000)
    partners = [
        hit("long_left", "GVOGm0760", 0, 100500),
        hit("left_edge", "GVOGm0054", 94000, 95000),
        hit("right_edge", "GVOGm0013", 106000, 107000),
        hit("outside", "GVOGm0461", 106001, 107000),
    ]
    witnesses = find_spatial_witnesses([capsid, *partners])
    pairs = [w for w in witnesses if w.rule == "capsid_partner_5kb"]
    assert {w.evidence_ids for w in pairs} == {
        tuple(h.evidence_id for h in sorted([capsid, partner], key=lambda h: h.gene)) for partner in partners[:3]
    }
    assert any(w.start == 0 and w.end == 101000 for w in pairs)
    assert not any(w.rule == "three_core_15kb" and w.start == 0 for w in witnesses)


def test_partial_overlap_adds_distinct_seed_and_containment_does_not():
    evidence = [hit("a", "GVOGm0003", 100, 1000)]
    witnesses = find_spatial_witnesses(evidence)
    old = MergedSeed("s", 0, 500, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", 0, 500),))
    seeds = add_rescue_seeds([old], witnesses)
    assert seeds[0] == old
    assert len(seeds) == 2
    assert seeds[1].start == 100 and seeds[1].end == 1000
    assert seeds[1].seed_id != old.seed_id
    containing = replace(old, end=1500, confirmed_seeds=(ConfirmedSeed("original", "s", 0, 1500),))
    assert add_rescue_seeds([containing], witnesses) == [containing]


def test_partial_rescue_does_not_merge_across_host_stop_or_split_its_core():
    old = MergedSeed("s", 0, 1500, seed_id="original", confirmed_seeds=(ConfirmedSeed("original", "s", 0, 1500),))
    evidence = [hit("a", "GVOGm0003", 1000, 2000)]
    seeds = add_rescue_seeds([old], find_spatial_witnesses(evidence))
    boundaries = [extend_confirmed_seed(s, [], {}, (0, 5000), None, 0.35) for s in seeds]
    boundaries[0].right_stop, boundaries[0].right_stop_reason = 1500, "host_stop"
    assert len(merge_reachable_boundaries(boundaries)) == 2
    assert boundaries[1].confirmed_seeds == seeds[1].confirmed_seeds
    assert accept_seed_partition(boundaries[1], [replace(boundaries[1], end=1600)]) is None


def test_nested_new_hulls_keep_all_witnesses_but_one_confirmed_seed():
    evidence = [hit("a", "GVOGm0003", 100, 1000), hit("b", "GVOGm0760", 1100, 2000)]
    witnesses = find_spatial_witnesses(evidence)
    assert len(witnesses) == 2
    (seed,) = add_rescue_seeds([], witnesses)
    assert (seed.start, seed.end) == (100, 2000)


def test_raw_rescue_survives_normal_extension_scoring_acceptance_and_ani(tmp_path):
    evidence = [hit("a", "GVOGm0003")]
    witnesses = find_spatial_witnesses(evidence)
    (seed,) = add_rescue_seeds([], witnesses)
    boundary = extend_confirmed_seed(seed, [], {}, (0, 5000), None, 0.35)
    synthesizer = EvidenceSynthesizer(
        config=EvidenceSynthesizerConfig(use_boltz=False, use_phylogenetic_validation=False, use_tmvec_database=False),
        work_dir=tmp_path,
    )
    result = synthesizer.verify_eve(boundary, [], hallmark_evidence=evidence, hallmark_witnesses=witnesses)
    assert (result.start, result.end) == (100, 1000)
    assert result.confirmed_seeds == seed.confirmed_seeds
    assert result.confidence_tier == "MEDIUM"
    assert result.taxonomy_class == "UNKNOWN"
    assert result.hallmark_evidence[0].validation_status == "unvalidated"
    selection = select_phase3_acceptance([result], AblationID.A0)
    assert selection.canonical_results == (result,)
    assert selection.normal_gate_decisions[0].reason == "hallmark_spatial_rescue_pass"
    assert unsupported_eve_ids([result]) == set()
    export_path = OutputGenerator(tmp_path / "export").write_evidence_profiles([result])
    exported = json.loads(export_path.read_text())[result.eve_id]
    exported_by_id = {record["evidence_id"]: record for record in exported["hallmark_evidence"]}
    member_ids = {member for witness in exported["hallmark_witnesses"] for member in witness["evidence_ids"]}
    assert member_ids == set(exported_by_id)
    assert set(exported["hallmark_assessment"]["summary"]["qualifying_evidence_ids"]) <= set(exported_by_id)
    assert exported_by_id[evidence[0].evidence_id]["evalue"] == 0.001
    assert exported_by_id[evidence[0].evidence_id]["validation_status"] == "unvalidated"
    result.end -= 1
    assert not select_phase3_acceptance([result], AblationID.A0).canonical_results
    assert unsupported_eve_ids([result]) == {result.eve_id}


def test_json_exports_catalogue_roles_without_turning_portal_into_capsid(tmp_path):
    evidence = [hit("mcp", "Mirus_MCP"), hit("portal", "Mirus_Portal")]
    result = VerificationResult("eve", "s", 0, 2000, hallmark_evidence=evidence)
    export_path = OutputGenerator(tmp_path).write_evidence_profiles([result])
    exported = json.loads(export_path.read_text())["eve"]["hallmark_evidence"]
    by_model = {record["model"]: record for record in exported}
    mcp_roles = by_model["Mirus_MCP"]["model_roles"]
    portal_roles = by_model["Mirus_Portal"]["model_roles"]
    assert mcp_roles == [
        {
            "panel": "Mirusviricota",
            "marker": "MCP",
            "module": "capsid_packaging",
            "role": "core_marker",
            "capsid": True,
            "structural_anchor": True,
            "excluded_from_scoring": False,
        }
    ]
    assert portal_roles == [
        {
            "panel": "Mirusviricota",
            "marker": "Portal",
            "module": "capsid_packaging",
            "role": "core_marker",
            "capsid": False,
            "structural_anchor": True,
            "excluded_from_scoring": False,
        }
    ]
    assert {record["validation_status"] for record in exported} == {"unvalidated"}
    assert not result.hallmark_assessment


@pytest.mark.parametrize(
    "evalues,expected",
    [
        ([0.001], 0.01875),
        ([1e-8], 0.075),
        ([0.001, 0.001], 0.1),
        ([1e-8, 0.001], 0.15),
        ([1e-8, 1e-8], 0.2),
        ([0.001, 0.001, 0.001], 0.325),
    ],
)
def test_fixed_nonspatial_scores_match_selected_prototype(evalues, expected):
    models = ["GVOGm0760", "GVOGm0054", "GVOGm0013"]
    evidence = [hit(f"g{i}", models[i], evalue=value, identity=30) for i, value in enumerate(evalues)]
    result = VerificationResult("eve", "s", 0, 2000, gene_count=5, hallmark_evidence=evidence)
    assert calculate_eve_confidence(result, 0) == pytest.approx(expected)
    result.hallmark_genes = ["GVOGm0003", "GVOGm0760", "GVOGm0054", "GVOGm0013"]
    result.ncldv_completeness_ratio = 1.0
    result.jelly_roll_confidence_bonus = 0.5
    assert calculate_eve_confidence(result, 0) == pytest.approx(expected)


@pytest.mark.parametrize("mutation", ["wrong_panel", "wrong_rule", "outside", "duplicate_family"])
def test_invalid_witness_does_not_gain_an_acceptance_exception(mutation):
    evidence = [hit("a", "GVOGm0760", 100, 1000), hit("b", "GVOGm0054", 1100, 2000), hit("c", "GVOGm0013", 2100, 3000)]
    (witness,) = find_spatial_witnesses(evidence)
    if mutation == "wrong_panel":
        witness = replace(witness, panel="CRESS")
    elif mutation == "wrong_rule":
        witness = replace(witness, rule="capsid_partner_5kb", tier="HIGH")
    elif mutation == "outside":
        evidence[2] = replace(evidence[2], end=4000)
    else:
        evidence[2] = replace(evidence[2], model="GVOGm0760")
    result = VerificationResult(
        "eve", "s", 0, 3500, gene_count=5, hallmark_evidence=evidence, hallmark_witnesses=[witness]
    )
    result.final_confidence = calculate_eve_confidence(result, 0)
    result.confidence_tier = "MEDIUM"
    assert not result.hallmark_assessment.witnesses
    assert not select_phase3_acceptance([result], AblationID.A0).canonical_results


def test_context_scans_run_once_before_scoring_and_exclude_partial_genes(tmp_path, monkeypatch):
    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome import phase3
    from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary, boundary_candidate_id
    from virosync.pipeline.phase3.gene_functions import GeneFunction, GeneFunctionResult

    boundary = RefinedBoundary("s", 100, 3000, confidence=0)
    identity = boundary_candidate_id(boundary)
    taxonomy = {
        identity: (
            [
                {"porf_id": "inside", "start": 100, "end": 300, "is_flanking": False},
                {"porf_id": "partial", "start": 2999, "end": 3100, "is_flanking": False},
            ],
            {},
        )
    }
    integration = {
        identity: {
            "integration_gene_hits": [
                {"protein_id": "inside", "profile": "rve", "evalue": 1e-8, "location": "interior"},
                {"protein_id": "flank", "profile": "rve", "evalue": 1e-8, "location": "downstream"},
            ]
        }
    }
    monkeypatch.setattr(phase3, "_scan_candidate_integration", lambda *args, **kwargs: integration)
    monkeypatch.setattr(phase3, "find_gene_function_resources", lambda path: object())
    calls = []

    def label(genes, *args, **kwargs):
        calls.append(genes)
        return GeneFunctionResult([GeneFunction(identity, "inside", "ATPase", "capsid", "Pfam:Pox_A32", 1e-8)], [])

    monkeypatch.setattr(phase3, "label_core_genes", label)
    config = PipelineConfig()
    config.databases.hmm_database = tmp_path / "combined.hmm"
    annotations = phase3._prepare_candidate_annotations(
        [boundary],
        taxonomy,
        [],
        proteome_path=tmp_path / "p.faa",
        validated_markers=[],
        validated_marker_index={},
        output_dir=tmp_path,
        config=config,
    )
    assert calls == [[(identity, "inside")]]
    assert {h.gene for h in annotations[identity]["hallmark_support"]} == {"inside"}
    assert len(annotations[identity]["gene_function_labels"]) == 1


def test_batched_verification_routes_contained_evidence_and_witnesses_in_input_order(tmp_path, monkeypatch):
    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome import phase3
    from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary, boundary_candidate_id

    class CountedRecords(list):
        iterations = 0

        def __iter__(self):
            self.iterations += 1
            return super().__iter__()

    inside_later = hit("later", "GVOGm0003", 700, 900)
    inside_earlier = hit("earlier", "GVOGm0760", 200, 400)
    partial = hit("partial", "GVOGm0003", 0, 300)
    remote = replace(inside_later, gene="remote", scaffold="elsewhere")
    evidence = CountedRecords([inside_later, remote, partial, inside_earlier])
    witnesses = CountedRecords(find_spatial_witnesses(evidence))
    evidence.iterations = 0
    boundaries = [RefinedBoundary("s", 100, 1000), RefinedBoundary("s", 2000, 3000)]
    captured = {}

    def capture(_task, **kwargs):
        captured.update(kwargs)
        return []

    monkeypatch.setattr(phase3, "call_task", capture)
    phase3._run_verification(
        boundaries=boundaries,
        hallmark_hits_map={},
        hallmark_evidence=evidence,
        hallmark_witnesses=witnesses,
        candidate_annotations={},
        gene_taxonomy_map={},
        interproscan_map={},
        jelly_roll_map={},
        masked_path=tmp_path / "g.fna",
        proteome_path=tmp_path / "p.faa",
        output_dir=tmp_path,
        host_signatures=set(),
        host_signature_model_payload=None,
        precomputed_tmvec=None,
        tmvec_device="cpu",
        config=PipelineConfig(),
    )
    first, second = map(boundary_candidate_id, boundaries)
    assert captured["hallmark_evidence_map"] == {first: [inside_later, inside_earlier], second: []}
    assert captured["hallmark_witnesses_map"][first] == [
        w for w in witnesses if w.scaffold == "s" and 100 <= w.start and w.end <= 1000
    ]
    assert captured["hallmark_witnesses_map"][second] == []
    assert evidence.iterations == 1
    assert witnesses.iterations == 2  # Once for indexing, once for the expected-value comparison.


def test_candidate_annotations_reuse_marker_index_and_scan_evidence_once(tmp_path, monkeypatch):
    from types import SimpleNamespace

    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome import phase3
    from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary, boundary_candidate_id
    from virosync.pipeline.phase3.gene_functions import GeneFunctionResult

    class CountedRecords(list):
        iterations = 0

        def __iter__(self):
            self.iterations += 1
            return super().__iter__()

    evidence = CountedRecords([hit("raw", "GVOGm0003", 300, 500), hit("clipped", "GVOGm0003", 0, 200)])
    markers = CountedRecords(
        [
            SimpleNamespace(query_porf="frameshift|aa1-20", scaffold="s", start=600, end=900),
            SimpleNamespace(query_porf="outside|aa1-20", scaffold="other", start=600, end=900),
        ]
    )
    marker_index = phase3._build_scaffold_start_index(markers)
    boundaries = [RefinedBoundary("s", 100, 1000), RefinedBoundary("s", 2000, 3000)]
    monkeypatch.setattr(
        phase3,
        "_scan_candidate_integration",
        lambda *args, **kwargs: {boundary_candidate_id(b): {"integration_gene_hits": []} for b in boundaries},
    )
    monkeypatch.setattr(phase3, "find_gene_function_resources", lambda _: object())
    labelled = []

    def label(genes, *args, **kwargs):
        labelled.extend(genes)
        return GeneFunctionResult([], [])

    monkeypatch.setattr(phase3, "label_core_genes", label)
    config = PipelineConfig()
    config.databases.hmm_database = tmp_path / "combined.hmm"
    phase3._prepare_candidate_annotations(
        boundaries,
        {},
        evidence,
        proteome_path=tmp_path / "p.faa",
        validated_markers=markers,
        validated_marker_index=marker_index,
        output_dir=tmp_path,
        config=config,
    )
    identity = boundary_candidate_id(boundaries[0])
    assert labelled == [(identity, "frameshift"), (identity, "raw")]
    assert evidence.iterations == markers.iterations == 1
