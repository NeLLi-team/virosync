"""Native retained-hallmark scoring and protected spatial rescue contracts."""

import csv
import json
import logging
from dataclasses import replace
from pathlib import Path

import pytest

from virosync.ablation import AblationID
from virosync.pipeline.hallmark_evidence import (
    CATALOGUE_PATH,
    HallmarkEvidence,
    HallmarkSupport,
    add_rescue_seeds,
    assess_hallmarks,
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


def hit(gene, model, start=100, end=1000, *, evalue=0.001, status="validated", identity=30):
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


@pytest.mark.parametrize("model", ["VP_MCP_6", "VS000792", "VS000652", "VS000745", "VS000783", "VS000788", "VS000791"])
@pytest.mark.parametrize("evalue", [0.0021, 1e-8])
@pytest.mark.parametrize("status,identity", [("unvalidated", 0), ("supported", 20), ("unvalidated", 30)])
def test_unvalidated_capsid_cannot_supply_core_floor_or_rescue_seed(model, evalue, status, identity):
    evidence = [hit("capsid", model, evalue=evalue, status=status, identity=identity)]
    original = evidence[0]

    summary = summarize_hallmarks(evidence)
    witnesses = find_spatial_witnesses(evidence)

    assert summary.core_count == 0
    assert summary.qualifying_evidence_ids == ()
    assert summary.floor_tier == ""
    assert witnesses == []
    assert add_rescue_seeds([], witnesses) == []
    assert evidence == [original]


def test_stale_unvalidated_witness_cannot_rescue_or_appear_in_hallmark_column(tmp_path):
    original = hit("capsid", "VP_MCP_6", evalue=0.0021)
    witnesses = find_spatial_witnesses([original])
    (seed,) = add_rescue_seeds([], witnesses)
    boundary = extend_confirmed_seed(seed, [], {}, (0, 5000), None, 0.35)
    evidence = [replace(original, validation_status="unvalidated", viral_top10_max_pident=0)]
    synthesizer = EvidenceSynthesizer(
        config=EvidenceSynthesizerConfig(use_boltz=False, use_phylogenetic_validation=False, use_tmvec_database=False),
        work_dir=tmp_path,
    )

    result = synthesizer.verify_eve(boundary, [], hallmark_evidence=evidence, hallmark_witnesses=witnesses)
    selection = select_phase3_acceptance([result], AblationID.A0)
    generator = OutputGenerator(tmp_path / "export")
    detailed_path = generator.write_predictions_detailed_tsv([result])
    profiles_path = generator.write_evidence_profiles([result])

    assert result.hallmark_assessment.summary.core_count == 0
    assert result.hallmark_assessment.summary.floor_tier == ""
    assert result.hallmark_assessment.witnesses == ()
    assert selection.canonical_results == ()
    assert selection.normal_gate_decisions[0].reason != "hallmark_spatial_rescue_pass"
    with detailed_path.open() as stream:
        assert next(csv.DictReader(stream, delimiter="\t"))["qualifying_hallmark_genes"] == "."
    profile = json.loads(profiles_path.read_text())[result.eve_id]
    assert profile["hallmark_evidence"][0]["validation_status"] == "unvalidated"


def test_validated_frameshift_twin_supplies_core_without_unvalidated_ordinary_copy():
    raw = hit("ordinary", "GVOGm0003", evalue=1e-20, status="unvalidated", identity=0)
    rescued = replace(hit("rescued", "GVOGm0003", evalue=1e-10), source="frameshift")

    summary = summarize_hallmarks([raw, rescued], scoring_locus_ids={"rescued": "ordinary"})
    (witness,) = find_spatial_witnesses([raw, rescued])

    assert summary.core_count == summary.strong_count == 1
    assert summary.qualifying_evidence_ids == (rescued.evidence_id,)
    assert witness.evidence_ids == (rescued.evidence_id,)


@pytest.mark.parametrize("extended", [False, True])
@pytest.mark.parametrize("published_class", ["NCLDV", "MIRUS"])
def test_qualifying_hallmark_genes_export_uses_scored_panel(tmp_path, extended, published_class):
    result = VerificationResult(
        eve_id="eve",
        scaffold="s",
        start=0,
        end=10000,
        taxonomy_class=published_class,
        confidence_tier="LOW",
        final_confidence=0.2,
        hallmark_evidence=[
            hit("capsid", "GVOGm0003", evalue=1e-8),
            hit("atpase", "GVOGm0760", evalue=1e-8, identity=30),
            hit("atpase_copy", "GVOGm0760", evalue=1e-8, identity=30),
            hit("context", "VS000151", evalue=1e-8, identity=30),
            hit("outside", "GVOGm0054", start=9990, end=11000, identity=30),
        ],
    )
    result.hallmark_assessment = assess_hallmarks(result)
    assessment = result.hallmark_assessment
    generator = OutputGenerator(tmp_path, extended_output=extended)

    canonical_path = generator.write_predictions_tsv([result])
    detailed_path = generator.write_predictions_detailed_tsv([result])

    with canonical_path.open() as stream:
        canonical = next(csv.DictReader(stream, delimiter="\t"))
    with detailed_path.open() as stream:
        detailed = next(csv.DictReader(stream, delimiter="\t"))
    assert canonical["qualifying_hallmark_genes"] == "A32|MCP"
    assert detailed["qualifying_hallmark_genes"] == "A32|MCP"
    assert canonical["effective_eve_class"] == published_class
    assert result.hallmark_assessment is assessment
    assert result.final_confidence == 0.2


@pytest.mark.parametrize("assessed", [False, True])
def test_qualifying_hallmark_genes_export_does_not_infer_from_raw_hits(tmp_path, assessed):
    result = VerificationResult(
        eve_id="eve",
        scaffold="s",
        start=0,
        end=10000,
        hallmark_genes=["GVOGm0760", "VS000151"],
        hallmark_evidence=[
            hit("unvalidated_atpase", "GVOGm0760", evalue=1e-8, status="unvalidated", identity=0),
            hit("context", "VS000151", evalue=1e-8, identity=30),
        ],
    )
    if assessed:
        result.hallmark_assessment = assess_hallmarks(result)
    generator = OutputGenerator(tmp_path)

    canonical_path = generator.write_predictions_tsv([result])
    detailed_path = generator.write_predictions_detailed_tsv([result])

    with canonical_path.open() as stream:
        canonical = next(csv.DictReader(stream, delimiter="\t"))
    with detailed_path.open() as stream:
        detailed = next(csv.DictReader(stream, delimiter="\t"))
    assert canonical["qualifying_hallmark_genes"] == "."
    assert detailed["qualifying_hallmark_genes"] == "."


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
    assert (
        accept_seed_partition(
            boundaries[1], [replace(boundaries[1], end=1600)], allowed_interval=(boundaries[1].start, boundaries[1].end)
        )
        is None
    )


def test_nested_new_hulls_keep_all_witnesses_but_one_confirmed_seed():
    evidence = [hit("a", "GVOGm0003", 100, 1000), hit("b", "GVOGm0760", 1100, 2000)]
    witnesses = find_spatial_witnesses(evidence)
    assert len(witnesses) == 2
    (seed,) = add_rescue_seeds([], witnesses)
    assert (seed.start, seed.end) == (100, 2000)


@pytest.mark.parametrize("status,identity", [("validated", 30), ("validated_novel", 0)])
def test_validated_rescue_survives_normal_extension_scoring_acceptance_and_ani(tmp_path, status, identity):
    evidence = [hit("a", "GVOGm0003", status=status, identity=identity)]
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
    assert result.hallmark_evidence[0].validation_status == status
    assert result.hallmark_assessment.summary.weighted_core_count == 0.25
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
    assert exported_by_id[evidence[0].evidence_id]["validation_status"] == status
    result.end -= 1
    with pytest.raises(ValueError, match="supplementary.*interval"):
        select_phase3_acceptance([result], AblationID.A0)
    assert unsupported_eve_ids([result]) == {result.eve_id}


def test_json_exports_catalogue_roles_without_turning_portal_into_capsid(tmp_path):
    evidence = [
        hit("mcp", "Mirus_MCP", status="unvalidated", identity=0),
        hit("portal", "Mirus_Portal", status="unvalidated", identity=0),
    ]
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


@pytest.mark.parametrize("use_vog", [False, True])
def test_context_scans_run_once_before_scoring_and_exclude_partial_genes(tmp_path, monkeypatch, use_vog):
    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome import phase3
    from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary, boundary_candidate_id
    from virosync.pipeline.phase3.gene_functions import (
        GeneFunction,
        GeneFunctionResult,
        Profile,
        ProfileHit,
        label_gene,
    )

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
        if use_vog:
            labelled = label_gene(
                identity,
                "inside",
                model_hits=[],
                profile_hits=[ProfileHit("VOG03169", 1e-10), ProfileHit("VOG05151", 1e-9)],
                model_annotations={},
                profiles={
                    "VOG03169": Profile(
                        "VOG03169", "VOGDB latest-2025-09-30", "MCP", "capsid", 1e-9, "Major capsid protein"
                    ),
                    "VOG05151": Profile(
                        "VOG05151", "VOGDB latest-2025-09-30", "MCP", "capsid", 1e-9, "Major capsid protein"
                    ),
                },
            )
            assert labelled is not None
            return GeneFunctionResult([labelled], [])
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
    if use_vog:
        support = annotations[identity]["hallmark_support"]
        assert {(row.model, row.stage) for row in support} == {
            ("VOG03169", "phase3_function"),
            ("rve", "phase3_integration"),
        }
        summary = summarize_hallmarks([hit("core", "GVOGm0760", evalue=1e-8)], support)
        assert summary.core_count == 1
        assert summary.compatible_count == 2
        assert summary.capsid_fraction == 0
        assert summary.floor_tier == ""


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


@pytest.mark.parametrize("ordinary_score", [65.0, 200.0])
def test_frameshift_alternative_does_not_add_gene_taxonomy_or_family_credit(
    tmp_path: Path, ordinary_score: float
) -> None:
    """One physical MCP keeps both records but contributes one scoring locus."""
    from virosync.config import PipelineConfig
    from virosync.orchestration._flows.single_genome import phase3
    from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy
    from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
    from virosync.pipeline.phase3.evidence_synthesizer import taxonomy_class_votes

    ordinary = ValidatedMarkerHit(
        query_porf="s_9|aa15-237",
        scaffold="s",
        start=10774,
        end=11929,
        strand="+",
        hmm_target="VS000309",
        hmm_score=ordinary_score,
        hmm_evalue=2e-21,
        validation_status="validated",
        top10_prefixes="NCLDV__",
        top10_targets="NCLDV__reference|protein",
        top10_pidents="87.1",
        best_hit_target="NCLDV__reference|protein",
        best_hit_pident=87.1,
        best_hit_bits=465,
        has_ncldv=1,
        has_mirus=0,
        has_plv=0,
        has_vp=0,
        has_viral=1,
    )
    rescued_id = "s_VSR0123456789abcdef"
    rescued = replace(
        ordinary,
        query_porf=f"{rescued_id}|aa1-343",
        start=10873,
        end=11898,
        hmm_score=144,
        hmm_evalue=2e-46,
    )
    marker_taxonomy = GeneTaxonomy(
        "s_9",
        "s",
        10774,
        11929,
        top1_target="NCLDV__reference|protein",
        top1_prefix="NCLDV__",
        top1_pident=87.1,
        top10_prefixes=["NCLDV__"],
        top10_targets=["NCLDV__reference|protein"],
        top10_pidents=[87.1],
        top10_bits=[100.0],
        top10_evalues=[1e-10],
        top10_qcovhsp=[80.0],
        has_ncldv_mirus=True,
        has_viral=True,
        has_hit=True,
    )
    taxonomies = [
        GeneTaxonomy(f"s_{index}", "s", start, start + 100)
        for index, start in enumerate((6000, 6500, 7000, 7500, 8000))
    ] + [marker_taxonomy, replace(marker_taxonomy, porf_id=rescued_id, start=10873, end=11898)]
    core = ConfirmedSeed("original", "s", 5774, 16929)
    boundary = RefinedBoundary("s", 5774, 16929, confirmed_seeds=(core,))
    native = collect_hallmark_evidence(
        [ordinary, rescued],
        [HMMHit("s_9", "VS000309", ordinary_score, 2e-21, 65, 15, 237)],
        frameshift_ids={rescued_id},
    )
    built = phase3._build_boundary_evidence(
        boundary=boundary,
        taxonomy_index=phase3._build_scaffold_start_index(taxonomies),
        marker_index=phase3._build_scaffold_start_index([ordinary, rescued]),
        boundary_taxonomy_map={},
        proteome_index={},
        boundary_diamond_query=None,
        interproscan_map=None,
        use_precomputed_taxonomy=True,
        proteome_path=tmp_path / "unused.faa",
        config=PipelineConfig(),
        logger=logging.getLogger(__name__),
    )
    records, summary = built.gene_taxonomy
    result = EvidenceSynthesizer(
        config=EvidenceSynthesizerConfig(use_boltz=False, use_phylogenetic_validation=False, use_tmvec_database=False),
        work_dir=tmp_path,
    ).verify_eve(
        boundary,
        [],
        hallmark_hits=built.hallmarks,
        hallmark_evidence=native,
        hallmark_witnesses=find_spatial_witnesses(native),
        gene_taxonomy_records=records,
        gene_taxonomy_summary=summary,
    )

    assert result.gene_count == result.gene_taxonomy_total == 6
    assert result.gene_taxonomy_viral_top10 == result.gene_taxonomy_viral_interior == 1
    assert result.gene_taxonomy_dominant_fraction == pytest.approx(1 / 6)
    assert result.score_components["scores"]["gene_taxonomy"] == pytest.approx(1 / 6)
    assert result.family_consistency_score == 0
    assert result.hallmark_count == 1
    assert taxonomy_class_votes(built.hallmarks, gene_taxonomy_records=records) == {"NCLDV": 5}
    assert set(result.mcp_gene_ids) == {"s_9", rescued_id}
    assert result.frameshift_rescue_marker_ids == [rescued_id]
    assert result.hallmark_evidence == native
    assert result.hallmark_assessment.summary.core_count == 1
    assert result.confirmed_seeds == (core,)
    assert (result.start, result.end) == (5774, 16929)


@pytest.mark.parametrize(
    "different",
    [
        {"strand": "-"},
        {"scaffold": "other"},
        {"hallmark_gene": "GVOGm0003"},
        {"start": 400, "end": 600},
        {"start": 401, "end": 600},
        {"porf_id": "ordinary_second"},
    ],
)
def test_rescued_locus_count_keeps_independent_genes(different):
    from virosync.orchestration._flows.single_genome.phase3 import _scoring_locus_ids

    ordinary = {
        "porf_id": "ordinary",
        "scaffold": "s",
        "strand": "+",
        "start": 100,
        "end": 400,
        "hallmark_gene": "VS000309",
        "hmm_score": 65,
    }
    other = {**ordinary, "porf_id": "s_VSR0123456789abcdef", "start": 200, "end": 350, **different}
    loci = _scoring_locus_ids([ordinary, other])
    for marker in (ordinary, other):
        marker["scoring_locus_id"] = loci.get(marker["porf_id"], marker["porf_id"])
    retained, _, _ = EvidenceSynthesizer._deduplicate_hallmark_hits([ordinary, other])
    assert len(retained) == 2


def test_rescue_bridges_fragments_without_losing_raw_marker_records():
    from virosync.orchestration._flows.single_genome.phase3 import _scoring_locus_ids

    base = {"scaffold": "s", "strand": "+", "hallmark_gene": "VS000309", "hmm_score": 65}
    markers = [
        {**base, "porf_id": "left", "start": 100, "end": 300},
        {**base, "porf_id": "right", "start": 400, "end": 600},
        {**base, "porf_id": "s_VSR0123456789abcdef", "start": 200, "end": 500, "hmm_score": 144},
    ]
    for ordered in (markers, list(reversed(markers))):
        loci = _scoring_locus_ids(ordered)
        assert set(loci.values()) == {"s_VSR0123456789abcdef"}
        for marker in ordered:
            marker["scoring_locus_id"] = loci[marker["porf_id"]]
        retained, _, _ = EvidenceSynthesizer._deduplicate_hallmark_hits(ordered)
        assert len(retained) == 1
        assert len(ordered) == 3


def test_rescued_multimodel_protein_cannot_bridge_unrelated_marker_loci():
    from virosync.orchestration._flows.single_genome.phase3 import _scoring_locus_ids

    base = {"scaffold": "s", "strand": "+", "hmm_score": 65}
    rescued = "s_VSR0123456789abcdef"
    markers = [
        {**base, "porf_id": "capsid", "start": 100, "end": 300, "hallmark_gene": "VS000309"},
        {**base, "porf_id": "atpase", "start": 400, "end": 600, "hallmark_gene": "GVOGm0760"},
        {**base, "porf_id": rescued, "start": 200, "end": 500, "hallmark_gene": "VS000309"},
        {**base, "porf_id": rescued, "start": 200, "end": 500, "hallmark_gene": "GVOGm0760"},
    ]
    assert _scoring_locus_ids(markers) == {}
    assert _scoring_locus_ids(list(reversed(markers))) == {}


def test_rescue_representative_uses_shared_model_score_only():
    from virosync.orchestration._flows.single_genome.phase3 import _scoring_locus_ids

    base = {"scaffold": "s", "strand": "+", "start": 100, "end": 400, "hallmark_gene": "VS000309"}
    rescued = "s_VSR0123456789abcdef"
    markers = [
        {**base, "porf_id": "ordinary", "hmm_score": 100},
        {**base, "porf_id": rescued, "hmm_score": 80},
        {**base, "porf_id": rescued, "hmm_score": 200, "hallmark_gene": "GVOGm0760"},
    ]
    assert _scoring_locus_ids(markers) == {"ordinary": "ordinary", rescued: "ordinary"}


@pytest.mark.parametrize("batched", [False, True])
@pytest.mark.parametrize("independent_partner", [False, True])
def test_linked_translations_cannot_supply_two_independent_catalogue_cores(
    tmp_path, monkeypatch, batched, independent_partner
):
    from virosync.orchestration import tasks, utils
    from virosync.pipeline.hallmark_evidence import assess_hallmarks
    from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary, boundary_candidate_id

    rescued = "s_VSR0123456789abcdef"
    evidence = [
        hit(gene, model, 100, 400, evalue=1e-20, status="validated", identity=80)
        for gene, model in [("ordinary", "VS000309"), (rescued, "VS000309"), (rescued, "GVOGm0760")]
    ]
    if independent_partner:
        evidence.append(hit("independent", "GVOGm0760", 700, 800, evalue=1e-20, status="validated", identity=80))
    markers = [
        {
            "porf_id": item.query_porf,
            "hallmark_gene": item.model,
            "hmm_score": item.score,
            **({"scoring_locus_id": "ordinary"} if item.gene != "independent" else {}),
            "validation_status": "validated",
            "top10_prefixes": "NCLDV__",
            "top10_pidents": "80",
        }
        for item in evidence
    ]
    seed = ConfirmedSeed("seed", "s", 0, 500)
    boundary = RefinedBoundary("s", 0, 1000, confirmed_seeds=(seed,))
    raw_witnesses = find_spatial_witnesses(evidence)
    if batched:
        genome = tmp_path / "genome.fna"
        genome.write_text(">s\n" + "ACGT" * 250 + "\n")
        monkeypatch.setattr(utils, "get_genes_for_boundary", lambda **kwargs: [])
        identity = boundary_candidate_id(boundary)
        [result] = tasks.verify_eve_candidates_batched_task(
            boundaries=[boundary],
            genome_path=genome,
            work_dir=tmp_path / "work",
            proteome_path=tmp_path / "unused.faa",
            hallmark_hits_map={identity: markers},
            hallmark_evidence_map={identity: evidence},
            hallmark_witnesses_map={identity: raw_witnesses},
            candidate_annotations={},
            max_workers=1,
        )
    else:
        result = EvidenceSynthesizer(
            config=EvidenceSynthesizerConfig(
                use_boltz=False, use_phylogenetic_validation=False, use_tmvec_database=False
            ),
            work_dir=tmp_path,
        ).verify_eve(boundary, [], hallmark_hits=markers, hallmark_evidence=evidence, hallmark_witnesses=raw_witnesses)
    expected_cores = 2 if independent_partner else 1
    expected_tier = "HIGH" if independent_partner else "MEDIUM"
    assert result.hallmark_assessment.summary.core_count == expected_cores
    assert result.hallmark_assessment.summary.strong_count == expected_cores
    assert result.hallmark_assessment.summary.floor_tier == result.confidence_tier == expected_tier
    assert result.hallmark_assessment.summary.coherence_bonus == (0.05 if independent_partner else 0)
    assert any(
        witness.rule == "capsid" and witness.tier == "MEDIUM" for witness in result.hallmark_assessment.witnesses
    )
    assert any(witness.tier == "HIGH" for witness in result.hallmark_assessment.witnesses) == independent_partner
    if independent_partner:
        assert any(evidence[-1].evidence_id in witness.evidence_ids for witness in result.hallmark_assessment.witnesses)
    assert set(result.hallmark_assessment.summary.qualifying_evidence_ids) == {item.evidence_id for item in evidence}
    assert result.hallmark_evidence == evidence
    assert result.hallmark_witnesses == raw_witnesses
    assert result.confirmed_seeds == (seed,)
    assert result.to_dict()["scoring_locus_ids"] == {"ordinary": "ordinary", rescued: "ordinary"}
    result.end = 900
    reassessed = assess_hallmarks(result)
    assert reassessed.summary.core_count == expected_cores
    assert reassessed.summary.floor_tier == expected_tier
    assert any(witness.tier == "HIGH" for witness in reassessed.witnesses) == independent_partner
    assert result.start <= seed.start < seed.end <= result.end


@pytest.mark.parametrize("model", [f"EVE_CAPSID_NOVEL_{index:04d}" for index in range(1, 34)] + ["VS000382"])
def test_completed_catalogue_resolves_exact_models_and_excludes_noncapsid_context(model: str) -> None:
    """Every retained model resolves, including explicit non-scoring records."""
    evidence = hit("g", model, evalue=1e-8)
    assert evidence.model_roles
    excluded = {
        "EVE_CAPSID_NOVEL_0018",
        "EVE_CAPSID_NOVEL_0026",
    }
    assert all(role["excluded_from_scoring"] for role in evidence.model_roles) == (model in excluded)
    summary = summarize_hallmarks([evidence])
    noncore = excluded | {
        "EVE_CAPSID_NOVEL_0001",
        "EVE_CAPSID_NOVEL_0006",
        "EVE_CAPSID_NOVEL_0015",
        "EVE_CAPSID_NOVEL_0020",
    }
    assert summary.core_count == (0 if model in noncore else 1)
    assert summary.compatible_count == (0 if model in noncore else 1)
    assert summary.coherence_bonus == 0
    assert summary.floor_tier == ""
    assert summary.capsid_fraction == 0
    assert find_spatial_witnesses([evidence]) == []


@pytest.mark.parametrize(
    "models,loci",
    [
        (["EVE_CAPSID_NOVEL_0016", "VP_Penton_1", "EVE_CAPSID_NOVEL_0032"], {}),
        (["EVE_CAPSID_NOVEL_0031", "EVE_CAPSID_NOVEL_0028"], {}),
        (["EVE_CAPSID_NOVEL_0027", "EVE_CAPSID_NOVEL_0028"], {}),
        (["VS000382", "VS000106"], {}),
        (["EVE_CAPSID_NOVEL_0029", "GVOGm0760"], {"g0": "one", "g1": "one"}),
    ],
)
def test_completed_catalogue_deduplicates_families_and_scoring_loci(models: list[str], loci: dict[str, str]) -> None:
    evidence = [hit(f"g{index}", model, evalue=1e-8) for index, model in enumerate(models)]
    summary = summarize_hallmarks(evidence, scoring_locus_ids=loci)
    assert summary.core_count == summary.compatible_count == 1
    assert summary.coherence_bonus == 0
    assert summary.floor_tier == ""
    assert summary.capsid_fraction == 0


@pytest.mark.parametrize("model", ["EVE_CAPSID_NOVEL_0031", "EVE_CAPSID_NOVEL_0029"])
@pytest.mark.parametrize(
    "status,evalue,weighted",
    [
        ("unvalidated", 1e-8, 0),
        ("supported", 1e-8, 0),
        ("validated", -1, 0),
        ("validated", 0.010001, 0),
        ("validated", 0.01, 0.25),
        ("validated_novel", 0.01, 0.25),
    ],
)
def test_completed_catalogue_keeps_validation_and_evalue_gates(
    model: str, status: str, evalue: float, weighted: float
) -> None:
    evidence = [hit("g", model, status=status, evalue=evalue)]
    summary = summarize_hallmarks(evidence)
    assert summary.weighted_core_count == weighted
    if weighted == 0:
        assert summary.qualifying_evidence_ids == ()
        assert summary.floor_tier == ""
        assert find_spatial_witnesses(evidence) == []
    elif model == "EVE_CAPSID_NOVEL_0031":
        assert summary.marker_score == pytest.approx(0.05)
        assert summary.capsid_fraction == 0
        assert summary.floor_tier == ""
        assert find_spatial_witnesses(evidence) == []


def test_completed_catalogue_unresolved_panel_cannot_combine_with_established_cores() -> None:
    evidence = [
        hit("penton", "EVE_CAPSID_NOVEL_0002", evalue=1e-8),
        hit("atpase", "GVOGm0760", evalue=1e-8),
        hit("protease", "EVE_CAPSID_NOVEL_0001", evalue=1e-8),
    ]
    result = VerificationResult("eve", "s", 0, 2000, hallmark_evidence=evidence)
    result.hallmark_assessment = assess_hallmarks(result)
    summary = result.hallmark_assessment.summary
    assert summary.core_count == summary.compatible_count == 1
    assert summary.marker_score == pytest.approx(0.2)
    assert summary.coherence_bonus == summary.capsid_fraction == 0
    assert summary.floor_tier == ""
    assert find_spatial_witnesses(evidence) == []


@pytest.mark.parametrize("status", ["validated", "validated_novel", "unvalidated"])
def test_completed_catalogue_putative_mcp_keeps_core_credit_without_catalogue_capsid_effects(status: str) -> None:
    capsid = hit("capsid", "EVE_CAPSID_NOVEL_0029", evalue=0.01, status=status)
    partner = hit("atpase", "GVOGm0760", 1100, 1400, evalue=0.01)
    singleton = summarize_hallmarks([capsid])
    pair = summarize_hallmarks([capsid, partner])
    witnesses = find_spatial_witnesses([capsid, partner])
    assert singleton.floor_tier == pair.floor_tier == ""
    assert singleton.capsid_fraction == pair.capsid_fraction == 0
    assert witnesses == []
    if status == "unvalidated":
        assert singleton.core_count == 0
        assert pair.core_count == 1
    else:
        assert singleton.core_count == 1
        assert singleton.weighted_core_count == 0.25
        assert singleton.marker_score == pytest.approx(0.05)
        assert pair.panel == "Nucleocytoviricota"
        assert pair.core_count == 2
        assert pair.weighted_core_count == 1
        assert pair.marker_score == pytest.approx(0.2)


def test_completed_catalogue_penton_writer_formats_outputs_and_preserves_seed(tmp_path: Path) -> None:
    """Injected validation tests output formatting and seed retention, not acceptance."""
    from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary

    marker = ValidatedMarkerHit(
        query_porf="penton|aa1-20",
        scaffold="s",
        start=100,
        end=1000,
        strand="+",
        hmm_target="EVE_CAPSID_NOVEL_0031",
        hmm_score=70,
        hmm_evalue=1e-8,
        validation_status="validated_novel",
        top10_prefixes="",
        best_hit_target="",
        best_hit_pident=0,
        best_hit_bits=0,
        has_ncldv=0,
        has_mirus=0,
        has_plv=0,
        has_vp=0,
        has_viral=0,
    )
    evidence = collect_hallmark_evidence([marker], [HMMHit("penton", marker.hmm_target, 70, 1e-8, 60, 1, 20)])
    seed = ConfirmedSeed("original", "s", 50, 1500)
    boundary = RefinedBoundary("s", 50, 1500, confirmed_seeds=(seed,))
    result = EvidenceSynthesizer(
        config=EvidenceSynthesizerConfig(use_boltz=False, use_phylogenetic_validation=False, use_tmvec_database=False),
        work_dir=tmp_path,
    ).verify_eve(boundary, [], hallmark_evidence=evidence, hallmark_witnesses=find_spatial_witnesses(evidence))
    generator = OutputGenerator(tmp_path / "export")
    canonical_path = generator.write_predictions_tsv([result])
    detailed_path = generator.write_predictions_detailed_tsv([result])
    profiles_path = generator.write_evidence_profiles([result])

    assert len(evidence) == 1
    assert result.hallmark_assessment.summary.core_count == 1
    assert result.hallmark_assessment.summary.marker_score == pytest.approx(0.2)
    assert result.hallmark_assessment.summary.floor_tier == ""
    assert result.hallmark_assessment.witnesses == ()
    assert result.confirmed_seeds == (seed,)
    assert (result.start, result.end) == (50, 1500)
    with canonical_path.open() as stream:
        assert next(csv.DictReader(stream, delimiter="\t"))["qualifying_hallmark_genes"] == "Penton"
    with detailed_path.open() as stream:
        assert next(csv.DictReader(stream, delimiter="\t"))["qualifying_hallmark_genes"] == "Penton"
    profile = json.loads(profiles_path.read_text())[result.eve_id]
    assert profile["hallmark_evidence"][0]["model_roles"][0]["panel"] == "Virophages"
    assert profile["hallmark_evidence"][0]["validation_status"] == "validated_novel"


def test_completed_catalogue_has_unique_consistent_noncapsid_links() -> None:
    links = json.loads(CATALOGUE_PATH.read_text(encoding="utf-8"))["links"]
    stage_keys = [(link["stage"], link["model"], link["panel"]) for link in links]
    unresolved = {link["model"] for link in links if link["panel"].startswith("Unresolved_")}
    established = {
        link["model"] for link in links if not link["panel"].startswith("Unresolved_") and link["role"] == "core_marker"
    }
    added = [link for link in links if link["model"].startswith("EVE_CAPSID_NOVEL_") or link["model"] == "VS000382"]
    pairs = [(link["model"], link["panel"]) for link in added]

    assert len(stage_keys) == len(set(stage_keys))
    assert len(pairs) == len(set(pairs))
    assert unresolved.isdisjoint(established)
    assert len({link["model"] for link in added}) == 34
    assert all(link["capsid"] is False for link in added)


@pytest.mark.parametrize(
    "model,panels",
    [
        ("VS000652", {"Nucleocytoviricota", "Polinton_like_viruses", "Virophages"}),
        ("VS000783", {"Nucleocytoviricota", "Polinton_like_viruses", "Virophages"}),
        ("VS000745", {"Nucleocytoviricota", "Polinton_like_viruses"}),
        ("VS000788", {"Nucleocytoviricota", "Polinton_like_viruses"}),
        ("VS000791", {"Nucleocytoviricota", "Polinton_like_viruses"}),
    ],
)
def test_shared_mcp_alias_supplies_existing_capsid_permissions(model, panels):
    evidence = hit("capsid", model, evalue=1e-8)
    summary = summarize_hallmarks([evidence])
    witnesses = find_spatial_witnesses([evidence])

    assert {role["panel"] for role in evidence.model_roles} == panels
    assert summary.core_count == 1
    assert summary.capsid_fraction == 1
    assert summary.floor_tier == "MEDIUM"
    assert {w.rule for w in witnesses} == {"capsid"}
    assert add_rescue_seeds([], witnesses)


@pytest.mark.parametrize("same_locus", [False, True])
def test_shared_mcp_aliases_do_not_multiply_family_or_locus(same_locus):
    evidence = [hit("a", "VS000652", evalue=1e-8), hit("b", "VS000783", evalue=1e-8)]
    summary = summarize_hallmarks(evidence, scoring_locus_ids={"b": "a"} if same_locus else {})
    assert summary.core_count == summary.compatible_count == 1
    assert summary.coherence_bonus == 0


@pytest.mark.parametrize("model", ["EVE_CAPSID_NOVEL_0019", "EVE_CAPSID_NOVEL_0027", "EVE_CAPSID_NOVEL_0031"])
def test_penton_alias_joins_virophage_core_without_extra_minor_capsid_family(model):
    evidence = [hit("penton", model, evalue=1e-8), hit("other", "VP_Penton_1", evalue=1e-8)]
    assert summarize_hallmarks(evidence[:1]).panel == "Virophages"
    summary = summarize_hallmarks(evidence)
    assert summary.panel == "Virophages"
    assert summary.core_count == summary.compatible_count == 1
    assert summary.capsid_fraction == 0
    assert find_spatial_witnesses(evidence) == []


@pytest.mark.parametrize(
    "models,core",
    [
        (["EVE_CAPSID_NOVEL_0001", "VS000374"], "VP_ATPase_1"),
        (["VS000151", "VS000375"], "GVOGm0760"),
        (["VS000146", "VS000149", "VS000072", "VS000164"], "GVOGm0760"),
    ],
)
def test_real_context_families_count_once_without_creating_core_or_witness(models, core):
    context = [hit(f"g{i}", model, evalue=1e-8) for i, model in enumerate(models)]
    alone = summarize_hallmarks(context)
    combined = summarize_hallmarks([hit("core", core, evalue=1e-8), *context])

    assert alone.core_count == alone.compatible_count == 0
    assert alone.marker_score == 0
    assert find_spatial_witnesses(context) == []
    assert add_rescue_seeds([], find_spatial_witnesses(context)) == []
    assert combined.core_count == 1
    assert combined.compatible_count == 2
    assert combined.coherence_bonus == pytest.approx(0.05)
    assert combined.floor_tier == ""
    assert find_spatial_witnesses([hit("core", core, evalue=1e-8), *context]) == []


@pytest.mark.parametrize(
    "evalue,status,loci,compatible",
    [
        (1e-5, "validated", {}, 2),
        (1.1e-5, "validated", {}, 1),
        (1e-8, "unvalidated", {}, 1),
        (1e-8, "validated", {"context": "core"}, 1),
    ],
)
def test_real_context_keeps_strength_validation_and_locus_gates(evalue, status, loci, compatible):
    evidence = [hit("core", "GVOGm0760", evalue=1e-8), hit("context", "VS000375", evalue=evalue, status=status)]
    summary = summarize_hallmarks(evidence, scoring_locus_ids=loci)
    assert summary.core_count == 1
    assert summary.compatible_count == compatible
    assert summary.floor_tier == ""


def test_ncldv_context_cannot_support_virophage_core():
    evidence = [hit("core", "VP_ATPase_1", evalue=1e-8), hit("context", "VS000375", evalue=1e-8)]
    summary = summarize_hallmarks(evidence)
    assert summary.panel == "Virophages"
    assert summary.compatible_count == summary.core_count == 1
    assert summary.coherence_bonus == 0


@pytest.mark.parametrize(
    "models,panel",
    [
        (["EVE_CAPSID_NOVEL_0029", "GVOGm0760", "VS000106"], "Nucleocytoviricota"),
        (["EVE_CAPSID_NOVEL_0016", "VP_ATPase_1", "VP_PRO_1"], "Virophages"),
    ],
)
def test_completed_catalogue_noncapsid_cores_join_three_core_medium_witnesses(models: list[str], panel: str) -> None:
    evidence = [hit(f"g{index}", model, index * 1100, index * 1100 + 1000) for index, model in enumerate(models)]
    witnesses = find_spatial_witnesses(evidence)
    result = VerificationResult(
        "eve", "s", 0, 4000, gene_count=5, hallmark_evidence=evidence, hallmark_witnesses=witnesses
    )

    calculate_eve_confidence(result, 0)

    assert result.hallmark_assessment.summary.core_count == 3
    assert result.hallmark_assessment.summary.capsid_fraction == 0
    assert result.hallmark_assessment.summary.floor_tier == ""
    assert {(w.panel, w.rule, w.tier) for w in result.hallmark_assessment.witnesses} == {
        (panel, "three_core_15kb", "MEDIUM")
    }
    assert result.score_components["hallmark_floor"] == pytest.approx(0.2)


def test_completed_catalogue_virophage_penton_partners_real_mcp_for_high_floor() -> None:
    evidence = [hit("penton", "EVE_CAPSID_NOVEL_0016"), hit("mcp", "VP_MCP_6", 1100, 2000)]
    witnesses = find_spatial_witnesses(evidence)
    result = VerificationResult(
        "eve", "s", 0, 3000, gene_count=5, hallmark_evidence=evidence, hallmark_witnesses=witnesses
    )

    confidence = calculate_eve_confidence(result, 0)

    assert result.hallmark_assessment.summary.core_count == 2
    assert {(w.rule, w.tier) for w in result.hallmark_assessment.witnesses} == {
        ("capsid", "MEDIUM"),
        ("capsid_partner_5kb", "HIGH"),
    }
    assert result.score_components["hallmark_floor"] == pytest.approx(0.7)
    assert confidence >= 0.7


@pytest.fixture
def explicit_support_catalogue(monkeypatch: pytest.MonkeyPatch) -> None:
    """Add synthetic exact links without assigning a biological model function."""
    from virosync.pipeline import hallmark_evidence

    catalogue = dict(hallmark_evidence._catalogue())
    link = {
        "panel": "Nucleocytoviricota",
        "marker": "Synthetic_context",
        "module": "context",
        "role": "support_only",
        "capsid": False,
        "structural_anchor": False,
        "excluded_from_scoring": False,
    }
    catalogue["phase1_marker", "TEST_SUPPORT_A"] = [link]
    catalogue["phase1_marker", "TEST_SUPPORT_B"] = [link]
    catalogue["gene_function", "TEST_SUPPORT_A"] = [link]
    monkeypatch.setattr(hallmark_evidence, "_catalogue", lambda: catalogue)


def test_explicit_support_alone_cannot_supply_core_score_witness_or_seed(explicit_support_catalogue: None) -> None:
    evidence = [hit("context", "TEST_SUPPORT_A", evalue=1e-8)]

    summary = summarize_hallmarks(evidence)
    witnesses = find_spatial_witnesses(evidence)

    assert summary.core_count == summary.compatible_count == 0
    assert summary.marker_score == summary.coherence_bonus == 0
    assert summary.qualifying_evidence_ids == ()
    assert summary.floor_tier == ""
    assert witnesses == []
    assert add_rescue_seeds([], witnesses) == []


@pytest.mark.parametrize("evalue,compatible", [(1e-5, 2), (1.0001e-5, 1), (0.01, 1)])
@pytest.mark.parametrize("status", ["validated", "validated_novel"])
def test_explicit_support_requires_same_strong_threshold_in_both_stages(
    explicit_support_catalogue: None, evalue: float, compatible: int, status: str
) -> None:
    core = hit("core", "GVOGm0760", evalue=0.01)
    context = hit("context", "TEST_SUPPORT_A", evalue=evalue, status=status)

    phase1 = summarize_hallmarks([core, context])
    phase3 = summarize_hallmarks([core], [HallmarkSupport("context", "TEST_SUPPORT_A", evalue, "gene_function")])

    assert phase1 == phase3
    assert phase1.compatible_count == compatible
    assert phase1.core_count == 1
    assert phase1.weighted_core_count == 0.25


@pytest.mark.parametrize("loci,compatible,bonus", [({}, 2, 0.05), ({"core": "one", "a": "one", "b": "one"}, 1, 0)])
def test_explicit_support_deduplicates_function_and_physical_locus(
    explicit_support_catalogue: None, loci: dict[str, str], compatible: int, bonus: float
) -> None:
    core = hit("core", "GVOGm0760", evalue=1e-8)
    evidence = [core, hit("a", "TEST_SUPPORT_A", evalue=1e-8), hit("b", "TEST_SUPPORT_B", evalue=1e-8)]

    summary = summarize_hallmarks(evidence, scoring_locus_ids=loci)

    assert summary.core_count == 1
    assert summary.compatible_count == compatible
    assert summary.coherence_bonus == bonus
    assert summary.qualifying_evidence_ids == (core.evidence_id,)
    assert find_spatial_witnesses(evidence) == []


@pytest.mark.parametrize("other_bonus,base_score,supported_score", [(0.1, 0.175, 0.225), (0.9, 0.699, 0.699)])
def test_explicit_support_changes_single_core_score_without_bypassing_high_cap(
    explicit_support_catalogue: None, other_bonus: float, base_score: float, supported_score: float
) -> None:
    core = hit("core", "GVOGm0760", evalue=1e-8)
    baseline = VerificationResult(
        "eve", "s", 0, 2000, gene_count=5, hallmark_evidence=[core], clustering_bonus=other_bonus
    )
    supported = replace(baseline, hallmark_evidence=[core, hit("context", "TEST_SUPPORT_A", evalue=1e-8)])

    assert calculate_eve_confidence(baseline, 0) == pytest.approx(base_score)
    assert calculate_eve_confidence(supported, 0) == pytest.approx(supported_score)
    assert supported.hallmark_assessment.summary.core_count == 1
    assert supported.hallmark_assessment.summary.coherence_bonus == pytest.approx(0.05)
    assert supported.hallmark_assessment.summary.floor_tier == ""
    assert supported.hallmark_assessment.witnesses == ()
