"""The reference class EUK-pEVE counts as viral and never as cellular or host.

Resource bundle v1.2.0 labels viral-like proteins of eukaryotic assemblies
``EUK-pEVE__<assembly>``. ``"EUK-pEVE__".startswith("EUK__")`` is False, so
before this rule every vocabulary read the hit as unknown. These tests pin the
decisions: the class is viral in every viral/cellular vocabulary, never
cellular or host, and carries no viral family.
"""

from __future__ import annotations

import csv
from pathlib import Path
from types import SimpleNamespace

import pytest

from virosync.config.thresholds import get_config
from virosync.orchestration._flows.single_genome.phase3 import _summarize_boundary_taxonomy
from virosync.output_contract import (
    DETAILED_TAXONOMY_PARTITION,
    EFFECTIVE_EVE_CLASSES,
    LINEAGE_EVE_CLASSES,
)
from virosync.pipeline.host_signatures import (
    DOMAIN_TOKENS,
    HostSignatureModel,
    TaxonomyLabelLookup,
    build_host_signature_model,
    dominant_lineage_tokens,
    host_lineage_share,
    host_signature_density_evalue_weighted,
    set_taxonomy_lookup,
)
from virosync.pipeline.phase1.frameshift_screening import diamond_query_coverages
from virosync.pipeline.phase1.hhg_seeding import Anchor, HMMHit
from virosync.pipeline.phase1.marker_validation import (
    VALIDATED_PREFIXES,
    NovelMarkerCriteria,
    ValidationStatus,
    collect_host_signatures,
    filter_validated_markers,
    parse_taxonomy_prefix,
    validate_hmm_hit,
)
from virosync.pipeline.phase1.seed_merger import MergedSeed
from virosync.pipeline.phase1.taxonomy_expansion import is_gene_viral_positive
from virosync.pipeline.phase2.boundary_diamond import (
    DiamondHit,
    build_gene_taxonomy_record,
    classify_all_porfs,
    has_identity_qualified_viral_hit,
    pORF,
)
from virosync.pipeline.phase2.boundary_diamond import (
    extract_prefix as boundary_extract_prefix,
)
from virosync.pipeline.phase2.boundary_refiner import (
    RefinedBoundary,
    _classify_gene_label,
    _has_viral_taxonomy,
    should_trim_gene_as_host,
)
from virosync.pipeline.phase3.evidence_synthesizer import (
    EvidenceSynthesizer,
    EvidenceSynthesizerConfig,
    VerificationResult,
    _high_confidence_euk_gene_count,
    consensus_taxonomy_class,
    infer_likely_family,
    marker_taxonomy_category,
    mcp_vote_classes,
    taxonomy_class_votes,
)
from virosync.pipeline.phase3.gene_taxonomy import (
    CELLULAR_NAMESPACES,
    VIRAL_PREFIXES,
    GeneDiamondHit,
    GeneTaxonomy,
    classify_gene_taxonomy,
    extract_prefix,
    extract_raw_prefix,
    has_only_peve_viral_hits,
    qualified_viral_hits,
    summarize_dominant_family,
    viral_hit_categories,
)
from virosync.pipeline.phase3.output_generator import (
    OutputGenerator,
    _taxonomy_partition_bucket,
    evaluate_v2_quality_gate,
)
from virosync.pipeline.phase3.phase1_surface import build_phase1_seed_surface
from virosync.pipeline.phase3.phylogenetic_validation import (
    GVClassValidation,
    PhylogeneticValidator,
)
from virosync.pipeline.phase3.structural_homology import BoltzFoldSeekAnalyzer
from virosync.pipeline.search_backend import select_target_hits
from virosync.validation.tsv_invariants import run_tsv_invariant_checks

PEVE_TARGET = "EUK-pEVE__EP00022|prot1"


def _kept(subject_ids: list[str]) -> list[str]:
    return [subject_ids[idx] for idx in select_target_hits(subject_ids, 10)]


def _marker(prefixes: str, pidents: str, *, gene: str = "gvogm0100", porf_id: str = "") -> dict:
    return {
        "hallmark_gene": gene,
        "porf_id": porf_id,
        "top10_prefixes": prefixes,
        "top10_pidents": pidents,
        "top10_targets": "",
        "validation_status": "validated",
    }


def _gene(porf_id: str, prefixes: str, pidents: str) -> dict:
    return {
        "porf_id": porf_id,
        "top10_prefixes": prefixes,
        "top10_pidents": pidents,
        "top10_targets": "",
        "is_flanking": False,
    }


# --- Phase 3 gene taxonomy vocabulary ---


def test_gene_taxonomy_reports_peve_as_its_own_viral_token() -> None:
    assert extract_prefix(PEVE_TARGET) == "EUK-pEVE"
    assert extract_raw_prefix(PEVE_TARGET) == "EUK-pEVE"
    assert "EUK-pEVE" in VIRAL_PREFIXES
    assert "EUK-pEVE" not in CELLULAR_NAMESPACES
    assert extract_prefix("EUK__EP00022|prot1") == "EUK"


def test_peve_hit_makes_a_gene_viral_but_sets_no_family_flag() -> None:
    gene = classify_gene_taxonomy(
        "porf1",
        [GeneDiamondHit("porf1", PEVE_TARGET, 1e-30, 200.0, 80.0, 90.0)],
        start=1,
        end=100,
    )

    assert gene.top1_prefix == "EUK-pEVE"
    assert gene.has_viral is True
    assert gene.has_ncldv_mirus is False
    assert gene.has_vp_plv is False
    assert gene.is_high_pident_euk is False


def test_peve_hit_below_the_identity_cutoff_is_not_viral() -> None:
    gene = classify_gene_taxonomy(
        "porf1",
        [GeneDiamondHit("porf1", PEVE_TARGET, 1e-30, 200.0, 24.9, 90.0)],
        start=1,
        end=100,
    )

    assert gene.has_viral is False


def test_peve_hits_give_no_dominant_family() -> None:
    taxonomies = [
        GeneTaxonomy("p1", 1, 100, top10_prefixes=["EUK-pEVE"], top10_pidents=[90.0]),
        GeneTaxonomy("p2", 101, 200, top10_prefixes=["EUK-pEVE__"], top10_pidents=[80.0]),
    ]

    assert summarize_dominant_family(taxonomies) == ("UNKNOWN", 0.0)


def test_qualified_viral_hits_keep_the_peve_spelling() -> None:
    record = {
        "top10_prefixes": "EUK-pEVE__,PPV__,EUK-pEVE__",
        "top10_pidents": "80,60,20",
        "top10_targets": f"{PEVE_TARGET},PPV__A|p,EUK-pEVE__B|q",
    }

    assert qualified_viral_hits(record) == [("EUK-pEVE", PEVE_TARGET), ("PPV", "PPV__A|p")]


@pytest.mark.parametrize("prefixes", ["EUK-pEVE__", "EUK-pEVE", "EUK-PEVE"])
def test_peve_only_hits_vote_viral_unknown_in_every_spelling(prefixes: str) -> None:
    record = {"top10_prefixes": prefixes, "top10_pidents": "80"}

    assert viral_hit_categories(record) == ["VIRAL_UNKNOWN"]
    assert has_only_peve_viral_hits([record]) is True


def test_peve_hit_never_overrides_a_lineage_hit() -> None:
    record = {"top10_prefixes": "EUK-pEVE__,PPV__", "top10_pidents": "90,40"}

    assert viral_hit_categories(record) == ["PPV"]
    assert has_only_peve_viral_hits([record]) is False


def test_low_identity_peve_hit_casts_no_vote() -> None:
    record = {"top10_prefixes": "EUK-pEVE__", "top10_pidents": "24.9"}

    assert viral_hit_categories(record) == []
    assert has_only_peve_viral_hits([record]) is False


# --- Phase 2 boundary vocabulary ---


def test_boundary_prefix_reads_peve_as_viral() -> None:
    assert boundary_extract_prefix(PEVE_TARGET) == "EUK-pEVE__"
    assert boundary_extract_prefix("EUK__EP00022|prot1") == "EUK__"
    assert has_identity_qualified_viral_hit(["EUK-pEVE__"], [80.0])
    assert has_identity_qualified_viral_hit(["EUK-pEVE"], [25.0])
    assert not has_identity_qualified_viral_hit(["EUK-pEVE__"], [24.9])


def test_boundary_classification_counts_peve_gene_as_viral_without_family_flags() -> None:
    proteome_index = {"ctg": [pORF("g1", "ctg", 0, 300)]}
    hits = {"g1": [DiamondHit("g1", PEVE_TARGET, 1e-30, 200.0, 80.0, 90.0)]}

    taxonomy = classify_all_porfs(["g1"], hits, proteome_index, host_prefix="EUK__")["g1"]

    assert taxonomy.top1_prefix == "EUK-pEVE__"
    assert taxonomy.has_viral is True
    assert taxonomy.has_ncldv_mirus is False
    assert taxonomy.has_vp_plv is False
    record = build_gene_taxonomy_record(taxonomy)
    assert record["top10_prefixes"] == ["EUK-pEVE"]
    assert record["has_viral"] is True


def _tax(**overrides):
    values = {
        "has_hit": True,
        "has_ncldv_mirus": False,
        "has_vp_plv": False,
        "has_viral": True,
        "top1_prefix": "EUK-pEVE__",
        "top10_prefixes": ["EUK-pEVE__"],
        "top10_pidents": [80.0],
        "top10_bits": [200.0],
        "taxonomy_fingerprint": None,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def test_peve_gene_counts_as_viral_for_gap_merging() -> None:
    assert _has_viral_taxonomy(_tax()) is True
    assert _has_viral_taxonomy(_tax(top10_pidents=[24.9], has_viral=False)) is False


def test_peve_gene_is_labelled_viral_and_never_trimmed_as_host() -> None:
    tax = _tax(top1_prefix="EUK__", top10_prefixes=["EUK__", "EUK-pEVE__"], top10_pidents=[92.0, 45.0])

    trim, score = should_trim_gene_as_host(
        gene_taxonomy=tax,
        control_stats=None,
        host_prefix="EUK__",
        host_baseline_fingerprint={},
    )
    label = _classify_gene_label(
        tax=_tax(),
        control_stats=None,
        host_prefix="EUK__",
        host_baseline_fingerprint={},
        min_overlap_score=0.40,
        neighbor_context=None,
        unknown_host_penalty=2.0,
        unknown_viral_bonus=2.0,
        host_signature_model=None,
        host_signature_threshold=0.5,
    )

    assert trim is False
    assert score == 0.0
    assert label == "V"


# --- Phase 1 marker validation ---


def _hmm_hit(model: str = "VS000803") -> HMMHit:
    return HMMHit(
        query_name="contig_1_1",
        target_name=model,
        score=120.0,
        evalue=1e-40,
        domain_score=120.0,
        query_start=1,
        query_end=100,
    )


def _validate(tmp_path: Path, rows: list[tuple[str, float]], model: str = "VS000803"):
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(f">contig_1_1 # 1 # 300 # + # ID=1_1;partial=00\n{'M' * 100}\n")
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("".join(f"contig_1_1\t{t}\t1e-40\t{200 - n}\t{p}\t95.0\n" for n, (t, p) in enumerate(rows)))
    return filter_validated_markers([_hmm_hit(model)], diamond, proteome)


def test_marker_prefix_parser_reads_peve_as_its_own_prefix() -> None:
    assert parse_taxonomy_prefix(PEVE_TARGET) == "EUK-pEVE__"
    assert parse_taxonomy_prefix("EUK__EP00022|prot1") == "EUK__"


@pytest.mark.parametrize(
    ("pident", "expected"),
    [(80.0, ValidationStatus.VALIDATED), (24.9, ValidationStatus.UNVALIDATED)],
)
def test_peve_reference_validates_a_marker_at_the_viral_identity_cutoff(
    pident: float,
    expected: ValidationStatus,
) -> None:
    status = validate_hmm_hit(
        hmm_score=100.0,
        hmm_coverage=0.9,
        diamond_hits=[(PEVE_TARGET, 200.0, pident, 1e-40)],
        novel_criteria=NovelMarkerCriteria(),
    )

    assert status == expected


def test_peve_only_top10_validates_a_marker_that_no_host_signature_collects(tmp_path: Path) -> None:
    markers = _validate(tmp_path, [(PEVE_TARGET, 80.0), ("EUK-pEVE__EP00023|prot2", 70.0)])

    assert markers[0].validation_status == "validated"
    assert markers[0].has_viral == 1
    assert markers[0].top10_prefixes.split(",")[:2] == ["EUK-pEVE__", "EUK-pEVE__"]
    assert collect_host_signatures(markers) == set()


def test_unvalidated_peve_best_hit_is_never_a_host_signature(tmp_path: Path) -> None:
    markers = _validate(tmp_path, [(PEVE_TARGET, 20.0)])

    assert markers[0].validation_status == "unvalidated"
    assert collect_host_signatures(markers) == set()


def test_cellular_best_hit_still_gives_a_host_signature(tmp_path: Path) -> None:
    markers = _validate(tmp_path, [("EUK__EP00022|prot1", 80.0)])

    assert markers[0].validation_status == "unvalidated"
    assert collect_host_signatures(markers) == {"EUK__EP00022"}


def test_taxonomy_expansion_reads_peve_as_viral_positive() -> None:
    assert is_gene_viral_positive([("EUK__host|p1", 95.0), (PEVE_TARGET, 25.0)])
    assert not is_gene_viral_positive([(PEVE_TARGET, 24.9)])


# --- Twin rule ---


def test_twin_rule_keeps_the_peve_row_when_the_cellular_row_comes_first() -> None:
    subjects = ["EUK__G1|p1", "EUK__G2|p5", "EUK-pEVE__G1|p1"]

    assert _kept(subjects) == ["EUK-pEVE__G1|p1", "EUK__G2|p5"]


def test_twin_rule_keeps_the_peve_row_when_it_comes_first() -> None:
    subjects = ["EUK-pEVE__G1|p1", "EUK__G2|p5", "EUK__G1|p1"]

    assert _kept(subjects) == ["EUK-pEVE__G1|p1", "EUK__G2|p5"]


@pytest.mark.parametrize("cellular", ["BAC", "ARC", "MITO", "PLASTID"])
def test_twin_rule_drops_every_cellular_twin_of_a_peve_row(cellular: str) -> None:
    assert _kept([f"{cellular}__G1|p1", "EUK-pEVE__G1|p1"]) == ["EUK-pEVE__G1|p1"]


def test_peve_row_alone_and_other_proteins_of_its_assembly_survive() -> None:
    assert _kept(["EUK-pEVE__G1|p1"]) == ["EUK-pEVE__G1|p1"]
    assert _kept(["EUK-pEVE__G1|p1", "EUK__G1|p2"]) == ["EUK-pEVE__G1|p1", "EUK__G1|p2"]


# --- Host signatures and lineage ---


HOST = "EUK|NA|Bigyra|Thraustochytrida|Thraustochytriidae|Aurantiochytrium|Aurantiochytrium limacinum"
PEVE_SUFFIXED = "EUK-pEVE|NA|Bigyra-pEVE|Thraustochytrida-pEVE|NA|NA|Aurantiochytrium limacinum-pEVE"
PEVE_UNSUFFIXED = "EUK-pEVE|NA|Bigyra|Fakeclade|Thraustochytriidae|Aurantiochytrium|Aurantiochytrium limacinum"
LABELS = {
    "EUK__SELF": HOST,
    "EUK-pEVE__SUFFIXED": PEVE_SUFFIXED,
    "EUK-pEVE__UNSUFFIXED": PEVE_UNSUFFIXED,
}


def _host_model() -> HostSignatureModel:
    return HostSignatureModel(
        token_weights={
            "bigyra": 90.0,
            "thraustochytrida": 80.0,
            "thraustochytriidae": 70.0,
            "aurantiochytrium": 60.0,
            "aurantiochytrium limacinum": 50.0,
        },
        host_prefixes=["EUK__"],
    )


def _record(targets: list[str]) -> dict:
    return {
        "top10_targets": targets,
        "top10_bitscores": [100.0] * len(targets),
        "top10_pidents": [90.0] * len(targets),
        "top10_evalues": [1e-40] * len(targets),
    }


@pytest.fixture
def lookup():
    set_taxonomy_lookup(TaxonomyLabelLookup(LABELS))
    yield LABELS
    set_taxonomy_lookup(None)


def test_peve_hits_have_zero_host_lineage_share(lookup) -> None:
    model = _host_model()

    assert host_lineage_share(_record(["EUK__SELF"]), model) == 1.0
    assert host_lineage_share(_record(["EUK-pEVE__SUFFIXED"] * 3), model) == 0.0
    # A pEVE lineage that kept the host's own rank names is still not host.
    assert host_lineage_share(_record(["EUK-pEVE__UNSUFFIXED"] * 3), model) == 0.0


def test_peve_hit_dilutes_the_host_share_of_a_mixed_gene(lookup) -> None:
    share = host_lineage_share(_record(["EUK__SELF", "EUK-pEVE__UNSUFFIXED"]), _host_model())

    # Rank weights are 10 and 9; only the host hit carries host lineage.
    assert share == pytest.approx(10 / 19)


def test_peve_lineage_never_enters_the_dominant_host_lineage(lookup) -> None:
    model = _host_model()
    lineages = [PEVE_SUFFIXED, PEVE_UNSUFFIXED, HOST]

    tokens = dominant_lineage_tokens(model, lineages)

    assert not any(token.endswith("-peve") for token in tokens)
    # The unsuffixed pEVE lineage lists the host species first; it must not win.
    assert "fakeclade" not in tokens
    assert "aurantiochytrium" in tokens


def test_dominant_lineage_drops_tokens_that_end_in_peve() -> None:
    model = HostSignatureModel(token_weights={"aurantiochytrium limacinum-peve": 50.0, "bigyra-peve": 90.0})
    host_rooted = "EUK|NA|Bigyra-pEVE|Thraustochytrida-pEVE|NA|NA|Aurantiochytrium limacinum-pEVE"

    assert dominant_lineage_tokens(model, [host_rooted]) == frozenset()
    assert dominant_lineage_tokens(model, [PEVE_SUFFIXED]) == frozenset()


def test_euk_peve_is_a_domain_token() -> None:
    assert "euk-peve" in DOMAIN_TOKENS


def test_marker_with_a_peve_hit_is_left_out_of_the_host_model(lookup) -> None:
    marker = {
        "hmm_target": "COG0085",
        "top10_prefixes": "EUK__,EUK-pEVE__",
        "top10_targets": "EUK__SELF|p1,EUK-pEVE__SUFFIXED|p2",
        "top10_bitscores": "100,90",
        "top10_pidents": "90,80",
        "top10_evalues": "1e-40,1e-30",
    }
    host_marker = {**marker, "top10_prefixes": "EUK__", "top10_targets": "EUK__SELF|p1"}

    skipped = build_host_signature_model([marker], host_prefixes={"EUK__"})
    kept = build_host_signature_model([host_marker], host_prefixes={"EUK__"})

    assert skipped.token_weights == {}
    assert kept.token_weights
    assert not any("peve" in token for token in kept.token_weights)


def test_pure_peve_gene_is_not_host_like(lookup) -> None:
    model = _host_model()
    records = [{**_record(["EUK-pEVE__UNSUFFIXED"]), "top1_target": "EUK-pEVE__UNSUFFIXED|p1"}]

    host_like, _mean, _weighted = host_signature_density_evalue_weighted(records, model, host_prefixes=["EUK__"])

    assert host_like == 0


# --- Phase 3 scoring helpers ---


@pytest.mark.parametrize("has_viral", [False, True])
def test_three_high_identity_peve_hits_are_never_a_host_gene(has_viral: bool) -> None:
    record = {
        "is_flanking": False,
        "has_hit": True,
        "has_viral": has_viral,
        "top10_prefixes": ["EUK-pEVE", "EUK-pEVE", "EUK-pEVE"],
        "top10_pidents": [90.0, 85.0, 80.0],
    }
    host_record = {**record, "has_viral": False, "top10_prefixes": ["EUK", "EUK", "EUK"]}

    assert _high_confidence_euk_gene_count([record]) == 0
    assert _high_confidence_euk_gene_count([host_record]) == 1


def test_marker_with_only_peve_hits_votes_viral_unknown() -> None:
    assert marker_taxonomy_category(_marker("EUK-pEVE__,EUK-pEVE__", "80,70")) == "VIRAL_UNKNOWN"
    assert marker_taxonomy_category(_marker("EUK-pEVE__,PPV__", "80,40")) == "PPV"
    assert marker_taxonomy_category(_marker("EUK-pEVE__", "24.9")) == ""


def test_region_with_only_peve_hits_is_viral_unknown_never_a_peve_class() -> None:
    markers = [_marker("EUK-pEVE__", "80")]
    genes = [_gene("g1", "EUK-pEVE", "70"), _gene("g2", "EUK-pEVE", "60")]

    taxonomy_class = consensus_taxonomy_class(markers, gene_taxonomy_records=genes)

    assert taxonomy_class == "VIRAL_UNKNOWN"
    assert taxonomy_class in EFFECTIVE_EVE_CLASSES
    assert taxonomy_class not in LINEAGE_EVE_CLASSES


def test_peve_gene_does_not_change_the_lineage_of_a_ppv_region() -> None:
    markers = [_marker("PPV__", "80", porf_id="g1")]
    genes = [_gene("g1", "PPV,EUK-pEVE", "80,70"), _gene("g2", "PPV", "60"), _gene("g3", "PPV,EUK-pEVE", "50,90")]

    assert consensus_taxonomy_class(markers, gene_taxonomy_records=genes) == "PPV"


def test_mcp_marker_with_only_peve_hits_does_not_override_lineage_votes() -> None:
    mcp = _marker("EUK-pEVE__", "80", gene="EVE_CAPSID_NOVEL_0029", porf_id="gm")
    genes = [_gene("g1", "NCLDV", "80"), _gene("g2", "NCLDV", "70"), _gene("g3", "NCLDV", "60")]

    assert mcp_vote_classes([mcp]) == set()
    assert consensus_taxonomy_class([mcp], gene_taxonomy_records=genes) == "NCLDV"


def _lineage_markers_and_classless_peve_markers() -> tuple[list[dict], list[dict]]:
    """Two NCLDV markers and two markers with no profile class whose only hits are EUK-pEVE."""
    markers = [
        _marker("NCLDV", "40", gene="GVOGm0013", porf_id="sc|1"),
        _marker("NCLDV", "45", gene="GVOGm0022", porf_id="sc|2"),
        _marker("EUK-pEVE", "60", gene="GVOGm0100", porf_id="sc|3"),
        _marker("EUK-pEVE", "55", gene="EVE_CAPSID_NOVEL_0002", porf_id="sc|4"),
    ]
    genes = [
        _gene("sc|1", "NCLDV", "40"),
        _gene("sc|2", "NCLDV", "45"),
        _gene("sc|3", "EUK-pEVE", "60"),
        _gene("sc|4", "EUK-pEVE", "55"),
    ]
    return markers, genes


def test_peve_only_markers_cast_no_vote_beside_lineage_markers() -> None:
    markers, genes = _lineage_markers_and_classless_peve_markers()

    assert taxonomy_class_votes(markers, gene_taxonomy_records=genes) == {"NCLDV": 6}
    assert consensus_taxonomy_class(markers, gene_taxonomy_records=genes) == "NCLDV"


@pytest.mark.parametrize("lineage", ["NCLDV", "MIRUS", "PPV", "CRESS", "PHAGE"])
def test_peve_only_genes_never_dilute_a_lineage_vote(lineage: str) -> None:
    genes = [_gene(f"l{index}", lineage, "80") for index in range(3)]
    genes += [_gene(f"e{index}", "EUK-pEVE", "70") for index in range(7)]

    assert taxonomy_class_votes([], gene_taxonomy_records=genes) == {lineage: 3}
    assert consensus_taxonomy_class([], gene_taxonomy_records=genes) == lineage


def test_peve_only_votes_stand_when_no_gene_votes_for_a_lineage() -> None:
    markers = [_marker("EUK-pEVE", "80", porf_id="g1")]
    genes = [_gene("g1", "EUK-pEVE", "80"), _gene("g2", "EUK-pEVE", "70")]

    assert taxonomy_class_votes(markers, gene_taxonomy_records=genes) == {"VIRAL_UNKNOWN": 4}
    assert consensus_taxonomy_class(markers, gene_taxonomy_records=genes) == "VIRAL_UNKNOWN"


def test_conflicting_lineage_votes_stay_while_peve_only_genes_drop_out() -> None:
    markers = [
        _marker("NCLDV,PPV", "80,70", porf_id="g1"),
        _marker("NCLDV", "80", porf_id="g2"),
    ]
    genes = [
        _gene("g1", "NCLDV,PPV", "80,70"),
        _gene("g2", "NCLDV", "80"),
        _gene("g3", "EUK-pEVE", "70"),
        _gene("g4", "NCLDV", "60"),
    ]

    assert taxonomy_class_votes(markers, gene_taxonomy_records=genes) == {"VIRAL_UNKNOWN": 3, "NCLDV": 4}
    assert consensus_taxonomy_class(markers, gene_taxonomy_records=genes) == "NCLDV"


def test_conflicting_lineage_votes_keep_the_peve_only_votes_beside_them() -> None:
    markers = [_marker("NCLDV,PPV", "80,70", porf_id="g1")]
    genes = [_gene("g1", "NCLDV,PPV", "80,70"), _gene("g2", "EUK-pEVE", "70")]

    assert taxonomy_class_votes(markers, gene_taxonomy_records=genes) == {"VIRAL_UNKNOWN": 4}
    assert consensus_taxonomy_class(markers, gene_taxonomy_records=genes) == "VIRAL_UNKNOWN"


def test_region_with_only_peve_viral_evidence_takes_the_mixed_gate_path() -> None:
    peve_genes = [_gene("g1", "EUK-pEVE", "80"), _gene("g2", "EUK-pEVE", "70")]
    peve_only = VerificationResult(
        eve_id="e",
        scaffold="s",
        start=0,
        end=3000,
        taxonomy_class="VIRAL_UNKNOWN",
        gene_taxonomy_records=peve_genes,
    )
    mixed_votes = VerificationResult(
        eve_id="e",
        scaffold="s",
        start=0,
        end=3000,
        taxonomy_class="VIRAL_UNKNOWN",
        gene_taxonomy_records=[_gene("g1", "NCLDV,PPV", "80,70")],
    )
    flank_only = VerificationResult(
        eve_id="e",
        scaffold="s",
        start=0,
        end=3000,
        taxonomy_class="VIRAL_UNKNOWN",
        gene_taxonomy_records=[{**peve_genes[0], "is_flanking": True}],
    )
    not_viral = VerificationResult(
        eve_id="e",
        scaffold="s",
        start=0,
        end=3000,
        taxonomy_class="UNKNOWN",
        gene_taxonomy_records=peve_genes,
    )

    assert infer_likely_family(peve_only) == "MIXED"
    assert infer_likely_family(mixed_votes) == "UNKNOWN"
    assert infer_likely_family(flank_only) == "UNKNOWN"
    assert infer_likely_family(not_viral) == "UNKNOWN"


def _verify_markers_and_genes(markers: list[dict], genes: list[dict]) -> VerificationResult:
    """Run the evidence steps that decide the class and the family of a 4 kb region."""
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    result = VerificationResult(eve_id="e", scaffold="sc", start=0, end=4000, confidence_tier="MEDIUM")
    synthesizer._process_hallmark_hits(result, markers)
    synthesizer._process_gene_taxonomy(result, genes, None)
    synthesizer._assign_taxonomy_class(result, markers, genes)
    result.likely_family = infer_likely_family(result)
    return result


def test_marker_database_peve_evidence_alone_takes_the_mixed_gate_path() -> None:
    markers = [_marker("EUK-pEVE", "60", gene="EVE_CAPSID_NOVEL_0002", porf_id="sc|3")]
    genes = [_gene("sc|3", "EUK", "55")]

    result = _verify_markers_and_genes(markers, genes)
    decision = evaluate_v2_quality_gate(result)

    assert result.taxonomy_class == "VIRAL_UNKNOWN"
    assert result.has_capsid_evidence is True
    assert result.likely_family == "MIXED"
    assert (decision.kept, decision.reason) == (True, "mixed_high_medium_pass")


def test_marker_and_gene_database_peve_evidence_together_take_the_mixed_gate_path() -> None:
    markers = [_marker("EUK-pEVE", "60", gene="EVE_CAPSID_NOVEL_0002", porf_id="sc|3")]
    genes = [_gene("sc|3", "EUK", "55"), _gene("sc|5", "EUK-pEVE", "50")]

    assert _verify_markers_and_genes(markers, genes).likely_family == "MIXED"


def test_peve_marker_beside_conflicting_lineage_markers_stays_without_a_family() -> None:
    markers = [
        _marker("NCLDV,PPV", "80,70", gene="EVE_CAPSID_NOVEL_0008", porf_id="sc|1"),
        _marker("EUK-pEVE", "60", gene="EVE_CAPSID_NOVEL_0002", porf_id="sc|3"),
    ]
    genes = [_gene("sc|1", "EUK", "30"), _gene("sc|3", "EUK", "55")]

    result = _verify_markers_and_genes(markers, genes)
    decision = evaluate_v2_quality_gate(result)

    assert result.taxonomy_class == "VIRAL_UNKNOWN"
    assert result.likely_family == "UNKNOWN"
    assert (decision.kept, decision.reason) == (False, "unsupported_class")


def test_peve_marker_below_the_identity_cutoff_gives_no_family() -> None:
    markers = [_marker("EUK-pEVE", "20", gene="EVE_CAPSID_NOVEL_0008", porf_id="sc|3")]
    genes = [_gene("sc|3", "EUK", "55")]

    result = _verify_markers_and_genes(markers, genes)

    assert result.likely_family == "UNKNOWN"


# Family-associated major capsid HMMs: the Phase 1 region class names a lineage
# for each, so the family rules would accept a short region without the size gate.
FAMILY_MCP_MODELS = [
    ("GVOGm0003", "NCLDV"),
    ("mirus_mcp_1", "MIRUS"),
    ("plv_mcp_1", "PPV"),
    ("vp_mcp_1", "PPV"),
]


def _verify_one_marker_region(
    monkeypatch: pytest.MonkeyPatch,
    *,
    gene: str,
    region_class: str,
    end: int,
    prefixes: str = "EUK-pEVE",
    pidents: str = "80",
) -> VerificationResult:
    """Run verify_eve on one validated marker whose only reference is given."""
    synthesizer = EvidenceSynthesizer(config=EvidenceSynthesizerConfig())
    monkeypatch.setattr(synthesizer, "_run_tiebreakers", lambda *args, **kwargs: None)
    hit = {**_marker(prefixes, pidents, gene=gene, porf_id="sc|1"), "hmm_score": 100.0}
    boundary = RefinedBoundary(scaffold="sc", start=0, end=end, confidence=0.95, predicted_family=region_class)
    return synthesizer.verify_eve(
        boundary,
        window_features=[],
        hallmark_hits=[hit],
        gene_taxonomy_records=[_gene("sc|1", prefixes, pidents)],
        gene_taxonomy_summary={"total": 1},
    )


@pytest.mark.parametrize(("gene", "region_class"), FAMILY_MCP_MODELS)
def test_peve_only_region_is_mixed_before_any_marker_or_region_family_rule(
    monkeypatch: pytest.MonkeyPatch, gene: str, region_class: str
) -> None:
    by_marker = _verify_one_marker_region(monkeypatch, gene=gene, region_class="", end=1000)
    by_region = _verify_one_marker_region(monkeypatch, gene=gene, region_class=region_class, end=1000)

    assert by_marker.taxonomy_class == by_region.taxonomy_class == "VIRAL_UNKNOWN"
    assert by_marker.likely_family == "MIXED"
    assert by_region.likely_family == "MIXED"


@pytest.mark.parametrize(("gene", "region_class"), FAMILY_MCP_MODELS)
def test_peve_only_region_with_a_family_marker_takes_the_mixed_size_gate_at_one_kb(
    monkeypatch: pytest.MonkeyPatch, gene: str, region_class: str
) -> None:
    result = _verify_one_marker_region(monkeypatch, gene=gene, region_class=region_class, end=1000)
    decision = evaluate_v2_quality_gate(result)

    assert (decision.kept, decision.effective_class, decision.reason) == (False, "MIXED", "mixed_high_medium_gate")


@pytest.mark.parametrize(("gene", "region_class"), FAMILY_MCP_MODELS)
def test_peve_only_region_with_a_family_marker_passes_the_mixed_gate_at_three_kb(
    monkeypatch: pytest.MonkeyPatch, gene: str, region_class: str
) -> None:
    result = _verify_one_marker_region(monkeypatch, gene=gene, region_class=region_class, end=3000)
    decision = evaluate_v2_quality_gate(result)

    assert (decision.kept, decision.effective_class, decision.reason) == (True, "MIXED", "mixed_high_medium_pass")


@pytest.mark.parametrize(("gene", "region_class"), FAMILY_MCP_MODELS)
def test_family_marker_region_with_a_lineage_hit_keeps_its_family_logic(
    monkeypatch: pytest.MonkeyPatch, gene: str, region_class: str
) -> None:
    result = _verify_one_marker_region(
        monkeypatch,
        gene=gene,
        region_class=region_class,
        end=1000,
        prefixes=f"EUK-pEVE,{region_class}",
        pidents="80,70",
    )

    assert result.taxonomy_class == region_class
    assert result.likely_family == region_class


# --- A1 seed surface ---


def _a1_anchor(
    gene: str,
    *,
    porf_id: str = "gene-a",
    status: str = "validated",
    prefixes: str = "EUK-pEVE__",
    pidents: str = "80",
) -> Anchor:
    return Anchor(
        porf_id=porf_id,
        scaffold="ctg",
        start=10,
        end=40,
        strand="+",
        hallmark_gene=gene,
        score=50.0,
        evalue=1e-20,
        validation_status=status,
        top10_prefixes=prefixes,
        top10_pidents=pidents,
    )


def _a1_result(region_class: str, *anchors: Anchor) -> VerificationResult:
    seed = MergedSeed(
        scaffold="ctg",
        start=0,
        end=3000,
        seed_id="seed",
        sources=["marker_validation"],
        hhg_score=0.0,
        novelty_score=0.0,
        compositional_score=0.0,
        max_kfd=0.0,
        gc_deviation=0.0,
        cub_deviation=0.0,
        confidence="medium",
        predicted_family=region_class,
        region_classification_ncldv_markers=0,
        region_classification_vp_plv_markers=0,
        region_classification_mirus_markers=0,
        anchors=list(anchors),
        hhg_anchors=[],
    )
    return build_phase1_seed_surface([seed]).results[0]


@pytest.mark.parametrize(("gene", "region_class"), FAMILY_MCP_MODELS)
def test_a1_publishes_viral_unknown_and_mixed_for_peve_only_evidence(gene: str, region_class: str) -> None:
    result = _a1_result(region_class, _a1_anchor(gene))

    assert result.taxonomy_class == "VIRAL_UNKNOWN"
    assert result.likely_family == "MIXED"
    assert result.region_classification == region_class


@pytest.mark.parametrize(("gene", "region_class"), FAMILY_MCP_MODELS)
def test_a1_keeps_the_family_when_a_lineage_hit_joins_the_peve_hit(gene: str, region_class: str) -> None:
    result = _a1_result(region_class, _a1_anchor(gene, prefixes=f"EUK-pEVE__,{region_class}__", pidents="80,70"))

    assert result.taxonomy_class == region_class
    assert result.likely_family == region_class


def test_a1_keeps_the_family_when_another_marker_has_the_lineage_hit() -> None:
    peve_only = _a1_anchor("GVOGm0003", porf_id="gene-a")
    lineage = _a1_anchor("GVOGm0022", porf_id="gene-b", prefixes="NCLDV__", pidents="60")

    result = _a1_result("NCLDV", peve_only, lineage)

    assert (result.taxonomy_class, result.likely_family) == ("NCLDV", "NCLDV")


@pytest.mark.parametrize(
    ("status", "prefixes", "pidents"),
    [
        ("validated_novel", "", ""),
        ("validated", "EUK-pEVE__", "20"),
        ("supported", "EUK-pEVE__", "80"),
        ("unvalidated", "EUK-pEVE__", "80"),
        ("validated", "", ""),
    ],
)
def test_a1_keeps_the_family_without_qualified_peve_evidence_from_validated_markers(
    status: str, prefixes: str, pidents: str
) -> None:
    result = _a1_result("NCLDV", _a1_anchor("GVOGm0003", status=status, prefixes=prefixes, pidents=pidents))

    assert (result.taxonomy_class, result.likely_family) == ("NCLDV", "NCLDV")


def test_a1_ignores_an_unvalidated_lineage_hit_beside_validated_peve_evidence() -> None:
    peve_only = _a1_anchor("GVOGm0003", porf_id="gene-a")
    unvalidated_lineage = _a1_anchor(
        "GVOGm0022", porf_id="gene-b", status="supported", prefixes="NCLDV__", pidents="60"
    )

    result = _a1_result("NCLDV", peve_only, unvalidated_lineage)

    assert (result.taxonomy_class, result.likely_family) == ("VIRAL_UNKNOWN", "MIXED")


# --- Detailed output ---


def test_taxonomy_partition_has_a_peve_bucket_after_phage() -> None:
    assert DETAILED_TAXONOMY_PARTITION[-2:] == ("PHAGE", "EUK-pEVE")


def test_partition_bucket_counts_peve_genes_in_their_own_bucket() -> None:
    peve = {
        "top1_target": PEVE_TARGET,
        "top1_prefix": "EUK-pEVE",
        "top10_prefixes": ["EUK-pEVE"],
        "top10_pidents": [80.0],
        "top10_targets": [PEVE_TARGET],
    }
    lineage_first = {**peve, "top10_prefixes": ["PPV", "EUK-pEVE"], "top10_pidents": [50.0, 80.0]}
    weak = {**peve, "top10_pidents": [20.0]}

    assert _taxonomy_partition_bucket(peve) == "EUK-pEVE"
    assert _taxonomy_partition_bucket(lineage_first) == "PPV"
    assert _taxonomy_partition_bucket(weak) == "UNK"


def test_detailed_tsv_counts_peve_genes_and_passes_the_invariants(tmp_path: Path) -> None:
    output_dir = tmp_path / "phase3_synthesis"
    output_dir.mkdir()
    result = VerificationResult(
        eve_id="EVE_scaf1_100-400",
        scaffold="scaf1",
        start=100,
        end=400,
        confidence_tier="HIGH",
        final_confidence=0.9,
        likely_family="MIXED",
        taxonomy_class="VIRAL_UNKNOWN",
        gene_count=3,
        gene_taxonomy_total=3,
        gene_taxonomy_records=[
            {
                "top1_target": PEVE_TARGET,
                "top1_prefix": "EUK-pEVE",
                "top10_prefixes": ["EUK-pEVE"],
                "top10_pidents": [80.0],
                "is_flanking": False,
            },
            {
                "top1_target": PEVE_TARGET,
                "top1_prefix": "EUK-pEVE",
                "top10_prefixes": ["EUK-pEVE"],
                "top10_pidents": [70.0],
                "is_flanking": False,
            },
            {
                "top1_target": "NCLDV__hit",
                "top1_prefix": "NCLDV",
                "top10_prefixes": ["NCLDV"],
                "top10_pidents": [50.0],
                "is_flanking": False,
            },
        ],
    )

    detailed_tsv = OutputGenerator(output_dir=output_dir, extended_output=True).write_predictions_detailed_tsv([result])

    with detailed_tsv.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert "EUK-pEVE:2" in row["taxonomy_best_hits"]
    assert row["taxonomy_best_hits"].endswith("EUK-pEVE:2")
    assert row["effective_eve_class"] == "VIRAL_UNKNOWN"
    assert row["ppv_top10_proteins"] == "0"
    report = run_tsv_invariant_checks(detailed_tsv)
    assert report.passed, report.issues


# --- Other viral vocabularies ---


def test_boundary_summary_counts_a_peve_gene_as_viral_but_not_host_or_family() -> None:
    taxonomy = classify_all_porfs(
        ["g1"],
        {"g1": [DiamondHit("g1", PEVE_TARGET, 1e-30, 200.0, 80.0, 90.0)]},
        {"ctg": [pORF("g1", "ctg", 0, 300)]},
        host_prefix="EUK__",
    )["g1"]
    boundary = RefinedBoundary(
        scaffold="ctg",
        start=0,
        end=300,
        seed_id="seed",
        original_start=0,
        original_end=300,
        confidence=0.5,
        posterior_probability=0.5,
    )
    config = SimpleNamespace(host=SimpleNamespace(label="EUK"), phase2=SimpleNamespace(diamond_flank_genes=0))

    _records, summary = _summarize_boundary_taxonomy(
        boundary=boundary,
        filtered_taxonomy=[taxonomy],
        boundary_taxonomy_map={},
        proteome_index={},
        boundary_diamond_query=None,
        config=config,
        logger=SimpleNamespace(info=lambda *args, **kwargs: None),
    )

    assert summary["viral_top10"] == 1
    assert summary["high_pident_euk"] == 0
    assert summary["ncldv_mirus"] == summary["vp_plv"] == 0
    assert summary["dominant_family"] == "UNKNOWN"
    assert summary["dominant_fraction"] == 0.0


def test_peve_is_a_validated_prefix_and_frameshift_rescue_confirms_it(tmp_path: Path) -> None:
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text(f"q1\t{PEVE_TARGET}\t1e-40\t200\t80.0\t95.0\nq2\t{PEVE_TARGET}\t1e-40\t200\t20.0\t95.0\n")

    coverages = diamond_query_coverages(diamond, min_pident=25.0, validated_prefixes=VALIDATED_PREFIXES)

    assert coverages == {"q1": 95.0}


def test_phylogenetic_diamond_check_reads_peve_as_viral(tmp_path: Path) -> None:
    validator = PhylogeneticValidator(genome_path=tmp_path / "genome.fna", work_dir=tmp_path)
    hits = tmp_path / "hits.tsv"
    row = ["q1", PEVE_TARGET, *["0"] * 9, "200"]
    hits.write_text("\t".join(row) + "\n")

    diamond = validator._parse_diamond_output("eve", hits, 1)

    assert diamond.best_domain == "EUK-pEVE"
    assert diamond.viral_protein_count == 1
    assert diamond.non_viral_protein_count == 0
    assert diamond.is_viral is True
    assert diamond.has_viral_in_neighbors is True


def test_phylogenetic_gvclass_domain_peve_is_viral_never_non_viral() -> None:
    gvclass = GVClassValidation(
        eve_id="eve",
        domain="EUK-pEVE",
        domain_percent=90.0,
        mcp_count=0,
        mirus_count=0,
        gvog_count=0,
        total_markers=0,
    )

    assert gvclass.is_viral is True
    assert gvclass.is_non_viral is False


def test_structural_hits_and_default_config_read_peve_as_viral() -> None:
    assert "EUK-pEVE__" in get_config().database.viral_prefixes
    assert BoltzFoldSeekAnalyzer._is_viral_target(object(), "EUK-pEVE__EP00022|p1") is True
    assert BoltzFoldSeekAnalyzer._is_viral_target(object(), "EUK__EP00022|p1") is False
