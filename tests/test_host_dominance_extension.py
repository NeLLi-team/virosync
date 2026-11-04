"""Three-host dominance preserves positive support, confirmed cores, and audit evidence."""

from dataclasses import replace
from types import SimpleNamespace

import pytest
from test_joint_seed_rescue import child, rejected
from test_joint_seed_rescue import pair_inputs as pair_inputs

from virosync.orchestration._flows.single_genome.phase2_resume_state import _gene_from_document, _gene_to_document
from virosync.pipeline.phase1.seed_merger import ConfirmedSeed, MergedSeed
from virosync.pipeline.phase2.boundary_diamond import (
    BoundaryDiamondConfig,
    DiamondHit,
    GeneTaxonomy,
    GenomeDiamondQuery,
    SeedGeneMapping,
    build_gene_taxonomy_record,
    classify_cached_diamond_query,
    pORF,
)
from virosync.pipeline.phase2.seed_extension import extend_confirmed_seed, merge_reachable_boundaries


def classified(name, start, end, *, viral=80.0, shallow_viral=False, nohit=False):
    """Use separate actual cached search lanes, never infer deep rank from shallow."""
    gene = pORF(name, "ctg", start, end)
    cell = DiamondHit(name, "EUK__host|p", 1e-20, 100.0, 70.0, 90.0)
    virus = DiamondHit(name, "NCLDV__virus|p", 1e-20, viral, 30.0, 90.0)
    shallow = [] if nohit else (sorted([virus, cell], key=lambda hit: -hit.bits) if shallow_viral else [cell])
    deep = [] if nohit else sorted([virus, cell], key=lambda hit: -hit.bits)
    tax = classify_cached_diamond_query(
        GenomeDiamondQuery(all_porf_ids=[name]),
        {name: shallow},
        {"ctg": [gene]},
        BoundaryDiamondConfig(viral_hit_depth=50),
        viral_veto_hits={name: deep},
    )[name]
    return gene, tax


def grow(
    monkeypatch,
    *,
    direction=1,
    viral=80.0,
    model=True,
    middle="hit",
    disabled=False,
    protected=(),
    marker=False,
    shallow_hosts=False,
):
    """Keep native classification/extension; isolate the configured host-share input."""
    monkeypatch.setattr(
        "virosync.pipeline.phase2.seed_extension.host_lineage_share", lambda tax, model: 0.6 if tax.has_hit else 0.0
    )
    coordinates = [(2100, 2200), (2300, 2400), (2500, 2600), (2700, 2800), (2900, 3100)]
    if direction == -1:
        coordinates = [(10000 - end, 10000 - start) for start, end in coordinates]
    records = [
        classified(
            name, start, end, viral=150.0 if index in (0, 4) else viral, shallow_viral=shallow_hosts or index in (0, 4)
        )
        for index, (name, (start, end)) in enumerate(
            zip(["near", "h1", "h2", "h3", "remote"], coordinates, strict=True)
        )
    ]
    genes = [gene for gene, _ in records]
    taxonomy = {gene.id: tax for gene, tax in records}
    if marker:
        taxonomy["h2"].phase1_validated_viral_marker = True
    if middle == "nohit":
        _, taxonomy["h2"] = classified("h2", genes[2].start, genes[2].end, nohit=True)
    elif middle == "missing":
        taxonomy.pop("h2")
    start, end = (1000, 2000) if direction == 1 else (8000, 9000)
    seed = MergedSeed("ctg", start, end, seed_id="core", confirmed_seeds=(ConfirmedSeed("core", "ctg", start, end),))
    host_model = SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0}) if model else None
    protection = {"protected_seed_intervals": protected} if protected else {}
    boundary = extend_confirmed_seed(
        seed, genes, taxonomy, (0, 10000), host_model, 0.35, disable_host_stop=disabled, **protection
    )
    return boundary, genes, taxonomy


@pytest.mark.parametrize("direction", [1, -1])
def test_three_cellular_dominant_weak_viral_genes_stop_both_directions(monkeypatch, direction):
    boundary, _, _ = grow(monkeypatch, direction=direction)
    assert (boundary.start, boundary.end) == ((1000, 2200) if direction == 1 else (7800, 9000))
    assert (boundary.right_stop if direction == 1 else boundary.left_stop) == (2300 if direction == 1 else 7700)
    assert (boundary.right_stop_reason if direction == 1 else boundary.left_stop_reason) == "host_stop"
    assert boundary.confirmed_seeds[0].seed_id == "core"


@pytest.mark.parametrize("viral,model", [(100.0, True), (150.0, True), (80.0, False)])
def test_tie_stronger_deep_viral_or_no_model_does_not_stop(monkeypatch, viral, model):
    boundary, _, _ = grow(monkeypatch, viral=viral, model=model)
    assert boundary.end == 3100
    assert boundary.right_stop is None


@pytest.mark.parametrize("middle", ["nohit", "missing"])
def test_nohit_and_missing_break_host_run_without_inventing_support(monkeypatch, middle):
    boundary, _, _ = grow(monkeypatch, middle=middle)
    assert boundary.right_stop_reason != "host_stop"
    assert boundary.right_stop is None
    if middle == "missing":
        assert boundary.end <= 2500  # Existing unavailable-taxonomy boundary remains authoritative.
    else:
        assert boundary.end == 3100


def test_a4_disables_host_stop_but_isolated_gene_support_is_unchanged(monkeypatch):
    boundary, _, _ = grow(monkeypatch, disabled=True)
    assert boundary.end == 3100
    assert boundary.right_stop is None
    # The owner approved a three-gene stop, not suppression of isolated support.
    gene, tax = classified("weak", 2100, 2200, viral=80.0, shallow_viral=True)
    seed = MergedSeed("ctg", 1000, 2000, seed_id="core", confirmed_seeds=(ConfirmedSeed("core", "ctg", 1000, 2000),))
    model = SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0})
    normal = extend_confirmed_seed(seed, [gene], {"weak": tax}, (0, 10000), model, 0.35)
    disabled = extend_confirmed_seed(seed, [gene], {"weak": tax}, (0, 10000), model, 0.35, disable_host_stop=True)
    assert normal.end == 2200
    assert disabled.end == 2200


def test_competition_uses_union_and_preserves_original10_fields():
    _, weak = classified("g", 100, 400, viral=80.0)
    _, strong = classified("g", 100, 400, viral=150.0)
    for field in ("top10_targets", "top10_bits", "top10_pidents", "top10_evalues", "top10_qcovhsp"):
        assert getattr(weak, field) == getattr(strong, field)
    assert weak.host_competition_assessed is True
    assert (weak.host_competition_cellular_bits, weak.host_competition_viral_bits) == (100.0, 80.0)
    assert (strong.host_competition_cellular_bits, strong.host_competition_viral_bits) == (100.0, 150.0)
    _, empty = classified("empty", 100, 400, nohit=True)
    assert empty.host_competition_assessed is True
    assert (empty.host_competition_cellular_bits, empty.host_competition_viral_bits) == (None, None)
    unavailable = GeneTaxonomy("missing", "ctg", 100, 400)
    assert unavailable.host_competition_assessed is False


def test_validated_phase1_marker_retains_its_independent_veto(monkeypatch):
    boundary, _, _ = grow(monkeypatch, marker=True)
    assert boundary.right_stop is None
    assert boundary.end == 3100


@pytest.mark.parametrize("protected", [((2250, 2850),), ((2500, 2600),)])
def test_new_barrier_does_not_form_inside_another_confirmed_core(monkeypatch, protected):
    boundary, _, _ = grow(monkeypatch, protected=protected)
    assert boundary.right_stop is None
    assert boundary.end == 3100


def test_typed_checkpoint_and_record_preserve_competition():
    _, tax = classified("g", 100, 400)
    document = _gene_to_document(tax, "gene")
    assert _gene_from_document(document, "gene") == tax
    record = build_gene_taxonomy_record(tax)
    for field in ("host_competition_assessed", "host_competition_cellular_bits", "host_competition_viral_bits"):
        assert record[field] == getattr(tax, field)
    malformed = dict(document, host_competition_assessed="false")
    with pytest.raises(ValueError):
        _gene_from_document(malformed, "gene")


def test_inside_host_genes_and_separate_confirmed_seed_remain(monkeypatch):
    boundary, _, _ = grow(monkeypatch)
    genes_and_tax = [classified(f"inside{i}", 1100 + i * 200, 1200 + i * 200) for i in range(3)]
    seed = MergedSeed(
        "ctg", 1000, 2000, seed_id="inside", confirmed_seeds=(ConfirmedSeed("inside", "ctg", 1000, 2000),)
    )
    model = SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0})
    inside = extend_confirmed_seed(
        seed, [g for g, _ in genes_and_tax], {g.id: t for g, t in genes_and_tax}, (0, 10000), model, 0.35
    )
    assert (inside.start, inside.end) == (1000, 2000)
    separate = replace(
        boundary,
        start=2900,
        end=3100,
        seed_id="other",
        confirmed_seeds=(ConfirmedSeed("other", "ctg", 2900, 3100),),
        right_stop=None,
    )
    merged = merge_reachable_boundaries([boundary, separate])
    assert {core.seed_id for item in merged for core in item.confirmed_seeds} == {"core", "other"}


@pytest.mark.parametrize("viral,expected", [(80.0, True), (100.0, False), (150.0, False)])
def test_joint_gap_uses_same_three_gene_competition(monkeypatch, viral, expected):
    from virosync.pipeline.phase3.joint_seed_rescue import _gap_has_host_barrier

    monkeypatch.setattr("virosync.pipeline.phase3.joint_seed_rescue.host_lineage_share", lambda tax, model: 0.6)
    pairs = [classified(f"h{i}", 2200 + i * 200, 2300 + i * 200, viral=viral) for i in range(3)]
    children = (SimpleNamespace(end=2000), SimpleNamespace(start=3000))
    model = SimpleNamespace(token_weights={"host": 1}, lineage_weights={"euk|host": 1.0})
    assert (
        _gap_has_host_barrier(children, [g for g, _ in pairs], {g.id: t for g, t in pairs}, model, 0.35, frozenset())
        is expected
    )


@pytest.mark.parametrize("prefix", ["PHAGE", "EUK-pEVE", "CRESS", "GVMAG"])
@pytest.mark.parametrize("identity,expected", [(24.9, None), (25.0, 80.0)])
def test_host_competition_keeps_native_viral_namespaces_and_identity(prefix, identity, expected):
    cell = DiamondHit("g", "EUK__host|p", 1e-20, 100.0, 70.0, 90.0)
    virus = DiamondHit("g", prefix + "__virus|p", 1e-20, 80.0, identity, 90.0)
    taxonomy = classify_cached_diamond_query(
        GenomeDiamondQuery(all_porf_ids=["g"]),
        {"g": [cell]},
        {"ctg": [pORF("g", "ctg", 100, 400)]},
        BoundaryDiamondConfig(viral_hit_depth=50),
        viral_veto_hits={"g": [cell, virus]},
    )["g"]
    assert taxonomy.host_competition_viral_bits == expected
    assert taxonomy.host_competition_cellular_target == "EUK__host|p"
    assert taxonomy.viral_hit_within_depth is (expected is not None)


def test_depth_ten_ignores_an_unrequested_deeper_search():
    cell = DiamondHit("g", "EUK__host|p", 1e-20, 100.0, 70.0, 90.0)
    virus = DiamondHit("g", "MIRUS__virus|p", 1e-20, 150.0, 80.0, 90.0)
    taxonomy = classify_cached_diamond_query(
        GenomeDiamondQuery(all_porf_ids=["g"]),
        {"g": [cell]},
        {"ctg": [pORF("g", "ctg", 100, 400)]},
        BoundaryDiamondConfig(viral_hit_depth=10),
        viral_veto_hits={"g": [virus]},
    )["g"]
    assert taxonomy.host_competition_assessed
    assert taxonomy.host_competition_viral_bits is None
    assert taxonomy.host_competition_viral_target == ""


@pytest.mark.parametrize(
    "change",
    [
        {"host_competition_viral_bits": float("nan")},
        {"host_competition_viral_bits": -1.0},
        {"host_competition_assessed": False},
        {"host_competition_viral_target": ""},
        {"viral_hit_within_depth": False},
    ],
)
def test_checkpoint_rejects_inconsistent_competition(change):
    _, taxonomy = classified("g", 100, 400)
    document = _gene_to_document(taxonomy, "gene")
    document.update(change)
    with pytest.raises(ValueError):
        _gene_from_document(document, "gene")


def test_assessed_nohit_and_unavailable_roundtrip_as_distinct_states():
    _, nohit = classified("g", 100, 400, nohit=True)
    unavailable = GeneTaxonomy("g", "ctg", 100, 400)
    assert _gene_from_document(_gene_to_document(nohit, "gene"), "gene") == nohit
    assert _gene_from_document(_gene_to_document(unavailable, "gene"), "gene") == unavailable
    assert nohit != unavailable


def test_top10_positive_host_run_stops_only_when_three_gene_barrier_completes(monkeypatch):
    boundary, _, taxonomy = grow(monkeypatch, shallow_hosts=True)
    for name in ("h1", "h2", "h3"):
        assert taxonomy[name].has_viral
        assert taxonomy[name].host_competition_cellular_bits > taxonomy[name].host_competition_viral_bits
    assert boundary.right_stop_reason == "host_stop"
    assert boundary.right_stop == 2300
    assert boundary.end == 2200
    without_stop, _, _ = grow(monkeypatch, shallow_hosts=True, disabled=True)
    assert without_stop.end == 3100
    # The same positive genes extend when no third host completes the barrier.
    interrupted, _, _ = grow(monkeypatch, shallow_hosts=True, protected=((2500, 2600),))
    assert interrupted.right_stop is None
    assert interrupted.end == 3100


def test_joint_pair_passes_separate_confirmed_core_to_real_gap_guard(pair_inputs, monkeypatch):
    from virosync.pipeline.phase3 import joint_seed_rescue

    third = child("third", 2700, 3500)
    pair_inputs["boundaries"].append(third)
    pair_inputs["results"].append(rejected(third))
    pair_inputs["query"].seed_gene_mappings["third"] = SeedGeneMapping(
        "third", "s", 2700, 3500, flank_start_bp=2000, flank_end_bp=4000
    )
    native_gap_guard = joint_seed_rescue._gap_has_host_barrier
    checked = []

    def observe_gap_guard(children, *args, **kwargs):
        checked.append((tuple(child.seed_id for child in children), kwargs["protected_seed_intervals"]))
        return native_gap_guard(children, *args, **kwargs)

    monkeypatch.setattr(joint_seed_rescue, "_gap_has_host_barrier", observe_gap_guard)
    pairs = joint_seed_rescue.eligible_joint_pairs(**pair_inputs)
    assert len(pairs) == 2
    assert checked == [
        (("left", "right"), ((100, 1000), (1500, 2500), (2700, 3500))),
        (("right", "third"), ((100, 1000), (1500, 2500), (2700, 3500))),
    ]
