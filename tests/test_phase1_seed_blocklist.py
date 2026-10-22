"""Kinase marker models never seed a region but still count inside one."""

from pathlib import Path

import pytest

from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit
from virosync.pipeline.phase1.region_assembly import assemble_candidate_regions
from virosync.pipeline.phase1.viral_markers import SEED_BLOCKED_MARKER_MODELS, is_seed_blocked_model

BLOCKED_MODEL = "VS000087"
SEED_MODEL = "GVOGm0003"


def _marker(porf: str, scaffold: str, start: int, end: int, model: str) -> ValidatedMarkerHit:
    return ValidatedMarkerHit(
        query_porf=porf,
        scaffold=scaffold,
        start=start,
        end=end,
        strand="+",
        hmm_target=model,
        hmm_score=120.0,
        hmm_evalue=1e-30,
        validation_status="validated",
        top10_prefixes="NCLDV__",
        best_hit_target="NCLDV__virus|p1",
        best_hit_pident=40.0,
        best_hit_bits=100.0,
        has_ncldv=1,
        has_mirus=0,
        has_plv=0,
        has_vp=0,
        has_viral=1,
    )


def _inputs(tmp_path: Path) -> tuple[Path, Path]:
    genome = tmp_path / "genome.fna"
    contigs = ("ctgA", "ctgB", "ctgC", "ctgD", "ctgE")
    genome.write_text("".join(f">{name}\n{'ACGT' * 5000}\n" for name in contigs), encoding="utf-8")
    proteome = tmp_path / "proteome.faa"
    genes = [f">ctgE_{n} # {1000 + n * 500} # {1300 + n * 500} # 1 # ID=1_{n};partial=00\nMKK\n" for n in range(1, 10)]
    proteome.write_text("".join(genes), encoding="utf-8")
    return genome, proteome


def _assemble(tmp_path: Path, markers: list[ValidatedMarkerHit]):
    genome, proteome = _inputs(tmp_path)
    return assemble_candidate_regions(
        markers,
        genome,
        proteome,
        tmp_path / "out",
        initial_window_bp=100,
        initial_window_genes=1,
        min_markers_initial=1,
        extension_kb=5,
        merge_distance=0,
    )


@pytest.mark.parametrize("model", ["VS000084", "VS000087", "VS000094", "VS000191", "VS000319"])
def test_seed_blocklist_names_the_kinase_models(model: str) -> None:
    assert model in SEED_BLOCKED_MARKER_MODELS


def test_blocked_marker_alone_does_not_seed_a_region(tmp_path: Path) -> None:
    regions = _assemble(tmp_path, [_marker("ctgA_1", "ctgA", 1000, 1300, BLOCKED_MODEL)])

    assert regions == []


def test_two_blocked_markers_do_not_seed_a_region(tmp_path: Path) -> None:
    markers = [
        _marker("ctgB_1", "ctgB", 1000, 1300, BLOCKED_MODEL),
        _marker("ctgB_2", "ctgB", 1500, 1800, "VS000084"),
    ]

    regions = _assemble(tmp_path, markers)

    assert regions == []


def test_other_marker_alone_still_seeds_a_region(tmp_path: Path) -> None:
    regions = _assemble(tmp_path, [_marker("ctgC_1", "ctgC", 1000, 1300, SEED_MODEL)])

    assert len(regions) == 1


def test_blocked_marker_next_to_a_seed_marker_counts_in_the_region(tmp_path: Path) -> None:
    markers = [
        _marker("ctgD_1", "ctgD", 1000, 1300, SEED_MODEL),
        _marker("ctgD_2", "ctgD", 1500, 1800, BLOCKED_MODEL),
    ]

    regions = _assemble(tmp_path, markers)

    assert [sorted(m.hmm_target for m in region.markers) for region in regions] == [[SEED_MODEL, BLOCKED_MODEL]]


def test_blocked_marker_inside_the_extension_of_a_seed_counts_in_the_region(tmp_path: Path) -> None:
    markers = [
        _marker("ctgE_1", "ctgE", 1500, 1800, SEED_MODEL),
        _marker("ctgE_9", "ctgE", 5500, 5800, BLOCKED_MODEL),
    ]

    regions = _assemble(tmp_path, markers)

    assert [sorted(m.hmm_target for m in region.markers) for region in regions] == [[SEED_MODEL, BLOCKED_MODEL]]


@pytest.mark.parametrize(
    "name",
    ["VS000087", "vs000087", "VS000087_Pkinase", "vs000087_Pkinase", "VS000087.hmm"],
)
def test_blocklist_matches_suffixed_and_lower_case_model_names(name: str) -> None:
    assert is_seed_blocked_model(name)


@pytest.mark.parametrize("name", ["GVOGm0003", "VS0000871", "VS000088_Pkinase", ""])
def test_blocklist_ignores_other_model_names(name: str) -> None:
    assert not is_seed_blocked_model(name)


@pytest.mark.parametrize("name", ["VS000087_Pkinase", "vs000087_Pkinase", "VS000087.hmm"])
def test_suffixed_blocked_marker_does_not_seed_a_region(tmp_path: Path, name: str) -> None:
    regions = _assemble(tmp_path, [_marker("ctgA_1", "ctgA", 1000, 1300, name)])

    assert regions == []
