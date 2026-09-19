"""Tests for bounded integration-gene evidence collection."""

from __future__ import annotations

import random
from dataclasses import asdict
from pathlib import Path

import pyhmmer
import pytest

from virosync.pipeline.phase1.marker_validation import ValidatedMarkerHit
from virosync.pipeline.phase3 import integration_genes
from virosync.pipeline.phase3.integration_genes import (
    INTEGRATION_PROFILE_PATH,
    IntegrationRegion,
    collect_annotated_integration_genes,
    scan_integration_genes,
)


def _profile_consensuses() -> dict[str, str]:
    with pyhmmer.plan7.HMMFile(INTEGRATION_PROFILE_PATH) as handle:
        return {profile.name.decode(): profile.consensus for profile in handle}


def _random_protein(length: int, *, seed: int) -> str:
    generator = random.Random(seed)
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(generator.choice(amino_acids) for _position in range(length))


def _shuffled(sequence: str, *, seed: int) -> str:
    residues = list(sequence)
    random.Random(seed).shuffle(residues)
    return "".join(residues)


def _write_proteome(path: Path, records: list[tuple[str, int, str, str]]) -> None:
    lines: list[str] = []
    for gene_index, (protein_id, start, strand, sequence) in enumerate(records, start=1):
        end = start + len(sequence) * 3
        strand_token = "1" if strand == "+" else "-1"
        lines.extend(
            (
                f">{protein_id} # {start + 1} # {end} # {strand_token} # ID=1_{gene_index};partial=00",
                sequence,
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _marker(
    protein_id: str,
    model: str,
    start: int,
    end: int,
    *,
    score: float,
    status: str,
) -> ValidatedMarkerHit:
    return ValidatedMarkerHit(
        query_porf=f"{protein_id}|aa1-20",
        scaffold="ctg",
        start=start,
        end=end,
        strand="+",
        hmm_target=model,
        hmm_score=score,
        hmm_evalue=1e-20,
        validation_status=status,
        top10_prefixes="",
        best_hit_target="",
        best_hit_pident=0.0,
        best_hit_bits=0.0,
        has_ncldv=0,
        has_mirus=0,
        has_plv=0,
        has_vp=0,
        has_viral=0,
    )


def test_profile_screen_batches_selected_genes_with_profile_ga_cutoffs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    consensuses = _profile_consensuses()
    proteome_path = tmp_path / "proteome.faa"
    _write_proteome(
        proteome_path,
        [
            ("ctg_1", 100, "+", consensuses["Phage_integrase"]),
            ("ctg_2", 3000, "-", consensuses["Resolvase"]),
            ("ctg_3", 1900, "+", consensuses["rve"]),
            ("ctg_4", 5000, "+", consensuses["rve_2"]),
            ("ctg_5", 8500, "-", consensuses["rve_3"]),
            ("ctg_6", 6000, "+", _random_protein(220, seed=7)),
            ("ctg_7", 6800, "+", _shuffled(consensuses["Phage_integrase"], seed=11)),
            ("ctg_8", 10_500, "+", consensuses["Phage_integrase"]),
        ],
    )
    region = IntegrationRegion(
        region_id="eve-1",
        scaffold="ctg",
        scan_start=0,
        scan_end=10_000,
        interior_start=2000,
        interior_end=8000,
    )
    real_hmmsearch = pyhmmer.hmmsearch
    hmmsearch_calls = 0

    def counting_hmmsearch(*args: object, **kwargs: object) -> object:
        nonlocal hmmsearch_calls
        hmmsearch_calls += 1
        assert "Z" not in kwargs
        return real_hmmsearch(*args, **kwargs)

    monkeypatch.setattr(integration_genes.pyhmmer, "hmmsearch", counting_hmmsearch)

    scan = scan_integration_genes(proteome_path, [region], threads=1)
    hits = scan.hits

    by_protein_profile = {(hit.protein_id, hit.profile): hit for hit in hits}
    assert hmmsearch_calls == 1
    assert scan.unsearched == []
    assert set(by_protein_profile) == {
        ("ctg_1", "Phage_integrase"),
        ("ctg_2", "Resolvase"),
        ("ctg_3", "rve"),
        ("ctg_4", "rve_2"),
        ("ctg_4", "rve_3"),
        ("ctg_5", "rve_2"),
        ("ctg_5", "rve_3"),
    }
    assert by_protein_profile[("ctg_1", "Phage_integrase")].start == 100
    assert by_protein_profile[("ctg_1", "Phage_integrase")].end == 616
    assert by_protein_profile[("ctg_1", "Phage_integrase")].location == "upstream"
    assert by_protein_profile[("ctg_2", "Resolvase")].strand == "-"
    assert by_protein_profile[("ctg_2", "Resolvase")].location == "interior"
    assert by_protein_profile[("ctg_3", "rve")].location == "boundary_overlap"
    assert by_protein_profile[("ctg_5", "rve_3")].location == "downstream"
    assert by_protein_profile[("ctg_1", "Phage_integrase")].sequence_score > 27.1
    assert by_protein_profile[("ctg_1", "Phage_integrase")].domain_score > 27.1
    assert {hit.mechanism for hit in hits} == {
        "tyrosine_recombinase",
        "serine_recombinase",
        "dde_integrase",
    }


def test_mixed_screen_preserves_hits_and_records_every_unsearched_context(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Oversized proteins retain context while sequence E-values use all targets."""
    proteome = tmp_path / "proteome.faa"
    integrase = _profile_consensuses()["Phage_integrase"]
    region = IntegrationRegion("eve-1", "ctg", 0, 400_000, 0, 400_000)
    other_region = IntegrationRegion("eve-2", "ctg", 500, 500_000, 400_000, 500_000)
    _write_proteome(proteome, [("ctg_1", 100, "+", integrase), ("ctg_2", 1_000, "-", "M" * 100_001)])
    original_bytes = proteome.read_bytes()
    real_hmmsearch = pyhmmer.hmmsearch
    search_options: list[dict[str, object]] = []

    def recording_search(profiles: object, sequences: object, **kwargs: object) -> object:
        search_options.append(kwargs)
        return real_hmmsearch(profiles, sequences, **kwargs)

    monkeypatch.setattr(integration_genes.pyhmmer, "hmmsearch", recording_search)
    scan = scan_integration_genes(proteome, [region, other_region], threads=1)

    assert search_options == [
        {
            "cpus": 1,
            "parallel": "targets",
            "E": float("inf"),
            "domE": float("inf"),
            "incE": float("inf"),
            "incdomE": float("inf"),
            "Z": 2,
        }
    ]
    supported_proteome = tmp_path / "supported.faa"
    _write_proteome(supported_proteome, [("ctg_1", 100, "+", integrase)])
    baseline = scan_integration_genes(supported_proteome, [region], threads=1).hits[0]
    hit = scan.hits[0]
    assert hit.evalue == pytest.approx(baseline.evalue * 2, abs=0.0)
    assert asdict(hit) | {"evalue": baseline.evalue} == asdict(baseline)
    assert [asdict(record) for record in scan.unsearched] == [
        {
            "region_id": "eve-1",
            "protein_id": "ctg_2",
            "length_aa": 100_001,
            "scaffold": "ctg",
            "start": 1_000,
            "end": 301_003,
            "strand": "-",
            "location": "interior",
            "reason": "sequence_length_limit",
            "limit_aa": 100_000,
        },
        {
            "region_id": "eve-2",
            "protein_id": "ctg_2",
            "length_aa": 100_001,
            "scaffold": "ctg",
            "start": 1_000,
            "end": 301_003,
            "strand": "-",
            "location": "upstream",
            "reason": "sequence_length_limit",
            "limit_aa": 100_000,
        },
    ]
    assert proteome.read_bytes() == original_bytes


@pytest.mark.parametrize("length_aa, expected_searches", [(100_000, 1), (100_001, 0)])
def test_screen_length_boundary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, length_aa: int, expected_searches: int
) -> None:
    """The engine limit is inclusive and an entirely excluded batch is not searched."""
    proteome = tmp_path / "proteome.faa"
    _write_proteome(proteome, [("ctg_1", 0, "+", "M" * length_aa)])
    searched: list[bytes] = []

    def empty_search(profiles: object, sequences: object, **kwargs: object) -> object:
        assert "Z" not in kwargs
        searched.extend(sequence.name for sequence in sequences)
        return [[] for _profile in profiles]

    monkeypatch.setattr(integration_genes.pyhmmer, "hmmsearch", empty_search)
    scan = scan_integration_genes(proteome, [IntegrationRegion("eve-1", "ctg", 0, 400_000, 0, 400_000)], threads=1)

    assert len(searched) == expected_searches
    assert scan.hits == []
    assert len(scan.unsearched) == 1 - expected_searches


@pytest.mark.parametrize("regions", [[], [IntegrationRegion("eve-1", "other", 0, 100, 0, 100)]])
def test_screen_empty_selection_has_no_exclusions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, regions: list[IntegrationRegion]
) -> None:
    """No regions and regions without proteins both complete without search."""
    proteome = tmp_path / "proteome.faa"
    _write_proteome(proteome, [("ctg_1", 0, "+", "M" * 30)])

    def unexpected_search(*args: object, **kwargs: object) -> object:
        pytest.fail("empty selections must not invoke HMMER")

    monkeypatch.setattr(integration_genes.pyhmmer, "hmmsearch", unexpected_search)
    scan = scan_integration_genes(proteome, regions, threads=1)

    assert scan.hits == []
    assert scan.unsearched == []


def test_unrelated_search_failure_propagates(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Length handling must not suppress unrelated HMMER failures."""
    proteome = tmp_path / "proteome.faa"
    _write_proteome(proteome, [("ctg_1", 0, "+", "M" * 30)])

    def failed_search(*args: object, **kwargs: object) -> object:
        raise ValueError("unrelated profile failure")

    monkeypatch.setattr(integration_genes.pyhmmer, "hmmsearch", failed_search)

    with pytest.raises(ValueError, match="unrelated profile failure"):
        scan_integration_genes(proteome, [IntegrationRegion("eve-1", "ctg", 0, 100, 0, 100)], threads=1)


def test_missing_profiles_do_not_count_as_a_completed_screen(tmp_path: Path) -> None:
    """A supported selection still requires the bundled HMM resource."""
    proteome = tmp_path / "proteome.faa"
    _write_proteome(proteome, [("ctg_1", 0, "+", "M" * 30)])

    with pytest.raises(FileNotFoundError):
        scan_integration_genes(
            proteome,
            [IntegrationRegion("eve-1", "ctg", 0, 100, 0, 100)],
            threads=1,
            hmm_path=tmp_path / "missing.hmm",
        )


def test_existing_annotations_keep_source_status_and_ignore_domain_counts(tmp_path: Path) -> None:
    proteome_path = tmp_path / "proteome.faa"
    _write_proteome(
        proteome_path,
        [
            ("ctg_1", 0, "+", "M" * 30),
            ("ctg_2", 120, "+", "M" * 30),
            ("ctg_3", 240, "+", "M" * 30),
            ("ctg_4", 360, "+", "M" * 30),
            ("ctg_5", 480, "+", "M" * 30),
            ("ctg_6", 600, "+", "M" * 30),
            ("ctg_7", 720, "-", "M" * 30),
        ],
    )
    annotations_path = tmp_path / "model_annotations_with_interpro.tsv"
    annotations_path.write_text(
        "model_name\tdescription\tmajority_annotation\tpfam_signature\tpfam_top_domains\n"
        "TYR\tfamily protein\tTyrosine recombinase\t\t.\n"
        "SER\tSerine-recombinase catalytic subunit\tunknown\t\t.\n"
        "INT\tintegrase catalytic protein\tunknown\t\t.\n"
        "REC\trecombinase family protein\tunknown\t\t.\n"
        "YQAJ\tphage recombinase\tunknown\tYqaJ\t.\n"
        "FALSE\tnonintegrase family protein\tunknown\t\tPhage_integrase:99\n",
        encoding="utf-8",
    )
    interproscan_path = tmp_path / "interproscan_batch.tsv"
    interproscan_path.write_text(
        "eve-1|ctg_7\tmd5\t30\tPfam\tPF00665\tIntegrase core domain\t3\t27\t3.4E-5\tT\t2026-01-01\t"
        "IPR001584\tIntegrase catalytic core\n"
        "eve-1|ctg_7\tmd5\t30\tPfam\tPF00665\tIntegrase core domain\t2\t28\t1.2E-30\tT\t2026-01-01\t"
        "IPR001584\tIntegrase catalytic core\n",
        encoding="utf-8",
    )
    region = IntegrationRegion(
        region_id="eve-1",
        scaffold="ctg",
        scan_start=0,
        scan_end=900,
        interior_start=100,
        interior_end=700,
    )
    markers = [
        _marker("ctg_1", "TYR", 0, 90, score=71.0, status="validated"),
        _marker("ctg_2", "SER", 120, 210, score=62.0, status="validated_novel"),
        _marker("ctg_3", "INT", 240, 330, score=53.0, status="supported"),
        _marker("ctg_4", "REC", 360, 450, score=44.0, status="unvalidated"),
        _marker("ctg_5", "YQAJ", 480, 570, score=35.0, status="validated"),
        _marker("ctg_6", "FALSE", 600, 690, score=99.0, status="validated"),
    ]

    hits = collect_annotated_integration_genes(
        proteome_path,
        markers,
        [region],
        model_annotations_path=annotations_path,
        interproscan_path=interproscan_path,
    )

    by_source_and_protein = {(hit.source, hit.protein_id): hit for hit in hits}
    assert set(by_source_and_protein) == {
        ("marker_annotation", "ctg_1"),
        ("marker_annotation", "ctg_2"),
        ("marker_annotation", "ctg_3"),
        ("marker_annotation", "ctg_4"),
        ("marker_annotation", "ctg_5"),
        ("interproscan", "ctg_7"),
    }
    assert by_source_and_protein[("marker_annotation", "ctg_1")].mechanism == "tyrosine_recombinase"
    assert by_source_and_protein[("marker_annotation", "ctg_2")].mechanism == "serine_recombinase"
    assert by_source_and_protein[("marker_annotation", "ctg_3")].mechanism == "integrase"
    assert by_source_and_protein[("marker_annotation", "ctg_4")].mechanism == "recombinase"
    assert by_source_and_protein[("marker_annotation", "ctg_5")].mechanism == "repair_recombinase"
    assert by_source_and_protein[("marker_annotation", "ctg_4")].validation_status == "unvalidated"
    assert by_source_and_protein[("interproscan", "ctg_7")].mechanism == "dde_integrase"
    assert by_source_and_protein[("interproscan", "ctg_7")].start == 720
    assert by_source_and_protein[("interproscan", "ctg_7")].end == 810
    assert by_source_and_protein[("interproscan", "ctg_7")].strand == "-"
    assert by_source_and_protein[("interproscan", "ctg_7")].location == "downstream"
    assert by_source_and_protein[("interproscan", "ctg_7")].domain_score is None
    assert by_source_and_protein[("interproscan", "ctg_7")].evalue == 1.2e-30
    assert by_source_and_protein[("interproscan", "ctg_7")].domain_start == 2
    assert by_source_and_protein[("interproscan", "ctg_7")].domain_end == 28
