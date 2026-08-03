"""Phase 3 per-genome ANI clustering and MCP taxonomy-class propagation."""

import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from virosync.pipeline.phase3.eve_ani_clustering import (
    MIN_CLUSTER_ALIGNED_FRACTION,
    MIN_CLUSTER_ANI,
    cluster_accepted_eves,
    recluster_survivors,
    unsupported_eve_ids,
)
from virosync.pipeline.phase3.evidence_synthesizer import VerificationResult

_SKANI_HEADER = (
    "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\n"
)


def _result(eve_id: str, *, taxonomy_class: str, has_mcp: bool) -> VerificationResult:
    """One accepted result. ``has_mcp`` means an MCP marker decided its class.

    Both flags are set together because that is the shape Phase 3 produces for
    an MCP-marker call. ``test_structural_mcp_alone_does_not_donate`` covers the
    case where they diverge.
    """
    return VerificationResult(
        eve_id=eve_id,
        scaffold="ctg1",
        start=0,
        end=600,
        taxonomy_class=taxonomy_class,
        has_mcp=has_mcp,
        taxonomy_class_from_mcp=has_mcp,
    )


def _genome(tmp_path: Path) -> Path:
    genome = tmp_path / "genome.fna"
    genome.write_text(">ctg1\n" + "ACGT" * 200 + "\n")
    return genome


def _fake_skani(monkeypatch, edges: list[tuple[str, str, float, float, float]]):
    """Run skani's contract without skani: EVE-keyed edges, file-keyed output."""
    monkeypatch.setattr(shutil, "which", lambda name: "/fake/skani")

    def run(command, capture_output=True, text=True):
        listed = Path(command[command.index("-l") + 1]).read_text().split()
        file_by_eve = {
            Path(path).read_text().splitlines()[0][1:]: path for path in listed
        }
        rows = "".join(
            f"{file_by_eve[eve_a]}\t{file_by_eve[eve_b]}\t"
            f"{ani}\t{af_a}\t{af_b}\n"
            for eve_a, eve_b, ani, af_a, af_b in edges
        )
        Path(command[command.index("-o") + 1]).write_text(_SKANI_HEADER + rows)
        return SimpleNamespace(returncode=0, stderr="")

    monkeypatch.setattr(subprocess, "run", run)


def _edge_rows(output_dir: Path) -> list[str]:
    lines = (
        (output_dir / "phase3_synthesis" / "eve_ani_edges.tsv")
        .read_text()
        .splitlines()
    )
    assert lines[0].split("\t") == ["eve_a", "eve_b", "ani", "af_a", "af_b"]
    return lines[1:]


def test_mcp_class_propagates_to_mcp_free_cluster_members(
    monkeypatch,
    tmp_path: Path,
) -> None:
    results = [
        _result("EVE_1", taxonomy_class="NCLDV", has_mcp=True),
        _result("EVE_2", taxonomy_class="VIRAL_UNKNOWN", has_mcp=False),
        _result("EVE_3", taxonomy_class="VIRAL_UNKNOWN", has_mcp=False),
    ]
    _fake_skani(
        monkeypatch,
        [
            ("EVE_1", "EVE_2", 99.1, 88.0, 90.0),
            ("EVE_1", "EVE_3", 97.4, 71.0, 73.0),
        ],
    )

    cluster_accepted_eves(
        results,
        genome_fasta=_genome(tmp_path),
        output_dir=tmp_path,
        threads=2,
    )

    assert [r.taxonomy_class for r in results] == ["NCLDV", "NCLDV", "NCLDV"]
    assert [r.taxonomy_class_before_ani for r in results] == [
        "",
        "VIRAL_UNKNOWN",
        "VIRAL_UNKNOWN",
    ]
    assert [r.taxonomy_class_propagated_from for r in results] == [
        "",
        "EVE_1",
        "EVE_1",
    ]
    assert [r.cluster_id for r in results] == [0, 0, 0]
    assert [r.cluster_size for r in results] == [3, 3, 3]
    assert [r.max_cluster_ani for r in results] == [99.1, 99.1, 97.4]
    assert len(_edge_rows(tmp_path)) == 2
    # Confidence scoring ran before acceptance, so clustering must not feed it.
    assert all(r.clustering_bonus == 0.0 for r in results)


def test_structural_mcp_alone_does_not_donate(monkeypatch, tmp_path: Path) -> None:
    """has_mcp is broader than "an MCP marker decided this class".

    A capsid annotation, a structural jelly-roll call, or phylogenetic evidence
    all raise has_mcp without casting a taxonomy vote. Such an EVE must not hand
    its consensus-derived class to a relative.
    """
    donor = _result("EVE_1", taxonomy_class="NCLDV", has_mcp=True)
    donor.taxonomy_class_from_mcp = False
    results = [donor, _result("EVE_2", taxonomy_class="VIRAL_UNKNOWN", has_mcp=False)]
    _fake_skani(monkeypatch, [("EVE_1", "EVE_2", 99.1, 88.0, 90.0)])

    cluster_accepted_eves(
        results,
        genome_fasta=_genome(tmp_path),
        output_dir=tmp_path,
        threads=1,
    )

    assert [r.taxonomy_class for r in results] == ["NCLDV", "VIRAL_UNKNOWN"]
    assert [r.taxonomy_class_propagated_from for r in results] == ["", ""]
    # They still cluster; only the class handover is withheld.
    assert [r.cluster_size for r in results] == [2, 2]


def test_disagreeing_mcp_members_propagate_nothing(
    monkeypatch,
    tmp_path: Path,
) -> None:
    results = [
        _result("EVE_1", taxonomy_class="NCLDV", has_mcp=True),
        _result("EVE_2", taxonomy_class="MIRUS", has_mcp=True),
        _result("EVE_3", taxonomy_class="VIRAL_UNKNOWN", has_mcp=False),
    ]
    _fake_skani(
        monkeypatch,
        [
            ("EVE_1", "EVE_2", 99.1, 88.0, 90.0),
            ("EVE_1", "EVE_3", 98.2, 84.0, 85.0),
        ],
    )

    cluster_accepted_eves(
        results,
        genome_fasta=_genome(tmp_path),
        output_dir=tmp_path,
        threads=1,
    )

    assert [r.taxonomy_class for r in results] == [
        "NCLDV",
        "MIRUS",
        "VIRAL_UNKNOWN",
    ]
    assert all(r.taxonomy_class_propagated_from == "" for r in results)
    # The cluster is still recorded; only the relabelling is withheld.
    assert [r.cluster_size for r in results] == [3, 3, 3]


def test_viral_unknown_mcp_member_propagates_nothing(
    monkeypatch,
    tmp_path: Path,
) -> None:
    results = [
        _result("EVE_1", taxonomy_class="VIRAL_UNKNOWN", has_mcp=True),
        _result("EVE_2", taxonomy_class="UNKNOWN", has_mcp=False),
    ]
    _fake_skani(monkeypatch, [("EVE_1", "EVE_2", 99.9, 95.0, 96.0)])

    cluster_accepted_eves(
        results,
        genome_fasta=_genome(tmp_path),
        output_dir=tmp_path,
        threads=1,
    )

    assert [r.taxonomy_class for r in results] == ["VIRAL_UNKNOWN", "UNKNOWN"]
    assert all(r.taxonomy_class_before_ani == "" for r in results)


def test_below_threshold_pairs_stay_separate_components(
    monkeypatch,
    tmp_path: Path,
) -> None:
    ani, af_a, af_b = MIN_CLUSTER_ANI - 0.1, 90.0, 90.0
    results = [
        _result("EVE_1", taxonomy_class="NCLDV", has_mcp=True),
        _result("EVE_2", taxonomy_class="VIRAL_UNKNOWN", has_mcp=False),
    ]
    _fake_skani(monkeypatch, [("EVE_1", "EVE_2", ani, af_a, af_b)])

    cluster_accepted_eves(
        results,
        genome_fasta=_genome(tmp_path),
        output_dir=tmp_path,
        threads=1,
    )

    assert [r.taxonomy_class for r in results] == ["NCLDV", "VIRAL_UNKNOWN"]
    assert [r.cluster_id for r in results] == [-1, -1]
    assert [r.cluster_size for r in results] == [1, 1]
    # The sub-threshold pair is still published, so the notebook can filter it.
    assert len(_edge_rows(tmp_path)) == 1


def test_uniform_mcp_state_clusters_leave_classes_untouched(
    monkeypatch,
    tmp_path: Path,
) -> None:
    results = [
        _result("EVE_1", taxonomy_class="NCLDV", has_mcp=True),
        _result("EVE_2", taxonomy_class="MIRUS", has_mcp=True),
        _result("EVE_3", taxonomy_class="VIRAL_UNKNOWN", has_mcp=False),
        _result("EVE_4", taxonomy_class="UNKNOWN", has_mcp=False),
    ]
    _fake_skani(
        monkeypatch,
        [
            ("EVE_1", "EVE_2", 99.1, 90.0, 91.0),
            ("EVE_3", "EVE_4", 98.5, 80.0, 81.0),
        ],
    )

    cluster_accepted_eves(
        results,
        genome_fasta=_genome(tmp_path),
        output_dir=tmp_path,
        threads=1,
    )

    assert [r.taxonomy_class for r in results] == [
        "NCLDV",
        "MIRUS",
        "VIRAL_UNKNOWN",
        "UNKNOWN",
    ]
    assert all(r.taxonomy_class_propagated_from == "" for r in results)
    assert [r.cluster_id for r in results] == [0, 0, 1, 1]


def _no_skani_output(monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda name: "/fake/skani")
    monkeypatch.setattr(
        subprocess, "run",
        lambda *a, **k: SimpleNamespace(
            returncode=1, stderr="ERROR No genomes/sketches found.\n"
        ),
    )


def _no_skani_binary(monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda name: None)


@pytest.mark.parametrize(
    "arrange,eve_count",
    [
        (_no_skani_output, 2),   # skani ran and could not sketch the sequences
        (_no_skani_binary, 2),   # skani is not installed
        (lambda monkeypatch: None, 1),  # only one EVE, nothing to compare
    ],
)
def test_every_no_op_path_writes_a_header_only_table_and_changes_nothing(
    monkeypatch, tmp_path: Path, arrange, eve_count: int
) -> None:
    """Three different reasons to skip clustering, one required outcome.

    The header-only table matters: the notebook reads this file, and its absence
    is not the same statement as "no pairs".
    """
    results = [
        _result("EVE_1", taxonomy_class="NCLDV", has_mcp=True),
        _result("EVE_2", taxonomy_class="VIRAL_UNKNOWN", has_mcp=False),
    ][:eve_count]
    arrange(monkeypatch)

    cluster_accepted_eves(
        results, genome_fasta=_genome(tmp_path), output_dir=tmp_path, threads=1
    )

    assert _edge_rows(tmp_path) == []
    assert [r.cluster_id for r in results] == [-1] * eve_count
    assert [r.cluster_size for r in results] == [1] * eve_count
    assert [r.taxonomy_class for r in results] == [
        "NCLDV", "VIRAL_UNKNOWN"
    ][:eve_count]


def _accepted(eve_id: str, taxonomy_class: str, *, hallmarks: int, cluster_id: int):
    result = VerificationResult(
        eve_id=eve_id, scaffold="ctg1", start=0, end=600,
        taxonomy_class=taxonomy_class, hallmark_count=hallmarks,
    )
    result.cluster_id = cluster_id
    return result


def test_unsupported_eves_are_dropped_only_without_a_marker_bearing_relative() -> None:
    """UNKNOWN means no validated marker and no viral gene hit: nothing viral."""
    # Alone and unsupported: dropped.
    lone = _accepted("EVE_lone", "UNKNOWN", hallmarks=0, cluster_id=-1)
    assert unsupported_eve_ids([lone]) == {"EVE_lone"}

    # Clustered with a marker-bearing EVE: that relative vouches for it.
    rescued = _accepted("EVE_rescued", "UNKNOWN", hallmarks=0, cluster_id=0)
    donor = _accepted("EVE_donor", "NCLDV", hallmarks=3, cluster_id=0)
    assert unsupported_eve_ids([rescued, donor]) == set()

    # Clustered only with other unsupported EVEs: still nothing viral anywhere.
    pair = [
        _accepted("EVE_a", "UNKNOWN", hallmarks=0, cluster_id=1),
        _accepted("EVE_b", "UNKNOWN", hallmarks=0, cluster_id=1),
    ]
    assert unsupported_eve_ids(pair) == {"EVE_a", "EVE_b"}


def test_viral_unknown_and_classified_eves_are_never_dropped() -> None:
    # VIRAL_UNKNOWN carries a validated marker; it is viral, just unresolved.
    kept = [
        _accepted("EVE_vu", "VIRAL_UNKNOWN", hallmarks=2, cluster_id=-1),
        _accepted("EVE_ncldv", "NCLDV", hallmarks=2, cluster_id=-1),
        _accepted("EVE_ppv", "PPV", hallmarks=0, cluster_id=-1),
    ]
    assert unsupported_eve_ids(kept) == set()


def test_a_weight_settled_class_does_not_make_an_eve_a_donor() -> None:
    """Donating needs an MCP to have decided the class, not merely to be present.

    Two MCP markers that disagree hand the decision to the weighted vote. A class
    the weights chose is not one a capsid decided, so it must not propagate.
    """
    from virosync.pipeline.phase3.evidence_synthesizer import EvidenceSynthesizer

    marker = lambda prefixes, gene, porf: {  # noqa: E731
        "hallmark_gene": gene, "porf_id": porf, "top10_prefixes": prefixes,
        "top10_pidents": "80", "top10_targets": "", "validation_status": "validated",
        "score": 100.0,
    }
    result = VerificationResult(eve_id="EVE_x", scaffold="ctg1", start=0, end=600)
    synthesizer = EvidenceSynthesizer.__new__(EvidenceSynthesizer)
    hits = [
        marker("NCLDV", "plv_mcp_1", "p1"),
        marker("MIRUS", "mirus_mcp_2", "p2"),
        marker("NCLDV", "gvogm0100", "p3"),
    ]
    synthesizer._assign_taxonomy_class(result, hits, None)

    # NCLDV 5 + 2 against MIRUS 5, so the weights pick NCLDV.
    assert result.taxonomy_class == "NCLDV"
    # But no single MCP decided it, so this EVE cannot donate.
    assert result.taxonomy_class_from_mcp is False


def test_cluster_sizes_are_recounted_after_a_drop(monkeypatch, tmp_path: Path) -> None:
    """A published row must not claim a relative that is no longer published."""
    survivor = _result("EVE_1", taxonomy_class="PPV", has_mcp=False)
    dropped = _result("EVE_2", taxonomy_class="UNKNOWN", has_mcp=False)
    results = [survivor, dropped]
    _fake_skani(monkeypatch, [("EVE_1", "EVE_2", 99.1, 88.0, 90.0)])

    _edges, pairs = cluster_accepted_eves(
        results, genome_fasta=_genome(tmp_path), output_dir=tmp_path, threads=1
    )
    assert [r.cluster_size for r in results] == [2, 2]

    # Neither carries a marker, so nothing rescues the UNKNOWN member.
    assert unsupported_eve_ids(results) == {"EVE_2"}
    remaining = [r for r in results if r.eve_id != "EVE_2"]
    recluster_survivors(remaining, pairs)

    assert survivor.cluster_id == -1
    assert survivor.cluster_size == 1
    assert survivor.max_cluster_ani == 0.0
