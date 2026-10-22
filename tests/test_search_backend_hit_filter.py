"""Tests for the top-hit filter in the shared DIAMOND search wrapper."""

from pathlib import Path

import pytest

from virosync.pipeline import search_backend
from virosync.pipeline.phase2.boundary_diamond import parse_diamond_output
from virosync.pipeline.phase3.gene_taxonomy import parse_diamond_top10
from virosync.pipeline.search_backend import (
    TARGET_FETCH_MULTIPLIER,
    filter_hits_file,
    run_sequence_search,
    select_target_hits,
)


def _kept(subject_ids: list[str], max_targets: int = 10) -> list[str]:
    return [subject_ids[idx] for idx in select_target_hits(subject_ids, max_targets)]


def test_select_target_hits_drops_repeated_subject_ids() -> None:
    subjects = ["EUK__G1|p1", "EUK__G1|p1", "EUK__G1|p1", "NCLDV__V1|p9"]

    assert _kept(subjects) == ["EUK__G1|p1", "NCLDV__V1|p9"]


def test_select_target_hits_keeps_best_hit_of_each_genome() -> None:
    subjects = ["EUK__G1|p1", "EUK__G1|p2", "EUK__G2|p1", "EUK__G1|p3"]

    assert _kept(subjects) == ["EUK__G1|p1", "EUK__G2|p1"]


def test_select_target_hits_treats_labels_as_part_of_the_genome() -> None:
    subjects = ["BAC__G1|p1", "PHAGE__G1|p2"]

    assert _kept(subjects) == subjects


def test_select_target_hits_twin_keeps_viral_row_when_cellular_comes_first() -> None:
    subjects = ["EUK__G1|p1", "EUK__G2|p5", "PPV__G1|p1"]

    assert _kept(subjects) == ["PPV__G1|p1", "EUK__G2|p5"]


def test_select_target_hits_twin_keeps_viral_row_when_viral_comes_first() -> None:
    subjects = ["PPV__G1|p1", "EUK__G2|p5", "EUK__G1|p1"]

    assert _kept(subjects) == ["PPV__G1|p1", "EUK__G2|p5"]


def test_select_target_hits_twin_moves_viral_row_to_first_occurrence() -> None:
    subjects = ["EUK__G1|p1", "EUK__G2|p5", "BAC__G3|p6", "NCLDV__G1|p1"]

    assert select_target_hits(subjects, 10) == [3, 1, 2]


@pytest.mark.parametrize("unrelated", ["OTHER__G|p", "G|p"], ids=["unknown-label", "no-label"])
def test_select_target_hits_twin_keeps_an_unrelated_row_of_the_same_protein(unrelated: str) -> None:
    subjects = [unrelated, "EUK__G|p", "NCLDV__G|p"]

    assert _kept(subjects) == [unrelated, "NCLDV__G|p"]


@pytest.mark.parametrize("cellular", ["EUK", "BAC", "ARC", "MITO", "PLASTID"])
def test_select_target_hits_twin_rule_covers_every_cellular_label(cellular: str) -> None:
    subjects = [f"{cellular}__G1|p1", "MIRUS__G1|p1"]

    assert _kept(subjects) == ["MIRUS__G1|p1"]


def test_select_target_hits_keeps_two_cellular_copies_of_one_protein() -> None:
    subjects = ["EUK__G1|p1", "MITO__G1|p1"]

    assert _kept(subjects) == subjects


def test_select_target_hits_truncates_after_filtering() -> None:
    subjects = [f"EUK__G{n // 5}|p{n}" for n in range(50)] + ["NCLDV__V1|p1"]

    kept = _kept(subjects, max_targets=3)

    assert kept == ["EUK__G0|p0", "EUK__G1|p5", "EUK__G2|p10"]


def test_select_target_hits_reaches_viral_hit_buried_under_host_copies() -> None:
    subjects = ["EUK__G1|p1"] * 7 + [f"EUK__G{n}|p1" for n in range(2, 10)]
    subjects.append("NCLDV__V1|p1")

    assert _kept(subjects, max_targets=10)[-1] == "NCLDV__V1|p1"


def test_select_target_hits_keeps_ids_without_a_label() -> None:
    subjects = ["plainid", "plainid", "other|p1", "other|p2"]

    assert _kept(subjects) == ["plainid", "other|p1"]


def _write_hits(path: Path, rows: list[tuple[str, str]]) -> None:
    lines = [f"{query}\t{subject}\t1e-30\t{200 - n}\t55.0\t90\n" for n, (query, subject) in enumerate(rows)]
    path.write_text("".join(lines), encoding="utf-8")


def test_filter_hits_file_filters_each_query_on_its_own(tmp_path: Path) -> None:
    hits = tmp_path / "hits.tsv"
    _write_hits(
        hits,
        [
            ("q1", "EUK__G1|p1"),
            ("q1", "EUK__G1|p2"),
            ("q1", "NCLDV__V1|p1"),
            ("q2", "EUK__G1|p1"),
            ("q2", "EUK__G2|p1"),
        ],
    )

    filter_hits_file(hits, ["qseqid", "sseqid", "evalue", "bitscore", "pident", "qcovhsp"], 2)

    rows = [line.split("\t")[:2] for line in hits.read_text(encoding="utf-8").splitlines()]
    assert rows == [
        ["q1", "EUK__G1|p1"],
        ["q1", "NCLDV__V1|p1"],
        ["q2", "EUK__G1|p1"],
        ["q2", "EUK__G2|p1"],
    ]


def test_filter_hits_file_finds_columns_by_name(tmp_path: Path) -> None:
    hits = tmp_path / "hits.tsv"
    hits.write_text(
        "q1\tEUK__G1|p1\t99\nq1\tEUK__G1|p2\t98\nq1\tBAC__G2|p1\t97\n",
        encoding="utf-8",
    )

    filter_hits_file(hits, ["qseqid", "sseqid", "bitscore"], 10)

    assert hits.read_text(encoding="utf-8") == "q1\tEUK__G1|p1\t99\nq1\tBAC__G2|p1\t97\n"


def test_filter_hits_file_rejects_columns_without_the_ids(tmp_path: Path) -> None:
    hits = tmp_path / "hits.tsv"
    hits.write_text("q1\t99\n", encoding="utf-8")

    with pytest.raises(ValueError, match="sseqid"):
        filter_hits_file(hits, ["qseqid", "bitscore"], 10)


class _FakeDiamond:
    """Stand-in for subprocess.Popen that records the command and writes hits."""

    commands: list[list[str]] = []
    rows: list[tuple[str, str]] = []

    def __init__(self, cmd: list[str], **_kwargs: object) -> None:
        self.pid = 0
        type(self).commands.append(cmd)
        output = Path(cmd[cmd.index("--out") + 1])
        _write_hits(output, type(self).rows)

    def poll(self) -> int:
        return 0


@pytest.fixture
def fake_diamond(monkeypatch: pytest.MonkeyPatch) -> type[_FakeDiamond]:
    _FakeDiamond.commands = []
    _FakeDiamond.rows = []
    monkeypatch.setattr(search_backend.subprocess, "Popen", _FakeDiamond)
    return _FakeDiamond


def _search(tmp_path: Path, **kwargs: object) -> Path:
    return run_sequence_search(
        query_fasta=tmp_path / "q.faa",
        db_path=tmp_path / "db.dmnd",
        output_tsv=tmp_path / "out.tsv",
        threads=1,
        **kwargs,
    )


def test_run_sequence_search_asks_diamond_for_a_multiple_of_the_targets(
    tmp_path: Path, fake_diamond: type[_FakeDiamond]
) -> None:
    _search(tmp_path, max_target_seqs=10)

    cmd = fake_diamond.commands[0]
    assert cmd[cmd.index("--max-target-seqs") + 1] == str(10 * TARGET_FETCH_MULTIPLIER)
    assert TARGET_FETCH_MULTIPLIER == 10


def test_run_sequence_search_truncates_the_filtered_hits_to_the_targets(
    tmp_path: Path, fake_diamond: type[_FakeDiamond]
) -> None:
    fake_diamond.rows = [("q1", f"EUK__G{n}|p{n}") for n in range(30)]

    output = _search(tmp_path, max_target_seqs=10)

    assert len(output.read_text(encoding="utf-8").splitlines()) == 10


def test_run_sequence_search_leaves_unlimited_search_unfiltered(
    tmp_path: Path, fake_diamond: type[_FakeDiamond]
) -> None:
    fake_diamond.rows = [("q1", "EUK__G1|p1")] * 3

    output = _search(tmp_path, max_target_seqs=0)

    cmd = fake_diamond.commands[0]
    assert cmd[cmd.index("--max-target-seqs") + 1] == "0"
    assert len(output.read_text(encoding="utf-8").splitlines()) == 3


def test_run_sequence_search_output_still_parses_in_consumers(tmp_path: Path, fake_diamond: type[_FakeDiamond]) -> None:
    fake_diamond.rows = [
        ("q1", "EUK__G1|p1"),
        ("q1", "EUK__G1|p1"),
        ("q1", "EUK__G1|p2"),
        ("q1", "NCLDV__V1|p3"),
    ]

    output = _search(tmp_path, max_target_seqs=10)

    boundary_hits = parse_diamond_output(output)
    gene_hits = parse_diamond_top10(output)
    assert [hit.target for hit in boundary_hits["q1"]] == ["EUK__G1|p1", "NCLDV__V1|p3"]
    assert [hit.target for hit in gene_hits["q1"]] == ["EUK__G1|p1", "NCLDV__V1|p3"]
