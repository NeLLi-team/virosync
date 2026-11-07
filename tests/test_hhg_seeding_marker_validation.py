"""Regression test for NameError on `genome_fasta` in hhg_seeding.

Past bug: `validate_hmm_hits_with_combined_db` referenced `genome_fasta`
inside its body, but that name was not a parameter of the function. Any
happy-path run that reached the `filter_validated_markers` call raised
`NameError: name 'genome_fasta' is not defined`. The only caller
(`hhg_seeding_pipeline`) already carried `genome_fasta` in its own scope,
so the fix is to thread it through as an optional parameter.

This test exercises the full happy path (marker_db exists, Diamond
succeeds) with stubbed Diamond and validator calls, and asserts the
function returns cleanly instead of NameError-ing.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from virosync.pipeline.phase1.hhg_seeding import (
    HMMHit,
    validate_hmm_hits_with_combined_db,
)
from virosync.pipeline.phase1.marker_validation import (
    NovelMarkerCriteria,
    filter_validated_markers,
    validate_hmm_hit,
)


def test_marker_validation_uses_viral_reference_at_rank_twenty():
    hits = [(f"EUK__G{n}|p1", 100.0, 60.0, 1e-20) for n in range(19)]
    hits.append(("NCLDV__virus|p1", 50.0, 30.0, 1e-10))
    assert validate_hmm_hit(50, 1, hits, NovelMarkerCriteria()).value == "validated"


def test_production_marker_validation_and_class_vote_use_rank_twenty(tmp_path):
    from virosync.pipeline.phase3.evidence_synthesizer import marker_taxonomy_category

    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">scaffold_1 # 101 # 400 # 1 # ID=1_1\n" + "M" * 100 + "\n")
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text(
        "".join(f"scaffold_1\tEUK__G{n}|p1\t1e-30\t{100 - n}\t60\t100\n" for n in range(19))
        + "scaffold_1\tNCLDV__virus|p1\t1e-10\t50\t30\t100\n"
    )
    (marker,) = filter_validated_markers([HMMHit("scaffold_1", "VS000185", 50, 1e-8, 50, 1, 50)], diamond, proteome)
    assert marker.validation_status == "validated"
    assert len(marker.top10_targets.split(",")) == 20
    assert marker_taxonomy_category(marker) == "NCLDV"


def test_capsid_winner_does_not_restore_losing_non_capsid_on_same_segment(tmp_path):
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">scaffold_1 # 101 # 400 # 1 # ID=1_1\n" + "M" * 100 + "\n")
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("scaffold_1\tNCLDV__virus|p1\t1e-30\t100\t90\t100\n")
    markers = filter_validated_markers(
        [
            HMMHit("scaffold_1", "EVE_CAPSID_NOVEL_0002", 200, 1e-30, 200, 10, 50),
            HMMHit("scaffold_1", "VS000185", 50, 1e-8, 50, 10, 50),
            HMMHit("scaffold_1", "EVE_CAPSID_NOVEL_0012", 20, 0.1, 20, 10, 50),
        ],
        diamond,
        proteome,
    )
    assert [marker.hmm_target for marker in markers] == ["EVE_CAPSID_NOVEL_0002", "EVE_CAPSID_NOVEL_0012"]
    assert all(marker.validation_status == "validated" for marker in markers)


def test_production_marker_validation_does_not_use_rank_twenty_one(tmp_path):
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">scaffold_1 # 101 # 400 # 1 # ID=1_1\n" + "M" * 100 + "\n")
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text(
        "".join(f"scaffold_1\tEUK__G{n}|p1\t1e-30\t{100 - n}\t60\t100\n" for n in range(20))
        + "scaffold_1\tNCLDV__virus|p1\t1e-10\t50\t30\t100\n"
    )
    (marker,) = filter_validated_markers(
        [HMMHit("scaffold_1", "VS000185", 50, 1e-8, 50, 1, 50)], diamond, proteome, max_seqs=20
    )
    assert marker.validation_status == "unvalidated"
    assert len(marker.top10_targets.split(",")) == 20


def test_competing_segment_retains_unvalidated_capsid_and_full_gene(tmp_path):
    from virosync.ablation import AblationID
    from virosync.orchestration._flows.single_genome.phase1 import _candidate_regions_to_seeds
    from virosync.pipeline.phase1.marker_roles import decide_marker_hit_role
    from virosync.pipeline.phase1.region_assembly import assemble_candidate_regions

    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">scaffold_1 # 101 # 400 # 1 # ID=1_1\n" + "M" * 100 + "\n")
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("scaffold_1\tEUK__host|p1\t1e-30\t100\t90\t100\n")
    hits = [
        HMMHit("scaffold_1", "VS000087", 200, 1e-30, 200, 10, 50),
        HMMHit("scaffold_1", "EVE_CAPSID_NOVEL_0002", -5, 1000, -5, 10, 50),
    ]
    markers = filter_validated_markers(hits, diamond, proteome)
    assert {marker.hmm_target for marker in markers} == {"VS000087", "EVE_CAPSID_NOVEL_0002"}
    capsid = next(marker for marker in markers if marker.hmm_target == "EVE_CAPSID_NOVEL_0002")
    assert (capsid.start, capsid.end, capsid.validation_status) == (100, 400, "unvalidated")
    retained = [
        marker for marker in markers if decide_marker_hit_role(marker, ablation_id=AblationID.A0).is_retained_evidence
    ]
    assert retained == [capsid]
    genome = tmp_path / "genome.fna"
    genome.write_text(">scaffold\n" + "A" * 1000 + "\n")
    regions = assemble_candidate_regions(retained, genome, proteome, tmp_path / "regions")
    seeds = _candidate_regions_to_seeds(regions, markers, [], [])
    assert len(seeds) == 1
    assert (seeds[0].confirmed_seeds[0].start, seeds[0].confirmed_seeds[0].end) == (100, 400)


def _fake_extract(hits, proteome_fasta, output_fasta):
    """Stub for extract_hmm_hit_sequences: create an empty FASTA and claim 1 seq."""
    Path(output_fasta).write_text(">q1\nM\n")
    return 1


def _fake_run_diamond(hmm_hit_fasta, diamond_db, output_tsv, threads, evalue, max_seqs, masked_genomes):
    """Stub for run_diamond_on_hmm_hits: create an empty output TSV."""
    Path(output_tsv).write_text("")


def test_validate_hmm_hits_does_not_nameerror_on_happy_path(tmp_path: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">contig_1_1\nMGKPSALVR\n")

    marker_db = tmp_path / "marker.dmnd"
    marker_db.write_bytes(b"")  # presence is enough; Diamond is mocked

    hits = [
        HMMHit(
            query_name="contig_1_1",
            target_name="mcp",
            score=120.0,
            evalue=1e-40,
            domain_score=120.0,
            query_start=1,
            query_end=9,
        )
    ]

    with (
        patch(
            "virosync.pipeline.phase1.hhg_seeding.extract_hmm_hit_sequences",
            side_effect=_fake_extract,
        ),
        patch(
            "virosync.pipeline.phase1.hhg_seeding.run_diamond_on_hmm_hits",
            side_effect=_fake_run_diamond,
        ),
        patch(
            "virosync.pipeline.phase1.hhg_seeding.filter_validated_markers",
            return_value=[],
        ),
    ):
        validated_hits, validated_markers = validate_hmm_hits_with_combined_db(
            hits=hits,
            proteome_fasta=proteome,
            marker_db=marker_db,
            threads=1,
            top_k=10,
            output_dir=tmp_path / "out",
        )

    # No validated markers were returned (our mock returned []),
    # so validated_hits should be empty, not a NameError traceback.
    assert validated_hits == []
    assert validated_markers == []


def test_validate_hmm_hits_forwards_optional_genome_fasta(tmp_path: Path) -> None:
    """After the fix, `genome_fasta` is an optional kwarg and is forwarded to
    `filter_validated_markers`. Verify the forwarding without assuming the
    caller ever sets it."""
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">contig_1_1\nMGKPSALVR\n")
    genome = tmp_path / "genome.fna"
    genome.write_text(">contig_1\nACGT\n")
    marker_db = tmp_path / "marker.dmnd"
    marker_db.write_bytes(b"")

    hits = [
        HMMHit(
            query_name="contig_1_1",
            target_name="mcp",
            score=120.0,
            evalue=1e-40,
            domain_score=120.0,
            query_start=1,
            query_end=9,
        )
    ]

    seen_kwargs: dict = {}

    def _capture_filter(**kwargs):
        seen_kwargs.update(kwargs)
        return []

    with (
        patch(
            "virosync.pipeline.phase1.hhg_seeding.extract_hmm_hit_sequences",
            side_effect=_fake_extract,
        ),
        patch(
            "virosync.pipeline.phase1.hhg_seeding.run_diamond_on_hmm_hits",
            side_effect=_fake_run_diamond,
        ),
        patch(
            "virosync.pipeline.phase1.hhg_seeding.filter_validated_markers",
            side_effect=_capture_filter,
        ),
    ):
        validate_hmm_hits_with_combined_db(
            hits=hits,
            proteome_fasta=proteome,
            marker_db=marker_db,
            threads=1,
            top_k=10,
            output_dir=tmp_path / "out",
            genome_fasta=genome,
        )

    assert seen_kwargs.get("genome_fasta") == genome


def test_validate_hmm_hits_gives_the_reference_mask_to_the_marker_search(tmp_path: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">contig_1_1\nMGKPSALVR\n")
    marker_db = tmp_path / "marker.dmnd"
    marker_db.write_bytes(b"")
    hits = [HMMHit("contig_1_1", "mcp", 120.0, 1e-40, 120.0, 1, 9)]
    masks_seen: list[frozenset[str]] = []

    def _capture_search(**kwargs):
        masks_seen.append(kwargs["masked_genomes"])
        Path(kwargs["output_tsv"]).write_text("")

    with (
        patch(
            "virosync.pipeline.phase1.hhg_seeding.extract_hmm_hit_sequences",
            side_effect=_fake_extract,
        ),
        patch(
            "virosync.pipeline.phase1.hhg_seeding.run_diamond_on_hmm_hits",
            side_effect=_capture_search,
        ),
        patch(
            "virosync.pipeline.phase1.hhg_seeding.filter_validated_markers",
            return_value=[],
        ),
    ):
        validate_hmm_hits_with_combined_db(
            hits=hits,
            proteome_fasta=proteome,
            marker_db=marker_db,
            masked_genomes=frozenset({"GCA-000001-1"}),
        )

    assert masks_seen == [frozenset({"GCA-000001-1"})]


def test_filter_uses_prodigal_coordinates_without_reading_genome(tmp_path: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">contig_1_1 # 1 # 300 # + # ID=1_1;partial=00\n"
        f"{'M' * 100}\n"
        ">contig_1_2 # 301 # 600 # + # ID=1_2;partial=00\n"
        f"{'M' * 100}\n"
        ">contig_2_1 # 1 # 300 # + # ID=2_1;partial=00\n"
        f"{'M' * 100}\n"
    )
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("")
    missing_genome = tmp_path / "not-read.fna"
    hits = [HMMHit(query, "mcp", 120.0, 1e-40, 120.0, 1, 100) for query in ("contig_1_1", "contig_1_2", "contig_2_1")]

    markers = filter_validated_markers(
        hits,
        diamond,
        proteome,
        genome_fasta=missing_genome,
        novel_criteria=NovelMarkerCriteria(
            min_hmm_coverage=0.3,
            require_cluster=False,
        ),
    )

    assert [marker.validation_status for marker in markers] == [
        "validated_novel",
        "validated_novel",
        "validated_novel",
    ]
    assert [(marker.scaffold, marker.start, marker.end) for marker in markers] == [
        ("contig_1", 0, 300),
        ("contig_1", 300, 600),
        ("contig_2", 0, 300),
    ]


def test_novel_coverage_uses_full_protein_length(tmp_path: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">contig_1_1 # 1 # 300 # + # ID=1_1;partial=00\n"
        f"{'M' * 100}\n"
        ">contig_1_2 # 301 # 3300 # + # ID=1_2;partial=00\n"
        f"{'M' * 1000}\n"
    )
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("")
    hits = [
        HMMHit("contig_1_1", "mcp", 80.0, 1e-20, 80.0, 1, 100),
        HMMHit("contig_1_2", "mcp", 80.0, 1e-20, 80.0, 1, 150),
    ]

    markers = filter_validated_markers(
        hits,
        diamond,
        proteome,
        novel_criteria=NovelMarkerCriteria(require_cluster=False),
    )

    assert [marker.validation_status for marker in markers] == [
        "validated_novel",
        "unvalidated",
    ]


def test_novel_marker_with_unknown_protein_length_is_unvalidated(tmp_path: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">other\nMMMM\n")
    genome = tmp_path / "genome.fna"
    genome.write_text(f">contig_1\n{'A' * 300}\n")
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("")
    hit = HMMHit("contig_1_frame=1", "mcp", 80.0, 1e-20, 80.0, 1, 100)

    markers = filter_validated_markers(
        [hit],
        diamond,
        proteome,
        genome_fasta=genome,
        novel_criteria=NovelMarkerCriteria(
            min_hmm_coverage=0.0,
            require_cluster=False,
        ),
    )

    assert markers[0].validation_status == "unvalidated"


def test_nearby_diamond_validated_protein_supports_novel_marker_in_any_order(
    tmp_path: Path,
) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">contig_1_1 # 1 # 300 # + # ID=1_1;partial=00\n"
        f"{'M' * 100}\n"
        ">contig_1_2 # 10001 # 10300 # + # ID=1_2;partial=00\n"
        f"{'M' * 100}\n"
    )
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text(
        "".join(f"contig_1_2\tEUK__host_{index}\t1e-40\t{300 - index}\t80\t100\n" for index in range(10))
        + "contig_1_2\tNCLDV__reference\t1e-40\t200\t35\t100\n"
    )
    novel = HMMHit("contig_1_1", "mcp", 80.0, 1e-20, 80.0, 1, 100)
    validated = HMMHit("contig_1_2", "polb", 80.0, 1e-20, 80.0, 1, 100)

    forward = filter_validated_markers([novel, validated], diamond, proteome, max_seqs=11)
    reverse = filter_validated_markers([validated, novel], diamond, proteome, max_seqs=11)

    expected = {
        "contig_1_1|aa1-100": "validated_novel",
        "contig_1_2|aa1-100": "validated",
    }
    assert {marker.query_porf: marker.validation_status for marker in forward} == expected
    assert {marker.query_porf: marker.validation_status for marker in reverse} == expected


@pytest.mark.parametrize(
    ("support_query", "support_start", "support_pident"),
    [
        # A 20,001 bp gap exceeds the default 20 kb cluster-support window.
        ("contig_1_2", 20302, 35.0),
        ("contig_1_2", 10001, 24.9),
        ("contig_2_1", 10001, 35.0),
    ],
    ids=["distant", "weak", "cross-scaffold"],
)
def test_invalid_cluster_support_does_not_validate_novel_marker(
    tmp_path: Path,
    support_query: str,
    support_start: int,
    support_pident: float,
) -> None:
    support_index = support_query.rsplit("_", 1)[1]
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">contig_1_1 # 1 # 300 # + # ID=1_1;partial=00\n"
        f"{'M' * 100}\n"
        f">{support_query} # {support_start} # {support_start + 299} # + # "
        f"ID=2_{support_index};partial=00\n"
        f"{'M' * 100}\n"
    )
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text(f"{support_query}\tNCLDV__reference\t1e-40\t200\t{support_pident}\t100\n")
    hits = [
        HMMHit("contig_1_1", "mcp", 80.0, 1e-20, 80.0, 1, 100),
        HMMHit(support_query, "polb", 80.0, 1e-20, 80.0, 1, 100),
    ]

    markers = filter_validated_markers(hits, diamond, proteome)

    assert markers[0].validation_status == "unvalidated"


def test_nearby_hmm_only_hits_do_not_support_each_other(tmp_path: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">contig_1_1 # 1 # 300 # + # ID=1_1;partial=00\n"
        f"{'M' * 100}\n"
        ">contig_1_2 # 301 # 600 # + # ID=1_2;partial=00\n"
        f"{'M' * 100}\n"
    )
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("")
    hits = [HMMHit(query, "mcp", 80.0, 1e-20, 80.0, 1, 100) for query in ("contig_1_1", "contig_1_2")]

    markers = filter_validated_markers(hits, diamond, proteome)

    assert [marker.validation_status for marker in markers] == [
        "unvalidated",
        "unvalidated",
    ]


def test_second_hmm_segment_on_same_protein_cannot_support_novel_segment(
    tmp_path: Path,
) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(f">contig_1_1 # 1 # 300 # + # ID=1_1;partial=00\n{'M' * 100}\n")
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("contig_1_1|aa61-100\tNCLDV__reference\t1e-40\t200\t35\t100\n")
    hits = [
        HMMHit("contig_1_1", "mcp", 80.0, 1e-20, 80.0, 1, 60),
        HMMHit("contig_1_1", "polb", 80.0, 1e-20, 80.0, 61, 100),
    ]

    markers = filter_validated_markers(hits, diamond, proteome)

    assert [marker.validation_status for marker in markers] == [
        "unvalidated",
        "validated",
    ]


def test_filter_reads_genome_for_frame_coordinates(tmp_path: Path) -> None:
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">contig_1_frame=1\nMMMMMMMMMM\n")
    genome = tmp_path / "genome.fna"
    genome.write_text(f">contig_1\n{'A' * 300}\n")
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("")
    hit = HMMHit("contig_1_frame=1", "mcp", 120.0, 1e-40, 120.0, 1, 10)

    markers = filter_validated_markers(
        [hit],
        diamond,
        proteome,
        genome_fasta=genome,
        novel_criteria=NovelMarkerCriteria(require_cluster=False),
    )

    assert len(markers) == 1
    assert markers[0].validation_status == "validated_novel"
    assert (markers[0].scaffold, markers[0].start, markers[0].end, markers[0].strand) == (
        "contig_1",
        0,
        30,
        "+",
    )


@pytest.mark.parametrize("override, expected", [("contradicted", "unvalidated"), ("", "validated_novel")])
def test_only_ordinary_capsid_can_supply_novel_marker_neighbor_support(tmp_path, override, expected):
    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        ">scaffold_1 # 101 # 400 # 1 # ID=1_1\n" + "M" * 100 + "\n"
        ">scaffold_2 # 501 # 800 # 1 # ID=1_2\n" + "M" * 100 + "\n"
    )
    diamond = tmp_path / "diamond.tsv"
    diamond.write_text("scaffold_1\tNCLDV__virus|p1\t1e-30\t100\t90\t100\n")
    markers = filter_validated_markers(
        [
            HMMHit("scaffold_1", "GVOGm0003", 100, 1e-20, 100, 1, 100, pfam_override=override),
            HMMHit("scaffold_2", "GVOGm0054", 100, 1e-20, 100, 1, 100),
        ],
        diamond,
        proteome,
    )
    assert next(marker for marker in markers if marker.hmm_target == "GVOGm0054").validation_status == expected
    assert next(marker for marker in markers if marker.hmm_target == "GVOGm0003").validation_status == "validated"


def test_marker_wrapper_preserves_qualifying_viral_member_of_selected_genome(tmp_path, monkeypatch):
    from virosync.pipeline.phase1.marker_validation import run_diamond_on_hmm_hits
    from virosync.pipeline.phase1.region_assembly import can_seed_marker

    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">ctg_1 # 1 # 300 # 1 # ID=1_1\n" + "M" * 100 + "\n")

    def diamond(**kwargs):
        kwargs["output_tsv"].write_text(
            "ctg_1\tNCLDV__v|weak\t1e-20\t300\t24\t100\nctg_1\tNCLDV__v|qualifying\t1e-20\t290\t40\t100\n"
        )

    monkeypatch.setattr("virosync.pipeline.search_backend._run_diamond", diamond)
    output = tmp_path / "selected.tsv"
    run_diamond_on_hmm_hits(proteome, tmp_path / "db", output, threads=1)
    markers = filter_validated_markers([HMMHit("ctg_1", "GVOGm0054", 100, 1e-20, 100, 1, 100)], output, proteome)
    assert len(markers) == 1
    assert markers[0].validation_status == "validated"
    assert "qualifying" in markers[0].top10_targets
    assert can_seed_marker(markers[0].hmm_target, markers[0].validation_status)


def test_marker_preference_keeps_mask_twins_and_first20_genome_slots(tmp_path, monkeypatch):
    from virosync.pipeline import search_backend
    from virosync.pipeline.phase1.marker_validation import run_diamond_on_hmm_hits

    rows = [
        "q\tNCLDV__masked|p\t1e-20\t400\t90\t100\n",
        "q\tEUK__v|twin\t1e-20\t350\t80\t100\n",
        "q\tNCLDV__v|twin\t1e-20\t340\t24\t100\n",
    ]
    rows.extend(f"q\tEUK__g{i}|p\t1e-20\t{330 - i}\t80\t100\n" for i in range(19))
    rows.extend(["q\tNCLDV__v|qualifying\t1e-20\t290\t40\t100\n", "q\tNCLDV__late|p\t1e-20\t280\t90\t100\n"])

    def diamond(**kwargs):
        kwargs["output_tsv"].write_text("".join(rows))

    monkeypatch.setattr(search_backend, "_run_diamond", diamond)
    selected = tmp_path / "marker.tsv"
    run_diamond_on_hmm_hits(
        tmp_path / "q.faa", tmp_path / "db", selected, threads=1, masked_genomes=frozenset({"masked"})
    )
    targets = [line.split("\t")[1] for line in selected.read_text().splitlines()]
    assert targets == ["NCLDV__v|qualifying", *[f"EUK__g{i}|p" for i in range(19)]]
    ordinary = tmp_path / "ordinary.tsv"
    search_backend.run_genome_filtered_search(
        tmp_path / "q.faa",
        tmp_path / "db",
        ordinary,
        threads=1,
        max_target_seqs=20,
        masked_genomes=frozenset({"masked"}),
    )
    assert ordinary.read_text().splitlines()[0].split("\t")[1] == "NCLDV__v|twin"


def test_marker_representative_restores_nearby_three_core_rescue(tmp_path, monkeypatch):
    from virosync.orchestration.tasks import marker_validation_task
    from virosync.pipeline.hallmark_evidence import add_rescue_seeds, collect_hallmark_evidence, find_spatial_witnesses

    proteome = tmp_path / "proteome.faa"
    proteome.write_text(
        "".join(
            f">ctg_{i + 1} # {i * 1000 + 1} # {i * 1000 + 300} # 1 # ID=1_{i + 1}\n" + "M" * 100 + "\n"
            for i in range(4)
        )
    )
    genome = tmp_path / "genome.fna"
    genome.write_text(">ctg\n" + "A" * 5000 + "\n")
    hits = [
        HMMHit(f"ctg_{i + 1}", model, 100, 1e-20, 100, 1, 100)
        for i, model in enumerate(["VS000375", "GVOGm0054", "GVOGm0760", "VS000106"])
    ]

    def diamond(**kwargs):
        kwargs["output_tsv"].write_text(
            "ctg_1\tNCLDV__v|weak\t1e-20\t300\t24\t100\nctg_1\tNCLDV__v|qualifying\t1e-20\t290\t40\t100\n"
        )

    monkeypatch.setattr("virosync.pipeline.search_backend._run_diamond", diamond)
    markers, _ = marker_validation_task(
        hits, proteome, tmp_path / "db", tmp_path / "out", genome_path=genome, threads=1
    )
    assert [m.validation_status for m in markers] == [
        "validated",
        "validated_novel",
        "validated_novel",
        "validated_novel",
    ]
    witnesses = find_spatial_witnesses(collect_hallmark_evidence(markers, hits))
    assert [(w.rule, w.start, w.end) for w in witnesses] == [("three_core_15kb", 1000, 3300)]
    assert [(seed.start, seed.end) for seed in add_rescue_seeds([], witnesses)] == [(1000, 3300)]
