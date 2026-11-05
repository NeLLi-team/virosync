"""One retained search supplies distinct positive and host-veto slices."""

import pytest

from virosync.pipeline import search_backend
from virosync.pipeline.phase2 import boundary_diamond as bd


@pytest.mark.parametrize("path", ["batched", "chunked", "superset"])
def test_one_native_search_supplies_positive_and_veto_depths(tmp_path, monkeypatch, path):
    fasta = tmp_path / "proteins.faa"
    fasta.write_text(">g1\nMAAA\n>g2\nMAAA\n")
    genes = [bd.pORF("g1", "s", 0, 100), bd.pORF("g2", "s", 100, 200)]
    query = bd.GenomeDiamondQuery(all_porf_ids=[g.id for g in genes])
    calls = []

    def native(**kwargs):
        calls.append((kwargs["max_target_seqs"], kwargs["threads"]))
        names = [line[1:] for line in kwargs["query_fasta"].read_text().splitlines() if line.startswith(">")]
        subjects = [f"EUK__G{i}|p" for i in range(14)] + ["NCLDV__virus|p"]
        kwargs["output_tsv"].write_text(
            "".join(
                f"{name}\t{target}\t1e-20\t{1000 - i}\t80\t90\n" for name in names for i, target in enumerate(subjects)
            )
        )

    monkeypatch.setattr(search_backend, "_run_diamond", native)
    config = bd.BoundaryDiamondConfig(viral_hit_depth=50, threads=32, chunk_size=1 if path == "chunked" else 10)
    if path == "superset":
        hits = bd.run_full_proteome_diamond(fasta, tmp_path / "db", tmp_path / "out", max_target_seqs=50, threads=32)
        taxonomy = bd.classify_cached_diamond_query(query, hits, {"s": genes}, config)
    else:
        taxonomy = bd.run_batched_diamond(query, fasta, tmp_path / "db", tmp_path / "out", {"s": genes}, config)

    assert calls == [(500, 8)] * (2 if path == "chunked" else 1)
    assert taxonomy["g1"].top10_targets == [f"EUK__G{i}|p" for i in range(10)]
    assert taxonomy["g1"].has_viral is False
    assert taxonomy["g1"].viral_hit_within_depth is True
    assert taxonomy["g1"].host_competition_viral_bits == 986
