"""Keep the original native search independent from deeper veto retrieval."""

from dataclasses import replace

import pytest

from virosync.pipeline import search_backend
from virosync.pipeline.phase2 import boundary_diamond as bd


@pytest.mark.parametrize("path", ["batched", "chunked", "superset"])
@pytest.mark.parametrize("tail", ["NCLDV__new|p", "EUK-pEVE__G0|p0"])
def test_native_fetch_window_mask_duplicates_and_late_twin_do_not_change_top10(tmp_path, monkeypatch, path, tail):
    fasta = tmp_path / "proteins.faa"
    fasta.write_text(">g1\nMAAA\n>g2\nMAAA\n")
    genes = [bd.pORF("g1", "s", 0, 100), bd.pORF("g2", "s", 100, 200)]
    query = bd.GenomeDiamondQuery(all_porf_ids=[g.id for g in genes])
    subjects = [f"EUK__self|p{i}" for i in range(20)] + [f"EUK__G{i % 6}|p{i}" for i in range(80)] + [tail]
    calls = []

    def native(**kwargs):
        limit = kwargs["max_target_seqs"]
        calls.append(limit)
        names = [line[1:].split()[0] for line in kwargs["query_fasta"].read_text().splitlines() if line.startswith(">")]
        # Native -k can change even the scores inside the apparent common window.
        base = 2000 if limit > 100 else 1000
        kwargs["output_tsv"].write_text(
            "".join(
                f"{name}\t{target}\t1e-20\t{base - i}\t80\t90\n"
                for name in names
                for i, target in enumerate(subjects[:limit])
            )
        )

    monkeypatch.setattr(search_backend, "_run_diamond", native)

    def run(depth):
        config = bd.BoundaryDiamondConfig(
            viral_hit_depth=depth, masked_genomes=frozenset({"self"}), chunk_size=1 if path == "chunked" else 10
        )
        directory = tmp_path / f"depth-{depth}"
        if path == "superset":
            common = dict(
                proteome_fasta=fasta, diamond_db=tmp_path / "db", threads=1, masked_genomes=config.masked_genomes
            )
            hits = bd.run_full_proteome_diamond(**common, output_dir=directory / "original", max_target_seqs=10)
            veto = (
                bd.run_full_proteome_diamond(**common, output_dir=directory / "veto", max_target_seqs=depth)
                if depth > 10
                else None
            )
            return bd.classify_cached_diamond_query(query, hits, {"s": genes}, config, viral_veto_hits=veto)
        return bd.run_batched_diamond(query, fasta, tmp_path / "db", directory, {"s": genes}, config)

    original = run(10)
    assert calls == [100] * (2 if path == "chunked" else 1)
    calls.clear()
    deep = run(50)
    assert calls == [100] * (2 if path == "chunked" else 1) + [500] * (2 if path == "chunked" else 1)
    assert all(row.viral_hit_within_depth and not row.has_viral for row in deep.values())
    # Scores in the independent native deep search differ by design. Only the
    # new competition diagnostics and deep veto may change, never original10.
    assert {
        key: replace(
            row,
            viral_hit_within_depth=original[key].viral_hit_within_depth,
            host_competition_assessed=original[key].host_competition_assessed,
            host_competition_viral_bits=original[key].host_competition_viral_bits,
            host_competition_cellular_bits=original[key].host_competition_cellular_bits,
            host_competition_viral_target=original[key].host_competition_viral_target,
            host_competition_cellular_target=original[key].host_competition_cellular_target,
        )
        for key, row in deep.items()
    } == original
    assert all(row.host_competition_assessed for row in deep.values())
    assert all(row.host_competition_viral_bits is not None for row in deep.values())
    assert all(
        row.host_competition_cellular_bits > original[key].host_competition_cellular_bits for key, row in deep.items()
    )
    assert all(len(row.top10_targets) == 6 for row in deep.values())
