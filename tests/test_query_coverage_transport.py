"""Retained raw coverage survives native classification and resume transport."""

from virosync.orchestration._flows.single_genome.phase2_resume_state import _gene_from_document, _gene_to_document
from virosync.pipeline.phase2.boundary_diamond import (
    build_gene_taxonomy_record,
    classify_all_porfs,
    parse_diamond_output,
    pORF,
)


def test_raw_coverage_rank_and_no_hit_roundtrip(tmp_path):
    raw = tmp_path / "hits.tsv"
    raw.write_text("g\tNCLDV__a\t1e-20\t80\t35\t72.4\ng\tEUK__b\t1e-10\t60\t40\t91.2\n")
    genes = classify_all_porfs(
        ["g", "empty"],
        parse_diamond_output(raw),
        {"s": [pORF("g", "s", 10, 100), pORF("empty", "s", 110, 200)]},
        "EUK__",
    )
    assert genes["g"].top10_qcovhsp == [72.4, 91.2]
    assert genes["empty"].top10_qcovhsp == []
    for gene in genes.values():
        restored = _gene_from_document(_gene_to_document(gene, "gene"), "gene")
        assert restored.top10_qcovhsp == gene.top10_qcovhsp
        assert build_gene_taxonomy_record(restored)["top10_qcovhsp"] == gene.top10_qcovhsp


def test_cached_producer_transports_actual_physical_and_alignment_fields(tmp_path):
    from virosync.pipeline.phase3.gene_taxonomy import GeneDiamondHit, materialize_gene_taxonomy_batch_from_cached_hits

    proteome = tmp_path / "proteome.faa"
    proteome.write_text(">scaffold_1 # 11 # 100 # 1 # ID=1_1\nMKK\n>scaffold_2 # 111 # 200 # 1 # ID=1_2\nMKK\n")
    hits = {
        "scaffold_1": [
            GeneDiamondHit("scaffold_1", "NCLDV__a", 1e-20, 80, 35, 72.4),
            GeneDiamondHit("scaffold_1", "EUK__b", 1e-10, 60, 40, 91.2),
        ]
    }
    records, _ = materialize_gene_taxonomy_batch_from_cached_hits(
        [{"eve_id": "eve", "scaffold": "scaffold", "start": 0, "end": 250}], proteome, hits, tmp_path / "out"
    )["eve"]
    hit, nohit = records
    assert (hit.scaffold, hit.start, hit.end) == ("scaffold", 10, 100)
    assert hit.top10_targets == ["NCLDV__a", "EUK__b"]
    assert hit.top10_bitscores == [80, 60]
    assert hit.top10_pidents == [35, 40]
    assert hit.top10_evalues == [1e-20, 1e-10]
    assert hit.top10_qcovhsp == [72.4, 91.2]
    assert (nohit.scaffold, nohit.start, nohit.end) == ("scaffold", 110, 200)
    assert nohit.top10_qcovhsp == []
