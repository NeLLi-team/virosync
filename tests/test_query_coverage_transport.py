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
