"""A run labels the core genes of the EVEs it outputs and changes nothing else.

The resource bundle is a tiny stand-in: one Pfam-style profile built from a seed
protein, one curated model, and the two gene-function files of a v1.2.0 bundle.
"""

from __future__ import annotations

import logging
from pathlib import Path
from unittest.mock import Mock

import pyhmmer
import pytest

from virosync.orchestration._flows.single_genome import orchestrator
from virosync.orchestration.tasks import generate_outputs_task
from virosync.pipeline.phase3 import gene_functions
from virosync.pipeline.phase3.evidence_synthesizer import (
    VerificationResult,
    VerificationStatus,
)
from virosync.pipeline.phase3.gene_functions import (
    GeneFunction,
    GeneFunctionInputs,
    find_gene_function_resources,
    label_core_genes,
    read_gene_functions,
)
from virosync.pipeline.phase3.output_generator import OutputGenerator

SEED_POLYMERASE = (
    "RVRRTWHGTSYGERLFDVCYPRYGYATDCCHIYARMRWHTILSADRKQVDKMITLADWEPELPDAAHHCSPPQDWHKMDLMAQE"
    "FIEACRSGVHRTHFQPEPQHAKWLAHGPYWECFHRKAYMLPDDDHWIAYNNYRFWSWFPGFLIYIHGVHPSYDQCE"
)
UNRELATED = (
    "SKTNFPANSKRYIVAYFRNGMHCWHDTMPDACTIDQREQFVMYVGCVGTDPYQYSSYPVADHKNNPLEKIMNNTWTGAPQCTAIQ"
    "CPHYEVIGDKCQKSNYCTRNHMLRSSIGRVNGHIYAKMGIATCIYENVTFAMKEITTSSNKRGRLYCKMTRHLTK"
)
LABELS_HEADER = "profile\tsource\tlabel\tmodule\tmax_evalue\tdescription\n"
ANNOTATIONS_HEADER = (
    "model_name\tsource\tdescription\tmajority_annotation\tpfam_signature\tpfam_top_domains\tpfam_n_sampled\n"
)
LENGTH_EXCLUSIONS_HEADER = "eve_id\tporf_id\tstatus\treason\tlength_aa\tlimit_aa\n"


def _write_profile_library(path: Path, name: str, seed: str) -> None:
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequence = pyhmmer.easel.TextSequence(name=b"seed", sequence=seed).digitize(alphabet)
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(sequence, pyhmmer.plan7.Background(alphabet))
    hmm.name = name.encode()
    with path.open("wb") as handle:
        hmm.write(handle, binary=False)


def _write_bundle(root: Path, *, gene_function_files: bool = True) -> Path:
    models = root / "models"
    models.mkdir(parents=True)
    hmm_database = models / "combined.hmm"
    hmm_database.write_text("HMMER3/f\nNAME  GVOGm0054\n//\n", encoding="utf-8")
    (models / "model_annotations_with_interpro.tsv").write_text(
        ANNOTATIONS_HEADER + "GVOGm0054\tgvclass_existing\tExisting gvclass model\t\t\t\t0\n",
        encoding="utf-8",
    )
    if gene_function_files:
        _write_profile_library(models / "gene_function.hmm", "DNA_pol_B", SEED_POLYMERASE)
        (models / "gene_function_labels.tsv").write_text(
            LABELS_HEADER + "DNA_pol_B\tPfam 38.2\tPolB\treplication\t1e-05\tDNA polymerase family B\n",
            encoding="utf-8",
        )
    return hmm_database


def _add_decoy_profiles(models: Path, names: list[str]) -> None:
    """Append profiles that match nothing to the gene-function library and its labels table."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    with (models / "gene_function.hmm").open("ab") as handle:
        for name in names:
            sequence = pyhmmer.easel.TextSequence(name=b"seed", sequence=UNRELATED).digitize(alphabet)
            hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(sequence, pyhmmer.plan7.Background(alphabet))
            hmm.name = name.encode()
            hmm.write(handle, binary=False)
    with (models / "gene_function_labels.tsv").open("a", encoding="utf-8") as handle:
        for name in names:
            handle.write(f"{name}\tPfam 38.2\tMTase\tother\t1e-05\tdecoy\n")


def _write_proteome(path: Path) -> Path:
    proteins = {
        "core_pol": SEED_POLYMERASE,
        "core_other": UNRELATED,
        "core_curated": UNRELATED,
        "flank_pol": SEED_POLYMERASE,
    }
    path.write_text(
        "".join(f">{name} # 1 # 480 # 1 # ID=1_1;partial=00\n{sequence}*\n" for name, sequence in proteins.items()),
        encoding="utf-8",
    )
    return path


def _write_model_hits(path: Path) -> Path:
    path.write_text(
        "query\ttarget\tscore\tevalue\ncore_curated\tGVOGm0054\t80.0\t1e-20\ncore_other\tGVOGm0054\t12.0\t1e-2\n",
        encoding="utf-8",
    )
    return path


def _result(eve_id: str, genes: list[tuple[str, bool]]) -> VerificationResult:
    result = VerificationResult(
        eve_id=eve_id,
        scaffold="contig_1",
        start=0,
        end=6001,
        confidence_tier="HIGH",
        status=VerificationStatus.HIGH_CONFIDENCE,
        final_confidence=0.9,
        region_classification="NCLDV",
        likely_family="NCLDV",
        taxonomy_class="NCLDV",
    )
    # The pipeline stores each record as a dict, although the field is annotated with a dataclass.
    records = [
        {
            "porf_id": porf_id,
            "scaffold": "contig_1",
            "porf_start": index * 600,
            "porf_end": index * 600 + 480,
            "top1_prefix": "EUK",
            "top1_target": "EUK__host",
            "top10_prefixes": ["EUK"],
            "is_flanking": is_flanking,
        }
        for index, (porf_id, is_flanking) in enumerate(genes)
    ]
    setattr(result, "gene_taxonomy_records", records)
    return result


def _eve_with_one_flank() -> VerificationResult:
    return _result(
        "EVE_1",
        [("core_pol", False), ("core_other", False), ("core_curated", False), ("flank_pol", True)],
    )


def _inputs(hmm_database: Path, model_hits: Path) -> GeneFunctionInputs:
    return GeneFunctionInputs(hmm_database=hmm_database, model_hits=model_hits, threads=1)


def test_find_gene_function_resources_returns_the_files_next_to_the_hmm_database(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")

    resources = find_gene_function_resources(hmm_database)

    assert resources is not None
    assert resources.profiles_hmm == hmm_database.parent / "gene_function.hmm"
    assert resources.profiles_table == hmm_database.parent / "gene_function_labels.tsv"
    assert resources.model_annotations == hmm_database.parent / "model_annotations_with_interpro.tsv"


def test_find_gene_function_resources_returns_none_for_a_bundle_before_v1_2_0(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync", gene_function_files=False)

    assert find_gene_function_resources(hmm_database) is None


@pytest.mark.parametrize(
    ("present", "missing"),
    [("gene_function.hmm", "gene_function_labels.tsv"), ("gene_function_labels.tsv", "gene_function.hmm")],
)
def test_find_gene_function_resources_rejects_a_bundle_with_only_one_of_the_two_files(
    tmp_path: Path, present: str, missing: str
) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")
    (hmm_database.parent / missing).unlink()
    assert (hmm_database.parent / present).is_file()

    with pytest.raises(FileNotFoundError, match=missing):
        find_gene_function_resources(hmm_database)


def test_write_gene_functions_fails_on_a_bundle_with_only_one_of_the_two_files(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")
    (hmm_database.parent / "gene_function_labels.tsv").unlink()
    generator = _generator(tmp_path, hmm_database)

    with pytest.raises(FileNotFoundError, match="gene_function_labels.tsv"):
        generator.write_gene_functions([_eve_with_one_flank()])


def test_label_core_genes_labels_genes_from_pfam_hits_and_curated_models(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")
    resources = find_gene_function_resources(hmm_database)
    assert resources is not None
    genes = [("EVE_1", "core_pol"), ("EVE_1", "core_other"), ("EVE_1", "core_curated")]

    rows = label_core_genes(
        genes,
        resources,
        model_hits=_write_model_hits(tmp_path / "hits.tsv"),
        proteomes=[_write_proteome(tmp_path / "proteome.fasta")],
        threads=1,
    ).labels

    assert [(row.porf_id, row.label, row.module, row.evidence) for row in rows] == [
        ("core_pol", "PolB", "replication", "Pfam:DNA_pol_B"),
        ("core_curated", "PolB", "replication", "GVOGm0054"),
    ]
    assert rows[0].evalue <= 1e-5
    assert rows[1].evalue == 1e-20


def test_label_core_genes_works_without_a_model_hit_table(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")
    resources = find_gene_function_resources(hmm_database)
    assert resources is not None

    rows = label_core_genes(
        [("EVE_1", "core_pol"), ("EVE_1", "core_curated")],
        resources,
        model_hits=tmp_path / "absent.tsv",
        proteomes=[_write_proteome(tmp_path / "proteome.fasta")],
        threads=1,
    ).labels

    assert [row.porf_id for row in rows] == ["core_pol"]


def test_label_core_genes_keeps_one_row_per_eve_for_a_gene_that_two_eves_share(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")
    resources = find_gene_function_resources(hmm_database)
    assert resources is not None

    rows = label_core_genes(
        [("EVE_1", "core_pol"), ("EVE_2", "core_pol")],
        resources,
        model_hits=tmp_path / "absent.tsv",
        proteomes=[_write_proteome(tmp_path / "proteome.fasta")],
        threads=1,
    ).labels

    assert [(row.eve_id, row.porf_id) for row in rows] == [("EVE_1", "core_pol"), ("EVE_2", "core_pol")]


def test_label_core_genes_leaves_a_gene_missing_from_the_proteome_unlabelled_and_warns(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")
    resources = find_gene_function_resources(hmm_database)
    assert resources is not None

    with caplog.at_level(logging.WARNING, logger="virosync.pipeline.phase3.gene_functions"):
        rows = label_core_genes(
            [("EVE_1", "ghost"), ("EVE_1", "core_pol")],
            resources,
            model_hits=tmp_path / "absent.tsv",
            proteomes=[_write_proteome(tmp_path / "proteome.fasta")],
            threads=1,
        ).labels

    assert [row.porf_id for row in rows] == ["core_pol"]
    assert "ghost" in caplog.text


def _write_crowded_proteome(path: Path, n_other_genes: int) -> Path:
    """Write core_pol and n unrelated proteins, each a rotation of one sequence."""
    records = [("core_pol", SEED_POLYMERASE)]
    records += [(f"other_{index}", UNRELATED[index:] + UNRELATED[:index]) for index in range(n_other_genes)]
    path.write_text(
        "".join(f">{name} # 1 # 480 # 1 # ID=1_1;partial=00\n{sequence}*\n" for name, sequence in records),
        encoding="utf-8",
    )
    return path


def _label_core_pol_with_other_genes(tmp_path: Path, n_other_genes: int) -> GeneFunction:
    """Label core_pol while the run also searches n other core genes."""
    hmm_database = _write_bundle(tmp_path / "virosync")
    resources = find_gene_function_resources(hmm_database)
    assert resources is not None
    genes = [("EVE_1", "core_pol")] + [("EVE_1", f"other_{index}") for index in range(n_other_genes)]
    rows = label_core_genes(
        genes,
        resources,
        model_hits=tmp_path / "absent.tsv",
        proteomes=[_write_crowded_proteome(tmp_path / "proteome.fasta", n_other_genes)],
        threads=1,
    ).labels
    return next(row for row in rows if row.porf_id == "core_pol")


def test_label_core_genes_gives_a_gene_the_same_evalue_whether_it_is_searched_alone_or_with_many_others(
    tmp_path: Path,
) -> None:
    alone = _label_core_pol_with_other_genes(tmp_path / "alone", 0)
    crowded = _label_core_pol_with_other_genes(tmp_path / "crowded", 150)

    assert (alone.label, alone.evalue) == (crowded.label, crowded.evalue)
    assert alone.label == "PolB"


def test_label_core_genes_scales_evalues_by_the_number_of_profiles_in_the_library(tmp_path: Path) -> None:
    proteome = _write_crowded_proteome(tmp_path / "proteome.fasta", 0)
    one_profile = _label_core_pol_with_other_genes(tmp_path / "one", 0)
    bundle = _write_bundle(tmp_path / "three" / "virosync")
    _add_decoy_profiles(bundle.parent, ["Decoy_A", "Decoy_B"])
    resources = find_gene_function_resources(bundle)
    assert resources is not None

    three_profiles = label_core_genes(
        [("EVE_1", "core_pol")],
        resources,
        model_hits=tmp_path / "absent.tsv",
        proteomes=[proteome],
        threads=1,
    ).labels

    assert three_profiles[0].evalue / one_profile.evalue == pytest.approx(3)


def _generator(tmp_path: Path, hmm_database: Path) -> OutputGenerator:
    return OutputGenerator(
        output_dir=tmp_path / "phase3_synthesis",
        proteome_fasta=_write_proteome(tmp_path / "proteome.fasta"),
        gene_function_inputs=_inputs(hmm_database, _write_model_hits(tmp_path / "hits.tsv")),
    )


def test_write_gene_functions_labels_the_core_genes_of_the_output_eves_only(tmp_path: Path) -> None:
    generator = _generator(tmp_path, _write_bundle(tmp_path / "virosync"))

    output_files = generator.write_gene_functions([_eve_with_one_flank()])

    table = output_files["gene_functions"]
    assert table == tmp_path / "phase3_synthesis" / "gene_functions.tsv"
    assert [(row.eve_id, row.porf_id, row.label) for row in read_gene_functions(table)] == [
        ("EVE_1", "core_pol", "PolB"),
        ("EVE_1", "core_curated", "PolB"),
    ]


def test_gene_function_length_exclusions_preserve_labels_and_first_file_sequences(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    generator = _generator(tmp_path, _write_bundle(tmp_path / "virosync"))
    baseline_files = generator.write_gene_functions([_eve_with_one_flank()])
    baseline = {row.porf_id: row for row in read_gene_functions(baseline_files["gene_functions"])}
    proteins = {
        "core_pol": SEED_POLYMERASE,
        "core_other": UNRELATED,
        "core_curated": "M" * 100_001,
        "long_unknown": "M" * 100_002,
    }
    proteome = tmp_path / "proteome.fasta"
    proteome.write_text(
        "".join(f">{name} # 1 # 480 # 1 # ID=1_1;partial=00\n{sequence}*\n" for name, sequence in proteins.items())
    )
    rescue = tmp_path / "phase1" / "frameshift_screening" / "confirmed_frameshift_proteins.faa"
    rescue.parent.mkdir(parents=True)
    rescue.write_text(
        f">core_curated # 1 # 480 # 1 # ID=1_1\n{SEED_POLYMERASE}\n"
        f">long_unknown # 1 # 480 # 1 # ID=1_1\n{SEED_POLYMERASE}\n"
    )
    original_proteome, original_rescue = proteome.read_bytes(), rescue.read_bytes()
    eves = [
        _result("EVE_1", [(name, False) for name in proteins]),
        _result("EVE_2", [("long_unknown", False)]),
    ]

    with caplog.at_level(logging.WARNING, logger="virosync.pipeline.phase3.gene_functions"):
        output_files = generator.generate_all(eves)

    assert read_gene_functions(Path(output_files["gene_functions"])) == [
        baseline["core_pol"],
        baseline["core_curated"],
    ]
    assert Path(output_files["gene_function_length_exclusions"]).read_text() == (
        LENGTH_EXCLUSIONS_HEADER
        + "EVE_1\tcore_curated\tunsearched\tsequence_length_limit\t100001\t100000\n"
        + "EVE_1\tlong_unknown\tunsearched\tsequence_length_limit\t100002\t100000\n"
        + "EVE_2\tlong_unknown\tunsearched\tsequence_length_limit\t100002\t100000\n"
    )
    assert "have no sequence" not in caplog.text
    assert proteome.read_bytes() == original_proteome
    assert rescue.read_bytes() == original_rescue


@pytest.mark.parametrize("length_aa", [100_000, 100_001])
@pytest.mark.parametrize("terminal_stop", ["", "*"])
def test_gene_function_length_limit_follows_stop_removal(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, length_aa: int, terminal_stop: str
) -> None:
    generator = _generator(tmp_path, _write_bundle(tmp_path / "virosync"))
    (tmp_path / "proteome.fasta").write_text(f">core_other\n{'M' * length_aa}{terminal_stop}\n")
    search = Mock(wraps=pyhmmer.hmmsearch)
    monkeypatch.setattr(gene_functions.pyhmmer, "hmmsearch", search)

    output_files = generator.write_gene_functions([_result("EVE_1", [("core_other", False)])])

    assert read_gene_functions(output_files["gene_functions"]) == []
    excluded = "EVE_1\tcore_other\tunsearched\tsequence_length_limit\t100001\t100000\n" if length_aa > 100_000 else ""
    assert output_files["gene_function_length_exclusions"].read_text() == LENGTH_EXCLUSIONS_HEADER + excluded
    assert search.call_count == (0 if excluded else 1)
    if not excluded:
        assert len(search.call_args.args[1][0]) == 100_000
        assert search.call_args.kwargs == {"cpus": 1, "E": 1e-5, "Z": 1, "domZ": 1}


def test_gene_function_length_handling_propagates_unrelated_search_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    generator = _generator(tmp_path, _write_bundle(tmp_path / "virosync"))
    monkeypatch.setattr(gene_functions.pyhmmer, "hmmsearch", Mock(side_effect=ValueError("unrelated profile failure")))

    with pytest.raises(ValueError, match="unrelated profile failure"):
        generator.write_gene_functions([_eve_with_one_flank()])


def _write_rescue_protein(tmp_path: Path, start: int, end: int) -> Path:
    rescue = tmp_path / "phase1" / "frameshift_screening" / "confirmed_frameshift_proteins.faa"
    rescue.parent.mkdir(parents=True)
    rescue.write_text(
        f">contig_1_VSR0123456789abcdef # {start} # {end} # 1 # "
        "ID=0_VSR0123456789abcdef;annotation=frameshift_rescued_domain;model=VS000001;shifts=1;stops=0;literal_stops=0\n"
        f"{SEED_POLYMERASE}\n",
        encoding="utf-8",
    )
    return rescue


def test_write_gene_functions_labels_a_frameshift_rescued_protein_inside_an_eve(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")
    _write_rescue_protein(tmp_path, 2000, 2480)
    generator = _generator(tmp_path, hmm_database)

    output_files = generator.write_gene_functions([_eve_with_one_flank()])

    rows = {(row.eve_id, row.porf_id): row for row in read_gene_functions(output_files["gene_functions"])}
    rescued = rows[("EVE_1", "contig_1_VSR0123456789abcdef")]
    assert (rescued.label, rescued.evidence) == ("PolB", "Pfam:DNA_pol_B")


def test_write_gene_functions_leaves_a_frameshift_rescued_protein_outside_the_eve_unlabelled(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")
    _write_rescue_protein(tmp_path, 7000, 7480)
    generator = _generator(tmp_path, hmm_database)

    output_files = generator.write_gene_functions([_eve_with_one_flank()])

    assert "contig_1_VSR0123456789abcdef" not in {
        row.porf_id for row in read_gene_functions(output_files["gene_functions"])
    }


def _write_empty_rescue_file(tmp_path: Path) -> Path:
    """Phase 1 writes an empty confirmed-rescue file when it confirms no rescue."""
    rescue = tmp_path / "phase1" / "frameshift_screening" / "confirmed_frameshift_proteins.faa"
    rescue.parent.mkdir(parents=True)
    rescue.write_bytes(b"")
    return rescue


def test_label_core_genes_skips_a_proteome_file_without_sequences(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")
    resources = find_gene_function_resources(hmm_database)
    assert resources is not None
    empty = tmp_path / "empty.faa"
    empty.write_bytes(b"")

    rows = label_core_genes(
        [("EVE_1", "core_pol")],
        resources,
        model_hits=tmp_path / "absent.tsv",
        proteomes=[_write_proteome(tmp_path / "proteome.fasta"), empty],
        threads=1,
    ).labels

    assert [row.porf_id for row in rows] == ["core_pol"]


def test_label_core_genes_still_fails_on_a_missing_proteome_file(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")
    resources = find_gene_function_resources(hmm_database)
    assert resources is not None

    with pytest.raises(FileNotFoundError):
        label_core_genes(
            [("EVE_1", "core_pol")],
            resources,
            model_hits=tmp_path / "absent.tsv",
            proteomes=[_write_proteome(tmp_path / "proteome.fasta"), tmp_path / "missing.faa"],
            threads=1,
        )


def test_write_gene_functions_labels_ordinary_genes_when_the_rescue_file_is_empty(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")
    _write_empty_rescue_file(tmp_path)
    generator = _generator(tmp_path, hmm_database)

    output_files = generator.write_gene_functions([_eve_with_one_flank()])

    assert [row.porf_id for row in read_gene_functions(output_files["gene_functions"])] == [
        "core_pol",
        "core_curated",
    ]


def test_generate_all_completes_and_labels_genes_when_the_rescue_file_is_empty(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")
    _write_empty_rescue_file(tmp_path)

    output_files = _generator(tmp_path, hmm_database).generate_all([_eve_with_one_flank()])

    assert [row.porf_id for row in read_gene_functions(Path(output_files["gene_functions"]))] == [
        "core_pol",
        "core_curated",
    ]


def _add_annotation_row(hmm_database: Path, row: str) -> None:
    with (hmm_database.parent / "model_annotations_with_interpro.tsv").open("a", encoding="utf-8") as handle:
        handle.write(row)


def test_label_core_genes_ignores_an_unused_annotation_row_with_a_malformed_sample_count(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")
    _add_annotation_row(hmm_database, "VS000999\tgvclass_existing_og_to_vs\tOG\t\t\tDNA_pol_B:100\tNA\n")
    resources = find_gene_function_resources(hmm_database)
    assert resources is not None

    rows = label_core_genes(
        [("EVE_1", "core_pol"), ("EVE_1", "core_curated")],
        resources,
        model_hits=_write_model_hits(tmp_path / "hits.tsv"),
        proteomes=[_write_proteome(tmp_path / "proteome.fasta")],
        threads=1,
    ).labels

    assert [row.porf_id for row in rows] == ["core_pol", "core_curated"]


def test_label_core_genes_ignores_an_unused_annotation_row_without_a_sample_count_cell(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")
    _add_annotation_row(hmm_database, "VS000999\tgvclass_existing_og_to_vs\tOG\t\t\tDNA_pol_B:100\n")
    resources = find_gene_function_resources(hmm_database)
    assert resources is not None

    rows = label_core_genes(
        [("EVE_1", "core_pol"), ("EVE_1", "core_curated")],
        resources,
        model_hits=_write_model_hits(tmp_path / "hits.tsv"),
        proteomes=[_write_proteome(tmp_path / "proteome.fasta")],
        threads=1,
    ).labels

    assert [row.porf_id for row in rows] == ["core_pol", "core_curated"]


def test_label_core_genes_still_fails_on_a_corrupt_nonempty_proteome_file(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")
    resources = find_gene_function_resources(hmm_database)
    assert resources is not None
    corrupt = tmp_path / "rescue.faa.gz"
    corrupt.write_bytes(b"\x1f\x8b\x08\x00not a gzip stream")

    with pytest.raises((EOFError, OSError, ValueError)):
        label_core_genes(
            [("EVE_1", "core_pol")],
            resources,
            model_hits=tmp_path / "absent.tsv",
            proteomes=[_write_proteome(tmp_path / "proteome.fasta"), corrupt],
            threads=1,
        )


def test_label_core_genes_falls_back_to_profile_evidence_for_an_og_model_with_a_malformed_sample_count(
    tmp_path: Path,
) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")
    _add_annotation_row(hmm_database, "VS000999\tgvclass_existing_og_to_vs\tOG\t\t\tDNA_pol_B:100\tNA\n")
    resources = find_gene_function_resources(hmm_database)
    assert resources is not None
    hits = tmp_path / "hits.tsv"
    hits.write_text("query\ttarget\tscore\tevalue\ncore_pol\tVS000999\t90.0\t1e-30\n", encoding="utf-8")

    rows = label_core_genes(
        [("EVE_1", "core_pol")],
        resources,
        model_hits=hits,
        proteomes=[_write_proteome(tmp_path / "proteome.fasta")],
        threads=1,
    ).labels

    assert [(row.label, row.evidence) for row in rows] == [("PolB", "Pfam:DNA_pol_B")]


def test_write_gene_functions_writes_a_header_when_no_core_gene_gets_a_label(tmp_path: Path) -> None:
    generator = _generator(tmp_path, _write_bundle(tmp_path / "virosync"))

    output_files = generator.write_gene_functions([_result("EVE_2", [("core_other", False)])])

    assert read_gene_functions(output_files["gene_functions"]) == []
    assert output_files["gene_functions"].read_text(encoding="utf-8").startswith("eve_id\tporf_id\t")


def test_write_gene_functions_skips_a_bundle_before_v1_2_0_with_one_info_line(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    generator = _generator(tmp_path, _write_bundle(tmp_path / "virosync", gene_function_files=False))

    with caplog.at_level(logging.INFO, logger="virosync.pipeline.phase3.output_generator"):
        output_files = generator.write_gene_functions([_eve_with_one_flank()])

    assert output_files == {}
    assert not (tmp_path / "phase3_synthesis" / "gene_functions.tsv").exists()
    assert not (tmp_path / "phase3_synthesis" / "gene_function_length_exclusions.tsv").exists()
    messages = [record.getMessage() for record in caplog.records if "gene-function" in record.getMessage()]
    assert len(messages) == 1
    assert caplog.records[-1].levelno == logging.INFO


def test_write_gene_functions_does_nothing_without_gene_function_inputs(tmp_path: Path) -> None:
    generator = OutputGenerator(output_dir=tmp_path / "out", proteome_fasta=_write_proteome(tmp_path / "p.fasta"))

    assert generator.write_gene_functions([_eve_with_one_flank()]) == {}
    assert not (tmp_path / "out" / "gene_functions.tsv").exists()
    assert not (tmp_path / "out" / "gene_function_length_exclusions.tsv").exists()


@pytest.mark.parametrize("empty", [False, True])
def test_generate_all_registers_the_gene_function_table_among_the_output_files(tmp_path: Path, empty: bool) -> None:
    generator = _generator(tmp_path, _write_bundle(tmp_path / "virosync"))

    output_files = generator.generate_all([] if empty else [_eve_with_one_flank()])

    expected = set() if empty else {"gene_functions", "gene_function_length_exclusions"}
    assert set(output_files) & {"gene_functions", "gene_function_length_exclusions"} == expected
    assert (tmp_path / "phase3_synthesis" / "gene_functions.tsv").exists() is not empty
    assert (tmp_path / "phase3_synthesis" / "gene_function_length_exclusions.tsv").exists() is not empty
    if not empty:
        assert Path(output_files["gene_function_length_exclusions"]).read_text() == LENGTH_EXCLUSIONS_HEADER
        assert "gene_taxonomy_all" in output_files


def test_generate_all_keeps_every_other_output_byte_identical(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")

    labelled = _generator(tmp_path, hmm_database).generate_all([_eve_with_one_flank()])
    unlabelled = OutputGenerator(output_dir=tmp_path / "plain").generate_all([_eve_with_one_flank()])

    assert set(labelled) - set(unlabelled) == {"gene_functions", "gene_function_length_exclusions"}
    for key in ("predictions_tsv", "predictions_detailed_tsv", "gene_taxonomy_all", "evidence_json"):
        assert Path(labelled[key]).read_bytes() == Path(unlabelled[key]).read_bytes()


def test_generate_outputs_task_passes_the_gene_function_inputs_to_the_generator(tmp_path: Path) -> None:
    hmm_database = _write_bundle(tmp_path / "virosync")

    output_files = generate_outputs_task(
        [_eve_with_one_flank()],
        tmp_path / "phase3_synthesis",
        proteome_path=_write_proteome(tmp_path / "proteome.fasta"),
        gene_function_inputs=_inputs(hmm_database, _write_model_hits(tmp_path / "hits.tsv")),
    )

    assert [
        (row.eve_id, row.porf_id, row.label) for row in read_gene_functions(Path(output_files["gene_functions"]))
    ] == [
        ("EVE_1", "core_pol", "PolB"),
        ("EVE_1", "core_curated", "PolB"),
    ]


def test_the_run_reads_the_phase1_validated_hits_beside_its_hmm_database(tmp_path: Path) -> None:
    hmm_database = tmp_path / "virosync" / "models" / "combined.hmm"

    inputs = orchestrator._gene_function_inputs(hmm_database, tmp_path / "run", threads=8)

    assert inputs == GeneFunctionInputs(
        hmm_database=hmm_database,
        model_hits=tmp_path / "run" / "phase1" / "hmm" / "hmm_hits_validated.tsv",
        threads=8,
    )


def test_a_run_without_an_hmm_database_labels_no_genes(tmp_path: Path) -> None:
    assert orchestrator._gene_function_inputs(None, tmp_path / "run", threads=8) is None
