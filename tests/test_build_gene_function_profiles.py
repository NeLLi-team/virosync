"""The profile builder keeps the Pfam and VOG profiles that give a gene a label.

The inputs are tiny profile libraries built in each test with pyhmmer: a Pfam-A
stand-in that holds every Pfam family of the vocabulary plus decoys, and a VOGDB
stand-in with an annotation table that mixes real label hits with real look-alikes.
"""

from __future__ import annotations

from pathlib import Path

import pyhmmer
import pytest

from scripts import build_gene_function_profiles as builder
from virosync.pipeline.phase3.gene_functions import (
    DOMAIN_LABELS,
    read_gene_function_profiles,
)

VOG_ANNOTATIONS_HEADER = "#GroupName\tProteinCount\tSpeciesCount\tFunctionalCategory\tConsensusFunctionalDescription\n"
VOG_DESCRIPTIONS = {
    "VOG00001": "REFSEQ Major capsid protein",
    "VOG00002": "sp|P03709|FI_LAMBD DNA-packaging protein FI",
    "VOG00003": "REFSEQ hypothetical protein",
    "VOG00004": "sp|Q84513|PCNA1_PBCV1 Probable DNA polymerase sliding clamp 1",
    "VOG00005": "sp|P20214|VINT_SSV1 Integrase",
    "VOG00006": "TFIIB: glycosyltransferase family A domain",
}


def _hmm(name: str, description: str = "") -> pyhmmer.plan7.HMM:
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequence = pyhmmer.easel.TextSequence(name=b"seed", sequence="MKVLAAGIVALLAAGCSSQWERTY")
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(sequence.digitize(alphabet), pyhmmer.plan7.Background(alphabet))
    hmm.name = name.encode()
    hmm.description = description.encode() or None
    return hmm


def _write_library(path: Path, hmms: list[pyhmmer.plan7.HMM]) -> Path:
    with path.open("wb") as handle:
        for hmm in hmms:
            hmm.write(handle, binary=False)
    return path


def _hmm_names(path: Path) -> list[str]:
    with pyhmmer.plan7.HMMFile(path) as handle:
        names = [hmm.name.decode() for hmm in handle]
    return names


def _write_pfam(path: Path, *, without: str = "") -> Path:
    names = [name for name in DOMAIN_LABELS if name != without]
    hmms = [_hmm(name, f"{name} family") for name in ["Peptidase_S8", *names, "Pkinase"]]
    return _write_library(path, hmms)


def _write_vog(path: Path, annotations: Path, *, without: str = "") -> Path:
    annotations.write_text(
        VOG_ANNOTATIONS_HEADER + "".join(f"{vog_id}\t5\t5\tXu\t{text}\n" for vog_id, text in VOG_DESCRIPTIONS.items()),
        encoding="utf-8",
    )
    return _write_library(path, [_hmm(vog_id) for vog_id in VOG_DESCRIPTIONS if vog_id != without])


def _request(tmp_path: Path, *, pfam_without: str = "", vog_without: str = "") -> builder.BuildRequest:
    annotations = tmp_path / "vog.annotations.tsv"
    return builder.BuildRequest(
        pfam=_write_pfam(tmp_path / "Pfam-A.hmm", without=pfam_without),
        vog=_write_vog(tmp_path / "vogdb.hmm", annotations, without=vog_without),
        vog_annotations=annotations,
        pfam_release="38.2",
        vog_release="latest-2025-09-30",
        output_dir=tmp_path / "out",
    )


def _build(tmp_path: Path) -> dict[str, builder.Profile]:
    builder.build_gene_function_profiles(_request(tmp_path))
    return read_gene_function_profiles(tmp_path / "out" / "gene_function_labels.tsv")


def test_builder_keeps_the_pfam_families_of_the_vocabulary_and_drops_decoys(tmp_path: Path) -> None:
    profiles = _build(tmp_path)

    pfam = [name for name, profile in profiles.items() if profile.database == "Pfam"]
    assert len(pfam) == 38
    assert "Peptidase_S8" not in pfam
    assert "Pkinase" not in pfam


def test_builder_labels_a_pfam_profile_with_its_release_cutoff_and_description(tmp_path: Path) -> None:
    profile = _build(tmp_path)["DNA_pol_B_2"]

    assert profile == builder.Profile(
        "DNA_pol_B_2",
        "Pfam 38.2",
        "pPolB",
        "replication",
        1e-5,
        "DNA_pol_B_2 family",
    )


def test_builder_keeps_the_vogs_whose_description_names_a_vocabulary_gene(tmp_path: Path) -> None:
    profiles = _build(tmp_path)

    assert [name for name, profile in profiles.items() if profile.database == "VOG"] == [
        "VOG00001",
        "VOG00005",
    ]


@pytest.mark.parametrize("look_alike", ["VOG00002", "VOG00003", "VOG00004", "VOG00006"])
def test_builder_drops_a_vog_look_alike_and_a_hypothetical_protein(tmp_path: Path, look_alike: str) -> None:
    assert look_alike not in _build(tmp_path)


def test_builder_labels_a_vog_profile_with_its_release_cutoff_and_clean_description(tmp_path: Path) -> None:
    profile = _build(tmp_path)["VOG00001"]

    assert profile == builder.Profile(
        "VOG00001",
        "VOGDB latest-2025-09-30",
        "MCP",
        "capsid",
        1e-9,
        "Major capsid protein",
    )


def test_builder_writes_the_hmm_library_in_the_order_of_the_labels_table(tmp_path: Path) -> None:
    profiles = _build(tmp_path)

    assert _hmm_names(tmp_path / "out" / "gene_function.hmm") == list(profiles)


def test_builder_writes_a_hmmer3_text_library(tmp_path: Path) -> None:
    _build(tmp_path)

    first_line = (tmp_path / "out" / "gene_function.hmm").read_text(encoding="utf-8").splitlines()[0]
    assert first_line.startswith("HMMER3/f")


def test_builder_fails_when_a_vocabulary_family_is_missing_from_the_pfam_release(tmp_path: Path) -> None:
    request = _request(tmp_path, pfam_without="Tet_JBP")

    with pytest.raises(ValueError, match="Tet_JBP"):
        builder.build_gene_function_profiles(request)


def test_builder_fails_when_a_selected_vog_has_no_model(tmp_path: Path) -> None:
    request = _request(tmp_path, vog_without="VOG00005")

    with pytest.raises(ValueError, match="VOG00005"):
        builder.build_gene_function_profiles(request)


def test_builder_fails_on_a_vog_annotation_table_with_another_header(tmp_path: Path) -> None:
    request = _request(tmp_path)
    request.vog_annotations.write_text("GroupName\tDescription\nVOG00001\tMajor capsid protein\n", encoding="utf-8")

    with pytest.raises(ValueError, match="annotation"):
        builder.build_gene_function_profiles(request)


def test_main_reports_the_profile_counts_per_database_and_module(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    request = _request(tmp_path)

    exit_code = builder.main(
        [
            "--pfam",
            str(request.pfam),
            "--vog",
            str(request.vog),
            "--vog-annotations",
            str(request.vog_annotations),
            "--pfam-release",
            "38.2",
            "--vog-release",
            "latest-2025-09-30",
            "--output-dir",
            str(request.output_dir),
        ]
    )

    output = capsys.readouterr().out
    assert exit_code == 0
    assert "Pfam profiles: 38" in output
    assert "VOG profiles: 2" in output
    assert "capsid: " in output


def test_main_returns_two_and_names_the_problem_when_the_build_fails(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    request = _request(tmp_path, pfam_without="Tet_JBP")

    exit_code = builder.main(
        [
            "--pfam",
            str(request.pfam),
            "--vog",
            str(request.vog),
            "--vog-annotations",
            str(request.vog_annotations),
            "--pfam-release",
            "38.2",
            "--vog-release",
            "latest-2025-09-30",
            "--output-dir",
            str(request.output_dir),
        ]
    )

    assert exit_code == 2
    assert "Tet_JBP" in capsys.readouterr().err
