"""Gene-function labels name the genes inside accepted EVEs without changing a call.

A label uses the gene vocabulary of Bellas and Sommaruga (2026, Fig. 4). It comes
from the gene's best-scoring curated ViroSync model, or from a renamed OG-backed
model whose top Pfam domain holds in at least half of its sampled members, otherwise
from the gene's strongest Pfam hit, otherwise from the gene's strongest VOG hit.
The rules and the VOG look-alikes that must stay unlabelled come from the benchmark
repository's real-synteny figure.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from virosync.pipeline.phase3.gene_functions import (
    DESCRIPTION_RULES,
    DOMAIN_LABELS,
    LABEL_MODULES,
    MODEL_LABELS,
    VOG_DESCRIPTION_RULES,
    GeneFunction,
    ModelAnnotation,
    ModelHit,
    Profile,
    ProfileHit,
    label_gene,
    read_gene_function_profiles,
    read_gene_functions,
    read_model_annotations,
    vog_label,
    write_gene_function_profiles,
    write_gene_functions,
)

PFAM = "Pfam 38.2"
VOG = "VOGDB latest-2025-09-30"


def _pfam_profile(name: str, label: str, module: str) -> Profile:
    return Profile(name, PFAM, label, module, 1e-5, f"{name} family")


def _vog_profile(vog_id: str, label: str, module: str, description: str) -> Profile:
    return Profile(vog_id, VOG, label, module, 1e-9, description)


def _profiles() -> dict[str, Profile]:
    profiles = [
        _pfam_profile("DNA_methylase", "MTase", "other"),
        _pfam_profile("DNA_pol_B_2", "pPolB", "replication"),
        _pfam_profile("Phage_integrase", "YR", "integration"),
        _pfam_profile("rve", "RVE", "integration"),
        _pfam_profile("Pox_VLTF3", "VLTF3", "mriya"),
        _pfam_profile("TFIIB", "TFIIB", "hallmark"),
        _pfam_profile("Patatin", "Lipase", "other"),
        _vog_profile("VOG00010", "MCP", "capsid", "Major capsid protein"),
        _vog_profile("VOG00011", "mCP", "capsid", "Minor capsid protein P2"),
        _vog_profile("VOG00012", "PolB", "replication", "DNA polymerase type-B"),
    ]
    return {profile.profile: profile for profile in profiles}


def _curated(description: str = "", majority: str = "") -> ModelAnnotation:
    return ModelAnnotation("curated", description, majority, "", 0)


def _og(source: str, top_domains: str, n_sampled: int) -> ModelAnnotation:
    return ModelAnnotation(source, "Existing NCLDV gvclass OG marker", "", top_domains, n_sampled)


def _annotations() -> dict[str, ModelAnnotation]:
    return {
        "VP_Penton_1": _curated(
            "Virophage Penton Protein",
            "Minor capsid protein V18/19-like, jelly roll domain superfamily",
        ),
        "EVE_CAPSID_NOVEL_0024": _curated("Putative penton protein"),
        "GVOGm0003": _curated("Existing gvclass model"),
        "Mirus_Portal": _curated("Portal Protein", "consensus disorder prediction"),
        "Mirus_Triplex1": _curated("Triplex Protein 1"),
        "ATPase": _curated("Existing gvclass model"),
        "GVOGm0054": _curated("Existing gvclass model"),
        "GVOGm0023": _curated("Existing gvclass model"),
        "VLTF3": _curated("Existing gvclass model"),
        "VS000355": _og("gvclass_existing_og_to_vs", "Peptidase_S74:150", 150),
        "VS000100": _og("gvclass_existing_og_to_vs", "Pox_VLTF3:90;AAA:10", 150),
        "VS000101": _og("OG_backed_relaxed_low_overlap_NCLDV", "DNA_pol_B:75", 150),
        "VS000102": _og("OG_backed_relaxed_low_overlap_VP", "DNA_pol_B:74", 150),
        "VS000103": _og("gvclass_existing", "DNA_pol_B:150", 150),
        "VS000104": _og("gvclass_existing_og_to_vs", "", 150),
        "VS000105": _og("gvclass_existing_og_to_vs", "rve:1", 0),
        "VS000106": _og("OG_backed_relaxed_low_overlap_MIRUS", "AAA:80;Pox_VLTF3:70", 150),
        "VS000107": _og("OG_backed_relaxed_low_overlap_PLV", "Pox_A32:80;AAA:70", 150),
        "VS000108": _og("gvclass_existing_og_to_vs", "DNA_pol_B:NA", 0),
        "VS000109": _og("gvclass_existing_og_to_vs", "DNA_pol_B:NA", 100),
        "VS000110": _og("gvclass_existing_og_to_vs", "DNA_pol_B:101", 100),
        "VS000111": _og("gvclass_existing_og_to_vs", "AAA:80;Pox_VLTF3:NA", 100),
        "VS000112": _og("gvclass_existing_og_to_vs", "DNA_pol_B:-5", 100),
        "VS000113": _og("gvclass_existing_og_to_vs", "DNA_pol_B", 100),
        "VS000114": _og("gvclass_existing_og_to_vs", "DNA_pol_B:100", 100),
    }


def _label(
    model_hits: list[ModelHit],
    profile_hits: list[ProfileHit],
) -> GeneFunction | None:
    return label_gene(
        "EVE_1",
        "EVE_1_3",
        model_hits=model_hits,
        profile_hits=profile_hits,
        model_annotations=_annotations(),
        profiles=_profiles(),
    )


def _model(name: str, evalue: float = 1e-30, score: float = 100.0) -> list[ModelHit]:
    return [ModelHit(name, score, evalue)]


@pytest.mark.parametrize(
    ("model", "label", "module"),
    [
        ("VP_Penton_1", "mCP", "capsid"),
        ("EVE_CAPSID_NOVEL_0024", "mCP", "capsid"),
        ("GVOGm0003", "MCP", "capsid"),
        ("Mirus_Portal", "Portal", "capsid"),
        ("Mirus_Triplex1", "Triplex", "capsid"),
        ("ATPase", "ATPase", "capsid"),
        ("GVOGm0054", "PolB", "replication"),
        ("GVOGm0023", "RNAPL", "hallmark"),
        ("VLTF3", "VLTF3", "mriya"),
    ],
)
def test_label_gene_names_a_curated_model_with_the_model_as_evidence(model: str, label: str, module: str) -> None:
    assert _label(_model(model), []) == GeneFunction("EVE_1", "EVE_1_3", label, module, model, 1e-30)


def test_label_gene_takes_the_best_scoring_model_hit_below_the_evalue_cutoff() -> None:
    hits = [
        ModelHit("GVOGm0003", 40.0, 1e-30),
        ModelHit("GVOGm0054", 60.0, 1e-8),
        ModelHit("ATPase", 99.0, 1e-2),
    ]

    assert _label(hits, []) == GeneFunction("EVE_1", "EVE_1_3", "PolB", "replication", "GVOGm0054", 1e-8)


def test_label_gene_ignores_a_model_hit_above_the_evalue_cutoff() -> None:
    assert _label(_model("GVOGm0003", evalue=1e-4), []) is None


def test_label_gene_uses_only_the_best_model_even_if_a_weaker_one_is_curated() -> None:
    hits = [ModelHit("VS000355", 90.0, 1e-30), ModelHit("GVOGm0003", 50.0, 1e-30)]

    assert _label(hits, []) is None


def test_label_gene_ignores_a_curated_model_missing_from_the_annotation_table() -> None:
    assert _label(_model("VP_Unknown_9"), []) is None


def test_label_gene_ignores_the_top_domain_of_a_renamed_og_model_below_the_floor() -> None:
    assert _label(_model("VS000102"), []) is None


def test_label_gene_ignores_an_og_model_whose_top_domain_is_outside_the_vocabulary() -> None:
    assert _label(_model("VS000355"), []) is None


def test_label_gene_ignores_a_vocabulary_domain_below_the_top_domain_of_an_og_model() -> None:
    assert _label(_model("VS000106"), []) is None


def test_label_gene_ignores_an_og_model_without_sampled_members() -> None:
    assert _label(_model("VS000105"), []) is None


def test_label_gene_ignores_an_og_model_without_pfam_domains() -> None:
    assert _label(_model("VS000104"), []) is None


def test_label_gene_ignores_the_pfam_signature_of_a_model_that_is_not_og_backed() -> None:
    assert _label(_model("VS000103"), []) is None


def test_label_gene_trusts_an_og_model_whose_top_domain_holds_in_most_members() -> None:
    assert _label(_model("VS000100", evalue=1e-12), []) == GeneFunction(
        "EVE_1",
        "EVE_1_3",
        "VLTF3",
        "mriya",
        "VS000100 (OG; Pox_VLTF3 in 90/150 members)",
        1e-12,
    )


def test_label_gene_trusts_an_og_model_at_exactly_half_of_its_members() -> None:
    gene_function = _label(_model("VS000101"), [])

    assert gene_function is not None
    assert gene_function.label == "PolB"
    assert gene_function.evidence == "VS000101 (OG; DNA_pol_B in 75/150 members)"


def test_label_gene_trusts_a_relaxed_og_model_for_its_vocabulary_top_domain() -> None:
    gene_function = _label(_model("VS000107"), [])

    assert gene_function is not None
    assert gene_function.label == "ATPase"
    assert gene_function.module == "capsid"


@pytest.mark.parametrize(
    ("domain", "label", "module"),
    [
        ("DNA_methylase", "MTase", "other"),
        ("DNA_pol_B_2", "pPolB", "replication"),
        ("Phage_integrase", "YR", "integration"),
        ("rve", "RVE", "integration"),
        ("Pox_VLTF3", "VLTF3", "mriya"),
        ("TFIIB", "TFIIB", "hallmark"),
        ("Patatin", "Lipase", "other"),
    ],
)
def test_label_gene_falls_back_to_a_pfam_hit_of_the_profile_set(domain: str, label: str, module: str) -> None:
    assert _label(_model("VS000355"), [ProfileHit(domain, 1e-20)]) == GeneFunction(
        "EVE_1", "EVE_1_3", label, module, f"Pfam:{domain}", 1e-20
    )


def test_label_gene_labels_a_gene_without_a_model_hit_from_its_pfam_hit() -> None:
    gene_function = _label([], [ProfileHit("DNA_methylase", 1e-6)])

    assert gene_function is not None
    assert gene_function.evidence == "Pfam:DNA_methylase"


def test_label_gene_ignores_a_pfam_hit_above_the_profile_cutoff() -> None:
    assert _label([], [ProfileHit("DNA_methylase", 1e-4)]) is None


def test_label_gene_accepts_a_pfam_hit_at_the_profile_cutoff() -> None:
    assert _label([], [ProfileHit("DNA_methylase", 1e-5)]) is not None


def test_label_gene_takes_the_strongest_pfam_hit_first() -> None:
    hits = [ProfileHit("Patatin", 1e-8), ProfileHit("DNA_pol_B_2", 1e-40)]

    gene_function = _label([], hits)

    assert gene_function is not None
    assert gene_function.label == "pPolB"
    assert gene_function.evalue == 1e-40


def test_label_gene_skips_a_stronger_hit_without_a_profile_label() -> None:
    hits = [ProfileHit("Peptidase_S8", 1e-80), ProfileHit("Patatin", 1e-8)]

    gene_function = _label([], hits)

    assert gene_function is not None
    assert gene_function.evidence == "Pfam:Patatin"


def test_label_gene_labels_a_gene_from_a_vog_hit_with_the_vog_as_evidence() -> None:
    assert _label([], [ProfileHit("VOG00011", 1e-30)]) == GeneFunction(
        "EVE_1", "EVE_1_3", "mCP", "capsid", "VOG:VOG00011", 1e-30
    )


def test_label_gene_ignores_a_vog_hit_above_the_vog_cutoff() -> None:
    assert _label([], [ProfileHit("VOG00010", 1e-8)]) is None


def test_label_gene_accepts_a_vog_hit_at_the_vog_cutoff() -> None:
    assert _label([], [ProfileHit("VOG00010", 1e-9)]) is not None


def test_label_gene_takes_the_strongest_vog_hit_first() -> None:
    hits = [ProfileHit("VOG00010", 1e-12), ProfileHit("VOG00012", 1e-90)]

    gene_function = _label([], hits)

    assert gene_function is not None
    assert gene_function.evidence == "VOG:VOG00012"


def test_label_gene_prefers_a_pfam_hit_to_a_stronger_vog_hit() -> None:
    hits = [ProfileHit("VOG00012", 1e-90), ProfileHit("DNA_pol_B_2", 1e-6)]

    gene_function = _label([], hits)

    assert gene_function is not None
    assert gene_function.evidence == "Pfam:DNA_pol_B_2"


def test_label_gene_prefers_a_curated_model_to_a_pfam_hit() -> None:
    gene_function = _label(_model("GVOGm0054"), [ProfileHit("DNA_pol_B_2", 1e-60)])

    assert gene_function is not None
    assert gene_function.evidence == "GVOGm0054"


def test_label_gene_prefers_a_trusted_og_model_to_a_pfam_hit() -> None:
    gene_function = _label(_model("VS000100"), [ProfileHit("DNA_methylase", 1e-60)])

    assert gene_function is not None
    assert gene_function.label == "VLTF3"


def test_label_gene_falls_back_from_an_untrusted_og_model_to_a_vog_hit() -> None:
    gene_function = _label(_model("VS000102"), [ProfileHit("VOG00010", 1e-50)])

    assert gene_function is not None
    assert gene_function.evidence == "VOG:VOG00010"


def test_label_gene_leaves_a_gene_without_evidence_unlabelled() -> None:
    assert _label([], [ProfileHit("Peptidase_S8", 1e-50)]) is None


@pytest.mark.parametrize(
    ("description", "expected"),
    [
        ("sp|Q84656|P2_PBCV1 Minor capsid protein P2", "mCP"),
        ("sp|A7U6E9|MCP_POV01 Major capsid protein", "MCP"),
        ("A32-like packaging ATPase", "ATPase"),
        ("sp|A0A7H0DND1|PG160_MONPV DNA packaging protein OPG160", "ATPase"),
        ("sp|P03709|FI_LAMBD DNA-packaging protein FI", ""),
        ("sp|Q38087|DPOL_BPR69 DNA-directed DNA polymerase", "PolB"),
        ("highly derived D5-like helicase-primase", "SF3"),
        ("helicase-primase primase subunit", ""),
        ("tyrosine-type recombinase/integrase", "YR"),
        ("site specific recombinase, tyrosine", "YR"),
        ("sp|P20214|VINT_SSV1 Integrase", "Int"),
        ("putative recombination protein Bet, RecT recombinase", ""),
        ("Viral late gene transcription factor 2", "VLTF2"),
        ("Viral late gene transcription factor 3", "VLTF3"),
        ("transcription initiation factor TFIIB", "TFIIB"),
        ("sp|Q65160|TFIIB_ASFB7 Initiation factor TFIIB homolog", "TFIIB"),
        ("TFIIB: glycosyltransferase family A domain", ""),
        ("zf-TFIIB domain-containing protein", ""),
        ("sp|P34877|MTSA_LACLC Type II methyltransferase M1.ScrFI", "MTase"),
        ("sp|Q6QGD5|PORTL_BPT5 Portal protein", "Portal"),
        ("terminase large subunit", "Terminase"),
        ("sp|Q84513|PCNA1_PBCV1 Probable DNA polymerase sliding clamp 1", ""),
        ("REFSEQ hypothetical protein", ""),
    ],
)
def test_vog_label_follows_the_consensus_description(description: str, expected: str) -> None:
    assert vog_label(description) == expected


def test_label_modules_group_the_vocabulary_into_six_modules() -> None:
    assert set(LABEL_MODULES.values()) == {
        "capsid",
        "replication",
        "integration",
        "mriya",
        "other",
        "hallmark",
    }


@pytest.mark.parametrize(
    ("label", "module"),
    [
        ("MCP", "capsid"),
        ("mCP", "capsid"),
        ("ATPase", "capsid"),
        ("Protease", "capsid"),
        ("Portal", "capsid"),
        ("Terminase", "capsid"),
        ("Triplex", "capsid"),
        ("PolB", "replication"),
        ("pPolB", "replication"),
        ("PolA", "replication"),
        ("SF1", "replication"),
        ("SF3", "replication"),
        ("AEP", "replication"),
        ("YR", "integration"),
        ("RVE", "integration"),
        ("Int", "integration"),
        ("VLTF2", "mriya"),
        ("VLTF3", "mriya"),
        ("HUH", "mriya"),
        ("MTase", "other"),
        ("Tet", "other"),
        ("SLATT", "other"),
        ("Lipase", "other"),
        ("RNAPL", "hallmark"),
        ("TFIIB", "hallmark"),
    ],
)
def test_label_modules_place_each_vocabulary_label_in_its_module(label: str, module: str) -> None:
    assert LABEL_MODULES[label] == module


def test_every_rule_table_label_has_a_module() -> None:
    rule_labels = {
        *MODEL_LABELS.values(),
        *DOMAIN_LABELS.values(),
        *(label for _, label in DESCRIPTION_RULES),
        *(label for _, label in VOG_DESCRIPTION_RULES),
    }

    assert rule_labels <= set(LABEL_MODULES)


def test_profiles_round_trip_through_the_labels_table(tmp_path: Path) -> None:
    path = tmp_path / "gene_function_labels.tsv"
    profiles = list(_profiles().values())

    write_gene_function_profiles(profiles, path)

    assert list(read_gene_function_profiles(path).values()) == profiles


def test_profiles_are_written_with_the_documented_columns(tmp_path: Path) -> None:
    path = tmp_path / "gene_function_labels.tsv"

    write_gene_function_profiles([_vog_profile("VOG00010", "MCP", "capsid", "Major capsid protein")], path)

    assert path.read_text(encoding="utf-8") == (
        "profile\tsource\tlabel\tmodule\tmax_evalue\tdescription\n"
        "VOG00010\tVOGDB latest-2025-09-30\tMCP\tcapsid\t1e-09\tMajor capsid protein\n"
    )


def test_read_profiles_rejects_a_label_outside_the_vocabulary(tmp_path: Path) -> None:
    path = tmp_path / "gene_function_labels.tsv"
    path.write_text(
        "profile\tsource\tlabel\tmodule\tmax_evalue\tdescription\nrve\tPfam 38.2\tIntegrase\tintegration\t1e-05\tx\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="label 'Integrase'"):
        read_gene_function_profiles(path)


def test_read_profiles_rejects_a_module_that_disagrees_with_the_label(
    tmp_path: Path,
) -> None:
    path = tmp_path / "gene_function_labels.tsv"
    path.write_text(
        "profile\tsource\tlabel\tmodule\tmax_evalue\tdescription\nrve\tPfam 38.2\tRVE\tcapsid\t1e-05\tx\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="module 'capsid'"):
        read_gene_function_profiles(path)


def test_read_profiles_rejects_a_source_that_is_neither_pfam_nor_vogdb(
    tmp_path: Path,
) -> None:
    path = tmp_path / "gene_function_labels.tsv"
    path.write_text(
        "profile\tsource\tlabel\tmodule\tmax_evalue\tdescription\nrve\tInterPro 100\tRVE\tintegration\t1e-05\tx\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="source 'InterPro 100'"):
        read_gene_function_profiles(path)


def test_read_profiles_rejects_a_table_with_other_columns(tmp_path: Path) -> None:
    path = tmp_path / "gene_function_labels.tsv"
    path.write_text("profile\tlabel\nrve\tRVE\n", encoding="utf-8")

    with pytest.raises(ValueError, match="columns"):
        read_gene_function_profiles(path)


def test_read_profiles_rejects_a_repeated_profile(tmp_path: Path) -> None:
    path = tmp_path / "gene_function_labels.tsv"
    row = "rve\tPfam 38.2\tRVE\tintegration\t1e-05\tx\n"
    path.write_text(
        "profile\tsource\tlabel\tmodule\tmax_evalue\tdescription\n" + row + row,
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="repeats profile 'rve'"):
        read_gene_function_profiles(path)


def test_gene_functions_are_written_with_the_documented_columns(tmp_path: Path) -> None:
    path = tmp_path / "gene_functions.tsv"

    write_gene_functions(
        [GeneFunction("EVE_1", "EVE_1_3", "PolB", "replication", "GVOGm0054", 2.5e-12)],
        path,
    )

    assert path.read_text(encoding="utf-8") == (
        "eve_id\tporf_id\tlabel\tmodule\tevidence\tevalue\nEVE_1\tEVE_1_3\tPolB\treplication\tGVOGm0054\t2.5e-12\n"
    )


def test_gene_functions_round_trip_through_the_output_table(tmp_path: Path) -> None:
    path = tmp_path / "gene_functions.tsv"
    rows = [
        GeneFunction("EVE_1", "EVE_1_3", "PolB", "replication", "GVOGm0054", 2.5e-12),
        GeneFunction(
            "EVE_2",
            "EVE_2_9",
            "VLTF3",
            "mriya",
            "VS000100 (OG; Pox_VLTF3 in 90/150 members)",
            1e-30,
        ),
    ]

    write_gene_functions(rows, path)

    assert read_gene_functions(path) == rows


def test_an_empty_output_table_keeps_its_header(tmp_path: Path) -> None:
    path = tmp_path / "gene_functions.tsv"

    write_gene_functions([], path)

    assert read_gene_functions(path) == []
    assert path.read_text(encoding="utf-8").startswith("eve_id\tporf_id\t")


def test_read_model_annotations_keeps_the_fields_that_decide_a_label(
    tmp_path: Path,
) -> None:
    path = tmp_path / "model_annotations_with_interpro.tsv"
    path.write_text(
        "model_name\tsource\tdescription\tmajority_annotation\tpfam_signature\t"
        "pfam_top_domains\tpfam_n_sampled\tpfam_n_with_pfam\tsource_marker\t"
        "source_scope\traw_name\n"
        "VS000100\tgvclass_existing_og_to_vs\tOG marker\tmajority\tPox_VLTF3\t"
        "Pox_VLTF3:90;AAA:10\t150\t100\t\t\t\n"
        "Mirus_Portal\tMirusviricota_Delmont2025\tPortal Protein\t\t\t\t150\t0\t\t\t\n",
        encoding="utf-8",
    )

    annotations = read_model_annotations(path)

    assert annotations == {
        "VS000100": ModelAnnotation(
            "gvclass_existing_og_to_vs",
            "OG marker",
            "majority",
            "Pox_VLTF3:90;AAA:10",
            150,
        ),
        "Mirus_Portal": ModelAnnotation("Mirusviricota_Delmont2025", "Portal Protein", "", "", 150),
    }


@pytest.mark.parametrize(
    "model",
    [
        "VS000108",  # malformed count and no sampled members
        "VS000109",  # malformed count
        "VS000110",  # more members with the domain than were sampled
        "VS000111",  # one malformed entry spoils the whole listing
        "VS000112",  # negative count
        "VS000113",  # entry without a count
    ],
)
def test_label_gene_distrusts_an_og_model_with_invalid_domain_counts_and_uses_profile_evidence(model: str) -> None:
    gene_function = _label(_model(model), [ProfileHit("DNA_methylase", 1e-20)])

    assert gene_function is not None
    assert gene_function.evidence == "Pfam:DNA_methylase"


def test_label_gene_trusts_an_og_model_when_every_sampled_member_has_the_domain() -> None:
    gene_function = _label(_model("VS000114"), [])

    assert gene_function is not None
    assert gene_function.evidence == "VS000114 (OG; DNA_pol_B in 100/100 members)"


def test_read_model_annotations_treats_a_malformed_sample_count_as_no_sampled_members(tmp_path: Path) -> None:
    path = tmp_path / "model_annotations_with_interpro.tsv"
    path.write_text(
        "model_name\tsource\tdescription\tmajority_annotation\tpfam_signature\t"
        "pfam_top_domains\tpfam_n_sampled\n"
        "VS000001\tgvclass_existing_og_to_vs\tOG\t\t\tDNA_pol_B:100\tNA\n",
        encoding="utf-8",
    )

    assert read_model_annotations(path)["VS000001"].pfam_n_sampled == 0


@pytest.mark.parametrize(
    ("row_end", "expected"),
    [("DNA_pol_B:100\n", 0), ("DNA_pol_B:100\t 100 \n", 100)],
    ids=["missing-cell", "padded-count"],
)
def test_read_model_annotations_reads_a_missing_or_padded_sample_count(
    tmp_path: Path, row_end: str, expected: int
) -> None:
    path = tmp_path / "model_annotations_with_interpro.tsv"
    path.write_text(
        "model_name\tsource\tdescription\tmajority_annotation\tpfam_signature\t"
        "pfam_top_domains\tpfam_n_sampled\n"
        f"VS000001\tgvclass_existing_og_to_vs\tOG\t\t\t{row_end}",
        encoding="utf-8",
    )

    assert read_model_annotations(path)["VS000001"].pfam_n_sampled == expected
