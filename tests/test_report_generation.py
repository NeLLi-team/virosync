from __future__ import annotations

import csv
import json
import os
import re
import struct
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from types import SimpleNamespace

import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from matplotlib.lines import Line2D

from virosync.pipeline.phase3.eve_ani_clustering import (
    MIN_CLUSTER_ALIGNED_FRACTION,
    MIN_CLUSTER_ANI,
)
from virosync.report import generate as report_generate
from virosync.report.generate import generate_eve_report
from virosync.report.graphviz_runtime import (
    ANI_NETWORK_ENGINE,
    ANI_NETWORK_RENDER_ATTR,
)


def test_generate_eve_report_writes_jupyter_notebook(
    monkeypatch,
    tmp_path: Path,
) -> None:
    captured: dict[str, object] = {}

    def _execute_notebook(
        input_path: str,
        output_path: str,
        parameters: dict[str, str],
        cwd: str,
        kernel_name: str,
        progress_bar: bool,
    ) -> None:
        # Source is a real notebook materialized from the jupytext .py source.
        assert input_path.endswith(".ipynb")
        assert Path(cwd) == tmp_path
        captured["parameters"] = parameters
        captured["kernel_name"] = kernel_name
        captured["progress_bar"] = progress_bar
        Path(output_path).write_text("{}")

    monkeypatch.setitem(
        sys.modules,
        "papermill",
        SimpleNamespace(execute_notebook=_execute_notebook),
    )
    monkeypatch.setattr(report_generate.logger, "isEnabledFor", lambda _level: False)

    report_paths = generate_eve_report(
        output_dir=tmp_path,
        genome_id="demo",
    )

    assert report_paths.jupyter == tmp_path / "notebooks" / "jupyter" / "eve_analysis.ipynb"
    assert report_paths.jupyter.exists()
    assert captured["parameters"]["GENOME_ID"] == "demo"
    assert captured["parameters"]["RESULTS_DIR"] == str(tmp_path)
    assert captured["kernel_name"] == "python3"
    assert captured["progress_bar"] is False
    # marimo output must no longer be produced
    assert not (tmp_path / "notebooks" / "marimo").exists()
    assert not hasattr(report_paths, "marimo")


def test_generate_eve_report_requires_papermill_success(
    monkeypatch,
    tmp_path: Path,
) -> None:
    def _execute_notebook(*args, **kwargs) -> None:
        raise RuntimeError("papermill failure")

    monkeypatch.setitem(
        sys.modules,
        "papermill",
        SimpleNamespace(execute_notebook=_execute_notebook),
    )

    with pytest.raises(RuntimeError, match="papermill failure"):
        generate_eve_report(output_dir=tmp_path, genome_id="demo")

    assert not (tmp_path / "notebooks" / "marimo").exists()


def test_notebook_source_is_valid_jupytext_with_parameters_cell() -> None:
    """The single .py source must round-trip through jupytext and expose the
    papermill 'parameters' cell, otherwise parameter injection silently breaks."""
    jupytext = pytest.importorskip("jupytext")

    assert report_generate._SOURCE.exists()
    notebook = jupytext.read(report_generate._SOURCE)
    tagged = [cell for cell in notebook.cells if "parameters" in cell.get("metadata", {}).get("tags", [])]
    assert tagged, "notebook source is missing a 'parameters'-tagged cell"


def _notebook_cell(needle: str) -> str:
    jupytext = pytest.importorskip("jupytext")
    notebook = jupytext.read(report_generate._SOURCE)
    cells = [cell.source for cell in notebook.cells if needle in cell.source]
    assert len(cells) == 1, f"{needle!r} matched {len(cells)} cells"
    # Keep report-only DIAMOND searches inside the test runner's single-CPU budget.
    return cells[0].replace("SYNTVIZ_DIAMOND_THREADS = 4", "SYNTVIZ_DIAMOND_THREADS = 1")


def _notebook_code() -> str:
    jupytext = pytest.importorskip("jupytext")
    notebook = jupytext.read(report_generate._SOURCE)
    return "\n".join(cell.source for cell in notebook.cells if cell.cell_type == "code")


def test_notebook_repeat_evidence_displays_assessment_without_promoting_candidates() -> None:
    """The report retains a negative result and positive incomplete evidence separately."""
    predictions = pd.DataFrame(
        [
            {
                "eve_id": "EVE_negative",
                "repeat_assessed_start": 0,
                "repeat_assessed_end": 1000,
                "repeat_left_status": "completed",
                "repeat_right_status": "completed",
                "direct_repeat_candidate_count": 0,
                "direct_repeat_display_id": ".",
                "tsd_assessment_status": "not_assessed_no_anchor",
                "final_confidence": 0.9,
            },
            {
                "eve_id": "EVE_incomplete",
                "repeat_assessed_start": 1200,
                "repeat_assessed_end": 2500,
                "repeat_left_status": "incomplete",
                "repeat_right_status": "completed",
                "direct_repeat_candidate_count": 2,
                "direct_repeat_display_id": "pair-display",
                "tsd_assessment_status": "incomplete",
                "final_confidence": 0.4,
            },
            {
                "eve_id": "EVE_unassessed",
                "repeat_assessed_start": ".",
                "repeat_assessed_end": ".",
                "repeat_left_status": "not_assessed",
                "repeat_right_status": "not_assessed",
                "direct_repeat_candidate_count": 0,
                "direct_repeat_display_id": ".",
                "tsd_assessment_status": "not_assessed",
                "final_confidence": 0.3,
            },
        ]
    )
    before = predictions.copy(deep=True)
    displayed: list[pd.DataFrame] = []

    exec(_notebook_cell("_repeat_columns = ["), {"predictions": predictions, "display": displayed.append})

    assert len(displayed) == 1
    assert displayed[0]["eve_id"].tolist() == ["EVE_negative", "EVE_incomplete", "EVE_unassessed"]
    assert displayed[0]["repeat_assessed_start"].tolist() == [0, 1200, "."]
    assert displayed[0]["direct_repeat_candidate_count"].tolist() == [0, 2, 0]
    assert displayed[0]["repeat_left_status"].tolist() == ["completed", "incomplete", "not_assessed"]
    assert displayed[0]["tsd_assessment_status"].tolist() == ["not_assessed_no_anchor", "incomplete", "not_assessed"]
    pd.testing.assert_frame_equal(predictions, before)


@pytest.mark.parametrize("rows", [[], [{"eve_id": "EVE_legacy", "start": 0, "end": 1000}]])
def test_notebook_repeat_evidence_handles_empty_and_legacy_outputs(
    rows: list[dict[str, object]], capsys: pytest.CaptureFixture[str]
) -> None:
    """Older reports without repeat fields and empty canonical calls still render."""
    displayed: list[pd.DataFrame] = []

    exec(_notebook_cell("_repeat_columns = ["), {"predictions": pd.DataFrame(rows), "display": displayed.append})

    assert displayed == []
    assert capsys.readouterr().out == "No repeat assessment records in this result.\n"


def _standalone_helper_cell() -> str:
    """Return the notebook's self-contained output readers and plot helpers."""
    return _notebook_cell("def _report_profile_hallmarks(")


def _parameter_defaults() -> dict[str, object]:
    jupytext = pytest.importorskip("jupytext")
    notebook = jupytext.read(report_generate._SOURCE)
    cells = [cell.source for cell in notebook.cells if "parameters" in cell.get("metadata", {}).get("tags", [])]
    assert len(cells) == 1
    namespace: dict[str, object] = {}
    exec(cells[0], namespace)
    return namespace


def test_notebook_no_longer_runs_skani_itself() -> None:
    """Clustering moved into Phase 3; the report must only read its results."""
    code = _notebook_code()

    # Only the invariant: one clustering, computed by Phase 3 and read here. A
    # second skani run in the notebook would let its clusters disagree with the
    # ones the pipeline published in the detailed TSV.
    assert "skani" not in code
    assert "eve_ani_edges.tsv" in code


def test_notebook_ani_parameters_match_the_pipeline_clustering() -> None:
    defaults = _parameter_defaults()

    assert defaults["ANI_THRESHOLD"] == MIN_CLUSTER_ANI
    assert defaults["MIN_ALIGNED_FRACTION"] == MIN_CLUSTER_ALIGNED_FRACTION


def test_notebook_reads_published_taxonomy_class_not_likely_family() -> None:
    """canon_family owns the only remaining read of the retired label."""
    jupytext = pytest.importorskip("jupytext")
    notebook = jupytext.read(report_generate._SOURCE)
    cells = [cell.source for cell in notebook.cells if cell.cell_type == "code" and "likely_family" in cell.source]

    assert len(cells) == 1, "likely_family is read outside canon_family"
    assert "def canon_family(" in cells[0]
    assert 'profile.get("taxonomy_class") or profile.get("likely_family")' in cells[0]
    assert "canon_family(" in _notebook_code()


def test_notebook_canon_family_folds_retired_labels_and_falls_back() -> None:
    namespace: dict[str, object] = {}
    exec(_standalone_helper_cell(), namespace)
    canon_family = namespace["canon_family"]

    assert canon_family({"taxonomy_class": "MIXED"}) == "VIRAL_UNKNOWN"
    assert canon_family({"taxonomy_class": "mixed"}) == "VIRAL_UNKNOWN"
    assert canon_family({"taxonomy_class": "PLV"}) == "PPV"
    assert canon_family({"taxonomy_class": "PHAGE"}) == "PHAGE"
    # Result files written before taxonomy_class existed still render.
    assert canon_family({"likely_family": "MIXED"}) == "VIRAL_UNKNOWN"
    assert canon_family({"likely_family": "VP"}) == "PPV"
    assert canon_family({"likely_family": "NCLDV"}) == "NCLDV"
    assert canon_family({"taxonomy_class": "", "likely_family": "NCLDV"}) == "NCLDV"
    # A published class always wins over the legacy one.
    assert canon_family({"taxonomy_class": "PHAGE", "likely_family": "NCLDV"}) == "PHAGE"
    assert canon_family({}) == "UNKNOWN"


def _ani_namespace(tmp_path: Path) -> dict[str, object]:
    """Namespace for the cluster and network cells, with two published clusters."""
    (tmp_path / "phase3_synthesis").mkdir()
    profiles = {
        "EVE_demo|contig_1_10-90": {
            "taxonomy_class": "NCLDV",
            "has_mcp": True,
            "cluster_id": 0,
            "cluster_size": 2,
            "confidence_tier": "HIGH",
        },
        "EVE_demo|contig_2_10-90": {
            "taxonomy_class": "NCLDV",
            "has_mcp": False,
            "cluster_id": 0,
            "cluster_size": 2,
            "confidence_tier": "HIGH",
        },
        "EVE_demo|contig_3_10-90": {
            "taxonomy_class": "PPV",
            "has_mcp": False,
            "cluster_id": -1,
            "cluster_size": 1,
            "confidence_tier": "LOW",
        },
        # Pre-clustering result file: no cluster_id, no taxonomy_class.
        "EVE_demo|contig_4_10-90": {"likely_family": "VP", "confidence_tier": "LOW"},
    }
    defaults = _parameter_defaults()
    return {
        "profiles": profiles,
        "BASE": tmp_path,
        "SYNTHESIS": tmp_path / "phase3_synthesis",
        "EVE_PREFIX": "EVE_demo",
        "GENOME_ID": "demo",
        "Path": Path,
        "plt": plt,
        "defaultdict": defaultdict,
        "FAMILY_COLORS": {"NCLDV": "#009988", "PPV": "#EE7733", "UNKNOWN": "#BBBBBB"},
        "canon_family": _canon_family(),
        "ANI_THRESHOLD": defaults["ANI_THRESHOLD"],
        "MIN_ALIGNED_FRACTION": defaults["MIN_ALIGNED_FRACTION"],
    }


def _canon_family():
    namespace: dict[str, object] = {}
    exec(_standalone_helper_cell(), namespace)
    return namespace["canon_family"]


def _write_edges(tmp_path: Path, rows: list[tuple[str, str, float, float, float]]) -> None:
    (tmp_path / "phase3_synthesis" / "eve_ani_edges.tsv").write_text(
        "eve_a\teve_b\tani\taf_a\taf_b\n"
        + "".join(f"{a}\t{b}\t{ani}\t{af_a}\t{af_b}\n" for a, b, ani, af_a, af_b in rows)
    )


def test_notebook_clusters_come_from_the_pipeline(tmp_path: Path) -> None:
    namespace = _ani_namespace(tmp_path)

    exec(_notebook_cell("eve_cluster_map = {eid:"), namespace)

    assert namespace["eve_cluster_map"] == {
        "EVE_demo|contig_1_10-90": "Cluster_1",
        "EVE_demo|contig_2_10-90": "Cluster_1",
        "EVE_demo|contig_3_10-90": "singleton",
        "EVE_demo|contig_4_10-90": "singleton",
    }
    assert [len(members) for _, members in namespace["multi_clusters"]] == [2]
    assert (tmp_path / "eve_cluster_composition.png").is_file()


def test_notebook_network_filters_edges_at_95_ani_and_50_aligned_fraction(
    tmp_path: Path,
) -> None:
    namespace = _ani_namespace(tmp_path)
    exec(_notebook_cell("eve_cluster_map = {eid:"), namespace)
    _write_edges(
        tmp_path,
        [
            # kept: at both thresholds exactly
            ("EVE_demo|contig_1_10-90", "EVE_demo|contig_2_10-90", 95.0, 50.0, 12.0),
            # dropped: below the ANI threshold
            ("EVE_demo|contig_1_10-90", "EVE_demo|contig_3_10-90", 94.9, 99.0, 99.0),
            # dropped: neither sequence reaches the aligned fraction
            ("EVE_demo|contig_2_10-90", "EVE_demo|contig_3_10-90", 99.9, 49.9, 49.9),
            # dropped: an endpoint outside the displayed tiers
            ("EVE_demo|contig_1_10-90", "EVE_demo|contig_9_10-90", 99.9, 99.0, 99.0),
        ],
    )

    exec(_notebook_cell("import graphviz"), namespace)

    assert [(a, b) for a, b, _ in namespace["network_edges"]] == [
        ("EVE_demo|contig_1_10-90", "EVE_demo|contig_2_10-90")
    ]
    # Only the connected pair is drawn; the two unconnected EVEs are counted out.
    assert namespace["network_eves"] == [
        "EVE_demo|contig_1_10-90",
        "EVE_demo|contig_2_10-90",
    ]
    assert namespace["n_omitted"] == 2
    assert (tmp_path / "eve_ani_network.png").is_file()


def test_notebook_dense_network_renders_with_bounded_dimensions(tmp_path: Path) -> None:
    namespace = _ani_namespace(tmp_path)
    node_ids = [f"EVE_demo|contig_{index}_10-90" for index in range(48)]
    namespace["profiles"] = {
        eve_id: {
            "taxonomy_class": "NCLDV" if index % 2 else "PPV",
            "cluster_id": 0,
            "cluster_size": len(node_ids),
            "confidence_tier": "HIGH",
        }
        for index, eve_id in enumerate(node_ids)
    }
    exec(_notebook_cell("eve_cluster_map = {eid:"), namespace)
    _write_edges(
        tmp_path,
        [
            (node_ids[index], node_ids[(index + offset) % len(node_ids)], 99.0, 90.0, 90.0)
            for index in range(len(node_ids))
            for offset in (1, 5)
        ],
    )

    exec(_notebook_cell("import graphviz"), namespace)

    png = (tmp_path / "eve_ani_network.png").read_bytes()
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    width, height = struct.unpack(">II", png[16:24])
    assert width <= 3_840
    assert height <= 3_840
    assert namespace["graph"].engine == ANI_NETWORK_ENGINE
    assert {name: namespace["graph"].graph_attr[name] for name in ANI_NETWORK_RENDER_ATTR} == ANI_NETWORK_RENDER_ATTR
    assert len(namespace["network_eves"]) == len(node_ids)
    assert len(namespace["network_edges"]) == 96


@pytest.mark.parametrize("rows", [None, []])
def test_notebook_network_skips_cleanly_without_qualifying_edges(
    tmp_path: Path,
    rows: list | None,
) -> None:
    """A missing edge table and a header-only one must not write a figure."""
    namespace = _ani_namespace(tmp_path)
    exec(_notebook_cell("eve_cluster_map = {eid:"), namespace)
    if rows is not None:
        _write_edges(tmp_path, rows)

    exec(_notebook_cell("import graphviz"), namespace)

    assert namespace["network_edges"] == []
    assert not (tmp_path / "eve_ani_network.png").exists()


def test_notebook_path_helpers_run_from_result_dir_without_virosync(
    tmp_path: Path,
) -> None:
    helper_source = _standalone_helper_cell()
    assert "from virosync" not in helper_source
    script = helper_source + "\nassert _report_profile_hallmarks({}) == []\n"
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)

    result = subprocess.run(
        [sys.executable, "-I", "-c", script],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr


def test_notebook_safe_ratio_handles_zero_denominators() -> None:
    namespace: dict[str, object] = {}
    exec(_standalone_helper_cell(), namespace)

    safe_ratio = namespace["safe_ratio"]
    assert safe_ratio(3, 2) == 1.5
    assert safe_ratio(7, 0) == 0.0
    assert safe_ratio(7, None) == 0.0


def test_notebook_reads_accepted_phylogenetic_rejection_without_assessment() -> None:
    from types import SimpleNamespace

    from virosync.pipeline.phase1.seed_merger import ConfirmedSeed
    from virosync.pipeline.phase2.boundary_refiner import RefinedBoundary
    from virosync.pipeline.phase3.eve_ani_clustering import unsupported_eve_ids
    from virosync.pipeline.phase3.evidence_synthesizer import EvidenceSynthesizer, EvidenceSynthesizerConfig
    from virosync.pipeline.phase3.output_generator import evaluate_v2_quality_gate

    class RejectedPhylogeny(EvidenceSynthesizer):
        def _run_phylogenetic_validation(self, result, boundary):
            result.phylogenetic_result = SimpleNamespace(has_viral_in_any_neighbor=False, rejects_viral=True)
            result.phylogenetic_score = 0.0

    hits = [
        {
            "porf_id": "p1|aa1-20",
            "hallmark_gene": "GVOGm0054",
            "hmm_score": 70,
            "validation_status": "validated",
            "top10_prefixes": "NCLDV__",
            "top10_pidents": "40",
            "top10_targets": "NCLDV__reference|p1",
        },
        {
            "porf_id": "p2|aa1-20",
            "hallmark_gene": "GVOGm0760",
            "hmm_score": 70,
            "validation_status": "validated",
            "top10_prefixes": "NCLDV__",
            "top10_pidents": "40",
            "top10_targets": "NCLDV__reference|p2",
        },
    ]
    boundary = RefinedBoundary(
        "s",
        0,
        6000,
        predicted_family="NCLDV",
        confirmed_seeds=(ConfirmedSeed("seed", "s", 0, 6000),),
    )
    result = RejectedPhylogeny(config=EvidenceSynthesizerConfig()).verify_eve(boundary, [], hallmark_hits=hits)
    assert result.confidence_tier == "LOW"
    assert evaluate_v2_quality_gate(result).kept
    assert unsupported_eve_ids([result]) == set()
    profile = result.to_dict()
    assert profile["hallmark_assessment"] is None

    namespace = {}
    exec(_standalone_helper_cell(), namespace)
    assert namespace["_report_profile_hallmarks"](profile) == []
    assert profile["hallmark_assessment"] is None
    assert profile["taxonomy_class"] == "NCLDV"


def test_notebook_displays_contained_profile_rescue_without_revalidating_it(tmp_path: Path) -> None:
    synthesis = tmp_path / "phase3_synthesis"
    synthesis.mkdir()
    phase0 = tmp_path / "phase0"
    phase0.mkdir()
    phase0.joinpath("proteome.fasta").write_text(
        ">s_1 # 101 # 200 # 1 # ID=s_1;partial=00\nMAAA\n>s_2 # 301 # 400 # 1 # ID=s_2;partial=00\nMAAA\n"
    )
    role = {
        "panel": "Mirusviricota",
        "marker": "MCP",
        "module": "capsid_packaging",
        "role": "core_marker",
        "capsid": True,
        "structural_anchor": True,
        "excluded_from_scoring": False,
    }
    mcp = {
        "evidence_id": "mcp",
        "gene": "s_1",
        "scaffold": "s",
        "start": 100,
        "end": 200,
        "model": "Mirus_MCP",
        "evalue": 0.001,
        "validation_status": "unvalidated",
        "model_roles": [role],
    }
    portal = {
        **mcp,
        "evidence_id": "portal",
        "gene": "s_2",
        "start": 300,
        "end": 400,
        "model": "Mirus_Portal",
        "model_roles": [{**role, "marker": "Portal", "capsid": False}],
    }
    profile = {
        "scaffold": "s",
        "start": 0,
        "end": 1000,
        "status": "accepted",
        "final_confidence": 0.7,
        "confidence_tier": "HIGH",
        "taxonomy_class": "UNKNOWN",
        "hallmark_evidence": [
            mcp,
            portal,
            {**mcp, "evidence_id": "duplicate", "evalue": 0.002},
            {**mcp, "evidence_id": "outside", "gene": "s_3", "end": 1001},
            {**mcp, "evidence_id": "unqualified", "gene": "s_4", "evalue": 2},
        ],
        "hallmark_assessment": {
            "summary": {"qualifying_evidence_ids": ["mcp", "duplicate", "outside"]},
            "witnesses": [{"evidence_ids": ["mcp", "portal"]}],
        },
    }
    synthesis.joinpath("evidence_profiles.json").write_text(json.dumps({"EVE_s_0-1000": profile}))
    namespace = {
        "BASE": tmp_path,
        "SYNTHESIS": synthesis,
        "GENOME_ID": "demo",
        "SHOW_TIERS": ["HIGH"],
        "json": json,
        "pd": pd,
        "defaultdict": defaultdict,
        "re": re,
        "csv": csv,
    }
    exec(_standalone_helper_cell(), namespace)
    exec(_notebook_cell("# ---- Load data ----"), namespace)
    assert namespace["mcp_porf_ids"] == {"s_1"}
    assert namespace["capsid_porf_ids"] == {"s_1"}
    retained = namespace["profile_hallmarks"]["EVE_s_0-1000"]
    assert {(h["gene"], h["model"]) for h in retained} == {("s_1", "Mirus_MCP"), ("s_2", "Mirus_Portal")}
    assert len(retained) == 2
    assert all(h["validation_status"] == "unvalidated" for h in retained)
    assert namespace["profiles"]["EVE_s_0-1000"] == profile
    labels = namespace["gene_function_by_gene"]
    assert labels["EVE_s_0-1000", "s_1"]["label"] == "MCP"
    assert "unvalidated" in labels["EVE_s_0-1000", "s_1"]["evidence"]
    assert "0.001" in labels["EVE_s_0-1000", "s_1"]["evidence"]


def test_notebook_separates_mcp_support_fold_and_gene_taxonomy(
    monkeypatch,
    tmp_path: Path,
) -> None:
    synthesis = tmp_path / "phase3_synthesis"
    synthesis.mkdir()
    (synthesis / "evidence_profiles.json").write_text(
        json.dumps(
            {
                "EVE_contig_0-150": {
                    "scaffold": "contig",
                    "start": 0,
                    "end": 150,
                    "status": "accepted",
                    "final_confidence": 0.9,
                    "confidence_tier": "HIGH",
                    "taxonomy_class": "NCLDV",
                    "hallmark_evidence": [
                        {
                            "evidence_id": "mcp",
                            "gene": "contig_4",
                            "scaffold": "contig",
                            "start": 90,
                            "end": 120,
                            "model": "GVOGm0003",
                            "evalue": 1e-8,
                            "validation_status": "validated_novel",
                            "model_roles": [
                                {
                                    "panel": "NCLDV",
                                    "marker": "MCP",
                                    "module": "capsid_packaging",
                                    "role": "core_marker",
                                    "capsid": True,
                                }
                            ],
                        }
                    ],
                    "hallmark_assessment": {"summary": {"qualifying_evidence_ids": ["mcp"]}, "witnesses": []},
                }
            }
        )
    )
    phase0 = tmp_path / "phase0"
    phase0.mkdir()
    phase0.joinpath("proteome.fasta").write_text(
        "".join(
            f">contig_{index} # {start} # {start + 29} # 1 # ID=contig_{index};partial=00\nM\n"
            for index, start in enumerate(range(1, 151, 30), start=1)
        )
    )
    synthesis.joinpath("virosync_jelly_roll_proteins.tsv").write_text(
        "protein_id\ttype\tconfidence\tlength\tmarker\tevidence\t"
        "mcp_support\tvalidation_status\n"
        "contig_1|aa1-200\tHK97\t0.900\t200\tMirus_MCP\tmarker_family\t"
        "structure_supported\tunvalidated\n"
        "contig_2\tUNKNOWN\t0.000\t300\tMCP\tinsufficient_evidence\t"
        "sequence_supported\tvalidated\n"
        "contig_3\tDJR\t0.500\t410\tMCP\tlength\tcandidate\tunvalidated\n"
        "contig_5\tSJR\t0.700\t300\tMCP\ttmvec:reference\n"
    )
    helper_namespace: dict[str, object] = {"pd": pd}
    exec(_standalone_helper_cell(), helper_namespace)
    namespace: dict[str, object] = {
        "BASE": tmp_path,
        "SYNTHESIS": synthesis,
        "GENOME_ID": "demo",
        "SHOW_TIERS": ["HIGH"],
        "json": json,
        "pd": pd,
        "defaultdict": defaultdict,
        "re": re,
        "csv": csv,
        "_report_profile_hallmarks": helper_namespace["_report_profile_hallmarks"],
        "_report_prodigal_gene": helper_namespace["_report_prodigal_gene"],
        "_function_label_rows": helper_namespace["_function_label_rows"],
    }
    exec(_notebook_cell("# ---- Load data ----"), namespace)

    assert namespace["mcp_porf_ids"] == {"contig_1", "contig_2", "contig_4"}
    classifications = namespace["mcp_classification_map"]
    assert classifications["contig_1"] == {
        "type": "HK97",
        "confidence": 0.9,
        "mcp_support": "structure_supported",
        "validation_status": "unvalidated",
    }
    assert classifications["contig_5"]["mcp_support"] == "candidate"

    namespace.update(
        {
            "BOUNDARY_TAXONOMY": {
                "contig_1": {"top1_prefix": "BAC", "top1_target": "bacterium"},
                "contig_2": {"top1_prefix": "ARC", "top1_target": "archaeon"},
                "contig_3": {"top1_prefix": "BAC", "top1_target": "bacterium"},
                "contig_4": {"top1_prefix": "EUK", "top1_target": "other-eukaryote"},
                "contig_5": {"top1_prefix": "ARC", "top1_target": "archaeon"},
            },
            "_report_viral_category": lambda _record: None,
            "_canonical_viral_category": lambda _origin, _target: None,
            "is_host_gene": lambda _target: False,
            "FLANK_SIZE": 0,
            "canon_family": lambda profile: profile["taxonomy_class"],
        }
    )
    exec(_notebook_cell("def _make_gene_record(gene, region):"), namespace)

    genes = {gene["porf_id"]: gene for gene in namespace["eve_genes_extended"]["EVE_contig_0-150"]["eve_genes"]}
    assert (genes["contig_1"]["is_mcp"], genes["contig_1"]["fold_type"]) == (
        True,
        "HK97",
    )
    assert (genes["contig_2"]["is_mcp"], genes["contig_2"]["fold_type"]) == (
        True,
        "UNKNOWN",
    )
    assert (genes["contig_3"]["is_mcp"], genes["contig_3"]["fold_type"]) == (
        False,
        "DJR",
    )
    assert (genes["contig_5"]["is_mcp"], genes["contig_5"]["fold_type"]) == (
        False,
        "SJR",
    )
    assert [genes[f"contig_{index}"]["category"] for index in range(1, 6)] == [
        "BAC",
        "ARC",
        "BAC",
        "EUK",
        "ARC",
    ]

    # A gene the Phase 2 search never covered has no taxonomy record at all,
    # which the map must not show as a search that returned nothing.
    namespace["BOUNDARY_TAXONOMY"] = {
        key: value for key, value in namespace["BOUNDARY_TAXONOMY"].items() if key != "contig_4"
    }
    exec(_notebook_cell("# ---- Build extended gene map data ----"), namespace)
    unsearched = {
        gene["porf_id"]: gene["category"]
        for eve in namespace["eve_genes_extended"].values()
        for gene in eve["left_flank_genes"] + eve["eve_genes"] + eve["right_flank_genes"]
    }
    assert unsearched["contig_4"] == "UNSEARCHED"
    assert unsearched["contig_1"] == "BAC"

    namespace.update(
        {
            "plt": plt,
            "mpatches": mpatches,
            "Line2D": Line2D,
            "Counter": Counter,
            "COLORS": {
                category: "#BBBBBB"
                for category in (
                    "PPV",
                    "MIRUS",
                    "NCLDV",
                    "GVMAG",
                    "PHAGE",
                    "EUK-pEVE",
                    "HOST",
                    "EUK",
                    "BAC",
                    "ARC",
                    "UNKNOWN",
                    "UNSEARCHED",
                )
            },
            "UNSEARCHED_EDGE": "#777777",
            "UNSEARCHED_LINESTYLE": (0, (2.0, 1.2)),
            "HOST_GENUS": "",
        }
    )
    monkeypatch.setattr(plt, "savefig", lambda *args, **kwargs: None)
    monkeypatch.setattr(plt, "show", lambda: None)
    monkeypatch.setattr(plt, "tight_layout", lambda *args, **kwargs: None)
    exec(
        _notebook_cell("# ---- Gene Map: all canonical EVEs with colored flanking genes ----"),
        namespace,
    )

    marker_lines = [line for line in namespace["ax"].lines if line.get_marker() in {"o", "*"}]
    circles = {float(line.get_xdata()[0]) for line in marker_lines if line.get_marker() == "o"}
    stars = {
        float(line.get_xdata()[0]): line.get_markeredgecolor() for line in marker_lines if line.get_marker() == "*"
    }
    assert circles == {15.0, 45.0, 105.0}
    assert stars == {15.0: "#CC3311", 75.0: "#000000", 135.0: "#0077BB"}
    legend_labels = [text.get_text() for text in namespace["ax"].get_legend().get_texts()]
    assert "Supported MCP" in legend_labels
    assert "HK97 fold classification" in legend_labels
    assert "DJR fold classification" in legend_labels
    assert "SJR fold classification" in legend_labels
    assert not any("structural" in label.lower() for label in legend_labels)
    plt.close(namespace["fig"])


def test_notebook_taxonomy_resolver_matches_pipeline() -> None:
    from virosync.pipeline.taxonomy_utils import resolve_org_id

    taxonomy_lookup = {
        "EUK__EP00224": "Eukaryota|Chlorophyta",
        "PHAGE__GCA-000906975-1": ("Viruses|Varidnaviria|Bamfordvirae|Preplasmiviricota"),
    }
    targets = [
        "EUK__EP00224|protein",
        "EUK__EP00224_Organism_Name|protein",
        "PHAGE__VARDNA__GCA-000906975-1_23|protein",
        "UNKNOWN__record_1|protein",
        "plain-target|protein",
        "",
    ]
    namespace: dict[str, object] = {"pd": pd}
    exec(_standalone_helper_cell(), namespace)
    report_resolve = namespace["_report_resolve_taxonomy_target"]

    assert [report_resolve(target, taxonomy_lookup) for target in targets] == [
        resolve_org_id(target, taxonomy_lookup) for target in targets
    ]


def test_notebook_taxonomy_csv_and_display_helpers_cover_supported_values() -> None:
    from virosync.pipeline.phase2.boundary_diamond import MIN_VIRAL_HIT_PIDENT

    taxonomy_lookup = {
        "EUK__host": "Eukaryota|Chlorophyta",
        "PHAGE__plv": "Viruses|Varidnaviria|Preplasmiviricota",
        "GVMAG__example": "Viruses|Varidnaviria|Nucleocytoviricota",
    }
    namespace: dict[str, object] = {"pd": pd}
    exec(_standalone_helper_cell(), namespace)

    csv_values = namespace["_csv_values"]
    canonical = namespace["_report_canonical_viral_category"]
    report_viral = namespace["_report_identity_qualified_viral_category"]
    assert namespace["_REPORT_MIN_VIRAL_HIT_PIDENT"] == MIN_VIRAL_HIT_PIDENT
    assert csv_values(None) == []
    assert csv_values(float("nan")) == []
    assert csv_values(("NCLDV__", "GVMAG__")) == ["NCLDV__", "GVMAG__"]
    assert csv_values("NCLDV__,GVMAG__") == ["NCLDV__", "GVMAG__"]
    assert canonical("PHAGE", "PHAGE__plv|protein", taxonomy_lookup) == "PPV"
    assert canonical("GVMAG", "GVMAG__example|protein", taxonomy_lookup) == "GVMAG"
    assert (
        report_viral(
            {
                "top10_prefixes": "EUK__,GVMAG__",
                "top10_targets": "EUK__host|protein,GVMAG__example|protein",
                "top10_pidents": "99.0,31.0",
            },
            taxonomy_lookup,
        )
        == "GVMAG"
    )
    assert (
        report_viral(
            {
                "top10_prefixes": "GVMAG__,EUK__,NCLDV__",
                "top10_targets": "GVMAG__example|protein",
                "top10_pidents": "not-a-number,99.0",
            },
            taxonomy_lookup,
        )
        is None
    )
    assert (
        report_viral(
            {
                "top10_prefixes": "GVMAG__",
                "top10_targets": "GVMAG__example|protein",
                "top10_pidents": str(MIN_VIRAL_HIT_PIDENT - 0.1),
            },
            taxonomy_lookup,
        )
        is None
    )


PEVE_TARGET = "EUK-pEVE__EP00022|XP_1"
NCLDV_TARGET = "NCLDV__GCA-000001234-1|YP_2"


def test_notebook_euk_peve_label_matches_the_pipeline_in_any_case() -> None:
    from virosync.pipeline.taxonomy_utils import EUK_PEVE_LABEL

    namespace: dict[str, object] = {"pd": pd}
    exec(_standalone_helper_cell(), namespace)
    canonical = namespace["_report_canonical_viral_category"]

    assert namespace["_REPORT_EUK_PEVE"] == EUK_PEVE_LABEL
    # The taxonomy map stores the prefix with its separator; upper case is the token other code sees.
    assert canonical("EUK-pEVE__", PEVE_TARGET, {}) == "EUK-pEVE"
    assert canonical("EUK-PEVE", PEVE_TARGET, {}) == "EUK-pEVE"
    assert canonical("EUK__", "EUK__EP00224|protein", {}) is None


@pytest.mark.parametrize(
    ("record", "expected"),
    [
        ({"top10_prefixes": "EUK-pEVE__", "top10_targets": PEVE_TARGET, "top10_pidents": "60.0"}, "EUK-pEVE"),
        (
            {
                "top10_prefixes": "EUK-pEVE__,NCLDV__",
                "top10_targets": f"{PEVE_TARGET},{NCLDV_TARGET}",
                "top10_pidents": "90.0,30.0",
            },
            "NCLDV",
        ),
        (
            {
                "top10_prefixes": "NCLDV__,EUK-pEVE__",
                "top10_targets": f"{NCLDV_TARGET},{PEVE_TARGET}",
                "top10_pidents": "20.0,30.0",
            },
            "EUK-pEVE",
        ),
        ({"top10_prefixes": "EUK-pEVE__", "top10_targets": PEVE_TARGET, "top10_pidents": "24.9"}, None),
    ],
    ids=["only-euk-peve", "lineage-hit-wins", "lineage-hit-below-cutoff", "euk-peve-below-cutoff"],
)
def test_notebook_viral_category_takes_euk_peve_only_without_a_qualifying_lineage_hit(
    record: dict[str, str], expected: str | None
) -> None:
    namespace: dict[str, object] = {"pd": pd}
    exec(_standalone_helper_cell(), namespace)

    assert namespace["_report_identity_qualified_viral_category"](record, {}) == expected


CAPSID = (
    "MSNKVQLTPEFAAAIRKGLLDNGYTVEEIAKALGVSRQTVYNWLKGEGSLPDADLLLKLASIFGVSVDYLLGLESEEVRSHAQSLAEKLNRLSPEQRNAVIQLVEALAQTSVSH"
)
PRIMASE = (
    "MKKLLPTAAAGLLLLAAQPAMAQVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVS"
    "AISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK"
)
FLANK = "MDEHWQTRPLAVFGGSIKNCYEQPMLRDTSAHGVIFNKWEYCPLQDRTSMAGHVEKNFWLYTRCPDSIQAMNE"
SYNTVIZ_CELL = "# ---- Gene neighborhoods in Syntviz ----"


def _syntviz_gene(
    porf_id: str,
    start: int,
    end: int,
    category: str,
    strand: str = "+",
    *,
    is_mcp: bool = False,
    is_capsid: bool = False,
    fold_type: str = "",
    function_label: str = "",
    function_evidence: str = "",
):
    return {
        "porf_id": porf_id,
        "start": start,
        "end": end,
        "strand": strand,
        "category": category,
        "is_mcp": is_mcp,
        "is_capsid": is_capsid,
        "fold_type": fold_type,
        "function_label": function_label,
        "function_evidence": function_evidence,
    }


def _syntviz_eve(contig: str, tier: str, confidence: float, core: list[dict]) -> dict:
    return {
        "scaffold": contig,
        "eve_start": 1_000,
        "eve_end": 2_000,
        "ext_start": 0,
        "ext_end": 3_000,
        "length": 3_000,
        "eve_length": 1_000,
        "confidence": confidence,
        "confidence_tier": tier,
        "left_flank_genes": [_syntviz_gene(f"{contig}_1", 100, 400, "HOST")],
        "eve_genes": core,
        "right_flank_genes": [],
    }


def _syntviz_namespace(tmp_path: Path, eves: dict[str, dict], proteins: dict[str, str]) -> dict[str, object]:
    (tmp_path / "phase0").mkdir(exist_ok=True)
    flanks = {gene["porf_id"]: FLANK for info in eves.values() for gene in info["left_flank_genes"]}
    (tmp_path / "phase0" / "proteome.fasta").write_text(
        "".join(f">{porf_id} # 1 # 2 # 1 # ID=x\n{seq}\n" for porf_id, seq in (flanks | proteins).items())
    )
    displayed: list[object] = []
    return {
        "Path": Path,
        "display": displayed.append,
        "displayed": displayed,
        "SYNTVIZ_INTERACTIVE": False,
        "FLANK_SIZE": 5_000,
        "BASE": tmp_path,
        "GENOME_ID": "demo",
        "SHOW_TIERS": ["HIGH", "MEDIUM", "LOW"],
        "COLORS": {
            "PPV": "#EE7733",
            "NCLDV": "#009988",
            "EUK-pEVE": "#882255",
            "HOST": "#0077BB",
            "UNKNOWN": "#BBBBBB",
        },
        "eve_genes_extended": eves,
    }


def _three_eves(tmp_path: Path) -> dict[str, object]:
    """Two HIGH EVEs share a capsid protein; one MEDIUM EVE shares nothing."""
    eves = {
        "EVE_demo|contig_1_1000-2000": _syntviz_eve(
            "contig_1",
            "HIGH",
            0.9,
            [_syntviz_gene("contig_1_2", 1_100, 1_463, "PPV", is_mcp=True, fold_type="DJR")],
        ),
        "EVE_demo|contig_2_1000-2000": _syntviz_eve(
            "contig_2",
            "HIGH",
            0.8,
            [_syntviz_gene("contig_2_2", 1_200, 1_563, "PPV", is_mcp=True)],
        ),
        "EVE_demo|contig_3_1000-2000": _syntviz_eve(
            "contig_3",
            "MEDIUM",
            0.5,
            [_syntviz_gene("contig_3_2", 1_100, 1_463, "NCLDV")],
        ),
    }
    return _syntviz_namespace(
        tmp_path,
        eves,
        {"contig_1_2": CAPSID, "contig_2_2": CAPSID, "contig_3_2": PRIMASE},
    )


def _syntviz_session(tmp_path: Path):
    from syntviz.models import Session

    return Session.model_validate_json((tmp_path / "eve_gene_map.syntviz.json").read_text())


def _displayed_types(namespace: dict[str, object]) -> list[str]:
    return [type(item).__name__ for item in namespace["displayed"]]


def test_notebook_syntviz_figure_links_shared_proteins_and_groups_tiers(
    tmp_path: Path,
) -> None:
    namespace = _three_eves(tmp_path)

    exec(_notebook_cell(SYNTVIZ_CELL), namespace)

    session = _syntviz_session(tmp_path)
    rows = {row.id: row for row in session.sequences}
    assert [rows[row_id].name.split()[0] for row_id in session.view.row_order] == ["H", "H", "M"]
    assert rows["EVE_demo|contig_3_1000-2000:window"].name == "M 0.50 c3 1–2 kb"
    assert {row.length_bp for row in rows.values()} == {3_000}
    assert len(session.features) == 6
    assert {annotation.value for annotation in session.annotations} == {
        "PPV",
        "NCLDV",
        "HOST",
    }
    groups = {group.id: set(group.feature_ids) for group in session.view.style_groups}
    capsids = {
        "EVE_demo|contig_1_1000-2000:contig_1_2",
        "EVE_demo|contig_2_1000-2000:contig_2_2",
    }
    assert groups["mcp"] == capsids
    assert groups["flanks"] == {f"EVE_demo|contig_{n}_1000-2000:contig_{n}_1" for n in (1, 2, 3)}
    assert len(session.homology_hits) == 1
    hit = session.homology_hits[0]
    assert {hit.query_feature_id, hit.target_feature_id} == capsids
    assert hit.identity_fraction == 1.0
    assert len(session.runs[0].searched_pairs) == 4
    assert session.view.link_filters.hit_color == "identity"
    # The source range stays in the tooltip of each row.
    assert session.view.row_coordinates == "none"
    gene_markers = {
        (marker.anchor.feature_id, marker.glyph)
        for marker in session.view.markers
        if marker.anchor.kind == "feature_segment"
    }
    assert gene_markers == {
        ("EVE_demo|contig_1_1000-2000:contig_1_2", "circle"),
        ("EVE_demo|contig_1_1000-2000:contig_1_2", "star"),
        ("EVE_demo|contig_2_1000-2000:contig_2_2", "circle"),
    }
    boundaries = {
        (marker.anchor.sequence_id, marker.anchor.position_bp, marker.glyph)
        for marker in session.view.markers
        if marker.anchor.kind == "sequence_position"
    }
    assert boundaries == {(row_id, bp, "diamond") for row_id in rows for bp in (1_000, 1_999)}
    assert (tmp_path / "eve_gene_map_syntviz.svg").read_text().startswith("<svg")
    assert (tmp_path / "eve_gene_map_syntviz.png").read_bytes().startswith(b"\x89PNG")
    # The shipped report embeds the static PNG, not the multi-megabyte viewer.
    assert _displayed_types(namespace) == ["Image"]


def test_notebook_syntviz_figure_orders_similar_rows_and_flips_reverse_capsids(
    tmp_path: Path,
) -> None:
    """Confidence order A, B, C, D; A and C share a capsid, B and D a primase."""
    eves = {
        f"EVE_demo|contig_{name}_1000-2000": _syntviz_eve(
            f"contig_{name}",
            "HIGH",
            confidence,
            [
                _syntviz_gene(
                    f"contig_{name}_2",
                    1_100,
                    1_463,
                    "PPV",
                    strand,
                    is_mcp=protein == CAPSID,
                )
            ],
        )
        for name, confidence, protein, strand in [
            ("a", 0.9, CAPSID, "+"),
            ("b", 0.8, PRIMASE, "+"),
            ("c", 0.7, CAPSID, "-"),
            ("d", 0.6, PRIMASE, "+"),
        ]
    }
    proteins = {
        "contig_a_2": CAPSID,
        "contig_b_2": PRIMASE,
        "contig_c_2": CAPSID,
        "contig_d_2": PRIMASE,
    }
    namespace = _syntviz_namespace(tmp_path, eves, proteins)

    exec(_notebook_cell(SYNTVIZ_CELL), namespace)

    session = _syntviz_session(tmp_path)
    order = [row_id.split("_")[2] for row_id in session.view.row_order]
    assert abs(order.index("a") - order.index("c")) == 1
    assert abs(order.index("b") - order.index("d")) == 1
    assert session.view.orientations["EVE_demo|contig_c_1000-2000:window"] == -1
    assert session.view.orientations["EVE_demo|contig_a_1000-2000:window"] == 1


@pytest.mark.parametrize("reverse_hits", [False, True])
def test_notebook_syntviz_figure_bounds_hits_ordering_and_final_coverage(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    reverse_hits: bool,
) -> None:
    from syntviz import runtime
    from syntviz.models import Session

    eves = {
        f"{tier}_{index:03}": _syntviz_eve(
            f"{tier}_{index:03}",
            tier,
            1 - index / 1_000,
            [
                _syntviz_gene(
                    f"{tier}_{index:03}_2",
                    1_100,
                    1_463,
                    "PPV",
                    "-" if index == 0 else "+",
                    is_mcp=True,
                    function_label="PolB",
                    function_evidence="Pfam:DNA_pol_B",
                )
            ],
        )
        for tier, count in (("HIGH", 501), ("MEDIUM", 500), ("LOW", 1))
        for index in range(count)
    }
    proteins = {f"{eve_id}_2": CAPSID for eve_id in eves if eve_id != "HIGH_500"}
    namespace = _syntviz_namespace(tmp_path, eves, proteins)
    namespace["SYNTVIZ_INTERACTIVE"] = True
    # Score, E-value, stable target ID and alignment coordinates each break a tie.
    # Rows 1/0 and 3/2 also test reciprocal collapse with equal and unequal scores.
    hit_rows = [
        "0\t1\t80\t5\t104\t5\t104\t1e-30\t100",
        "0\t1\t80\t1\t99\t1\t99\t1e-20\t150",
        "0\t1\t80\t3\t102\t3\t102\t1e-25\t150",
        "0\t1\t80\t1\t100\t1\t100\t1e-25\t150",
        "0\t2\t80\t1\t100\t1\t100\t1e-25\t150",
        "1\t0\t80\t1\t100\t1\t100\t1e-25\t150",
        "2\t4\t70\t1\t100\t1\t100\t1e-20\t90",
        "2\t3\t90\t1\t100\t1\t100\t1e-20\t100",
        "3\t2\t90\t1\t100\t1\t100\t1e-20\t120",
        "4\t0\t19\t1\t100\t1\t100\t1e-20\t200",
        "0\t0\t100\t1\t100\t1\t100\t1e-30\t200",
    ]
    hit_text = "\n".join(reversed(hit_rows) if reverse_hits else hit_rows) + "\n"
    diamond = tmp_path / "diamond"
    diamond.write_text(
        f"#!{sys.executable}\n"
        "import pathlib\n"
        "import sys\n"
        "assert sys.argv[sys.argv.index('--threads') + 1] == '1'\n"
        "if sys.argv[1] == 'blastp':\n"
        f"    pathlib.Path(sys.argv[sys.argv.index('-o') + 1]).write_text({hit_text!r})\n"
    )
    diamond.chmod(0o755)
    node = tmp_path / "node"
    node.touch()
    node.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    high_rows = [f"HIGH_{index:03}:window" for index in range(501)]
    medium_rows = [f"MEDIUM_{index:03}:window" for index in range(500)]
    reordered_medium = medium_rows[-1:] + medium_rows[:-1]
    ordered_inputs = []
    rendered = []

    def order_rows(session: Session, operations: list[dict]) -> Session:
        assert operations == [{"op": "orderRowsBySimilarity", "sequence_ids": medium_rows}]
        ordered_inputs.append(session.model_dump())
        updated = session.model_dump()
        updated["view"]["row_order"] = high_rows + reordered_medium + ["LOW_000:window"]
        return Session.model_validate(updated)

    def render_figure(session: Session, output: Path) -> None:
        # Attribute access requires validated GenomePair models at both renders.
        assert session.runs[0].searched_pairs[0].query_genome_id == "HIGH_000"
        rendered.append((output.suffix, session.model_dump()))
        output.write_text("render spy")

    monkeypatch.setattr(runtime, "execute", order_rows)
    monkeypatch.setattr(runtime, "render", render_figure)

    exec(_notebook_cell(SYNTVIZ_CELL), namespace)

    session = _syntviz_session(tmp_path)
    assert len(ordered_inputs) == 1
    assert session.view.row_order == high_rows + reordered_medium + ["LOW_000:window"]
    assert {row.id for row in session.sequences} == {f"{eve_id}:window" for eve_id in eves}
    assert {feature.id for feature in session.features} == {
        f"{eve_id}:{eve_id}_{suffix}" for eve_id in eves for suffix in (1, 2)
    }
    assert {annotation.id for annotation in session.annotations} == {
        f"{eve_id}:{eve_id}_{suffix}:{kind}"
        for eve_id in eves
        for suffix, kind in ((1, "category"), (2, "category"), (2, "function"))
    }
    assert session.view.orientations == {f"{eve_id}:window": -1 if eve_id.endswith("_000") else 1 for eve_id in eves}
    saved = session.model_dump()
    assert {key: saved[key] for key in ("genomes", "sequences", "features", "annotations")} == {
        key: ordered_inputs[0][key] for key in ("genomes", "sequences", "features", "annotations")
    }
    assert rendered == [(".svg", saved), (".png", saved)]
    assert session.runs[0].parameters == {
        "link_evalue": 1e-5,
        "link_min_pident": 20.0,
        "search_ran": True,
        "searched_eve_count": 1_000,
        "core_protein_count": 1_000,
        "qualifying_rows": 9,
        "retained_hits": 2,
        "retention_policy": "best cross-EVE hit per query protein, then reciprocal collapse",
    }
    assert [
        (hit.query_feature_id, hit.target_feature_id, hit.score, hit.evalue, hit.query_interval.start_aa)
        for hit in session.homology_hits
    ] == [
        ("HIGH_000:HIGH_000_2", "HIGH_001:HIGH_001_2", 150.0, 1e-25, 0),
        ("HIGH_003:HIGH_003_2", "HIGH_002:HIGH_002_2", 120.0, 1e-20, 0),
    ]
    pairs = {(pair.query_genome_id, pair.target_genome_id) for pair in session.runs[0].searched_pairs}
    # Include searched neighbours with zero hits, but neither LOW nor absent-protein rows.
    assert ("MEDIUM_499", "MEDIUM_000") in pairs
    assert ("MEDIUM_499", "MEDIUM_498") not in pairs
    assert ("HIGH_499", "MEDIUM_499") not in pairs
    assert not any("HIGH_500" in pair or "LOW_000" in pair for pair in pairs)
    expected_pairs = {
        (f"{tier}_{index:03}", f"{tier}_{index + 1:03}")
        for tier, count in (("HIGH", 499), ("MEDIUM", 498))
        for index in range(count)
    } | {("MEDIUM_499", "MEDIUM_000")}
    assert pairs == expected_pairs | {(target, query) for query, target in expected_pairs}
    assert len(session.runs[0].searched_pairs) == 1_996 <= 2 * (len(eves) - 1)


def test_notebook_syntviz_figure_clips_genes_at_the_window_edges(
    tmp_path: Path,
) -> None:
    eve = _syntviz_eve("contig_1", "HIGH", 0.9, [_syntviz_gene("contig_1_2", 1_100, 3_300, "PPV")])
    eve["left_flank_genes"] = [_syntviz_gene("contig_1_1", -60, 400, "HOST")]
    namespace = _syntviz_namespace(tmp_path, {"EVE_demo|contig_1_1000-2000": eve}, {"contig_1_2": CAPSID})

    exec(_notebook_cell(SYNTVIZ_CELL), namespace)

    features = {feature.id.split(":")[1]: feature for feature in _syntviz_session(tmp_path).features}
    left, core = features["contig_1_1"], features["contig_1_2"]
    assert (left.segments[0].start_bp, left.partial_start) == (0, True)
    assert (core.segments[0].end_bp, core.partial_end) == (3_000, True)
    assert not core.partial_start


def test_notebook_syntviz_figure_links_no_low_confidence_eves(tmp_path: Path) -> None:
    eves = {
        f"EVE_demo|contig_{n}_1000-2000": _syntviz_eve(
            f"contig_{n}", "LOW", 0.1, [_syntviz_gene(f"contig_{n}_2", 1_100, 1_463, "PPV")]
        )
        for n in (1, 2)
    }
    namespace = _syntviz_namespace(tmp_path, eves, {"contig_1_2": CAPSID, "contig_2_2": CAPSID})

    exec(_notebook_cell(SYNTVIZ_CELL), namespace)

    session = _syntviz_session(tmp_path)
    assert session.homology_hits == []
    assert session.runs[0].searched_pairs == []
    assert session.runs[0].parameters["search_ran"] is False
    assert session.runs[0].parameters["searched_eve_count"] == 0
    assert session.runs[0].parameters["core_protein_count"] == 0


def test_notebook_syntviz_figure_keeps_proteins_that_a_low_eve_shares(tmp_path: Path) -> None:
    """A LOW EVE overlapping a HIGH one shares its capsid gene."""
    eves = {
        f"EVE_demo|contig_{n}_1000-2000": _syntviz_eve(
            f"contig_{n}", "HIGH", 0.9, [_syntviz_gene(f"contig_{n}_2", 1_100, 1_463, "PPV")]
        )
        for n in (1, 2)
    }
    eves["EVE_demo|contig_2_1050-2050"] = _syntviz_eve(
        "contig_2", "LOW", 0.1, [_syntviz_gene("contig_2_2", 1_100, 1_463, "PPV")]
    )
    namespace = _syntviz_namespace(tmp_path, eves, {"contig_1_2": CAPSID, "contig_2_2": CAPSID})

    exec(_notebook_cell(SYNTVIZ_CELL), namespace)

    session = _syntviz_session(tmp_path)
    (hit,) = session.homology_hits
    assert {hit.query_feature_id, hit.target_feature_id} == {
        "EVE_demo|contig_1_1000-2000:contig_1_2",
        "EVE_demo|contig_2_1000-2000:contig_2_2",
    }
    assert len(session.runs[0].searched_pairs) == 2


def test_notebook_syntviz_figure_records_searches_without_hits(tmp_path: Path) -> None:
    eves = {
        "EVE_demo|contig_1_1000-2000": _syntviz_eve("contig_1", "HIGH", 0.9, [_syntviz_gene("contig_1_2", 1_100, 1_463, "PPV")]),
        "EVE_demo|contig_2_1000-2000": _syntviz_eve("contig_2", "HIGH", 0.8, [_syntviz_gene("contig_2_2", 1_100, 1_463, "PPV")]),
    }  # fmt: skip
    namespace = _syntviz_namespace(tmp_path, eves, {"contig_1_2": CAPSID, "contig_2_2": PRIMASE})

    exec(_notebook_cell(SYNTVIZ_CELL), namespace)

    session = _syntviz_session(tmp_path)
    assert session.homology_hits == []
    assert len(session.runs[0].searched_pairs) == 2
    assert session.runs[0].parameters["search_ran"] is True
    assert session.runs[0].parameters["qualifying_rows"] == 0
    assert session.runs[0].parameters["retained_hits"] == 0


def test_notebook_syntviz_figure_fits_many_rows_in_the_renderer_limit(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    eves = {
        f"EVE_demo|contig_{index}_1000-2000": _syntviz_eve(
            f"contig_{index}",
            "LOW",
            0.1,
            [_syntviz_gene(f"contig_{index}_2", 1_100, 1_463, "PPV", is_mcp=True)],
        )
        for index in range(350)
    }
    namespace = _syntviz_namespace(tmp_path, eves, {})

    exec(_notebook_cell(SYNTVIZ_CELL), namespace)

    session = _syntviz_session(tmp_path)
    assert len(session.sequences) == 350
    assert session.view.height_px <= 16_384
    assert (tmp_path / "eve_gene_map_syntviz.png").read_bytes().startswith(b"\x89PNG")
    # Syntviz draws at most 1,000 markers; each EVE keeps both diamonds and its circle or none.
    markers_per_row = Counter(marker.id.split(":")[1] for marker in session.view.markers)
    assert set(markers_per_row.values()) == {3}
    assert len(markers_per_row) == 333
    assert "markers omitted for 17 EVEs" in capsys.readouterr().out


def test_notebook_syntviz_figure_is_interactive_on_request(tmp_path: Path) -> None:
    namespace = _three_eves(tmp_path)
    namespace["SYNTVIZ_INTERACTIVE"] = True

    exec(_notebook_cell(SYNTVIZ_CELL), namespace)

    assert _displayed_types(namespace) == ["Figure"]


def test_notebook_syntviz_figure_needs_neither_diamond_nor_node(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    namespace = _three_eves(tmp_path)
    monkeypatch.setenv("PATH", str(tmp_path))

    exec(_notebook_cell(SYNTVIZ_CELL), namespace)

    session = _syntviz_session(tmp_path)
    assert len(session.sequences) == 3
    assert session.homology_hits == []
    assert session.runs[0].searched_pairs == []
    assert not (tmp_path / "eve_gene_map_syntviz.svg").exists()
    assert session.runs[0].parameters["search_ran"] is False
    assert session.runs[0].parameters["searched_eve_count"] == 0
    assert session.runs[0].parameters["core_protein_count"] == 3
    assert _displayed_types(namespace) == ["Figure"]


def test_notebook_syntviz_figure_names_a_gene_by_its_function_and_annotates_the_evidence(tmp_path: Path) -> None:
    eves = {
        "EVE_demo|contig_1_1000-2000": _syntviz_eve(
            "contig_1",
            "HIGH",
            0.9,
            [
                _syntviz_gene(
                    "contig_1_2", 1_100, 1_463, "PPV", function_label="PolB", function_evidence="Pfam:DNA_pol_B"
                ),
                _syntviz_gene("contig_1_3", 1_500, 1_800, "PPV", is_mcp=True),
                _syntviz_gene("contig_1_4", 1_850, 1_990, "PPV"),
            ],
        )
    }
    namespace = _syntviz_namespace(tmp_path, eves, {"contig_1_2": CAPSID, "contig_1_3": CAPSID, "contig_1_4": PRIMASE})

    exec(_notebook_cell(SYNTVIZ_CELL), namespace)

    session = _syntviz_session(tmp_path)
    names = {feature.id.split(":")[1]: feature.name for feature in session.features}
    assert names["contig_1_2"] == "PolB"
    assert names["contig_1_3"] == "MCP"
    assert names["contig_1_4"] == "contig_1_4"
    function = [annotation for annotation in session.annotations if annotation.namespace == "function"]
    assert [(a.target.feature_id, a.value) for a in function] == [("EVE_demo|contig_1_1000-2000:contig_1_2", "PolB")]
    assert "Pfam:DNA_pol_B" in function[0].evidence_source


def test_notebook_syntviz_figure_keeps_the_capsid_outline_when_a_capsid_has_a_function_label(tmp_path: Path) -> None:
    eves = {
        "EVE_demo|contig_1_1000-2000": _syntviz_eve(
            "contig_1",
            "HIGH",
            0.9,
            [
                _syntviz_gene(
                    "contig_1_2",
                    1_100,
                    1_463,
                    "PPV",
                    is_mcp=True,
                    function_label="MCP",
                    function_evidence="GVOGm0003",
                ),
                _syntviz_gene(
                    "contig_1_3", 1_500, 1_800, "PPV", function_label="MCP", function_evidence="Pfam:Capsid_NCLDV"
                ),
            ],
        )
    }
    namespace = _syntviz_namespace(tmp_path, eves, {"contig_1_2": CAPSID, "contig_1_3": CAPSID})

    exec(_notebook_cell(SYNTVIZ_CELL), namespace)

    groups = {group.id: set(group.feature_ids) for group in _syntviz_session(tmp_path).view.style_groups}
    assert groups["mcp"] == {"EVE_demo|contig_1_1000-2000:contig_1_2"}


def test_notebook_syntviz_figure_outlines_capsid_evidence_with_a_style_group_not_markers(tmp_path: Path) -> None:
    """An MCP keeps its group and circle; other capsid evidence, core or flank, gets the same outline."""
    eve = _syntviz_eve(
        "contig_1",
        "HIGH",
        0.9,
        [
            _syntviz_gene("contig_1_2", 1_100, 1_463, "EUK-pEVE", is_mcp=True, is_capsid=True),
            _syntviz_gene("contig_1_3", 1_500, 1_800, "EUK-pEVE", is_capsid=True),
            _syntviz_gene("contig_1_4", 1_850, 1_990, "PPV"),
        ],
    )
    eve["right_flank_genes"] = [_syntviz_gene("contig_1_5", 2_100, 2_400, "HOST", is_capsid=True)]
    proteins = {"contig_1_2": CAPSID, "contig_1_3": CAPSID, "contig_1_4": PRIMASE, "contig_1_5": FLANK}
    namespace = _syntviz_namespace(tmp_path, {"EVE_demo|contig_1_1000-2000": eve}, proteins)

    exec(_notebook_cell(SYNTVIZ_CELL), namespace)

    view = _syntviz_session(tmp_path).view
    groups = {group.id: group for group in view.style_groups}
    assert set(groups["mcp"].feature_ids) == {"EVE_demo|contig_1_1000-2000:contig_1_2"}
    assert set(groups["capsid"].feature_ids) == {
        "EVE_demo|contig_1_1000-2000:contig_1_3",
        "EVE_demo|contig_1_1000-2000:contig_1_5",
    }
    assert groups["capsid"].style == groups["mcp"].style
    assert groups["capsid"].style.stroke == "#000000"
    gene_markers = {
        (marker.anchor.feature_id, marker.glyph) for marker in view.markers if marker.anchor.kind == "feature_segment"
    }
    assert gene_markers == {("EVE_demo|contig_1_1000-2000:contig_1_2", "circle")}
    assert view.colors.palette["EUK-pEVE"] == "#882255"


def test_notebook_syntviz_palette_colors_categories_only(tmp_path: Path) -> None:
    eves = {
        "EVE_demo|contig_1_1000-2000": _syntviz_eve(
            "contig_1",
            "HIGH",
            0.9,
            [_syntviz_gene("contig_1_2", 1_100, 1_463, "PPV", function_label="PolB", function_evidence="x")],
        )
    }
    namespace = _syntviz_namespace(tmp_path, eves, {"contig_1_2": CAPSID})

    exec(_notebook_cell(SYNTVIZ_CELL), namespace)

    palette = _syntviz_session(tmp_path).view.colors.palette
    assert set(palette) == {"PPV", "HOST"}


@pytest.mark.parametrize(
    ("spans", "expected"),
    [
        ([(0, 10), (12, 20), (30, 40)], [0, 0, 0]),
        ([(0, 10), (5, 15), (30, 40)], [0, 1, 0]),
        ([(0, 10), (10.5, 20)], [0, 1]),
        ([(0, 10), (5, 15), (8, 18)], [0, 1, None]),
        ([], []),
    ],
)
def test_notebook_function_labels_move_to_a_free_row_or_are_dropped(
    spans: list[tuple[float, float]], expected: list[int | None]
) -> None:
    namespace: dict[str, object] = {"pd": pd}
    exec(_standalone_helper_cell(), namespace)

    assert namespace["_function_label_rows"](spans, 1.0, 2) == expected


GENE_FUNCTIONS_HEADER = "eve_id\tporf_id\tlabel\tmodule\tevidence\tevalue\n"


def _load_data_namespace(tmp_path: Path, gene_functions: str | None) -> dict[str, object]:
    """Run the notebook's load cell on one EVE whose five genes are 30 bp apart."""
    synthesis = tmp_path / "phase3_synthesis"
    synthesis.mkdir()
    (synthesis / "evidence_profiles.json").write_text(
        json.dumps(
            {
                "EVE_contig_30-120": {
                    "scaffold": "contig",
                    "start": 30,
                    "end": 120,
                    "status": "accepted",
                    "final_confidence": 0.9,
                    "confidence_tier": "HIGH",
                    "taxonomy_class": "NCLDV",
                }
            }
        )
    )
    (tmp_path / "phase0").mkdir()
    (tmp_path / "phase0" / "proteome.fasta").write_text(
        "".join(
            f">contig_{index} # {start} # {start + 29} # 1 # ID=contig_{index};partial=00\nM\n"
            for index, start in enumerate(range(1, 151, 30), start=1)
        )
    )
    if gene_functions is not None:
        (synthesis / "gene_functions.tsv").write_text(GENE_FUNCTIONS_HEADER + gene_functions)
    helpers: dict[str, object] = {"pd": pd}
    exec(_standalone_helper_cell(), helpers)
    namespace: dict[str, object] = {
        "BASE": tmp_path,
        "SYNTHESIS": synthesis,
        "GENOME_ID": "demo",
        "SHOW_TIERS": ["HIGH"],
        "json": json,
        "pd": pd,
        "defaultdict": defaultdict,
        "re": re,
        "csv": csv,
        "_report_profile_hallmarks": helpers["_report_profile_hallmarks"],
        "_report_prodigal_gene": helpers["_report_prodigal_gene"],
        "_function_label_rows": helpers["_function_label_rows"],
    }
    exec(_notebook_cell("# ---- Load data ----"), namespace)
    return namespace


def _gene_record_namespace(tmp_path: Path, gene_functions: str | None) -> dict[str, object]:
    namespace = _load_data_namespace(tmp_path, gene_functions)
    namespace.update(
        {
            "BOUNDARY_TAXONOMY": {},
            "_report_viral_category": lambda _record: None,
            "_canonical_viral_category": lambda _origin, _target: None,
            "is_host_gene": lambda _target: False,
            "FLANK_SIZE": 0,
            "canon_family": lambda profile: profile["taxonomy_class"],
        }
    )
    exec(_notebook_cell("def _make_gene_record(gene, region):"), namespace)
    return namespace


def test_notebook_loads_gene_functions_by_eve_and_gene(tmp_path: Path) -> None:
    namespace = _load_data_namespace(tmp_path, "EVE_contig_30-120\tcontig_2\tPolB\treplication\tGVOGm0054\t1e-20\n")

    assert namespace["gene_function_by_gene"] == {
        ("EVE_contig_30-120", "contig_2"): {"label": "PolB", "module": "replication", "evidence": "GVOGm0054"}
    }


def test_notebook_loads_no_gene_functions_when_the_table_is_absent(tmp_path: Path) -> None:
    assert _load_data_namespace(tmp_path, None)["gene_function_by_gene"] == {}


def test_notebook_gives_core_genes_their_function_and_flanking_genes_none(tmp_path: Path) -> None:
    namespace = _gene_record_namespace(
        tmp_path,
        "EVE_contig_30-120\tcontig_2\tPolB\treplication\tGVOGm0054\t1e-20\n"
        "EVE_contig_30-120\tcontig_1\tMCP\tcapsid\tGVOGm0003\t1e-20\n",
    )

    eve = namespace["eve_genes_extended"]["EVE_contig_30-120"]

    assert [(g["porf_id"], g["function_label"], g["function_evidence"]) for g in eve["eve_genes"]] == [
        ("contig_2", "PolB", "GVOGm0054"),
        ("contig_3", "", ""),
        ("contig_4", "", ""),
    ]
    assert all(not gene.get("function_label") for gene in eve["left_flank_genes"] + eve["right_flank_genes"])


def _draw_gene_map(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, gene_functions: str | None) -> dict[str, object]:
    namespace = _gene_record_namespace(tmp_path, gene_functions)
    namespace.update(
        {
            "plt": plt,
            "mpatches": mpatches,
            "Line2D": Line2D,
            "Counter": Counter,
            "COLORS": defaultdict(lambda: "#BBBBBB"),
            "UNSEARCHED_EDGE": "#777777",
            "UNSEARCHED_LINESTYLE": (0, (2.0, 1.2)),
            "HOST_GENUS": "",
        }
    )
    monkeypatch.setattr(plt, "savefig", lambda *args, **kwargs: None)
    monkeypatch.setattr(plt, "show", lambda: None)
    monkeypatch.setattr(plt, "tight_layout", lambda *args, **kwargs: None)
    exec(_notebook_cell("# ---- Gene Map: all canonical EVEs with colored flanking genes ----"), namespace)
    return namespace


def test_notebook_gene_map_prints_the_function_under_each_labelled_core_gene(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    namespace = _draw_gene_map(
        tmp_path,
        monkeypatch,
        "EVE_contig_30-120\tcontig_2\tPolB\treplication\tGVOGm0054\t1e-20\n"
        "EVE_contig_30-120\tcontig_4\tPortal\tcapsid\tMirus_Portal\t1e-30\n",
    )

    ax = namespace["ax"]
    labels = {text.get_text(): text for text in ax.texts if text.get_text() in {"PolB", "Portal"}}
    assert set(labels) == {"PolB", "Portal"}
    assert labels["PolB"].get_position()[0] == pytest.approx(45.0)
    assert labels["Portal"].get_position()[0] == pytest.approx(105.0)
    assert all(text.get_position()[1] < 0 for text in labels.values())
    plt.close(namespace["fig"])


def test_notebook_gene_map_labels_never_overlap(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    long_label = "Terminase" * 7
    crowded = "".join(f"EVE_contig_30-120\tcontig_{index}\t{long_label}\tcapsid\tx\t1e-20\n" for index in (2, 3, 4))
    namespace = _draw_gene_map(tmp_path, monkeypatch, crowded)

    ax = namespace["ax"]
    renderer = namespace["fig"].canvas.get_renderer()
    boxes = [text.get_window_extent(renderer) for text in ax.texts if text.get_text() == long_label]
    assert len(boxes) == 3
    assert not any(a.overlaps(b) for i, a in enumerate(boxes) for b in boxes[i + 1 :])
    assert len({round(text.get_position()[1], 3) for text in ax.texts if text.get_text() == long_label}) == 2
    plt.close(namespace["fig"])


def test_notebook_gene_map_draws_no_function_labels_without_the_table(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    namespace = _draw_gene_map(tmp_path, monkeypatch, None)

    assert not [text for text in namespace["ax"].texts if text.get_text() in {"PolB", "Portal", "MCP"}]
    plt.close(namespace["fig"])


def test_notebook_gene_map_ignores_a_rescued_protein_that_no_drawn_gene_carries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    namespace = _draw_gene_map(
        tmp_path,
        monkeypatch,
        "EVE_contig_30-120\tcontig_2\tPolB\treplication\tGVOGm0054\t1e-20\n"
        "EVE_contig_30-120\tcontig_VSR0123456789abcdef\tTerminase\tcapsid\tPfam:Terminase_6N\t1e-30\n",
    )

    eve = namespace["eve_genes_extended"]["EVE_contig_30-120"]
    assert "contig_VSR0123456789abcdef" not in {gene["porf_id"] for gene in eve["eve_genes"]}
    assert [gene["function_label"] for gene in eve["eve_genes"] if gene["function_label"]] == ["PolB"]
    assert {text.get_text() for text in namespace["ax"].texts} & {"PolB", "Terminase"} == {"PolB"}
    plt.close(namespace["fig"])


GENE_MAP_CELL = "# ---- Gene Map: all canonical EVEs with colored flanking genes ----"
DETAIL_CELL = "DETAIL_MAX_EVES = 6"
DISTRIBUTION_CELL = "# ---- Stacked bar: core vs flank gene composition ----"
NOTEBOOK_CELLS_TO_GENE_RECORDS = (
    "from IPython import get_ipython",
    "def _report_profile_hallmarks(",
    "# --- Style setup (beautiful-data-viz) ---",
    "# ---- Load data ----",
    "# ---- Build extended gene map data ----",
)
HOST_LINEAGE = "EUK|Amorphea|Amoebozoa|Evosea|Eumycetozoa|Dictyostelium|Dictyostelium discoideum"
HOST_TOKENS = {"amorphea", "amoebozoa", "evosea", "eumycetozoa", "dictyostelium", "dictyostelium discoideum"}
PEVE_LINEAGE = (
    "EUK-pEVE|Amorphea-pEVE|Amoebozoa-pEVE|Evosea-pEVE|Eumycetozoa-pEVE|Cavenderia-pEVE|Cavenderia fasciculata-pEVE"
)
PEVE_TOKENS = {
    "euk-peve",
    "amorphea-peve",
    "amoebozoa-peve",
    "evosea-peve",
    "eumycetozoa-peve",
    "cavenderia-peve",
    "cavenderia fasciculata-peve",
}
FIVE_GENES_30_BP_APART = "".join(
    f">contig_{index} # {start} # {start + 29} # 1 # ID=contig_{index};partial=00\nM\n"
    for index, start in enumerate(range(1, 151, 30), start=1)
)


@pytest.fixture
def euk_peve_report(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Run the notebook through the gene records of one EVE with host, EUK-pEVE and capsid genes.

    contig_1 is the left flank and contig_5 the right flank; contig_2 to contig_4 form the core.
    The host model weighs every -pEVE rank above the host's own ranks.
    """
    labels = tmp_path / "labels.tsv"
    labels.write_text(
        f"EUK__EP00001\t{HOST_LINEAGE}\n"
        "EUK__EP00099\tEUK|Sar|Alveolata|Ciliophora|Oligohymenophorea|Tetrahymena|Tetrahymena thermophila\n"
        f"EUK-pEVE__EP00022\t{PEVE_LINEAGE}\n"
        "NCLDV__GCA-000001234-1\tNCLDV|Nucleocytoviricota|Megaviricetes|Imitervirales|Mimiviridae\n"
    )
    marker_dir = tmp_path / "phase1" / "marker_validation"
    marker_dir.mkdir(parents=True)
    token_weights = dict.fromkeys(HOST_TOKENS, 50.0) | dict.fromkeys(PEVE_TOKENS, 90.0) | {"sar": 10.0}
    (marker_dir / "host_signature_model.json").write_text(json.dumps({"token_weights": token_weights}))
    boundary_dir = tmp_path / "phase2" / "boundary_diamond"
    boundary_dir.mkdir(parents=True)
    (boundary_dir / "taxonomy_map.tsv").write_text(
        "porf_id\ttop1_target\ttop1_prefix\ttop10_prefixes\ttop10_targets\ttop10_pidents\n"
        "contig_1\tEUK__EP00001|XP_10\tEUK__\tEUK__\tEUK__EP00001|XP_10\t80.00\n"
        f"contig_2\t{PEVE_TARGET}\tEUK-pEVE__\tEUK-pEVE__\t{PEVE_TARGET}\t60.00\n"
        f"contig_3\t{PEVE_TARGET}\tEUK-pEVE__\tEUK-pEVE__,NCLDV__\t{PEVE_TARGET},{NCLDV_TARGET}\t90.00,30.00\n"
        f"contig_4\t{PEVE_TARGET}\tEUK-pEVE__\tEUK-pEVE__\t{PEVE_TARGET}\t20.00\n"
        "contig_5\tEUK__EP00099|XP_5\tEUK__\tEUK__\tEUK__EP00099|XP_5\t70.00\n"
    )
    synthesis = tmp_path / "phase3_synthesis"
    synthesis.mkdir()
    (synthesis / "evidence_profiles.json").write_text(
        json.dumps(
            {
                "EVE_contig_30-120": {
                    "scaffold": "contig",
                    "start": 30,
                    "end": 120,
                    "status": "accepted",
                    "final_confidence": 0.9,
                    "confidence_tier": "HIGH",
                    "taxonomy_class": "VIRAL_UNKNOWN",
                    "mcp_gene_ids": ["contig_4"],
                    "capsid_gene_ids": ["contig_2", "contig_4"],
                }
            }
        )
    )
    (synthesis / "virosync_jelly_roll_proteins.tsv").write_text(
        "protein_id\ttype\tconfidence\tlength\tmarker\tevidence\tmcp_support\tvalidation_status\n"
        "contig_5\tDJR\t0.900\t300\tMCP\tlength\tsequence_supported\tvalidated\n"
    )
    (tmp_path / "phase0").mkdir()
    (tmp_path / "phase0" / "proteome.fasta").write_text(FIVE_GENES_30_BP_APART)
    monkeypatch.setattr(plt, "savefig", lambda *args, **kwargs: None)
    monkeypatch.setattr(plt, "show", lambda: None)
    monkeypatch.setattr(plt, "tight_layout", lambda *args, **kwargs: None)
    namespace = _parameter_defaults() | {
        "RESULTS_DIR": str(tmp_path),
        "TAX_LABELS_PATH": str(labels),
        "SHOW_TIERS": ["HIGH"],
        "FLANK_SIZE": 0,
    }
    # The style cell sets global rcParams; keep them from leaking into other tests.
    with plt.rc_context():
        for needle in NOTEBOOK_CELLS_TO_GENE_RECORDS:
            exec(_notebook_cell(needle), namespace)
        yield namespace
    plt.close("all")


def _report_genes(namespace: dict[str, object]) -> dict[str, dict]:
    eve = namespace["eve_genes_extended"]["EVE_contig_30-120"]
    regions = ("left_flank_genes", "eve_genes", "right_flank_genes")
    return {gene["porf_id"]: gene for region in regions for gene in eve[region]}


def test_notebook_counts_euk_peve_genes_as_viral_and_never_as_host(euk_peve_report: dict[str, object]) -> None:
    genes = _report_genes(euk_peve_report)

    assert {porf_id: gene["category"] for porf_id, gene in genes.items()} == {
        "contig_1": "HOST",
        "contig_2": "EUK-pEVE",
        "contig_3": "NCLDV",
        "contig_4": "EUK-pEVE",
        "contig_5": "EUK",
    }
    assert {porf_id: gene["host_lineage_match"] for porf_id, gene in genes.items()} == {
        "contig_1": True,
        "contig_2": False,
        "contig_3": False,
        "contig_4": False,
        "contig_5": False,
    }


def test_notebook_host_lineage_drops_euk_peve_ranks(euk_peve_report: dict[str, object]) -> None:
    is_host_gene = euk_peve_report["is_host_gene"]

    assert euk_peve_report["HOST_TAXA"] == HOST_TOKENS
    assert euk_peve_report["HOST_GENUS"] == "dictyostelium"
    assert is_host_gene(PEVE_TARGET) is False
    assert is_host_gene("EUK__EP00001|XP_10") is True


def test_notebook_gives_euk_peve_a_color_no_other_category_uses(euk_peve_report: dict[str, object]) -> None:
    colors = euk_peve_report["COLORS"]

    assert list(colors.values()).count(colors["EUK-pEVE"]) == 1


def test_notebook_marks_capsid_evidence_apart_from_mcp(euk_peve_report: dict[str, object]) -> None:
    genes = _report_genes(euk_peve_report)

    assert euk_peve_report["capsid_porf_ids"] == {"contig_2", "contig_4", "contig_5"}
    assert euk_peve_report["mcp_porf_ids"] == {"contig_4", "contig_5"}
    assert {porf_id: (gene["is_capsid"], gene["is_mcp"]) for porf_id, gene in genes.items()} == {
        "contig_1": (False, False),
        "contig_2": (True, False),
        "contig_3": (False, False),
        "contig_4": (True, True),
        "contig_5": (True, True),
    }


def test_notebook_gene_map_outlines_capsid_evidence_and_keeps_the_mcp_circle(
    euk_peve_report: dict[str, object],
) -> None:
    colors = euk_peve_report["COLORS"]

    exec(_notebook_cell(GENE_MAP_CELL), euk_peve_report)

    ax = euk_peve_report["ax"]
    arrows = {
        round(min(x for x, _ in patch.get_xy())): patch for patch in ax.patches if isinstance(patch, mpatches.Polygon)
    }
    outlined = {start for start, arrow in arrows.items() if arrow.get_linewidth() > 0}
    assert outlined == {30, 90, 120}
    assert {mcolors.to_hex(arrows[start].get_edgecolor()) for start in outlined} == {"#000000"}
    assert mcolors.to_hex(arrows[30].get_facecolor()) == colors["EUK-pEVE"]
    assert {float(line.get_xdata()[0]) for line in ax.lines if line.get_marker() == "o"} == {105.0, 135.0}
    legend = ax.get_legend()
    handles = {text.get_text(): handle for text, handle in zip(legend.get_texts(), legend.legend_handles, strict=True)}
    assert mcolors.to_hex(handles["EUK-pEVE (viral-like, eukaryotic genome)"].get_facecolor()) == colors["EUK-pEVE"]
    assert mcolors.to_hex(handles["Thick black outline: capsid evidence"].get_edgecolor()) == "#000000"
    assert "Supported MCP" in handles


def test_notebook_detail_map_outlines_capsid_evidence_and_lists_euk_peve(
    euk_peve_report: dict[str, object],
) -> None:
    exec(_notebook_cell(DETAIL_CELL), euk_peve_report)

    fig = euk_peve_report["fig"]
    arrows = {
        round(min(x for x, _ in patch.get_xy()), 2): patch
        for patch in fig.axes[0].patches
        if isinstance(patch, mpatches.FancyArrow)
    }
    black = {start for start, arrow in arrows.items() if mcolors.to_hex(arrow.get_edgecolor()) == "#000000"}
    assert black == {0.0, 0.06, 0.09}
    assert min(arrows[start].get_linewidth() for start in black) > max(
        arrows[start].get_linewidth() for start in arrows.keys() - black
    )
    labels = [text.get_text() for text in fig.legends[0].get_texts()]
    assert "EUK-pEVE (viral-like, eukaryotic genome)" in labels
    assert "Thick black outline: capsid evidence" in labels
    assert "Supported MCP" in labels


def test_notebook_category_distribution_counts_euk_peve_and_sums_to_one(
    euk_peve_report: dict[str, object],
) -> None:
    exec(_notebook_cell(GENE_MAP_CELL), euk_peve_report)
    exec(_notebook_cell(DISTRIBUTION_CELL), euk_peve_report)

    containers = {container.get_label(): container for container in euk_peve_report["ax"].containers}
    heights = {label: [bar.get_height() for bar in container] for label, container in containers.items()}
    assert heights["EUK-pEVE"] == pytest.approx([0.0, 2 / 3, 0.0])
    assert [sum(region) for region in zip(*heights.values(), strict=True)] == pytest.approx([1.0, 1.0, 1.0])
    assert mcolors.to_hex(containers["EUK-pEVE"][1].get_facecolor()) == euk_peve_report["COLORS"]["EUK-pEVE"]
