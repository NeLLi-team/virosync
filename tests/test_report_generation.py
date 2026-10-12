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
    return cells[0]


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
    """Return the notebook's self-contained helper cell.

    Anchored on the MCP detector because that helper mirrors a pipeline
    predicate and is parity-tested below, so it cannot quietly leave the cell.
    """
    return _notebook_cell("def _report_is_mcp_gene(")


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
    script = helper_source + "\nassert _report_is_mcp_gene('vp_mcp_1')\n"
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


def test_notebook_mcp_helper_matches_canonical_detector() -> None:
    from virosync.pipeline.phase3.mcp_detection import (
        MCP_HMM_EXACT,
        MCP_HMM_PREFIXES,
        is_mcp_gene,
    )

    namespace: dict[str, object] = {}
    exec(_standalone_helper_cell(), namespace)

    assert namespace["_REPORT_MCP_HMM_EXACT"] == MCP_HMM_EXACT
    assert namespace["_REPORT_MCP_HMM_PREFIXES"] == MCP_HMM_PREFIXES
    report_is_mcp_gene = namespace["_report_is_mcp_gene"]
    corpus = [
        "OG1352",
        "gamadvirusMCP",
        "plv_MCP_1",
        "vp_MCP_3",
        "mcp_mirus",
        "dmcp",
        "mcp_lookalike",
        "ncmcp_pseudoprotein",
        "polb",
        "",
        None,
    ]
    assert [report_is_mcp_gene(name) for name in corpus] == [is_mcp_gene(name) for name in corpus]


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
    marker_dir = tmp_path / "phase1" / "marker_validation"
    marker_dir.mkdir(parents=True)
    marker_dir.joinpath("validated_marker_hits.tsv").write_text(
        "query_porf\thmm_target\tvalidation_status\ncontig_4\tMCP\tvalidated_novel\ncontig_3\tMCP\tunvalidated\n"
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
        "_report_is_mcp_gene": helper_namespace["_report_is_mcp_gene"],
        "_report_prodigal_gene": helper_namespace["_report_prodigal_gene"],
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
                    "HOST",
                    "EUK",
                    "BAC",
                    "ARC",
                    "UNKNOWN",
                )
            },
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
    fold_type: str = "",
):
    return {
        "porf_id": porf_id,
        "start": start,
        "end": end,
        "strand": strand,
        "category": category,
        "is_mcp": is_mcp,
        "fold_type": fold_type,
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
    assert len(session.runs[0].searched_pairs) == 6
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
    assert _displayed_types(namespace) == ["Figure"]
