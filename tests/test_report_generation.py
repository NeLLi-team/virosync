from __future__ import annotations

from collections import defaultdict
import os
import struct
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import matplotlib.pyplot as plt
import pandas as pd
import pytest

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
    tagged = [
        cell
        for cell in notebook.cells
        if "parameters" in cell.get("metadata", {}).get("tags", [])
    ]
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
    return "\n".join(
        cell.source for cell in notebook.cells if cell.cell_type == "code"
    )


def _standalone_helper_cell() -> str:
    """Return the notebook's self-contained helper cell.

    Anchored on the MCP detector because that helper mirrors a pipeline
    predicate and is parity-tested below, so it cannot quietly leave the cell.
    """
    return _notebook_cell("def _report_is_mcp_gene(name):")


def _parameter_defaults() -> dict[str, object]:
    jupytext = pytest.importorskip("jupytext")
    notebook = jupytext.read(report_generate._SOURCE)
    cells = [
        cell.source
        for cell in notebook.cells
        if "parameters" in cell.get("metadata", {}).get("tags", [])
    ]
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
    cells = [
        cell.source
        for cell in notebook.cells
        if cell.cell_type == "code" and "likely_family" in cell.source
    ]

    assert len(cells) == 1, "likely_family is read outside canon_family"
    assert "def canon_family(profile):" in cells[0]
    assert "profile.get('taxonomy_class') or profile.get('likely_family')" in cells[0]
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
            "taxonomy_class": "NCLDV", "has_mcp": True,
            "cluster_id": 0, "cluster_size": 2, "confidence_tier": "HIGH",
        },
        "EVE_demo|contig_2_10-90": {
            "taxonomy_class": "NCLDV", "has_mcp": False,
            "cluster_id": 0, "cluster_size": 2, "confidence_tier": "HIGH",
        },
        "EVE_demo|contig_3_10-90": {
            "taxonomy_class": "PPV", "has_mcp": False,
            "cluster_id": -1, "cluster_size": 1, "confidence_tier": "LOW",
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
    assert {
        name: namespace["graph"].graph_attr[name]
        for name in ANI_NETWORK_RENDER_ATTR
    } == ANI_NETWORK_RENDER_ATTR
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
    assert [report_is_mcp_gene(name) for name in corpus] == [
        is_mcp_gene(name) for name in corpus
    ]


def test_notebook_taxonomy_resolver_matches_pipeline() -> None:
    from virosync.pipeline.taxonomy_utils import resolve_org_id

    taxonomy_lookup = {
        "EUK__EP00224": "Eukaryota|Chlorophyta",
        "PHAGE__GCA-000906975-1": (
            "Viruses|Varidnaviria|Bamfordvirae|Preplasmiviricota"
        ),
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
    assert report_viral(
        {
            "top10_prefixes": "EUK__,GVMAG__",
            "top10_targets": "EUK__host|protein,GVMAG__example|protein",
            "top10_pidents": "99.0,31.0",
        },
        taxonomy_lookup,
    ) == "GVMAG"
    assert report_viral(
        {
            "top10_prefixes": "GVMAG__,EUK__,NCLDV__",
            "top10_targets": "GVMAG__example|protein",
            "top10_pidents": "not-a-number,99.0",
        },
        taxonomy_lookup,
    ) is None
    assert report_viral(
        {
            "top10_prefixes": "GVMAG__",
            "top10_targets": "GVMAG__example|protein",
            "top10_pidents": str(MIN_VIRAL_HIT_PIDENT - 0.1),
        },
        taxonomy_lookup,
    ) is None
