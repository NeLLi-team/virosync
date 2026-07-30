from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from virosync.report import generate as report_generate
from virosync.report.generate import generate_eve_report
from virosync.utils.path_safety import safe_filename_component, safe_filename_components


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


def _standalone_path_safety_cell() -> str:
    jupytext = pytest.importorskip("jupytext")
    notebook = jupytext.read(report_generate._SOURCE)
    cells = [
        cell.source
        for cell in notebook.cells
        if "def safe_filename_component(value):" in cell.source
    ]
    assert len(cells) == 1
    return cells[0]


def test_notebook_path_helpers_run_from_result_dir_without_virosync(
    tmp_path: Path,
) -> None:
    helper_source = _standalone_path_safety_cell()
    assert "from virosync" not in helper_source
    script = helper_source + "\nassert safe_filename_component('../EVE/1').isascii()\n"
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
    exec(_standalone_path_safety_cell(), namespace)

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
    exec(_standalone_path_safety_cell(), namespace)

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
    exec(_standalone_path_safety_cell(), namespace)
    report_resolve = namespace["_report_resolve_taxonomy_target"]

    assert [report_resolve(target, taxonomy_lookup) for target in targets] == [
        resolve_org_id(target, taxonomy_lookup) for target in targets
    ]


def test_notebook_taxonomy_csv_and_display_helpers_cover_supported_values() -> None:
    taxonomy_lookup = {
        "PHAGE__plv": "Viruses|Varidnaviria|Preplasmiviricota",
        "GVMAG__example": "Viruses|Varidnaviria|Nucleocytoviricota",
    }
    namespace: dict[str, object] = {"pd": pd}
    exec(_standalone_path_safety_cell(), namespace)

    csv_values = namespace["_csv_values"]
    canonical = namespace["_report_canonical_viral_category"]
    assert csv_values(None) == []
    assert csv_values(float("nan")) == []
    assert csv_values(("NCLDV__", "GVMAG__")) == ["NCLDV__", "GVMAG__"]
    assert csv_values("NCLDV__,GVMAG__") == ["NCLDV__", "GVMAG__"]
    assert canonical("PHAGE", "PHAGE__plv|protein", taxonomy_lookup) == "PPV"
    assert canonical("GVMAG", "GVMAG__example|protein", taxonomy_lookup) == "GVMAG"


def test_notebook_path_helpers_match_runtime_encoder() -> None:
    namespace: dict[str, object] = {}
    exec(_standalone_path_safety_cell(), namespace)
    raw_ids = [
        "EVE_ok",
        "a/b",
        "a|b",
        "../EVE",
        "EVE λ",
        "EVE\n",
        "x" * 300,
        "λ" * 100,
    ]

    notebook_components = namespace["safe_filename_components"](
        raw_ids,
        label="EVE ID",
    )

    assert notebook_components == safe_filename_components(raw_ids, label="EVE ID")
    assert all(
        namespace["safe_filename_component"](raw_id) == safe_filename_component(raw_id)
        for raw_id in raw_ids
    )
