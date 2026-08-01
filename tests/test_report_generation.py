from __future__ import annotations

from collections import defaultdict
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from virosync.report import generate as report_generate
from virosync.report.generate import generate_eve_report
from virosync.utils.path_safety import (
    require_strict_child,
    safe_filename_component,
    safe_filename_components,
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


def _ani_clustering_cell() -> str:
    jupytext = pytest.importorskip("jupytext")
    notebook = jupytext.read(report_generate._SOURCE)
    cells = [
        cell.source
        for cell in notebook.cells
        if "if len(profiles) < 3:" in cell.source
    ]
    assert len(cells) == 1
    return cells[0]


def _ani_namespace(tmp_path: Path) -> dict[str, object]:
    eve_ids = ["eve1", "eve2", "eve3"]
    fasta = tmp_path / "demo_eves.fna"
    fasta.write_text(
        "".join(f">{eve_id} tier=LOW\nACGTACGT\n" for eve_id in eve_ids)
    )
    return {
        "profiles": {eve_id: {} for eve_id in eve_ids},
        "SHOW_TIERS": ("HIGH", "MEDIUM", "LOW"),
        "BASE": tmp_path,
        "tempfile": tempfile,
        "Path": Path,
        "subprocess": subprocess,
        "safe_filename_components": safe_filename_components,
        "require_strict_child": require_strict_child,
        "ANI_THRESHOLD": 90.0,
        "defaultdict": defaultdict,
        "plt": plt,
    }


def test_notebook_ani_clustering_keeps_unsketchable_eves_as_singletons(
    monkeypatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(shutil, "which", lambda name: "/fake/skani")
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(
            returncode=1,
            stderr="ERROR No genomes/sketches found.\n",
        ),
    )
    namespace = _ani_namespace(tmp_path)

    exec(_ani_clustering_cell(), namespace)

    assert namespace["eve_cluster_map"] == {
        "eve1": "singleton",
        "eve2": "singleton",
        "eve3": "singleton",
    }
    assert namespace["multi_clusters"] == []
    assert (tmp_path / "eve_cluster_composition.png").is_file()


def test_notebook_ani_clustering_rejects_missing_filtered_eves(
    monkeypatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(shutil, "which", lambda name: "/fake/skani")
    namespace = _ani_namespace(tmp_path)
    (tmp_path / "demo_eves.fna").write_text(
        ">eve1 tier=REJECTED\nACGTACGT\n"
    )

    with pytest.raises(RuntimeError, match="EVE FASTA/profile mismatch"):
        exec(_ani_clustering_cell(), namespace)


def test_notebook_ani_clustering_keeps_other_skani_errors_fatal(
    monkeypatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(shutil, "which", lambda name: "/fake/skani")
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(
            returncode=1,
            stderr="ERROR corrupted sketch database\n",
        ),
    )

    with pytest.raises(RuntimeError, match="skani triangle failed"):
        exec(_ani_clustering_cell(), _ani_namespace(tmp_path))


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
    from virosync.pipeline.phase2.boundary_diamond import MIN_VIRAL_HIT_PIDENT

    taxonomy_lookup = {
        "EUK__host": "Eukaryota|Chlorophyta",
        "PHAGE__plv": "Viruses|Varidnaviria|Preplasmiviricota",
        "GVMAG__example": "Viruses|Varidnaviria|Nucleocytoviricota",
    }
    namespace: dict[str, object] = {"pd": pd}
    exec(_standalone_path_safety_cell(), namespace)

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
