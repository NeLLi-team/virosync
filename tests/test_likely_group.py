from __future__ import annotations

from virosync.pipeline.phase3.evidence_synthesizer import (
    infer_likely_group,
    load_marker_annotation_index,
)


def test_annotation_index_parses_capscan_group(tmp_path) -> None:
    ann = tmp_path / "ann.tsv"
    ann.write_text(
        "model_name\tsource\tdescription\tmajority_annotation\n"
        "plv_mcp_caps_PgVV_Aquinto\tcapscan_Bellas2026\tcapscan PgVV Major Capsid Protein\tGroup II dsDNA virus coat\n"
        "plv_mcp_caps_Trimcap_B1_Aquinto\tcapscan_Bellas2026\tcapscan Trimcap_cluster_2 Major Capsid Protein\tx\n"
        "PLV_MCP_3\tVP_PLV_update_Dec25\tPLV Major Capsid Protein\ty\n"
    )
    idx = load_marker_annotation_index(ann)
    assert idx["plv_mcp_caps_pgvv_aquinto"]["capscan_group"] == "PgVV"
    assert idx["plv_mcp_caps_trimcap_b1_aquinto"]["capscan_group"] == "Trimcap_cluster_2"
    # non-capscan markers carry no group
    assert idx["plv_mcp_3"]["capscan_group"] == ""


def test_infer_likely_group_picks_best_hit() -> None:
    idx = {
        "plv_mcp_caps_pgvv_aquinto": {"capscan_group": "PgVV"},
        "plv_mcp_caps_trimcap_b1_aquinto": {"capscan_group": "Trimcap_cluster_2"},
    }
    # highest-scoring capscan hit wins, even when another group has more hits
    hits = [
        ("plv_mcp_caps_Trimcap_B1_Aquinto", 80.0),
        ("plv_mcp_caps_Trimcap_B1_Aquinto", 90.0),
        ("plv_mcp_caps_PgVV_Aquinto", 150.0),
    ]
    assert infer_likely_group(hits, idx) == "PgVV"
    # a higher-scoring non-group marker does not suppress the best group-bearing hit
    assert (
        infer_likely_group([("VP_ATPase_1", 300.0), ("plv_mcp_caps_PgVV_Aquinto", 120.0)], idx)
        == "PgVV"
    )
    # no capscan group-defining marker -> empty
    assert infer_likely_group([("PLV_MCP_3", 200.0), ("VP_ATPase_1", 100.0)], idx) == ""
    assert infer_likely_group([], idx) == ""
    assert infer_likely_group([("plv_mcp_caps_PgVV_Aquinto", 100.0)], None) == ""


def test_infer_likely_group_best_score_decides_between_groups() -> None:
    idx = {
        "plv_mcp_caps_a": {"capscan_group": "GroupA"},
        "plv_mcp_caps_b": {"capscan_group": "GroupB"},
    }
    # the stronger hit decides, even with one marker per group
    assert infer_likely_group([("plv_mcp_caps_A", 120.0), ("plv_mcp_caps_B", 90.0)], idx) == "GroupA"
    assert infer_likely_group([("plv_mcp_caps_B", 110.0), ("plv_mcp_caps_A", 80.0)], idx) == "GroupB"
