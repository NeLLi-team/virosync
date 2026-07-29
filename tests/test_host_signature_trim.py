from __future__ import annotations

from types import SimpleNamespace

from virosync.pipeline.host_signatures import HostSignatureModel
from virosync.pipeline.phase1.seed_merger import MergedSeed
from virosync.pipeline.phase2.host_signature_trim import (
    HostTrimParams,
    trim_seed_by_host_signature,
    trim_seeds_by_host_signature,
)


def _host_model() -> HostSignatureModel:
    return HostSignatureModel(
        token_weights={"euk": 1.0},
        token_counts={"euk": 1},
        max_weight=1.0,
        min_token_length=3,
        host_prefixes=["EUK__"],
        weight_mode="rank",
    )


def test_trim_seed_summary_uses_placeholder_consensus() -> None:
    seed = MergedSeed(scaffold="ctg1", start=100, end=500)
    records = [
        {
            "porf_id": "g1",
            "porf_start": 150,
            "porf_end": 260,
            "top10_targets": "",
            "top10_bitscores": "",
            "top10_pidents": "",
            "top10_evalues": "",
            "has_ncldv_mirus": True,
            "has_vp_plv": False,
            "has_viral": True,
        }
    ]

    _trimmed_seed, summary = trim_seed_by_host_signature(
        seed=seed,
        gene_records=records,
        host_model=_host_model(),
        validated_markers=[],
        params=HostTrimParams(window_bp=200, step_bp=200),
    )

    assert summary["host_consensus_taxonomy"] == "."


def test_trim_seeds_summary_uses_placeholder_consensus() -> None:
    seed = MergedSeed(scaffold="ctg2", start=100, end=500)
    eve_id = f"EVE_{seed.scaffold}_{seed.start}-{seed.end}"
    gene_taxonomy_map = {
        eve_id: (
            [
                {
                    "porf_id": "g2",
                    "porf_start": 160,
                    "porf_end": 300,
                    "top10_targets": "",
                    "top10_bitscores": "",
                    "top10_pidents": "",
                    "top10_evalues": "",
                    "has_ncldv_mirus": True,
                    "has_vp_plv": False,
                    "has_viral": True,
                }
            ],
            {},
        )
    }

    _trimmed_seeds, summaries = trim_seeds_by_host_signature(
        seeds=[seed],
        gene_taxonomy_map=gene_taxonomy_map,
        host_model=_host_model(),
        validated_markers=[],
        params=HostTrimParams(window_bp=200, step_bp=200),
    )

    assert summaries[0]["host_consensus_taxonomy"] == "."


def _viral_record(porf_id: str, start: int, end: int) -> dict:
    return {
        "porf_id": porf_id,
        "porf_start": start,
        "porf_end": end,
        "top10_targets": "",
        "top10_bitscores": "",
        "top10_pidents": "",
        "top10_evalues": "",
        "has_ncldv_mirus": True,
        "has_vp_plv": False,
        "has_viral": True,
    }


def test_marker_at_window_right_edge_belongs_only_to_next_window() -> None:
    seed = MergedSeed(scaffold="ctg", start=0, end=20)
    records = [_viral_record("left", 0, 10), _viral_record("right", 10, 20)]
    marker = SimpleNamespace(
        validation_status="validated",
        scaffold="ctg",
        start=10,
        end=11,
    )

    trimmed, summary = trim_seed_by_host_signature(
        seed=seed,
        gene_records=records,
        host_model=_host_model(),
        validated_markers=[marker],
        params=HostTrimParams(window_bp=10, step_bp=10, buffer_kb=0),
    )

    assert summary["reason"] == "marker_windows"
    assert (trimmed.start, trimmed.end) == (10, 20)


def test_seed_midpoint_at_window_right_edge_belongs_only_to_next_window() -> None:
    seed = MergedSeed(scaffold="ctg", start=0, end=20)
    records = [_viral_record("left", 0, 10), _viral_record("right", 10, 20)]

    trimmed, summary = trim_seed_by_host_signature(
        seed=seed,
        gene_records=records,
        host_model=_host_model(),
        validated_markers=[],
        params=HostTrimParams(window_bp=10, step_bp=10, buffer_kb=0),
    )

    assert summary["reason"] == "midpoint_windows"
    assert (trimmed.start, trimmed.end) == (10, 20)
