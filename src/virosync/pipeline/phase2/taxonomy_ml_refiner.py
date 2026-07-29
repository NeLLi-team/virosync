"""Phase 2 taxonomy-based ML boundary refinement."""

from __future__ import annotations

import csv
import json
import logging
import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, Literal, Optional

from virosync.pipeline.host_signatures import HostSignatureModel, score_host_signature_record
from virosync.pipeline.phase1.seed_merger import MergedSeed
from virosync.pipeline.phase2.boundary_diamond import (
    MIN_VIRAL_HIT_PIDENT,
    VIRAL_PREFIXES,
)
from virosync.pipeline.taxonomy_utils import compute_hit_weight

if TYPE_CHECKING:
    from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy, GenomeDiamondQuery, SeedGeneMapping

logger = logging.getLogger(__name__)

TaxonomyMLModelType = Literal["logreg", "gbdt", "xgboost"]
FEATURE_COLUMNS = [
    "has_hit", "has_viral", "has_ncldv_mirus", "has_vp_plv", "top1_pident",
    "top1_evalue_log10", "top1_is_host", "host_rank_score", "viral_rank_score",
    "host_signature_score", "top10_hit_count", "mean_top10_bits",
    "viral_top10_fraction", "host_top10_fraction",
    "neighbor_viral_fraction_mean", "neighbor_host_score_mean", "neighbor_host_top1_fraction",
]


@dataclass(frozen=True)
class GeneFeatureRow:
    seed_id: str
    porf_id: str
    role: str
    label: int
    scaffold: str
    start: int
    end: int
    gene_order_index: int
    features: dict[str, float]


def refine_seeds_by_taxonomy_ml(
    merged_seeds: list[MergedSeed],
    taxonomy_map: dict[str, "GeneTaxonomy"],
    boundary_query: "GenomeDiamondQuery",
    host_signature_model: Optional[HostSignatureModel] = None,
    model_type: TaxonomyMLModelType = "logreg",
    host_prefix: str = "EUK__",
    probability_threshold: float = 0.5,
    neighbor_window: int = 3,
    output_dir: Optional[Path] = None,
    random_state: int = 42,
    taxonomy_weight_mode: str = "rank",
    host_signature_threshold: float = 0.5,
    host_guard_viral_rank_min: float = 8.0,
) -> list[MergedSeed]:
    """Train per-genome taxonomy model and refine seeds within flank envelopes."""
    if not merged_seeds or not taxonomy_map or not boundary_query.seed_gene_mappings:
        return merged_seeds

    rows = _build_rows(
        boundary_query,
        taxonomy_map,
        host_signature_model,
        host_prefix,
        taxonomy_weight_mode,
        neighbor_window,
    )
    train_rows = _select_training_rows(rows)
    model, model_summary = _fit_model(train_rows, model_type, random_state)

    prediction_rows = [row for row in rows if row.role != "control"]
    scores = _predict_scores(model, prediction_rows)
    refined_seeds, refine_summary = _refine_seed_boundaries(
        merged_seeds,
        boundary_query.seed_gene_mappings,
        prediction_rows,
        scores,
        probability_threshold,
        host_signature_threshold=host_signature_threshold,
        host_guard_viral_rank_min=host_guard_viral_rank_min,
    )

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        _write_feature_table(out / "feature_table.tsv", rows)
        _write_predictions(out / "gene_predictions.tsv", prediction_rows, scores, probability_threshold)
        summary = {
            "requested_model_type": model_type,
            "used_model_type": model_summary["used_model_type"],
            "fallback_reason": model_summary["fallback_reason"],
            "feature_columns": FEATURE_COLUMNS,
            "threshold": probability_threshold,
            "neighbor_window": neighbor_window,
            "host_signature_threshold": host_signature_threshold,
            "host_guard_viral_rank_min": host_guard_viral_rank_min,
            "training": model_summary,
            "refinement": refine_summary,
            "n_feature_rows": len(rows),
            "n_prediction_rows": len(prediction_rows),
        }
        (out / "model_summary.json").write_text(json.dumps(summary, indent=2))

    logger.info(
        "Taxonomy ML refinement updated %d/%d seeds (%s)",
        refine_summary["changed"],
        len(merged_seeds),
        model_summary["used_model_type"],
    )
    return refined_seeds


def _build_rows(
    boundary_query: "GenomeDiamondQuery",
    taxonomy_map: dict[str, "GeneTaxonomy"],
    host_model: Optional[HostSignatureModel],
    host_prefix: str,
    taxonomy_weight_mode: str,
    neighbor_window: int,
) -> list[GeneFeatureRow]:
    rows: list[GeneFeatureRow] = []
    for seed_id, mapping in boundary_query.seed_gene_mappings.items():
        seen: set[str] = set()
        seed_rows: list[GeneFeatureRow] = []
        _add_seed_role_rows(
            seed_rows,
            seed_id,
            mapping.eve_porf_ids,
            "seed_interior",
            1,
            seen,
            taxonomy_map,
            host_model,
            host_prefix,
            taxonomy_weight_mode,
        )
        _add_seed_role_rows(
            seed_rows,
            seed_id,
            mapping.upstream_porf_ids,
            "upstream",
            -1,
            seen,
            taxonomy_map,
            host_model,
            host_prefix,
            taxonomy_weight_mode,
        )
        _add_seed_role_rows(
            seed_rows,
            seed_id,
            mapping.downstream_porf_ids,
            "downstream",
            -1,
            seen,
            taxonomy_map,
            host_model,
            host_prefix,
            taxonomy_weight_mode,
        )
        seed_rows.sort(key=lambda row: (row.start, row.end, row.porf_id))
        seed_rows = _annotate_neighbor_features(seed_rows, neighbor_window)
        for order_index, row in enumerate(seed_rows):
            rows.append(replace(row, gene_order_index=order_index))

    seen_control: set[str] = set()
    for porf_id in boundary_query.control_porf_ids:
        if porf_id in seen_control or porf_id not in taxonomy_map:
            continue
        seen_control.add(porf_id)
        tax = taxonomy_map[porf_id]
        rows.append(
            _build_row(
                "__control__",
                porf_id,
                "control",
                0,
                tax,
                host_model,
                host_prefix,
                taxonomy_weight_mode,
            )
        )
    return rows


def _add_seed_role_rows(
    rows: list[GeneFeatureRow],
    seed_id: str,
    porf_ids: list[str],
    role: str,
    label: int,
    seen: set[str],
    taxonomy_map: dict[str, "GeneTaxonomy"],
    host_model: Optional[HostSignatureModel],
    host_prefix: str,
    taxonomy_weight_mode: str,
) -> None:
    for porf_id in porf_ids:
        if porf_id in seen or porf_id not in taxonomy_map:
            continue
        seen.add(porf_id)
        rows.append(
            _build_row(
                seed_id, porf_id, role, label, taxonomy_map[porf_id],
                host_model, host_prefix, taxonomy_weight_mode,
            )
        )


def _build_row(
    seed_id: str,
    porf_id: str,
    role: str,
    label: int,
    tax: "GeneTaxonomy",
    host_model: Optional[HostSignatureModel],
    host_prefix: str,
    taxonomy_weight_mode: str,
) -> GeneFeatureRow:
    host_rank, viral_rank = _ranked_prefix_scores(tax, host_prefix, taxonomy_weight_mode)
    evalue = max(float(getattr(tax, "top1_evalue", 1.0) or 1.0), 1e-300)
    bits = [float(v) for v in (getattr(tax, "top10_bits", []) or [])]
    prefixes = [str(p) for p in (getattr(tax, "top10_prefixes", []) or [])]
    pidents = [float(v) for v in (getattr(tax, "top10_pidents", []) or [])]
    host_score = score_host_signature_record(tax, host_model) if host_model is not None else 0.0
    top10_count = float(len(prefixes))
    viral_fraction = (
        sum(
            1
            for idx, prefix in enumerate(prefixes)
            if prefix in VIRAL_PREFIXES
            and idx < len(pidents)
            and pidents[idx] >= MIN_VIRAL_HIT_PIDENT
        ) / top10_count
        if top10_count > 0
        else 0.0
    )
    host_fraction = (
        sum(1 for prefix in prefixes if prefix == host_prefix) / top10_count
        if top10_count > 0
        else 0.0
    )
    features = {
        "has_hit": float(bool(getattr(tax, "has_hit", False))),
        "has_viral": float(bool(getattr(tax, "has_viral", False))),
        "has_ncldv_mirus": float(bool(getattr(tax, "has_ncldv_mirus", False))),
        "has_vp_plv": float(bool(getattr(tax, "has_vp_plv", False))),
        "top1_pident": float(getattr(tax, "top1_pident", 0.0) or 0.0),
        "top1_evalue_log10": -math.log10(evalue),
        "top1_is_host": float(getattr(tax, "top1_prefix", "") == host_prefix),
        "host_rank_score": host_rank,
        "viral_rank_score": viral_rank,
        "host_signature_score": float(host_score),
        "top10_hit_count": top10_count,
        "mean_top10_bits": (sum(bits) / len(bits)) if bits else 0.0,
        "viral_top10_fraction": viral_fraction,
        "host_top10_fraction": host_fraction,
        "neighbor_viral_fraction_mean": 0.0,
        "neighbor_host_score_mean": 0.0,
        "neighbor_host_top1_fraction": 0.0,
    }
    return GeneFeatureRow(
        seed_id,
        porf_id,
        role,
        label,
        tax.scaffold,
        tax.start,
        tax.end,
        -1,
        features,
    )


def _ranked_prefix_scores(tax: "GeneTaxonomy", host_prefix: str, mode: str) -> tuple[float, float]:
    prefixes = getattr(tax, "top10_prefixes", []) or []
    if not prefixes:
        host = 1.0 if getattr(tax, "top1_prefix", "") == host_prefix else 0.0
        viral = 1.0 if bool(getattr(tax, "has_viral", False)) else 0.0
        return host, viral

    bits = getattr(tax, "top10_bits", []) or []
    pidents = getattr(tax, "top10_pidents", []) or []
    host_score = 0.0
    viral_score = 0.0
    for rank, prefix in enumerate(prefixes):
        bit = float(bits[rank]) if rank < len(bits) else 0.0
        weight = compute_hit_weight(rank, bit, mode)
        if prefix == host_prefix:
            host_score += weight
        if prefix in VIRAL_PREFIXES:
            try:
                pident = float(pidents[rank])
            except (IndexError, TypeError, ValueError):
                pident = 0.0
            if pident >= MIN_VIRAL_HIT_PIDENT:
                viral_score += weight
    return host_score, viral_score


def _annotate_neighbor_features(
    rows: list[GeneFeatureRow],
    neighbor_window: int,
) -> list[GeneFeatureRow]:
    if not rows:
        return rows
    window = max(1, int(neighbor_window))
    updated_rows: list[GeneFeatureRow] = []
    n_rows = len(rows)
    for idx, row in enumerate(rows):
        left = max(0, idx - window)
        right = min(n_rows, idx + window + 1)
        neighborhood = rows[left:right]
        size = float(len(neighborhood)) or 1.0
        viral_mean = sum(r.features["viral_top10_fraction"] for r in neighborhood) / size
        host_score_mean = sum(r.features["host_signature_score"] for r in neighborhood) / size
        host_top1_mean = sum(r.features["top1_is_host"] for r in neighborhood) / size
        features = dict(row.features)
        features["neighbor_viral_fraction_mean"] = viral_mean
        features["neighbor_host_score_mean"] = host_score_mean
        features["neighbor_host_top1_fraction"] = host_top1_mean
        updated_rows.append(replace(row, features=features))
    return updated_rows


def _select_training_rows(rows: list[GeneFeatureRow]) -> list[GeneFeatureRow]:
    dedup: dict[str, GeneFeatureRow] = {}
    for row in rows:
        if row.label < 0:
            continue
        previous = dedup.get(row.porf_id)
        if previous is None or row.label > previous.label:
            dedup[row.porf_id] = row
    return list(dedup.values())


def _fit_model(
    training_rows: list[GeneFeatureRow],
    model_type: TaxonomyMLModelType,
    random_state: int,
) -> tuple[Optional[object], dict[str, object]]:
    labels = [row.label for row in training_rows]
    n_pos = sum(1 for v in labels if v == 1)
    n_neg = sum(1 for v in labels if v == 0)
    summary: dict[str, object] = {
        "n_training_genes": len(training_rows),
        "n_positive_genes": n_pos,
        "n_negative_genes": n_neg,
        "used_model_type": "heuristic",
        "fallback_reason": "",
        "metrics": {},
    }
    if n_pos == 0 or n_neg == 0:
        summary["fallback_reason"] = "insufficient_classes"
        return None, summary

    matrix = [[row.features[name] for name in FEATURE_COLUMNS] for row in training_rows]
    model, used, fallback = _build_estimator(model_type, random_state, n_pos, n_neg)
    model.fit(matrix, labels)
    summary["used_model_type"] = used
    summary["fallback_reason"] = fallback

    try:
        from sklearn.metrics import accuracy_score, roc_auc_score

        probs = model.predict_proba(matrix)[:, 1]
        preds = [1 if p >= 0.5 else 0 for p in probs]
        summary["metrics"] = {
            "train_accuracy": float(accuracy_score(labels, preds)),
            "train_auc": float(roc_auc_score(labels, probs)),
        }
    except Exception:
        summary["metrics"] = {}
    return model, summary


def _build_estimator(
    model_type: TaxonomyMLModelType,
    random_state: int,
    n_pos: int,
    n_neg: int,
) -> tuple[object, str, str]:
    if model_type == "logreg":
        from sklearn.linear_model import LogisticRegression
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler

        model = make_pipeline(
            StandardScaler(),
            LogisticRegression(class_weight="balanced", max_iter=1000, solver="liblinear", random_state=random_state),
        )
        return model, "logreg", ""

    if model_type == "xgboost":
        try:
            from xgboost import XGBClassifier

            ratio = float(n_neg) / float(max(1, n_pos))
            model = XGBClassifier(
                n_estimators=200,
                max_depth=3,
                learning_rate=0.1,
                subsample=0.8,
                colsample_bytree=0.8,
                objective="binary:logistic",
                eval_metric="logloss",
                scale_pos_weight=ratio,
                random_state=random_state,
                n_jobs=1,
            )
            return model, "xgboost", ""
        except Exception:
            logger.warning("xgboost unavailable; using GradientBoostingClassifier")

    from sklearn.ensemble import GradientBoostingClassifier

    fallback = "xgboost_unavailable" if model_type == "xgboost" else ""
    return GradientBoostingClassifier(random_state=random_state), "gbdt", fallback


def _predict_scores(model: Optional[object], rows: list[GeneFeatureRow]) -> dict[tuple[str, str], float]:
    if model is None:
        return {(row.seed_id, row.porf_id): _heuristic_probability(row) for row in rows}
    matrix = [[row.features[name] for name in FEATURE_COLUMNS] for row in rows]
    probs = model.predict_proba(matrix)[:, 1]
    return {(row.seed_id, row.porf_id): float(prob) for row, prob in zip(rows, probs)}


def _heuristic_probability(row: GeneFeatureRow) -> float:
    viral = row.features["viral_rank_score"] + row.features["has_viral"]
    host = row.features["host_rank_score"] + row.features["host_signature_score"]
    score = viral / (viral + host) if (viral + host) > 0 else row.features["has_viral"]
    if row.features["top1_is_host"] > 0 and row.features["has_viral"] <= 0:
        score *= 0.5
    return max(0.0, min(1.0, score))


def _refine_seed_boundaries(
    merged_seeds: list[MergedSeed],
    seed_mappings: dict[str, "SeedGeneMapping"],
    prediction_rows: list[GeneFeatureRow],
    scores: dict[tuple[str, str], float],
    threshold: float,
    host_signature_threshold: float = 0.5,
    host_guard_viral_rank_min: float = 8.0,
) -> tuple[list[MergedSeed], dict[str, int]]:
    mapping_by_coords = {
        (m.scaffold, int(m.seed_start), int(m.seed_end)): key
        for key, m in seed_mappings.items()
    }

    rows_by_seed: dict[str, list[GeneFeatureRow]] = {}
    for row in prediction_rows:
        rows_by_seed.setdefault(row.seed_id, []).append(row)

    refined: list[MergedSeed] = []
    changed = 0
    extended = 0
    contracted = 0
    for seed in merged_seeds:
        seed_key = seed.seed_id or mapping_by_coords.get((seed.scaffold, seed.start, seed.end), "")
        mapping = seed_mappings.get(seed_key)
        rows = sorted(rows_by_seed.get(seed_key, []), key=lambda r: (r.start, r.end))
        if mapping is None or not rows:
            refined.append(seed)
            continue

        positive_candidates = [r for r in rows if scores.get((r.seed_id, r.porf_id), 0.0) >= threshold]
        positives = [
            row
            for row in positive_candidates
            if _passes_host_guard(row, host_signature_threshold, host_guard_viral_rank_min)
        ]
        if not positives:
            refined.append(seed)
            continue

        new_start = max(mapping.flank_start_bp, min(r.start for r in positives))
        new_end = min(mapping.flank_end_bp, max(r.end for r in positives))
        if new_start >= new_end:
            refined.append(seed)
            continue

        if new_start != seed.start or new_end != seed.end:
            changed += 1
            if new_start < seed.start or new_end > seed.end:
                extended += 1
            if new_start > seed.start or new_end < seed.end:
                contracted += 1
        refined.append(replace(seed, start=new_start, end=new_end))

    return refined, {
        "changed": changed,
        "extended": extended,
        "contracted": contracted,
    }


def _passes_host_guard(
    row: GeneFeatureRow,
    host_signature_threshold: float,
    host_guard_viral_rank_min: float,
) -> bool:
    """
    Block host-like positives unless they also have strong viral evidence.
    """
    host_score = float(row.features.get("host_signature_score", 0.0))
    if host_score < host_signature_threshold:
        return True
    if float(row.features.get("has_viral", 0.0)) > 0.0:
        return True
    viral_rank = float(row.features.get("viral_rank_score", 0.0))
    return viral_rank >= host_guard_viral_rank_min


def _write_feature_table(path: Path, rows: list[GeneFeatureRow]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "seed_id",
                "porf_id",
                "scaffold",
                "start",
                "end",
                "gene_order_index",
                "is_seed_interior",
                "is_flank_upstream",
                "is_flank_downstream",
                "is_control_gene",
                "train_label",
                *FEATURE_COLUMNS,
            ]
        )
        for row in sorted(rows, key=lambda r: (r.seed_id, r.start, r.end, r.porf_id)):
            values = [
                row.seed_id,
                row.porf_id,
                row.scaffold,
                row.start,
                row.end,
                row.gene_order_index,
                int(row.role == "seed_interior"),
                int(row.role == "upstream"),
                int(row.role == "downstream"),
                int(row.role == "control"),
                row.label,
            ]
            writer.writerow([*values, *[row.features[name] for name in FEATURE_COLUMNS]])


def _write_predictions(
    path: Path,
    rows: list[GeneFeatureRow],
    scores: dict[tuple[str, str], float],
    threshold: float,
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "seed_id",
                "porf_id",
                "scaffold",
                "start",
                "end",
                "gene_order_index",
                "is_seed_interior",
                "is_flank_upstream",
                "is_flank_downstream",
                "is_control_gene",
                "train_label",
                "pred_proba",
                "pred_label",
            ]
        )
        for row in sorted(rows, key=lambda r: (r.seed_id, r.start, r.end, r.porf_id)):
            probability = scores.get((row.seed_id, row.porf_id), 0.0)
            values = [
                row.seed_id,
                row.porf_id,
                row.scaffold,
                row.start,
                row.end,
                row.gene_order_index,
                int(row.role == "seed_interior"),
                int(row.role == "upstream"),
                int(row.role == "downstream"),
                int(row.role == "control"),
                row.label,
                probability,
            ]
            writer.writerow([*values, int(probability >= threshold)])


__all__ = ["TaxonomyMLModelType", "refine_seeds_by_taxonomy_ml"]
