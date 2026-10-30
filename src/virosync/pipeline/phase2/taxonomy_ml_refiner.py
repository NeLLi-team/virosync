"""Host-model-independent per-gene ML evidence for the extension controller."""

from __future__ import annotations

import csv
import json
import logging
import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, Literal, Protocol

from virosync.pipeline.phase2.boundary_diamond import (
    MIN_VIRAL_HIT_PIDENT,
    VIRAL_PREFIXES,
)
from virosync.pipeline.taxonomy_utils import compute_hit_weight

if TYPE_CHECKING:
    import numpy as np
    from numpy.typing import NDArray

    from virosync.pipeline.phase2.boundary_diamond import GeneTaxonomy, GenomeDiamondQuery

logger = logging.getLogger(__name__)

TaxonomyMLModelType = Literal["logreg", "gbdt", "xgboost"]
FEATURE_COLUMNS = [
    "has_hit",
    "has_viral",
    "has_ncldv_mirus",
    "has_vp_plv",
    "top1_pident",
    "top1_evalue_log10",
    "top1_is_host",
    "host_rank_score",
    "viral_rank_score",
    "top10_hit_count",
    "mean_top10_bits",
    "viral_top10_fraction",
    "host_top10_fraction",
    "neighbor_viral_fraction_mean",
    "neighbor_host_top1_fraction",
]


class TaxonomyClassifier(Protocol):
    """Operations needed from the selected scikit-learn-compatible model."""

    def fit(self, matrix: list[list[float]], labels: list[int]) -> object:
        """Fit the training matrix."""
        ...

    def predict_proba(self, matrix: list[list[float]]) -> NDArray[np.float64]:
        """Return negative and positive class probabilities."""
        ...


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


def predict_taxonomy_support(
    taxonomy_map: dict[str, GeneTaxonomy],
    boundary_query: GenomeDiamondQuery,
    model_type: TaxonomyMLModelType = "logreg",
    host_prefix: str = "EUK__",
    probability_threshold: float = 0.5,
    neighbor_window: int = 3,
    output_dir: Path | None = None,
    random_state: int = 42,
    taxonomy_weight_mode: str = "rank",
) -> dict[str, frozenset[str]]:
    """Return host-model-independent per-gene support, never coordinates."""
    rows = _build_rows(boundary_query, taxonomy_map, host_prefix, taxonomy_weight_mode, neighbor_window)
    model, summary = _fit_model(_select_training_rows(rows), model_type, random_state)
    if model is None:
        logger.info("Taxonomy ML unavailable (%s); using deterministic extension", summary["fallback_reason"])
        return {}
    prediction_rows = [row for row in rows if row.role != "control"]
    if not prediction_rows:
        return {}
    scores = _predict_scores(model, prediction_rows)
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        _write_feature_table(output_dir / "feature_table.tsv", rows)
        _write_predictions(output_dir / "gene_predictions.tsv", prediction_rows, scores, probability_threshold)
        (output_dir / "model_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return {
        seed_id: frozenset(
            porf_id for (sid, porf_id), score in scores.items() if sid == seed_id and score >= probability_threshold
        )
        for seed_id in boundary_query.seed_gene_mappings
    }


def _build_rows(
    boundary_query: GenomeDiamondQuery,
    taxonomy_map: dict[str, GeneTaxonomy],
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
    taxonomy_map: dict[str, GeneTaxonomy],
    host_prefix: str,
    taxonomy_weight_mode: str,
) -> None:
    for porf_id in porf_ids:
        if porf_id in seen or porf_id not in taxonomy_map:
            continue
        seen.add(porf_id)
        rows.append(
            _build_row(
                seed_id,
                porf_id,
                role,
                label,
                taxonomy_map[porf_id],
                host_prefix,
                taxonomy_weight_mode,
            )
        )


def _build_row(
    seed_id: str,
    porf_id: str,
    role: str,
    label: int,
    tax: GeneTaxonomy,
    host_prefix: str,
    taxonomy_weight_mode: str,
) -> GeneFeatureRow:
    host_rank, viral_rank = _ranked_prefix_scores(tax, host_prefix, taxonomy_weight_mode)
    evalue = max(float(getattr(tax, "top1_evalue", 1.0) or 1.0), 1e-300)
    bits = [float(v) for v in (getattr(tax, "top10_bits", []) or [])]
    prefixes = [str(p) for p in (getattr(tax, "top10_prefixes", []) or [])]
    pidents = [float(v) for v in (getattr(tax, "top10_pidents", []) or [])]
    top10_count = float(len(prefixes))
    viral_fraction = (
        sum(
            1
            for idx, prefix in enumerate(prefixes)
            if prefix in VIRAL_PREFIXES and idx < len(pidents) and pidents[idx] >= MIN_VIRAL_HIT_PIDENT
        )
        / top10_count
        if top10_count > 0
        else 0.0
    )
    host_fraction = sum(1 for prefix in prefixes if prefix == host_prefix) / top10_count if top10_count > 0 else 0.0
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
        "top10_hit_count": top10_count,
        "mean_top10_bits": (sum(bits) / len(bits)) if bits else 0.0,
        "viral_top10_fraction": viral_fraction,
        "host_top10_fraction": host_fraction,
        "neighbor_viral_fraction_mean": 0.0,
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


def _ranked_prefix_scores(tax: GeneTaxonomy, host_prefix: str, mode: str) -> tuple[float, float]:
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
        host_top1_mean = sum(r.features["top1_is_host"] for r in neighborhood) / size
        features = dict(row.features)
        features["neighbor_viral_fraction_mean"] = viral_mean
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
) -> tuple[TaxonomyClassifier | None, dict[str, object]]:
    labels = [row.label for row in training_rows]
    n_pos = sum(1 for v in labels if v == 1)
    n_neg = sum(1 for v in labels if v == 0)
    summary: dict[str, object] = {
        "n_training_genes": len(training_rows),
        "n_positive_genes": n_pos,
        "n_negative_genes": n_neg,
        "used_model_type": "unavailable",
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
) -> tuple[TaxonomyClassifier, str, str]:
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

    from sklearn.ensemble import GradientBoostingClassifier

    return GradientBoostingClassifier(random_state=random_state), "gbdt", ""


def _predict_scores(model: TaxonomyClassifier, rows: list[GeneFeatureRow]) -> dict[tuple[str, str], float]:
    matrix = [[row.features[name] for name in FEATURE_COLUMNS] for row in rows]
    probs = model.predict_proba(matrix)[:, 1]
    return {(row.seed_id, row.porf_id): float(prob) for row, prob in zip(rows, probs)}


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


__all__ = ["TaxonomyMLModelType", "predict_taxonomy_support"]
