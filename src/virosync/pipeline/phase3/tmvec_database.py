"""
TMVec database search for structural similarity.

Loads manifest-bound BFVD TMVec2 embeddings and runs cosine-similarity search
against query embeddings.
"""

from __future__ import annotations

import csv
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from virosync.utils.database_manager import (
    TMVEC_EMBEDDING_WIDTH,
    ViroSyncDatabaseManager,
)

from .tmvec_predictor import get_tmvec_predictor

logger = logging.getLogger(__name__)


def build_db_paths(database_root: Path, manifest: dict) -> dict[str, dict[str, Path]]:
    """Resolve hash-bound BFVD paths from the TMVec2 manifest."""
    bfvd = manifest["databases"]["bfvd"]
    return {
        "bfvd": {
            "embeddings": database_root / bfvd["embeddings"]["path"],
            "metadata": database_root / bfvd["metadata"]["path"],
        }
    }


@dataclass
class TMVecHit:
    target_id: str
    tm_score: float
    database: str
    protein_name: Optional[str] = None
    organism: Optional[str] = None
    lineage: Optional[str] = None
    keywords: Optional[str] = None


class TMVecDatabaseSearch:
    """Search precomputed TMVec databases for structural similarity."""

    def __init__(
        self,
        device: str = "cuda",
        databases: Optional[list[str]] = None,
        min_tm: float = 0.0,
        database_root: Optional[Path] = None,
        require_gpu: bool = False,
        fail_on_unavailable: bool = False,
    ) -> None:
        self.device = device
        self.databases = databases or ["bfvd"]
        unsupported = sorted(set(self.databases) - {"bfvd"})
        if unsupported:
            raise ValueError(
                "Unsupported TMVec2 database key(s): " + ", ".join(unsupported)
            )
        self.min_tm = min_tm
        if database_root is None:
            raise ValueError(
                "TMVec database_root is required; resolve phase3.tmvec_database_dir "
                "before constructing TMVecDatabaseSearch"
            )
        self.database_root = Path(database_root)
        self._manifest = None
        self._db_paths = None
        self._db_cache: dict[str, dict] = {}
        self._predictor = None
        self.require_gpu = require_gpu
        self.fail_on_unavailable = fail_on_unavailable

    @property
    def _must_fail(self) -> bool:
        return self.require_gpu or self.fail_on_unavailable

    def _load_manifest(self) -> dict:
        if self._manifest is None:
            self._manifest = ViroSyncDatabaseManager.load_tmvec_manifest(
                self.database_root,
                verify_hashes=True,
                databases=self.databases,
            )
            self._db_paths = build_db_paths(self.database_root, self._manifest)
        return self._manifest

    @staticmethod
    def _deduplicate_proteins(
        proteins: list[tuple[str, str]],
    ) -> tuple[list[tuple[str, str]], int]:
        """Deduplicate proteins by pORF ID while preserving input order."""
        unique: list[tuple[str, str]] = []
        seen_sequences: dict[str, str] = {}
        conflicting = 0
        for porf_id, sequence in proteins:
            prior = seen_sequences.get(porf_id)
            if prior is None:
                seen_sequences[porf_id] = sequence
                unique.append((porf_id, sequence))
            elif prior != sequence:
                conflicting += 1
        return unique, conflicting

    @property
    def predictor(self):
        if self._predictor is None:
            try:
                self._load_manifest()
                self._predictor = get_tmvec_predictor(
                    device=self.device,
                    model_root=self.database_root / "models",
                    require_gpu=self.require_gpu,
                    fail_on_unavailable=self.fail_on_unavailable,
                )
                if not self._predictor.available:
                    logger.warning("TMVec predictor not available on device=%s", self.device)
                    if self._must_fail:
                        raise RuntimeError(
                            f"TMVec predictor unavailable on device={self.device}."
                        )
                    self._predictor = None
            except RuntimeError as e:
                logger.warning("TMVec predictor initialization failed: %s", e)
                if self._must_fail:
                    raise
                self._predictor = None
            except Exception as e:
                logger.error("Unexpected error initializing TMVec predictor: %s", e)
                if self._must_fail:
                    raise
                self._predictor = None
        return self._predictor

    def _load_db(self, name: str) -> Optional[dict]:
        if name in self._db_cache:
            return self._db_cache[name]

        try:
            self._load_manifest()
        except Exception as exc:
            if self._must_fail:
                raise
            logger.warning("TMVec2 resource manifest is invalid: %s", exc)
            return None
        assert self._db_paths is not None
        paths = self._db_paths.get(name)
        if not paths:
            if self._must_fail:
                raise RuntimeError(f"TMVec database {name} is not configured")
            logger.warning("TMVec database %s not configured", name)
            return None

        emb_path = paths.get("embeddings")
        if not emb_path or not emb_path.is_file():
            if self._must_fail:
                raise RuntimeError(
                    f"TMVec database {name} embeddings not found: {emb_path}"
                )
            logger.warning("TMVec database %s embeddings not found: %s", name, emb_path)
            return None

        embeddings = np.load(emb_path, mmap_mode="r")
        if (
            embeddings.ndim != 2
            or any(size < 1 for size in embeddings.shape)
            or embeddings.shape[1] != TMVEC_EMBEDDING_WIDTH
            or embeddings.dtype.kind not in "iuf"
        ):
            message = (
                f"TMVec database {name} has invalid embeddings "
                f"shape={embeddings.shape} dtype={embeddings.dtype}"
            )
            if self._must_fail:
                raise RuntimeError(message)
            logger.warning("%s", message)
            return None
        logger.info("TMVec database %s: embeddings shape=%s", name, embeddings.shape)

        metadata_path = paths.get("metadata")
        if not metadata_path or not metadata_path.is_file():
            if self._must_fail:
                raise RuntimeError(f"TMVec BFVD metadata missing: {metadata_path}")
            logger.warning("TMVec BFVD metadata missing: %s", metadata_path)
            return None
        annotations: list[dict] = []
        with metadata_path.open("r", encoding="utf-8", newline="") as handle:
            if metadata_path.suffix == ".tsv":
                annotations = list(csv.DictReader(handle, delimiter="\t"))
            elif metadata_path.suffix == ".jsonl":
                annotations = [json.loads(line) for line in handle if line.strip()]
            else:
                raise RuntimeError(
                    "TMVec BFVD metadata must use TSV or JSONL, not pickle"
                )
        ids = [str(item["id"]) for item in annotations]

        if len(ids) != embeddings.shape[0]:
            message = (
                f"TMVec database {name} embedding/metadata row count mismatch: "
                f"{embeddings.shape[0]} != {len(ids)}"
            )
            if self._must_fail:
                raise RuntimeError(message)
            logger.warning("%s", message)
            return None

        db = {
            "embeddings": embeddings,
            "ids": ids,
            "annotations": annotations,
        }
        self._db_cache[name] = db
        return db

    def search_sequence(
        self,
        sequence: str,
        databases: Optional[list[str]] = None,
        top_k: int = 1,
    ) -> dict[str, TMVecHit | None]:
        predictor = self.predictor
        if predictor is None or not predictor.available:
            if self._must_fail:
                raise RuntimeError("TMVec predictor unavailable after preflight")
            return {name: None for name in (databases or self.databases)}

        db_names = databases or self.databases
        hits: dict[str, TMVecHit | None] = {}

        for name in db_names:
            db = self._load_db(name)
            if not db:
                hits[name] = None
                continue

            try:
                results = predictor.search_database(
                    sequence,
                    db["embeddings"],
                    db["ids"],
                    top_k=top_k,
                    min_tm=self.min_tm,
                )
            except Exception as exc:
                logger.warning("TMVec search failed for %s: %s", name, exc)
                if self._must_fail:
                    raise
                hits[name] = None
                continue

            if not results:
                hits[name] = None
                continue

            target_id, score = results[0]
            annotations = db.get("annotations")
            protein_name = organism = lineage = keywords = None
            if annotations is not None:
                try:
                    idx = db["ids"].index(target_id)
                except ValueError:
                    idx = -1
                if idx >= 0:
                    info = annotations[idx]
                    protein_name = info.get("protein_name") or None
                    organism = info.get("organism") or None
                    lineage = info.get("lineage") or None
                    keywords = info.get("keywords") or None

            hits[name] = TMVecHit(
                target_id=target_id,
                tm_score=float(score),
                database=name,
                protein_name=protein_name,
                organism=organism,
                lineage=lineage,
                keywords=keywords,
            )

        return hits

    def search_batch(
        self,
        proteins: list[tuple[str, str]],
        databases: Optional[list[str]] = None,
    ) -> dict[str, dict[str, TMVecHit | None]]:
        """
        Batch search for multiple proteins at once.

        This is much more efficient than calling search_sequence for each protein
        because it:
        1. Loads the predictor model only once
        2. Batch embeds all proteins together (GPU batch size managed internally)
        3. Searches each database via efficient matrix multiplication

        Args:
            proteins: List of (porf_id, sequence) tuples
            databases: Database names to search (default: self.databases)

        Returns:
            Dict mapping porf_id -> {db_name -> TMVecHit | None}
        """
        if not proteins:
            return {}

        raw_count = len(proteins)
        proteins, conflicting_ids = self._deduplicate_proteins(proteins)
        if len(proteins) != raw_count:
            logger.info(
                "TMVec batch: deduplicated proteins by pORF ID (%d -> %d)",
                raw_count,
                len(proteins),
            )
        if conflicting_ids:
            logger.warning(
                "TMVec batch: %d duplicate pORF IDs had conflicting sequences; using first occurrence",
                conflicting_ids,
            )

        predictor = self.predictor
        if predictor is None or not predictor.available:
            if self._must_fail:
                raise RuntimeError("TMVec predictor unavailable after preflight")
            db_names = databases or self.databases
            return {pid: {name: None for name in db_names} for pid, _ in proteins}

        db_names = databases or self.databases
        porf_ids = [pid for pid, _ in proteins]
        sequences = [seq for _, seq in proteins]

        logger.info("TMVec batch: embedding %d proteins with %s...",
                   len(sequences), type(predictor).__name__)

        # Embed all proteins with the predictor's token-budgeted batches.
        try:
            query_embeddings = np.asarray(predictor.embed_batch(sequences))
            logger.info("TMVec batch: query embeddings shape=%s", query_embeddings.shape)
        except Exception as exc:
            logger.error("TMVec batch embedding failed: %s", exc)
            if self._must_fail:
                raise
            return {pid: {name: None for name in db_names} for pid in porf_ids}

        # Initialize results
        results: dict[str, dict[str, TMVecHit | None]] = {
            pid: {} for pid in porf_ids
        }

        if query_embeddings.ndim != 2 or query_embeddings.shape[0] != len(porf_ids):
            msg = (
                "TMVec batch embedding returned invalid shape "
                f"{query_embeddings.shape}; expected ({len(porf_ids)}, n_dims)"
            )
            if self._must_fail:
                raise RuntimeError(msg)
            logger.warning("%s; disabling TMVec hits for this batch", msg)
            return {pid: {name: None for name in db_names} for pid in porf_ids}

        embedding_norms = np.linalg.norm(query_embeddings, axis=1)
        valid_embeddings = (
            np.isfinite(query_embeddings).all(axis=1)
            & np.isfinite(embedding_norms)
            & (embedding_norms > 1e-8)
        )
        invalid_porf_ids = [
            porf_id
            for porf_id, is_valid in zip(porf_ids, valid_embeddings)
            if not bool(is_valid)
        ]
        if invalid_porf_ids:
            msg = (
                "TMVec produced invalid or zero embeddings for "
                f"{len(invalid_porf_ids)} protein(s)"
            )
            if self._must_fail:
                raise RuntimeError(msg)
            logger.warning("%s; those proteins will have no TMVec hits", msg)

        if not bool(valid_embeddings.any()):
            return {pid: {name: None for name in db_names} for pid in porf_ids}

        # Normalize query embeddings for cosine similarity. Invalid rows stay
        # zero and are skipped below.
        query_norms = np.zeros_like(query_embeddings, dtype=np.float32)
        query_norms[valid_embeddings] = (
            query_embeddings[valid_embeddings]
            / embedding_norms[valid_embeddings, None]
        )

        # Search each database sequentially
        for db_name in db_names:
            logger.info("TMVec batch: searching %s database...", db_name)
            db = self._load_db(db_name)
            if not db:
                for pid in porf_ids:
                    results[pid][db_name] = None
                continue

            try:
                # Normalize database embeddings
                db_embeddings = db["embeddings"]
                db_norms = db_embeddings / (
                    np.linalg.norm(db_embeddings, axis=1, keepdims=True) + 1e-8
                )

                # Compute all similarities: [n_queries, n_db]
                similarities = query_norms @ db_norms.T

                # Get annotations if available
                annotations = db.get("annotations")
                db_ids = db["ids"]

                # Extract top hit for each protein
                for i, porf_id in enumerate(porf_ids):
                    if not bool(valid_embeddings[i]):
                        results[porf_id][db_name] = None
                        continue

                    top_idx = int(np.argmax(similarities[i]))
                    top_score = float(similarities[i, top_idx])

                    if top_score < self.min_tm:
                        results[porf_id][db_name] = None
                        continue

                    target_id = db_ids[top_idx]
                    protein_name = organism = lineage = keywords = None

                    if annotations is not None and top_idx < len(annotations):
                        info = annotations[top_idx]
                        protein_name = info.get("protein_name") or None
                        organism = info.get("organism") or None
                        lineage = info.get("lineage") or None
                        keywords = info.get("keywords") or None

                    results[porf_id][db_name] = TMVecHit(
                        target_id=target_id,
                        tm_score=top_score,
                        database=db_name,
                        protein_name=protein_name,
                        organism=organism,
                        lineage=lineage,
                        keywords=keywords,
                    )

            except Exception as exc:
                logger.warning("TMVec batch search failed for %s: %s", db_name, exc)
                if self._must_fail:
                    raise
                for pid in porf_ids:
                    results[pid][db_name] = None

        logger.info("TMVec batch: completed %d proteins across %d databases",
                   len(proteins), len(db_names))
        return results
