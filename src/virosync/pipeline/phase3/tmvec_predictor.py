"""TMVec2 sequence embeddings and cosine-similarity search.

TMVec2 uses residue features from Lobster-24M and maps them to a 512-value
structure-aware vector. Model files must be installed below the selected
TMVec resource root. This module does not download model files at runtime.

The model architecture and inference path match paarth-b/tmvec-bench commit
59eb75eb75fa1a7524eec59227a83155baadb9b0 and the pinned TMVec-2 model revision.
"""

from __future__ import annotations

import json
import logging
import threading
import weakref
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import torch
import torch.nn as nn

logger = logging.getLogger(__name__)

LOBSTER_MODEL_ID = "asalam91/lobster_24M"
LOBSTER_MODEL_REVISION = "9c36ae05d277e312ac319cbc41b5759472f5bd90"
TMVEC2_MODEL_ID = "scikit-bio/TMVec-2"
TMVEC2_MODEL_REVISION = "91fbaaefbacd72ff6bc2f2126e8a0c165b2a9d92"

_MODEL_LOADING_LOCK = threading.Lock()
_MAX_TOKEN_LENGTH = 512


@dataclass(frozen=True)
class TMvecConfig:
    """TMVec2 inference architecture."""

    d_model: int = 408
    nhead: int = 8
    num_layers: int = 4
    dim_feedforward: int = 2048
    out_dim: int = 512
    dropout: float = 0.2
    activation: str = "gelu"
    max_length: int = _MAX_TOKEN_LENGTH
    projection_hidden_dim: int = 1024

    @classmethod
    def from_json(cls, path: Path) -> "TMvecConfig":
        """Read the architecture fields from the installed TMVec2 parameters."""
        with path.open(encoding="utf-8") as handle:
            data = json.load(handle)
        if not isinstance(data, dict):
            raise RuntimeError("Installed TMVec2 parameters must contain a JSON object")
        fields = set(cls.__dataclass_fields__)
        missing = sorted(fields - set(data))
        if missing:
            raise RuntimeError(
                "Installed TMVec2 parameters are missing: " + ", ".join(missing)
            )
        return cls(**{key: data[key] for key in fields})


class TMvecModel(nn.Module):
    """TMVec2 transformer and projection head."""

    def __init__(self, config: TMvecConfig):
        super().__init__()
        self.config = config
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=config.d_model,
            nhead=config.nhead,
            dim_feedforward=config.dim_feedforward,
            dropout=config.dropout,
            activation=config.activation,
            batch_first=True,
            norm_first=False,
        )
        self.encoder = nn.TransformerEncoder(
            encoder_layer,
            num_layers=config.num_layers,
        )
        self.dropout = nn.Dropout(config.dropout)
        self.projection = nn.Sequential(
            nn.Linear(config.d_model, config.projection_hidden_dim),
            nn.ReLU(),
            nn.Dropout(config.dropout),
            nn.Linear(config.projection_hidden_dim, config.out_dim),
        )

    def forward(
        self,
        embeddings: torch.Tensor,
        padding_mask: torch.Tensor,
    ) -> torch.Tensor:
        """Map Lobster residue features to one TMVec2 vector per sequence."""
        hidden = self.encoder(
            embeddings,
            src_key_padding_mask=padding_mask,
        )
        lengths = (~padding_mask).sum(dim=1, keepdim=True).float().clamp(min=1e-9)
        mask = (~padding_mask).unsqueeze(-1).float()
        pooled = (hidden * mask).sum(dim=1) / lengths
        return self.projection(self.dropout(pooled))


class TMvecPredictor:
    """Create TMVec2 embeddings from local Lobster-24M and TMVec2 files."""

    def __init__(
        self,
        device: str = "cuda",
        model_root: Optional[Path] = None,
        require_gpu: bool = False,
        fail_on_unavailable: bool = False,
    ) -> None:
        self.require_gpu = require_gpu
        self.fail_on_unavailable = fail_on_unavailable
        if device not in {"cpu", "cuda"}:
            raise ValueError(f"Unsupported TMVec device: {device}")
        if device == "cuda" and not torch.cuda.is_available():
            if self._must_fail:
                raise RuntimeError(
                    "TMvecPredictor requested CUDA but no CUDA device is available."
                )
            logger.warning("CUDA is not available. TMVec2 will use the CPU.")
            device = "cpu"

        self.device = device
        self.model_root = Path(model_root).expanduser() if model_root else None
        self._lobster_model = None
        self._tokenizer = None
        self._tmvec_model = None
        self._config: Optional[TMvecConfig] = None
        self._initialized = False
        self._available: Optional[bool] = None
        self._batch_oom_fallbacks = 0
        self._per_seq_fallback_failures = 0

    @property
    def _must_fail(self) -> bool:
        return self.require_gpu or self.fail_on_unavailable

    @property
    def _lobster_path(self) -> Optional[Path]:
        return self.model_root / "lobster_24M" if self.model_root else None

    @property
    def _tmvec_path(self) -> Optional[Path]:
        return self.model_root / "tmvec-2" if self.model_root else None

    @property
    def available(self) -> bool:
        """Return true when dependencies and local model files are present."""
        if self._available is None:
            try:
                from lobster.model._mlm import LobsterPMLM  # noqa: F401

                required = (
                    self._lobster_path / "config.json" if self._lobster_path else None,
                    self._lobster_path / "pytorch_model.bin"
                    if self._lobster_path
                    else None,
                    self._lobster_path / "vocab.txt" if self._lobster_path else None,
                    self._lobster_path / "tokenizer_config.json"
                    if self._lobster_path
                    else None,
                    self._lobster_path / "special_tokens_map.json"
                    if self._lobster_path
                    else None,
                    self._tmvec_path / "params.json" if self._tmvec_path else None,
                    self._tmvec_path / "tmvec-2.ckpt" if self._tmvec_path else None,
                )
                self._available = all(
                    path is not None and path.is_file() for path in required
                )
            except ImportError as exc:
                logger.warning("TMVec2 dependency is not available: %s", exc)
                self._available = False
        return self._available

    def _lazy_init(self) -> None:
        if self._initialized:
            return
        with _MODEL_LOADING_LOCK:
            if self._initialized:
                return
            if not self.available:
                raise RuntimeError(
                    "TMVec2 model files or the Lobster runtime are not available. "
                    "Run ViroSync optional-resource setup for TMVec2."
                )
            self._load_tmvec2()
            self._load_lobster()
            self._initialized = True
            logger.info("TMVec2 loaded on %s", self.device)

    def _load_tmvec2(self) -> None:
        if self._tmvec_path is None:
            raise RuntimeError("TMVec2 model root is not configured")
        config = TMvecConfig.from_json(self._tmvec_path / "params.json")
        expected = TMvecConfig()
        if config != expected:
            raise RuntimeError(
                "Installed TMVec2 parameters do not match the supported architecture"
            )

        model = TMvecModel(config)
        checkpoint = torch.load(
            self._tmvec_path / "tmvec-2.ckpt",
            map_location="cpu",
            weights_only=True,
        )
        state_dict = checkpoint.get("state_dict", checkpoint)
        model.load_state_dict(state_dict, strict=True)
        self._tmvec_model = model.to(self.device).eval()
        self._config = config

    def _load_lobster(self) -> None:
        if self._lobster_path is None:
            raise RuntimeError("Lobster model root is not configured")
        from lobster.model._mlm import LobsterPMLM

        model = LobsterPMLM(str(self._lobster_path))
        if int(model.model.config.hidden_size) != TMvecConfig().d_model:
            raise RuntimeError(
                "Installed Lobster model does not produce 408-value residue features"
            )
        self._tokenizer = model.tokenizer
        self._lobster_model = model.to(self.device).eval()

    @staticmethod
    def _prepare_sequence(sequence: str) -> str:
        prepared = "".join(sequence.split()).upper().rstrip("*")
        if not prepared:
            raise ValueError("TMVec2 sequence is empty after terminal stop removal")
        if "*" in prepared:
            raise ValueError("TMVec2 sequence contains an internal stop symbol")
        return prepared

    def _embed_prepared_batch(self, sequences: list[str]) -> np.ndarray:
        encoded = self._tokenizer(
            sequences,
            padding=True,
            truncation=True,
            max_length=_MAX_TOKEN_LENGTH,
            return_tensors="pt",
        )
        input_ids = encoded["input_ids"].to(self.device)
        attention_mask = encoded["attention_mask"].to(self.device)
        with torch.inference_mode():
            outputs = self._lobster_model.model(
                input_ids=input_ids,
                attention_mask=attention_mask,
                output_hidden_states=True,
            )
            residue_features = outputs.hidden_states[-1]
            if residue_features.shape[-1] != TMvecConfig().d_model:
                raise RuntimeError(
                    "Lobster returned residue features with an invalid width: "
                    f"{residue_features.shape[-1]}"
                )
            embeddings = self._tmvec_model(
                residue_features,
                padding_mask=attention_mask.eq(0),
            )
        result = embeddings.detach().cpu().numpy().astype(np.float32, copy=False)
        norms = np.linalg.norm(result, axis=1)
        if not np.isfinite(result).all() or not np.all(np.isfinite(norms) & (norms > 1e-8)):
            raise RuntimeError("TMVec2 produced an invalid or zero embedding")
        return result

    def release(self) -> None:
        """Release model memory. The predictor can load the models again later."""
        self._lobster_model = None
        self._tmvec_model = None
        self._tokenizer = None
        self._config = None
        self._initialized = False
        try:
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        except Exception:
            pass

    def embed(self, sequence: str) -> np.ndarray:
        """Return one 512-value TMVec2 vector."""
        self._lazy_init()
        return self._embed_prepared_batch([self._prepare_sequence(sequence)])[0]

    def embed_batch(
        self,
        sequences: list[str],
        max_residues: int = 4000,
        max_batch: int = 100,
    ) -> np.ndarray:
        """Return TMVec2 vectors with a bounded token count per model call."""
        self._lazy_init()
        if max_residues < 1 or max_batch < 1:
            raise ValueError("max_residues and max_batch must be positive")
        if self._config is None:
            raise RuntimeError("TMVec2 is not initialized")

        prepared: list[tuple[int, str]] = []
        results = np.zeros((len(sequences), self._config.out_dim), dtype=np.float32)
        for index, sequence in enumerate(sequences):
            prepared.append((index, self._prepare_sequence(sequence)))
        prepared.sort(key=lambda item: len(item[1]), reverse=True)

        batches: list[list[tuple[int, str]]] = []
        current: list[tuple[int, str]] = []
        residue_count = 0
        for item in prepared:
            item_length = max(1, min(len(item[1]), _MAX_TOKEN_LENGTH))
            if current and (
                residue_count + item_length > max_residues
                or len(current) >= max_batch
            ):
                batches.append(current)
                current = []
                residue_count = 0
            current.append(item)
            residue_count += item_length
        if current:
            batches.append(current)

        for batch_index, batch in enumerate(batches):
            indices = [item[0] for item in batch]
            batch_sequences = [item[1] for item in batch]
            try:
                batch_embeddings = self._embed_prepared_batch(batch_sequences)
                results[indices] = batch_embeddings
            except Exception as exc:
                if "out of memory" in str(exc).lower():
                    self._batch_oom_fallbacks += 1
                    if torch.cuda.is_available():
                        torch.cuda.empty_cache()
                if self._must_fail:
                    raise
                logger.warning(
                    "TMVec2 batch %d failed; retrying one sequence at a time: %s",
                    batch_index,
                    exc,
                )
                for index, sequence in batch:
                    try:
                        results[index] = self._embed_prepared_batch([sequence])[0]
                    except Exception as sequence_exc:
                        self._per_seq_fallback_failures += 1
                        logger.warning(
                            "TMVec2 sequence %d failed: %s",
                            index,
                            sequence_exc,
                        )
        return results

    def search_database(
        self,
        query_sequence: str,
        database_embeddings: np.ndarray,
        database_ids: list[str],
        top_k: int = 10,
        min_tm: float = 0.3,
    ) -> list[tuple[str, float]]:
        """Return database IDs with the highest cosine similarity."""
        query = self.embed(query_sequence)
        if database_embeddings.ndim != 2 or database_embeddings.shape[1] != query.shape[0]:
            raise ValueError(
                "TMVec2 database embeddings must be a 2-D array with width "
                f"{query.shape[0]}"
            )
        if len(database_ids) != database_embeddings.shape[0]:
            raise ValueError("TMVec2 database ID and embedding row counts differ")
        query_normalized = query / np.linalg.norm(query)
        database_normalized = database_embeddings / (
            np.linalg.norm(database_embeddings, axis=1, keepdims=True) + 1e-8
        )
        similarities = database_normalized @ query_normalized
        hits: list[tuple[str, float]] = []
        for index in np.argsort(similarities)[::-1]:
            score = float(similarities[index])
            if score < min_tm or len(hits) >= top_k:
                break
            hits.append((database_ids[index], score))
        return hits


_live_predictors: weakref.WeakSet = weakref.WeakSet()


def get_tmvec_predictor(
    device: str = "cuda",
    model_root: Optional[Path] = None,
    require_gpu: bool = False,
    fail_on_unavailable: bool = False,
) -> TMvecPredictor:
    """Create a TMVec2 predictor from local model files."""
    try:
        predictor = TMvecPredictor(
            device=device,
            model_root=model_root,
            require_gpu=require_gpu,
            fail_on_unavailable=fail_on_unavailable,
        )
        if predictor.available:
            _live_predictors.add(predictor)
            return predictor
    except Exception as exc:
        logger.warning("TMVec2 is not available: %s", exc)
        if require_gpu or fail_on_unavailable:
            raise
    raise RuntimeError("No TMVec2 predictor is available")


def release_tmvec_predictor() -> None:
    """Release all predictor instances made by the factory."""
    for predictor in list(_live_predictors):
        try:
            predictor.release()
        except Exception as exc:
            logger.debug("TMVec2 predictor release failed: %s", exc)
    _live_predictors.clear()
