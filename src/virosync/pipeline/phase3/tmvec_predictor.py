"""
TMvec Structural Similarity Predictor.

This module provides GPU-accelerated TM-score prediction using the TMvec models.
TMvec predicts structural similarity (TM-score) directly from protein sequences
without requiring actual structure prediction.

Key features:
- Uses ProtT5-XL for per-residue embeddings
- TMvec transformer encoder learns structural representations
- Cosine similarity of embeddings ≈ TM-score
- Supports batch processing for efficient database searches

Model variants:
- tmvec_swiss_model: Base model (sequences ≤300 residues)
- tmvec_swiss_model_large: Large model (sequences ≤1000 residues)
- tm_vec_cath_model: CATH-trained base model
- tm_vec_cath_model_large: CATH-trained large model

References:
- Hamamsy et al., Nature Biotechnology 2023
  https://www.nature.com/articles/s41587-023-01917-2
- GitHub: https://github.com/tymor22/tm-vec
"""

# Set CUDA_VISIBLE_DEVICES before torch import if specified in environment
# This must happen before any torch imports to take effect
import os

if os.environ.get("CUDA_VISIBLE_DEVICES"):
    # Already set by parent process, ensure it persists
    pass
elif os.environ.get("VIROSYNC_GPU"):
    # Alternative: use VIROSYNC_GPU config to set device
    os.environ["CUDA_VISIBLE_DEVICES"] = os.environ["VIROSYNC_GPU"]

import logging
import weakref
import re
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import json
import hashlib

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

logger = logging.getLogger(__name__)

PROTTRANS_MODEL_ID = "Rostlab/prot_t5_xl_uniref50"
PROTTRANS_MODEL_REVISION = "973be27c52ee6474de9c945952a8008aeb2a1a73"
PROTTRANS_PROXY_MODEL_ID = "Rostlab/prot_t5_xl_half_uniref50-enc"
PROTTRANS_PROXY_MODEL_REVISION = "94a6abc029ae13029317b140b7424e012bf8dfbf"

# Global lock to prevent concurrent model loading (meta tensor issues)
_model_loading_lock = threading.Lock()
_PROTTRANS_EMBED_DIM = 1024


@dataclass
class TMvecConfig:
    """Configuration for TMvec model."""

    d_model: int = 1024
    nhead: int = 4
    num_layers: int = 2
    dim_feedforward: int = 2048
    out_dim: int = 512
    dropout: float = 0.1
    activation: str = "relu"
    max_length: int = 1000

    projection_hidden_dim: Optional[int] = None  # If set, use two-layer projection

    @classmethod
    def from_json(cls, path: Path) -> "TMvecConfig":
        """Load config from JSON file."""
        with open(path) as f:
            data = json.load(f)
        # Filter to known fields
        known_fields = {
            "d_model", "nhead", "num_layers", "dim_feedforward",
            "out_dim", "dropout", "activation", "max_length",
            "projection_hidden_dim"
        }
        filtered = {k: v for k, v in data.items() if k in known_fields}
        return cls(**filtered)


class TMvecModel(nn.Module):
    """
    TMvec Transformer Encoder for structural similarity.

    Architecture:
    1. TransformerEncoder processes per-residue embeddings
    2. Global average pooling over sequence length
    3. Dropout + projection to output dimension
    4. Cosine similarity of output vectors ≈ TM-score
    """

    def __init__(self, config: TMvecConfig):
        super().__init__()
        self.config = config

        # Transformer encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=config.d_model,
            nhead=config.nhead,
            dim_feedforward=config.dim_feedforward,
            dropout=config.dropout,
            activation=config.activation,
            batch_first=True,
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=config.num_layers)

        self.dropout = nn.Dropout(config.dropout)

        # Projection layer(s)
        if config.projection_hidden_dim:
            # Two-layer projection for checkpoints that use it.
            self.projection = nn.Sequential(
                nn.Linear(config.d_model, config.projection_hidden_dim),
                nn.Dropout(config.dropout),
                nn.GELU() if config.activation == "gelu" else nn.ReLU(),
                nn.Linear(config.projection_hidden_dim, config.out_dim),
            )
        else:
            # Single linear projection (original tm-vec)
            self.projection = nn.Linear(config.d_model, config.out_dim)

    def forward(
        self,
        x: torch.Tensor,
        padding_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Forward pass through TMvec encoder.

        Args:
            x: Input embeddings [batch, seq_len, d_model]
            padding_mask: Padding mask [batch, seq_len] (True = padded)

        Returns:
            Output embeddings [batch, out_dim]
        """
        # Transformer encoder
        x = self.encoder(x, src_key_padding_mask=padding_mask)

        # Global average pooling (excluding padding)
        if padding_mask is not None:
            mask = ~padding_mask  # True for real tokens
            lens = mask.sum(dim=1, keepdim=True).float().clamp(min=1)
            x = (x * mask.unsqueeze(-1).float()).sum(dim=1) / lens
        else:
            x = x.mean(dim=1)

        x = self.dropout(x)

        # Project to output dimension
        if isinstance(self.projection, nn.Sequential):
            x = self.projection(x)
        else:
            x = self.projection(x)

        return x


class TMvecPredictor:
    """
    GPU-accelerated TM-score prediction using TMvec.

    This class provides the proper TMvec implementation that:
    1. Uses ProtT5-XL for per-residue protein embeddings
    2. Processes through trained TMvec transformer encoder
    3. Computes cosine similarity as TM-score prediction

    For remote homologs (sequence similarity ≤10%) with high structural
    similarity (TM-score ≥0.6), TMvec predicts TM-scores within 0.026
    of TM-align computed values.

    Example:
        predictor = TMvecPredictor(device="cuda")
        tm_score = predictor.predict_tm_score(seq1, seq2)
    """

    # Model download locations
    FIGSHARE_MODELS = {
        "tmvec_swiss_model": "https://ndownloader.figshare.com/files/47502971",
        "tmvec_swiss_model_large": "https://ndownloader.figshare.com/files/49181530",
        "tm_vec_cath_model": "https://ndownloader.figshare.com/files/46296322",
        "tm_vec_cath_model_large": "https://ndownloader.figshare.com/files/49181533",
    }

    FIGSHARE_CONFIGS = {
        "tmvec_swiss_model": "https://ndownloader.figshare.com/files/47502968",
        "tmvec_swiss_model_large": "https://ndownloader.figshare.com/files/49181515",
        "tm_vec_cath_model": "https://ndownloader.figshare.com/files/46296310",
        "tm_vec_cath_model_large": "https://ndownloader.figshare.com/files/49181518",
    }

    FIGSHARE_MD5 = {
        "tm_vec_swiss_model_params.json": "ca3135729eae51d3fe7e462202f537e7",
        "tm_vec_swiss_model.ckpt": "316572493c242a6f6140c9ed87428ce1",
        "tm_vec_swiss_model_large_params.json": "fbb1f2288be74ad6c5ac1c05a19f876d",
        "tm_vec_swiss_model_large.ckpt": "69d8ef7a3286b8f6077fb89310ac19dd",
        "tm_vec_cath_model_params.json": "0fcbc713c2c55a49b7c825cb1b9c3ea5",
        "tm_vec_cath_model.ckpt": "e2ed6e09d6fd9b578018a8cd66bc07ff",
        "tm_vec_cath_model_large_params.json": "fbb1f2288be74ad6c5ac1c05a19f876d",
        "tm_vec_cath_model_large.ckpt": "150849899c58cd06cab03401c30cc589",
    }

    def __init__(
        self,
        device: str = "cuda",
        model_name: str = "tmvec_swiss_model_large",
        cache_dir: Optional[Path] = None,
        require_gpu: bool = False,
        fail_on_unavailable: bool = False,
    ):
        """
        Initialize TMvec predictor.

        Args:
            device: Torch device ("cuda" or "cpu")
            model_name: Which TMvec model to use
            cache_dir: Directory to cache downloaded models
            require_gpu: If True, raise RuntimeError when the requested
                device is "cuda" but no CUDA device is available, and
                disallow silent CPU fallback from meta-tensor / OOM paths.
            fail_on_unavailable: If True, propagate unforeseen initialization,
                embedding, and fallback failures after preflight enabled TMVec.
        """
        self.require_gpu = require_gpu
        self.fail_on_unavailable = fail_on_unavailable
        if device not in {"cpu", "cuda"}:
            raise ValueError(f"Unsupported TMVec device: {device}")
        if device == "cuda" and not torch.cuda.is_available():
            if self._must_fail:
                raise RuntimeError(
                    "TMvecPredictor requested CUDA but no CUDA device is available."
                )
            logger.warning(
                "TMvecPredictor: CUDA requested but not available; degrading to CPU. "
                "Set tmvec_require_gpu=True to fail fast instead of silently using CPU."
            )
            self.device = "cpu"
        else:
            self.device = device
        self.model_name = model_name
        self.cache_dir = cache_dir or Path.home() / ".cache" / "virosync" / "tmvec"
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Lazy-loaded models
        self._prottrans_model = None
        self._tokenizer = None
        self._tmvec_model = None
        self._config = None
        self._initialized = False
        self._available = None
        # Counters surfaced in logs so the pipeline can record GPU OOM
        # fallback severity instead of silently degrading scoring.
        self._batch_oom_fallbacks = 0
        self._per_seq_fallback_failures = 0

    @property
    def _must_fail(self) -> bool:
        return self.require_gpu or self.fail_on_unavailable

    @property
    def available(self) -> bool:
        """Check if TMvec dependencies are available."""
        if self._available is None:
            try:
                # Direct imports avoid transformers' lazy import issues in forked processes
                from transformers.models.t5.modeling_t5 import T5EncoderModel
                from transformers.models.t5.tokenization_t5 import T5Tokenizer
                import pytorch_lightning  # Required for loading checkpoints
                self._available = True
            except ImportError as e:
                logger.warning(f"TMvec not available: {e}")
                self._available = False
        return self._available

    @staticmethod
    def _expected_filenames(model_name: str) -> tuple[str, str]:
        """Return Figshare config/checkpoint filenames for a TM-Vec model key."""
        if model_name == "tmvec_swiss_model":
            stem = "tm_vec_swiss_model"
        elif model_name == "tmvec_swiss_model_large":
            stem = "tm_vec_swiss_model_large"
        elif model_name in {"tm_vec_cath_model", "tm_vec_cath_model_large"}:
            stem = model_name
        else:
            raise ValueError(f"Unknown TMVec model: {model_name}")
        return f"{stem}_params.json", f"{stem}.ckpt"

    @staticmethod
    def _md5(path: Path) -> str:
        digest = hashlib.md5(usedforsecurity=False)
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    def _download_model(
        self,
        url: str,
        filename: str,
        expected_md5: Optional[str] = None,
    ) -> Path:
        """Download model file from URL."""
        import urllib.request

        filepath = self.cache_dir / filename
        if filepath.exists():
            if expected_md5 and self._md5(filepath) != expected_md5:
                logger.warning("Cached TMVec file failed checksum; re-downloading %s", filename)
                filepath.unlink()
            else:
                return filepath

        logger.info(f"Downloading {filename}...")
        urllib.request.urlretrieve(url, filepath)

        if expected_md5:
            actual_md5 = self._md5(filepath)
            if actual_md5 != expected_md5:
                filepath.unlink(missing_ok=True)
                raise RuntimeError(
                    f"Checksum mismatch for {filename}: expected {expected_md5}, got {actual_md5}"
                )

            return filepath
        return filepath

    def _lazy_init(self) -> None:
        """Lazy initialization of models.

        Uses a global lock to prevent concurrent model loading which can
        cause meta tensor issues when multiple workers load simultaneously.
        """
        if self._initialized:
            return

        # Use lock to prevent concurrent model loading
        with _model_loading_lock:
            # Double-check after acquiring lock
            if self._initialized:
                return

            if not self.available:
                raise ImportError(
                    "TMvec dependencies not installed. Install with: "
                    "pip install transformers pytorch-lightning"
                )

            logger.info(f"Loading TMvec model: {self.model_name}...")

            # Load TMvec first so incompatible checkpoint/config pairs fail
            # before the much larger ProtT5 model is loaded.
            self._load_tmvec()

            # Load ProtT5 for embeddings
            self._load_prottrans()

            self._initialized = True
            logger.info(f"TMvec loaded on {self.device}")

    def _load_prottrans(self) -> None:
        """Load ProtT5-XL model for embeddings.

        Note: Uses low_cpu_mem_usage=False to avoid meta tensor issues when
        multiple workers try to load the model simultaneously.
        Uses direct submodule imports to avoid lazy import issues in forked processes.
        """
        # Direct imports avoid transformers' lazy import mechanism which breaks in forked processes
        from transformers.models.t5.modeling_t5 import T5EncoderModel
        from transformers.models.t5.tokenization_t5 import T5Tokenizer

        logger.info("Loading ProtT5-XL-UniRef50...")
        self._tokenizer = T5Tokenizer.from_pretrained(
            PROTTRANS_MODEL_ID,
            revision=PROTTRANS_MODEL_REVISION,
            do_lower_case=False,
        )

        # Load model with low_cpu_mem_usage=False to avoid meta tensor issues
        # This ensures the model is fully materialized before moving to GPU
        try:
            self._prottrans_model = T5EncoderModel.from_pretrained(
                PROTTRANS_MODEL_ID,
                revision=PROTTRANS_MODEL_REVISION,
                low_cpu_mem_usage=False,
            )
            self._prottrans_model = self._prottrans_model.to(self.device)
        except NotImplementedError as e:
            # Meta tensor error - fall back to CPU only (unless require_gpu)
            if "meta tensor" in str(e).lower() or "Cannot copy" in str(e):
                if self._must_fail:
                    logger.error(
                        "Meta tensor error loading ProtT5 to GPU and "
                        "tmvec_require_gpu=True; refusing silent CPU fallback."
                    )
                    raise
                logger.warning(
                    "Meta tensor error loading to GPU, falling back to CPU: %s", e
                )
                self.device = "cpu"
                self._prottrans_model = T5EncoderModel.from_pretrained(
                    PROTTRANS_MODEL_ID,
                    revision=PROTTRANS_MODEL_REVISION,
                    low_cpu_mem_usage=False,
                )
            else:
                raise
        self._prottrans_model.eval()

    def _load_tmvec(self) -> None:
        """Load trained original TM-Vec weights from Figshare."""
        if self.model_name in {"scikit-bio/tmvec-2", "TMVec-2", "tmvec-2"}:
            raise RuntimeError(
                "scikit-bio/tmvec-2 is not compatible with ViroSync's current "
                "TMVec database path: that checkpoint expects 408-dimensional "
                "input features, while ViroSync and the supported TM-Vec "
                "databases use 1024-dimensional ProtT5 residue embeddings. "
                "Use tmvec_swiss_model_large or tm_vec_cath_model_large."
            )

        if self.model_name not in self.FIGSHARE_MODELS:
            raise RuntimeError(
                f"Unsupported TMVec model '{self.model_name}'. Supported models: "
                f"{', '.join(sorted(self.FIGSHARE_MODELS))}"
            )

        config_filename, ckpt_filename = self._expected_filenames(self.model_name)
        config_path = self._download_model(
            self.FIGSHARE_CONFIGS[self.model_name],
            config_filename,
            self.FIGSHARE_MD5.get(config_filename),
        )
        ckpt_path = self._download_model(
            self.FIGSHARE_MODELS[self.model_name],
            ckpt_filename,
            self.FIGSHARE_MD5.get(ckpt_filename),
        )
        self._config = TMvecConfig.from_json(config_path)

        if self._config.d_model != _PROTTRANS_EMBED_DIM:
            raise RuntimeError(
                "TMVec model/input mismatch: trained checkpoint expects "
                f"{self._config.d_model}-dimensional residue embeddings, but "
                f"ViroSync's current TMVec path provides {_PROTTRANS_EMBED_DIM}-"
                "dimensional ProtT5 embeddings."
            )

        # Create model and load checkpoint
        self._tmvec_model = TMvecModel(self._config)

        ckpt = torch.load(ckpt_path, map_location="cpu")
        state_dict = ckpt.get("state_dict", ckpt)
        remapped = {
            key.replace("mlp.", "projection."): value
            for key, value in state_dict.items()
        }
        self._tmvec_model.load_state_dict(remapped, strict=True)

        self._tmvec_model = self._tmvec_model.to(self.device)
        self._tmvec_model.eval()
        logger.info("Loaded trained TMVec model %s", self.model_name)

    def _get_prottrans_embedding(self, sequence: str) -> torch.Tensor:
        """
        Get ProtT5 per-residue embeddings for a sequence.

        Args:
            sequence: Amino acid sequence

        Returns:
            Embeddings tensor [seq_len, 1024]

        Raises:
            RuntimeError: If meta tensor error occurs (GPU memory issue)
        """
        # Format sequence with spaces between amino acids
        seq_spaced = " ".join(list(sequence))
        seq_spaced = re.sub(r"[UZOB]", "X", seq_spaced)

        # Tokenize
        ids = self._tokenizer.batch_encode_plus(
            [seq_spaced], add_special_tokens=True, padding=True
        )
        input_ids = torch.tensor(ids["input_ids"]).to(self.device)
        attention_mask = torch.tensor(ids["attention_mask"]).to(self.device)

        try:
            with torch.no_grad():
                outputs = self._prottrans_model(
                    input_ids=input_ids,
                    attention_mask=attention_mask,
                )
        except NotImplementedError as e:
            # Meta tensor error - model not properly loaded
            if "meta tensor" in str(e).lower() or "Cannot copy" in str(e):
                raise RuntimeError(
                    f"Meta tensor error during inference: {e}. "
                    "This typically occurs when multiple workers try to load the model simultaneously. "
                    "Try running with fewer workers or use CPU-only mode."
                ) from e
            raise

        # Get embeddings excluding special tokens
        seq_len = (attention_mask[0] == 1).sum()
        embedding = outputs.last_hidden_state[0, :seq_len - 1]

        return embedding

    def _get_prottrans_embeddings_batch(
        self, sequences: list[str]
    ) -> list[torch.Tensor]:
        """
        Get ProtT5 per-residue embeddings for multiple sequences in ONE forward pass.

        Args:
            sequences: List of amino acid sequences (already cleaned/truncated)

        Returns:
            List of embedding tensors [seq_len_i, 1024] for each sequence
        """
        # Format sequences with spaces between amino acids
        seq_spaced = []
        for seq in sequences:
            s = " ".join(list(seq))
            s = re.sub(r"[UZOB]", "X", s)
            seq_spaced.append(s)

        # Tokenize ALL sequences together with padding
        ids = self._tokenizer.batch_encode_plus(
            seq_spaced,
            add_special_tokens=True,
            padding=True,
            return_tensors="pt",
        )
        input_ids = ids["input_ids"].to(self.device)
        attention_mask = ids["attention_mask"].to(self.device)

        try:
            with torch.no_grad():
                outputs = self._prottrans_model(
                    input_ids=input_ids,
                    attention_mask=attention_mask,
                )
        except NotImplementedError as e:
            if "meta tensor" in str(e).lower() or "Cannot copy" in str(e):
                raise RuntimeError(
                    f"Meta tensor error during batch inference: {e}"
                ) from e
            raise

        # Extract per-sequence embeddings (excluding special tokens)
        # Use .item() for safe tensor-to-int conversion (Codex review)
        embeddings = []
        for hidden_states, mask in zip(outputs.last_hidden_state, attention_mask):
            seq_len = int((mask == 1).sum().item())  # Safe conversion
            embedding = hidden_states[: seq_len - 1]
            embeddings.append(embedding)

        return embeddings

    def release(self) -> None:
        """Release GPU memory held by loaded models.

        Sets all model references to None and empties the CUDA cache.
        The predictor can be re-initialized later via ``_lazy_init()``.
        """
        self._prottrans_model = None
        self._tmvec_model = None
        self._tokenizer = None
        self._initialized = False

        try:
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        except Exception:
            pass

        logger.info("TMvecPredictor released GPU resources")

    def embed(self, sequence: str) -> np.ndarray:
        """
        Get TMvec embedding for a sequence.

        Args:
            sequence: Amino acid sequence

        Returns:
            Embedding vector (out_dim dimensional, typically 512)
        """
        self._lazy_init()

        # Clean sequence
        sequence = "".join(c for c in sequence.upper() if c in "ACDEFGHIKLMNPQRSTVWY")

        if len(sequence) < 5:
            raise ValueError(f"Sequence too short: {len(sequence)}")

        # Truncate if too long
        max_len = self._config.max_length
        if len(sequence) > max_len:
            logger.warning(f"Truncating sequence from {len(sequence)} to {max_len}")
            sequence = sequence[:max_len]

        # Get ProtT5 embeddings
        prottrans_emb = self._get_prottrans_embedding(sequence)

        # Process through TMvec
        with torch.no_grad():
            # Add batch dimension
            prottrans_emb = prottrans_emb.unsqueeze(0)
            tmvec_emb = self._tmvec_model(prottrans_emb)

        embedding = tmvec_emb.squeeze(0).cpu().numpy()
        norm = float(np.linalg.norm(embedding))
        if not np.isfinite(embedding).all() or norm <= 1e-8:
            raise RuntimeError("TMVec produced an invalid or zero embedding")
        return embedding

    def embed_batch(
        self,
        sequences: list[str],
        max_residues: int = 4000,
        max_batch: int = 100,
    ) -> np.ndarray:
        """
        Get TMvec embeddings for multiple sequences with token-budgeted batching.

        Uses ProtTrans-style batching that limits total residues per batch rather
        than a fixed sequence count. This prevents OOM with long sequences while
        maximizing throughput for short sequences.

        Sequences are sorted by length (descending) for optimal padding efficiency,
        then accumulated into batches based on total residue count.

        Args:
            sequences: List of amino acid sequences
            max_residues: Maximum total residues per batch (default 4000, per ProtTrans)
            max_batch: Maximum sequences per batch regardless of length (default 100)

        Returns:
            Embeddings array [n_sequences, out_dim]
        """
        self._lazy_init()

        max_len = self._config.max_length

        # Step 1: Clean all sequences and track original indices
        cleaned_with_idx: list[tuple[int, str, int]] = []  # (orig_idx, sequence, length)
        for i, seq in enumerate(sequences):
            clean = "".join(c for c in seq.upper() if c in "ACDEFGHIKLMNPQRSTVWY")
            if len(clean) < 5:
                logger.warning(f"Sequence {i} too short ({len(clean)} aa), will return zeros")
                continue
            if len(clean) > max_len:
                logger.debug(f"Truncating sequence {i} from {len(clean)} to {max_len}")
                clean = clean[:max_len]
            cleaned_with_idx.append((i, clean, len(clean)))

        # Initialize results array with zeros (float32 to match model output)
        results = np.zeros((len(sequences), self._config.out_dim), dtype=np.float32)

        if not cleaned_with_idx:
            return results

        # Step 2: Sort by length descending for optimal padding (ProtTrans strategy)
        cleaned_with_idx.sort(key=lambda x: x[2], reverse=True)

        # Step 3: Accumulate into token-budgeted batches
        batches: list[list[tuple[int, str]]] = []
        current_batch: list[tuple[int, str]] = []
        current_residues = 0

        for orig_idx, seq, seq_len in cleaned_with_idx:
            # Check if adding this sequence would exceed limits
            if current_batch and (
                current_residues + seq_len > max_residues
                or len(current_batch) >= max_batch
            ):
                # Process current batch first
                batches.append(current_batch)
                current_batch = []
                current_residues = 0

            current_batch.append((orig_idx, seq))
            current_residues += seq_len

        # Don't forget the last batch
        if current_batch:
            batches.append(current_batch)

        logger.info(
            f"Token-budgeted batching: {len(cleaned_with_idx)} sequences -> "
            f"{len(batches)} batches (max_residues={max_residues}, max_batch={max_batch})"
        )

        # Step 4: Process each batch
        for batch_idx, batch in enumerate(batches):
            orig_indices = [item[0] for item in batch]
            batch_seqs = [item[1] for item in batch]

            try:
                # Batched ProtT5 inference
                prottrans_embs = self._get_prottrans_embeddings_batch(batch_seqs)

                # Pad embeddings to same length for batched TMVec
                max_emb_len = max(e.shape[0] for e in prottrans_embs)
                padded = []
                masks = []

                for emb in prottrans_embs:
                    pad_len = max_emb_len - emb.shape[0]
                    if pad_len > 0:
                        emb_padded = F.pad(emb, (0, 0, 0, pad_len))
                        mask = torch.cat([
                            torch.zeros(emb.shape[0], dtype=torch.bool, device=self.device),
                            torch.ones(pad_len, dtype=torch.bool, device=self.device),
                        ])
                    else:
                        emb_padded = emb
                        mask = torch.zeros(emb.shape[0], dtype=torch.bool, device=self.device)
                    padded.append(emb_padded)
                    masks.append(mask)

                # Batched TMVec forward pass
                with torch.no_grad():
                    batch_tensor = torch.stack(padded)
                    mask_tensor = torch.stack(masks)
                    tmvec_embs = self._tmvec_model(batch_tensor, padding_mask=mask_tensor)

                # Assign results to original indices
                for orig_idx, emb in zip(orig_indices, tmvec_embs):
                    results[orig_idx] = emb.cpu().numpy()

            except RuntimeError as e:
                # Handle CUDA OOM specifically
                if "out of memory" in str(e).lower():
                    self._batch_oom_fallbacks += 1
                    logger.warning(
                        f"Batch {batch_idx} OOM ({len(batch_seqs)} seqs, "
                        f"~{sum(len(s) for s in batch_seqs)} residues), "
                        "falling back to per-sequence"
                    )
                    if torch.cuda.is_available():
                        torch.cuda.empty_cache()
                    if self._must_fail:
                        raise
                else:
                    logger.warning(f"Batch {batch_idx} failed: {e}, falling back to per-sequence")
                    if self._must_fail:
                        raise

                # Fallback to per-sequence
                for orig_idx, seq in batch:
                    try:
                        results[orig_idx] = self.embed(seq)
                    except Exception as seq_e:
                        self._per_seq_fallback_failures += 1
                        logger.warning(f"Per-sequence fallback failed for idx {orig_idx}: {seq_e}")

            except Exception as e:
                logger.warning(f"Batch {batch_idx} failed: {e}, falling back to per-sequence")
                if self._must_fail:
                    raise
                for orig_idx, seq in batch:
                    try:
                        results[orig_idx] = self.embed(seq)
                    except Exception as seq_e:
                        self._per_seq_fallback_failures += 1
                        logger.warning(f"Per-sequence fallback failed for idx {orig_idx}: {seq_e}")

        valid_indices = [orig_idx for orig_idx, _, _ in cleaned_with_idx]
        invalid_indices = [
            idx for idx in valid_indices
            if not np.isfinite(results[idx]).all()
            or float(np.linalg.norm(results[idx])) <= 1e-8
        ]
        if invalid_indices:
            msg = (
                "TMVec produced invalid or zero embeddings for "
                f"{len(invalid_indices)} valid sequence(s)"
            )
            if self._must_fail:
                raise RuntimeError(msg)
            logger.warning("%s; those sequences will have no TMVec hits", msg)
            for idx in invalid_indices:
                results[idx] = 0.0

        return results

    def predict_tm_score(self, seq1: str, seq2: str) -> float:
        """
        Predict TM-score between two sequences.

        Args:
            seq1: First amino acid sequence
            seq2: Second amino acid sequence

        Returns:
            Predicted TM-score (0-1)
        """
        self._lazy_init()

        # Clean sequences
        seq1 = "".join(c for c in seq1.upper() if c in "ACDEFGHIKLMNPQRSTVWY")
        seq2 = "".join(c for c in seq2.upper() if c in "ACDEFGHIKLMNPQRSTVWY")

        if len(seq1) < 5 or len(seq2) < 5:
            return 0.0

        # Get embeddings
        emb1 = self.embed(seq1)
        emb2 = self.embed(seq2)

        # Cosine similarity = TM-score prediction
        emb1_norm = emb1 / (np.linalg.norm(emb1) + 1e-8)
        emb2_norm = emb2 / (np.linalg.norm(emb2) + 1e-8)
        tm_score = np.dot(emb1_norm, emb2_norm)

        # TM-score is [0, 1]
        return float(np.clip(tm_score, 0.0, 1.0))

    def search_database(
        self,
        query_sequence: str,
        database_embeddings: np.ndarray,
        database_ids: list[str],
        top_k: int = 10,
        min_tm: float = 0.3,
    ) -> list[tuple[str, float]]:
        """
        Search for structurally similar proteins in a database.

        Args:
            query_sequence: Query amino acid sequence
            database_embeddings: Pre-computed embeddings [n, out_dim]
            database_ids: IDs corresponding to embeddings
            top_k: Number of top hits to return
            min_tm: Minimum TM-score threshold

        Returns:
            List of (id, tm_score) tuples sorted by score
        """
        self._lazy_init()

        query_emb = self.embed(query_sequence)

        query_norm_value = float(np.linalg.norm(query_emb))
        if not np.isfinite(query_emb).all() or query_norm_value <= 1e-8:
            msg = "TMVec produced an invalid or zero query embedding"
            if self._must_fail:
                raise RuntimeError(msg)
            logger.warning("%s; returning no database hits", msg)
            return []

        # Normalize for cosine similarity
        query_norm = query_emb / (np.linalg.norm(query_emb) + 1e-8)
        db_norms = database_embeddings / (
            np.linalg.norm(database_embeddings, axis=1, keepdims=True) + 1e-8
        )

        # Compute all similarities
        similarities = db_norms @ query_norm

        # Filter and sort
        hits = []
        for idx in np.argsort(similarities)[::-1]:
            score = float(similarities[idx])
            if score < min_tm:
                break
            if len(hits) >= top_k:
                break
            hits.append((database_ids[idx], score))

        return hits


class TMvecProxyPredictor:
    """
    Fallback predictor using ProtTrans embeddings directly.

    This is a simplified version that uses raw ProtT5 embeddings
    with mean pooling and cosine similarity. Less accurate than
    proper TMvec but works without the TMvec model weights.

    Use this when TMvec models are unavailable.
    """

    def __init__(self, device: str = "cuda"):
        self.device = device if torch.cuda.is_available() else "cpu"
        self._model = None
        self._tokenizer = None
        self._initialized = False

    @property
    def available(self) -> bool:
        try:
            # Direct imports avoid transformers' lazy import issues in forked processes
            from transformers.models.t5.modeling_t5 import T5EncoderModel
            return True
        except ImportError:
            return False

    def _lazy_init(self) -> None:
        """Lazy initialization of ProtTrans model.

        Uses a global lock to prevent concurrent model loading which can
        cause meta tensor issues when multiple workers load simultaneously.
        Uses direct submodule imports to avoid lazy import issues in forked processes.
        """
        if self._initialized:
            return

        # Use lock to prevent concurrent model loading
        with _model_loading_lock:
            # Double-check after acquiring lock
            if self._initialized:
                return

            # Direct imports avoid transformers' lazy import mechanism
            from transformers.models.t5.modeling_t5 import T5EncoderModel
            from transformers.models.t5.tokenization_t5 import T5Tokenizer

            logger.info("Loading ProtTrans T5 for similarity estimation...")
            self._tokenizer = T5Tokenizer.from_pretrained(
                PROTTRANS_PROXY_MODEL_ID,
                revision=PROTTRANS_PROXY_MODEL_REVISION,
                do_lower_case=False,
            )
            # Disable low_cpu_mem_usage to avoid meta tensor issues with .to()
            self._model = T5EncoderModel.from_pretrained(
                PROTTRANS_PROXY_MODEL_ID,
                revision=PROTTRANS_PROXY_MODEL_REVISION,
                low_cpu_mem_usage=False,
            )
            self._model = self._model.to(self.device)
            self._model.eval()
            self._initialized = True

    def embed(self, sequence: str) -> np.ndarray:
        """Get mean-pooled ProtTrans embedding."""
        self._lazy_init()

        sequence = "".join(c for c in sequence.upper() if c in "ACDEFGHIKLMNPQRSTVWY")
        seq_spaced = " ".join(list(sequence))
        seq_spaced = re.sub(r"[UZOB]", "X", seq_spaced)

        ids = self._tokenizer.batch_encode_plus(
            [seq_spaced], add_special_tokens=True, padding=True
        )
        input_ids = torch.tensor(ids["input_ids"]).to(self.device)
        attention_mask = torch.tensor(ids["attention_mask"]).to(self.device)

        with torch.no_grad():
            outputs = self._model(input_ids=input_ids, attention_mask=attention_mask)

        seq_len = (attention_mask[0] == 1).sum()
        embedding = outputs.last_hidden_state[0, :seq_len - 1]
        mean_emb = embedding.mean(dim=0)

        return mean_emb.cpu().numpy()

    def predict_tm_score(self, seq1: str, seq2: str) -> float:
        """Predict TM-score using cosine similarity of ProtTrans embeddings."""
        emb1 = self.embed(seq1)
        emb2 = self.embed(seq2)

        emb1_norm = emb1 / (np.linalg.norm(emb1) + 1e-8)
        emb2_norm = emb2 / (np.linalg.norm(emb2) + 1e-8)

        return float(max(0.0, np.dot(emb1_norm, emb2_norm)))

    def search_database(
        self,
        query_sequence: str,
        database_embeddings: np.ndarray,
        database_ids: list[str],
        top_k: int = 10,
        min_tm: float = 0.3,
    ) -> list[tuple[str, float]]:
        """
        Search for structurally similar proteins in a database.

        Uses cosine similarity of ProtTrans embeddings as a proxy for TM-score.
        Less accurate than proper TMvec but provides interface compatibility.

        Args:
            query_sequence: Query amino acid sequence
            database_embeddings: Pre-computed embeddings [n, out_dim]
            database_ids: IDs corresponding to embeddings
            top_k: Number of top hits to return
            min_tm: Minimum similarity threshold

        Returns:
            List of (id, similarity_score) tuples sorted by score
        """
        query_emb = self.embed(query_sequence)

        # Normalize for cosine similarity
        query_norm = query_emb / (np.linalg.norm(query_emb) + 1e-8)
        db_norms = database_embeddings / (
            np.linalg.norm(database_embeddings, axis=1, keepdims=True) + 1e-8
        )

        # Compute all similarities
        similarities = db_norms @ query_norm

        # Filter and sort
        hits = []
        for idx in np.argsort(similarities)[::-1]:
            score = float(similarities[idx])
            if score < min_tm:
                break
            if len(hits) >= top_k:
                break
            hits.append((database_ids[idx], score))

        return hits


# Every predictor handed out by get_tmvec_predictor() is tracked here so
# release_tmvec_predictor() can free GPU memory for the instances actually in
# use. Callers keep their own references (TMVecDatabaseSearch._predictor holds
# one), so a weak set frees nothing on its own and only observes.
_live_predictors: "weakref.WeakSet" = weakref.WeakSet()


def get_tmvec_predictor(
    device: str = "cuda",
    model_name: str = "tmvec_swiss_model_large",
    fallback_to_proxy: bool = False,
    require_gpu: bool = False,
    fail_on_unavailable: bool = False,
) -> TMvecPredictor | TMvecProxyPredictor:
    """
    Get the best available TMvec predictor.

    Args:
        device: Torch device
        model_name: TMvec model variant
        fallback_to_proxy: If True, fall back to proxy if TMvec unavailable
        require_gpu: If True, forward to TMvecPredictor so that a missing
            CUDA device or meta-tensor fallback raises instead of silently
            demoting to CPU.
        fail_on_unavailable: Propagate runtime failures after preflight enabled
            TMVec, independently of whether GPU use is mandatory.

    Returns:
        TMvec predictor instance
    """
    try:
        predictor = TMvecPredictor(
            device=device,
            model_name=model_name,
            require_gpu=require_gpu,
            fail_on_unavailable=fail_on_unavailable,
        )
        if predictor.available:
            _live_predictors.add(predictor)
            return predictor
    except Exception as e:
        logger.warning(f"TMvec not available: {e}")
        if require_gpu or fail_on_unavailable:
            raise

    if fallback_to_proxy:
        logger.info("Using ProtTrans proxy for structural similarity")
        proxy = TMvecProxyPredictor(device=device)
        _live_predictors.add(proxy)
        return proxy

    raise RuntimeError("No TMvec predictor available")


# Module-level predictor cache for reuse across calls within a worker.
_cached_predictor: TMvecPredictor | TMvecProxyPredictor | None = None


def get_cached_tmvec_predictor(
    device: str = "cuda",
    model_name: str = "tmvec_swiss_model_large",
    fallback_to_proxy: bool = False,
    require_gpu: bool = False,
    fail_on_unavailable: bool = False,
) -> TMvecPredictor | TMvecProxyPredictor:
    """Get or create a cached TMvec predictor.

    Identical to :func:`get_tmvec_predictor` but reuses a module-level
    singleton so that repeated calls within the same worker do not
    reload models.
    """
    global _cached_predictor
    if _cached_predictor is None:
        _cached_predictor = get_tmvec_predictor(
            device,
            model_name,
            fallback_to_proxy,
            require_gpu=require_gpu,
            fail_on_unavailable=fail_on_unavailable,
        )
    return _cached_predictor


def release_tmvec_predictor() -> None:
    """Release every live TMvec predictor and free GPU memory.

    Call this after Phase 3 completes to ensure the large ProtT5 + TMvec
    models (~45 GiB on GPU) are freed before other phases run on the same
    worker.

    This used to clear only the module cache written by
    get_cached_tmvec_predictor(), which nothing calls, so it freed nothing:
    the predictor actually in use is created through get_tmvec_predictor()
    and held by TMVecDatabaseSearch. Every predictor that factory hands out
    is now released here.
    """
    for predictor in list(_live_predictors):
        release = getattr(predictor, "release", None)
        if release is None:
            continue
        try:
            release()
        except Exception as exc:  # a stuck predictor must not block the rest
            logger.debug("TMVec predictor release failed: %s", exc)
    _live_predictors.clear()

    global _cached_predictor
    if _cached_predictor is not None:
        if hasattr(_cached_predictor, "release"):
            _cached_predictor.release()
        _cached_predictor = None
        logger.info("Cached TMvec predictor released")
