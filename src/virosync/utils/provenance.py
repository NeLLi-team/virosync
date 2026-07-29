"""Runtime provenance capture for reproducibility.

This module captures tool versions, database checksums, configuration parameters,
and execution metadata to enable exact reproduction of ViroSync results.
"""

from __future__ import annotations

import hashlib
import json
import logging
import subprocess
import sys
from dataclasses import fields, is_dataclass
from enum import Enum
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def _json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Enum):
        return value.value
    if is_dataclass(value):
        return {
            item.name: _json_safe(getattr(value, item.name))
            for item in fields(value)
        }
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def capture_tool_versions() -> dict[str, str]:
    """
    Capture versions of all tools used in the pipeline.

    Returns:
        Dictionary mapping tool names to version strings.
    """
    versions = {}

    # Python version
    versions["python"] = sys.version.split()[0]

    # ViroSync package version
    try:
        import virosync
        versions["virosync"] = getattr(virosync, "__version__", "unknown")
    except Exception:
        versions["virosync"] = "unknown"

    # PyTorch and PyTorch Geometric
    try:
        import torch
        versions["pytorch"] = torch.__version__
    except Exception:
        versions["pytorch"] = "not installed"

    try:
        import torch_geometric
        versions["pytorch_geometric"] = torch_geometric.__version__
    except Exception:
        versions["pytorch_geometric"] = "not installed"

    # pyhmmer
    try:
        import pyhmmer
        versions["pyhmmer"] = pyhmmer.__version__
    except Exception:
        versions["pyhmmer"] = "not installed"

    # External tools (via subprocess)
    tools = {
        "diamond": ["diamond", "version"],
        "prodigal-gv": ["prodigal-gv", "-v"],
        "foldseek": ["foldseek", "version"],
    }

    for tool, cmd in tools.items():
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=5,
                check=False,
            )
            # Parse first line of output (most tools print version there)
            if result.stdout:
                versions[tool] = result.stdout.strip().split('\n')[0]
            elif result.stderr:
                versions[tool] = result.stderr.strip().split('\n')[0]
            else:
                versions[tool] = "version unavailable"
        except FileNotFoundError:
            versions[tool] = "not found"
        except subprocess.TimeoutExpired:
            versions[tool] = "timeout"
        except Exception as e:
            versions[tool] = f"error: {e}"

    return versions


def compute_file_checksum(file_path: Path, algorithm: str = "sha256") -> str:
    """
    Compute checksum of a file.

    Args:
        file_path: Path to file.
        algorithm: Hash algorithm (md5, sha256).

    Returns:
        Hexadecimal checksum string.
    """
    if not file_path.exists():
        return "file not found"

    try:
        hasher = hashlib.new(algorithm)
        with open(file_path, 'rb') as f:
            # Read in chunks to handle large files
            for chunk in iter(lambda: f.read(65536), b''):
                hasher.update(chunk)
        return hasher.hexdigest()
    except Exception as e:
        logger.warning(f"Failed to compute checksum for {file_path}: {e}")
        return f"error: {e}"


def capture_database_info(config: dict[str, Any]) -> dict[str, dict[str, str]]:
    """
    Capture database paths and checksums.

    Args:
        config: Pipeline configuration dictionary.

    Returns:
        Dictionary mapping database names to path and checksum info.
    """
    databases = {}

    # HMM database
    hmm_db_path = config.get("hmm_database")
    if hmm_db_path:
        hmm_path = Path(hmm_db_path)
        databases["hmm_database"] = {
            "path": str(hmm_path.absolute()),
            "exists": hmm_path.exists(),
            "size_bytes": hmm_path.stat().st_size if hmm_path.exists() else 0,
            # HMM files are large, skip checksum for performance
            "checksum": "not computed (large file)",
        }

    # Marker Diamond DB
    marker_db_path = config.get("marker_db")
    if marker_db_path:
        marker_path = Path(marker_db_path)
        databases["marker_db"] = {
            "path": str(marker_path.absolute()),
            "exists": marker_path.exists(),
            "size_bytes": marker_path.stat().st_size if marker_path.exists() else 0,
            "checksum": "not computed (binary DB)",
        }

    # Gene taxonomy Diamond DB
    gene_tax_db_path = config.get("gene_taxonomy_faa_db")
    if gene_tax_db_path:
        gene_tax_path = Path(gene_tax_db_path)
        databases["gene_taxonomy_db"] = {
            "path": str(gene_tax_path.absolute()),
            "exists": gene_tax_path.exists(),
            "size_bytes": gene_tax_path.stat().st_size if gene_tax_path.exists() else 0,
            "checksum": "not computed (large binary DB)",
        }

    # Taxonomy labels
    tax_labels_path = config.get("taxonomy_labels")
    if tax_labels_path:
        tax_path = Path(tax_labels_path)
        databases["taxonomy_labels"] = {
            "path": str(tax_path.absolute()),
            "exists": tax_path.exists(),
            "size_bytes": tax_path.stat().st_size if tax_path.exists() else 0,
            # Taxonomy labels are small TSV, compute checksum
            "checksum": compute_file_checksum(tax_path) if tax_path.exists() else "N/A",
        }

    return databases


def write_provenance(
    output_dir: Path,
    config: dict[str, Any],
    input_genome: Path | None = None,
) -> None:
    """
    Write complete provenance information to run directory.

    Args:
        output_dir: Output directory for provenance.json.
        config: Pipeline configuration dictionary.
        input_genome: Optional path to input genome FASTA.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Capture once and reuse: version probing spawns external subprocesses.
    tool_versions = capture_tool_versions()
    masking_status = None
    masking_status_path = config.get("masking_status_path")
    if masking_status_path:
        from virosync.config import MaskingConfig
        from virosync.pipeline.phase0.masking import load_masking_result

        expected_masking = config.get("masking")
        if isinstance(expected_masking, dict):
            expected_masking = MaskingConfig(**expected_masking)
        masking_result = load_masking_result(
            Path(masking_status_path),
            expected_config=expected_masking,
        )
        expected_status_sha256 = config.get("masking_status_sha256")
        if (
            expected_status_sha256 is not None
            and masking_result.status_sha256 != expected_status_sha256
        ):
            raise ValueError("provenance masking status SHA256 mismatch")
        status_payload = json.loads(Path(masking_status_path).read_text(encoding="utf-8"))
        masking_status = {
            "path": str(Path(masking_status_path)),
            "sha256": masking_result.status_sha256,
            "result_fingerprint": status_payload["result_fingerprint"],
            "status": status_payload["status"],
            "benchmark_eligible": status_payload["benchmark_eligible"],
            "backend_versions": status_payload["backend_versions"],
            "requested_backend": status_payload["requested_backend"],
            "effective_backend": status_payload["effective_backend"],
            "failure_policy": status_payload["failure_policy"],
        }
    provenance = {
        "virosync_version": {
            "package": "virosync",
            "version": tool_versions.get("virosync", "unknown"),
        },
        "tool_versions": tool_versions,
        "databases": capture_database_info(config),
        "input_genome": {
            "path": str(input_genome.absolute()) if input_genome else "N/A",
            "checksum": compute_file_checksum(input_genome) if input_genome and input_genome.exists() else "N/A",
        },
        "configuration": {
            # Store subset of config to avoid leaking sensitive paths
            "use_tmvec_database": config.get("use_tmvec_database", False),
            "use_interproscan": config.get("use_interproscan", False),
            "masking": _json_safe(config.get("masking")),
            "marker_validation_top_k": config.get("marker_validation_top_k", 10),
        },
        "masking_status": masking_status,
    }

    output_file = output_dir / "provenance.json"
    try:
        with open(output_file, 'w') as f:
            json.dump(provenance, f, indent=2)
        logger.info(f"Provenance information written to {output_file}")
    except Exception as e:
        logger.warning(f"Failed to write provenance file: {e}")


