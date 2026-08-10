"""
Structural Homology Module for EVE Verification.

Primary pipeline path uses Boltz-2 structure prediction (optional),
followed by FoldSeek searches against viral protein structure databases.

This module provides a "tie-breaker" for ambiguous predictions where
sequence homology is insufficient but structural similarity can reveal
viral ancestry.
"""

import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from virosync.config import get_config
from virosync.utils.executables import resolve_boltz_executable
from virosync.utils.path_safety import (
    require_strict_child,
    safe_filename_component,
    safe_filename_components,
)

logger = logging.getLogger(__name__)


@dataclass
class StructurePrediction:
    """Result of structure prediction."""

    porf_id: str
    sequence: str
    pdb_string: str
    plddt_scores: np.ndarray
    mean_plddt: float
    ptm_score: float  # Predicted TM-score

    # Issue #3 fix: Explicit flags for prediction source
    is_predicted: bool = True  # True if real structure prediction, False if fallback/unavailable
    method: str = "predictor"  # "predictor", "esm2_fallback", "tmvec", "none"

    @property
    def is_confident(self) -> bool:
        """High-confidence structure (uses config.structural.plddt_high_confidence).

        Returns False if this is not a real structure prediction.
        """
        if not self.is_predicted:
            return False
        return self.mean_plddt > get_config().structural.plddt_high_confidence

    @property
    def is_very_confident(self) -> bool:
        """Very high-confidence structure (uses config.structural.plddt_very_high_confidence).

        Returns False if this is not a real structure prediction.
        """
        if not self.is_predicted:
            return False
        return self.mean_plddt > get_config().structural.plddt_very_high_confidence

    @property
    def has_structure(self) -> bool:
        """Whether this prediction has an actual 3D structure."""
        return bool(self.pdb_string)


@dataclass
class FoldSeekHit:
    """Result from FoldSeek structure search."""

    query_id: str
    target_id: str
    target_description: str
    tm_score: float  # Always available - primary metric

    # Issue #4 fix: These are Optional - only set when from real structural alignment
    evalue: Optional[float] = None  # Only from FoldSeek, not TMvec
    score: Optional[float] = None  # Bitscore from alignment
    qcov: Optional[float] = None  # Query coverage from alignment
    tcov: Optional[float] = None  # Target coverage from alignment
    lddt: Optional[float] = None  # Local distance difference test
    target_taxonomy: str = ""

    # Issue #4 fix: Flag indicating source of hit
    is_from_alignment: bool = True  # True if from FoldSeek, False if from TMvec


@dataclass
class StructuralHomologyResult:
    """Combined result of structure prediction and database search."""

    porf_id: str
    prediction: Optional[StructurePrediction]
    foldseek_hits: list[FoldSeekHit] = field(default_factory=list)

    # Classification
    has_viral_hit: bool = False
    best_viral_hit: Optional[FoldSeekHit] = None
    viral_taxonomy: str = ""

    # Confidence
    structural_evidence_score: float = 0.0

    @property
    def supports_viral_origin(self) -> bool:
        """Structure supports viral origin."""
        return self.has_viral_hit and self.structural_evidence_score > 0.5


class FoldSeekSearcher:
    """
    Structure-based homology search using FoldSeek.

    Searches predicted structures against viral protein structure databases.
    """

    def __init__(
        self,
        database_path: Optional[Path] = None,
        threads: int = 8,
        sensitivity: float = 7.5,
    ):
        """
        Initialize FoldSeek searcher.

        Args:
            database_path: Path to FoldSeek database (created if not exists)
            threads: Number of threads for search
            sensitivity: Search sensitivity (1-8, higher = more sensitive)
        """
        self.database_path = database_path
        self.threads = threads
        self.sensitivity = sensitivity
        self._check_foldseek()

    def _check_foldseek(self) -> None:
        """Check if FoldSeek is available."""
        try:
            result = subprocess.run(
                ["foldseek", "--version"],
                capture_output=True,
                text=True,
            )
            logger.info(f"FoldSeek version: {result.stdout.strip()}")
        except FileNotFoundError:
            logger.error("FoldSeek not found. Install via: pixi install")
            raise

    def search(
        self,
        query_pdb: Path,
        output_prefix: Optional[str] = None,
    ) -> list[FoldSeekHit]:
        """
        Search a query structure against the database.

        Args:
            query_pdb: Path to query PDB file
            output_prefix: Prefix for output files

        Returns:
            List of FoldSeekHit objects
        """
        if self.database_path is None:
            raise ValueError("No database path configured")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create query database
            query_db = tmpdir / "query"
            subprocess.run(
                ["foldseek", "createdb", str(query_pdb), str(query_db)],
                check=True,
                capture_output=True,
            )

            # Search
            result_db = tmpdir / "result"
            cmd = [
                "foldseek",
                "search",
                str(query_db),
                str(self.database_path),
                str(result_db),
                str(tmpdir / "tmp"),
                "--threads",
                str(self.threads),
                "-s",
                str(self.sensitivity),
                "-e",
                "10",  # E-value cutoff
                "-a",  # Required for convertalis structural alignment fields
            ]

            subprocess.run(cmd, check=True, capture_output=True)

            # Convert to tabular output
            result_tsv = tmpdir / "result.tsv"
            subprocess.run(
                [
                    "foldseek",
                    "convertalis",
                    str(query_db),
                    str(self.database_path),
                    str(result_db),
                    str(result_tsv),
                    "--format-output",
                    "query,target,evalue,bits,qcov,tcov,lddt,alntmscore",
                ],
                check=True,
                capture_output=True,
            )

            # Parse results
            return self._parse_results(result_tsv)

    def search_batch(
        self,
        pdb_files: list[Path],
        output_dir: Optional[Path] = None,
    ) -> dict[str, list[FoldSeekHit]]:
        """
        Search multiple structures against the database.

        Args:
            pdb_files: List of PDB file paths
            output_dir: Optional directory for results

        Returns:
            Dictionary mapping query IDs to hit lists
        """
        results = {}

        for pdb_file in pdb_files:
            query_id = pdb_file.stem
            try:
                hits = self.search(pdb_file)
                results[query_id] = hits
            except Exception as e:
                raise RuntimeError(
                    f"FoldSeek search failed for {query_id}; an empty hit list "
                    f"would be indistinguishable from a real absence of "
                    f"structural homology: {e}"
                ) from e

        return results

    def _parse_results(self, result_file: Path) -> list[FoldSeekHit]:
        """Parse FoldSeek tabular output."""
        hits = []

        if not result_file.exists():
            return hits

        with open(result_file) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 8:
                    try:
                        hit = FoldSeekHit(
                            query_id=parts[0],
                            target_id=parts[1],
                            target_description=parts[1],  # Would need DB info
                            tm_score=float(parts[7]),  # Required field first
                            # Real alignment statistics from FoldSeek
                            evalue=float(parts[2]),
                            score=float(parts[3]),
                            qcov=float(parts[4]),
                            tcov=float(parts[5]),
                            lddt=float(parts[6]),
                            is_from_alignment=True,  # Issue #4: Real FoldSeek hit
                        )
                        hits.append(hit)
                    except (ValueError, IndexError) as e:
                        logger.debug(f"Failed to parse hit: {e}")

        # Sort by TM-score descending
        hits.sort(key=lambda x: x.tm_score, reverse=True)

        return hits


class BoltzFoldSeekAnalyzer:
    """
    Structural analysis using Boltz-2 predictions + FoldSeek search.

    This is optional and should be used only for MCP candidates due to runtime cost.
    """

    def __init__(
        self,
        viral_db_path: Optional[Path],
        device: str = "cuda",
        threads: int = 8,
        use_msa_server: bool = False,
        min_seq_len: int = 100,
        max_seq_len: int = 1000,
        no_kernels: bool = True,
    ) -> None:
        self.viral_db_path = viral_db_path
        self.device = device
        self.threads = threads
        self.use_msa_server = use_msa_server
        self.min_seq_len = min_seq_len
        self.max_seq_len = max_seq_len
        self.no_kernels = no_kernels
        self._searcher: Optional[FoldSeekSearcher] = None
        self._boltz_executable: Optional[str] = None

    def available(self) -> bool:
        if self.viral_db_path is None:
            logger.warning("Boltz/FoldSeek disabled: viral_structure_db not set")
            return False
        if not self.use_msa_server:
            logger.warning(
                "Boltz/FoldSeek disabled: ViroSync's Boltz integration writes "
                "single-sequence YAML without local MSA files. Set "
                "boltz_use_msa_server=true to let Boltz generate MSAs."
            )
            return False
        self._boltz_executable = resolve_boltz_executable()
        if self._boltz_executable is None:
            logger.warning("Boltz not found (PATH or scripts/boltz); skipping structural analysis")
            return False
        if shutil.which("foldseek") is None:
            logger.warning("FoldSeek not found on PATH; skipping structural analysis")
            return False
        return True

    @property
    def searcher(self) -> FoldSeekSearcher:
        if self._searcher is None:
            self._searcher = FoldSeekSearcher(
                database_path=self.viral_db_path,
                threads=self.threads,
            )
        return self._searcher

    @staticmethod
    def _sanitize_id(porf_id: str) -> str:
        return safe_filename_component(porf_id)

    def _write_boltz_yaml(self, protein_id: str, sequence: str, yaml_dir: Path) -> Path:
        yaml_path = yaml_dir / f"{protein_id}.yaml"
        require_strict_child(yaml_dir, yaml_path)
        yaml_content = (
            "version: 1\n"
            "sequences:\n"
            "  - protein:\n"
            "      id: A\n"
            f"      sequence: {sequence}\n"
        )
        yaml_path.write_text(yaml_content)
        return yaml_path

    def _run_boltz(self, yaml_dir: Path, output_dir: Path) -> bool:
        boltz_executable = self._boltz_executable or resolve_boltz_executable()
        if boltz_executable is None:
            logger.warning("Boltz executable unavailable during run")
            return False

        cmd = [
            boltz_executable,
            "predict",
            str(yaml_dir),
            "--out_dir",
            str(output_dir),
        ]
        if self.use_msa_server:
            cmd.append("--use_msa_server")
        if self.no_kernels:
            cmd.append("--no_kernels")
        if self.device == "cpu":
            cmd.extend(["--accelerator", "cpu"])
        try:
            proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
        except Exception as exc:
            logger.warning("Boltz prediction failed: %s", exc)
            return False

        combined_output = "\n".join(part for part in (proc.stdout, proc.stderr) if part)
        if proc.returncode != 0:
            logger.warning("Boltz prediction failed with exit code %s: %s", proc.returncode, combined_output[-2000:])
            return False
        if "Failed to process" in combined_output or "Missing MSA" in combined_output:
            logger.warning("Boltz prediction reported failed inputs: %s", combined_output[-2000:])
            return False
        return True

    def _is_viral_target(self, target_id: str) -> bool:
        db_config = get_config().database
        if any(target_id.startswith(prefix) for prefix in db_config.viral_prefixes):
            return True
        viral_keywords = [
            "virus",
            "viral",
            "phage",
            "ncldv",
            "mimivirus",
            "pandoravirus",
            "megavirus",
        ]
        target_lower = target_id.lower()
        return any(kw in target_lower for kw in viral_keywords)

    def _calculate_evidence_score(self, result: StructuralHomologyResult) -> float:
        score = 0.0
        if result.has_viral_hit and result.best_viral_hit:
            hit = result.best_viral_hit
            tm_contribution = min(1.0, hit.tm_score) * 0.4
            if hit.evalue and hit.evalue > 0:
                evalue_contribution = min(1.0, -np.log10(hit.evalue) / 20) * 0.3
            else:
                evalue_contribution = 0.3
            score += tm_contribution + evalue_contribution
        return score

    def analyze_batch(
        self,
        porfs: list[tuple[str, str]],
        work_dir: Path,
    ) -> list[StructuralHomologyResult]:
        if not porfs or not self.available():
            return []

        eligible_porfs = [
            (porf_id, sequence)
            for porf_id, sequence in porfs
            if self.min_seq_len <= len(sequence) <= self.max_seq_len
        ]
        filename_components = safe_filename_components(
            (porf_id for porf_id, _sequence in eligible_porfs),
            label="protein ID",
        )

        work_dir.mkdir(parents=True, exist_ok=True)
        yaml_dir = work_dir / "boltz_yaml"
        pred_dir = work_dir / "boltz_predictions"
        require_strict_child(work_dir, yaml_dir)
        require_strict_child(work_dir, pred_dir)
        yaml_dir.mkdir(exist_ok=True)
        pred_dir.mkdir(exist_ok=True)

        id_map: dict[str, str] = {}
        sequences: list[tuple[str, str]] = []
        for porf_id, sequence in eligible_porfs:
            safe_id = filename_components[porf_id]
            id_map[safe_id] = porf_id
            sequences.append((safe_id, sequence))
            self._write_boltz_yaml(safe_id, sequence, yaml_dir)

        if not sequences:
            return []

        if not self._run_boltz(yaml_dir, pred_dir):
            return []

        structure_files = list(pred_dir.rglob("*.cif")) + list(pred_dir.rglob("*.pdb"))
        if not structure_files:
            logger.warning("Boltz produced no structures in %s", pred_dir)
            return []

        hit_map = self.searcher.search_batch(structure_files)
        results: list[StructuralHomologyResult] = []
        for query_id, hits in hit_map.items():
            porf_id = id_map.get(query_id, query_id)
            result = StructuralHomologyResult(porf_id=porf_id, prediction=None)
            result.foldseek_hits = hits or []
            viral_hits = [h for h in result.foldseek_hits if self._is_viral_target(h.target_id)]
            if viral_hits:
                best = max(viral_hits, key=lambda h: h.tm_score)
                result.has_viral_hit = True
                result.best_viral_hit = best
            result.structural_evidence_score = self._calculate_evidence_score(result)
            results.append(result)

        return results
