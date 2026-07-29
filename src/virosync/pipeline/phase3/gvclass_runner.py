"""GVClass batch runner for EVE classification."""

import csv
import logging
import subprocess
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def run_gvclass_batch(
    eve_fasta_dir: Path,
    output_dir: Path,
    gvclass_path: Path,
    threads: int = 8,
    gvclass_db: Optional[Path] = None,
) -> Optional[Path]:
    """
    Run GVClass on EVE nucleotide sequences.

    Args:
        eve_fasta_dir: Directory with EVE .fna files (phase3_synthesis/gvclass_input/nucleotide/)
        output_dir: Output directory for GVClass results
        gvclass_path: Path to the GVClass installation directory
        threads: Number of threads
        gvclass_db: Optional path to GVClass database directory

    Returns:
        Path to summary TSV or None if failed
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check if there are any input files
    fasta_files = list(eve_fasta_dir.glob("*.fna"))
    if not fasta_files:
        logger.warning(f"No FASTA files found in {eve_fasta_dir}")
        return None

    cmd = [
        str(gvclass_path / "gvclass"),
        str(eve_fasta_dir),
        "-o", str(output_dir),
        "-t", str(threads),
        "--mode-fast",
    ]

    # Add database path if provided
    if gvclass_db is not None:
        cmd.extend(["-d", str(gvclass_db)])

    logger.info(f"Running GVClass on {len(fasta_files)} EVE sequences")
    logger.debug(f"GVClass command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            check=True,
            timeout=1800,  # 30 minute timeout
            capture_output=True,
            text=True,
        )
        logger.debug(f"GVClass stdout: {result.stdout[:500] if result.stdout else ''}")

        # GVClass outputs: gvclass_summary.tsv or *.summary.tab
        for pattern in ["gvclass_summary.tsv", "*.summary.tab", "*_summary.tsv"]:
            matches = list(output_dir.glob(pattern))
            if matches:
                logger.info(f"GVClass completed: {matches[0]}")
                return matches[0]

        logger.warning(f"GVClass completed but no summary file found in {output_dir}")
        return None

    except subprocess.TimeoutExpired:
        logger.error("GVClass timed out after 30 minutes")
        return None
    except subprocess.CalledProcessError as e:
        logger.error(f"GVClass failed with exit code {e.returncode}")
        logger.error(f"GVClass stderr: {e.stderr[:500] if e.stderr else ''}")
        return None
    except FileNotFoundError:
        logger.error(f"GVClass executable not found at {gvclass_path / 'gvclass'}")
        return None
    except Exception as e:
        logger.error(f"GVClass failed with unexpected error: {e}")
        return None


def _is_header_line(parts: list[str]) -> bool:
    """Check if a line is a header by looking for known column names."""
    known_headers = {"file", "genome", "domain", "classification", "gvog", "mcp", "mirus"}
    for part in parts:
        if part.lower() in known_headers or any(h in part.lower() for h in known_headers):
            return True
    return False


def _strip_fna_suffix(filename: str) -> str:
    """Strip only .fna suffix from filename, preserving dots in ID."""
    if filename.endswith(".fna"):
        return filename[:-4]
    return filename


def load_gvclass_id_map(manifest_path: Path) -> dict[str, str]:
    """Load encoded nucleotide FASTA stems mapped to raw EVE identifiers."""
    if not manifest_path.exists():
        return {}

    id_map: dict[str, str] = {}
    with manifest_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            raw_id = row.get("eve_id", "")
            relative_path = row.get("nucleotide_fasta", "")
            if not raw_id or not relative_path:
                continue
            encoded_stem = _strip_fna_suffix(Path(relative_path).name)
            prior = id_map.get(encoded_stem)
            if prior is not None and prior != raw_id:
                raise ValueError(
                    f"GVClass manifest ID collision for {encoded_stem!r}: "
                    f"{prior!r} and {raw_id!r}"
                )
            id_map[encoded_stem] = raw_id
    return id_map


def parse_gvclass_results(
    summary_path: Path,
    id_map: Optional[dict[str, str]] = None,
) -> dict[str, dict]:
    """
    Parse GVClass summary output.

    Args:
        summary_path: Path to GVClass summary TSV
        id_map: Optional encoded input stem to raw EVE ID mapping

    Returns:
        Dictionary mapping eve_id to classification results:
        {eve_id: {domain, gvog_count, mcp_count, mirus_count}}
    """
    results = {}

    if not summary_path.exists():
        logger.warning(f"GVClass summary file not found: {summary_path}")
        return results

    with open(summary_path) as f:
        header = None
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if header is None and _is_header_line(parts):
                # Detect header by known column names
                header = parts
                continue

            if len(parts) < 2:
                continue

            # If we haven't found a header yet, skip data lines
            if header is None:
                logger.warning("GVClass output missing header, skipping")
                continue

            # First column is typically the input file name or path. Normalize to
            # its encoded stem, then recover the raw biological ID when a manifest
            # mapping is available. Without a manifest, preserve legacy behavior.
            legacy_id = _strip_fna_suffix(parts[0])
            encoded_stem = _strip_fna_suffix(Path(parts[0]).name)
            eve_id = (
                id_map.get(encoded_stem, legacy_id)
                if id_map is not None
                else legacy_id
            )

            # Build result dict based on available columns
            result = {
                "domain": "",
                "gvog_count": 0,
                "mcp_count": 0,
                "mirus_count": 0,
            }

            # Parse based on header if available
            for i, col in enumerate(header):
                if i >= len(parts):
                    break
                col_lower = col.lower()
                val = parts[i]

                if "domain" in col_lower or "classification" in col_lower:
                    result["domain"] = val
                elif "gvog" in col_lower and "count" in col_lower:
                    try:
                        result["gvog_count"] = int(val) if val.isdigit() else 0
                    except ValueError:
                        pass
                elif "mcp" in col_lower and "count" in col_lower:
                    try:
                        result["mcp_count"] = int(val) if val.isdigit() else 0
                    except ValueError:
                        pass
                elif "mirus" in col_lower and "count" in col_lower:
                    try:
                        result["mirus_count"] = int(val) if val.isdigit() else 0
                    except ValueError:
                        pass

            results[eve_id] = result

    logger.info(f"Parsed GVClass results for {len(results)} EVEs")
    return results


def write_gvclass_results_tsv(results: dict[str, dict], output_path: Path) -> Path:
    """
    Write parsed GVClass results to a simplified TSV file.

    Args:
        results: Dictionary from parse_gvclass_results()
        output_path: Output TSV path

    Returns:
        Path to written file
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["eve_id", "gvclass_domain", "gvog_count", "mcp_count", "mirus_count"]
        )
        for eve_id, data in sorted(results.items()):
            writer.writerow(
                [
                    eve_id,
                    data["domain"],
                    data["gvog_count"],
                    data["mcp_count"],
                    data["mirus_count"],
                ]
            )

    logger.info(f"Wrote GVClass results TSV: {output_path}")
    return output_path
