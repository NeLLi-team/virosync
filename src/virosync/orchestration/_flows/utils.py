"""Shared utilities for ViroSync orchestration functions.

Internal module - import from virosync.orchestration.flows instead.
"""

import inspect
from dataclasses import replace
from pathlib import Path
from typing import Any

from virosync.config import PipelineConfig
from virosync.config.pipeline_config import FIELD_SPECS


def _detect_explicit_overrides(
    signature: inspect.Signature,
    passed_kwargs: dict[str, Any],
    exclude_keys: set[str],
) -> dict[str, Any]:
    """Detect which parameters were explicitly overridden vs using defaults.

    Strategy: Compare passed values to signature defaults.
    If value != default, treat as explicit override.

    Args:
        signature: Function signature to inspect
        passed_kwargs: Dict of passed keyword arguments (from locals())
        exclude_keys: Keys to exclude from detection (e.g., 'config', 'genome_path')

    Returns:
        Dict containing only explicitly overridden parameters
    """
    defaults = {
        param.name: param.default
        for param in signature.parameters.values()
        if param.default is not inspect.Parameter.empty
    }

    explicit = {}
    for key, value in passed_kwargs.items():
        if key in exclude_keys:
            continue

        # If value differs from default, treat as explicit
        default_value = defaults.get(key, inspect.Parameter.empty)
        if default_value is inspect.Parameter.empty:
            # No default - always explicit
            explicit[key] = value
        elif value != default_value:
            # Value differs from default - explicit override
            explicit[key] = value
        # else: value == default, so config can override

    return explicit


def _merge_config_with_kwargs(
    config: PipelineConfig,
    explicit_kwargs: dict[str, Any],
    exclude_keys: set[str],
    *,
    include_none: bool = False,
) -> PipelineConfig:
    """Return nested pipeline config with flat public-flow values applied.

    The canonical field mapping keeps the public flat signature at the boundary
    while internal orchestration consumes the existing nested sections. Supplied
    configs ignore ``None`` overrides, matching historical precedence.

    Args:
        config: Base pipeline configuration.
        explicit_kwargs: Dict of explicitly provided kwargs
        exclude_keys: Keys to exclude from merging (e.g., 'config', 'genome_path')
        include_none: Apply ``None`` values when materializing no-config defaults.

    Returns:
        A copied configuration with updated nested sections.
    """
    updates_by_section: dict[str, dict[str, Any]] = {}
    for spec in FIELD_SPECS:
        if spec.flat in exclude_keys or spec.flat not in explicit_kwargs:
            continue
        value = explicit_kwargs[spec.flat]
        if value is None and not include_none:
            continue
        updates_by_section.setdefault(spec.section, {})[spec.field] = value

    merged = config
    for section_name, updates in updates_by_section.items():
        section = replace(getattr(merged, section_name), **updates)
        merged = replace(merged, **{section_name: section})
    return merged


def log_region_statistics(
    regions: list,
    logger,
    label: str = "Candidate regions",
) -> None:
    """Calculate and log summary statistics for candidate regions.

    Computes region count, length range/mean, marker counts,
    lineage support (NCLDV, MIRUS, PPV), MCP presence.

    Args:
        regions: List of candidate region objects (Anchor or similar)
        logger: Logger instance
        label: Descriptive label for log output
    """
    if not regions:
        logger.info(f"{label}: 0 regions")
        return

    # Compute statistics
    lengths = [r.length for r in regions]
    marker_counts = [r.marker_count for r in regions]
    ncldv_counts = [sum(1 for m in r.markers if getattr(m, "has_ncldv", 0)) for r in regions]
    mirus_counts = [sum(1 for m in r.markers if getattr(m, "has_mirus", 0)) for r in regions]
    # ``has_plv`` covers legacy PLV__ and current PPV__ labels.
    ppv_counts = [sum(1 for m in r.markers if getattr(m, "has_plv", 0)) for r in regions]
    mcp_counts = [sum(1 for m in r.markers if getattr(m, "is_mcp", False)) for r in regions]

    # Log summary
    logger.info(f"{label} summary:")
    logger.info(f"  Regions: {len(regions)}")
    logger.info(f"  Length: {min(lengths)}-{max(lengths)} bp (mean={sum(lengths) // len(lengths)})")
    logger.info(
        f"  Markers per region: {min(marker_counts)}-{max(marker_counts)} "
        f"(mean={sum(marker_counts) / len(marker_counts):.1f})"
    )
    logger.info("  Lineage support (top-10 validated):")
    ncldv_count = sum(1 for c in ncldv_counts if c > 0)
    logger.info(f"    NCLDV: {ncldv_count} regions ({ncldv_count * 100 / len(regions):.0f}%)")
    mirus_count = sum(1 for c in mirus_counts if c > 0)
    logger.info(f"    MIRUS: {mirus_count} regions ({mirus_count * 100 / len(regions):.0f}%)")
    ppv_count = sum(1 for c in ppv_counts if c > 0)
    logger.info(f"    PPV: {ppv_count} regions ({ppv_count * 100 / len(regions):.0f}%)")
    mcp_count = sum(1 for c in mcp_counts if c > 0)
    logger.info(f"  Regions with MCP markers: {mcp_count} ({mcp_count * 100 / len(regions):.0f}%)")


def build_marker_faa(
    marker_faa_dir: Path,
    output_path: Path,
    logger,
    rebuild: bool = False,
) -> Path | None:
    """Build combined marker.faa from individual marker FAA files."""
    marker_faa_dir = Path(marker_faa_dir)
    if not marker_faa_dir.exists():
        logger.error("Marker FAA directory not found: %s", marker_faa_dir)
        return None

    marker_files = sorted(marker_faa_dir.glob("*.faa"))
    if not marker_files:
        logger.error("No marker FAA files found in %s", marker_faa_dir)
        return None

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists() and not rebuild:
        output_mtime = output_path.stat().st_mtime
        newest_input = max((p.stat().st_mtime for p in marker_files), default=0)
        if output_mtime >= newest_input:
            logger.info("marker.faa is up to date; skipping rebuild")
            return output_path

    logger.info("Building marker.faa from %d marker FAA files", len(marker_files))
    with output_path.open("w") as out_handle:
        for faa in marker_files:
            with faa.open() as in_handle:
                last_char = "\n"
                for line in in_handle:
                    out_handle.write(line)
                    if line:
                        last_char = line[-1]
                if last_char != "\n":
                    out_handle.write("\n")
    logger.info("Wrote marker.faa: %s", output_path)
    return output_path


def ensure_combined_faa(
    faa_dir: Path,
    logger,
    marker_faa: Path | None = None,
    rebuild: bool = False,
    output_path: Path | None = None,
) -> Path | None:
    """Ensure combined.faa exists with all FAA files merged."""
    combined = Path(output_path) if output_path is not None else Path(faa_dir) / "combined.faa"
    combined.parent.mkdir(parents=True, exist_ok=True)
    if combined.exists() and not rebuild:
        combined_mtime = combined.stat().st_mtime
        inputs = [p for p in sorted(Path(faa_dir).glob("*.faa")) if p.name not in {"combined.faa", "marker.faa"}]
        marker_mtime = Path(marker_faa).stat().st_mtime if marker_faa and Path(marker_faa).exists() else 0
        newest_input = max([p.stat().st_mtime for p in inputs], default=0)
        if combined_mtime >= max(newest_input, marker_mtime):
            logger.info("combined.faa is up to date; skipping rebuild")
            return combined

    inputs = [p for p in sorted(Path(faa_dir).glob("*.faa")) if p.name not in {"combined.faa", "marker.faa"}]
    if not inputs:
        logger.error("No FAA files found to build combined.faa in %s", faa_dir)
        return None

    logger.info("Building combined.faa from %d FAA files", len(inputs))
    header_ids: set[str] = set()
    with combined.open("w") as out_handle:
        for faa in inputs:
            with faa.open() as in_handle:
                last_char = "\n"
                for line in in_handle:
                    out_handle.write(line)
                    if line:
                        last_char = line[-1]
                    if line.startswith(">"):
                        header_ids.add(line[1:].split()[0])
                if last_char != "\n":
                    out_handle.write("\n")

        if marker_faa and Path(marker_faa).exists():
            logger.info(
                "Appending marker.faa entries not present in combined.faa (%s)",
                Path(marker_faa).name,
            )
            with Path(marker_faa).open() as marker_handle:
                write_entry = True
                last_line = ""
                for line in marker_handle:
                    if line.startswith(">"):
                        entry_id = line[1:].split()[0]
                        write_entry = entry_id not in header_ids
                        if write_entry:
                            header_ids.add(entry_id)
                    if write_entry:
                        out_handle.write(line)
                    last_line = line
                if last_line and not last_line.endswith("\n"):
                    out_handle.write("\n")
    logger.info("Wrote combined.faa: %s", combined)
    return combined
