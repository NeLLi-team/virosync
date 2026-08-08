"""Strict, side-effect-free decoding for the public ViroSync configuration."""

from __future__ import annotations

import copy
import hashlib
import json
import re
from dataclasses import dataclass, fields, is_dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Optional

import yaml
from yaml.nodes import MappingNode

from .pipeline_config import (
    ConfigError,
    PipelineConfig,
    _decode_dataclass,
    _unknown_key,
)


class _UniqueKeyLoader(yaml.SafeLoader):
    """Safe YAML loader that rejects duplicate mapping keys at every depth."""


def _construct_unique_mapping(loader, node, deep=False):
    if not isinstance(node, MappingNode):
        raise ConfigError("Expected a YAML mapping")
    seen: set[Any] = set()
    for key_node, _ in node.value:
        key = loader.construct_object(key_node, deep=deep)
        try:
            duplicate = key in seen
        except TypeError as exc:
            raise ConfigError(
                f"YAML mapping key on line {key_node.start_mark.line + 1} is not scalar"
            ) from exc
        if duplicate:
            raise ConfigError(
                f"Duplicate YAML key {key!r} on line {key_node.start_mark.line + 1}"
            )
        seen.add(key)
    return yaml.SafeLoader.construct_mapping(loader, node, deep=deep)


_UniqueKeyLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
    _construct_unique_mapping,
)


@dataclass
class OrchestrationConfig:
    """Process-level settings that do not belong to a single genome flow."""

    database_root: Optional[Path] = None
    core_resources_url: Optional[str] = None
    core_resources_version: Optional[str] = None
    core_resources_sha256: Optional[str] = None
    core_resources_manifest_sha256: Optional[str] = None
    tmvec_resources_url: Optional[str] = None
    interproscan_resources_url: Optional[str] = None
    max_concurrent_genomes: int = 4
    retries: int = 1
    retry_delay_seconds: int = 60
    gpu_id: Optional[int] = None

    def validate_semantics(self) -> list[str]:
        errors = []
        identity = (
            self.core_resources_version,
            self.core_resources_sha256,
            self.core_resources_manifest_sha256,
        )
        if any(value is not None for value in identity) and not all(
            value is not None for value in identity
        ):
            errors.append(
                "orchestration core resource version, archive SHA-256, and manifest "
                "SHA-256 must be configured together"
            )
        if self.core_resources_url is not None and not all(
            value is not None for value in identity
        ):
            errors.append(
                "orchestration.core_resources_url requires core_resources_version, "
                "core_resources_sha256, and core_resources_manifest_sha256"
            )
        if self.core_resources_url is None and any(value is not None for value in identity):
            errors.append(
                "orchestration core resource identity requires core_resources_url"
            )
        if self.core_resources_version is not None and not re.fullmatch(
            r"v[0-9]+\.[0-9]+\.[0-9]+", self.core_resources_version
        ):
            errors.append(
                "orchestration.core_resources_version must have the form vMAJOR.MINOR.PATCH"
            )
        for name, value in (
            ("core_resources_sha256", self.core_resources_sha256),
            (
                "core_resources_manifest_sha256",
                self.core_resources_manifest_sha256,
            ),
        ):
            if value is not None and re.fullmatch(r"[0-9a-f]{64}", value) is None:
                errors.append(
                    f"orchestration.{name} must be a lowercase 64-character SHA-256"
                )
        if self.max_concurrent_genomes < 1:
            errors.append("orchestration.max_concurrent_genomes must be >= 1")
        if self.retries < 0:
            errors.append("orchestration.retries must be >= 0")
        if self.retry_delay_seconds < 0:
            errors.append("orchestration.retry_delay_seconds must be >= 0")
        if self.gpu_id is not None and self.gpu_id < 0:
            errors.append("orchestration.gpu_id must be >= 0 when set")
        return errors


@dataclass(frozen=True)
class FeatureResolution:
    """Requested and effective state for an optional runtime feature."""

    requested: bool
    required: bool
    enabled: bool
    reason_code: Optional[str] = None
    details: tuple[str, ...] = ()


_PIPELINE_SECTIONS = {
    "ablation",
    "databases",
    "compute",
    "host",
    "phase1",
    "phase2",
    "phase3",
    "execution",
}
_APPLICATION_KEYS = {"schema_version", "orchestration", *_PIPELINE_SECTIONS}
_ORCHESTRATION_FIELDS = {item.name for item in fields(OrchestrationConfig)}


# Compatibility aliases are deliberately centralized here. Each source key maps
# to exactly one canonical key and a collision with the canonical spelling is an
# error rather than an order-dependent override.
_LEGACY_ORCHESTRATION_ALIASES = {
    "hmm_db": "databases.hmm_database",
    "hmm_database": "databases.hmm_database",
    "hmm_allowlist": "databases.hmm_allowlist",
    "marker_faa_dir": "databases.marker_faa_dir",
    "marker_faa_db": "databases.marker_faa_db",
    "marker_db": "databases.marker_db",
    "seed_marker_allowlist": "databases.seed_marker_allowlist",
    "gene_taxonomy_faa_db": "databases.gene_taxonomy_faa_db",
    "faa_dir": "databases.faa_dir",
    "gvclass_db": "databases.gvclass_db",
    "diamond_db": "databases.diamond_db",
    "taxonomy_labels_file": "databases.taxonomy_labels_file",
    "threads": "compute.threads",
    "threads_per_worker": "compute.threads",
    "threads_per_genome": "compute.threads",
    "max_threads": "compute.max_threads",
    "device": "compute.device",
    "search_backend": "compute.search_backend",
    "gene_taxonomy_threads": "compute.gene_taxonomy_threads",
    "interproscan_threads": "compute.interproscan_threads",
    "rebuild_db": "phase1.rebuild_db",
    "hmm_chunk_size": "phase1.hmm_chunk_size",
    "assembly_mode": "phase1.assembly_mode",
    "enable_phylogenetic": "phase3.enable_phylogenetic",
    "skip_structural": "phase3.skip_structural",
    "use_structural_homology": "phase3.use_boltz",
    "use_boltz": "phase3.use_boltz",
    "boltz_mcp_only": "phase3.boltz_mcp_only",
    "boltz_use_msa_server": "phase3.boltz_use_msa_server",
    "boltz_min_seq_len": "phase3.boltz_min_seq_len",
    "boltz_max_seq_len": "phase3.boltz_max_seq_len",
    "boltz_no_kernels": "phase3.boltz_no_kernels",
    "use_tmvec_database": "phase3.use_tmvec_database",
    "tmvec_require_gpu": "phase3.tmvec_require_gpu",
    "tmvec_databases": "phase3.tmvec_databases",
    "tmvec_database_dir": "phase3.tmvec_database_dir",
    "tmvec_min_score": "phase3.tmvec_min_score",
    "viral_structure_db": "phase3.viral_structure_db",
    "export_all_eve_sequences": "phase3.export_all_eve_sequences",
    "run_gvclass": "phase3.run_gvclass",
    "gvclass_path": "phase3.gvclass_path",
    "interproscan_enabled": "phase3.interproscan_enabled",
    "interproscan_dir": "phase3.interproscan_dir",
    "interproscan_keywords": "phase3.interproscan_keywords",
    "interproscan_applications": "phase3.interproscan_applications",
    "resume": "execution.resume",
    "workers": "orchestration.max_concurrent_genomes",
}

_RETIRED_ORCHESTRATION_ALIASES = {
    "tier1_model_path": "databases.tier1_model",
    "tier2_model_path": "databases.tier2_model",
}


def _assign_alias(target: dict, dotted: str, value: Any, source: str) -> None:
    parts = dotted.split(".")
    current = target
    for part in parts[:-1]:
        existing = current.get(part)
        if existing is None:
            existing = {}
            current[part] = existing
        if not isinstance(existing, dict):
            raise ConfigError(
                f"Configuration alias '{source}' conflicts with section '{part}'"
            )
        current = existing
    leaf = parts[-1]
    if leaf in current:
        raise ConfigError(
            f"Configuration alias '{source}' conflicts with canonical key '{dotted}'"
        )
    current[leaf] = copy.deepcopy(value)


def _require_mapping(value: Any, path: str) -> dict:
    if not isinstance(value, dict):
        raise ConfigError(f"Configuration section '{path}' must be a mapping")
    return value


def _normalize_nested_aliases(normalized: dict) -> None:
    phase1 = _require_mapping(normalized.get("phase1", {}), "phase1")
    if "composition_clusterer" in phase1:
        raise ConfigError(
            "Unknown configuration key 'phase1.composition_clusterer': this "
            "retired control never affected runtime behavior; remove it "
            "(no replacement)."
        )
    hhg = phase1.pop("hhg_seeding", None)
    if hhg is not None:
        hhg = _require_mapping(hhg, "phase1.hhg_seeding")
        allowed = {
            "marker_top_k",
            "novel_marker_min_score",
            "novel_marker_min_coverage",
            "novel_marker_require_cluster",
            "validate_hmm_only",
        }
        for key in hhg:
            if key not in allowed:
                raise _unknown_key("phase1.hhg_seeding", str(key), allowed)
        if "validate_hmm_only" in hhg:
            raise ConfigError(
                "Configuration key 'phase1.hhg_seeding.validate_hmm_only' was "
                "removed because HMM-only validation is not a supported mode"
            )
        hhg_aliases = {
            "marker_top_k": "phase1.marker_validation_top_k",
            "novel_marker_min_score": "phase1.novel_marker_min_score",
            "novel_marker_min_coverage": "phase1.novel_marker_min_coverage",
            "novel_marker_require_cluster": "phase1.novel_marker_require_cluster",
        }
        for key, dotted in hhg_aliases.items():
            if key in hhg:
                _assign_alias(normalized, dotted, hhg[key], f"phase1.hhg_seeding.{key}")

    phase2 = _require_mapping(normalized.get("phase2", {}), "phase2")
    diamond = phase2.pop("boundary_diamond", None)
    if diamond is not None:
        diamond = _require_mapping(diamond, "phase2.boundary_diamond")
        aliases = {
            "flank_genes": "phase2.diamond_flank_genes",
            "control_sample_size": "phase2.diamond_control_sample_size",
            "control_min_distance": "phase2.diamond_control_min_distance",
            "top_k": "phase2.diamond_top_k",
            "chunk_size": "phase2.diamond_chunk_size",
            "random_seed": "phase2.diamond_random_seed",
            "superset_prototype_enabled": (
                "phase2.diamond_superset_prototype_enabled"
            ),
        }
        removed = {
            "threads": (
                "phase2.boundary_diamond.threads was removed; use "
                "compute.gene_taxonomy_threads"
            ),
            "host_prefix": (
                "phase2.boundary_diamond.host_prefix was removed; the primary "
                "prefix is derived from host.label"
            ),
        }
        for key in diamond:
            if key in removed:
                raise ConfigError(removed[key])
            if key not in aliases:
                raise _unknown_key(
                    "phase2.boundary_diamond",
                    str(key),
                    set(aliases) | set(removed),
                )
        for key, dotted in aliases.items():
            if key in diamond:
                _assign_alias(
                    normalized,
                    dotted,
                    diamond[key],
                    f"phase2.boundary_diamond.{key}",
                )
    if "host_trim" in phase2:
        host_trim = _require_mapping(phase2.pop("host_trim"), "phase2.host_trim")
        if "use_control_baseline" in host_trim:
            raise ConfigError(
                "Configuration key 'phase2.host_trim.use_control_baseline' was "
                "removed because control-baseline adjustment is not implemented"
            )
        if host_trim:
            key = next(iter(host_trim))
            raise _unknown_key("phase2.host_trim", str(key), set())

    phase3 = _require_mapping(normalized.get("phase3", {}), "phase3")
    if "use_structural_homology" in phase3:
        _assign_alias(
            normalized,
            "phase3.use_boltz",
            phase3.pop("use_structural_homology"),
            "phase3.use_structural_homology",
        )
    if "high_pident_euk_threshold" in phase3:
        _assign_alias(
            normalized,
            "host.high_pident_threshold",
            phase3.pop("high_pident_euk_threshold"),
            "phase3.high_pident_euk_threshold",
        )

    host = _require_mapping(normalized.get("host", {}), "host")
    if "high_pident_euk_threshold" in host:
        _assign_alias(
            normalized,
            "host.high_pident_threshold",
            host.pop("high_pident_euk_threshold"),
            "host.high_pident_euk_threshold",
        )


def _normalize_application_dict(data: dict) -> dict:
    if not isinstance(data, dict):
        raise ConfigError("Application configuration must be a mapping")
    for key in data:
        if not isinstance(key, str):
            raise ConfigError(f"Configuration keys must be strings, got {key!r}")
        if key not in _APPLICATION_KEYS:
            raise _unknown_key("", key, _APPLICATION_KEYS)

    normalized: dict[str, Any] = {}
    if "schema_version" in data:
        normalized["schema_version"] = copy.deepcopy(data["schema_version"])
    for section in _PIPELINE_SECTIONS:
        if section in data:
            normalized[section] = copy.deepcopy(data[section])

    raw_orchestration = _require_mapping(data.get("orchestration", {}), "orchestration")
    normalized["orchestration"] = {}
    for key, value in raw_orchestration.items():
        if not isinstance(key, str):
            raise ConfigError(
                f"Configuration keys must be strings, got orchestration.{key!r}"
            )
        if key in _ORCHESTRATION_FIELDS:
            _assign_alias(
                normalized,
                f"orchestration.{key}",
                value,
                f"orchestration.{key}",
            )
        elif key in _RETIRED_ORCHESTRATION_ALIASES:
            canonical = _RETIRED_ORCHESTRATION_ALIASES[key]
            raise ConfigError(
                f"Unknown configuration key 'orchestration.{key}': retired alias "
                f"for '{canonical}'; remove it (no replacement)."
            )
        elif key in _LEGACY_ORCHESTRATION_ALIASES:
            _assign_alias(
                normalized,
                _LEGACY_ORCHESTRATION_ALIASES[key],
                value,
                f"orchestration.{key}",
            )
        else:
            raise _unknown_key(
                "orchestration",
                key,
                _ORCHESTRATION_FIELDS | set(_LEGACY_ORCHESTRATION_ALIASES),
            )

    _normalize_nested_aliases(normalized)
    return normalized


def _serialize(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Enum):
        return value.value
    if is_dataclass(value):
        return {
            item.name: _serialize(getattr(value, item.name)) for item in fields(value)
        }
    if isinstance(value, dict):
        return {str(key): _serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serialize(item) for item in value]
    return value


@dataclass
class ApplicationConfig:
    """Versioned process and pipeline configuration loaded without side effects."""

    schema_version: int
    orchestration: OrchestrationConfig
    pipeline: PipelineConfig

    @classmethod
    def from_yaml(cls, path: Path) -> "ApplicationConfig":
        path = Path(path)
        try:
            with path.open(encoding="utf-8") as handle:
                raw = yaml.load(handle, Loader=_UniqueKeyLoader)
        except OSError as exc:
            raise ConfigError(
                f"Cannot read configuration file '{path}': {exc}"
            ) from exc
        except yaml.YAMLError as exc:
            raise ConfigError(
                f"Invalid YAML in configuration file '{path}': {exc}"
            ) from exc
        if raw is None:
            raise ConfigError(f"Configuration file '{path}' is empty")
        return cls.from_dict(raw, base_dir=path.parent.resolve())

    @classmethod
    def from_dict(
        cls,
        data: dict,
        *,
        base_dir: Optional[Path] = None,
    ) -> "ApplicationConfig":
        normalized = _normalize_application_dict(data)
        # Pre-v1 files had no version marker. Treat only that absence as v1 so
        # documented legacy aliases remain migratable; explicit versions stay strict.
        schema_version = normalized.get("schema_version", 1)
        if type(schema_version) is not int:
            raise ConfigError("Configuration key 'schema_version' must be an integer")
        if schema_version != 1:
            raise ConfigError(
                f"Unsupported schema_version {schema_version!r}; expected 1"
            )

        config_dir = Path(base_dir) if base_dir is not None else None
        orchestration = _decode_dataclass(
            OrchestrationConfig,
            normalized["orchestration"],
            "orchestration",
            config_dir,
        )
        orchestration_errors = orchestration.validate_semantics()
        if orchestration_errors:
            raise ConfigError(
                "Invalid orchestration configuration: "
                + "; ".join(orchestration_errors)
            )
        pipeline = PipelineConfig.from_dict(
            {section: normalized.get(section, {}) for section in _PIPELINE_SECTIONS},
            base_dir=config_dir,
        )
        return cls(
            schema_version=schema_version,
            orchestration=orchestration,
            pipeline=pipeline,
        )

    def to_dict(self) -> dict[str, Any]:
        """Return the canonical, JSON/YAML-safe representation including nulls."""
        pipeline = _serialize(self.pipeline)
        return {
            "schema_version": self.schema_version,
            "orchestration": _serialize(self.orchestration),
            **pipeline,
        }

    def effective_payload(
        self,
        optional_features: Optional[dict[str, FeatureResolution]] = None,
    ) -> dict[str, Any]:
        """Return deterministic runtime configuration for provenance output."""
        payload = self.to_dict()
        payload["optional_features"] = _serialize(optional_features or {})
        canonical = json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        payload["effective_config_sha256"] = hashlib.sha256(canonical).hexdigest()
        return payload

    def to_yaml(self, path: Path) -> None:
        """Write the canonical schema without invoking resource resolution."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(self.to_dict(), handle, sort_keys=False)
