"""
ViroSync Pipeline Configuration Dataclasses.

Provides a unified configuration system for the ViroSync orchestration pipeline.
Consolidates 50+ parameters into nested dataclasses with validation.
"""

import types
from dataclasses import dataclass, field, fields, is_dataclass, replace
from difflib import get_close_matches
from enum import Enum
from pathlib import Path
from typing import Any, Optional, Union, get_args, get_origin, get_type_hints

from virosync.ablation import ABLATION_CONTRACT_SHA256, AblationID


class ConfigError(ValueError):
    """A source configuration value does not satisfy the public schema."""


def _unknown_key(path: str, key: str, allowed: set[str]) -> ConfigError:
    dotted = f"{path}.{key}" if path else key
    matches = get_close_matches(key, sorted(allowed), n=1, cutoff=0.6)
    hint = (
        f" Did you mean '{path + '.' if path else ''}{matches[0]}'?" if matches else ""
    )
    return ConfigError(f"Unknown configuration key '{dotted}'.{hint}")


def _decode_value(annotation, value, path: str, base_dir: Optional[Path]):
    origin = get_origin(annotation)
    args = get_args(annotation)

    if origin in (Union, types.UnionType):
        if value is None and type(None) in args:
            return None
        candidates = [arg for arg in args if arg is not type(None)]
        failures = []
        for candidate in candidates:
            try:
                return _decode_value(candidate, value, path, base_dir)
            except ConfigError as exc:
                failures.append(str(exc))
        raise ConfigError(
            f"Invalid value for '{path}': {value!r}. {'; '.join(failures)}"
        )

    if value is None:
        raise ConfigError(f"Configuration key '{path}' may not be null")

    if origin is list:
        if not isinstance(value, list):
            raise ConfigError(f"Configuration key '{path}' must be a list")
        item_type = args[0] if args else Any
        return [
            _decode_value(item_type, item, f"{path}[{index}]", base_dir)
            for index, item in enumerate(value)
        ]

    if origin is tuple:
        if not isinstance(value, (list, tuple)):
            raise ConfigError(f"Configuration key '{path}' must be a list")
        item_type = args[0] if args else Any
        return tuple(
            _decode_value(item_type, item, f"{path}[{index}]", base_dir)
            for index, item in enumerate(value)
        )

    if annotation is Any:
        return value
    if annotation is bool:
        if type(value) is not bool:
            raise ConfigError(f"Configuration key '{path}' must be a boolean")
        return value
    if annotation is int:
        if type(value) is not int:
            raise ConfigError(f"Configuration key '{path}' must be an integer")
        return value
    if annotation is float:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ConfigError(f"Configuration key '{path}' must be a number")
        return float(value)
    if annotation is str:
        if not isinstance(value, str):
            raise ConfigError(f"Configuration key '{path}' must be a string")
        return value
    if annotation is Path:
        if not isinstance(value, (str, Path)):
            raise ConfigError(f"Configuration key '{path}' must be a path string")
        decoded = Path(value).expanduser()
        if base_dir is not None and not decoded.is_absolute():
            decoded = base_dir / decoded
        return decoded
    if isinstance(annotation, type) and issubclass(annotation, Enum):
        if not isinstance(value, str):
            raise ConfigError(f"Configuration key '{path}' must be a string enum value")
        try:
            return annotation(value)
        except ValueError as exc:
            choices = ", ".join(repr(member.value) for member in annotation)
            raise ConfigError(
                f"Configuration key '{path}' must be one of: {choices}"
            ) from exc
    if isinstance(annotation, type) and is_dataclass(annotation):
        return _decode_dataclass(annotation, value, path, base_dir)
    if isinstance(annotation, type) and isinstance(value, annotation):
        return value
    raise ConfigError(f"Unsupported schema type for '{path}': {annotation!r}")


def _decode_dataclass(cls, data, path: str, base_dir: Optional[Path] = None):
    if not isinstance(data, dict):
        label = path or cls.__name__
        raise ConfigError(f"Configuration section '{label}' must be a mapping")
    type_hints = get_type_hints(cls)
    allowed = {item.name for item in fields(cls)}
    for key in data:
        if key not in allowed:
            raise _unknown_key(path, str(key), allowed)
    kwargs = {
        key: _decode_value(
            type_hints[key],
            value,
            f"{path}.{key}" if path else key,
            base_dir,
        )
        for key, value in data.items()
    }
    try:
        return cls(**kwargs)
    except (TypeError, ValueError) as exc:
        label = path or cls.__name__
        raise ConfigError(f"Invalid configuration section '{label}': {exc}") from exc


class AssemblyMode(str, Enum):
    """HMM assembly mode for seed detection."""

    DEFAULT = "default"
    FRAGMENTED = "fragmented"
    RELAXED = "relaxed"
    STRICT = "strict"


class Device(str, Enum):
    """Compute device for GPU-accelerated tasks."""

    CUDA = "cuda"
    CPU = "cpu"


class SearchBackend(str, Enum):
    """Sequence search backend. Diamond is the sole supported backend (since v1.0.5)."""

    DIAMOND = "diamond"


class MaskingBackend(str, Enum):
    """Requested repeat-masking backend."""

    OFF = "off"
    TRF = "trf"
    REPEATMASKER = "repeatmasker"
    TRF_REPEATMASKER = "trf_repeatmasker"


class MaskingFailurePolicy(str, Enum):
    """Behavior when a requested masking backend fails."""

    STRICT = "strict"
    FALLBACK = "fallback"


@dataclass
class DatabasePaths:
    """Database and reference file paths for pipeline execution."""

    # HMM Detection
    hmm_database: Optional[Path] = None
    hmm_allowlist: Optional[Path] = None

    # Marker Validation
    marker_faa_dir: Optional[Path] = None
    marker_faa_db: Optional[Path] = None
    marker_db: Optional[Path] = None
    seed_marker_allowlist: Optional[list[str]] = None

    # Gene Taxonomy (Phase 3)
    gene_taxonomy_faa_db: Optional[Path] = None

    # Novelty/Taxonomy
    faa_dir: Optional[Path] = None

    # Phylogenetic Validation
    gvclass_db: Optional[Path] = None
    diamond_db: Optional[Path] = None

    # Taxonomy Labels (for host signature comparison)
    taxonomy_labels_file: Optional[Path] = None

    def __post_init__(self):
        """Convert string paths to Path objects."""
        for attr in [
            "hmm_database",
            "hmm_allowlist",
            "marker_faa_dir",
            "marker_faa_db",
            "marker_db",
            "gene_taxonomy_faa_db",
            "faa_dir",
            "gvclass_db",
            "diamond_db",
            "taxonomy_labels_file",
        ]:
            value = getattr(self, attr)
            if value is not None and not isinstance(value, Path):
                setattr(self, attr, Path(value))


@dataclass
class ComputeConfig:
    """Compute resource settings."""

    threads: int = 8
    max_threads: Optional[int] = None
    device: Device = Device.CPU
    search_backend: SearchBackend = SearchBackend.DIAMOND
    gene_taxonomy_threads: Optional[int] = None
    interproscan_threads: Optional[int] = None

    def __post_init__(self):
        """Convert device/search_backend strings to enum if needed."""
        if isinstance(self.device, str):
            self.device = Device(self.device)
        if isinstance(self.search_backend, str):
            self.search_backend = SearchBackend(self.search_backend)

    def effective_threads(self) -> int:
        """Return threads capped by max_threads if set."""
        if self.max_threads:
            return min(self.threads, self.max_threads)
        return self.threads


@dataclass
class HostConfig:
    """Host taxonomy configuration for host-specific heuristics."""

    prefixes: list[str] = field(default_factory=lambda: ["EUK__", "MITO__", "PLASTID__"])
    label: str = "EUK"
    high_pident_threshold: float = 70.0


@dataclass
class AblationConfig:
    """Closed benchmark intervention selected for this pipeline run."""

    id: AblationID = AblationID.A0
    contract_sha256: str = ABLATION_CONTRACT_SHA256

    def __post_init__(self) -> None:
        if isinstance(self.id, str):
            self.id = AblationID(self.id)
        if not isinstance(self.id, AblationID):
            raise ValueError("ablation.id must be one of A0-A6")
        if self.contract_sha256 != ABLATION_CONTRACT_SHA256:
            raise ValueError(
                "ablation.contract_sha256 does not match this ViroSync build"
            )


@dataclass
class Phase1Config:
    """Phase 1: Seeding parameters."""

    rebuild_db: bool = False

    # HMM settings
    assembly_mode: AssemblyMode = AssemblyMode.DEFAULT
    hmm_chunk_size: Optional[int] = None
    frameshift_screening_enabled: bool = False

    # Region assembly
    initial_window_bp: int = 10000
    initial_window_genes: int = 5
    min_markers_initial: int = 1
    extension_kb: int = 5
    merge_distance: int = 1000

    # Host taxonomy deviation (Phase 1 optional extension signal)
    host_taxonomy_deviation_enabled: bool = False
    host_taxonomy_deviation_allow_seeds: bool = False
    host_taxonomy_deviation_min_token_len: int = 6
    host_taxonomy_deviation_min_tokens: int = 3
    host_taxonomy_deviation_overlap_threshold: float = 0.3
    host_taxonomy_deviation_max_pident: float = 70.0
    host_taxonomy_deviation_max_hits: int = 5
    host_taxonomy_deviation_window_bp: int = 5000
    host_taxonomy_deviation_window_count: int = 25
    host_taxonomy_deviation_window_seed: int = 13
    host_taxonomy_deviation_window_min_markers: int = 1
    host_taxonomy_deviation_seed_window_bp: int = 10000
    host_taxonomy_deviation_seed_min_markers: int = 3

    # HMM-gated marker validation
    marker_validation_top_k: int = 10
    novel_marker_min_score: float = 30.0
    novel_marker_min_coverage: float = 0.5
    novel_marker_require_cluster: bool = True

    def __post_init__(self):
        """Convert assembly_mode string to enum if needed."""
        if isinstance(self.assembly_mode, str):
            self.assembly_mode = AssemblyMode(self.assembly_mode)


@dataclass
class Phase2Config:
    """Phase 2: Boundary refinement parameters."""

    taxonomy_weight_mode: str = "rank"  # "rank" (10..1) or "bitscore" for taxonomy weights
    # Taxonomy boundary ML refinement (Phase 2d/2e classifier backend)
    taxonomy_ml_enabled: bool = False
    taxonomy_ml_model: str = "logreg"  # "logreg", "gbdt", or "xgboost"
    taxonomy_ml_threshold: float = 0.5
    taxonomy_ml_neighbor_window: int = 3
    # Phase 2a: host-signature trimming
    host_trim_enabled: bool = True
    host_trim_window_bp: int = 5000
    host_trim_step_bp: int = 1000
    host_trim_max_host_fraction: float = 0.3
    host_trim_min_viral_fraction: float = 0.05
    host_trim_score_threshold: float = 0.3
    host_trim_buffer_kb: int = 5
    host_signature_min_token_len: int = 3
    host_trim_min_overlap_score: float = 0.40  # Min fingerprint overlap for host classification
    # Phase 2b: batched Diamond for boundary refinement (always enabled as of Jan 2026)
    diamond_flank_genes: int = 10  # Seeds already extended by ±5 genes; ±10 flank gives ±15 total
    diamond_control_sample_size: int = 100
    diamond_control_min_distance: int = 30
    diamond_top_k: int = 10
    diamond_chunk_size: int = 10000
    diamond_random_seed: int = 42
    diamond_superset_prototype_enabled: bool = False

@dataclass
class Phase3Config:
    """Phase 3: Evidence synthesis parameters."""

    host_signature_evidence_threshold: float = 0.3

    # Confidence tier thresholds (for output classification)
    # HIGH: confidence >= high_tier_threshold
    # MEDIUM: low_tier_threshold <= confidence < high_tier_threshold
    # LOW: confidence < low_tier_threshold
    high_tier_threshold: float = 0.7
    low_tier_threshold: float = 0.2
    use_crf_in_final_score: bool = False
    priority_marker_list: list[str] = field(default_factory=lambda: ["mcp"])
    marker_floor_priority_only: float = 0.55
    marker_floor_priority_plus_family: float = 0.70
    marker_floor_priority_multi_family: float = 0.80
    marker_family_bonus_per_family: float = 0.06
    marker_multi_family_bonus: float = 0.08
    enable_phylogenetic: bool = False
    skip_structural: bool = True
    use_boltz: bool = False
    boltz_mcp_only: bool = True
    boltz_use_msa_server: bool = False
    boltz_min_seq_len: int = 100
    boltz_max_seq_len: int = 1000
    boltz_no_kernels: bool = True  # Use --no_kernels flag for safer Boltz execution
    use_tmvec_database: bool = False
    tmvec_require_gpu: bool = False
    tmvec_databases: Optional[list[str]] = None
    tmvec_database_dir: Optional[Path] = None
    tmvec_min_score: float = 0.5
    viral_structure_db: Optional[Path] = None
    extended_output: bool = True
    export_all_eve_sequences: bool = True

    # InterProScan
    interproscan_enabled: bool = False
    interproscan_dir: Optional[Path] = None
    interproscan_keywords: Optional[list[str]] = None
    interproscan_applications: Optional[list[str]] = None

    # GVClass batch classification
    run_gvclass: bool = False
    gvclass_path: Optional[Path] = None

    def __post_init__(self):
        """Convert path fields to Path if needed."""
        if self.interproscan_dir is not None and not isinstance(
            self.interproscan_dir, Path
        ):
            self.interproscan_dir = Path(self.interproscan_dir)
        if self.viral_structure_db is not None and not isinstance(
            self.viral_structure_db, Path
        ):
            self.viral_structure_db = Path(self.viral_structure_db)
        if self.tmvec_database_dir is not None and not isinstance(
            self.tmvec_database_dir, Path
        ):
            self.tmvec_database_dir = Path(self.tmvec_database_dir)
        if self.gvclass_path is not None and not isinstance(
            self.gvclass_path, Path
        ):
            self.gvclass_path = Path(self.gvclass_path)


@dataclass(frozen=True)
class MaskingConfig:
    """Validated repeat-masking request and explicit fallback policy."""

    backend: MaskingBackend = MaskingBackend.OFF
    failure_policy: MaskingFailurePolicy = MaskingFailurePolicy.STRICT
    fallback_backend: Optional[MaskingBackend] = None
    repeatmasker_species: Optional[str] = None
    repeatmasker_library: Optional[Path] = None

    def __post_init__(self) -> None:
        if isinstance(self.backend, str):
            object.__setattr__(self, "backend", MaskingBackend(self.backend))
        if isinstance(self.failure_policy, str):
            object.__setattr__(
                self,
                "failure_policy",
                MaskingFailurePolicy(self.failure_policy),
            )
        if isinstance(self.fallback_backend, str):
            object.__setattr__(
                self,
                "fallback_backend",
                MaskingBackend(self.fallback_backend),
            )
        if self.repeatmasker_library is not None and not isinstance(
            self.repeatmasker_library,
            Path,
        ):
            object.__setattr__(
                self,
                "repeatmasker_library",
                Path(self.repeatmasker_library),
            )
        if self.repeatmasker_species is not None:
            species = self.repeatmasker_species.strip()
            object.__setattr__(self, "repeatmasker_species", species or None)

    @property
    def uses_repeatmasker(self) -> bool:
        return self.backend in {
            MaskingBackend.REPEATMASKER,
            MaskingBackend.TRF_REPEATMASKER,
        }

    def with_backend(self, backend: MaskingBackend | str) -> "MaskingConfig":
        """Return this request with a different backend and the same latent target."""
        selected = MaskingBackend(backend)
        if selected is MaskingBackend.OFF:
            return replace(
                self,
                backend=selected,
                failure_policy=MaskingFailurePolicy.STRICT,
                fallback_backend=None,
            )
        return replace(self, backend=selected)

    def validate(self) -> list[str]:
        """Return semantic errors without touching configured files."""
        errors: list[str] = []
        has_species = self.repeatmasker_species is not None
        has_library = self.repeatmasker_library is not None
        if has_species and has_library:
            errors.append(
                "execution.masking must set exactly one of repeatmasker_species "
                "or repeatmasker_library, not both"
            )
        if self.uses_repeatmasker and not (has_species or has_library):
            errors.append(
                "execution.masking backend requires exactly one of "
                "repeatmasker_species or repeatmasker_library"
            )
        if self.failure_policy is MaskingFailurePolicy.STRICT:
            if self.fallback_backend is not None:
                errors.append(
                    "execution.masking.fallback_backend must be null under strict policy"
                )
        else:
            if self.fallback_backend not in {
                MaskingBackend.OFF,
                MaskingBackend.TRF,
            }:
                errors.append(
                    "execution.masking fallback policy requires fallback_backend "
                    "to be 'off' or 'trf'"
                )
            elif self.fallback_backend is self.backend:
                errors.append(
                    "execution.masking.fallback_backend must differ from backend"
                )
            if self.backend is MaskingBackend.OFF:
                errors.append(
                    "execution.masking fallback policy is invalid when backend is 'off'"
                )
        return errors


@dataclass
class ExecutionConfig:
    """Execution control settings."""

    masking: MaskingConfig = field(default_factory=MaskingConfig)
    resume: bool = True

    @property
    def skip_masking(self) -> bool:
        """Derived compatibility view; canonical state lives in ``masking``."""
        return self.masking.backend is MaskingBackend.OFF

_UNSET = object()


@dataclass(frozen=True)
class FieldSpec:
    """Single source of truth mapping one flat flow-kwarg to its nested config field.

    Drives all three config maps:
    - ``with_overrides`` (flat kwarg -> nested field), incl. enum/Path coercion + aliases
    - ``to_flow_kwargs`` (nested field -> flat kwarg), incl. enum ``.value`` emission
    - ``_from_flat_dict`` (YAML/flat dict -> nested field) for the table-driven sections;
      genuinely special parsing (databases legacy aliases, host block, high_pident
      coupling, use_boltz cascade) stays explicit in that method.
    """

    flat: str                       # canonical flat name (with_overrides key + to_flow_kwargs output key)
    section: str                    # nested config attribute on PipelineConfig
    field: str                      # dataclass field name on the section
    enum: Optional[type] = None     # enum class: with_overrides coerces str; to_flow_kwargs emits .value
    path: bool = False              # with_overrides coerces str -> Path
    overridable: bool = True        # appears in with_overrides
    emit: bool = True               # appears in to_flow_kwargs
    wo_aliases: tuple = ()          # extra override-only flat names accepted by with_overrides
    ff_block: Optional[str] = None  # _from_flat_dict source block: None|phase1|phase2|phase3|boundary_diamond
    ff_section_key: Optional[str] = None  # key within ff_block (defaults to field)
    ff_flat_aliases: tuple = ()     # extra flat keys tried (after `flat`) in _from_flat_dict
    ff_default: object = _UNSET     # _from_flat_dict default override (else dataclass default)
    ff_skip: bool = False           # not table-driven in _from_flat_dict (handled explicitly)


def _spec(flat, section, field, **kw):
    return FieldSpec(flat=flat, section=section, field=field, **kw)


_HTD = "host_taxonomy_deviation_"

FIELD_SPECS: list[FieldSpec] = [
    # --- ablation (mutually exclusive benchmark arm) ---
    _spec("ablation_id", "ablation", "id", enum=AblationID),
    _spec(
        "ablation_contract_sha256",
        "ablation",
        "contract_sha256",
        overridable=False,
    ),
    # --- databases (legacy 'or' aliases handled explicitly in _from_flat_dict) ---
    _spec("hmm_database", "databases", "hmm_database", path=True, ff_skip=True),
    _spec("hmm_allowlist", "databases", "hmm_allowlist", path=True, ff_skip=True),
    _spec("marker_faa_dir", "databases", "marker_faa_dir", path=True, ff_skip=True),
    _spec("marker_faa_db", "databases", "marker_faa_db", path=True, ff_skip=True),
    _spec("marker_db", "databases", "marker_db", path=True, ff_skip=True),
    _spec("gene_taxonomy_faa_db", "databases", "gene_taxonomy_faa_db", path=True, ff_skip=True),
    _spec("seed_marker_allowlist", "databases", "seed_marker_allowlist", ff_skip=True),
    _spec("faa_dir", "databases", "faa_dir", path=True, ff_skip=True),
    _spec("gvclass_db", "databases", "gvclass_db", path=True, ff_skip=True),
    _spec("diamond_db", "databases", "diamond_db", path=True, ff_skip=True),
    _spec(
        "taxonomy_labels_file",
        "databases",
        "taxonomy_labels_file",
        path=True,
        ff_skip=True,
    ),
    # --- compute (flat-only in _from_flat_dict) ---
    _spec("threads", "compute", "threads"),
    _spec("max_threads", "compute", "max_threads"),
    _spec("device", "compute", "device", enum=Device),
    _spec("search_backend", "compute", "search_backend", enum=SearchBackend),
    _spec("gene_taxonomy_threads", "compute", "gene_taxonomy_threads"),
    _spec("interproscan_threads", "compute", "interproscan_threads"),
    # --- host (legacy aliases + high_pident coupling handled explicitly) ---
    _spec("host_prefixes", "host", "prefixes", ff_skip=True),
    _spec("host_label", "host", "label", ff_skip=True),
    _spec("high_pident_host_threshold", "host", "high_pident_threshold", ff_skip=True),
    # --- phase1 (flat names == field names) ---
    _spec("rebuild_db", "phase1", "rebuild_db", ff_block="phase1"),
    _spec("assembly_mode", "phase1", "assembly_mode", enum=AssemblyMode, ff_block="phase1"),
    _spec("hmm_chunk_size", "phase1", "hmm_chunk_size", ff_block="phase1"),
    _spec(
        "frameshift_screening_enabled",
        "phase1",
        "frameshift_screening_enabled",
        ff_block="phase1",
    ),
    _spec("initial_window_bp", "phase1", "initial_window_bp", ff_block="phase1"),
    _spec("initial_window_genes", "phase1", "initial_window_genes", ff_block="phase1"),
    _spec("min_markers_initial", "phase1", "min_markers_initial", ff_block="phase1"),
    _spec("extension_kb", "phase1", "extension_kb", ff_block="phase1"),
    _spec("merge_distance", "phase1", "merge_distance", ff_block="phase1"),
    _spec(_HTD + "enabled", "phase1", _HTD + "enabled", ff_block="phase1"),
    _spec(_HTD + "allow_seeds", "phase1", _HTD + "allow_seeds", ff_block="phase1"),
    _spec(_HTD + "min_token_len", "phase1", _HTD + "min_token_len", ff_block="phase1"),
    _spec(_HTD + "min_tokens", "phase1", _HTD + "min_tokens", ff_block="phase1"),
    _spec(_HTD + "overlap_threshold", "phase1", _HTD + "overlap_threshold", ff_block="phase1"),
    _spec(_HTD + "max_pident", "phase1", _HTD + "max_pident", ff_block="phase1"),
    _spec(_HTD + "max_hits", "phase1", _HTD + "max_hits", ff_block="phase1"),
    _spec(_HTD + "window_bp", "phase1", _HTD + "window_bp", ff_block="phase1"),
    _spec(_HTD + "window_count", "phase1", _HTD + "window_count", ff_block="phase1"),
    _spec(_HTD + "window_seed", "phase1", _HTD + "window_seed", ff_block="phase1"),
    _spec(
        _HTD + "window_min_markers",
        "phase1",
        _HTD + "window_min_markers",
        ff_block="phase1",
    ),
    _spec(
        _HTD + "seed_window_bp", "phase1", _HTD + "seed_window_bp", ff_block="phase1"
    ),
    _spec(
        _HTD + "seed_min_markers",
        "phase1",
        _HTD + "seed_min_markers",
        ff_block="phase1",
    ),
    _spec(
        "marker_validation_top_k",
        "phase1",
        "marker_validation_top_k",
        ff_block="phase1",
    ),
    _spec(
        "novel_marker_min_score", "phase1", "novel_marker_min_score", ff_block="phase1"
    ),
    _spec(
        "novel_marker_min_coverage",
        "phase1",
        "novel_marker_min_coverage",
        ff_block="phase1",
    ),
    _spec(
        "novel_marker_require_cluster",
        "phase1",
        "novel_marker_require_cluster",
        ff_block="phase1",
    ),
    # --- phase2 (flat names are boundary_-prefixed; section keys are the bare field names) ---
    _spec("taxonomy_weight_mode", "phase2", "taxonomy_weight_mode", ff_block="phase2"),
    _spec("boundary_taxonomy_ml_enabled", "phase2", "taxonomy_ml_enabled",
          wo_aliases=("taxonomy_ml_enabled",), ff_block="phase2", ff_flat_aliases=("taxonomy_ml_enabled",)),
    _spec("boundary_taxonomy_ml_model", "phase2", "taxonomy_ml_model",
          wo_aliases=("taxonomy_ml_model",), ff_block="phase2", ff_flat_aliases=("taxonomy_ml_model",)),
    _spec("boundary_taxonomy_ml_threshold", "phase2", "taxonomy_ml_threshold",
          wo_aliases=("taxonomy_ml_threshold",), ff_block="phase2", ff_flat_aliases=("taxonomy_ml_threshold",)),
    _spec("boundary_taxonomy_ml_neighbor_window", "phase2", "taxonomy_ml_neighbor_window",
          wo_aliases=("taxonomy_ml_neighbor_window",), ff_block="phase2",
          ff_flat_aliases=("taxonomy_ml_neighbor_window",)),
    _spec("boundary_host_trim_enabled", "phase2", "host_trim_enabled", ff_block="phase2"),
    _spec("boundary_host_trim_window_bp", "phase2", "host_trim_window_bp", ff_block="phase2"),
    _spec("boundary_host_trim_step_bp", "phase2", "host_trim_step_bp", ff_block="phase2"),
    _spec("boundary_host_trim_max_host_fraction", "phase2", "host_trim_max_host_fraction", ff_block="phase2"),
    _spec("boundary_host_trim_min_viral_fraction", "phase2", "host_trim_min_viral_fraction", ff_block="phase2"),
    _spec("boundary_host_trim_score_threshold", "phase2", "host_trim_score_threshold", ff_block="phase2"),
    _spec("boundary_host_trim_buffer_kb", "phase2", "host_trim_buffer_kb", ff_block="phase2"),
    _spec("boundary_host_trim_min_overlap_score", "phase2", "host_trim_min_overlap_score", ff_block="phase2"),
    # NOTE _from_flat_dict default for this is 6, diverging from the dataclass default (3).
    _spec("boundary_host_signature_min_token_len", "phase2", "host_signature_min_token_len",
          ff_block="phase2", ff_default=6),
    # phase2b Diamond + ANI: emitted but NOT overridable; parsed from the nested boundary_diamond block.
    _spec(
        "boundary_diamond_flank_genes",
        "phase2",
        "diamond_flank_genes",
        overridable=False,
        ff_block="boundary_diamond",
        ff_section_key="flank_genes",
    ),
    _spec(
        "boundary_diamond_control_sample_size",
        "phase2",
        "diamond_control_sample_size",
        overridable=False,
        ff_block="boundary_diamond",
        ff_section_key="control_sample_size",
    ),
    _spec(
        "boundary_diamond_control_min_distance",
        "phase2",
        "diamond_control_min_distance",
        overridable=False,
        ff_block="boundary_diamond",
        ff_section_key="control_min_distance",
    ),
    _spec(
        "boundary_diamond_top_k",
        "phase2",
        "diamond_top_k",
        overridable=False,
        ff_block="boundary_diamond",
        ff_section_key="top_k",
    ),
    _spec(
        "boundary_diamond_chunk_size",
        "phase2",
        "diamond_chunk_size",
        overridable=False,
        ff_block="boundary_diamond",
        ff_section_key="chunk_size",
    ),
    _spec(
        "boundary_diamond_random_seed",
        "phase2",
        "diamond_random_seed",
        overridable=False,
        ff_block="boundary_diamond",
        ff_section_key="random_seed",
    ),
    _spec(
        "boundary_diamond_superset_prototype_enabled",
        "phase2",
        "diamond_superset_prototype_enabled",
        overridable=False,
        ff_block="boundary_diamond",
        ff_section_key="superset_prototype_enabled",
    ),
    # --- phase3 (flat names == field names) ---
    _spec(
        "host_signature_evidence_threshold",
        "phase3",
        "host_signature_evidence_threshold",
        ff_block="phase3",
    ),
    _spec("high_tier_threshold", "phase3", "high_tier_threshold", ff_block="phase3"),
    _spec("low_tier_threshold", "phase3", "low_tier_threshold", ff_block="phase3"),
    _spec("use_crf_in_final_score", "phase3", "use_crf_in_final_score", ff_block="phase3"),
    _spec("priority_marker_list", "phase3", "priority_marker_list", ff_block="phase3"),
    _spec("marker_floor_priority_only", "phase3", "marker_floor_priority_only", ff_block="phase3"),
    _spec("marker_floor_priority_plus_family", "phase3", "marker_floor_priority_plus_family", ff_block="phase3"),
    _spec("marker_floor_priority_multi_family", "phase3", "marker_floor_priority_multi_family", ff_block="phase3"),
    _spec("marker_family_bonus_per_family", "phase3", "marker_family_bonus_per_family", ff_block="phase3"),
    _spec("marker_multi_family_bonus", "phase3", "marker_multi_family_bonus", ff_block="phase3"),
    _spec("enable_phylogenetic", "phase3", "enable_phylogenetic", ff_block="phase3"),
    _spec("skip_structural", "phase3", "skip_structural", ff_block="phase3"),
    _spec("use_boltz", "phase3", "use_boltz", ff_skip=True),
    _spec("boltz_mcp_only", "phase3", "boltz_mcp_only", ff_block="phase3"),
    _spec("boltz_use_msa_server", "phase3", "boltz_use_msa_server", ff_block="phase3"),
    _spec("boltz_min_seq_len", "phase3", "boltz_min_seq_len", ff_block="phase3"),
    _spec("boltz_max_seq_len", "phase3", "boltz_max_seq_len", ff_block="phase3"),
    _spec("boltz_no_kernels", "phase3", "boltz_no_kernels", ff_block="phase3"),
    _spec("use_tmvec_database", "phase3", "use_tmvec_database", ff_block="phase3"),
    _spec("tmvec_require_gpu", "phase3", "tmvec_require_gpu", ff_block="phase3"),
    _spec("tmvec_databases", "phase3", "tmvec_databases", ff_block="phase3"),
    _spec("tmvec_database_dir", "phase3", "tmvec_database_dir", path=True, ff_block="phase3"),
    _spec("tmvec_min_score", "phase3", "tmvec_min_score", ff_block="phase3"),
    _spec("viral_structure_db", "phase3", "viral_structure_db", path=True, ff_block="phase3"),
    _spec("extended_output", "phase3", "extended_output", ff_block="phase3"),
    _spec("export_all_eve_sequences", "phase3", "export_all_eve_sequences", ff_block="phase3"),
    _spec("interproscan_enabled", "phase3", "interproscan_enabled", ff_block="phase3"),
    _spec("interproscan_dir", "phase3", "interproscan_dir", path=True, ff_block="phase3"),
    _spec("interproscan_keywords", "phase3", "interproscan_keywords", ff_block="phase3"),
    _spec("interproscan_applications", "phase3", "interproscan_applications", ff_block="phase3"),
    _spec("run_gvclass", "phase3", "run_gvclass", ff_block="phase3"),
    _spec("gvclass_path", "phase3", "gvclass_path", path=True, ff_block="phase3"),
    # --- execution (flat-only) ---
    _spec("masking", "execution", "masking"),
    _spec("resume", "execution", "resume"),
]

# Lookups derived once from FIELD_SPECS.
_WO_LOOKUP: dict[str, FieldSpec] = {}
for _s in FIELD_SPECS:
    if _s.overridable:
        _WO_LOOKUP[_s.flat] = _s
        for _alias in _s.wo_aliases:
            _WO_LOOKUP[_alias] = _s
_EMIT_SPECS: list[FieldSpec] = [_s for _s in FIELD_SPECS if _s.emit]


_RETIRED_PIPELINE_KEYS = frozenset(
    {
        "databases.tier1_model",
        "databases.tier2_model",
        "phase1.composition_window_bp",
        "phase1.composition_step_bp",
        "phase1.composition_classifier_threshold",
        "phase1.composition_training_step_bp",
        "phase1.composition_expansion_step_bp",
        "phase1.composition_max_gap_without_marker_kb",
        "phase1.composition_max_region_bp",
        "phase2.mode",
        "phase2.window_bp",
        "phase2.step_bp",
        "phase2.classifier_threshold",
        "phase2.ensemble_strategy",
        "phase2.refinement_margin_kb",
        "phase2.ani_threshold",
        "phase2.min_cluster_bp",
    }
)


def _reject_retired_pipeline_keys(data: dict) -> None:
    """Reject accepted-but-inert v1 controls with a migration diagnostic."""
    for section, section_data in data.items():
        if not isinstance(section_data, dict):
            continue
        for key in section_data:
            dotted = f"{section}.{key}"
            if dotted in _RETIRED_PIPELINE_KEYS:
                raise ConfigError(
                    f"Unknown configuration key '{dotted}': this retired control "
                    "never affected runtime behavior; remove it (no replacement)."
                )


@dataclass
class PipelineConfig:
    """
    Complete ViroSync pipeline configuration.

    Usage:
        # From YAML file
        config = PipelineConfig.from_yaml("config.yaml")

        # Programmatic
        config = PipelineConfig(
            databases=DatabaseConfig(hmm_database=Path("markers.hmm")),
            phase1=Phase1Config(assembly_mode=AssemblyMode.RELAXED),
        )

        # In flow
        results = single_genome_flow(
            genome_path=genome,
            output_dir=output / "my_genome",
            genome_id="my_genome",
            config=config,
        )
    """

    # Nested configs with defaults
    ablation: AblationConfig = field(default_factory=AblationConfig)
    databases: DatabasePaths = field(default_factory=DatabasePaths)
    compute: ComputeConfig = field(default_factory=ComputeConfig)
    host: HostConfig = field(default_factory=HostConfig)
    phase1: Phase1Config = field(default_factory=Phase1Config)
    phase2: Phase2Config = field(default_factory=Phase2Config)
    phase3: Phase3Config = field(default_factory=Phase3Config)
    execution: ExecutionConfig = field(default_factory=ExecutionConfig)

    def validate_database_paths(self, check_files: bool = True) -> list[str]:
        """
        Validate that required database files exist on disk.

        Args:
            check_files: If True, verify files exist on disk. Set to False for tests.

        Returns:
            List of validation error messages.
        """
        errors = []

        if not check_files:
            return errors

        selected_paths: list[tuple[str, Path, str]] = []

        def selected(label: str, path: Optional[Path], kind: str = "path") -> None:
            if path is not None:
                selected_paths.append((label, Path(path), kind))

        selected("HMM database", self.databases.hmm_database, "file")
        selected("HMM allowlist", self.databases.hmm_allowlist, "file")
        selected(
            "Gene taxonomy Diamond database",
            self.databases.gene_taxonomy_faa_db,
            "file",
        )
        selected(
            "Taxonomy labels file",
            self.databases.taxonomy_labels_file,
            "file",
        )

        # A prebuilt marker DB supersedes the marker-FAA build inputs.
        if self.databases.marker_db is not None:
            selected(
                "Marker Diamond database",
                self.databases.marker_db,
                "file",
            )
        else:
            selected("Marker FAA database", self.databases.marker_faa_db, "file")
            selected("Marker FAA directory", self.databases.marker_faa_dir, "directory")
            selected("FAA build directory", self.databases.faa_dir, "directory")

        if self.phase3.enable_phylogenetic:
            selected("GVClass database", self.databases.gvclass_db)
            selected("Phylogenetic Diamond database", self.databases.diamond_db)
        elif self.phase3.run_gvclass:
            selected("GVClass database", self.databases.gvclass_db)

        masking_library = self.execution.masking.repeatmasker_library
        if self.execution.masking.uses_repeatmasker and masking_library is not None:
            if not masking_library.is_file():
                errors.append(
                    "RepeatMasker library is not a regular file: "
                    f"{masking_library}"
                )
            elif masking_library.stat().st_size == 0:
                errors.append(f"RepeatMasker library is empty: {masking_library}")

        seen: set[tuple[str, Path]] = set()
        for label, path, kind in selected_paths:
            identity = (label, path)
            if identity in seen:
                continue
            seen.add(identity)
            if kind == "file":
                exists = path.is_file()
            elif kind == "directory":
                exists = path.is_dir()
            else:
                exists = path.exists()
            if not exists:
                errors.append(f"{label} not found: {path}")

        return errors

    def validate(self) -> list[str]:
        """Validate semantic and configured resource requirements."""
        errors = self.validate_semantics()

        # HMM-gated workflow requires HMM database (only path supported)
        if not self.databases.hmm_database:
            errors.append("HMM-gated workflow requires databases.hmm_database")
        if self.databases.marker_db is None and (
            self.databases.faa_dir is None
            or (
                self.databases.marker_faa_dir is None
                and self.databases.marker_faa_db is None
            )
        ):
            errors.append(
                "HMM-gated workflow requires databases.marker_db, or "
                "databases.faa_dir plus marker_faa_db/marker_faa_dir"
            )
        if self.phase3.enable_phylogenetic:
            if not self.databases.gvclass_db and not self.databases.diamond_db:
                errors.append(
                    "Phylogenetic validation requires gvclass_db or diamond_db"
                )
        if self.phase3.interproscan_enabled and not self.phase3.interproscan_dir:
            errors.append("InterProScan enabled but interproscan_dir not set")
        return errors

    def validate_semantics(self) -> list[str]:
        """Validate types-independent ranges and cross-field invariants."""
        errors = list(self.execution.masking.validate())
        if not isinstance(self.ablation.id, AblationID):
            errors.append("ablation.id must be one of A0-A6")
        if self.ablation.contract_sha256 != ABLATION_CONTRACT_SHA256:
            errors.append(
                "ablation.contract_sha256 does not match this ViroSync build"
            )
        if self.ablation.id is AblationID.A4 and self.phase2.taxonomy_ml_enabled:
            errors.append(
                "ablation A4 cannot use phase2.taxonomy_ml_enabled because its "
                "model contains host-derived features"
            )
        if (
            self.phase2.diamond_superset_prototype_enabled
            and self.phase2.diamond_top_k != 10
        ):
            errors.append(
                "phase2.diamond_superset_prototype_enabled requires "
                "phase2.diamond_top_k=10 so both cached consumers use the same "
                "DIAMOND hit limit"
            )

        positive_ints = {
            "compute.threads": self.compute.threads,
            "phase1.initial_window_bp": self.phase1.initial_window_bp,
            "phase1.initial_window_genes": self.phase1.initial_window_genes,
            "phase1.min_markers_initial": self.phase1.min_markers_initial,
            "phase1.host_taxonomy_deviation_min_token_len": self.phase1.host_taxonomy_deviation_min_token_len,
            "phase1.host_taxonomy_deviation_min_tokens": self.phase1.host_taxonomy_deviation_min_tokens,
            "phase1.host_taxonomy_deviation_max_hits": self.phase1.host_taxonomy_deviation_max_hits,
            "phase1.host_taxonomy_deviation_window_bp": self.phase1.host_taxonomy_deviation_window_bp,
            "phase1.host_taxonomy_deviation_window_count": self.phase1.host_taxonomy_deviation_window_count,
            "phase1.host_taxonomy_deviation_window_min_markers": self.phase1.host_taxonomy_deviation_window_min_markers,
            "phase1.host_taxonomy_deviation_seed_window_bp": self.phase1.host_taxonomy_deviation_seed_window_bp,
            "phase1.host_taxonomy_deviation_seed_min_markers": self.phase1.host_taxonomy_deviation_seed_min_markers,
            "phase1.marker_validation_top_k": self.phase1.marker_validation_top_k,
            "phase2.host_trim_window_bp": self.phase2.host_trim_window_bp,
            "phase2.host_trim_step_bp": self.phase2.host_trim_step_bp,
            "phase2.host_signature_min_token_len": self.phase2.host_signature_min_token_len,
            "phase2.taxonomy_ml_neighbor_window": self.phase2.taxonomy_ml_neighbor_window,
            "phase2.diamond_control_sample_size": self.phase2.diamond_control_sample_size,
            "phase2.diamond_top_k": self.phase2.diamond_top_k,
            "phase2.diamond_chunk_size": self.phase2.diamond_chunk_size,
            "phase3.boltz_min_seq_len": self.phase3.boltz_min_seq_len,
            "phase3.boltz_max_seq_len": self.phase3.boltz_max_seq_len,
        }
        optional_positive = {
            "compute.max_threads": self.compute.max_threads,
            "compute.gene_taxonomy_threads": self.compute.gene_taxonomy_threads,
            "compute.interproscan_threads": self.compute.interproscan_threads,
            "phase1.hmm_chunk_size": self.phase1.hmm_chunk_size,
        }
        for name, value in positive_ints.items():
            if value < 1:
                errors.append(f"{name} must be >= 1")
        for name, value in optional_positive.items():
            if value is not None and value < 1:
                errors.append(f"{name} must be >= 1 when set")

        nonnegative_ints = {
            "phase1.extension_kb": self.phase1.extension_kb,
            "phase1.merge_distance": self.phase1.merge_distance,
            "phase1.host_taxonomy_deviation_window_seed": self.phase1.host_taxonomy_deviation_window_seed,
            "phase2.host_trim_buffer_kb": self.phase2.host_trim_buffer_kb,
            "phase2.diamond_flank_genes": self.phase2.diamond_flank_genes,
            "phase2.diamond_control_min_distance": self.phase2.diamond_control_min_distance,
            "phase2.diamond_random_seed": self.phase2.diamond_random_seed,
        }
        for name, value in nonnegative_ints.items():
            if value < 0:
                errors.append(f"{name} must be >= 0")

        nonnegative_numbers = {
            "phase1.novel_marker_min_score": self.phase1.novel_marker_min_score,
        }
        for name, value in nonnegative_numbers.items():
            if value < 0.0:
                errors.append(f"{name} must be >= 0.0")

        probabilities = {
            "phase1.host_taxonomy_deviation_overlap_threshold": self.phase1.host_taxonomy_deviation_overlap_threshold,
            "phase1.novel_marker_min_coverage": self.phase1.novel_marker_min_coverage,
            "phase2.taxonomy_ml_threshold": self.phase2.taxonomy_ml_threshold,
            "phase2.host_trim_max_host_fraction": self.phase2.host_trim_max_host_fraction,
            "phase2.host_trim_min_viral_fraction": self.phase2.host_trim_min_viral_fraction,
            "phase2.host_trim_score_threshold": self.phase2.host_trim_score_threshold,
            "phase2.host_trim_min_overlap_score": self.phase2.host_trim_min_overlap_score,
            "phase3.host_signature_evidence_threshold": self.phase3.host_signature_evidence_threshold,
            "phase3.tmvec_min_score": self.phase3.tmvec_min_score,
        }
        for name, value in probabilities.items():
            if not 0.0 <= value <= 1.0:
                errors.append(f"{name} must be between 0.0 and 1.0")

        pidents = {
            "host.high_pident_threshold": self.host.high_pident_threshold,
            "phase1.host_taxonomy_deviation_max_pident": self.phase1.host_taxonomy_deviation_max_pident,
        }
        for name, value in pidents.items():
            if not 0.0 <= value <= 100.0:
                errors.append(f"{name} must be between 0.0 and 100.0")

        # Tier threshold ranges
        if not 0.0 <= self.phase3.high_tier_threshold <= 1.0:
            errors.append("high_tier_threshold must be between 0.0 and 1.0")
        if not 0.0 <= self.phase3.low_tier_threshold <= 1.0:
            errors.append("low_tier_threshold must be between 0.0 and 1.0")
        if self.phase3.high_tier_threshold <= self.phase3.low_tier_threshold:
            errors.append("high_tier_threshold must be > low_tier_threshold")
        floor_values = [
            self.phase3.marker_floor_priority_only,
            self.phase3.marker_floor_priority_plus_family,
            self.phase3.marker_floor_priority_multi_family,
            self.phase3.marker_family_bonus_per_family,
            self.phase3.marker_multi_family_bonus,
        ]
        if any(v < 0.0 or v > 1.0 for v in floor_values):
            errors.append("phase3 marker floors/bonuses must be between 0.0 and 1.0")
        if self.phase3.marker_floor_priority_plus_family < self.phase3.marker_floor_priority_only:
            errors.append("marker_floor_priority_plus_family must be >= marker_floor_priority_only")
        if self.phase3.marker_floor_priority_multi_family < self.phase3.marker_floor_priority_plus_family:
            errors.append("marker_floor_priority_multi_family must be >= marker_floor_priority_plus_family")
        if not self.phase3.priority_marker_list:
            errors.append("phase3.priority_marker_list must include at least one marker token")

        if self.phase2.taxonomy_ml_model not in {"logreg", "gbdt", "xgboost"}:
            errors.append(
                "phase2.taxonomy_ml_model must be one of: logreg, gbdt, xgboost"
            )
        if self.phase2.taxonomy_weight_mode not in {"rank", "bitscore"}:
            errors.append("phase2.taxonomy_weight_mode must be one of: rank, bitscore")
        if self.phase3.boltz_min_seq_len > self.phase3.boltz_max_seq_len:
            errors.append("phase3.boltz_min_seq_len must be <= boltz_max_seq_len")
        if self.phase3.tmvec_require_gpu and not self.phase3.use_tmvec_database:
            errors.append(
                "phase3.tmvec_require_gpu requires phase3.use_tmvec_database=true"
            )
        if self.phase3.tmvec_databases is not None:
            unsupported = sorted(
                set(self.phase3.tmvec_databases) - {"bfvd", "cath", "swissprot", "pdb"}
            )
            if unsupported:
                errors.append(
                    f"phase3.tmvec_databases contains unsupported values: {unsupported}"
                )

        if not self.host.prefixes:
            errors.append("host.prefixes must include at least one prefix")
        elif len(set(self.host.prefixes)) != len(self.host.prefixes):
            errors.append("host.prefixes must not contain duplicates")
        elif any(not prefix.endswith("__") for prefix in self.host.prefixes):
            errors.append("host.prefixes values must end with '__'")
        if not self.host.label.strip():
            errors.append("host.label must not be empty")
        elif f"{self.host.label}__" not in self.host.prefixes:
            errors.append(
                f"host.prefixes must contain the primary prefix '{self.host.label}__'"
            )

        return errors

    def with_overrides(self, **kwargs) -> "PipelineConfig":
        """
        Return a new config with specified overrides applied.

        Maps flat kwargs (from flow function signatures) to nested config
        structure via FIELD_SPECS. Only non-None values are applied.

        Args:
            **kwargs: Flat kwargs matching flow function parameter names.

        Returns:
            New PipelineConfig with overrides applied.
        """
        import copy

        new_config = copy.deepcopy(self)
        for key, value in kwargs.items():
            if value is None:
                continue
            if key == "skip_masking":
                if type(value) is not bool:
                    raise ConfigError(
                        "Pipeline override 'skip_masking' must be a boolean"
                    )
                current = new_config.execution.masking
                new_config.execution.masking = (
                    current.with_backend(MaskingBackend.OFF)
                    if value
                    else current.with_backend(MaskingBackend.TRF_REPEATMASKER)
                    if current.backend is MaskingBackend.OFF
                    else current
                )
                continue
            spec = _WO_LOOKUP.get(key)
            if spec is None:
                raise ConfigError(f"Unknown pipeline override '{key}'")
            if key == "masking" and isinstance(value, dict):
                value = _decode_dataclass(
                    MaskingConfig,
                    value,
                    "execution.masking",
                )
            if spec.enum is not None and isinstance(value, str):
                try:
                    value = spec.enum(value)
                except ValueError as exc:
                    choices = ", ".join(repr(item.value) for item in spec.enum)
                    raise ConfigError(
                        f"Pipeline override '{key}' must be one of: {choices}"
                    ) from exc
            elif spec.path and isinstance(value, str):
                value = Path(value)
            setattr(getattr(new_config, spec.section), spec.field, value)
        errors = new_config.validate_semantics()
        if errors:
            raise ConfigError("Invalid pipeline overrides: " + "; ".join(errors))
        return new_config

    def to_flow_kwargs(self) -> dict[str, Any]:
        """Convert to kwargs dict for backward compatibility with flows."""
        out: dict[str, Any] = {}
        for spec in _EMIT_SPECS:
            value = getattr(getattr(self, spec.section), spec.field)
            if spec.enum is not None and value is not None:
                value = value.value
            out[spec.flat] = value
        return out

    @classmethod
    def from_yaml(cls, path: Path) -> "PipelineConfig":
        """Load the pipeline section of a strict application configuration."""
        from .application_config import ApplicationConfig

        return ApplicationConfig.from_yaml(path).pipeline

    @classmethod
    def from_dict(
        cls,
        data: dict,
        *,
        base_dir: Optional[Path] = None,
    ) -> "PipelineConfig":
        """Decode canonical nested pipeline configuration without side effects."""
        if not isinstance(data, dict):
            raise ConfigError("Pipeline configuration must be a mapping")
        _reject_retired_pipeline_keys(data)
        allowed = {
            "ablation",
            "databases",
            "compute",
            "host",
            "phase1",
            "phase2",
            "phase3",
            "execution",
        }
        for key in data:
            if key not in allowed:
                raise _unknown_key("", str(key), allowed)

        section_types = {
            "ablation": AblationConfig,
            "databases": DatabasePaths,
            "compute": ComputeConfig,
            "host": HostConfig,
            "phase1": Phase1Config,
            "phase2": Phase2Config,
            "phase3": Phase3Config,
            "execution": ExecutionConfig,
        }
        kwargs = {
            section: _decode_dataclass(
                section_type,
                data.get(section, {}),
                section,
                Path(base_dir) if base_dir is not None else None,
            )
            for section, section_type in section_types.items()
        }
        config = cls(**kwargs)
        errors = config.validate_semantics()
        if errors:
            raise ConfigError("Invalid pipeline configuration: " + "; ".join(errors))
        return config

    @classmethod
    def _from_flat_dict(cls, data: dict) -> "PipelineConfig":
        """Reject the retired ambiguous flat parser with a migration pointer."""
        raise ConfigError(
            "Flat pipeline configuration is no longer accepted. Use canonical "
            "ablation/databases/compute/host/phase1/phase2/phase3/execution sections."
        )

    def to_yaml(self, path: Path) -> None:
        """Write a complete v1 application file that ``from_yaml`` can read."""
        from .application_config import ApplicationConfig, OrchestrationConfig

        ApplicationConfig(
            schema_version=1,
            orchestration=OrchestrationConfig(),
            pipeline=self,
        ).to_yaml(path)
