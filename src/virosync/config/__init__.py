"""ViroSync Pipeline Configuration."""

# New unified pipeline configuration system
from .application_config import (
    ApplicationConfig,
    FeatureResolution,
    OrchestrationConfig,
)
from .pipeline_config import (
    AblationConfig,
    AssemblyMode,
    ComputeConfig,
    ConfigError,
    DatabasePaths,
    Device,
    ExecutionConfig,
    MaskingBackend,
    MaskingConfig,
    MaskingFailurePolicy,
    Phase1Config,
    Phase2Config,
    Phase3Config,
    PipelineConfig,
    SearchBackend,
)

# Legacy threshold configuration (from thresholds.py)
from .thresholds import (
    DatabaseConfig,
    EvidenceThresholds,
    StructuralThresholds,
    ViroSyncConfig,
    get_config,
)

__all__ = [
    # New pipeline config
    "ConfigError",
    "ApplicationConfig",
    "FeatureResolution",
    "OrchestrationConfig",
    "PipelineConfig",
    "AblationConfig",
    "DatabasePaths",
    "ComputeConfig",
    "Phase1Config",
    "Phase2Config",
    "Phase3Config",
    "ExecutionConfig",
    "MaskingConfig",
    "MaskingBackend",
    "MaskingFailurePolicy",
    "AssemblyMode",
    "Device",
    "SearchBackend",
    # Legacy threshold config
    "StructuralThresholds",
    "EvidenceThresholds",
    "DatabaseConfig",
    "ViroSyncConfig",
    "get_config",
]
