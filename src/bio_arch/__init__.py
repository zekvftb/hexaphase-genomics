"""Biology as Information Architecture.

An open-source, reproducible framework examining biological systems
through computational and information-theoretic perspectives.
"""

__version__ = "1.0.0"

from bio_arch.contracts import (
    AnalysisRun,
    DatasetManifest,
    EvidenceClass,
    Finding,
    InterpretationRecord,
    ModuleResult,
    ValidationStatus,
)
from bio_arch.provenance import (
    SeedManager,
    compute_sha256,
    get_system_environment,
    now_iso,
)

__all__ = [
    "AnalysisRun",
    "DatasetManifest",
    "EvidenceClass",
    "Finding",
    "InterpretationRecord",
    "ModuleResult",
    "ValidationStatus",
    "SeedManager",
    "compute_sha256",
    "get_system_environment",
    "now_iso",
]
