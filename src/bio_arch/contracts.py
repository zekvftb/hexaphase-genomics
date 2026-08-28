"""Shared data contracts and evidence classification models.

Defines typed dataclasses for data exchange between modules, provenance recording,
and scientific evidence categorization.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime
from enum import Enum
import json
import re
from typing import Any


class EvidenceClass(str, Enum):
    """Scientific evidence classification with required epistemological wording."""

    MEASUREMENT = "measurement"
    SIMULATION = "simulation"
    INTERPRETATION = "interpretation"
    HYPOTHESIS = "hypothesis"
    CROSS_VERIFICATION = "cross_verification"

    @property
    def required_prefix(self) -> str:
        """Return the required phrasing prefix for claims of this class."""
        prefixes = {
            EvidenceClass.MEASUREMENT: "The analysis measured",
            EvidenceClass.SIMULATION: "Under this model",
            EvidenceClass.INTERPRETATION: "One interpretation is",
            EvidenceClass.HYPOTHESIS: "This predicts",
            EvidenceClass.CROSS_VERIFICATION: "Cross-engine verification established that",
        }
        return prefixes[self]


class ValidationStatus(str, Enum):
    """Validation outcome for ingested datasets and artifacts."""

    VALID = "valid"
    INVALID = "invalid"
    WARNING = "warning"


_SEMVER_REGEX = re.compile(r"^\d+\.\d+\.\d+(-[0-9A-Za-z.-]+)?(\+[0-9A-Za-z.-]+)?$")
_SHA256_REGEX = re.compile(r"^[0-9a-fA-F]{64}$")


def _validate_iso_timestamp(timestamp_str: str) -> None:
    """Validate that a string conforms to ISO 8601."""
    try:
        datetime.fromisoformat(timestamp_str.replace("Z", "+00:00"))
    except (ValueError, TypeError) as exc:
        raise ValueError(f"Invalid ISO 8601 timestamp: '{timestamp_str}'") from exc


def _validate_sha256(checksum: str) -> None:
    """Validate 64-character hexadecimal SHA-256 string."""
    if not isinstance(checksum, str) or not _SHA256_REGEX.match(checksum):
        raise ValueError(f"Invalid SHA-256 checksum (expected 64 hex characters): '{checksum}'")


def _validate_semver(version: str) -> None:
    """Validate semantic version format."""
    if not isinstance(version, str) or not _SEMVER_REGEX.match(version):
        raise ValueError(f"Invalid semantic version (e.g. '0.1.0'): '{version}'")


@dataclass
class DatasetManifest:
    """Provenance and metadata contract for ingested datasets."""

    dataset_id: str
    source: str
    license: str
    retrieval_date: str
    checksum: str
    organism: str
    url: str | None = None
    assembly: str | None = None
    sequence_type: str | None = None
    annotations: dict[str, Any] = field(default_factory=dict)
    validation_status: ValidationStatus = ValidationStatus.VALID

    def __post_init__(self) -> None:
        if not self.dataset_id or not isinstance(self.dataset_id, str):
            raise ValueError("dataset_id must be a non-empty string.")
        if not self.source:
            raise ValueError("source must be specified.")
        if not self.license:
            raise ValueError("license must be specified.")
        _validate_iso_timestamp(self.retrieval_date)
        _validate_sha256(self.checksum)
        if isinstance(self.validation_status, str):
            self.validation_status = ValidationStatus(self.validation_status)

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["validation_status"] = self.validation_status.value
        return data

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> DatasetManifest:
        clean_data = dict(data)
        if "validation_status" in clean_data:
            clean_data["validation_status"] = ValidationStatus(clean_data["validation_status"])
        return cls(**clean_data)

    @classmethod
    def from_json(cls, json_str: str) -> DatasetManifest:
        return cls.from_dict(json.loads(json_str))


@dataclass
class AnalysisRun:
    """Metadata tracking a specific module analysis execution."""

    run_id: str
    timestamp: str
    module: str
    version: str
    input_ids: list[str] = field(default_factory=list)
    parameters: dict[str, Any] = field(default_factory=dict)
    seed: int | None = None
    environment: dict[str, Any] = field(default_factory=dict)
    status: str = "success"
    warnings: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if not self.run_id:
            raise ValueError("run_id must be a non-empty string.")
        _validate_iso_timestamp(self.timestamp)
        _validate_semver(self.version)
        if self.status not in ("success", "failed", "running", "warning"):
            raise ValueError(f"Invalid status: '{self.status}'. Must be success/failed/running/warning.")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> AnalysisRun:
        return cls(**data)

    @classmethod
    def from_json(cls, json_str: str) -> AnalysisRun:
        return cls.from_dict(json.loads(json_str))


@dataclass
class Finding:
    """A structured, empirical finding with observed values and controls."""

    finding_id: str
    metric: str
    observed_value: Any
    control_distribution: dict[str, Any] | None = None
    effect_size: float | None = None
    uncertainty: dict[str, float] | None = None
    adjusted_p_value: float | None = None
    biological_context: str = ""

    def __post_init__(self) -> None:
        if not self.finding_id:
            raise ValueError("finding_id must be a non-empty string.")
        if not self.metric:
            raise ValueError("metric must be specified.")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Finding:
        return cls(**data)

    @classmethod
    def from_json(cls, json_str: str) -> Finding:
        return cls.from_dict(json.loads(json_str))


@dataclass
class InterpretationRecord:
    """An explicit interpretation, simulation, or hypothesis derived from findings."""

    finding_ids: list[str]
    classification: EvidenceClass
    claim: str
    alternatives: list[str] = field(default_factory=list)
    limitations: list[str] = field(default_factory=list)
    proposed_test: str = ""

    def __post_init__(self) -> None:
        if not self.finding_ids:
            raise ValueError("finding_ids must contain at least one linked finding ID.")
        if isinstance(self.classification, str):
            self.classification = EvidenceClass(self.classification)
        if not self.claim or not isinstance(self.claim, str):
            raise ValueError("claim must be a non-empty string.")

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["classification"] = self.classification.value
        return data

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> InterpretationRecord:
        clean_data = dict(data)
        if "classification" in clean_data:
            clean_data["classification"] = EvidenceClass(clean_data["classification"])
        return cls(**clean_data)

    @classmethod
    def from_json(cls, json_str: str) -> InterpretationRecord:
        return cls.from_dict(json.loads(json_str))


@dataclass
class ModuleResult:
    """Standardized output container emitted by all modules."""

    run_metadata: AnalysisRun
    outputs: dict[str, Any] = field(default_factory=dict)
    findings: list[Finding] = field(default_factory=list)
    interpretations: list[InterpretationRecord] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    artifact_paths: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "run_metadata": self.run_metadata.to_dict(),
            "outputs": self.outputs,
            "findings": [f.to_dict() for f in self.findings],
            "interpretations": [i.to_dict() for i in self.interpretations],
            "warnings": self.warnings,
            "errors": self.errors,
            "artifact_paths": self.artifact_paths,
        }

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ModuleResult:
        run_metadata = AnalysisRun.from_dict(data["run_metadata"])
        findings = [Finding.from_dict(f) for f in data.get("findings", [])]
        interpretations = [InterpretationRecord.from_dict(i) for i in data.get("interpretations", [])]
        return cls(
            run_metadata=run_metadata,
            outputs=data.get("outputs", {}),
            findings=findings,
            interpretations=interpretations,
            warnings=data.get("warnings", []),
            errors=data.get("errors", []),
            artifact_paths=data.get("artifact_paths", []),
        )

    @classmethod
    def from_json(cls, json_str: str) -> ModuleResult:
        return cls.from_dict(json.loads(json_str))


class LogicGateType(str, Enum):
    """Biological logic gate classification."""

    FRAMESHIFT_BRANCH = "frameshift_branch"
    G4_CIRCUIT_BREAKER = "g4_circuit_breaker"
    READTHROUGH_OVERFLOW = "readthrough_overflow"
    XOR_COLLISION = "xor_collision"
    ROLLING_CIRCLE_LOOP = "rolling_circle_loop"


@dataclass
class BiologicalLogicGate:
    """Represents an identified biological logic gate / execution switch."""

    gate_id: str
    gate_type: LogicGateType
    start_pos: int
    end_pos: int
    strand: str
    trigger_motif: str
    downstream_barrier_energy: float
    predicted_efficiency: float
    target_subroutine_id: str
    description: str
    metrics: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["gate_type"] = self.gate_type.value if isinstance(self.gate_type, LogicGateType) else self.gate_type
        return data

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> BiologicalLogicGate:
        clean = dict(data)
        if "gate_type" in clean:
            clean["gate_type"] = LogicGateType(clean["gate_type"])
        return cls(**clean)


@dataclass
class LogicGateScanReport:
    """Container for biological logic gate scan outputs across a genome."""

    genome_id: str
    genome_length: int
    gates_found: list[BiologicalLogicGate] = field(default_factory=list)
    gate_counts_by_type: dict[str, int] = field(default_factory=dict)
    summary: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "genome_id": self.genome_id,
            "genome_length": self.genome_length,
            "gates_found": [g.to_dict() for g in self.gates_found],
            "gate_counts_by_type": self.gate_counts_by_type,
            "summary": self.summary,
        }

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)


@dataclass
class DualEngineEvidenceRecord:
    """Official cryptographic & mathematical congruence certificate between Python and SMC."""

    record_id: str
    genome_id: str
    sequence_length_bp: int
    python_gc_pct: float
    smc_gc_pct: float
    gc_congruence: bool
    python_phase0_codons: int
    smc_phase0_codons: int
    codon_congruence: bool
    python_execution_time_ms: float
    smc_execution_time_ms: float
    speedup_ratio: float
    sha256_checksum: str
    timestamp_iso: str
    status: str = "CERTIFIED_100_PERCENT"
    evidence_class: str = EvidenceClass.CROSS_VERIFICATION.value

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)


@dataclass
class CrisprArray:
    """Decompiled CRISPR append-only security log and repeat-spacer array."""

    array_id: str
    start_pos: int
    end_pos: int
    repeat_length_bp: int
    repeat_consensus: str
    repeats_count: int
    spacers_count: int
    spacers: list[str]
    description: str = "CRISPR Append-Only Hardware Security Array"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class CpgMemoryIsland:
    """Epigenetic 1-bit non-volatile memory register (Gardiner-Garden & Frommer criteria)."""

    island_id: str
    start_pos: int
    end_pos: int
    length_bp: int
    gc_content_pct: float
    cpg_obs_exp_ratio: float
    cpg_count: int
    is_promoter_latch: bool
    description: str = "Epigenetic Non-Volatile Memory (NVRAM) Latch"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class RiboswitchAdc:
    """RNA Analog-to-Digital chemical sensor and transcriptional switch."""

    switch_id: str
    start_pos: int
    end_pos: int
    ligand_class: str
    aptamer_motif: str
    terminator_mfe_dG: float
    switching_delta_dG: float
    predicted_state: str
    description: str = "RNA Riboswitch Analog-to-Digital Chemical Converter (ADC)"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class BiologicalCircuitReport:
    """Container for biological circuit discoveries across a genome."""

    genome_id: str
    genome_length_bp: int
    crispr_arrays: list[CrisprArray] = field(default_factory=list)
    cpg_memory_islands: list[CpgMemoryIsland] = field(default_factory=list)
    riboswitch_adcs: list[RiboswitchAdc] = field(default_factory=list)
    summary: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "genome_id": self.genome_id,
            "genome_length_bp": self.genome_length_bp,
            "crispr_arrays": [a.to_dict() for a in self.crispr_arrays],
            "cpg_memory_islands": [c.to_dict() for c in self.cpg_memory_islands],
            "riboswitch_adcs": [r.to_dict() for r in self.riboswitch_adcs],
            "summary": self.summary,
        }

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)


@dataclass
class CryptochromeRadicalPairRecord:
    """Quantum-entangled radical pair magnetoreception parameters."""

    gene_id: str
    organism: str
    tryptophan_chain_count: int
    tryptophan_positions: list[int]
    electron_hopping_pathway_valid: bool
    estimated_spin_coherence_time_ns: float
    singlet_triplet_yield_phi_s: float
    magnetoreception_viable: bool
    description: str = "Cryptochrome Radical Pair Quantum Magnetosensor"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class QuantumTunnelingRecord:
    """Löwdin quantum proton tunneling mutation vulnerability."""

    target_id: str
    codon_number: int
    codon_triplet: str
    tunneling_barrier_height_ev: float
    wkb_transmission_probability: float
    quantum_mutation_rate_per_sec: float
    is_quantum_hotspot: bool
    description: str = "Löwdin Quantum Proton Tunneling Potential"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class QuantumConductanceRecord:
    r"""1D Tight-Binding DNA $\pi$-Stack Quantum Conductance & Telemetry."""

    sequence_id: str
    length_bp: int
    guanine_trap_count: int
    mean_transmission_coefficient: float
    electronic_bandgap_ev: float
    damage_telemetry_intact: bool
    description: str = r"1D DNA \pi-Stack Quantum Telemetry Circuit"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class QuantumBiologyReport:
    """Container for quantum biological audits across genomes and sensors."""

    genome_id: str
    radical_pair_records: list[CryptochromeRadicalPairRecord] = field(default_factory=list)
    tunneling_records: list[QuantumTunnelingRecord] = field(default_factory=list)
    conductance_records: list[QuantumConductanceRecord] = field(default_factory=list)
    summary: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "genome_id": self.genome_id,
            "radical_pair_records": [r.to_dict() for r in self.radical_pair_records],
            "tunneling_records": [t.to_dict() for t in self.tunneling_records],
            "conductance_records": [c.to_dict() for c in self.conductance_records],
            "summary": self.summary,
        }

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)
