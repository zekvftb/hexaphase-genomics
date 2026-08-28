"""Scaffold and contracts demonstration.

Demonstrates creating typed manifests, analysis runs, findings, and interpretation records
completely offline without external services or paid LLMs.
"""

from pathlib import Path
import sys

# Add src to sys.path for direct execution without prior pip install
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from bio_arch.contracts import (
    AnalysisRun,
    DatasetManifest,
    EvidenceClass,
    Finding,
    InterpretationRecord,
    ModuleResult,
    ValidationStatus,
)
from bio_arch.logger import setup_logger
from bio_arch.provenance import (
    SeedManager,
    compute_sha256,
    get_system_environment,
    now_iso,
)


def main() -> None:
    logger = setup_logger("scaffold_demo")
    logger.info("Initializing Biology as Information Architecture Scaffold Demo...")

    # 1. Simulate an in-memory synthetic sequence dataset
    synthetic_fasta_content = (
        ">synth_seq_1 description=lac_promoter_model\n"
        "AATTGTGAGCGGATAACAATTTCACACAGGAAACAGCTATGACCATGATTACGGATTCACT\n"
    ).encode("utf-8")

    # Compute exact checksum
    checksum = compute_sha256(synthetic_fasta_content)
    logger.info("Computed SHA-256 for synthetic input: %s", checksum)

    # 2. Build DatasetManifest
    manifest = DatasetManifest(
        dataset_id="ds_synth_001",
        source="In-silico Synthetic Benchmark",
        license="MIT",
        retrieval_date=now_iso(),
        checksum=checksum,
        organism="Escherichia coli model fragment",
        sequence_type="DNA",
        annotations={"length_bp": 60, "region": "promoter_operator"},
        validation_status=ValidationStatus.VALID,
    )
    logger.info("Created DatasetManifest: %s", manifest.dataset_id)

    # 3. Seed and Provenance Management
    seed_mgr = SeedManager(master_seed=42)
    module_seed = seed_mgr.derive_seed("module_1_information")
    env_metadata = get_system_environment()

    # 4. Record AnalysisRun
    run_meta = AnalysisRun(
        run_id="run_demo_20260827",
        timestamp=now_iso(),
        module="module_1_information",
        version="0.1.0",
        input_ids=[manifest.dataset_id],
        parameters={"k_mer": 2, "null_model": "dinucleotide_shuffle", "shuffles": 100},
        seed=module_seed,
        environment=env_metadata,
        status="success",
        warnings=[],
    )

    # 5. Create a Finding (Measurement)
    finding = Finding(
        finding_id="f_001_gc_skew",
        metric="gc_skew",
        observed_value=0.0833,
        control_distribution={"mean": 0.001, "std": 0.045, "null_model": "dinucleotide_shuffle"},
        effect_size=1.82,
        uncertainty={"ci_95_lower": 0.052, "ci_95_upper": 0.114},
        adjusted_p_value=0.034,
        biological_context="Lac operator region sequence asymmetry",
    )

    # 6. Create an InterpretationRecord (strictly distinguishing observation from claim)
    interpretation = InterpretationRecord(
        finding_ids=[finding.finding_id],
        classification=EvidenceClass.INTERPRETATION,
        claim=f"{EvidenceClass.INTERPRETATION.required_prefix} that local base skew reflects transcriptional directionality.",
        alternatives=["Biased mutational pressure during leading vs lagging strand replication."],
        limitations=["Small sample length (60 bp); statistical confidence requires wider genomic windows."],
        proposed_test="Benchmark against 1000 bp flanking regions and replication origin distance.",
    )

    # 7. Aggregate into ModuleResult
    result = ModuleResult(
        run_metadata=run_meta,
        outputs={"sequence_length": 60, "gc_content": 0.40},
        findings=[finding],
        interpretations=[interpretation],
        artifact_paths=["outputs/demo_finding.json"],
    )

    logger.info("ModuleResult constructed successfully.")
    print("\n--- Standardized ModuleResult JSON Output ---")
    print(result.to_json(indent=2))
    print("---------------------------------------------\n")
    logger.info("Demo complete. 100% offline, reproducible, and zero external costs.")


if __name__ == "__main__":
    main()
