"""Unit tests for shared data contracts and validation."""

import json
import pytest

from bio_arch.contracts import (
    AnalysisRun,
    DatasetManifest,
    EvidenceClass,
    Finding,
    InterpretationRecord,
    ModuleResult,
    ValidationStatus,
)


def test_evidence_class_prefixes():
    """Verify evidence classes provide their standardized epistemological prefixes."""
    assert EvidenceClass.MEASUREMENT.required_prefix == "The analysis measured"
    assert EvidenceClass.SIMULATION.required_prefix == "Under this model"
    assert EvidenceClass.INTERPRETATION.required_prefix == "One interpretation is"
    assert EvidenceClass.HYPOTHESIS.required_prefix == "This predicts"


def test_dataset_manifest_roundtrip():
    """Test DatasetManifest creation, validation, and JSON serialization."""
    manifest = DatasetManifest(
        dataset_id="ds_001_synthetic",
        source="Synthetic Control Generator",
        license="MIT",
        retrieval_date="2026-08-27T12:00:00+00:00",
        checksum="e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
        organism="Synthetic Lambda",
        url="https://example.org/datasets/synth.fasta",
        sequence_type="DNA",
        annotations={"feature_count": 12},
        validation_status=ValidationStatus.VALID,
    )

    data = manifest.to_dict()
    assert data["dataset_id"] == "ds_001_synthetic"
    assert data["validation_status"] == "valid"

    json_str = manifest.to_json()
    reconstructed = DatasetManifest.from_json(json_str)
    assert reconstructed.dataset_id == manifest.dataset_id
    assert reconstructed.checksum == manifest.checksum
    assert reconstructed.validation_status == ValidationStatus.VALID


def test_dataset_manifest_validation_failures():
    """Verify DatasetManifest catches invalid checksums and timestamps."""
    valid_sha = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"

    # Invalid SHA-256 (wrong length or characters)
    with pytest.raises(ValueError, match="Invalid SHA-256 checksum"):
        DatasetManifest(
            dataset_id="ds_bad",
            source="Test",
            license="MIT",
            retrieval_date="2026-08-27T12:00:00+00:00",
            checksum="short_hash_123",
            organism="Test Org",
        )

    # Invalid timestamp
    with pytest.raises(ValueError, match="Invalid ISO 8601 timestamp"):
        DatasetManifest(
            dataset_id="ds_bad",
            source="Test",
            license="MIT",
            retrieval_date="not-a-timestamp",
            checksum=valid_sha,
            organism="Test Org",
        )

    # Empty dataset_id
    with pytest.raises(ValueError, match="dataset_id must be a non-empty string"):
        DatasetManifest(
            dataset_id="",
            source="Test",
            license="MIT",
            retrieval_date="2026-08-27T12:00:00+00:00",
            checksum=valid_sha,
            organism="Test Org",
        )


def test_analysis_run_roundtrip():
    """Test AnalysisRun creation and serialization."""
    run = AnalysisRun(
        run_id="run_20260827_001",
        timestamp="2026-08-27T14:30:00+00:00",
        module="module_1_information",
        version="0.1.0",
        input_ids=["ds_001_synthetic"],
        parameters={"k_mer": 3, "shuffles": 100},
        seed=42,
        environment={"python": "3.11"},
        status="success",
        warnings=[],
    )

    json_str = run.to_json()
    reconstructed = AnalysisRun.from_json(json_str)
    assert reconstructed.run_id == run.run_id
    assert reconstructed.seed == 42
    assert reconstructed.parameters["k_mer"] == 3


def test_analysis_run_validation():
    """Test AnalysisRun semver and status validation."""
    with pytest.raises(ValueError, match="Invalid semantic version"):
        AnalysisRun(
            run_id="run_bad",
            timestamp="2026-08-27T14:30:00+00:00",
            module="mod1",
            version="version_1",  # Invalid semver
        )

    with pytest.raises(ValueError, match="Invalid status"):
        AnalysisRun(
            run_id="run_bad",
            timestamp="2026-08-27T14:30:00+00:00",
            module="mod1",
            version="1.0.0",
            status="aborted",  # Not an allowed status
        )


def test_finding_and_interpretation_linkage():
    """Verify findings and interpretations can be linked and serialized properly."""
    finding = Finding(
        finding_id="f_001_shannon_entropy",
        metric="shannon_entropy_k1",
        observed_value=1.982,
        control_distribution={"mean": 1.995, "std": 0.002, "null_model": "dinucleotide_shuffle"},
        effect_size=-0.013,
        uncertainty={"ci_lower": -0.015, "ci_upper": -0.011},
        adjusted_p_value=0.0004,
        biological_context="Promoter region of lac operon",
    )

    interp = InterpretationRecord(
        finding_ids=[finding.finding_id],
        classification=EvidenceClass.INTERPRETATION,
        claim="One interpretation is that the sequence constraint reduces local Shannon entropy.",
        alternatives=["Biased GC content naturally reduces entropy under a random background model."],
        limitations=["Analysis restricted to single-nucleotide scale."],
        proposed_test="Test with higher k-mer orders (k=2 to 4) against dinucleotide-preserving nulls.",
    )

    run_meta = AnalysisRun(
        run_id="run_test",
        timestamp="2026-08-27T15:00:00+00:00",
        module="module_1_information",
        version="0.1.0",
    )

    mod_result = ModuleResult(
        run_metadata=run_meta,
        outputs={"mean_entropy": 1.982},
        findings=[finding],
        interpretations=[interp],
        artifact_paths=["outputs/runs/test/entropy.json"],
    )

    json_str = mod_result.to_json()
    parsed = json.loads(json_str)
    assert len(parsed["findings"]) == 1
    assert parsed["findings"][0]["metric"] == "shannon_entropy_k1"
    assert parsed["interpretations"][0]["classification"] == "interpretation"

    restored = ModuleResult.from_json(json_str)
    assert len(restored.findings) == 1
    assert restored.findings[0].observed_value == 1.982
    assert restored.interpretations[0].classification == EvidenceClass.INTERPRETATION
