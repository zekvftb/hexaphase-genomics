"""Unit tests for bio_arch.modules.cross_verification (Module 7)."""

import pytest
from bio_arch.contracts import EvidenceClass, DualEngineEvidenceRecord
from bio_arch.modules.cross_verification import (
    execute_python_pipeline,
    execute_smc_bytecode_pipeline,
    run_dual_engine_audit,
)


def test_dual_engine_congruence_synthetic():
    """Verify that Python and SMC Bytecode VM produce 100% identical outputs."""
    seq = "ATGCGTACGTTAGCTAGCTAGCTA" * 5
    record = run_dual_engine_audit(seq, genome_id="SYNTHETIC_BENCHMARK")

    assert record.gc_congruence is True
    assert record.codon_congruence is True
    assert record.status == "CERTIFIED_100_PERCENT"
    assert record.evidence_class == EvidenceClass.CROSS_VERIFICATION.value
    assert record.python_gc_pct == record.smc_gc_pct
    assert record.python_phase0_codons == record.smc_phase0_codons
    assert record.python_execution_time_ms >= 0.0
    assert record.smc_execution_time_ms >= 0.0


def test_dual_engine_empty_sequence_guard():
    """Verify that empty sequence raises ValueError."""
    with pytest.raises(ValueError, match="Sequence cannot be empty"):
        run_dual_engine_audit("", genome_id="EMPTY_TEST")
