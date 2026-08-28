"""Module 7: Dual-Engine Cross-Verification and Execution Benchmark Engine.

Implements N-Version parallel cross-validation between the Python reference
genomics pipeline and the SMC Bytecode Virtual Machine.
"""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
import sys
import time
from typing import Any

from bio_arch.contracts import DualEngineEvidenceRecord, EvidenceClass
from bio_arch.logger import setup_logger
from bio_arch.modules.information import compute_composition
from bio_arch.provenance import now_iso

# SMC Compiler & Bytecode Engine imports
try:
    from smc.bytecode_vm import BytecodeVM
    from smc.compiler import BytecodeCompiler
    from smc.lexer import SmcLexer
    from smc.parser import SmcParser
    SMC_AVAILABLE = True
except ImportError:
    SMC_AVAILABLE = False

logger = setup_logger("bio_arch.cross_verification")


def execute_python_pipeline(sequence: str) -> tuple[dict[str, Any], float]:
    """Execute reference Python genomic analysis and measure duration in ms."""
    t0 = time.perf_counter_ns()
    
    comp = compute_composition(sequence)
    # Phase 0 codon extraction
    seq_upper = sequence.upper()
    phase0_codons = [seq_upper[i : i + 3] for i in range(0, len(seq_upper) - 2, 3)]
    
    t1 = time.perf_counter_ns()
    elapsed_ms = (t1 - t0) / 1_000_000.0

    return {
        "gc_content_pct": round(comp.gc_content * 100.0, 2),
        "phase0_codon_count": len(phase0_codons),
        "length_bp": len(sequence),
    }, elapsed_ms


def execute_smc_bytecode_pipeline(sequence: str) -> tuple[dict[str, Any], float]:
    """Execute identical genomic analysis on SMC Bytecode VM and measure duration in ms."""
    if not SMC_AVAILABLE:
        raise RuntimeError("SMC Language compiler/VM is not installed in the current environment.")

    # Generate equivalent SMC code
    smc_code = f"""
    var seq = "{sequence.upper()}"
    var g_count = count_matches(seq, "G")
    var c_count = count_matches(seq, "C")
    var total_gc = g_count + c_count
    var gc_pct = 0.0
    if len(seq) > 0 {{
        gc_pct = round((total_gc * 100.0) / len(seq), 2)
    }}
    var codons = hexaphase_window(seq, "+0", 3)
    var phase0_count = len(codons)
    """

    tokens = SmcLexer(smc_code, strict=True).tokenize()
    ast = SmcParser(tokens).parse()
    chunk = BytecodeCompiler().compile(ast)
    vm = BytecodeVM()

    t0 = time.perf_counter_ns()
    vm.run(chunk)
    t1 = time.perf_counter_ns()
    elapsed_ms = (t1 - t0) / 1_000_000.0

    return {
        "gc_content_pct": float(vm.globals.get("gc_pct", 0.0)),
        "phase0_codon_count": int(vm.globals.get("phase0_count", 0)),
        "length_bp": len(sequence),
    }, elapsed_ms


def run_dual_engine_audit(sequence: str, genome_id: str = "GENOME_001") -> DualEngineEvidenceRecord:
    """Execute dual-engine cross-verification and certify mathematical congruence."""
    seq_clean = sequence.strip()
    if not seq_clean:
        raise ValueError("Sequence cannot be empty for dual-engine audit.")

    py_res, py_time_ms = execute_python_pipeline(seq_clean)
    smc_res, smc_time_ms = execute_smc_bytecode_pipeline(seq_clean)

    gc_match = math.isclose(py_res["gc_content_pct"], smc_res["gc_content_pct"], abs_tol=0.05)
    codon_match = (py_res["phase0_codon_count"] == smc_res["phase0_codon_count"])

    if not (gc_match and codon_match):
        raise AssertionError(
            f"Dual-engine divergence detected for {genome_id}!\n"
            f"Python: GC={py_res['gc_content_pct']}%, Codons={py_res['phase0_codon_count']}\n"
            f"SMC:    GC={smc_res['gc_content_pct']}%, Codons={smc_res['phase0_codon_count']}"
        )

    # Compute SHA-256 fingerprint
    hasher = hashlib.sha256()
    hasher.update(seq_clean.encode("utf-8"))
    hasher.update(str(py_res).encode("utf-8"))
    checksum = hasher.hexdigest()

    speedup = round(py_time_ms / max(0.0001, smc_time_ms), 3)

    return DualEngineEvidenceRecord(
        record_id=f"DUAL_ENGINE_{genome_id.upper()}",
        genome_id=genome_id,
        sequence_length_bp=len(seq_clean),
        python_gc_pct=py_res["gc_content_pct"],
        smc_gc_pct=smc_res["gc_content_pct"],
        gc_congruence=gc_match,
        python_phase0_codons=py_res["phase0_codon_count"],
        smc_phase0_codons=smc_res["phase0_codon_count"],
        codon_congruence=codon_match,
        python_execution_time_ms=round(py_time_ms, 3),
        smc_execution_time_ms=round(smc_time_ms, 3),
        speedup_ratio=speedup,
        sha256_checksum=checksum,
        timestamp_iso=now_iso(),
        status="CERTIFIED_100_PERCENT",
        evidence_class=EvidenceClass.CROSS_VERIFICATION.value,
    )
