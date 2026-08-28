"""Unit tests for bio_arch.modules.logic_gates."""

import pytest
from bio_arch.contracts import LogicGateType, LogicGateScanReport
from bio_arch.modules.logic_gates import (
    calculate_stem_loop_mfe,
    scan_frameshift_branches,
    scan_g_quadruplexes,
    scan_readthrough_gates,
    scan_xor_collisions,
    scan_all_logic_gates,
)


def test_calculate_stem_loop_mfe():
    """Verify minimum free energy estimation for a GC-rich hairpin."""
    # Stable stem: GGGCCC (stem) AAAA (loop) GGGCCC (complement)
    hairpin = "GGGCCCAAAAGGGCCC"
    dG, stem_len, loop_len = calculate_stem_loop_mfe(hairpin)
    assert dG < -8.0, f"Expected strong negative dG, got {dG}"
    assert stem_len >= 4

    # Unstructured poly-A
    poly_a = "AAAAAAAAAAAAAAAAAAAA"
    dG_null, _, _ = calculate_stem_loop_mfe(poly_a)
    assert dG_null == 0.0


def test_scan_frameshift_branches_positive():
    """Verify detection of a programmed -1 ribosomal frameshift multiplexer."""
    # Construct synthetic coronavirus-like slippery site: UUUAAAC + 6nt spacer + GC hairpin
    synthetic_fs = (
        "ATGCGTACGTTTAAAC"  # contains TTTAAAC (start at 9)
        "ACGTAC"            # 6nt spacer
        "GCGCGCGTAAAGCGCGC" # Stable GC stem loop (dG < -10)
        "TAGCTA"
    )
    gates = scan_frameshift_branches(synthetic_fs, min_barrier_energy=-7.0)
    assert len(gates) >= 1
    g = gates[0]
    assert g.gate_type == LogicGateType.FRAMESHIFT_BRANCH
    assert "TTTAAAC" in g.trigger_motif
    assert g.downstream_barrier_energy <= -7.0
    assert 0.05 <= g.predicted_efficiency <= 0.85


def test_scan_g_quadruplexes_positive():
    """Verify detection of forward G4 and reverse i-Motif circuit breakers."""
    seq_with_g4 = "NNNN" + "GGGAATGGGAATGGGAATGGG" + "NNNN" + "CCCAATCCCAATCCCAATCCC" + "NNNN"
    gates = scan_g_quadruplexes(seq_with_g4, min_hunter_score=0.8)
    
    types = [g.gate_type for g in gates]
    assert len(gates) == 2
    assert all(t == LogicGateType.G4_CIRCUIT_BREAKER for t in types)
    
    strands = [g.strand for g in gates]
    assert "+" in strands
    assert "-" in strands


def test_scan_readthrough_gates():
    """Verify detection of leaky stop-codon bypass gates."""
    seq_readthrough = "ATGCGATCG" + "TGACAA" + "GCGCGCGAAAAACGCGCGC" + "TAA"
    gates = scan_readthrough_gates(seq_readthrough)
    assert len(gates) >= 1
    g = gates[0]
    assert g.gate_type == LogicGateType.READTHROUGH_OVERFLOW
    assert g.trigger_motif.startswith("TGAC")
    assert g.predicted_efficiency > 0.05


def test_scan_xor_collisions():
    """Verify detection of convergent overlapping gene collision XOR switches."""
    genes = [
        {"name": "Gene_Sense", "start": 100, "end": 600, "strand": "+"},
        {"name": "Gene_Antisense", "start": 450, "end": 900, "strand": "-"},
    ]
    gates = scan_xor_collisions(genes)
    assert len(gates) == 1
    g = gates[0]
    assert g.gate_type == LogicGateType.XOR_COLLISION
    assert g.start_pos == 450
    assert g.end_pos == 600
    assert g.metrics["overlap_length_bp"] == 151


def test_scan_all_logic_gates_e2e():
    """Verify end-to-end multi-gate scan report generation and JSON serialization."""
    test_genome = (
        "ATGCGTACGTTTAAACACGTACGCGCGCGTAAAGCGCGCTAGCTA"
        "NNNNGGGAATGGGAATGGGAATGGGNNNN"
        "TGACAAGCGCGCGAAAAACGCGCGCTAA"
    )
    report = scan_all_logic_gates(test_genome, genome_id="TEST_VIRUS_01")
    assert isinstance(report, LogicGateScanReport)
    assert report.genome_id == "TEST_VIRUS_01"
    assert len(report.gates_found) >= 3
    assert report.summary["total_logic_gates_detected"] >= 3
    
    # Test JSON roundtrip
    json_str = report.to_json()
    assert "TEST_VIRUS_01" in json_str
    assert "gates_found" in json_str
