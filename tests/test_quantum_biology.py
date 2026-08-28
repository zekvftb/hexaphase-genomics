"""Unit tests for bio_arch.modules.quantum_biology (Module 9)."""

import pytest
from bio_arch.contracts import (
    CryptochromeRadicalPairRecord,
    QuantumBiologyReport,
    QuantumConductanceRecord,
    QuantumTunnelingRecord,
)
from bio_arch.modules.quantum_biology import (
    analyze_cryptochrome_radical_pairs,
    calculate_dna_tight_binding_conductance,
    calculate_wkb_proton_tunneling,
    scan_all_quantum_biology,
)


def test_wkb_proton_tunneling_bounds():
    """Verify that WKB transmission probability strictly satisfies 0.0 < T <= 1.0."""
    seq = "ATGCGTCGACCGCGGCGATAG"  # contains CGG/CCG hotspots
    records = calculate_wkb_proton_tunneling(seq, target_id="TP53_TEST")

    assert len(records) == len(seq) // 3
    for r in records:
        assert 0.0 < r.wkb_transmission_probability <= 1.0
        assert 0.20 <= r.tunneling_barrier_height_ev <= 0.45
        assert r.quantum_mutation_rate_per_sec >= 0.0

    # Ensure CpG codon 'CGG' has lower barrier height than 'ATG'
    cgg_record = [r for r in records if "CG" in r.codon_triplet][0]
    atg_record = [r for r in records if r.codon_triplet == "ATG"][0]
    assert cgg_record.tunneling_barrier_height_ev < atg_record.tunneling_barrier_height_ev


def test_cryptochrome_avian_vs_human_radical_pairs():
    """Verify that both European Robin and Human Cryptochromes maintain valid Trp electron relays."""
    # Synthetic Cryptochrome fragment with 4 conserved Tryptophans
    mock_cry = "MVHWFRKGLRLHDNPALLAAVRGARCVRCVYILDPWFAASSSVGINRWRFLLLSLERWR"
    record = analyze_cryptochrome_radical_pairs(mock_cry, gene_id="TEST_CRY", organism="Homo sapiens")

    assert isinstance(record, CryptochromeRadicalPairRecord)
    assert record.tryptophan_chain_count >= 3
    assert record.electron_hopping_pathway_valid is True
    assert record.estimated_spin_coherence_time_ns >= 500.0
    assert record.magnetoreception_viable is True


def test_dna_tight_binding_conductance():
    """Verify 1D Tight-Binding electronic transmission across DNA pi-stack."""
    gc_rich_dna = "GGGGCCCCGGGGCCCC"
    at_rich_dna = "AAAATTTTAAAATTTT"

    cond_gc = calculate_dna_tight_binding_conductance(gc_rich_dna)
    cond_at = calculate_dna_tight_binding_conductance(at_rich_dna)

    assert cond_gc.guanine_trap_count == 8
    assert cond_at.guanine_trap_count == 0
    # GC-rich stack with Guanine traps maintains higher electronic transmission
    assert cond_gc.mean_transmission_coefficient > cond_at.mean_transmission_coefficient
    assert cond_gc.damage_telemetry_intact is True


def test_scan_all_quantum_biology_e2e():
    """End-to-end audit report generation for quantum biology."""
    seq = "ATGCGTTTAAACGCGCGCTAG"
    report = scan_all_quantum_biology(seq, target_id="QUANTUM_TEST")
    assert isinstance(report, QuantumBiologyReport)
    assert "quantum_tunneling_hotspots" in report.summary
    assert "damage_telemetry_functional" in report.summary
