"""Unit tests for Module 5: Universal Biological Compiler & Decompiler."""

import pytest
from bio_arch.modules.universal_compiler import (
    compute_shannon_entropy,
    analyze_wobble_synchronization,
    decompile_genomic_bytecode,
    WobbleSyncReport,
    AssemblyInstruction,
)


def test_compute_shannon_entropy():
    # 4 equally probable bases = exactly 2.0 bits of entropy
    symbols = ["A", "C", "G", "T"] * 100
    h = compute_shannon_entropy(symbols)
    assert 1.99 <= h <= 2.01

    # Monomer = 0.0 bits
    assert compute_shannon_entropy(["A"] * 50) == 0.0


def test_analyze_wobble_synchronization():
    # Synthetic multi-codon sequence
    seq = "ATGCGTACGTTAGCCTAG" * 20
    report = analyze_wobble_synchronization(seq)
    assert isinstance(report, WobbleSyncReport)
    assert report.total_codons_analyzed == 120
    assert report.pos3_wobble_entropy_bits > 0.0
    assert report.pos3_information_capacity_ratio > 0.5


def test_decompile_genomic_bytecode():
    # Construct sequence with Start (ATG), Frameshift (TTTAAAC + GC stem), and Stop (TAA)
    synthetic_code = (
        "ATG"               # Start Methionine (0x0000)
        "CGT"               # Arg
        "TTTAAAC"           # Slippery heptamer
        "ACGTAC"            # Spacer
        "GCGCGCGTAAAGCGCGC" # GC Hairpin barrier
        "TAA"               # Stop Codon
    )
    instructions, asm_text = decompile_genomic_bytecode(synthetic_code, genome_id="SYNTH_VIRUS_01")
    assert len(instructions) >= 5
    assert "BIOLOGICAL MACHINE CODE DISASSEMBLY LISTING" in asm_text
    assert "START_SUBROUTINE" in asm_text
    assert "BRANCH_SLIP_MUX" in asm_text
    assert "TERM_STACK_POP" in asm_text

    # Verify first instruction is SYS_INIT_PROMOTER
    assert instructions[0].opcode == "SYS_INIT_PROMOTER"
