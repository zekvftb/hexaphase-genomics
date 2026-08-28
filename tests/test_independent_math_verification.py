"""Independent Mathematical Sanity Check and Ground-Truth Verification.

Zero-dependency verification suite testing core mathematical bounds,
known theoretical limits, and NCBI canonical translation tables.
"""

import math
import pytest
from bio_arch.modules.information import shannon_entropy
from bio_arch.modules.recompiler import translate_sequence, recompile_dual_protein_dna


def test_independent_shannon_entropy_theoretical_bounds():
    """Verify Shannon entropy calculation against known analytical mathematical limits."""
    
    # 1. Zero Entropy Limit (Single character = 0.0 bits)
    zero_seq = "AAAAAAAAAAAAAAAAAAAA"
    h_zero = shannon_entropy(zero_seq, k=1)
    assert h_zero == 0.0, f"Expected 0.0 bits for pure sequence, got {h_zero}"

    # 2. 1-Bit Entropy Limit (Equal 50/50 2-character split = 1.0 bit)
    binary_seq = "AC" * 50
    h_binary = shannon_entropy(binary_seq, k=1)
    assert math.isclose(h_binary, 1.0, rel_tol=1e-5), f"Expected 1.0 bit for 50/50 split, got {h_binary}"

    # 3. Maximum Theoretical 4-Base Limit (Equal 25% distribution = log2(4) = 2.0 bits)
    max_seq = "ACGT" * 25
    h_max = shannon_entropy(max_seq, k=1)
    assert math.isclose(h_max, 2.0, rel_tol=1e-5), f"Expected 2.0 bits for 4-base uniform, got {h_max}"


def test_independent_ncbi_standard_genetic_code():
    """Verify that canonical start, stop, and standard codons match NCBI Genetic Code Table 1."""
    # Standard start codon
    assert translate_sequence("ATG", offset=0) == "M"
    # Standard stop codons
    assert translate_sequence("TAA", offset=0) == "*"
    assert translate_sequence("TAG", offset=0) == "*"
    assert translate_sequence("TGA", offset=0) == "*"
    # Distinct amino acids
    assert translate_sequence("TGG", offset=0) == "W"  # Tryptophan
    assert translate_sequence("TTT", offset=0) == "F"  # Phenylalanine
    assert translate_sequence("CCG", offset=0) == "P"  # Proline


def test_independent_recompiler_dual_channel_exactness():
    """Verify that dual-phase compilation produces 100% exact primary frame translation."""
    p0 = "MCLV"
    p1 = "CLWV"
    res = recompile_dual_protein_dna(p0, p1)
    
    # Independently re-translate the synthesized physical DNA
    independent_f0 = translate_sequence(res.synthesized_dna, offset=0)
    assert independent_f0 == p0, f"Primary frame mismatch: expected {p0}, got {independent_f0}"
    assert res.f0_identity_pct == 100.0
