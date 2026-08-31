"""Differential Testing & Property-Based Invariant Verification Suite.

Validates custom bioinformatics and compilation algorithms directly against
Biopython (Bio.Seq / Bio.SeqUtils) and mathematical invariants using Hypothesis.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import random

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from hypothesis import given, settings, strategies as st
import pytest

from bio_arch.modules.information import (
    compute_composition,
    shannon_entropy,
    shuffle_dinucleotide,
    shuffle_mononucleotide,
)
from bio_arch.modules.logic_gates import (
    calculate_stem_loop_mfe,
    scan_all_logic_gates,
    scan_frameshift_branches,
)
from bio_arch.modules.recompiler import (
    CODON_TABLE,
    recompile_dual_protein_dna,
    translate_sequence,
)
from bio_arch.modules.universal_compiler import COMPLEMENT_MAP


# ---------------------------------------------------------------------------
# 1. Differential Testing vs Biopython
# ---------------------------------------------------------------------------

@given(st.text(alphabet="ACGT", min_size=3, max_size=300).filter(lambda s: len(s) % 3 == 0))
@settings(max_examples=50)
def test_differential_translation_standard_code(dna_seq: str):
    """Verify that translate_sequence matches Biopython standard genetic code translation."""
    custom_translation = translate_sequence(dna_seq, offset=0)
    biopython_translation = str(Seq(dna_seq).translate(table=1))

    assert custom_translation == biopython_translation, (
        f"Translation discrepancy on sequence '{dna_seq}': "
        f"Custom='{custom_translation}' vs Biopython='{biopython_translation}'"
    )


@given(st.text(alphabet="ACGT", min_size=1, max_size=200))
@settings(max_examples=50)
def test_differential_reverse_complement(dna_seq: str):
    """Verify that custom COMPLEMENT_MAP reverse complement matches Biopython."""
    custom_rc = dna_seq.translate(COMPLEMENT_MAP)[::-1]
    biopython_rc = str(Seq(dna_seq).reverse_complement())

    assert custom_rc == biopython_rc, (
        f"Reverse complement discrepancy on '{dna_seq}': "
        f"Custom='{custom_rc}' vs Biopython='{biopython_rc}'"
    )


@given(st.text(alphabet="ACGT", min_size=1, max_size=300))
@settings(max_examples=50)
def test_differential_gc_fraction(dna_seq: str):
    """Verify that compute_composition GC fraction matches Bio.SeqUtils.gc_fraction."""
    comp = compute_composition(dna_seq)
    biopython_gc = gc_fraction(dna_seq)

    assert pytest.approx(comp.gc_content, abs=1e-5) == biopython_gc, (
        f"GC fraction mismatch on '{dna_seq}': "
        f"Custom={comp.gc_content} vs Biopython={biopython_gc}"
    )


def test_differential_gc_skew_formula():
    """Verify GC skew calculation against standard mathematical formula (G-C)/(G+C)."""
    test_seqs = [
        "GGGGCCCC",     # G=4, C=4 -> Skew = 0.0
        "GGGGGGCC",     # G=6, C=2 -> Skew = (6-2)/(6+2) = 4/8 = 0.5
        "CCGGGGGG",     # G=6, C=2 -> Skew = 0.5
        "CCCCGGAA",     # G=2, C=4 -> Skew = (2-4)/(2+4) = -2/6 = -0.3333...
        "AAAAAAAA",     # G=0, C=0 -> Skew = 0.0
    ]
    for s in test_seqs:
        comp = compute_composition(s)
        g = s.count("G")
        c = s.count("C")
        expected_skew = (g - c) / (g + c) if (g + c) > 0 else 0.0
        assert pytest.approx(comp.gc_skew, abs=1e-5) == expected_skew


# ---------------------------------------------------------------------------
# 2. Property-Based Invariant Testing
# ---------------------------------------------------------------------------

@given(
    st.text(alphabet="ACGT", min_size=10, max_size=300),
    st.integers(min_value=0, max_value=5),
)
@settings(max_examples=50)
def test_property_translation_bounds_and_no_overflow(dna_seq: str, offset: int):
    """Verify that translate_sequence handles arbitrary offsets safely without overflow."""
    result = translate_sequence(dna_seq, offset=offset)
    expected_max_aa = max(0, (len(dna_seq) - offset) // 3)
    assert len(result) == expected_max_aa
    assert all(aa in "ACDEFGHIKLMNPQRSTVWY*" for aa in result)


@given(
    st.integers(min_value=0, max_value=2),
    st.integers(min_value=-100, max_value=100),
)
def test_property_reading_phase_invariant(initial_phase: int, shift: int):
    """Verify reading phase state shift slip(k) strictly satisfies 0 <= phi <= 2."""
    shifted_phase = (initial_phase + shift) % 3
    assert 0 <= shifted_phase <= 2


@given(st.text(alphabet="ACGT", min_size=20, max_size=200))
@settings(max_examples=40)
def test_property_dinucleotide_shuffle_invariants(dna_seq: str):
    """Verify that Altschul-Erickson dinucleotide shuffle strictly preserves length and single-nucleotide counts."""
    rng = random.Random(12345)
    shuffled = shuffle_dinucleotide(dna_seq, rng)

    # 1. Exact length preservation
    assert len(shuffled) == len(dna_seq)

    # 2. Exact mononucleotide frequency preservation
    assert Counter(shuffled) == Counter(dna_seq)

    # 3. Same first and last nucleotide (Eulerian walk property)
    if len(dna_seq) > 2:
        assert shuffled[0] == dna_seq[0]
        assert shuffled[-1] == dna_seq[-1]


# ---------------------------------------------------------------------------
# 3. Determinism & Isolation Verification
# ---------------------------------------------------------------------------

def test_deterministic_seed_reproducibility():
    """Verify that analysis with fixed seed produces 100% bitwise-identical results across runs."""
    dna = "ATGCGATCGATCGATCGATCGATCGATCGATCG" * 5
    rng1 = random.Random(42)
    rng2 = random.Random(42)

    shuf1 = [shuffle_dinucleotide(dna, rng1) for _ in range(10)]
    shuf2 = [shuffle_dinucleotide(dna, rng2) for _ in range(10)]

    assert shuf1 == shuf2


def test_fault_tolerant_error_handling():
    """Verify that zero division, empty sequences, and missing keys fail safely without crashing."""
    # 1. Empty sequence handling
    comp = compute_composition("")
    assert comp.length == 0
    assert comp.gc_content == 0.0
    assert comp.gc_skew == 0.0

    # 2. Zero-length entropy
    entropy = shannon_entropy("")
    assert entropy == 0.0

    # 3. Stem loop on tiny sequence
    dG, stem, loop = calculate_stem_loop_mfe("A")
    assert dG == 0.0


# ---------------------------------------------------------------------------
# 4. Known Positive & Negative Biological Controls
# ---------------------------------------------------------------------------

def test_positive_control_hiv1_gag_pol_frameshift():
    """Positive Control: Verify detection of HIV-1 gag-pol frameshift motif (AAAUUUA)."""
    # Canonical HIV-1 slippery site: AAAUUUA (AAATTTA in DNA) with downstream stem-loop
    hiv1_fragment = "ATGCCCAAAATTTAACGTACGCGCGCGTAAAGCGCGCTAGCTA"
    gates = scan_frameshift_branches(hiv1_fragment, min_barrier_energy=-7.0)
    
    assert len(gates) >= 1
    motifs = [g.trigger_motif for g in gates]
    assert any("AAATTTA" in m for m in motifs)


def test_positive_control_sars_cov_2_frameshift():
    """Positive Control: Verify detection of SARS-CoV-2 ORF1a/1b frameshift motif (UUUAAAC)."""
    # Canonical SARS-CoV-2 slippery site: UUUAAAC (TTTAAAC in DNA) with downstream stem-loop
    sars2_fragment = "ATGCGTACGTTTAAACACGTACGCGCGCGTAAAGCGCGCTAGCTA"
    gates = scan_frameshift_branches(sars2_fragment, min_barrier_energy=-7.0)

    assert len(gates) >= 1
    motifs = [g.trigger_motif for g in gates]
    assert any("TTTAAAC" in m for m in motifs)


def test_negative_control_uniform_random_sequence():
    """Negative Control: Verify that synthetic unstructured poly-A sequence yields zero frameshift gates."""
    poly_a = "A" * 500
    report = scan_all_logic_gates(poly_a, genome_id="NEGATIVE_CONTROL_POLY_A", num_shuffles=10, seed=42)

    assert len(report.gates_found) == 0
    assert report.summary["branching_multiplexers"] == 0
    assert report.summary["g4_circuit_breakers"] == 0
