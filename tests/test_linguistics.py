"""Unit tests for DNA Linguistics and Bytecode Profiling."""

import pytest

from bio_arch.contracts import EvidenceClass
from bio_arch.modules.linguistics import (
    analyze_linguistic_architecture,
    extract_codons,
    extract_lz_subroutines,
    fit_zipfs_law,
)


def test_extract_codons():
    """Verify codon extraction across reading frames."""
    seq = "ATGCAGTAA"
    # Frame 0: ATG, CAG, TAA
    assert extract_codons(seq, frame=0) == ["ATG", "CAG", "TAA"]
    # Frame 1: TGC, AGT
    assert extract_codons(seq, frame=1) == ["TGC", "AGT"]


def test_fit_zipfs_law():
    """Verify power-law linear regression in log-log space."""
    # Synthetic token counts following power law: rank 1 has 100, rank 2 has 50, rank 4 has 25
    counts = {"A": 100, "B": 50, "C": 33, "D": 25}
    fit = fit_zipfs_law(counts)

    assert fit.vocabulary_size == 4
    assert fit.total_tokens == 208
    # alpha should be near 1.0 (classic Zipf's law)
    assert 0.8 <= fit.alpha <= 1.2
    assert fit.r_squared > 0.95


def test_extract_lz_subroutines():
    """Verify detection of repeated modular sequence blocks."""
    # "CAG" * 10 has repeated subroutines
    seq = "ATGC" + "CAG" * 10 + "ATGC"
    subroutines = extract_lz_subroutines(seq, min_len=3)

    assert len(subroutines) > 0
    # CAG should appear with high repetition count
    assert any("CAG" in phrase for phrase in subroutines)


def test_analyze_linguistic_architecture():
    """Test full linguistic profiling integration."""
    seq = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA" * 5
    out, findings, interps = analyze_linguistic_architecture(seq, record_id="test_seq", seed=42)

    assert out["sequence_length"] == len(seq)
    assert out["codon_vocabulary_size"] > 0
    assert len(findings) >= 2
    assert len(interps) >= 1
    assert interps[0].classification == EvidenceClass.INTERPRETATION
