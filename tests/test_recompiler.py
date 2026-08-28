"""Unit tests for Module 6: Dual-Phase Biological Recompiler."""

import pytest
from bio_arch.modules.recompiler import (
    recompile_dual_protein_dna,
    translate_sequence,
    RecompilationResult,
)


def test_recompile_dual_protein_exact():
    # Test compatible peptide pair (M + C + L in Frame 0 and C + L in Frame +1)
    p0 = "MCLV"
    p1 = "CLWV"
    res = recompile_dual_protein_dna(p0, p1)
    assert isinstance(res, RecompilationResult)
    assert len(res.synthesized_dna) == 13  # 4 codons * 3 + 1
    assert res.f0_identity_pct == 100.0
    assert res.compression_ratio > 1.8


def test_translate_sequence_offsets():
    dna = "ATGCGATCGTAA"
    f0 = translate_sequence(dna, offset=0)
    f1 = translate_sequence(dna, offset=1)
    assert f0.startswith("M")
    assert len(f1) >= 3
