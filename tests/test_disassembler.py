"""Unit tests for Biological Disassembler and Decompiler."""

import pytest

from bio_arch.modules.disassembler import (
    decompile_to_pseudocode,
    disassemble_sequence,
    generate_assembly_listing,
    scan_tokens,
)


def test_scan_tokens_promoter_and_rbs():
    """Verify detection of Pribnow box and Shine-Dalgarno tokens."""
    # Contains Pribnow box TATAAT and RBS AGGAGG
    seq = "GGCTATAATCCAGGAGGTTCATGAAACCCTAA"
    tokens = scan_tokens(seq)

    types = [t.token_type for t in tokens]
    assert "PROMOTER_MINUS10" in types
    assert "RBS_SHINE_DALGARNO" in types


def test_scan_tokens_orf():
    """Verify open reading frame identification."""
    # ATG + 10 codons + TAA = 36 bp
    orf_seq = "ATG" + "GTC" * 10 + "TAA"
    tokens = scan_tokens(orf_seq)

    types = [t.token_type for t in tokens]
    assert "START_CODON" in types
    assert "STOP_CODON" in types
    assert "OPEN_READING_FRAME" in types


def test_disassemble_sequence_full():
    """Test end-to-end biological sequence disassembly and decompilation."""
    # Full synthetic operon: promoter + RBS + start + codons + stop
    operon_seq = (
        "TTGACACTGATATAATCGT"  # -35 and -10 promoter
        "AGGAGGCC"              # RBS
        "ATG" + "GCC" * 12 + "TAA"  # CDS
    )

    result = disassemble_sequence(operon_seq, routine_name="test_operon")

    assert result.routine_id == "test_operon"
    assert len(result.tokens) > 0
    assert len(result.assembly_listing) > 0
    assert "def test_operon" in result.decompiled_pseudocode
    assert "return halt" in result.decompiled_pseudocode
