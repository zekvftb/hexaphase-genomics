"""Test Suite for Synthetic Dual-Coding, AAV Packaging Optimizer, and DNA Storage Parity."""

from __future__ import annotations

import pytest

from bio_arch.modules.dna_storage import (
    decode_storage_payload,
    encode_storage_payload,
)
from bio_arch.modules.recompiler import (
    calculate_aav_packaging_savings,
    compute_cai,
    count_homopolymer_runs,
    recompile_dual_protein_dna,
    scan_restriction_sites,
    translate_sequence,
)


def test_recompile_dual_protein_with_cai_and_quality():
    """Verify dual-protein recompilation generates high-identity overlapping DNA with valid CAI."""
    # Target proteins: Green Fluorescent Protein short epitope vs RFP short epitope
    p0 = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYG"
    p1 = "MASSEDVIKEFMRFKVRMEGSVNGHEFEIEGEGEGRPYEG"

    res = recompile_dual_protein_dna(p0, p1, optimize_cai=True, filter_restriction=True)

    assert res.total_length_bp == len(p0) * 3 + 1
    assert res.f0_identity_pct == 100.0
    assert res.f1_identity_pct >= 25.0
    assert res.compression_ratio >= 1.85
    assert 0.0 < res.codon_adaptation_index <= 1.0

    # Ensure translation matches synthesized DNA
    trans_0 = translate_sequence(res.synthesized_dna, offset=0)
    trans_1 = translate_sequence(res.synthesized_dna, offset=1)
    assert trans_0[:len(p0)] == res.translated_f0
    assert trans_1[:len(p1)] == res.translated_f1


def test_recompile_dual_protein_designed_pair():
    """Verify dual-protein recompilation on a co-designed dual-coding peptide pair."""
    # Co-designed dual-coding pair (where Frame +1 amino acids are reachable via Frame 0 codons)
    # E.g. "AAAA" in F0 can yield "LLLL" or "QQQQ" in F1
    # GCT (A), GCA (A), GCC (A), GCG (A) -> Suffixes CT, CA, CC, CG + next G -> CTG (L), CAG (Q), CCG (P), CGG (R)
    p0 = "AAAA" * 5
    p1 = "LLLL" * 5

    res = recompile_dual_protein_dna(p0, p1, optimize_cai=True)
    assert res.f0_identity_pct == 100.0
    assert res.f1_identity_pct == 100.0


def test_aav_vector_packaging_savings():
    """Verify AAV vector capacity calculations and footprint reduction."""
    # Two standard 400 aa therapeutic proteins
    len_f0 = 400
    len_f1 = 400
    synthesized_bp = (400 * 3) + 1  # 1201 bp

    savings = calculate_aav_packaging_savings(
        protein_len_f0_aa=len_f0,
        protein_len_f1_aa=len_f1,
        synthesized_dna_bp=synthesized_bp,
        aav_limit_bp=4700,
    )

    # Separate genes: 2*(600+250) + 2400 + 290 = 4390 bp
    # Compact construct: (600+250) + 1201 + 290 = 2341 bp
    assert savings["separate_dual_cassette_total_bp"] == 4390
    assert savings["compact_overlapping_cassette_total_bp"] == 2341
    assert savings["bp_saved"] == 2049
    assert savings["percent_footprint_reduction"] > 45.0
    assert savings["fits_standard_aav"] is True
    assert savings["headroom_remaining_bp"] > 2000


def test_restriction_site_and_homopolymer_detection():
    """Verify detection of restriction enzyme recognition sequences and homopolymers."""
    seq_with_sites = "ATGCGAATTCGGATCCTAAGCTTAAAAAAAAGGGGGG"
    
    sites = scan_restriction_sites(seq_with_sites)
    assert any("EcoRI" in s for s in sites)
    assert any("BamHI" in s for s in sites)
    assert any("HindIII" in s for s in sites)

    homopolymer_count = count_homopolymer_runs(seq_with_sites, max_run=5)
    assert homopolymer_count >= 2


def test_dna_storage_encoding_and_parity_integrity():
    """Verify encoding digital binary data into DNA with frame-offset parity checks."""
    test_payload = "HexaPhase DNA Storage v2.0 - Parity Check Verified!"
    packet = encode_storage_payload(test_payload)

    assert packet.original_bytes_count == len(test_payload.encode("utf-8"))
    assert packet.total_length_bp > 0
    assert packet.information_density_bits_per_base > 1.0

    # Decode and verify exact roundtrip
    decoded = decode_storage_payload(packet.interleaved_dna)
    assert decoded["parity_valid"] is True
    assert decoded["recovered_bytes"].decode("utf-8") == test_payload
    assert decoded["sha256_checksum"] == packet.sha256_checksum


def test_dna_storage_error_detection_on_corruption():
    """Verify parity check fails when DNA sequence is corrupted."""
    test_payload = "Sensitive Genomic Cryptographic Key"
    packet = encode_storage_payload(test_payload)

    # Inject substitution corruption in primary DNA (after 16-base header)
    corrupted_dna = packet.interleaved_dna[:20] + ("T" if packet.interleaved_dna[20] != "T" else "A") + packet.interleaved_dna[21:]
    decoded_corrupted = decode_storage_payload(corrupted_dna)

    # Parity check catches the alteration
    assert decoded_corrupted["parity_valid"] is False
