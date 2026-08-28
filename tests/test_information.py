"""Unit tests for Module 1: DNA/RNA Information Architecture."""

from collections import Counter
import math
from pathlib import Path
import random
import pytest

from bio_arch.contracts import EvidenceClass
from bio_arch.modules.information import (
    analyze_sequence,
    autocorrelation,
    benjamini_hochberg_correction,
    compare_against_null,
    compression_ratio,
    compute_composition,
    conditional_entropy,
    detect_tandem_repeats,
    run_module_1,
    shannon_entropy,
    shuffle_dinucleotide,
    shuffle_mononucleotide,
)


def test_compute_composition():
    """Verify composition, GC content, and GC/AT skew formulas."""
    # Balanced ACGT
    comp = compute_composition("ACGT")
    assert comp.length == 4
    assert comp.gc_content == 0.5
    assert comp.gc_skew == 0.0
    assert comp.at_skew == 0.0

    # High GC with G bias: GGGCCC -> 3 G, 3 C
    comp_gc = compute_composition("GGGCCC")
    assert comp_gc.gc_content == 1.0
    assert comp_gc.gc_skew == 0.0

    # G-skew: 3 G, 1 C -> (3 - 1) / (3 + 1) = 0.5
    comp_skew = compute_composition("GGGC")
    assert comp_skew.gc_skew == 0.5

    # Empty string
    empty = compute_composition("")
    assert empty.length == 0
    assert empty.gc_content == 0.0


def test_shannon_entropy():
    """Verify Shannon entropy computation in bits."""
    # Homopolymer has 0 entropy
    assert shannon_entropy("AAAAAAA", k=1) == 0.0

    # Perfectly uniform 4-base sequence: max entropy is log2(4) = 2.0 bits
    uniform = "ACGT" * 100
    h1 = shannon_entropy(uniform, k=1)
    assert math.isclose(h1, 2.0, abs_tol=1e-3)

    # Bounded k-mer check (k > 6 must raise ValueError)
    with pytest.raises(ValueError, match="k-mer size must be between 1 and 6"):
        shannon_entropy("ACGT", k=7)


def test_conditional_entropy():
    """Verify conditional entropy H(X_t | X_{t-1})."""
    # Deterministic sequence: knowing previous base completely predicts next
    repeating = "AC" * 100
    # In repeating "AC", after 'A' always 'C', after 'C' always 'A' -> H(X_t|X_{t-1}) is near 0
    cond = conditional_entropy(repeating, k=1)
    assert cond < 0.1


def test_compression_ratio_complexity():
    """Verify compression ratio reflects descriptive sequence complexity."""
    repetitive = "A" * 1000
    # Highly repetitive sequence compresses down significantly
    rep_ratio = compression_ratio(repetitive)

    rng = random.Random(42)
    random_seq = "".join(rng.choices("ACGT", k=1000))
    rand_ratio = compression_ratio(random_seq)

    assert rep_ratio < rand_ratio


def test_tandem_repeats_detection():
    """Verify exact tandem repeat identification."""
    seq = "ATGC" + "CAG" * 5 + "TTTT"
    repeats = detect_tandem_repeats(seq, min_unit=3, max_unit=3, min_repeats=3)
    assert len(repeats) >= 1
    assert repeats[0]["unit"] == "CAG"
    assert repeats[0]["count"] == 5


def test_autocorrelation():
    """Verify lag autocorrelation calculation."""
    # Period-3 repeat: ATGATGATGATG
    seq = "ATG" * 20
    corr = autocorrelation(seq, max_lag=6)
    # Lag 3 should have 100% match
    assert corr[3] == 1.0
    # Lag 1 should have 0% match (A!=T, T!=G, G!=A)
    assert corr[1] == 0.0


def test_shuffle_mononucleotide():
    """Verify single-nucleotide frequencies are exactly preserved."""
    seq = "AAACCCGGGTTT"
    original_counts = Counter(seq)

    rng = random.Random(42)
    shuffled = shuffle_mononucleotide(seq, rng)

    assert len(shuffled) == len(seq)
    assert Counter(shuffled) == original_counts


def test_shuffle_dinucleotide():
    """Verify dinucleotide transition frequencies are exactly preserved."""
    seq = "ATGCAATTTGCTAGCTAGCATGCA"
    rng = random.Random(42)

    def get_dinuc_counts(s: str) -> Counter:
        return Counter(s[i : i + 2] for i in range(len(s) - 1))

    original_dinucs = get_dinuc_counts(seq)
    shuffled = shuffle_dinucleotide(seq, rng)

    assert len(shuffled) == len(seq)
    assert Counter(shuffled) == Counter(seq)
    assert get_dinuc_counts(shuffled) == original_dinucs


def test_benjamini_hochberg_correction():
    """Verify FDR p-value correction."""
    p_vals = [0.01, 0.04, 0.03, 0.20]
    adj = benjamini_hochberg_correction(p_vals)

    assert len(adj) == len(p_vals)
    # Adjusted p-values must be >= raw p-values
    for raw, corrected in zip(p_vals, adj):
        assert corrected >= raw
        assert 0.0 <= corrected <= 1.0


def test_run_module_1(tmp_path: Path):
    """Integration test for Module 1 execution and contract output."""
    sample_fasta = tmp_path / "test.fasta"
    sample_fasta.write_text(">test_gene\nATGCATGCATGCATGCATGCATGCATGCATGC\n")

    result = run_module_1(
        input_file=sample_fasta,
        output_dir=tmp_path,
        num_shuffles=20,
        seed=123,
    )

    assert result.run_metadata.status == "success"
    assert len(result.findings) >= 3
    assert len(result.interpretations) >= 1
    assert result.interpretations[0].classification == EvidenceClass.INTERPRETATION

    # Verify output JSON artifact exists
    out_file = tmp_path / "test_information_metrics.json"
    assert out_file.is_file()
