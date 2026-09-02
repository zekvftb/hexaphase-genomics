"""Automated Verification Suite for Pre-Registered Empirical Hypothesis Testing."""

from __future__ import annotations

import json
from pathlib import Path
import pytest

from scripts.empirical.test_kolmogorov_compression import (
    lz76_complexity,
    run_compression_experiment,
    sliding_window_shannon_entropy,
    zlib_compression_ratio,
)
from scripts.empirical.test_genomic_grammar_zipf import (
    fit_zipf_power_law,
    run_zipf_experiment,
)
from scripts.empirical.test_automaton_stability import (
    execute_pda,
    run_automaton_experiment,
)
from scripts.empirical.run_all_hypothesis_tests import run_full_empirical_evaluation


def test_lz76_and_compression_metrics():
    """Verify Lempel-Ziv complexity and Shannon entropy calculations."""
    seq_rep = "ATGC" * 100  # 400 bp highly repetitive
    seq_rand = ("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC" * 6)[:400]

    lz_rep = lz76_complexity(seq_rep)
    lz_rand = lz76_complexity(seq_rand)
    comp_rep = zlib_compression_ratio(seq_rep)
    comp_rand = zlib_compression_ratio(seq_rand)

    # For equal length (400 bp), repetitive sequence has lower compression ratio
    assert comp_rep < comp_rand
    assert lz_rep < lz_rand

    ent_64 = sliding_window_shannon_entropy(seq_rand, window=32)
    assert 1.0 < ent_64 <= 2.0


def test_zipf_power_law_fitting():
    """Verify linear regression on log-log k-mer rank-frequency distributions."""
    seq = "ATGCGATCGATCGATCGATCGATCGATCGATC" * 10
    alpha, r2, unique_kmers = fit_zipf_power_law(seq, k=4)

    assert unique_kmers > 0
    assert 0.0 <= r2 <= 1.0
    assert alpha >= 0.0


def test_push_down_automaton_execution():
    """Verify deterministic PDA execution, stack bounds, and halt detection."""
    # Sequence with clean stop codon at end
    dna = "ATGCTGGTTAAATAG"  # M, L, V, K, *
    res = execute_pda(dna, frame_offset=0, max_steps=50, max_stack=10)

    assert res["steps_executed"] > 0
    assert res["terminated_cleanly"] == 1.0
    assert res["underflow_crashes"] == 0
    assert res["max_stack_depth"] >= 0


def test_raw_hypothesis_ledger_and_report_artifacts():
    """Verify raw json ledger and markdown report artifacts exist and contain valid results."""
    root = Path(__file__).parent.parent
    raw_json = root / "outputs" / "empirical_hypothesis_testing_raw.json"
    report_file = root / "EMPIRICAL_COMPUTATIONAL_HYPOTHESIS_REPORT.md"

    assert raw_json.is_file(), "Raw hypothesis ledger missing in outputs/"
    assert report_file.is_file(), "Master empirical hypothesis report missing"

    data = json.loads(raw_json.read_text(encoding="utf-8"))
    assert "experiment_a_compression" in data
    assert "experiment_b_zipf_grammar" in data
    assert "experiment_c_automaton_stability" in data
    assert data["significance_threshold_alpha"] == 0.01


def test_quick_empirical_runner_integration():
    """Verify end-to-end hypothesis evaluation executes cleanly on small shuffle counts."""
    root = Path(__file__).parent.parent
    ledger = run_full_empirical_evaluation(root, n_shuffles=10)

    assert "experiment_a_compression" in ledger
    assert "NC_003977.2" in ledger["experiment_a_compression"]
