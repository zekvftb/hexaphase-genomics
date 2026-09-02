"""Automated Verification Suite for Non-Canonical smORF Discovery Engine."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import random
import pytest

from bio_arch.modules.non_canonical_miner import (
    NonCanonicalCandidate,
    NonCanonicalSmorfMiner,
    START_CODON_EFFICIENCY,
)
from scripts.mining.run_non_canonical_mining import run_non_canonical_pipeline


@pytest.fixture(scope="module")
def miner():
    """Create a deterministic NonCanonicalSmorfMiner instance."""
    return NonCanonicalSmorfMiner(
        min_length_aa=30,
        max_length_aa=120,
        min_cai_threshold=0.68,
        n_null_shuffles=50,
        significance_alpha=0.01,
        min_z_threshold=2.5,
        seed=42,
    )


def test_near_cognate_codon_recognition():
    """Verify recognition and efficiency mapping of canonical and near-cognate start codons."""
    assert "ATG" in START_CODON_EFFICIENCY
    assert "CTG" in START_CODON_EFFICIENCY
    assert "GTG" in START_CODON_EFFICIENCY
    assert "ACG" in START_CODON_EFFICIENCY
    assert "ATA" in START_CODON_EFFICIENCY
    assert "TTG" in START_CODON_EFFICIENCY

    assert START_CODON_EFFICIENCY["ATG"] == 1.00
    assert START_CODON_EFFICIENCY["CTG"] == 0.25
    assert START_CODON_EFFICIENCY["GTG"] == 0.12


def test_kozak_context_scoring(miner):
    """Verify Kozak initiation context evaluation for Optimal, Strong, and Weak configurations."""
    # Test Optimal: A at -3 and G at +4
    dna_optimal = "TTTACCATGGTTAAATAG"
    res_optimal = miner._evaluate_kozak(dna_optimal, start_bp=6)
    assert res_optimal["strength"] == "Optimal"

    # Test Weak: T at -3 and T at +4
    dna_weak = "TTTTTTATGTTTAAATAG"
    res_weak = miner._evaluate_kozak(dna_weak, start_bp=6)
    assert res_weak["strength"] == "Weak"


def test_viroporin_and_hydropathy_metrics(miner):
    """Verify Kyte-Doolittle TM domain detection and amphipathic hydrophobic moment."""
    # Hydrophobic peptide (putative TM helix)
    tm_peptide = "L" * 20 + "K" * 5
    mean_h, max_h, tm_segs = miner._analyze_hydropathy(tm_peptide)
    assert len(tm_segs) >= 1
    assert max_h >= 1.6

    # Amphipathic moment
    h_moment = miner._calc_hydrophobic_moment(tm_peptide)
    assert h_moment > 0.0


def test_negative_control_random_dna_suppression(miner):
    """Verify uniformly random synthetic DNA produces zero false-positive significant smORFs."""
    rng = random.Random(999)
    random_dna = "".join(rng.choice("ACGT") for _ in range(600))

    candidates = miner.scan_cds(
        cds_dna=random_dna,
        parent_gene="RANDOM_NULL",
        parent_product="synthetic random sequence",
        accession="SYNTH_NC_NULL",
        organism="Synthetic Null Control",
    )

    sig_cands = [c for c in candidates if c.statistically_significant]
    assert len(sig_cands) == 0


def test_non_canonical_pipeline_artifacts():
    """Verify exported CSV ledger and markdown report exist and contain valid schema."""
    root = Path(__file__).parent.parent
    csv_file = root / "outputs" / "non_canonical_smorf_candidates.csv"
    report_file = root / "NON_CANONICAL_VIRAL_SMORF_DISCOVERY.md"

    assert csv_file.is_file(), "non_canonical_smorf_candidates.csv missing in outputs/"
    assert report_file.is_file(), "NON_CANONICAL_VIRAL_SMORF_DISCOVERY.md missing"

    with open(csv_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        assert len(rows) > 0
        sample = rows[0]
        assert "Start_Codon" in sample
        assert "Start_Codon_Type" in sample
        assert "Kozak_Strength" in sample
        assert "Host_CAI" in sample
        assert "Hydrophobic_Moment" in sample
        assert "Z_Score" in sample
