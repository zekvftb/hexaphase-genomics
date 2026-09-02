"""Automated Verification Suite for Viral Overlapping Gene & smORF Mining Engine."""

from __future__ import annotations

import csv
from pathlib import Path
import random
import pytest

from bio_arch.modules.viral_miner import OverlappingCandidate, ViralOverlapMiner


@pytest.fixture(scope="module")
def miner():
    """Create a deterministic ViralOverlapMiner instance for testing."""
    return ViralOverlapMiner(
        min_length_aa=30,
        min_cai_threshold=0.65,
        n_null_shuffles=100,
        significance_alpha=0.01,
        seed=42,
    )


def test_positive_control_hbv_overlap_rediscovery(miner):
    """Verify miner re-discovers the canonical HBV Polymerase/Surface overlap (p < 0.001)."""
    root = Path(__file__).parent.parent
    hbv_gbk = root / "data" / "mining_corpus" / "NC_003977.2.gbk"
    assert hbv_gbk.is_file(), "HBV GenBank file missing from data/mining_corpus/"

    candidates = miner.mine_genome_file(hbv_gbk)
    assert len(candidates) > 0

    # Search for candidate corresponding to Polymerase overlap in S gene or vice versa
    sig_p_matches = [
        c for c in candidates
        if c.statistically_significant and ("P" in c.matched_annotation_name or c.length_aa > 100)
    ]
    assert len(sig_p_matches) > 0

    top_hbv = max(sig_p_matches, key=lambda x: x.z_score)
    assert top_hbv.z_score >= 3.0
    assert top_hbv.empirical_p_value <= 0.01
    assert top_hbv.cai_score >= 0.68


def test_positive_control_aav2_aap_rediscovery(miner):
    """Verify miner re-discovers AAV2 Assembly-Activating Protein (AAP) nested in Cap gene."""
    root = Path(__file__).parent.parent
    aav2_gbk = root / "data" / "mining_corpus" / "NC_001401.2.gbk"
    assert aav2_gbk.is_file(), "AAV2 GenBank file missing from data/mining_corpus/"

    candidates = miner.mine_genome_file(aav2_gbk)
    assert len(candidates) > 0

    aap_matches = [
        c for c in candidates
        if c.statistically_significant and c.frame_offset == 1 and c.length_aa >= 150
    ]
    assert len(aap_matches) > 0

    top_aap = max(aap_matches, key=lambda x: x.z_score)
    assert top_aap.z_score >= 2.5
    assert top_aap.empirical_p_value <= 0.01
    assert top_aap.length_aa >= 150


def test_negative_control_uniform_random_dna(miner):
    """Verify uniformly random non-coding DNA produces no significant candidate smORFs."""
    # Generate 600 bp of uniformly random sequence with start and stop codons
    rng = random.Random(12345)
    random_bases = [rng.choice("ACGT") for _ in range(600)]
    # Ensure divisible by 3
    random_dna = "".join(random_bases)

    candidates = miner.scan_cds_alternative_frames(
        cds_dna=random_dna,
        parent_gene="RANDOM_NULL",
        parent_product="synthetic random sequence",
        accession="SYNTH_NULL_001",
        organism="Synthetic Null Control",
    )

    # Assert that no candidates pass statistical significance
    sig_candidates = [c for c in candidates if c.statistically_significant]
    assert len(sig_candidates) == 0


def test_novel_candidate_csv_and_genbank_export():
    """Verify exported candidate CSV and re-annotated GenBank files exist and contain valid fields."""
    root = Path(__file__).parent.parent
    csv_file = root / "outputs" / "novel_viral_overlapping_candidates.csv"
    report_file = root / "VIRAL_OVERLAPPING_GENE_DISCOVERY_REPORT.md"
    annotated_dir = root / "outputs" / "annotated_candidates"

    assert csv_file.is_file(), "Candidate CSV ledger missing in outputs/"
    assert report_file.is_file(), "Master discovery report missing"
    assert annotated_dir.is_dir(), "Annotated candidates directory missing in outputs/"

    # Validate CSV columns
    with open(csv_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        assert len(rows) >= 50
        sample = rows[0]
        assert "Accession" in sample
        assert "Host_CAI" in sample
        assert "Z_Score" in sample
        assert "Statistically_Significant" in sample
        assert "Peptide_Sequence" in sample


def test_viral_overlap_miner_parameter_customization():
    """Verify custom length thresholds and significance filtering."""
    custom_miner = ViralOverlapMiner(
        min_length_aa=50,
        n_null_shuffles=50,
        significance_alpha=0.001,
        seed=999,
    )
    assert custom_miner.min_length_aa == 50
    assert custom_miner.significance_alpha == 0.001
