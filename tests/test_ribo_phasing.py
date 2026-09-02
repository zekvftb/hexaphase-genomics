"""Automated Verification Suite for Ribo-seq Phasing & Translation Evidence Engine."""

from __future__ import annotations

from pathlib import Path
import pytest

from bio_arch.modules.ribo_phasing import (
    analyze_ribo_phasing,
    calculate_chi_square_uniformity,
    generate_synthetic_ribo_reads,
    get_psite_offset,
)
from bio_arch.modules.ribo_registry import PUBLIC_RIBO_DATASETS, get_dataset_by_target
from scripts.validation.verify_candidates_riboseq import run_riboseq_verification_pipeline


def test_psite_offset_calibration():
    """Verify length-dependent P-site offsets (28nt -> 12, 30nt -> 13)."""
    assert get_psite_offset(28) == 12
    assert get_psite_offset(29) == 12
    assert get_psite_offset(30) == 13
    assert get_psite_offset(31) == 13
    assert get_psite_offset(35) == 12  # fallback default


def test_chi_square_uniformity_calculation():
    """Verify Chi-square statistic and analytic p-value calculation (df=2)."""
    # Uniform counts -> Chi2 = 0, p = 1.0
    chi2_uni, p_uni = calculate_chi_square_uniformity([100, 100, 100])
    assert chi2_uni == 0.0
    assert p_uni == 1.0

    # Highly biased counts -> high Chi2, p < 0.001
    chi2_bias, p_bias = calculate_chi_square_uniformity([150, 10, 10])
    assert chi2_bias > 50.0
    assert p_bias < 0.001


def test_reading_frame_phasing_analysis():
    """Verify phasing analysis identifies dominant target frame with high TPI."""
    reads = generate_synthetic_ribo_reads(
        start_nt=1000,
        end_nt=1300,
        target_frame=1,
        signal_to_noise=0.85,
        total_reads=100,
        seed=42,
    )

    res = analyze_ribo_phasing(
        read_alignments=reads,
        candidate_id="TEST_CAND_01",
        accession="NC_000001.1",
        parent_gene="TEST_GENE",
        start_nt=1000,
        end_nt=1300,
        target_frame=1,
        min_coverage_threshold=15,
    )

    assert res.total_psites == 100
    assert res.triplet_periodicity_index >= 0.70
    assert res.frame_counts[1] > res.frame_counts[0]
    assert res.frame_counts[1] > res.frame_counts[2]
    assert res.chi_square_p_value < 0.01
    assert res.classification == "Confirmed Active Translation"


def test_public_ribo_registry():
    """Verify curated dataset registry contains required viral host models."""
    assert len(PUBLIC_RIBO_DATASETS) >= 3
    sv40_ds = get_dataset_by_target("SV40")
    assert sv40_ds is not None
    assert "PRJNA639148" in sv40_ds.bioproject_id

    hbv_ds = get_dataset_by_target("HBV")
    assert hbv_ds is not None
    assert "PRJNA516397" in hbv_ds.bioproject_id


def test_riboseq_verification_report_artifact():
    """Verify end-to-end report generation and markdown contents."""
    root = Path(__file__).parent.parent
    report_file = root / "outputs" / "RIBOSEQ_CANDIDATE_TRANSLATION_REPORT.md"

    assert report_file.is_file(), "RIBOSEQ_CANDIDATE_TRANSLATION_REPORT.md missing"
    content = report_file.read_text(encoding="utf-8")
    assert "Public Ribo-seq Mining & Triplet Phasing Verification Report" in content
    assert "Candidate_1_SV40" in content
    assert "Candidate_2_PhiX174" in content
    assert "Candidate_3_HBV" in content
    assert "Confirmed Active Translation" in content
