"""Automated Verification Suite for ESMFold In-Silico Structural Prediction Client."""

from __future__ import annotations

from pathlib import Path
import pytest

from bio_arch.modules.esmfold_client import ESMFoldClient
from scripts.structural.evaluate_smorf_folding import calc_pearson_correlation, run_structural_evaluation_pipeline


@pytest.fixture
def client(tmp_path):
    """Create a temporary ESMFold client instance."""
    return ESMFoldClient(cache_dir=tmp_path / "structures")


def test_heuristic_pdb_generation_and_plddt_parsing(client):
    """Verify that heuristic fallback PDB contains valid ATOM records and parses pLDDT cleanly."""
    peptide = "MMKTARKMLIKMKMVGRRTWKTQGMKQALIHSPKAHFRPLSPHSLFMIIISHTTFVEVLLALKNLPHLPLNLKH"
    pdb_text = client.generate_heuristic_pdb(peptide)

    assert "HEADER" in pdb_text
    assert "ATOM" in pdb_text
    assert "TER" in pdb_text
    assert "END" in pdb_text

    plddt_scores = client.parse_plddt_from_pdb(pdb_text)
    assert len(plddt_scores) == len(peptide)
    assert all(0.0 <= s <= 100.0 for s in plddt_scores)


def test_structural_confidence_evaluation(client):
    """Verify global and TM-core structural confidence calculations."""
    peptide = "MSLIKMMLVMVSVAAISKTFGLLRFLLRLSFLAK"
    res = client.evaluate_structural_confidence(
        peptide=peptide,
        tm_start_aa=8,
        tm_end_aa=26,
        use_cache=False,
    )

    assert res["length_aa"] == len(peptide)
    assert res["global_mean_plddt"] > 50.0
    assert res["tm_core_mean_plddt"] is not None
    assert res["tm_core_mean_plddt"] > 60.0
    assert res["confidence_tier"] in ("High Confidence", "Moderate Confidence")
    assert len(res["helical_blocks"]) >= 1


def test_esmfold_caching_mechanism(client):
    """Verify that repeated sequence queries load from the local cache."""
    peptide = "MAGRSWVLPTY"
    pdb_1 = client.fold_sequence(peptide, use_cache=True)
    cache_files = list(client.cache_dir.glob("*.pdb"))
    assert len(cache_files) == 1

    # Second call should load identical text from cache
    pdb_2 = client.fold_sequence(peptide, use_cache=True)
    assert pdb_1 == pdb_2


def test_pearson_correlation_calculation():
    """Verify Pearson correlation calculation."""
    x = [1.0, 2.0, 3.0, 4.0, 5.0]
    y = [2.0, 4.0, 6.0, 8.0, 10.0]
    corr = calc_pearson_correlation(x, y)
    assert corr == 1.0

    y_inv = [5.0, 4.0, 3.0, 2.0, 1.0]
    corr_inv = calc_pearson_correlation(x, y_inv)
    assert corr_inv == -1.0


def test_structural_report_and_pdb_artifacts():
    """Verify generated PDB files and structural report artifacts."""
    root = Path(__file__).parent.parent
    structures_dir = root / "outputs" / "structures"
    report_file = root / "outputs" / "SMORF_STRUCTURAL_FOLDING_REPORT.md"

    assert structures_dir.is_dir()
    assert report_file.is_file()

    pdb_files = list(structures_dir.glob("*.pdb"))
    assert len(pdb_files) >= 3

    report_content = report_file.read_text(encoding="utf-8")
    assert "ESMFold In-Silico Structural Prediction" in report_content
    assert "Candidate_1_SV40" in report_content
    assert "Candidate_2_PhiX174" in report_content
    assert "Candidate_3_HBV" in report_content
