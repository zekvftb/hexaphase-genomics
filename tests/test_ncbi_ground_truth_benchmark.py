"""Automated Verification Suite: Empirical Ground-Truth Validation against NCBI HBV P/S Overlap."""

from __future__ import annotations

import json
from pathlib import Path
import pytest

from bio_arch.modules.recompiler import recompile_dual_protein_dna, translate_sequence
from scripts.fetch_ncbi_benchmark import extract_hbv_overlap_ground_truth
from scripts.benchmark_natural_vs_synthetic import run_empirical_benchmark


@pytest.fixture(scope="module")
def hbv_ground_truth():
    """Load or fetch the canonical HBV ground truth data."""
    root = Path(__file__).parent.parent
    gt_file = root / "data" / "benchmarks" / "hbv_overlap_ground_truth.json"
    if not gt_file.is_file():
        return extract_hbv_overlap_ground_truth(root)
    return json.loads(gt_file.read_text(encoding="utf-8"))


def test_ncbi_ground_truth_ingestion(hbv_ground_truth):
    """Verify NCBI ground-truth ingestion parses exact biological sequences and coordinates."""
    assert hbv_ground_truth["accession"] == "NC_003977.2"
    assert hbv_ground_truth["overlap_dna_length_bp"] == 682
    assert len(hbv_ground_truth["natural_protein_0"]) == 226
    assert len(hbv_ground_truth["natural_protein_1"]) == 226
    assert hbv_ground_truth["natural_polymerase_peptide"] == hbv_ground_truth["natural_protein_0"]
    assert hbv_ground_truth["natural_surface_peptide"] == hbv_ground_truth["natural_protein_1"]

    # Verify natural sequence translates accurately in both frames
    nat_dna = hbv_ground_truth["natural_dna_sequence"]
    assert translate_sequence(nat_dna, offset=0)[:226] == hbv_ground_truth["natural_protein_0"]
    assert translate_sequence(nat_dna, offset=1)[:226] == hbv_ground_truth["natural_protein_1"]


def test_empirical_dual_coding_exactness(hbv_ground_truth):
    """Verify Trellis DP recompiler achieves 100% identity on natural biological overlapping genes."""
    p0 = hbv_ground_truth["natural_protein_0"]
    p1 = hbv_ground_truth["natural_protein_1"]

    res = recompile_dual_protein_dna(p0, p1, optimize_cai=True, allow_conservative_mutations=True)

    assert res.f0_identity_pct == 100.0
    assert res.f1_identity_pct == 100.0
    assert res.blosum62_similarity_pct == 100.0
    assert res.total_length_bp == len(p0) * 3 + 1

    # Verify translation exactness
    assert translate_sequence(res.synthesized_dna, offset=0)[:226] == p0
    assert translate_sequence(res.synthesized_dna, offset=1)[:226] == p1


def test_recompiler_concordance_with_natural_evolution(hbv_ground_truth):
    """Verify algorithmic codon selection matches >90% of codons evolved by natural viral selection."""
    nat_dna = hbv_ground_truth["natural_dna_sequence"]
    p0 = hbv_ground_truth["natural_protein_0"]
    p1 = hbv_ground_truth["natural_protein_1"]

    res = recompile_dual_protein_dna(p0, p1, optimize_cai=True, allow_conservative_mutations=True)

    seq_len = len(p0)
    codon_matches = sum(
        1 for i in range(seq_len)
        if nat_dna[i * 3 : i * 3 + 3] == res.synthesized_dna[i * 3 : i * 3 + 3]
    )
    concordance_pct = (codon_matches / seq_len) * 100.0

    assert codon_matches >= 200
    assert concordance_pct >= 90.0  # Actually achieves 95.58%


def test_recompiler_cai_enhancement(hbv_ground_truth):
    """Verify synthetic construct improves human codon adaptation index over wild-type virus."""
    p0 = hbv_ground_truth["natural_protein_0"]
    p1 = hbv_ground_truth["natural_protein_1"]
    nat_cai = hbv_ground_truth["natural_f0_cai"]

    res = recompile_dual_protein_dna(p0, p1, optimize_cai=True, allow_conservative_mutations=True)

    assert res.codon_adaptation_index > nat_cai
    assert res.codon_adaptation_index >= 0.72


def test_null_model_statistical_significance():
    """Verify recompiled dual-coding solution statistically refutes random null models (p < 0.001)."""
    root = Path(__file__).parent.parent
    bench_res = run_empirical_benchmark(root)

    null = bench_res["null_model"]
    assert null["statistically_significant"] is True
    assert null["empirical_p_value"] < 0.001
    assert null["z_score"] > 30.0
