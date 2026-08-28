"""Unit tests for bio_arch.modules.biological_circuits (Module 8)."""

import pytest
from bio_arch.contracts import (
    BiologicalCircuitReport,
    CpgMemoryIsland,
    CrisprArray,
    RiboswitchAdc,
)
from bio_arch.modules.biological_circuits import (
    scan_all_biological_circuits,
    scan_cpg_islands,
    scan_crispr_arrays,
    scan_riboswitch_adcs,
)


def test_scan_cpg_islands_positive():
    """Verify detection of a bona fide CpG island memory register."""
    # Construct synthetic CpG island with high GC (60%) and high CpG frequency
    cpg_window = ("CGCG" * 25) + ("ATAT" * 25)  # 200 bp
    flank = "A" * 100
    full_seq = flank + cpg_window + flank

    islands = scan_cpg_islands(full_seq, window_size=200, min_gc=0.45, min_obs_exp=0.55)
    assert len(islands) >= 1
    isl = islands[0]
    assert isl.length_bp >= 200
    assert isl.cpg_obs_exp_ratio >= 0.55
    assert isl.gc_content_pct >= 45.0


def test_scan_crispr_arrays_synthetic():
    """Verify identification and extraction of direct repeats and viral spacers."""
    repeat = "GTTTTCCCCGCACCACGCGGGGA"  # 23-24 bp palindromic repeat
    spacer1 = "ATCGATCGATCGATCGATCGATCGATCG"  # 28 bp spacer
    spacer2 = "TGCATGCATGCATGCATGCATGCATGCA"  # 28 bp spacer
    spacer3 = "GCTAGCTAGCTAGCTAGCTAGCTAGCTA"  # 28 bp spacer

    crispr_locus = repeat + spacer1 + repeat + spacer2 + repeat + spacer3 + repeat
    arrays = scan_crispr_arrays(crispr_locus, min_repeat_len=len(repeat), min_repeats=3)

    assert len(arrays) >= 1
    arr = arrays[0]
    assert arr.repeats_count >= 3
    assert arr.spacers_count >= 2
    assert spacer1 in arr.spacers or spacer2 in arr.spacers


def test_scan_riboswitch_adcs_positive():
    """Verify detection of a metabolite-sensing RNA Riboswitch ADC with hairpin."""
    # Purine aptamer motif + 6bp GC stem-loop terminator
    aptamer = "ACTCC" + "AAAAAAAAAAAA" + "GGAGT"  # 22 nt
    terminator = "GCGCGC" + "AAAA" + "GCGCGC"  # Hairpin dG < -10
    seq = "NNNN" + aptamer + terminator + "NNNN"

    switches = scan_riboswitch_adcs(seq, min_switching_mfe=-7.0)
    assert len(switches) >= 1
    sw = switches[0]
    assert sw.ligand_class == "PURINE_GUANINE"
    assert sw.terminator_mfe_dG <= -7.0
    assert sw.predicted_state == "OFF_IN_HIGH_LIGAND"


def test_scan_all_biological_circuits_e2e():
    """End-to-end audit report generation."""
    seq = "ATGCGTACGT" * 50
    report = scan_all_biological_circuits(seq, genome_id="TEST_CIRCUITS")
    assert isinstance(report, BiologicalCircuitReport)
    assert report.genome_id == "TEST_CIRCUITS"
    assert "crispr_security_arrays_found" in report.summary
