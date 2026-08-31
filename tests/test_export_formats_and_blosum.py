"""Test Suite for BLOSUM62 Conservation Relaxation and Vector File Exporters."""

from __future__ import annotations

import pytest

from bio_arch.modules.export_formats import (
    export_to_csv_summary,
    export_to_fasta,
    export_to_genbank,
    export_to_sbol3,
)
from bio_arch.modules.recompiler import (
    blosum62_score,
    recompile_dual_protein_dna,
)


def test_blosum62_scoring_matrix():
    """Verify BLOSUM62 matrix lookups for identical, conservative, and non-conservative pairs."""
    # Identical
    assert blosum62_score("W", "W") == 11
    assert blosum62_score("C", "C") == 9

    # Conservative (favorable positive score)
    assert blosum62_score("I", "V") == 3
    assert blosum62_score("E", "D") == 2
    assert blosum62_score("K", "R") == 2

    # Non-conservative (negative score)
    assert blosum62_score("W", "D") == -4
    assert blosum62_score("P", "W") == -4


def test_recompile_with_blosum62_relaxation():
    """Verify dual-protein recompilation records BLOSUM62 similarity and substituted positions."""
    p0 = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYG"
    p1 = "MASSEDVIKEFMRFKVRMEGSVNGHEFEIEGEGEGRPYEG"

    res = recompile_dual_protein_dna(p0, p1, allow_conservative_mutations=True)

    assert res.blosum62_similarity_pct >= 60.0
    assert isinstance(res.substituted_positions, list)
    assert len(res.substituted_positions) > 0
    # Check substituted position structure
    first_sub = res.substituted_positions[0]
    assert "position" in first_sub
    assert "target_f1" in first_sub
    assert "actual_f1" in first_sub
    assert "blosum62_score" in first_sub


def test_genbank_export_format():
    """Verify that GenBank exporter produces valid .gbk text with CDS annotations."""
    p0 = "MSKGEELFTG"
    p1 = "MASSEDVIKE"
    res = recompile_dual_protein_dna(p0, p1)

    gbk_text = export_to_genbank(res, locus_name="TEST_DUAL_AAV")

    assert "LOCUS       TEST_DUAL_AAV" in gbk_text
    assert "CDS             1..30" in gbk_text
    assert "CDS             2..31" in gbk_text
    assert "/gene=\"Target_Protein_F0\"" in gbk_text
    assert "/gene=\"Target_Protein_F1\"" in gbk_text
    assert "ORIGIN" in gbk_text
    assert "//" in gbk_text


def test_sbol3_xml_export_format():
    """Verify that SBOL3 exporter produces valid XML RDF component structure."""
    p0 = "MSKGEELFTG"
    p1 = "MASSEDVIKE"
    res = recompile_dual_protein_dna(p0, p1)

    sbol_xml = export_to_sbol3(res, display_id="my_aav_construct")

    assert '<?xml version="1.0" encoding="UTF-8"?>' in sbol_xml
    assert "<sbol:Component" in sbol_xml
    assert "<sbol:displayId>my_aav_construct</sbol:displayId>" in sbol_xml
    assert "<sbol:Sequence" in sbol_xml
    assert res.synthesized_dna.upper() in sbol_xml


def test_fasta_and_csv_exports():
    """Verify FASTA and CSV summary reports."""
    p0 = "MSKGEE"
    p1 = "MASSED"
    res = recompile_dual_protein_dna(p0, p1)

    fasta_text = export_to_fasta(res, name="vector_demo")
    assert ">vector_demo_DNA" in fasta_text
    assert ">vector_demo_Protein_Frame_0" in fasta_text
    assert ">vector_demo_Protein_Frame_1" in fasta_text

    csv_text = export_to_csv_summary(res)
    assert "Residue_Index,Frame0_Target_AA" in csv_text
    assert "Frame1_Match_Status" in csv_text
    lines = csv_text.strip().splitlines()
    assert len(lines) == len(p0) + 1
