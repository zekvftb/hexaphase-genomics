"""Industry-Standard Vector Exporters for Synthetic Biology Constructs.

Exports dual-coding recompiled genetic constructs to:
- Annotated GenBank format (.gbk) with CDS features for Frame 0 and Frame +1
- Synthetic Biology Open Language (SBOL v3) XML/RDF format
- Standard multi-FASTA format
- CSV summary alignment report
"""

from __future__ import annotations

import csv
import datetime
import io
from typing import Any

from bio_arch.modules.recompiler import RecompilationResult


def export_to_genbank(
    res: RecompilationResult,
    locus_name: str = "SYNTH_DUAL_VECTOR",
    definition: str = "Synthetic dual-coding compact expression vector",
) -> str:
    """Generate an annotated GenBank (.gbk) file representation."""
    dna_seq = res.synthesized_dna.lower()
    seq_len = len(dna_seq)
    today = datetime.datetime.now().strftime("%d-%b-%Y").upper()

    f0_end = len(res.target_protein_f0) * 3
    f1_end = len(res.target_protein_f1) * 3 + 1

    lines = [
        f"LOCUS       {locus_name:<16} {seq_len:>6} bp    DNA     linear   SYN {today}",
        f"DEFINITION  {definition}.",
        f"ACCESSION   {locus_name}",
        "VERSION     1.0",
        "KEYWORDS    synthetic biology; dual-coding; HexaPhase; AAV optimization.",
        "SOURCE      synthetic construct",
        "  ORGANISM  synthetic construct",
        "            other sequences; artificial sequences.",
        "FEATURES             Location/Qualifiers",
        f"     source          1..{seq_len}",
        '                     /organism="synthetic construct"',
        '                     /mol_type="other DNA"',
        f"     CDS             1..{f0_end}",
        '                     /gene="Target_Protein_F0"',
        '                     /codon_start=1',
        f'                     /translation="{res.translated_f0}"',
        f'                     /note="Frame 0 primary CDS (Identity: {res.f0_identity_pct}%)"',
        f"     CDS             2..{f1_end}",
        '                     /gene="Target_Protein_F1"',
        '                     /codon_start=1',
        f'                     /translation="{res.translated_f1}"',
        f'                     /note="Frame +1 overlapping CDS (Identity: {res.f1_identity_pct}%, BLOSUM62: {res.blosum62_similarity_pct}%)"',
        '     misc_feature    1..' + str(seq_len),
        f'                     /note="Human CAI: {res.codon_adaptation_index}; Compression: {res.compression_ratio}x"',
    ]

    # Add restriction site notes if any
    for rs in res.restriction_sites_detected:
        lines.append(f'     misc_feature    1..{seq_len}')
        lines.append(f'                     /note="Restriction Site Detected: {rs}"')

    lines.append("ORIGIN")
    # Format sequence in standard 60-bp numbered lines
    for i in range(0, seq_len, 60):
        chunk = dna_seq[i : i + 60]
        # Split into 10-bp blocks
        blocks = " ".join(chunk[j : j + 10] for j in range(0, len(chunk), 10))
        lines.append(f"{i + 1:>9} {blocks}")
    lines.append("//")
    lines.append("")

    return "\n".join(lines)


def export_to_sbol3(
    res: RecompilationResult,
    display_id: str = "synth_dual_vector",
    namespace: str = "https://hexaphase.org/constructs/",
) -> str:
    """Generate Synthetic Biology Open Language (SBOL v3) XML/RDF specification."""
    dna_seq = res.synthesized_dna.upper()
    seq_len = len(dna_seq)
    uri = f"{namespace}{display_id}"

    xml_lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"',
        '         xmlns:sbol="http://sbols.org/v3#"',
        '         xmlns:dcterms="http://purl.org/dc/terms/"',
        '         xmlns:so="https://identifiers.org/SO:">',
        f'  <sbol:Component rdf:about="{uri}">',
        f'    <sbol:displayId>{display_id}</sbol:displayId>',
        f'    <sbol:name>{display_id}</sbol:name>',
        f'    <dcterms:description>Dual-coding synthetic construct (Compression: {res.compression_ratio}x, CAI: {res.codon_adaptation_index})</dcterms:description>',
        '    <sbol:type rdf:resource="https://identifiers.org/SBO:0000251"/>',
        f'    <sbol:hasSequence rdf:resource="{uri}_seq"/>',
        '    <sbol:hasFeature>',
        f'      <sbol:SubComponent rdf:about="{uri}_f0">',
        '        <sbol:displayId>CDS_Frame_0</sbol:displayId>',
        '        <sbol:role rdf:resource="https://identifiers.org/SO:0000316"/>',
        '      </sbol:SubComponent>',
        '    </sbol:hasFeature>',
        '    <sbol:hasFeature>',
        f'      <sbol:SubComponent rdf:about="{uri}_f1">',
        '        <sbol:displayId>CDS_Frame_1</sbol:displayId>',
        '        <sbol:role rdf:resource="https://identifiers.org/SO:0000316"/>',
        '      </sbol:SubComponent>',
        '    </sbol:hasFeature>',
        '  </sbol:Component>',
        f'  <sbol:Sequence rdf:about="{uri}_seq">',
        f'    <sbol:displayId>{display_id}_seq</sbol:displayId>',
        f'    <sbol:elements>{dna_seq}</sbol:elements>',
        '    <sbol:encoding rdf:resource="https://identifiers.org/edam:format_1207"/>',
        '  </sbol:Sequence>',
        '</rdf:RDF>',
    ]
    return "\n".join(xml_lines)


def export_to_fasta(
    res: RecompilationResult,
    name: str = "synth_dual_vector",
) -> str:
    """Export DNA and translated proteins into a multi-FASTA document."""
    lines = [
        f">{name}_DNA length={res.total_length_bp}bp compression={res.compression_ratio}x cai={res.codon_adaptation_index}",
        res.synthesized_dna,
        f">{name}_Protein_Frame_0 identity={res.f0_identity_pct}%",
        res.translated_f0,
        f">{name}_Protein_Frame_1 identity={res.f1_identity_pct}% blosum62={res.blosum62_similarity_pct}%",
        res.translated_f1,
    ]
    return "\n".join(lines)


def export_to_csv_summary(res: RecompilationResult) -> str:
    """Export a residue-by-residue alignment table in CSV format."""
    output = io.StringIO()
    writer = csv.writer(output)

    writer.writerow([
        "Residue_Index",
        "Frame0_Target_AA",
        "Frame0_Actual_AA",
        "Frame0_Codon",
        "Frame1_Target_AA",
        "Frame1_Actual_AA",
        "Frame1_Codon",
        "Frame1_Match_Status",
    ])

    min_len = len(res.target_protein_f0)
    dna = res.synthesized_dna

    for i in range(min_len):
        f0_tgt = res.target_protein_f0[i]
        f0_act = res.translated_f0[i] if i < len(res.translated_f0) else "-"
        f0_codon = dna[i * 3 : (i + 1) * 3]

        f1_tgt = res.target_protein_f1[i]
        f1_act = res.translated_f1[i] if i < len(res.translated_f1) else "-"
        f1_codon = dna[i * 3 + 1 : (i + 1) * 3 + 1] if (i * 3 + 4) <= len(dna) else dna[i * 3 + 1 :]

        status = "EXACT_MATCH" if f1_tgt == f1_act else "CONSERVATIVE_MUT"

        writer.writerow([
            i + 1,
            f0_tgt,
            f0_act,
            f0_codon,
            f1_tgt,
            f1_act,
            f1_codon,
            status,
        ])

    return output.getvalue()
