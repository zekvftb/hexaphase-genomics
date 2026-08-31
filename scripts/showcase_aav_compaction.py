"""End-to-End Benchmark Showcase: Dual-Coding Vector Compaction for AAV Packaging.

Demonstrates compaction of two canonical fluorescent reporter proteins (GFP in Frame 0,
mCherry in Frame +1) into a single overlapping DNA construct, evaluates CAI, BLOSUM62
conservation, restriction site avoidance, and exports to annotated GenBank and SBOL3.
"""

from __future__ import annotations

from pathlib import Path
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from bio_arch.modules.export_formats import (
    export_to_fasta,
    export_to_genbank,
    export_to_sbol3,
)
from bio_arch.modules.recompiler import (
    calculate_aav_packaging_savings,
    recompile_dual_protein_dna,
)

# Canonical GFP (238 aa)
GFP_SEQUENCE = (
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCF"
    "SRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEY"
    "NYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKR"
    "DHMVLLEFVTAAGITHGMDELYK"
)

# Canonical mCherry (236 aa)
MCHERRY_SEQUENCE = (
    "MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQF"
    "MYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVM"
    "QKKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNED"
    "YTIVEQYERAEGRHSTGGMDELYK"
)


def run_showcase():
    print("=" * 80)
    print("🧬 HEXAPHASE GENOMICS: DUAL-CODING AAV VECTOR COMPACTION SHOWCASE")
    print("=" * 80)

    print(f"\n[1] Input Target Peptides:")
    print(f"    - Frame 0: GFP ({len(GFP_SEQUENCE)} aa)")
    print(f"    - Frame +1: mCherry ({len(MCHERRY_SEQUENCE)} aa)")

    print("\n[2] Executing Trellis Dynamic Programming Recompiler...")
    res = recompile_dual_protein_dna(
        protein_f0=GFP_SEQUENCE,
        protein_f1=MCHERRY_SEQUENCE,
        optimize_cai=True,
        filter_restriction=True,
        allow_conservative_mutations=True,
    )

    print("\n[3] Synthesis Results:")
    print(f"    - Overlapping DNA Construct Length: {res.total_length_bp} bp")
    print(f"    - Frame 0 (GFP) Identity:           {res.f0_identity_pct}%")
    print(f"    - Frame +1 (mCherry) Identity:      {res.f1_identity_pct}%")
    print(f"    - BLOSUM62 Conservation Score:      {res.blosum62_similarity_pct}%")
    print(f"    - Mammalian Codon Adaptation (CAI): {res.codon_adaptation_index}")
    print(f"    - Physical Footprint Compression:   {res.compression_ratio}x")
    print(f"    - Restriction Sites Detected:       {res.restriction_sites_detected or 'None (Clean)'}")
    print(f"    - Homopolymer Runs:                 {res.homopolymer_runs_detected}")

    print("\n[4] AAV Capsid Packaging Capacity Analysis (4.7 kb Limit):")
    savings = calculate_aav_packaging_savings(
        protein_len_f0_aa=len(GFP_SEQUENCE),
        protein_len_f1_aa=len(MCHERRY_SEQUENCE),
        synthesized_dna_bp=res.total_length_bp,
        aav_limit_bp=4700,
    )
    print(f"    - Separate Dual-Gene Footprint:     {savings['separate_dual_cassette_total_bp']:,} bp")
    print(f"    - Dual-Phase Compact Footprint:     {savings['compact_overlapping_cassette_total_bp']:,} bp")
    print(f"    - Physical Space Saved:             {savings['bp_saved']:,} bp ({savings['percent_footprint_reduction']}%)")
    print(f"    - Fits Standard AAV Packaging:      {savings['fits_standard_aav']}")
    print(f"    - Extra AAV Capsid Headroom:        {savings['headroom_remaining_bp']:,} bp")

    # Exporting files to examples/
    examples_dir = Path("examples")
    examples_dir.mkdir(parents=True, exist_ok=True)

    gbk_path = examples_dir / "gfp_mcherry_dual_coding.gbk"
    gbk_content = export_to_genbank(
        res,
        locus_name="GFP_MCHERRY_AAV",
        definition="Dual-coding GFP (F0) and mCherry (F1) compact AAV expression vector",
    )
    gbk_path.write_text(gbk_content, encoding="utf-8")
    print(f"\n[5] Exported Annotated GenBank: {gbk_path.resolve()}")

    sbol_path = examples_dir / "gfp_mcherry_dual_coding.sbol3"
    sbol_content = export_to_sbol3(res, display_id="gfp_mcherry_dual_coding")
    sbol_path.write_text(sbol_content, encoding="utf-8")
    print(f"    Exported SBOL3 XML:         {sbol_path.resolve()}")

    fasta_path = examples_dir / "gfp_mcherry_dual_coding.fasta"
    fasta_content = export_to_fasta(res, name="gfp_mcherry_dual_coding")
    fasta_path.write_text(fasta_content, encoding="utf-8")
    print(f"    Exported Multi-FASTA:       {fasta_path.resolve()}")

    print("\n" + "=" * 80)
    print("✅ SHOWCASE RUN COMPLETED SUCCESSFULLY")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    run_showcase()
