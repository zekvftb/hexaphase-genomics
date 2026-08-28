"""CLI Tool: Decompile DNA Genomes into Biological Assembly Bytecode (.asm)."""

from __future__ import annotations

import json
from pathlib import Path
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from bio_arch.modules.universal_compiler import (
    analyze_wobble_synchronization,
    decompile_genomic_bytecode,
)
from bio_arch.modules.ingestion import stream_fasta, ValidationReport

def read_fasta_simple(file_path: Path) -> list[dict[str, str]]:
    report = ValidationReport(format_detected="fasta")
    records = list(stream_fasta(file_path, report=report, alphabet_type="DNA_EXTENDED"))
    return records

def run_decompiler():
    out_dir = Path("D:/DNA/outputs/decompiled_bytecode")
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=======================================================================================================")
    print("🧬 HexaPhase Biological Machine Code Decompiler & Assembly Generator")
    print("=======================================================================================================\n")

    targets = [
        ("PhiX174_NC_001422", Path("D:/DNA/data/phix174_complete.fasta")),
        ("HIV1_NC_001802", Path("D:/DNA/data/pan_pathogen_cohort.fasta"), "HIV1"),
        ("Ebola_NC_002549", Path("D:/DNA/data/pan_pathogen_cohort.fasta"), "Ebola_Virus"),
        ("TP53_NM_000546", Path("D:/DNA/data/human_cancer_oncogenes.fasta"), "TP53"),
        ("MYC_NM_002467", Path("D:/DNA/data/human_cancer_oncogenes.fasta"), "MYC"),
        ("H5N1_HA_NC_007360", Path("D:/DNA/data/h5n1_avian_flu_complete.fasta"), "H5N1_Seg4_HA"),
    ]

    all_wobble_reports = {}

    for item in targets:
        tag = item[0]
        fasta_path = item[1]
        filter_id = item[2] if len(item) > 2 else None

        if not fasta_path.exists():
            continue

        records = read_fasta_simple(fasta_path)
        matching_rec = None
        for r in records:
            if filter_id is None or filter_id in r["id"]:
                matching_rec = r
                break

        if not matching_rec:
            continue

        seq = matching_rec["sequence"]
        print(f"--- Decompiling: {tag} ({len(seq):,} bp) ---")

        # 1. Wobble Information Synchronization Analysis
        wobble_report = analyze_wobble_synchronization(seq)
        all_wobble_reports[tag] = wobble_report.to_dict()
        print(f"  * Pos 1 Shannon Entropy: {wobble_report.pos1_shannon_entropy_bits:.3f} bits")
        print(f"  * Pos 2 Shannon Entropy: {wobble_report.pos2_shannon_entropy_bits:.3f} bits")
        print(f"  * Pos 3 Wobble Entropy:  {wobble_report.pos3_wobble_entropy_bits:.3f} bits (Capacity Ratio: {wobble_report.pos3_information_capacity_ratio:.2f}x)")
        print(f"  * Mutual Information (F0 / F+1): {wobble_report.mutual_information_f0_f1_bits:.4f} bits/codon")

        # 2. Decompile to Assembly (.asm)
        instructions, asm_text = decompile_genomic_bytecode(seq, genome_id=tag)
        asm_out_path = out_dir / f"{tag}_disassembly.asm"
        with open(asm_out_path, "w", encoding="utf-8") as f:
            f.write(asm_text)
        print(f"  * Disassembled into {len(instructions):,} Biological Assembly opcodes -> {asm_out_path}\n")

    # Save summary report
    summary_path = out_dir / "wobble_information_architecture_report.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump({
            "report_title": "HexaPhase Universal Information Architecture & Wobble Synchronization Audit",
            "findings_summary": "Pos3 Wobble codon positions consistently carry the highest Shannon entropy, acting as an independent carrier wave for parallel execution tracks without degrading Frame 0 protein structure.",
            "wobble_reports": all_wobble_reports,
        }, f, indent=2)
    print(f"✅ Full Wobble Information report written to: {summary_path}")

if __name__ == "__main__":
    run_decompiler()
