"""CLI Runner: Scan and characterize Biological Logic Gates across viral and host genomes."""

from __future__ import annotations

import json
from pathlib import Path
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from bio_arch.modules.logic_gates import scan_all_logic_gates
from bio_arch.modules.ingestion import stream_fasta, ValidationReport

def read_fasta_simple(file_path: Path) -> list[dict[str, str]]:
    report = ValidationReport(format_detected="fasta")
    records = list(stream_fasta(file_path, report=report, alphabet_type="DNA_EXTENDED"))
    return records

def run_logic_gates_scan():
    outputs_dir = Path("D:/DNA/outputs/logic_gates")
    outputs_dir.mkdir(parents=True, exist_ok=True)
    
    cohort_files = [
        ("Reference Cohort", Path("D:/DNA/data/study_cohort/study_cohort.fasta")),
        ("Global Pathogen Panel", Path("D:/DNA/data/pan_pathogen_cohort.fasta")),
    ]

    print("=======================================================================================================")
    print("🧬 HexaPhase Biological Logic Gate & Hardware Switch Global Audit")
    print("=======================================================================================================\n")
    
    all_reports = []
    total_fs = 0
    total_g4 = 0
    total_rt = 0
    total_bp = 0
    
    for cohort_name, fasta_path in cohort_files:
        if not fasta_path.exists():
            continue
            
        records = read_fasta_simple(fasta_path)
        print(f"--- Cohort: {cohort_name} ({len(records)} genomes) ---")
        print(f"{'Genome ID':<24} | {'Length (bp)':<11} | {'Frameshift':<12} | {'G4 Transistor':<14} | {'Readthrough':<11} | {'Density (/kb)'}")
        print("-" * 96)
        
        for rec in records:
            genome_id = rec["id"]
            seq = rec["sequence"]
            total_bp += len(seq)
            
            report = scan_all_logic_gates(seq, genome_id=genome_id)
            all_reports.append(report.to_dict())
            
            fs_count = report.gate_counts_by_type.get("frameshift_branch", 0)
            g4_count = report.gate_counts_by_type.get("g4_circuit_breaker", 0)
            rt_count = report.gate_counts_by_type.get("readthrough_overflow", 0)
            density = report.summary.get("computational_density_per_kb", 0.0)
            
            total_fs += fs_count
            total_g4 += g4_count
            total_rt += rt_count
            
            print(f"{genome_id:<24} | {len(seq):<11} | {fs_count:<12} | {g4_count:<14} | {rt_count:<11} | {density:.3f}")
        print("\n")

    print("=" * 96)
    print(f"GRAND TOTALS across {len(all_reports)} Genomes ({total_bp:,} bp analyzed):")
    print(f"  * Total Programmed Frameshifting Multiplexers: {total_fs:,}")
    print(f"  * Total G4 Molecular Transistors:            {total_g4:,}")
    print(f"  * Total Leaky Stop Readthrough Gates:        {total_rt:,}")
    print(f"  * Total Biological Hardware Gates Discovered: {total_fs + total_g4 + total_rt:,}")
    print(f"  * Mean Genomic Logic Density:                {((total_fs + total_g4 + total_rt) / total_bp) * 1000.0:.3f} gates per kb")
    print("=" * 96)
    
    # Save combined report
    combined_output_path = outputs_dir / "pan_pathogen_logic_gates_audit.json"
    with open(combined_output_path, "w", encoding="utf-8") as f:
        json.dump({
            "audit_title": "HexaPhase Pan-Pathogen & Reference Biological Logic Gate Audit",
            "total_genomes_scanned": len(all_reports),
            "total_base_pairs_scanned": total_bp,
            "summary_totals": {
                "total_frameshift_multiplexers": total_fs,
                "total_g4_circuit_breakers": total_g4,
                "total_readthrough_gates": total_rt,
                "grand_total_logic_switches": total_fs + total_g4 + total_rt,
                "mean_computational_density_per_kb": round(((total_fs + total_g4 + total_rt) / total_bp) * 1000.0, 3),
            },
            "genome_reports": all_reports,
        }, f, indent=2)
        
    print(f"\n✅ Full certified Pan-Pathogen logic gate audit saved to: {combined_output_path}")

if __name__ == "__main__":
    run_logic_gates_scan()
