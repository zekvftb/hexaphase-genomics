"""CLI Runner: Scan and characterize Biological Logic Gates across viral and host genomes."""

from __future__ import annotations

import json
from pathlib import Path
import sys

# Ensure UTF-8 stdout for Windows consoles
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
    
    fasta_path = Path("D:/DNA/data/study_cohort/study_cohort.fasta")
    if not fasta_path.exists():
        print(f"File not found: {fasta_path}")
        return

    print("==========================================================================================")
    print("HexaPhase Biological Logic Gate & Hardware Switch Scanner")
    print("==========================================================================================")
    
    # Read multi-fasta
    records = read_fasta_simple(fasta_path)
    print(f"Loaded {len(records)} reference genomes for hardware gate mining.\n")
    
    all_reports = []
    total_fs = 0
    total_g4 = 0
    total_rt = 0
    
    print(f"{'Genome ID':<22} | {'Length (bp)':<11} | {'Frameshift Gates':<16} | {'G4 Transistors':<14} | {'Readthrough':<11} | {'Density (/kb)'}")
    print("-" * 95)
    
    for rec in records:
        genome_id = rec["id"]
        seq = rec["sequence"]
        
        report = scan_all_logic_gates(seq, genome_id=genome_id)
        all_reports.append(report.to_dict())
        
        fs_count = report.gate_counts_by_type.get("frameshift_branch", 0)
        g4_count = report.gate_counts_by_type.get("g4_circuit_breaker", 0)
        rt_count = report.gate_counts_by_type.get("readthrough_overflow", 0)
        density = report.summary.get("computational_density_per_kb", 0.0)
        
        total_fs += fs_count
        total_g4 += g4_count
        total_rt += rt_count
        
        print(f"{genome_id:<22} | {len(seq):<11} | {fs_count:<16} | {g4_count:<14} | {rt_count:<11} | {density:.3f}")

    print("-" * 95)
    print(f"{'TOTALS':<22} | {'--':<11} | {total_fs:<16} | {total_g4:<14} | {total_rt:<11} | --\n")
    
    # Save combined report
    combined_output_path = outputs_dir / "biological_logic_gates_audit.json"
    with open(combined_output_path, "w", encoding="utf-8") as f:
        json.dump({
            "audit_title": "HexaPhase Biological Logic Gates and Hardware Execution Audit",
            "total_genomes_scanned": len(records),
            "summary_totals": {
                "total_frameshift_multiplexers": total_fs,
                "total_g4_circuit_breakers": total_g4,
                "total_readthrough_gates": total_rt,
                "grand_total_logic_switches": total_fs + total_g4 + total_rt,
            },
            "genome_reports": all_reports,
        }, f, indent=2)
        
    print(f"Full certified logic gate audit written to: {combined_output_path}")

if __name__ == "__main__":
    run_logic_gates_scan()
