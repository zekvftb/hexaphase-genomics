"""CLI Tool: Analyze 3.8-Billion-Year-Old Primordial Archaea Genomes."""

from __future__ import annotations

import json
from pathlib import Path
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from bio_arch.modules.logic_gates import scan_all_logic_gates
from bio_arch.modules.universal_compiler import analyze_wobble_synchronization
from bio_arch.modules.ingestion import stream_fasta, ValidationReport

def read_fasta_simple(file_path: Path) -> list[dict[str, str]]:
    report = ValidationReport(format_detected="fasta")
    records = list(stream_fasta(file_path, report=report, alphabet_type="DNA_EXTENDED"))
    return records

def scan_primordial_archaea():
    out_dir = Path("D:/DNA/outputs/primordial_kernel")
    out_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = Path("D:/DNA/data/primordial_archaea/primordial_archaea_cohort.fasta")
    if not fasta_path.exists():
        print(f"File not found: {fasta_path}")
        return

    print("=======================================================================================================")
    print("🧬 Primordial Archaea OS Kernel Census: 3.8-Billion-Year Invariance Audit")
    print("=======================================================================================================\n")

    records = read_fasta_simple(fasta_path)
    print(f"Loaded {len(records)} ancient Archaea genomes ({sum(len(r['sequence']) for r in records):,} bp total).\n")

    all_reports = []
    total_fs = 0
    total_g4 = 0
    total_rt = 0
    total_bp = 0

    print(f"{'Archaea Species':<32} | {'Length (bp)':<11} | {'Frameshift':<12} | {'G4 Transistors':<14} | {'Readthrough':<11} | {'Density (/kb)'}")
    print("-" * 105)

    for rec in records:
        species_id = rec["id"]
        seq = rec["sequence"]
        total_bp += len(seq)

        # Logic Gate audit
        report = scan_all_logic_gates(seq, genome_id=species_id)
        
        # Wobble entropy analysis
        wobble_report = analyze_wobble_synchronization(seq)

        fs_count = report.gate_counts_by_type.get("frameshift_branch", 0)
        g4_count = report.gate_counts_by_type.get("g4_circuit_breaker", 0)
        rt_count = report.gate_counts_by_type.get("readthrough_overflow", 0)
        density = report.summary.get("computational_density_per_kb", 0.0)

        total_fs += fs_count
        total_g4 += g4_count
        total_rt += rt_count

        all_reports.append({
            "species_id": species_id,
            "length_bp": len(seq),
            "gate_counts": report.gate_counts_by_type,
            "computational_density_per_kb": density,
            "wobble_information": wobble_report.to_dict(),
        })

        print(f"{species_id:<32} | {len(seq):<11} | {fs_count:<12} | {g4_count:<14} | {rt_count:<11} | {density:.3f}")

    print("-" * 105)
    mean_density = ((total_fs + total_g4 + total_rt) / total_bp) * 1000.0
    print(f"{'PRIMORDIAL TOTALS':<32} | {total_bp:<11} | {total_fs:<12} | {total_g4:<14} | {total_rt:<11} | {mean_density:.3f}\n")

    output_data = {
        "audit_title": "HexaPhase Primordial Archaea OS Kernel & Evolutionary Invariance Audit",
        "evolutionary_era": "Archean Eon (~3.8 to ~3.5 Billion Years Ago)",
        "total_archaea_genomes_scanned": len(records),
        "total_base_pairs_scanned": total_bp,
        "summary_totals": {
            "total_frameshift_multiplexers": total_fs,
            "total_g4_circuit_breakers": total_g4,
            "total_readthrough_gates": total_rt,
            "grand_total_hardware_switches": total_fs + total_g4 + total_rt,
            "primordial_computational_density_per_kb": round(mean_density, 3),
        },
        "species_reports": all_reports,
    }

    out_file = out_dir / "primordial_archaea_audit.json"
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(output_data, f, indent=2)

    print(f"✅ Certified 3.8-Billion-Year Primordial Audit saved to: {out_file}")

if __name__ == "__main__":
    scan_primordial_archaea()
