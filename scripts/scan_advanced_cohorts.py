"""Advanced Batch Scanner: Analyze Human Cancer Oncogenes and H5N1 Avian Flu."""

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

def scan_cohort_file(title: str, fasta_path: Path, output_json_path: Path) -> dict:
    if not fasta_path.exists():
        print(f"Error: {fasta_path} does not exist.")
        return {}

    records = read_fasta_simple(fasta_path)
    print(f"\n=======================================================================================================")
    print(f"🧬 {title} ({len(records)} sequences)")
    print(f"=======================================================================================================")
    print(f"{'Target ID':<20} | {'Length (bp)':<11} | {'Frameshift':<12} | {'G4 Transistor':<14} | {'Readthrough':<11} | {'Density (/kb)'}")
    print("-" * 96)

    all_reports = []
    total_fs = 0
    total_g4 = 0
    total_rt = 0
    total_bp = 0

    for rec in records:
        target_id = rec["id"]
        seq = rec["sequence"]
        total_bp += len(seq)

        report = scan_all_logic_gates(seq, genome_id=target_id)
        all_reports.append(report.to_dict())

        fs_count = report.gate_counts_by_type.get("frameshift_branch", 0)
        g4_count = report.gate_counts_by_type.get("g4_circuit_breaker", 0)
        rt_count = report.gate_counts_by_type.get("readthrough_overflow", 0)
        density = report.summary.get("computational_density_per_kb", 0.0)

        total_fs += fs_count
        total_g4 += g4_count
        total_rt += rt_count

        print(f"{target_id:<20} | {len(seq):<11} | {fs_count:<12} | {g4_count:<14} | {rt_count:<11} | {density:.3f}")

    print("-" * 96)
    print(f"{'TOTALS':<20} | {total_bp:<11} | {total_fs:<12} | {total_g4:<14} | {total_rt:<11} | {((total_fs + total_g4 + total_rt) / max(1, total_bp)) * 1000.0:.3f}")

    output_data = {
        "cohort_title": title,
        "total_sequences_scanned": len(records),
        "total_base_pairs": total_bp,
        "summary_totals": {
            "total_frameshift_multiplexers": total_fs,
            "total_g4_circuit_breakers": total_g4,
            "total_readthrough_gates": total_rt,
            "grand_total_hardware_switches": total_fs + total_g4 + total_rt,
            "mean_density_per_kb": round(((total_fs + total_g4 + total_rt) / max(1, total_bp)) * 1000.0, 3),
        },
        "reports": all_reports,
    }

    output_json_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json_path, "w", encoding="utf-8") as f:
        json.dump(output_data, f, indent=2)

    print(f"✅ Certified audit written to: {output_json_path}")
    return output_data

def main():
    outputs_dir = Path("D:/DNA/outputs/logic_gates")
    
    # 1. Human Cancer Oncogene Screen
    cancer_fasta = Path("D:/DNA/data/human_cancer_oncogenes.fasta")
    cancer_out = outputs_dir / "cancer_oncogenes_logic_audit.json"
    scan_cohort_file("Human Cancer Oncogene & Tumor Suppressor Panel", cancer_fasta, cancer_out)

    # 2. H5N1 Avian Flu Screen
    h5n1_fasta = Path("D:/DNA/data/h5n1_avian_flu_complete.fasta")
    h5n1_out = outputs_dir / "h5n1_avian_flu_audit.json"
    scan_cohort_file("Highly Pathogenic Avian Influenza A (H5N1) 8-Segment Genome", h5n1_fasta, h5n1_out)

if __name__ == "__main__":
    main()
