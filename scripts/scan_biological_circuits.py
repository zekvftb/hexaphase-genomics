"""Global Census: Biological Memory & Security Architecture (Module 8).

Audits:
1. Human Cancer Oncogenes for CpG Epigenetic Memory Registers (NVRAM)
2. Primordial Archaea Cohort for CRISPR Append-Only Security Arrays
3. Pan-Pathogen Cohort for RNA Riboswitch Chemical ADCs
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from bio_arch.modules.biological_circuits import scan_all_biological_circuits
from bio_arch.provenance import now_iso


def parse_fasta_records(filepath: Path) -> list[tuple[str, str]]:
    """Parse multi-FASTA file into list of (header, sequence) tuples."""
    if not filepath.exists():
        return []
    records = []
    current_header = ""
    current_seq: list[str] = []

    for line in filepath.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header and current_seq:
                records.append((current_header, "".join(current_seq)))
            current_header = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line)

    if current_header and current_seq:
        records.append((current_header, "".join(current_seq)))

    return records


def run_biological_circuits_census():
    print("=" * 80)
    print("💾 MODULE 8 CENSUS: BIOLOGICAL MEMORY & SECURITY ARCHITECTURE")
    print("   Auditing: CRISPR Security Logs, CpG NVRAM Registers, and Riboswitch ADCs")
    print("=" * 80)

    data_dir = Path(__file__).parent.parent / "data"
    outputs_dir = Path(__file__).parent.parent / "outputs"
    outputs_dir.mkdir(parents=True, exist_ok=True)

    all_reports = []
    total_crispr_arrays = 0
    total_spacers = 0
    total_cpg_islands = 0
    total_promoter_latches = 0
    total_riboswitches = 0
    total_bp_scanned = 0

    # 1. Human Cancer Oncogenes (CpG Epigenetic Registers)
    oncogene_path = data_dir / "human_cancer_oncogenes.fasta"
    if oncogene_path.exists():
        records = parse_fasta_records(oncogene_path)
        print(f"\n🧬 [Cohort 1] Human Cancer Oncogenes ({len(records)} targets)...")
        for header, seq in records:
            name = header.split()[0]
            rep = scan_all_biological_circuits(seq, genome_id=name)
            all_reports.append(rep.to_dict())
            total_bp_scanned += len(seq)
            total_cpg_islands += len(rep.cpg_memory_islands)
            total_promoter_latches += rep.summary.get("promoter_memory_latches", 0)
            total_riboswitches += len(rep.riboswitch_adcs)
            total_crispr_arrays += len(rep.crispr_arrays)
            total_spacers += rep.summary.get("total_spacers_archived", 0)

            print(f"   • {name:<12} ({len(seq):>6,} bp): "
                  f"{len(rep.cpg_memory_islands)} CpG NVRAM Registers "
                  f"({rep.summary.get('promoter_memory_latches', 0)} Promoter Latches)")

    # 2. Primordial Archaea (CRISPR Arrays & Security Logs)
    archaea_dir = data_dir / "primordial_archaea"
    if archaea_dir.exists():
        archaea_files = list(archaea_dir.glob("*.fasta"))
        print(f"\n🌋 [Cohort 2] Primordial Archaea ({len(archaea_files)} ancient genomes)...")
        for fpath in archaea_files:
            records = parse_fasta_records(fpath)
            for header, seq in records:
                name = fpath.stem.upper()
                rep = scan_all_biological_circuits(seq, genome_id=name)
                all_reports.append(rep.to_dict())
                total_bp_scanned += len(seq)
                total_cpg_islands += len(rep.cpg_memory_islands)
                total_promoter_latches += rep.summary.get("promoter_memory_latches", 0)
                total_riboswitches += len(rep.riboswitch_adcs)
                total_crispr_arrays += len(rep.crispr_arrays)
                total_spacers += rep.summary.get("total_spacers_archived", 0)

                print(f"   • {name:<20} ({len(seq):>10,} bp): "
                      f"{len(rep.crispr_arrays)} CRISPR Security Arrays "
                      f"({rep.summary.get('total_spacers_archived', 0)} Viral Spacers Logged), "
                      f"{len(rep.riboswitch_adcs)} Riboswitch ADCs")

    # 3. Pan-Pathogen Cohort
    pathogen_path = data_dir / "pan_pathogen_cohort.fasta"
    if pathogen_path.exists():
        records = parse_fasta_records(pathogen_path)
        print(f"\n🦠 [Cohort 3] Pan-Pathogen Reference Cohort ({len(records)} genomes)...")
        for header, seq in records:
            name = header.split()[0]
            rep = scan_all_biological_circuits(seq, genome_id=name)
            all_reports.append(rep.to_dict())
            total_bp_scanned += len(seq)
            total_cpg_islands += len(rep.cpg_memory_islands)
            total_promoter_latches += rep.summary.get("promoter_memory_latches", 0)
            total_riboswitches += len(rep.riboswitch_adcs)
            total_crispr_arrays += len(rep.crispr_arrays)
            total_spacers += rep.summary.get("total_spacers_archived", 0)

    # Export comprehensive audit
    out_file = outputs_dir / "biological_circuits_audit.json"
    audit_data = {
        "status": "PASSED",
        "timestamp": now_iso(),
        "audit_version": "1.0.0-MODULE-8",
        "total_bp_scanned": total_bp_scanned,
        "global_summary": {
            "total_cpg_memory_islands": total_cpg_islands,
            "total_promoter_memory_latches": total_promoter_latches,
            "total_crispr_security_arrays": total_crispr_arrays,
            "total_viral_spacers_archived": total_spacers,
            "total_riboswitch_adcs": total_riboswitches,
        },
        "reports": all_reports,
    }
    out_file.write_text(json.dumps(audit_data, indent=2), encoding="utf-8")
    print(f"\n📁 Census report saved to {out_file}")

    print("=" * 80)
    print(f"🎯 CENSUS SUMMARY ACROSS {total_bp_scanned:,} BP:")
    print(f"   💾 CpG Epigenetic Memory Registers: {total_cpg_islands:,} ({total_promoter_latches:,} Active Latches)")
    print(f"   🛡️ CRISPR Security Arrays:          {total_crispr_arrays:,} ({total_spacers:,} Viral Spacers Logged)")
    print(f"   🎛️ RNA Riboswitch Chemical ADCs:    {total_riboswitches:,}")
    print("=" * 80)


if __name__ == "__main__":
    run_biological_circuits_census()
