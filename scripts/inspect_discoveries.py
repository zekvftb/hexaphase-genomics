import json
from pathlib import Path

findings_path = Path("D:/DNA/outputs/unaccounted_subroutines_findings.json")
data = json.loads(findings_path.read_text(encoding="utf-8"))

print("==========================================================================")
print("[SMC LAB] DEEP ANOMALY & NOVEL BIOLOGICAL DISCOVERY REPORT")
print("==========================================================================")

for genome_id, record in data["records"].items():
    if genome_id == "NC_001422.1":
        continue  # Skip duplicate
    print(f"\n[GENOME] {genome_id.upper()} (Length: {record['genome_length_bp']} bp)")
    print(f"   Canonical Major Genes: {record['canonical_major_genes']} | Discovered Subroutines: {record['unaccounted_subroutines_count']}")
    
    rbs_subs = [s for s in record["subroutines"] if s["upstream_shine_dalgarno"]]
    print(f"   Hardware RBS-Verified Subroutines: {len(rbs_subs)}")
    
    # 1. Look for high-confidence novel subroutines
    for s in rbs_subs[:5]:
        print(f"\n   [*] ID: {s['subroutine_id']}")
        print(f"       Category: {s['category']} | Strand: {s['strand']} | Frame: {s['frame']}")
        print(f"       Coordinates: {s['start']} -> {s['end']} (Length: {s['length_bp']} bp / {s['length_aa']} AA)")
        print(f"       Hardware RBS Motif: '{s['upstream_rbs_motif']}' (4-15bp upstream of ATG)")
        print(f"       Overlapping Parent Gene: {s['overlapping_parent_id']}")
        print(f"       Downstream Hairpin Gate: {s['hairpin_terminator_downstream']}")
        print(f"       Peptide: {s['protein_sequence']}")
        print(f"       Notes: {s['notes']}")
