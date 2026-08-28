"""Fetch 4 ancient extremophile Archaea genomes (3.8-billion-year-old Primordial OS Kernel)."""

import urllib.request
import ssl
import time
from pathlib import Path

data_dir = Path("D:/DNA/data/primordial_archaea")
data_dir.mkdir(parents=True, exist_ok=True)

ARCHAEA_GENOMES = [
    ("Methanocaldococcus_jannaschii", "NC_000909.1", "Methanocaldococcus jannaschii DSM 2661 complete genome (Hydrothermal Vent)"),
    ("Pyrococcus_furiosus", "NC_003413.1", "Pyrococcus furiosus DSM 3638 complete genome (100°C Boiling Vents)"),
    ("Sulfolobus_solfataricus", "NC_002754.1", "Sulfolobus solfataricus P2 complete genome (Volcanic Hot Springs)"),
    ("Haloferax_volcanii", "NC_013967.1", "Haloferax volcanii DS2 complete chromosome (Dead Sea Hypersaline)"),
]

def fetch_ncbi_fasta(accession: str) -> str:
    ctx = ssl._create_unverified_context()
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"})
    with urllib.request.urlopen(req, context=ctx, timeout=30) as resp:
        return resp.read().decode("utf-8")

def build_archaea_panel():
    out_file = data_dir / "primordial_archaea_cohort.fasta"
    print(f"Fetching {len(ARCHAEA_GENOMES)} ancient Archaea genomes from NCBI RefSeq...")
    
    with open(out_file, "w", encoding="utf-8") as out:
        for tag, acc, desc in ARCHAEA_GENOMES:
            try:
                fasta_data = fetch_ncbi_fasta(acc)
                if fasta_data.startswith(">"):
                    lines = fasta_data.strip().split("\n")
                    clean_header = f">{tag} {desc} [{acc}]"
                    seq_content = "\n".join(lines[1:])
                    out.write(f"{clean_header}\n{seq_content}\n")
                    seq_len = len("".join(lines[1:]))
                    print(f"  [DOWNLOADED] {tag:<30} ({acc}) -> {seq_len:,} bp")
                else:
                    print(f"  [FAILED] {tag}: Invalid response")
                time.sleep(0.5)
            except Exception as e:
                print(f"  [ERROR] {tag} ({acc}): {e}")

    print(f"\nSaved Primordial Archaea Cohort to: {out_file}")

if __name__ == "__main__":
    build_archaea_panel()
