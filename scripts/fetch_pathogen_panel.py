"""Fetch major viral pathogens from NCBI for Logic Gate analysis."""

import urllib.request
import ssl
import time
from pathlib import Path

data_dir = Path("D:/DNA/data")
data_dir.mkdir(parents=True, exist_ok=True)

TARGETS = [
    ("HIV1", "NC_001802.1", "Human Immunodeficiency Virus 1 (HIV-1)"),
    ("Ebola_Virus", "NC_002549.1", "Zaire ebolavirus complete genome"),
    ("RSV", "NC_001803.1", "Human respiratory syncytial virus"),
    ("Dengue_Virus_2", "NC_001474.2", "Dengue virus type 2 complete genome"),
    ("Zika_Virus", "NC_012532.1", "Zika virus isolate complete genome"),
    ("Hepatitis_C", "NC_004102.1", "Hepatitis C virus genotype 1"),
    ("Mpox_Virus", "NC_063383.1", "Monkeypox virus complete genome"),
    ("Influenza_A_PB2", "NC_026433.1", "Influenza A virus PB2 segment"),
    ("Influenza_A_PB1", "NC_026434.1", "Influenza A virus PB1 segment"),
    ("Influenza_A_HA", "NC_026436.1", "Influenza A virus Hemagglutinin"),
]

def fetch_ncbi_fasta(accession: str) -> str:
    ctx = ssl._create_unverified_context()
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"})
    with urllib.request.urlopen(req, context=ctx, timeout=15) as resp:
        return resp.read().decode("utf-8")

def build_pathogen_panel():
    out_file = data_dir / "pan_pathogen_cohort.fasta"
    print(f"Fetching {len(TARGETS)} major viral pathogens from NCBI...")
    
    with open(out_file, "w", encoding="utf-8") as out:
        for tag, acc, desc in TARGETS:
            try:
                fasta_data = fetch_ncbi_fasta(acc)
                if fasta_data.startswith(">"):
                    lines = fasta_data.strip().split("\n")
                    clean_header = f">{tag} {desc} [{acc}]"
                    seq_content = "\n".join(lines[1:])
                    out.write(f"{clean_header}\n{seq_content}\n")
                    seq_len = len("".join(lines[1:]))
                    print(f"  [DOWNLOADED] {tag:<18} ({acc}) -> {seq_len:,} bp")
                else:
                    print(f"  [FAILED] {tag}: Invalid response")
                time.sleep(0.3) # Friendly rate limit
            except Exception as e:
                print(f"  [ERROR] {tag} ({acc}): {e}")

    print(f"\nSaved complete Pan-Pathogen panel to: {out_file}")

if __name__ == "__main__":
    build_pathogen_panel()
