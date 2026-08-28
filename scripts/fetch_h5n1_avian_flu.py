"""Fetch complete 8-segment genome of Highly Pathogenic Avian Influenza A (H5N1)."""

import urllib.request
import ssl
import time
from pathlib import Path

data_dir = Path("D:/DNA/data")
data_dir.mkdir(parents=True, exist_ok=True)

# Highly Pathogenic Avian Influenza A (H5N1) reference segments (A/Goose/Guangdong/1/96)
H5N1_SEGMENTS = [
    ("H5N1_Seg1_PB2", "NC_007357.1", "Influenza A virus (A/Goose/Guangdong/1/1996(H5N1)) segment 1 PB2"),
    ("H5N1_Seg2_PB1", "NC_007358.1", "Influenza A virus (A/Goose/Guangdong/1/1996(H5N1)) segment 2 PB1"),
    ("H5N1_Seg3_PA",  "NC_007359.1", "Influenza A virus (A/Goose/Guangdong/1/1996(H5N1)) segment 3 PA"),
    ("H5N1_Seg4_HA",  "NC_007360.1", "Influenza A virus (A/Goose/Guangdong/1/1996(H5N1)) segment 4 HA"),
    ("H5N1_Seg5_NP",  "NC_007361.1", "Influenza A virus (A/Goose/Guangdong/1/1996(H5N1)) segment 5 NP"),
    ("H5N1_Seg6_NA",  "NC_007362.1", "Influenza A virus (A/Goose/Guangdong/1/1996(H5N1)) segment 6 NA"),
    ("H5N1_Seg7_M",   "NC_007363.1", "Influenza A virus (A/Goose/Guangdong/1/1996(H5N1)) segment 7 M1/M2"),
    ("H5N1_Seg8_NS",  "NC_007364.1", "Influenza A virus (A/Goose/Guangdong/1/1996(H5N1)) segment 8 NS1/NEP"),
]

def fetch_ncbi_fasta(accession: str) -> str:
    ctx = ssl._create_unverified_context()
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"})
    with urllib.request.urlopen(req, context=ctx, timeout=15) as resp:
        return resp.read().decode("utf-8")

def build_h5n1_panel():
    out_file = data_dir / "h5n1_avian_flu_complete.fasta"
    print(f"Fetching all {len(H5N1_SEGMENTS)} segments of Avian Influenza A (H5N1) from NCBI...")
    
    with open(out_file, "w", encoding="utf-8") as out:
        for tag, acc, desc in H5N1_SEGMENTS:
            try:
                fasta_data = fetch_ncbi_fasta(acc)
                if fasta_data.startswith(">"):
                    lines = fasta_data.strip().split("\n")
                    clean_header = f">{tag} {desc} [{acc}]"
                    seq_content = "\n".join(lines[1:])
                    out.write(f"{clean_header}\n{seq_content}\n")
                    seq_len = len("".join(lines[1:]))
                    print(f"  [DOWNLOADED] {tag:<15} ({acc}) -> {seq_len:,} bp")
                else:
                    print(f"  [FAILED] {tag}: Invalid response")
                time.sleep(0.3)
            except Exception as e:
                print(f"  [ERROR] {tag} ({acc}): {e}")

    print(f"\nSaved complete H5N1 Bird Flu Genome to: {out_file}")

if __name__ == "__main__":
    build_h5n1_panel()
