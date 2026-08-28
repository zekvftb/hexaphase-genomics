"""Fetch master human cancer oncogenes and tumor suppressors from NCBI RefSeq."""

import urllib.request
import ssl
import time
from pathlib import Path

data_dir = Path("D:/DNA/data")
data_dir.mkdir(parents=True, exist_ok=True)

CANCER_GENES = [
    ("TP53", "NM_000546.6", "Homo sapiens tumor protein p53 (TP53), transcript variant 1, mRNA"),
    ("MYC", "NM_002467.6", "Homo sapiens MYC proto-oncogene, bHLH transcription factor (MYC), mRNA"),
    ("KRAS", "NM_004985.5", "Homo sapiens KRAS proto-oncogene, GTPase (KRAS), transcript variant a, mRNA"),
    ("EGFR", "NM_005228.5", "Homo sapiens epidermal growth factor receptor (EGFR), transcript variant 1, mRNA"),
    ("BRCA1", "NM_007294.4", "Homo sapiens BRCA1 DNA repair associated (BRCA1), transcript variant 1, mRNA"),
    ("PTEN", "NM_000314.8", "Homo sapiens phosphatase and tensin homolog (PTEN), transcript variant 1, mRNA"),
    ("BRAF", "NM_004333.6", "Homo sapiens B-Raf proto-oncogene, serine/threonine kinase (BRAF), transcript variant 1, mRNA"),
]

def fetch_ncbi_fasta(accession: str) -> str:
    ctx = ssl._create_unverified_context()
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"})
    with urllib.request.urlopen(req, context=ctx, timeout=15) as resp:
        return resp.read().decode("utf-8")

def build_cancer_panel():
    out_file = data_dir / "human_cancer_oncogenes.fasta"
    print(f"Fetching {len(CANCER_GENES)} master human cancer genes from NCBI RefSeq...")
    
    with open(out_file, "w", encoding="utf-8") as out:
        for tag, acc, desc in CANCER_GENES:
            try:
                fasta_data = fetch_ncbi_fasta(acc)
                if fasta_data.startswith(">"):
                    lines = fasta_data.strip().split("\n")
                    clean_header = f">{tag} {desc} [{acc}]"
                    seq_content = "\n".join(lines[1:])
                    out.write(f"{clean_header}\n{seq_content}\n")
                    seq_len = len("".join(lines[1:]))
                    print(f"  [DOWNLOADED] {tag:<10} ({acc}) -> {seq_len:,} bp")
                else:
                    print(f"  [FAILED] {tag}: Invalid response")
                time.sleep(0.3)
            except Exception as e:
                print(f"  [ERROR] {tag} ({acc}): {e}")

    print(f"\nSaved Human Cancer Oncogene Panel to: {out_file}")

if __name__ == "__main__":
    build_cancer_panel()
