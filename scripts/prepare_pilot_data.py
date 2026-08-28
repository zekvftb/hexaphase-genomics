"""Fetch Bacteriophage PhiX174 (NC_001422.1) and extract functional partitions."""

from pathlib import Path
import urllib.request

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)


def fetch_phix174() -> str:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text"
    req = urllib.request.Request(url, headers={"User-Agent": "BioArchResearch/1.0"})
    with urllib.request.urlopen(req) as resp:
        content = resp.read().decode("utf-8")

    lines = content.strip().splitlines()
    seq_lines = [l.strip() for l in lines if not l.startswith(">")]
    genome_seq = "".join(seq_lines).upper()

    if len(genome_seq) != 5386:
        raise ValueError(f"Expected PhiX174 genome length 5386, got {len(genome_seq)}")

    return genome_seq


def main() -> None:
    print("Fetching Bacteriophage PhiX174 (NC_001422.1) from NCBI...")
    genome = fetch_phix174()
    print(f"Downloaded complete circular genome: {len(genome)} bp")

    # 1. Save complete genome
    genome_file = DATA_DIR / "phix174_complete.fasta"
    genome_file.write_text(f">NC_001422.1 Bacteriophage phiX174 complete genome\n{genome}\n", encoding="utf-8")

    # 2. Extract Preregistered Partitions (1-based to 0-based slices)
    # Partition 1: Overlapping Gene E (inside Gene D: coords 568..843)
    overlap_e = genome[567:843]

    # Partition 2: Overlapping Gene K (overlaps Gene A and C: coords 51..221)
    overlap_k = genome[50:221]

    # Partition 3: Non-overlapping Gene F (Major coat protein: coords 1001..2284)
    non_overlap_f = genome[1000:2284]

    # Partition 4: Non-overlapping Gene G (Major spike protein: coords 2395..2922)
    non_overlap_g = genome[2394:2922]

    # Partition 5: Non-overlapping Gene H (Pilot protein: coords 2931..3917)
    non_overlap_h = genome[2930:3917]

    # Partition 6: Intergenic Spacer (coords 879..1000 between Gene J and F)
    intergenic_jf = genome[878:1000]

    partitions = [
        ("phix174_whole_genome", "category=whole_genome length=5386", genome),
        ("phix174_overlap_gene_e", "category=overlapping_cds coords=568-843 gene=E_inside_D", overlap_e),
        ("phix174_overlap_gene_k", "category=overlapping_cds coords=51-221 gene=K_overlaps_A_C", overlap_k),
        ("phix174_nonoverlap_gene_f", "category=non_overlapping_cds coords=1001-2284 gene=F_coat", non_overlap_f),
        ("phix174_nonoverlap_gene_g", "category=non_overlapping_cds coords=2395-2922 gene=G_spike", non_overlap_g),
        ("phix174_nonoverlap_gene_h", "category=non_overlapping_cds coords=2931-3917 gene=H_pilot", non_overlap_h),
        ("phix174_intergenic_jf", "category=intergenic coords=879-1000 between=J_and_F", intergenic_jf),
    ]

    part_file = DATA_DIR / "phix174_partitions.fasta"
    with part_file.open("w", encoding="utf-8") as f:
        for rec_id, desc, s in partitions:
            f.write(f">{rec_id} {desc}\n{s}\n")

    print(f"Generated {len(partitions)} preregistered functional partitions in: {part_file}")


if __name__ == "__main__":
    main()
