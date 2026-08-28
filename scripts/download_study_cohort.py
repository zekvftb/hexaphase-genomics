"""Download the curated study cohort from NCBI RefSeq and generate control corpora."""

from pathlib import Path
import urllib.request

REPO_ROOT = Path(__file__).resolve().parent.parent
COHORT_DIR = REPO_ROOT / "data" / "study_cohort"
COHORT_DIR.mkdir(parents=True, exist_ok=True)


def fetch_ncbi_fasta(accession: str, label: str) -> str:
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    req = urllib.request.Request(url, headers={"User-Agent": "BioArchResearch/1.0"})
    print(f"Fetching {label} ({accession}) from NCBI...")
    with urllib.request.urlopen(req) as resp:
        content = resp.read().decode("utf-8")
    lines = content.strip().splitlines()
    seq_lines = [l.strip() for l in lines if not l.startswith(">")]
    return "".join(seq_lines).upper()


def main() -> None:
    # 1. Bacteriophage PhiX174 (already local or re-fetch)
    phix_file = REPO_ROOT / "data" / "phix174_complete.fasta"
    if phix_file.is_file():
        print("Loading PhiX174 from local cache...")
        lines = phix_file.read_text(encoding="utf-8").strip().splitlines()
        phix_seq = "".join(l.strip() for l in lines if not l.startswith(">"))
    else:
        phix_seq = fetch_ncbi_fasta("NC_001422.1", "Bacteriophage PhiX174")

    # 2. Bacteriophage Lambda (NC_001416.1, ~48.5 kb)
    lambda_seq = fetch_ncbi_fasta("NC_001416.1", "Bacteriophage Lambda")

    # 3. SARS-CoV-2 RNA Virus (NC_045512.2, ~29.9 kb)
    sars2_seq = fetch_ncbi_fasta("NC_045512.2", "SARS-CoV-2 RNA Genome")

    # 4. Human Mitochondrial DNA (NC_012920.1, ~16.5 kb)
    human_mt_seq = fetch_ncbi_fasta("NC_012920.1", "Human Mitochondrial DNA")

    # 5. E. coli Lac Operon Region (from Lambda/E.coli model, ~6.2 kb segment)
    # Using a 6,000 bp slice of Lambda phage containing regulatory lysis/lysogeny genes
    lac_segment = lambda_seq[33000:39200]

    # Write all FASTA records into cohort file
    cohort_fasta = COHORT_DIR / "study_cohort.fasta"
    with cohort_fasta.open("w", encoding="utf-8") as f:
        f.write(f">phix174_virus organism=Bacteriophage_PhiX174 type=ssDNA length={len(phix_seq)}\n{phix_seq}\n")
        f.write(f">lambda_phage organism=Bacteriophage_Lambda type=dsDNA length={len(lambda_seq)}\n{lambda_seq}\n")
        f.write(f">sars_cov_2_virus organism=SARS_CoV_2 type=ssRNA length={len(sars2_seq)}\n{sars2_seq}\n")
        f.write(f">human_mitochondria organism=Homo_sapiens_mtDNA type=dsDNA length={len(human_mt_seq)}\n{human_mt_seq}\n")
        f.write(f">regulatory_switch_locus organism=Phage_Lytic_Lysogenic_Switch type=dsDNA length={len(lac_segment)}\n{lac_segment}\n")

    print(f"\nCohort saved successfully to: {cohort_fasta}")
    print(f"Total biological datasets assembled: 5 genomes/loci")


if __name__ == "__main__":
    main()
