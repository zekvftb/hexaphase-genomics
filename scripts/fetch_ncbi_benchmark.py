"""Ground-Truth Ingestion: Fetch and parse canonical HBV overlapping genes from NCBI.

Extracts the experimentally validated dual-coding overlap between Polymerase (Frame 0)
and Small Surface Antigen (Frame +1) from Hepatitis B Virus (NCBI RefSeq: NC_003977.2).
"""

from __future__ import annotations

import io
import json
from pathlib import Path
import sys
import urllib.request

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from Bio import SeqIO
from Bio.Seq import Seq

from bio_arch.modules.recompiler import compute_cai

NCBI_ACCESSION = "NC_003977.2"
EUTILS_URL = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={NCBI_ACCESSION}&rettype=gbwithparts&retmode=text"


def fetch_ncbi_record(cache_dir: Path) -> str:
    """Fetch GenBank record from NCBI E-utilities or read from local cache."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{NCBI_ACCESSION}.gbk"

    if cache_file.is_file() and cache_file.stat().st_size > 1000:
        return cache_file.read_text(encoding="utf-8")

    req = urllib.request.Request(
        EUTILS_URL,
        headers={"User-Agent": "BioArch-NCBI-GroundTruth-Validator/1.0"},
    )
    with urllib.request.urlopen(req, timeout=20) as response:
        content = response.read().decode("utf-8")

    cache_file.write_text(content, encoding="utf-8")
    return content


def extract_hbv_overlap_ground_truth(base_dir: Path | None = None) -> dict:
    """Extract and validate the exact overlapping coding sequence between P and S genes."""
    root = base_dir if base_dir is not None else Path(__file__).parent.parent
    benchmarks_dir = root / "data" / "benchmarks"
    benchmarks_dir.mkdir(parents=True, exist_ok=True)

    gbk_text = fetch_ncbi_record(benchmarks_dir)
    record = SeqIO.read(io.StringIO(gbk_text), "genbank")

    # Locate Polymerase (P) and Small Surface Antigen (S) CDS features
    p_cds = [f for f in record.features if f.type == "CDS" and f.qualifiers.get("gene") == ["P"]][0]
    s_cds = [f for f in record.features if f.type == "CDS" and f.qualifiers.get("product") == ["small envelope protein"]][0]

    p_full_dna = str(p_cds.extract(record.seq))
    s_full_dna = str(s_cds.extract(record.seq))

    # The S gene begins at nucleotide offset 1030 in the Polymerase CDS (+1 phase offset)
    s_start_in_p = p_full_dna.find(s_full_dna[:30])
    assert s_start_in_p % 3 == 1, f"Expected offset +1 (+1 phase), got {s_start_in_p % 3}"

    # Extract overlapping segment: Frame 0 starts at codon boundary (s_start_in_p - 1 = 1029)
    overlap_start = s_start_in_p - 1
    overlap_len_bp = len(s_full_dna) + 1  # 682 bp total
    overlap_dna = p_full_dna[overlap_start : overlap_start + overlap_len_bp]

    # Translate both reading frames
    # Frame 0: Polymerase segment (226 aa)
    p_trans = str(Seq(overlap_dna).translate())[:226]
    # Frame +1: Small Surface Antigen segment (226 aa + stop)
    s_trans = str(Seq(overlap_dna[1:]).translate())[:226]

    # Calculate wild-type human CAI
    cai_f0 = round(compute_cai(overlap_dna), 4)
    cai_f1 = round(compute_cai(overlap_dna[1:]), 4)

    payload = {
        "accession": NCBI_ACCESSION,
        "organism": record.annotations.get("organism", "Hepatitis B virus"),
        "total_genome_length_bp": len(record.seq),
        "overlap_dna_length_bp": len(overlap_dna),
        "natural_dna_sequence": overlap_dna,
        "natural_protein_0_name": "HBV Polymerase (overlapping domain)",
        "natural_protein_0": p_trans,
        "natural_polymerase_peptide": p_trans,
        "natural_protein_1_name": "HBV Small Surface Antigen (HBsAg)",
        "natural_protein_1": s_trans,
        "natural_surface_peptide": s_trans,
        "natural_f0_cai": cai_f0,
        "natural_cai_f0": cai_f0,
        "natural_f1_cai": cai_f1,
        "natural_cai_f1": cai_f1,
        "overlap_start_in_p_cds": overlap_start,
        "overlap_frame_offset": 1,
    }

    out_file = benchmarks_dir / "hbv_overlap_ground_truth.json"
    out_file.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


if __name__ == "__main__":
    res = extract_hbv_overlap_ground_truth()
    print(f"[NCBI INGESTION] Fetched {res['accession']} ({res['organism']})")
    print(f"  Overlap Length: {res['overlap_dna_length_bp']} bp ({len(res['natural_protein_0'])} aa per frame)")
    print(f"  Frame 0 (Polymerase): {res['natural_protein_0'][:30]}...")
    print(f"  Frame +1 (HBsAg):     {res['natural_protein_1'][:30]}...")
    print(f"  Natural CAI: F0={res['natural_f0_cai']}, F1={res['natural_f1_cai']}")
    print(f"  Saved ground-truth to: data/benchmarks/hbv_overlap_ground_truth.json")
