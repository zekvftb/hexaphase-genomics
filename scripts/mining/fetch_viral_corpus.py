"""Automated Ingestion of RefSeq Viral Genome Corpus for Overlapping Gene Mining.

Downloads and parses complete genomes across target compact viral families:
- Circoviridae (PCV2, BFDV)
- Parvoviridae (AAV2, B19)
- Anelloviridae (TTV)
- Polyomaviridae (SV40)
- Microviridae / Leviviridae (phiX174, MS2)
- Hepadnaviridae (HBV)
"""

from __future__ import annotations

import io
import json
from pathlib import Path
import sys
import time
import urllib.request

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from Bio import SeqIO

TARGET_VIRAL_ACCESSIONS = [
    {"accession": "NC_003977.2", "family": "Hepadnaviridae", "name": "Hepatitis B virus"},
    {"accession": "NC_001401.2", "family": "Parvoviridae", "name": "Adeno-associated virus 2 (AAV2)"},
    {"accession": "NC_005148.1", "family": "Circoviridae", "name": "Porcine circovirus 2"},
    {"accession": "NC_001456.1", "family": "Circoviridae", "name": "Beak and feather disease virus"},
    {"accession": "NC_000883.2", "family": "Parvoviridae", "name": "Human parvovirus B19"},
    {"accession": "NC_015783.1", "family": "Anelloviridae", "name": "Torque teno virus 1"},
    {"accession": "NC_001669.1", "family": "Polyomaviridae", "name": "Simian virus 40 (SV40)"},
    {"accession": "NC_001422.1", "family": "Microviridae", "name": "Bacteriophage phiX174"},
    {"accession": "NC_001417.2", "family": "Leviviridae", "name": "Bacteriophage MS2"},
]


def fetch_gbk_record(accession: str, cache_dir: Path) -> str:
    """Fetch GenBank record from NCBI E-utilities or return cached copy."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{accession}.gbk"

    if cache_file.is_file() and cache_file.stat().st_size > 500:
        return cache_file.read_text(encoding="utf-8")

    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=gbwithparts&retmode=text"
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "BioArch-ViralMiner/1.0 (Research Pipeline)"},
    )

    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=20) as resp:
                content = resp.read().decode("utf-8")
            if len(content) > 500:
                cache_file.write_text(content, encoding="utf-8")
                return content
        except Exception as e:
            if attempt == 2:
                raise RuntimeError(f"Failed to fetch {accession} after 3 attempts: {e}")
            time.sleep(1.0)

    raise RuntimeError(f"Failed to retrieve valid record for {accession}")


def ingest_viral_corpus(base_dir: Path | None = None) -> list[dict]:
    """Ingest and structure viral genomes for candidate smORF mining."""
    root = base_dir if base_dir is not None else Path(__file__).parent.parent.parent
    corpus_dir = root / "data" / "mining_corpus"
    corpus_dir.mkdir(parents=True, exist_ok=True)

    records_index = []

    for item in TARGET_VIRAL_ACCESSIONS:
        acc = item["accession"]
        try:
            gbk_text = fetch_gbk_record(acc, corpus_dir)
            rec = SeqIO.read(io.StringIO(gbk_text), "genbank")

            full_dna = str(rec.seq).upper()
            cds_features = []

            for f_idx, feat in enumerate(rec.features):
                if feat.type == "CDS":
                    gene_name = (
                        feat.qualifiers.get("gene", [""])[0] or
                        feat.qualifiers.get("product", [""])[0] or
                        f"CDS_{f_idx}"
                    )
                    product_desc = feat.qualifiers.get("product", [""])[0]
                    protein_id = feat.qualifiers.get("protein_id", [""])[0]
                    strand = feat.location.strand or 1
                    cds_dna = str(feat.extract(rec.seq)).upper()
                    translation = feat.qualifiers.get("translation", [""])[0]

                    cds_features.append({
                        "feature_index": f_idx,
                        "gene_name": gene_name,
                        "product": product_desc,
                        "protein_id": protein_id,
                        "strand": strand,
                        "start": int(feat.location.start),
                        "end": int(feat.location.end),
                        "dna_length": len(cds_dna),
                        "cds_dna": cds_dna,
                        "translation": translation,
                    })

            record_meta = {
                "accession": acc,
                "family": item["family"],
                "name": item["name"],
                "organism": rec.annotations.get("organism", item["name"]),
                "topology": rec.annotations.get("topology", "linear"),
                "genome_length_bp": len(full_dna),
                "num_annotated_cds": len(cds_features),
                "cds_features": cds_features,
            }

            json_file = corpus_dir / f"{acc}_metadata.json"
            json_file.write_text(json.dumps(record_meta, indent=2), encoding="utf-8")
            records_index.append(record_meta)

            print(f"✅ Ingested {acc}: {item['name']} ({len(full_dna)} bp, {len(cds_features)} CDS)")

        except Exception as e:
            print(f"⚠️ Warning: Could not ingest {acc}: {e}")

    # Write master corpus index
    index_file = corpus_dir / "viral_corpus_manifest.json"
    index_file.write_text(json.dumps(records_index, indent=2), encoding="utf-8")
    print(f"\n[CORPUS READY] Ingested {len(records_index)} viral genomes in data/mining_corpus/\n")
    return records_index


if __name__ == "__main__":
    ingest_viral_corpus()
