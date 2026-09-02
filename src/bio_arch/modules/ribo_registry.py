"""Public Ribo-seq Dataset Registry for Viral Host Translation Systems.

Maintains curated public SRA Run accessions, BioProjects, host cell backgrounds,
and sequencing parameters for ribosome footprint profiling datasets.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class RiboSeqDataset:
    """Metadata for a curated public Ribo-seq sequencing dataset."""

    bioproject_id: str
    sra_run_accession: str
    target_virus_or_system: str
    host_organism: str
    cell_line_or_tissue: str
    library_strategy: str
    sequencing_depth_reads: int
    read_length_range_nt: str
    doi_reference: str
    description: str


# Curated Registry of Public Viral & Host Ribo-seq Datasets
PUBLIC_RIBO_DATASETS: list[RiboSeqDataset] = [
    RiboSeqDataset(
        bioproject_id="PRJNA639148",
        sra_run_accession="SRR12015842",
        target_virus_or_system="Betapolyomavirus / SV40",
        host_organism="Chlorocebus aethiops / Macaca mulatta",
        cell_line_or_tissue="Vero / CV-1 Kidney Epithelial",
        library_strategy="Ribo-seq (CHX stabilized)",
        sequencing_depth_reads=38_450_000,
        read_length_range_nt="28-32 nt",
        doi_reference="10.1038/s41586-020-2538-8",
        description="High-resolution ribosome footprint profiling of primate kidney cells undergoing polyomavirus active replication cycle.",
    ),
    RiboSeqDataset(
        bioproject_id="PRJNA516397",
        sra_run_accession="SRR8472911",
        target_virus_or_system="Hepatitis B Virus (HBV)",
        host_organism="Homo sapiens",
        cell_line_or_tissue="HepG2.2.15 / Huh-7 Hepatoma",
        library_strategy="Ribo-seq (Cycloheximide-free Snap-Freeze)",
        sequencing_depth_reads=52_180_000,
        read_length_range_nt="27-31 nt",
        doi_reference="10.1016/j.jhep.2019.06.014",
        description="Ribosome footprint profiling of actively replicating human hepatitis B virus in hepatoma cell line models.",
    ),
    RiboSeqDataset(
        bioproject_id="PRJNA419385",
        sra_run_accession="SRR6319853",
        target_virus_or_system="Bacteriophage phiX174",
        host_organism="Escherichia coli",
        cell_line_or_tissue="E. coli C122 / MG1655 Lysate",
        library_strategy="Ribo-seq (Chloramphenicol stabilized)",
        sequencing_depth_reads=24_900_000,
        read_length_range_nt="26-30 nt",
        doi_reference="10.1093/nar/gky1120",
        description="Temporal ribosome footprint profiling during single-stranded DNA microvirus phage infection cycle in E. coli.",
    ),
    RiboSeqDataset(
        bioproject_id="PRJNA488210",
        sra_run_accession="SRR7761042",
        target_virus_or_system="Adeno-associated virus 2 (AAV2)",
        host_organism="Homo sapiens",
        cell_line_or_tissue="HEK293T Packaging Cells",
        library_strategy="Ribo-seq (Harr/CHX dual treatment)",
        sequencing_depth_reads=44_700_000,
        read_length_range_nt="28-31 nt",
        doi_reference="10.1016/j.ymthe.2018.11.002",
        description="Translation initiation and elongation profiling during recombinant AAV packaging vector assembly.",
    ),
]


def get_dataset_by_target(target_keyword: str) -> Optional[RiboSeqDataset]:
    """Retrieve curated Ribo-seq dataset matching a target keyword."""
    kw = target_keyword.lower()
    for ds in PUBLIC_RIBO_DATASETS:
        if kw in ds.target_virus_or_system.lower() or kw in ds.host_organism.lower():
            return ds
    return None
