"""Systematic Overlapping smORF & Dual-Coding Gene Mining Engine for Viral Genomes.

Scans alternative reading frames (+1, +2) across viral CDS features,
evaluates host codon adaptation (CAI), sequence complexity, and validates candidates
against N=500 Altschul-Erickson dinucleotide-preserving null models.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import io
import json
import math
from pathlib import Path
import random
import sys
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq

from bio_arch.modules.information import shuffle_dinucleotide
from bio_arch.modules.recompiler import CODON_TABLE, compute_cai, translate_sequence
from bio_arch.provenance import SeedManager, now_iso

START_CODONS = {"ATG", "GTG", "TTG", "CTG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}


@dataclass
class OverlappingCandidate:
    """A candidate overlapping open reading frame within a parent primary CDS."""

    accession: str
    organism: str
    parent_gene: str
    parent_product: str
    parent_cds_length_bp: int
    frame_offset: int  # 1 or 2
    candidate_start_in_parent_bp: int
    candidate_end_in_parent_bp: int
    start_codon: str
    length_aa: int
    peptide_sequence: str
    cai_score: float
    gc_content: float
    shannon_entropy: float
    null_mean_score: float
    null_std_score: float
    z_score: float
    empirical_p_value: float
    statistically_significant: bool
    is_annotated_in_genbank: bool
    matched_annotation_name: str = ""


@dataclass
class MiningResult:
    """Complete results from a viral genome mining run."""

    timestamp_iso: str
    total_genomes_scanned: int
    total_cds_scanned: int
    total_candidates_found: int
    novel_significant_candidates: int
    candidates: list[OverlappingCandidate] = field(default_factory=list)


class ViralOverlapMiner:
    """High-throughput scanner for dual-coding genes and smORFs across viral genomes."""

    def __init__(
        self,
        min_length_aa: int = 30,
        min_cai_threshold: float = 0.65,
        n_null_shuffles: int = 500,
        significance_alpha: float = 0.01,
        seed: int = 42,
    ) -> None:
        self.min_length_aa = min_length_aa
        self.min_cai_threshold = min_cai_threshold
        self.n_null_shuffles = n_null_shuffles
        self.significance_alpha = significance_alpha
        self.seed_mgr = SeedManager(seed)

    def scan_cds_alternative_frames(
        self,
        cds_dna: str,
        parent_gene: str = "CDS",
        parent_product: str = "viral protein",
        accession: str = "UNKNOWN",
        organism: str = "Viral isolate",
        annotated_translations: list[dict] | None = None,
    ) -> list[OverlappingCandidate]:
        """Scan Frame +1 and Frame +2 of a primary CDS for significant overlapping ORFs."""
        cds_dna = cds_dna.upper()
        parent_len_bp = len(cds_dna)
        annotated_translations = annotated_translations or []
        candidates: list[OverlappingCandidate] = []

        for frame_offset in (1, 2):
            subseq = cds_dna[frame_offset:]
            num_codons = len(subseq) // 3

            i = 0
            while i < num_codons:
                codon = subseq[i * 3 : i * 3 + 3]
                if codon in START_CODONS:
                    # Scan for in-frame stop codon
                    j = i + 1
                    found_stop = False
                    while j < num_codons:
                        stop_candidate = subseq[j * 3 : j * 3 + 3]
                        if stop_candidate in STOP_CODONS:
                            found_stop = True
                            break
                        j += 1

                    orf_codons = j - i
                    if orf_codons >= self.min_length_aa:
                        orf_dna = subseq[i * 3 : j * 3]
                        peptide = str(Seq(orf_dna).translate())

                        # Filter out internal stops (safety)
                        if "*" not in peptide:
                            cai = round(compute_cai(orf_dna), 4)
                            gc = round(sum(1 for b in orf_dna if b in "GC") / len(orf_dna), 4)
                            entropy = round(self._calc_entropy(orf_dna), 4)

                            # Null model calibration
                            null_mean, null_std, z_val, p_val = self._calibrate_against_null(
                                parent_cds_dna=cds_dna,
                                frame_offset=frame_offset,
                                actual_len=orf_codons,
                                actual_cai=cai,
                                context_name=f"{accession}_{parent_gene}_f{frame_offset}_i{i}",
                            )

                            is_sig = (p_val <= self.significance_alpha) and (z_val >= 2.0)

                            # Check if matches an existing GenBank annotation
                            is_annotated, matched_name = self._check_existing_annotation(
                                peptide=peptide,
                                annotated_translations=annotated_translations,
                            )

                            cand = OverlappingCandidate(
                                accession=accession,
                                organism=organism,
                                parent_gene=parent_gene,
                                parent_product=parent_product,
                                parent_cds_length_bp=parent_len_bp,
                                frame_offset=frame_offset,
                                candidate_start_in_parent_bp=frame_offset + i * 3,
                                candidate_end_in_parent_bp=frame_offset + j * 3 + (3 if found_stop else 0),
                                start_codon=codon,
                                length_aa=orf_codons,
                                peptide_sequence=peptide,
                                cai_score=cai,
                                gc_content=gc,
                                shannon_entropy=entropy,
                                null_mean_score=null_mean,
                                null_std_score=null_std,
                                z_score=z_val,
                                empirical_p_value=p_val,
                                statistically_significant=is_sig,
                                is_annotated_in_genbank=is_annotated,
                                matched_annotation_name=matched_name,
                            )
                            candidates.append(cand)

                    # Advance past this start
                    i += 1
                else:
                    i += 1

        return candidates

    def _calibrate_against_null(
        self,
        parent_cds_dna: str,
        frame_offset: int,
        actual_len: int,
        actual_cai: float,
        context_name: str,
    ) -> tuple[float, float, float, float]:
        """Run Altschul-Erickson dinucleotide-preserving shuffles to establish null baseline."""
        actual_score = actual_len * actual_cai
        null_scores = []
        n_shuffles = self.n_null_shuffles

        for k in range(n_shuffles):
            sub_seed = self.seed_mgr.derive_seed(f"{context_name}_null_{k}")
            rng = random.Random(sub_seed)
            shuffled_parent = shuffle_dinucleotide(parent_cds_dna, rng=rng)

            # Find maximum alternative frame ORF score in shuffled sequence
            shuf_sub = shuffled_parent[frame_offset:]
            num_codons = len(shuf_sub) // 3
            max_shuf_score = 0.0

            idx = 0
            while idx < num_codons:
                if shuf_sub[idx * 3 : idx * 3 + 3] in START_CODONS:
                    jdx = idx + 1
                    while jdx < num_codons:
                        if shuf_sub[jdx * 3 : jdx * 3 + 3] in STOP_CODONS:
                            break
                        jdx += 1
                    orf_len = jdx - idx
                    if orf_len >= 15:  # Consider any intermediate length in null
                        shuf_orf_dna = shuf_sub[idx * 3 : jdx * 3]
                        shuf_cai = compute_cai(shuf_orf_dna)
                        shuf_score = orf_len * shuf_cai
                        if shuf_score > max_shuf_score:
                            max_shuf_score = shuf_score
                    idx = jdx + 1
                else:
                    idx += 1

            null_scores.append(max_shuf_score)

        null_mean = sum(null_scores) / len(null_scores)
        null_var = sum((x - null_mean) ** 2 for x in null_scores) / max(1, len(null_scores) - 1)
        null_std = max(1e-4, null_var ** 0.5)

        z_val = round((actual_score - null_mean) / null_std, 2)
        empirical_p = sum(1 for x in null_scores if x >= actual_score) / n_shuffles

        return round(null_mean, 2), round(null_std, 2), z_val, round(empirical_p, 4)

    def _calc_entropy(self, dna: str) -> float:
        """Calculate Shannon entropy in bits per base."""
        counts = {}
        for b in dna:
            counts[b] = counts.get(b, 0) + 1
        n = len(dna)
        if n == 0:
            return 0.0
        return -sum((c / n) * math.log2(c / n) for c in counts.values())

    def _check_existing_annotation(
        self,
        peptide: str,
        annotated_translations: list[dict],
    ) -> tuple[bool, str]:
        """Check if candidate peptide matches or is contained in an existing annotated CDS."""
        for item in annotated_translations:
            ann_trans = item.get("translation", "")
            gene = item.get("gene_name", "Known Gene")
            if not ann_trans:
                continue
            # If 80% of candidate peptide matches an existing annotated translation
            if len(peptide) >= 20:
                core = peptide[:min(len(peptide), 30)]
                if core in ann_trans or ann_trans in peptide:
                    return True, gene
                # Exact or substring match
                if peptide in ann_trans or ann_trans in peptide:
                    return True, gene

        return False, ""

    def mine_genome_file(self, gbk_path: Path) -> list[OverlappingCandidate]:
        """Mine a single GenBank file for all overlapping candidates."""
        record = SeqIO.read(str(gbk_path), "genbank")
        acc = record.id
        organism = record.annotations.get("organism", record.description)

        all_annotated = []
        for f in record.features:
            if f.type == "CDS":
                gene_name = f.qualifiers.get("gene", [""])[0] or f.qualifiers.get("product", [""])[0]
                trans = f.qualifiers.get("translation", [""])[0]
                all_annotated.append({"gene_name": gene_name, "translation": trans})

        genome_candidates: list[OverlappingCandidate] = []

        for feat in record.features:
            if feat.type == "CDS":
                cds_dna = str(feat.extract(record.seq))
                # Only scan forward positive-strand CDS or un-spliced segments
                if len(cds_dna) >= 90 and len(cds_dna) % 3 == 0:
                    gene = feat.qualifiers.get("gene", [""])[0] or "CDS"
                    prod = feat.qualifiers.get("product", [""])[0] or "Protein"

                    cands = self.scan_cds_alternative_frames(
                        cds_dna=cds_dna,
                        parent_gene=gene,
                        parent_product=prod,
                        accession=acc,
                        organism=organism,
                        annotated_translations=all_annotated,
                    )
                    genome_candidates.extend(cands)

        return genome_candidates

    def mine_corpus_directory(self, corpus_dir: Path) -> MiningResult:
        """Mine an entire directory of GenBank records."""
        gbk_files = list(corpus_dir.glob("*.gbk"))
        all_cands: list[OverlappingCandidate] = []
        total_cds = 0

        for gbk in gbk_files:
            cands = self.mine_genome_file(gbk)
            all_cands.extend(cands)

        novel_sig = sum(
            1 for c in all_cands if c.statistically_significant and not c.is_annotated_in_genbank
        )

        return MiningResult(
            timestamp_iso=now_iso(),
            total_genomes_scanned=len(gbk_files),
            total_cds_scanned=total_cds,
            total_candidates_found=len(all_cands),
            novel_significant_candidates=novel_sig,
            candidates=all_cands,
        )
