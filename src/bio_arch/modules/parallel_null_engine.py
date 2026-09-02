"""Multi-Core Parallelized Null-Model Screening Engine.

Executes high-throughput scanning of alternative reading frames (+1, +2) across viral CDS features,
distributes N=500 Eulerian-walk dinucleotide shuffles across CPU cores, and enforces a hard
barrier filter (z <= 3.0, p >= 0.001 dropped immediately) to preserve downstream compute.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass, field
import io
import math
import os
from pathlib import Path
import random
import sys
import time
from typing import Any

from Bio.Seq import Seq

_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

from bio_arch.modules.information import shuffle_dinucleotide
from bio_arch.modules.non_canonical_miner import START_CODON_EFFICIENCY, STOP_CODONS
from bio_arch.modules.recompiler import compute_cai
from bio_arch.provenance import SeedManager


@dataclass
class FastCandidate:
    """Lightweight representation of an alternative-frame candidate smORF."""

    accession: str
    organism: str
    parent_gene: str
    parent_product: str
    frame_offset: int
    start_in_parent_bp: int
    end_in_parent_bp: int
    start_codon: str
    length_aa: int
    peptide_sequence: str
    host_cai: float
    parent_cds_seq: str
    null_mean: float = 0.0
    null_std: float = 0.0
    z_score: float = 0.0
    empirical_p_value: float = 1.0
    passes_hard_barrier: bool = False


def _evaluate_null_worker(cand: FastCandidate, n_shuffles: int, seed: int) -> FastCandidate:
    """Worker function to run N=500 Eulerian dinucleotide null shuffles on a single candidate."""
    actual_score = cand.length_aa * cand.host_cai
    parent_seq = cand.parent_cds_seq
    frame_offset = cand.frame_offset
    seed_mgr = SeedManager(seed)

    null_scores = []
    for k in range(n_shuffles):
        sub_seed = seed_mgr.derive_seed(f"{cand.accession}_{cand.parent_gene}_f{frame_offset}_null_{k}")
        rng = random.Random(sub_seed)
        shuf_parent = shuffle_dinucleotide(parent_seq, rng=rng)

        shuf_sub = shuf_parent[frame_offset:]
        num_codons = len(shuf_sub) // 3
        max_s = 0.0

        idx = 0
        while idx < num_codons:
            codon = shuf_sub[idx * 3 : idx * 3 + 3]
            if codon in START_CODON_EFFICIENCY:
                jdx = idx + 1
                while jdx < num_codons:
                    if shuf_sub[jdx * 3 : jdx * 3 + 3] in STOP_CODONS:
                        break
                    jdx += 1
                orf_len = jdx - idx
                if orf_len >= 15:
                    orf_dna = shuf_sub[idx * 3 : jdx * 3]
                    cai = compute_cai(orf_dna)
                    s = orf_len * cai
                    if s > max_s:
                        max_s = s
                idx = jdx + 1
            else:
                idx += 1
        null_scores.append(max_s)

    mean = sum(null_scores) / len(null_scores)
    var = sum((x - mean) ** 2 for x in null_scores) / max(1, len(null_scores) - 1)
    std = max(1e-4, var ** 0.5)

    z = round((actual_score - mean) / std, 2)
    p = round(sum(1 for x in null_scores if x >= actual_score) / n_shuffles, 4)

    cand.null_mean = round(mean, 2)
    cand.null_std = round(std, 2)
    cand.z_score = z
    cand.empirical_p_value = p
    cand.passes_hard_barrier = (z > 3.0) and (p <= 0.001)

    return cand


class ParallelNullEngine:
    """High-throughput multi-core scanner and null-model calibration engine."""

    def __init__(
        self,
        max_workers: int | None = None,
        n_shuffles: int = 500,
        min_length_aa: int = 30,
        max_length_aa: int = 250,
        min_cai_threshold: float = 0.70,
        z_barrier_threshold: float = 3.0,
        p_barrier_threshold: float = 0.001,
        seed: int = 42,
    ) -> None:
        self.max_workers = max_workers or min(8, os.cpu_count() or 4)
        self.n_shuffles = n_shuffles
        self.min_length_aa = min_length_aa
        self.max_length_aa = max_length_aa
        self.min_cai_threshold = min_cai_threshold
        self.z_barrier_threshold = z_barrier_threshold
        self.p_barrier_threshold = p_barrier_threshold
        self.seed = seed

    def scan_cds_fast(self, cds_row: dict[str, Any], genome_meta: dict[str, Any]) -> list[FastCandidate]:
        """Stage 1: Fast sequence scan for candidate smORFs in Frame +1 and Frame +2."""
        cds_dna = cds_row.get("sequence", "").upper()
        if len(cds_dna) < 90 or len(cds_dna) % 3 != 0:
            return []

        acc = genome_meta.get("accession", cds_row.get("accession", "UNKNOWN"))
        org = genome_meta.get("organism", "Viral Genome")
        gene = cds_row.get("gene", "CDS")
        prod = cds_row.get("product", "viral protein")

        initial_cands = []

        for frame_offset in (1, 2):
            subseq = cds_dna[frame_offset:]
            num_codons = len(subseq) // 3

            i = 0
            while i < num_codons:
                codon = subseq[i * 3 : i * 3 + 3]
                if codon in START_CODON_EFFICIENCY:
                    j = i + 1
                    found_stop = False
                    while j < num_codons:
                        if subseq[j * 3 : j * 3 + 3] in STOP_CODONS:
                            found_stop = True
                            break
                        j += 1

                    orf_len = j - i
                    if self.min_length_aa <= orf_len <= self.max_length_aa:
                        orf_dna = subseq[i * 3 : j * 3]
                        peptide = str(Seq(orf_dna).translate())

                        if "*" not in peptide:
                            cai = round(compute_cai(orf_dna), 4)
                            if cai >= self.min_cai_threshold:
                                abs_start = frame_offset + i * 3
                                abs_end = frame_offset + j * 3 + (3 if found_stop else 0)

                                cand = FastCandidate(
                                    accession=acc,
                                    organism=org,
                                    parent_gene=gene,
                                    parent_product=prod,
                                    frame_offset=frame_offset,
                                    start_in_parent_bp=abs_start,
                                    end_in_parent_bp=abs_end,
                                    start_codon=codon,
                                    length_aa=orf_len,
                                    peptide_sequence=peptide,
                                    host_cai=cai,
                                    parent_cds_seq=cds_dna,
                                )
                                initial_cands.append(cand)
                    i += 1
                else:
                    i += 1

        return initial_cands

    def run_screening_pipeline(
        self,
        genomes_with_cds: list[tuple[dict[str, Any], list[dict[str, Any]]]],
    ) -> tuple[list[FastCandidate], dict[str, Any]]:
        """Stage 1 + Stage 2: Parallel screening with hard barrier filtering."""
        start_time = time.time()

        # Stage 1: Collect all candidate smORFs across all genomes
        all_stage1_candidates: list[FastCandidate] = []
        total_cds_count = 0

        for genome_meta, cds_list in genomes_with_cds:
            total_cds_count += len(cds_list)
            for cds_row in cds_list:
                cands = self.scan_cds_fast(cds_row, genome_meta)
                all_stage1_candidates.extend(cands)

        # Stage 2: Parallelized N=500 null shuffle evaluation across CPU cores
        processed_candidates: list[FastCandidate] = []
        if all_stage1_candidates:
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                futures = [
                    executor.submit(_evaluate_null_worker, cand, self.n_shuffles, self.seed + idx)
                    for idx, cand in enumerate(all_stage1_candidates)
                ]
                for f in futures:
                    processed_candidates.append(f.result())

        # Hard Barrier Filter: Drop z <= 3.0 immediately
        surviving_candidates = [c for c in processed_candidates if c.passes_hard_barrier]

        elapsed_sec = max(0.001, time.time() - start_time)
        n_genomes = len(genomes_with_cds)
        n_initial = len(all_stage1_candidates)
        n_surviving = len(surviving_candidates)

        telemetry = {
            "total_genomes_screened": n_genomes,
            "total_cds_screened": total_cds_count,
            "stage1_candidates_found": n_initial,
            "surviving_candidates_z_gt_3": n_surviving,
            "hard_barrier_rejection_rate_pct": round(((n_initial - n_surviving) / max(1, n_initial)) * 100.0, 2),
            "execution_time_seconds": round(elapsed_sec, 3),
            "throughput_genomes_per_sec": round(n_genomes / elapsed_sec, 2),
            "throughput_candidates_per_sec": round(n_initial / elapsed_sec, 2),
            "parallel_workers": self.max_workers,
            "null_shuffles_per_candidate": self.n_shuffles,
        }

        return surviving_candidates, telemetry
