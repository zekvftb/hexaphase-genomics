"""P-Site Extraction & Triplet-Periodicity Phasing Engine for Ribo-seq Data.

Extracts ribosome P-sites from protected footprint alignments using length-dependent
offset calibration (28nt -> +12, 29nt -> +12, 30nt -> +13), computes 3-frame phasing distributions,
calculates the Triplet Periodicity Index (TPI), performs Chi-square goodness-of-fit testing,
and classifies active in-vivo translation evidence.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import io
import math
from pathlib import Path
import random
import sys
from typing import Any

# Length-dependent P-site offset calibration (nt from 5' end of read)
P_SITE_OFFSETS: dict[int, int] = {
    26: 12,
    27: 12,
    28: 12,
    29: 12,
    30: 13,
    31: 13,
    32: 13,
}
DEFAULT_P_SITE_OFFSET = 12


@dataclass
class RiboPhasingResult:
    """Quantitative results of reading-frame triplet periodicity evaluation."""

    candidate_id: str
    accession: str
    parent_gene: str
    target_frame: int
    start_nt: int
    end_nt: int
    total_psites: int
    frame_counts: dict[int, int]
    frame_fractions: dict[int, float]
    triplet_periodicity_index: float
    chi_square_stat: float
    chi_square_p_value: float
    translation_evidence_score: float
    classification: str
    bioproject_source: str = ""
    sra_run_source: str = ""


def get_psite_offset(read_length_nt: int) -> int:
    """Get length-calibrated offset from 5' end to ribosome P-site."""
    return P_SITE_OFFSETS.get(read_length_nt, DEFAULT_P_SITE_OFFSET)


def calculate_chi_square_uniformity(counts: list[int]) -> tuple[float, float]:
    """Calculate Chi-square statistic and exact p-value against uniform distribution (df=2).

    For df = 2, the Chi-square survival function (1 - CDF) is exact: P(X >= x) = exp(-x / 2).
    """
    total = sum(counts)
    if total < 3:
        return 0.0, 1.0

    expected = total / 3.0
    chi2 = sum(((obs - expected) ** 2) / expected for obs in counts)

    # Exact survival function for df = 2 Chi-square distribution
    p_val = math.exp(-chi2 / 2.0)
    return round(chi2, 3), max(1e-12, min(1.0, round(p_val, 6)))


def analyze_ribo_phasing(
    read_alignments: list[dict],
    candidate_id: str,
    accession: str,
    parent_gene: str,
    start_nt: int,
    end_nt: int,
    target_frame: int,
    min_coverage_threshold: int = 15,
    bioproject_id: str = "",
    sra_run_id: str = "",
) -> RiboPhasingResult:
    """Analyze reading frame phasing and triplet periodicity for a candidate smORF window."""
    frame_counts = {0: 0, 1: 0, 2: 0}
    total_psites = 0

    for read in read_alignments:
        read_start = read.get("start_nt", 0)
        read_len = read.get("length_nt", 28)
        offset = get_psite_offset(read_len)
        psite_pos = read_start + offset

        # Check if P-site falls within the candidate genomic window
        if start_nt <= psite_pos <= end_nt:
            # Phase relative to start of primary CDS / window
            rel_pos = psite_pos - start_nt
            frame = rel_pos % 3
            frame_counts[frame] += 1
            total_psites += 1

    obs_list = [frame_counts[0], frame_counts[1], frame_counts[2]]
    chi2_stat, p_val = calculate_chi_square_uniformity(obs_list)

    if total_psites > 0:
        fractions = {f: round(cnt / total_psites, 4) for f, cnt in frame_counts.items()}
        tpi = fractions.get(target_frame, 0.0)
    else:
        fractions = {0: 0.0, 1: 0.0, 2: 0.0}
        tpi = 0.0

    # Translation Evidence Score: TPI * -log10(p_value)
    neg_log10_p = -math.log10(max(1e-12, p_val))
    evidence_score = round(tpi * neg_log10_p, 3)

    # Classification
    if total_psites < min_coverage_threshold:
        classification = "Ambiguous / Low Coverage"
    elif tpi >= 0.60 and p_val < 0.01:
        classification = "Confirmed Active Translation"
    elif tpi >= 0.50 and p_val < 0.05:
        classification = "Moderate Translation Evidence"
    else:
        classification = "Non-Translating / Background Noise"

    return RiboPhasingResult(
        candidate_id=candidate_id,
        accession=accession,
        parent_gene=parent_gene,
        target_frame=target_frame,
        start_nt=start_nt,
        end_nt=end_nt,
        total_psites=total_psites,
        frame_counts=frame_counts,
        frame_fractions=fractions,
        triplet_periodicity_index=tpi,
        chi_square_stat=chi2_stat,
        chi_square_p_value=p_val,
        translation_evidence_score=evidence_score,
        classification=classification,
        bioproject_source=bioproject_id,
        sra_run_source=sra_run_id,
    )


def generate_synthetic_ribo_reads(
    start_nt: int,
    end_nt: int,
    target_frame: int,
    signal_to_noise: float = 0.75,
    total_reads: int = 150,
    seed: int = 42,
) -> list[dict]:
    """Generate deterministic synthetic Ribo-seq reads with specified frame phasing signal."""
    rng = random.Random(seed)
    reads = []
    orf_len_nt = end_nt - start_nt

    for _ in range(total_reads):
        # Choose length from realistic Ribo-seq distribution
        read_len = rng.choice([28, 28, 29, 29, 30, 31])
        offset = get_psite_offset(read_len)

        # Decide whether this read aligns with the target frame signal or background noise
        if rng.random() < signal_to_noise:
            frame = target_frame
        else:
            frame = rng.choice([0, 1, 2])

        # Pick a codon position in the ORF
        codon_idx = rng.randint(0, max(0, (orf_len_nt // 3) - 1))
        psite_pos = start_nt + (codon_idx * 3) + frame
        read_start = psite_pos - offset

        reads.append({
            "start_nt": read_start,
            "length_nt": read_len,
        })

    return reads
