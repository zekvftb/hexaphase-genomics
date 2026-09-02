"""Pre-Registered Hypothesis Test B: Genomic Grammar & Zipf's Law Power-Law Scaling.

Evaluates whether k-mer rank-frequency distributions in natural genomes follow a true
power law distinct from dinucleotide-preserved randomized sequences.
"""

from __future__ import annotations

from collections import Counter
import json
import math
from pathlib import Path
import random
import sys

_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from Bio import SeqIO

from bio_arch.modules.information import shuffle_dinucleotide
from bio_arch.provenance import SeedManager


def fit_zipf_power_law(sequence: str, k: int = 4) -> tuple[float, float, float]:
    """Fit rank-frequency distribution of k-mers to log10(freq) = log10(C) - alpha * log10(rank).

    Returns (alpha_slope, R2_goodness_of_fit, total_unique_kmers).
    """
    n = len(sequence)
    if n < k:
        return 0.0, 0.0, 0

    kmers = [sequence[i : i + k] for i in range(n - k + 1)]
    counts = Counter(kmers)
    sorted_freqs = sorted(counts.values(), reverse=True)
    num_unique = len(sorted_freqs)

    if num_unique < 3:
        return 0.0, 0.0, num_unique

    # Log10-transform ranks and frequencies
    x = [math.log10(r) for r in range(1, num_unique + 1)]
    y = [math.log10(f) for f in sorted_freqs]

    # Linear regression: y = slope * x + intercept
    n_pts = len(x)
    sum_x = sum(x)
    sum_y = sum(y)
    sum_xx = sum(xi * xi for xi in x)
    sum_yy = sum(yi * yi for yi in y)
    sum_xy = sum(xi * yi for xi, yi in zip(x, y))

    denom = (n_pts * sum_xx - sum_x * sum_x)
    if abs(denom) < 1e-12:
        return 0.0, 0.0, num_unique

    slope = (n_pts * sum_xy - sum_x * sum_y) / denom
    intercept = (sum_y - slope * sum_x) / n_pts

    # R^2 calculation
    y_mean = sum_y / n_pts
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)
    ss_res = sum((yi - (slope * xi + intercept)) ** 2 for xi, yi in zip(x, y))

    r2 = 1.0 - (ss_res / max(1e-12, ss_tot)) if ss_tot > 0 else 0.0
    alpha = -slope  # Standard Zipf exponent is positive

    return round(alpha, 4), round(r2, 4), num_unique


def run_zipf_experiment(
    corpus_dir: Path,
    n_shuffles: int = 500,
    seed: int = 42,
) -> dict:
    seed_mgr = SeedManager(seed)
    target_accessions = ["NC_003977.2", "NC_005148.1", "NC_001401.2", "NC_001422.1", "NC_001669.1"]

    experiment_results = {}

    for acc in target_accessions:
        gbk_file = corpus_dir / f"{acc}.gbk"
        if not gbk_file.is_file():
            continue

        rec = SeqIO.read(str(gbk_file), "genbank")
        seq = str(rec.seq).upper()
        organism = rec.annotations.get("organism", rec.description)

        k_results = {}
        for k in (3, 4, 5, 6):
            nat_alpha, nat_r2, nat_unique = fit_zipf_power_law(seq, k=k)

            null_alphas = []
            null_r2s = []

            for s_idx in range(n_shuffles):
                sub_seed = seed_mgr.derive_seed(f"zipf_null_{acc}_k{k}_{s_idx}")
                rng = random.Random(sub_seed)
                shuffled = shuffle_dinucleotide(seq, rng=rng)

                alpha, r2, _ = fit_zipf_power_law(shuffled, k=k)
                null_alphas.append(alpha)
                null_r2s.append(r2)

            def calc_metric_stats(actual: float, null_dist: list[float]) -> dict:
                mean = sum(null_dist) / len(null_dist)
                var = sum((x - mean) ** 2 for x in null_dist) / max(1, len(null_dist) - 1)
                std = max(1e-6, var ** 0.5)
                z = round((actual - mean) / std, 2)
                diff = abs(actual - mean)
                p = round(sum(1 for x in null_dist if abs(x - mean) >= diff) / len(null_dist), 4)
                decision = "Reject H0 (Deviates from Markov-1)" if p < 0.01 else "Accept H0 (Consistent with Markov-1)"
                return {
                    "natural_value": actual,
                    "null_mean": round(mean, 4),
                    "null_std": round(std, 4),
                    "z_score": z,
                    "empirical_p_value": p,
                    "decision": decision,
                }

            k_results[f"k_{k}"] = {
                "k": k,
                "unique_kmers": nat_unique,
                "zipf_alpha": calc_metric_stats(nat_alpha, null_alphas),
                "zipf_r2": calc_metric_stats(nat_r2, null_r2s),
            }

        experiment_results[acc] = {
            "accession": acc,
            "organism": organism,
            "genome_length_bp": len(seq),
            "kmer_fits": k_results,
        }

    return experiment_results


if __name__ == "__main__":
    root = Path(__file__).parent.parent.parent
    c_dir = root / "data" / "mining_corpus"
    res = run_zipf_experiment(c_dir, n_shuffles=200)
    print(json.dumps(res, indent=2))
