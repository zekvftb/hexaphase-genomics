"""Pre-Registered Hypothesis Test A: Kolmogorov & Lempel-Ziv Compressibility.

Evaluates whether real viral genomes exhibit non-random algorithmic compressibility
compared against N=500 dinucleotide-preserving shuffled controls.
"""

from __future__ import annotations

from collections import Counter
import io
import json
import math
from pathlib import Path
import random
import sys
import zlib

_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from Bio import SeqIO

from bio_arch.modules.information import shuffle_dinucleotide
from bio_arch.provenance import SeedManager, now_iso


def lz76_complexity(sequence: str) -> int:
    """Compute the standard Lempel-Ziv (1976) production complexity.

    Counts the number of unique substrings generated sequentially from left to right.
    """
    n = len(sequence)
    if n == 0:
        return 0

    history = set()
    phrases = 0
    curr = ""

    for char in sequence:
        curr += char
        if curr not in history:
            history.add(curr)
            phrases += 1
            curr = ""

    if curr:
        phrases += 1

    return phrases


def zlib_compression_ratio(sequence: str) -> float:
    """Compute DEFLATE/LZW compression ratio (compressed_bytes / raw_bytes)."""
    raw_bytes = sequence.encode("ascii")
    comp_bytes = zlib.compress(raw_bytes, level=9)
    return round(len(comp_bytes) / max(1, len(raw_bytes)), 4)


def sliding_window_shannon_entropy(sequence: str, window: int = 64) -> float:
    """Compute mean Shannon entropy across sliding windows of length w."""
    n = len(sequence)
    if n < window:
        return 0.0

    entropies = []
    for i in range(0, n - window + 1, window // 2):
        win = sequence[i : i + window]
        counts = Counter(win)
        w_len = len(win)
        ent = -sum((cnt / w_len) * math.log2(cnt / w_len) for cnt in counts.values())
        entropies.append(ent)

    return round(sum(entropies) / len(entropies), 4)


def run_compression_experiment(
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

        # 1. Natural metrics
        nat_lz76 = lz76_complexity(seq)
        nat_comp_ratio = zlib_compression_ratio(seq)
        nat_ent_64 = sliding_window_shannon_entropy(seq, window=64)
        nat_ent_128 = sliding_window_shannon_entropy(seq, window=128)

        # 2. Null distribution via N=500 Altschul-Erickson shuffles
        null_lz76 = []
        null_comp = []
        null_ent_64 = []
        null_ent_128 = []

        for k in range(n_shuffles):
            sub_seed = seed_mgr.derive_seed(f"lz_null_{acc}_{k}")
            rng = random.Random(sub_seed)
            shuffled = shuffle_dinucleotide(seq, rng=rng)

            null_lz76.append(lz76_complexity(shuffled))
            null_comp.append(zlib_compression_ratio(shuffled))
            null_ent_64.append(sliding_window_shannon_entropy(shuffled, window=64))
            null_ent_128.append(sliding_window_shannon_entropy(shuffled, window=128))

        def calc_stats(actual: float, null_dist: list[float]) -> dict:
            mean = sum(null_dist) / len(null_dist)
            var = sum((x - mean) ** 2 for x in null_dist) / max(1, len(null_dist) - 1)
            std = max(1e-6, var ** 0.5)
            z = round((actual - mean) / std, 2)
            # Two-tailed empirical p-value
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

        experiment_results[acc] = {
            "accession": acc,
            "organism": organism,
            "genome_length_bp": len(seq),
            "lz76_complexity": calc_stats(nat_lz76, null_lz76),
            "compression_ratio_zlib": calc_stats(nat_comp_ratio, null_comp),
            "entropy_w64": calc_stats(nat_ent_64, null_ent_64),
            "entropy_w128": calc_stats(nat_ent_128, null_ent_128),
        }

    return experiment_results


if __name__ == "__main__":
    root = Path(__file__).parent.parent.parent
    c_dir = root / "data" / "mining_corpus"
    res = run_compression_experiment(c_dir, n_shuffles=200)
    print(json.dumps(res, indent=2))
