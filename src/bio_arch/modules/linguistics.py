"""Module: DNA Linguistics & Bytecode Profiling.

Analyzes DNA/RNA as a formal information language and bytecode:
- Zipf's Law rank-frequency analysis (power-law exponent alpha).
- Lempel-Ziv modular dictionary phrase extraction (reusable subroutines).
- Vocabulary growth dynamics (Heaps' Law).
- Evaluates whether coding sequences mathematically resemble compiled machine code,
  human natural language, or stochastic noise.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import random
import sys
from typing import Any

from bio_arch.contracts import (
    AnalysisRun,
    EvidenceClass,
    Finding,
    InterpretationRecord,
    ModuleResult,
)
from bio_arch.logger import setup_logger
from bio_arch.provenance import (
    SeedManager,
    get_system_environment,
    now_iso,
)

logger = setup_logger("bio_arch.linguistics")


@dataclass
class ZipfFitResult:
    """Outcome of fitting Zipf's Law (f(r) ~ r^(-alpha)) to token frequencies."""

    alpha: float  # Power law exponent (English ~ 1.0, random noise ~ 0.0)
    r_squared: float  # Goodness of fit in log-log space
    vocabulary_size: int
    total_tokens: int
    ranked_tokens: list[dict[str, Any]]  # [{rank, token, count, freq}]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def extract_codons(sequence: str, frame: int = 0) -> list[str]:
    """Extract 3-letter codons from a specified reading frame (0, 1, or 2)."""
    seq = sequence.upper()
    codons = []
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        if len(codon) == 3:
            codons.append(codon)
    return codons


def fit_zipfs_law(token_counts: dict[str, int]) -> ZipfFitResult:
    """Fit Zipf's power law log(freq) = C - alpha * log(rank) using linear regression."""
    if not token_counts:
        return ZipfFitResult(0.0, 0.0, 0, 0, [])

    sorted_items = sorted(token_counts.items(), key=lambda x: x[1], reverse=True)
    total_tokens = sum(token_counts.values())
    vocab_size = len(sorted_items)

    ranked_tokens = []
    log_ranks = []
    log_freqs = []

    for rank, (token, count) in enumerate(sorted_items, start=1):
        freq = count / total_tokens
        log_r = math.log(rank)
        log_f = math.log(freq)
        log_ranks.append(log_r)
        log_freqs.append(log_f)
        ranked_tokens.append({
            "rank": rank,
            "token": token,
            "count": count,
            "frequency": round(freq, 6),
        })

    n = len(log_ranks)
    if n < 2:
        return ZipfFitResult(0.0, 1.0, vocab_size, total_tokens, ranked_tokens)

    # Linear regression: log_f = intercept - alpha * log_r
    mean_x = sum(log_ranks) / n
    mean_y = sum(log_freqs) / n

    ss_xx = sum((x - mean_x) ** 2 for x in log_ranks)
    ss_yy = sum((y - mean_y) ** 2 for y in log_freqs)
    ss_xy = sum((x - mean_x) * (y - mean_y) for x, y in zip(log_ranks, log_freqs))

    slope = ss_xy / ss_xx if ss_xx > 0 else 0.0
    alpha = -slope  # By convention, f(r) ~ r^(-alpha)

    # R-squared
    r_squared = (ss_xy ** 2) / (ss_xx * ss_yy) if (ss_xx * ss_yy) > 0 else 0.0

    return ZipfFitResult(
        alpha=round(alpha, 4),
        r_squared=round(r_squared, 4),
        vocabulary_size=vocab_size,
        total_tokens=total_tokens,
        ranked_tokens=ranked_tokens,
    )


def extract_lz_subroutines(sequence: str, min_len: int = 3) -> dict[str, int]:
    """Lempel-Ziv-78 phrase extraction identifying recurring sequence motifs ('subroutines')."""
    seq = sequence.upper()
    n = len(seq)
    dictionary: set[str] = set()
    phrase_counts: dict[str, int] = defaultdict(int)

    i = 0
    while i < n:
        j = i + 1
        while j <= n and seq[i:j] in dictionary:
            j += 1
        phrase = seq[i:j]
        dictionary.add(phrase)
        if len(phrase) >= min_len:
            phrase_counts[phrase] += 1
        i = j

    # Count actual total occurrences in the entire sequence for high-frequency phrases
    reusable_subroutines: dict[str, int] = {}
    for phrase in phrase_counts:
        # Find true total occurrences
        count = seq.count(phrase)
        if count >= 2:
            reusable_subroutines[phrase] = count

    # Sort descending by count * length (information weight)
    sorted_subroutines = dict(
        sorted(reusable_subroutines.items(), key=lambda x: x[1] * len(x[0]), reverse=True)
    )
    return sorted_subroutines


def analyze_linguistic_architecture(
    sequence: str,
    record_id: str = "seq_01",
    seed: int = 42,
) -> tuple[dict[str, Any], list[Finding], list[InterpretationRecord]]:
    """Perform linguistic and bytecode profiling on a biological sequence."""
    seq = sequence.upper()

    # 1. Codon opcodes (frame 0)
    codons = extract_codons(seq, frame=0)
    codon_counts = Counter(codons)
    codon_zipf = fit_zipfs_law(dict(codon_counts))

    # 2. k-mer vocabulary at k=4 (short instruction blocks)
    k4_counts: dict[str, int] = defaultdict(int)
    for i in range(len(seq) - 3):
        k4_counts[seq[i : i + 4]] += 1
    k4_zipf = fit_zipfs_law(k4_counts)

    # 3. Lempel-Ziv modular subroutines
    subroutines = extract_lz_subroutines(seq, min_len=4)
    top_subroutines = list(subroutines.items())[:10]

    # 4. Randomized null control for Zipf exponent
    rng = random.Random(seed)
    shuffled_seq = list(seq)
    rng.shuffle(shuffled_seq)
    shuffled_str = "".join(shuffled_seq)
    shuffled_codons = Counter(extract_codons(shuffled_str, frame=0))
    null_zipf = fit_zipfs_law(dict(shuffled_codons))

    findings: list[Finding] = [
        Finding(
            finding_id=f"f_{record_id}_codon_zipf",
            metric="codon_zipf_exponent_alpha",
            observed_value=codon_zipf.alpha,
            control_distribution={
                "null_model": "mononucleotide_shuffle",
                "null_alpha": null_zipf.alpha,
                "r_squared": codon_zipf.r_squared,
            },
            effect_size=round(codon_zipf.alpha - null_zipf.alpha, 4),
            biological_context=f"Codon rank-frequency power law exponent for {record_id}.",
        ),
        Finding(
            finding_id=f"f_{record_id}_modular_subroutines",
            metric="reusable_subroutine_count",
            observed_value=len(subroutines),
            biological_context=f"Count of repeating LZ-78 multi-nucleotide subroutines (length >= 4) in {record_id}.",
        ),
    ]

    interpretations: list[InterpretationRecord] = [
        InterpretationRecord(
            finding_ids=[f"f_{record_id}_codon_zipf"],
            classification=EvidenceClass.INTERPRETATION,
            claim=f"{EvidenceClass.INTERPRETATION.required_prefix} that a codon Zipf exponent of alpha={codon_zipf.alpha} (vs null alpha={null_zipf.alpha}) reflects non-random codon usage bias and translational efficiency constraints.",
            alternatives=[
                "Biased GC composition naturally skews codon availability without linguistic intent.",
                "tRNA abundance in host organism driving codon selection pressure.",
            ],
            limitations=[
                "Zipf's law is a necessary but not sufficient condition for formal human-like or programming language.",
            ],
            proposed_test="Benchmark against highly expressed ribosomal genes vs rarely expressed viral regulators.",
        )
    ]

    outputs = {
        "record_id": record_id,
        "sequence_length": len(seq),
        "total_codons": codon_zipf.total_tokens,
        "codon_vocabulary_size": codon_zipf.vocabulary_size,
        "codon_zipf_alpha": codon_zipf.alpha,
        "codon_zipf_r2": codon_zipf.r_squared,
        "k4_zipf_alpha": k4_zipf.alpha,
        "null_shuffled_alpha": null_zipf.alpha,
        "reusable_subroutines_found": len(subroutines),
        "top_subroutines": top_subroutines,
        "codon_rank_distribution": codon_zipf.ranked_tokens[:20],
    }

    return outputs, findings, interpretations


def run_linguistics(
    input_file: Path | str,
    output_dir: Path | str | None = None,
    seed: int = 42,
) -> ModuleResult:
    """Run linguistic and bytecode profiling on an input FASTA file."""
    in_path = Path(input_file)
    if not in_path.is_file():
        raise FileNotFoundError(f"Input file not found: {in_path}")

    seed_mgr = SeedManager(master_seed=seed)
    mod_seed = seed_mgr.derive_seed("linguistics_profiler")

    run_meta = AnalysisRun(
        run_id=f"run_ling_{in_path.stem}_{mod_seed}",
        timestamp=now_iso(),
        module="linguistics_profiler",
        version="0.1.0",
        input_ids=[in_path.name],
        parameters={"seed": mod_seed},
        seed=mod_seed,
        environment=get_system_environment(),
        status="running",
    )

    # Read fasta records
    records: list[tuple[str, str]] = []
    curr_id = None
    curr_seq: list[str] = []
    with in_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if curr_id and curr_seq:
                    records.append((curr_id, "".join(curr_seq)))
                curr_id = line[1:].split()[0]
                curr_seq = []
            elif curr_id:
                curr_seq.append(line)
        if curr_id and curr_seq:
            records.append((curr_id, "".join(curr_seq)))

    all_findings: list[Finding] = []
    all_interpretations: list[InterpretationRecord] = []
    results_list: list[dict[str, Any]] = []

    for rec_id, seq in records:
        out, findings, interps = analyze_linguistic_architecture(seq, record_id=rec_id, seed=mod_seed)
        results_list.append(out)
        all_findings.extend(findings)
        all_interpretations.extend(interps)

    run_meta.status = "success"

    artifact_paths: list[str] = []
    if output_dir:
        out_p = Path(output_dir)
        out_p.mkdir(parents=True, exist_ok=True)
        out_file = out_p / f"{in_path.stem}_linguistics_profile.json"
        payload = {
            "run_metadata": run_meta.to_dict(),
            "results": results_list,
            "findings": [f.to_dict() for f in all_findings],
            "interpretations": [i.to_dict() for i in all_interpretations],
        }
        out_file.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        artifact_paths.append(str(out_file))

    return ModuleResult(
        run_metadata=run_meta,
        outputs={"records_analyzed": len(records), "profiles": results_list},
        findings=all_findings,
        interpretations=all_interpretations,
        warnings=[],
        errors=[],
        artifact_paths=artifact_paths,
    )


def main() -> None:
    """CLI entrypoint for linguistics profiler."""
    parser = argparse.ArgumentParser(
        description="Bio-Arch DNA Linguistics & Bytecode Profiler",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("file", type=str, help="Path to input FASTA file")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for null comparison")
    parser.add_argument("--outdir", type=str, default="outputs/linguistics", help="Output directory")

    args = parser.parse_args()

    result = run_linguistics(input_file=args.file, output_dir=args.outdir, seed=args.seed)

    print("\n--- DNA Linguistics & Bytecode Results ---")
    for profile in result.outputs["profiles"]:
        print(f"Record: {profile['record_id']} ({profile['sequence_length']} bp)")
        print(f"  Codon Zipf Exponent alpha: {profile['codon_zipf_alpha']} (R^2={profile['codon_zipf_r2']})")
        print(f"  Null Shuffled alpha:      {profile['null_shuffled_alpha']}")
        print(f"  Reusable Subroutines:     {profile['reusable_subroutines_found']}")
        if profile["top_subroutines"]:
            top_str = ", ".join(f"{phrase} ({cnt}x)" for phrase, cnt in profile["top_subroutines"][:5])
            print(f"  Top Subroutines:          {top_str}")
    print("------------------------------------------\n")


if __name__ == "__main__":
    main()
