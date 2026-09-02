"""Pre-Registered Hypothesis Test C: Push-Down Automaton Stability & Execution Dynamics.

Maps DNA codons to formal push-down stack automaton instructions to test whether
natural reading frames exhibit higher automaton stability or lower crash rates
than dinucleotide-preserving shuffled null models.
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
from bio_arch.modules.recompiler import CODON_TABLE
from bio_arch.provenance import SeedManager

# Canonical, unbiased mapping of 20 amino acids + stops to 5 PDA operations
OP_PUSH = "PUSH"
OP_POP = "POP"
OP_BRANCH = "BRANCH"
OP_NOP = "NOP"
OP_HALT = "HALT"

AA_TO_PDA_OP = {
    # Hydrophobic aliphatics -> PUSH
    "L": OP_PUSH, "I": OP_PUSH, "V": OP_PUSH, "F": OP_PUSH, "M": OP_PUSH,
    # Charged residues -> POP
    "K": OP_POP, "R": OP_POP, "D": OP_POP, "E": OP_POP, "H": OP_POP,
    # Polar & turns -> BRANCH
    "S": OP_BRANCH, "T": OP_BRANCH, "N": OP_BRANCH, "Q": OP_BRANCH, "P": OP_BRANCH, "G": OP_BRANCH,
    # Structural & aromatic -> NOP
    "A": OP_NOP, "C": OP_NOP, "W": OP_NOP, "Y": OP_NOP,
    # Stop codons -> HALT
    "*": OP_HALT,
}


def codon_to_op(codon: str) -> str:
    aa = CODON_TABLE.get(codon, "*")
    return AA_TO_PDA_OP.get(aa, OP_NOP)


def execute_pda(sequence: str, frame_offset: int = 0, max_steps: int = 2000, max_stack: int = 100) -> dict:
    """Execute sequence as a linear push-down automaton program."""
    subseq = sequence[frame_offset:]
    codons = [subseq[i : i + 3] for i in range(0, len(subseq) - 2, 3)]
    num_codons = len(codons)

    stack = []
    pc = 0
    steps = 0
    stack_depths = [0]
    op_history = []
    terminated_cleanly = False
    underflow_crashes = 0
    overflow_crashes = 0

    while pc < num_codons and steps < max_steps:
        steps += 1
        codon = codons[pc]
        op = codon_to_op(codon)
        op_history.append(op)

        if op == OP_HALT:
            terminated_cleanly = True
            break

        elif op == OP_PUSH:
            if len(stack) < max_stack:
                stack.append(pc % 16)
            else:
                overflow_crashes += 1
            pc += 1

        elif op == OP_POP:
            if stack:
                stack.pop()
            else:
                underflow_crashes += 1
            pc += 1

        elif op == OP_BRANCH:
            # Branch forward or backward based on top of stack
            top_val = stack[-1] if stack else 0
            offset = (top_val % 7) - 3  # offset between -3 and +3 codons
            if offset == 0:
                pc += 1
            else:
                pc = max(0, pc + offset)

        elif op == OP_NOP:
            pc += 1

        stack_depths.append(len(stack))

    # Operation entropy
    counts = Counter(op_history)
    n_ops = len(op_history)
    op_entropy = -sum((c / n_ops) * math.log2(c / n_ops) for c in counts.values()) if n_ops > 0 else 0.0

    return {
        "steps_executed": steps,
        "terminated_cleanly": 1.0 if terminated_cleanly else 0.0,
        "underflow_crashes": underflow_crashes,
        "overflow_crashes": overflow_crashes,
        "max_stack_depth": max(stack_depths),
        "mean_stack_depth": round(sum(stack_depths) / len(stack_depths), 2),
        "operation_entropy": round(op_entropy, 4),
    }


def run_automaton_experiment(
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

        # Measure natural across 3 frames
        frame_metrics = {}
        for frame in (0, 1, 2):
            nat_res = execute_pda(seq, frame_offset=frame)

            null_term = []
            null_max_stack = []
            null_mean_stack = []
            null_entropy = []
            null_underflows = []

            for s_idx in range(n_shuffles):
                sub_seed = seed_mgr.derive_seed(f"pda_null_{acc}_f{frame}_{s_idx}")
                rng = random.Random(sub_seed)
                shuffled = shuffle_dinucleotide(seq, rng=rng)

                shuf_res = execute_pda(shuffled, frame_offset=frame)
                null_term.append(shuf_res["terminated_cleanly"])
                null_max_stack.append(shuf_res["max_stack_depth"])
                null_mean_stack.append(shuf_res["mean_stack_depth"])
                null_entropy.append(shuf_res["operation_entropy"])
                null_underflows.append(shuf_res["underflow_crashes"])

            def calc_metric(actual: float, null_dist: list[float]) -> dict:
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

            frame_metrics[f"frame_{frame}"] = {
                "frame": frame,
                "clean_termination": calc_metric(nat_res["terminated_cleanly"], null_term),
                "max_stack_depth": calc_metric(nat_res["max_stack_depth"], null_max_stack),
                "mean_stack_depth": calc_metric(nat_res["mean_stack_depth"], null_mean_stack),
                "operation_entropy": calc_metric(nat_res["operation_entropy"], null_entropy),
                "underflow_crashes": calc_metric(nat_res["underflow_crashes"], null_underflows),
            }

        experiment_results[acc] = {
            "accession": acc,
            "organism": organism,
            "genome_length_bp": len(seq),
            "frame_execution": frame_metrics,
        }

    return experiment_results


if __name__ == "__main__":
    root = Path(__file__).parent.parent.parent
    c_dir = root / "data" / "mining_corpus"
    res = run_automaton_experiment(c_dir, n_shuffles=200)
    print(json.dumps(res, indent=2))
