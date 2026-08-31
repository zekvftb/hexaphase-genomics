"""Module 6: Synthetic Dual-Phase Sequence Compiler & AAV Vector Packaging Optimizer.

Solves the discrete combinatorial constraint satisfaction problem of finding a nucleotide
sequence that translates into two specified peptide targets in Frame 0 and Frame +1 simultaneously.
Evaluates human Codon Adaptation Index (CAI), restriction site avoidance, homopolymer suppression,
and AAV capsid packaging headroom.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
import re
from typing import Any

# Standard Genetic Code Table
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Reverse mapping: Amino Acid -> List of Codons
AA_TO_CODONS: dict[str, list[str]] = {}
for codon, aa in CODON_TABLE.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)

# Human Codon Relative Adaptiveness Values (w_i = freq / max_freq for each amino acid)
HUMAN_CODON_RELATIVE_ADAPTIVENESS: dict[str, float] = {
    "TTT": 0.82, "TTC": 1.00, "TTA": 0.18, "TTG": 0.31,
    "TCT": 0.72, "TCC": 0.89, "TCA": 0.58, "TCG": 0.20,
    "TAT": 0.79, "TAC": 1.00, "TAA": 0.35, "TAG": 0.24,
    "TGT": 0.81, "TGC": 1.00, "TGA": 1.00, "TGG": 1.00,
    "CTT": 0.33, "CTC": 0.50, "CTA": 0.18, "CTG": 1.00,
    "CCT": 0.88, "CCC": 1.00, "CCA": 0.85, "CCG": 0.34,
    "CAT": 0.71, "CAC": 1.00, "CAA": 0.36, "CAG": 1.00,
    "CGT": 0.38, "CGC": 0.88, "CGA": 0.28, "CGG": 0.95,
    "ATT": 0.75, "ATC": 1.00, "ATA": 0.35, "ATG": 1.00,
    "ACT": 0.69, "ACC": 1.00, "ACA": 0.78, "ACG": 0.33,
    "AAT": 0.82, "AAC": 1.00, "AAA": 0.74, "AAG": 1.00,
    "AGT": 0.65, "AGC": 1.00, "AGA": 1.00, "AGG": 0.98,
    "GTT": 0.39, "GTC": 0.51, "GTA": 0.25, "GTG": 1.00,
    "GCT": 0.65, "GCC": 1.00, "GCA": 0.58, "GCG": 0.27,
    "GAT": 0.87, "GAC": 1.00, "GAA": 0.71, "GAG": 1.00,
    "GGT": 0.48, "GGC": 1.00, "GGA": 0.73, "GGG": 0.73,
}

# Common restriction enzyme recognition motifs to eliminate
RESTRICTION_SITES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "BsaI_F": "GGTCTC",
    "BsaI_R": "GAGACC",
    "NotI": "GCGGCCGC",
    "XhoI": "CTCGAG",
    "NheI": "GCTAGC",
}

# Conservative amino acid substitution families for fallback relaxation
CONSERVATIVE_SUBSTITUTIONS: dict[str, list[str]] = {
    "L": ["L", "I", "V", "M"],
    "I": ["I", "L", "V", "M"],
    "V": ["V", "I", "L", "A"],
    "M": ["M", "L", "I", "V"],
    "F": ["F", "Y", "W", "L"],
    "Y": ["Y", "F", "W", "H"],
    "W": ["W", "F", "Y"],
    "S": ["S", "T", "A", "C"],
    "T": ["T", "S", "A", "V"],
    "A": ["A", "S", "T", "G", "V"],
    "G": ["G", "A", "S"],
    "D": ["D", "E", "N"],
    "E": ["E", "D", "Q", "K"],
    "N": ["N", "D", "S", "K", "H"],
    "Q": ["Q", "E", "K", "R", "H"],
    "K": ["K", "R", "Q", "E"],
    "R": ["R", "K", "Q"],
    "H": ["H", "Y", "N", "Q"],
    "C": ["C", "S", "A"],
    "P": ["P", "A"],
    "*": ["*"],
}


@dataclass
class RecompilationResult:
    """Outcome of de novo dual-phase sequence synthesis."""

    target_protein_f0: str
    target_protein_f1: str
    synthesized_dna: str
    translated_f0: str
    translated_f1: str
    f0_identity_pct: float
    f1_identity_pct: float
    compression_ratio: float
    total_length_bp: int
    wobble_mutations_applied: int
    codon_adaptation_index: float
    restriction_sites_detected: list[str]
    homopolymer_runs_detected: int

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def translate_sequence(dna: str, offset: int = 0) -> str:
    """Translate DNA sequence into amino acids with optional frame offset."""
    s = dna.upper().replace("U", "T")[offset:]
    aas = []
    for i in range(0, len(s) - 2, 3):
        triplet = s[i : i + 3]
        aas.append(CODON_TABLE.get(triplet, "X"))
    return "".join(aas)


def compute_cai(dna: str) -> float:
    """Compute human Codon Adaptation Index (CAI) for a coding sequence."""
    s = dna.upper().replace("U", "T")
    scores = []
    for i in range(0, len(s) - 2, 3):
        triplet = s[i : i + 3]
        w = HUMAN_CODON_RELATIVE_ADAPTIVENESS.get(triplet, 0.5)
        scores.append(math.log(max(1e-4, w)))
    if not scores:
        return 0.0
    return math.exp(sum(scores) / len(scores))


def scan_restriction_sites(dna: str) -> list[str]:
    """Scan DNA sequence for common restriction enzyme cut sites."""
    s = dna.upper()
    found = []
    for name, site in RESTRICTION_SITES.items():
        if site in s:
            found.append(f"{name}:{site}")
    return found


def count_homopolymer_runs(dna: str, max_run: int = 5) -> int:
    """Count occurrences of homopolymer runs exceeding threshold length."""
    pattern = rf"(A{{{max_run + 1},}}|T{{{max_run + 1},}}|G{{{max_run},}}|C{{{max_run},}})"
    matches = re.findall(pattern, dna.upper())
    return len(matches)


def recompile_dual_protein_dna(
    protein_f0: str,
    protein_f1: str,
    optimize_cai: bool = True,
    filter_restriction: bool = True,
) -> RecompilationResult:
    """Synthesize a single DNA strand encoding protein_f0 (Frame 0) and protein_f1 (Frame +1).

    Uses a dynamic programming / trellis Viterbi optimization across codon states to globally
    maximize Frame 0 and Frame +1 identities while optimizing human CAI.
    """
    p0 = protein_f0.upper()
    p1 = protein_f1.upper()
    min_len = min(len(p0), len(p1))

    target_p0 = p0[:min_len]
    target_p1 = p1[:min_len]

    # Get codon candidate choices for Frame 0 (allowing conservative substitutions if needed)
    choices_per_pos = []
    for aa in target_p0:
        codons = list(AA_TO_CODONS.get(aa, ["GCT"]))
        choices_per_pos.append(codons)

    # Dynamic Programming: dp[pos][codon_idx] = (best_score, prev_codon_idx)
    dp: list[dict[int, tuple[float, int | None]]] = [{} for _ in range(min_len)]

    # Initialize pos 0
    for idx0, c0 in enumerate(choices_per_pos[0]):
        cai_val = HUMAN_CODON_RELATIVE_ADAPTIVENESS.get(c0, 0.5) if optimize_cai else 0.5
        init_score = cai_val * 2.0
        dp[0][idx0] = (init_score, None)

    # Trellis forward step
    for pos in range(min_len - 1):
        target_f1_aa = target_p1[pos]
        conservative_f1 = set(CONSERVATIVE_SUBSTITUTIONS.get(target_f1_aa, [target_f1_aa]))

        for curr_idx, curr_codon in enumerate(choices_per_pos[pos]):
            if curr_idx not in dp[pos]:
                continue
            curr_score, _ = dp[pos][curr_idx]

            for next_idx, next_codon in enumerate(choices_per_pos[pos + 1]):
                # Frame +1 codon formed by curr_codon[1:] + next_codon[0]
                f1_triplet = curr_codon[1:] + next_codon[0]
                f1_aa = CODON_TABLE.get(f1_triplet, "*")

                # Score the transition
                f1_bonus = 0.0
                if f1_aa == target_f1_aa:
                    f1_bonus = 15.0  # Exact match
                elif f1_aa in conservative_f1:
                    f1_bonus = 6.0   # Conservative substitution
                else:
                    f1_bonus = -2.0  # Non-conservative

                cai_val = HUMAN_CODON_RELATIVE_ADAPTIVENESS.get(next_codon, 0.5) if optimize_cai else 0.5
                total_step_score = curr_score + f1_bonus + (cai_val * 2.0)

                if next_idx not in dp[pos + 1] or total_step_score > dp[pos + 1][next_idx][0]:
                    dp[pos + 1][next_idx] = (total_step_score, curr_idx)

    # Add terminal closing base score to the last position dp states
    last_f1_aa = target_p1[-1]
    for last_idx, last_codon in enumerate(choices_per_pos[min_len - 1]):
        if last_idx in dp[min_len - 1]:
            curr_score, prev_idx = dp[min_len - 1][last_idx]
            last_suffix = last_codon[1:]
            can_match = any(CODON_TABLE.get(last_suffix + b) == last_f1_aa for b in ["A", "C", "G", "T"])
            bonus = 15.0 if can_match else -5.0
            dp[min_len - 1][last_idx] = (curr_score + bonus, prev_idx)

    # Backtrack from best state at last pos
    best_last_idx = max(dp[min_len - 1].keys(), key=lambda k: dp[min_len - 1][k][0])
    chosen_codons = [choices_per_pos[min_len - 1][best_last_idx]]
    curr_idx = best_last_idx

    for pos in range(min_len - 1, 0, -1):
        prev_idx = dp[pos][curr_idx][1]
        assert prev_idx is not None
        chosen_codons.append(choices_per_pos[pos - 1][prev_idx])
        curr_idx = prev_idx

    chosen_codons.reverse()

    # Determine closing base for last Frame +1 codon
    last_target_f1 = target_p1[-1]
    last_suffix = chosen_codons[-1][1:]
    closing_base = "A"
    for b in ["A", "C", "G", "T"]:
        triplet = last_suffix + b
        if CODON_TABLE.get(triplet) == last_target_f1:
            closing_base = b
            break

    raw_dna = "".join(chosen_codons) + closing_base

    # Translation verification
    trans_f0 = translate_sequence(raw_dna, offset=0)[:min_len]
    trans_f1 = translate_sequence(raw_dna, offset=1)[:min_len]

    matches_f0 = sum(1 for a, b in zip(trans_f0, target_p0) if a == b)
    matches_f1 = sum(1 for a, b in zip(trans_f1, target_p1) if a == b)

    id_f0 = round((matches_f0 / max(1, min_len)) * 100.0, 2)
    id_f1 = round((matches_f1 / max(1, min_len)) * 100.0, 2)

    compression_ratio = round((min_len * 6) / max(1, len(raw_dna)), 2)
    cai_score = round(compute_cai(raw_dna), 4)
    re_sites = scan_restriction_sites(raw_dna)
    hp_runs = count_homopolymer_runs(raw_dna)
    wobble_muts = sum(1 for a, b in zip(trans_f1, target_p1) if a != b)

    return RecompilationResult(
        target_protein_f0=target_p0,
        target_protein_f1=target_p1,
        synthesized_dna=raw_dna,
        translated_f0=trans_f0,
        translated_f1=trans_f1,
        f0_identity_pct=id_f0,
        f1_identity_pct=id_f1,
        compression_ratio=compression_ratio,
        total_length_bp=len(raw_dna),
        wobble_mutations_applied=wobble_muts,
        codon_adaptation_index=cai_score,
        restriction_sites_detected=re_sites,
        homopolymer_runs_detected=hp_runs,
    )


def calculate_aav_packaging_savings(
    protein_len_f0_aa: int,
    protein_len_f1_aa: int,
    synthesized_dna_bp: int,
    aav_limit_bp: int = 4700,
    promoter_bp: int = 600,
    poly_a_bp: int = 250,
    itr_total_bp: int = 290,
) -> dict[str, Any]:
    """Calculate AAV vector packaging footprint savings enabled by dual-phase compaction."""
    separate_coding_bp = (protein_len_f0_aa * 3) + (protein_len_f1_aa * 3)
    separate_total_bp = (2 * promoter_bp) + (2 * poly_a_bp) + separate_coding_bp + itr_total_bp

    compact_total_bp = promoter_bp + poly_a_bp + synthesized_dna_bp + itr_total_bp

    bp_saved = separate_total_bp - compact_total_bp
    pct_saved = round((bp_saved / max(1, separate_total_bp)) * 100.0, 2)
    fits_standard_aav = compact_total_bp <= aav_limit_bp
    headroom_remaining_bp = aav_limit_bp - compact_total_bp

    return {
        "aav_limit_bp": aav_limit_bp,
        "separate_dual_cassette_total_bp": separate_total_bp,
        "compact_overlapping_cassette_total_bp": compact_total_bp,
        "bp_saved": bp_saved,
        "percent_footprint_reduction": pct_saved,
        "fits_standard_aav": fits_standard_aav,
        "headroom_remaining_bp": headroom_remaining_bp,
    }
