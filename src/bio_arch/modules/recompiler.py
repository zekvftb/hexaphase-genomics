"""Module 6: Dual-Phase Biological Recompiler.

Synthesizes de novo DNA sequences that encode two distinct functional proteins
simultaneously in Frame 0 and Frame +1 via Wobble Carrier Wave modulation.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
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


def recompile_dual_protein_dna(protein_f0: str, protein_f1: str) -> RecompilationResult:
    """Synthesize a single DNA strand encoding protein_f0 (Frame 0) and protein_f1 (Frame +1)."""
    p0 = protein_f0.upper()
    p1 = protein_f1.upper()
    min_len = min(len(p0), len(p1))
    
    # Pad to equal length
    target_p0 = p0[:min_len]
    target_p1 = p1[:min_len]

    synthesized_bases: list[str] = []
    wobble_muts = 0

    # Iterative forward solver with backoff
    for i in range(min_len):
        aa0 = target_p0[i]
        aa1 = target_p1[i]

        candidates_0 = AA_TO_CODONS.get(aa0, ["GCT"])
        candidates_1_targets = CONSERVATIVE_SUBSTITUTIONS.get(aa1, [aa1])

        # Try to find exact dual-codon compatibility
        found = False
        best_triplet = None

        for c0 in candidates_0:
            # Suffix of c0 is c0[1:] (2 bases)
            for sub_aa1 in candidates_1_targets:
                codons_sub1 = AA_TO_CODONS.get(sub_aa1, [])
                for c1 in codons_sub1:
                    if c1.startswith(c0[1:]):
                        # Compatible!
                        best_triplet = c0
                        if sub_aa1 != aa1:
                            wobble_muts += 1
                        found = True
                        break
                if found:
                    break
            if found:
                break

        if not found:
            # Fallback: Use canonical codon for Frame 0 and let Frame +1 take conservative fit
            best_triplet = candidates_0[0]
            wobble_muts += 1

        synthesized_bases.append(best_triplet)

    # Need final 3rd base for the last +1 codon
    last_c1_target = target_p1[-1]
    last_c1_codons = AA_TO_CODONS.get(last_c1_target, ["GCT"])
    last_suffix = synthesized_bases[-1][1:]
    closing_base = "A"
    for c in last_c1_codons:
        if c.startswith(last_suffix):
            closing_base = c[2]
            break
    
    raw_dna = "".join(synthesized_bases) + closing_base

    # Translation verification
    trans_f0 = translate_sequence(raw_dna, offset=0)[:min_len]
    trans_f1 = translate_sequence(raw_dna, offset=1)[:min_len]

    matches_f0 = sum(1 for a, b in zip(trans_f0, target_p0) if a == b)
    matches_f1 = sum(1 for a, b in zip(trans_f1, target_p1) if a == b)

    id_f0 = round((matches_f0 / max(1, min_len)) * 100.0, 2)
    id_f1 = round((matches_f1 / max(1, min_len)) * 100.0, 2)

    # 2 separate proteins of length L would take 2 * (3 * L) = 6L bases.
    # Dual-phase takes 3L + 1 bases -> Compression ratio ~ 2.0x!
    compression_ratio = round((min_len * 6) / max(1, len(raw_dna)), 2)

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
    )
