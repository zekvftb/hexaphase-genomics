"""Module: Symbolic Genomic Disassembler & Sequence Annotation Parser.

Maps canonical biological sequence features to symbolic intermediate representation (IR) opcodes:
- Promoters (-35 / -10 Pribnow consensus): Transcription initiation sites.
- Operators / Transcription Factor Binding Sites: Regulatory control regions.
- Ribosome Binding Sites (Shine-Dalgarno): Translation initiation signals.
- Open Reading Frames (ORFs): Protein-coding sequence regions.
- Stop Codons & Intrinsic Terminators: Translation and transcription termination signals.

Emits human-readable symbolic assembly listings (.asm) for structural visualization.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
import re
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
    get_system_environment,
    now_iso,
)

logger = setup_logger("bio_arch.disassembler")

# Genetic Code translation table (Standard Code 11)
GENETIC_CODE = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "STOP", "TAG": "STOP",
    "TGT": "Cys", "TGC": "Cys", "TGA": "STOP", "TGG": "Trp",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}


@dataclass
class BiologicalToken:
    """A recognized syntactic opcode or token in the genetic program."""

    token_type: str  # PROMOTER, OPERATOR, RBS, START, CODON, STOP, TERMINATOR
    start: int
    end: int
    sequence: str
    label: str
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class DisassembledRoutine:
    """A decompiled biological function / operon routine."""

    routine_id: str
    start_offset: int
    end_offset: int
    length_bp: int
    tokens: list[BiologicalToken]
    assembly_listing: list[str]
    decompiled_pseudocode: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "routine_id": self.routine_id,
            "start_offset": self.start_offset,
            "end_offset": self.end_offset,
            "length_bp": self.length_bp,
            "tokens": [t.to_dict() for t in self.tokens],
            "assembly_listing": self.assembly_listing,
            "decompiled_pseudocode": self.decompiled_pseudocode,
        }


def hamming_dist(s1: str, s2: str) -> int:
    """Hamming distance between two equal-length strings."""
    return sum(1 for a, b in zip(s1, s2) if a != b)


def scan_tokens(sequence: str) -> list[BiologicalToken]:
    """Scan raw sequence for biological compiler opcodes and functional tokens."""
    seq = sequence.upper()
    tokens: list[BiologicalToken] = []

    # 1. Scan for -10 Pribnow Box consensus: TATAAT (allow 1 mismatch)
    for i in range(len(seq) - 5):
        sub = seq[i : i + 6]
        if hamming_dist(sub, "TATAAT") <= 1:
            tokens.append(
                BiologicalToken(
                    token_type="PROMOTER_MINUS10",
                    start=i,
                    end=i + 6,
                    sequence=sub,
                    label="Pribnow_box_-10",
                    metadata={"consensus": "TATAAT", "role": "Function entry point"},
                )
            )

    # 2. Scan for -35 Box consensus: TTGACA (allow 1 mismatch)
    for i in range(len(seq) - 5):
        sub = seq[i : i + 6]
        if hamming_dist(sub, "TTGACA") <= 1:
            tokens.append(
                BiologicalToken(
                    token_type="PROMOTER_MINUS35",
                    start=i,
                    end=i + 6,
                    sequence=sub,
                    label="Sigma70_box_-35",
                    metadata={"consensus": "TTGACA", "role": "Recognition anchor"},
                )
            )

    # 3. Scan for Shine-Dalgarno Ribosome Binding Sites (RBS): AGGAGG, GGAGG, AAGG
    rbs_patterns = [("AGGAGG", 0), ("GGAGG", 0), ("AGGA", 0)]
    for pat, max_mismatch in rbs_patterns:
        pat_len = len(pat)
        for i in range(len(seq) - pat_len):
            sub = seq[i : i + pat_len]
            if hamming_dist(sub, pat) <= max_mismatch:
                tokens.append(
                    BiologicalToken(
                        token_type="RBS_SHINE_DALGARNO",
                        start=i,
                        end=i + pat_len,
                        sequence=sub,
                        label=f"RBS_{pat}",
                        metadata={"consensus": pat, "role": "Hardware stack pointer setup"},
                    )
                )

    # 4. Scan for Open Reading Frames (ORFs): START (ATG) -> STOP (TAA, TAG, TGA)
    stop_codons = {"TAA", "TAG", "TGA"}
    for frame in (0, 1, 2):
        i = frame
        while i < len(seq) - 2:
            if seq[i : i + 3] == "ATG":
                start_pos = i
                # Follow codon stream until stop codon
                j = start_pos + 3
                found_stop = False
                while j < len(seq) - 2:
                    triplet = seq[j : j + 3]
                    if triplet in stop_codons:
                        found_stop = True
                        orf_len = j + 3 - start_pos
                        if orf_len >= 30:  # Minimum ORF threshold: 10 codons
                            tokens.append(
                                BiologicalToken(
                                    token_type="START_CODON",
                                    start=start_pos,
                                    end=start_pos + 3,
                                    sequence="ATG",
                                    label="INIT_Met",
                                    metadata={"frame": frame, "role": "Execution start"},
                                )
                            )
                            tokens.append(
                                BiologicalToken(
                                    token_type="STOP_CODON",
                                    start=j,
                                    end=j + 3,
                                    sequence=triplet,
                                    label=f"HALT_{triplet}",
                                    metadata={"frame": frame, "role": "Execution halt / return"},
                                )
                            )
                            tokens.append(
                                BiologicalToken(
                                    token_type="OPEN_READING_FRAME",
                                    start=start_pos,
                                    end=j + 3,
                                    sequence=seq[start_pos : j + 3],
                                    label=f"CDS_Frame_{frame}",
                                    metadata={
                                        "codons": orf_len // 3,
                                        "length_bp": orf_len,
                                        "role": "Bytecode instruction stream",
                                    },
                                )
                            )
                        break
                    j += 3
                if found_stop:
                    i = j + 3
                else:
                    i += 3
            else:
                i += 3

    # Sort tokens by starting offset
    tokens.sort(key=lambda t: t.start)
    return tokens


def generate_assembly_listing(sequence: str, tokens: list[BiologicalToken]) -> list[str]:
    """Emit an assembly-style disassembly listing mapping offsets to biological opcodes."""
    lines: list[str] = []
    lines.append("; ==========================================================================")
    lines.append("; BIOLOGICAL GENOMIC DISASSEMBLER - BYTECODE LISTING")
    lines.append("; ==========================================================================")

    for tok in tokens:
        addr = f"{tok.start:08X}"
        if tok.token_type.startswith("PROMOTER"):
            lines.append(f"{addr}:  ENTRY_CALL    {tok.token_type:<18} [{tok.sequence}] ; {tok.metadata.get('role', '')}")
        elif tok.token_type == "RBS_SHINE_DALGARNO":
            lines.append(f"{addr}:  STACK_INIT    {tok.token_type:<18} [{tok.sequence}] ; Ribosome translation setup")
        elif tok.token_type == "START_CODON":
            lines.append(f"{addr}:  EXEC_START    ALLOC_MET          [{tok.sequence}] ; Begin polypeptide execution")
        elif tok.token_type == "STOP_CODON":
            lines.append(f"{addr}:  RETURN_HALT   HALT_OPCODE        [{tok.sequence}] ; Release chain and disassemble complex")
        elif tok.token_type == "OPEN_READING_FRAME":
            orf_seq = tok.sequence
            amino_acids = [GENETIC_CODE.get(orf_seq[c : c + 3], "X") for c in range(0, min(len(orf_seq), 30), 3)]
            sample_aa = "-".join(amino_acids)
            lines.append(f"{addr}:  CODE_BLOCK    STREAM_{tok.metadata.get('codons')}d       [{len(tok.sequence)} bp] ; AA stream: {sample_aa}...")

    return lines


def decompile_to_pseudocode(routine_name: str, tokens: list[BiologicalToken]) -> str:
    """Decompile parsed biological tokens into high-level Python pseudocode."""
    code_lines = [
        f"def {routine_name}(runtime_environment):",
        '    """Decompiled Biological Routine (Operon / Transcript Architecture)."""',
    ]

    has_promoter = any("PROMOTER" in t.token_type for t in tokens)
    has_rbs = any(t.token_type == "RBS_SHINE_DALGARNO" for t in tokens)
    orfs = [t for t in tokens if t.token_type == "OPEN_READING_FRAME"]

    if has_promoter:
        prom = next(t for t in tokens if "PROMOTER" in t.token_type)
        code_lines.append(f"    # Function Entry Address @ offset {prom.start} ({prom.sequence})")
        code_lines.append("    if not runtime_environment.repressor_bound:")

    indent = "        " if has_promoter else "    "

    if has_rbs:
        rbs = next(t for t in tokens if t.token_type == "RBS_SHINE_DALGARNO")
        code_lines.append(f"{indent}# Bind translation hardware: RBS sequence {rbs.sequence} @ offset {rbs.start}")
        code_lines.append(f"{indent}hardware.allocate_ribosome_frame()")

    if orfs:
        for idx, orf in enumerate(orfs, start=1):
            codons_count = orf.metadata.get("codons", len(orf.sequence) // 3)
            code_lines.append(f"{indent}# Subroutine {idx}: {orf.label} (Length: {codons_count} instructions)")
            code_lines.append(f"{indent}yield execute_polypeptide_{idx}(offset={orf.start}, instructions={codons_count})")
    else:
        code_lines.append(f"{indent}# Non-coding functional transcript block")
        code_lines.append(f"{indent}execute_regulatory_rna()")

    code_lines.append(f"{indent}return halt(status=SUCCESS)")

    if has_promoter:
        code_lines.append("    else:")
        code_lines.append("        return halt(status=REPRESSED_GATE_LOCKED)")

    return "\n".join(code_lines)


def disassemble_sequence(
    sequence: str,
    routine_name: str = "biological_routine",
) -> DisassembledRoutine:
    """Disassemble raw biological sequence into tokens, assembly, and pseudocode."""
    tokens = scan_tokens(sequence)
    assembly = generate_assembly_listing(sequence, tokens)
    pseudocode = decompile_to_pseudocode(routine_name, tokens)

    return DisassembledRoutine(
        routine_id=routine_name,
        start_offset=0,
        end_offset=len(sequence),
        length_bp=len(sequence),
        tokens=tokens,
        assembly_listing=assembly,
        decompiled_pseudocode=pseudocode,
    )


def main() -> None:
    """CLI entrypoint for biological disassembler."""
    parser = argparse.ArgumentParser(
        description="Bio-Arch Biological Disassembler & Decompiler",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("file", type=str, help="Path to input FASTA file")
    parser.add_argument("--name", type=str, default="operon_routine", help="Routine name")

    args = parser.parse_args()

    seq_lines = []
    with open(args.file, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith(">"):
                seq_lines.append(line.strip())
    full_seq = "".join(seq_lines)

    res = disassemble_sequence(full_seq, routine_name=args.name)

    print("\n--- ASSEMBLY LISTING ---")
    for line in res.assembly_listing:
        print(line)

    print("\n--- DECOMPILED PSEUDOCODE ---")
    print(res.decompiled_pseudocode)
    print("-----------------------------\n")


if __name__ == "__main__":
    main()
