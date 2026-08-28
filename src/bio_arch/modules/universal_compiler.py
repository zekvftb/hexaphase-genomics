"""Module 5: Universal Biological Compiler & Genomic Machine Code Decompiler.

Decodes the fundamental assembly instruction set architecture (ISA) of biological genomes.
Provides:
1. Wobble-position entropy and dual-phase information synchronization analysis.
2. Genomic Bytecode Decompiler: Converts raw nucleotide sequences into cycle-accurate
   Biological Assembly (.asm) and executable SMC machine code.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Any

from bio_arch.modules.logic_gates import (
    calculate_stem_loop_mfe,
    scan_all_logic_gates,
)

# Standard Genetic Code Table 1 (Universal)
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

COMPLEMENT_MAP = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn", "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")


@dataclass
class WobbleSyncReport:
    """Mathematical analysis of multi-track information transmission via wobble positions."""

    total_codons_analyzed: int
    pos1_shannon_entropy_bits: float
    pos2_shannon_entropy_bits: float
    pos3_wobble_entropy_bits: float
    pos3_information_capacity_ratio: float
    mutual_information_f0_f1_bits: float
    dual_channel_capacity_bits_per_base: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class AssemblyInstruction:
    """A single disassembled biological machine code instruction."""

    address_hex: str
    offset_bp: int
    raw_triplet: str
    opcode: str
    operand_primary: str
    parallel_tracks: dict[str, str] = field(default_factory=dict)
    hardware_event: str | None = None
    comment: str = ""

    def to_asm_line(self) -> str:
        hw_str = f" ; [{self.hardware_event}]" if self.hardware_event else ""
        comm_str = f" -- {self.comment}" if self.comment else ""
        tracks_str = f"|| +1:{self.parallel_tracks.get('+1', '-')} -0:{self.parallel_tracks.get('-0', '-')}"
        return f"{self.address_hex:<8} {self.raw_triplet:<5} {self.opcode:<16} {self.operand_primary:<12} {tracks_str:<20}{hw_str}{comm_str}"


def compute_shannon_entropy(symbols: list[str]) -> float:
    """Compute Shannon entropy H(X) in bits."""
    if not symbols:
        return 0.0
    counts: dict[str, int] = {}
    for s in symbols:
        counts[s] = counts.get(s, 0) + 1
    total = len(symbols)
    entropy = 0.0
    for cnt in counts.values():
        p = cnt / total
        if p > 0:
            entropy -= p * math.log2(p)
    return round(entropy, 4)


def analyze_wobble_synchronization(sequence: str) -> WobbleSyncReport:
    """Analyze how the 3rd wobble position serves as an independent parallel data carrier."""
    clean_seq = sequence.upper().replace("U", "T")
    num_codons = len(clean_seq) // 3
    if num_codons == 0:
        return WobbleSyncReport(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    pos1 = [clean_seq[i * 3] for i in range(num_codons)]
    pos2 = [clean_seq[i * 3 + 1] for i in range(num_codons)]
    pos3 = [clean_seq[i * 3 + 2] for i in range(num_codons)]

    h1 = compute_shannon_entropy(pos1)
    h2 = compute_shannon_entropy(pos2)
    h3 = compute_shannon_entropy(pos3)

    # In dual coding, Pos3 in Frame 0 serves as Pos2 in Frame +1
    # Mutual information estimation I(F0; F+1)
    joint_pairs = [f"{pos3[i]}_{pos1[(i + 1) % num_codons]}" for i in range(num_codons)]
    h_joint = compute_shannon_entropy(joint_pairs)
    mi = max(0.0, round(h3 + h1 - h_joint, 4))

    capacity_ratio = round(h3 / max(0.001, (h1 + h2) / 2.0), 3)
    total_capacity = round((h1 + h2 + h3 + mi) / 3.0, 3)

    return WobbleSyncReport(
        total_codons_analyzed=num_codons,
        pos1_shannon_entropy_bits=h1,
        pos2_shannon_entropy_bits=h2,
        pos3_wobble_entropy_bits=h3,
        pos3_information_capacity_ratio=capacity_ratio,
        mutual_information_f0_f1_bits=mi,
        dual_channel_capacity_bits_per_base=total_capacity,
    )


def decompile_genomic_bytecode(sequence: str, genome_id: str = "GENOME_01") -> tuple[list[AssemblyInstruction], str]:
    """Decompile raw nucleotide DNA into structured Biological Assembly code (.asm)."""
    clean_seq = sequence.upper().replace("U", "T")
    rev_seq = clean_seq.translate(COMPLEMENT_MAP)[::-1]
    n = len(clean_seq)
    
    # Run logic gate audit to detect hardware interrupts
    gate_report = scan_all_logic_gates(clean_seq, genome_id=genome_id)
    gates_by_pos: dict[int, list[Any]] = {}
    for gate in gate_report.gates_found:
        gates_by_pos.setdefault(gate.start_pos, []).append(gate)

    instructions: list[AssemblyInstruction] = []
    
    # Header entry
    instructions.append(
        AssemblyInstruction(
            address_hex="0x0000",
            offset_bp=0,
            raw_triplet="---",
            opcode="SYS_INIT_PROMOTER",
            operand_primary=f"TARGET:{genome_id}",
            parallel_tracks={"+1": "---", "-0": "---"},
            comment="Entry point into biological execution tape",
        )
    )

    for i in range(0, n - 2, 3):
        triplet_f0 = clean_seq[i : i + 3]
        triplet_f1 = clean_seq[i + 1 : i + 4] if i + 4 <= n else clean_seq[i + 1 :] + "N" * (3 - (n - i - 1))
        
        # Reverse complement coordinate for antisense track -0
        rev_idx = max(0, n - (i + 3))
        triplet_rev0 = rev_seq[rev_idx : rev_idx + 3] if rev_idx + 3 <= len(rev_seq) else "NNN"

        aa_f0 = CODON_TABLE.get(triplet_f0, "X")
        aa_f1 = CODON_TABLE.get(triplet_f1, "X")
        aa_rev0 = CODON_TABLE.get(triplet_rev0, "X")

        addr_hex = f"0x{i:04X}"
        opcode = "PUSH_PEPTIDE"
        operand = f"AA:{aa_f0} ({triplet_f0})"
        hw_event = None
        comment = ""

        # Check for start/stop states
        if triplet_f0 == "ATG":
            opcode = "START_SUBROUTINE"
            comment = "Initiation Methionine Codon"
        elif triplet_f0 in ("TAA", "TAG", "TGA"):
            opcode = "TERM_STACK_POP"
            comment = "Termination Stop Codon"

        # Check for hardware interrupts at this coordinate
        if i in gates_by_pos or (i - 1) in gates_by_pos or (i - 2) in gates_by_pos:
            active_gates = gates_by_pos.get(i, []) + gates_by_pos.get(i - 1, []) + gates_by_pos.get(i - 2, [])
            for g in active_gates:
                if g.gate_type.value == "frameshift_branch":
                    opcode = "BRANCH_SLIP_MUX"
                    hw_event = f"SLIP -1 ({int(g.predicted_efficiency * 100)}% Divert)"
                    comment = f"Hairpin barrier dG: {g.downstream_barrier_energy} kcal/mol"
                elif g.gate_type.value == "g4_circuit_breaker":
                    opcode = "CIRCUIT_LATCH_G4"
                    hw_event = "G4_TRANSISTOR_TRIP"
                    comment = f"G4 Stack Hunter Score: {g.metrics.get('g4_score', 1.0)}"
                elif g.gate_type.value == "readthrough_overflow":
                    opcode = "OVERFLOW_BYPASS"
                    hw_event = f"LEAKY_STOP ({int(g.predicted_efficiency * 100)}% Readthrough)"

        inst = AssemblyInstruction(
            address_hex=addr_hex,
            offset_bp=i,
            raw_triplet=triplet_f0,
            opcode=opcode,
            operand_primary=operand,
            parallel_tracks={"+1": f"{aa_f1} ({triplet_f1})", "-0": f"{aa_rev0} ({triplet_rev0})"},
            hardware_event=hw_event,
            comment=comment,
        )
        instructions.append(inst)

    # Format into complete Assembly listing
    asm_lines = [
        "; ===========================================================================================",
        f"; BIOLOGICAL MACHINE CODE DISASSEMBLY LISTING: {genome_id}",
        f"; Total Genome Length: {n:,} base pairs | Total Assembly Opcodes: {len(instructions):,}",
        "; Architecture: HexaPhase 6-Track Concurrent Molecular Virtual Machine",
        "; ===========================================================================================",
        f"{'OFFSET':<8} {'RAW':<5} {'OPCODE':<16} {'PRIMARY TRACK':<12} {'PARALLEL TRACKS':<20} {'INTERRUPTS / COMMENTS'}",
        "-" * 105,
    ]
    for inst in instructions:
        asm_lines.append(inst.to_asm_line())

    asm_lines.append("-" * 105)
    asm_lines.append(f"; [HALT] End of biological execution tape for {genome_id}.\n")
    
    return instructions, "\n".join(asm_lines)
