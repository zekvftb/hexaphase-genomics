"""CLI Tool: De Novo Dual-Phase DNA Recompiler."""

from __future__ import annotations

import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from bio_arch.modules.recompiler import recompile_dual_protein_dna
from bio_arch.modules.universal_compiler import decompile_genomic_bytecode

def main():
    print("=======================================================================================================")
    print("🧬 HexaPhase De Novo Dual-Protein Recompiler Engine")
    print("=======================================================================================================\n")

    # Example 1: Synthesize a therapeutic signaling peptide + antiviral viroporin in ONE sequence
    p0_therapeutic = "MSQYLSLIPASSYYLSHLRSMLQANMLTKVC" # 31 aa (Mitochondrial D-loop signaling peptide)
    p1_viroporin   = "MLLTLLCTLLLVIYYMLLTLLCTLLLVIYYA" # 31 aa (Viroporin transmembrane core)

    print(f"Target Protein Frame 0 (Therapeutic Peptide): {p0_therapeutic}")
    print(f"Target Protein Frame +1 (Antiviral Viroporin): {p1_viroporin}\n")

    res = recompile_dual_protein_dna(p0_therapeutic, p1_viroporin)

    print("--- Recompilation & Synthesis Results ---")
    print(f"  * Synthesized Dual-Phase DNA ({res.total_length_bp} bp):")
    print(f"    {res.synthesized_dna}")
    print(f"  * Frame 0 Translation: {res.translated_f0} (Fidelity: {res.f0_identity_pct}%)")
    print(f"  * Frame +1 Translation: {res.translated_f1} (Fidelity: {res.f1_identity_pct}%)")
    print(f"  * Physical Information Compression Ratio: {res.compression_ratio}x (Saved {(31*6) - res.total_length_bp} bp)")
    print(f"  * Wobble Carrier Adjustments Applied: {res.wobble_mutations_applied}\n")

    # Disassemble synthesized sequence
    instructions, asm_text = decompile_genomic_bytecode(res.synthesized_dna, genome_id="SYNTHETIC_DUAL_PHASE_01")
    print(f"✅ Successfully compiled & disassembled into {len(instructions)} biological machine code opcodes!")

if __name__ == "__main__":
    main()
