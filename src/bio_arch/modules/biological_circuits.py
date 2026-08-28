"""Module 8: Biological Memory & Security Architecture.

Decompiles advanced cellular computing primitives:
1. CRISPR Direct-Repeat & Spacer Arrays (Append-Only Hardware Malware Logs)
2. CpG Epigenetic Memory Islands (Non-Volatile 1-Bit Latches / NVRAM)
3. RNA Riboswitch Transducers (Analog-to-Digital Chemical Converters / ADCs)
"""

from __future__ import annotations

from collections import defaultdict
import math
import re
from typing import Any

from bio_arch.contracts import (
    BiologicalCircuitReport,
    CpgMemoryIsland,
    CrisprArray,
    RiboswitchAdc,
)
from bio_arch.logger import setup_logger
from bio_arch.modules.logic_gates import calculate_stem_loop_mfe

logger = setup_logger("bio_arch.biological_circuits")

# Canonical Riboswitch Aptamer Core Signatures (Rfam consensus approximations)
RIBOSWITCH_SIGNATURES = [
    (
        re.compile(r"ACTCC[ACGT]{10,25}GGAGT", re.IGNORECASE),
        "PURINE_GUANINE",
        "Purine / Guanine Metabolite Sensing Aptamer",
    ),
    (
        re.compile(r"GAGTC[ACGT]{8,20}GACTC", re.IGNORECASE),
        "SAM_I",
        "S-Adenosylmethionine (SAM-I) Methylation Sensor",
    ),
    (
        re.compile(r"AGGAA[ACGT]{10,25}TTCCT", re.IGNORECASE),
        "FMN_FLAVIN",
        "Flavin Mononucleotide (FMN) Redox State Sensor",
    ),
    (
        re.compile(r"CTGAGA[ACGT]{8,22}TCTCAG", re.IGNORECASE),
        "TPP_THIAMINE",
        "Thiamine Pyrophosphate (TPP) Vitamin B1 Sensor",
    ),
]


def scan_cpg_islands(
    sequence: str,
    window_size: int = 200,
    step: int = 50,
    min_gc: float = 0.50,
    min_obs_exp: float = 0.60,
) -> list[CpgMemoryIsland]:
    """Scan sequence for epigenetic 1-bit non-volatile memory registers (Gardiner-Garden & Frommer criteria).

    Criteria:
    - Window length >= 200 bp
    - GC content >= 50%
    - Observed-to-Expected CpG ratio >= 0.60
      Obs/Exp = (Count(CG) * Length) / (Count(C) * Count(G))
    """
    seq_upper = sequence.upper()
    n = len(seq_upper)
    if n < window_size:
        return []

    raw_islands: list[tuple[int, int, float, float, int]] = []

    for i in range(0, n - window_size + 1, step):
        window = seq_upper[i : i + window_size]
        c_count = window.count("C")
        g_count = window.count("G")
        cg_count = window.count("CG")
        gc_content = (c_count + g_count) / window_size

        if c_count == 0 or g_count == 0:
            obs_exp = 0.0
        else:
            obs_exp = (cg_count * window_size) / (c_count * g_count)

        if gc_content >= min_gc and obs_exp >= min_obs_exp:
            raw_islands.append((i + 1, i + window_size, gc_content, obs_exp, cg_count))

    if not raw_islands:
        return []

    # Merge overlapping qualifying windows
    merged: list[CpgMemoryIsland] = []
    curr_start, curr_end, curr_gc, curr_oe, curr_cg = raw_islands[0]
    island_idx = 1

    for start, end, gc, oe, cg in raw_islands[1:]:
        if start <= curr_end:
            # Overlap: extend
            curr_end = max(curr_end, end)
            curr_gc = max(curr_gc, gc)
            curr_oe = max(curr_oe, oe)
            curr_cg += cg
        else:
            length = curr_end - curr_start + 1
            merged.append(
                CpgMemoryIsland(
                    island_id=f"CPG_MEM_{island_idx:03d}",
                    start_pos=curr_start,
                    end_pos=curr_end,
                    length_bp=length,
                    gc_content_pct=round(curr_gc * 100.0, 2),
                    cpg_obs_exp_ratio=round(curr_oe, 3),
                    cpg_count=curr_cg,
                    is_promoter_latch=(curr_oe >= 0.65),
                )
            )
            island_idx += 1
            curr_start, curr_end, curr_gc, curr_oe, curr_cg = start, end, gc, oe, cg

    # Append final merged island
    length = curr_end - curr_start + 1
    merged.append(
        CpgMemoryIsland(
            island_id=f"CPG_MEM_{island_idx:03d}",
            start_pos=curr_start,
            end_pos=curr_end,
            length_bp=length,
            gc_content_pct=round(curr_gc * 100.0, 2),
            cpg_obs_exp_ratio=round(curr_oe, 3),
            cpg_count=curr_cg,
            is_promoter_latch=(curr_oe >= 0.65),
        )
    )

    return merged


def scan_crispr_arrays(
    sequence: str,
    min_repeat_len: int = 24,
    max_repeat_len: int = 40,
    min_spacer_len: int = 20,
    max_spacer_len: int = 75,
    min_repeats: int = 3,
) -> list[CrisprArray]:
    """Identify and decompile CRISPR direct-repeat / spacer arrays (Append-Only Malware Logs)."""
    seq_upper = sequence.upper()
    n = len(seq_upper)
    if n < 100:
        return []

    arrays: list[CrisprArray] = []
    kmer_positions: dict[str, list[int]] = defaultdict(list)

    # Index candidates of repeat size k
    k = min_repeat_len
    for i in range(0, n - k + 1):
        kmer = seq_upper[i : i + k]
        kmer_positions[kmer].append(i)

    visited_starts: set[int] = set()

    for kmer, positions in kmer_positions.items():
        if len(positions) < min_repeats:
            continue

        # Check if positions have regular periodic spacer spacing (20-75 bp)
        valid_chain = [positions[0]]
        for pos in positions[1:]:
            prev = valid_chain[-1]
            dist = pos - (prev + k)
            if min_spacer_len <= dist <= max_spacer_len:
                valid_chain.append(pos)
            elif dist > max_spacer_len:
                if len(valid_chain) >= min_repeats:
                    break
                valid_chain = [pos]

        if len(valid_chain) >= min_repeats:
            start_idx = valid_chain[0]
            if start_idx in visited_starts:
                continue
            visited_starts.add(start_idx)

            # Extract spacers
            spacers: list[str] = []
            for idx in range(len(valid_chain) - 1):
                s_start = valid_chain[idx] + k
                s_end = valid_chain[idx + 1]
                spacers.append(seq_upper[s_start:s_end])

            end_idx = valid_chain[-1] + k
            array_id = f"CRISPR_ARRAY_{start_idx+1}_{end_idx}"

            arrays.append(
                CrisprArray(
                    array_id=array_id,
                    start_pos=start_idx + 1,
                    end_pos=end_idx,
                    repeat_length_bp=k,
                    repeat_consensus=kmer,
                    repeats_count=len(valid_chain),
                    spacers_count=len(spacers),
                    spacers=spacers,
                )
            )

    return arrays


def scan_riboswitch_adcs(sequence: str, min_switching_mfe: float = -6.0) -> list[RiboswitchAdc]:
    """Scan sequence for RNA Riboswitch Analog-to-Digital chemical sensors & switches."""
    seq_upper = sequence.upper()
    n = len(seq_upper)
    switches: list[RiboswitchAdc] = []

    switch_idx = 1
    for pat, ligand_class, desc in RIBOSWITCH_SIGNATURES:
        for match in pat.finditer(seq_upper):
            start = match.start()
            end = match.end()
            motif = match.group(0)

            # Check downstream regulatory expression platform (30-60 nt downstream)
            downstream_start = end
            downstream_end = min(n, end + 50)
            if downstream_end - downstream_start < 15:
                continue

            downstream_seq = seq_upper[downstream_start:downstream_end]
            dG, stem_len, loop_len = calculate_stem_loop_mfe(downstream_seq)

            if dG <= min_switching_mfe:
                switches.append(
                    RiboswitchAdc(
                        switch_id=f"RIBOSWITCH_ADC_{switch_idx:03d}",
                        start_pos=start + 1,
                        end_pos=downstream_end,
                        ligand_class=ligand_class,
                        aptamer_motif=motif,
                        terminator_mfe_dG=dG,
                        switching_delta_dG=round(abs(dG) * 0.65, 2),
                        predicted_state="OFF_IN_HIGH_LIGAND",
                        description=f"{desc} (Downstream Hairpin dG: {dG} kcal/mol)",
                    )
                )
                switch_idx += 1

    return switches


def scan_all_biological_circuits(
    sequence: str,
    genome_id: str = "unknown",
) -> BiologicalCircuitReport:
    """Execute complete biological memory & security architecture audit across a genome."""
    crispr = scan_crispr_arrays(sequence)
    cpg = scan_cpg_islands(sequence)
    riboswitches = scan_riboswitch_adcs(sequence)

    summary = {
        "crispr_security_arrays_found": len(crispr),
        "total_spacers_archived": sum(a.spacers_count for a in crispr),
        "cpg_memory_islands_mapped": len(cpg),
        "promoter_memory_latches": sum(1 for c in cpg if c.is_promoter_latch),
        "riboswitch_adcs_detected": len(riboswitches),
    }

    return BiologicalCircuitReport(
        genome_id=genome_id,
        genome_length_bp=len(sequence),
        crispr_arrays=crispr,
        cpg_memory_islands=cpg,
        riboswitch_adcs=riboswitches,
        summary=summary,
    )
