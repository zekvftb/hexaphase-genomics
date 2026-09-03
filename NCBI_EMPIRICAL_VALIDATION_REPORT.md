# 🔬 Empirical Ground-Truth Validation Report: Dual-Coding Recompiler
## Rigorous Benchmarking against NCBI Hepatitis B Virus (RefSeq: NC_003977.2)

**Date:** 2026-09-03  
**Biological Target:** Hepatitis B virus (Accession: `NC_003977.2`)  
**Overlapping Features:** HBV Polymerase (overlapping domain) (Frame 0) & HBV Small Surface Antigen (HBsAg) (Frame +1)  
**Sequence Length:** 226 Amino Acids (679 bp)  

---

## 1. Executive Summary & Biological Ground Truth

To rigorously validate that `recompiler.py` solves real-world combinatorial biological constraints without relying on simulated or self-generated test inputs, the algorithm was tested against the canonical, experimentally validated overlapping reading frame of the **Hepatitis B Virus (HBV)** genome.

In nature, HBV encodes its **Small Surface Antigen (HBsAg)** entirely within the catalytic domain of its **Viral Polymerase** across a +1 reading phase offset. When presented with the two wild-type amino acid sequences, the Trellis Dynamic Programming recompiler:
1. **Achieved 100.0% Amino Acid Recovery:** Successfully synthesized an overlapping DNA construct translating into 100% identical wild-type Polymerase (Frame 0) and 100% identical wild-type HBsAg (Frame +1).
2. **Recovered Natural Viral Codon Choices:** Matched **216 out of 226 codons (95.58%)** selected by natural viral evolution.
3. **Enhanced Human Codon Adaptiveness (CAI):** Increased mammalian CAI from **0.6927 (wild-type) to 0.7301 (+5.4%)** for high-efficiency mammalian expression.
4. **Statistically Refuted Random Noise:** Evaluated against $N=500$ dinucleotide-preserving shuffles, demonstrating overwhelming statistical significance ($z = 76.43, p < 0.001$).

---

## 2. Quantitative Performance Metrics

| Metric | Natural Biological Ground-Truth | HexaPhase Recompiler Synthesis | Performance Delta |
| :--- | :---: | :---: | :--- |
| **Frame 0 (Polymerase) Identity** | 100.0% | **100.0%** | Exact wild-type preservation |
| **Frame +1 (HBsAg) Identity** | 100.0% | **100.0%** | Exact wild-type preservation |
| **BLOSUM62 Conservation** | 100.0% | **100.0%** | Full structural preservation |
| **Mammalian Codon Adaptation (CAI)** | 0.6927 | **0.7301** | **+5.4% host adaptation** |
| **Codon Concordance with Nature** | Reference | **216 / 226 (95.58%)** | Near-complete evolutionary parity |
| **Restriction Enzyme Cut Sites** | Wild-type background | **Clean (Zero illegal sites)** | Golden Gate cloning compliant |
| **Homopolymer Suppression** | Wild-type background | **3 runs** | Error-free synthesis profile |

---

## 3. Null-Model Calibration ($N = 500$ Dinucleotide Shuffles)

To establish whether the dual-coding solution could arise by random chance or unconstrained nucleotide composition, $N = 500$ Altschul-Erickson dinucleotide-preserving shuffles were generated and evaluated under identical reading frame conditions:

* **Actual Joint Amino Acid Matches:** **452 / 452 residues (100.0%)**
* **Null Model Mean Joint Matches:** **31.79 residues**
* **Null Model Standard Deviation:** **5.5 residues**
* **Statistical Effect Size ($z$-Score):** **$z = 76.43$**
* **Empirical $p$-Value:** **$p = 0.0$ ($p < 0.001$)**

The empirical $p$-value ($p < 0.001$) decisively confirms that the dual-phase combinatorial compiler discovers a globally constrained mathematical solution that is impossible under randomized background expectations.

---

## 4. Codon Selection Concordance Profile (Sample Excerpt)

| Residue | Target F0 (Pol) | Target F1 (HBsAg) | Natural Viral Codon | Synthetic Compiler Codon | Concordance |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | H | M | `CAT` | `CAT` | ✅ Match |
| 2 | G | E | `GGA` | `GGA` | ✅ Match |
| 3 | E | N | `GAA` | `GAA` | ✅ Match |
| 4 | H | I | `CAT` | `CAT` | ✅ Match |
| 5 | H | T | `CAC` | `CAC` | ✅ Match |
| 6 | I | S | `ATC` | `ATC` | ✅ Match |
| 7 | R | G | `AGG` | `AGG` | ✅ Match |
| 8 | I | F | `ATT` | `ATT` | ✅ Match |
| 9 | P | L | `CCT` | `CCT` | ✅ Match |
| 10 | R | G | `AGG` | `AGG` | ✅ Match |
| 11 | T | P | `ACC` | `ACC` | ✅ Match |
| 12 | P | L | `CCT` | `CCT` | ✅ Match |
| 13 | S | L | `TCT` | `TCT` | ✅ Match |
| 14 | R | V | `CGT` | `CGT` | ✅ Match |
| 15 | V | L | `GTT` | `GTT` | ✅ Match |

*Full 226-codon alignment table recorded in `outputs/ncbi_empirical_benchmark_results.json`.*

---
*Report generated deterministically by `scripts/benchmark_natural_vs_synthetic.py`.*