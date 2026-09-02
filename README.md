# 🧬 HexaPhase Genomics: Multi-Phase Sequence Analysis, Information Theory & Synthetic Optimization

[![DOI](https://zenodo.org/badge/1349560488.svg)](https://doi.org/10.5281/zenodo.22147682)
[![License: CC-BY-4.0](https://img.shields.io/badge/License-CC--BY--4.0-blue.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Tests: 121 Passed](https://img.shields.io/badge/tests-121%20passed-brightgreen.svg)](#reproducibility--testing)
[![Evidence: 4-Tier Categorized](https://img.shields.io/badge/Evidence-4--Tier%20Standard-blue.svg)](docs/scientific_limitations.md)

An open-source computational biology platform providing descriptive information-theoretic sequence analysis, structural regulatory element scanning, and combinatorial optimization for compact synthetic DNA constructs.

---

## 📖 Manuscripts & Working Papers

Full-text working manuscripts and journal submission packages are compiled in Markdown and PDF format under [`manuscripts/`](manuscripts/):

### 📄 1. Information Density & Multi-Phase Analysis
* 📜 **Manuscript:** [`HexaPhase_Universal_Multiplexed_Subroutines_Manuscript.md`](manuscripts/HexaPhase_Universal_Multiplexed_Subroutines_Manuscript.md) | 📥 **PDF:** [`HexaPhase_Universal_Multiplexed_Subroutines_Manuscript.pdf`](manuscripts/HexaPhase_Universal_Multiplexed_Subroutines_Manuscript.pdf)
* ✉️ **Cover Letter:** [`HexaPhase_Universal_Multiplexed_Subroutines_Cover_Letter.md`](manuscripts/HexaPhase_Universal_Multiplexed_Subroutines_Cover_Letter.md)
* *Preprint:* bioRxiv / Zenodo DOI: `10.5281/zenodo.22147682`

### 📄 2. Mitochondrial D-Loop Open Reading Frames & Evolutionary Conservation
* 📜 **Manuscript:** [`Mitochondrial_DLoop_Signaling_Peptide_Oncology_Manuscript.md`](manuscripts/Mitochondrial_DLoop_Signaling_Peptide_Oncology_Manuscript.md) | 📥 **PDF:** [`Mitochondrial_DLoop_Signaling_Peptide_Oncology_Manuscript.pdf`](manuscripts/Mitochondrial_DLoop_Signaling_Peptide_Oncology_Manuscript.pdf)
* ✉️ **Cover Letter:** [`Mitochondrial_DLoop_Signaling_Peptide_Oncology_Cover_Letter.md`](manuscripts/Mitochondrial_DLoop_Signaling_Peptide_Oncology_Cover_Letter.md)

### 📄 3. Coronavirus Viroporin Conserved Motifs & Sequence Analysis
* 📜 **Manuscript:** [`Pan_Coronavirus_Frozen_Viroporin_Antiviral_Manuscript.md`](manuscripts/Pan_Coronavirus_Frozen_Viroporin_Antiviral_Manuscript.md) | 📥 **PDF:** [`Pan_Coronavirus_Frozen_Viroporin_Antiviral_Manuscript.pdf`](manuscripts/Pan_Coronavirus_Frozen_Viroporin_Antiviral_Manuscript.pdf)
* ✉️ **Cover Letter:** [`Pan_Coronavirus_Frozen_Viroporin_Antiviral_Cover_Letter.md`](manuscripts/Pan_Coronavirus_Frozen_Viroporin_Antiviral_Cover_Letter.md)

### 📄 4. Descriptive Symbolic Intermediate Representations for Genomes
* 📜 **Manuscript:** [`Instruction_Set_Architecture_Living_Genomes_Manuscript.md`](manuscripts/Instruction_Set_Architecture_Living_Genomes_Manuscript.md) | 📥 **PDF:** [`Instruction_Set_Architecture_Living_Genomes_Manuscript.pdf`](manuscripts/Instruction_Set_Architecture_Living_Genomes_Manuscript.pdf)
* ✉️ **Cover Letter:** [`Instruction_Set_Architecture_Living_Genomes_Cover_Letter.md`](manuscripts/Instruction_Set_Architecture_Living_Genomes_Cover_Letter.md)

### 📑 Project Roadmap & Documentation:
* 🗺️ **[`PROJECT_DISCOVERY_JOURNAL.md`](PROJECT_DISCOVERY_JOURNAL.md)** *(Chronological Project Journal & Evidence Classification)*
* ⚠️ **[`docs/scientific_limitations.md`](docs/scientific_limitations.md)** *(Scientific Limitations, Epistemological Boundaries & Null Controls)*
* 📋 **[`docs/reproducibility.md`](docs/reproducibility.md)** *(Reproducibility Protocol & Environment Metadata)*

---

## ⚡ Computational Architecture & Modules

The platform is structured into modular Python packages with typed data contracts, deterministic random seeds, and composition-preserving null models:

| Module | Core Functionality | Evidence Classification |
| :--- | :--- | :--- |
| **Module 0: Ingestion** | Schema validation, checksum verification, and case normalization | `[Measurement]` |
| **Module 1: Information** | Shannon entropy, autocorrelation, and dinucleotide-shuffled null models | `[Measurement]` |
| **Module 2: Regulation** | Feed-forward loops, feedback cycles, and degree-preserving randomized graphs | `[Measurement]` / `[Simulation]` |
| **Module 3: Simulation** | Discrete Boolean network attractor dynamics and basin stability | `[Simulation]` |
| **Module 4: Structural Features** | Programmed frameshift motifs, G4 quadruplexes, and readthrough contexts | `[Measurement]` / `[Interpretation]` |
| **Module 5: Symbolic IR** | Codon position entropy and symbolic assembly disassembly (`.asm` IR) | `[Measurement]` |
| **Module 6: Recompiler** | Trellis DP sequence optimization & BLOSUM62 conservation for dual-coding AAV constructs | `[Simulation]` / Optimization |
| **Module 7: Parity Verification** | Implementation cross-verification between Python reference and Bytecode VM | `[Measurement]` |
| **Module 8: Feature Extraction** | CRISPR direct-repeat arrays, Gardiner-Garden CpG islands, and riboswitches | `[Measurement]` / `[Interpretation]` |
| **Module 9: Quantum Biophysics** | 1D WKB proton tunneling barrier simulations and radical pair dynamics | `[Simulation]` / `[Hypothesis]` |

---

## 🧬 Real-World Showcase: AAV Vector Compaction (GFP + mCherry)

Dual-coding sequence compilation interleaves two distinct protein sequences into a single overlapping DNA construct across Frame 0 and Frame +1, dramatically increasing payload headroom for size-constrained viral vectors such as **Adeno-Associated Virus (AAV, 4.7 kb capacity limit)**:

| Metric | Separate Dual Monomers | HexaPhase Dual-Phase Compaction | Practical Impact |
| :--- | :---: | :---: | :--- |
| **Target Reporter Proteins** | GFP (238 aa) + mCherry (236 aa) | GFP (F0) + mCherry (F1) | Dual-reporter construct |
| **Promoters & PolyA Signals** | 2 Promoters + 2 PolyA | 1 Single Overlapping Cassette | Saves 850 bp non-coding overhead |
| **Total Cassette Footprint** | **3,412 bp** | **1,849 bp** | **45.8% Physical Size Reduction** |
| **AAV 4.7 kb Capsid Headroom** | 1,288 bp remaining | **2,851 bp remaining** | **+1,563 bp Extra Payload Space** |
| **Human Codon Adaptation (CAI)** | 0.7500 | **0.7687** | High mammalian expression efficiency |
| **BLOSUM62 Conservation** | N/A | **65.68%** | Preserves charge & hydrophobicity |
| **Cloning Exclusions** | Standard | Zero `EcoRI`, `BamHI`, `NotI`, `BsaI` | Ready for Golden Gate assembly |

### 🚀 Python SDK Quickstart: Dual-Protein Compaction & GenBank Export
```python
from bio_arch.modules.recompiler import recompile_dual_protein_dna, calculate_aav_packaging_savings
from bio_arch.modules.export_formats import export_to_genbank, export_to_sbol3

# 1. Define target amino acid sequences
gfp = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYG..."
mcherry = "MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEG..."

# 2. Synthesize dual-coding DNA with mammalian CAI and BLOSUM62 relaxation
res = recompile_dual_protein_dna(gfp, mcherry, optimize_cai=True, allow_conservative_mutations=True)

# 3. Calculate AAV packaging headroom savings
savings = calculate_aav_packaging_savings(len(gfp), len(mcherry), res.total_length_bp)
print(f"Saved {savings['bp_saved']} bp ({savings['percent_footprint_reduction']}%) in AAV vector!")

# 4. Export annotated GenBank and SBOL v3 files
gbk_content = export_to_genbank(res, locus_name="GFP_MCHERRY_AAV")
with open("vector_construct.gbk", "w") as f:
    f.write(gbk_content)
```

---

## 🔬 Scientific Data & Audit Outputs

* 📊 [`outputs/master_cohort_integrity_ledger.json`](outputs/master_cohort_integrity_ledger.json) — 54-sequence cohort registry with cryptographic SHA-256 digests.
* 📊 [`outputs/dual_engine_verification_audit.json`](outputs/dual_engine_verification_audit.json) — Python vs SMC VM execution parity benchmark.
* 📊 [`outputs/biological_circuits_audit.json`](outputs/biological_circuits_audit.json) — CRISPR arrays, CpG islands, and riboswitch census.
* 📊 [`outputs/quantum_biology_audit.json`](outputs/quantum_biology_audit.json) — Theoretical radical pair and tunneling simulation outputs.
* 📊 [`outputs/synbio_and_neurogenomics_audit.json`](outputs/synbio_and_neurogenomics_audit.json) — Dual-coding synthetic construct optimization data.

---

## 🔍 Reproducibility & Testing

All analysis runs, null models, and tests can be executed locally:

```bash
# 1. Run the complete automated test suite (81 tests)
python -m pytest tests/ -v

# 2. Run the Universal Master Cohort Integrity Audit
python scripts/run_universal_cohort_audit.py

# 3. Run the Dual-Engine Cross-Verification Benchmark
python scripts/run_dual_engine_audit.py

# 4. Run the Structural Feature & Circuit Scanner
python scripts/scan_biological_circuits.py
```

---

## 📑 Citation

```bibtex
@software{Rezek_HexaPhase_Genomics_2026,
  author       = {Rezek, Jason},
  title        = {{HexaPhase Genomics: Multi-Phase Sequence Analysis, Information Theory & Synthetic Optimization}},
  month        = aug,
  year         = 2026,
  publisher    = {GitHub / Zenodo},
  version      = {1.0.0},
  doi          = {10.5281/zenodo.22147682},
  url          = {https://github.com/zekvftb/hexaphase-genomics}
}
```
