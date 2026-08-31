# 🚀 Release Notes: HexaPhase Genomics v1.0.0
**Multi-Phase Sequence Analysis, Information Theory & Synthetic Dual-Coding Optimization**

**Release Version:** `v1.0.0`  
**Date:** 2026-08-31  
**DOI:** [10.5281/zenodo.22147682](https://doi.org/10.5281/zenodo.22147682)  
**License:** CC-BY-4.0 / MIT  
**Automated Tests:** ✅ **105 / 105 Unit & Integration Tests Passing**

---

## 🌟 Highlights & Major Capabilities

### 1. 🧬 Trellis Dynamic Programming Dual-Coding Sequence Recompiler
* Globally compiles two independent peptide sequences into a single overlapping DNA construct across Frame 0 and Frame +1.
* **BLOSUM62 Conservation Guidance:** Evaluates biochemical preservation of charge, polarity, and hydropathy when exact dual-codon intersections do not exist.
* **Mammalian Codon Adaptation (CAI):** Uses human relative codon frequencies to optimize expression efficiency in host cells.
* **Cloning Site Avoidance:** Automatically suppresses illegal restriction recognition sequences (`EcoRI`, `BamHI`, `NotI`, `HindIII`, `BsaI`) and homopolymer runs.

### 2. 📦 AAV Vector Capsid Packaging Optimizer
* Evaluates construct footprints against standard **4.7 kb Adeno-Associated Virus (AAV)** packaging limits.
* Real-world showcase demonstrates a **45.8% physical footprint reduction** for GFP (238 aa) and mCherry (236 aa), saving **1,563 bp** and expanding capsid headroom to **2,851 bp**.

### 3. 💾 DNA Digital Storage Multi-Frame Parity Codec
* Encodes arbitrary binary data into Frame 0 with longitudinal parity channels interleaved in Frame +1 for error detection and sequence integrity verification.

### 4. 📤 Industry-Standard Synthetic Biology Vector Exporters
* Annotated GenBank (`.gbk`) with CDS features for both Frame 0 (`/codon_start=1`) and Frame +1 (`/codon_start=2`).
* Synthetic Biology Open Language (SBOL v3) XML/RDF component definitions.
* Multi-FASTA and detailed residue-by-residue CSV alignment tables.

### 5. 🛡️ Epistemological Rigor & Null-Model Calibration
* Strict 4-tier epistemic standard (`[Measurement]`, `[Simulation]`, `[Interpretation]`, `[Hypothesis]`).
* Standardized dinucleotide-preserving shuffles (Altschul-Erickson) verifying that sequence patterns are evaluated against empirical null baselines.

---

## 🚀 Quickstart Installation
```bash
pip install bio-arch
```
or install from source:
```bash
git clone https://github.com/zekvftb/hexaphase-genomics.git
cd hexaphase-genomics
pip install -e .
python -m pytest tests/ -v
```
