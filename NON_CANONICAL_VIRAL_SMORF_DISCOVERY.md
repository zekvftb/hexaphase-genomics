# 🔬 Non-Canonical smORF Discovery & Translational Anomaly Report
## Systematic Mining of Near-Cognate Initiation Codons & Viroporin Motifs Across Viral Genomes

**Date:** 2026-09-02  
**Genomes Screened:** 9 RefSeq Viral Genomes  
**smORF Window:** 30 to 120 Amino Acids  
**Host CAI Threshold:** $\ge 0.70$  
**Null Controls:** $N=200$ Altschul-Erickson Dinucleotide Shuffles per Candidate  
**Falsification Standard:** Significant only if $z > 3.0$ and $p < 0.001$  

---

## 1. Executive Discovery Summary

Standard genome annotation pipelines rely strictly on canonical AUG start codons, systematically missing functional microproteins (smORFs) initiated by **near-cognate start codons (CUG, GUG, ACG, AUA, UUG)**. Using Kozak initiation context scoring, host codon adaptation modeling, and $N=500$ dinucleotide null model calibrations, this engine systematically screened the RefSeq viral corpus.

### Key Findings:
- **Total Candidates Evaluated (CAI $\ge 0.70$, 30–120 aa):** 142
- **Statistically Significant ($z > 3.0, p < 0.001$):** 1
- **Novel Unannotated Candidates:** **0**

---

## 2. Top Non-Canonical smORF Candidates ($z > 3.0, p < 0.001$)

| Organism (Accession) | Parent CDS | Frame | Start Codon (Type) | Kozak | Length | Host CAI | TM Helices | $\mu_H$ | $z$-Score | Viroporin Profile |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| Escherichia phage ph (`NC_001422.1`) | `CDS` | +2 | `ATG` (CANONICAL_AUG) | Optimal | **120 aa** | **0.7019** | 0 TM | 0.044 | **$z = 4.24$** | Soluble / Globular |

---

## 3. High-Priority Viroporin & Transmembrane smORF Profiles

---

## 4. Methodological Controls & Experimental Next Steps

1. **Ribo-seq P-Site Offset Mapping:** Align public ribo-seq datasets targeting non-canonical start sites (CUG/GUG) to verify translation initiation.
2. **Chemical Crosslinking & Patch-Clamp:** For putative viroporin candidates with high amphipathic moments, perform planar lipid bilayer electrophysiology to test ion channel conductance.
3. **Tagged Expression:** Express with C-terminal FLAG tag in HEK293T / host cell lines.

---
*Report generated deterministically by `scripts/mining/run_non_canonical_mining.py`.*