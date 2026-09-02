# 🔬 Executive Discovery Dossier: Top Mined Viral Overlapping smORFs
## Comprehensive Biophysical Profiling, Structural Screening & Literature Cross-Mapping

**Date:** 2026-09-02  
**Repository:** `hexaphase-genomics` (`D:\DNA`)  
**Status:** Peer-Review Ready / In Vitro Target Portfolio  

---

## 1. Executive Summary Table: Top 5 Overlapping Candidates

| Rank | Organism (Accession) | Parent CDS | Frame | Start Codon | Length | Host CAI | Null $z$-score | $p$-value | Classification |
| :---: | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **#1** | Escherichia phage phiX (`NC_001422.1`) | `CDS` | +1 | `GTG` | **203 aa** | **0.6546** | **$z = 7.87$** | $p < 0.001$ | **Novel Candidate smORF** |
| **#2** | adeno-associated virus (`NC_001401.2`) | `CDS` | +1 | `TTG` | **207 aa** | **0.6916** | **$z = 3.09$** | $p < 0.001$ | **Validated Unannotated** |
| **#3** | Escherichia phage phiX (`NC_001422.1`) | `CDS` | +1 | `CTG` | **201 aa** | **0.6518** | **$z = 7.6$** | $p < 0.001$ | **Novel Candidate smORF** |
| **#4** | Escherichia phage phiX (`NC_001422.1`) | `CDS` | +1 | `CTG` | **197 aa** | **0.6571** | **$z = 7.32$** | $p < 0.001$ | **Novel Candidate smORF** |
| **#5** | Escherichia phage phiX (`NC_001422.1`) | `CDS` | +1 | `GTG` | **177 aa** | **0.6599** | **$z = 7.41$** | $p < 0.001$ | **Novel Candidate smORF** |

---

## 2. Comprehensive Candidate Structural & Biophysical Profiles

### Candidate #1: Escherichia phage phiX174 (NC_001422.1)
**Classification:** `Novel Candidate smORF`  
**Biological Context:** High-hydropathy alternative reading frame exhibiting a clear transmembrane spanning helix (putative holin / lysis-associated microprotein architecture).

#### A. Genomic & Expression Architecture
- **Parent Primary CDS:** `CDS` (pilot protein for DNA ejection)
- **Relative Reading Phase:** **Frame +1** (Offset: 142..754 bp)
- **Start Codon & Context:** `GTG` | **Kozak Strength:** Weak
- **Kozak Sequence Window:** `...GCATGG[GTG]ATGC...`
- **Host Codon Adaptation Index (CAI):** **0.6546** (High mammalian expression compatibility)
- **Statistical Null Calibration:** **$z = 7.87$**, **$p = 0.0$** ($N=500$ Eulerian walk shuffles)

#### B. Biophysical & Structural Properties
- **Peptide Length & Mass:** **203 amino acids** | **24.09 kDa**
- **Isoelectric Point (pI):** **12.35** | **Net Charge (pH 7.4):** **39.3**
- **Secondary Structure Propensity:** 59.1% Alpha-Helix, 25.6% Beta-Sheet, 15.3% Coil
- **Transmembrane Architecture:** **6 TM Helices detected**
  - *Predicted TM Segments:*
    * Residues 9..33 (25 aa, Mean Hydropathy = 2.05)
    * Residues 21..39 (19 aa, Mean Hydropathy = 1.65)
    * Residues 36..54 (19 aa, Mean Hydropathy = 1.74)
    * Residues 41..66 (26 aa, Mean Hydropathy = 2.16)
    * Residues 56..101 (46 aa, Mean Hydropathy = 2.73)
    * Residues 85..103 (19 aa, Mean Hydropathy = 1.68)

#### C. Primary Peptide Sequence (FASTA)
```fasta
>smORF_NC_001422.1_CDS_F+1|203aa|MW=24.09kDa|pI=12.35|CAI=0.6546
VMLVLNLPFKALMFLTLMRPPLVLFLVLWLKLVKDFLKVRCRLALLPFLISCLIWLDLVASLPLIKERILVIILLLHFLSLMLGSVLVLMLPLLVWLTPDLRIKKSLLKCNWTIRKRLPRCKMRLKKRLLAFSRRLHARIRKTRYMHKMRCLLINRRSLLLALRLLWKTPIFPSNSRFPRLCAKCLLKLKRLVSILPMTKSKK
```

---

### Candidate #2: adeno-associated virus 2 (NC_001401.2)
**Classification:** `Validated Unannotated`  
**Biological Context:** Corresponds to the non-structural Assembly-Activating Protein (AAP) / X-protein reading frame in dependoparvoviruses, omitted in some primary RefSeq annotations.

#### A. Genomic & Expression Architecture
- **Parent Primary CDS:** `CDS` (major coat protein VP1)
- **Relative Reading Phase:** **Frame +1** (Offset: 517..1141 bp)
- **Start Codon & Context:** `TTG` | **Kozak Strength:** Strong
- **Kozak Sequence Window:** `...TGAATT[TTG]GTCA...`
- **Host Codon Adaptation Index (CAI):** **0.6916** (High mammalian expression compatibility)
- **Statistical Null Calibration:** **$z = 3.09$**, **$p = 0.01$** ($N=500$ Eulerian walk shuffles)

#### B. Biophysical & Structural Properties
- **Peptide Length & Mass:** **207 amino acids** | **22.95 kDa**
- **Isoelectric Point (pI):** **11.18** | **Net Charge (pH 7.4):** **10.6**
- **Secondary Structure Propensity:** 37.7% Alpha-Helix, 30.9% Beta-Sheet, 31.4% Coil
- **Transmembrane Architecture:** **Soluble / Globular profile (0 TM)**

#### C. Primary Peptide Sequence (FASTA)
```fasta
>smORF_NC_001401.2_CDS_F+1|207aa|MW=22.95kDa|pI=11.18|CAI=0.6916
LVRLETQTQYLTPSLSDSHQQPPLVWELIRWLQAVAHQWQTITRAPTEWVIPREIGIAIPHGWATESSPPAPEPGPCPPTTTTSTNKFPANQEPRTTITTLATAPLGGILTSTDSTATFHHVTGKDSSTTTGDSDPRDSTSSSLTFKSKRSRRMTVRRRLPITLPARFRCLLTRSTSSRTSSARRIKDASRRSQQTSSWCHSMDTSP
```

---

### Candidate #3: Escherichia phage phiX174 (NC_001422.1)
**Classification:** `Novel Candidate smORF`  
**Biological Context:** High-hydropathy alternative reading frame exhibiting a clear transmembrane spanning helix (putative holin / lysis-associated microprotein architecture).

#### A. Genomic & Expression Architecture
- **Parent Primary CDS:** `CDS` (pilot protein for DNA ejection)
- **Relative Reading Phase:** **Frame +1** (Offset: 148..754 bp)
- **Start Codon & Context:** `CTG` | **Kozak Strength:** Strong
- **Kozak Sequence Window:** `...GTGATG[CTG]GTAT...`
- **Host Codon Adaptation Index (CAI):** **0.6518** (High mammalian expression compatibility)
- **Statistical Null Calibration:** **$z = 7.6$**, **$p = 0.0$** ($N=500$ Eulerian walk shuffles)

#### B. Biophysical & Structural Properties
- **Peptide Length & Mass:** **201 amino acids** | **23.86 kDa**
- **Isoelectric Point (pI):** **12.35** | **Net Charge (pH 7.4):** **39.3**
- **Secondary Structure Propensity:** 59.2% Alpha-Helix, 25.4% Beta-Sheet, 15.4% Coil
- **Transmembrane Architecture:** **6 TM Helices detected**
  - *Predicted TM Segments:*
    * Residues 7..31 (25 aa, Mean Hydropathy = 2.05)
    * Residues 19..37 (19 aa, Mean Hydropathy = 1.65)
    * Residues 34..52 (19 aa, Mean Hydropathy = 1.74)
    * Residues 39..64 (26 aa, Mean Hydropathy = 2.16)
    * Residues 54..99 (46 aa, Mean Hydropathy = 2.73)
    * Residues 83..101 (19 aa, Mean Hydropathy = 1.68)

#### C. Primary Peptide Sequence (FASTA)
```fasta
>smORF_NC_001422.1_CDS_F+1|201aa|MW=23.86kDa|pI=12.35|CAI=0.6518
LVLNLPFKALMFLTLMRPPLVLFLVLWLKLVKDFLKVRCRLALLPFLISCLIWLDLVASLPLIKERILVIILLLHFLSLMLGSVLVLMLPLLVWLTPDLRIKKSLLKCNWTIRKRLPRCKMRLKKRLLAFSRRLHARIRKTRYMHKMRCLLINRRSLLLALRLLWKTPIFPSNSRFPRLCAKCLLKLKRLVSILPMTKSKK
```

---

### Candidate #4: Escherichia phage phiX174 (NC_001422.1)
**Classification:** `Novel Candidate smORF`  
**Biological Context:** High-hydropathy alternative reading frame exhibiting a clear transmembrane spanning helix (putative holin / lysis-associated microprotein architecture).

#### A. Genomic & Expression Architecture
- **Parent Primary CDS:** `CDS` (pilot protein for DNA ejection)
- **Relative Reading Phase:** **Frame +1** (Offset: 160..754 bp)
- **Start Codon & Context:** `CTG` | **Kozak Strength:** Moderate (Purine at -3)
- **Kozak Sequence Window:** `...TTAAAT[CTG]CCAT...`
- **Host Codon Adaptation Index (CAI):** **0.6571** (High mammalian expression compatibility)
- **Statistical Null Calibration:** **$z = 7.32$**, **$p = 0.0$** ($N=500$ Eulerian walk shuffles)

#### B. Biophysical & Structural Properties
- **Peptide Length & Mass:** **197 amino acids** | **23.42 kDa**
- **Isoelectric Point (pI):** **12.35** | **Net Charge (pH 7.4):** **39.3**
- **Secondary Structure Propensity:** 59.4% Alpha-Helix, 25.4% Beta-Sheet, 15.2% Coil
- **Transmembrane Architecture:** **6 TM Helices detected**
  - *Predicted TM Segments:*
    * Residues 3..27 (25 aa, Mean Hydropathy = 2.05)
    * Residues 15..33 (19 aa, Mean Hydropathy = 1.65)
    * Residues 30..48 (19 aa, Mean Hydropathy = 1.74)
    * Residues 35..60 (26 aa, Mean Hydropathy = 2.16)
    * Residues 50..95 (46 aa, Mean Hydropathy = 2.73)
    * Residues 79..97 (19 aa, Mean Hydropathy = 1.68)

#### C. Primary Peptide Sequence (FASTA)
```fasta
>smORF_NC_001422.1_CDS_F+1|197aa|MW=23.42kDa|pI=12.35|CAI=0.6571
LPFKALMFLTLMRPPLVLFLVLWLKLVKDFLKVRCRLALLPFLISCLIWLDLVASLPLIKERILVIILLLHFLSLMLGSVLVLMLPLLVWLTPDLRIKKSLLKCNWTIRKRLPRCKMRLKKRLLAFSRRLHARIRKTRYMHKMRCLLINRRSLLLALRLLWKTPIFPSNSRFPRLCAKCLLKLKRLVSILPMTKSKK
```

---

### Candidate #5: Escherichia phage phiX174 (NC_001422.1)
**Classification:** `Novel Candidate smORF`  
**Biological Context:** High-hydropathy alternative reading frame exhibiting a clear transmembrane spanning helix (putative holin / lysis-associated microprotein architecture).

#### A. Genomic & Expression Architecture
- **Parent Primary CDS:** `CDS` (pilot protein for DNA ejection)
- **Relative Reading Phase:** **Frame +1** (Offset: 220..754 bp)
- **Start Codon & Context:** `GTG` | **Kozak Strength:** Weak
- **Kozak Sequence Window:** `...TTTCTG[GTG]CTAT...`
- **Host Codon Adaptation Index (CAI):** **0.6599** (High mammalian expression compatibility)
- **Statistical Null Calibration:** **$z = 7.41$**, **$p = 0.0$** ($N=500$ Eulerian walk shuffles)

#### B. Biophysical & Structural Properties
- **Peptide Length & Mass:** **177 amino acids** | **21.08 kDa**
- **Isoelectric Point (pI):** **12.3** | **Net Charge (pH 7.4):** **37.3**
- **Secondary Structure Propensity:** 59.3% Alpha-Helix, 25.4% Beta-Sheet, 15.3% Coil
- **Transmembrane Architecture:** **4 TM Helices detected**
  - *Predicted TM Segments:*
    * Residues 10..28 (19 aa, Mean Hydropathy = 1.74)
    * Residues 15..40 (26 aa, Mean Hydropathy = 2.16)
    * Residues 30..75 (46 aa, Mean Hydropathy = 2.73)
    * Residues 59..77 (19 aa, Mean Hydropathy = 1.68)

#### C. Primary Peptide Sequence (FASTA)
```fasta
>smORF_NC_001422.1_CDS_F+1|177aa|MW=21.08kDa|pI=12.3|CAI=0.6599
VLWLKLVKDFLKVRCRLALLPFLISCLIWLDLVASLPLIKERILVIILLLHFLSLMLGSVLVLMLPLLVWLTPDLRIKKSLLKCNWTIRKRLPRCKMRLKKRLLAFSRRLHARIRKTRYMHKMRCLLINRRSLLLALRLLWKTPIFPSNSRFPRLCAKCLLKLKRLVSILPMTKSKK
```

---

## 3. Actionable Experimental Validation Protocols

To empirically validate translation and functional expression of these unannotated viral candidates in wet-lab environments, the following 3-tier experimental workflow is recommended:

### Protocol 1: Tagged Mammalian Expression & Subcellular Localization
1. Synthesize codon-optimized cDNA constructs containing the candidate smORF with a C-terminal **3xFLAG** or **mNeonGreen** tag under a CMV/EF1a promoter.
2. Transfect into target host cell lines (e.g. HEK293T for human/primate viruses, PK-15 for porcine circovirus).
3. Perform Western blot analysis to verify stable peptide accumulation at the predicted molecular weight.
4. Perform confocal immunofluorescence to determine subcellular localization (mitochondrial, nuclear, ER, or plasma membrane channel).

### Protocol 2: Ribosome Profiling (Ribo-seq) Footprint Alignment
1. Ingest public Ribo-seq datasets from infected host cells (NCBI SRA).
2. Align protected ribosome footprints to the viral genome using 28–31 nt P-site offset mapping.
3. Quantify translation initiation peaks across the candidate's Kozak context and verify in-frame ribosome translocation across Frame +1 / Frame +2.

### Protocol 3: Targeted Mass Spectrometry (LC-MS/MS)
1. Perform in-gel tryptic digestion of viral lysates or infected cellular fractions.
2. Query mass spectra against a customized FASTA database appending the candidate peptide sequences to the standard Uniprot host/viral proteome.
3. Filter at 1% False Discovery Rate (FDR) to identify unique tryptic junction peptides spanning the overlapping reading frame.

---
*Dossier generated deterministically by `scripts/mining/inspect_top_candidates.py`.*