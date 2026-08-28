# Stage 7 Pilot Preregistration

## 1. Study Title
**Information Architecture of Dual-Frame Overlapping Genes in Bacteriophage $\Phi\text{X174}$**

## 2. Primary Research Question
Does the higher-order information structure (Shannon entropy, sequence compressibility, and $k$-mer distribution) of dual-frame overlapping coding regions in Bacteriophage $\Phi\text{X174}$ exhibit greater statistical constraint than single-frame coding regions and intergenic regions when benchmarked against dinucleotide-preserving null controls?

## 3. Scope Rule Compliance
- **One Primary Question**: Do overlapping coding frames constrain sequence information architecture beyond local dinucleotide composition?
- **One Primary Outcome**: Empirical effect size ($z$-score) and Benjamini-Hochberg FDR-adjusted $p$-value for sequence compressibility (zlib proxy) and 2-mer Shannon entropy.
- **One Null Model**: Dinucleotide-preserving Eulerian walk sequence shuffles (Altschul-Erickson algorithm, 100 permutations per region).
- **One Manageable Dataset**: Bacteriophage $\Phi\text{X174}$ complete genome (NCBI RefSeq: `NC_001422.1`, 5,386 bp, single-stranded circular DNA).

## 4. Hypotheses & Evidence Classification
- **$H_0$ (Null Hypothesis)**: The sequence compressibility and 2-mer Shannon entropy of overlapping genes are fully explained by single-nucleotide and dinucleotide transition frequencies ($z \approx 0, p_{\text{adj}} > 0.05$).
- **$H_1$ (Alternative Hypothesis)**: Simultaneous translation in two reading frames imposes multi-nucleotide coding constraints that reduce sequence compressibility and alter higher-order entropy relative to dinucleotide-matched controls ($|z| > 2.0, p_{\text{adj}} \le 0.05$).
- **Evidence Classification Standard**:
  - Raw computed values will be reported as `MEASUREMENT` ("The analysis measured...").
  - Comparisons against shuffled controls will be reported as `INTERPRETATION` ("One interpretation is...").
  - Predictions regarding selective pressure will be reported as `HYPOTHESIS` ("This predicts...").

## 5. Partitions of Analysis
The 5,386 bp genome will be segmented into 4 functional partitions:
1. **Whole Genome**: Complete 5,386 bp sequence.
2. **Overlapping Coding Regions**:
   - Gene B inside Gene A (coords: 5075..5386 / 1..51)
   - Gene E inside Gene D (coords: 568..843)
   - Gene K overlapping Genes A and C (coords: 51..221)
3. **Non-Overlapping Coding Regions**:
   - Gene F (Major coat protein, coords: 1001..2284)
   - Gene G (Major spike protein, coords: 2395..2922)
   - Gene H (DNA pilot protein, coords: 2931..3917)
4. **Intergenic / Regulatory Spacers**: Regions between defined CDS features.
