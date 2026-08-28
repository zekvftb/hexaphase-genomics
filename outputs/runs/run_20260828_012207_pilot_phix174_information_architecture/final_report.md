# Scientific Analysis Report: pilot_phix174_information_architecture
**Run ID**: `run_20260828_012207_pilot_phix174_information_architecture` | **Status**: `COMPLETED` | **Date**: 2026-08-28T01:22:08.939662+00:00

## 1. Run Summary
- **Execution Status**: completed
- **Global Seed**: 2026
- **Modules Executed**: module_0_ingestion, module_1_information, module_2_regulation, module_3_simulation
- **Total Duration**: 1.15s

## 2. Datasets and Provenance
- **ds_phix174_refseq** (Bacteriophage phiX174 (NC_001422.1))
  - Checksum (SHA-256): `5af131154db6ed660b9a54da006c10eb3d2187d038982ad98b2411578c26d810`
  - Source: NCBI Entrez / RefSeq NC_001422.1 | License: NCBI RefSeq Public Domain
  - Validation: `valid`

## 3. Methods and Parameters
- **Config snapshot**: Saved at `frozen_config.yaml`
- **Resource Limits**: {"max_kmer_length": 6, "max_sequence_length_mb": 25, "max_permutations": 500}

## 4. Empirical Measurements
- **ingestion_record_count**: `{'total': 7, 'valid': 7}` — *Ingested fasta records for Bacteriophage phiX174 (NC_001422.1).*
- **gc_content**: `0.447642` — *GC content and base composition for phix174_whole_genome.*
- **gc_content**: `0.478261` — *GC content and base composition for phix174_overlap_gene_e.*
- **gc_content**: `0.403509` — *GC content and base composition for phix174_overlap_gene_k.*
- **gc_content**: `0.449377` — *GC content and base composition for phix174_nonoverlap_gene_f.*
- **gc_content**: `0.420455` — *GC content and base composition for phix174_nonoverlap_gene_g.*
- **gc_content**: `0.450861` — *GC content and base composition for phix174_nonoverlap_gene_h.*
- **gc_content**: `0.516393` — *GC content and base composition for phix174_intergenic_jf.*
- **evidence_type_ratio**: `{'experimental': 3}` — *Proportion of curated vs correlation-derived edges in valid_network.*
- **max_knockout_hamming_distance**: `3` — *Maximum state shift resulting from single-node in-silico knockout.*

## 5. Control Comparisons (Null Models)
- **shannon_entropy_k2**: observed=3.951712 (z=0.995, adj-p=1.0)
  - Context: 2-mer Shannon entropy compared to dinucleotide-preserving controls.
  - Null: dinucleotide_shuffle (N=100)
- **conditional_entropy_k1**: observed=1.967133 (z=0.995, adj-p=1.0)
  - Context: First-order conditional entropy H(X_t | X_{t-1}).
  - Null: dinucleotide_shuffle (N=100)
- **compression_ratio_zlib**: observed=0.315262 (z=-5.6677, adj-p=0.029703)
  - Context: Descriptive compression ratio proxy for sequence redundancy.
  - Null: dinucleotide_shuffle (N=100)
- **shannon_entropy_k2**: observed=3.941424 (z=0.995, adj-p=1.0)
  - Context: 2-mer Shannon entropy compared to dinucleotide-preserving controls.
  - Null: dinucleotide_shuffle (N=100)
- **conditional_entropy_k1**: observed=1.947064 (z=-0.995, adj-p=1.0)
  - Context: First-order conditional entropy H(X_t | X_{t-1}).
  - Null: dinucleotide_shuffle (N=100)
- **compression_ratio_zlib**: observed=0.427536 (z=-1.2988, adj-p=1.0)
  - Context: Descriptive compression ratio proxy for sequence redundancy.
  - Null: dinucleotide_shuffle (N=100)
- **shannon_entropy_k2**: observed=3.82979 (z=-0.995, adj-p=1.0)
  - Context: 2-mer Shannon entropy compared to dinucleotide-preserving controls.
  - Null: dinucleotide_shuffle (N=100)
- **conditional_entropy_k1**: observed=1.859535 (z=0.995, adj-p=1.0)
  - Context: First-order conditional entropy H(X_t | X_{t-1}).
  - Null: dinucleotide_shuffle (N=100)
- **compression_ratio_zlib**: observed=0.48538 (z=-0.1788, adj-p=1.0)
  - Context: Descriptive compression ratio proxy for sequence redundancy.
  - Null: dinucleotide_shuffle (N=100)
- **shannon_entropy_k2**: observed=3.929822 (z=-0.995, adj-p=1.0)
  - Context: 2-mer Shannon entropy compared to dinucleotide-preserving controls.
  - Null: dinucleotide_shuffle (N=100)
- **conditional_entropy_k1**: observed=1.957952 (z=-0.995, adj-p=1.0)
  - Context: First-order conditional entropy H(X_t | X_{t-1}).
  - Null: dinucleotide_shuffle (N=100)
- **compression_ratio_zlib**: observed=0.355919 (z=1.4819, adj-p=0.475248)
  - Context: Descriptive compression ratio proxy for sequence redundancy.
  - Null: dinucleotide_shuffle (N=100)
- **shannon_entropy_k2**: observed=3.871008 (z=-0.995, adj-p=1.0)
  - Context: 2-mer Shannon entropy compared to dinucleotide-preserving controls.
  - Null: dinucleotide_shuffle (N=100)
- **conditional_entropy_k1**: observed=1.92966 (z=-0.995, adj-p=1.0)
  - Context: First-order conditional entropy H(X_t | X_{t-1}).
  - Null: dinucleotide_shuffle (N=100)
- **compression_ratio_zlib**: observed=0.371212 (z=-2.7958, adj-p=0.059406)
  - Context: Descriptive compression ratio proxy for sequence redundancy.
  - Null: dinucleotide_shuffle (N=100)
- **shannon_entropy_k2**: observed=3.915709 (z=0.995, adj-p=1.0)
  - Context: 2-mer Shannon entropy compared to dinucleotide-preserving controls.
  - Null: dinucleotide_shuffle (N=100)
- **conditional_entropy_k1**: observed=1.929787 (z=-0.995, adj-p=1.0)
  - Context: First-order conditional entropy H(X_t | X_{t-1}).
  - Null: dinucleotide_shuffle (N=100)
- **compression_ratio_zlib**: observed=0.345491 (z=-2.7513, adj-p=0.029703)
  - Context: Descriptive compression ratio proxy for sequence redundancy.
  - Null: dinucleotide_shuffle (N=100)
- **shannon_entropy_k2**: observed=3.904327 (z=-0.995, adj-p=1.0)
  - Context: 2-mer Shannon entropy compared to dinucleotide-preserving controls.
  - Null: dinucleotide_shuffle (N=100)
- **conditional_entropy_k1**: observed=1.917653 (z=-0.995, adj-p=1.0)
  - Context: First-order conditional entropy H(X_t | X_{t-1}).
  - Null: dinucleotide_shuffle (N=100)
- **compression_ratio_zlib**: observed=0.557377 (z=1.0461, adj-p=1.0)
  - Context: Descriptive compression ratio proxy for sequence redundancy.
  - Null: dinucleotide_shuffle (N=100)
- **feed_forward_loop_count**: observed=0 (z=0.0, adj-p=1.0)
  - Context: FFL motif census for network valid_network across 3 nodes.
  - Null: degree_preserving_directed_swap (N=50)
- **unique_attractors_found**: observed=5 (N/A, N/A)
  - Context: Attractor landscape under synchronous Boolean rules for valid_network.
  - Null: None (N=None)

## 6. Simulations
- Under this model, the network converges to 5 distinct attractor states from 10 random starts.
  - *Limitations*: Boolean assumptions discretize expression levels into binary 0/1 states.

## 7. Interpretations and Alternatives
- **Claim**: One interpretation is that deviations in compressibility from null models indicate higher-order repeat structure or motifs.
  - **Alternative Explanations**: Local GC skews or replication origins causing periodic structural motifs without selective pressure.
  - **Limitations**: Compression algorithms (zlib) are heuristics, not exact Kolmogorov complexity calculations.
- **Claim**: One interpretation is that deviations in compressibility from null models indicate higher-order repeat structure or motifs.
  - **Alternative Explanations**: Local GC skews or replication origins causing periodic structural motifs without selective pressure.
  - **Limitations**: Compression algorithms (zlib) are heuristics, not exact Kolmogorov complexity calculations.
- **Claim**: One interpretation is that deviations in compressibility from null models indicate higher-order repeat structure or motifs.
  - **Alternative Explanations**: Local GC skews or replication origins causing periodic structural motifs without selective pressure.
  - **Limitations**: Compression algorithms (zlib) are heuristics, not exact Kolmogorov complexity calculations.
- **Claim**: One interpretation is that deviations in compressibility from null models indicate higher-order repeat structure or motifs.
  - **Alternative Explanations**: Local GC skews or replication origins causing periodic structural motifs without selective pressure.
  - **Limitations**: Compression algorithms (zlib) are heuristics, not exact Kolmogorov complexity calculations.
- **Claim**: One interpretation is that deviations in compressibility from null models indicate higher-order repeat structure or motifs.
  - **Alternative Explanations**: Local GC skews or replication origins causing periodic structural motifs without selective pressure.
  - **Limitations**: Compression algorithms (zlib) are heuristics, not exact Kolmogorov complexity calculations.
- **Claim**: One interpretation is that deviations in compressibility from null models indicate higher-order repeat structure or motifs.
  - **Alternative Explanations**: Local GC skews or replication origins causing periodic structural motifs without selective pressure.
  - **Limitations**: Compression algorithms (zlib) are heuristics, not exact Kolmogorov complexity calculations.
- **Claim**: One interpretation is that deviations in compressibility from null models indicate higher-order repeat structure or motifs.
  - **Alternative Explanations**: Local GC skews or replication origins causing periodic structural motifs without selective pressure.
  - **Limitations**: Compression algorithms (zlib) are heuristics, not exact Kolmogorov complexity calculations.
- **Claim**: One interpretation is that statistically enriched feed-forward loops provide signal filtering or pulse generation.
  - **Alternative Explanations**: Topological by-product of gene duplication and preferential attachment during evolutionary history.
  - **Limitations**: 0 of 3 edges are correlation-derived and must not be interpreted as confirmed causal mechanisms.

## 8. Testable Hypotheses
- Derived from `f_ingest_ds_phix174_refseq`: Verify against reference assembly or independent biological assay annotations.
- Derived from `f_phix174_whole_genome_compression_ratio, f_phix174_whole_genome_entropy_k2`: Perform comparative analysis against known non-coding neutral regions or simulated genomes.
- Derived from `f_phix174_overlap_gene_e_compression_ratio, f_phix174_overlap_gene_e_entropy_k2`: Perform comparative analysis against known non-coding neutral regions or simulated genomes.
- Derived from `f_phix174_overlap_gene_k_compression_ratio, f_phix174_overlap_gene_k_entropy_k2`: Perform comparative analysis against known non-coding neutral regions or simulated genomes.
- Derived from `f_phix174_nonoverlap_gene_f_compression_ratio, f_phix174_nonoverlap_gene_f_entropy_k2`: Perform comparative analysis against known non-coding neutral regions or simulated genomes.
- Derived from `f_phix174_nonoverlap_gene_g_compression_ratio, f_phix174_nonoverlap_gene_g_entropy_k2`: Perform comparative analysis against known non-coding neutral regions or simulated genomes.
- Derived from `f_phix174_nonoverlap_gene_h_compression_ratio, f_phix174_nonoverlap_gene_h_entropy_k2`: Perform comparative analysis against known non-coding neutral regions or simulated genomes.
- Derived from `f_phix174_intergenic_jf_compression_ratio, f_phix174_intergenic_jf_entropy_k2`: Perform comparative analysis against known non-coding neutral regions or simulated genomes.
- Derived from `f_valid_network_ffl_enrichment`: Perform targeted gene knockdown/CRISPR interference to experimentally verify edge directionality.
- Derived from `f_valid_network_attractor_count`: Compare predicted attractors with single-cell RNA-seq cluster expression profiles.

## 9. Limitations & Guardrails
- Statistical patterns and computational models must not be treated as literal biological mechanisms.
- Correlation does not imply regulatory or causal direction.
- Measurements are restricted to documented datasets and parameter bounds.

## 10. Next Experiments
- Expand testing to diverse homologous regions across related prokaryotic and eukaryotic clades.
- Benchmark Boolean state attractor predictions against single-cell dynamic assays.