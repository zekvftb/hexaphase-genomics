# Data Directory

This directory holds local input datasets for the **Biology as Information Architecture** research pipeline.

## Guidelines for Laptop Environments

To ensure fast analysis and prevent disk storage exhaustion:
1. **Lightweight Files**: For initial explorations and testing, use small sequences, targeted regions, viral/bacterial genomes, or synthetic benchmark sequences (< 10 MB).
2. **Never Silently Repair Data**: Module 0 will validate files as-is. If a sequence is malformed, the pipeline will flag or reject it, documenting the issue in a validation report.
3. **Licensing and Provenance**: Always record the source, license, and organism when adding new datasets.
4. **Git Tracking**: Raw sequence files (`*.fasta`, `*.fastq`, `*.gff`, etc.) should generally be kept local or referenced via URLs rather than committed into version control if large.
