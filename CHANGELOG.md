# Changelog

All notable changes to the "Biology as Information Architecture" project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-28

### Added
- **Repository Scaffold & Contracts**:
  - Typed models for `DatasetManifest`, `AnalysisRun`, `Finding`, `InterpretationRecord`, and `ModuleResult`.
  - Strict evidence classification: `measurement`, `simulation`, `interpretation`, and `hypothesis`.
  - Provenance utilities: streaming 64 KB chunked SHA-256 calculation, deterministic `SeedManager`, and runtime environment logging.
- **Module 0 (Data Ingestion & Validation)**:
  - Streaming parsers for FASTA, FASTQ, BED, and Network Edge Lists (TSV/CSV).
  - Explicit normalization tracking (zero silent data repair).
  - Validation reports with error and warning severity classifications.
- **Module 1 (DNA/RNA Information Architecture)**:
  - Nucleotide composition, GC content, GC skew, and AT skew.
  - Shannon entropy and conditional entropy ($k \le 6$).
  - Compression ratio (zlib/DEFLATE) complexity proxy.
  - Palindromes, tandem repeats, and match autocorrelation.
  - Altschul-Erickson dinucleotide-preserving Eulerian walk shuffles.
  - Benjamini-Hochberg False Discovery Rate (FDR) multiple-testing correction.
- **Module 2 (Regulatory Network Discovery)**:
  - NetworkX attributed directed graph modeling.
  - Separation of evidence types: curated, experimental, predicted, and correlation.
  - Topological metrics: reciprocity, density, strongly connected components, feedback cycles.
  - 3-node Feed-Forward Loop (FFL) motif census.
  - Degree-preserving randomized null network generation.
- **Module 3 (Emergent Behavior & Simulation)**:
  - Transparent Boolean network dynamic simulator.
  - Synchronous and asynchronous update modes.
  - Single-node in-silico knockouts and forced activation sensitivity analysis.
  - Attractor detection (fixed points and limit cycles) and Hamming distance robustness tracking.
- **Orchestrator**:
  - YAML configuration loader with resource bounds.
  - Modes: `full`, `single`, `dry-run`, `resume`, `validate-only`.
  - Frozen configuration and structured execution logging.
  - Automated generation of machine-readable `final_summary.json` and human-readable `final_report.md`.
- **Test Suite**:
  - 44 comprehensive unit and integration tests passing in < 1 second.
