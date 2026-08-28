# Reproducibility Statement

The *Biology as Information Architecture* framework is designed for end-to-end, deterministic scientific reproducibility.

## How Reproducibility is Guaranteed

1. **Deterministic Random Seeds**:
   - Every randomized process (sequence shuffling, network degree-preserving swaps, initial simulation state sampling) utilizes an isolated random number generator initialized by `SeedManager`.
   - The master seed derives context-specific sub-seeds deterministically via SHA-256 digests.

2. **Complete Provenance Tracking**:
   - Ingested files have their exact SHA-256 content checksum recorded in a `DatasetManifest`.
   - Normalizations (such as case folding or coordinate adjustments) are explicitly recorded; raw data is never silently repaired.
   - Every module execution stores an `AnalysisRun` record documenting the software version, parameters, seed, and complete Python environment metadata (platform, OS, installed package versions).

3. **Frozen Configurations**:
   - When the orchestrator executes, it immediately snapshots the input YAML configuration into `frozen_config.yaml` within the unique run directory.

4. **Self-Contained Local Execution**:
   - The framework has zero dependencies on external cloud APIs or online services. Any run can be repeated identically on an air-gapped machine using the same seed and input datasets.
