# Biology as Information Architecture

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests: 44 Passed](https://img.shields.io/badge/tests-44%20passed-brightgreen.svg)](#testing)

An open-source, reproducible Python research framework examining biological information systems (DNA, RNA, regulation, proteins, cellular dynamics) through computational lenses—including information theory, software architecture, network science, and control theory.

---

## 🔬 Core Scientific Philosophy & Guardrails

Computational concepts are applied strictly as **measurable features and formal models**, never as assumed biological equivalences.

### 1. Explicit Evidence Classification
Every finding or claim emitted by this pipeline is categorized into one of four evidence classes:
- **`measurement`**: Directly computed from an identified empirical dataset. *(Phrased: "The analysis measured...")*
- **`simulation`**: Produced by a defined computational model. *(Phrased: "Under this model...")*
- **`interpretation`**: A plausible explanation connecting observations. *(Phrased: "One interpretation is...")*
- **`hypothesis`**: A falsifiable claim requiring experimental testing. *(Phrased: "This predicts...")*

### 2. Mandatory Controls & Null Models
- **Sequences**: Mono- and dinucleotide-preserving shuffled controls via Eulerian walks (Altschul-Erickson).
- **Networks**: Degree-preserving randomized network topologies via double-edge swaps.
- **Hypothesis Testing**: Multiple-testing corrections (Benjamini-Hochberg FDR) and effect sizes ($z$-scores) with empirical uncertainty intervals.
- **Causality Guardrail**: Never infer function, mechanism, or causality from statistical or topological correlation alone.

---

## 💻 Designed for Local Laptop Execution

- **100% Free & Local**: Runs entirely on your local CPU without requiring external paid LLM APIs or cloud accounts.
- **Storage-Conscious**: Streaming I/O and bounded parameter limits (e.g., $k$-mer limits $k \le 6$) prevent out-of-memory errors or excessive disk usage.
- **Full Reproducibility**: Tracks SHA-256 checksums, deterministic random seeds, execution timestamps, and exact environment package versions.

---

## 📁 Repository Structure

```
d:\DNA\
├── pyproject.toml                     # Package dependencies & configuration
├── config\
│   └── example.yaml                   # Analysis configuration template
├── src\bio_arch\
│   ├── contracts.py                   # Typed data models & JSON serialization
│   ├── provenance.py                  # SHA-256 hashing, seed derivation, environment capture
│   ├── logger.py                      # Structured logging utilities
│   ├── orchestrator.py                # Pipeline orchestrator and report builder
│   └── modules\
│       ├── ingestion.py               # Module 0: Data ingestion and validation
│       ├── information.py             # Module 1: DNA/RNA information architecture
│       ├── regulation.py              # Module 2: Regulatory network discovery
│       └── simulation.py              # Module 3: Emergent behavior & Boolean simulation
├── tests\                             # 44 unit and integration tests
├── examples\                          # Runnable demos and synthetic examples
├── data\                              # Local input datasets
├── docs\                              # Scientific limitations, reproducibility, and guides
└── outputs\                           # Generated run artifacts, reports, and manifests
```

---

## 🚀 Quickstart & Usage

### 1. Verification & Tests
Verify that all 44 unit tests pass on your machine:
```powershell
python -m pytest tests/
```

### 2. Run the Full Pipeline
Execute all modules end-to-end using the pipeline orchestrator:
```powershell
python -m bio_arch.orchestrator --config config/example.yaml --mode full
```
This automatically produces a unique run folder in `outputs/runs/` containing:
- `frozen_config.yaml`: Exact configuration snapshot.
- `execution.log`: Complete timestamped execution log.
- `final_summary.json`: Machine-readable results and manifest checksums.
- `final_report.md`: Human-readable scientific report with all 10 standard report sections.

### 3. Run Individual Modules via CLI

#### Module 0: Data Ingestion & Validation
```powershell
python -m bio_arch.modules.ingestion data/synthetic_sample.fasta --id lambda_demo --organism "Lambda Phage"
```

#### Module 1: DNA/RNA Information Architecture
```powershell
python -m bio_arch.modules.information data/synthetic_sample.fasta --shuffles 50 --seed 42
```

#### Module 2: Regulatory Network Discovery
```powershell
python -m bio_arch.modules.regulation tests/fixtures/valid_network.tsv --null-graphs 20
```

#### Module 3: Emergent Behavior & Simulation
```powershell
python -m bio_arch.modules.simulation tests/fixtures/valid_network.tsv --steps 20 --trials 10
```

---

## 📚 Documentation
- [Scientific Limitations & Boundaries](docs/scientific_limitations.md)
- [Reproducibility Statement](docs/reproducibility.md)
- [Contributing Guidelines](CONTRIBUTING.md)
- [Code of Conduct](CODE_OF_CONDUCT.md)
- [Security Policy](SECURITY.md)
- [Changelog](CHANGELOG.md)
