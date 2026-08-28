# Contributing to Biology as Information Architecture

Thank you for your interest in contributing! This project is an open, reproducible research framework. We welcome contributions that adhere to our scientific guardrails and engineering standards.

## Scientific Principles & Guardrails
All contributions introducing new analyses or interpretations must follow our core standards:
1. **Evidence Classification**:
   - `measurement`: Directly computed from an identified dataset. Must use: *"The analysis measured..."*
   - `simulation`: Produced by a defined computational model. Must use: *"Under this model..."*
   - `interpretation`: A plausible explanation connecting findings. Must use: *"One interpretation is..."*
   - `hypothesis`: A testable, falsifiable claim. Must use: *"This predicts..."*
2. **Mandatory Null Controls**:
   - Sequence analyses must compare against mono- and dinucleotide-preserving shuffled controls.
   - Network topology metrics must benchmark against degree-preserving randomized graphs.
   - All hypothesis tests must report effect sizes, uncertainty intervals, and multiple-testing corrected $p$-values.
3. **No Silent Repair**:
   - Ingestion parsers must never silently alter or repair data without recording normalizations in `ValidationReport`.
4. **No Causal Overstatement**:
   - Correlation or topological association alone must never be described as established biological function or causal mechanism.

## Development Workflow
1. Use Python 3.11+.
2. Follow strict type hints and docstrings.
3. Run tests using `pytest -v tests/`. All tests must pass before opening a Pull Request.
4. Keep memory and CPU footprints bounded so that research remains accessible to researchers using personal laptops.
