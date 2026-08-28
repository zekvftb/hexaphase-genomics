# Scientific Limitations & Boundaries

The *Biology as Information Architecture* project uses computational lenses (software engineering, information theory, network science, dynamical systems) to extract measurable patterns from biological data.

To prevent overreach and maintain scientific integrity, researchers using this platform must adhere to the following limitations:

## 1. Analogy vs. Mechanism
- **Metaphor is not Mechanism**: A sequence region exhibiting high compressibility or a network motif resembling a logic gate does **not** prove that the biological cell operates like digital hardware.
- Computational concepts are descriptive lenses and predictive models, not literal biological equivalences.

## 2. Statistical Significance vs. Biological Function
- Statistical rarity relative to a null model (such as a dinucleotide-preserving sequence shuffle or degree-preserving graph swap) does not inherently imply biological function, adaptive selection, or evolutionary design.
- Mutational biases, structural physics, replication timing, and neutral drift can produce complex non-random patterns naturally.

## 3. Correlation vs. Causality
- Edge associations derived from high-throughput co-expression or proximity ligation (Hi-C) must be categorized strictly as `correlation` or `predicted`, never as confirmed causal regulatory mechanisms.
- Directionality and causal interaction require targeted perturbation assays (e.g. CRISPR interference, gene knockout).

## 4. Modeling Simplifications
- **Boolean Discrete States**: Boolean network simulations discretize continuous biochemical concentrations into binary $\{0, 1\}$ states. Attractors identified in Boolean models represent qualitative dynamical possibilities and must be verified against continuous differential equation models and experimental assays.
- **Compression Heuristics**: Algorithms like `zlib` (LZ77/Huffman) are practical proxies for descriptive sequence complexity, not mathematically uncomputable Kolmogorov complexity.
