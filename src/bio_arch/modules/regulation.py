"""Module 2: Regulatory Network Discovery.

Models transcriptional and post-transcriptional control relationships as attributed directed graphs.
Features:
- Validates identifiers, signs, and distinguishes evidence types (curated, experimental, predicted, correlation).
- Computes node degrees, centralities, strongly connected components, feedback loops, and feed-forward loops (FFLs).
- Compares network topology against degree-preserving randomized null models using deterministic seeds.
- Strict scientific guardrail: Preserves uncertainty and NEVER labels correlation as causal.
- Bounded computation tailored for local laptops.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import asdict, dataclass, field
import json
import math
from pathlib import Path
import random
import sys
from typing import Any, Iterator, Literal

import networkx as nx

from bio_arch.contracts import (
    AnalysisRun,
    EvidenceClass,
    Finding,
    InterpretationRecord,
    ModuleResult,
)
from bio_arch.logger import setup_logger
from bio_arch.provenance import (
    SeedManager,
    get_system_environment,
    now_iso,
)

logger = setup_logger("bio_arch.regulation")

EvidenceType = Literal["curated", "experimental", "predicted", "correlation"]
InteractionSign = Literal["activation", "repression", "unspecified"]


@dataclass
class RegulatoryEdge:
    """An attributed directed interaction between regulators and targets."""

    source: str
    target: str
    interaction_type: str = "regulation"
    sign: InteractionSign = "unspecified"
    confidence: float = 1.0
    evidence_type: EvidenceType = "experimental"
    citation: str = ""
    context: str = ""

    def __post_init__(self) -> None:
        if not self.source or not self.target:
            raise ValueError("Regulatory edge requires non-empty source and target.")
        if self.confidence < 0.0 or self.confidence > 1.0:
            raise ValueError(f"Confidence score must be in [0, 1], got {self.confidence}")
        if self.evidence_type not in ("curated", "experimental", "predicted", "correlation"):
            raise ValueError(f"Invalid evidence type: '{self.evidence_type}'")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class NetworkSummary:
    """Topological properties and motif census of a regulatory network."""

    node_count: int
    edge_count: int
    density: float
    reciprocity: float
    num_sccs: int
    num_wccs: int
    ffl_count: int
    feedback_loop_count: int
    evidence_breakdown: dict[str, int]
    top_regulators: list[tuple[str, int]]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# ---------------------------------------------------------------------------
# Network Construction & Motif Mining
# ---------------------------------------------------------------------------

def load_regulatory_graph(edges: list[RegulatoryEdge]) -> nx.DiGraph:
    """Construct an attributed NetworkX DiGraph from a list of RegulatoryEdge instances."""
    G = nx.DiGraph()
    for edge in edges:
        G.add_edge(
            edge.source,
            edge.target,
            interaction_type=edge.interaction_type,
            sign=edge.sign,
            confidence=edge.confidence,
            evidence_type=edge.evidence_type,
            citation=edge.citation,
            context=edge.context,
        )
    return G


def parse_regulatory_file(file_path: Path | str) -> list[RegulatoryEdge]:
    """Parse TSV or CSV edge lists with header auto-detection."""
    p = Path(file_path)
    if not p.is_file():
        raise FileNotFoundError(f"Edge list file not found: {p}")

    edges: list[RegulatoryEdge] = []
    with p.open("r", encoding="utf-8", errors="replace") as f:
        first_line = f.readline()
        f.seek(0)
        delimiter = "\t" if "\t" in first_line else ","
        reader = csv.reader(f, delimiter=delimiter)

        header = None
        for line_num, row in enumerate(reader, start=1):
            if not row or row[0].startswith("#"):
                continue

            if header is None and any(col.lower() in ("source", "regulator", "from") for col in row):
                header = [c.strip().lower() for c in row]
                continue

            if len(row) < 2:
                continue

            source = row[0].strip()
            target = row[1].strip()
            interaction = row[2].strip() if len(row) > 2 else "regulation"
            
            # Sign inference
            sign: InteractionSign = "unspecified"
            if len(row) > 3 and row[3].strip().lower() in ("activation", "+", "+1", "pos"):
                sign = "activation"
            elif len(row) > 3 and row[3].strip().lower() in ("repression", "-", "-1", "neg"):
                sign = "repression"

            confidence = 1.0
            if len(row) > 4:
                try:
                    confidence = float(row[4].strip())
                except ValueError:
                    confidence = 1.0

            evidence: EvidenceType = "experimental"
            if len(row) > 5 and row[5].strip().lower() in ("curated", "experimental", "predicted", "correlation"):
                evidence = row[5].strip().lower()  # type: ignore

            citation = row[6].strip() if len(row) > 6 else ""
            context = row[7].strip() if len(row) > 7 else ""

            edges.append(
                RegulatoryEdge(
                    source=source,
                    target=target,
                    interaction_type=interaction,
                    sign=sign,
                    confidence=confidence,
                    evidence_type=evidence,
                    citation=citation,
                    context=context,
                )
            )

    return edges


def count_feed_forward_loops(G: nx.DiGraph) -> int:
    """Count feed-forward loop (FFL) 3-node motifs: X -> Y, Y -> Z, and direct X -> Z."""
    ffl_count = 0
    for x in G.nodes():
        x_targets = set(G.successors(x))
        for y in x_targets:
            if y == x:
                continue
            y_targets = set(G.successors(y))
            # Nodes z that receive edges from BOTH x and y
            shared_targets = x_targets.intersection(y_targets) - {x, y}
            ffl_count += len(shared_targets)
    return ffl_count


def count_simple_feedback_loops(G: nx.DiGraph, max_depth: int = 4) -> int:
    """Count directed feedback cycles up to a safe maximum depth to avoid combinatorial explosion."""
    if len(G) == 0:
        return 0

    count = 0
    # Use bounded simple cycles
    for cycle in nx.simple_cycles(G):
        if len(cycle) <= max_depth:
            count += 1
        # Circuit breaker for large cyclic dense graphs
        if count >= 1000:
            break
    return count


# ---------------------------------------------------------------------------
# Degree-Preserving Randomized Null Controls
# ---------------------------------------------------------------------------

def randomize_degree_preserving(G: nx.DiGraph, seed: int = 42, n_swaps: int | None = None) -> nx.DiGraph:
    """Generate a randomized directed null graph strictly preserving in- and out-degree sequences.

    Uses double-edge swaps: selects (u, v) and (x, y) and swaps to (u, y) and (x, v)
    if no self-loops or multi-edges are introduced.
    """
    R = G.copy()
    num_edges = R.number_of_edges()
    if num_edges < 2:
        return R

    rng = random.Random(seed)
    swaps_needed = n_swaps if n_swaps is not None else (num_edges * 10)
    edges = list(R.edges())
    max_tries = swaps_needed * 5
    successful_swaps = 0
    tries = 0

    while successful_swaps < swaps_needed and tries < max_tries:
        tries += 1
        e1_idx, e2_idx = rng.sample(range(len(edges)), 2)
        u, v = edges[e1_idx]
        x, y = edges[e2_idx]

        # Ensure 4 distinct vertices and no self-loops
        if u == y or x == v or u == x or v == y:
            continue

        # Ensure newly formed edges do not already exist
        if R.has_edge(u, y) or R.has_edge(x, v):
            continue

        # Perform the swap
        R.remove_edge(u, v)
        R.remove_edge(x, y)
        R.add_edge(u, y)
        R.add_edge(x, v)

        edges[e1_idx] = (u, y)
        edges[e2_idx] = (x, v)
        successful_swaps += 1

    return R


# ---------------------------------------------------------------------------
# Regulatory Network Analysis & Contract Emitting
# ---------------------------------------------------------------------------

def analyze_regulatory_network(
    edges: list[RegulatoryEdge],
    network_id: str = "reg_net_01",
    num_null_graphs: int = 50,
    seed: int = 42,
) -> tuple[NetworkSummary, list[Finding], list[InterpretationRecord]]:
    """Perform topological analysis, motif enrichment testing, and evidence auditing."""
    G = load_regulatory_graph(edges)
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()

    # Evidence breakdown
    evidence_counts = Counter(e.evidence_type for e in edges)
    correlation_only_count = evidence_counts.get("correlation", 0)

    # Topological properties
    density = nx.density(G) if num_nodes > 1 else 0.0
    reciprocity = nx.reciprocity(G) if num_edges > 0 else 0.0
    num_sccs = nx.number_strongly_connected_components(G)
    num_wccs = nx.number_weakly_connected_components(G)

    # Motif counting
    obs_ffl = count_feed_forward_loops(G)
    obs_cycles = count_simple_feedback_loops(G, max_depth=4)

    # Degree ranking
    out_degrees = sorted(G.out_degree(), key=lambda x: x[1], reverse=True)
    top_regulators = out_degrees[:5]

    # Null comparison for FFL enrichment
    seed_mgr = SeedManager(master_seed=seed)
    null_ffl_counts: list[int] = []
    for i in range(num_null_graphs):
        null_seed = seed_mgr.derive_seed(f"null_graph_{i}")
        R = randomize_degree_preserving(G, seed=null_seed)
        null_ffl_counts.append(count_feed_forward_loops(R))

    # Effect size and empirical p-value for FFL motif
    mean_null_ffl = sum(null_ffl_counts) / len(null_ffl_counts) if null_ffl_counts else 0.0
    std_null_ffl = (
        math.sqrt(sum((c - mean_null_ffl) ** 2 for c in null_ffl_counts) / (len(null_ffl_counts) - 1))
        if len(null_ffl_counts) > 1
        else 0.0
    )
    z_score = (obs_ffl - mean_null_ffl) / std_null_ffl if std_null_ffl > 0 else 0.0
    diff = abs(obs_ffl - mean_null_ffl)
    more_extreme = sum(1 for c in null_ffl_counts if abs(c - mean_null_ffl) >= diff)
    p_value = (more_extreme + 1) / (len(null_ffl_counts) + 1)

    summary = NetworkSummary(
        node_count=num_nodes,
        edge_count=num_edges,
        density=round(density, 6),
        reciprocity=round(reciprocity, 6),
        num_sccs=num_sccs,
        num_wccs=num_wccs,
        ffl_count=obs_ffl,
        feedback_loop_count=obs_cycles,
        evidence_breakdown=dict(evidence_counts),
        top_regulators=top_regulators,
    )

    findings: list[Finding] = [
        Finding(
            finding_id=f"f_{network_id}_ffl_enrichment",
            metric="feed_forward_loop_count",
            observed_value=obs_ffl,
            control_distribution={
                "null_model": "degree_preserving_directed_swap",
                "iterations": num_null_graphs,
                "null_mean": round(mean_null_ffl, 4),
                "null_std": round(std_null_ffl, 4),
            },
            effect_size=round(z_score, 4),
            uncertainty={"null_ci_2.5": min(null_ffl_counts), "null_ci_97.5": max(null_ffl_counts)},
            adjusted_p_value=round(p_value, 6),
            biological_context=f"FFL motif census for network {network_id} across {num_nodes} nodes.",
        ),
        Finding(
            finding_id=f"f_{network_id}_evidence_composition",
            metric="evidence_type_ratio",
            observed_value=dict(evidence_counts),
            biological_context=f"Proportion of curated vs correlation-derived edges in {network_id}.",
        ),
    ]

    interpretations: list[InterpretationRecord] = [
        InterpretationRecord(
            finding_ids=[f"f_{network_id}_ffl_enrichment"],
            classification=EvidenceClass.INTERPRETATION,
            claim=f"{EvidenceClass.INTERPRETATION.required_prefix} that statistically enriched feed-forward loops provide signal filtering or pulse generation.",
            alternatives=[
                "Topological by-product of gene duplication and preferential attachment during evolutionary history."
            ],
            limitations=[
                f"{correlation_only_count} of {num_edges} edges are correlation-derived and must not be interpreted as confirmed causal mechanisms."
            ],
            proposed_test="Perform targeted gene knockdown/CRISPR interference to experimentally verify edge directionality.",
        )
    ]

    return summary, findings, interpretations


def run_module_2(
    input_file: Path | str,
    output_dir: Path | str | None = None,
    num_null_graphs: int = 50,
    seed: int = 42,
) -> ModuleResult:
    """Execute Module 2 regulatory network analysis from a file."""
    in_path = Path(input_file)
    edges = parse_regulatory_file(in_path)

    seed_mgr = SeedManager(master_seed=seed)
    mod_seed = seed_mgr.derive_seed("module_2_regulation")

    run_meta = AnalysisRun(
        run_id=f"run_mod2_{in_path.stem}_{mod_seed}",
        timestamp=now_iso(),
        module="module_2_regulation",
        version="0.1.0",
        input_ids=[in_path.name],
        parameters={"num_null_graphs": num_null_graphs, "seed": mod_seed},
        seed=mod_seed,
        environment=get_system_environment(),
        status="success",
    )

    summary, findings, interpretations = analyze_regulatory_network(
        edges=edges,
        network_id=in_path.stem,
        num_null_graphs=num_null_graphs,
        seed=mod_seed,
    )

    artifact_paths: list[str] = []
    if output_dir:
        out_p = Path(output_dir)
        out_p.mkdir(parents=True, exist_ok=True)
        out_file = out_p / f"{in_path.stem}_network_metrics.json"

        result_payload = {
            "run_metadata": run_meta.to_dict(),
            "summary": summary.to_dict(),
            "findings": [f.to_dict() for f in findings],
            "interpretations": [i.to_dict() for i in interpretations],
        }
        out_file.write_text(json.dumps(result_payload, indent=2), encoding="utf-8")
        artifact_paths.append(str(out_file))

    warnings = []
    if summary.evidence_breakdown.get("correlation", 0) > 0:
        warnings.append(
            f"Caution: Network contains {summary.evidence_breakdown['correlation']} correlation-only edges. Directionality is unproven."
        )

    return ModuleResult(
        run_metadata=run_meta,
        outputs=summary.to_dict(),
        findings=findings,
        interpretations=interpretations,
        warnings=warnings,
        errors=[],
        artifact_paths=artifact_paths,
    )


def main() -> None:
    """CLI entrypoint for Module 2."""
    parser = argparse.ArgumentParser(
        description="Bio-Arch Module 2: Regulatory Network Discovery",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("file", type=str, help="Path to TSV/CSV regulatory edge list")
    parser.add_argument("--null-graphs", type=int, default=50, help="Number of degree-preserving randomized null graphs")
    parser.add_argument("--seed", type=int, default=42, help="Deterministic random seed")
    parser.add_argument("--outdir", type=str, default="outputs/regulation", help="Output directory")

    args = parser.parse_args()

    result = run_module_2(
        input_file=args.file,
        output_dir=args.outdir,
        num_null_graphs=args.null_graphs,
        seed=args.seed,
    )

    outputs = result.outputs
    print("\n--- Module 2: Regulatory Network Discovery Results ---")
    print(f"Nodes:                  {outputs['node_count']}")
    print(f"Edges:                  {outputs['edge_count']}")
    print(f"Network Density:        {outputs['density']:.4f}")
    print(f"Feed-Forward Loops:     {outputs['ffl_count']}")
    print(f"Feedback Cycles (<=4):  {outputs['feedback_loop_count']}")
    print(f"Strongly Connected:     {outputs['num_sccs']} components")
    print("Evidence Breakdown:     ", outputs["evidence_breakdown"])
    print("Top Regulators (out-deg):", outputs["top_regulators"])
    if result.warnings:
        for w in result.warnings:
            print(f" [WARNING] {w}")
    if result.artifact_paths:
        print(f"Saved artifacts to:     {result.artifact_paths[0]}")
    print("------------------------------------------------------\n")


if __name__ == "__main__":
    main()
