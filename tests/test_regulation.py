"""Unit tests for Module 2: Regulatory Network Discovery."""

from pathlib import Path
import networkx as nx
import pytest

from bio_arch.contracts import EvidenceClass
from bio_arch.modules.regulation import (
    count_feed_forward_loops,
    count_simple_feedback_loops,
    load_regulatory_graph,
    parse_regulatory_file,
    randomize_degree_preserving,
    RegulatoryEdge,
    run_module_2,
)


def test_regulatory_edge_validation():
    """Verify field constraints and validations on RegulatoryEdge."""
    # Valid edge
    e = RegulatoryEdge(
        source="TF1",
        target="GeneA",
        interaction_type="activation",
        sign="activation",
        confidence=0.9,
        evidence_type="curated",
    )
    assert e.source == "TF1"
    assert e.confidence == 0.9

    # Invalid confidence (< 0 or > 1)
    with pytest.raises(ValueError, match="Confidence score must be in"):
        RegulatoryEdge(source="TF1", target="GeneA", confidence=1.5)

    # Invalid evidence type
    with pytest.raises(ValueError, match="Invalid evidence type"):
        RegulatoryEdge(source="TF1", target="GeneA", evidence_type="hearsay")  # type: ignore


def test_feed_forward_loop_counting():
    """Verify exact identification of 3-node Feed-Forward Loop (FFL) motifs."""
    # Build canonical FFL: X -> Y, Y -> Z, X -> Z
    G = nx.DiGraph()
    G.add_edge("X", "Y")
    G.add_edge("Y", "Z")
    G.add_edge("X", "Z")

    assert count_feed_forward_loops(G) == 1

    # Linear cascade without shortcut: X -> Y, Y -> Z
    G_linear = nx.DiGraph()
    G_linear.add_edge("X", "Y")
    G_linear.add_edge("Y", "Z")
    assert count_feed_forward_loops(G_linear) == 0


def test_simple_feedback_cycle_counting():
    """Verify directed feedback cycle counting."""
    # 3-node cycle: A -> B -> C -> A
    G = nx.DiGraph()
    G.add_edge("A", "B")
    G.add_edge("B", "C")
    G.add_edge("C", "A")

    assert count_simple_feedback_loops(G, max_depth=4) == 1


def test_degree_preserving_randomization():
    """Verify that randomize_degree_preserving preserves exact degree sequences."""
    # Build small bipartite-like network
    edges = [
        RegulatoryEdge("TF1", "G1"),
        RegulatoryEdge("TF1", "G2"),
        RegulatoryEdge("TF2", "G2"),
        RegulatoryEdge("TF2", "G3"),
        RegulatoryEdge("TF3", "G1"),
        RegulatoryEdge("TF3", "G3"),
    ]
    G = load_regulatory_graph(edges)

    in_degrees_orig = dict(G.in_degree())
    out_degrees_orig = dict(G.out_degree())

    R = randomize_degree_preserving(G, seed=123, n_swaps=20)

    # Node count and edge count must be identical
    assert R.number_of_nodes() == G.number_of_nodes()
    assert R.number_of_edges() == G.number_of_edges()

    # Every node's in-degree and out-degree must match exactly
    for node in G.nodes():
        assert R.in_degree(node) == in_degrees_orig[node]
        assert R.out_degree(node) == out_degrees_orig[node]


def test_run_module_2_integration(tmp_path: Path):
    """Test full Module 2 execution and warning generation for correlation data."""
    # Create test regulatory network with mixed evidence
    tsv_content = (
        "source\ttarget\tinteraction\tsign\tconfidence\tevidence\n"
        "X\tY\tactivation\t+\t0.9\texperimental\n"
        "Y\tZ\tactivation\t+\t0.85\texperimental\n"
        "X\tZ\tactivation\t+\t0.95\texperimental\n"
        "W\tZ\tassociation\tunspecified\t0.5\tcorrelation\n"
    )
    tsv_file = tmp_path / "test_net.tsv"
    tsv_file.write_text(tsv_content)

    result = run_module_2(
        input_file=tsv_file,
        output_dir=tmp_path,
        num_null_graphs=10,
        seed=42,
    )

    assert result.run_metadata.status == "success"
    assert result.outputs["node_count"] == 4
    assert result.outputs["edge_count"] == 4
    assert result.outputs["ffl_count"] == 1

    # Ensure correlation warning is raised
    assert any("correlation-only" in w for w in result.warnings)

    # Ensure interpretation record exists and is labeled INTERPRETATION
    assert len(result.interpretations) >= 1
    assert result.interpretations[0].classification == EvidenceClass.INTERPRETATION

    # Ensure output file was written
    out_file = tmp_path / "test_net_network_metrics.json"
    assert out_file.is_file()
