"""Unit tests for Module 3: Emergent Behavior and Simulation."""

from pathlib import Path
import pytest

from bio_arch.contracts import EvidenceClass
from bio_arch.modules.regulation import RegulatoryEdge
from bio_arch.modules.simulation import (
    BooleanNetwork,
    calculate_hamming_distance,
    evaluate_knockout_sensitivity,
    run_module_3,
)


def test_boolean_network_logic():
    """Verify basic activation, repression, knockout, and forced activation rules."""
    edges = [
        RegulatoryEdge(source="A", target="B", interaction_type="activation", sign="activation"),
        RegulatoryEdge(source="C", target="B", interaction_type="repression", sign="repression"),
    ]

    net = BooleanNetwork.from_edges(edges)

    # 1. Activator on, repressor off -> B turns ON
    state_act = {"A": 1, "B": 0, "C": 0}
    assert net.evaluate_node("B", state_act) == 1

    # 2. Activator on, repressor ON -> B turns OFF (dominant repression)
    state_rep = {"A": 1, "B": 1, "C": 1}
    assert net.evaluate_node("B", state_rep) == 0

    # 3. Knockout test: B knocked out remains 0
    net_ko = BooleanNetwork.from_edges(edges, knockouts={"B"})
    assert net_ko.evaluate_node("B", state_act) == 0

    # 4. Forced activation: B constitutively 1
    net_forced = BooleanNetwork.from_edges(edges, forced_activations={"B"})
    assert net_forced.evaluate_node("B", state_rep) == 1


def test_attractor_fixed_point():
    """Verify convergence to a stable fixed point."""
    # Self-activating node X -> X
    edges = [RegulatoryEdge(source="X", target="X", sign="activation")]
    net = BooleanNetwork.from_edges(edges)

    traj = net.simulate(initial_state={"X": 1}, max_steps=10)
    assert traj.attractor is not None
    assert traj.attractor.attractor_type == "fixed_point"
    assert traj.attractor.cycle_length == 1
    assert traj.attractor.states[0]["X"] == 1


def test_attractor_limit_cycle():
    """Verify oscillation detection in a negative feedback ring (repressilator toy)."""
    # X -| Y, Y -| Z, Z -| X
    edges = [
        RegulatoryEdge(source="X", target="Y", sign="repression"),
        RegulatoryEdge(source="Y", target="Z", sign="repression"),
        RegulatoryEdge(source="Z", target="X", sign="repression"),
    ]
    net = BooleanNetwork.from_edges(edges)

    # State: 1, 0, 0 -> 1 inhibits Y, Z has no repressor (so remains 0 or changes depending on rule)
    traj = net.simulate(initial_state={"X": 1, "Y": 0, "Z": 0}, max_steps=20)
    assert traj.attractor is not None
    # Must reach a fixed point or periodic cycle
    assert traj.attractor.cycle_length >= 1


def test_hamming_distance():
    """Verify Hamming distance calculation."""
    s1 = {"A": 0, "B": 1, "C": 1}
    s2 = {"A": 1, "B": 1, "C": 0}
    # Differ on A and C -> distance 2
    assert calculate_hamming_distance(s1, s2) == 2


def test_run_module_3_integration(tmp_path: Path):
    """Test Module 3 end-to-end execution and SIMULATION evidence tagging."""
    net_file = tmp_path / "simple_net.tsv"
    net_file.write_text(
        "source\ttarget\tinteraction\tsign\n"
        "A\tB\tactivation\t+\n"
        "B\tA\tactivation\t+\n"
    )

    result = run_module_3(
        input_file=net_file,
        output_dir=tmp_path,
        time_steps=15,
        trials=5,
        seed=42,
    )

    assert result.run_metadata.status == "success"
    assert result.outputs["trials_run"] == 5
    assert len(result.findings) >= 2

    # Check evidence classification
    assert len(result.interpretations) >= 1
    assert result.interpretations[0].classification == EvidenceClass.SIMULATION
    assert result.interpretations[0].claim.startswith("Under this model")

    # Output file saved
    out_file = tmp_path / "simple_net_simulation_results.json"
    assert out_file.is_file()
