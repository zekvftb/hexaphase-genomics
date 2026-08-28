"""Module 3: Emergent Behavior and Simulation.

A transparent, discrete-state Boolean regulatory network simulator.
Features:
- Synchronous and asynchronous state updates.
- Node state definitions (0: inactive, 1: active).
- Threshold and signed regulatory logic rules (activators vs repressors).
- Perturbation and sensitivity: single-gene knockouts, forced constitutive activation, stochastic noise.
- Deterministic attractor detection (fixed points and limit cycles).
- Robustness evaluation (Hamming distance across trajectory perturbations).
- Strict scientific guardrail: results are clearly labeled SIMULATION ('Under this model...').
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
import random
import sys
from typing import Any, Callable, Literal

from bio_arch.contracts import (
    AnalysisRun,
    EvidenceClass,
    Finding,
    InterpretationRecord,
    ModuleResult,
)
from bio_arch.logger import setup_logger
from bio_arch.modules.regulation import RegulatoryEdge, parse_regulatory_file
from bio_arch.provenance import (
    SeedManager,
    get_system_environment,
    now_iso,
)

logger = setup_logger("bio_arch.simulation")

UpdateMode = Literal["synchronous", "asynchronous"]


@dataclass
class Attractor:
    """A dynamic attractor (fixed point or periodic limit cycle) in state space."""

    attractor_type: Literal["fixed_point", "limit_cycle"]
    cycle_length: int
    states: list[dict[str, int]]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class SimulationTrajectory:
    """Full time-course record of a simulation trial."""

    time_steps: int
    states: list[dict[str, int]]  # states[t] = {node_name: 0/1}
    attractor: Attractor | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "time_steps": self.time_steps,
            "states": self.states,
            "attractor": self.attractor.to_dict() if self.attractor else None,
        }


class BooleanNetwork:
    """A discrete Boolean network driven by regulatory interactions."""

    def __init__(
        self,
        nodes: list[str],
        activators: dict[str, list[str]],
        repressors: dict[str, list[str]],
        knockouts: set[str] | None = None,
        forced_activations: set[str] | None = None,
    ) -> None:
        self.nodes = sorted(nodes)
        self.activators = activators
        self.repressors = repressors
        self.knockouts = knockouts or set()
        self.forced_activations = forced_activations or set()

    @classmethod
    def from_edges(
        cls,
        edges: list[RegulatoryEdge],
        knockouts: set[str] | None = None,
        forced_activations: set[str] | None = None,
    ) -> BooleanNetwork:
        """Construct network from RegulatoryEdge list."""
        all_nodes = set()
        activators: dict[str, list[str]] = defaultdict(list)
        repressors: dict[str, list[str]] = defaultdict(list)

        for e in edges:
            all_nodes.add(e.source)
            all_nodes.add(e.target)
            if e.sign == "repression":
                repressors[e.target].append(e.source)
            else:
                activators[e.target].append(e.source)

        return cls(
            nodes=list(all_nodes),
            activators=dict(activators),
            repressors=dict(repressors),
            knockouts=knockouts,
            forced_activations=forced_activations,
        )

    def evaluate_node(self, node: str, current_state: dict[str, int]) -> int:
        """Compute the next state of a node according to standard biological logic:

        - If knocked out -> 0
        - If constitutively forced -> 1
        - If any active repressor is present -> 0 (dominant inhibition)
        - Otherwise, active if any activator is active, or if node has no inputs retain self.
        """
        if node in self.knockouts:
            return 0
        if node in self.forced_activations:
            return 1

        active_repressors = any(current_state.get(rep, 0) == 1 for rep in self.repressors.get(node, []))
        if active_repressors:
            return 0

        node_activators = self.activators.get(node, [])
        if node_activators:
            return 1 if any(current_state.get(act, 0) == 1 for act in node_activators) else 0

        # Self-maintaining if no regulatory inputs
        return current_state.get(node, 0)

    def step(
        self,
        current_state: dict[str, int],
        mode: UpdateMode = "synchronous",
        noise: float = 0.0,
        rng: random.Random | None = None,
    ) -> dict[str, int]:
        """Advance the network by one time step."""
        next_state = dict(current_state)

        if mode == "synchronous":
            for node in self.nodes:
                next_state[node] = self.evaluate_node(node, current_state)
        else:
            # Asynchronous update: random node ordering
            shuffled_nodes = list(self.nodes)
            if rng:
                rng.shuffle(shuffled_nodes)
            for node in shuffled_nodes:
                next_state[node] = self.evaluate_node(node, next_state)

        # Inject noise if requested
        if noise > 0.0 and rng:
            for node in self.nodes:
                if rng.random() < noise:
                    next_state[node] = 1 - next_state[node]

        return next_state

    def state_to_key(self, state: dict[str, int]) -> tuple[int, ...]:
        """Convert state dictionary to hashable tuple for fast attractor detection."""
        return tuple(state[node] for node in self.nodes)

    def simulate(
        self,
        initial_state: dict[str, int],
        max_steps: int = 50,
        mode: UpdateMode = "synchronous",
        noise: float = 0.0,
        rng: random.Random | None = None,
    ) -> SimulationTrajectory:
        """Run trajectory from initial state and identify attractors."""
        history = [dict(initial_state)]
        seen_states: dict[tuple[int, ...], int] = {self.state_to_key(initial_state): 0}
        attractor = None

        current = dict(initial_state)
        for t in range(1, max_steps + 1):
            next_s = self.step(current, mode=mode, noise=noise, rng=rng)
            history.append(dict(next_s))

            key = self.state_to_key(next_s)
            if noise == 0.0 and mode == "synchronous" and key in seen_states:
                # Cycle found
                start_step = seen_states[key]
                cycle_states = history[start_step : len(history) - 1]
                cycle_len = len(cycle_states)
                attractor = Attractor(
                    attractor_type="fixed_point" if cycle_len == 1 else "limit_cycle",
                    cycle_length=cycle_len,
                    states=cycle_states,
                )
                break

            seen_states[key] = t
            current = next_s

        return SimulationTrajectory(
            time_steps=len(history) - 1,
            states=history,
            attractor=attractor,
        )


# ---------------------------------------------------------------------------
# Sensitivity and Perturbation Analysis
# ---------------------------------------------------------------------------

def calculate_hamming_distance(s1: dict[str, int], s2: dict[str, int]) -> int:
    """Number of state differences between two states."""
    return sum(1 for k in s1 if s1[k] != s2.get(k, 0))


def evaluate_knockout_sensitivity(
    edges: list[RegulatoryEdge],
    initial_state: dict[str, int],
    steps: int = 30,
    seed: int = 42,
) -> dict[str, Any]:
    """Simulate single-gene knockouts and measure trajectory deviation from baseline."""
    rng = random.Random(seed)
    net_baseline = BooleanNetwork.from_edges(edges)
    traj_baseline = net_baseline.simulate(initial_state, max_steps=steps, rng=rng)
    baseline_final = traj_baseline.states[-1]

    knockout_impacts: dict[str, int] = {}
    for node in net_baseline.nodes:
        net_ko = BooleanNetwork.from_edges(edges, knockouts={node})
        traj_ko = net_ko.simulate(initial_state, max_steps=steps, rng=rng)
        ko_final = traj_ko.states[-1]
        dist = calculate_hamming_distance(baseline_final, ko_final)
        knockout_impacts[node] = dist

    return {
        "baseline_attractor": traj_baseline.attractor.to_dict() if traj_baseline.attractor else None,
        "knockout_hamming_distances": knockout_impacts,
    }


# ---------------------------------------------------------------------------
# Module 3 Execution & Contracts
# ---------------------------------------------------------------------------

def run_module_3(
    input_file: Path | str,
    output_dir: Path | str | None = None,
    time_steps: int = 30,
    trials: int = 10,
    seed: int = 42,
) -> ModuleResult:
    """Execute Module 3 discrete simulation from an edge list file."""
    in_path = Path(input_file)
    edges = parse_regulatory_file(in_path)

    seed_mgr = SeedManager(master_seed=seed)
    mod_seed = seed_mgr.derive_seed("module_3_simulation")
    rng = random.Random(mod_seed)

    run_meta = AnalysisRun(
        run_id=f"run_mod3_{in_path.stem}_{mod_seed}",
        timestamp=now_iso(),
        module="module_3_simulation",
        version="0.1.0",
        input_ids=[in_path.name],
        parameters={"time_steps": time_steps, "trials": trials, "seed": mod_seed},
        seed=mod_seed,
        environment=get_system_environment(),
        status="success",
    )

    net = BooleanNetwork.from_edges(edges)

    # Run multi-trial attractor discovery from random initial states
    discovered_attractors: list[dict[str, Any]] = []
    fixed_points = 0
    limit_cycles = 0

    for trial_idx in range(trials):
        init_state = {node: rng.choice([0, 1]) for node in net.nodes}
        traj = net.simulate(init_state, max_steps=time_steps, mode="synchronous", rng=rng)
        if traj.attractor:
            att_dict = traj.attractor.to_dict()
            if traj.attractor.attractor_type == "fixed_point":
                fixed_points += 1
            else:
                limit_cycles += 1
            if att_dict not in discovered_attractors:
                discovered_attractors.append(att_dict)

    # Run knockout sensitivity
    default_init = {node: 1 for node in net.nodes}
    ko_sensitivity = evaluate_knockout_sensitivity(edges, default_init, steps=time_steps, seed=mod_seed)

    # Findings with SIMULATION evidence classification
    findings: list[Finding] = [
        Finding(
            finding_id=f"f_{in_path.stem}_attractor_count",
            metric="unique_attractors_found",
            observed_value=len(discovered_attractors),
            control_distribution={"trials": trials, "fixed_points": fixed_points, "limit_cycles": limit_cycles},
            biological_context=f"Attractor landscape under synchronous Boolean rules for {in_path.stem}.",
        ),
        Finding(
            finding_id=f"f_{in_path.stem}_knockout_impact",
            metric="max_knockout_hamming_distance",
            observed_value=max(ko_sensitivity["knockout_hamming_distances"].values(), default=0),
            biological_context="Maximum state shift resulting from single-node in-silico knockout.",
        ),
    ]

    # Explicit simulation claim
    interpretations: list[InterpretationRecord] = [
        InterpretationRecord(
            finding_ids=[f"f_{in_path.stem}_attractor_count"],
            classification=EvidenceClass.SIMULATION,
            claim=f"{EvidenceClass.SIMULATION.required_prefix}, the network converges to {len(discovered_attractors)} distinct attractor states from {trials} random starts.",
            alternatives=["Continuous kinetic models (ODEs) or stochastic chemical master equations may exhibit different multistability."],
            limitations=["Boolean assumptions discretize expression levels into binary 0/1 states."],
            proposed_test="Compare predicted attractors with single-cell RNA-seq cluster expression profiles.",
        )
    ]

    outputs = {
        "nodes": net.nodes,
        "trials_run": trials,
        "unique_attractors": len(discovered_attractors),
        "fixed_points": fixed_points,
        "limit_cycles": limit_cycles,
        "knockout_sensitivity": ko_sensitivity["knockout_hamming_distances"],
    }

    artifact_paths: list[str] = []
    if output_dir:
        out_p = Path(output_dir)
        out_p.mkdir(parents=True, exist_ok=True)
        out_file = out_p / f"{in_path.stem}_simulation_results.json"

        result_payload = {
            "run_metadata": run_meta.to_dict(),
            "outputs": outputs,
            "findings": [f.to_dict() for f in findings],
            "interpretations": [i.to_dict() for i in interpretations],
            "attractors": discovered_attractors,
        }
        out_file.write_text(json.dumps(result_payload, indent=2), encoding="utf-8")
        artifact_paths.append(str(out_file))

    return ModuleResult(
        run_metadata=run_meta,
        outputs=outputs,
        findings=findings,
        interpretations=interpretations,
        warnings=[],
        errors=[],
        artifact_paths=artifact_paths,
    )


def main() -> None:
    """CLI entrypoint for Module 3."""
    parser = argparse.ArgumentParser(
        description="Bio-Arch Module 3: Emergent Behavior & Simulation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("file", type=str, help="Path to TSV/CSV regulatory edge list")
    parser.add_argument("--steps", type=int, default=30, help="Maximum simulation time steps")
    parser.add_argument("--trials", type=int, default=20, help="Number of random initial trials")
    parser.add_argument("--seed", type=int, default=42, help="Deterministic random seed")
    parser.add_argument("--outdir", type=str, default="outputs/simulation", help="Output directory")

    args = parser.parse_args()

    result = run_module_3(
        input_file=args.file,
        output_dir=args.outdir,
        time_steps=args.steps,
        trials=args.trials,
        seed=args.seed,
    )

    outputs = result.outputs
    print("\n--- Module 3: Simulation Results ---")
    print(f"Nodes:               {len(outputs['nodes'])}")
    print(f"Trials:              {outputs['trials_run']}")
    print(f"Unique Attractors:   {outputs['unique_attractors']}")
    print(f"Fixed Points:        {outputs['fixed_points']}")
    print(f"Limit Cycles:        {outputs['limit_cycles']}")
    print(f"Knockout Impacts:    {outputs['knockout_sensitivity']}")
    if result.artifact_paths:
        print(f"Saved artifacts to:  {result.artifact_paths[0]}")
    print("------------------------------------\n")


if __name__ == "__main__":
    main()
