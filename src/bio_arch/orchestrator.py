"""Pipeline Orchestrator.

Coordinates sequential execution of bio_arch modules using shared data contracts.
Features:
- Validates YAML configuration.
- Creates unique run directories with frozen configuration and run logs.
- Deterministic module sub-seeds derived via SeedManager.
- Execution modes: full run, single module, dry-run, resume, validation-only.
- Graceful failure handling: completed artifacts are preserved.
- Final reporting: machine-readable JSON + human-readable Markdown adhering to Section 7.3.
- Pure coordination: strictly keeps scientific logic inside modules.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import time
from typing import Any, Literal
import yaml

from bio_arch.contracts import (
    AnalysisRun,
    DatasetManifest,
    EvidenceClass,
    Finding,
    InterpretationRecord,
    ModuleResult,
    ValidationStatus,
)
from bio_arch.logger import setup_logger
from bio_arch.modules.ingestion import ingest_file
from bio_arch.modules.information import run_module_1
from bio_arch.modules.regulation import run_module_2
from bio_arch.modules.simulation import run_module_3
from bio_arch.provenance import (
    SeedManager,
    compute_sha256,
    get_system_environment,
    now_iso,
)


@dataclass
class PipelineConfig:
    """Validated configuration for pipeline execution."""

    project_name: str
    version: str
    global_seed: int = 42
    dry_run: bool = False
    output_dir: str = "outputs/runs"
    resource_limits: dict[str, Any] = field(default_factory=dict)
    modules: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_yaml(cls, config_path: Path | str) -> PipelineConfig:
        p = Path(config_path)
        if not p.is_file():
            raise FileNotFoundError(f"Configuration file not found: {p}")

        raw = yaml.safe_load(p.read_text(encoding="utf-8"))
        if not isinstance(raw, dict):
            raise ValueError(f"Invalid YAML config structure in {p}")

        return cls(
            project_name=raw.get("project_name", "bio_information_architecture"),
            version=raw.get("version", "0.1.0"),
            global_seed=raw.get("global_seed", 42),
            dry_run=raw.get("dry_run", False),
            output_dir=raw.get("output_dir", "outputs/runs"),
            resource_limits=raw.get("resource_limits", {}),
            modules=raw.get("modules", {}),
        )


class Orchestrator:
    """Coordinates pipeline execution, artifact handoff, and report generation."""

    def __init__(
        self,
        config: PipelineConfig,
        run_dir: Path | None = None,
        resume: bool = False,
    ) -> None:
        self.config = config
        self.resume = resume
        self.seed_mgr = SeedManager(master_seed=config.global_seed)

        # Setup run directory
        if run_dir:
            self.run_dir = Path(run_dir)
        else:
            timestamp_str = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
            self.run_dir = Path(config.output_dir) / f"run_{timestamp_str}_{config.project_name}"

        self.run_dir.mkdir(parents=True, exist_ok=True)

        # Structured logger targeting both console and run directory
        self.logger = setup_logger(
            name="bio_arch.orchestrator",
            log_file=self.run_dir / "execution.log",
        )

        # Freeze config snapshot for reproducibility
        self.frozen_config_path = self.run_dir / "frozen_config.yaml"
        if not self.frozen_config_path.exists():
            self.frozen_config_path.write_text(yaml.dump(asdict(self.config)), encoding="utf-8")

        self.manifests: list[DatasetManifest] = []
        self.module_results: dict[str, ModuleResult] = {}
        self.timing_records: dict[str, float] = {}

    def run(
        self,
        mode: Literal["full", "single", "dry-run", "resume", "validate-only"] = "full",
        target_module: str | None = None,
    ) -> dict[str, Any]:
        """Execute the pipeline in the specified mode."""
        self.logger.info("Starting orchestrator in mode '%s' (Run dir: %s)", mode, self.run_dir)

        if mode == "dry-run" or self.config.dry_run:
            return self._execute_dry_run()

        modules_to_run = self._determine_execution_queue(mode, target_module)
        self.logger.info("Execution queue: %s", ", ".join(modules_to_run))

        try:
            # 1. Module 0: Ingestion
            if "module_0_ingestion" in modules_to_run:
                self._run_module_0()

            # 2. Module 1: Sequence Information
            if "module_1_information" in modules_to_run:
                self._run_module_1()

            # 3. Module 2: Regulatory Networks
            if "module_2_regulation" in modules_to_run:
                self._run_module_2()

            # 4. Module 3: Simulation
            if "module_3_simulation" in modules_to_run:
                self._run_module_3()

        except Exception as exc:
            self.logger.error("Pipeline failure: %s", exc, exc_info=True)
            # Ensure partial results and reports are saved gracefully
            self._write_reports(status="failed", error_message=str(exc))
            raise

        # Generate final outputs
        summary = self._write_reports(status="completed")
        self.logger.info("Pipeline completed successfully.")
        return summary

    def _determine_execution_queue(
        self,
        mode: str,
        target_module: str | None,
    ) -> list[str]:
        """Determine which modules should be executed based on mode and resume status."""
        if mode == "validate-only":
            return ["module_0_ingestion"]

        if mode == "single":
            if not target_module:
                raise ValueError("Single module mode requires a target module name.")
            return [target_module]

        all_modules = [
            "module_0_ingestion",
            "module_1_information",
            "module_2_regulation",
            "module_3_simulation",
        ]

        active = [m for m in all_modules if self.config.modules.get(m, {}).get("enabled", True)]

        if self.resume:
            # Skip modules whose result files already exist
            filtered = []
            for m in active:
                res_file = self.run_dir / f"{m}_result.json"
                if res_file.exists():
                    self.logger.info("Skipping already completed module '%s' (resume mode)", m)
                else:
                    filtered.append(m)
            return filtered

        return active

    def _execute_dry_run(self) -> dict[str, Any]:
        """Validate input existence and parameters without expensive computation."""
        self.logger.info("Performing DRY-RUN validation...")
        issues = []

        mod0_cfg = self.config.modules.get("module_0_ingestion", {})
        for item in mod0_cfg.get("input_files", []):
            p = Path(item.get("path", ""))
            if not p.is_file():
                issues.append(f"Input file not found: '{p}'")

        status = "dry_run_passed" if not issues else "dry_run_failed"
        self.logger.info("Dry run status: %s (Issues: %s)", status, issues)
        return {"status": status, "issues": issues}

    def _run_module_0(self) -> None:
        """Execute Module 0 Ingestion."""
        start_t = time.perf_counter()
        self.logger.info("Executing Module 0 (Data Ingestion & Validation)...")
        mod_dir = self.run_dir / "module_0_ingestion"
        mod_dir.mkdir(parents=True, exist_ok=True)

        cfg = self.config.modules.get("module_0_ingestion", {})
        input_files = cfg.get("input_files", [])

        if not input_files:
            self.logger.warning("No input files configured for Module 0. Skipping.")
            return

        for idx, item in enumerate(input_files):
            p = Path(item["path"])
            dataset_id = item.get("dataset_id", f"dataset_{p.stem}")
            organism = item.get("organism", "Unknown")
            lic = item.get("license", "Unspecified")
            src = item.get("source", "Local Ingestion")
            seq_type = item.get("sequence_type", "DNA")

            manifest, report, result = ingest_file(
                file_path=p,
                dataset_id=dataset_id,
                organism=organism,
                license=lic,
                source=src,
                sequence_type=seq_type,
                output_dir=mod_dir,
                seed=self.seed_mgr.derive_seed(f"ingest_{dataset_id}"),
            )
            self.manifests.append(manifest)
            self.module_results[f"module_0_{dataset_id}"] = result

        # Save module result checkpoint
        (self.run_dir / "module_0_ingestion_result.json").write_text(
            json.dumps([r.to_dict() for r in self.module_results.values() if r.run_metadata.module == "module_0_ingestion"], indent=2)
        )
        self.timing_records["module_0_ingestion"] = time.perf_counter() - start_t

    def _run_module_1(self) -> None:
        """Execute Module 1 Sequence Information Architecture."""
        start_t = time.perf_counter()
        self.logger.info("Executing Module 1 (DNA/RNA Information Architecture)...")
        mod_dir = self.run_dir / "module_1_information"
        mod_dir.mkdir(parents=True, exist_ok=True)

        cfg = self.config.modules.get("module_1_information", {})
        shuffles = cfg.get("num_shuffles", 100)

        # Look for normalized outputs from module 0 or fall back to sample fasta
        norm_files = list((self.run_dir / "module_0_ingestion").glob("*_normalized.json"))
        target_inputs = norm_files or [Path("data/synthetic_sample.fasta")]

        for in_file in target_inputs:
            if not in_file.is_file():
                continue
            res = run_module_1(
                input_file=in_file,
                output_dir=mod_dir,
                num_shuffles=shuffles,
                seed=self.seed_mgr.derive_seed(f"mod1_{in_file.stem}"),
            )
            self.module_results[f"module_1_{in_file.stem}"] = res

        (self.run_dir / "module_1_information_result.json").write_text(
            json.dumps([r.to_dict() for r in self.module_results.values() if r.run_metadata.module == "module_1_information"], indent=2)
        )
        self.timing_records["module_1_information"] = time.perf_counter() - start_t

    def _run_module_2(self) -> None:
        """Execute Module 2 Regulatory Networks."""
        start_t = time.perf_counter()
        self.logger.info("Executing Module 2 (Regulatory Network Discovery)...")
        mod_dir = self.run_dir / "module_2_regulation"
        mod_dir.mkdir(parents=True, exist_ok=True)

        cfg = self.config.modules.get("module_2_regulation", {})
        null_graphs = cfg.get("num_null_graphs", 50)
        edge_file = Path(cfg.get("input_file", "tests/fixtures/valid_network.tsv"))

        if edge_file.is_file():
            res = run_module_2(
                input_file=edge_file,
                output_dir=mod_dir,
                num_null_graphs=null_graphs,
                seed=self.seed_mgr.derive_seed("mod2_regulation"),
            )
            self.module_results["module_2_regulation"] = res
            (self.run_dir / "module_2_regulation_result.json").write_text(res.to_json())

        self.timing_records["module_2_regulation"] = time.perf_counter() - start_t

    def _run_module_3(self) -> None:
        """Execute Module 3 Simulation."""
        start_t = time.perf_counter()
        self.logger.info("Executing Module 3 (Emergent Behavior & Simulation)...")
        mod_dir = self.run_dir / "module_3_simulation"
        mod_dir.mkdir(parents=True, exist_ok=True)

        cfg = self.config.modules.get("module_3_simulation", {})
        steps = cfg.get("time_steps", 30)
        trials = cfg.get("trials", 10)
        edge_file = Path(cfg.get("input_file", "tests/fixtures/valid_network.tsv"))

        if edge_file.is_file():
            res = run_module_3(
                input_file=edge_file,
                output_dir=mod_dir,
                time_steps=steps,
                trials=trials,
                seed=self.seed_mgr.derive_seed("mod3_simulation"),
            )
            self.module_results["module_3_simulation"] = res
            (self.run_dir / "module_3_simulation_result.json").write_text(res.to_json())

        self.timing_records["module_3_simulation"] = time.perf_counter() - start_t

    def _write_reports(self, status: str, error_message: str | None = None) -> dict[str, Any]:
        """Generate final machine-readable JSON and human-readable Markdown reports."""
        all_findings: list[Finding] = []
        all_interpretations: list[InterpretationRecord] = []
        all_warnings: list[str] = []

        for res in self.module_results.values():
            all_findings.extend(res.findings)
            all_interpretations.extend(res.interpretations)
            all_warnings.extend(res.warnings)

        summary_payload = {
            "run_id": self.run_dir.name,
            "status": status,
            "error_message": error_message,
            "timestamp": now_iso(),
            "global_seed": self.config.global_seed,
            "timings_seconds": self.timing_records,
            "manifests": [m.to_dict() for m in self.manifests],
            "total_findings": len(all_findings),
            "total_interpretations": len(all_interpretations),
            "warnings": all_warnings,
        }

        # 1. JSON Report
        json_report_path = self.run_dir / "final_summary.json"
        json_report_path.write_text(json.dumps(summary_payload, indent=2), encoding="utf-8")

        # 2. Markdown Report conforming to Section 7.3 Final Report Structure
        md_lines = [
            f"# Scientific Analysis Report: {self.config.project_name}",
            f"**Run ID**: `{self.run_dir.name}` | **Status**: `{status.upper()}` | **Date**: {now_iso()}\n",
            "## 1. Run Summary",
            f"- **Execution Status**: {status}",
            f"- **Global Seed**: {self.config.global_seed}",
            f"- **Modules Executed**: {', '.join(self.timing_records.keys()) or 'None'}",
            f"- **Total Duration**: {sum(self.timing_records.values()):.2f}s",
            "",
            "## 2. Datasets and Provenance",
        ]

        if self.manifests:
            for m in self.manifests:
                md_lines.append(f"- **{m.dataset_id}** ({m.organism})")
                md_lines.append(f"  - Checksum (SHA-256): `{m.checksum}`")
                md_lines.append(f"  - Source: {m.source} | License: {m.license}")
                md_lines.append(f"  - Validation: `{m.validation_status.value}`")
        else:
            md_lines.append("No independent datasets registered in this run.")

        md_lines.extend([
            "",
            "## 3. Methods and Parameters",
            f"- **Config snapshot**: Saved at `frozen_config.yaml`",
            f"- **Resource Limits**: {json.dumps(self.config.resource_limits)}",
            "",
            "## 4. Empirical Measurements",
        ])

        measurements = [f for f in all_findings if f.control_distribution is None]
        for m in measurements:
            md_lines.append(f"- **{m.metric}**: `{m.observed_value}` — *{m.biological_context}*")

        md_lines.extend([
            "",
            "## 5. Control Comparisons (Null Models)",
        ])

        controlled = [f for f in all_findings if f.control_distribution is not None]
        for c in controlled:
            z_str = f"z={c.effect_size}" if c.effect_size is not None else "N/A"
            p_str = f"adj-p={c.adjusted_p_value}" if c.adjusted_p_value is not None else "N/A"
            md_lines.append(f"- **{c.metric}**: observed={c.observed_value} ({z_str}, {p_str})")
            md_lines.append(f"  - Context: {c.biological_context}")
            md_lines.append(f"  - Null: {c.control_distribution.get('null_model')} (N={c.control_distribution.get('iterations')})")

        md_lines.extend([
            "",
            "## 6. Simulations",
        ])

        sim_interps = [i for i in all_interpretations if i.classification == EvidenceClass.SIMULATION]
        if sim_interps:
            for s in sim_interps:
                md_lines.append(f"- {s.claim}")
                md_lines.append(f"  - *Limitations*: {', '.join(s.limitations)}")
        else:
            md_lines.append("No active simulations recorded.")

        md_lines.extend([
            "",
            "## 7. Interpretations and Alternatives",
        ])

        norm_interps = [i for i in all_interpretations if i.classification == EvidenceClass.INTERPRETATION]
        for i in norm_interps:
            md_lines.append(f"- **Claim**: {i.claim}")
            if i.alternatives:
                md_lines.append(f"  - **Alternative Explanations**: {', '.join(i.alternatives)}")
            if i.limitations:
                md_lines.append(f"  - **Limitations**: {', '.join(i.limitations)}")

        md_lines.extend([
            "",
            "## 8. Testable Hypotheses",
        ])
        for i in all_interpretations:
            if i.proposed_test:
                md_lines.append(f"- Derived from `{', '.join(i.finding_ids)}`: {i.proposed_test}")

        md_lines.extend([
            "",
            "## 9. Limitations & Guardrails",
            "- Statistical patterns and computational models must not be treated as literal biological mechanisms.",
            "- Correlation does not imply regulatory or causal direction.",
            "- Measurements are restricted to documented datasets and parameter bounds.",
            "",
            "## 10. Next Experiments",
            "- Expand testing to diverse homologous regions across related prokaryotic and eukaryotic clades.",
            "- Benchmark Boolean state attractor predictions against single-cell dynamic assays.",
        ])

        md_report_path = self.run_dir / "final_report.md"
        md_report_path.write_text("\n".join(md_lines), encoding="utf-8")

        return summary_payload


def main() -> None:
    """CLI interface for the bio_arch orchestrator."""
    parser = argparse.ArgumentParser(
        description="Bio-Arch Pipeline Orchestrator",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--config", type=str, default="config/example.yaml", help="Path to YAML config")
    parser.add_argument(
        "--mode",
        type=str,
        choices=["full", "single", "dry-run", "resume", "validate-only"],
        default="full",
        help="Execution mode",
    )
    parser.add_argument("--module", type=str, default=None, help="Target module for 'single' mode")
    parser.add_argument("--run-dir", type=str, default=None, help="Custom run directory or resume target")

    args = parser.parse_args()

    config = PipelineConfig.from_yaml(args.config)
    orchestrator = Orchestrator(
        config=config,
        run_dir=Path(args.run_dir) if args.run_dir else None,
        resume=(args.mode == "resume"),
    )

    summary = orchestrator.run(mode=args.mode, target_module=args.module)
    print("\n--- Pipeline Run Summary ---")
    print(f"Status:          {summary['status']}")
    if "run_id" in summary:
        print(f"Run Directory:   {summary['run_id']}")
        print(f"Total Findings:  {summary['total_findings']}")
        print(f"Interpretations: {summary['total_interpretations']}")
    print("----------------------------\n")


if __name__ == "__main__":
    main()
