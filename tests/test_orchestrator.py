"""Unit and integration tests for Pipeline Orchestrator."""

import json
from pathlib import Path
import pytest
import yaml

from bio_arch.orchestrator import Orchestrator, PipelineConfig


def test_pipeline_config_from_yaml(tmp_path: Path):
    """Test loading and validating PipelineConfig from YAML."""
    yaml_file = tmp_path / "test_config.yaml"
    yaml_file.write_text(
        "project_name: 'test_project'\n"
        "version: '0.1.0'\n"
        "global_seed: 99\n"
        "dry_run: false\n"
        "output_dir: 'outputs/test_runs'\n"
    )

    cfg = PipelineConfig.from_yaml(yaml_file)
    assert cfg.project_name == "test_project"
    assert cfg.global_seed == 99
    assert cfg.dry_run is False


def test_orchestrator_dry_run(tmp_path: Path):
    """Test dry-run mode without executing computation."""
    cfg = PipelineConfig(
        project_name="dry_test",
        version="0.1.0",
        output_dir=str(tmp_path),
        modules={
            "module_0_ingestion": {
                "input_files": [{"path": str(tmp_path / "nonexistent.fasta")}]
            }
        },
    )

    orch = Orchestrator(config=cfg, run_dir=tmp_path / "run_dry")
    res = orch.run(mode="dry-run")

    assert res["status"] == "dry_run_failed"
    assert len(res["issues"]) == 1
    assert "not found" in res["issues"][0]


def test_orchestrator_full_e2e_run(tmp_path: Path):
    """Test end-to-end execution of all modules with synthetic data."""
    # Create tiny synthetic fasta
    fasta = tmp_path / "sample.fasta"
    fasta.write_text(">seq1\nATGCATGCATGCATGCATGCATGCATGCATGC\n")

    # Create tiny synthetic network
    net = tmp_path / "sample_net.tsv"
    net.write_text("source\ttarget\tinteraction\tsign\nA\tB\tactivation\t+\nB\tC\trepression\t-\n")

    cfg = PipelineConfig(
        project_name="e2e_test",
        version="0.1.0",
        global_seed=42,
        output_dir=str(tmp_path / "runs"),
        modules={
            "module_0_ingestion": {
                "enabled": True,
                "input_files": [{
                    "path": str(fasta),
                    "dataset_id": "ds_synth",
                    "organism": "Synthetic",
                    "license": "CC0",
                    "source": "E2E Test",
                }],
            },
            "module_1_information": {
                "enabled": True,
                "num_shuffles": 10,
            },
            "module_2_regulation": {
                "enabled": True,
                "input_file": str(net),
                "num_null_graphs": 10,
            },
            "module_3_simulation": {
                "enabled": True,
                "input_file": str(net),
                "time_steps": 10,
                "trials": 5,
            },
        },
    )

    run_dir = tmp_path / "runs" / "test_run"
    orch = Orchestrator(config=cfg, run_dir=run_dir)
    summary = orch.run(mode="full")

    assert summary["status"] == "completed"
    assert summary["total_findings"] > 0
    assert summary["total_interpretations"] > 0

    # Verify generated artifacts
    assert (run_dir / "frozen_config.yaml").is_file()
    assert (run_dir / "execution.log").is_file()
    assert (run_dir / "final_summary.json").is_file()
    assert (run_dir / "final_report.md").is_file()

    # Verify Markdown report contains the 10 required sections from Section 7.3
    report_text = (run_dir / "final_report.md").read_text(encoding="utf-8")
    assert "## 1. Run Summary" in report_text
    assert "## 2. Datasets and Provenance" in report_text
    assert "## 3. Methods and Parameters" in report_text
    assert "## 4. Empirical Measurements" in report_text
    assert "## 5. Control Comparisons (Null Models)" in report_text
    assert "## 6. Simulations" in report_text
    assert "## 7. Interpretations and Alternatives" in report_text
    assert "## 8. Testable Hypotheses" in report_text
    assert "## 9. Limitations & Guardrails" in report_text
    assert "## 10. Next Experiments" in report_text


def test_orchestrator_resume_mode(tmp_path: Path):
    """Verify resume mode skips existing completed modules."""
    run_dir = tmp_path / "resumable_run"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Pre-populate a completed module 2 result
    dummy_result = run_dir / "module_2_regulation_result.json"
    dummy_result.write_text('{"status": "pre_computed"}')

    cfg = PipelineConfig(
        project_name="resume_test",
        version="0.1.0",
        modules={
            "module_2_regulation": {"enabled": True},
        },
    )

    orch = Orchestrator(config=cfg, run_dir=run_dir, resume=True)
    queue = orch._determine_execution_queue(mode="resume", target_module=None)

    # Module 2 must be filtered out since result file already exists
    assert "module_2_regulation" not in queue
