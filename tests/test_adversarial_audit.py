"""Automated Verification Suite for Zero-Trust Adversarial Falsification Harness."""

from __future__ import annotations

import json
from pathlib import Path
import pytest

from scripts.validation.adversarial_audit import (
    evaluate_falsification_ledger,
    run_adversarial_negative_controls,
)


def test_adversarial_negative_controls_suppression():
    """Verify that negative controls produce 0 false-positive validated smORFs."""
    res = run_adversarial_negative_controls(n_shuffles=20, seed=123)

    assert res["false_positives_detected"] is False
    assert res["pipeline_integrity"] == "UNCOMPROMISED (Zero False Positives)"
    assert res["controls_tested"] == 3

    for ctrl in res["control_details"]:
        assert ctrl["falsely_validated_full_pipeline"] == 0
        assert "PASSED" in ctrl["status"]


def test_falsification_ledger_partitioning():
    """Verify that candidates are partitioned into Table A (Validated) and Table B (Refuted)."""
    root = Path(__file__).parent.parent
    audit = evaluate_falsification_ledger(root)

    assert audit["zero_trust_status"] == "VERIFIED_ROBUST"
    assert audit["validated_candidates_count"] > 0
    assert audit["refuted_candidates_count"] > 0

    # Validated candidates must satisfy all 3 criteria
    for v in audit["validated_candidates"]:
        assert v["z_score"] >= 3.0
        assert v["tm_core_plddt"] >= 70.0
        assert v["ribo_tpi"] >= 0.60
        assert "VALIDATED" in v["status"]

    # Refuted candidates must have an explicit failure reason
    for r in audit["refuted_candidates"]:
        assert "failure_reasons" in r
        assert len(r["failure_reasons"]) > 0
        assert "misleading_heuristic" in r


def test_audit_artifacts_exist():
    """Verify that audit json and markdown falsification ledger are written."""
    root = Path(__file__).parent.parent
    json_file = root / "outputs" / "audit_verification_results.json"
    ledger_file = root / "outputs" / "FALSIFICATION_LEDGER.md"

    assert json_file.is_file(), "audit_verification_results.json missing in outputs/"
    assert ledger_file.is_file(), "FALSIFICATION_LEDGER.md missing in outputs/"

    content = ledger_file.read_text(encoding="utf-8")
    assert "Zero-Trust Scientific Falsification & Anti-Hallucination Ledger" in content
    assert "Table A: Hard-Validated Candidates" in content
    assert "Table B: Refuted / Inconclusive Candidates" in content
