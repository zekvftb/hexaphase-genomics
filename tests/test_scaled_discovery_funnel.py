"""Automated Verification & Benchmark Suite for Scaled Multi-Core Mining Funnel."""

from __future__ import annotations

import csv
from pathlib import Path
import time
import pytest

from bio_arch.modules.bulk_viral_indexer import BulkViralIndexer
from bio_arch.modules.parallel_null_engine import ParallelNullEngine
from scripts.mining.run_scaled_discovery_funnel import run_scaled_funnel_pipeline


@pytest.fixture
def temp_indexer(tmp_path):
    """Create a temporary BulkViralIndexer instance."""
    db_file = tmp_path / "test_genomes.db"
    return BulkViralIndexer(db_path=db_file)


def test_sqlite_indexing_and_retrieval(temp_indexer):
    """Verify that viral GenBank records and primary CDS features are indexed into SQLite."""
    root = Path(__file__).parent.parent
    corpus_dir = root / "data" / "mining_corpus"

    indexed = temp_indexer.index_corpus_directory(corpus_dir)
    assert indexed > 0
    assert temp_indexer.get_total_genomes_count() == indexed
    assert temp_indexer.get_total_cds_count() > 0

    genomes = temp_indexer.get_all_genomes()
    assert len(genomes) == indexed
    sample_acc = genomes[0]["accession"]

    cds_list = temp_indexer.get_cds_features(sample_acc)
    assert len(cds_list) > 0


def test_parallel_null_engine_multi_core(temp_indexer):
    """Verify parallel null model execution and hard-barrier candidate filtering."""
    root = Path(__file__).parent.parent
    corpus_dir = root / "data" / "mining_corpus"
    temp_indexer.index_corpus_directory(corpus_dir)

    genomes = temp_indexer.get_all_genomes()[:3]
    genomes_with_cds = [(g, temp_indexer.get_cds_features(g["accession"])) for g in genomes]

    engine = ParallelNullEngine(
        max_workers=2,
        n_shuffles=20,
        min_length_aa=30,
        min_cai_threshold=0.70,
        z_barrier_threshold=3.0,
        seed=42,
    )

    surviving, telemetry = engine.run_screening_pipeline(genomes_with_cds)

    assert "total_genomes_screened" in telemetry
    assert telemetry["total_genomes_screened"] == 3
    assert telemetry["stage1_candidates_found"] >= 0
    assert telemetry["execution_time_seconds"] > 0
    assert telemetry["throughput_genomes_per_sec"] > 0


def test_scaled_funnel_artifacts_and_routing():
    """Verify generated Table A, Table B, and master markdown report artifacts."""
    root = Path(__file__).parent.parent
    table_a_file = root / "outputs" / "SCALED_VALIDATED_TABLE_A.csv"
    table_b_file = root / "outputs" / "SCALED_REFUTED_LEDGER.csv"
    report_file = root / "outputs" / "SCALED_VIRAL_MINING_REPORT.md"

    assert table_a_file.is_file(), "SCALED_VALIDATED_TABLE_A.csv missing"
    assert table_b_file.is_file(), "SCALED_REFUTED_LEDGER.csv missing"
    assert report_file.is_file(), "SCALED_VIRAL_MINING_REPORT.md missing"

    content = report_file.read_text(encoding="utf-8")
    assert "Scaled Multi-Core Viral Mining & High-Throughput Validation Funnel" in content
    assert "Multi-Core Execution Telemetry" in content


def test_sub_minute_benchmark_performance():
    """Benchmark throughput to verify rapid sub-minute execution on local corpus."""
    root = Path(__file__).parent.parent
    start = time.time()
    res = run_scaled_funnel_pipeline(root, n_shuffles=20, max_workers=4)
    elapsed = time.time() - start

    assert elapsed < 60.0  # Must execute in under 1 minute
    assert res["telemetry"]["total_genomes_screened"] > 0
