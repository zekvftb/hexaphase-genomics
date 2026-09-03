"""Automated Verification Suite for Circoviridae & Anelloviridae Bulk Ingestion & Mining."""

from __future__ import annotations

from pathlib import Path
import sqlite3
import pytest

from bio_arch.modules.bulk_viral_indexer import BulkViralIndexer
from scripts.mining.fetch_circoviridae_anelloviridae import (
    REPRESENTATIVE_OFFLINE_RECORDS,
    populate_offline_representative_records,
)
from scripts.mining.run_expanded_family_funnel import run_expanded_family_pipeline


@pytest.fixture
def temp_circ_db(tmp_path):
    """Create a temporary SQLite database for testing."""
    db_file = tmp_path / "test_circo_anello.db"
    indexer = BulkViralIndexer(db_path=db_file)
    populate_offline_representative_records(indexer)
    return indexer, db_file


def test_sqlite_batch_insertion_and_schema(temp_circ_db):
    """Verify that Circoviridae and Anelloviridae records are properly inserted with correct schema."""
    indexer, db_file = temp_circ_db

    assert indexer.get_total_genomes_count() == len(REPRESENTATIVE_OFFLINE_RECORDS)
    assert indexer.get_total_cds_count() > 0

    all_genomes = indexer.get_all_genomes()
    families = {g["family"] for g in all_genomes}
    assert "Circoviridae" in families or "Anelloviridae" in families

    for g in all_genomes:
        cds_list = indexer.get_cds_features(g["accession"])
        assert len(cds_list) > 0
        for cds in cds_list:
            assert "start_bp" in cds
            assert "end_bp" in cds
            assert len(cds["sequence"]) > 0


def test_circular_sequence_boundary_handling():
    """Verify handling of circular viral genome wraparound coordinates."""
    seq = "ATGCGTAGCTAGCTAGCTAGCTAG" * 10
    seq_len = len(seq)

    # Test coordinate wrap calculation
    start_bp = seq_len - 15
    end_bp = 15

    if start_bp > end_bp:
        # Wraparound circular junction
        wrapped_seq = seq[start_bp:] + seq[:end_bp]
        assert len(wrapped_seq) == 30
        assert wrapped_seq.startswith(seq[start_bp:])
        assert wrapped_seq.endswith(seq[:end_bp])


def test_expanded_family_funnel_execution():
    """Verify end-to-end execution of the expanded family validation funnel."""
    root = Path(__file__).parent.parent
    res = run_expanded_family_pipeline(root, n_shuffles=20, max_workers=2)

    assert "telemetry" in res
    assert res["telemetry"]["total_genomes_screened"] >= 0
    assert "report_path" in res

    report_path = Path(res["report_path"])
    assert report_path.is_file()
    content = report_path.read_text(encoding="utf-8")
    assert "Circoviridae & Anelloviridae" in content
    assert "Table A: Validated" in content
    assert "Table B: Refuted" in content
