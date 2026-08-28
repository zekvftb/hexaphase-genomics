"""Unit and integration tests for Module 0 (Ingestion & Validation)."""

import json
from pathlib import Path
import pytest

from bio_arch.contracts import ValidationStatus
from bio_arch.modules.ingestion import (
    detect_format,
    ingest_file,
    stream_bed,
    stream_edge_list,
    stream_fasta,
    stream_fastq,
    ValidationReport,
)

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


def test_detect_format():
    """Verify format detection works from extension and file contents."""
    assert detect_format(FIXTURES_DIR / "valid_dna.fasta") == "fasta"
    assert detect_format(FIXTURES_DIR / "valid_reads.fastq") == "fastq"
    assert detect_format(FIXTURES_DIR / "valid_features.bed") == "bed"
    assert detect_format(FIXTURES_DIR / "valid_network.tsv") == "tsv"


def test_ingest_valid_fasta(tmp_path: Path):
    """Test ingestion of valid FASTA sequence file."""
    manifest, report, result = ingest_file(
        file_path=FIXTURES_DIR / "valid_dna.fasta",
        dataset_id="test_dna_01",
        organism="Synthetic Test",
        license="CC0",
        source="Unit Test Fixture",
        output_dir=tmp_path,
    )

    assert report.is_valid is True
    assert report.total_records == 2
    assert report.valid_records == 2
    assert manifest.validation_status == ValidationStatus.VALID
    assert result.run_metadata.status == "success"

    # Verify emitted files
    manifest_file = tmp_path / "test_dna_01_manifest.json"
    assert manifest_file.is_file()
    assert json.loads(manifest_file.read_text())["checksum"] == manifest.checksum


def test_ingest_invalid_alphabet_fasta(tmp_path: Path):
    """Verify invalid alphabets are flagged and not silently accepted."""
    manifest, report, result = ingest_file(
        file_path=FIXTURES_DIR / "invalid_alphabet.fasta",
        dataset_id="test_corrupt_dna",
        organism="Synthetic Corrupt",
        license="CC0",
        source="Unit Test Fixture",
        output_dir=tmp_path,
    )

    assert report.is_valid is False
    assert manifest.validation_status == ValidationStatus.INVALID
    assert result.run_metadata.status == "failed"
    assert len(report.issues) > 0
    assert any("Illegal characters" in issue.message for issue in report.issues)


def test_ingest_valid_fastq(tmp_path: Path):
    """Test streaming and validation of FASTQ file."""
    manifest, report, result = ingest_file(
        file_path=FIXTURES_DIR / "valid_reads.fastq",
        dataset_id="test_reads_01",
        organism="Sequencing Sample",
        license="MIT",
        source="Unit Test Fixture",
        output_dir=tmp_path,
    )

    assert report.is_valid is True
    assert report.total_records == 2
    assert report.valid_records == 2
    assert manifest.validation_status == ValidationStatus.VALID


def test_ingest_mismatched_length_fastq():
    """Verify FASTQ with mismatched sequence/quality length fails validation."""
    report = ValidationReport(format_detected="fastq")
    records = list(stream_fastq(FIXTURES_DIR / "mismatched_length.fastq", report))

    assert len(records) == 0
    assert report.is_valid is False
    assert any("does not match quality score length" in issue.message for issue in report.issues)


def test_ingest_bed_coordinates():
    """Test BED file coordinate checking (negative coordinates, inverted ranges)."""
    # Valid BED
    report_valid = ValidationReport(format_detected="bed")
    valid_records = list(stream_bed(FIXTURES_DIR / "valid_features.bed", report_valid))
    assert report_valid.is_valid is True
    assert len(valid_records) == 3
    assert valid_records[0]["chrom"] == "chr1"
    assert valid_records[0]["start"] == 100
    assert valid_records[0]["end"] == 200

    # Invalid coordinates
    report_invalid = ValidationReport(format_detected="bed")
    invalid_records = list(stream_bed(FIXTURES_DIR / "invalid_coordinates.bed", report_invalid))
    assert report_invalid.is_valid is False
    assert any("less than start" in issue.message for issue in report_invalid.issues)
    assert any("cannot be negative" in issue.message for issue in report_invalid.issues)


def test_ingest_edge_list():
    """Test edge list network file parsing."""
    report = ValidationReport(format_detected="tsv")
    edges = list(stream_edge_list(FIXTURES_DIR / "valid_network.tsv", report))

    assert report.is_valid is True
    assert len(edges) == 3
    assert edges[0]["source"] == "GeneA"
    assert edges[0]["target"] == "GeneB"
    assert edges[0]["confidence"] == 0.95


def test_normalization_tracking(tmp_path: Path):
    """Ensure casing normalizations are explicitly recorded, never silent."""
    test_fasta = tmp_path / "mixed_case.fasta"
    test_fasta.write_text(">mixed\natgcATGC\n")

    report = ValidationReport(format_detected="fasta")
    records = list(stream_fasta(test_fasta, report, normalize_case=True))

    assert len(records) == 1
    assert records[0]["sequence"] == "ATGCATGC"
    # Must record normalization
    assert "uppercase_normalization" in report.normalizations
