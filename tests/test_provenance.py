"""Unit tests for provenance, checksums, environment logging, and seed management."""

from datetime import datetime
from pathlib import Path
import pytest

from bio_arch.provenance import (
    SeedManager,
    compute_sha256,
    get_system_environment,
    now_iso,
)


def test_now_iso():
    """Verify now_iso returns valid ISO-8601 UTC timestamp."""
    iso_str = now_iso()
    assert isinstance(iso_str, str)
    dt = datetime.fromisoformat(iso_str)
    assert dt.tzinfo is not None


def test_compute_sha256_bytes():
    """Verify SHA-256 computation on raw bytes."""
    # SHA-256 of empty string is e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
    empty_hash = compute_sha256(b"")
    assert empty_hash == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"

    hello_hash = compute_sha256(b"hello world")
    assert hello_hash == "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9"


def test_compute_sha256_file(tmp_path: Path):
    """Verify SHA-256 computation on files using chunked streaming."""
    test_file = tmp_path / "sample.txt"
    test_file.write_bytes(b"hello world")

    computed = compute_sha256(test_file)
    assert computed == "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9"


def test_compute_sha256_file_not_found():
    """Verify FileNotFoundError raised when target path does not exist."""
    with pytest.raises(FileNotFoundError):
        compute_sha256("non_existent_file_xyz_123.bin")


def test_get_system_environment():
    """Verify runtime environment metadata dictionary has required fields."""
    env = get_system_environment()
    assert "python_version" in env
    assert "platform" in env
    assert "installed_packages" in env
    assert isinstance(env["installed_packages"], dict)


def test_seed_manager_determinism():
    """Verify SeedManager derives reproducible and distinct sub-seeds."""
    sm1 = SeedManager(master_seed=12345)
    sm2 = SeedManager(master_seed=12345)

    # Identical contexts must yield identical seeds
    seed_a1 = sm1.derive_seed("module_1_shuffles")
    seed_a2 = sm2.derive_seed("module_1_shuffles")
    assert seed_a1 == seed_a2

    # Different contexts must yield distinct seeds
    seed_b1 = sm1.derive_seed("module_2_graphs")
    assert seed_a1 != seed_b1

    # Seeds must be non-negative integers
    assert 0 <= seed_a1 < (2**31 - 1)
