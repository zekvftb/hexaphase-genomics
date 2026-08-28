"""Provenance, hashing, environment capture, and deterministic seed management.

Ensures strict reproducibility, auditability, and data integrity across runs.
"""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import importlib.metadata
import os
from pathlib import Path
import platform
import sys
from typing import Any, Union


def now_iso() -> str:
    """Return the current UTC timestamp formatted as ISO 8601."""
    return datetime.now(timezone.utc).isoformat()


def compute_sha256(target: Union[str, Path, bytes]) -> str:
    """Compute the SHA-256 hexadecimal digest for a file path or raw bytes.

    Uses 64 KB memory chunks to safely hash files of any size without RAM exhaustion.
    """
    hasher = hashlib.sha256()

    if isinstance(target, bytes):
        hasher.update(target)
        return hasher.hexdigest()

    file_path = Path(target)
    if not file_path.is_file():
        raise FileNotFoundError(f"Cannot compute checksum: file not found at '{file_path}'")

    chunk_size = 65536  # 64 KB chunks
    with file_path.open("rb") as f:
        while chunk := f.read(chunk_size):
            hasher.update(chunk)

    return hasher.hexdigest()


def get_system_environment() -> dict[str, Any]:
    """Capture runtime environment metadata for full scientific auditability."""
    # List of relevant packages to record versions for if installed
    packages_to_check = [
        "bio-arch",
        "numpy",
        "scipy",
        "networkx",
        "pandas",
        "biopython",
        "matplotlib",
        "pytest",
        "pyyaml",
    ]

    installed_versions: dict[str, str] = {}
    for pkg in packages_to_check:
        try:
            installed_versions[pkg] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            pass

    return {
        "python_version": sys.version,
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        "installed_packages": installed_versions,
    }


class SeedManager:
    """Manages deterministic random seeds for reproducible analysis and simulations."""

    def __init__(self, master_seed: int | None = 42) -> None:
        self.master_seed = master_seed

    def derive_seed(self, context_name: str) -> int:
        """Derive a deterministic sub-seed for a specific component or module.

        Combines the master seed and context name using SHA-256 to ensure uncorrelated
        yet fully deterministic sub-seeds across different modules.
        """
        if self.master_seed is None:
            import random
            return random.randint(0, 2**31 - 1)

        seed_string = f"{self.master_seed}:{context_name}".encode("utf-8")
        hash_digest = hashlib.sha256(seed_string).hexdigest()
        # Take first 8 hex characters (32 bits) modulo 2^31-1 for safe signed 32-bit int
        return int(hash_digest[:8], 16) % (2**31 - 1)
