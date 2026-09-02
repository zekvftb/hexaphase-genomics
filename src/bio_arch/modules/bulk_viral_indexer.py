"""Bulk Viral Genome Registry & SQLite Caching Module.

Ingests, parses, and indexes complete viral RefSeq genomes and their primary CDS features
into an idempotent SQLite database (data/genome_cache/genomes.db) for high-throughput mining.
"""

from __future__ import annotations

from pathlib import Path
import sqlite3
import sys
from typing import Any, Optional

from Bio import SeqIO

_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))


class BulkViralIndexer:
    """Manages idempotent SQLite indexing and retrieval for bulk viral genomes."""

    DEFAULT_DB_PATH = _ROOT / "data" / "genome_cache" / "genomes.db"

    def __init__(self, db_path: Path | None = None) -> None:
        self.db_path = db_path if db_path is not None else self.DEFAULT_DB_PATH
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.init_db()

    def _get_connection(self) -> sqlite3.Connection:
        conn = sqlite3.connect(str(self.db_path))
        conn.row_factory = sqlite3.Row
        return conn

    def init_db(self) -> None:
        """Create database tables if they do not exist."""
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS genomes (
                    accession TEXT PRIMARY KEY,
                    taxid TEXT,
                    organism TEXT,
                    family TEXT,
                    genome_length_bp INTEGER,
                    sequence TEXT,
                    gbk_raw TEXT
                )
            """)
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS cds_features (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    accession TEXT,
                    gene TEXT,
                    product TEXT,
                    start_bp INTEGER,
                    end_bp INTEGER,
                    strand INTEGER,
                    sequence TEXT,
                    translation TEXT,
                    FOREIGN KEY (accession) REFERENCES genomes(accession)
                )
            """)
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_cds_acc ON cds_features(accession)")
            conn.commit()

    def index_genbank_file(self, gbk_path: Path) -> bool:
        """Parse and index a single GenBank file idempotently."""
        if not gbk_path.is_file():
            return False

        try:
            record = SeqIO.read(str(gbk_path), "genbank")
            raw_text = gbk_path.read_text(encoding="utf-8", errors="replace")
        except Exception:
            return False

        acc = record.id
        organism = record.annotations.get("organism", record.description)
        taxonomy = record.annotations.get("taxonomy", [])
        family = taxonomy[-1] if taxonomy else "Unclassified Virus"
        seq_str = str(record.seq).upper()
        length_bp = len(seq_str)
        taxid = ""

        # Extract TaxID if present
        for feat in record.features:
            if feat.type == "source":
                db_xrefs = feat.qualifiers.get("db_xref", [])
                for xref in db_xrefs:
                    if xref.startswith("taxon:"):
                        taxid = xref.split(":")[-1]
                        break

        with self._get_connection() as conn:
            cursor = conn.cursor()
            # Check if already indexed
            cursor.execute("SELECT accession FROM genomes WHERE accession = ?", (acc,))
            if cursor.fetchone() is not None:
                return True  # Already indexed

            # Insert genome
            cursor.execute(
                """
                INSERT OR REPLACE INTO genomes (accession, taxid, organism, family, genome_length_bp, sequence, gbk_raw)
                VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                (acc, taxid, organism, family, length_bp, seq_str, raw_text),
            )

            # Insert CDS features
            for feat in record.features:
                if feat.type == "CDS":
                    gene = feat.qualifiers.get("gene", [""])[0] or "CDS"
                    prod = feat.qualifiers.get("product", [""])[0] or "viral protein"
                    trans = feat.qualifiers.get("translation", [""])[0]
                    strand = feat.location.strand or 1
                    cds_seq = str(feat.extract(record.seq)).upper()
                    start_bp = int(feat.location.start)
                    end_bp = int(feat.location.end)

                    cursor.execute(
                        """
                        INSERT INTO cds_features (accession, gene, product, start_bp, end_bp, strand, sequence, translation)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        (acc, gene, prod, start_bp, end_bp, strand, cds_seq, trans),
                    )

            conn.commit()
        return True

    def index_corpus_directory(self, corpus_dir: Path) -> int:
        """Index all GenBank files in a directory."""
        if not corpus_dir.is_dir():
            return 0

        indexed_count = 0
        for gbk_path in corpus_dir.glob("*.gbk"):
            if self.index_genbank_file(gbk_path):
                indexed_count += 1
        return indexed_count

    def get_all_genomes(self) -> list[dict[str, Any]]:
        """Retrieve all indexed viral genomes."""
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT accession, taxid, organism, family, genome_length_bp, sequence FROM genomes")
            return [dict(row) for row in cursor.fetchall()]

    def get_cds_features(self, accession: str) -> list[dict[str, Any]]:
        """Retrieve all primary CDS features for an accession."""
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT * FROM cds_features WHERE accession = ?", (accession,))
            return [dict(row) for row in cursor.fetchall()]

    def get_total_genomes_count(self) -> int:
        """Count total genomes stored in SQLite."""
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT COUNT(*) FROM genomes")
            return cursor.fetchone()[0]

    def get_total_cds_count(self) -> int:
        """Count total CDS features stored in SQLite."""
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT COUNT(*) FROM cds_features")
            return cursor.fetchone()[0]
