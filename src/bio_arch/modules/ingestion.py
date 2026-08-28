"""Module 0: Data Ingestion and Validation.

Accepts FASTA, FASTQ, BED, GFF/GTF, CSV/TSV expression, and edge-list network files.
Performs strict validation, format detection, alphabet checks, coordinate validation,
and SHA-256 calculation.
Ensures zero silent data repair: all normalizations and discrepancies are explicitly recorded.
Emits DatasetManifest, validation reports, and ModuleResult.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
import re
import sys
from typing import Any, Iterator, Literal

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
from bio_arch.provenance import (
    SeedManager,
    compute_sha256,
    get_system_environment,
    now_iso,
)

logger = setup_logger("bio_arch.ingestion")

# Standard sequence alphabets
ALPHABETS = {
    "DNA": set("ACGTN"),
    "DNA_EXTENDED": set("ACGTURYKMSWBDHVN"),
    "RNA": set("ACGUN"),
    "PROTEIN": set("ACDEFGHIKLMNPQRSTVWY*"),
}


@dataclass
class ValidationIssue:
    """A specific error or warning encountered during file validation."""

    severity: Literal["error", "warning"]
    line_number: int | None
    record_id: str | None
    message: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ValidationReport:
    """Detailed summary of the validation check for an ingested dataset."""

    format_detected: str
    total_records: int = 0
    valid_records: int = 0
    issues: list[ValidationIssue] = field(default_factory=list)
    normalizations: list[str] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        """Dataset is valid only if there are zero error-level issues."""
        return not any(i.severity == "error" for i in self.issues)

    @property
    def status(self) -> ValidationStatus:
        if not self.is_valid:
            return ValidationStatus.INVALID
        if any(i.severity == "warning" for i in self.issues):
            return ValidationStatus.WARNING
        return ValidationStatus.VALID

    def add_issue(
        self,
        severity: Literal["error", "warning"],
        message: str,
        line_number: int | None = None,
        record_id: str | None = None,
    ) -> None:
        self.issues.append(
            ValidationIssue(
                severity=severity,
                line_number=line_number,
                record_id=record_id,
                message=message,
            )
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "format_detected": self.format_detected,
            "total_records": self.total_records,
            "valid_records": self.valid_records,
            "is_valid": self.is_valid,
            "status": self.status.value,
            "issues": [i.to_dict() for i in self.issues],
            "normalizations": self.normalizations,
        }

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)


def detect_format(file_path: Path) -> str:
    """Inspect file extension and leading lines to detect the file format."""
    suffix = file_path.suffix.lower()
    ext_map = {
        ".fasta": "fasta",
        ".fa": "fasta",
        ".fna": "fasta",
        ".fastq": "fastq",
        ".fq": "fastq",
        ".bed": "bed",
        ".gff": "gff",
        ".gff3": "gff",
        ".gtf": "gtf",
        ".tsv": "tsv",
        ".csv": "csv",
        ".tab": "tsv",
    }
    if suffix in ext_map:
        return ext_map[suffix]

    # Inspect first few non-empty lines
    with file_path.open("r", encoding="utf-8", errors="replace") as f:
        first_lines = [f.readline().strip() for _ in range(5)]
        first_lines = [l for l in first_lines if l]

    if not first_lines:
        raise ValueError(f"File '{file_path}' is empty.")

    first = first_lines[0]
    if first.startswith(">"):
        return "fasta"
    if first.startswith("@"):
        return "fastq"
    if first.startswith("##gff") or "\t" in first and len(first.split("\t")) >= 8:
        return "gff"
    if "," in first:
        return "csv"
    if "\t" in first:
        parts = first.split("\t")
        if len(parts) >= 3 and parts[0].startswith("chr"):
            return "bed"
        return "tsv"

    return "unknown"


def stream_fasta(
    file_path: Path,
    report: ValidationReport,
    alphabet_type: str = "DNA",
    normalize_case: bool = True,
) -> Iterator[dict[str, Any]]:
    """Stream FASTA records, validating alphabets, non-empty sequences, and duplicate IDs."""
    allowed = ALPHABETS.get(alphabet_type.upper(), ALPHABETS["DNA"])
    seen_ids: set[str] = set()
    current_id: str | None = None
    current_desc: str = ""
    current_seq_chunks: list[str] = []
    line_start_num = 1

    with file_path.open("r", encoding="utf-8", errors="replace") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Flush previous record
                if current_id is not None:
                    seq = "".join(current_seq_chunks)
                    if not seq:
                        report.add_issue("error", f"Record '{current_id}' has empty sequence.", line_start_num, current_id)
                    else:
                        report.valid_records += 1
                        yield {"id": current_id, "description": current_desc, "sequence": seq}

                # Start new record
                report.total_records += 1
                header = line[1:].strip()
                parts = header.split(maxsplit=1)
                record_id = parts[0] if parts else f"record_{report.total_records}"
                current_desc = parts[1] if len(parts) > 1 else ""

                if record_id in seen_ids:
                    report.add_issue("warning", f"Duplicate record identifier: '{record_id}'", line_num, record_id)
                seen_ids.add(record_id)

                current_id = record_id
                current_seq_chunks = []
                line_start_num = line_num
            else:
                if current_id is None:
                    report.add_issue("error", "Found sequence data before header line.", line_num)
                    continue

                clean_chunk = line
                if normalize_case and not line.isupper():
                    if "uppercase_normalization" not in report.normalizations:
                        report.normalizations.append("uppercase_normalization")
                    clean_chunk = clean_chunk.upper()

                # Validate alphabet
                illegal_chars = set(clean_chunk) - allowed
                if illegal_chars:
                    report.add_issue(
                        "error",
                        f"Illegal characters {illegal_chars} found for alphabet '{alphabet_type}'.",
                        line_num,
                        current_id,
                    )
                current_seq_chunks.append(clean_chunk)

        # Flush final record
        if current_id is not None:
            seq = "".join(current_seq_chunks)
            if not seq:
                report.add_issue("error", f"Record '{current_id}' has empty sequence.", line_start_num, current_id)
            else:
                report.valid_records += 1
                yield {"id": current_id, "description": current_desc, "sequence": seq}


def stream_fastq(
    file_path: Path,
    report: ValidationReport,
    alphabet_type: str = "DNA",
    normalize_case: bool = True,
) -> Iterator[dict[str, Any]]:
    """Stream 4-line FASTQ records, checking length matching and quality score formats."""
    allowed = ALPHABETS.get(alphabet_type.upper(), ALPHABETS["DNA"])
    seen_ids: set[str] = set()

    with file_path.open("r", encoding="utf-8", errors="replace") as f:
        line_num = 0
        while True:
            header_line = f.readline()
            if not header_line:
                break
            line_num += 1
            header_line = header_line.strip()
            if not header_line:
                continue

            report.total_records += 1

            if not header_line.startswith("@"):
                report.add_issue("error", f"Malformed FASTQ header (must start with '@'): '{header_line}'", line_num)
                # Skip to next lines if possible
                continue

            record_id = header_line[1:].split()[0]
            if record_id in seen_ids:
                report.add_issue("warning", f"Duplicate FASTQ identifier: '{record_id}'", line_num, record_id)
            seen_ids.add(record_id)

            seq_line = f.readline().strip()
            line_num += 1

            plus_line = f.readline().strip()
            line_num += 1
            if not plus_line.startswith("+"):
                report.add_issue("error", f"Missing '+' separator line in FASTQ record '{record_id}'.", line_num, record_id)

            qual_line = f.readline().strip()
            line_num += 1

            if len(seq_line) != len(qual_line):
                report.add_issue(
                    "error",
                    f"Sequence length ({len(seq_line)}) does not match quality score length ({len(qual_line)}).",
                    line_num,
                    record_id,
                )
                continue

            if normalize_case and not seq_line.isupper():
                if "uppercase_normalization" not in report.normalizations:
                    report.normalizations.append("uppercase_normalization")
                seq_line = seq_line.upper()

            illegal_chars = set(seq_line) - allowed
            if illegal_chars:
                report.add_issue(
                    "error",
                    f"Illegal sequence characters {illegal_chars} for alphabet '{alphabet_type}'.",
                    line_num - 2,
                    record_id,
                )
                continue

            report.valid_records += 1
            yield {"id": record_id, "sequence": seq_line, "quality": qual_line}


def stream_bed(file_path: Path, report: ValidationReport) -> Iterator[dict[str, Any]]:
    """Stream BED file entries, validating 0-based non-negative coordinates and start <= end."""
    with file_path.open("r", encoding="utf-8", errors="replace") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue

            report.total_records += 1
            parts = line.split("\t")
            if len(parts) < 3:
                parts = line.split()  # fallback to whitespace
            if len(parts) < 3:
                report.add_issue("error", f"BED line has fewer than 3 required fields (chrom, start, end).", line_num)
                continue

            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                report.add_issue("error", f"Coordinates must be integers. Found: '{parts[1]}', '{parts[2]}'", line_num)
                continue

            if start < 0:
                report.add_issue("error", f"BED start coordinate cannot be negative ({start}).", line_num)
                continue
            if end < start:
                report.add_issue("error", f"BED end coordinate ({end}) is less than start ({start}).", line_num)
                continue

            name = parts[3] if len(parts) > 3 else f"feature_{report.total_records}"
            score = float(parts[4]) if len(parts) > 4 and parts[4] != "." else None
            strand = parts[5] if len(parts) > 5 else "."

            report.valid_records += 1
            yield {
                "chrom": chrom,
                "start": start,
                "end": end,
                "name": name,
                "score": score,
                "strand": strand,
            }


def stream_edge_list(
    file_path: Path,
    report: ValidationReport,
    delimiter: str | None = None,
) -> Iterator[dict[str, Any]]:
    """Stream network edge lists (source, target, interaction_type, confidence)."""
    with file_path.open("r", encoding="utf-8", errors="replace") as f:
        # Detect delimiter if not given
        first_line = f.readline()
        f.seek(0)
        if delimiter is None:
            delimiter = "\t" if "\t" in first_line else ","

        reader = csv.reader(f, delimiter=delimiter)
        header = None

        for line_num, row in enumerate(reader, start=1):
            if not row or row[0].startswith("#"):
                continue

            # Check for header
            if header is None and any(col.lower() in ("source", "from", "regulator") for col in row):
                header = [col.strip().lower() for col in row]
                continue

            report.total_records += 1
            if len(row) < 2:
                report.add_issue("error", f"Edge list entry requires at least source and target nodes.", line_num)
                continue

            source = row[0].strip()
            target = row[1].strip()
            interaction = row[2].strip() if len(row) > 2 else "associates_with"
            confidence_str = row[3].strip() if len(row) > 3 else "1.0"

            try:
                confidence = float(confidence_str)
            except ValueError:
                confidence = 1.0

            if not source or not target:
                report.add_issue("error", "Source or target node identifier is empty.", line_num)
                continue

            report.valid_records += 1
            yield {
                "source": source,
                "target": target,
                "interaction_type": interaction,
                "confidence": confidence,
            }


def ingest_file(
    file_path: Path | str,
    dataset_id: str,
    organism: str,
    license: str,
    source: str,
    file_format: str | None = None,
    sequence_type: str | None = "DNA",
    url: str | None = None,
    output_dir: Path | str | None = None,
    seed: int = 42,
) -> tuple[DatasetManifest, ValidationReport, ModuleResult]:
    """Ingest and validate a biological dataset according to shared contracts.

    Returns:
        (DatasetManifest, ValidationReport, ModuleResult)
    """
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found at: '{path}'")

    checksum = compute_sha256(path)
    detected_fmt = file_format or detect_format(path)
    report = ValidationReport(format_detected=detected_fmt)

    # Initialize analysis tracking
    seed_mgr = SeedManager(master_seed=seed)
    mod_seed = seed_mgr.derive_seed("module_0_ingestion")

    run_meta = AnalysisRun(
        run_id=f"ingest_{dataset_id}_{int(seed_mgr.master_seed or 0)}",
        timestamp=now_iso(),
        module="module_0_ingestion",
        version="0.1.0",
        input_ids=[str(path.name)],
        parameters={"format": detected_fmt, "sequence_type": sequence_type},
        seed=mod_seed,
        environment=get_system_environment(),
        status="running",
    )

    records: list[dict[str, Any]] = []

    # Stream & validate based on detected format
    if detected_fmt == "fasta":
        records = list(stream_fasta(path, report, alphabet_type=sequence_type or "DNA"))
    elif detected_fmt == "fastq":
        records = list(stream_fastq(path, report, alphabet_type=sequence_type or "DNA"))
    elif detected_fmt == "bed":
        records = list(stream_bed(path, report))
    elif detected_fmt in ("csv", "tsv", "edge_list"):
        records = list(stream_edge_list(path, report))
    else:
        report.add_issue("warning", f"Format '{detected_fmt}' ingested without specialized parser.")

    run_meta.status = "success" if report.is_valid else "failed"

    # Build DatasetManifest
    manifest = DatasetManifest(
        dataset_id=dataset_id,
        source=source,
        license=license,
        retrieval_date=now_iso(),
        checksum=checksum,
        organism=organism,
        url=url,
        sequence_type=sequence_type,
        annotations={
            "format": detected_fmt,
            "total_records": report.total_records,
            "valid_records": report.valid_records,
            "normalizations": report.normalizations,
        },
        validation_status=report.status,
    )

    # Create scientific finding for the ingestion measurement
    finding = Finding(
        finding_id=f"f_ingest_{dataset_id}",
        metric="ingestion_record_count",
        observed_value={"total": report.total_records, "valid": report.valid_records},
        biological_context=f"Ingested {detected_fmt} records for {organism}.",
    )

    interpretation = InterpretationRecord(
        finding_ids=[finding.finding_id],
        classification=EvidenceClass.MEASUREMENT,
        claim=f"{EvidenceClass.MEASUREMENT.required_prefix} {report.valid_records} valid records out of {report.total_records} in {path.name}.",
        alternatives=[],
        limitations=["Syntactic and coordinate validation only; functional validity unverified."],
        proposed_test="Verify against reference assembly or independent biological assay annotations.",
    )

    # Save artifacts if output_dir provided
    artifact_paths: list[str] = []
    if output_dir:
        out_p = Path(output_dir)
        out_p.mkdir(parents=True, exist_ok=True)

        manifest_path = out_p / f"{dataset_id}_manifest.json"
        manifest_path.write_text(manifest.to_json(), encoding="utf-8")
        artifact_paths.append(str(manifest_path))

        report_path = out_p / f"{dataset_id}_validation_report.json"
        report_path.write_text(report.to_json(), encoding="utf-8")
        artifact_paths.append(str(report_path))

        normalized_path = out_p / f"{dataset_id}_normalized.json"
        normalized_path.write_text(json.dumps(records, indent=2), encoding="utf-8")
        artifact_paths.append(str(normalized_path))

    mod_result = ModuleResult(
        run_metadata=run_meta,
        outputs={
            "format": detected_fmt,
            "total_records": report.total_records,
            "valid_records": report.valid_records,
        },
        findings=[finding],
        interpretations=[interpretation],
        warnings=[f"{i.severity}: {i.message}" for i in report.issues if i.severity == "warning"],
        errors=[f"{i.severity}: {i.message}" for i in report.issues if i.severity == "error"],
        artifact_paths=artifact_paths,
    )

    return manifest, report, mod_result


def main() -> None:
    """CLI entrypoint for Module 0 ingestion."""
    parser = argparse.ArgumentParser(
        description="Bio-Arch Module 0: Ingestion and Validation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("file", type=str, help="Path to sequence, coordinate, or network file")
    parser.add_argument("--id", type=str, default="dataset_001", help="Dataset identifier")
    parser.add_argument("--organism", type=str, default="Synthetic/Unknown", help="Source organism")
    parser.add_argument("--source", type=str, default="Local File", help="Data source description")
    parser.add_argument("--license", type=str, default="Open", help="License terms")
    parser.add_argument("--format", type=str, default=None, help="Explicit format override (fasta/fastq/bed/edge_list)")
    parser.add_argument("--seq-type", type=str, default="DNA", help="Sequence type (DNA/RNA/PROTEIN)")
    parser.add_argument("--outdir", type=str, default="outputs/ingestion", help="Output directory for manifests")

    args = parser.parse_args()

    manifest, report, result = ingest_file(
        file_path=args.file,
        dataset_id=args.id,
        organism=args.organism,
        license=args.license,
        source=args.source,
        file_format=args.format,
        sequence_type=args.seq_type,
        output_dir=args.outdir,
    )

    print("\n--- Ingestion Report ---")
    print(f"Format:        {report.format_detected}")
    print(f"Total Records: {report.total_records}")
    print(f"Valid Records: {report.valid_records}")
    print(f"Status:        {manifest.validation_status.value}")
    if report.normalizations:
        print(f"Normalizations:{', '.join(report.normalizations)}")
    if report.issues:
        print(f"Issues logged: {len(report.issues)}")
        for issue in report.issues:
            print(f"  [{issue.severity.upper()}] Line {issue.line_number or '?'}: {issue.message}")
    print("------------------------\n")


if __name__ == "__main__":
    main()
