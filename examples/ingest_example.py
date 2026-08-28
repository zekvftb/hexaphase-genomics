"""Example demonstrating Module 0 (Data Ingestion & Validation) on synthetic data."""

from pathlib import Path
import sys

# Ensure src is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from bio_arch.modules.ingestion import ingest_file


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    sample_file = repo_root / "data" / "synthetic_sample.fasta"
    output_dir = repo_root / "outputs" / "sample_ingest"

    print("Running Module 0 ingestion on:", sample_file)

    manifest, report, result = ingest_file(
        file_path=sample_file,
        dataset_id="ds_sample_lambda",
        organism="Synthetic Lambda/Lac model",
        license="CC0-1.0",
        source="In-silico benchmark",
        sequence_type="DNA",
        output_dir=output_dir,
    )

    print("\n--- Ingestion Summary ---")
    print(f"Format:        {report.format_detected}")
    print(f"Total records: {report.total_records}")
    print(f"Valid records: {report.valid_records}")
    print(f"Validation:    {manifest.validation_status.value}")
    print(f"SHA-256:       {manifest.checksum}")
    print(f"Artifacts saved to: {output_dir}")
    for p in result.artifact_paths:
        print(f"  -> {p}")


if __name__ == "__main__":
    main()
