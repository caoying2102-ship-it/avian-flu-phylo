#!/usr/bin/env python3

import argparse
import shutil
from pathlib import Path


SEGMENTS = (
    "MP",
    "HA",
    "NA",
    "PB2",
    "PB1",
    "PA",
    "NP",
    "NS",
)


def get_project_root():
    """
    This script is expected to live at:
    project_root/02_msa/scripts/prepare_msa_inputs.py
    """
    return Path(__file__).resolve().parents[2]


def get_file_pairs(project_root, segments=SEGMENTS):
    """Create the source-to-destination mapping for all files."""
    genes_dir = (
        project_root
        / "01_sample_preparation"
        / "output"
        / "genes"
    )

    msa_input_dir = (
        project_root
        / "02_msa"
        / "input"
    )

    file_pairs = []

    for segment in segments:
        source_dir = genes_dir / segment
        destination_dir = msa_input_dir / segment

        source_fasta = (
            source_dir
            / f"{segment}_reference_dedup.fasta"
        )

        source_metadata = (
            source_dir
            / f"{segment}_metadata_reference_dedup.xlsx"
        )

        destination_fasta = (
            destination_dir
            / source_fasta.name
        )

        destination_metadata = (
            destination_dir
            / source_metadata.name
        )

        file_pairs.extend(
            [
                (
                    segment,
                    "FASTA",
                    source_fasta,
                    destination_fasta,
                ),
                (
                    segment,
                    "metadata",
                    source_metadata,
                    destination_metadata,
                ),
            ]
        )

    return file_pairs


def validate_source_files(file_pairs):
    """
    Check all source files before copying starts.

    If any required file is missing, the script stops so incomplete MSA inputs are not produced.
    """
    missing_files = []

    for _, _, source_file, _ in file_pairs:
        if not source_file.exists():
            missing_files.append(source_file)
        elif not source_file.is_file():
            missing_files.append(source_file)

    if missing_files:
        missing_text = "\n".join(
            f"  - {path}"
            for path in missing_files
        )

        raise FileNotFoundError(
            "The following required input files are missing:\n"
            f"{missing_text}\n\n"
            "Run the extraction and SeqKit deduplication "
            "steps first."
        )


def copy_file(source_file, destination_file):
    """
    Copy files while preserving metadata such as modification time.

    If the destination file already exists, it is updated to match the latest source file.
    """
    destination_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    shutil.copy2(
        source_file,
        destination_file,
    )

    if not destination_file.exists():
        raise RuntimeError(
            f"Copy failed: {destination_file}"
        )

    if source_file.stat().st_size != (
        destination_file.stat().st_size
    ):
        raise RuntimeError(
            "File-size verification failed:\n"
            f"Source: {source_file}\n"
            f"Destination: {destination_file}"
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--segment", required=True, type=str.upper, choices=SEGMENTS)
    segment = parser.parse_args().segment
    project_root = get_project_root()

    sample_output_dir = (
        project_root
        / "01_sample_preparation"
        / "output"
        / "genes"
    )

    msa_dir = project_root / "02_msa"

    if not sample_output_dir.exists():
        raise FileNotFoundError(
            "Gene output directory does not exist: "
            f"{sample_output_dir}"
        )

    if not msa_dir.exists():
        raise FileNotFoundError(
            f"02_msa directory does not exist: {msa_dir}"
        )

    file_pairs = get_file_pairs(project_root, (segment,))

    # Check all files before copying
    validate_source_files(file_pairs)

    copied_files = 0

    for (
        segment,
        file_type,
        source_file,
        destination_file,
    ) in file_pairs:
        copy_file(
            source_file=source_file,
            destination_file=destination_file,
        )

        copied_files += 1

        print(
            f"{segment} {file_type}: "
            f"{destination_file}"
        )

    print()
    print(
        f"Completed: {copied_files} files copied "
        "for 1 segment."
    )
    print(
        f"MSA input directory: "
        f"{project_root / '02_msa' / 'input'}"
    )


if __name__ == "__main__":
    main()
