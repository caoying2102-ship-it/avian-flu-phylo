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
    Script expected location:

    project_root/
    └── 03_trimal/
        └── scripts/
            └── prepare_trimal_inputs.py
    """
    return Path(__file__).resolve().parents[2]


def build_file_pairs(project_root, segments=SEGMENTS):
    """Create the mapping between MSA outputs and TrimAl inputs."""
    msa_output_dir = (
        project_root
        / "02_msa"
        / "output"
    )

    trimal_input_dir = (
        project_root
        / "03_trimal"
        / "input"
    )

    file_pairs = []

    for segment in segments:
        source_dir = (
            msa_output_dir / segment
        )

        destination_dir = (
            trimal_input_dir / segment
        )

        source_fasta = (
            source_dir
            / f"{segment}_reference_dedup_aligned.fasta"
        )

        source_metadata = (
            source_dir
            / f"{segment}_metadata_reference_dedup_aligned.xlsx"
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
    Check all MSA results before copying.

    If any file is missing or empty, the script stops to avoid generating incomplete TrimAl inputs.
    """
    missing_files = []
    empty_files = []

    for _, _, source_file, _ in file_pairs:
        if not source_file.exists():
            missing_files.append(source_file)
        elif not source_file.is_file():
            missing_files.append(source_file)
        elif source_file.stat().st_size == 0:
            empty_files.append(source_file)

    messages = []

    if missing_files:
        messages.append(
            "Missing files:\n"
            + "\n".join(
                f"  - {path}"
                for path in missing_files
            )
        )

    if empty_files:
        messages.append(
            "Empty files:\n"
            + "\n".join(
                f"  - {path}"
                for path in empty_files
            )
        )

    if messages:
        raise FileNotFoundError(
            "\n\n".join(messages)
            + "\n\nRun msa.py and matchmetadata.py first."
        )


def copy_and_verify(
    source_file,
    destination_file,
):
    """Copy a file and verify the destination size."""
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

    source_size = source_file.stat().st_size
    destination_size = (
        destination_file.stat().st_size
    )

    if source_size != destination_size:
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

    msa_output_dir = (
        project_root
        / "02_msa"
        / "output"
    )

    trimal_dir = (
        project_root
        / "03_trimal"
    )

    if not msa_output_dir.exists():
        raise FileNotFoundError(
            f"MSA output directory does not exist: "
            f"{msa_output_dir}"
        )

    if not trimal_dir.exists():
        raise FileNotFoundError(
            f"TrimAl directory does not exist: "
            f"{trimal_dir}"
        )

    file_pairs = build_file_pairs(project_root, (segment,))

    # Check all files before copying
    validate_source_files(file_pairs)

    copied_count = 0

    for (
        segment,
        file_type,
        source_file,
        destination_file,
    ) in file_pairs:
        copy_and_verify(
            source_file=source_file,
            destination_file=destination_file,
        )

        copied_count += 1

        print(
            f"{segment} {file_type}: "
            f"{destination_file}"
        )

    print()
    print(
        f"Completed: {copied_count} files copied "
        "for 1 segment."
    )
    print(
        f"TrimAl input directory: "
        f"{trimal_dir / 'input'}"
    )


if __name__ == "__main__":
    main()
