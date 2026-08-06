#!/usr/bin/env python3

import argparse
import re
import csv
import os
import shutil
import subprocess
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

ISOLATE_IDS_TO_REMOVE = {
    "EPI_ISL_20075028": "MANUAL_EXCLUSION",
    "EPI_ISL_20075029": "MANUAL_EXCLUSION",
    "EPI_ISL_20075031": "MANUAL_EXCLUSION",
    "EPI_ISL_20075045": "MANUAL_EXCLUSION",
    "EPI_ISL_20075046": "MANUAL_EXCLUSION",
}

# Use all available CPUs. You can also set a fixed number such as 8.
MAFFT_THREADS = -1

def extract_isolate_id(header):
    """Extract the EPI_ISL identifier from a FASTA header."""
    match = re.search(
        r"(?<![A-Za-z0-9_])EPI_ISL_\d+"
        r"(?![A-Za-z0-9_])",
        header,
    )

    if match is None:
        return None

    return match.group(0)
    
def get_project_root():
    """
    Script expected location:
    project_root/02_msa/scripts/msa.py
    """
    return Path(__file__).resolve().parents[2]


def find_mafft():
    """Check whether MAFFT is installed."""
    mafft_path = shutil.which("mafft")

    if mafft_path is None:
        raise FileNotFoundError(
            "MAFFT was not found.\n"
            "Install it with:\n"
            "  brew install mafft"
        )

    result = subprocess.run(
        [mafft_path, "--version"],
        capture_output=True,
        text=True,
        check=False,
    )

    version = (
        result.stdout.strip()
        or result.stderr.strip()
        or "version unknown"
    )

    print(f"MAFFT found: {version}")

    return mafft_path


def count_fasta_sequences(fasta_file):
    """Count the number of sequences in a FASTA file."""
    count = 0

    with fasta_file.open(
        "r",
        encoding="utf-8",
    ) as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1

    return count
def write_fasta_record(handle, header, sequence):
    """Write a FASTA record with 80-character line wrapping."""
    handle.write(f">{header}\n")

    line_width = 80

    for position in range(
        0,
        len(sequence),
        line_width,
    ):
        handle.write(
            sequence[position:position + line_width]
            + "\n"
        )
def filter_valid_dna_sequences(
    input_fasta,
    valid_fasta,
    removed_report,
):
    """
    Perform two filtering steps:

    1. Remove sequences containing non-ATGC characters.
    2. Deduplicate by Isolate_Id, keeping the longest valid sequence for each ID.
    """
    valid_characters = set("ATGC")

    records = []
    current_header = None
    sequence_lines = []

    def save_current_record():
        if current_header is None:
            return

        sequence = "".join(sequence_lines)
        sequence = "".join(
            sequence.split()
        ).upper()

        records.append(
            {
                "header": current_header,
                "sequence": sequence,
            }
        )

    # Read all FASTA records
    with input_fasta.open(
        "r",
        encoding="utf-8",
    ) as input_handle:
        for raw_line in input_handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                save_current_record()

                current_header = line[1:].strip()
                sequence_lines = []
            else:
                if current_header is None:
                    raise ValueError(
                        "Invalid FASTA format: sequence "
                        "appears before the first header in "
                        f"{input_fasta}"
                    )

                sequence_lines.append(line)

        save_current_record()

    removed_records = []
    records_by_id = {}

    for record in records:
        header = record["header"]
        sequence = record["sequence"]

        isolate_id = extract_isolate_id(header)

        if isolate_id is None:
            removed_records.append(
                {
                    "Header": header,
                    "Isolate_Id": "",
                    "Reason": "MISSING_ISOLATE_ID",
                    "Invalid_Characters": "",
                    "Sequence_Length": len(sequence),
                }
            )
            continue
                    # Remove manually excluded Isolate_Id values
        if isolate_id in ISOLATE_IDS_TO_REMOVE:
            removed_records.append(
                {
                    "Header": header,
                    "Isolate_Id": isolate_id,
                    "Reason": ISOLATE_IDS_TO_REMOVE[
                        isolate_id
                    ],
                    "Invalid_Characters": "",
                    "Sequence_Length": len(sequence),
                }
            )

            print(
                f"  Removing manually excluded sample: "
                f"{isolate_id}"
            )

            continue

        if not sequence:
            removed_records.append(
                {
                    "Header": header,
                    "Isolate_Id": isolate_id,
                    "Reason": "EMPTY_SEQUENCE",
                    "Invalid_Characters": "",
                    "Sequence_Length": 0,
                }
            )
            continue

        invalid_characters = sorted(
            set(sequence) - valid_characters
        )

        if invalid_characters:
            removed_records.append(
                {
                    "Header": header,
                    "Isolate_Id": isolate_id,
                    "Reason": "INVALID_CHARACTERS",
                    "Invalid_Characters": "".join(
                        invalid_characters
                    ),
                    "Sequence_Length": len(sequence),
                }
            )
            continue

        # First time we see this ID, keep it directly
        if isolate_id not in records_by_id:
            records_by_id[isolate_id] = record
            continue

        existing_record = records_by_id[
            isolate_id
        ]

        # For duplicate IDs, keep the longer sequence
        if len(sequence) > len(
            existing_record["sequence"]
        ):
            removed_records.append(
                {
                    "Header": existing_record["header"],
                    "Isolate_Id": isolate_id,
                    "Reason": "DUPLICATE_ID_SHORTER_SEQUENCE",
                    "Invalid_Characters": "",
                    "Sequence_Length": len(
                        existing_record["sequence"]
                    ),
                }
            )

            records_by_id[isolate_id] = record
        else:
            removed_records.append(
                {
                    "Header": header,
                    "Isolate_Id": isolate_id,
                    "Reason": "DUPLICATE_ID_SHORTER_OR_EQUAL",
                    "Invalid_Characters": "",
                    "Sequence_Length": len(sequence),
                }
            )

    valid_fasta.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    # Write the filtered FASTA with unique IDs
    with valid_fasta.open(
        "w",
        encoding="utf-8",
    ) as output_handle:
        for record in records_by_id.values():
            write_fasta_record(
                handle=output_handle,
                header=record["header"],
                sequence=record["sequence"],
            )

    # Write the removal report
    with removed_report.open(
        "w",
        newline="",
        encoding="utf-8-sig",
    ) as report_handle:
        writer = csv.DictWriter(
            report_handle,
            fieldnames=[
                "Header",
                "Isolate_Id",
                "Reason",
                "Invalid_Characters",
                "Sequence_Length",
            ],
        )

        writer.writeheader()
        writer.writerows(removed_records)

    valid_count = len(records_by_id)
    removed_count = len(removed_records)

    if valid_count == 0:
        raise ValueError(
            "No valid unique sequences remain after "
            f"filtering: {input_fasta}"
        )

    return valid_count, removed_count

    with input_fasta.open(
        "r",
        encoding="utf-8",
    ) as input_handle, valid_fasta.open(
        "w",
        encoding="utf-8",
    ) as output_handle:
        for raw_line in input_handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                kept, removed = process_record(
                    header=current_header,
                    lines=sequence_lines,
                    output_handle=output_handle,
                )

                valid_count += kept
                removed_count += removed

                current_header = line[1:].strip()
                sequence_lines = []
            else:
                if current_header is None:
                    raise ValueError(
                        "Invalid FASTA format: sequence "
                        "appears before the first header in "
                        f"{input_fasta}"
                    )

                sequence_lines.append(line)

        # Process the last sequence
        kept, removed = process_record(
            header=current_header,
            lines=sequence_lines,
            output_handle=output_handle,
        )

        valid_count += kept
        removed_count += removed

    with removed_report.open(
        "w",
        newline="",
        encoding="utf-8-sig",
    ) as report_handle:
        writer = csv.DictWriter(
            report_handle,
            fieldnames=[
                "Header",
                "Invalid_Characters",
                "Sequence_Length",
            ],
        )

        writer.writeheader()
        writer.writerows(removed_records)

    if valid_count == 0:
        raise ValueError(
            "No valid ATGC-only sequences remain after "
            f"filtering: {input_fasta}"
        )

    return valid_count, removed_count

def validate_inputs(input_dir, segments=SEGMENTS):
    """Check all segment inputs before running MAFFT."""
    missing_files = []
    empty_files = []

    for segment in segments:
        input_fasta = (
            input_dir
            / segment
            / f"{segment}_reference_dedup.fasta"
        )

        if not input_fasta.exists():
            missing_files.append(input_fasta)
        elif input_fasta.stat().st_size == 0:
            empty_files.append(input_fasta)

    messages = []

    if missing_files:
        messages.append(
            "Missing FASTA files:\n"
            + "\n".join(
                f"  - {path}"
                for path in missing_files
            )
        )

    if empty_files:
        messages.append(
            "Empty FASTA files:\n"
            + "\n".join(
                f"  - {path}"
                for path in empty_files
            )
        )

    if messages:
        raise FileNotFoundError(
            "\n\n".join(messages)
        )


def run_mafft(
    mafft_path,
    input_fasta,
    output_fasta,
    log_file,
):
    """
    Perform multiple sequence alignment with MAFFT.

    --auto:
        Let MAFFT choose an algorithm based on dataset size.

    --thread:
        Set the CPU thread count.

    --inputorder:
        Preserve input sequence order as much as possible.
    """
    output_fasta.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary_output = output_fasta.with_suffix(
        output_fasta.suffix + ".tmp"
    )

    command = [
        mafft_path,
        "--auto",
        "--thread",
        str(MAFFT_THREADS),
        "--inputorder",
        str(input_fasta),
    ]

    try:
        with temporary_output.open(
            "w",
            encoding="utf-8",
        ) as output_handle, log_file.open(
            "w",
            encoding="utf-8",
        ) as log_handle:
            result = subprocess.run(
                command,
                stdout=output_handle,
                stderr=log_handle,
                text=True,
                check=False,
            )

        if result.returncode != 0:
            temporary_output.unlink(
                missing_ok=True
            )

            raise RuntimeError(
                "MAFFT failed.\n"
                f"Input: {input_fasta}\n"
                f"Log: {log_file}\n"
                f"Exit code: {result.returncode}"
            )

        if not temporary_output.exists():
            raise RuntimeError(
                f"MAFFT did not create output: "
                f"{temporary_output}"
            )

        if temporary_output.stat().st_size == 0:
            temporary_output.unlink(
                missing_ok=True
            )

            raise RuntimeError(
                f"MAFFT output is empty: {input_fasta}"
            )

        # Replace the final output only after MAFFT succeeds to avoid leaving partial files
        os.replace(
            temporary_output,
            output_fasta,
        )

    except Exception:
        temporary_output.unlink(
            missing_ok=True
        )
        raise


def process_segment(
    segment,
    mafft_path,
    input_dir,
    output_dir,
):
    """Clean invalid sequences and run MAFFT for one gene segment."""
    input_fasta = (
        input_dir
        / segment
        / f"{segment}_reference_dedup.fasta"
    )

    segment_output_dir = (
        output_dir / segment
    )

    segment_output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # Intermediate FASTA containing only ATGC-compliant sequences
    valid_fasta = (
        segment_output_dir
        / f"{segment}_reference_dedup_valid.fasta"
    )

    # Report of removed sequences
    removed_report = (
        segment_output_dir
        / f"{segment}_removed_invalid_sequences.csv"
    )

    # Final MAFFT output
    output_fasta = (
        segment_output_dir
        / f"{segment}_reference_dedup_aligned.fasta"
    )

    log_file = (
        segment_output_dir
        / f"{segment}_mafft.log"
    )

    original_count = count_fasta_sequences(
        input_fasta
    )

    print()
    print(f"Processing {segment}")
    print(f"  Original sequences: {original_count}")
    print(f"  Input: {input_fasta}")

    # Step 1: remove entire sequences containing non-ATGC characters
    valid_count, removed_count = (
        filter_valid_dna_sequences(
            input_fasta=input_fasta,
            valid_fasta=valid_fasta,
            removed_report=removed_report,
        )
    )

    print(f"  Valid ATGC sequences: {valid_count}")
    print(
    f"  Invalid or duplicate-ID sequences removed: "
    f"{removed_count}"
)
    print(f"  Cleaned FASTA: {valid_fasta}")
    print(f"  Removal report: {removed_report}")

    if original_count != valid_count + removed_count:
        raise RuntimeError(
            f"{segment}: filtering count verification "
            "failed.\n"
            f"Original: {original_count}\n"
            f"Valid: {valid_count}\n"
            f"Removed: {removed_count}"
        )

    # Step 2: run MAFFT on the cleaned FASTA
    run_mafft(
        mafft_path=mafft_path,
        input_fasta=valid_fasta,
        output_fasta=output_fasta,
        log_file=log_file,
    )

    aligned_count = count_fasta_sequences(
        output_fasta
    )

    if valid_count != aligned_count:
        raise RuntimeError(
            f"{segment}: sequence count changed during "
            "MAFFT alignment.\n"
            f"Valid input sequences: {valid_count}\n"
            f"Aligned sequences: {aligned_count}"
        )

    print(f"  Aligned sequences: {aligned_count}")
    print(f"  Alignment output: {output_fasta}")
    print(f"  MAFFT log: {log_file}")
    print(f"  {segment} completed")
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--segment", required=True, type=str.upper, choices=SEGMENTS)
    segment = parser.parse_args().segment
    project_root = get_project_root()

    input_dir = (
        project_root
        / "02_msa"
        / "input"
    )

    output_dir = (
        project_root
        / "02_msa"
        / "output"
    )

    if not input_dir.exists():
        raise FileNotFoundError(
            f"MSA input directory does not exist: "
            f"{input_dir}"
        )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    mafft_path = find_mafft()

    # Check all inputs once before execution
    validate_inputs(input_dir, (segment,))

    completed = 0

    process_segment(segment=segment, mafft_path=mafft_path, input_dir=input_dir, output_dir=output_dir)
    completed += 1

    print()
    print(
        f"MAFFT completed successfully for "
        f"{completed}/1 segment."
    )
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
