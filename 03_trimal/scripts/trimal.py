#!/usr/bin/env python3

import argparse
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


def get_project_root():
    """
    Script expected location:

    project_root/
    └── 03_trimal/
        └── scripts/
            └── trimal.py
    """
    return Path(__file__).resolve().parents[2]


def find_trimal():
    """Check whether TrimAl is installed."""
    trimal_path = shutil.which("trimal")

    if trimal_path is None:
        raise FileNotFoundError(
            "TrimAl was not found.\n"
            "Install it with:\n"
            "  conda install -c bioconda trimal"
        )

    result = subprocess.run(
        [trimal_path, "--version"],
        capture_output=True,
        text=True,
        check=False,
    )

    version = (
        result.stdout.strip()
        or result.stderr.strip()
        or "version unknown"
    )

    print(f"TrimAl found: {version}")

    return trimal_path


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


def read_alignment_length(fasta_file):
    """
    Check that all sequences have the same length and return the alignment length.
    """
    sequence_lengths = []
    current_sequence = []

    with fasta_file.open(
        "r",
        encoding="utf-8",
    ) as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_sequence:
                    sequence_lengths.append(
                        len("".join(current_sequence))
                    )

                current_sequence = []
            else:
                current_sequence.append(line)

        if current_sequence:
            sequence_lengths.append(
                len("".join(current_sequence))
            )

    if not sequence_lengths:
        raise ValueError(
            f"No sequences found in: {fasta_file}"
        )

    unique_lengths = set(sequence_lengths)

    if len(unique_lengths) != 1:
        raise ValueError(
            "Input FASTA is not a valid alignment because "
            "the sequences have different lengths:\n"
            f"{fasta_file}"
        )

    return sequence_lengths[0]


def validate_inputs(input_dir, segments=SEGMENTS):
    """
    Check all input files once before running TrimAl.
    """
    missing_files = []
    empty_files = []

    for segment in segments:
        input_fasta = (
            input_dir
            / segment
            / f"{segment}_reference_dedup_aligned.fasta"
        )

        input_metadata = (
            input_dir
            / segment
            / f"{segment}_metadata_reference_dedup_aligned.xlsx"
        )

        for input_file in (
            input_fasta,
            input_metadata,
        ):
            if not input_file.exists():
                missing_files.append(input_file)
            elif input_file.stat().st_size == 0:
                empty_files.append(input_file)

    messages = []

    if missing_files:
        messages.append(
            "Missing input files:\n"
            + "\n".join(
                f"  - {path}"
                for path in missing_files
            )
        )

    if empty_files:
        messages.append(
            "Empty input files:\n"
            + "\n".join(
                f"  - {path}"
                for path in empty_files
            )
        )

    if messages:
        raise FileNotFoundError(
            "\n\n".join(messages)
        )


def run_trimal(
    trimal_path,
    input_fasta,
    output_fasta,
    log_file,
):
    """
    Trim the alignment using TrimAl's automated1 mode.
    """
    output_fasta.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary_output = output_fasta.with_suffix(
        output_fasta.suffix + ".tmp"
    )

    command = [
        trimal_path,
        "-in",
        str(input_fasta),
        "-out",
        str(temporary_output),
        "-automated1",
    ]

    print(
        "  Running: "
        + " ".join(command)
    )

    try:
        with log_file.open(
            "w",
            encoding="utf-8",
        ) as log_handle:
            result = subprocess.run(
                command,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )

        if result.returncode != 0:
            temporary_output.unlink(
                missing_ok=True
            )

            raise RuntimeError(
                "TrimAl failed.\n"
                f"Input: {input_fasta}\n"
                f"Log: {log_file}\n"
                f"Exit code: {result.returncode}"
            )

        if not temporary_output.exists():
            raise RuntimeError(
                "TrimAl did not create an output file:\n"
                f"{temporary_output}"
            )

        if temporary_output.stat().st_size == 0:
            temporary_output.unlink(
                missing_ok=True
            )

            raise RuntimeError(
                f"TrimAl output is empty: {input_fasta}"
            )

        # Replace the final output only after success
        os.replace(
            temporary_output,
            output_fasta,
        )

    except Exception:
        temporary_output.unlink(
            missing_ok=True
        )
        raise


def validate_output_consistency(
    segment,
    input_sequence_count,
    output_sequence_count,
    input_alignment_length,
    output_alignment_length,
):
    """Validate trimAl output while allowing removal of gap-only sequences."""
    if output_alignment_length > input_alignment_length:
        raise RuntimeError(
            f"{segment}: alignment length increased after TrimAl.\n"
            f"Before: {input_alignment_length}\n"
            f"After: {output_alignment_length}"
        )

    if input_sequence_count != output_sequence_count:
        print(
            f"Warning: {segment} sequence count changed after TrimAl.\n"
            f"Before: {input_sequence_count}\n"
            f"After: {output_sequence_count}\n"
            "This can happen when trimAl removes gap-only sequences; proceeding anyway.",
            flush=True,
        )


def process_segment(
    segment,
    trimal_path,
    input_dir,
    output_dir,
):
    """Process one gene segment."""
    segment_input_dir = (
        input_dir / segment
    )

    segment_output_dir = (
        output_dir / segment
    )

    input_fasta = (
        segment_input_dir
        / f"{segment}_reference_dedup_aligned.fasta"
    )

    input_metadata = (
        segment_input_dir
        / f"{segment}_metadata_reference_dedup_aligned.xlsx"
    )

    output_fasta = (
        segment_output_dir
        / f"{segment}_reference_dedup_aligned_trimmed.fasta"
    )

    output_metadata = (
        segment_output_dir
        / f"{segment}_metadata_reference_dedup_aligned_trimmed.xlsx"
    )

    log_file = (
        segment_output_dir
        / f"{segment}_trimal.log"
    )

    segment_output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    input_sequence_count = count_fasta_sequences(
        input_fasta
    )

    input_alignment_length = read_alignment_length(
        input_fasta
    )

    print()
    print(f"Processing {segment}")
    print(
        f"  Input sequences: "
        f"{input_sequence_count}"
    )
    print(
        f"  Input alignment length: "
        f"{input_alignment_length}"
    )
    print(f"  Input: {input_fasta}")

    run_trimal(
        trimal_path=trimal_path,
        input_fasta=input_fasta,
        output_fasta=output_fasta,
        log_file=log_file,
    )

    output_sequence_count = count_fasta_sequences(
        output_fasta
    )

    output_alignment_length = read_alignment_length(
        output_fasta
    )

    validate_output_consistency(
        segment=segment,
        input_sequence_count=input_sequence_count,
        output_sequence_count=output_sequence_count,
        input_alignment_length=input_alignment_length,
        output_alignment_length=output_alignment_length,
    )

    # TrimAl trims columns but does not remove sequences, so metadata can be copied directly
    shutil.copy2(
        input_metadata,
        output_metadata,
    )

    print(
        f"  Output sequences: "
        f"{output_sequence_count}"
    )
    print(
        f"  Output alignment length: "
        f"{output_alignment_length}"
    )
    print(
        f"  Columns removed: "
        f"{input_alignment_length - output_alignment_length}"
    )
    print(f"  FASTA output: {output_fasta}")
    print(f"  Metadata output: {output_metadata}")
    print(f"  Log: {log_file}")
    print(f"  {segment} completed")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--segment", required=True, type=str.upper, choices=SEGMENTS)
    segment = parser.parse_args().segment
    project_root = get_project_root()

    input_dir = (
        project_root
        / "03_trimal"
        / "input"
    )

    output_dir = (
        project_root
        / "03_trimal"
        / "output"
    )

    if not input_dir.exists():
        raise FileNotFoundError(
            f"TrimAl input directory does not exist: "
            f"{input_dir}"
        )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    trimal_path = find_trimal()

    # Check all inputs before starting
    validate_inputs(input_dir, (segment,))

    completed = 0

    process_segment(segment=segment, trimal_path=trimal_path, input_dir=input_dir, output_dir=output_dir)
    completed += 1

    print()
    print(
        f"TrimAl completed successfully for "
        f"{completed}/1 segment."
    )
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
