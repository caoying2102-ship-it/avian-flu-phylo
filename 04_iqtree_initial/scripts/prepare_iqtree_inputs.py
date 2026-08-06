#!/usr/bin/env python3

import argparse
import shutil
from pathlib import Path


FASTA_SUFFIXES = {".fa", ".fas", ".fasta", ".fna", ".faa", ".aln"}


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
    Script location:

    project_root/
    └── 04_iqtree_initial/
        └── scripts/
            └── prepare_iqtree_inputs.py
    """
    return Path(__file__).resolve().parents[2]


def resolve_segments(requested_segment):
    """Resolve a CLI segment selection to a tuple of segments."""
    if requested_segment is None:
        return SEGMENTS
    normalized = requested_segment.upper()
    if normalized == "ALL":
        return SEGMENTS
    if normalized in SEGMENTS:
        return (normalized,)
    raise ValueError(f"Unsupported segment: {requested_segment}")


def build_file_pairs(project_root, segments=SEGMENTS):
    """Create the mapping between TrimAl outputs and IQ-TREE inputs."""
    trimal_output_dir = (
        project_root
        / "03_trimal"
        / "output"
    )

    iqtree_input_dir = (
        project_root
        / "04_iqtree_initial"
        / "input"
    )

    file_pairs = []

    for segment in segments:
        source_dir = (
            trimal_output_dir / segment
        )

        destination_dir = (
            iqtree_input_dir / segment
        )

        source_fasta = (
            source_dir
            / f"{segment}_reference_dedup_aligned_trimmed.fasta"
        )

        source_metadata = (
            source_dir
            / f"{segment}_metadata_reference_dedup_aligned_trimmed.xlsx"
        )

        destination_fasta = (
            destination_dir / source_fasta.name
        )

        destination_metadata = (
            destination_dir / source_metadata.name
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


def validate_source_files(file_pairs, strict=True):
    """Check all TrimAl output files before copying."""
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

    if messages and strict:
        raise FileNotFoundError(
            "\n\n".join(messages)
            + "\n\nRun the TrimAl step first."
        )

    return messages


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

    if (
        source_file.stat().st_size
        != destination_file.stat().st_size
    ):
        raise RuntimeError(
            "File-size verification failed:\n"
            f"Source: {source_file}\n"
            f"Destination: {destination_file}"
        )


def deduplicate_fasta(source_file, destination_file):
    """Write a FASTA file that keeps only the first occurrence of each header."""
    seen_headers = set()
    kept_records = []

    with source_file.open("r", encoding="utf-8") as infile:
        current_header = None
        current_sequence = []

        def flush_current_record():
            nonlocal current_header, current_sequence
            if current_header is None:
                return
            if current_header not in seen_headers:
                seen_headers.add(current_header)
                kept_records.append(f">{current_header}\n")
                kept_records.extend(current_sequence)
                if kept_records[-1] != "\n":
                    kept_records.append("\n")
            current_header = None
            current_sequence = []

        for raw_line in infile:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush_current_record()
                current_header = line[1:].strip()
                current_sequence = []
            else:
                if current_header is None:
                    continue
                current_sequence.append(line + "\n")

        flush_current_record()

    destination_file.parent.mkdir(parents=True, exist_ok=True)
    with destination_file.open("w", encoding="utf-8") as outfile:
        outfile.write("".join(kept_records))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--segment",
        default="ALL",
        help="Segment to prepare; use ALL (default) to process every segment.",
    )
    requested_segment = parser.parse_args().segment
    segments = resolve_segments(requested_segment)
    project_root = get_project_root()

    trimal_output_dir = (
        project_root
        / "03_trimal"
        / "output"
    )

    iqtree_dir = (
        project_root
        / "04_iqtree_initial"
    )

    if not trimal_output_dir.exists():
        raise FileNotFoundError(
            f"TrimAl output directory does not exist: "
            f"{trimal_output_dir}"
        )

    if not iqtree_dir.exists():
        raise FileNotFoundError(
            f"IQ-TREE directory does not exist: "
            f"{iqtree_dir}"
        )

    file_pairs = build_file_pairs(project_root, segments)

    # Begin copying only after all files pass validation; in ALL mode, missing segments are skipped.
    strict = len(segments) == 1
    validate_source_files(file_pairs, strict=strict)

    copied_count = 0
    skipped_count = 0

    for (
        segment,
        file_type,
        source_file,
        destination_file,
    ) in file_pairs:
        if not source_file.exists() or not source_file.is_file():
            skipped_count += 1
            print(
                f"Skipping {segment} {file_type}: source missing: {source_file}"
            )
            continue

        if file_type == "FASTA":
            deduplicate_fasta(
                source_file=source_file,
                destination_file=destination_file,
            )
        else:
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
        f"for {len(segments)} segment(s)."
    )
    if skipped_count:
        print(f"Skipped: {skipped_count} file(s) due to missing source inputs.")
    print(
        f"IQ-TREE input directory: "
        f"{iqtree_dir / 'input'}"
    )


if __name__ == "__main__":
    main()
