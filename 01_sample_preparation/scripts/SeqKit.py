#!/usr/bin/env python3

import argparse
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

GENES_DIR = (
    PROJECT_ROOT
    / "01_sample_preparation"
    / "output"
    / "genes"
)

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


def check_seqkit():
    """Check whether SeqKit is installed."""
    seqkit_path = shutil.which("seqkit")

    if seqkit_path is None:
        raise FileNotFoundError(
            "SeqKit is not installed or is not available "
            "in PATH.\n"
            "Install it with:\n"
            "brew install seqkit"
        )

    result = subprocess.run(
        [seqkit_path, "version"],
        capture_output=True,
        text=True,
        check=True,
    )

    version_text = (
        result.stdout.strip()
        or result.stderr.strip()
    )

    print(f"SeqKit found: {version_text}")

    return seqkit_path


def extract_isolate_id(header):
    """
    Extract the EPI_ISL identifier from a FASTA header.

    Example:
    A/fox/England/015850/2022|2022-08-23|
    EPI_ISL_17072388|Europe
    """
    match = re.search(
        r"(?<![A-Za-z0-9_])EPI_ISL_\d+"
        r"(?![A-Za-z0-9_])",
        header,
    )

    if match is None:
        return None

    return match.group(0)


def read_fasta_ids(fasta_file):
    """
    Read Isolate_Id values in FASTA order.

    Returns:
        FASTA sequence count
        list of Isolate_Id values
        list of headers that could not be recognized
    """
    isolate_ids = []
    missing_id_headers = []
    sequence_count = 0

    with fasta_file.open(
        "r",
        encoding="utf-8",
    ) as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line.startswith(">"):
                continue

            sequence_count += 1
            header = line[1:].strip()

            isolate_id = extract_isolate_id(header)

            if isolate_id is None:
                missing_id_headers.append(header)
            else:
                isolate_ids.append(isolate_id)

    return (
        sequence_count,
        isolate_ids,
        missing_id_headers,
    )


def run_seqkit_rmdup(
    seqkit_path,
    input_fasta,
    output_fasta,
):
    """
    Deduplicate based on sequence content.

    -s: deduplicate based on sequence content rather than header names
    -i: ignore sequence case
    """
    command = [
        seqkit_path,
        "rmdup",
        "-s",
        "-i",
        str(input_fasta),
        "-o",
        str(output_fasta),
    ]

    print("  Running: " + " ".join(command))

    subprocess.run(
        command,
        check=True,
    )


def generate_deduplicated_metadata(
    input_metadata,
    output_metadata,
    retained_ids,
):
    """
    Extract and reorder the corresponding metadata based on the Isolate_Id values retained in the deduplicated FASTA.
    """
    metadata = pd.read_excel(
        input_metadata,
        dtype={"Isolate_Id": str},
    )

    if "Isolate_Id" not in metadata.columns:
        raise ValueError(
            f"Isolate_Id column is missing from: "
            f"{input_metadata}"
        )

    metadata = metadata.copy()

    metadata["_Clean_Isolate_Id"] = (
        metadata["Isolate_Id"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    # Save the order from the deduplicated FASTA
    retained_order = {
        isolate_id: position
        for position, isolate_id
        in enumerate(retained_ids)
    }

    selected_metadata = metadata[
        metadata["_Clean_Isolate_Id"].isin(
            retained_order
        )
    ].copy()

    selected_metadata["_FASTA_Order"] = (
        selected_metadata["_Clean_Isolate_Id"].map(
            retained_order
        )
    )

    selected_metadata = selected_metadata.sort_values(
        by="_FASTA_Order",
        kind="stable",
    )

    metadata_ids = set(
        selected_metadata["_Clean_Isolate_Id"]
    )

    missing_metadata_ids = [
        isolate_id
        for isolate_id in retained_ids
        if isolate_id not in metadata_ids
    ]

    selected_metadata = selected_metadata.drop(
        columns=[
            "_Clean_Isolate_Id",
            "_FASTA_Order",
        ],
        errors="ignore",
    )

    selected_metadata.to_excel(
        output_metadata,
        index=False,
    )

    return len(selected_metadata), missing_metadata_ids


def process_segment(segment, seqkit_path):
    """Process one gene segment."""
    segment_dir = GENES_DIR / segment

    input_fasta = (
        segment_dir
        / f"{segment}_reference.fasta"
    )

    input_metadata = (
        segment_dir
        / f"{segment}_metadata_reference.xlsx"
    )

    output_fasta = (
        segment_dir
        / f"{segment}_reference_dedup.fasta"
    )

    output_metadata = (
        segment_dir
        / f"{segment}_metadata_reference_dedup.xlsx"
    )

    print(f"\nProcessing {segment}...")

    if not input_fasta.exists():
        print(
            f"  WARNING: FASTA file not found: "
            f"{input_fasta}"
        )
        return

    if not input_metadata.exists():
        print(
            f"  WARNING: Metadata file not found: "
            f"{input_metadata}"
        )
        return

    (
        original_sequence_count,
        _,
        original_missing_headers,
    ) = read_fasta_ids(input_fasta)

    run_seqkit_rmdup(
        seqkit_path=seqkit_path,
        input_fasta=input_fasta,
        output_fasta=output_fasta,
    )

    (
        deduplicated_sequence_count,
        retained_ids,
        deduplicated_missing_headers,
    ) = read_fasta_ids(output_fasta)

    metadata_count, missing_metadata_ids = (
        generate_deduplicated_metadata(
            input_metadata=input_metadata,
            output_metadata=output_metadata,
            retained_ids=retained_ids,
        )
    )

    removed_count = (
        original_sequence_count
        - deduplicated_sequence_count
    )

    print(
        f"  Original FASTA sequences: "
        f"{original_sequence_count}"
    )
    print(
        f"  Deduplicated FASTA sequences: "
        f"{deduplicated_sequence_count}"
    )
    print(
        f"  Duplicate sequences removed: "
        f"{removed_count}"
    )
    print(
        f"  Deduplicated metadata rows: "
        f"{metadata_count}"
    )
    print(f"  FASTA output: {output_fasta}")
    print(f"  Metadata output: {output_metadata}")

    if original_missing_headers:
        print(
            "  WARNING: Headers without EPI_ISL ID "
            "in the input FASTA: "
            f"{len(original_missing_headers)}"
        )

    if deduplicated_missing_headers:
        print(
            "  WARNING: Headers without EPI_ISL ID "
            "in the deduplicated FASTA: "
            f"{len(deduplicated_missing_headers)}"
        )

    if missing_metadata_ids:
        print(
            "  WARNING: Retained FASTA IDs missing "
            "from metadata: "
            f"{len(missing_metadata_ids)}"
        )

        for isolate_id in missing_metadata_ids[:10]:
            print(f"    {isolate_id}")

        if len(missing_metadata_ids) > 10:
            print("    ...")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--segment", required=True, type=str.upper, choices=SEGMENTS)
    segment = parser.parse_args().segment
    if not GENES_DIR.exists():
        raise FileNotFoundError(
            f"Genes directory does not exist: "
            f"{GENES_DIR}"
        )

    seqkit_path = check_seqkit()

    process_segment(segment=segment, seqkit_path=seqkit_path)

    print("\nSeqKit deduplication completed.")


if __name__ == "__main__":
    main()
