#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

REFERENCE_ID_FILE = (
    PROJECT_ROOT
    / "Data"
    / "reference ID.csv"
)

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


def clean_isolate_id(value):
    """Clean spaces and Excel-empty values from Isolate_Id."""
    if pd.isna(value):
        return ""

    return str(value).strip()


def read_reference_ids(reference_file):
    """Read Isolate_Id values from reference ID.csv."""
    if not reference_file.exists():
        raise FileNotFoundError(
            f"Reference ID file does not exist: "
            f"{reference_file}"
        )

    reference_table = pd.read_csv(
        reference_file,
        dtype=str,
        encoding="utf-8-sig",
    )

    if reference_table.empty:
        raise ValueError(
            f"Reference ID file is empty: {reference_file}"
        )

    # Prefer the Isolate_Id column; otherwise fall back to the first column
    if "Isolate_Id" in reference_table.columns:
        id_column = "Isolate_Id"
    else:
        id_column = reference_table.columns[0]

        print(
            "WARNING: Isolate_Id column was not found in "
            "reference ID.csv."
        )
        print(
            f"The first column will be used: {id_column}"
        )

    reference_ids = (
        reference_table[id_column]
        .map(clean_isolate_id)
        .loc[lambda values: values != ""]
        .drop_duplicates()
        .tolist()
    )

    if not reference_ids:
        raise ValueError(
            "No valid Isolate_Id values were found in "
            f"{reference_file}"
        )

    print(
        f"Loaded {len(reference_ids)} unique reference IDs."
    )

    return reference_ids


def extract_metadata(
    metadata_file,
    output_file,
    reference_ids,
):
    """Extract metadata based on Isolate_Id."""
    if not metadata_file.exists():
        print(
            f"WARNING: Metadata file not found: "
            f"{metadata_file}"
        )
        return set()

    metadata = pd.read_excel(
        metadata_file,
        dtype={"Isolate_Id": str},
    )

    if "Isolate_Id" not in metadata.columns:
        raise ValueError(
            f"Isolate_Id column is missing from: "
            f"{metadata_file}"
        )

    metadata = metadata.copy()

    metadata["_Clean_Isolate_Id"] = (
        metadata["Isolate_Id"].map(clean_isolate_id)
    )

    # Used to order rows according to reference ID.csv
    reference_order = {
        isolate_id: index
        for index, isolate_id in enumerate(reference_ids)
    }

    selected_metadata = metadata[
        metadata["_Clean_Isolate_Id"].isin(
            reference_order
        )
    ].copy()

    selected_metadata["_Reference_Order"] = (
        selected_metadata["_Clean_Isolate_Id"].map(
            reference_order
        )
    )

    selected_metadata = selected_metadata.sort_values(
        by="_Reference_Order",
        kind="stable",
    )

    found_ids = set(
        selected_metadata["_Clean_Isolate_Id"]
    )

    selected_metadata = selected_metadata.drop(
        columns=[
            "_Clean_Isolate_Id",
            "_Reference_Order",
        ],
        errors="ignore",
    )

    selected_metadata.to_excel(
        output_file,
        index=False,
    )

    return found_ids


def read_fasta_records(fasta_file):
    """
    Read a FASTA file.

    Returns:
        [(header, sequence), ...]
    """
    records = []
    current_header = None
    sequence_lines = []

    with fasta_file.open(
        "r",
        encoding="utf-8",
    ) as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    records.append(
                        (
                            current_header,
                            "".join(sequence_lines),
                        )
                    )

                current_header = line[1:].strip()
                sequence_lines = []
            else:
                if current_header is None:
                    raise ValueError(
                        "Invalid FASTA file: sequence found "
                        "before the first header in "
                        f"{fasta_file}"
                    )

                sequence_lines.append(line)

    if current_header is not None:
        records.append(
            (
                current_header,
                "".join(sequence_lines),
            )
        )

    return records


def get_isolate_id_from_header(header):
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


def extract_fasta(
    fasta_file,
    output_file,
    reference_ids,
):
    """Extract sequences from the FASTA headers based on Isolate_Id."""
    if not fasta_file.exists():
        print(
            f"WARNING: FASTA file not found: {fasta_file}"
        )
        return set()

    reference_id_set = set(reference_ids)
    records = read_fasta_records(fasta_file)

    # In theory, one Isolate_Id should map to a single deduplicated sequence.
    # A list is used to tolerate the case where the same ID appears multiple times unexpectedly.
    records_by_id = {}

    headers_without_id = 0

    for header, sequence in records:
        isolate_id = get_isolate_id_from_header(header)

        if isolate_id is None:
            headers_without_id += 1
            continue

        if isolate_id not in reference_id_set:
            continue

        records_by_id.setdefault(
            isolate_id,
            [],
        ).append(
            (header, sequence)
        )

    found_ids = set(records_by_id)

    # Write the FASTA in the order specified by reference ID.csv
    with output_file.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for isolate_id in reference_ids:
            for header, sequence in records_by_id.get(
                isolate_id,
                [],
            ):
                write_fasta_record(
                    handle=handle,
                    header=header,
                    sequence=sequence,
                )

    if headers_without_id:
        print(
            f"  FASTA headers without EPI_ISL ID: "
            f"{headers_without_id}"
        )

    return found_ids


def write_missing_ids(
    segment,
    reference_ids,
    metadata_ids,
    fasta_ids,
    output_file,
):
    """Write unmatched IDs for each gene to aid inspection."""
    rows = []

    for isolate_id in reference_ids:
        in_metadata = isolate_id in metadata_ids
        in_fasta = isolate_id in fasta_ids

        if not in_metadata or not in_fasta:
            rows.append(
                {
                    "Segment": segment,
                    "Isolate_Id": isolate_id,
                    "Found_in_Metadata": in_metadata,
                    "Found_in_FASTA": in_fasta,
                }
            )

    missing_table = pd.DataFrame(
        rows,
        columns=[
            "Segment",
            "Isolate_Id",
            "Found_in_Metadata",
            "Found_in_FASTA",
        ],
    )

    missing_table.to_csv(
        output_file,
        index=False,
        encoding="utf-8-sig",
    )

    return len(missing_table)


def process_segment(segment, reference_ids):
    """Process one gene segment."""
    segment_dir = GENES_DIR / segment

    metadata_file = (
        segment_dir
        / f"{segment}_metadata_simplified.xlsx"
    )

    fasta_file = (
        segment_dir
        / f"{segment}_dedup_renamed.fasta"
    )

    output_metadata_file = (
        segment_dir
        / f"{segment}_metadata_reference.xlsx"
    )

    output_fasta_file = (
        segment_dir
        / f"{segment}_reference.fasta"
    )

    missing_ids_file = (
        segment_dir
        / f"{segment}_missing_reference_IDs.csv"
    )

    if not segment_dir.exists():
        print(
            f"WARNING: Segment folder not found: "
            f"{segment_dir}"
        )
        return

    print(f"\nProcessing {segment}...")

    metadata_ids = extract_metadata(
        metadata_file=metadata_file,
        output_file=output_metadata_file,
        reference_ids=reference_ids,
    )

    fasta_ids = extract_fasta(
        fasta_file=fasta_file,
        output_file=output_fasta_file,
        reference_ids=reference_ids,
    )

    missing_count = write_missing_ids(
        segment=segment,
        reference_ids=reference_ids,
        metadata_ids=metadata_ids,
        fasta_ids=fasta_ids,
        output_file=missing_ids_file,
    )

    print(
        f"  Metadata samples extracted: "
        f"{len(metadata_ids)}"
    )
    print(
        f"  FASTA samples extracted: "
        f"{len(fasta_ids)}"
    )
    print(
        f"  Reference IDs with missing data: "
        f"{missing_count}"
    )
    print(f"  Metadata output: {output_metadata_file}")
    print(f"  FASTA output: {output_fasta_file}")
    print(f"  Missing-ID report: {missing_ids_file}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--segment", required=True, type=str.upper, choices=SEGMENTS)
    segment = parser.parse_args().segment
    if not GENES_DIR.exists():
        raise FileNotFoundError(
            f"Genes directory does not exist: {GENES_DIR}"
        )

    reference_ids = read_reference_ids(
        REFERENCE_ID_FILE
    )

    process_segment(segment=segment, reference_ids=reference_ids)

    print("\nExtraction completed.")


if __name__ == "__main__":
    main()
