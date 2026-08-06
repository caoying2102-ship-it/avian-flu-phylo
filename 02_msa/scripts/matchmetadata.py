#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd


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
    project_root/02_msa/scripts/matchmetadata.py
    """
    return Path(__file__).resolve().parents[2]


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
    """
    fasta_ids = []
    headers_without_id = []

    with fasta_file.open(
        "r",
        encoding="utf-8",
    ) as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line.startswith(">"):
                continue

            header = line[1:].strip()
            isolate_id = extract_isolate_id(header)

            if isolate_id is None:
                headers_without_id.append(header)
            else:
                fasta_ids.append(isolate_id)

    if headers_without_id:
        examples = "\n".join(
            f"  - {header}"
            for header in headers_without_id[:10]
        )

        raise ValueError(
            f"{len(headers_without_id)} FASTA headers "
            "do not contain an EPI_ISL ID:\n"
            f"{examples}"
        )

    if not fasta_ids:
        raise ValueError(
            f"No EPI_ISL IDs were found in: {fasta_file}"
        )

    duplicated_ids = pd.Series(fasta_ids)[
        pd.Series(fasta_ids).duplicated()
    ].unique()

    if len(duplicated_ids) > 0:
        examples = ", ".join(
            duplicated_ids[:10]
        )

        raise ValueError(
            "Duplicate Isolate_Id values were found in "
            f"{fasta_file}: {examples}"
        )

    return fasta_ids


def match_metadata(
    fasta_file,
    metadata_file,
    output_metadata,
    missing_report,
):
    """Extract matching metadata rows according to Isolate_Id values in the FASTA headers."""
    fasta_ids = read_fasta_ids(fasta_file)

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
        metadata["Isolate_Id"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    # Check whether metadata contains duplicate IDs
    duplicated_metadata = metadata[
        metadata["_Clean_Isolate_Id"].duplicated(
            keep=False
        )
        & (metadata["_Clean_Isolate_Id"] != "")
    ]

    if not duplicated_metadata.empty:
        duplicate_count = (
            duplicated_metadata["_Clean_Isolate_Id"]
            .nunique()
        )

        print(
            f"  WARNING: Metadata contains "
            f"{duplicate_count} duplicated Isolate_Id "
            "values. The first row will be retained."
        )

        metadata = metadata.drop_duplicates(
            subset="_Clean_Isolate_Id",
            keep="first",
        )

    # Build a FASTA-order table
    fasta_order = pd.DataFrame(
        {
            "_Clean_Isolate_Id": fasta_ids,
            "_FASTA_Order": range(len(fasta_ids)),
        }
    )

    # Use a left merge so every FASTA ID can be checked
    matched_metadata = fasta_order.merge(
        metadata,
        on="_Clean_Isolate_Id",
        how="left",
        indicator=True,
        sort=False,
    )

    missing_ids = matched_metadata.loc[
        matched_metadata["_merge"] == "left_only",
        "_Clean_Isolate_Id",
    ].tolist()

    # Save FASTA IDs that are missing from metadata
    pd.DataFrame(
        {
            "Isolate_Id": missing_ids,
        }
    ).to_csv(
        missing_report,
        index=False,
        encoding="utf-8-sig",
    )

    if missing_ids:
        examples = "\n".join(
            f"  - {isolate_id}"
            for isolate_id in missing_ids[:10]
        )

        raise ValueError(
            f"{len(missing_ids)} FASTA Isolate_Id values "
            "were not found in metadata.\n"
            f"Examples:\n{examples}\n"
            f"Full report: {missing_report}"
        )

    matched_metadata = matched_metadata.sort_values(
        by="_FASTA_Order",
        kind="stable",
    )

    matched_metadata = matched_metadata.drop(
        columns=[
            "_Clean_Isolate_Id",
            "_FASTA_Order",
            "_merge",
        ],
        errors="ignore",
    )

    matched_metadata.to_excel(
        output_metadata,
        index=False,
    )

    if len(matched_metadata) != len(fasta_ids):
        raise RuntimeError(
            "FASTA and metadata counts do not match.\n"
            f"FASTA: {len(fasta_ids)}\n"
            f"Metadata: {len(matched_metadata)}"
        )

    return len(fasta_ids), len(matched_metadata)


def process_segment(
    segment,
    msa_input_dir,
    msa_output_dir,
):
    """Process one gene segment."""
    fasta_file = (
        msa_output_dir
        / segment
        / f"{segment}_reference_dedup_aligned.fasta"
    )

    metadata_file = (
        msa_input_dir
        / segment
        / f"{segment}_metadata_reference_dedup.xlsx"
    )

    output_metadata = (
        msa_output_dir
        / segment
        / f"{segment}_metadata_reference_dedup_aligned.xlsx"
    )

    missing_report = (
        msa_output_dir
        / segment
        / f"{segment}_missing_metadata_IDs.csv"
    )

    print()
    print(f"Processing {segment}")

    if not fasta_file.exists():
        raise FileNotFoundError(
            f"Aligned FASTA file not found: {fasta_file}"
        )

    if not metadata_file.exists():
        raise FileNotFoundError(
            f"Metadata file not found: {metadata_file}"
        )

    fasta_count, metadata_count = match_metadata(
        fasta_file=fasta_file,
        metadata_file=metadata_file,
        output_metadata=output_metadata,
        missing_report=missing_report,
    )

    print(f"  FASTA sequences: {fasta_count}")
    print(f"  Matched metadata rows: {metadata_count}")
    print(f"  Output: {output_metadata}")
    print(f"  {segment} completed")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--segment", required=True, type=str.upper, choices=SEGMENTS)
    segment = parser.parse_args().segment
    project_root = get_project_root()

    msa_input_dir = (
        project_root
        / "02_msa"
        / "input"
    )

    msa_output_dir = (
        project_root
        / "02_msa"
        / "output"
    )

    if not msa_input_dir.exists():
        raise FileNotFoundError(
            f"MSA input directory does not exist: "
            f"{msa_input_dir}"
        )

    msa_output_dir.mkdir(parents=True, exist_ok=True)

    completed = 0

    process_segment(segment=segment, msa_input_dir=msa_input_dir, msa_output_dir=msa_output_dir)
    completed += 1

    print()
    print(
        f"Metadata matching completed for "
        f"{completed}/1 segment."
    )


if __name__ == "__main__":
    main()
