#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def reclassify_host(host):
    host = str(host).lower().strip()

    if host == "human":
        return "Human"
    if host == "dairy cow":
        return "Dairy cow"
    if host == "chicken":
        return "Chicken"
    if host == "duck":
        return "Duck"
    if host == "goose":
        return "Goose"

    mammals = [
        "bovine",
        "canine",
        "equine",
        "feline",
        "ferret",
        "mink",
        "mouse",
        "rodent",
        "seal",
        "swine",
        "mammals",
        "other mammals",
    ]

    if host in mammals:
        return "Other Animal"

    return "Wild Bird"


SEGMENTS = ("HA", "NA", "MP", "PB2", "PB1", "PA", "NP", "NS")


def generate_summary(input_file, output_file):
    metadata = pd.read_excel(input_file)

    required_columns = [
        "Subcontinent",
        "Collection_Date",
        "Host",
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in metadata.columns
    ]

    if missing_columns:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(missing_columns)
        )

    # Convert dates such as 2022-06-28 into the year 2022
    collection_dates = pd.to_datetime(
        metadata["Collection_Date"],
        errors="coerce",
    )

    metadata["Year"] = collection_dates.dt.year.astype("Int64")

    # Unrecognized dates are classified as Unknown
    metadata["Year"] = metadata["Year"].astype("string")
    metadata["Year"] = metadata["Year"].fillna("Unknown")

    # Reclassify host values
    metadata["Host_Category"] = metadata["Host"].map(
        reclassify_host
    )

    # Empty Subcontinent values are classified as Other/Unknown
    metadata["Subcontinent"] = (
        metadata["Subcontinent"]
        .fillna("Other/Unknown")
        .astype(str)
        .str.strip()
        .replace("", "Other/Unknown")
    )

    # Count rows by Subcontinent, year, and Host
    summary = (
        metadata.groupby(
            [
                "Subcontinent",
                "Year",
                "Host_Category",
            ],
            dropna=False,
        )
        .size()
        .reset_index(name="Count")
    )

    summary = summary.sort_values(
        by=[
            "Subcontinent",
            "Year",
            "Host_Category",
        ]
    )

    summary = summary.rename(
        columns={
            "Year": "Collection_Year",
            "Host_Category": "Host",
        }
    )

    summary.to_excel(
        output_file,
        index=False,
    )

    print(f"Input rows: {len(metadata)}")
    print(f"Summary rows: {len(summary)}")
    print(f"Output file: {output_file}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--segment", required=True, type=str.upper, choices=SEGMENTS)
    segment = parser.parse_args().segment
    # Assume this script is located at:
    # 01_sample_preparation/scripts/summary_count.py
    project_root = Path(__file__).resolve().parents[2]

    segment_dir = (
        project_root
        / "01_sample_preparation"
        / "output"
        / "genes"
        / segment
    )

    input_file = segment_dir / f"{segment}_metadata_simplified.xlsx"
    output_file = segment_dir / f"{segment}_summary.xlsx"

    if not input_file.exists():
        raise FileNotFoundError(
            f"Input file does not exist: {input_file}"
        )

    generate_summary(
        input_file=input_file,
        output_file=output_file,
    )


if __name__ == "__main__":
    main()
