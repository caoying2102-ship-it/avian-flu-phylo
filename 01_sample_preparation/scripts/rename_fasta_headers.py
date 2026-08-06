#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO


SEGMENTS = ("MP", "HA", "NA", "PB2", "PB1", "PA", "NP", "NS")


def remove_all_whitespace(value, default):
    if pd.isna(value) or str(value).strip() == "":
        return default

    return re.sub(r"\s+", "", str(value))


def rename_fasta_headers(segment_dir, segment):
    metadata_path = (
        segment_dir
        / f"{segment}_metadata_simplified.xlsx"
    )

    input_fasta = (
        segment_dir
        / f"{segment}_dedup.fasta"
    )

    output_fasta = (
        segment_dir
        / f"{segment}_dedup_renamed.fasta"
    )

    if not metadata_path.exists():
        print(f"WARNING: metadata not found: {metadata_path}")
        return

    if not input_fasta.exists():
        print(f"WARNING: FASTA not found: {input_fasta}")
        return

    df = pd.read_excel(
        metadata_path,
        dtype={"Isolate_Id": str},
    )

    df["Isolate_Id"] = df["Isolate_Id"].map(
        lambda value: remove_all_whitespace(
            value,
            "",
        )
    )

    meta_dict = (
        df.drop_duplicates(
            subset="Isolate_Id",
            keep="first",
        )
        .set_index("Isolate_Id")
        .to_dict("index")
    )

    renamed_count = 0
    skipped_count = 0

    with open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            parts = record.description.split("|")

            isolate_id = (
                remove_all_whitespace(parts[3], "")
                if len(parts) > 3
                else None
            )

            if isolate_id and isolate_id in meta_dict:
                info = meta_dict[isolate_id]

                isolate_name = remove_all_whitespace(
                    info.get("Isolate_Name"),
                    "Unknown_Isolate",
                )

                collection_date = remove_all_whitespace(
                    info.get("Collection_Date"),
                    "Unknown_Date",
                )

                subcontinent = remove_all_whitespace(
                    info.get("Subcontinent"),
                    "Other/Unknown",
                )

                new_header = (
                    f"{isolate_name}|"
                    f"{collection_date}|"
                    f"{isolate_id}|"
                    f"{subcontinent}"
                )

                new_header = re.sub(
                    r"\s+",
                    "",
                    new_header,
                )

                record.id = new_header
                record.name = new_header
                record.description = ""

                SeqIO.write(
                    record,
                    output_handle,
                    "fasta",
                )

                renamed_count += 1

            else:
                skipped_count += 1

    print(
        f"{segment}: renamed={renamed_count}, "
        f"skipped={skipped_count}, "
        f"output={output_fasta}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--segment", required=True, type=str.upper, choices=SEGMENTS)
    segment = parser.parse_args().segment
    project_root = Path(__file__).resolve().parents[2]

    genes_dir = (
        project_root
        / "01_sample_preparation"
        / "output"
        / "genes"
    )

    rename_fasta_headers(segment_dir=genes_dir / segment, segment=segment)


if __name__ == "__main__":
    main()
