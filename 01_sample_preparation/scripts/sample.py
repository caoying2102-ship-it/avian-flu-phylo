#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


TARGET_SAMPLE_SIZE = 5200
RANDOM_SEED = 42

REQUIRED_ISOLATE_IDS = {
    "EPI_ISL_1254",
}

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


def prepare_metadata(metadata):
    required_columns = [
        "Collection_Date",
        "Subcontinent",
        "Host",
        "Isolate_Id",
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in metadata.columns
    ]

    if missing_columns:
        raise ValueError(
            "Simplified metadata is missing columns: "
            + ", ".join(missing_columns)
        )

    metadata = metadata.copy()

    # Convert dates such as 2022-06-28 into the year 2022
    parsed_dates = pd.to_datetime(
        metadata["Collection_Date"],
        errors="coerce",
    )

    metadata["_Collection_Year"] = (
        parsed_dates.dt.year.astype("Int64").astype("string")
    )

    metadata["_Collection_Year"] = metadata[
        "_Collection_Year"
    ].fillna("Unknown")

    metadata["_Subcontinent"] = (
        metadata["Subcontinent"]
        .fillna("Other/Unknown")
        .astype(str)
        .str.strip()
        .replace("", "Other/Unknown")
    )

    metadata["_Host_Category"] = metadata["Host"].map(
        reclassify_host
    )

    return metadata


def prepare_summary(summary):
    summary = summary.copy()

    # Support both Collection_Year and Collection_Date column names
    if "Collection_Year" in summary.columns:
        summary["_Collection_Year"] = (
            summary["Collection_Year"]
            .astype("string")
            .str.replace(r"\.0$", "", regex=True)
            .fillna("Unknown")
        )
    elif "Collection_Date" in summary.columns:
        parsed_dates = pd.to_datetime(
            summary["Collection_Date"],
            errors="coerce",
        )

        summary["_Collection_Year"] = (
            parsed_dates.dt.year.astype("Int64").astype("string")
        )

        summary["_Collection_Year"] = summary[
            "_Collection_Year"
        ].fillna("Unknown")
    else:
        raise ValueError(
            "Summary file must contain Collection_Year "
            "or Collection_Date"
        )

    required_columns = [
        "Subcontinent",
        "Host",
        "Count",
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in summary.columns
    ]

    if missing_columns:
        raise ValueError(
            "Summary file is missing columns: "
            + ", ".join(missing_columns)
        )

    summary["_Subcontinent"] = (
        summary["Subcontinent"]
        .fillna("Other/Unknown")
        .astype(str)
        .str.strip()
        .replace("", "Other/Unknown")
    )

    summary["_Host_Category"] = (
        summary["Host"]
        .fillna("Wild Bird")
        .astype(str)
        .str.strip()
    )

    summary["Count"] = pd.to_numeric(
        summary["Count"],
        errors="coerce",
    ).fillna(0).astype(int)

    return summary


def allocate_sample_counts(summary, target_size):
    """
    Allocate sample counts for each cell.

    Rule A:
        If Count <= 5, keep all samples.

    Rule B:
        If Count > 5, take at least 5 samples and distribute the remaining slots proportionally to cell capacity.
    """
    counts = summary["Count"].to_numpy(dtype=int)

    # Guaranteed minimum number of samples per cell
    minimum_samples = np.where(
        counts <= 5,
        counts,
        5,
    )

    total_available = int(counts.sum())
    minimum_total = int(minimum_samples.sum())

    actual_target = min(target_size, total_available)

    # If the guaranteed minimum for all cells already exceeds 5200, the floor rule must be applied first
    if minimum_total >= actual_target:
        sample_counts = minimum_samples.copy()

        print(
            "WARNING: Mandatory minimum samples exceed "
            f"the requested target ({target_size})."
        )
        print(
            "The minimum-rule sample size will be used: "
            f"{sample_counts.sum()}"
        )

        return sample_counts

    sample_counts = minimum_samples.copy()

    remaining_slots = actual_target - minimum_total
    remaining_capacity = counts - minimum_samples

    total_remaining_capacity = int(remaining_capacity.sum())

    if total_remaining_capacity == 0:
        return sample_counts

    # Distribute the remaining slots proportionally to each cell's remaining capacity
    proportional = (
        remaining_capacity
        / total_remaining_capacity
        * remaining_slots
    )

    extra_samples = np.floor(proportional).astype(int)

    # Prevent allocations from exceeding the cell's own sample count
    extra_samples = np.minimum(
        extra_samples,
        remaining_capacity,
    )

    sample_counts += extra_samples

    slots_left = actual_target - int(sample_counts.sum())

    # Use the largest-remainder method for the remaining slots after flooring
    remainders = proportional - np.floor(proportional)
    order = np.argsort(-remainders)

    for index in order:
        if slots_left <= 0:
            break

        if sample_counts[index] < counts[index]:
            sample_counts[index] += 1
            slots_left -= 1

    # If there are still remaining slots, assign them to cells that are not yet full
    while slots_left > 0:
        available_indexes = np.where(
            sample_counts < counts
        )[0]

        if len(available_indexes) == 0:
            break

        for index in available_indexes:
            if slots_left <= 0:
                break

            sample_counts[index] += 1
            slots_left -= 1

    return sample_counts


def sample_metadata(metadata, summary):
    sampled_parts = []

    metadata = metadata.copy()

    # Clean Isolate_Id values to avoid matching failures caused by whitespace
    metadata["_Clean_Isolate_Id"] = (
        metadata["Isolate_Id"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    available_ids = set(metadata["_Clean_Isolate_Id"])

    missing_required_ids = (
        REQUIRED_ISOLATE_IDS - available_ids
    )

    if missing_required_ids:
        raise ValueError(
            "The following required Isolate_Id values "
            "were not found in the simplified metadata: "
            + ", ".join(sorted(missing_required_ids))
        )

    for _, row in summary.iterrows():
        year = row["_Collection_Year"]
        subcontinent = row["_Subcontinent"]
        host_category = row["_Host_Category"]
        sample_count = int(row["sample_count"])

        grid = metadata[
            (metadata["_Collection_Year"] == year)
            & (metadata["_Subcontinent"] == subcontinent)
            & (metadata["_Host_Category"] == host_category)
        ].copy()

        if sample_count > len(grid):
            print(
                "WARNING: Requested sample count exceeds "
                "available metadata rows for grid: "
                f"{year}, {subcontinent}, {host_category}. "
                f"Requested={sample_count}, "
                f"Available={len(grid)}"
            )

            sample_count = len(grid)

        if sample_count == 0:
            continue

        # Identify the samples that must be preserved for this cell
        required_rows = grid[
            grid["_Clean_Isolate_Id"].isin(
                REQUIRED_ISOLATE_IDS
            )
        ].copy()

        # Remove the mandatory samples first, then randomly sample from the remaining ones
        remaining_grid = grid[
            ~grid["_Clean_Isolate_Id"].isin(
                REQUIRED_ISOLATE_IDS
            )
        ].copy()

        required_count = len(required_rows)

        if required_count > sample_count:
            raise ValueError(
                "The number of required samples in grid "
                f"{year}, {subcontinent}, {host_category} "
                f"is greater than sample_count."
            )

        random_sample_count = (
            sample_count - required_count
        )

        if random_sample_count >= len(remaining_grid):
            randomly_sampled_rows = remaining_grid.copy()
        elif random_sample_count > 0:
            randomly_sampled_rows = remaining_grid.sample(
                n=random_sample_count,
                random_state=RANDOM_SEED,
                replace=False,
            )
        else:
            randomly_sampled_rows = remaining_grid.iloc[
                0:0
            ].copy()

        sampled_grid = pd.concat(
            [
                required_rows,
                randomly_sampled_rows,
            ],
            ignore_index=False,
        )

        sampled_parts.append(sampled_grid)

    if not sampled_parts:
        raise ValueError("No metadata rows were sampled.")

    sampled_metadata = pd.concat(
        sampled_parts,
        ignore_index=True,
    )

    # Verify once more that the mandatory samples still exist at the end
    sampled_ids = set(
        sampled_metadata["_Clean_Isolate_Id"]
    )

    missing_after_sampling = (
        REQUIRED_ISOLATE_IDS - sampled_ids
    )

    if missing_after_sampling:
        raise ValueError(
            "Required Isolate_Id was not included after "
            "sampling: "
            + ", ".join(sorted(missing_after_sampling))
        )

    # Randomly shuffle the final sample order
    sampled_metadata = sampled_metadata.sample(
        frac=1,
        random_state=RANDOM_SEED,
    ).reset_index(drop=True)

    # Remove temporary internal ID columns
    sampled_metadata = sampled_metadata.drop(
        columns=["_Clean_Isolate_Id"],
        errors="ignore",
    )

    print(
        "Required Isolate_Id retained: "
        + ", ".join(sorted(REQUIRED_ISOLATE_IDS))
    )

    return sampled_metadata

def main():
    global REQUIRED_ISOLATE_IDS
    parser = argparse.ArgumentParser()
    parser.add_argument("--segment", required=True, type=str.upper, choices=("HA", "NA", "MP", "PB2", "PB1", "PA", "NP", "NS"))
    segment = parser.parse_args().segment
    # This isolate must remain in the sampled set for all gene segments.
    REQUIRED_ISOLATE_IDS = REQUIRED_ISOLATE_IDS.copy()
    # Assume sample.py is located at:
    # 01_sample_preparation/scripts/sample.py
    project_root = Path(__file__).resolve().parents[2]

    segment_dir = (
        project_root
        / "01_sample_preparation"
        / "output"
        / "genes"
        / segment
    )

    data_dir = project_root / "Data"
    data_dir.mkdir(parents=True, exist_ok=True)

    metadata_file = segment_dir / f"{segment}_metadata_simplified.xlsx"
    summary_file = segment_dir / f"{segment}_summary.xlsx"

    sampled_metadata_file = (
        segment_dir / f"{segment}_metadata_5000sample.xlsx"
    )

    new_summary_file = (
        segment_dir / f"{segment}_New_Summary_Table.xlsx"
    )

    reference_id_file = (
        data_dir / "reference ID.csv"
    )

    if not metadata_file.exists():
        raise FileNotFoundError(
            f"Metadata file does not exist: {metadata_file}"
        )

    if not summary_file.exists():
        raise FileNotFoundError(
            f"Summary file does not exist: {summary_file}"
        )

    metadata = pd.read_excel(metadata_file)
    summary = pd.read_excel(summary_file)

    metadata = prepare_metadata(metadata)
    summary = prepare_summary(summary)

    # Use the actual metadata count to correct the Count values in the summary
    metadata_counts = (
        metadata.groupby(
            [
                "_Collection_Year",
                "_Subcontinent",
                "_Host_Category",
            ],
            dropna=False,
        )
        .size()
        .reset_index(name="_Actual_Count")
    )

    summary = summary.merge(
        metadata_counts,
        on=[
            "_Collection_Year",
            "_Subcontinent",
            "_Host_Category",
        ],
        how="left",
    )

    summary["_Actual_Count"] = (
        summary["_Actual_Count"]
        .fillna(0)
        .astype(int)
    )

    count_mismatch = (
        summary["Count"] != summary["_Actual_Count"]
    )

    if count_mismatch.any():
        print(
            "WARNING: Some Count values in the summary file "
            "do not match the metadata."
        )
        print(
            "Actual metadata counts will be used for sampling."
        )

    summary["Count"] = summary["_Actual_Count"]

    summary["sample_count"] = allocate_sample_counts(
        summary=summary,
        target_size=TARGET_SAMPLE_SIZE,
    )

    sampled_metadata = sample_metadata(
        metadata=metadata,
        summary=summary,
    )

    # Remove temporary classification columns used internally by the script
    temporary_metadata_columns = [
        "_Collection_Year",
        "_Subcontinent",
        "_Host_Category",
    ]

    sampled_metadata = sampled_metadata.drop(
        columns=temporary_metadata_columns,
        errors="ignore",
    )

    temporary_summary_columns = [
        "_Collection_Year",
        "_Subcontinent",
        "_Host_Category",
        "_Actual_Count",
    ]

    output_summary = summary.drop(
        columns=temporary_summary_columns,
        errors="ignore",
    )

    # Save the sampled metadata
    sampled_metadata.to_excel(
        sampled_metadata_file,
        index=False,
    )

    # Save the new summary including sample_count
    output_summary.to_excel(
        new_summary_file,
        index=False,
    )

    # Generate reference ID.csv
    reference_ids = sampled_metadata[
        ["Isolate_Id"]
    ].copy()

    reference_ids.to_csv(
        reference_id_file,
        index=False,
        encoding="utf-8-sig",
    )

    print()
    print(f"Target sample size: {TARGET_SAMPLE_SIZE}")
    print(f"Actual sample size: {len(sampled_metadata)}")
    print(f"Sampled metadata: {sampled_metadata_file}")
    print(f"New summary table: {new_summary_file}")
    print(f"Reference ID file: {reference_id_file}")


if __name__ == "__main__":
    main()
