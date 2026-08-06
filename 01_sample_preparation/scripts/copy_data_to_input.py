#!/usr/bin/env python3
"""Copy raw data files from Data/ into 01_sample_preparation/input/."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = PROJECT_ROOT / "Data"
INPUT_DIR = PROJECT_ROOT / "01_sample_preparation" / "input"
SUPPORTED_EXTENSIONS = {".fa", ".fas", ".fasta", ".fna", ".xls", ".xlsx"}


def copy_data(overwrite: bool = False) -> int:
    if not DATA_DIR.exists():
        raise FileNotFoundError(f"Data directory not found: {DATA_DIR}")
    INPUT_DIR.mkdir(parents=True, exist_ok=True)

    files = sorted(
        path
        for path in DATA_DIR.rglob("*")
        if path.is_file() and path.suffix.lower() in SUPPORTED_EXTENSIONS
    )
    if not files:
        raise FileNotFoundError(
            f"No supported files found in Data/. Place metadata and FASTA files there first."
        )

    copied = 0
    for source in files:
        destination = INPUT_DIR / source.name
        if destination.exists() and not overwrite:
            print(f"Skipping existing file: {destination}")
            continue
        shutil.copy2(source, destination)
        print(f"Copied: {source.name}")
        copied += 1

    if copied == 0:
        print("No new files were copied. Use --overwrite to replace existing files.")
    else:
        print(f"Copied {copied} file(s) to: {INPUT_DIR}")
    return copied


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Copy raw FASTA and metadata files from Data/ to 01_sample_preparation/input/."
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing files in 01_sample_preparation/input/.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        copy_data(overwrite=args.overwrite)
    except Exception as exc:
        print(f"Error: {exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
