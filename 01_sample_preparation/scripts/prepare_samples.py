#!/usr/bin/env python3
"""Merge, split, deduplicate, and validate metadata and FASTA inputs."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
STAGE_DIR = PROJECT_ROOT / "01_sample_preparation"
DEFAULT_GENES = ("PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS")
EXCEL_SUFFIXES = {".xls", ".xlsx"}
FASTA_SUFFIXES = {".fasta", ".fa", ".fas", ".fna"}
SIMPLIFIED_COLUMNS = (
    "strain",
    "Isolate_Id",
    "Continent",
    "Country",
    "Region",
    "Location",
    "Isolate_Name",
    "Collection_Date",
    "Clade",
    "Source",
    "Segment",
    "Subcontinent",
)


def get_subcontinent_map() -> dict[str, str]:
    return {
        "China": "EastAsia", "Japan": "EastAsia", "Korea,Republicof": "EastAsia",
        "HongKong(SAR)": "EastAsia", "Taiwan": "EastAsia", "Mongolia": "EastAsia",
        "Korea,DemocraticPeople'sRepublicof": "EastAsia", "Tokoname": "EastAsia",
        "Vietnam": "SoutheastAsia", "Thailand": "SoutheastAsia",
        "Indonesia": "SoutheastAsia", "Singapore": "SoutheastAsia",
        "Cambodia": "SoutheastAsia", "Malaysia": "SoutheastAsia",
        "Laos": "SoutheastAsia", "Myanmar": "SoutheastAsia",
        "Philippines": "SoutheastAsia",
        "Lao,People'sDemocraticRepublic": "SoutheastAsia",
        "PreyVeng": "SoutheastAsia", "India": "SouthAsia",
        "Pakistan": "SouthAsia", "Bangladesh": "SouthAsia",
        "SriLanka": "SouthAsia", "Bhutan": "SouthAsia", "Nepal": "SouthAsia",
        "Afghanistan": "SouthAsia", "Israel": "WestAsia", "Lebanon": "WestAsia",
        "Iraq": "WestAsia", "Turkey": "WestAsia", "Azerbaijan": "WestAsia",
        "PalestinianTerritory": "WestAsia", "Georgia": "WestAsia",
        "Iran,IslamicRepublicof": "WestAsia", "SaudiArabia": "MiddleEastAsia",
        "UnitedArabEmirates": "MiddleEastAsia", "Kuwait": "MiddleEastAsia",
        "Kazakhstan": "CentralAsia",
    }


def remove_all_whitespace(value: object) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", "", str(value))


def classify_region(location: object) -> str:
    loc = remove_all_whitespace(location).casefold()
    usa_terms = (
        "usa", "unitedstates", "northamerica", "california", "newyork",
        "ohio", "atlanta", "swiftwater", "texas", "michigan", "colorado",
    )
    if any(term in loc for term in usa_terms):
        return "North America"
    if "sudan" in loc:
        return "Africa"
    asia_map = {
        "Asia_Central": ("kazakhstan", "uzbekistan"),
        "Asia_EastAsia": ("china", "japan", "korea", "taiwan", "hongkong"),
        "Asia_MiddleEast": ("saudiarabia", "iran", "iraq", "israel", "turkey"),
        "Asia_SEAsia": ("vietnam", "thailand", "indonesia", "cambodia", "laos", "malaysia"),
        "Asia_SouthAsia": ("india", "bangladesh", "pakistan"),
    }
    for region, countries in asia_map.items():
        if any(country in loc for country in countries):
            return region
    for continent in ("Europe", "South America", "Oceania", "Africa"):
        if remove_all_whitespace(continent).casefold() in loc:
            return continent
    return "Other/Unknown"


def normalize_collection_date(value: object) -> str:
    if pd.isna(value) or not str(value).strip():
        return "Unknown_Date"
    parsed = pd.to_datetime(value, errors="coerce")
    if pd.isna(parsed):
        return "Unknown_Date"
    return parsed.strftime("%Y-%m-%d")


def make_simplified_metadata(metadata: pd.DataFrame, segment: str) -> pd.DataFrame:
    core_columns = (
        "Isolate_Id", "Continent", "Country", "Region", "Location",
        "Isolate_Name", "Collection_Date", "Clade",
    )
    simplified = pd.DataFrame(index=metadata.index)
    for column in core_columns:
        if column in metadata.columns:
            simplified[column] = metadata[column]
        else:
            simplified[column] = ""

    # Remove every whitespace character from metadata values before standardization.
    for column in core_columns:
        simplified[column] = simplified[column].map(remove_all_whitespace)

    # Location has priority for geographic fields when it contains slash-separated values.
    location_parts = simplified["Location"].str.split("/", n=2, expand=True)
    for index, column in enumerate(("Continent", "Country")):
        if index in location_parts.columns:
            extracted = location_parts[index].fillna("")
            simplified[column] = extracted.where(extracted.ne(""), simplified[column])

    simplified["Region"] = simplified["Location"].map(classify_region)
    simplified["Collection_Date"] = metadata.get(
        "Collection_Date", pd.Series("", index=metadata.index)
    ).map(normalize_collection_date)
    simplified["Source"] = "GISAID"
    simplified["Segment"] = segment
    subcontinent_map = get_subcontinent_map()
    simplified["Subcontinent"] = simplified["Country"].map(subcontinent_map).fillna(
        "Other/Unknown"
    )
    simplified["strain"] = (
        simplified["Continent"].replace("", "Unknown_Continent") + "|"
        + simplified["Isolate_Name"].replace("", "Unknown_Isolate") + "|"
        + simplified["Collection_Date"] + "|"
        + simplified["Isolate_Id"].replace("", "Unknown_Id")
    )
    return simplified.loc[:, SIMPLIFIED_COLUMNS]


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str
    source_file: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge any number of Excel and FASTA files, split sequences by gene, "
            "remove duplicate sequences, and match sequence IDs to metadata."
        )
    )
    parser.add_argument("--input-dir", type=Path, default=STAGE_DIR / "input")
    parser.add_argument("--output-dir", type=Path, default=STAGE_DIR / "output")
    parser.add_argument("--genes", nargs="+", default=list(DEFAULT_GENES))
    parser.add_argument("--id-column", default="Isolate_Id")
    parser.add_argument(
        "--header-id-index",
        type=int,
        default=3,
        help="Zero-based index of the isolate ID in pipe-delimited FASTA headers.",
    )
    parser.add_argument(
        "--all-sheets",
        action="store_true",
        help="Read every Excel worksheet instead of only the first worksheet.",
    )
    parser.add_argument(
        "--no-recursive",
        action="store_true",
        help="Do not search for inputs in subdirectories.",
    )
    return parser.parse_args()


def discover_files(input_dir: Path, suffixes: set[str], recursive: bool) -> list[Path]:
    paths = input_dir.rglob("*") if recursive else input_dir.glob("*")
    return sorted(p for p in paths if p.is_file() and p.suffix.lower() in suffixes)


def read_excel_files(files: Sequence[Path], all_sheets: bool) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in files:
        try:
            if all_sheets:
                sheet_items = pd.read_excel(path, sheet_name=None).items()
            else:
                sheet_items = [("first_sheet", pd.read_excel(path, sheet_name=0))]
        except Exception as exc:
            raise RuntimeError(f"Could not read Excel file {path}: {exc}") from exc

        for sheet_name, frame in sheet_items:
            frame = frame.copy()
            frame.columns = [str(column).strip() for column in frame.columns]
            frame["Source_File"] = path.name
            frame["Source_Sheet"] = str(sheet_name)
            frames.append(frame)
    return pd.concat(frames, ignore_index=True, sort=False)


def make_record(header: str, lines: Sequence[str], path: Path) -> FastaRecord:
    sequence = "".join(lines).upper()
    if not sequence:
        raise ValueError(f"Empty FASTA record in {path}: >{header}")
    return FastaRecord(header, sequence, path.name)


def read_fasta(path: Path) -> Iterator[FastaRecord]:
    header: str | None = None
    sequence_lines: list[str] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield make_record(header, sequence_lines, path)
                header = line[1:].strip()
                sequence_lines = []
            else:
                if header is None:
                    raise ValueError(
                        f"Sequence before first FASTA header in {path}, line {line_number}."
                    )
                sequence_lines.append(re.sub(r"\s+", "", line))
    if header is not None:
        yield make_record(header, sequence_lines, path)


def write_fasta(records: Sequence[FastaRecord], path: Path, width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(f">{record.header}\n")
            for start in range(0, len(record.sequence), width):
                handle.write(record.sequence[start : start + width] + "\n")


def detect_gene(header: str, genes: Sequence[str]) -> str | None:
    # Exact token matching prevents short names such as NA from matching words such as China.
    tokens = {token.upper() for token in re.findall(r"[A-Za-z0-9]+", header)}
    matches = [gene for gene in genes if gene.upper() in tokens]
    if not matches:
        return None
    return sorted(matches, key=lambda gene: (-len(gene), genes.index(gene)))[0]


def deduplicate_by_sequence(records: Sequence[FastaRecord]) -> list[FastaRecord]:
    seen: set[str] = set()
    unique: list[FastaRecord] = []
    for record in records:
        if record.sequence not in seen:
            seen.add(record.sequence)
            unique.append(record)
    return unique


def extract_isolate_id(header: str, index: int) -> str | None:
    fields = [field.strip() for field in header.split("|")]
    if -len(fields) <= index < len(fields) and fields[index]:
        return fields[index]
    return None


def save_metadata(frame: pd.DataFrame, output_dir: Path) -> None:
    frame.to_csv(output_dir / "combined_metadata.csv", index=False, encoding="utf-8-sig")
    frame.to_excel(output_dir / "combined_metadata.xlsx", index=False)


def run(args: argparse.Namespace) -> int:
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    recursive = not args.no_recursive
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    excel_files = discover_files(input_dir, EXCEL_SUFFIXES, recursive)
    fasta_files = discover_files(input_dir, FASTA_SUFFIXES, recursive)
    if not excel_files:
        raise FileNotFoundError(f"No .xls or .xlsx files found in {input_dir}")
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in {input_dir}")
    print(f"Found {len(excel_files)} Excel files and {len(fasta_files)} FASTA files.")

    metadata = read_excel_files(excel_files, args.all_sheets)
    if args.id_column not in metadata.columns:
        columns = ", ".join(map(str, metadata.columns))
        raise KeyError(f"Required column {args.id_column!r} is absent. Columns: {columns}")
    metadata[args.id_column] = metadata[args.id_column].astype("string").str.strip()
    save_metadata(metadata, output_dir)

    all_records = [record for path in fasta_files for record in read_fasta(path)]
    write_fasta(all_records, output_dir / "combined_sequences.fasta")
    by_gene: dict[str, list[FastaRecord]] = {gene: [] for gene in args.genes}
    unassigned: list[FastaRecord] = []
    for record in all_records:
        gene = detect_gene(record.header, args.genes)
        (by_gene[gene] if gene else unassigned).append(record)

    qc_rows: list[dict[str, object]] = []
    for gene, records in by_gene.items():
        gene_dir = output_dir / "genes" / gene
        unique = deduplicate_by_sequence(records)
        write_fasta(records, gene_dir / f"{gene}.fasta")
        write_fasta(unique, gene_dir / f"{gene}_dedup.fasta")
        isolate_ids = {
            isolate_id
            for record in unique
            if (isolate_id := extract_isolate_id(record.header, args.header_id_index))
        }
        matched = metadata[metadata[args.id_column].isin(isolate_ids)].copy()
        matched.to_csv(gene_dir / f"{gene}_metadata.csv", index=False, encoding="utf-8-sig")
        matched.to_excel(gene_dir / f"{gene}_metadata.xlsx", index=False)
        simplified = make_simplified_metadata(matched, gene)
        simplified.to_excel(
            gene_dir / f"{gene}_metadata_simplified.xlsx", index=False
        )
        qc_rows.append(
            {
                "Gene": gene,
                "Original_Sequences": len(records),
                "After_Deduplication": len(unique),
                "Removed_Duplicates": len(records) - len(unique),
                "FASTA_Isolate_IDs": len(isolate_ids),
                "Metadata_Matches": len(matched),
            }
        )

    if unassigned:
        write_fasta(unassigned, output_dir / "unassigned_sequences.fasta")
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(output_dir / "qc_summary.csv", index=False, encoding="utf-8-sig")
    qc.to_excel(output_dir / "qc_summary.xlsx", index=False)

    fasta_counts = Counter(record.source_file for record in all_records)
    with (output_dir / "input_manifest.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["File_Type", "File_Name", "Record_Count"])
        writer.writerows(["Excel", path.name, ""] for path in excel_files)
        writer.writerows(["FASTA", path.name, fasta_counts[path.name]] for path in fasta_files)

    print(qc.to_string(index=False))
    print(f"Unassigned sequences: {len(unassigned)}")
    print(f"Sample preparation completed: {output_dir}")
    return 0


def main() -> int:
    try:
        return run(parse_args())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
