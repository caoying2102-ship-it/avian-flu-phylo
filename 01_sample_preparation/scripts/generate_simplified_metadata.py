#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd


SEGMENTS = ("MP", "HA", "NA", "PB2", "PB1", "PA", "NP", "NS")

PROJECT_ROOT = Path(__file__).resolve().parents[2]

OUTPUT_COLUMNS = [
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
    "Host",
]


def remove_whitespace(value):
    if pd.isna(value):
        return ""

    return re.sub(r"\s+", "", str(value))


def normalize_location(value):
    if pd.isna(value):
        return ""

    return re.sub(r"\s+", "", str(value)).casefold()


def classify_subcontinent(location):
    loc = normalize_location(location)

    if not loc:
        return "Other/Unknown"

    region_map = {
        "North America": [
            "northamerica",
            "usa",
            "unitedstates",
            "canada",
            "mexico",
            "california",
            "newyork",
            "ohio",
            "atlanta",
            "swiftwater",
            "texas",
            "michigan",
            "colorado",
            "alaska",
        ],
        "Central America": [
            "centralamerica",
            "guatemala",
            "belize",
            "honduras",
            "elsalvador",
            "nicaragua",
            "costarica",
            "panama",
        ],
        "Caribbean": [
            "caribbean",
            "cuba",
            "jamaica",
            "haiti",
            "dominicanrepublic",
            "puertorico",
            "bahamas",
            "barbados",
            "trinidadandtobago",
        ],
        "South America": [
            "southamerica",
            "argentina",
            "bolivia",
            "brazil",
            "chile",
            "colombia",
            "ecuador",
            "guyana",
            "paraguay",
            "peru",
            "suriname",
            "uruguay",
            "venezuela",
        ],
        "Asia_Central": [
            "centralasia",
            "kazakhstan",
            "uzbekistan",
            "kyrgyzstan",
            "tajikistan",
            "turkmenistan",
        ],
        "Asia_EastAsia": [
            "eastasia",
            "china",
            "japan",
            "korea",
            "taiwan",
            "hongkong",
            "mongolia",
            "tokoname",
        ],
        "Asia_MiddleEast": [
            "middleeast",
            "westasia",
            "saudiarabia",
            "iran",
            "iraq",
            "israel",
            "turkey",
            "türkiye",
            "lebanon",
            "jordan",
            "syria",
            "yemen",
            "oman",
            "qatar",
            "kuwait",
            "unitedarabemirates",
            "bahrain",
            "palestine",
            "palestinianterritory",
            "azerbaijan",
            "georgia",
            "armenia",
        ],
        "Asia_SEAsia": [
            "southeastasia",
            "vietnam",
            "thailand",
            "indonesia",
            "cambodia",
            "laos",
            "laopeople'sdemocraticrepublic",
            "malaysia",
            "singapore",
            "myanmar",
            "philippines",
            "brunei",
            "timorleste",
            "preyveng",
        ],
        "Asia_SouthAsia": [
            "southasia",
            "india",
            "bangladesh",
            "pakistan",
            "srilanka",
            "nepal",
            "bhutan",
            "afghanistan",
            "maldives",
        ],
        "Europe": [
            "europe",
            "unitedkingdom",
            "england",
            "scotland",
            "wales",
            "ireland",
            "france",
            "germany",
            "italy",
            "spain",
            "portugal",
            "netherlands",
            "belgium",
            "switzerland",
            "austria",
            "poland",
            "denmark",
            "sweden",
            "norway",
            "finland",
            "iceland",
            "greece",
            "romania",
            "bulgaria",
            "hungary",
            "czechia",
            "czechrepublic",
            "slovakia",
            "slovenia",
            "croatia",
            "serbia",
            "ukraine",
            "russia",
        ],
        "Africa": [
            "africa",
            "sudan",
            "southsudan",
            "egypt",
            "libya",
            "algeria",
            "morocco",
            "tunisia",
            "ethiopia",
            "kenya",
            "uganda",
            "tanzania",
            "nigeria",
            "ghana",
            "senegal",
            "cameroon",
            "congo",
            "southafrica",
            "zambia",
            "zimbabwe",
            "mozambique",
            "botswana",
            "namibia",
            "madagascar",
        ],
        "Oceania": [
            "oceania",
            "australia",
            "newzealand",
            "papuanewguinea",
            "fiji",
            "samoa",
            "tonga",
            "vanuatu",
            "solomonislands",
            "micronesia",
        ],
    }

    for region, terms in region_map.items():
        if any(term in loc for term in terms):
            return region

    return "Other/Unknown"


def normalize_date(value):
    if pd.isna(value) or str(value).strip() == "":
        return "Unknown_Date"

    if isinstance(value, (int, float)) and not isinstance(value, bool):
        parsed = pd.to_datetime(
            value,
            unit="D",
            origin="1899-12-30",
            errors="coerce",
        )
    else:
        parsed = pd.to_datetime(
            value,
            errors="coerce",
        )

    if pd.isna(parsed):
        return "Unknown_Date"

    return parsed.strftime("%Y-%m-%d")


def get_column(df, column):
    if column in df.columns:
        return df[column].copy()

    return pd.Series(
        "",
        index=df.index,
        dtype="object",
    )


def extract_geography(simplified):
    location_parts = simplified["Location"].str.split(
        "/",
        n=2,
        expand=True,
    )

    geographic_columns = {
        0: "Continent",
        1: "Country",
        2: "Region",
    }

    for position, column in geographic_columns.items():
        if position not in location_parts.columns:
            continue

        extracted = location_parts[position].fillna("")

        simplified[column] = extracted.where(
            extracted != "",
            simplified[column],
        )

    return simplified


def generate_simplified_metadata(input_file, segment):
    df = pd.read_excel(input_file)

    df.columns = [
        remove_whitespace(column)
        for column in df.columns
    ]

    simplified = pd.DataFrame(index=df.index)

    core_columns = [
        "Isolate_Id",
        "Continent",
        "Country",
        "Region",
        "Location",
        "Isolate_Name",
        "Collection_Date",
        "Clade",
        "Host",
    ]

    for column in core_columns:
        simplified[column] = get_column(
            df,
            column,
        ).map(remove_whitespace)

    simplified = extract_geography(simplified)

    simplified["Collection_Date"] = get_column(
        df,
        "Collection_Date",
    ).map(normalize_date)

    simplified["Source"] = "GISAID"
    simplified["Segment"] = segment

    simplified["Subcontinent"] = simplified["Location"].map(
        classify_subcontinent
    )

    simplified["strain"] = (
        simplified["Continent"].replace(
            "",
            "Unknown_Continent",
        )
        + "|"
        + simplified["Isolate_Name"].replace(
            "",
            "Unknown_Isolate",
        )
        + "|"
        + simplified["Collection_Date"]
        + "|"
        + simplified["Isolate_Id"].replace(
            "",
            "Unknown_Id",
        )
    )

    simplified = simplified[OUTPUT_COLUMNS]

    output_file = (
        input_file.parent
        / f"{segment}_metadata_simplified.xlsx"
    )

    simplified.to_excel(
        output_file,
        index=False,
    )

    print(
        f"{segment}: {len(simplified)} rows -> {output_file}"
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

    if not genes_dir.exists():
        raise FileNotFoundError(
            f"Gene output directory does not exist: {genes_dir}\n"
            "Run 01_sample_preparation/scripts/run.sh first."
        )

    processed = 0

    input_file = genes_dir / segment / f"{segment}_metadata.xlsx"

    if not input_file.exists():
        raise FileNotFoundError(f"file not found: {input_file}")

    generate_simplified_metadata(input_file=input_file, segment=segment)
    processed += 1

    print(
        f"Completed simplified metadata for "
        f"{processed}/1 segment."
    )


if __name__ == "__main__":
    main()
