#!/usr/bin/env python3
"""Run one influenza segment through all four pipeline stages."""

import argparse
import importlib.util
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import yaml


PROJECT_ROOT = Path(__file__).resolve().parent
CONFIG_FILE = PROJECT_ROOT / "config.yaml"


def load_config():
    with CONFIG_FILE.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    segments = tuple(config.get("segments", ()))
    if not segments:
        raise ValueError("config.yaml must define a non-empty segments list")
    return config, segments


def run_script(stage_name, relative_path, *arguments):
    script_path = PROJECT_ROOT / relative_path
    if not script_path.is_file():
        raise FileNotFoundError(f"Script not found: {script_path}")
    command = [sys.executable, str(script_path), *map(str, arguments)]
    print(f"\n{'=' * 60}\nRunning step: {stage_name}\nCommand: {' '.join(command)}\n{'=' * 60}")
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)


def check_requirements(config):
    missing = []
    for module in config.get("python_modules", []):
        if importlib.util.find_spec(module) is None:
            missing.append(f"Python module: {module}")

    resolved = {}
    for name, candidates in config.get("executables", {}).items():
        if isinstance(candidates, str):
            candidates = [candidates]
        executable = next((shutil.which(item) for item in candidates if shutil.which(item)), None)
        if executable is None:
            missing.append(f"Software: {name} (command: {' / '.join(candidates)})")
        else:
            resolved[name] = executable

    if missing:
        raise RuntimeError("Missing runtime dependencies:\n  - " + "\n  - ".join(missing))
    print("Dependency check passed: " + ", ".join(sorted(resolved)))


def copy_data_to_input(overwrite: bool = False) -> bool:
    """Copy supported raw data from Data/ into 01_sample_preparation/input/."""
    data_dir = PROJECT_ROOT / "Data"
    input_dir = PROJECT_ROOT / "01_sample_preparation" / "input"
    input_dir.mkdir(parents=True, exist_ok=True)

    if not data_dir.exists():
        return False

    supported_exts = {".fa", ".fas", ".fasta", ".fna", ".xls", ".xlsx"}
    files = sorted(
        path
        for path in data_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in supported_exts
    )
    if not files:
        return False

    copied = 0
    for source in files:
        destination = input_dir / source.name
        if destination.exists() and not overwrite:
            print(f"Skipping existing file: {destination}")
            continue
        shutil.copy2(source, destination)
        print(f"Copied: {source.name} -> {destination}")
        copied += 1

    if copied == 0:
        print("No new files were copied; input directory already contains the required files.")
    else:
        print(f"Copied {copied} file(s) from Data/ to: {input_dir}")
    return True


def check_data_inputs():
    data_dir = PROJECT_ROOT / "Data"
    if not data_dir.is_dir():
        raise FileNotFoundError(f"Data directory does not exist: {data_dir}")
    has_metadata = any(data_dir.rglob("*.xlsx")) or any(data_dir.rglob("*.xls"))
    has_fasta = (
        any(data_dir.rglob("*.fasta"))
        or any(data_dir.rglob("*.fa"))
        or any(data_dir.rglob("*.fas"))
        or any(data_dir.rglob("*.fna"))
    )
    if not has_metadata or not has_fasta:
        raise FileNotFoundError("Data requires at least one metadata file and one FASTA file")


def check_initial_inputs():
    input_dir = PROJECT_ROOT / "01_sample_preparation" / "input"
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")
    has_excel = any(input_dir.glob("*.xlsx")) or any(input_dir.glob("*.xls"))
    has_fasta = any(input_dir.glob("*.fasta")) or any(input_dir.glob("*.fa")) or any(input_dir.glob("*.fas")) or any(input_dir.glob("*.fna"))
    if not has_excel or not has_fasta:
        raise FileNotFoundError("01_sample_preparation/input requires both Excel and FASTA files")


def archive_previous_outputs(segment: str) -> None:
    """Archive any existing output directories for the selected segment before a fresh run."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archive_root = PROJECT_ROOT / "logs" / "archive" / timestamp
    paths = [
        PROJECT_ROOT / "01_sample_preparation" / "output",
        PROJECT_ROOT / "02_msa" / "output" / segment,
        PROJECT_ROOT / "03_trimal" / "output" / segment,
        PROJECT_ROOT / "04_iqtree_initial" / "output" / segment,
    ]

    moved_any = False
    for path in paths:
        if path.exists():
            target = archive_root / path.relative_to(PROJECT_ROOT)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(path), str(target))
            print(f"Archived existing output: {path} -> {target}")
            moved_any = True

    if moved_any:
        print(f"Archived previous outputs to: {archive_root}")
    else:
        print("No previous segment outputs were found to archive.")


def run_pipeline(segment, config):
    check_requirements(config)
    check_data_inputs()
    copy_data_to_input()
    check_initial_inputs()

    scripts_01 = [
        ("prepare_samples.py", "--genes", segment),
        ("generate_simplified_metadata.py", "--segment", segment),
        ("rename_fasta_headers.py", "--segment", segment),
        ("summary_count.py", "--segment", segment),
        ("sample.py", "--segment", segment),
        ("extract.py", "--segment", segment),
        ("SeqKit.py", "--segment", segment),
    ]
    for script, *args in scripts_01:
        run_script(f"01 sample preparation: {script}", Path("01_sample_preparation/scripts") / script, *args)

    # Metadata matching depends on the aligned FASTA, so it must follow MAFFT.
    for label, script in [
        ("02 prepare MAFFT inputs", "prepare_msa_inputs.py"),
        ("02 run MAFFT", "msa.py"),
        ("02 match aligned metadata", "matchmetadata.py"),
    ]:
        run_script(label, Path("02_msa/scripts") / script, "--segment", segment)

    for label, script in [
        ("03 prepare trimAl inputs", "prepare_trimal_inputs.py"),
        ("03 run trimAl", "trimal.py"),
        ("04 prepare IQ-TREE inputs", "prepare_iqtree_inputs.py"),
    ]:
        stage_dir = "03_trimal/scripts" if label.startswith("03") else "04_iqtree_initial/scripts"
        run_script(label, Path(stage_dir) / script, "--segment", segment)

    fasta = f"input/{segment}/{segment}_reference_dedup_aligned_trimmed.fasta"
    outdir = f"output/{segment}"
    run_script(f"04 run IQ-TREE: {segment}", Path("04_iqtree_initial/scripts/run_iqtree.py"), fasta, "--outdir", outdir)
    print(f"\nPipeline completed successfully: {segment}")


def main():
    config, valid_segments = load_config()
    parser = argparse.ArgumentParser(description="Run one segment through stages 01–04.")
    parser.add_argument("segment_positional", nargs="?", type=str.upper, choices=valid_segments, metavar="SEGMENT")
    parser.add_argument("--segment", dest="segment_option", type=str.upper, choices=valid_segments)
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Archive and remove existing outputs before starting a fresh run.",
    )
    args = parser.parse_args()
    if args.segment_positional and args.segment_option:
        parser.error("Use either SEGMENT or --segment SEGMENT, not both")
    segment = args.segment_option or args.segment_positional
    if not segment:
        parser.error("A segment must be specified, for example: --segment HA")

    if args.clean:
        archive_previous_outputs(segment)

    try:
        run_pipeline(segment, config)
    except subprocess.CalledProcessError as error:
        print(f"\nPipeline failed with exit code: {error.returncode}", file=sys.stderr)
        return error.returncode
    except Exception as error:
        print(f"\nError: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
