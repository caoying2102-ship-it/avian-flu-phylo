#!/usr/bin/env python3
"""Run IQ-TREE 2 using the aligned FASTA in the project's input folder."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


# Expected portable folder layout:
# 04_iqtree_initial/
# |-- input/HA/HA_reference_dedup_aligned_trimmed.fasta
# |-- scripts/run_iqtree.py
# `parent` moves from scripts/ back to the project folder.
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_FASTA = (
    PROJECT_DIR / "input" / "HA" / "HA_reference_dedup_aligned_trimmed.fasta"
)
DEFAULT_OUTDIR = PROJECT_DIR / "output"
FASTA_SUFFIXES = {".fa", ".fas", ".fasta", ".fna", ".faa", ".aln"}


def copy_metadata_to_output(fasta: Path, outdir: Path) -> Path | None:
    """Copy the paired metadata workbook from the input folder to the output folder."""
    if not fasta.exists():
        return None

    segment = fasta.stem.split("_")[0]
    metadata_candidates = [
        path
        for path in sorted(fasta.parent.iterdir())
        if path.is_file()
        and path.suffix.lower() in {".xlsx", ".xls"}
        and "metadata" in path.name.lower()
        and segment.lower() in path.name.lower()
    ]

    if not metadata_candidates:
        return None

    metadata_source = metadata_candidates[0]
    outdir.mkdir(parents=True, exist_ok=True)
    destination = outdir / metadata_source.name
    shutil.copy2(metadata_source, destination)
    return destination


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Use IQ-TREE 2 to infer a maximum-likelihood phylogenetic tree."
    )
    parser.add_argument(
        "fasta",
        nargs="?",
        type=Path,
        default=DEFAULT_FASTA,
        help=(
            "Aligned FASTA (default: "
            "../input/HA/HA_reference_dedup_aligned_trimmed.fasta)"
        ),
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help="Output directory (default: ../output)",
    )
    parser.add_argument(
        "-B",
        "--bootstrap",
        type=int,
        default=1000,
        help="Number of ultrafast bootstrap replicates (default: 1000)",
    )
    parser.add_argument(
        "-T",
        "--threads",
        default="AUTO",
        help="Number of threads, or AUTO (default: AUTO)",
    )
    parser.add_argument(
        "-m",
        "--model",
        default="MFP",
        help="Substitution model; MFP selects it automatically (default: MFP)",
    )
    parser.add_argument(
        "--redo",
        action="store_true",
        help="Overwrite results from an earlier run with the same prefix",
    )
    parser.add_argument(
        "--background",
        action="store_true",
        help=(
            "Run IQ-TREE in a detached background process; it continues after "
            "the terminal closes or the macOS user is switched"
        ),
    )
    return parser.parse_args()


def project_relative(path: Path) -> Path:
    """Resolve a user-supplied relative path from the portable project folder."""
    path = path.expanduser()
    return path.resolve() if path.is_absolute() else (PROJECT_DIR / path).resolve()


def find_default_fasta(expected: Path) -> Path | None:
    """Find the input FASTA even if Finder has hidden an extra extension."""
    if expected.is_file():
        return expected

    input_dir = PROJECT_DIR / "input"
    if not input_dir.is_dir():
        return None

    candidates = sorted(
        path
        for path in input_dir.rglob("*")
        if path.is_file()
        and (
            path.suffix.lower() in FASTA_SUFFIXES
            or ".fasta" in path.name.lower()
            or ".fa." in path.name.lower()
        )
    )
    if len(candidates) == 1:
        print(f"Expected filename was not found; automatically using: {candidates[0].name}")
        return candidates[0].resolve()
    if len(candidates) > 1:
        print("Error: multiple possible FASTA files were found:", file=sys.stderr)
        for candidate in candidates:
            print(f"  - {candidate.name}", file=sys.stderr)
        print("Pass the desired filename to the script explicitly.", file=sys.stderr)
    return None


def main() -> int:
    args = parse_args()
    fasta = project_relative(args.fasta)
    if not fasta.is_file():
        # Automatic discovery is only used for the built-in default input.
        if args.fasta == DEFAULT_FASTA:
            discovered = find_default_fasta(fasta)
            if discovered is not None:
                fasta = discovered
            else:
                print(f"Error: FASTA file not found: {fasta}", file=sys.stderr)
                input_dir = PROJECT_DIR / "input"
                if input_dir.is_dir():
                    names = sorted(path.name for path in input_dir.iterdir())
                    print("Files currently present in input/:", file=sys.stderr)
                    for name in names:
                        print(f"  - {name}", file=sys.stderr)
                else:
                    print(f"The input folder does not exist: {input_dir}", file=sys.stderr)
                return 1
        else:
            print(f"Error: FASTA file not found: {fasta}", file=sys.stderr)
            return 1

    iqtree = shutil.which("iqtree2") or shutil.which("iqtree")
    if iqtree is None:
        print(
            "Error: IQ-TREE was not found. Install it first, for example:\n"
            "  conda install -c bioconda iqtree",
            file=sys.stderr,
        )
        return 1

    if args.bootstrap < 1000:
        print("Warning: at least 1000 ultrafast bootstrap replicates are recommended.")

    outdir = project_relative(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outdir / fasta.stem

    command = [
        iqtree,
        "-s",
        str(fasta),
        "-m",
        args.model,
        "-B",
        str(args.bootstrap),
        "-T",
        str(args.threads),
        "--prefix",
        str(prefix),
    ]
    if args.redo:
        command.append("-redo")

    print("Input: ", fasta)
    print("Output:", outdir)
    print("Running:", " ".join(command))

    if args.background:
        background_log = Path(f"{prefix}.background.log")
        pid_file = Path(f"{prefix}.pid")
        try:
            with background_log.open("a", encoding="utf-8") as log_handle:
                process = subprocess.Popen(
                    command,
                    cwd=PROJECT_DIR,
                    stdin=subprocess.DEVNULL,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                    close_fds=True,
                )
            pid_file.write_text(f"{process.pid}\n", encoding="utf-8")
        except OSError as exc:
            print(f"Could not start IQ-TREE in the background: {exc}", file=sys.stderr)
            return 1

        print("\nIQ-TREE is now running in the background.")
        print(f"Process ID: {process.pid}")
        print(f"PID file:   {pid_file}")
        print(f"Progress:   {background_log}")
        print(f"Best tree (when finished): {prefix}.treefile")
        return 0

    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"IQ-TREE failed with exit code {exc.returncode}.", file=sys.stderr)
        return exc.returncode

    copied_metadata = copy_metadata_to_output(fasta, outdir)
    if copied_metadata is not None:
        print(f"Metadata:    {copied_metadata}")

    print(f"\nFinished. Best tree: {prefix}.treefile")
    print(f"Full report: {prefix}.iqtree")
    print(f"Run log:     {prefix}.log")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
