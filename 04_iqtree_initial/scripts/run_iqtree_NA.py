#!/usr/bin/env python3
"""Run IQ-TREE for the NA segment."""
import runpy
import sys
from pathlib import Path

GENE = "NA"
SCRIPT_DIR = Path(__file__).resolve().parent
sys.argv = [
    str(SCRIPT_DIR / "run_iqtree.py"),
    f"input/{GENE}/{GENE}_reference_dedup_aligned_trimmed.fasta",
    "-o", f"output/{GENE}",
    *sys.argv[1:],
]
runpy.run_path(SCRIPT_DIR / "run_iqtree.py", run_name="__main__")
