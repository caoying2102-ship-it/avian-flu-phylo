from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent
INPUT_DIR = ROOT / "01_sample_preparation" / "input"
SEGMENTS = ["HA", "NA", "MP", "PB2", "PB1", "PA", "NP", "NS"]

def clear_dir(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for child in directory.iterdir():
        if child.is_file() or child.is_symlink():
            child.unlink()
        elif child.is_dir():
            shutil.rmtree(child)

def save_uploaded_files(uploaded_files, target_dir: Path) -> None:
    clear_dir(target_dir)
    for uploaded in uploaded_files:
        target = target_dir / uploaded.name
        with open(target, "wb") as f:
            f.write(uploaded.getbuffer())

def detect_segment_from_names(file_names):
    tokens = {
        "HA": ["ha", "hemagglutinin"],
        "NA": ["na", "neuraminidase"],
        "MP": ["mp", "matrix"],
        "PB2": ["pb2", "polymerase basic 2"],
        "PB1": ["pb1", "polymerase basic 1"],
        "PA": ["pa", "polymerase acidic"],
        "NP": ["np", "nucleoprotein"],
        "NS": ["ns", "nonstructural"],
    }
    name_text = " ".join(str(n).lower() for n in file_names)
    for segment, keywords in tokens.items():
        if any(keyword in name_text for keyword in keywords):
            return segment
    return None

def detect_segment_from_fasta(fasta_path: Path):
    try:
        with open(fasta_path, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    header = line[1:].strip().lower()
                    for segment, keywords in {
                        "HA": ["ha", "hemagglutinin"],
                        "NA": ["na", "neuraminidase"],
                        "MP": ["mp", "matrix"],
                        "PB2": ["pb2", "polymerase basic 2"],
                        "PB1": ["pb1", "polymerase basic 1"],
                        "PA": ["pa", "polymerase acidic"],
                        "NP": ["np", "nucleoprotein"],
                        "NS": ["ns", "nonstructural"],
                    }.items():
                        if any(keyword in header for keyword in keywords):
                            return segment
                    break
    except Exception:
        pass
    return None

def run_pipeline(segment: str, clean: bool = True):
    cmd = [sys.executable, str(ROOT / "run_pipeline.py"), "--segment", segment]
    if clean:
        cmd.append("--clean")
    return subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True)

def find_output_dir(segment: str) -> Path:
    return ROOT / "04_iqtree_initial" / "output" / segment

def download_button_for_file(path: Path):
    if not path.exists():
        return
    with open(path, "rb") as f:
        data = f.read()
    st.download_button(
        label=f"Download {path.name}",
        data=data,
        file_name=path.name,
        mime="application/octet-stream",
    )

def read_metadata_preview(path: Path):
    if not path.exists():
        return None
    try:
        if path.suffix.lower() in {".xlsx", ".xls"}:
            return pd.read_excel(path)
        if path.suffix.lower() == ".csv":
            return pd.read_csv(path)
    except Exception:
        return None
    return None

st.set_page_config(page_title="Influenza IQ-TREE Pipeline", layout="wide")
st.title("Influenza IQ-TREE Pipeline")
st.caption("Upload FASTA and metadata files, auto-detect the segment, run the pipeline, and download the result files.")

segment_option = st.selectbox("Select segment", ["AUTO"] + SEGMENTS)

fasta_files = st.file_uploader(
    "Upload one or more FASTA files",
    type=["fa", "fas", "fasta", "fna"],
    accept_multiple_files=True,
)
metadata_files = st.file_uploader(
    "Upload one or more metadata files (.xlsx, .xls, .csv)",
    type=["xlsx", "xls", "csv"],
    accept_multiple_files=True,
)
clean_run = st.checkbox("Start from a clean state", value=True)

if st.button("Run pipeline", type="primary"):
    if not fasta_files:
        st.error("Please upload at least one FASTA file.")
        st.stop()
    if not metadata_files:
        st.error("Please upload at least one metadata file.")
        st.stop()

    detected = detect_segment_from_names([f.name for f in fasta_files + metadata_files])
    if detected is None and fasta_files:
        for fasta in fasta_files:
            tmp = ROOT / "tmp_detect.fasta"
            with open(tmp, "wb") as f:
                f.write(fasta.getbuffer())
            detected = detect_segment_from_fasta(tmp)
            tmp.unlink(missing_ok=True)
            if detected:
                break

    chosen_segment = segment_option
    if chosen_segment == "AUTO":
        if detected is None:
            st.error("Auto-detect failed. Please choose the segment manually.")
            st.stop()
        chosen_segment = detected

    chosen_segment = chosen_segment.upper()
    if chosen_segment not in SEGMENTS:
        st.error(f"Unsupported segment: {chosen_segment}")
        st.stop()

    st.info(f"Using segment: {chosen_segment}")

    save_uploaded_files(fasta_files, INPUT_DIR)
    for metadata in metadata_files:
        target = INPUT_DIR / metadata.name
        with open(target, "wb") as f:
            f.write(metadata.getbuffer())

    with st.spinner(f"Running pipeline for {chosen_segment}..."):
        result = run_pipeline(chosen_segment, clean=clean_run)

    if result.returncode != 0:
        st.error("Pipeline execution failed.")
        st.code(result.stdout)
        st.code(result.stderr)
        st.stop()

    st.success("Pipeline finished successfully.")

    output_dir = find_output_dir(chosen_segment)
    st.subheader("Output directory")
    st.write(output_dir)

    if output_dir.exists():
        files = sorted([p for p in output_dir.iterdir() if p.is_file()])
        cols = st.columns(3)
        for i, path in enumerate(files):
            with cols[i % 3]:
                download_button_for_file(path)

    tree_path = None
    for candidate in output_dir.glob("*.treefile"):
        tree_path = candidate
        break
    if tree_path is not None:
        st.subheader("Tree preview")
        tree_text = tree_path.read_text(encoding="utf-8", errors="ignore")
        st.code(tree_text[:4000])
    else:
        st.warning("No tree file found in the final output directory.")

    iqtree_path = None
    for candidate in output_dir.glob("*.iqtree"):
        iqtree_path = candidate
        break
    if iqtree_path is not None:
        st.subheader("IQ-TREE summary")
        text = iqtree_path.read_text(encoding="utf-8", errors="ignore")
        st.code(text[:4000])

    metadata_preview = None
    for candidate in sorted(output_dir.rglob("*")):
        if candidate.is_file() and candidate.suffix.lower() in {".xlsx", ".xls", ".csv"} and "metadata" in candidate.name.lower():
            df = read_metadata_preview(candidate)
            if df is not None:
                metadata_preview = df
                st.subheader(f"Metadata preview: {candidate.name}")
                st.dataframe(df.head(20), use_container_width=True)
                break

    if metadata_preview is None:
        st.warning("No metadata file was found in the final output directory.")