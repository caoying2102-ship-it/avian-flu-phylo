# 04 IQ-TREE Initial Analysis

This stage takes the trimmed per-gene alignments from stage 03 and runs IQ-TREE 2 to infer a maximum-likelihood phylogenetic tree.

The actual workflow implemented by this project is:

1. Copy the trimmed FASTA and paired metadata from the TrimAl output into the IQ-TREE input directory.
2. Remove duplicate FASTA headers if present.
3. Validate that the source files exist and are non-empty.
4. Run IQ-TREE 2 on the selected alignment.
5. Save the output tree, model report, and log files into the output directory.
6. Copy the matching metadata workbook into the output folder when available.

---

## Scripts in this folder

- `prepare_iqtree_inputs.py`
  - Copies aligned and trimmed FASTA files from stage 03 into the IQ-TREE input folder.
  - Deduplicates FASTA headers so each header appears only once.

- `run_iqtree.py`
  - Main tree-inference script.
  - Finds IQ-TREE, runs the selected FASTA, and writes the resulting phylogenetic output files.

---

## Inputs

This stage expects the trimmed alignment files from stage 03:

```text
03_trimal/output/{SEGMENT}/
    {SEGMENT}_reference_dedup_aligned_trimmed.fasta
    {SEGMENT}_metadata_reference_dedup_aligned_trimmed.xlsx
```

These are copied into:

```text
04_iqtree_initial/input/{SEGMENT}/
```

Supported genes:

```text
MP, HA, NA, PB2, PB1, PA, NP, NS
```

The default behavior is to use HA as the main target input.

---

## Step 1: Prepare IQ-TREE inputs

```bash
python3 04_iqtree_initial/scripts/prepare_iqtree_inputs.py --segment HA
```

This script:

- builds the source-to-destination file mapping from TrimAl outputs to IQ-TREE inputs
- checks that all files exist and are non-empty
- copies the metadata workbook and FASTA file
- removes repeated FASTA headers while keeping the first occurrence

Example output:

```text
04_iqtree_initial/input/HA/
    HA_reference_dedup_aligned_trimmed.fasta
    HA_metadata_reference_dedup_aligned_trimmed.xlsx
```

In `ALL` mode, the script can prepare all segments. In single-segment mode, it fails if the requested source file is missing.

---

## Step 2: Run IQ-TREE

```bash
python3 04_iqtree_initial/scripts/run_iqtree.py \
  04_iqtree_initial/input/HA/HA_reference_dedup_aligned_trimmed.fasta \
  -o 04_iqtree_initial/output
```

The script also supports the default HA path, so a simple command is often enough:

```bash
python3 04_iqtree_initial/scripts/run_iqtree.py
```

This default resolves approximately to:

```text
04_iqtree_initial/input/HA/HA_reference_dedup_aligned_trimmed.fasta
```

### Default options used by the script

- substitution model: `MFP`
- bootstrap replicates: `1000`
- threads: `AUTO`
- prefix: output directory + FASTA stem

So the command generated is effectively equivalent to:

```bash
iqtree2 -s input.fasta -m MFP -B 1000 -T AUTO --prefix output/HA_reference_dedup_aligned_trimmed
```

The script allows overriding these values with flags such as:

```bash
python3 04_iqtree_initial/scripts/run_iqtree.py \
  04_iqtree_initial/input/HA/HA_reference_dedup_aligned_trimmed.fasta \
  -B 2000 \
  -T 8 \
  -m MFP \
  -o 04_iqtree_initial/output
```

---

## Background execution

The script supports background execution:

```bash
python3 04_iqtree_initial/scripts/run_iqtree.py --background
```

This starts the IQ-TREE process in a detached session and writes:

- a PID file
- a background log file
- the final tree when the run completes

This is useful when running a long tree-building job and continuing with other work.

---

## Output files

IQ-TREE writes a standard set of result files into the selected output directory, for example:

```text
04_iqtree_initial/output/
    HA_reference_dedup_aligned_trimmed.treefile
    HA_reference_dedup_aligned_trimmed.iqtree
    HA_reference_dedup_aligned_trimmed.log
    HA_metadata_reference_dedup_aligned_trimmed.xlsx   # copied if present
```

The key files are:

- `.treefile` — the inferred best tree
- `.iqtree` — IQ-TREE summary report and model information
- `.log` — the run log

---

## Metadata copy behavior

When the FASTA file has a paired metadata workbook in the same input directory, the script automatically copies that workbook into the output directory.

This is handled by `copy_metadata_to_output()` in `run_iqtree.py`, which searches for an `.xlsx` or `.xls` file containing `metadata` in the filename and matching the same gene segment.

This is useful for preserving metadata alongside the final tree output for downstream interpretation.

---

## Important implementation details

- IQ-TREE is searched by command name using `iqtree2` or `iqtree`.
- If the default input FASTA is missing, the script attempts a limited automatic discovery under the input directory.
- If multiple candidate FASTA files are found, the script stops and asks the user to specify the exact file.
- The script resolves relative paths from the project root rather than the current shell directory.
- The default output directory is:

```text
04_iqtree_initial/output/
```

---

## Summary

Stage 04 runs the initial phylogenetic inference step. It prepares the trimmed alignments and metadata, runs IQ-TREE 2 for the selected segment, and writes the usual tree and report outputs needed for downstream evolutionary analysis.

The default project behavior is HA-oriented, but the scripts can be pointed at any segment or at all segments as needed.
