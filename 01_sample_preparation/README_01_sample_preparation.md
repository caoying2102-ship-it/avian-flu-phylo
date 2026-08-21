# 01 Sample Preparation

This stage prepares raw influenza metadata and FASTA inputs so that downstream steps can perform multiple-sequence alignment, trimming, and phylogenetic reconstruction on a consistent, cleaned, and sample-matched dataset.

This folder does not build trees by itself. Its job is to convert raw input files into a standardized set of per-gene files, perform initial QC, and create the final reference sample list used by the other stages.

## What this stage does

The actual script workflow in this project is:

1. Read all Excel metadata files and all FASTA files from the input directory.
2. Merge them into a combined metadata table and a combined FASTA file.
3. Detect the gene segment from each FASTA header (PB2, PB1, PA, HA, NP, NA, MP, NS).
4. Split sequences by gene and write per-segment files.
5. Deduplicate exact sequence duplicates within each segment.
6. Match FASTA isolate IDs to metadata by Isolate_Id.
7. Create a quality-control summary for each gene.
8. Summarize HA sample counts by year, subcontinent, and host category.
9. Perform stratified sampling on HA using that summary.
10. Save the selected HA isolate IDs to Data/reference ID.csv.
11. Extract the same isolates from each other segment.
12. Deduplicate the extracted reference sequences using SeqKit.

This means the project is designed around a single HA-based sampling reference set, and the other segments are then filtered to the same sampled isolates.

---

## Scripts in this folder

- `prepare_samples.py`
  - Main preprocessing script.
  - Merges metadata, merges FASTA records, identifies segments, deduplicates, matches IDs, and writes QC outputs.

- `run.sh`
  - Thin wrapper around `prepare_samples.py`.
  - This is the default first command for stage 01.

- `summary_count.py`
  - Counts HA samples by subcontinent, collection year, and host category.

- `sample.py`
  - Performs stratified HA sampling.
  - Writes the final `Data/reference ID.csv` isolate list.

- `extract.py`
  - Reads the reference isolate list and extracts matching metadata and sequences for each gene segment.

- `SeqKit.py`
  - Runs sequence deduplication on extracted segment FASTA files.

- `generate_simplified_metadata.py`
  - Produces a cleaned metadata table with standardized fields.

- `rename_fasta_headers.py`
  - Renames headers to match metadata after ID matching.

---

## Expected input

Place all raw inputs in:

```text
01_sample_preparation/input/
```

Supported metadata formats:

```text
.xls
.xlsx
```

Supported FASTA formats:

```text
.fasta
.fa
.fas
.fna
```

Required metadata columns include at least:

```text
Isolate_Id
Collection_Date
Subcontinent
Host
```

The pipeline also expects other common fields such as:

```text
Continent
Country
Region
Location
Isolate_Name
Clade
```

The FASTA headers need to include a gene name and an isolate ID in a parseable format. The script searches for gene names such as `PB2`, `PB1`, `PA`, `HA`, `NP`, `NA`, `MP`, `NS` and uses pipe-delimited headers with an isolate ID field.

---

## Main outputs

After running the default stage script, the main outputs are generated under:

```text
01_sample_preparation/output/
```

### Top-level outputs

```text
combined_metadata.csv
combined_metadata.xlsx
combined_sequences.fasta
qc_summary.csv
qc_summary.xlsx
input_manifest.csv
unassigned_sequences.fasta   # only if there are unrecognized FASTA headers
```

These files represent the merged and cleaned starting dataset.

### Per-gene outputs

For each segment, the script writes files under:

```text
01_sample_preparation/output/genes/{SEGMENT}/
```

Example:

```text
01_sample_preparation/output/genes/HA/
    HA.fasta
    HA_dedup.fasta
    HA_metadata.csv
    HA_metadata.xlsx
    HA_metadata_simplified.xlsx
```

The key point is that the project splits the data by gene before downstream sampling and selection.

---

## Stage 01 execution flow

### 1. Run the initial merge and split step

From the project root:

```bash
bash 01_sample_preparation/scripts/run.sh
```

This calls:

```bash
python3 01_sample_preparation/scripts/prepare_samples.py \
  --input-dir 01_sample_preparation/input \
  --output-dir 01_sample_preparation/output
```

This script does the initial work:

- reads metadata and FASTA inputs
- merges them into combined files
- detects gene segments from FASTA headers
- writes gene-specific FASTA files
- deduplicates exact sequence duplicates
- matches isolate IDs to metadata
- prints QC summary counts

The script is the backbone of this stage.

### 2. Prepare simplified metadata

```bash
python3 01_sample_preparation/scripts/generate_simplified_metadata.py
```

This script standardizes metadata columns and creates cleaned tables used later in sampling and header renaming.

### 3. Rename deduplicated FASTA headers

```bash
python3 01_sample_preparation/scripts/rename_fasta_headers.py
```

This script rewrites FASTA headers to a standardized ID format used for matching to metadata records.

### 4. Summarize HA sample structure

```bash
python3 01_sample_preparation/scripts/summary_count.py
```

This step works on HA only and groups samples by:

- subcontinent
- collection year
- host category

This summary is the basis for the sampling design.

### 5. Perform HA stratified sampling

```bash
python3 01_sample_preparation/scripts/sample.py
```

The sampling rules are implemented as follows:

- target size is 5200 samples
- groups with 5 or fewer samples are kept fully
- groups with more than 5 samples receive at least 5 samples
- remaining slots are allocated proportionally to remaining capacity
- random seed is fixed at 42 for reproducibility
- `EPI_ISL_1254` is retained by default

The output is:

```text
Data/reference ID.csv
```

This file becomes the master list of sampled isolates.

### 6. Extract the same isolates from all segments

```bash
python3 01_sample_preparation/scripts/extract.py
```

This script reads the reference IDs and extracts matching metadata and sequences from every segment. It writes files such as:

```text
{SEGMENT}_metadata_reference.xlsx
{SEGMENT}_reference.fasta
{SEGMENT}_missing_reference_IDs.csv
```

This is the point where the project enforces the same sampled isolate set across all gene segments.

### 7. Final deduplication with SeqKit

```bash
python3 01_sample_preparation/scripts/SeqKit.py
```

For each segment, the script removes duplicate sequences by sequence content and keeps the metadata aligned to the retained records. Final expected files are:

```text
{SEGMENT}_reference_dedup.fasta
{SEGMENT}_metadata_reference_dedup.xlsx
```

These are the cleaned reference files passed to downstream alignment and tree-building steps.

---

## Important logic in the project

### HA is the sampling anchor

The actual design in this repository is HA-first:

- HA is summarized
- HA is stratified and sampled
- the resulting HA isolate list is saved to `Data/reference ID.csv`
- each other segment is then filtered to the same set of isolate IDs

This means the phylogenetic analyses for other gene segments are not independent samplings. They use the same HA-derived sample set, only pulling the matching sequences for that gene.

### Sequence deduplication is strict

Both `prepare_samples.py` and the final SeqKit step deduplicate by sequence content. If multiple isolates share the same sequence, only the first occurrence is kept in the deduplicated result.

### Matching relies on isolate IDs

The pipeline links FASTA sequences to metadata by isolate IDs. If the ID format is inconsistent or missing, the later matching and extraction steps can fail or produce missing-ID reports.

### QC is explicit

The stage reports how many records were removed, how many were matched, and whether selected reference IDs are missing in each gene segment. These reports are essential for validating that the reference set is usable before later stages run.

---

## Typical workflow sequence

```bash
cd /Users/YingCao/Desktop/RecentWork/avian-flu-phylo

bash 01_sample_preparation/scripts/run.sh
python3 01_sample_preparation/scripts/generate_simplified_metadata.py
python3 01_sample_preparation/scripts/rename_fasta_headers.py
python3 01_sample_preparation/scripts/summary_count.py
python3 01_sample_preparation/scripts/sample.py
python3 01_sample_preparation/scripts/extract.py
python3 01_sample_preparation/scripts/SeqKit.py
```

The script set is intentionally staged so that quality issues can be inspected between steps.

---

## Summary

Stage 01 is the data-preparation and sampling stage. It turns raw sequence and metadata files into a cleaned, gene-split, deduplicated, and standardized reference dataset, with HA used as the sampling anchor for the final isolate set that is propagated across the other segments.

This stage does not perform alignment or phylogeny itself; it prepares the valid inputs for the later MSA, trimming, and IQ-TREE steps.
