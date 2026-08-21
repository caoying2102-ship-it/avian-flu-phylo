# 02 MSA

This stage prepares and aligns the cleaned reference sequences from each gene segment using MAFFT. It does not perform phylogenetic tree inference; it converts the per-segment reference FASTA files into aligned FASTA files that are ready for trimming and downstream IQ-TREE analysis.

The real workflow in this project is:

1. Copy the filtered per-segment reference files from the 01 stage into 02_msa/input.
2. Validate that the expected FASTA files exist.
3. Remove sequences containing non-ATGC characters.
4. Remove manually excluded isolate IDs.
5. Keep only the longest valid sequence for each Isolate_Id.
6. Run MAFFT alignment for the cleaned sequence set.
7. Match aligned FASTA results back to metadata and save aligned outputs.

---

## Scripts in this folder

- `prepare_msa_inputs.py`
  - Copies the per-segment FASTA and metadata files from stage 01 into the MSA input directory.

- `msa.py`
  - Main MSA script.
  - Validates gene inputs, cleans invalid sequences, runs MAFFT, and writes aligned FASTA files.

- `matchmetadata.py`
  - Matches aligned FASTA IDs back to metadata and writes metadata aligned to the FASTA order.

---

## Expected input files

This stage expects the stage-01 reference files under:

```text
01_sample_preparation/output/genes/{SEGMENT}/
```

The copied MSA inputs are placed in:

```text
02_msa/input/{SEGMENT}/
```

Expected source files include:

```text
{SEGMENT}_reference_dedup.fasta
{SEGMENT}_metadata_reference_dedup.xlsx
```

The files are copied per segment for:

```text
MP, HA, NA, PB2, PB1, PA, NP, NS
```

---

## Step 1: Prepare MSA inputs

```bash
python3 02_msa/scripts/prepare_msa_inputs.py --segment HA
```

This script copies the selected segment files into the MSA input directory. It validates the source files first and raises an error if the expected files are missing.

The destination structure is:

```text
02_msa/input/HA/
    HA_reference_dedup.fasta
    HA_metadata_reference_dedup.xlsx
```

This is the stage used to make each segment ready for independent processing.

---

## Step 2: Run MAFFT per segment

```bash
python3 02_msa/scripts/msa.py --segment HA
```

This script performs the real MSA processing for a single segment.

### What it does

- checks for MAFFT installation
- reads the input FASTA file
- strips whitespace in sequences
- extracts `EPI_ISL_...` IDs from FASTA headers
- removes manually excluded IDs from a hard-coded set:

```text
EPI_ISL_20075028
EPI_ISL_20075029
EPI_ISL_20075031
EPI_ISL_20075045
EPI_ISL_20075046
```

- removes records with:
  - missing isolate IDs
  - empty sequences
  - invalid characters not in `ATGC`
- keeps only the longest valid sequence for each isolate ID when duplicates exist
- writes a cleaned FASTA file
- writes a CSV report of removed records
- runs MAFFT with:

```bash
mafft --auto --thread -1 --inputorder input.fasta
```

The final output is:

```text
02_msa/output/HA/HA_reference_dedup_aligned.fasta
```

Additional QC files are also produced:

```text
02_msa/output/HA/HA_reference_dedup_valid.fasta
02_msa/output/HA/HA_removed_invalid_sequences.csv
02_msa/output/HA/HA_mafft.log
```

---

## Important filtering rules

The actual MSA script enforces several strict rules:

1. Only sequences composed entirely of `A`, `T`, `G`, `C` characters are kept.
2. Any record without an `EPI_ISL_` ID is removed.
3. Manually excluded isolate IDs are removed before alignment.
4. If the same isolate ID appears more than once, the longest valid sequence is retained.
5. If the output alignment changes the number of records, the script raises an error.

This is a strong QC stage, not just a raw MAFFT wrapper.

---

## Step 3: Match aligned sequences to metadata

```bash
python3 02_msa/scripts/matchmetadata.py --segment HA
```

This script takes the aligned FASTA and the metadata file from the MSA input directory, then checks that every FASTA isolate ID exists in the metadata.

### It does the following

- reads `EPI_ISL_...` IDs from the aligned FASTA headers
- reads the metadata Excel file
- validates the `Isolate_Id` column
- removes duplicate metadata IDs by keeping the first row
- merges FASTA and metadata by isolate ID in FASTA order
- raises an error if any FASTA ID is missing from metadata
- writes a metadata table matching the aligned FASTA order

Output:

```text
02_msa/output/HA/HA_metadata_reference_dedup_aligned.xlsx
```

This ensures the sequence order and metadata row order stay synchronized.

---

## Output structure

Typical output directory layout:

```text
02_msa/
├── input/
│   ├── HA/
│   │   ├── HA_reference_dedup.fasta
│   │   └── HA_metadata_reference_dedup.xlsx
│   └── ...
├── output/
│   ├── HA/
│   │   ├── HA_reference_dedup_valid.fasta
│   │   ├── HA_reference_dedup_aligned.fasta
│   │   ├── HA_removed_invalid_sequences.csv
│   │   ├── HA_mafft.log
│   │   └── HA_metadata_reference_dedup_aligned.xlsx
│   └── ...
└── scripts/
    ├── prepare_msa_inputs.py
    ├── msa.py
    └── matchmetadata.py
```

---

## Notes about the implementation

- The project uses MAFFT with `--auto` and all available CPU threads (`MAFFT_THREADS = -1`).
- The script expects MAFFT to be installed in the environment.
- The validation logic is intentionally strict: it fails early if inputs are missing, headers are malformed, or IDs do not match metadata.
- The output FASTA is sequence-order-preserving and the metadata file is reordered to match the aligned FASTA order.

---

## Summary

Stage 02 is the per-gene MSA preparation stage. It takes the cleaned reference files selected in stage 01, removes invalid or duplicate records, aligns them with MAFFT, and ensures the aligned FASTA and metadata remain matched and ready for the trimming stage.

At the end of this stage, each gene segment has an aligned FASTA file and a corresponding aligned metadata table suitable for subsequent trimming and phylogenetic analysis.
