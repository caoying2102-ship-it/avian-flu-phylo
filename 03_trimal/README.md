# 03 TrimAl

This stage trims the MSA outputs produced in stage 02 using trimAl. It does not build trees; it removes poorly aligned or ambiguous columns from each per-gene alignment while preserving the sequence set and copying the same metadata into the trimmed output.

The actual workflow implemented by the project is:

1. Copy the aligned FASTA and metadata files from 02_msa/output into 03_trimal/input.
2. Validate that all required files exist and are non-empty.
3. Confirm the alignment has equal sequence lengths.
4. Run trimAl with `-automated1`.
5. Validate that the trimmed alignment length is not longer than the original alignment.
6. Copy metadata output alongside the trimmed FASTA.

---

## Scripts in this folder

- `prepare_trimal_inputs.py`
  - Copies aligned FASTA and metadata from stage 02 into the TrimAl input directory.

- `trimal.py`
  - Validates input files, runs trimAl, checks alignment consistency, and writes trimmed outputs.

---

## Expected inputs

This stage expects the aligned outputs from stage 02, specifically:

```text
02_msa/output/{SEGMENT}/
    {SEGMENT}_reference_dedup_aligned.fasta
    {SEGMENT}_metadata_reference_dedup_aligned.xlsx
```

These are copied into:

```text
03_trimal/input/{SEGMENT}/
```

Supported segments:

```text
MP, HA, NA, PB2, PB1, PA, NP, NS
```

---

## Step 1: Prepare TrimAl inputs

```bash
python3 03_trimal/scripts/prepare_trimal_inputs.py --segment HA
```

This script copies the aligned FASTA and aligned metadata from the MSA stage into the TrimAl input directory. It checks whether the source files exist before copying, and it fails early if they are missing or empty.

Example destination:

```text
03_trimal/input/HA/
    HA_reference_dedup_aligned.fasta
    HA_metadata_reference_dedup_aligned.xlsx
```

---

## Step 2: Trim alignment with TrimAl

```bash
python3 03_trimal/scripts/trimal.py --segment HA
```

This script performs the actual trimming step for one gene segment.

### What it does

- checks whether TrimAl is installed
- validates input files exist
- verifies the alignment is a valid alignment with equal sequence lengths
- runs TrimAl using:

```bash
trimal -in input.fasta -out output.fasta -automated1
```

- validates the trimmed output is not longer than the input alignment
- warns if the sequence count changes, which is allowed when gap-only sequences are removed
- copies the original metadata table to the trimmed output directory unchanged

Output files:

```text
03_trimal/output/HA/
    HA_reference_dedup_aligned_trimmed.fasta
    HA_metadata_reference_dedup_aligned_trimmed.xlsx
    HA_trimal.log
```

---

## Alignment validation logic

The script is strict about the alignment structure before trimming:

- all sequences in the FASTA must have the same length
- the input FASTA must contain sequences
- the trimmed alignment must not expand in length
- if the sequence count changes after trimming, the code prints a warning and proceeds, because trimAl may remove gap-only sequences

This means the stage focuses on preserving valid sequence alignments while removing low-information positions.

---

## Metadata handling

TrimAl only trims columns in the alignment; it does not modify sample identity. Therefore, after trimming the FASTA, the script copies the metadata file directly to the output folder:

```text
{SEGMENT}_metadata_reference_dedup_aligned_trimmed.xlsx
```

This keeps metadata and sequence order aligned for the next stage.

---

## Output structure

Typical layout:

```text
03_trimal/
├── input/
│   ├── HA/
│   │   ├── HA_reference_dedup_aligned.fasta
│   │   └── HA_metadata_reference_dedup_aligned.xlsx
│   └── ...
├── output/
│   ├── HA/
│   │   ├── HA_reference_dedup_aligned_trimmed.fasta
│   │   ├── HA_metadata_reference_dedup_aligned_trimmed.xlsx
│   │   └── HA_trimal.log
│   └── ...
├── scripts/
│   ├── prepare_trimal_inputs.py
│   └── trimal.py
└── README.md
```

---

## Practical notes

- This stage assumes TrimAl is installed and available on the PATH.
- The pipeline runs segment-by-segment, not all segments in one call.
- The trimming approach is `-automated1`, which chooses a standardized trimming strategy automatically.
- The project does not re-filter the metadata during trimming; it simply preserves the metadata table alongside the trimmed sequence alignment.

---

## Summary

Stage 03 is the alignment-cleaning step. It takes the validated per-gene MSA outputs from stage 02, removes poorly aligned columns with TrimAl, preserves the same sample identities, and creates the trimmed alignment files that are ready for IQ-TREE tree inference.
