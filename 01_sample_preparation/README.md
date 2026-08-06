# 01 Sample Preparation

This stage prepares influenza sequence and metadata inputs for downstream analysis.

## Input

Place your metadata Excel file and FASTA sequence file in the input directory.

## Output

The stage generates:

- combined metadata CSV/XLSX
- combined sequence FASTA
- per-segment sequence files under output/genes
- QC summary reports

## Scripts

- prepare_samples.py
- generate_simplified_metadata.py
- rename_fasta_headers.py
- summary_count.py
- sample.py
- extract.py
- SeqKit.py
