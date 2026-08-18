# Avian Influenza Phylogenetic Tree Pipeline

This pipeline processes one influenza gene segment end-to-end:

1. `01_sample_preparation`: assemble raw inputs, clean metadata, and prepare segment-specific files
2. `02_msa`: align sequences with MAFFT and link metadata to aligned outputs
3. `03_trimal`: trim the alignment with trimAl
4. `04_iqtree_initial`: build the phylogenetic tree with IQ-TREE

Supported segments: `HA`, `NA`, `MP`, `PB2`, `PB1`, `PA`, `NP`, and `NS`.

## 1. Install Conda and create the environment

A Miniconda installer is included in the project root:

- `Miniconda3-latest-MacOSX-arm64.sh`
- `install_miniconda.sh`

On macOS, install Miniconda with:

```bash
chmod +x install_miniconda.sh
./install_miniconda.sh
```

After installation, reopen the terminal and activate Conda:

```bash
source ~/miniconda3/bin/activate
conda init zsh
```

Then start a new shell and create the project environment:

```bash
conda activate base
conda env create -f environment.yml
conda activate iqtree-pipeline
```

If the environment already exists, update it instead:

```bash
conda activate iqtree-pipeline
conda env update -f environment.yml
```

The environment defined in [environment.yml](environment.yml) includes the required Python packages and external tools:

- Python 3.11
- pandas
- openpyxl
- Biopython
- MAFFT
- trimAl
- IQ-TREE
- SeqKit
- PyYAML

## 2. Prepare raw input files

This repository only contains the pipeline code. Raw data are not stored in GitHub and must be supplied locally.

Place your input FASTA files and metadata Excel files directly in:

```text
Data/
```

The project accepts:

- metadata: `.xlsx` or `.xls`
- FASTA: `.fa`, `.fas`, `.fasta`, `.fna`

The script will automatically copy supported files from `Data/` into `01_sample_preparation/input/` before running the pipeline.

Manual refresh is also available:

```bash
python 01_sample_preparation/scripts/copy_data_to_input.py
```

To overwrite existing files in `01_sample_preparation/input/`:

```bash
python 01_sample_preparation/scripts/copy_data_to_input.py --overwrite
```

Important input rules:

- `Data/` must contain at least one metadata file and at least one FASTA file
- FASTA headers must include a recognizable segment name and isolate ID
- metadata must include the required columns used by the scripts, especially `Isolate_Id`

## 3. Example data for first-time testing

The repository includes an example dataset for a quick test of the file format and workflow.

Example layout:

```text
Data/
├── example1.fasta
├── example2.fasta
├── example_metadata.xlsx
├── README.md
└── other local files
```

This example is intended only for verifying the pipeline format, not for full research use.

The key points are:

- the FASTA may be split across multiple files
- there is one paired metadata file
- each FASTA header should contain the segment name and isolate ID so metadata can be matched to sequence records correctly

## 4. Run the pipeline

From the project root:

```bash
python run_pipeline.py --segment HA
```

You may also use the positional form:

```bash
python run_pipeline.py HA
```

Lowercase segment names are accepted and normalized to uppercase:

```bash
python run_pipeline.py --segment pb2
```

You cannot use both positional and `--segment` at the same time.

If you want to archive prior outputs before a fresh run:

```bash
python run_pipeline.py --segment HA --clean
```

This moves any previous segment outputs under `logs/archive/<timestamp>/` before the new run starts.

## 5. What the pipeline produces

The pipeline writes intermediate and final outputs for each selected segment under:

```text
04_iqtree_initial/output/SEGMENT/
```

For each run, the output directory contains the IQ-TREE results, including:

- `.treefile`
- `.iqtree`
- `.log`
- the matched metadata workbook copied alongside the final tree outputs

The metadata copied into the tree output directory is the paired `*_metadata_reference_dedup_aligned_trimmed.xlsx` file that matches the FASTA used for tree inference.

## 6. Execution flow

The controller runs the stages in order:

```text
01_sample_preparation/
  prepare_samples.py
  generate_simplified_metadata.py
  rename_fasta_headers.py
  summary_count.py
  sample.py
  extract.py
  SeqKit.py

02_msa/
  prepare_msa_inputs.py
  msa.py
  matchmetadata.py

03_trimal/
  prepare_trimal_inputs.py
  trimal.py

04_iqtree_initial/
  prepare_iqtree_inputs.py
  run_iqtree.py
```

## 7. Troubleshooting

- Missing dependency: activate the environment and check `mafft`, `trimal`, `seqkit`, and `iqtree2`
- Raw files not found: confirm they are placed in `Data/`, not directly in the stage input directories
- Input copy step skipped: verify the Data files have supported file extensions
- Previous output conflicts: use `--clean` or remove old results before rerunning
- Segment not recognized: use one of the supported values from `config.yaml`

## 8. Static checks

You can validate the project without reading real data:

```bash
python -m compileall -q run_pipeline.py 01_sample_preparation/scripts 02_msa/scripts 03_trimal/scripts 04_iqtree_initial/scripts
python run_pipeline.py --help
```

These checks only verify syntax and command-line usage; they do not run the pipeline on actual data.
