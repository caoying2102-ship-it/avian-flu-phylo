# avian-flu-phylo

This project runs one influenza gene segment through the following stages:

1. `01_sample_preparation`: merge inputs, clean metadata, sample, extract, and deduplicate
2. `02_msa`: prepare inputs, run MAFFT alignment, and match metadata to the aligned sequences
3. `03_trimal`: trim alignments with trimAl
4. `04_iqtree_initial`: infer a phylogenetic tree with IQ-TREE 2

The pipeline supports `HA`, `NA`, `MP`, `PB2`, `PB1`, `PA`, `NP`, and `NS`. A single run only reads and writes the selected segment's intermediate files and does not require the other seven segments to exist.

## Fresh runs and existing outputs

If you want each run to start from a clean state, use the new `--clean` option. It archives any existing intermediate/output directories for the selected segment under `logs/archive/<timestamp>/` before running.

Example:

```bash
python run_pipeline.py --segment HA --clean
```

This preserves prior result files in `logs/archive/` while letting the current execution start without stale outputs.

## Installation

Conda/Mamba is recommended:

```bash
conda env create -f environment.yml
conda activate iqtree-pipeline
```

If someone else clones this repository from GitHub, they can use the same commands to set up the environment.

For macOS or Linux, the `environment.yml` installs both the required Python packages and the external bioinformatics tools:
- `mafft`
- `trimal`
- `iqtree`
- `seqkit`

If the user already has a Conda environment, they can activate it first and then install the project dependencies with:

```bash
conda env update -f environment.yml
```

### Setup after cloning from GitHub

After cloning the repository, use the same environment file to reproduce the required tools and Python packages:

```bash
git clone https://github.com/<your-username>/<repo>.git
cd <repo>
conda env create -f environment.yml
conda activate iqtree-pipeline
```


```bash
conda activate phylo_arm
conda env update -f environment.yml
```

The controller checks Python modules and the executables `mafft`, `trimal`, `seqkit`, and IQ-TREE before running. IQ-TREE supports both `iqtree2` and `iqtree`. Dependencies and supported segments are configured in `config.yaml`.

### Data policy

This repository does not include raw `Data/` or `01_sample_preparation/input/` files. Those directories are locally ignored in Git to avoid uploading large GISAID FASTA and metadata files.

## Input

This repository contains only pipeline code. Raw GISAID input files are not uploaded to GitHub and must be obtained separately.

Place the downloaded metadata Excel file(s) and FASTA file(s) in:

```text
Data/
```

Then copy them into the pipeline input directory before running the workflow:

```bash
python 01_sample_preparation/scripts/copy_data_to_input.py
```

If you need to replace old files in `01_sample_preparation/input/`, add `--overwrite`:

```bash
python 01_sample_preparation/scripts/copy_data_to_input.py --overwrite
```

The pipeline expects at least one Excel metadata file (`.xls` or `.xlsx`) and at least one FASTA file (`.fa`, `.fas`, `.fasta`, or `.fna`) in `01_sample_preparation/input/`.

>>>>>>> 4e80ef2 (Update data policy and add copy_data_to_input helper)
FASTA headers must contain a recognizable segment name and isolate ID; metadata must include the columns required by the scripts, especially `Isolate_Id`.

## Running

Recommended usage:

```bash
python run_pipeline.py --segment HA
```

Positional arguments and lowercase input are also supported:

```bash
python run_pipeline.py NA
python run_pipeline.py --segment pb2
```

You cannot use positional arguments and `--segment` at the same time. If no segment is specified or the segment is not in the supported list, the program exits before processing data.

## Execution order

```text
01 prepare_samples.py --genes SEGMENT
   generate_simplified_metadata.py --segment SEGMENT
   rename_fasta_headers.py --segment SEGMENT
   summary_count.py --segment SEGMENT
   sample.py --segment SEGMENT
   extract.py --segment SEGMENT
   SeqKit.py --segment SEGMENT

02 prepare_msa_inputs.py --segment SEGMENT
   msa.py --segment SEGMENT
   matchmetadata.py --segment SEGMENT

03 prepare_trimal_inputs.py --segment SEGMENT
   trimal.py --segment SEGMENT

04 prepare_iqtree_inputs.py --segment SEGMENT
   run_iqtree.py input/SEGMENT/... --outdir output/SEGMENT
```

Final results are written to:

```text
04_iqtree_initial/output/SEGMENT/
```

The main output files include `.treefile`, `.iqtree`, and `.log`.

## Configuration

`config.yaml` is the single source of segment and software configuration for the controller. It contains:

- `segments`: supported segments
- `executables`: candidate commands for external tools
- `python_modules`: Python modules that are checked before execution

When changing the segment list, make sure the stage scripts support the same segment names.

## Troubleshooting

- Missing runtime dependencies: after activating the environment, check `which mafft`, `which trimal`, `which seqkit`, and `which iqtree2`.
- Initial inputs not found: confirm that the files are placed directly in `01_sample_preparation/input/`.
- Previous-stage outputs not found: run the full controller from the project root and inspect the earliest failing stage.
- Existing IQ-TREE results with the same prefix: move the old outputs or call the IQ-TREE script directly with `--redo` as needed.

## Static checks

You can run the following without using real data:

```bash
python -m compileall -q run_pipeline.py 01_sample_preparation/scripts 02_msa/scripts 03_trimal/scripts 04_iqtree_initial/scripts
python run_pipeline.py --help
```

These commands only check syntax and the command-line interface; they do not run the pipeline or read input data.
