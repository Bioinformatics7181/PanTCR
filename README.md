# PanTCR: Local Pan-clonotype Graphs Enable Robust Inference and Discovery of TCR Alleles from Fragmented Transcriptomes

This repository contains the implementation and paper-reproduction workflows
for **PanTCR**, a method for mature TRB V/J coding-region allele inference from
fragmented transcriptomic TCR observations using population-aware local
pan-allele graph priors.

![PanTCR workflow](figs/PanTCR.png)

PanTCR starts from MiXCR-processed TRB rearrangement tables, converts clone
observations into V/J allele evidence, and infers sample-level alleles with a
local pan-allele graph prior. The repository also includes the scripts used to
reproduce the synthetic, semi-synthetic, pseudo-bulk, and public bulk analyses
reported in the manuscript and Supplementary Tables.

## Prerequisites

- Linux or macOS is recommended for the bash experiment drivers.
- Python 3.8 or later.
- Conda or another isolated Python environment manager is recommended.
- Node.js 18 or later for lightweight `.mjs` summary scripts in experiments
  04, 06, and 08. The tested major version is recorded in `.node-version`;
  these scripts use only the Node.js standard library, so no `npm install` or
  package lock is required.
- Java and MiXCR v4 or later for clone-table generation and baseline calls.
- ART for in silico Illumina read simulation.
- samtools and SRA Toolkit for pseudo-bulk and public bulk RNA-seq reruns.

## Installation

Create an isolated environment and install the Python requirements:

```bash
conda create -n pantcr python=3.9 -y
conda activate pantcr
pip install -r requirements.txt
```

Install the bioinformatics tools needed by the workflows you plan to
run. For example:

```bash
conda install -c bioconda art samtools sra-tools
```

Install MiXCR from the official distribution or Bioconda and make it
available on `PATH`. Check the runtime before running the paper workflows:

```bash
python --version
node --version
mixcr --version
art_illumina --help
samtools --version
```

## Data Preparation

The lightweight references, example input, and pretrained graph material are
bundled under `data/`:

- `data/ref/`: TRB reference files and the MiXCR-all custom library.
- `data/examples/sc1/`: a small MiXCR-processed TRB clone-table example.
- `data/pretrained_graphs/`: pretrained graph material for the example run.

Before running simulation workflows, unpack the CDR3 dictionary:

```bash
unzip data/ref/cdr3_dict.zip -d data/ref/
```

If `unzip` is unavailable:

```bash
python -m zipfile -e data/ref/cdr3_dict.zip data/ref/
```

Obtain large public datasets according to the manuscript Data
Availability statement and Supplementary Table S1. The main public sources are:

- AIRR-1: ArrayExpress `E-MTAB-13593`
- AIRR-2: ENA `PRJEB28370`
- SC-1 to SC-8: official 10x Genomics datasets listed in Supplementary Table S1
- BULK-1: GEO `GSE112087`
- BULK-2: GEO `GSE107011`

Experiment-specific README files describe the expected `input/` layout for
large raw reads, MiXCR intermediates, and generated benchmark result roots.

## Quick Start: Run PanTCR On One Sample

First generate a TRB clone table with MiXCR:

```bash
mixcr analyze -s hsa rna-seq \
  sample_R1.fastq.gz sample_R2.fastq.gz \
  sample_mixcr
```

The bundled SC-1 example already includes a processed TRB clone table. Run
PanTCR with the provided pretrained graph:

```bash
python code/pantcr/run_sample.py \
  --clones data/examples/sc1/SC1.clones_TRB.tsv \
  --sample-id SC1 \
  --out-dir data/examples/sc1/output \
  --graph-root data/pretrained_graphs/trb_airr1 \
  --reference-dir data/ref
```

The output directory contains coordinate-projected V/J observations, inferred
allele calls, and retained-evidence versus graph-imputation audit tables.

## Build A New PanTCR Graph

For a new cohort, first convert MiXCR clone tables into V/J observation CSVs
with `code/pantcr/collect_mutation.py` or `code/pantcr/run_sample.py`. Then
build separate V and J graphs:

```bash
python code/pantcr/build_pangenome_graph.py \
  --in_dir cohort_observations/V \
  --out_dir custom_graph/V \
  --metadata metadata.csv \
  --write_paths

python code/pantcr/build_pangenome_graph.py \
  --in_dir cohort_observations/J \
  --out_dir custom_graph/J \
  --metadata metadata.csv \
  --write_paths
```

Minimal metadata schema:

```csv
sample_id,population_id,filename
sample_001,AFR,sample_001.V.csv
sample_002,EUR,sample_002.V.csv
```

`sample_id` and `population_id` are required. `filename` is optional when the
observation-table filename differs from `sample_id`.

## Repository Layout

- `code/pantcr/`: core PanTCR components for mutation-evidence extraction,
  graph construction, allele inference, and evidence auditing.
- `code/simulation/`: genotype sampling and controlled repertoire simulation.
- `code/experiments/`: shared benchmark utilities and experiment-specific
  result builders.
- `data/`: lightweight references, example input, and pretrained graph files.
- `experiments/`: paper experiment drivers and experiment-specific README
  files.
- `figs/`: workflow and documentation figures.

Generated outputs are written under experiment-specific `generated/`, `results/`,
or `runs/` directories.

## Reproduce Paper Experiments

Each folder under `experiments/` contains a README with inputs, commands,
parameters, expected outputs, and the corresponding manuscript or
Supplementary Table role.

The main in silico benchmark driver covers Scenario A, Scenario B, Scenario C,
read-depth/RPC settings, alpha settings, fragmentation/truncation settings, and
the high-quality full-length simulation configuration:

```bash
bash experiments/01_in_silico_trbv_benchmarks/run_all_scenarios.sh
```

The upstream benchmark hub also contains pipelines reused by downstream
analyses:

```text
experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/02_semi_synthetic_airr/
experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/03_pseudo_bulk_rnaseq/
experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/04_public_bulk_rnaseq/
```

After the relevant upstream result roots have been generated, downstream
analyses can be rebuilt with their experiment-level drivers, for example:

```bash
bash experiments/02_prediction_and_paralogy_audits/run_rebuild.sh
bash experiments/03_leave_target_allele_out/run_leaveout_allele_specific.sh
bash experiments/04_naive_diversity_threshold_sensitivity/run_rebuild.sh
bash experiments/05_public_bulk_candidate_nonimgt/run_rebuild.sh
bash experiments/06_population_prior_robustness/run_rebuild.sh
bash experiments/07_semi_synthetic_airr_benchmark/run_rebuild.sh
bash experiments/08_pseudo_bulk_rnaseq_benchmark/run_rebuild.sh
```

## License

This project is released under the MIT License. See `LICENSE`.

## Contact

Xinyang Qian: qianxy@stu.xjtu.edu.cn
