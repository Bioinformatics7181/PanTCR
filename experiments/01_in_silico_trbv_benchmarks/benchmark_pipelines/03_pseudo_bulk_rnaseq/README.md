# Pseudo-Bulk RNA-seq Benchmark Pipeline

Purpose: generate pseudo-bulk/scBulk benchmark outputs consumed by `experiments/08_pseudo_bulk_rnaseq_benchmark`.

Canonical output path:

`experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/03_pseudo_bulk_rnaseq/generated/per_dataset_results`

Package-local drivers:

- `run_one_dataset.sh`: driver for one SC/dataset pair. It starts from either
  an existing MiXCR `clones_TRB.tsv` file or a 10x Gene Expression BAM, then
  writes the directory layout consumed by experiment 08.
- `run_all_datasets.sh`: aggregate driver for the SC-1 through SC-8 manuscript
  mapping. It calls `run_one_dataset.sh` for each dataset and writes
  `generated/per_dataset_results/scbulk_manuscript_mode_mapping.csv`.

Reusable implementation code is under `code/pantcr/` and
`code/experiments/00_benchmark_utils/`.

Workflow:

1. Download or place the paired 10x Gene Expression BAM and targeted V(D)J proxy-truth outputs for each dataset under this folder's `input/` tree.
2. Extract TRB-region reads from the GEX BAM with `samtools`, convert the extracted reads to FASTQ, and run `mixcr analyze rna-seq`.
3. Run MiXCR contig assembly and `findAlleles` to obtain MiXCR and FindAlleles baseline calls.
4. Convert MiXCR clone tables into PanTCR mutation-observation CSVs with `code/pantcr/collect_mutation.py`.
5. Run PanTCR inference with the selected graph roots. The Figure 4E/S18/S19
   branch uses the default pseudo-bulk PanTCR inference settings defined in
   `run_one_dataset.sh`.
6. Evaluate predictions against targeted V(D)J proxy truth and write per-dataset result folders.

Required tools:

- `samtools`
- `mixcr`
- Python dependencies from `requirements.txt`

Dataset source:

The SC-1 through SC-8 datasets are the paired Gene Expression and targeted
V(D)J 10x Genomics datasets listed in the manuscript Data Availability
statement and Supplementary Table S1. Use the official 10x Genomics download
links in that table to obtain the GEX BAM or equivalent FASTQ/clone-table input
and the matched targeted V(D)J files used as proxy truth.

Downstream use:

`experiments/08_pseudo_bulk_rnaseq_benchmark/run_rebuild.sh` uses `generated/per_dataset_results` from this pipeline as `SC_BASE`.

Example lightweight rerun from existing MiXCR outputs for the final pseudo-bulk branch:

```bash
bash run_one_dataset.sh \
  --sc-id SC-1 \
  --dataset-id dataset2 \
  --clones input/SC-1_dataset2/trb_all.clones_TRB.tsv \
  --alleles-tsv input/SC-1_dataset2/trb_all.alleles.tsv \
  --alleles-json input/SC-1_dataset2/trb_all.customAlleles.json \
  --genotype input/SC-1_dataset2/genotype.csv
```

To run the manuscript SC-1 through SC-8 mapping after placing inputs under `input/SC-X_datasetY/`:

```bash
bash run_all_datasets.sh
```

The aggregate driver writes `scbulk_manuscript_mode_mapping.csv` in
`generated/per_dataset_results/`, which is the mapping file consumed by
experiment 08.
