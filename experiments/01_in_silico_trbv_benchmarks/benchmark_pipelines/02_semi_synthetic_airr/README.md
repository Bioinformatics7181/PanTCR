# Semi-Synthetic AIRR Benchmark Pipeline

Purpose: generate the semi-synthetic AIRR benchmark result root consumed by `experiments/07_semi_synthetic_airr_benchmark`.

Canonical output path:

`experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/02_semi_synthetic_airr/generated/results`

Expected output layout:

- `results/labels/{expr}/MIX/genotype_*.csv`
- `results/mutations/{expr}/MIX/{V,J}/*_sequences.csv`
- `results/infer/{expr}/{MiXCR,FindAlleles,Bayes,BayesNoPrior}/...`
- `results/eval/{expr}/...`

Package-local pipeline:

- `run_pipeline.sh`: end-to-end semi-synthetic driver.
- `run_all_airr_benchmarks.sh`: batch wrapper that runs the benchmark settings
  consumed by downstream semi-synthetic analyses (`expr_AIRR1` and `expr_AIRR2` by
  default).
- `run_degraded_airr_batch.sh`: degrades paired AIRR FASTQs, derives proxy
  genotype labels from full-length MiXCR FindAlleles outputs, and reruns MiXCR
  on degraded reads.
- `collect_airr_mutations.sh`, `copy_airr_genotypes.sh`,
  `infer_airr_mixcr.sh`, `infer_airr_findalleles.sh`,
  `infer_airr_pantcr.sh`, and `eval_airr_alleles.sh`: sample-name-driven
  downstream steps for AIRR sample IDs such as `ERR*`.

Required input layout:

Place the full-length AIRR source material under:

`experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/02_semi_synthetic_airr/input/full_length_airr/{POP}/{SAMPLE}/`

Each sample directory includes paired FASTQs (`*_1.fastq.gz` and
`*_2.fastq.gz`, or uncompressed FASTQ) plus full-length MiXCR FindAlleles
outputs (`{SAMPLE}.alleles.tsv` and `{SAMPLE}.customAlleles.json`) used to
derive proxy genotype labels.

The AIRR source datasets are public datasets listed in the manuscript Data
Availability statement and Supplementary Table S1: AIRR-1 is ArrayExpress
`E-MTAB-13593`, and AIRR-2 is ENA `PRJEB28370`. Download the full-length AIRR
FASTQs from those repositories, run MiXCR/FindAlleles on the full-length reads
to obtain the proxy-truth allele files above, and place the files in the input
layout shown here.

Downstream use:

Run the complete semi-synthetic pipeline so the canonical `generated/results`
directory exists. Then `experiments/07_semi_synthetic_airr_benchmark/run_rebuild.sh`
will read that generated directory by default.

Example:

```bash
bash run_all_airr_benchmarks.sh
```

To run one setting, use `run_pipeline.sh --expr-id expr_AIRR1` or
`run_pipeline.sh --expr-id expr_AIRR2`. The default single-run expression ID is
`expr_AIRR1`.
