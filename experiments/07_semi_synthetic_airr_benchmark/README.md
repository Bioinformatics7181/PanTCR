# 07 Semi-Synthetic AIRR Benchmark

Purpose: provide the semi-synthetic AIRR benchmark workflow, Supplementary Table S16 metric rebuild, and retained-evidence diagnostics used by the final Figure 4 / Supplementary Table S18 provenance chain.

Implementation:

- `code/experiments/07_semi_synthetic_airr_benchmark/semi_simu_fine_evidence_analysis.py`: stratifies non-default truth allele recovery by observed default-relative defining changes and retained PanTCR sample-level observations.

Evaluation truth: this benchmark uses full-length AIRR-derived allele calls as proxy truth.

Key parameters: retained PanTCR sample-level diagnostic summaries were filtered using the semi/real-data `NaiveDiversityIndex >= 2` setting.

Supplementary Table S18 is rebuilt by `code/experiments/08_pseudo_bulk_rnaseq_benchmark/summarize_figure4_observed_splits.py`, which recomputes the method-comparison observed-region stratum from unfiltered MiXCR-derived observed ranges (`min_naive=None`) and combines the semi-synthetic and pseudo-bulk rows in `experiments/08_pseudo_bulk_rnaseq_benchmark/generated/final_two_table_audit/method_comparison_by_observed_defining_site.csv`.

Inputs: full reruns default to the expanded benchmark output root `experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/02_semi_synthetic_airr/generated/results`, with references from `data/ref`. Use `--results-root` for another expanded full-provenance result root.

The package-local semi-synthetic driver is under
`experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/02_semi_synthetic_airr`.
It derives proxy genotype labels from each sample's full-length MiXCR
FindAlleles outputs (`*.alleles.tsv` and `*.customAlleles.json`) before running
the degraded-read benchmark. The raw/full-length AIRR source files are expected
under that pipeline's `input/full_length_airr` tree and are obtained from the
AIRR-1/AIRR-2 public accessions reported in Supplementary Table S1.

## Run

First generate the upstream semi-synthetic result root:

```bash
cd experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/02_semi_synthetic_airr
bash run_all_airr_benchmarks.sh
cd ../../../..
```

Then rebuild the experiment 07 retained-evidence analysis:

```bash
bash experiments/07_semi_synthetic_airr_benchmark/run_rebuild.sh
```

## Outputs And Downstream Use

Outputs are written under `generated/evidence_analysis`. The key file reused by
experiment 08 is `generated/evidence_analysis/per_truth_call_status.csv`.
