# Experiment Analysis Code

This directory contains code used to rebuild manuscript tables, benchmark
summaries, and evidence-audit sources.

## Directory Guide

- `00_benchmark_utils/`: shared evaluation and baseline-conversion utilities
  used by multiple benchmark folders.
- `02_prediction_and_paralogy_audits/`: builders for count matrices,
  normalized provenance inputs, observed-region recovery summaries,
  sequence-only primary metrics, and gene-label/paralogy audits.
- `03_leave_target_allele_out/`: preparation and summary builders for the
  leave-target-allele-out stress test.
- `04_naive_diversity_threshold_sensitivity/`: inference-only threshold
  sensitivity builders.
- `05_public_bulk_candidate_nonimgt/`: public bulk candidate non-IMGT evidence
  and imputation annotation.
- `06_population_prior_robustness/`: graph-prior robustness builders.
- `07_semi_synthetic_airr_benchmark/`: retained-evidence diagnostics for
  semi-synthetic AIRR outputs.
- `08_pseudo_bulk_rnaseq_benchmark/`: pseudo-bulk/scBulk evidence, Figure 4E,
  and Supplementary Tables S18-S19 builders.
Experiment folders under `../../experiments/` provide the runnable bash drivers
and document the required inputs.
