# Experiments

This directory contains the manuscript reproduction workflows. Each numbered
folder corresponds to a paper-facing benchmark, audit, or supplementary result
family. The folders contain bash drivers and README files; reusable Python
implementations live under `../code/`.

## Folder Guide

- `01_in_silico_trbv_benchmarks/`: controlled truth-known synthetic TRBV/TRBJ
  simulations, including Scenario A/B/C, full-length simulation, MiXCR-all, and
  manuscript stress grids. It also contains benchmark pipelines that generate
  semi-synthetic AIRR, pseudo-bulk/scBulk, and public bulk result roots for
  downstream experiments.
- `02_prediction_and_paralogy_audits/`: count, recovery-status,
  observed-region, and gene-label/paralogy audit builders for Supplementary
  Tables S2-S10.
- `03_leave_target_allele_out/`: leave-target-allele-out stress test supporting
  Supplementary Tables S12-S13.
- `04_naive_diversity_threshold_sensitivity/`: inference-only support-threshold
  sensitivity analysis for Supplementary Table S15.
- `05_public_bulk_candidate_nonimgt/`: retained-evidence and graph-imputation
  audit for candidate non-IMGT alleles in public bulk RNA-seq, supporting
  Supplementary Table S20.
- `06_population_prior_robustness/`: no-prior, mismatched-prior, global-prior,
  and matched-prior robustness analysis for Supplementary Table S14.
- `07_semi_synthetic_airr_benchmark/`: downstream summaries and retained
  evidence diagnostics for the semi-synthetic AIRR benchmark, supporting
  Supplementary Table S16.
- `08_pseudo_bulk_rnaseq_benchmark/`: pseudo-bulk/scBulk benchmark summaries and
  observed/evidence-supported non-default allele analyses for Figure 4E and
  Supplementary Tables S17-S19.

## Generated Outputs And Dependencies

Fresh reruns write to `results/`, `runs/`, or `generated/` directories.

Some workflows depend on earlier generated results. For example, experiments
02, 04, and 06 consume `01_in_silico_trbv_benchmarks/results/`, while
experiments 05, 07, and 08 consume benchmark-pipeline outputs generated under
experiment 01. Experiment 02 also assembles the semi-synthetic and pseudo-bulk
per-truth/per-prediction records generated from those benchmark-pipeline
outputs. See each experiment README for the expected input path and rerun
command.

## Recommended Reproduction Order

Run commands from the package root after completing the data preparation steps
in the root README.

1. Generate the controlled in silico benchmark outputs:

```bash
bash experiments/01_in_silico_trbv_benchmarks/run_all_scenarios.sh
```

This creates `experiments/01_in_silico_trbv_benchmarks/results/`, which is the
in silico upstream source for experiments 02, 04, and 06.

2. Generate benchmark-pipeline outputs used by biological/proxy-truth analyses:

```bash
cd experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/02_semi_synthetic_airr
bash run_all_airr_benchmarks.sh

cd ../03_pseudo_bulk_rnaseq
bash run_all_datasets.sh

cd ../04_public_bulk_rnaseq
bash run_sra_to_clones.sh --bulk-label BULK-1 --acc-list input/BULK-1/SRR_Acc_List.txt
bash run_from_clones.sh --bulk-label BULK-1 --clones-dir generated/BULK-1/clones
bash run_sra_to_clones.sh --bulk-label BULK-2 --acc-list input/BULK-2/SRR_Acc_List.txt
bash run_from_clones.sh --bulk-label BULK-2 --clones-dir generated/BULK-2/clones
cd ../../../..
```

The semi-synthetic, pseudo-bulk/scBulk, and public-bulk README files describe
the required raw or proxy-truth inputs for these commands.

3. Rebuild downstream analyses from the generated result roots. Experiment 02
   assembles its normalized full-provenance input tree during this step:

```bash
bash experiments/02_prediction_and_paralogy_audits/run_rebuild.sh
bash experiments/03_leave_target_allele_out/run_leaveout_allele_specific.sh
bash experiments/04_naive_diversity_threshold_sensitivity/run_rebuild.sh
bash experiments/05_public_bulk_candidate_nonimgt/run_rebuild.sh
bash experiments/06_population_prior_robustness/run_rebuild.sh
bash experiments/07_semi_synthetic_airr_benchmark/run_rebuild.sh
bash experiments/08_pseudo_bulk_rnaseq_benchmark/run_rebuild.sh
```
