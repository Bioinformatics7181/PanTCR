# 06 Population Prior Robustness

Purpose: support Supplementary Table S14 and the matched-population-prior robustness analysis in the synthetic benchmark.

Implementation scripts:

- `code/experiments/06_population_prior_robustness/run_prior_robustness.py`: runs graph-prior inference modes for AFR test samples under no-prior, mismatched, global, and matched prior-weight settings.
- `code/experiments/00_benchmark_utils/summarize_eval_results.py` and `code/experiments/06_population_prior_robustness/build_prior_robustness_source.py`: convert generated evaluation CSVs into the S14 source CSV.
- `code/experiments/06_population_prior_robustness/summarize_prior_robustness_s14.mjs`: writes compact CSV/Markdown summaries for the Supplementary Table S14 source values.

Key parameters: AFR test samples were evaluated in Scenario A (`expr_ScenarioA`) and Scenario B (`expr_ScenarioB`) under no-prior, EUR-mismatched, global, and AFR-matched prior modes. The graph-prior modes excluded evaluated test-sample reads and recorded genotypes. The runner explicitly passes the manuscript-facing synthetic inference settings `min_naive=1`, `penalty_K=0`, and `pi_min=0.1` unless overridden.

Inputs: generated Scenario A/B validation, graph, label, and mutation outputs from `experiments/01_in_silico_trbv_benchmarks/results`. Use `RESULTS_ROOT` when pointing the analysis to another relative rerun directory.

## Run

After experiment 01 has generated Scenario A and Scenario B result trees, run:

```bash
bash experiments/06_population_prior_robustness/run_rebuild.sh
```

Use `RESULTS_ROOT` to point to a separate in silico rerun:

```bash
RESULTS_ROOT=experiments/01_in_silico_trbv_benchmarks/results \
  bash experiments/06_population_prior_robustness/run_rebuild.sh
```

## Outputs And Downstream Use

Rerun summaries are written under `generated/summary` and `generated/report`.
The report step writes compact CSV/Markdown files rather than recreating the
formal workbook layout.
