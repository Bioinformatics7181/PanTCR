# 01 In Silico TRBV Benchmarks

Purpose: reproduce the controlled synthetic TRBV benchmark workflows used for Figure 2 and the synthetic rows in Supplementary Tables S2-S11.

## Drivers In This Folder

- `simu_validation_pipeline.sh`: end-to-end driver for simulation, mutation extraction, k-fold split, graph construction, inference, baselines, and evaluation.
- `run_all_scenarios.sh`: aggregate driver that enumerates the reproduced Scenario A/B, Scenario C, full-length, read-depth/RPC, alpha, and fragmentation/truncation settings by calling `simu_validation_pipeline.sh` with explicit parameters. For Scenario A/B/C it also runs the MiXCR-all/custom-library baseline unless `--skip-mixcr-all` is set.
- `run_simu_batch.sh`: scenario-specific simulation driver.
- `analyze_simu_mutation.sh`: runs PanTCR observation extraction on MiXCR clone outputs.
- `run_mixcr_all_custom.sh`: runs the MiXCR-all expanded-reference baseline on generated samples and evaluates it with `evaluate_allele_calls.py --mode mixcr-all`.
- `copy_genotype.sh`, `split_kfold.sh`, `build_graph_batch.sh`, `infer_bayes.sh`, `infer_mixcr.sh`, `infer_findalleles.sh`, `eval_alleles.sh`: individual pipeline steps.

All Python implementations are centralized:

- simulation: `code/simulation/`
- PanTCR: `code/pantcr/`
- baselines: `code/experiments/00_benchmark_utils/`
- evaluation: `code/experiments/00_benchmark_utils/evaluate_allele_calls.py`

## Key Parameters

The manuscript-facing synthetic benchmark uses simulation-specific inference settings:

- `min_naive=1`
- `penalty_K=0`
- `pi_min=0.1`
- Scenario A (`expr_ScenarioA`) is the baseline fragmented-read simulation:
  `10,000` clones, `50,000` simulated molecules, alpha `1.5`, truncation
  probability `0.7`, mean 5' loss `70`, HiSeq 2500 `2 x 150 bp`, and ART
  sampling/copy setting `1`.
- Scenario B (`expr_ScenarioB`) uses the Scenario A settings with MiXCR-default
  reference alleles excluded during genotype sampling.
- Scenario C (`expr_ScenarioC`) uses the Scenario A settings plus mutation-enabled
  genotype seeding:
  `--seed-start 50 --mut-prob 1 --j-mut-prob 0 --mut-range-v 70 --mut-range-j 28`.
- The high-quality full-length setting (`expr_FullLength`) uses
  `art_ss=MSv3`, `art_l=250`, `art_c=5`, `art_m=400`, and
  `art_s=5`.
- The manuscript RPC grid is `0.1`, `0.25`, `0.5`, `0.75`, `1.0`, `2.5`,
  `5.0`, `7.5`, and `10.0`. In `run_all_scenarios.sh`, Scenario A provides
  the RPC `5.0` baseline and the non-baseline RPC outputs use readable IDs
  such as `rpc_0p1`, `rpc_1p0`, and `rpc_10p0`.
- The clone-abundance alpha grid is `1.0`, `1.5`, and `2.0`. Scenario A
  provides alpha `1.5`; the non-baseline settings use `alpha_1p0` and
  `alpha_2p0`.
- The physical-information-loss grid varies mean 5' truncation length from
  `30` to `160` bp and truncation probability across `0.5`, `0.7`, `0.9`,
  and `1.0`, one axis at a time. Scenario A provides probability `0.7` and
  mean `70`; the non-baseline settings use readable IDs such as
  `trunc_mean_30`, `trunc_mean_160`, and `trunc_prob_1p0`.
- MiXCR-all evaluation is allele-name based and folds only mature-region-equivalent
  `TRBV7-7*01/02` and `TRBV15*01/05`; this rule is implemented in the shared
  evaluator.

The bash drivers pass these synthetic-benchmark inference settings explicitly.

## Data Requirements

Small reference inputs are bundled under `data/ref/`. Full reruns also require
`cdr3_dict.pkl`; unpack the bundled archive from the package root:

```bash
unzip data/ref/cdr3_dict.zip -d data/ref/
```

or provide the extracted file as:

```text
data/ref/cdr3_dict.pkl
```

or set:

```bash
export PANTCR_SIM_REF_ROOT=data/ref
```

Install the external tools used by the workflow, including MiXCR, ART, and
FindAlleles-related MiXCR commands, before running full simulations.

The MiXCR-all baseline additionally uses:

```text
data/ref/hsa_custom_library.json
```

## Run

Run all manuscript-facing in silico scenarios and stress grids from the package
root:

```bash
bash experiments/01_in_silico_trbv_benchmarks/run_all_scenarios.sh
```

To run one scenario group, pass `--scenarios`, for example:

```bash
bash experiments/01_in_silico_trbv_benchmarks/run_all_scenarios.sh \
  --scenarios "scenario_a scenario_c full_length"
```

The aggregate driver calls `simu_validation_pipeline.sh`, which performs
simulation, MiXCR evidence extraction, genotype-label preparation, k-fold graph
construction, PanTCR/MiXCR/FindAlleles inference, and evaluation.

## Outputs And Downstream Use

Fresh reruns write generated samples under `samples/` and evaluated method
outputs under `results/`. The `results/` tree is the upstream source for:

- `experiments/02_prediction_and_paralogy_audits/`
- `experiments/04_naive_diversity_threshold_sensitivity/`
- `experiments/06_population_prior_robustness/`

The benchmark-pipeline subdirectories under `benchmark_pipelines/` generate
semi-synthetic AIRR, pseudo-bulk/scBulk, and public-bulk result roots used by
experiments 05, 07, and 08.

## Scenario Coverage

`run_all_scenarios.sh` provides explicit commands for Scenario A, Scenario B,
Scenario C, high-quality full-length simulation, read-depth/RPC settings, alpha
settings, and fragmentation/truncation settings. Reruns write generated sample
and inference trees to this folder's `samples/` and `results/` directories by default. User-supplied
relative `--sample-root`, `--results-root`, and MiXCR-all library paths are
resolved relative to the package root.

The default aggregate script covers the synthetic settings described in the
current manuscript text, figure legends, supplementary figures, and
supplementary tables. Generated sample and result trees are produced by the
scripts under the paths described above.
