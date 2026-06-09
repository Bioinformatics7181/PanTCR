# 02 Prediction and Paralogy Audits

Purpose: support prediction-level classification matrices, observed-region compatibility summaries, and gene-label/paralogy audits for Supplementary Tables S2-S10.

Implementation scripts:

- `code/experiments/02_prediction_and_paralogy_audits/make_count_matrix_and_coverage.py`: builds count-based recovery, prediction, and coverage matrices.
- `code/experiments/02_prediction_and_paralogy_audits/assemble_full_provenance_inputs.py`: assembles the normalized per-truth/per-prediction record tree from generated upstream experiment outputs.
- `code/experiments/02_prediction_and_paralogy_audits/build_gene_level_recovery_matrices.py`: builds gene-aware performance and truth-recovery diagnostic matrices from the normalized benchmark records.
- `code/experiments/02_prediction_and_paralogy_audits/summarize_observed_region_recovery_s5_s6.py`: summarizes the S5/S6 observed-region three-class result values from normalized benchmark records.
- `code/experiments/02_prediction_and_paralogy_audits/rebuild_sequence_only_primary_tables.py`: rebuilds the final overall S2/S5/S7 tables from normalized per-truth and per-prediction records. Overall S2 and S5 use one-to-one exact sequence matching across the evaluated mature TRBV region while ignoring the upstream gene label. Per-gene matrices and prediction/gene-label audit tables remain gene-aware.

Primary metric entry point: `rebuild_sequence_only_primary_tables.py` is the final script for the overall S2/S5/S7 primary outputs under the sequence-only matching rule. The source-record builders remain in this directory because they generate coverage strata, observed-region compatibility inputs, per-gene diagnostics, and provenance records.

Inputs: `run_rebuild.sh` consumes generated in silico outputs from `experiments/01_in_silico_trbv_benchmarks/results` together with references in `data/ref`. It then assembles normalized per-truth and per-prediction records from generated in silico, semi-synthetic, and pseudo-bulk benchmark outputs into `generated/full_provenance_inputs/`.

The assembled in silico record root uses manuscript-facing subdirectories under
`01_count_matrix_and_coverage_strata/`: `scenario_a`, `scenario_b`,
`scenario_c`, and `full_length`.

Final S5/S6 reported values use MiXCR-derived `ObservedRanges` before the PanTCR `NaiveDiversityIndex` inference filter. The category "Observed-region Compatible Partial Recovery" means that a non-exact prediction is concordant with the truth at every MiXCR-derived observed position for that sample and gene. The prefix/suffix-compatible partial status produced by `make_count_matrix_and_coverage.py` is a diagnostic category; the final S5/S6 definition is observed-region compatible recovery.

Final overall S2/S5 values are sequence-recovery metrics. Gene-label ambiguity is audited separately in S7/S8 and S9/S10, so sequence-correct but gene-label-discordant predictions are reported through the gene-label audit tables. Per-gene rows retain gene grouping as gene-aware diagnostics.

For pseudo-bulk rows, the final normalized records use the generated pseudo-bulk PanTCR result root, matching the Figure 4E/S18/S19 package chain.

## Run

After experiment 01 has generated `experiments/01_in_silico_trbv_benchmarks/results`
and the semi-synthetic and pseudo-bulk benchmark pipelines have generated their
result roots, run:

```bash
bash experiments/02_prediction_and_paralogy_audits/run_rebuild.sh
```

If the experiment 07 or 08 evidence-analysis CSVs are not present, the driver
generates them from the corresponding benchmark-pipeline outputs before
assembling `generated/full_provenance_inputs/`.

To use a separately organized provenance input tree, provide it explicitly:

```bash
FULL_PROVENANCE_ROOT=experiments/02_prediction_and_paralogy_audits/generated/full_provenance_inputs \
  bash experiments/02_prediction_and_paralogy_audits/run_rebuild.sh
```

`PANTCR_FULL_PROVENANCE_ROOT` is accepted as an equivalent environment variable
for shared execution environments. To rebuild only the in silico count/coverage
strata, run with `SKIP_FULL_PROVENANCE=1`.

## Outputs And Downstream Use

Rerun outputs are written under `generated/`. `run_rebuild.sh` first builds
count/coverage strata from the 01 generated result tree, assembles
`generated/full_provenance_inputs/`, rebuilds gene-level diagnostic matrices,
rebuilds the S5/S6 observed-region matrices, and then rebuilds the final
sequence-only S2/S5/S7 primary tables.
