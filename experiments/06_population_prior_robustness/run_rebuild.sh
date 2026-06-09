#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PYTHON="${PYTHON:-python}"
NODE="${NODE:-node}"

RESULTS_ROOT="${RESULTS_ROOT:-${PACKAGE_ROOT}/experiments/01_in_silico_trbv_benchmarks/results}"
EXPRS="${EXPRS:-expr_ScenarioA expr_ScenarioB}"
EXPR_LABELS="${EXPR_LABELS:-expr_ScenarioA:Scenario_A expr_ScenarioB:Scenario_B}"
GENES="${GENES:-V}"
POPS="${POPS:-AFR}"
FOLDS="${FOLDS:-5}"
PENALTY_K="${PENALTY_K:-0}"
PI_MIN="${PI_MIN:-0.1}"
SKIP_INFERENCE="${SKIP_INFERENCE:-0}"
SKIP_EVAL="${SKIP_EVAL:-0}"

METHODS="BayesNoPrior_no_prior Bayes_priorMismatch Bayes_priorGlobal Bayes_priorMatched"
SUMMARY_DIR="${SUMMARY_DIR:-${SCRIPT_DIR}/generated/summary}"
EVAL_SUMMARY_DIR="${EVAL_SUMMARY_DIR:-${SUMMARY_DIR}/eval_summary}"
SOURCE_CSV="${SOURCE_CSV:-${SUMMARY_DIR}/population_prior_robustness_s14_source.csv}"
REPORT_DIR="${REPORT_DIR:-${SCRIPT_DIR}/generated/report}"

if [[ "${SKIP_INFERENCE}" -eq 0 ]]; then
  for expr in ${EXPRS}; do
    "${PYTHON}" "${PACKAGE_ROOT}/code/experiments/06_population_prior_robustness/run_prior_robustness.py" \
      --results-root "${RESULTS_ROOT}" \
      --expr "${expr}" \
      --genes "${GENES}" \
      --folds "${FOLDS}" \
      --penalty-K "${PENALTY_K}" \
      --pi-min "${PI_MIN}" \
      --modes "matched mismatch global no_prior" \
      --python "${PYTHON}"
  done
fi

if [[ "${SKIP_EVAL}" -eq 0 ]]; then
  for expr in ${EXPRS}; do
    bash "${PACKAGE_ROOT}/experiments/01_in_silico_trbv_benchmarks/eval_alleles.sh" \
      --expr "${expr}" \
      --pop "${POPS}" \
      --gene "${GENES}" \
      --method "${METHODS}" \
      --results-root "${RESULTS_ROOT}" \
      --python "${PYTHON}"
  done
fi

"${PYTHON}" "${PACKAGE_ROOT}/code/experiments/00_benchmark_utils/summarize_eval_results.py" \
  --results-root "${RESULTS_ROOT}" \
  --exprs "${EXPRS}" \
  --methods "${METHODS}" \
  --genes "${GENES}" \
  --populations "${POPS}" \
  --out-dir "${EVAL_SUMMARY_DIR}"

"${PYTHON}" "${PACKAGE_ROOT}/code/experiments/06_population_prior_robustness/build_prior_robustness_source.py" \
  --eval-summary-dir "${EVAL_SUMMARY_DIR}" \
  --out-csv "${SOURCE_CSV}" \
  --expr-labels "${EXPR_LABELS}" \
  --gene "V" \
  --population "AFR"

"${NODE}" "${PACKAGE_ROOT}/code/experiments/06_population_prior_robustness/summarize_prior_robustness_s14.mjs" \
  --source-csv "${SOURCE_CSV}" \
  --out-dir "${REPORT_DIR}" \
  "$@"
