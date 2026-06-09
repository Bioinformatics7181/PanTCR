#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PYTHON="${PYTHON:-python}"
NODE="${NODE:-node}"

RESULTS_ROOT="${RESULTS_ROOT:-${PACKAGE_ROOT}/experiments/01_in_silico_trbv_benchmarks/results}"
EXPR="${EXPR:-expr_ScenarioA}"
GENES="${GENES:-V J}"
POPS="${POPS:-AFR EUR AMR EAS SAS}"
FOLDS="${FOLDS:-5}"
THRESHOLDS="${THRESHOLDS:-0 1 2 3}"
PENALTY_K="${PENALTY_K:-0}"
PI_MIN="${PI_MIN:-0.1}"
SKIP_INFERENCE="${SKIP_INFERENCE:-0}"
SKIP_EVAL="${SKIP_EVAL:-0}"
RESUME="${RESUME:-1}"

SUMMARY_DIR="${SUMMARY_DIR:-${SCRIPT_DIR}/generated/summary}"
EVAL_SUMMARY_DIR="${EVAL_SUMMARY_DIR:-${SUMMARY_DIR}/eval_summary}"
REPORT_DIR="${REPORT_DIR:-${SCRIPT_DIR}/generated/report}"

METHODS=""
for threshold in ${THRESHOLDS}; do
  METHODS="${METHODS} Bayes_fixedGraph_inferMinNaive${threshold} BayesNoPrior_inferMinNaive${threshold}"
done
METHODS="${METHODS# }"

if [[ "${SKIP_INFERENCE}" -eq 0 ]]; then
  infer_args=(
    "${PYTHON}" "${PACKAGE_ROOT}/code/experiments/04_naive_diversity_threshold_sensitivity/run_inference_only_threshold_sensitivity.py"
    --results-root "${RESULTS_ROOT}"
    --expr "${EXPR}"
    --genes "${GENES}"
    --folds "${FOLDS}"
    --thresholds "${THRESHOLDS}"
    --penalty-K "${PENALTY_K}"
    --pi-min "${PI_MIN}"
    --python "${PYTHON}"
  )
  if [[ "${RESUME}" -eq 1 ]]; then
    infer_args+=(--resume)
  fi
  "${infer_args[@]}"
fi

if [[ "${SKIP_EVAL}" -eq 0 ]]; then
  bash "${PACKAGE_ROOT}/experiments/01_in_silico_trbv_benchmarks/eval_alleles.sh" \
    --expr "${EXPR}" \
    --pop "${POPS}" \
    --gene "${GENES}" \
    --method "${METHODS}" \
    --results-root "${RESULTS_ROOT}" \
    --python "${PYTHON}"
fi

"${PYTHON}" "${PACKAGE_ROOT}/code/experiments/00_benchmark_utils/summarize_eval_results.py" \
  --results-root "${RESULTS_ROOT}" \
  --exprs "${EXPR}" \
  --methods "${METHODS}" \
  --genes "${GENES}" \
  --populations "${POPS}" \
  --out-dir "${EVAL_SUMMARY_DIR}"

"${PYTHON}" "${PACKAGE_ROOT}/code/experiments/04_naive_diversity_threshold_sensitivity/build_threshold_sensitivity_sources.py" \
  --eval-summary-dir "${EVAL_SUMMARY_DIR}" \
  --out-dir "${SUMMARY_DIR}"

"${NODE}" "${PACKAGE_ROOT}/code/experiments/04_naive_diversity_threshold_sensitivity/summarize_threshold_sensitivity_s15.mjs" \
  --summary-dir "${SUMMARY_DIR}" \
  --out-dir "${REPORT_DIR}" \
  "$@"
