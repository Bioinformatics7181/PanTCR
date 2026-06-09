#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PYTHON="${PYTHON:-python}"

RESULTS_ROOT="${RESULTS_ROOT:-${PACKAGE_ROOT}/experiments/01_in_silico_trbv_benchmarks/results}"
REF_ROOT="${REF_ROOT:-${PACKAGE_ROOT}/data/ref}"
COUNT_EXPRS="${COUNT_EXPRS:-expr_ScenarioA expr_ScenarioB expr_ScenarioC expr_FullLength}"
COUNT_METHODS="${COUNT_METHODS:-MiXCR FindAlleles Bayes BayesNoPrior}"
COUNT_OUT_ROOT="${COUNT_OUT_ROOT:-${SCRIPT_DIR}/generated/count_matrix_and_coverage_strata}"
SKIP_COUNT="${SKIP_COUNT:-0}"

FULL_PROVENANCE_ROOT="${FULL_PROVENANCE_ROOT:-${PANTCR_FULL_PROVENANCE_ROOT:-}}"
AUTO_FULL_PROVENANCE_ROOT="${AUTO_FULL_PROVENANCE_ROOT:-${SCRIPT_DIR}/generated/full_provenance_inputs}"
SKIP_FULL_PROVENANCE="${SKIP_FULL_PROVENANCE:-0}"

BENCHMARK_ROOT="${PACKAGE_ROOT}/experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines"
SEMI_RESULTS_ROOT="${SEMI_RESULTS_ROOT:-${BENCHMARK_ROOT}/02_semi_synthetic_airr/generated/results}"
SEMI_EVIDENCE_DIR="${SEMI_EVIDENCE_DIR:-${PACKAGE_ROOT}/experiments/07_semi_synthetic_airr_benchmark/generated/evidence_analysis}"
SEMI_WORKSPACE="${SEMI_WORKSPACE:-${PACKAGE_ROOT}/experiments/07_semi_synthetic_airr_benchmark/generated/workspace}"
SC_BASE="${SC_BASE:-${BENCHMARK_ROOT}/03_pseudo_bulk_rnaseq/generated/per_dataset_results}"
SC_EVIDENCE_DIR="${SC_EVIDENCE_DIR:-${PACKAGE_ROOT}/experiments/08_pseudo_bulk_rnaseq_benchmark/generated/evidence_analysis}"
REAL_MIN_NAIVE="${REAL_MIN_NAIVE:-2}"

OBSERVED_OUT_DIR="${OBSERVED_OUT_DIR:-${SCRIPT_DIR}/generated/observed_region_recovery}"
CLASSIFIED_TRUTH_CSV="${CLASSIFIED_TRUTH_CSV:-${OBSERVED_OUT_DIR}/s5_s6_observed_region_three_class_per_truth.csv}"
SEQUENCE_ONLY_OUT_DIR="${SEQUENCE_ONLY_OUT_DIR:-${SCRIPT_DIR}/generated/sequence_only_primary_metrics}"
SEQUENCE_ONLY_WORKBOOK="${SEQUENCE_ONLY_WORKBOOK:-}"
SEQUENCE_ONLY_OUTPUT_WORKBOOK="${SEQUENCE_ONLY_OUTPUT_WORKBOOK:-${SCRIPT_DIR}/generated/Supplementary Tables sequence-only primary metrics.xlsx}"
STRICT_FULL_PROVENANCE="${STRICT_FULL_PROVENANCE:-0}"

has_required_pair() {
  local dir="$1"
  [[ -f "${dir}/per_truth_call_status.csv" && -f "${dir}/per_prediction_status.csv" ]]
}

ensure_semi_evidence() {
  if has_required_pair "${SEMI_EVIDENCE_DIR}"; then
    return 0
  fi
  if [[ ! -d "${SEMI_RESULTS_ROOT}" ]]; then
    cat >&2 <<EOF
[ERROR] Missing semi-synthetic evidence inputs for experiment 02.
Expected either generated evidence files:
  ${SEMI_EVIDENCE_DIR}/per_truth_call_status.csv
  ${SEMI_EVIDENCE_DIR}/per_prediction_status.csv
or the upstream semi-synthetic benchmark results:
  ${SEMI_RESULTS_ROOT}

Run:
  cd experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/02_semi_synthetic_airr
  bash run_all_airr_benchmarks.sh
EOF
    exit 2
  fi
  "${PYTHON}" "${PACKAGE_ROOT}/code/experiments/07_semi_synthetic_airr_benchmark/semi_simu_fine_evidence_analysis.py" \
    --results-root "${SEMI_RESULTS_ROOT}" \
    --ref-root "${REF_ROOT}" \
    --workspace "${SEMI_WORKSPACE}" \
    --out-dir "${SEMI_EVIDENCE_DIR}" \
    --min-naive "${REAL_MIN_NAIVE}"
}

ensure_scbulk_evidence() {
  if has_required_pair "${SC_EVIDENCE_DIR}"; then
    return 0
  fi
  if [[ ! -d "${SC_BASE}" ]]; then
    cat >&2 <<EOF
[ERROR] Missing pseudo-bulk evidence inputs for experiment 02.
Expected either generated evidence files:
  ${SC_EVIDENCE_DIR}/per_truth_call_status.csv
  ${SC_EVIDENCE_DIR}/per_prediction_status.csv
or the upstream pseudo-bulk benchmark results:
  ${SC_BASE}

Run:
  cd experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/03_pseudo_bulk_rnaseq
  bash run_all_datasets.sh
EOF
    exit 2
  fi
  "${PYTHON}" "${PACKAGE_ROOT}/code/experiments/08_pseudo_bulk_rnaseq_benchmark/scbulk_fine_evidence_analysis.py" \
    --base "${SC_BASE}" \
    --ref-root "${REF_ROOT}" \
    --out-dir "${SC_EVIDENCE_DIR}" \
    --min-naive "${REAL_MIN_NAIVE}"
}

if [[ "${SKIP_COUNT}" -eq 0 ]]; then
  for expr in ${COUNT_EXPRS}; do
    "${PYTHON}" "${PACKAGE_ROOT}/code/experiments/02_prediction_and_paralogy_audits/make_count_matrix_and_coverage.py" \
      --results-root "${RESULTS_ROOT}" \
      --ref-root "${REF_ROOT}" \
      --expr "${expr}" \
      --methods "${COUNT_METHODS}" \
      --out-dir "${COUNT_OUT_ROOT}/${expr}"
  done
fi

if [[ "${SKIP_FULL_PROVENANCE}" -eq 1 ]]; then
  echo "[DONE] Count/coverage strata completed; full-provenance analyses were skipped because SKIP_FULL_PROVENANCE=1."
  exit 0
fi

if [[ -z "${FULL_PROVENANCE_ROOT}" && "${SKIP_FULL_PROVENANCE}" -eq 0 ]]; then
  ensure_semi_evidence
  ensure_scbulk_evidence
  "${PYTHON}" "${PACKAGE_ROOT}/code/experiments/02_prediction_and_paralogy_audits/assemble_full_provenance_inputs.py" \
    --count-root "${COUNT_OUT_ROOT}" \
    --semi-evidence-dir "${SEMI_EVIDENCE_DIR}" \
    --scbulk-evidence-dir "${SC_EVIDENCE_DIR}" \
    --out-root "${AUTO_FULL_PROVENANCE_ROOT}"
  FULL_PROVENANCE_ROOT="${AUTO_FULL_PROVENANCE_ROOT}"
fi

if [[ -n "${FULL_PROVENANCE_ROOT}" && -d "${FULL_PROVENANCE_ROOT}" ]]; then
  "${PYTHON}" "${PACKAGE_ROOT}/code/experiments/02_prediction_and_paralogy_audits/build_gene_level_recovery_matrices.py" \
    --source-root "${FULL_PROVENANCE_ROOT}"

  "${PYTHON}" "${PACKAGE_ROOT}/code/experiments/02_prediction_and_paralogy_audits/summarize_observed_region_recovery_s5_s6.py" \
    --source-root "${FULL_PROVENANCE_ROOT}" \
    --out-dir "${OBSERVED_OUT_DIR}" \
    --write-per-truth \
    "$@"

  if [[ -n "${SEQUENCE_ONLY_WORKBOOK}" ]]; then
    "${PYTHON}" "${PACKAGE_ROOT}/code/experiments/02_prediction_and_paralogy_audits/rebuild_sequence_only_primary_tables.py" \
      --source-root "${FULL_PROVENANCE_ROOT}" \
      --classified-truth "${CLASSIFIED_TRUTH_CSV}" \
      --workbook "${SEQUENCE_ONLY_WORKBOOK}" \
      --output-workbook "${SEQUENCE_ONLY_OUTPUT_WORKBOOK}" \
      --out-dir "${SEQUENCE_ONLY_OUT_DIR}"
  else
    "${PYTHON}" "${PACKAGE_ROOT}/code/experiments/02_prediction_and_paralogy_audits/rebuild_sequence_only_primary_tables.py" \
      --source-root "${FULL_PROVENANCE_ROOT}" \
      --classified-truth "${CLASSIFIED_TRUTH_CSV}" \
      --out-dir "${SEQUENCE_ONLY_OUT_DIR}"
  fi
else
  if [[ -z "${FULL_PROVENANCE_ROOT}" ]]; then
    missing_msg="not provided"
  else
    missing_msg="${FULL_PROVENANCE_ROOT}"
  fi
  cat >&2 <<EOF
[WARN] Full provenance input is not present, so S5/S6 observed-region matrix rebuild was skipped:
  ${missing_msg}

The count/coverage strata above are generated from experiments/01 results. To
rebuild only the count/coverage strata, keep SKIP_FULL_PROVENANCE=1. To use a
separately organized provenance tree, provide FULL_PROVENANCE_ROOT or
PANTCR_FULL_PROVENANCE_ROOT.
EOF
  if [[ "${STRICT_FULL_PROVENANCE}" -eq 1 ]]; then
    exit 2
  fi
fi
