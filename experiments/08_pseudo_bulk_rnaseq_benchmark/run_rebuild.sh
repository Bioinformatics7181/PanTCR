#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PYTHON="${PYTHON:-python}"
NODE="${NODE:-node}"

REF_ROOT="${REF_ROOT:-${PACKAGE_ROOT}/data/ref}"

BENCHMARK_ROOT="${PACKAGE_ROOT}/experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines"
SEMI_EXPERIMENT_DIR="${PACKAGE_ROOT}/experiments/07_semi_synthetic_airr_benchmark"
SEMI_RESULTS_ROOT="${SEMI_RESULTS_ROOT:-${BENCHMARK_ROOT}/02_semi_synthetic_airr/generated/results}"
SEMI_WORKSPACE="${SEMI_WORKSPACE:-${SEMI_EXPERIMENT_DIR}/generated/workspace}"
SEMI_EVIDENCE_DIR="${SEMI_EVIDENCE_DIR:-${SEMI_EXPERIMENT_DIR}/generated/evidence_analysis}"
SEMI_TRUTH_CSV="${SEMI_TRUTH_CSV:-${SEMI_EVIDENCE_DIR}/per_truth_call_status.csv}"

SC_BASE="${SC_BASE:-${BENCHMARK_ROOT}/03_pseudo_bulk_rnaseq/generated/per_dataset_results}"
SC_EVIDENCE_DIR="${SC_EVIDENCE_DIR:-${SCRIPT_DIR}/generated/evidence_analysis}"
SC_TRUTH_CSV="${SC_TRUTH_CSV:-${SC_EVIDENCE_DIR}/per_truth_call_status.csv}"

FIG4_AUDIT_DIR="${FIG4_AUDIT_DIR:-${SCRIPT_DIR}/generated/fig4e_nonreference_variant_audit}"
FINAL_AUDIT_DIR="${FINAL_AUDIT_DIR:-${SCRIPT_DIR}/generated/final_two_table_audit}"
REPORT_DIR="${REPORT_DIR:-${SCRIPT_DIR}/generated/report}"
MIN_NAIVE="${MIN_NAIVE:-2}"

if [[ ! -f "${SEMI_TRUTH_CSV}" ]]; then
  if [[ ! -d "${SEMI_RESULTS_ROOT}" ]]; then
    cat >&2 <<EOF
[ERROR] Missing semi-synthetic evidence input.
Expected generated truth-status CSV:
  ${SEMI_TRUTH_CSV}
or expanded benchmark output:
  ${SEMI_RESULTS_ROOT}

Run the semi-synthetic benchmark pipeline under
experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/02_semi_synthetic_airr.
EOF
    exit 2
  fi
  SEMI_CMD=("${PYTHON}" "${PACKAGE_ROOT}/code/experiments/07_semi_synthetic_airr_benchmark/semi_simu_fine_evidence_analysis.py")
  SEMI_CMD+=(--results-root "${SEMI_RESULTS_ROOT}")
  SEMI_CMD+=(--ref-root "${REF_ROOT}")
  SEMI_CMD+=(--workspace "${SEMI_WORKSPACE}")
  SEMI_CMD+=(--out-dir "${SEMI_EVIDENCE_DIR}")
  SEMI_CMD+=(--min-naive "${MIN_NAIVE}")
  "${SEMI_CMD[@]}"
fi

if [[ ! -d "${SC_BASE}" ]]; then
  cat >&2 <<EOF
[ERROR] Missing pseudo-bulk generated input.
Expected package-local generated per-dataset results:
  ${SC_BASE}

Run the pseudo-bulk benchmark driver first:
  bash ${BENCHMARK_ROOT}/03_pseudo_bulk_rnaseq/run_all_datasets.sh
EOF
  exit 2
fi

"${PYTHON}" "${PACKAGE_ROOT}/code/experiments/08_pseudo_bulk_rnaseq_benchmark/scbulk_fine_evidence_analysis.py" \
  --base "${SC_BASE}" \
  --ref-root "${REF_ROOT}" \
  --out-dir "${SC_EVIDENCE_DIR}" \
  --min-naive "${MIN_NAIVE}"

"${PYTHON}" "${PACKAGE_ROOT}/code/experiments/08_pseudo_bulk_rnaseq_benchmark/audit_fig4e_nonreference_variants.py" \
  --truth-csv "${SC_TRUTH_CSV}" \
  --semi-truth-csv "${SEMI_TRUTH_CSV}" \
  --ref-root "${REF_ROOT}" \
  --out-dir "${FIG4_AUDIT_DIR}" \
  --min-naive "${MIN_NAIVE}"

"${PYTHON}" "${PACKAGE_ROOT}/code/experiments/08_pseudo_bulk_rnaseq_benchmark/summarize_figure4_observed_splits.py" \
  --semi-truth-csv "${SEMI_TRUTH_CSV}" \
  --semi-ref-root "${REF_ROOT}" \
  --sc-truth-csv "${SC_TRUTH_CSV}" \
  --sc-ref-root "${REF_ROOT}" \
  --examples-csv "${FIG4_AUDIT_DIR}/fig4e_recovered_nonreference_alleles_for_table.csv" \
  --out-dir "${FINAL_AUDIT_DIR}"

"${NODE}" "${PACKAGE_ROOT}/code/experiments/08_pseudo_bulk_rnaseq_benchmark/summarize_figure4_evidence_s18_s19.mjs" \
  --input-dir "${FINAL_AUDIT_DIR}" \
  --out-dir "${REPORT_DIR}" \
  "$@"
