#!/usr/bin/env bash
set -euo pipefail

# Manuscript/stress synthetic scenario dispatcher.
# This script calls the single-scenario driver in this directory and writes the
# `results/` tree consumed by downstream experiment folders.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

SCENARIOS="scenario_a scenario_b full_length read_depth alpha fragmentation scenario_c"
SELECTED="${SCENARIOS}"
STEPS="1 2 3 4 5 6 7 8 9 10"
POPS="AFR EUR AMR EAS SAS"
GENES="V J"
SAMPLE_ROOT="${SCRIPT_DIR}/samples"
RESULTS_ROOT="${SCRIPT_DIR}/results"
DRY_RUN=0
RUN_MIXCR_ALL=1
MIXCR_ALL_LIBRARY="${PACKAGE_ROOT}/data/ref/hsa_custom_library.json"

usage() {
  cat <<EOF
Usage:
  bash $0 [options]

Options:
  --scenarios STR       Space/comma separated scenario groups to run.
                        Default: ${SCENARIOS}
  --steps STR           Steps passed to simu_validation_pipeline.sh.
                        Default: "${STEPS}"
  --pops STR            Populations to run. Default: "${POPS}"
  --genes STR           Gene segments to process. Default: "${GENES}"
  --sample-root STR     Generated sample root. Default: ${SAMPLE_ROOT}
  --results-root STR    Generated result root. Default: ${RESULTS_ROOT}
  --skip-mixcr-all      Do not run the MiXCR-all/custom-library baseline after
                        Scenario A/B/C sample generation.
  --mixcr-all-library STR
                        MiXCR custom library JSON. Default: data/ref/hsa_custom_library.json
  --dry-run             Print commands without executing them.
  -h, --help            Show help.

Scenario groups:
  scenario_a       expr_ScenarioA baseline fragmented-read Scenario A.
  scenario_b       expr_ScenarioB hidden/default-excluded Scenario B.
  full_length      expr_FullLength high-quality full-length setting.
  read_depth       manuscript RPC grid; Scenario A provides RPC=5 baseline.
  alpha            manuscript clone-abundance alpha grid; Scenario A provides alpha=1.5 baseline.
  fragmentation    manuscript 5' truncation mean/probability stress grids.
  scenario_c       expr_ScenarioC mutation-enabled hidden/novel allele Scenario C.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --scenarios) SELECTED="$2"; shift 2 ;;
    --steps) STEPS="$2"; shift 2 ;;
    --pops) POPS="$2"; shift 2 ;;
    --genes) GENES="$2"; shift 2 ;;
    --sample-root) SAMPLE_ROOT="$2"; shift 2 ;;
    --results-root) RESULTS_ROOT="$2"; shift 2 ;;
    --skip-mixcr-all) RUN_MIXCR_ALL=0; shift 1 ;;
    --mixcr-all-library) MIXCR_ALL_LIBRARY="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown option: $1"; usage; exit 1 ;;
  esac
done

is_abs_path() {
  local p="$1"
  [[ "$p" = /* || "$p" =~ ^[A-Za-z]:[\\/] ]]
}

resolve_package_path() {
  local p="$1"
  if is_abs_path "$p"; then
    echo "$p"
  else
    echo "${PACKAGE_ROOT}/${p}"
  fi
}

SAMPLE_ROOT="$(resolve_package_path "${SAMPLE_ROOT}")"
RESULTS_ROOT="$(resolve_package_path "${RESULTS_ROOT}")"
MIXCR_ALL_LIBRARY="$(resolve_package_path "${MIXCR_ALL_LIBRARY}")"

has_scenario() {
  local wanted="$1"
  for item in ${SELECTED//,/ }; do
    [[ "${item}" == "${wanted}" ]] && return 0
  done
  return 1
}

run_cmd() {
  echo "+ $*"
  if [[ "${DRY_RUN}" -eq 0 ]]; then
    "$@"
  fi
}

run_pipeline() {
  local expr="$1"
  shift
  run_cmd bash "${SCRIPT_DIR}/simu_validation_pipeline.sh" \
    --expr-id "${expr}" \
    --steps "${STEPS}" \
    --pops "${POPS}" \
    --genes "${GENES}" \
    --path-sample "${SAMPLE_ROOT}" \
    --path-results "${RESULTS_ROOT}" \
    "$@"
}

run_mixcr_all() {
  local expr="$1"
  if [[ "${RUN_MIXCR_ALL}" -eq 0 ]]; then
    return 0
  fi
  run_cmd bash "${SCRIPT_DIR}/run_mixcr_all_custom.sh" \
    --expr "${expr}" \
    --pop "${POPS}" \
    --gene "${GENES}" \
    --sample-root "${SAMPLE_ROOT}" \
    --results-root "${RESULTS_ROOT}" \
    --library "${MIXCR_ALL_LIBRARY}"
}

if has_scenario scenario_a; then
  run_pipeline expr_ScenarioA \
    --n 50 --seed-start 0 --exclude-default F \
    --nc 10000 --nr 50000 --alpha 1.5 \
    --p-deg 0.7 --p-deg-sd 0.08 \
    --len-mean 340 --len-sd 6 --len-min 330 --len-max 360 \
    --cut-intact-mean 0 --cut-intact-sd 0 --cut-max-intact 0 \
    --cut-deg-mean 70 --cut-deg-sd 20 --cut-max-deg 140 \
    --c-max 80 --min-keep 250 \
    --art-ss HS25 --art-l 150 --art-c 1 --art-m 320 --art-s 5 \
    --mixcr-preset rna-seq --mixcr-threads 8
  run_mixcr_all expr_ScenarioA
fi

if has_scenario scenario_b; then
  run_pipeline expr_ScenarioB \
    --n 50 --seed-start 0 --exclude-default T \
    --nc 10000 --nr 50000 --alpha 1.5 \
    --p-deg 0.7 --p-deg-sd 0.08 \
    --len-mean 340 --len-sd 6 --len-min 330 --len-max 360 \
    --cut-intact-mean 0 --cut-intact-sd 0 --cut-max-intact 0 \
    --cut-deg-mean 70 --cut-deg-sd 20 --cut-max-deg 140 \
    --c-max 80 --min-keep 250 \
    --art-ss HS25 --art-l 150 --art-c 1 --art-m 320 --art-s 5 \
    --mixcr-preset rna-seq --mixcr-threads 8
  run_mixcr_all expr_ScenarioB
fi

if has_scenario full_length; then
  run_pipeline expr_FullLength \
    --n 50 --seed-start 0 --exclude-default F \
    --nc 10000 --nr 100000 --alpha 1.5 \
    --p-deg 0 --p-deg-sd 0 \
    --len-mean 400 --len-sd 0 --len-min 390 --len-max 410 \
    --cut-intact-mean 0 --cut-intact-sd 0 --cut-max-intact 0 \
    --cut-deg-mean 70 --cut-deg-sd 20 --cut-max-deg 140 \
    --c-max 80 --min-keep 250 \
    --art-ss MSv3 --art-l 250 --art-c 5 --art-m 400 --art-s 5 \
    --mixcr-preset rna-seq --mixcr-threads 8
fi

if has_scenario read_depth; then
  run_pipeline rpc_0p1 --n 50 --seed-start 0 --exclude-default F --nr 1000
  run_pipeline rpc_0p25 --n 50 --seed-start 0 --exclude-default F --nr 2500
  run_pipeline rpc_0p5 --n 50 --seed-start 0 --exclude-default F --nr 5000
  run_pipeline rpc_0p75 --n 50 --seed-start 0 --exclude-default F --nr 7500
  run_pipeline rpc_1p0 --n 50 --seed-start 0 --exclude-default F --nr 10000
  run_pipeline rpc_2p5 --n 50 --seed-start 0 --exclude-default F --nr 25000
  run_pipeline rpc_7p5 --n 50 --seed-start 0 --exclude-default F --nr 75000
  run_pipeline rpc_10p0 --n 50 --seed-start 0 --exclude-default F --nr 100000
fi

if has_scenario alpha; then
  run_pipeline alpha_1p0 --n 50 --seed-start 0 --exclude-default F --alpha 1.0
  run_pipeline alpha_2p0 --n 50 --seed-start 0 --exclude-default F --alpha 2.0
fi

if has_scenario fragmentation; then
  run_pipeline trunc_mean_30 --n 50 --seed-start 0 --exclude-default F --cut-deg-mean 30
  run_pipeline trunc_mean_50 --n 50 --seed-start 0 --exclude-default F --cut-deg-mean 50
  run_pipeline trunc_mean_90 --n 50 --seed-start 0 --exclude-default F --cut-deg-mean 90 --cut-deg-sd 12 --c-max 140 --min-keep 200
  run_pipeline trunc_mean_110 --n 50 --seed-start 0 --exclude-default F --cut-deg-mean 110 --cut-deg-sd 12 --cut-max-deg 170 --c-max 170 --min-keep 170
  run_pipeline trunc_mean_140 --n 50 --seed-start 0 --exclude-default F --cut-deg-mean 140 --cut-deg-sd 12 --cut-max-deg 220 --c-max 220 --min-keep 120
  run_pipeline trunc_mean_160 --n 50 --seed-start 0 --exclude-default F --cut-deg-mean 160 --cut-deg-sd 12 --cut-max-deg 240 --c-max 240 --min-keep 100
  run_pipeline trunc_prob_0p5 --n 50 --seed-start 0 --exclude-default F --p-deg 0.5
  run_pipeline trunc_prob_0p9 --n 50 --seed-start 0 --exclude-default F --p-deg 0.9
  run_pipeline trunc_prob_1p0 --n 50 --seed-start 0 --exclude-default F --p-deg 1.0
fi

if has_scenario scenario_c; then
  run_pipeline expr_ScenarioC \
    --n 50 --seed-start 50 --exclude-default F \
    --mut-prob 1 --j-mut-prob 0 --mut-range-v 70 --mut-range-j 28 \
    --nc 10000 --nr 50000 --alpha 1.5 \
    --p-deg 0.7 --p-deg-sd 0.08 \
    --len-mean 340 --len-sd 6 --len-min 330 --len-max 360 \
    --cut-intact-mean 0 --cut-intact-sd 0 --cut-max-intact 0 \
    --cut-deg-mean 70 --cut-deg-sd 20 --cut-max-deg 140 \
    --c-max 80 --min-keep 250 \
    --art-ss HS25 --art-l 150 --art-c 1 --art-m 320 --art-s 5 \
    --mixcr-preset rna-seq --mixcr-threads 8
  run_mixcr_all expr_ScenarioC
fi

echo "[DONE] Requested scenario groups processed. Generated results root: ${RESULTS_ROOT}"
