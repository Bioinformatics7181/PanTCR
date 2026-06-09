#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

EXPRS="${EXPRS:-expr_AIRR1 expr_AIRR2}"
TARGET_STEPS="${TARGET_STEPS:-1 2 3 4 5 6 7 8 9 10}"
POPS="${POPS:-MIX}"

EXTRA_ARGS=()

usage() {
  cat <<EOF
Usage:
  bash run_all_airr_benchmarks.sh [options] [extra run_pipeline.sh options]

Options:
  --exprs STR        Space-separated expression IDs. Default: "${EXPRS}"
  --steps STR        Space-separated pipeline steps. Default: "${TARGET_STEPS}"
  --pops STR         Populations. Default: "${POPS}"
  -h, --help         Show help.

All remaining options are forwarded to run_pipeline.sh for each expression ID.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --exprs) EXPRS="$2"; shift 2 ;;
    --steps) TARGET_STEPS="$2"; shift 2 ;;
    --pops) POPS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) EXTRA_ARGS+=("$1"); shift 1 ;;
  esac
done

for expr in ${EXPRS}; do
  echo "[RUN] Semi-synthetic AIRR benchmark ${expr}"
  bash "${SCRIPT_DIR}/run_pipeline.sh" \
    --expr-id "${expr}" \
    --steps "${TARGET_STEPS}" \
    --pops "${POPS}" \
    "${EXTRA_ARGS[@]}"
done

echo "[DONE] Semi-synthetic AIRR benchmark set completed: ${EXPRS}"
