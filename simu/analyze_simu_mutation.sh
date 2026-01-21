#!/usr/bin/env bash
set -euo pipefail

# =========================================
# analyze_simu_mutation.sh
#
# Adapted for new layout:
#   Input:  sample/{EXPR_ID}/{POP}/seed{N}/...
#   Output: results/mutations/{EXPR_ID}/{POP}/{Gene}/
#   Log:    logs/{EXPR_ID}/...
#
# It iterates over provided POPS and finds all seeds within.
# =========================================

# ---------------- Defaults ----------------
POPS="AFR EUR AMR EAS SAS"  # Default to 5 major populations
EXPR_ID="expr_0"
SAMPLE_ROOT="samples"
RESULTS_BASE="results"

# Reference file paths (adjust if needed)
REF_V="ref/IMGT_TRBV_pro.tsv"
REF_J="ref/IMGT_TRBJ_pro.tsv"

# Python script
PY_SCRIPT="collect_mutationV6.py"

# ---------------- Arg parsing ----------------
usage() {
  cat <<EOF
Usage:
  bash $0 --pop "AFR EUR" --expr expr_0

Options:
  --pop           Population list (space separated) (default: "${POPS}")
  --expr          Experiment id folder (default: ${EXPR_ID})
  --sample-root   Root folder containing sample/ (default: ${SAMPLE_ROOT})
  --results-root  Results root folder (default: ${RESULTS_BASE})
  --ref-v         V reference tsv (default: ${REF_V})
  --ref-j         J reference tsv (default: ${REF_J})
  --py            Python script (default: ${PY_SCRIPT})
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pop) POPS="$2"; shift 2;;
    --expr) EXPR_ID="$2"; shift 2;;
    --sample-root) SAMPLE_ROOT="$2"; shift 2;;
    --results-root) RESULTS_BASE="$2"; shift 2;;
    --ref-v) REF_V="$2"; shift 2;;
    --ref-j) REF_J="$2"; shift 2;;
    --py) PY_SCRIPT="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown arg: $1"; usage; exit 1;;
  esac
done

# ---------------- Log setup ----------------
# Modified: Log directory now groups by EXPR_ID
LOG_DIR="logs/${EXPR_ID}"
mkdir -p "$LOG_DIR"

TIMESTAMP=$(date "+%Y%m%d_%H%M%S")
# Just use the first pop in filename for brevity, or generic name
LOG_FILE="${LOG_DIR}/analyze_mutation_${TIMESTAMP}.log"
exec > "$LOG_FILE" 2>&1

echo "======================================================="
echo " Mutation Analysis Task Started"
echo " POPS=${POPS}"
echo " EXPR_ID=${EXPR_ID}"
echo " INPUT_ROOT=${SAMPLE_ROOT}/${EXPR_ID}"
echo " RESULTS_ROOT=${RESULTS_BASE}/mutations/${EXPR_ID}"
echo " LOG=${LOG_FILE}"
echo "======================================================="

# ---------------- Helpers ----------------
check_report_integrity() {
  local report_file="$1"
  local gene_name="$2"
  local target_pattern="\s*- Mismatches/Errors \(excluded from output\): 0.*"

  if [[ -f "$report_file" ]]; then
    if grep -q -E "$target_pattern" "$report_file"; then
      echo "    [PASS] ${gene_name} Report Integrity: Mismatches/Errors is 0."
    else
      echo "    [WARN] ${gene_name} Report Integrity: Mismatches/Errors > 0."
    fi
  else
    echo "    [ERROR] ${gene_name} Report file not found: ${report_file}"
  fi
}

find_clones_file() {
  local dir="$1"
  # Prefer TRB-specific
  local f
  f=$(find "$dir" -maxdepth 2 -type f -name "*clones_TRB.tsv" 2>/dev/null | head -n 1 || true)
  if [[ -n "${f}" ]]; then echo "${f}"; return 0; fi

  # Next: standard MiXCR analyze/export clones name
  f=$(find "$dir" -maxdepth 2 -type f -name "*.clones.tsv" 2>/dev/null | head -n 1 || true)
  if [[ -n "${f}" ]]; then echo "${f}"; return 0; fi

  # Fallback: anything containing clones
  f=$(find "$dir" -maxdepth 2 -type f -name "*clones*.tsv" 2>/dev/null | head -n 1 || true)
  if [[ -n "${f}" ]]; then echo "${f}"; return 0; fi

  echo ""
  return 1
}

# ---------------- Sanity checks ----------------
if [[ ! -f "$REF_V" ]]; then
  echo "[ERROR] REF_V not found: ${REF_V}"; exit 1
fi
if [[ ! -f "$REF_J" ]]; then
  echo "[ERROR] REF_J not found: ${REF_J}"; exit 1
fi
if [[ ! -f "$PY_SCRIPT" ]]; then
  echo "[ERROR] Python script not found: ${PY_SCRIPT}"; exit 1
fi

# ---------------- Main Loop: Iterate POPS ----------------
echo "============================================================"
echo " Execution Log - Started at: $(date)"
echo "============================================================"

# Iterate over each population in the list
for pop in ${POPS//,/ }; do
  
  # Modified: Input path is now sample/expr/pop
  POP_INPUT_DIR="${SAMPLE_ROOT}/${EXPR_ID}/${pop}"
  
  # Modified: Output path is results/mutations/expr/pop
  POP_RESULT_DIR="${RESULTS_BASE}/mutations/${EXPR_ID}/${pop}"

  if [[ ! -d "$POP_INPUT_DIR" ]]; then
    echo "[WARN] Directory for POP=${pop} not found at ${POP_INPUT_DIR}. Skipping."
    continue
  fi

  echo ">>> Processing Population: ${pop}"
  echo "    Input:  ${POP_INPUT_DIR}"
  echo "    Output: ${POP_RESULT_DIR}"

  # Create output directories for V and J
  mkdir -p "${POP_RESULT_DIR}/V"
  mkdir -p "${POP_RESULT_DIR}/J"

  shopt -s nullglob
  # Iterate seeds inside the pop folder
  for seed_path in "${POP_INPUT_DIR}"/seed*; do
    [[ -d "$seed_path" ]] || continue

    seed_base=$(basename "$seed_path")       # e.g. seed0
    seed_num="${seed_base#seed}"             # e.g. 0
    if ! [[ "$seed_num" =~ ^[0-9]+$ ]]; then
      continue
    fi
    
    # In the new structure, seed_path IS the container for data files
    # (Previously it was seed/expr, now it is just seed folder containing files)
    input_dir="$seed_path"
    
    input_file="$(find_clones_file "$input_dir" || true)"
    if [[ -z "$input_file" ]]; then
      echo "  [SKIP] ${pop} - ${seed_base}: no clones tsv found"
      continue
    fi

    echo "  ----------------------------------------------------------"
    echo "  [RUN] seed=${seed_num}"
    echo "    File: ${input_file}"

    # Sample prefix for filenames: POP_seedX
    sample_tag="${pop}_seed${seed_num}"

    # ---- V gene ----
    # Output directly to mutations/expr/pop/V/
    V_out_prefix="${POP_RESULT_DIR}/V/${sample_tag}.V"
    V_report="${V_out_prefix}_report.txt"
    
    # Run Python
    python "$PY_SCRIPT" \
      --input "$input_file" \
      --gene V \
      --ref "$REF_V" \
      --prefix "$V_out_prefix" >/dev/null 2>&1 || echo "    [FAIL] Python script failed for V gene."

    check_report_integrity "$V_report" "V"

    # ---- J gene ----
    # Output directly to mutations/expr/pop/J/
    J_out_prefix="${POP_RESULT_DIR}/J/${sample_tag}.J"
    J_report="${J_out_prefix}_report.txt"

    python "$PY_SCRIPT" \
      --input "$input_file" \
      --gene J \
      --ref "$REF_J" \
      --prefix "$J_out_prefix" >/dev/null 2>&1 || echo "    [FAIL] Python script failed for J gene."

    check_report_integrity "$J_report" "J"

  done
  shopt -u nullglob
done

echo "============================================================"
echo " Analysis Completed at: $(date)"
echo " All results in: ${RESULTS_BASE}/mutations/${EXPR_ID}"
echo " Log: ${LOG_FILE}"
echo "============================================================"