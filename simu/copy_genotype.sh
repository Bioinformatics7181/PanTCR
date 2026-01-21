#!/usr/bin/env bash
set -euo pipefail

# ======================================================
# copy_genotype.sh (Adapted for new layout)
#
# Input:  sample/{EXPR_ID}/{POP}/seed{N}/genotype*.csv
# Output: results/labels/{EXPR_ID}/{POP}/genotype*.csv
# ======================================================

# ---------------- Defaults ----------------
POPS="AFR EUR AMR EAS SAS"
SAMPLE_ROOT="samples"
RESULTS_ROOT="results"
EXPR_ID="expr_0"

# ---------------- Arg parsing -------------
usage() {
  cat <<EOF
Usage:
  bash $0 --expr expr_0 --pop "AFR EUR"

Options:
  --expr          Experiment ID (required)
  --pop           Population list (space separated, default: all 5)
  --sample-root   Root folder containing sample/ (default: ${SAMPLE_ROOT})
  --results-root  Results root folder (default: ${RESULTS_ROOT})
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --expr) EXPR_ID="$2"; shift 2 ;;
    --pop) POPS="$2"; shift 2 ;;
    --sample-root) SAMPLE_ROOT="$2"; shift 2 ;;
    --results-root) RESULTS_ROOT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "${EXPR_ID}" ]]; then
  echo "[ERROR] --expr is required."
  exit 1
fi

# ---------------- Log setup ---------------
LOG_DIR="logs/${EXPR_ID}"
mkdir -p "${LOG_DIR}"

TIMESTAMP=$(date "+%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/copy_genotype_${TIMESTAMP}.log"

# Redirect stdout and stderr to log file AND console
exec > "$LOG_FILE" 2>&1

echo "======================================================="
echo " Batch Copy Genotype"
echo " POPS:          ${POPS}"
echo " EXPR_ID:       ${EXPR_ID}"
echo " INPUT:         ${SAMPLE_ROOT}/${EXPR_ID}/{POP}/..."
echo " OUTPUT:        ${RESULTS_ROOT}/labels/${EXPR_ID}/{POP}/..."
echo " LOG:           ${LOG_FILE}"
echo "======================================================="

# ---------------- Main --------------------
shopt -s nullglob

# Iterate over populations
for pop in ${POPS//,/ }; do
  echo ">>> Processing Population: ${pop} <<<"

  # Modified: Input path structure
  SRC_POP_DIR="${SAMPLE_ROOT}/${EXPR_ID}/${pop}"
  
  # Modified: Output path structure
  TARGET_DEST="${RESULTS_ROOT}/labels/${EXPR_ID}/${pop}"
  
  if [[ ! -d "$SRC_POP_DIR" ]]; then
    echo "  [WARN] Source directory not found: ${SRC_POP_DIR}. Skipping."
    continue
  fi

  mkdir -p "$TARGET_DEST"

  count=0
  skipped=0
  missing=0

  # Iterate over seed directories
  # Use array to safely handle paths
  seed_dirs=( "${SRC_POP_DIR}"/seed* )

  if [[ ${#seed_dirs[@]} -eq 0 ]]; then
      echo "  [WARN] No seed directories found in $SRC_POP_DIR"
      continue
  fi

  for s_dir in "${seed_dirs[@]}"; do
      seed_name=$(basename "$s_dir") # e.g., seed0
      
      # Look for genotype file inside the seed folder
      # Since run_simu script saves it as genotype_{POP}...csv inside seed folder
      files=( "${s_dir}/genotype"*.csv )

      if [[ ${#files[@]} -gt 0 ]]; then
          target_file="${files[0]}"
          filename=$(basename "$target_file")

          # Check collision
          if [[ -f "${TARGET_DEST}/${filename}" ]]; then
              # Optional: Check if content differs? For now just skip.
              # echo "  [SKIP] ${seed_name}: ${filename} exists."
              skipped=$((skipped + 1))
          else
              cp "$target_file" "$TARGET_DEST/"
              count=$((count + 1))
          fi
      else
          echo "  [WARN] ${seed_name}: No genotype CSV found."
          missing=$((missing + 1))
      fi
  done
  
  echo "  [SUMMARY] Copied: $count | Skipped: $skipped | Missing: $missing"
  echo "------------------------------------------------------------"
done

shopt -u nullglob
echo "======================================================="
echo " Task Finished."
echo " Log: ${LOG_FILE}"
echo "======================================================="