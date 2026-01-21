#!/usr/bin/env bash
set -euo pipefail

# ======================================================
# split_kfold.sh
#
# Splits mutation CSV files into K folds for cross-validation.
#
# Input:  results/mutations/{EXPR_ID}/{POP}/{GENE}/*.csv
# Output: results/validation/{EXPR_ID}/{GENE}/fold_{0..K-1}/
# ======================================================

# ---------------- Defaults ----------------
K=5
SEED=42
EXPR_ID="expr_0"
RESULTS_ROOT="results"

# ---------------- Arg parsing -------------
usage() {
  cat <<EOF
Usage:
  bash $0 --expr expr_0 --k 5 --seed 42

Options:
  --expr          Experiment ID (required)
  --k             Number of folds (default: ${K})
  --seed          Random seed (default: ${SEED})
  --results-root  Root folder (default: ${RESULTS_ROOT})
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --expr) EXPR_ID="$2"; shift 2 ;;
    --k) K="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
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
LOG_FILE="${LOG_DIR}/split_kfold_${EXPR_ID}_${TIMESTAMP}.log"

exec > "$LOG_FILE" 2>&1

echo "======================================================="
echo " K-Fold Split Task Started"
echo " EXPR_ID: ${EXPR_ID} | K=${K} | SEED=${SEED}"
echo " INPUT:   ${RESULTS_ROOT}/mutations/${EXPR_ID}"
echo " OUTPUT:  ${RESULTS_ROOT}/validation/${EXPR_ID}"
echo "======================================================="

# ---------------- Python Splitter Helper ----------------
# Reads file paths from stdin, shuffles them with SEED, 
# and assigns a fold index (i % k).
PYTHON_SPLITTER=$(cat <<END
import sys
import random

lines = [line.strip() for line in sys.stdin if line.strip()]
k = int(sys.argv[1])
seed = int(sys.argv[2])

random.seed(seed)
random.shuffle(lines)

for i, line in enumerate(lines):
    fold_idx = i % k
    print(f"{line}\t{fold_idx}")
END
)

# ---------------- Main Logic ----------------
INPUT_BASE="${RESULTS_ROOT}/mutations/${EXPR_ID}"
OUTPUT_BASE="${RESULTS_ROOT}/validation/${EXPR_ID}"

if [[ ! -d "$INPUT_BASE" ]]; then
  echo "[ERROR] Input directory not found: ${INPUT_BASE}"
  exit 1
fi

# Process V and J genes separately
for GENE in V J; do
  echo ">>> Processing Gene: ${GENE}"

  # 1. Prepare output directories (fold_0 ... fold_K-1)
  GENE_OUT_DIR="${OUTPUT_BASE}/${GENE}"
  
  # Clean previous runs to ensure valid validation sets
  if [[ -d "$GENE_OUT_DIR" ]]; then
    rm -rf "$GENE_OUT_DIR"
  fi

  for ((i=0; i<K; i++)); do
    mkdir -p "${GENE_OUT_DIR}/fold_${i}"
  done

  # 2. Identify all Population directories under input root
  shopt -s nullglob
  pop_dirs=( "${INPUT_BASE}"/* )
  shopt -u nullglob

  if [[ ${#pop_dirs[@]} -eq 0 ]]; then
    echo "  [ERROR] No population directories found."
    continue
  fi

  total_files_count=0

  # 3. Iterate over each population (Stratified Split)
  # We split EACH population independently to ensure balanced folds.
  for p_dir in "${pop_dirs[@]}"; do
    pop_name=$(basename "$p_dir")
    
    # Path: results/mutations/{EXPR}/{POP}/{GENE}
    gene_input_dir="${p_dir}/${GENE}"
    
    if [[ ! -d "$gene_input_dir" ]]; then
      continue
    fi

    # Collect all CSV files for this POP + GENE
    mapfile -t files < <(find "$gene_input_dir" -maxdepth 1 -name "*.csv" | sort)
    
    num_files=${#files[@]}
    if [[ "$num_files" -eq 0 ]]; then
      echo "    [SKIP] ${pop_name}: No csv files found."
      continue
    fi

    echo "    [SPLIT] ${pop_name}: Distributing ${num_files} files..."

    # 4. Use Python to shuffle and assign folds
    printf "%s\n" "${files[@]}" | python3 -c "$PYTHON_SPLITTER" "$K" "$SEED" | while read -r filepath fold_id; do
      
      dest_dir="${GENE_OUT_DIR}/fold_${fold_id}"
      cp "$filepath" "$dest_dir/"
      
    done
    
    total_files_count=$((total_files_count + num_files))
  done

  echo "  [DONE] ${GENE} Gene finished. Total files: ${total_files_count}"
  echo "------------------------------------------------------------"
done

echo "======================================================="
echo " Task Finished."
echo " Log: ${LOG_FILE}"
echo "======================================================="