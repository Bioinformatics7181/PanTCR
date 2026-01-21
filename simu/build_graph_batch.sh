#!/usr/bin/env bash
set -euo pipefail

# ======================================================
# build_graph_kfold_train.sh
#
# Logic:
#   1. Iterates through Genes (V, J) and Folds (0..K-1).
#   2. Construct Training Set:
#      For Fold i (Validation), aggregate all CSVs from Fold j where j != i
#      into a TEMP directory.
#   3. Run Python script using the TEMP directory as input.
#   4. Output: results/pang/{EXPR_ID}/{GENE}/fold_{i}
#      (Note: This graph represents the model trained on "All except i")
# ======================================================

# ---------------- Defaults ----------------
EXPR_ID="expr_0"
GENES="V J"
K=5
MIN_NAIVE=1
METADATA="ref/metadata.csv"

RESULTS_ROOT="results"
# Validation root contains the split CSVs
VALIDATION_ROOT="${RESULTS_ROOT}/validation"
# Output root for pangenome graphs
PANG_ROOT="${RESULTS_ROOT}/pang"
PY_SCRIPT="build_pangenome_graph_v7.py"

# ---------------- Arg Parsing ----------------
usage() {
  cat <<EOF
Usage:
  bash $0 --expr expr_0 --k 5 --min-naive 1 --metadata ref/metadata.csv

Options:
  --expr          Experiment ID (required)
  --gene          Gene list (default: "${GENES}")
  --k             Number of folds (default: ${K})
  --min-naive     Min naive count for graph pruning (default: ${MIN_NAIVE})
  --metadata      Path to metadata CSV (default: ${METADATA})
  --results-root  Results root folder (default: ${RESULTS_ROOT})
  --py            Python script path (default: ${PY_SCRIPT})
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --expr) EXPR_ID="$2"; shift 2 ;;
    --gene) GENES="$2"; shift 2 ;;
    --k) K="$2"; shift 2 ;;
    --min-naive) MIN_NAIVE="$2"; shift 2 ;;
    --metadata) METADATA="$2"; shift 2 ;;
    --results-root) RESULTS_ROOT="$2"; shift 2 ;;
    --py) PY_SCRIPT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "${EXPR_ID}" ]]; then
  echo "[ERROR] --expr is required."
  exit 1
fi

if [[ ! -f "${PY_SCRIPT}" ]]; then
  echo "[ERROR] Python script not found: ${PY_SCRIPT}"
  exit 1
fi

if [[ ! -f "${METADATA}" ]]; then
  echo "[WARN] Metadata file not found: ${METADATA}. Script might fail if required."
fi

# ---------------- Log Setup ----------------
LOG_DIR="logs/${EXPR_ID}"
mkdir -p "${LOG_DIR}"

TIMESTAMP=$(date "+%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/build_graph_kfold_${TIMESTAMP}.log"

# Log to file, keep console clean
exec > "$LOG_FILE" 2>&1

echo "======================================================="
echo " Pangenome Graph Construction (K-Fold Train) Started"
echo " EXPR_ID:    ${EXPR_ID}"
echo " K Folds:    ${K}"
echo " MIN_NAIVE:  ${MIN_NAIVE}"
echo " STRATEGY:   Train on [All - Fold_i], Validate on [Fold_i]"
echo " INPUT BASE: ${VALIDATION_ROOT}/${EXPR_ID}/..."
echo " OUTPUT:     ${PANG_ROOT}/${EXPR_ID}/..."
echo " LOG:        ${LOG_FILE}"
echo "======================================================="

# ---------------- Main Execution ----------------

# Define a temporary directory root for aggregating training data
TEMP_ROOT="${PANG_ROOT}/${EXPR_ID}/_temp_staging_${TIMESTAMP}"
mkdir -p "$TEMP_ROOT"

# Iterate over Genes
for gene in ${GENES//,/ }; do
  # [FIX] Ensure variable name matches (loop var 'gene', logic uses 'GENE' -> standardize to uppercase)
  GENE=$(echo "$gene" | tr '[:lower:]' '[:upper:]')
  
  echo ">>> Processing Gene: ${GENE}"

  # Iterate over Folds (0 to K-1). 
  # Here 'i' represents the VALIDATION fold (the one we leave out).
  for ((i=0; i<K; i++)); do
    
    # 1. Setup Output Directory
    OUT_DIR="${PANG_ROOT}/${EXPR_ID}/${GENE}/fold_${i}"
    
    # 2. Setup Temporary Input Directory for this Fold's Training Set
    # We will copy all files EXCEPT fold_${i} into this folder
    TEMP_TRAIN_DIR="${TEMP_ROOT}/${GENE}/train_for_fold_${i}"
    
    # Clean/Create temp dir
    rm -rf "$TEMP_TRAIN_DIR"
    mkdir -p "$TEMP_TRAIN_DIR"
    
    echo "  --------------------------------------------------"
    echo "  [PREP] Building Training Set for Fold ${i} (Excluding Fold ${i})"
    
    count_files=0
    
    # 3. Aggregate Data: Loop through ALL folds 'j'
    for ((j=0; j<K; j++)); do
        # If j == i, this is our validation set, SKIP it for training input
        if [[ "$j" -eq "$i" ]]; then
            continue
        fi
        
        SOURCE_FOLD_DIR="${VALIDATION_ROOT}/${EXPR_ID}/${GENE}/fold_${j}"
        
        if [[ -d "$SOURCE_FOLD_DIR" ]]; then
            # Symlink files to save space/time (using find to avoid errors if empty)
            # We use absolute paths for symlinks to be safe
            abs_source=$(readlink -f "$SOURCE_FOLD_DIR")
            
            # Link all CSVs from source fold to temp train dir
            find "$abs_source" -maxdepth 1 -name "*.csv" -print0 | xargs -0 -I {} ln -s {} "$TEMP_TRAIN_DIR/" 2>/dev/null || true
        fi
    done
    
    # Check if we actually collected any files
    count_files=$(find "$TEMP_TRAIN_DIR" -name "*.csv" | wc -l)
    
    if [[ "$count_files" -eq 0 ]]; then
        echo "  [SKIP] No training files found for Gene ${GENE} Fold ${i} (Check input paths?)"
        continue
    fi

    echo "  [RUN] Construction Graph..."
    echo "    Training Input (Temp): ${TEMP_TRAIN_DIR} (${count_files} files)"
    echo "    Output (Graph):        ${OUT_DIR}"

    # Ensure output directory exists
    mkdir -p "$OUT_DIR"

    # 4. Run Python Script
    # Input is the aggregated TEMP_TRAIN_DIR
    if python "${PY_SCRIPT}" \
      --in_dir "${TEMP_TRAIN_DIR}" \
      --out_dir "${OUT_DIR}" \
      --metadata "${METADATA}" \
      --min_naive "${MIN_NAIVE}"; then
      
      echo "    [DONE] Graph built successfully."
    else
      echo "    [FAIL] Python script failed for Gene ${GENE} / Fold ${i}."
    fi
    
    # Optional: Clean up this specific temp dir immediately to save inode space
    # rm -rf "$TEMP_TRAIN_DIR"

  done
  echo "------------------------------------------------------------"
done

# Final Cleanup of all temp dirs
echo "[CLEANUP] Removing temporary directories..."
rm -rf "$TEMP_ROOT"

echo "======================================================="
echo " Task Finished."
echo " Log: ${LOG_FILE}"
echo "======================================================="