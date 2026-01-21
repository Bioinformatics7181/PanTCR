#!/usr/bin/env bash
set -euo pipefail

# ======================================================
# infer_Bayes.sh (K-Fold + NoPrior Support)
#
# Logic:
#   1. Iterate Genes (V, J) and Folds (0..K-1).
#   2. Input: Validation Set -> results/validation/{EXPR_ID}/{GENE}/fold_{i}
#   3. Prior (Optional): Pangenome Graph -> results/pang/{EXPR_ID}/{GENE}/fold_{i}
#      - If --no-prior is set, Pangenome Graph is IGNORED.
#   4. Process: infer_genotype_bayes_v4.py
#   5. Output: 
#      - With Prior: results/infer/{EXPR_ID}/Bayes/{POP}/{GENE}/
#      - No Prior:   results/infer/{EXPR_ID}/BayesNoPrior/{POP}/{GENE}/
# ======================================================

# ---------------- Defaults ----------------
EXPR_ID="expr_0"
GENES="V J"
K=5
MIN_NAIVE=1
USE_PRIOR=1  # Default: Use Pangenome

RESULTS_ROOT="results"
VALIDATION_ROOT="${RESULTS_ROOT}/validation"
PANG_ROOT="${RESULTS_ROOT}/pang"
INFER_ROOT="${RESULTS_ROOT}/infer"

PY_SCRIPT="infer_genotype_bayes_v4.py"

# ---------------- Arg Parsing ----------------
usage() {
  cat <<EOF
Usage:
  bash $0 --expr expr_0 --k 5 --min-naive 1 [--no-prior]

Options:
  --expr          Experiment ID (required)
  --gene          Gene list (default: "${GENES}")
  --k             Number of folds (default: ${K})
  --min-naive     Min naive count (default: ${MIN_NAIVE})
  --no-prior      If set, do NOT use pangenome graph (save to BayesNoPrior)
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
    --no-prior) USE_PRIOR=0; shift 1 ;;
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

# ---------------- Config Setup ----------------
# Determine Output Folder Name based on Prior Usage
if [[ "$USE_PRIOR" -eq 1 ]]; then
    INFER_SUBDIR="Bayes"
else
    INFER_SUBDIR="BayesNoPrior"
fi

# ---------------- Log Setup ----------------
LOG_DIR="logs/${EXPR_ID}"
mkdir -p "${LOG_DIR}"

TIMESTAMP=$(date "+%Y%m%d_%H%M%S")
# Append suffix to log file to distinguish runs
LOG_FILE="${LOG_DIR}/infer_${INFER_SUBDIR}_${TIMESTAMP}.log"

# STRICT LOGGING: Only to file
exec > "$LOG_FILE" 2>&1

echo "======================================================="
echo " Bayesian Inference Task Started"
echo " EXPR_ID:    ${EXPR_ID}"
echo " K Folds:    ${K}"
echo " MIN_NAIVE:  ${MIN_NAIVE}"
echo " USE PRIOR:  $( [[ $USE_PRIOR -eq 1 ]] && echo "YES" || echo "NO" )"
echo " INPUT:      ${VALIDATION_ROOT}/${EXPR_ID}/..."
if [[ "$USE_PRIOR" -eq 1 ]]; then
  echo " PRIOR MAP:  ${PANG_ROOT}/${EXPR_ID}/..."
fi
echo " OUTPUT:     ${INFER_ROOT}/${EXPR_ID}/${INFER_SUBDIR}/..."
echo " LOG:        ${LOG_FILE}"
echo "======================================================="

# ---------------- Main Execution ----------------

for gene in ${GENES//,/ }; do
  GENE=$(echo "$gene" | tr '[:lower:]' '[:upper:]')
  echo ">>> Processing Gene: ${GENE}"

  # Iterate over Folds (0 to K-1)
  for ((i=0; i<K; i++)); do
    
    # 1. Define Paths for this Fold
    VAL_DIR="${VALIDATION_ROOT}/${EXPR_ID}/${GENE}/fold_${i}"
    GRAPH_DIR="${PANG_ROOT}/${EXPR_ID}/${GENE}/fold_${i}"
    
    echo "  --------------------------------------------------"
    echo "  [FOLD ${i}] Gene: ${GENE}"
    echo "    Samples: ${VAL_DIR}"
    
    # 2. Check Input Directory
    if [[ ! -d "$VAL_DIR" ]]; then
      echo "    [SKIP] Validation directory not found."
      continue
    fi
    
    # 3. Check Graph Directory (Only if USE_PRIOR=1)
    if [[ "$USE_PRIOR" -eq 1 ]]; then
        echo "    Graph:   ${GRAPH_DIR}"
        if [[ ! -d "$GRAPH_DIR" ]]; then
            echo "    [WARN] Graph directory not found. Skipping fold."
            continue
        fi
    else
        echo "    Graph:   [DISABLED by --no-prior]"
    fi

    # 4. Iterate over Samples
    shopt -s nullglob
    input_files=( "${VAL_DIR}"/*.csv )
    shopt -u nullglob

    if [[ ${#input_files[@]} -eq 0 ]]; then
        echo "    [SKIP] No CSV files found in validation fold."
        continue
    fi

    for sample_file in "${input_files[@]}"; do
        filename=$(basename "$sample_file")
        
        # Regex: Expects POP_seedX... or POP_...
        if [[ $filename =~ ^([A-Za-z]+)_(seed[0-9]+) ]]; then
            POP_NAME="${BASH_REMATCH[1]}"
            SEED_STR="${BASH_REMATCH[2]}" # e.g. seed0
            
            # 5. Define Output Directory
            # Structure: results/infer/expr/{Bayes|BayesNoPrior}/POP/Gene
            OUT_DIR="${INFER_ROOT}/${EXPR_ID}/${INFER_SUBDIR}/${POP_NAME}/${GENE}"
            mkdir -p "$OUT_DIR"
            
            # Output Filename
            OUT_FILE="${OUT_DIR}/infer_${POP_NAME}_${SEED_STR}.${GENE}.csv"
            
            # 6. Construct Python Command (Using Array for robustness)
            CMD=(python "${PY_SCRIPT}")
            CMD+=(--sample_csv "${sample_file}")
            CMD+=(--population_id "${POP_NAME}")
            CMD+=(--min_naive "${MIN_NAIVE}")
            CMD+=(--out "${OUT_FILE}")

            # Conditionally add pangenome argument
            if [[ "$USE_PRIOR" -eq 1 ]]; then
                CMD+=(--pangenome_dir "${GRAPH_DIR}")
            fi
            
            # 7. Execute
            if "${CMD[@]}"; then
                echo "    [OK] ${filename} -> Saved"
            else
                echo "    [FAIL] Inference failed for ${filename}"
            fi

        else
            echo "    [WARN] Could not parse filename format: ${filename}. Skipping."
        fi
    done

  done
  echo "------------------------------------------------------------"
done

echo "======================================================="
echo " Task Finished."
echo " Log: ${LOG_FILE}"
echo "======================================================="