#!/usr/bin/env bash
set -euo pipefail

# ======================================================
# eval_alleles.sh (Universal Batch Evaluator)
#
# Logic:
#   1. Iterate Methods, Pops, Genes.
#   2. Match Infer file (seedX) with Genotype file (seedX).
#   3. Run eval_alleles_without_CDR3.py
#
# Directory Mapping:
#   Labels: results/labels/{EXPR_ID}/{POP}/
#   Input:  results/infer/{EXPR_ID}/{METHOD}/{POP}/{GENE}/
#   Output: results/eval/{EXPR_ID}/{METHOD}/{POP}/{GENE}/
# ======================================================

# ---------------- Defaults ----------------
EXPR_ID="expr_0"
POPS="AFR EUR AMR EAS SAS"
GENES="V J"
METHODS="MiXCR FindAlleles Bayes BayesNoPrior" # Default methods to check

RESULTS_ROOT="results"
INSEC="F"  # Intersect mode (T/F)

# Reference Helper Files (Adjust paths if needed)
PMTR_NAME="ref/pmTR_TRB_V_J_trim.csv"
TRBV_INDEX="ref/TRB_index.csv"
PY_SCRIPT="eval_alleles_without_CDR3.py"

# ---------------- Arg Parsing ----------------
usage() {
  cat <<EOF
Usage:
  bash $0 --expr expr_0

Options:
  --expr          Experiment ID (required)
  --pop           Population list (default: "${POPS}")
  --gene          Gene list (default: "${GENES}")
  --method        Method list (folder names in results/infer/) (default: "${METHODS}")
  --results-root  Results root folder (default: ${RESULTS_ROOT})
  --intersect     Enable intersect mode (only evaluate shared alleles)
  --ref-pmtr      Path to pmTR reference (default: ${PMTR_NAME})
  --ref-index     Path to TRB index (default: ${TRBV_INDEX})
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --expr) EXPR_ID="$2"; shift 2 ;;
    --pop) POPS="$2"; shift 2 ;;
    --gene) GENES="$2"; shift 2 ;;
    --method) METHODS="$2"; shift 2 ;;
    --results-root) RESULTS_ROOT="$2"; shift 2 ;;
    --ref-pmtr) PMTR_NAME="$2"; shift 2 ;;
    --ref-index) TRBV_INDEX="$2"; shift 2 ;;
    --intersect) INSEC="T"; shift 1 ;;
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

# ---------------- Log Setup ----------------
LOG_DIR="logs/${EXPR_ID}"
mkdir -p "${LOG_DIR}"
TIMESTAMP=$(date "+%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/eval_alleles_${TIMESTAMP}.log"

# Log to file only
exec > "$LOG_FILE" 2>&1

echo "======================================================="
echo " Batch Evaluation Task Started"
echo " EXPR_ID:    ${EXPR_ID}"
echo " POPS:       ${POPS}"
echo " GENES:      ${GENES}"
echo " METHODS:    ${METHODS}"
echo " INTERSECT:  ${INSEC}"
echo "======================================================="

# ---------------- Main Execution ----------------

# 1. Iterate Methods (MiXCR, Bayes, etc.)
for method in ${METHODS//,/ }; do
  echo ">>> Processing Method: ${method}"

  # 2. Iterate Populations
  for pop in ${POPS//,/ }; do
    
    # Define Label Directory (Truth)
    LABEL_DIR="${RESULTS_ROOT}/labels/${EXPR_ID}/${pop}"

    if [[ ! -d "$LABEL_DIR" ]]; then
       echo "  [SKIP] Label directory not found for ${pop}"
       continue
    fi

    # 3. Iterate Genes
    for gene in ${GENES//,/ }; do
      GENE=$(echo "$gene" | tr '[:lower:]' '[:upper:]')
      
      # Define Inference Directory (Input)
      INFER_DIR="${RESULTS_ROOT}/infer/${EXPR_ID}/${method}/${pop}/${GENE}"
      
      # Define Evaluation Directory (Output)
      # results/eval/{EXPR_ID}/{METHOD}/{POP}/{GENE}
      OUT_DIR="${RESULTS_ROOT}/eval/${EXPR_ID}/${method}/${pop}/${GENE}"

      if [[ ! -d "$INFER_DIR" ]]; then
        # Silent skip if folder doesn't exist (maybe only V was run)
        # echo "  [DEBUG] Missing infer dir: ${INFER_DIR}"
        continue
      fi

      echo "  --------------------------------------------------"
      echo "  [RUN] ${method} | ${pop} | ${GENE}"
      echo "    Infer Dir: ${INFER_DIR}"
      echo "    Eval Dir:  ${OUT_DIR}"

      mkdir -p "$OUT_DIR"
      
      shopt -s nullglob
      infer_files=( "${INFER_DIR}"/*.csv )
      
      count=0
      
      for infer_file in "${infer_files[@]}"; do
          filename=$(basename "$infer_file")
          
          # Match Seed ID from filename (e.g., infer_AFR_seed0.V.csv)
          if [[ $filename =~ seed([0-9]+) ]]; then
              seed_id="${BASH_REMATCH[1]}"
              
              # Find corresponding Label File
              # Look for genotype_*seed{ID}.csv in LABEL_DIR
              label_match=( "${LABEL_DIR}/genotype_"*"seed${seed_id}.csv" )
              
              if [[ ${#label_match[@]} -gt 0 ]]; then
                  label_path="${label_match[0]}"
                  
                  # Define Output Prefix
                  if [[ "$INSEC" == "T" ]]; then
                      out_prefix="${OUT_DIR}/eval_${method}_${pop}_seed${seed_id}_${GENE}_intersect"
                  else
                      out_prefix="${OUT_DIR}/eval_${method}_${pop}_seed${seed_id}_${GENE}"
                  fi
                  
                  # Build Command
                  CMD=(python "${PY_SCRIPT}")
                  CMD+=(--gt "${label_path}")
                  CMD+=(--infer "${infer_file}")
                  CMD+=(--pmtr "${PMTR_NAME}")
                  CMD+=(--index "${TRBV_INDEX}")
                  CMD+=(--gene_type "${GENE}")
                  CMD+=(--out_prefix "${out_prefix}")
                  
                  if [[ "$INSEC" == "T" ]]; then
                      CMD+=(--intersect)
                  fi
                  
                  # Execute
                  if "${CMD[@]}"; then
                      echo "    [OK] seed${seed_id}"
                      count=$((count + 1))
                  else
                      echo "    [FAIL] Python script failed for seed${seed_id}"
                  fi
              else
                  echo "    [WARN] No matching genotype file for seed${seed_id} in ${LABEL_DIR}"
              fi
          fi
      done
      
      if [[ $count -eq 0 ]]; then
          echo "    [WARN] No files processed for this combination."
      else
          echo "    [DONE] Processed $count samples."
      fi
      shopt -u nullglob

    done # End Gene
  done # End Pop
done # End Method

echo "======================================================="
echo " Evaluation Task Finished."
echo " Log: ${LOG_FILE}"
echo "======================================================="