#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
# ======================================================
# infer_FindAlleles.sh
#
# Logic:
#   1. Input: sample/{EXPR_ID}/{POP}/ERR*
#             (Looks for *.alleles.tsv and *.customAlleles.json)
#   2. Process: Runs infer_for_findalleles.py
#   3. Output: results/infer/{EXPR_ID}/FindAlleles/{POP}/{GENE}/
# ======================================================

# ---------------- Defaults ----------------
POPS="MIX"
GENES="V J"
EXPR_ID="expr_AIRR1"

SAMPLE_ROOT="${SCRIPT_DIR}/generated/samples"
RESULTS_ROOT="${SCRIPT_DIR}/generated/results"

# Reference files (Adjust paths if needed)
REF_V="${PACKAGE_ROOT}/data/ref/IMGT_TRB_default.csv"
REF_J="${PACKAGE_ROOT}/data/ref/IMGT_TRB_default.csv"

# Python script
PYTHON="${PYTHON:-python}"
PY_SCRIPT="${PACKAGE_ROOT}/code/experiments/00_benchmark_utils/infer_for_findalleles.py"

# ---------------- Arg Parsing ----------------
usage() {
  cat <<EOF
Usage:
  bash $0 --expr expr_AIRR1 --pop "AFR EUR" --gene "V J"

Options:
  --expr          Experiment ID (required)
  --pop           Population list (default: "${POPS}")
  --gene          Gene list (default: "${GENES}")
  --sample-root   Root folder containing sample/ (default: ${SAMPLE_ROOT})
  --results-root  Results root folder (default: ${RESULTS_ROOT})
  --ref-v         Reference file for V gene (default: ${REF_V})
  --ref-j         Reference file for J gene (default: ${REF_J})
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --expr) EXPR_ID="$2"; shift 2 ;;
    --pop) POPS="$2"; shift 2 ;;
    --gene) GENES="$2"; shift 2 ;;
    --sample-root) SAMPLE_ROOT="$2"; shift 2 ;;
    --results-root) RESULTS_ROOT="$2"; shift 2 ;;
    --ref-v) REF_V="$2"; shift 2 ;;
    --ref-j) REF_J="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "${EXPR_ID}" ]]; then
  echo "[ERROR] --expr is required."
  exit 1
fi

# ---------------- Log Setup ----------------
LOG_DIR="logs/${EXPR_ID}"
mkdir -p "${LOG_DIR}"

TIMESTAMP=$(date "+%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/infer_FindAlleles_${TIMESTAMP}.log"

# STRICT LOGGING: Only to file, nothing to console
exec > "$LOG_FILE" 2>&1

echo "======================================================="
echo " FindAlleles Inference Task Started"
echo " EXPR_ID:    ${EXPR_ID}"
echo " POPS:       ${POPS}"
echo " GENES:      ${GENES}"
echo " INPUT:      ${SAMPLE_ROOT}/${EXPR_ID}/{POP}/ERR..."
echo " OUTPUT:     ${RESULTS_ROOT}/infer/${EXPR_ID}/FindAlleles/..."
echo " LOG:        ${LOG_FILE}"
echo "======================================================="

# ---------------- Helpers ----------------
find_files_in_sample() {
  local dir="$1"
  # Look for any file ending in .alleles.tsv and .customAlleles.json
  # This makes it robust against naming conventions
  
  local tsv
  tsv=$(find "$dir" -maxdepth 1 -name "*.alleles.tsv" | head -n 1 || true)
  
  local json
  json=$(find "$dir" -maxdepth 1 -name "*.customAlleles.json" | head -n 1 || true)
  
  if [[ -n "$tsv" && -n "$json" ]]; then
      echo "$tsv|$json"
      return 0
  fi
  echo ""
  return 1
}

# ---------------- Main Execution ----------------

# 1. Iterate Populations
for pop in ${POPS//,/ }; do
  
  # Input Path: sample/expr/pop
  INPUT_POP_DIR="${SAMPLE_ROOT}/${EXPR_ID}/${pop}"
  
  if [[ ! -d "$INPUT_POP_DIR" ]]; then
      echo "[WARN] Population directory not found: ${INPUT_POP_DIR}. Skipping."
      continue
  fi

  echo ">>> Processing Population: ${pop}"

  # 2. Iterate Genes (V / J)
  for gene in ${GENES//,/ }; do
    gene=$(echo "$gene" | tr '[:lower:]' '[:upper:]') # Ensure uppercase
    
    # Select Reference
    if [[ "$gene" == "V" ]]; then
       CUR_REF="$REF_V"
    elif [[ "$gene" == "J" ]]; then
       CUR_REF="$REF_J"
    else
       echo "  [WARN] Unknown gene type: $gene. Skipping."
       continue
    fi
    
    if [[ ! -f "$CUR_REF" ]]; then
       echo "  [ERROR] Reference file not found for ${gene}: ${CUR_REF}. Skipping."
       continue
    fi

    # Output Path: results/infer/expr/FindAlleles/pop/gene
    OUT_DIR="${RESULTS_ROOT}/infer/${EXPR_ID}/FindAlleles/${pop}/${gene}"
    mkdir -p "$OUT_DIR"

    echo "  -> Gene: ${gene} | Output: ${OUT_DIR}"

    # 3. Iterate samples
    shopt -s nullglob
    sample_dirs=( "${INPUT_POP_DIR}"/ERR* )
    
    if [[ ${#sample_dirs[@]} -eq 0 ]]; then
       echo "     [WARN] No sample directories found for ${pop}."
    fi

    for sample_path in "${sample_dirs[@]}"; do
        [[ -d "$sample_path" ]] || continue

        sample_base=$(basename "$sample_path")

        # Find input files (TSV + JSON)
        files_found=$(find_files_in_sample "$sample_path")
        
        if [[ -n "$files_found" ]]; then
            # Split the result string "tsv|json"
            IFS='|' read -r input_tsv input_json <<< "$files_found"
            
            # Define Output Filename
            # infer_{sample_base}.{GENE}.csv
            output_filename="infer_${sample_base}.${gene}.csv"
            output_path="${OUT_DIR}/${output_filename}"

            # echo "     Processing ${sample_base}..."

            # Call Python Script
            if "${PYTHON}" "$PY_SCRIPT" \
                --tsv "$input_tsv" \
                --json "$input_json" \
                --gene "$gene" \
                --ref "$CUR_REF" \
                --output "$output_path"; then
                
                echo "     [OK] ${sample_base} -> ${output_filename}"
            else
                echo "     [FAIL] Python script failed on ${sample_base}"
            fi
        else
            echo "     [SKIP] ${sample_base}: Missing alleles.tsv or customAlleles.json"
        fi
    done
    shopt -u nullglob
  done # End Gene
  echo "------------------------------------------------------------"
done # End Pop

echo "============================================================"
echo " FindAlleles Inference Completed at: $(date)"
echo " Log: ${LOG_FILE}"
echo "============================================================"

