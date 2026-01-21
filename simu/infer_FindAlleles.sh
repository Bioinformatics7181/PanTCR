#!/usr/bin/env bash
set -euo pipefail

# ======================================================
# infer_FindAlleles.sh
#
# Logic:
#   1. Input: sample/{EXPR_ID}/{POP}/seed{N}
#             (Looks for *.alleles.tsv and *.customAlleles.json)
#   2. Process: Runs infer_for_findalleles.py
#   3. Output: results/infer/{EXPR_ID}/FindAlleles/{POP}/{GENE}/
# ======================================================

# ---------------- Defaults ----------------
POPS="AFR EUR AMR EAS SAS"
GENES="V J"
EXPR_ID="expr_0"

SAMPLE_ROOT="samples"
RESULTS_ROOT="results"

# Reference files (Adjust paths if needed)
REF_V="ref/IMGT_TRBV_pro.tsv"
REF_J="ref/IMGT_TRBJ_pro.tsv"

# Python script
PY_SCRIPT="infer_for_findalleles.py"

# ---------------- Arg Parsing ----------------
usage() {
  cat <<EOF
Usage:
  bash $0 --expr expr_0 --pop "AFR EUR" --gene "V J"

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
echo " INPUT:      ${SAMPLE_ROOT}/${EXPR_ID}/{POP}/seed..."
echo " OUTPUT:     ${RESULTS_ROOT}/infer/${EXPR_ID}/FindAlleles/..."
echo " LOG:        ${LOG_FILE}"
echo "======================================================="

# ---------------- Helpers ----------------
find_files_in_seed() {
  local dir="$1"
  # Look for any file ending in .alleles.tsv and .customAlleles.json
  # This makes it robust against naming conventions (e.g. POP_seedX or POP_default_seedX)
  
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

    # 3. Iterate Seeds
    shopt -s nullglob
    seed_dirs=( "${INPUT_POP_DIR}"/seed* )
    
    if [[ ${#seed_dirs[@]} -eq 0 ]]; then
       echo "     [WARN] No seed directories found for ${pop}."
    fi

    for seed_path in "${seed_dirs[@]}"; do
        seed_base=$(basename "$seed_path")       # e.g. seed0
        seed_num="${seed_base#seed}"             # e.g. 0
        
        if ! [[ "$seed_num" =~ ^[0-9]+$ ]]; then continue; fi

        # Find input files (TSV + JSON)
        files_found=$(find_files_in_seed "$seed_path")
        
        if [[ -n "$files_found" ]]; then
            # Split the result string "tsv|json"
            IFS='|' read -r input_tsv input_json <<< "$files_found"
            
            # Define Output Filename
            # infer_{POP}_seed{seed}.{GENE}.csv
            output_filename="infer_${pop}_seed${seed_num}.${gene}.csv"
            output_path="${OUT_DIR}/${output_filename}"

            # echo "     Processing ${seed_base}..."

            # Call Python Script
            if python "$PY_SCRIPT" \
                --tsv "$input_tsv" \
                --json "$input_json" \
                --gene "$gene" \
                --ref "$CUR_REF" \
                --output "$output_path"; then
                
                echo "     [OK] ${seed_base} -> ${output_filename}"
            else
                echo "     [FAIL] Python script failed on ${seed_base}"
            fi
        else
            echo "     [SKIP] ${seed_base}: Missing alleles.tsv or customAlleles.json"
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