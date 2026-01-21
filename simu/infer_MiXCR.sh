#!/usr/bin/env bash
set -euo pipefail

# ======================================================
# infer_MiXCR.sh (v2 Multi-Pop/Multi-Gene)
#
# Logic:
#   1. Iterates over defined POPs and GENEs (V, J).
#   2. Scans input: sample/{EXPR_ID}/{POP}/seed*
#   3. Runs clones2alleles.py
#   4. Output: results/infer/{EXPR_ID}/MiXCR/{POP}/{GENE}/
# ======================================================

# ---------------- Defaults ----------------
POPS="AFR EUR AMR EAS SAS"
GENES="V J"
EXPR_ID="expr_0"

SAMPLE_ROOT="samples"
RESULTS_ROOT="results"

# Default Reference Files (Adjust paths if needed)
REF_V="ref/IMGT_TRB_default.csv"
REF_J="ref/IMGT_TRB_default.csv"

# Python Script
PY_SCRIPT="clones2alleles.py"

# ---------------- Arg Parsing ----------------
usage() {
  cat <<EOF
Usage:
  bash $0 --expr expr_0 --pop "AFR EUR" --gene "V J"

Options:
  --expr          Experiment ID (required)
  --pop           Population list (space separated) (default: "${POPS}")
  --gene          Gene list (V, J or "V J") (default: "${GENES}")
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
LOG_FILE="${LOG_DIR}/infer_MiXCR_${TIMESTAMP}.log"

# Redirect output to log and console
exec > "$LOG_FILE" 2>&1

echo "======================================================="
echo " MiXCR Inference Task Started"
echo " EXPR_ID:    ${EXPR_ID}"
echo " POPS:       ${POPS}"
echo " GENES:      ${GENES}"
echo " INPUT:      ${SAMPLE_ROOT}/${EXPR_ID}/{POP}/seed..."
echo " OUTPUT:     ${RESULTS_ROOT}/infer/${EXPR_ID}/MiXCR/..."
echo " LOG:        ${LOG_FILE}"
echo "======================================================="

# ---------------- Helpers ----------------
find_clones_file() {
  local dir="$1"
  # Priority 1: standard clones_TRB export
  local f
  f=$(find "$dir" -maxdepth 2 -type f -name "*clones_TRB.tsv" 2>/dev/null | head -n 1 || true)
  if [[ -n "${f}" ]]; then echo "${f}"; return 0; fi

  # Priority 2: standard MiXCR export
  f=$(find "$dir" -maxdepth 2 -type f -name "*.clones.tsv" 2>/dev/null | head -n 1 || true)
  if [[ -n "${f}" ]]; then echo "${f}"; return 0; fi
  
  # Priority 3: Loose match
  f=$(find "$dir" -maxdepth 2 -type f -name "*clones*.tsv" 2>/dev/null | head -n 1 || true)
  if [[ -n "${f}" ]]; then echo "${f}"; return 0; fi

  echo ""
  return 1
}

# ---------------- Main Execution ----------------

# 1. Iterate Populations
for pop in ${POPS//,/ }; do
  
  INPUT_POP_DIR="${SAMPLE_ROOT}/${EXPR_ID}/${pop}"
  
  if [[ ! -d "$INPUT_POP_DIR" ]]; then
      echo "[WARN] Population directory not found: ${INPUT_POP_DIR}. Skipping."
      continue
  fi

  echo ">>> Processing Population: ${pop}"

  # 2. Iterate Genes (V and/or J)
  for gene in ${GENES//,/ }; do
    gene=$(echo "$gene" | tr '[:lower:]' '[:upper:]') # Ensure uppercase
    
    # Select Reference File based on Gene
    if [[ "$gene" == "V" ]]; then
       CUR_REF="$REF_V"
    elif [[ "$gene" == "J" ]]; then
       CUR_REF="$REF_J"
    else
       echo "  [WARN] Unknown gene type: $gene. Skipping."
       continue
    fi
    
    if [[ ! -f "$CUR_REF" ]]; then
       echo "  [ERROR] Reference file not found for ${gene}: ${CUR_REF}. Skipping this gene."
       continue
    fi

    # Define Output Directory
    # results/infer/{EXPR_ID}/MiXCR/{POP}/{GENE}
    OUT_DIR="${RESULTS_ROOT}/infer/${EXPR_ID}/MiXCR/${pop}/${gene}"
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

        # Find input clones file in the seed directory
        input_file=$(find_clones_file "$seed_path")

        if [[ -n "$input_file" ]]; then
            # Output Filename: infer_{POP}_seed{seed}.{GENE}.csv
            output_filename="infer_${pop}_seed${seed_num}.${gene}.csv"
            output_path="${OUT_DIR}/${output_filename}"

            # Only run if output doesn't exist or we want to overwrite (optional logic)
            # Here we just run it.
            
            # echo "     Processing ${seed_base}..."

            # Call Python Script
            # Ensure python script handles errors gracefully
            if python "$PY_SCRIPT" \
                --clones "$input_file" \
                --gene "$gene" \
                --ref "$CUR_REF" \
                --output "$output_path" ; then
                
                # Verify output creation
                if [[ -f "$output_path" ]]; then
                     echo "     [OK] ${seed_base} -> ${output_filename}"
                else
                     echo "     [FAIL] ${seed_base}: Output file not created."
                fi
            else
                echo "     [FAIL] Python script error on ${seed_base}"
            fi
        else
            echo "     [SKIP] ${seed_base}: No clones file found."
        fi
    done
    shopt -u nullglob
  done # End Gene Loop
  echo "------------------------------------------------------------"
done # End Pop Loop

echo "============================================================"
echo " Inference Completed at: $(date)"
echo " Log: ${LOG_FILE}"
echo "============================================================"