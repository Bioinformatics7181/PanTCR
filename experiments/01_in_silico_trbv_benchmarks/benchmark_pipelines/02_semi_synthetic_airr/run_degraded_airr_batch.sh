#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
# =========================
# Default Parameters
# =========================
# Directory Settings
EXPR_ID="expr_AIRR1"
INPUT_ROOT="${SCRIPT_DIR}/input/full_length_airr"
OUTPUT_ROOT="${SCRIPT_DIR}/generated/samples"            # Base output (Final: samples/expr_AIRR1/POP/SAMPLE)

# Populations to process (Space separated string)
# Default: "MIX". Example for multi: "AFR EUR AMR EAS SAS"
POPS="MIX" 

# Script & Resource Paths
PY_DEGRADE="${PACKAGE_ROOT}/code/experiments/07_semi_synthetic_airr_benchmark/degrade_airr_fastq_pairs.py"
PY_INFER="${PACKAGE_ROOT}/code/experiments/00_benchmark_utils/infer_for_findalleles.py"
REF_FILE="${PACKAGE_ROOT}/data/ref/IMGT_TRB_default.csv" 

PYTHON_EXE="python"

# MiXCR Settings
THREADS=8
MIXCR_SPECIES="hsa"

# Degradation Settings
KEEP_PROB=1
P_DEGRADED=1
CUT_DEGRADED_MEAN=100
CUT_DEGRADED_SD=12
CUT_MAX=200

CUT_INTACT_MEAN=0
CUT_INTACT_SD=0
CUT_MAX_INTACT=0

MIN_READ_LEN=50
SEED=42

# Trimming Options
TRIM_R1=1
R1_DIR="5"
R1_FIXED_TRIM=0

TRIM_R2=0
R2_DIR="3"
R2_FIXED_TRIM=0

MAX_PAIRS=0

DO_CLEANUP=1

# =========================
# Args Parsing
# =========================
while [[ $# -gt 0 ]]; do
  case "$1" in
    --expr) EXPR_ID="$2"; shift 2 ;;
    --input_root) INPUT_ROOT="$2"; shift 2 ;;
    --output_root) OUTPUT_ROOT="$2"; shift 2 ;;
    --pops) POPS="$2"; shift 2 ;;  # Accepts string: "AFR EUR"
    
    --ref_file) REF_FILE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    
    --keep_prob) KEEP_PROB="$2"; shift 2 ;;
    --p_degraded) P_DEGRADED="$2"; shift 2 ;;
    --cut_degraded_mean) CUT_DEGRADED_MEAN="$2"; shift 2 ;;
    --cut_degraded_sd) CUT_DEGRADED_SD="$2"; shift 2 ;;
    --cut_max) CUT_MAX="$2"; shift 2 ;;
    
    --cut_intact_mean) CUT_INTACT_MEAN="$2"; shift 2 ;;
    --cut_intact_sd) CUT_INTACT_SD="$2"; shift 2 ;;
    --cut_max_intact) CUT_MAX_INTACT="$2"; shift 2 ;;
    
    --min_read_len) MIN_READ_LEN="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    
    --trim_r2_on) TRIM_R2=1; shift 1 ;;
    --r1_fixed_trim) R1_FIXED_TRIM="$2"; shift 2 ;;

    --no_cleanup) DO_CLEANUP=0; shift 1 ;;
    
    -h|--help)
      echo "Usage: $0 --expr expr_AIRR1 --pops 'AFR EUR' --input_root input/full_length_airr"
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown arg: $1"
      exit 1
      ;;
  esac
done

# =========================
# Setup Paths & Logs
# =========================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FULL_PY_DEGRADE="${PY_DEGRADE}"
FULL_PY_INFER="${PY_INFER}"

if [[ ! -f "$FULL_PY_DEGRADE" ]]; then echo "[ERROR] Missing: $FULL_PY_DEGRADE"; exit 1; fi
if [[ ! -f "$FULL_PY_INFER" ]]; then echo "[ERROR] Missing: $FULL_PY_INFER"; exit 1; fi

# Remove trailing slash
INPUT_ROOT="${INPUT_ROOT%/}"

# Define Base Output Directory
BASE_OUTPUT_DIR="${OUTPUT_ROOT}/${EXPR_ID}"
LOG_DIR="logs/${EXPR_ID}"

mkdir -p "$LOG_DIR"
mkdir -p "$BASE_OUTPUT_DIR"

TIMESTAMP=$(date "+%Y%m%d_%H%M%S")
MAIN_LOG="${LOG_DIR}/run_pipeline_${EXPR_ID}_${TIMESTAMP}.log"

exec > >(tee -i "$MAIN_LOG") 2>&1

echo "======================================================="
echo " Pipeline Started"
echo " EXPR_ID:    ${EXPR_ID}"
echo " POPS:       ${POPS}"
echo " INPUT:      ${INPUT_ROOT}"
echo " OUTPUT:     ${BASE_OUTPUT_DIR}"
echo "======================================================="

# Save Config
CONFIG_FILE="${BASE_OUTPUT_DIR}/${EXPR_ID}_config.tsv"
if [[ ! -f "$CONFIG_FILE" ]]; then
    printf "timestamp\texpr_id\tinput_root\tpops\tseed\n" > "$CONFIG_FILE"
fi
printf "%s\t%s\t%s\t%s\t%s\n" \
  "$(date -Iseconds)" "$EXPR_ID" "$INPUT_ROOT" "$POPS" "$SEED" >> "$CONFIG_FILE"

# =========================
# Main Processing Loop
# =========================

if [[ ! -d "$INPUT_ROOT" ]]; then
    echo "[ERROR] Input directory not found: $INPUT_ROOT"
    exit 1
fi

# Loop through each population defined in POPS
for POP in ${POPS}; do

    POP_INPUT_DIR="${INPUT_ROOT}/${POP}"
    
    if [[ ! -d "$POP_INPUT_DIR" ]]; then
        echo "[WARN] Population directory not found: $POP_INPUT_DIR. Skipping."
        continue
    fi

    echo ">>> Processing Population: $POP"
    echo "    Scanning: $POP_INPUT_DIR"

    # Find samples within specific population folder
    find "$POP_INPUT_DIR" -mindepth 1 -maxdepth 2 -type f \( -name "*_1.fastq.gz" -o -name "*_1.fastq" -o -name "*_1.fq.gz" -o -name "*_1.fq" \) | sort | while read -r R1_FILE; do
        
        # 1. Parse Paths
        SAMPLE_DIR=$(dirname "$R1_FILE")
        SAMPLE_NAME=$(basename "$SAMPLE_DIR")
        R1_FILENAME=$(basename "$R1_FILE")
        
        # Target: samples/expr_AIRR1/[POP]/[SAMPLE_NAME]
        TARGET_DIR="${BASE_OUTPUT_DIR}/${POP}/${SAMPLE_NAME}"
        mkdir -p "$TARGET_DIR"
        
        # Check R2
        R2_FILENAME="${R1_FILENAME/_1/_2}"
        R2_FILE="${SAMPLE_DIR}/${R2_FILENAME}"
        
        if [[ ! -f "$R2_FILE" ]]; then
            echo "[WARN] Skipping $SAMPLE_NAME: R2 not found"
            continue
        fi
        
        echo "-------------------------------------------------------"
        echo "Sample: $SAMPLE_NAME ($POP)"
        
        OUT_R1="${TARGET_DIR}/${R1_FILENAME}"
        OUT_R2="${TARGET_DIR}/${R2_FILENAME}"
        TRIM_LOG="${TARGET_DIR}/out.trim.csv"
        
        # 2. Degrade
        echo "[STEP 1] Degradation..."
        CMD_ARGS=(
            "-1" "$R1_FILE" "-2" "$R2_FILE"
            "-o1" "$OUT_R1" "-o2" "$OUT_R2"
            "--log" "$TRIM_LOG"
            "--p_degraded" "$P_DEGRADED"
            "--cut_degraded_mean" "$CUT_DEGRADED_MEAN"
            "--cut_degraded_sd" "$CUT_DEGRADED_SD"
            "--cut_max" "$CUT_MAX"
            "--cut_intact_mean" "$CUT_INTACT_MEAN"
            "--cut_intact_sd" "$CUT_INTACT_SD"
            "--cut_max_intact" "$CUT_MAX_INTACT"
            "--min_read_len" "$MIN_READ_LEN"
            "--r1_fixed_trim" "$R1_FIXED_TRIM"
            "--keep_prob" "$KEEP_PROB"
            "--seed" "$SEED"
        )
        
        if [[ "$TRIM_R1" -eq 1 ]]; then CMD_ARGS+=("--trim_r1" "--r1_dir" "$R1_DIR"); fi
        if [[ "$TRIM_R2" -eq 1 ]]; then CMD_ARGS+=("--trim_r2" "--r2_dir" "$R2_DIR"); fi
        if [[ "$MAX_PAIRS" -gt 0 ]]; then CMD_ARGS+=("--max_pairs" "$MAX_PAIRS"); fi

        $PYTHON_EXE "$FULL_PY_DEGRADE" "${CMD_ARGS[@]}"

        # 3. Inference
        echo "[STEP 2] Inference..."
        INFER_TSV="${SAMPLE_DIR}/${SAMPLE_NAME}.alleles.tsv"
        INFER_JSON="${SAMPLE_DIR}/${SAMPLE_NAME}.customAlleles.json"
        INFER_OUT="${TARGET_DIR}/genotype_${SAMPLE_NAME}.csv"
        
        if [[ -f "$INFER_TSV" && -f "$INFER_JSON" ]]; then
            $PYTHON_EXE "$FULL_PY_INFER" \
                --tsv "$INFER_TSV" --json "$INFER_JSON" \
                --ref "$REF_FILE" --output "$INFER_OUT"
        else
            echo "[WARN] Inference inputs missing. Skipping."
        fi

        # 4. MiXCR
        echo "[STEP 3] MiXCR Pipeline..."
        if ! command -v mixcr &> /dev/null; then echo "[ERROR] mixcr not found."; exit 1; fi
        
        pushd "$TARGET_DIR" > /dev/null
        
        mixcr analyze -s "$MIXCR_SPECIES" --threads "$THREADS" rna-seq \
            --force-overwrite "$R1_FILENAME" "$R2_FILENAME" "$SAMPLE_NAME"
            
        mixcr assembleContigs --force-overwrite \
            --assemble-contigs-by VDJRegion \
            "${SAMPLE_NAME}.clna" "${SAMPLE_NAME}.contigs.VDJRegion.clns"
            
        mixcr findAlleles --force-overwrite --no-clns-output \
            --report "${SAMPLE_NAME}.findAlleles.report.txt" \
            --export-library "${SAMPLE_NAME}.customAlleles.json" \
            --export-alleles-mutations "${SAMPLE_NAME}.alleles.tsv" \
            "${SAMPLE_NAME}.contigs.VDJRegion.clns"

        if [[ "$DO_CLEANUP" -eq 1 ]]; then
            rm -f ./*.vdjca ./*.clna ./*.clns || true
        fi
            
        popd > /dev/null
        
        echo "[DONE] Finished $SAMPLE_NAME."
        
    done
done

echo "======================================================="
echo " Pipeline Complete."
echo " Output: ${BASE_OUTPUT_DIR}"
echo "======================================================="

