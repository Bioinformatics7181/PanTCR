#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"

# ==============================================================================
# Simulation & Validation Pipeline
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Work Dir Setting
# ------------------------------------------------------------------------------

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

echo "[INIT] Script is located at: $SCRIPT_DIR"
echo "[INIT] Original working dir: $(pwd)"

if [[ "$(pwd)" != "$SCRIPT_DIR" ]]; then
    echo "[INIT] Switching working directory to: $SCRIPT_DIR"
    cd "$SCRIPT_DIR"
fi

# ------------------------------------------------------------------------------
# 1. Default Configuration
# ------------------------------------------------------------------------------

# Execution Control
TARGET_STEPS="1 2 3 4 5 6 7 8 9 10"  # Default: Run all steps

# Experiment Identity
EXPR_ID="expr_AIRR1"
POPS="MIX"
GENES="V J"

# Simulation: General
SIM_KEEP_PROB=1

# Simulation: Degradation
SIM_P_DEGRADED=0.7

# Simulation: Fragmentation / Cutting
SIM_CUT_INTACT_MEAN=0
SIM_CUT_INTACT_SD=0
SIM_CUT_MAX_INTACT=0
SIM_CUT_DEGRADED_MEAN=100
SIM_CUT_DEGRADED_SD=12
SIM_CUT_MAX_DEGRADED=200
SIM_MIN_KEEP_LEN=50

# Validation / Inference
VAL_K_FOLDS=5
VAL_SPLIT_SEED=42
VAL_MIN_NAIVE=2
VAL_PENALTY_K=3
VAL_PI_MIN=0.2
VAL_METHODS="MiXCR FindAlleles Bayes BayesNoPrior"

# Paths
PATH_INPUT="${SCRIPT_DIR}/input/full_length_airr"
PATH_SAMPLE="${SCRIPT_DIR}/generated/samples"
PATH_RESULTS="${SCRIPT_DIR}/generated/results"
PATH_METADATA="${PACKAGE_ROOT}/data/ref/metadata.csv"

# ------------------------------------------------------------------------------
# 2. Argument Parsing
# ------------------------------------------------------------------------------

show_help() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Execution Control:"
    echo "  --steps STR             Space-separated steps to run (e.g., \"1 3 10\")"
    echo "                          Default: \"${TARGET_STEPS}\""
    echo ""
    echo "Experiment Identity:"
    echo "  --expr-id STR           Experiment ID (default: ${EXPR_ID})"
    echo "  --pops STR              Populations (default: \"${POPS}\")"
    echo "  --genes STR             Genes (default: \"${GENES}\")"
    echo ""
    echo "Simulation - General:"
    echo "  --keep_prob FLOAT       Keep prob (default: ${SIM_KEEP_PROB})"
    echo ""
    echo "Simulation - Degradation:"
    echo "  --p-deg FLOAT           Degradation prob (default: ${SIM_P_DEGRADED})"
    echo ""
    echo "Simulation - Fragmentation:"
    echo "  --cut-intact-mean INT   Intact cut mean (default: ${SIM_CUT_INTACT_MEAN})"
    echo "  --cut-intact-sd INT     Intact cut SD (default: ${SIM_CUT_INTACT_SD})"
    echo "  --cut-max-intact INT    Intact max cut (default: ${SIM_CUT_MAX_INTACT})"
    echo "  --cut-deg-mean INT      Degraded cut mean (default: ${SIM_CUT_DEGRADED_MEAN})"
    echo "  --cut-deg-sd INT        Degraded cut SD (default: ${SIM_CUT_DEGRADED_SD})"
    echo "  --cut-max-deg INT       Degraded max cut (default: ${SIM_CUT_MAX_DEGRADED})"
    echo "  --min-keep INT          Min keep length (default: ${SIM_MIN_KEEP_LEN})"
    echo ""
    echo "Validation:"
    echo "  --k-folds INT           K-Folds (default: ${VAL_K_FOLDS})"
    echo "  --split-seed INT        Split seed (default: ${VAL_SPLIT_SEED})"
    echo "  --min-naive INT         Min naive threshold (default: ${VAL_MIN_NAIVE})"
    echo "  --methods STR           Methods list (default: \"${VAL_METHODS}\")"
    echo ""
    echo "Paths:"
    echo "  --path-input STR        Input dir (default: ${PATH_INPUT})"
    echo "  --path-sample STR       Sample dir (default: ${PATH_SAMPLE})"
    echo "  --path-results STR      Results dir (default: ${PATH_RESULTS})"
    echo "  --path-metadata STR     Metadata file (default: ${PATH_METADATA})"
    echo ""
    echo "Misc:"
    echo "  --help, -h              Show this help"
    exit 0
}

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --steps)            TARGET_STEPS="$2"; shift ;;
        --expr-id)          EXPR_ID="$2"; shift ;;
        --pops)             POPS="$2"; shift ;;
        --genes)            GENES="$2"; shift ;;
        --keep_prob)        SIM_KEEP_PROB="$2"; shift ;;
        --p-deg)            SIM_P_DEGRADED="$2"; shift ;;
        --cut-intact-mean)  SIM_CUT_INTACT_MEAN="$2"; shift ;;
        --cut-intact-sd)    SIM_CUT_INTACT_SD="$2"; shift ;;
        --cut-max-intact)   SIM_CUT_MAX_INTACT="$2"; shift ;;
        --cut-deg-mean)     SIM_CUT_DEGRADED_MEAN="$2"; shift ;;
        --cut-deg-sd)       SIM_CUT_DEGRADED_SD="$2"; shift ;;
        --cut-max-deg)      SIM_CUT_MAX_DEGRADED="$2"; shift ;;
        --min-keep)         SIM_MIN_KEEP_LEN="$2"; shift ;;
        --k-folds)          VAL_K_FOLDS="$2"; shift ;;
        --split-seed)       VAL_SPLIT_SEED="$2"; shift ;;
        --min-naive)        VAL_MIN_NAIVE="$2"; shift ;;
        --penalty-k)        VAL_PENALTY_K="$2"; shift ;;
        --pi-min)           VAL_PI_MIN="$2"; shift ;;
        --methods)          VAL_METHODS="$2"; shift ;;
        --path-input)       PATH_INPUT="$2"; shift ;;
        --path-sample)      PATH_SAMPLE="$2"; shift ;;
        --path-results)     PATH_RESULTS="$2"; shift ;;
        --path-metadata)    PATH_METADATA="$2"; shift ;;
        --help|-h)          show_help ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# ------------------------------------------------------------------------------
# 3. Helper Functions & Logging
# ------------------------------------------------------------------------------

should_run_step() {
    local step=$1
    for target in ${TARGET_STEPS}; do
        if [[ "$step" == "$target" ]]; then return 0; fi
    done
    return 1
}

mkdir -p logs
PIPELINE_LOG="logs/pipeline_${EXPR_ID}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$PIPELINE_LOG") 2>&1

echo "################################################################"
echo " Experiment Configuration Summary"
echo "################################################################"
echo " [0] Meta Information"
echo "   EXPR_ID        : ${EXPR_ID}"
echo "   Populations    : ${POPS}"
echo "   Start Time     : $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo " [1] Axis 1: Bias Intensity (Degradation)"
echo "   SIM_P_DEGRADED : ${SIM_P_DEGRADED}  (Higher = Flatter Likelihood)"
echo ""
echo " [2] Axis 2: Information Loss (Truncation)"
echo "   SIM_CUT_MEAN   : ${SIM_CUT_DEGRADED_MEAN} bp (Mean degraded cut length)"
echo "   SIM_CUT_MAX    : ${SIM_CUT_MAX_DEGRADED} bp"
echo ""
echo " [3] Axis 3: Effective Evidence (Depth & Diversity)"
echo "   KEEP_PROB      : ${SIM_KEEP_PROB} (Total Reads / Depth)"
echo "################################################################"

# ------------------------------------------------------------------------------
# 4. Pipeline Steps
# ------------------------------------------------------------------------------

# Step 1: Simulation
if should_run_step 1; then
    echo -e "\n[Step 1/10] Running Simulation..."
    bash "${SCRIPT_DIR}/run_degraded_airr_batch.sh" \
      --expr "${EXPR_ID}" \
      --input_root "${PATH_INPUT}" \
      --output_root "${PATH_SAMPLE}" \
      --pops "${POPS}" \
      --keep_prob "${SIM_KEEP_PROB}" \
      --p_degraded "${SIM_P_DEGRADED}" \
      --cut_intact_mean "${SIM_CUT_INTACT_MEAN}" \
      --cut_intact_sd "${SIM_CUT_INTACT_SD}" \
      --cut_max_intact "${SIM_CUT_MAX_INTACT}" \
      --cut_degraded_mean "${SIM_CUT_DEGRADED_MEAN}" \
      --cut_degraded_sd "${SIM_CUT_DEGRADED_SD}" \
      --cut_max "${SIM_CUT_MAX_DEGRADED}" \
      --min_read_len "${SIM_MIN_KEEP_LEN}" \
      --seed "${VAL_SPLIT_SEED}"
fi

# Step 2: Analyze Mutations
if should_run_step 2; then
    echo -e "\n[Step 2/10] Analyzing Mutations..."
    bash "${SCRIPT_DIR}/collect_airr_mutations.sh" \
      --expr "${EXPR_ID}" \
      --pop "${POPS}" \
      --sample-root "${PATH_SAMPLE}" \
      --results-root "${PATH_RESULTS}"
fi

# Step 3: Copy Genotypes
if should_run_step 3; then
    echo -e "\n[Step 3/10] Preparing Genotype Labels..."
    bash "${SCRIPT_DIR}/copy_airr_genotypes.sh" \
      --expr "${EXPR_ID}" \
      --pop "${POPS}" \
      --sample-root "${PATH_SAMPLE}" \
      --results-root "${PATH_RESULTS}"
fi

# Step 4: K-Fold Split
if should_run_step 4; then
    echo -e "\n[Step 4/10] Splitting Data..."
    bash "${PACKAGE_ROOT}/experiments/01_in_silico_trbv_benchmarks/split_kfold.sh" \
      --expr "${EXPR_ID}" \
      --k "${VAL_K_FOLDS}" \
      --seed "${VAL_SPLIT_SEED}" \
      --results-root "${PATH_RESULTS}"
fi

# Step 5: MiXCR Inference
if should_run_step 5; then
    echo -e "\n[Step 5/10] Inference: MiXCR..."
    bash "${SCRIPT_DIR}/infer_airr_mixcr.sh" \
      --expr "${EXPR_ID}" \
      --pop "${POPS}" \
      --gene "${GENES}" \
      --sample-root "${PATH_SAMPLE}" \
      --results-root "${PATH_RESULTS}"
fi

# Step 6: FindAlleles Inference
if should_run_step 6; then
    echo -e "\n[Step 6/10] Inference: FindAlleles..."
    bash "${SCRIPT_DIR}/infer_airr_findalleles.sh" \
      --expr "${EXPR_ID}" \
      --pop "${POPS}" \
      --gene "${GENES}" \
      --sample-root "${PATH_SAMPLE}" \
      --results-root "${PATH_RESULTS}"
fi

# Step 7: Build Pangenome Graph
if should_run_step 7; then
    echo -e "\n[Step 7/10] Building Pangenome Graphs..."
    bash "${PACKAGE_ROOT}/experiments/01_in_silico_trbv_benchmarks/build_graph_batch.sh" \
      --expr "${EXPR_ID}" \
      --gene "${GENES}" \
      --k "${VAL_K_FOLDS}" \
      --min-naive "${VAL_MIN_NAIVE}" \
      --metadata "${PATH_METADATA}" \
      --results-root "${PATH_RESULTS}"
fi

# Step 8: Bayes Inference (With Prior)
if should_run_step 8; then
    echo -e "\n[Step 8/10] Inference: Bayes (With Prior)..."
    bash "${SCRIPT_DIR}/infer_airr_pantcr.sh" \
      --expr "${EXPR_ID}" \
      --gene "${GENES}" \
      --k "${VAL_K_FOLDS}" \
      --min-naive "${VAL_MIN_NAIVE}" \
      --penalty-k "${VAL_PENALTY_K}" \
      --pi-min "${VAL_PI_MIN}" \
      --results-root "${PATH_RESULTS}"
fi

# Step 9: Bayes Inference (No Prior)
if should_run_step 9; then
    echo -e "\n[Step 9/10] Inference: Bayes (No Prior)..."
    bash "${SCRIPT_DIR}/infer_airr_pantcr.sh" \
      --expr "${EXPR_ID}" \
      --gene "${GENES}" \
      --k "${VAL_K_FOLDS}" \
      --min-naive "${VAL_MIN_NAIVE}" \
      --penalty-k "${VAL_PENALTY_K}" \
      --pi-min "${VAL_PI_MIN}" \
      --results-root "${PATH_RESULTS}" \
      --no-prior
fi

# Step 10: Evaluation
if should_run_step 10; then
    echo -e "\n[Step 10/10] Evaluation..."
    bash "${SCRIPT_DIR}/eval_airr_alleles.sh" \
      --expr "${EXPR_ID}" \
      --method "${VAL_METHODS}" \
      --pop "${POPS}" \
      --gene "${GENES}" \
      --results-root "${PATH_RESULTS}"
fi

echo "################################################################"
echo " Pipeline Finished Successfully"
echo " End: $(date)"
echo " Check ${PATH_RESULTS}/eval/${EXPR_ID} for results."
echo "################################################################"


