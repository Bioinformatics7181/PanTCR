#!/usr/bin/env bash
set -euo pipefail

# =========================
# Defaults
# =========================
# Modified: Support multiple populations by default
POPS="AFR EUR AMR EAS SAS" 
SEED_START=0
N=1
EXPR_ID="expr_0"
OUT_ROOT="samples"

# Step1: genotype
INPUT_REF="ref/pmTR_TRB_V_J_cleaned.csv"
EXCLUDE_DEFAULT=""     # passed to select_genotype.py -d
FUNCTIONALITY=""       # passed to select_genotype.py -f

# Step2: repertoire
CDR3_DICT="ref/cdr3_dict.pkl"
NC=10000
NR=50000
ALPHA=1.5

# Step3: bias simulation
LEADER_TRIM_LEN=0
P_DEGRADED=0.7
P_DEGRADED_SD=0.08

LEN_MEAN=340
LEN_SD=6
LEN_MIN=330
LEN_MAX=360

CUT_INTACT_MEAN=0
CUT_INTACT_SD=0
CUT_MAX_INTACT=0

CUT_DEGRADED_MEAN=70
CUT_DEGRADED_SD=20
CUT_MAX=140

C_MAX=80
MIN_KEEP_LEN=250

# Step4: ART
ART_SS="HS25"
ART_L=150
ART_C=1
ART_M=320
ART_S=5

# MiXCR
MIXCR_THREADS=8
MIXCR_SPECIES="hsa"
MIXCR_PRESET="rna-seq"

# Cleanup
DO_CLEANUP=1

# =========================
# Args
# =========================
while [[ $# -gt 0 ]]; do
  case "$1" in
    --pop) POPS="$2"; shift 2 ;; # Modified: Input can now be a list like "AFR EUR"
    --seed_start) SEED_START="$2"; shift 2 ;;
    --n) N="$2"; shift 2 ;;
    --expr) EXPR_ID="$2"; shift 2 ;;
    --out_root) OUT_ROOT="$2"; shift 2 ;;

    --input_ref) INPUT_REF="$2"; shift 2 ;;
    --exclude_default) EXCLUDE_DEFAULT="$2"; shift 2 ;;
    --functionality) FUNCTIONALITY="$2"; shift 2 ;;

    --cdr3_dict) CDR3_DICT="$2"; shift 2 ;;
    --nc) NC="$2"; shift 2 ;;
    --nr) NR="$2"; shift 2 ;;
    --alpha) ALPHA="$2"; shift 2 ;;

    --leader_trim_len) LEADER_TRIM_LEN="$2"; shift 2 ;;
    --p_degraded) P_DEGRADED="$2"; shift 2 ;;
    --p_degraded_sd) P_DEGRADED_SD="$2"; shift 2 ;;

    --len_mean) LEN_MEAN="$2"; shift 2 ;;
    --len_sd) LEN_SD="$2"; shift 2 ;;
    --len_min) LEN_MIN="$2"; shift 2 ;;
    --len_max) LEN_MAX="$2"; shift 2 ;;

    --cut_intact_mean) CUT_INTACT_MEAN="$2"; shift 2 ;;
    --cut_intact_sd) CUT_INTACT_SD="$2"; shift 2 ;;
    --cut_max_intact) CUT_MAX_INTACT="$2"; shift 2 ;;

    --cut_degraded_mean) CUT_DEGRADED_MEAN="$2"; shift 2 ;;
    --cut_degraded_sd) CUT_DEGRADED_SD="$2"; shift 2 ;;
    --cut_max) CUT_MAX="$2"; shift 2 ;;

    --c_max) C_MAX="$2"; shift 2 ;;
    --min_keep_len) MIN_KEEP_LEN="$2"; shift 2 ;;

    --art_ss) ART_SS="$2"; shift 2 ;;
    --art_l) ART_L="$2"; shift 2 ;;
    --art_c) ART_C="$2"; shift 2 ;;
    --art_m) ART_M="$2"; shift 2 ;;
    --art_s) ART_S="$2"; shift 2 ;;

    --mixcr_threads) MIXCR_THREADS="$2"; shift 2 ;;
    --mixcr_species) MIXCR_SPECIES="$2"; shift 2 ;;
    --mixcr_preset) MIXCR_PRESET="$2"; shift 2 ;;

    --no_cleanup) DO_CLEANUP=0; shift 1 ;;
    -h|--help)
      echo "Usage:"
      echo "  $0 --pop \"AFR EUR\" --seed_start 0 --n 5 --expr expr_0"
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown arg: $1"
      exit 1
      ;;
  esac
done

# =========================
# Resolve script dir (so python paths are stable)
# =========================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# =========================
# EXCLUDE_DEFAULT -> suffix used in *sample_prefix*
# =========================
EXCLUDE_DEFAULT_UP="$(echo "$EXCLUDE_DEFAULT" | tr '[:lower:]' '[:upper:]' | tr -d ' ')"
TAG_SUFFIX=""
if [[ "$EXCLUDE_DEFAULT_UP" == "T" || "$EXCLUDE_DEFAULT_UP" == "TRUE" ]]; then
  TAG_SUFFIX="_exdefaultT"
fi

EXCLUDE_DEFAULT_SUFFIX=""
if [[ "$EXCLUDE_DEFAULT_UP" == "T" || "$EXCLUDE_DEFAULT_UP" == "TRUE" ]]; then
  EXCLUDE_DEFAULT_SUFFIX="_exdefaultT"
fi

# =========================
# Log setup
# =========================
LOG_DIR="logs/${EXPR_ID}"
mkdir -p $LOG_DIR
# Use the first pop in filename just for naming if multiple exist
FIRST_POP=$(echo $POPS | awk '{print $1}')
TIMESTAMP=$(date "+%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/run_simulation_${EXPR_ID}_${TIMESTAMP}.log"
exec > "$LOG_FILE" 2>&1

echo "======================================================="
echo " Analysis task started."
echo " POPS=${POPS}"
echo " EXPR_ID=${EXPR_ID}"
echo " LOG=${LOG_FILE}"
echo "======================================================="

# =========================
# Experiment Config (Global for this EXPR_ID)
# =========================
# Modified: EXPR_DIR is now the parent container
EXPR_DIR="${OUT_ROOT}/${EXPR_ID}"
mkdir -p "$EXPR_DIR"

# Modified: Save config file under EXPR_DIR, parallel to population folders
EXPERIMENTS_TSV="${EXPR_DIR}/experiments_config.tsv"

if [[ ! -f "$EXPERIMENTS_TSV" ]]; then
  printf "timestamp\texpr_id\tpop\tseed\texclude_default\tfunctionality\tnc\tnr\talpha\tp_degraded\tp_degraded_sd\tlen_mean\tlen_sd\tlen_min\tlen_max\tcut_intact_mean\tcut_intact_sd\tcut_max_intact\tcut_degraded_mean\tcut_degraded_sd\tcut_max\tc_max\tmin_keep_len\tart_ss\tart_l\tart_c\tart_m\tart_s\tmixcr_preset\tmixcr_threads\n" > "$EXPERIMENTS_TSV"
fi

timestamp="$(date -Iseconds)"

log_one_seed () {
  local pop="$1"
  local seed="$2"
  # Append to the global experiment TSV
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$timestamp" "$EXPR_ID" "$pop" "$seed" "$EXCLUDE_DEFAULT" "$FUNCTIONALITY" \
    "$NC" "$NR" "$ALPHA" \
    "$P_DEGRADED" "$P_DEGRADED_SD" \
    "$LEN_MEAN" "$LEN_SD" "$LEN_MIN" "$LEN_MAX" \
    "$CUT_INTACT_MEAN" "$CUT_INTACT_SD" "$CUT_MAX_INTACT" \
    "$CUT_DEGRADED_MEAN" "$CUT_DEGRADED_SD" "$CUT_MAX" \
    "$C_MAX" "$MIN_KEEP_LEN" \
    "$ART_SS" "$ART_L" "$ART_C" "$ART_M" "$ART_S" \
    "$MIXCR_PRESET" "$MIXCR_THREADS" >> "$EXPERIMENTS_TSV"
}

run_one_seed () {
  local pop="$1"
  local seed="$2"

  # Modified: Directory structure changed to sample/expr/pop/seed
  local out_dir="${EXPR_DIR}/${pop}/seed${seed}"
  mkdir -p "$out_dir"

  # Sample prefix: {POP}{_exdefaultT}_seed{seed}
  local sample_prefix="${pop}${TAG_SUFFIX}_seed${seed}"

  echo "============================================================"
  echo "[RUN] POP=${pop} seed=${seed} expr=${EXPR_ID}"
  echo "[DIR] ${out_dir}"
  echo "============================================================"

  log_one_seed "$pop" "$seed"

  pushd "$out_dir" >/dev/null

  # 1) select genotype
  python "${SCRIPT_DIR}/select_genotype.py" \
    -i "${SCRIPT_DIR}/${INPUT_REF}" \
    -p "$pop" \
    -o "." \
    -s "$seed" \
    --prefix "genotype" \
    -f "$FUNCTIONALITY" \
    -d "$EXCLUDE_DEFAULT"

  local geno_csv="genotype_${pop}${EXCLUDE_DEFAULT_SUFFIX}_seed${seed}.csv"
  [[ -f "$geno_csv" ]] || { echo "[ERROR] Missing $geno_csv"; exit 1; }

  # 2) simulate repertoire
  python "${SCRIPT_DIR}/simulate_repertoire.py" \
    --genotype "$geno_csv" \
    --dict "${SCRIPT_DIR}/${CDR3_DICT}" \
    -nc "$NC" -nr "$NR" --alpha "$ALPHA" \
    -o "${sample_prefix}" \
    -s "$seed"

  local rep_csv="${sample_prefix}_repertoire.csv"
  [[ -f "$rep_csv" ]] || { echo "[ERROR] Missing $rep_csv"; exit 1; }

  # 3) bias -> molecules fasta (+ meta)
  python "${SCRIPT_DIR}/simulate_bias_advanced.py" \
    -i "$rep_csv" \
    -o "${sample_prefix}_molecules.fasta" \
    --log "${sample_prefix}_bias_metadata.csv" \
    --leader_trim_len "$LEADER_TRIM_LEN" \
    --p_degraded "$P_DEGRADED" \
    --p_degraded_sd "$P_DEGRADED_SD" \
    --len_mean "$LEN_MEAN" \
    --len_sd "$LEN_SD" \
    --len_min "$LEN_MIN" \
    --len_max "$LEN_MAX" \
    --cut_intact_mean "$CUT_INTACT_MEAN" \
    --cut_intact_sd "$CUT_INTACT_SD" \
    --cut_max_intact "$CUT_MAX_INTACT" \
    --cut_degraded_mean "$CUT_DEGRADED_MEAN" \
    --cut_degraded_sd "$CUT_DEGRADED_SD" \
    --cut_max "$CUT_MAX" \
    --c_max "$C_MAX" \
    --min_keep_len "$MIN_KEEP_LEN" \
    --seed "$seed"

  # 4) ART
  art_illumina \
    -ss "$ART_SS" \
    -i "${sample_prefix}_molecules.fasta" \
    -p \
    -l "$ART_L" \
    -c "$ART_C" \
    -m "$ART_M" \
    -s "$ART_S" \
    -rs "$seed" \
    -o "${sample_prefix}_"

  # 5) MiXCR analyze
  mixcr analyze -s "$MIXCR_SPECIES" --threads "$MIXCR_THREADS" "$MIXCR_PRESET" \
    --force-overwrite \
    "${sample_prefix}_1.fq" "${sample_prefix}_2.fq" \
    "$sample_prefix"

  # 6) Assemble contigs
  mixcr assembleContigs \
    --force-overwrite \
    --assemble-contigs-by VDJRegion \
    "${sample_prefix}.clna" \
    "${sample_prefix}.contigs.VDJRegion.clns"

  # 7) FindAlleles
  mixcr findAlleles \
    --force-overwrite \
    --no-clns-output \
    --report "${sample_prefix}.findAlleles.report.txt" \
    --export-library "${sample_prefix}.customAlleles.json" \
    --export-alleles-mutations "${sample_prefix}.alleles.tsv" \
    "${sample_prefix}.contigs.VDJRegion.clns"

  # 8) Cleanup
  if [[ "$DO_CLEANUP" -eq 1 ]]; then
    rm -f ./*.vdjca ./*.clna ./*.clns || true
  fi

  popd >/dev/null
  echo "[DONE] pop=${pop} seed=${seed} -> ${out_dir}"
}

# =========================
# Main Loop: Iterate POPS -> SEEDS
# =========================
# Modified: Split POPS string by space/comma and iterate
for pop_item in ${POPS//,/ }; do
  for ((i=0; i<N; i++)); do
    seed=$((SEED_START + i))
    run_one_seed "$pop_item" "$seed"
  done
done

echo "============================================================"
echo "[ALL DONE] POPS=${POPS} seed_start=${SEED_START} n=${N} expr=${EXPR_ID}"
echo "[CONFIG] ${EXPERIMENTS_TSV}"
echo "============================================================"