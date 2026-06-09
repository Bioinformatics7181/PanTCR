#!/usr/bin/env bash
set -euo pipefail

# Population-agnostic leave-allele-out allele-specific simulation.
# All generated files are written under RUN_ROOT, which defaults to this
# self-contained experiment directory.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

SIMU_ROOT="${PANTCR_SIM_REF_ROOT:-${PACKAGE_ROOT}/data/ref}"
RUN_ROOT="${SCRIPT_DIR}/runs"
EXPR="expr_leaveout_allele_specific"
CODE_EXPERIMENT_ROOT="${PACKAGE_ROOT}/code/experiments/03_leave_target_allele_out"
SIMULATION_CODE_ROOT="${PACKAGE_ROOT}/code/simulation"
PANTCR_CODE_ROOT="${PACKAGE_ROOT}/code/pantcr"
BASELINE_CODE_ROOT="${PACKAGE_ROOT}/code/experiments/00_benchmark_utils"
EVALUATION_CODE_ROOT="${PACKAGE_ROOT}/code/experiments/00_benchmark_utils"

TARGET_ALLELES="TRBV7-7*03 TRBV6-6*02 TRBV12-2*02 TRBV11-2*02 TRBV24-1*02 TRBV30*02 TRBV5-1*02 TRBV5-4*05 TRBV5-5*02 TRBV6-1*03 TRBV6-6*08 TRBV7-9*09 TRBV7-9*10 TRBV12-2*03 TRBV12-3*05 TRBV12-5*02"

N_TRAIN=40
REPLICATES=5
TRAIN_POP="TRAIN"
TEST_POP="PANEL"
SEED_START_TRAIN=0
SEED_START_TEST=10000
RANDOM_SEED=20260521
SAMPLING_MODE="pooled"

PYTHON="python"
MIXCR="mixcr"
ART="art_illumina"
THREADS=8

GENES="V"
METHODS="MiXCR FindAlleles PanTCRLeaveout BayesNoPrior"
MIN_NAIVE_GRAPH=1
MIN_NAIVE_INFER=1
PENALTY_K=0
PI_MIN=0.1
SUMMARY_MIN_NAIVE=2
CROSS_GENE_DELTA=0

NC=10000
NR=50000
ALPHA=1.5

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

ART_SS="HS25"
ART_L=150
ART_C=1
ART_M=320
ART_S=5

MIXCR_SPECIES="hsa"
MIXCR_PRESET="rna-seq"

PREPARE_ONLY=0
SKIP_PREPARE=0
SKIP_SIMULATION=0
SKIP_DOWNSTREAM=0
SKIP_QC_EVAL=0
RESUME=0
DRY_RUN=0
DO_CLEANUP=1

usage() {
  cat <<EOF
Usage:
  bash $0 [options]

Main options:
  --simu-root STR           Simulation reference directory containing the package reference files (default: ${SIMU_ROOT})
  --run-root STR            Generated output root (default: ${RUN_ROOT})
  --expr STR                Experiment id (default: ${EXPR})
  --target-alleles STR      Space/comma separated target allele list
  --target-file STR         Optional target allele text file
  --n-train INT             Number of target-free training samples (default: ${N_TRAIN})
  --replicates INT          Replicates per target allele (default: ${REPLICATES})
  --sampling-mode STR       pooled or uniform (default: ${SAMPLING_MODE})
  --seed-start-train INT    First training seed (default: ${SEED_START_TRAIN})
  --seed-start-test INT     First forced-test seed (default: ${SEED_START_TEST})
  --random-seed INT         RNG seed for genotype generation (default: ${RANDOM_SEED})

Simulation options:
  --nc INT                  Number of clones (default: ${NC})
  --nr INT                  Number of reads (default: ${NR})
  --alpha FLOAT             Zipf alpha (default: ${ALPHA})
  --p-degraded FLOAT        Degraded molecule probability (default: ${P_DEGRADED})
  --p-degraded-sd FLOAT     Degraded probability SD (default: ${P_DEGRADED_SD})
  --threads INT             MiXCR threads (default: ${THREADS})

Pipeline options:
  --python STR              Python executable (default: ${PYTHON})
  --mixcr STR               MiXCR executable (default: ${MIXCR})
  --art STR                 ART executable (default: ${ART})
  --min-naive-graph INT     Graph min naive threshold (default: ${MIN_NAIVE_GRAPH})
  --min-naive-infer INT     Inference min naive threshold (default: ${MIN_NAIVE_INFER})
  --penalty-k FLOAT         Synthetic-run K=2 penalty passed to PanTCR inference (default: ${PENALTY_K})
  --pi-min FLOAT            Synthetic-run minimum mixture proportion passed to PanTCR inference (default: ${PI_MIN})
  --summary-min-naive INT   Mutation-evidence filter for target-level summary (default: ${SUMMARY_MIN_NAIVE})
  --cross-gene-delta FLOAT  Cross-gene best-match identity margin for summary (default: ${CROSS_GENE_DELTA})
  --methods STR             Summary/eval methods (default: "${METHODS}")
  --genes STR               Gene segments to process: V, J, or V,J (default: ${GENES})

Control flags:
  --prepare-only            Stop after writing genotypes and manifests
  --skip-prepare            Reuse existing manifests
  --skip-simulation         Reuse existing MiXCR/FindAlleles sample outputs
  --skip-downstream         Stop before mutation extraction, graph, inference, and summary
  --skip-qc-eval            Skip whole-repertoire QC evaluation
  --resume                  Skip sample simulation if .alleles.tsv already exists
  --no-cleanup              Keep MiXCR intermediate .vdjca/.clna/.clns files
  --dry-run                 Print commands without running them
EOF
}

TARGET_FILE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --simu-root) SIMU_ROOT="$2"; shift 2 ;;
    --run-root) RUN_ROOT="$2"; shift 2 ;;
    --expr) EXPR="$2"; shift 2 ;;
    --target-alleles) TARGET_ALLELES="$2"; shift 2 ;;
    --target-file) TARGET_FILE="$2"; shift 2 ;;
    --n-train) N_TRAIN="$2"; shift 2 ;;
    --replicates) REPLICATES="$2"; shift 2 ;;
    --sampling-mode) SAMPLING_MODE="$2"; shift 2 ;;
    --seed-start-train) SEED_START_TRAIN="$2"; shift 2 ;;
    --seed-start-test) SEED_START_TEST="$2"; shift 2 ;;
    --random-seed) RANDOM_SEED="$2"; shift 2 ;;
    --nc) NC="$2"; shift 2 ;;
    --nr) NR="$2"; shift 2 ;;
    --alpha) ALPHA="$2"; shift 2 ;;
    --p-degraded) P_DEGRADED="$2"; shift 2 ;;
    --p-degraded-sd) P_DEGRADED_SD="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --python) PYTHON="$2"; shift 2 ;;
    --mixcr) MIXCR="$2"; shift 2 ;;
    --art) ART="$2"; shift 2 ;;
    --min-naive-graph) MIN_NAIVE_GRAPH="$2"; shift 2 ;;
    --min-naive-infer) MIN_NAIVE_INFER="$2"; shift 2 ;;
    --penalty-k) PENALTY_K="$2"; shift 2 ;;
    --pi-min) PI_MIN="$2"; shift 2 ;;
    --summary-min-naive) SUMMARY_MIN_NAIVE="$2"; shift 2 ;;
    --cross-gene-delta) CROSS_GENE_DELTA="$2"; shift 2 ;;
    --methods) METHODS="$2"; shift 2 ;;
    --genes) GENES="$2"; shift 2 ;;
    --prepare-only) PREPARE_ONLY=1; shift 1 ;;
    --skip-prepare) SKIP_PREPARE=1; shift 1 ;;
    --skip-simulation) SKIP_SIMULATION=1; shift 1 ;;
    --skip-downstream) SKIP_DOWNSTREAM=1; shift 1 ;;
    --skip-qc-eval) SKIP_QC_EVAL=1; shift 1 ;;
    --resume) RESUME=1; shift 1 ;;
    --no-cleanup) DO_CLEANUP=0; shift 1 ;;
    --dry-run) DRY_RUN=1; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1"; usage; exit 1 ;;
  esac
done

run_cmd() {
  echo "+ $*"
  if [[ "${DRY_RUN}" -eq 0 ]]; then
    "$@"
  fi
}

abs_path() {
  local p="$1"
  if [[ "$p" = /* ]]; then
    readlink -f "$p"
  else
    readlink -f "${SCRIPT_DIR}/${p}"
  fi
}

find_clones_file() {
  local dir="$1"
  local f
  f=$(find "$dir" -maxdepth 2 -type f -name "*clones_TRB.tsv" 2>/dev/null | head -n 1 || true)
  if [[ -n "${f}" ]]; then echo "${f}"; return 0; fi
  f=$(find "$dir" -maxdepth 2 -type f -name "*.clones.tsv" 2>/dev/null | head -n 1 || true)
  if [[ -n "${f}" ]]; then echo "${f}"; return 0; fi
  f=$(find "$dir" -maxdepth 2 -type f -name "*clones*.tsv" 2>/dev/null | head -n 1 || true)
  if [[ -n "${f}" ]]; then echo "${f}"; return 0; fi
  echo ""
}

find_findalleles_files() {
  local dir="$1"
  local tsv json
  tsv=$(find "$dir" -maxdepth 1 -type f -name "*.alleles.tsv" 2>/dev/null | head -n 1 || true)
  json=$(find "$dir" -maxdepth 1 -type f -name "*.customAlleles.json" 2>/dev/null | head -n 1 || true)
  if [[ -n "${tsv}" && -n "${json}" ]]; then
    echo "${tsv}|${json}"
  else
    echo ""
  fi
}

copy_ref_if_needed() {
  local src="$1"
  local dst="$2"
  if [[ ! -f "${dst}" ]]; then
    run_cmd cp "${src}" "${dst}"
  fi
}

ref_for_gene() {
  local gene_upper="$1"
  if [[ "${gene_upper}" == "J" ]]; then
    echo "${IMGT_TRBJ_PRO}"
  elif [[ "${gene_upper}" == "V" ]]; then
    echo "${IMGT_TRBV_PRO}"
  else
    echo "[ERROR] Unsupported gene segment in --genes: ${gene_upper}" >&2
    return 1
  fi
}

SIMU_ROOT="$(abs_path "${SIMU_ROOT}")"
RUN_ROOT="$(abs_path "${RUN_ROOT}")"
if [[ -d "${SIMU_ROOT}/ref" ]]; then
  SIM_REF_ROOT="${SIMU_ROOT}/ref"
else
  SIM_REF_ROOT="${SIMU_ROOT}"
fi

mkdir -p "${RUN_ROOT}/logs/${EXPR}" "${RUN_ROOT}/ref"
LOG_FILE="${RUN_ROOT}/logs/${EXPR}/run_leaveout_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "======================================================="
echo " Leave-allele-out allele-specific simulation"
echo " SIM_REF_ROOT(read-only): ${SIM_REF_ROOT}"
echo " RUN_ROOT(generated):  ${RUN_ROOT}"
echo " EXPR:                 ${EXPR}"
echo " TRAIN/TEST:           ${N_TRAIN} / target x ${REPLICATES}"
echo " LOG:                  ${LOG_FILE}"
echo "======================================================="

if [[ ! -d "${SIMU_ROOT}" ]]; then
  echo "[ERROR] SIMU_ROOT not found: ${SIMU_ROOT}"
  exit 1
fi

PMTR="${RUN_ROOT}/ref/pmTR_TRB_V_J_cleaned.csv"
IMGT_DEFAULT="${RUN_ROOT}/ref/IMGT_TRB_default.csv"
TRB_INDEX="${RUN_ROOT}/ref/TRB_index.csv"
IMGT_TRBV_PRO="${RUN_ROOT}/ref/IMGT_TRBV_pro.tsv"
IMGT_TRBJ_PRO="${RUN_ROOT}/ref/IMGT_TRBJ_pro.tsv"
CDR3_DICT="${RUN_ROOT}/ref/cdr3_dict.pkl"

copy_ref_if_needed "${SIM_REF_ROOT}/pmTR_TRB_V_J_cleaned.csv" "${PMTR}"
copy_ref_if_needed "${SIM_REF_ROOT}/IMGT_TRB_default.csv" "${IMGT_DEFAULT}"
copy_ref_if_needed "${SIM_REF_ROOT}/TRB_index.csv" "${TRB_INDEX}"
copy_ref_if_needed "${SIM_REF_ROOT}/IMGT_TRBV_pro.tsv" "${IMGT_TRBV_PRO}"
copy_ref_if_needed "${SIM_REF_ROOT}/IMGT_TRBJ_pro.tsv" "${IMGT_TRBJ_PRO}"

if [[ ! -f "${CDR3_DICT}" ]]; then
  if [[ -f "${SIM_REF_ROOT}/cdr3_dict.pkl" ]]; then
    run_cmd cp "${SIM_REF_ROOT}/cdr3_dict.pkl" "${CDR3_DICT}"
  else
    echo "[ERROR] Missing ${SIM_REF_ROOT}/cdr3_dict.pkl."
    echo "[ERROR] Use the unzipped cdr3_dict.pkl directly; do not rely on decompressing cdr3_dict.zip during the job."
    exit 1
  fi
fi

MANIFEST_DIR="${RUN_ROOT}/results/leaveout/${EXPR}"
SIM_MANIFEST="${MANIFEST_DIR}/simulation_manifest.tsv"
TARGET_MANIFEST="${MANIFEST_DIR}/target_panel_manifest.tsv"
TARGET_METADATA="${MANIFEST_DIR}/target_metadata.csv"

if [[ "${SKIP_PREPARE}" -eq 0 ]]; then
  PREP_CMD=("${PYTHON}" "${CODE_EXPERIMENT_ROOT}/prepare_leaveout_experiment.py"
    --simu-root "${SIM_REF_ROOT}"
    --run-root "${RUN_ROOT}"
    --expr "${EXPR}"
    --target-alleles "${TARGET_ALLELES}"
    --n-train "${N_TRAIN}"
    --replicates "${REPLICATES}"
    --seed-start-train "${SEED_START_TRAIN}"
    --seed-start-test "${SEED_START_TEST}"
    --train-pop "${TRAIN_POP}"
    --test-pop "${TEST_POP}"
    --sampling-mode "${SAMPLING_MODE}"
    --random-seed "${RANDOM_SEED}"
    --pmtr "${PMTR}"
    --pmtr-ref "${PMTR}"
    --default-ref "${IMGT_DEFAULT}"
    --index "${TRB_INDEX}")
  if [[ -n "${TARGET_FILE}" ]]; then
    PREP_CMD+=(--target-file "${TARGET_FILE}")
  fi
  run_cmd "${PREP_CMD[@]}"
fi

if [[ "${PREPARE_ONLY}" -eq 1 ]]; then
  echo "[DONE] prepare-only mode. Manifest: ${SIM_MANIFEST}"
  exit 0
fi

if [[ ! -f "${SIM_MANIFEST}" ]]; then
  echo "[ERROR] Missing simulation manifest: ${SIM_MANIFEST}"
  exit 1
fi

simulate_one_sample() {
  local role="$1"
  local target_allele="$2"
  local target_gene="$3"
  local seed="$4"
  local replicate="$5"
  local pop="$6"
  local genotype_csv="$7"
  local sample_dir="$8"
  local sample_prefix="$9"

  mkdir -p "${sample_dir}"
  if [[ "${RESUME}" -eq 1 && -f "${sample_dir}/${sample_prefix}.alleles.tsv" ]]; then
    echo "[SKIP] resume: ${sample_prefix} already has FindAlleles output."
    return
  fi

  echo "-------------------------------------------------------"
  echo "[SIM] role=${role} target=${target_allele} gene=${target_gene} seed=${seed} rep=${replicate}"
  echo "[DIR] ${sample_dir}"

  pushd "${sample_dir}" >/dev/null
  run_cmd "${PYTHON}" "${SIMULATION_CODE_ROOT}/simulate_repertoire.py" \
    --genotype "${genotype_csv}" \
    --dict "${CDR3_DICT}" \
    -nc "${NC}" \
    -nr "${NR}" \
    --alpha "${ALPHA}" \
    -o "${sample_prefix}" \
    -s "${seed}"

  run_cmd "${PYTHON}" "${SIMULATION_CODE_ROOT}/simulate_bias_advanced.py" \
    -i "${sample_prefix}_repertoire.csv" \
    -o "${sample_prefix}_molecules.fasta" \
    --log "${sample_prefix}_bias_metadata.csv" \
    --leader_trim_len "${LEADER_TRIM_LEN}" \
    --p_degraded "${P_DEGRADED}" \
    --p_degraded_sd "${P_DEGRADED_SD}" \
    --len_mean "${LEN_MEAN}" \
    --len_sd "${LEN_SD}" \
    --len_min "${LEN_MIN}" \
    --len_max "${LEN_MAX}" \
    --cut_intact_mean "${CUT_INTACT_MEAN}" \
    --cut_intact_sd "${CUT_INTACT_SD}" \
    --cut_max_intact "${CUT_MAX_INTACT}" \
    --cut_degraded_mean "${CUT_DEGRADED_MEAN}" \
    --cut_degraded_sd "${CUT_DEGRADED_SD}" \
    --cut_max "${CUT_MAX}" \
    --c_max "${C_MAX}" \
    --min_keep_len "${MIN_KEEP_LEN}" \
    --seed "${seed}"

  run_cmd "${ART}" \
    -ss "${ART_SS}" \
    -i "${sample_prefix}_molecules.fasta" \
    -p \
    -l "${ART_L}" \
    -c "${ART_C}" \
    -m "${ART_M}" \
    -s "${ART_S}" \
    -rs "${seed}" \
    -o "${sample_prefix}_"

  run_cmd "${MIXCR}" analyze -s "${MIXCR_SPECIES}" --threads "${THREADS}" "${MIXCR_PRESET}" \
    --force-overwrite \
    "${sample_prefix}_1.fq" "${sample_prefix}_2.fq" \
    "${sample_prefix}"

  run_cmd "${MIXCR}" assembleContigs \
    --force-overwrite \
    --assemble-contigs-by VDJRegion \
    "${sample_prefix}.clna" \
    "${sample_prefix}.contigs.VDJRegion.clns"

  run_cmd "${MIXCR}" findAlleles \
    --force-overwrite \
    --no-clns-output \
    --report "${sample_prefix}.findAlleles.report.txt" \
    --export-library "${sample_prefix}.customAlleles.json" \
    --export-alleles-mutations "${sample_prefix}.alleles.tsv" \
    "${sample_prefix}.contigs.VDJRegion.clns"

  if [[ "${DO_CLEANUP}" -eq 1 && "${DRY_RUN}" -eq 0 ]]; then
    rm -f ./*.vdjca ./*.clna ./*.clns || true
  fi
  popd >/dev/null
}

if [[ "${SKIP_SIMULATION}" -eq 0 ]]; then
  tail -n +2 "${SIM_MANIFEST}" | while IFS=$'\t' read -r role target_allele target_gene seed replicate pop genotype_csv sample_dir sample_prefix; do
    [[ -n "${role}" ]] || continue
    simulate_one_sample "${role}" "${target_allele}" "${target_gene}" "${seed}" "${replicate}" "${pop}" "${genotype_csv}" "${sample_dir}" "${sample_prefix}"
  done
fi

if [[ "${SKIP_DOWNSTREAM}" -eq 1 ]]; then
  echo "[DONE] stopped before downstream steps."
  exit 0
fi

extract_mutations_for_pop() {
  local pop="$1"
  local pop_dir="${RUN_ROOT}/samples/${EXPR}/${pop}"
  if [[ ! -d "${pop_dir}" ]]; then
    echo "[WARN] Missing sample pop directory: ${pop_dir}"
    return
  fi
  for gene in ${GENES//,/ }; do
    local gene_upper
    gene_upper="$(echo "${gene}" | tr '[:lower:]' '[:upper:]')"
    local ref_file
    ref_file="$(ref_for_gene "${gene_upper}")"
    local out_dir="${RUN_ROOT}/results/mutations/${EXPR}/${pop}/${gene_upper}"
    mkdir -p "${out_dir}"
    for seed_path in "${pop_dir}"/seed*; do
      [[ -d "${seed_path}" ]] || continue
      local seed_base seed_num clones sample_tag out_prefix
      seed_base="$(basename "${seed_path}")"
      seed_num="${seed_base#seed}"
      [[ "${seed_num}" =~ ^[0-9]+$ ]] || continue
      clones="$(find_clones_file "${seed_path}")"
      if [[ -z "${clones}" ]]; then
        echo "[WARN] No clones file for ${pop} seed${seed_num}"
        continue
      fi
      sample_tag="${pop}_seed${seed_num}"
      out_prefix="${out_dir}/${sample_tag}.${gene_upper}"
      run_cmd "${PYTHON}" "${PANTCR_CODE_ROOT}/collect_mutation.py" \
        --input "${clones}" \
        --gene "${gene_upper}" \
        --ref "${ref_file}" \
        --prefix "${out_prefix}"
    done
  done
}

copy_labels_for_pop() {
  local pop="$1"
  local out_dir="${RUN_ROOT}/results/labels/${EXPR}/${pop}"
  mkdir -p "${out_dir}"
  awk -F '\t' -v pop="${pop}" 'NR > 1 && $6 == pop {print $4 "\t" $7}' "${SIM_MANIFEST}" | while IFS=$'\t' read -r seed genotype_csv; do
    [[ -n "${seed}" ]] || continue
    run_cmd cp "${genotype_csv}" "${out_dir}/genotype_${pop}_seed${seed}.csv"
  done
}

run_mixcr_baseline() {
  local pop="${TEST_POP}"
  for gene in ${GENES//,/ }; do
    local gene_upper
    gene_upper="$(echo "${gene}" | tr '[:lower:]' '[:upper:]')"
    local out_dir="${RUN_ROOT}/results/infer/${EXPR}/MiXCR/${pop}/${gene_upper}"
    mkdir -p "${out_dir}"
    for seed_path in "${RUN_ROOT}/samples/${EXPR}/${pop}"/seed*; do
      [[ -d "${seed_path}" ]] || continue
      local seed_num clones output_path
      seed_num="$(basename "${seed_path}")"
      seed_num="${seed_num#seed}"
      clones="$(find_clones_file "${seed_path}")"
      if [[ -z "${clones}" ]]; then
        echo "[WARN] No clones file for MiXCR baseline seed${seed_num}"
        continue
      fi
      output_path="${out_dir}/infer_${pop}_seed${seed_num}.${gene_upper}.csv"
      run_cmd "${PYTHON}" "${BASELINE_CODE_ROOT}/clones2alleles.py" \
        --clones "${clones}" \
        --gene "${gene_upper}" \
        --ref "${IMGT_DEFAULT}" \
        --output "${output_path}"
    done
  done
}

run_findalleles_baseline() {
  local pop="${TEST_POP}"
  for gene in ${GENES//,/ }; do
    local gene_upper ref_file
    gene_upper="$(echo "${gene}" | tr '[:lower:]' '[:upper:]')"
    ref_file="$(ref_for_gene "${gene_upper}")"
    local out_dir="${RUN_ROOT}/results/infer/${EXPR}/FindAlleles/${pop}/${gene_upper}"
    mkdir -p "${out_dir}"
    for seed_path in "${RUN_ROOT}/samples/${EXPR}/${pop}"/seed*; do
      [[ -d "${seed_path}" ]] || continue
      local seed_num files input_tsv input_json output_path
      seed_num="$(basename "${seed_path}")"
      seed_num="${seed_num#seed}"
      files="$(find_findalleles_files "${seed_path}")"
      if [[ -z "${files}" ]]; then
        echo "[WARN] No FindAlleles files for seed${seed_num}"
        continue
      fi
      IFS='|' read -r input_tsv input_json <<< "${files}"
      output_path="${out_dir}/infer_${pop}_seed${seed_num}.${gene_upper}.csv"
      run_cmd "${PYTHON}" "${BASELINE_CODE_ROOT}/infer_for_findalleles.py" \
        --tsv "${input_tsv}" \
        --json "${input_json}" \
        --gene "${gene_upper}" \
        --ref "${ref_file}" \
        --output "${output_path}"
    done
  done
}

build_leaveout_graph() {
  for gene in ${GENES//,/ }; do
    local gene_upper train_mut_dir meta graph_dir
    gene_upper="$(echo "${gene}" | tr '[:lower:]' '[:upper:]')"
    train_mut_dir="${RUN_ROOT}/results/mutations/${EXPR}/${TRAIN_POP}/${gene_upper}"
    meta="${RUN_ROOT}/results/leaveout/${EXPR}/leaveout_graph_metadata_${gene_upper}.csv"
    graph_dir="${RUN_ROOT}/results/pang/${EXPR}/LeaveoutGlobal/${gene_upper}"
    mkdir -p "$(dirname "${meta}")" "${graph_dir}"
    printf "sample_id,population_id,filename\n" > "${meta}"
    for sample_csv in "${train_mut_dir}"/*.csv; do
      [[ -f "${sample_csv}" ]] || continue
      local base sample_id
      base="$(basename "${sample_csv}")"
      sample_id="${base%.${gene_upper}_sequences.csv}"
      printf "%s,GLOBAL,%s\n" "${sample_id}" "${base}" >> "${meta}"
    done
    run_cmd "${PYTHON}" "${PANTCR_CODE_ROOT}/build_pangenome_graph.py" \
      --in_dir "${train_mut_dir}" \
      --metadata "${meta}" \
      --out_dir "${graph_dir}" \
      --min_naive "${MIN_NAIVE_GRAPH}" \
      --write_paths
  done
}

run_pantcr_inference() {
  local pop="${TEST_POP}"
  for gene in ${GENES//,/ }; do
    local gene_upper mut_dir graph_dir graph_out np_out
    gene_upper="$(echo "${gene}" | tr '[:lower:]' '[:upper:]')"
    mut_dir="${RUN_ROOT}/results/mutations/${EXPR}/${pop}/${gene_upper}"
    graph_dir="${RUN_ROOT}/results/pang/${EXPR}/LeaveoutGlobal/${gene_upper}"
    graph_out="${RUN_ROOT}/results/infer/${EXPR}/PanTCRLeaveout/${pop}/${gene_upper}"
    np_out="${RUN_ROOT}/results/infer/${EXPR}/BayesNoPrior/${pop}/${gene_upper}"
    mkdir -p "${graph_out}" "${np_out}"
    for sample_csv in "${mut_dir}"/*.csv; do
      [[ -f "${sample_csv}" ]] || continue
      local base seed_id
      base="$(basename "${sample_csv}")"
      if [[ "${base}" =~ seed([0-9]+) ]]; then
        seed_id="seed${BASH_REMATCH[1]}"
      else
        echo "[WARN] Could not parse seed from ${base}"
        continue
      fi
      run_cmd "${PYTHON}" "${PANTCR_CODE_ROOT}/infer_genotype_bayes.py" \
        --sample_csv "${sample_csv}" \
        --min_naive "${MIN_NAIVE_INFER}" \
        --penalty_K "${PENALTY_K}" \
        --pi_min "${PI_MIN}" \
        --pangenome_dir "${graph_dir}" \
        --out "${graph_out}/infer_${pop}_${seed_id}.${gene_upper}.csv"
      run_cmd "${PYTHON}" "${PANTCR_CODE_ROOT}/infer_genotype_bayes.py" \
        --sample_csv "${sample_csv}" \
        --min_naive "${MIN_NAIVE_INFER}" \
        --penalty_K "${PENALTY_K}" \
        --pi_min "${PI_MIN}" \
        --out "${np_out}/infer_${pop}_${seed_id}.${gene_upper}.csv"
    done
  done
}

run_qc_eval() {
  local pop="${TEST_POP}"
  for gene in ${GENES//,/ }; do
    local gene_upper
    gene_upper="$(echo "${gene}" | tr '[:lower:]' '[:upper:]')"
    for method in ${METHODS//,/ }; do
      local infer_dir="${RUN_ROOT}/results/infer/${EXPR}/${method}/${pop}/${gene_upper}"
      [[ -d "${infer_dir}" ]] || continue
      local out_dir="${RUN_ROOT}/results/eval/${EXPR}/${method}/${pop}/${gene_upper}"
      mkdir -p "${out_dir}"
      for infer_file in "${infer_dir}"/*.csv; do
        [[ -f "${infer_file}" ]] || continue
        local filename seed_id label_path out_prefix
        filename="$(basename "${infer_file}")"
        if [[ "${filename}" =~ seed([0-9]+) ]]; then
          seed_id="${BASH_REMATCH[1]}"
        else
          continue
        fi
        label_path="${RUN_ROOT}/results/labels/${EXPR}/${pop}/genotype_${pop}_seed${seed_id}.csv"
        if [[ ! -f "${label_path}" ]]; then
          echo "[WARN] Missing label for ${method} seed${seed_id}: ${label_path}"
          continue
        fi
        out_prefix="${out_dir}/eval_${method}_${pop}_seed${seed_id}_${gene_upper}"
        run_cmd "${PYTHON}" "${EVALUATION_CODE_ROOT}/evaluate_allele_calls.py" \
          --gt "${label_path}" \
          --infer "${infer_file}" \
          --index "${TRB_INDEX}" \
          --gene_type "${gene_upper}" \
          --out_prefix "${out_prefix}"
      done
    done
  done
}

extract_mutations_for_pop "${TRAIN_POP}"
extract_mutations_for_pop "${TEST_POP}"
copy_labels_for_pop "${TRAIN_POP}"
copy_labels_for_pop "${TEST_POP}"
run_mixcr_baseline
run_findalleles_baseline
build_leaveout_graph
run_pantcr_inference

if [[ "${SKIP_QC_EVAL}" -eq 0 ]]; then
  run_qc_eval
fi

run_cmd "${PYTHON}" "${CODE_EXPERIMENT_ROOT}/summarize_leaveout_target_recovery.py" \
  --run-root "${RUN_ROOT}" \
  --expr "${EXPR}" \
  --test-pop "${TEST_POP}" \
  --methods "${METHODS}" \
  --target-manifest "${TARGET_MANIFEST}" \
  --target-metadata "${TARGET_METADATA}" \
  --pmtr-ref "${PMTR}" \
  --index "${TRB_INDEX}" \
  --out-dir "${RUN_ROOT}/results/leaveout/${EXPR}/summary" \
  --min-naive "${SUMMARY_MIN_NAIVE}" \
  --cross-gene-delta "${CROSS_GENE_DELTA}"

run_cmd "${PYTHON}" "${CODE_EXPERIMENT_ROOT}/summarize_leaveout_allele_recovery.py" \
  --summary-dir "${RUN_ROOT}/results/leaveout/${EXPR}/summary" \
  --metadata-csv "${TARGET_METADATA}" \
  --out-dir "${RUN_ROOT}/results/leaveout/${EXPR}/summary"

echo "======================================================="
echo " Done."
echo " Target-level summary:"
echo " ${RUN_ROOT}/results/leaveout/${EXPR}/summary/target_recovery_summary.csv"
echo " Allele-level recovery tables:"
echo " ${RUN_ROOT}/results/leaveout/${EXPR}/summary/leaveout_allele_recovery_overall.csv"
echo " ${RUN_ROOT}/results/leaveout/${EXPR}/summary/leaveout_allele_recovery_by_target.csv"
echo " Per-sample status:"
echo " ${RUN_ROOT}/results/leaveout/${EXPR}/summary/per_target_method_status.csv"
echo "======================================================="
