#!/usr/bin/env bash
set -euo pipefail

# Run the MiXCR-all/custom-library baseline on already generated synthetic
# samples, then evaluate allele-name hits with the unified evaluator.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

EXPR_ID="expr_ScenarioA"
POPS="AFR EUR AMR EAS SAS"
GENES="V J"
SAMPLE_ROOT="${SCRIPT_DIR}/samples"
RESULTS_ROOT="${SCRIPT_DIR}/results"
LIBRARY_JSON="${PACKAGE_ROOT}/data/ref/hsa_custom_library.json"
SPECIES_NAME="hsa_custom"
MIXCR_PRESET="rna-seq"
MIXCR_THREADS=8
PYTHON="${PYTHON:-python}"
EVAL_PY="${PACKAGE_ROOT}/code/experiments/00_benchmark_utils/evaluate_allele_calls.py"
FORCE=1
CLEANUP=1

usage() {
  cat <<EOF
Usage:
  bash $0 --expr expr_ScenarioA [options]

Options:
  --expr STR           Experiment ID. Default: ${EXPR_ID}
  --pop STR            Population list. Default: "${POPS}"
  --gene STR           Gene list. Default: "${GENES}"
  --sample-root STR    Synthetic sample root. Default: ${SAMPLE_ROOT}
  --results-root STR   Result root. Default: ${RESULTS_ROOT}
  --library STR        MiXCR custom library JSON. Default: data/ref/hsa_custom_library.json
  --species STR        MiXCR species/library name. Default: ${SPECIES_NAME}
  --mixcr-preset STR   MiXCR analyze preset. Default: ${MIXCR_PRESET}
  --mixcr-threads INT  MiXCR threads. Default: ${MIXCR_THREADS}
  --python STR         Python executable. Default: ${PYTHON}
  --no-force           Do not overwrite existing MiXCR outputs.
  --no-cleanup         Keep MiXCR intermediate vdjca/cln files.
  -h, --help           Show help.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --expr) EXPR_ID="$2"; shift 2 ;;
    --pop) POPS="$2"; shift 2 ;;
    --gene) GENES="$2"; shift 2 ;;
    --sample-root) SAMPLE_ROOT="$2"; shift 2 ;;
    --results-root) RESULTS_ROOT="$2"; shift 2 ;;
    --library) LIBRARY_JSON="$2"; shift 2 ;;
    --species) SPECIES_NAME="$2"; shift 2 ;;
    --mixcr-preset) MIXCR_PRESET="$2"; shift 2 ;;
    --mixcr-threads) MIXCR_THREADS="$2"; shift 2 ;;
    --python) PYTHON="$2"; shift 2 ;;
    --no-force) FORCE=0; shift 1 ;;
    --no-cleanup) CLEANUP=0; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1"; usage; exit 1 ;;
  esac
done

resolve_path() {
  local p="$1"
  if [[ "$p" = /* || "$p" =~ ^[A-Za-z]:[\\/] ]]; then
    echo "$p"
  else
    echo "${PACKAGE_ROOT}/${p}"
  fi
}

SAMPLE_ROOT_ABS="$(resolve_path "${SAMPLE_ROOT}")"
RESULTS_ROOT_ABS="$(resolve_path "${RESULTS_ROOT}")"
LIBRARY_JSON_ABS="$(resolve_path "${LIBRARY_JSON}")"

if [[ ! -f "${LIBRARY_JSON_ABS}" ]]; then
  echo "[ERROR] MiXCR custom library not found: ${LIBRARY_JSON_ABS}"
  exit 1
fi
if [[ ! -f "${EVAL_PY}" ]]; then
  echo "[ERROR] Evaluator not found: ${EVAL_PY}"
  exit 1
fi

LOG_DIR="${RESULTS_ROOT_ABS}/logs/${EXPR_ID}"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/mixcr_all_custom_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "======================================================="
echo " MiXCR-all custom-library baseline"
echo " EXPR_ID:      ${EXPR_ID}"
echo " POPS:         ${POPS}"
echo " GENES:        ${GENES}"
echo " SAMPLE_ROOT:  ${SAMPLE_ROOT_ABS}"
echo " RESULTS_ROOT: ${RESULTS_ROOT_ABS}"
echo " LIBRARY:      ${LIBRARY_JSON_ABS}"
echo "======================================================="

find_fastq_pair() {
  local sample_dir="$1"
  local r1 r2
  r1="$(find "${sample_dir}" -maxdepth 1 -type f \( -name "*_1.fastq.gz" -o -name "*_1.fastq" -o -name "*_1.fq.gz" -o -name "*_1.fq" \) | sort | head -n 1)"
  r2="$(find "${sample_dir}" -maxdepth 1 -type f \( -name "*_2.fastq.gz" -o -name "*_2.fastq" -o -name "*_2.fq.gz" -o -name "*_2.fq" \) | sort | head -n 1)"
  if [[ -n "${r1}" && -n "${r2}" ]]; then
    printf "%s\t%s\n" "${r1}" "${r2}"
  fi
}

for pop in ${POPS//,/ }; do
  POP_SAMPLE_DIR="${SAMPLE_ROOT_ABS}/${EXPR_ID}/${pop}"
  if [[ ! -d "${POP_SAMPLE_DIR}" ]]; then
    echo "[SKIP] Missing sample population folder: ${POP_SAMPLE_DIR}"
    continue
  fi

  shopt -s nullglob
  for sample_dir in "${POP_SAMPLE_DIR}"/seed*; do
    [[ -d "${sample_dir}" ]] || continue
    seed_name="$(basename "${sample_dir}")"
    if [[ ! "${seed_name}" =~ seed([0-9]+) ]]; then
      echo "[SKIP] Cannot parse seed from ${sample_dir}"
      continue
    fi
    seed_id="${BASH_REMATCH[1]}"

    pair="$(find_fastq_pair "${sample_dir}")"
    if [[ -z "${pair}" ]]; then
      echo "[SKIP] Missing paired FASTQ in ${sample_dir}"
      continue
    fi
    r1="${pair%%$'\t'*}"
    r2="${pair#*$'\t'}"
    r1_base="$(basename "${r1}")"
    sample_prefix="${r1_base%%_1*}"

    label_match=( "${sample_dir}"/genotype_*"seed${seed_id}".csv )
    if [[ ${#label_match[@]} -eq 0 ]]; then
      echo "[SKIP] Missing genotype label for ${sample_dir}"
      continue
    fi
    genotype_csv="${label_match[0]}"

    out_sample_dir="${RESULTS_ROOT_ABS}/infer/${EXPR_ID}/MiXCR-all/${pop}/${seed_name}"
    mkdir -p "${out_sample_dir}"
    mixcr_prefix="${out_sample_dir}/${sample_prefix}_mixcr_all"
    tsv_file="${mixcr_prefix}.clones_TRB.tsv"

    echo "[RUN] ${EXPR_ID}/${pop}/${seed_name}: ${sample_prefix}"
    mixcr_cmd=(mixcr analyze -s "${SPECIES_NAME}" --library "${LIBRARY_JSON_ABS}" --threads "${MIXCR_THREADS}" "${MIXCR_PRESET}")
    if [[ "${FORCE}" -eq 1 ]]; then
      mixcr_cmd+=(--force-overwrite)
    fi
    mixcr_cmd+=("${r1}" "${r2}" "${mixcr_prefix}")
    "${mixcr_cmd[@]}"

    if [[ ! -f "${tsv_file}" ]]; then
      echo "[ERROR] Missing MiXCR-all clone TSV: ${tsv_file}"
      exit 1
    fi

    for gene in ${GENES//,/ }; do
      GENE="$(echo "${gene}" | tr '[:lower:]' '[:upper:]')"
      eval_dir="${RESULTS_ROOT_ABS}/eval/${EXPR_ID}/MiXCR-all/${pop}/${GENE}"
      mkdir -p "${eval_dir}"
      out_prefix="${eval_dir}/eval_MiXCR-all_${sample_prefix}_${GENE}"
      "${PYTHON}" "${EVAL_PY}" \
        --mode mixcr-all \
        --tsv "${tsv_file}" \
        --gt "${genotype_csv}" \
        --gene_type "${GENE}" \
        --out_prefix "${out_prefix}"
    done

    if [[ "${CLEANUP}" -eq 1 ]]; then
      rm -f "${mixcr_prefix}".vdjca "${mixcr_prefix}".clns "${mixcr_prefix}".clna || true
    fi
  done
  shopt -u nullglob
done

echo "[DONE] MiXCR-all baseline evaluated. Log: ${LOG_FILE}"
