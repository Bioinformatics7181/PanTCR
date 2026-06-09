#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

OUT_BASE="${OUT_BASE:-${SCRIPT_DIR}/generated/per_dataset_results}"
INPUT_ROOT="${INPUT_ROOT:-${SCRIPT_DIR}/input}"
THREADS="${THREADS:-8}"
MIN_NAIVE="${MIN_NAIVE:-2}"
PENALTY_K="${PENALTY_K:-3}"
PI_MIN="${PI_MIN:-0.1}"

DATASETS="${DATASETS:-SC-1:dataset2 SC-2:dataset5 SC-3:dataset4 SC-4:dataset3 SC-5:dataset10 SC-6:dataset9 SC-7:dataset8 SC-8:dataset1}"

usage() {
  cat <<EOF
Usage:
  bash run_all_datasets.sh [options]

Options:
  --datasets STR     Space-separated SC:dataset pairs.
                     Default: "${DATASETS}"
  --input-root DIR   Input root. Default: ${INPUT_ROOT}
  --out-base DIR     Output root. Default: ${OUT_BASE}
  --threads INT      External-tool threads. Default: ${THREADS}
  --min-naive INT    PanTCR min_naive. Default: ${MIN_NAIVE}
  --penalty-k INT    PanTCR penalty_K. Default: ${PENALTY_K}
  -h, --help         Show help.

For each SC:dataset pair, place inputs under:
  input/SC-X_datasetY/

Expected files are genotype.csv plus either trb_all.clones_TRB.tsv or
trb_all.bam. Optional FindAlleles files are trb_all.alleles.tsv and
trb_all.customAlleles.json.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --datasets) DATASETS="$2"; shift 2 ;;
    --input-root) INPUT_ROOT="$2"; shift 2 ;;
    --out-base) OUT_BASE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --min-naive) MIN_NAIVE="$2"; shift 2 ;;
    --penalty-k) PENALTY_K="$2"; shift 2 ;;
    --pi-min) PI_MIN="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

mkdir -p "${OUT_BASE}"
manifest="${OUT_BASE}/scbulk_manuscript_mode_mapping.csv"
printf "SC,Dataset,ResultDirectory\n" > "${manifest}"

for pair in ${DATASETS}; do
  sc_id="${pair%%:*}"
  dataset_id="${pair#*:}"
  input_dir="${INPUT_ROOT}/${sc_id}_${dataset_id}"
  genotype="${input_dir}/genotype.csv"
  clones="${input_dir}/trb_all.clones_TRB.tsv"
  gex_bam="${input_dir}/trb_all.bam"
  alleles_tsv="${input_dir}/trb_all.alleles.tsv"
  alleles_json="${input_dir}/trb_all.customAlleles.json"

  if [[ ! -f "${genotype}" ]]; then
    echo "[ERROR] Missing genotype CSV for ${sc_id}_${dataset_id}: ${genotype}" >&2
    exit 1
  fi

  args=(
    bash "${SCRIPT_DIR}/run_one_dataset.sh"
    --sc-id "${sc_id}"
    --dataset-id "${dataset_id}"
    --genotype "${genotype}"
    --out-base "${OUT_BASE}"
    --threads "${THREADS}"
    --min-naive "${MIN_NAIVE}"
    --penalty-k "${PENALTY_K}"
    --pi-min "${PI_MIN}"
  )
  if [[ -f "${clones}" ]]; then
    args+=(--clones "${clones}")
  elif [[ -f "${gex_bam}" ]]; then
    args+=(--gex-bam "${gex_bam}")
  else
    echo "[ERROR] Missing clones TSV or GEX BAM for ${sc_id}_${dataset_id} under ${input_dir}" >&2
    exit 1
  fi
  if [[ -f "${alleles_tsv}" && -f "${alleles_json}" ]]; then
    args+=(--alleles-tsv "${alleles_tsv}" --alleles-json "${alleles_json}")
  fi

  echo "[RUN] ${sc_id}_${dataset_id}"
  "${args[@]}"
  printf "%s,%s,%s\n" "${sc_id}" "${dataset_id}" "${sc_id}_${dataset_id}" >> "${manifest}"
done

echo "[DONE] pseudo-bulk dataset set written to ${OUT_BASE}"
echo "[DONE] mapping manifest written to ${manifest}"
