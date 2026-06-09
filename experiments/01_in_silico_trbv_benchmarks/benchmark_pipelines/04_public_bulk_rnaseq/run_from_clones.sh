#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"

PYTHON="${PYTHON:-python}"
BULK_LABEL="BULK-1"
CLONES_DIR=""
OUT_ROOT=""
GRAPH_ROOT="${GRAPH_ROOT:-${PACKAGE_ROOT}/data/pretrained_graphs/trb_airr1}"
MIN_NAIVE="${MIN_NAIVE:-2}"
PENALTY_K="${PENALTY_K:-3}"
PI_MIN="${PI_MIN:-0.2}"
SUMMARY_OUTPUT=""

usage() {
  cat <<EOF
Usage:
  bash run_from_clones.sh --bulk-label BULK-1 --clones-dir generated/BULK-1/clones

This driver starts from MiXCR *.clones_TRB.tsv files and writes the public-bulk
mutation and PanTCR inference tree expected by experiment 05. Raw SRA download
and MiXCR clone generation are documented in README.md and should be run before
this script when clones are not already available.

Options:
  --summary-output PATH  Override candidate summary path. Defaults to
                         candidate_allele_summary.csv under the output root.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bulk-label) BULK_LABEL="$2"; shift 2 ;;
    --clones-dir) CLONES_DIR="$2"; shift 2 ;;
    --out-root) OUT_ROOT="$2"; shift 2 ;;
    --graph-root) GRAPH_ROOT="$2"; shift 2 ;;
    --min-naive) MIN_NAIVE="$2"; shift 2 ;;
    --penalty-k) PENALTY_K="$2"; shift 2 ;;
    --pi-min) PI_MIN="$2"; shift 2 ;;
    --summary-output) SUMMARY_OUTPUT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "${CLONES_DIR}" ]]; then
  echo "[ERROR] --clones-dir is required." >&2
  usage
  exit 1
fi
if [[ -z "${OUT_ROOT}" ]]; then
  OUT_ROOT="${SCRIPT_DIR}/generated/${BULK_LABEL}"
fi

mkdir -p "${OUT_ROOT}/mutations/V" "${OUT_ROOT}/mutations/J" "${OUT_ROOT}/infer/PanTCR/V" "${OUT_ROOT}/infer/PanTCR/J" "${OUT_ROOT}/pang/semi"

if [[ -d "${GRAPH_ROOT}/V" && ! -e "${OUT_ROOT}/pang/semi/V" ]]; then
  cp -R "${GRAPH_ROOT}/V" "${OUT_ROOT}/pang/semi/V"
fi
if [[ -d "${GRAPH_ROOT}/J" && ! -e "${OUT_ROOT}/pang/semi/J" ]]; then
  cp -R "${GRAPH_ROOT}/J" "${OUT_ROOT}/pang/semi/J"
fi

shopt -s nullglob
for clones in "${CLONES_DIR}"/*.clones_TRB.tsv; do
  sample_id="$(basename "${clones}" .clones_TRB.tsv)"
  echo "[RUN] ${BULK_LABEL} ${sample_id}"
  "${PYTHON}" "${PACKAGE_ROOT}/code/pantcr/collect_mutation.py" --input "${clones}" --gene V --ref "${PACKAGE_ROOT}/data/ref/IMGT_TRBV_pro.tsv" --prefix "${OUT_ROOT}/mutations/V/${sample_id}.V"
  "${PYTHON}" "${PACKAGE_ROOT}/code/pantcr/collect_mutation.py" --input "${clones}" --gene J --ref "${PACKAGE_ROOT}/data/ref/IMGT_TRBJ_pro.tsv" --prefix "${OUT_ROOT}/mutations/J/${sample_id}.J"
  "${PYTHON}" "${PACKAGE_ROOT}/code/pantcr/infer_genotype_bayes.py" --sample_csv "${OUT_ROOT}/mutations/V/${sample_id}.V_sequences.csv" --out "${OUT_ROOT}/infer/PanTCR/V/infer_${sample_id}.V.csv" --pangenome_dir "${GRAPH_ROOT}/V" --min_naive "${MIN_NAIVE}" --penalty_K "${PENALTY_K}" --pi_min "${PI_MIN}"
  "${PYTHON}" "${PACKAGE_ROOT}/code/pantcr/infer_genotype_bayes.py" --sample_csv "${OUT_ROOT}/mutations/J/${sample_id}.J_sequences.csv" --out "${OUT_ROOT}/infer/PanTCR/J/infer_${sample_id}.J.csv" --pangenome_dir "${GRAPH_ROOT}/J" --min_naive "${MIN_NAIVE}" --penalty_K "${PENALTY_K}" --pi_min "${PI_MIN}"
done
shopt -u nullglob

if [[ -z "${SUMMARY_OUTPUT}" ]]; then
  SUMMARY_OUTPUT="${OUT_ROOT}/candidate_allele_summary.csv"
fi

"${PYTHON}" "${PACKAGE_ROOT}/code/experiments/05_public_bulk_candidate_nonimgt/summarize_public_bulk_candidates.py" \
  --input_dir "${OUT_ROOT}/mutations" \
  --index "${PACKAGE_ROOT}/data/ref/TRB_index.csv" \
  --ref "${PACKAGE_ROOT}/data/ref/pmTR_TRB_V_J_cleaned.csv" \
  --output "${SUMMARY_OUTPUT}"

echo "[DONE] public-bulk mutation and inference tree written to ${OUT_ROOT}"
echo "[DONE] public-bulk candidate summary written to ${SUMMARY_OUTPUT}"
