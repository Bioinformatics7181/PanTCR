#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BULK_LABEL="BULK-1"
ACC_LIST=""
OUT_CLONES_DIR=""
WORK_DIR=""
THREADS="${THREADS:-8}"

usage() {
  cat <<EOF
Usage:
  bash run_sra_to_clones.sh --bulk-label BULK-1 --acc-list input/BULK-1/SRR_Acc_List.txt

Downloads public SRA accessions, converts them to FASTQ, runs MiXCR RNA-seq,
and writes *.clones_TRB.tsv files for run_from_clones.sh.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bulk-label) BULK_LABEL="$2"; shift 2 ;;
    --acc-list) ACC_LIST="$2"; shift 2 ;;
    --out-clones-dir) OUT_CLONES_DIR="$2"; shift 2 ;;
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "${ACC_LIST}" ]]; then
  echo "[ERROR] --acc-list is required." >&2
  usage
  exit 1
fi
if [[ -z "${OUT_CLONES_DIR}" ]]; then
  OUT_CLONES_DIR="${SCRIPT_DIR}/generated/${BULK_LABEL}/clones"
fi
if [[ -z "${WORK_DIR}" ]]; then
  WORK_DIR="${SCRIPT_DIR}/generated/${BULK_LABEL}/sra_work"
fi

mkdir -p "${OUT_CLONES_DIR}" "${WORK_DIR}"

while read -r accession || [[ -n "${accession}" ]]; do
  [[ -z "${accession}" ]] && continue
  echo "[RUN] ${BULK_LABEL} ${accession}"
  sample_work="${WORK_DIR}/${accession}"
  mkdir -p "${sample_work}"
  pushd "${sample_work}" >/dev/null

  prefetch "${accession}"
  fasterq-dump --split-files "${accession}"

  input_files=()
  if [[ -f "${accession}.fastq" ]]; then
    input_files=("${accession}.fastq")
  elif [[ -f "${accession}_1.fastq" && -f "${accession}_2.fastq" ]]; then
    input_files=("${accession}_1.fastq" "${accession}_2.fastq")
  else
    echo "[ERROR] No FASTQ files found for ${accession}" >&2
    exit 1
  fi

  mixcr analyze rna-seq --species hsa --threads "${THREADS}" "${input_files[@]}" "${accession}"

  if [[ -f "${accession}.clones_TRB.tsv" ]]; then
    cp "${accession}.clones_TRB.tsv" "${OUT_CLONES_DIR}/"
  else
    echo "[ERROR] Missing MiXCR TRB clones output for ${accession}" >&2
    exit 1
  fi

  popd >/dev/null
done < "${ACC_LIST}"

echo "[DONE] MiXCR clone tables written to ${OUT_CLONES_DIR}"
