#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"

PYTHON="${PYTHON:-python}"
THREADS="${THREADS:-8}"
SC_ID=""
DATASET_ID=""
CLONES_TSV=""
ALLELES_TSV=""
ALLELES_JSON=""
GENOTYPE_CSV=""
GEX_BAM=""
CHR_RANGE="${CHR_RANGE:-chr7:142000000-143000000}"
GRAPH_ROOT="${GRAPH_ROOT:-${PACKAGE_ROOT}/data/pretrained_graphs/trb_airr1}"
TRB_INDEX="${TRB_INDEX:-${PACKAGE_ROOT}/data/ref/TRB_index.csv}"
OUT_BASE="${OUT_BASE:-${SCRIPT_DIR}/generated/per_dataset_results}"
MIN_NAIVE="${MIN_NAIVE:-2}"
PENALTY_K="${PENALTY_K:-3}"
PI_MIN="${PI_MIN:-0.1}"

usage() {
  cat <<EOF
Usage:
  bash run_one_dataset.sh --sc-id SC-1 --dataset-id dataset2 --genotype input/.../genotype.csv [inputs]

Inputs:
  --clones        Existing MiXCR clones_TRB.tsv. Preferred for lightweight reruns.
  --alleles-tsv   MiXCR findAlleles mutation TSV for FindAlleles baseline.
  --alleles-json  MiXCR findAlleles custom allele JSON for FindAlleles baseline.
  --gex-bam       Optional 10x Gene Expression BAM. If --clones is absent, samtools+MiXCR are used to generate it.

PanTCR parameters:
  --min-naive INT  PanTCR min_naive. Default: ${MIN_NAIVE}
  --penalty-k INT  PanTCR penalty_K. Default: ${PENALTY_K}

Outputs:
  ${OUT_BASE}/{SC_ID}_{DATASET_ID}/
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sc-id) SC_ID="$2"; shift 2 ;;
    --dataset-id) DATASET_ID="$2"; shift 2 ;;
    --clones) CLONES_TSV="$2"; shift 2 ;;
    --alleles-tsv) ALLELES_TSV="$2"; shift 2 ;;
    --alleles-json) ALLELES_JSON="$2"; shift 2 ;;
    --genotype) GENOTYPE_CSV="$2"; shift 2 ;;
    --gex-bam) GEX_BAM="$2"; shift 2 ;;
    --chr-range) CHR_RANGE="$2"; shift 2 ;;
    --graph-root) GRAPH_ROOT="$2"; shift 2 ;;
    --index) TRB_INDEX="$2"; shift 2 ;;
    --out-base) OUT_BASE="$2"; shift 2 ;;
    --min-naive) MIN_NAIVE="$2"; shift 2 ;;
    --penalty-k) PENALTY_K="$2"; shift 2 ;;
    --pi-min) PI_MIN="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "${SC_ID}" || -z "${DATASET_ID}" || -z "${GENOTYPE_CSV}" ]]; then
  echo "[ERROR] --sc-id, --dataset-id, and --genotype are required." >&2
  usage
  exit 1
fi

DATASET_DIR="${OUT_BASE}/${SC_ID}_${DATASET_ID}"
mkdir -p "${DATASET_DIR}"
cp "${GENOTYPE_CSV}" "${DATASET_DIR}/genotype.csv"

if [[ -z "${CLONES_TSV}" ]]; then
  if [[ -z "${GEX_BAM}" ]]; then
    echo "[ERROR] Provide --clones or --gex-bam." >&2
    exit 1
  fi
  mkdir -p "${DATASET_DIR}/gex"
  samtools view -b "${GEX_BAM}" "${CHR_RANGE}" > "${DATASET_DIR}/gex/trb.region.bam"
  samtools view "${DATASET_DIR}/gex/trb.region.bam" | cut -f 1 | sort | uniq > "${DATASET_DIR}/gex/read_names.txt"
  samtools view -b -N "${DATASET_DIR}/gex/read_names.txt" "${GEX_BAM}" > "${DATASET_DIR}/gex/trb.complete.bam"
  samtools sort -n -@ "${THREADS}" "${DATASET_DIR}/gex/trb.complete.bam" -o "${DATASET_DIR}/gex/trb.complete.name_sorted.bam"
  samtools fastq -@ "${THREADS}" -n "${DATASET_DIR}/gex/trb.complete.name_sorted.bam" | gzip > "${DATASET_DIR}/gex/trb_all.fastq.gz"
  mixcr analyze -s hsa --threads "${THREADS}" rna-seq --force-overwrite "${DATASET_DIR}/gex/trb_all.fastq.gz" "${DATASET_DIR}/gex/trb_all"
  mixcr assembleContigs --force-overwrite --assemble-contigs-by VDJRegion "${DATASET_DIR}/gex/trb_all.clna" "${DATASET_DIR}/gex/trb_all.contigs.VDJRegion.clns"
  mixcr findAlleles --no-clns-output --export-alleles-mutations "${DATASET_DIR}/gex/trb_all.alleles.tsv" --report "${DATASET_DIR}/gex/trb_all.findAlleles.report.txt" --export-library "${DATASET_DIR}/gex/trb_all.customAlleles.json" --force-overwrite "${DATASET_DIR}/gex/trb_all.contigs.VDJRegion.clns"
  CLONES_TSV="${DATASET_DIR}/gex/trb_all.clones_TRB.tsv"
  ALLELES_TSV="${DATASET_DIR}/gex/trb_all.alleles.tsv"
  ALLELES_JSON="${DATASET_DIR}/gex/trb_all.customAlleles.json"
fi

"${PYTHON}" "${PACKAGE_ROOT}/code/pantcr/collect_mutation.py" --input "${CLONES_TSV}" --gene V --ref "${PACKAGE_ROOT}/data/ref/IMGT_TRBV_pro.tsv" --prefix "${DATASET_DIR}/${DATASET_ID}.V"

mkdir -p "${DATASET_DIR}/MiXCR-default" "${DATASET_DIR}/FindAlleles" "${DATASET_DIR}/PanTCR.semi"
"${PYTHON}" "${PACKAGE_ROOT}/code/experiments/00_benchmark_utils/clones2alleles.py" --clones "${CLONES_TSV}" --ref "${PACKAGE_ROOT}/data/ref/IMGT_TRB_default.csv" --output "${DATASET_DIR}/MiXCR-default/infer_MiXCR.csv"

if [[ -n "${ALLELES_TSV}" && -n "${ALLELES_JSON}" ]]; then
  "${PYTHON}" "${PACKAGE_ROOT}/code/experiments/00_benchmark_utils/infer_for_findalleles.py" --tsv "${ALLELES_TSV}" --json "${ALLELES_JSON}" --ref "${PACKAGE_ROOT}/data/ref/IMGT_TRB_default.csv" --output "${DATASET_DIR}/FindAlleles/infer_findalleles.csv"
else
  echo "[WARN] Missing FindAlleles TSV/JSON; FindAlleles baseline will not be generated." >&2
fi

"${PYTHON}" "${PACKAGE_ROOT}/code/pantcr/infer_genotype_bayes.py" --sample_csv "${DATASET_DIR}/${DATASET_ID}.V_sequences.csv" --out "${DATASET_DIR}/PanTCR.semi/infer_PanTCR.semi.V.csv" --pangenome_dir "${GRAPH_ROOT}/V" --min_naive "${MIN_NAIVE}" --penalty_K "${PENALTY_K}" --pi_min "${PI_MIN}"

"${PYTHON}" "${PACKAGE_ROOT}/code/experiments/00_benchmark_utils/evaluate_allele_calls.py" --gt "${DATASET_DIR}/genotype.csv" --infer "${DATASET_DIR}/MiXCR-default/infer_MiXCR.csv" --index "${TRB_INDEX}" --gene_type V --out_prefix "${DATASET_DIR}/MiXCR-default/eval_MiXCR.V"
if [[ -f "${DATASET_DIR}/FindAlleles/infer_findalleles.csv" ]]; then
  "${PYTHON}" "${PACKAGE_ROOT}/code/experiments/00_benchmark_utils/evaluate_allele_calls.py" --gt "${DATASET_DIR}/genotype.csv" --infer "${DATASET_DIR}/FindAlleles/infer_findalleles.csv" --index "${TRB_INDEX}" --gene_type V --out_prefix "${DATASET_DIR}/FindAlleles/eval_findalleles.V"
fi
"${PYTHON}" "${PACKAGE_ROOT}/code/experiments/00_benchmark_utils/evaluate_allele_calls.py" --gt "${DATASET_DIR}/genotype.csv" --infer "${DATASET_DIR}/PanTCR.semi/infer_PanTCR.semi.V.csv" --index "${TRB_INDEX}" --gene_type V --out_prefix "${DATASET_DIR}/PanTCR.semi/eval_PanTCR.semi.V"

echo "[DONE] pseudo-bulk dataset result written to ${DATASET_DIR}"
