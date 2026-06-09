#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PYTHON="${PYTHON:-python}"

"${PYTHON}" "${PACKAGE_ROOT}/code/experiments/07_semi_synthetic_airr_benchmark/semi_simu_fine_evidence_analysis.py" "$@"
