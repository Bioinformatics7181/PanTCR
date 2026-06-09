#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PYTHON="${PYTHON:-python}"

"${PYTHON}" "${PACKAGE_ROOT}/code/experiments/05_public_bulk_candidate_nonimgt/annotate_bulk_evidence_imputation.py" "$@"
