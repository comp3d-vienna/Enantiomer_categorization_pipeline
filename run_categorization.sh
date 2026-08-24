#!/bin/bash
# Run Categorization (step 4) on production data under Data/.
#
# 4a whole_process.py (closest-pair RMSD) + 4b export_feature_table.py (SILIRID).
# Both write Data/Categorization/feature_table.csv.

if [ -z "${BASH_VERSION:-}" ]; then
    exec /bin/bash "$0" "$@"
fi

set -euo pipefail

unset PIPELINE_TEST

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/Data"
CATEGORIZATION_DIR="${SCRIPT_DIR}/Categorization"

if [[ -n "${PYTHON_CMD:-}" ]]; then
    PYTHON="$PYTHON_CMD"
elif command -v python3 >/dev/null 2>&1; then
    PYTHON=python3
elif command -v python >/dev/null 2>&1; then
    PYTHON=python
else
    echo "Error: python3 or python not found in PATH." >&2
    exit 1
fi

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Compute closest-pair RMSD and SILIRID features for paired enantiomer groups (step 4).

Steps:
  4a  whole_process.py         Closest-pair RMSD → feature_table.csv
  4b  export_feature_table.py  Add SILIRID similarity and fingerprints to that table

Prerequisite: complete Pharmacophore generation (steps 3.1-3.2), especially:
  Data/Pharmacophore/paired_enantiomers_pocket_based/
  Data/Pharmacophore/Enantiomer_aligned_structure_ligandextract_canonical/

Delete Data/Categorization/ yourself if you need a full redo.

Options:
  -h, --help  Show this help message

Environment:
  PYTHON_CMD  Python executable (default: python3, then python)

Examples:
  bash run_categorization.sh
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage
            exit 1
            ;;
    esac
done

check_python_deps() {
    if ! "$PYTHON" -c "import pandas, pytz, rdkit" 2>/dev/null; then
        echo "Error: missing Python dependencies (pandas, pytz, and/or rdkit)." >&2
        echo "Python used: $PYTHON" >&2
        echo "Activate Python 3.13 pipeline env, e.g.: conda activate categorize_pipeline" >&2
        echo "Or install: $PYTHON -m pip install -r ${SCRIPT_DIR}/requirements.txt" >&2
        exit 1
    fi
}

require_pharmacophore_inputs() {
    local missing=0
    if [[ ! -d "${DATA_DIR}/Pharmacophore/paired_enantiomers_pocket_based" ]] \
        || [[ -z "$(find "${DATA_DIR}/Pharmacophore/paired_enantiomers_pocket_based" -name '*.csv' -print -quit 2>/dev/null)" ]]; then
        echo "Missing: Data/Pharmacophore/paired_enantiomers_pocket_based/*.csv (step 3.2)"
        missing=1
    fi
    if [[ ! -d "${DATA_DIR}/Pharmacophore/Enantiomer_aligned_structure_ligandextract_canonical" ]]; then
        echo "Missing: Data/Pharmacophore/Enantiomer_aligned_structure_ligandextract_canonical/"
        missing=1
    fi
    if [[ "$missing" -eq 1 ]]; then
        echo "Error: Pharmacophore outputs not found." >&2
        echo "Run: bash run_pharmacophore.sh" >&2
        exit 1
    fi
}

echo "Production categorization pipeline"
echo "Data directory: ${DATA_DIR}"

check_python_deps
require_pharmacophore_inputs

echo ""
echo "=== Step 4a: Closest-pair RMSD ==="
"$PYTHON" "${CATEGORIZATION_DIR}/whole_process.py"

echo ""
echo "=== Step 4b: Add SILIRID similarity and fingerprints ==="
"$PYTHON" "${CATEGORIZATION_DIR}/export_feature_table.py"

echo ""
echo "Categorization finished"
echo "Outputs: ${DATA_DIR}/Categorization/"
echo "  feature_table.csv"
echo "  silirid_fingerprint_slots.csv"
