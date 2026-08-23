#!/bin/bash
# Run Categorization (step 4) on production data under Data/.
#
# Default: 4a whole_process + 4b export training table + 4c custom RMSD thresholds.

if [ -z "${BASH_VERSION:-}" ]; then
    exec /bin/bash "$0" "$@"
fi

set -euo pipefail

unset PIPELINE_TEST

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/Data"
CATEGORIZATION_DIR="${SCRIPT_DIR}/Categorization"
MANUAL_CURATION_DIR="${DATA_DIR}/Manual_curation"

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

FEATURES_ONLY=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Run binding-mode categorization for paired enantiomer groups (step 4).

Steps (default: all):
  4a  whole_process.py              Manual training/invalid + closest-pair RMSD
  4b  export_threshold_training_table.py   RMSD + SILIRID similarity
  4c  train_category_thresholds.py  custom RMSD bins: purity and recall

Prerequisite: complete Pharmacophore generation (steps 3.1-3.2), especially:
  Data/Pharmacophore/paired_enantiomers_pocket_based/
  Data/Pharmacophore/Enantiomer_aligned_structure_ligandextract_canonical/

For step 4b-4c, manual labels under Data/Manual_curation/ are expected.

Delete Data/Categorization/ yourself if you need a full redo.

Options:
  --features-only   Run 4a only (skip export + custom RMSD scoring)
  -h, --help        Show this help message

Environment:
  PYTHON_CMD  Python executable (default: python3, then python)

Examples:
  bash run_categorization.sh
  bash run_categorization.sh --features-only
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --features-only)
            FEATURES_ONLY=1
            shift
            ;;
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
        echo "Activate Python 3.9 pipeline env, e.g.: conda activate categorize_pipeline" >&2
        echo "Or install: $PYTHON -m pip install -r ${SCRIPT_DIR}/requirements.txt" >&2
        exit 1
    fi
}

require_pharmacophore_inputs() {
    local missing=0
    if [[ ! -d "${DATA_DIR}/Pharmacophore/paired_enantiomers_pocket_based" ]]; then
        echo "Missing: Data/Pharmacophore/paired_enantiomers_pocket_based/"
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

if [[ ! -d "${MANUAL_CURATION_DIR}" ]]; then
    echo "Error: ${MANUAL_CURATION_DIR} not found — manual curation is required for categorization." >&2
    exit 1
fi

echo ""
echo "=== Step 4a: Manual training/invalid + closest-pair RMSD ==="
"$PYTHON" "${CATEGORIZATION_DIR}/whole_process.py"

if [[ "$FEATURES_ONLY" -eq 1 ]]; then
    echo ""
    echo "Categorization finished (features only; skipped 4b-4c)"
    echo "Outputs: ${DATA_DIR}/Categorization/"
    exit 0
fi

echo ""
echo "=== Step 4b: Export threshold training table (RMSD + SILIRID) ==="
"$PYTHON" "${CATEGORIZATION_DIR}/export_threshold_training_table.py"

echo ""
echo "=== Step 4c: Custom RMSD thresholds (purity and recall) ==="
"$PYTHON" "${CATEGORIZATION_DIR}/train_category_thresholds.py"

echo ""
echo "Categorization finished"
echo "Outputs: ${DATA_DIR}/Categorization/"
echo "  complex_categories.csv"
echo "  threshold_training_table.csv"
echo "  silirid_fingerprint_slots.csv"
echo "  custom_rmsd_thresholds.csv"
echo "  custom_rmsd_threshold_report.txt"
