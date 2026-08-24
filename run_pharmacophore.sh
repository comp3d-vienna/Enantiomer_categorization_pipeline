#!/bin/bash
# Run Pharmacophore generation steps 3.1-3.2.

if [ -z "${BASH_VERSION:-}" ]; then
    exec /bin/bash "$0" "$@"
fi

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHARMACOPHORE_DIR="${SCRIPT_DIR}/Pharmacophore_generation_after_alignment"
PROJECT_ROOT="${SCRIPT_DIR}"

unset PIPELINE_TEST
# shellcheck source=Data_preprocess/data_root.sh
source "${PROJECT_ROOT}/Data_preprocess/data_root.sh"
DATA_DIR="$(select_pipeline_data_dir)"
export PIPELINE_DATA_DIR="$DATA_DIR"
ALIGNMENT_DIR="${DATA_DIR}/Alignment"

if [[ -n "${PYTHON_CMD:-}" ]]; then
    PYTHON="$PYTHON_CMD"
elif command -v python3 >/dev/null 2>&1; then
    PYTHON=python3
elif command -v python >/dev/null 2>&1; then
    PYTHON=python
else
    echo "Error: python3 or python not found in PATH." >&2
    echo "Activate your Python 3.13 env or set PYTHON_CMD, e.g.:" >&2
    echo "  conda activate categorize_pipeline" >&2
    echo "  export PYTHON_CMD=\$(which python)" >&2
    exit 1
fi
export PYTHON
export PYTHON_CMD="${PYTHON_CMD:-$PYTHON}"

SKIP_GENERATION=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Run pharmacophore generation for paired enantiomers (steps 3.1-3.2).

Prerequisite: complete Alignment (steps 2.1-2.4), especially:
  Data/Alignment/paired_enantiomers/
  Data/Alignment/paired_enantiomers_aligned_structures/
  Data/Alignment/Enantiomer_aligned_structure_ligandextract/
  Data/db_ids_with_m_layer.txt

Existing canonical SDFs, PML files, and interaction TSVs are skipped.
Delete Data/Pharmacophore/ yourself if you need a full redo.

Options:
  --skip-generation   Skip step 3.1 (gen_ia_ph4s_fg / interaction TSV export)
  -h, --help          Show this help message

Environment:
  PYTHON_CMD          Python executable (default: python3, then python)

Examples:
  conda activate categorize_pipeline
  bash run_pharmacophore.sh
  bash run_pharmacophore.sh --skip-generation
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-generation)
            SKIP_GENERATION=1
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

step() {
    echo ""
    echo "============================================================"
    echo "$1"
    echo "============================================================"
}

require_alignment_inputs() {
    local missing=0
    if [[ ! -d "${ALIGNMENT_DIR}/paired_enantiomers" ]]; then
        echo "Missing: ${DATA_DIR}/Alignment/paired_enantiomers/"
        missing=1
    fi
    if [[ ! -d "${ALIGNMENT_DIR}/paired_enantiomers_aligned_structures" ]]; then
        echo "Missing: ${DATA_DIR}/Alignment/paired_enantiomers_aligned_structures/"
        missing=1
    fi
    if [[ ! -d "${ALIGNMENT_DIR}/Enantiomer_aligned_structure_ligandextract" ]]; then
        echo "Missing: ${DATA_DIR}/Alignment/Enantiomer_aligned_structure_ligandextract/"
        missing=1
    fi
    if [[ ! -f "${DATA_DIR}/db_ids_with_m_layer.txt" ]]; then
        echo "Missing: ${DATA_DIR}/db_ids_with_m_layer.txt"
        missing=1
    fi
    if [[ "$missing" -eq 1 ]]; then
        echo "Error: Alignment outputs not found." >&2
        echo "Run: bash run_alignment.sh" >&2
        exit 1
    fi
}

echo "Pharmacophore pipeline"
echo "Data directory: ${DATA_DIR}"

require_alignment_inputs
echo "Python: $PYTHON"

if [[ "$SKIP_GENERATION" -eq 0 ]]; then
    step "Step 3.1: Interaction pharmacophore generation"
    bash "${PHARMACOPHORE_DIR}/pharmacophore_generation_new.sh"
else
    step "Step 3.1: Skipped - pharmacophore generation disabled"
fi

step "Step 3.2: Process interaction data (paired enantiomers)"
"$PYTHON" "${PHARMACOPHORE_DIR}/process_interaction_data.py"

step "Pharmacophore pipeline finished"
echo "Outputs are under: ${DATA_DIR}/Pharmacophore/"
echo "  Interaction_data/                         step 3.1"
echo "  Pharmacophore_generation_results/         step 3.1"
echo "  paired_enantiomers_pocket_based/          step 3.2 (one CSV per group)"
echo ""
echo "Next: bash run_categorization.sh"
