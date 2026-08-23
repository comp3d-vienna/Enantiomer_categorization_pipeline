#!/bin/bash
# Run Pharmacophore generation steps 3.1-3.2 on production data under Data/.

if [ -z "${BASH_VERSION:-}" ]; then
    exec /bin/bash "$0" "$@"
fi

set -euo pipefail

unset PIPELINE_TEST

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/Data"
PHARMACOPHORE_DIR="${SCRIPT_DIR}/Pharmacophore_generation_after_alignment"
ALIGNMENT_DIR="${DATA_DIR}/Alignment"
RESOLVE_PYTHON="${PHARMACOPHORE_DIR}/resolve_python.sh"

SKIP_GENERATION=0

# shellcheck source=Pharmacophore_generation_after_alignment/resolve_python.sh
source "$RESOLVE_PYTHON"

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
  PYTHON_CMD          Python 3.9 for all pipeline steps (default: active env)
  PIPELINE_CONDA_ENV  Conda env name fallback (default: categorize_pipeline)

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
        echo "Missing: Data/Alignment/paired_enantiomers/"
        missing=1
    fi
    if [[ ! -d "${ALIGNMENT_DIR}/paired_enantiomers_aligned_structures" ]]; then
        echo "Missing: Data/Alignment/paired_enantiomers_aligned_structures/"
        missing=1
    fi
    if [[ ! -d "${ALIGNMENT_DIR}/Enantiomer_aligned_structure_ligandextract" ]]; then
        echo "Missing: Data/Alignment/Enantiomer_aligned_structure_ligandextract/"
        missing=1
    fi
    if [[ ! -f "${DATA_DIR}/db_ids_with_m_layer.txt" ]]; then
        echo "Missing: Data/db_ids_with_m_layer.txt"
        missing=1
    fi
    if [[ "$missing" -eq 1 ]]; then
        echo "Error: Alignment outputs not found." >&2
        echo "Run: bash run_alignment.sh" >&2
        exit 1
    fi
}

echo "Production pharmacophore pipeline"
echo "Data directory: ${DATA_DIR}"

require_alignment_inputs
resolve_python
echo "Python: $PIPELINE_PYTHON"
echo "PYTHONPATH (CDPL): ${PYTHONPATH:-<not set>}"

if [[ "$SKIP_GENERATION" -eq 0 ]]; then
    step "Step 3.1: Interaction pharmacophore generation"
    bash "${PHARMACOPHORE_DIR}/pharmacophore_generation_new.sh"
else
    step "Step 3.1: Skipped - pharmacophore generation disabled"
fi

step "Step 3.2: Process interaction data (paired enantiomers)"
"$PIPELINE_PYTHON" "${PHARMACOPHORE_DIR}/process_interaction_data.py"

step "Pharmacophore pipeline finished"
echo "Outputs are under: ${DATA_DIR}/Pharmacophore/"
echo "  Interaction_data/                         step 3.1"
echo "  Pharmacophore_generation_results/         step 3.1"
echo "  paired_enantiomers_pocket_based/          step 3.2 (one CSV per group)"
echo ""
echo "Next: bash run_categorization.sh"
