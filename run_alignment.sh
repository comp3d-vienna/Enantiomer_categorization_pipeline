#!/bin/bash
# Run Alignment steps 2.1-2.4 on production data under Data/.

if [ -z "${BASH_VERSION:-}" ]; then
    exec /bin/bash "$0" "$@"
fi

set -euo pipefail

# Ensure production paths (not test/)
unset PIPELINE_TEST

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/Data"
ALIGNMENT_DIR="${SCRIPT_DIR}/Alignment"

if [[ -n "${PYTHON_CMD:-}" ]]; then
    PYTHON="$PYTHON_CMD"
elif command -v python3 >/dev/null 2>&1; then
    PYTHON=python3
elif command -v python >/dev/null 2>&1; then
    PYTHON=python
else
    echo "Error: python3 or python not found in PATH." >&2
    echo "Activate your Python 3.9 env or set PYTHON_CMD, e.g.:" >&2
    echo "  conda activate categorize_pipeline" >&2
    echo "  export PYTHON_CMD=\$(which python)" >&2
    exit 1
fi

SKIP_MAE_CONVERT=0
SKIP_STRUCTALIGN=0
SKIP_LIGAND_EXTRACT=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Run Alignment steps 2.1-2.4 on production outputs under Data/.

Prerequisite: complete Data_preprocess (steps 1.1-1.4), especially:
  Data/chiral_classification_results/filtered_enantiomer_with_uniprot_new.txt
  Data/db_ids_with_m_layer.txt
  Data/Enantiomer_protein_preperation/Enantiomer_prepwizard_results/*_prepped.mae

Existing prep PDBs, parsed chains, and ligand SDFs are skipped. Delete
Data/Alignment/ yourself if you need a full redo.

Options:
  --skip-mae-convert    Skip Schrödinger MAE to PDB in step 2.1-2.2
  --skip-structalign    Skip Schrödinger structalign (step 2.3)
  --skip-ligand-extract Skip LigandExtractor (step 2.4)
  -h, --help            Show this help message

Environment:
  PYTHON_CMD            Python executable (default: python3, then python)
  SCHRODINGER           Required for MAE->PDB and structalign unless skipped
                        (e.g. export SCHRODINGER=/path/to/your/schrodinger)

Examples:
  bash run_alignment.sh
  bash run_alignment.sh --skip-structalign --skip-ligand-extract
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-mae-convert)
            SKIP_MAE_CONVERT=1
            shift
            ;;
        --skip-structalign)
            SKIP_STRUCTALIGN=1
            shift
            ;;
        --skip-ligand-extract)
            SKIP_LIGAND_EXTRACT=1
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

require_preprocess_inputs() {
    local missing=0
    if [[ ! -f "${DATA_DIR}/chiral_classification_results/filtered_enantiomer_with_uniprot_new.txt" ]]; then
        echo "Missing: Data/chiral_classification_results/filtered_enantiomer_with_uniprot_new.txt"
        missing=1
    fi
    if [[ ! -f "${DATA_DIR}/db_ids_with_m_layer.txt" ]]; then
        echo "Missing: Data/db_ids_with_m_layer.txt"
        missing=1
    fi
    local prep_mae_dir="${DATA_DIR}/Enantiomer_protein_preperation/Enantiomer_prepwizard_results"
    if [[ "$SKIP_MAE_CONVERT" -eq 0 && -n "${SCHRODINGER:-}" ]]; then
        if [[ ! -d "$prep_mae_dir" ]] || [[ -z "$(find "$prep_mae_dir" -name '*_prepped.mae' -print -quit 2>/dev/null)" ]]; then
            echo "Warning: no *_prepped.mae files under ${prep_mae_dir}"
            echo "         Run Data_preprocess step 1.4 (maestro_protein_preparation.sh) first."
        fi
    fi
    if [[ "$missing" -eq 1 ]]; then
        echo "Error: Data_preprocess outputs not found under Data/." >&2
        echo "Complete Data_preprocess before running alignment." >&2
        exit 1
    fi
}

echo "Production alignment pipeline"
echo "Data directory: ${DATA_DIR}"
echo "SCHRODINGER:    ${SCHRODINGER:-<not set>}"

require_preprocess_inputs

step "Steps 2.1-2.2: Group enantiomers, MAE->PDB, parse chains"
GROUP_ARGS=()
if [[ "$SKIP_MAE_CONVERT" -eq 1 ]]; then
    GROUP_ARGS+=(--skip-mae-convert)
elif [[ -z "${SCHRODINGER:-}" ]]; then
    echo "Error: SCHRODINGER is not set (required for MAE to PDB)." >&2
    echo "Set SCHRODINGER or use --skip-mae-convert if PDBs already exist." >&2
    exit 1
fi
"$PYTHON" "${ALIGNMENT_DIR}/group_and_parse_chains.py" "${GROUP_ARGS[@]}"

if [[ "$SKIP_STRUCTALIGN" -eq 1 ]]; then
    step "Step 2.3: Skipped - structalign disabled"
elif [[ -z "${SCHRODINGER:-}" ]]; then
    echo "Error: SCHRODINGER is not set (required for structalign)." >&2
    exit 1
else
    step "Step 2.3: Schrödinger structalign"
    bash "${ALIGNMENT_DIR}/schrodinger.sh"
fi

if [[ "$SKIP_LIGAND_EXTRACT" -eq 1 ]]; then
    step "Step 2.4: Skipped - ligand extraction disabled"
else
    step "Step 2.4: Ligand extraction"
    bash "${ALIGNMENT_DIR}/ligand_extractor_prep.sh"
fi

step "Alignment pipeline finished"
echo "Outputs are under: ${DATA_DIR}/Alignment/"
echo "  paired_enantiomers/                       step 2.1"
echo "  single_enantiomers/                       step 2.1"
echo "  group_manifest.tsv                        step 2.1"
echo "  paired_enantiomers_parsed_chain/          step 2.2"
echo "  paired_enantiomers_aligned_structures/    step 2.3"
echo "  Enantiomer_aligned_structure_ligandextract/ step 2.4"
echo "Prep PDB (MAE conversion): ${DATA_DIR}/Alignment/Enantiomer_prep_pdbformat/"
echo ""
echo "Next: bash run_pharmacophore.sh"
