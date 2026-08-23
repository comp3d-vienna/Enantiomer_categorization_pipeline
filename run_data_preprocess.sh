#!/bin/bash
# Run all Data_preprocess steps on production data under Data/.

if [ -z "${BASH_VERSION:-}" ]; then
    exec /bin/bash "$0" "$@"
fi

set -euo pipefail

# Ensure production paths (not test/)
unset PIPELINE_TEST

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/Data"
PREPROCESS_DIR="${SCRIPT_DIR}/Data_preprocess"
export MMCIF_SOURCE="${MMCIF_SOURCE:-${DATA_DIR}/mmCIF}"
export MMCIF_RENAME="${MMCIF_RENAME:-${DATA_DIR}/mmCIF_rename}"

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

SKIP_UNIPROT=0
SKIP_PREPWIZARD=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Run the full Data_preprocess pipeline (steps 1.1-1.4) on production data.

Set SCHRODINGER, MMCIF_SOURCE, and MMCIF_RENAME to your local paths
before running (see Data_preprocess/README.md). If MMCIF_* are unset,
mmCIF input is Data/mmCIF and decompressed CIFs go to Data/mmCIF_rename.
Pipeline outputs are written under Data/.

Options:
  --skip-uniprot      Skip RCSB UniProt API queries in step 1.3
  --skip-prepwizard   Skip Schrödinger PrepWizard (step 1.4)
  -h, --help          Show this help message

Environment:
  PYTHON_CMD          Python executable (default: python3, then python)
  MMCIF_SOURCE        Your PDB mmCIF gzip archive (default: Data/mmCIF)
  MMCIF_RENAME        Your decompressed mmCIF directory (default: Data/mmCIF_rename)
  SCHRODINGER         Your Schrödinger install root (required for step 1.4
                    unless --skip-prepwizard is used)

Examples:
  export SCHRODINGER=/path/to/your/schrodinger
  export MMCIF_SOURCE=/path/to/your/mmCIF
  export MMCIF_RENAME=/path/to/your/mmCIF_rename
  bash run_data_preprocess.sh
  bash run_data_preprocess.sh --skip-uniprot --skip-prepwizard
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-uniprot)
            SKIP_UNIPROT=1
            shift
            ;;
        --skip-prepwizard)
            SKIP_PREPWIZARD=1
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

require_mmcif_source() {
    if [[ ! -d "$MMCIF_SOURCE" ]]; then
        echo "Error: mmCIF source not found: $MMCIF_SOURCE" >&2
        echo "Place downloaded mmCIF files in Data/mmCIF, or set MMCIF_SOURCE." >&2
        exit 1
    fi
    if [[ -z "$(find "$MMCIF_SOURCE" -name '*.cif.gz' -print -quit 2>/dev/null)" ]]; then
        echo "Error: no .cif.gz files found under $MMCIF_SOURCE" >&2
        exit 1
    fi
}

echo "Production Data_preprocess pipeline"
echo "Data directory:  ${DATA_DIR}"
echo "mmCIF source:    ${MMCIF_SOURCE}"
echo "mmCIF rename:    ${MMCIF_RENAME}"
echo "SCHRODINGER:     ${SCHRODINGER:-<not set>}"

require_mmcif_source

cd "$PREPROCESS_DIR"

step "Step 1.1a: Extract and split ligands"
bash extract_split_ligands.sh

step "Step 1.1b: Convert SDF to InChI"
"$PYTHON" convert_sdf_to_inchi.py

step "Step 1.1c: Classify chiral compounds"
"$PYTHON" classify_chiral_compounds.py

step "Step 1.2-1.3b: Filter drug-like enantiomers and annotate"
if [[ "$SKIP_UNIPROT" -eq 1 ]]; then
    "$PYTHON" filter_and_annotate_enantiomers.py --skip-uniprot
else
    "$PYTHON" filter_and_annotate_enantiomers.py
fi

if [[ "$SKIP_PREPWIZARD" -eq 1 ]]; then
    step "Step 1.4: Skipped - prepwizard disabled"
elif [[ -z "${SCHRODINGER:-}" ]]; then
    echo "Error: SCHRODINGER is not set (required for step 1.4)." >&2
    echo "Set SCHRODINGER or use --skip-prepwizard." >&2
    exit 1
else
    step "Step 1.4: Prepare protein structures"
    bash maestro_protein_preparation.sh
fi

step "Data_preprocess pipeline finished"
echo "Outputs are under: ${DATA_DIR}/"
echo "  ligand_sdf_structure_from_mmCIF/          step 1.1a"
echo "  ligand_sdf_structure_from_mmCIF_separated/ step 1.1a"
echo "  ligand_sdf_structure_from_mmCIF_separated_InChI/ step 1.1b"
echo "  chiral_classification_results/            steps 1.1c-1.3"
echo "  db_ids_with_m_layer.txt                   step 1.3b"
if [[ "$SKIP_PREPWIZARD" -eq 0 && -n "${SCHRODINGER:-}" ]]; then
    echo "  Enantiomer_pdbstructure_maeformat/        step 1.4"
    echo "  Enantiomer_protein_preperation/           step 1.4"
fi
echo ""
echo "Next: bash run_alignment.sh"
