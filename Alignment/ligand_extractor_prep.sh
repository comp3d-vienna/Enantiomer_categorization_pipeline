#!/bin/bash
# Alignment step 2.4 — extract ligands from aligned structures with LigandExtractor.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$(cd "$SCRIPT_DIR/../Data" && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/../test" && pwd)"

if [[ "${PIPELINE_TEST:-}" == "1" ]]; then
    WORK_DIR="${TEST_DIR}"
else
    WORK_DIR="${DATA_DIR}"
fi

ALIGNMENT_DIR="${WORK_DIR}/Alignment"
ALIGNED_STRUCTURE_PDB_DIR="${ALIGNMENT_DIR}/paired_enantiomers_aligned_structures"
LIGAND_OUT_BASE_DIR="${ALIGNMENT_DIR}/Enantiomer_aligned_structure_ligandextract"
PAIRED_ENANTIOMERS_DIR="${ALIGNMENT_DIR}/paired_enantiomers"
LigandExtractor_newversion="${SCRIPT_DIR}/LigandExtractor_1.0.1/LigandExtractor"
LOG_FILE="${ALIGNMENT_DIR}/ligand_extractor_prep_errors.log"

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

parse_complex_pdb_chain() {
    local fname="$1"
    COMPLEX_PDB_ID=""
    COMPLEX_CHAIN_ID=""
    if [[ "$fname" =~ _([a-zA-Z0-9]+)_chain_([a-zA-Z0-9])(_copied)?\.pdb$ ]]; then
        COMPLEX_PDB_ID="${BASH_REMATCH[1]}"
        COMPLEX_CHAIN_ID="${BASH_REMATCH[2]}"
        return 0
    fi
    return 1
}

mkdir -p "$ALIGNED_STRUCTURE_PDB_DIR" "$LIGAND_OUT_BASE_DIR"

echo "=== Ligand Extractor Prep Log - $(date) ===" > "$LOG_FILE"
echo "=========================================="
echo "Extract ligands from aligned structures"
echo "=========================================="

if [ ! -x "$LigandExtractor_newversion" ]; then
    echo "Error: LigandExtractor not found or not executable at $LigandExtractor_newversion"
    echo "$(date): LigandExtractor missing at $LigandExtractor_newversion" >> "$LOG_FILE"
    exit 1
fi

if [[ ! -d "$PAIRED_ENANTIOMERS_DIR" ]] || [[ -z "$(find "$PAIRED_ENANTIOMERS_DIR" -name '*.csv' -print -quit 2>/dev/null)" ]]; then
    echo "Error: no paired enantiomer CSVs in $PAIRED_ENANTIOMERS_DIR" >&2
    echo "Run Alignment steps 2.1-2.2 first." >&2
    exit 1
fi

paired_ligand_ids=()
while IFS= read -r full_id; do
    [[ -z "$full_id" ]] && continue
    paired_ligand_ids+=("$full_id")
done < <(
    PAIRED_ENANTIOMERS_DIR="$PAIRED_ENANTIOMERS_DIR" "$PYTHON" - <<'PY'
import csv
import glob
import os

paired_dir = os.environ["PAIRED_ENANTIOMERS_DIR"]
ids = set()
for path in sorted(glob.glob(os.path.join(paired_dir, "*.csv"))):
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            ident = (row.get("identifier") or "").strip()
            m_layer = (row.get("m_layer") or "").strip()
            if ident and m_layer:
                ids.add(f"{ident}_{m_layer}")
for item in sorted(ids):
    print(item)
PY
)

if [[ ${#paired_ligand_ids[@]} -eq 0 ]]; then
    echo "Warning: no paired ligand IDs found in $PAIRED_ENANTIOMERS_DIR"
    echo "$(date): No paired ligand IDs in $PAIRED_ENANTIOMERS_DIR" >> "$LOG_FILE"
else
for full_id in "${paired_ligand_ids[@]}"; do
    IFS='_' read -r pdb_id het chain resseq rest <<< "$full_id"
    m_layer="${full_id##*_}"

    if [ -z "$pdb_id" ] || [ -z "$het" ] || [ -z "$chain" ] || [ -z "$resseq" ]; then
        msg="Warning: Could not parse ligand identifier from '$full_id'"
        echo "$(date): $msg" >> "$LOG_FILE"
        continue
    fi

    ligand_identifier="${het}_${chain}_${resseq}"

    candidate_complex_files=( )
    while IFS= read -r file; do
        candidate_complex_files+=( "$file" )
    done < <(find "$ALIGNED_STRUCTURE_PDB_DIR" -type f \( \
        -name "*_${pdb_id}_chain_${chain}.pdb" \
        -o -name "*_${pdb_id}_chain_${chain}_copied.pdb" \
        \) 2>/dev/null)

    matched_files=( )
    for candidate in "${candidate_complex_files[@]}"; do
        matched_name=$(basename "$candidate")
        if ! parse_complex_pdb_chain "$matched_name"; then
            msg="Warning: Could not parse pdb/chain from '$matched_name'"
            echo "$(date): $msg" >> "$LOG_FILE"
            continue
        fi
        if [ "$COMPLEX_PDB_ID" = "$pdb_id" ] && [ "$COMPLEX_CHAIN_ID" = "$chain" ]; then
            matched_files+=( "$candidate" )
        fi
    done

    if [ ${#matched_files[@]} -eq 0 ]; then
        msg="Warning: No aligned complex for pdb_id=$pdb_id chain=$chain (ligand $ligand_identifier)"
        echo "$(date): $msg" >> "$LOG_FILE"
        continue
    fi

    for matched_complex_file in "${matched_files[@]}"; do
        output_folder="${LIGAND_OUT_BASE_DIR}/${pdb_id}_${ligand_identifier}_${m_layer}"
        output_file="${output_folder}.sdf"
        if [ -e "$output_file" ]; then
            echo "Output $output_file exists; skipping."
            continue
        fi

        echo "Processing ligand '$ligand_identifier' from $matched_complex_file"
        if ! "$LigandExtractor_newversion" -c "$matched_complex_file" -l "$ligand_identifier" -o "$output_folder"; then
            msg="Error: LigandExtractor failed for $matched_complex_file -l $ligand_identifier"
            echo "$(date): $msg" >> "$LOG_FILE"
        fi
    done

done
fi

echo ""
echo "=========================================="
echo "Renaming .sdf files and removing folders"
echo "=========================================="

find "$LIGAND_OUT_BASE_DIR" -mindepth 2 -type f -name "*.sdf" | while read -r sdf_file; do
    folder_path=$(dirname "$sdf_file")
    parent_folder=$(basename "$folder_path")
    target_path="${LIGAND_OUT_BASE_DIR}/${parent_folder}.sdf"

    if [ -e "$target_path" ]; then
        rm -rf "$folder_path"
    else
        mv "$sdf_file" "$target_path"
        rm -rf "$folder_path"
    fi
done

find "$LIGAND_OUT_BASE_DIR" -mindepth 1 -type d -empty -delete 2>/dev/null || true

echo ""
echo "=========================================="
echo "Processing Complete"
echo "=========================================="
echo "Outputs: $LIGAND_OUT_BASE_DIR"
echo "Log: $LOG_FILE"
