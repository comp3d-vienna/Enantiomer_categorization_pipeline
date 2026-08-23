#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$(cd "$SCRIPT_DIR/../Data" && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/../test" && pwd)"

if [[ "${PIPELINE_TEST:-}" == "1" ]]; then
    WORK_DIR="${TEST_DIR}"
    CIF_DIR="${MMCIF_RENAME:-${WORK_DIR}/mmcif_rename}"
    Enantiomer_list_identifier="${WORK_DIR}/db_ids_with_m_layer.txt"
    PREPWIZARD_TIMEOUT_LOG="${WORK_DIR}/prepwizard_timeout_pdbids.txt"
    PREPWIZARD_FAILED_LOG="${WORK_DIR}/prepwizard_failed_pdbids.txt"
else
    WORK_DIR="${DATA_DIR}"
    CIF_DIR="${MMCIF_RENAME:-${DATA_DIR}/mmCIF_rename}"
    Enantiomer_list_identifier="${WORK_DIR}/db_ids_with_m_layer.txt"
    PREPWIZARD_TIMEOUT_LOG="${WORK_DIR}/prepwizard_timeout_pdbids.txt"
    PREPWIZARD_FAILED_LOG="${WORK_DIR}/prepwizard_failed_pdbids.txt"
fi

MAE_OUTPUT_DIR="${WORK_DIR}/Enantiomer_pdbstructure_maeformat"
PREPWIZARD_OUTPUT_DIR="${WORK_DIR}/Enantiomer_protein_preperation/Enantiomer_prepwizard_results"

# All required directories must exist
mkdir -p "$MAE_OUTPUT_DIR" "$PREPWIZARD_OUTPUT_DIR"
cd "$WORK_DIR"
echo "Schrödinger launch directory: $(pwd)"
echo "CIF directory: $CIF_DIR"

# Get a sorted, unique list of PDB IDs from the identifier file
mapfile -t unique_pdb_id_list < <(awk -F_ '{print $1}' "$Enantiomer_list_identifier" | sort -u)
echo "The number of unique PDB IDs is: ${#unique_pdb_id_list[@]}"

# Only convert .cif to .mae if .mae doesn't already exist, for each PDB ID
for pdb_id in "${unique_pdb_id_list[@]}"; do
    cif_file="${CIF_DIR}/${pdb_id}.cif"
    mae_file="${MAE_OUTPUT_DIR}/${pdb_id}.mae"
    if [[ ! -f "$mae_file" ]]; then
        if [[ -f "$cif_file" ]]; then
            echo "Converting $cif_file to $mae_file ..."
            "$SCHRODINGER/utilities/structconvert" "$cif_file" "$mae_file"
            if [[ $? -ne 0 ]]; then
                echo "Conversion failed for $cif_file"
            else
                echo "Successfully converted $cif_file to $mae_file"
            fi
        else
            echo "CIF file not found for PDB ID $pdb_id at $cif_file"
        fi
    else
        echo "Mae file already exists for PDB ID $pdb_id at $mae_file. Skipping conversion."
    fi
done

# PrepWizard requires output paths and input files as absolute paths.
for pdb_id in "${unique_pdb_id_list[@]}"; do
    mae_file="${MAE_OUTPUT_DIR}/${pdb_id}.mae"
    prepped_mae_file="${PREPWIZARD_OUTPUT_DIR}/${pdb_id}_prepped.mae"
    if [[ -f "$prepped_mae_file" ]]; then
        echo "Prepped MAE already exists for $pdb_id. Skipping."
        continue
    fi
    if [[ ! -f "$mae_file" ]]; then
        echo "MAE file not found for PDB ID $pdb_id at $mae_file"
        continue
    fi
    echo "Running PrepWizard for $pdb_id ..."
    timeout 1200 "$SCHRODINGER/utilities/prepwizard" \
        -WAIT \
        -JOBNAME "prep_${pdb_id}_$$" \
        "$mae_file" "$prepped_mae_file"
    exit_status=$?
    if [[ $exit_status -eq 124 ]]; then
        echo "PrepWizard timed out for $pdb_id (over 20 min)."
        echo "$pdb_id" >> "$PREPWIZARD_TIMEOUT_LOG"
        continue
    fi
    if [[ $exit_status -ne 0 ]]; then
        echo "PrepWizard failed for $pdb_id (exit code $exit_status)."
        echo "$pdb_id" >> "$PREPWIZARD_FAILED_LOG"
        continue
    fi
    if [[ ! -s "$prepped_mae_file" ]]; then
        echo "PrepWizard produced empty output for $pdb_id."
        echo "$pdb_id" >> "$PREPWIZARD_FAILED_LOG"
        continue
    fi
    echo "Successfully prepared $pdb_id"
done

# Collect Schrödinger PrepWizard job logs (prep_*.log) from the launch directory
PREPWIZARD_LOG_DIR="${WORK_DIR}/Enantiomer_protein_preperation/log"
mkdir -p "$PREPWIZARD_LOG_DIR"
moved_logs=0
shopt -s nullglob
for log_file in "${WORK_DIR}"/prep_*.log; do
    mv -v "$log_file" "$PREPWIZARD_LOG_DIR/"
    moved_logs=$((moved_logs + 1))
done
shopt -u nullglob

if [[ "$moved_logs" -gt 0 ]]; then
    echo "Moved $moved_logs prep_*.log file(s) to $PREPWIZARD_LOG_DIR"
else
    echo "No prep_*.log files found in $WORK_DIR"
fi
