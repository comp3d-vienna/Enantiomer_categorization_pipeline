#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=data_root.sh
source "${PROJECT_ROOT}/Data_preprocess/data_root.sh"
DATA_DIR="$(select_pipeline_data_dir)"
WORK_DIR="${DATA_DIR}"

if [[ "${PIPELINE_TEST:-}" == "1" ]]; then
    CIF_DIR="${MMCIF_RENAME:-${WORK_DIR}/mmcif_rename}"
else
    CIF_DIR="${MMCIF_RENAME:-${DATA_DIR}/mmCIF_rename}"
fi
Enantiomer_list_identifier="${WORK_DIR}/db_ids_with_m_layer.txt"
PREPWIZARD_TIMEOUT_LOG="${WORK_DIR}/prepwizard_timeout_pdbids.txt"
PREPWIZARD_FAILED_LOG="${WORK_DIR}/prepwizard_failed_pdbids.txt"

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

# parallelised protein prep
MAX_JOBS=$(( $(nproc) - 2 ))
# maximum no. of jobs is cpu dependent

run_prepwizard() {
    pdb_id="$1"
    mae_file="${MAE_OUTPUT_DIR}/${pdb_id}.mae"
    prepped_mae_file="${PREPWIZARD_OUTPUT_DIR}/${pdb_id}_prepped.mae"
    jobname="prep_${pdb_id}_$$"   # $$ is this worker's PID, so the name is unique per PDB

    if [[ -f "$prepped_mae_file" ]]; then
        echo "Prepped MAE already exists for $pdb_id. Skipping."
        return 0
    fi
    if [[ ! -f "$mae_file" ]]; then
        echo "MAE file not found for PDB ID $pdb_id at $mae_file"
        return 0
    fi

    echo "Running PrepWizard for $pdb_id ..."
    timeout 1200 "$SCHRODINGER/utilities/prepwizard" \
        -WAIT \
        -JOBNAME "$jobname" \
        "$mae_file" "$prepped_mae_file"
    exit_status=$?

    if [[ $exit_status -eq 124 ]]; then
        echo "PrepWizard timed out for $pdb_id (over 20 min)."
        echo "$pdb_id" >> "$PREPWIZARD_TIMEOUT_LOG"
        schroedinger_job_id=$("$SCHRODINGER/jobcontrol" -list | grep "$jobname" | awk '{ print $1 }')
        "$SCHRODINGER/jobcontrol" -kill "$schroedinger_job_id"
        sleep 5
        pkill -f "$jobname"          # fallback for stray child processes
        rm -f "$prepped_mae_file"    # don't leave a partial output
        return 0
    fi
    if [[ $exit_status -ne 0 ]]; then
        echo "PrepWizard failed for $pdb_id (exit code $exit_status)."
        echo "$pdb_id" >> "$PREPWIZARD_FAILED_LOG"
        return 0
    fi
    if [[ ! -s "$prepped_mae_file" ]]; then
        echo "PrepWizard produced empty output for $pdb_id."
        echo "$pdb_id" >> "$PREPWIZARD_FAILED_LOG"
        return 0
    fi
    echo "Successfully prepared $pdb_id"
}

# workers are fresh bash processes, so the function and variables must be exported
export -f run_prepwizard
export SCHRODINGER MAE_OUTPUT_DIR PREPWIZARD_OUTPUT_DIR PREPWIZARD_TIMEOUT_LOG PREPWIZARD_FAILED_LOG

printf '%s\n' "${unique_pdb_id_list[@]}" | xargs -P "$MAX_JOBS" -I{} bash -c 'run_prepwizard "$@"' _ {}

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
