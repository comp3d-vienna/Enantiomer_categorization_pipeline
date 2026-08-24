#!/bin/bash
# Pharmacophore generation (step 3.1) — interaction pharmacophores from aligned structures.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=../Data_preprocess/data_root.sh
source "${PROJECT_ROOT}/Data_preprocess/data_root.sh"
DATA_DIR="$(select_pipeline_data_dir)"
WORK_DIR="${DATA_DIR}"

ALIGNMENT_DIR="${WORK_DIR}/Alignment"
PHARMACOPHORE_DIR="${WORK_DIR}/Pharmacophore"
PAIRED_ENANTIOMERS_DIR="${ALIGNMENT_DIR}/paired_enantiomers"
ALIGNED_STRUCTURE_DIR="${ALIGNMENT_DIR}/paired_enantiomers_aligned_structures"
LIGAND_SDF_DIR="${ALIGNMENT_DIR}/Enantiomer_aligned_structure_ligandextract"
CANONICAL_SDF_DIR="${PHARMACOPHORE_DIR}/Enantiomer_aligned_structure_ligandextract_canonical"
OUTPUT_DIR="${PHARMACOPHORE_DIR}/Pharmacophore_generation_results"
INTERACTION_DATA_DIR="${PHARMACOPHORE_DIR}/Interaction_data"
LIGAND_LIST="${WORK_DIR}/db_ids_with_m_layer.txt"
PAIRED_LIGAND_IDS="${PHARMACOPHORE_DIR}/paired_ligand_ids.txt"

FAILED_JOBS_FILE="${PHARMACOPHORE_DIR}/generation_failed_jobs.txt"
FAILED_JOBS_RECEPTOR_NOT_FOUND="${PHARMACOPHORE_DIR}/generation_receptor_not_found.txt"
TIME_OUT_FILE="${PHARMACOPHORE_DIR}/generation_timeout.txt"

GEN_SCRIPT="${SCRIPT_DIR}/gen_ia_ph4s_fg.py"
CANON_SCRIPT="${SCRIPT_DIR}/canon_mols_1.py"
WATER_RESIDUES="HOH"
TIMEOUT_SEC=90

if [[ -n "${PYTHON_CMD:-}" ]]; then
    PYTHON="$PYTHON_CMD"
elif [[ -n "${PYTHON:-}" ]] && command -v "$PYTHON" >/dev/null 2>&1; then
    :
elif command -v python3 >/dev/null 2>&1; then
    PYTHON=python3
elif command -v python >/dev/null 2>&1; then
    PYTHON=python
else
    echo "Error: python3 or python not found in PATH." >&2
    echo "Activate categorize_pipeline or set PYTHON_CMD." >&2
    exit 1
fi
echo "Step 3.1 Python: $PYTHON"

parse_complex_pdb_chain() {
    local fname="$1"
    COMPLEX_PDB_ID=""
    COMPLEX_CHAIN_ID=""
    if [[ "$fname" =~ _([a-zA-Z0-9]+)_chain_([a-zA-Z0-9])(_copied)?\.pdb$ ]]; then
        COMPLEX_PDB_ID="${BASH_REMATCH[1]}"
        COMPLEX_CHAIN_ID="${BASH_REMATCH[2]}"
        return 0
    fi
    if [[ "$fname" =~ ^([a-zA-Z0-9]+)_prepped(-[0-9]+)?\.pdb$ ]]; then
        COMPLEX_PDB_ID="${BASH_REMATCH[1]}"
        return 0
    fi
    return 1
}

ligand_id_in_paired_set() {
    local line="$1"
    grep -qxF "$line" "$PAIRED_LIGAND_IDS"
}

build_paired_ligand_ids() {
  "$PYTHON" - <<'PY'
import csv
import glob
import os
import sys

paired_dir = os.environ["PAIRED_ENANTIOMERS_DIR"]
out_path = os.environ["PAIRED_LIGAND_IDS"]
ids = set()
for path in sorted(glob.glob(os.path.join(paired_dir, "*.csv"))):
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            ident = row["identifier"].strip()
            m_layer = row.get("m_layer", "").strip()
            ids.add(ident)
            if m_layer:
                ids.add(f"{ident}_{m_layer}")
with open(out_path, "w") as f:
    for item in sorted(ids):
        f.write(item + "\n")
print(f"Wrote {len(ids)} paired ligand IDs to {out_path}")
PY
}

if [[ ! -d "$PAIRED_ENANTIOMERS_DIR" ]]; then
    echo "Error: paired enantiomer groups not found: $PAIRED_ENANTIOMERS_DIR" >&2
    echo "Run Alignment steps 2.1-2.4 first." >&2
    exit 1
fi

if [[ ! -d "$ALIGNED_STRUCTURE_DIR" ]]; then
    echo "Error: aligned structures not found: $ALIGNED_STRUCTURE_DIR" >&2
    exit 1
fi

if [[ ! -f "$LIGAND_LIST" ]]; then
    echo "Error: ligand list not found: $LIGAND_LIST" >&2
    exit 1
fi

mkdir -p "$CANONICAL_SDF_DIR" "$OUTPUT_DIR" "$INTERACTION_DATA_DIR" "$PHARMACOPHORE_DIR"
: > "$FAILED_JOBS_FILE"
: > "$FAILED_JOBS_RECEPTOR_NOT_FOUND"
: > "$TIME_OUT_FILE"

export PAIRED_ENANTIOMERS_DIR PAIRED_LIGAND_IDS
build_paired_ligand_ids

echo "Aligned structures: $ALIGNED_STRUCTURE_DIR"
echo "Ligand SDF dir:       $LIGAND_SDF_DIR"
echo "Pharmacophore output: $PHARMACOPHORE_DIR"

echo "---------------------------------------------"
echo "Step 3.1a: Canonicalize ligands (paired only)"
echo "---------------------------------------------"
total_start_time=$(date +%s)
shopt -s nullglob
for ligand_file in "${LIGAND_SDF_DIR}"/*.sdf; do
    base=$(basename "$ligand_file" .sdf)
    if ! ligand_id_in_paired_set "$base"; then
        continue
    fi
    canonical_ligand_file="${CANONICAL_SDF_DIR}/${base}_canon.sdf"
    if [[ -e "$canonical_ligand_file" ]]; then
        continue
    fi
    "$PYTHON" "$CANON_SCRIPT" -i "$ligand_file" -o "$canonical_ligand_file"
done
shopt -u nullglob
total_end_time=$(date +%s)
echo "Canonicalization time: $((total_end_time - total_start_time)) seconds."

echo "---------------------------------------------"
echo "Step 3.1b: Generate interaction pharmacophores"
echo "---------------------------------------------"

while IFS= read -r line; do
    line=$(echo "$line" | tr -d '[:space:]')
    [[ -z "$line" ]] && continue

    if ! ligand_id_in_paired_set "$line"; then
        continue
    fi

    IFS='_' read -r pdbid ligand_id chain_id residue_number rest <<< "$line"
    if [[ -z "$pdbid" || -z "$ligand_id" || -z "$chain_id" || -z "$residue_number" ]]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') Could not parse '$line'" | tee -a "$FAILED_JOBS_RECEPTOR_NOT_FOUND"
        continue
    fi

    receptor_files=()
    while IFS= read -r file; do
        matched_name=$(basename "$file")
        if parse_complex_pdb_chain "$matched_name"; then
            if [[ "$COMPLEX_PDB_ID" == "$pdbid" && "$COMPLEX_CHAIN_ID" == "$chain_id" ]]; then
                receptor_files+=("$file")
            fi
        fi
    done < <(find "$ALIGNED_STRUCTURE_DIR" -type f \( \
        -name "*_${pdbid}_chain_${chain_id}.pdb" \
        -o -name "*_${pdbid}_chain_${chain_id}_copied.pdb" \
        \) 2>/dev/null)

    if [[ ${#receptor_files[@]} -eq 0 ]]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') No aligned receptor for pdb=$pdbid chain=$chain_id" \
            | tee -a "$FAILED_JOBS_RECEPTOR_NOT_FOUND"
        continue
    fi

    ligand_file="${CANONICAL_SDF_DIR}/${line}_canon.sdf"
    if [[ ! -f "$ligand_file" ]]; then
        shopt -s nullglob
        matches=( "${CANONICAL_SDF_DIR}/${line}"*.sdf )
        shopt -u nullglob
        if [[ ${#matches[@]} -eq 0 || ! -f "${matches[0]}" ]]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') Ligand SDF not found for $line (expected $ligand_file)" \
                | tee -a "$FAILED_JOBS_FILE"
            continue
        fi
        ligand_file="${matches[0]}"
    fi
    if [[ ! -s "$ligand_file" ]]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') Empty ligand SDF for $line: $ligand_file" \
            | tee -a "$FAILED_JOBS_FILE"
        continue
    fi

    ligand_basename=$(basename "$ligand_file" .sdf)
    ligand_basename="${ligand_basename/_canon/}"
    output_file="${OUTPUT_DIR}/ph4_${ligand_basename}.pml"
    interaction_data_file="${INTERACTION_DATA_DIR}/interaction_${ligand_basename}.tsv"
    strip_res="${chain_id}_${ligand_id}_${residue_number}"

    if [[ -e "$output_file" && -e "$interaction_data_file" ]]; then
        continue
    fi

    for receptor_file in "${receptor_files[@]}"; do
        if [[ ! -s "$receptor_file" ]]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') Empty receptor $receptor_file" \
                | tee -a "$FAILED_JOBS_RECEPTOR_NOT_FOUND"
            continue
        fi

        echo "Processing $line with receptor $(basename "$receptor_file")"
        start_time=$(date +%s)
        if timeout "$TIMEOUT_SEC" "$PYTHON" "$GEN_SCRIPT" \
            -r "$receptor_file" \
            -l "$ligand_file" \
            -o "$output_file" \
            -s "$strip_res" \
            "$WATER_RESIDUES" \
            -i "$interaction_data_file"; then
            echo "Finished $line in $(( $(date +%s) - start_time )) seconds."
        else
            exit_status=$?
            if [[ $exit_status -eq 124 ]]; then
                echo "$(date '+%Y-%m-%d %H:%M:%S') Timeout for $line receptor $receptor_file" \
                    | tee -a "$TIME_OUT_FILE"
            else
                echo "$(date '+%Y-%m-%d %H:%M:%S') Failed for $line receptor $receptor_file (exit $exit_status)" \
                    | tee -a "$FAILED_JOBS_FILE"
            fi
        fi
    done
done < "$LIGAND_LIST"

echo "---------------------------------------------"
echo "Pharmacophore generation complete"
echo "  Interaction data: $INTERACTION_DATA_DIR"
echo "  Pharmacophores:   $OUTPUT_DIR"
echo "---------------------------------------------"
