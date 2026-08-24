#!/bin/bash
# Alignment step 2.3 — Schrödinger structalign on parsed chain structures.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=../Data_preprocess/data_root.sh
source "${PROJECT_ROOT}/Data_preprocess/data_root.sh"
DATA_DIR="$(select_pipeline_data_dir)"
WORK_DIR="${DATA_DIR}"

ALIGNMENT_DIR="${WORK_DIR}/Alignment"
PREP_PDB_DIR="${ALIGNMENT_DIR}/Enantiomer_prep_pdbformat"
PARSED_CHAIN_DIR="${ALIGNMENT_DIR}/paired_enantiomers_parsed_chain"
ALIGNED_DIR="${ALIGNMENT_DIR}/paired_enantiomers_aligned_structures"
LIST_WITH_UNIPROT_ID="${PARSED_CHAIN_DIR}/list_with_uniprot_id.txt"
REFERENCE_OUTPUT_FILE="${ALIGNED_DIR}/selected_reference_structures.txt"
ERROR_LOG_FILE="${ALIGNED_DIR}/errors_alignment.log"

if [[ -z "${SCHROEDINGER_UTILITIES_DIR:-}" ]]; then
    if [[ -z "${SCHRODINGER:-}" ]]; then
        echo "Error: SCHRODINGER is not set." >&2
        echo "Set SCHRODINGER to your Schrödinger install root (the directory that contains utilities/structalign)." >&2
        exit 1
    fi
    SCHROEDINGER_UTILITIES_DIR="${SCHRODINGER}/utilities"
fi

STRUCTALIGN="${SCHROEDINGER_UTILITIES_DIR}/structalign"

pdb_id_from_stem() {
    local stem="$1"
    if [[ "$stem" =~ ^g[0-9]{6}_[a-f0-9]{10}_([a-zA-Z0-9]+)_chain_ ]]; then
        echo "${BASH_REMATCH[1]}"
    fi
}

find_prep_pdb() {
    local pdb_id="$1"
    local match
    match=$(ls "${PREP_PDB_DIR}/${pdb_id}_prepped"*.pdb 2>/dev/null | head -n 1 || true)
    echo "$match"
}

if [[ ! -x "$STRUCTALIGN" ]]; then
    echo "Error: structalign not found or not executable: $STRUCTALIGN" >&2
    echo "Set SCHRODINGER or SCHROEDINGER_UTILITIES_DIR." >&2
    exit 1
fi

if [[ ! -f "$LIST_WITH_UNIPROT_ID" ]]; then
    echo "Error: list file not found: $LIST_WITH_UNIPROT_ID" >&2
    echo "Run Alignment/group_and_parse_chains.py first (steps 2.1-2.2)." >&2
    exit 1
fi

mkdir -p "$ALIGNED_DIR"
cd "$ALIGNMENT_DIR"
: > "$REFERENCE_OUTPUT_FILE"
: > "$ERROR_LOG_FILE"

echo "Alignment working directory: $(pwd)"
echo "Parsed chains:  $PARSED_CHAIN_DIR"
echo "Prep PDB dir:   $PREP_PDB_DIR"
echo "Aligned output: $ALIGNED_DIR"

process_group() {
    local target_key="$1"
    shift
    local group_items=("$@")

    echo "[GROUP] $target_key"
    echo "  Members:"
    for item in "${group_items[@]}"; do
        echo "    - $item"
    done

    local existing_files=()
    local structalign_input=()

    for p_fname in "${group_items[@]}"; do
        local pdb_path="${PARSED_CHAIN_DIR}/${p_fname}.pdb"
        if [[ -f "$pdb_path" ]]; then
            existing_files+=("$p_fname")
            structalign_input+=("$pdb_path")
        else
            echo "[Warning] Missing parsed chain PDB: $pdb_path"
            echo "[Warning] Missing parsed chain PDB: $pdb_path" >> "$ERROR_LOG_FILE"
        fi
    done

    local count=${#existing_files[@]}

    if (( count == 0 )); then
        echo "[Info] No PDB files for group ${target_key}; skipping."
        echo "[Info] No PDB files for group ${target_key}; skipping." >> "$ERROR_LOG_FILE"
        return
    fi

    if (( count == 1 )); then
        local first_file="${existing_files[0]}"
        local pdb_id
        pdb_id=$(pdb_id_from_stem "$first_file")
        echo "[Info] Single structure in group ${target_key}; copying."
        if [[ -n "$pdb_id" ]]; then
            local prep_pdb
            prep_pdb=$(find_prep_pdb "$pdb_id")
            if [[ -n "$prep_pdb" && -f "$prep_pdb" ]]; then
                cp "$prep_pdb" "${ALIGNED_DIR}/${pdb_id}_prepped_copied.pdb" 2>>"$ERROR_LOG_FILE" || \
                    echo "[ERROR] Failed to copy prepped PDB for ${pdb_id}" >> "$ERROR_LOG_FILE"
            else
                echo "[Warning] Missing prepped PDB for ${pdb_id}" >> "$ERROR_LOG_FILE"
            fi
        fi
        cp "${PARSED_CHAIN_DIR}/${first_file}.pdb" "${ALIGNED_DIR}/${first_file}_copied.pdb" 2>>"$ERROR_LOG_FILE" || \
            echo "[ERROR] Failed to copy parsed chain PDB for ${first_file}" >> "$ERROR_LOG_FILE"
        return
    fi

    local ref_val="${existing_files[0]}"
    echo "$ref_val" >> "$REFERENCE_OUTPUT_FILE"
    echo "[Info] Aligning ${count} structures (reference: ${ref_val})"

    local unique_pdb_ids=()
    for fname in "${existing_files[@]}"; do
        local pdb_id
        pdb_id=$(pdb_id_from_stem "$fname")
        if [[ -n "$pdb_id" && ! " ${unique_pdb_ids[*]:-} " =~ " ${pdb_id} " ]]; then
            unique_pdb_ids+=("$pdb_id")
        fi
    done

    if (( ${#unique_pdb_ids[@]} == 1 )); then
        local pdb_id="${unique_pdb_ids[0]}"
        echo "[Info] Same PDB ID (${pdb_id}); copying reference prepped structure."
        local prep_pdb
        prep_pdb=$(find_prep_pdb "$pdb_id")
        if [[ -n "$prep_pdb" && -f "$prep_pdb" ]]; then
            cp "$prep_pdb" "${ALIGNED_DIR}/${pdb_id}_prepped_copied.pdb" 2>>"$ERROR_LOG_FILE" || \
                echo "[ERROR] Failed to copy prepped PDB for ${pdb_id}" >> "$ERROR_LOG_FILE"
        else
            echo "[Warning] Missing prepped PDB for ${pdb_id}" >> "$ERROR_LOG_FILE"
        fi
    else
        echo "[Info] Multiple PDB IDs (${unique_pdb_ids[*]}); aligning full prepped structures."
        local structalign_pdb_input=()
        local ref_pdb_id
        ref_pdb_id=$(pdb_id_from_stem "$ref_val")
        local ref_pdb_file
        ref_pdb_file=$(find_prep_pdb "$ref_pdb_id")

        for pdb_id in "${unique_pdb_ids[@]}"; do
            local this_pdb_file
            this_pdb_file=$(find_prep_pdb "$pdb_id")
            if [[ -n "$this_pdb_file" && -f "$this_pdb_file" ]]; then
                structalign_pdb_input+=("$this_pdb_file")
            else
                echo "[Warning] Missing prepped PDB for ${pdb_id}" >> "$ERROR_LOG_FILE"
            fi
        done

        if [[ -n "$ref_pdb_file" && -f "$ref_pdb_file" ]]; then
            structalign_pdb_input=("$ref_pdb_file" $(for f in "${structalign_pdb_input[@]}"; do [[ "$f" != "$ref_pdb_file" ]] && echo "$f"; done))
            if "$STRUCTALIGN" "${structalign_pdb_input[@]}"; then
                echo "[SUCCESS] Full-structure alignment completed for ${target_key}"
            else
                echo "[ERROR] Full-structure alignment failed for ${target_key}" >> "$ERROR_LOG_FILE"
            fi
        else
            echo "[ERROR] Reference prepped PDB not found for ${ref_pdb_id}" >> "$ERROR_LOG_FILE"
        fi
    fi

    if "$STRUCTALIGN" "${structalign_input[@]}"; then
        echo "[SUCCESS] Chain alignment completed for ${target_key}"
    else
        echo "[ERROR] Chain alignment failed for ${target_key}" >> "$ERROR_LOG_FILE"
    fi
}

tmp_group_file=$(mktemp)
while IFS=$'\t' read -r u_id i_hash p_id c_id g_id || [[ -n "$u_id" ]]; do
    [[ -z "$u_id" || "$u_id" =~ ^# ]] && continue
    u_id="${u_id//$'\r'/}"
    i_hash="${i_hash//$'\r'/}"
    p_id="${p_id//$'\r'/}"
    c_id="${c_id//$'\r'/}"
    g_id="${g_id//$'\r'/}"

    if [[ -z "$g_id" ]]; then
        echo "[Warning] Missing group_id for ${u_id} ${i_hash} ${p_id} ${c_id}; skipping" >> "$ERROR_LOG_FILE"
        continue
    fi

    f_name="${g_id}_${i_hash}_${p_id}_chain_${c_id}"
    g_key="${g_id}_${i_hash}"
    printf "%s\t%s\n" "$g_key" "$f_name" >> "$tmp_group_file"
done < "$LIST_WITH_UNIPROT_ID"

current_group_key=""
collected_files=()

while IFS=$'\t' read -r row_key row_fname; do
    if [[ -z "$current_group_key" ]]; then
        current_group_key="$row_key"
    fi

    if [[ "$row_key" != "$current_group_key" ]]; then
        process_group "$current_group_key" "${collected_files[@]}"
        collected_files=()
        current_group_key="$row_key"
    fi

    collected_files+=("$row_fname")
done < <(LC_ALL=C sort -t $'\t' -k1,1 "$tmp_group_file")

if [[ ${#collected_files[@]} -gt 0 ]]; then
    process_group "$current_group_key" "${collected_files[@]}"
fi

rm -f "$tmp_group_file"

find "$ALIGNMENT_DIR" -maxdepth 1 -type f -name "*.pdb" | while read -r pdbfile; do
    mv "$pdbfile" "$ALIGNED_DIR/"
done

echo "[Info] Alignment outputs collected in $ALIGNED_DIR"
