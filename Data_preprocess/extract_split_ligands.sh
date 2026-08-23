#!/bin/bash
# Combined script: Extract ligand structure from PDB, split multi-component SDF files
#
# Environment:
#   PIPELINE_TEST=1   Route outputs under ../test/
#   MMCIF_SOURCE      Gzip mmCIF/PDB archive
#                     (production default: ../Data/mmCIF;
#                      test default: ../test/mmcif)
#   MMCIF_RENAME      Decompressed structures
#                     (production default: ../Data/mmCIF_rename;
#                      test default: ../test/mmcif_rename)
#
# Resume: existing non-empty SDF and decompressed structure files are skipped
# per PDB. A PDB is not re-split if any separated SDF for that code already exists.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$(cd "$SCRIPT_DIR/../Data" && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/../test" && pwd)"

if [[ "${PIPELINE_TEST:-}" == "1" ]]; then
    # TEST (run test/run_data_preprocess_test.sh or export PIPELINE_TEST=1)
    INPUT_DIR="${MMCIF_SOURCE:-${TEST_DIR}/mmcif}"
    PDB_RENAME_DIR="${MMCIF_RENAME:-${TEST_DIR}/mmcif_rename}"
    LIGAND_SDF_DIR="${TEST_DIR}/ligand_sdf"
    SEPARATED_SDF_DIR="${TEST_DIR}/ligand_sdf_separated"
    TIMER_DIR="${TEST_DIR}/timer_mmCIF"
    log_file="${TEST_DIR}/extract_split_procedure_errors.log"
else
    INPUT_DIR="${MMCIF_SOURCE:-${DATA_DIR}/mmCIF}"
    PDB_RENAME_DIR="${MMCIF_RENAME:-${DATA_DIR}/mmCIF_rename}"
    LIGAND_SDF_DIR="${DATA_DIR}/ligand_sdf_structure_from_mmCIF"
    SEPARATED_SDF_DIR="${LIGAND_SDF_DIR}_separated"
    TIMER_DIR="${DATA_DIR}/timer_mmCIF"
    log_file="${DATA_DIR}/extract_split_procedure_errors.log"
fi

mkdir -p "$TIMER_DIR"
mkdir -p "$PDB_RENAME_DIR"
mkdir -p "$LIGAND_SDF_DIR"
mkdir -p "$SEPARATED_SDF_DIR"

echo "mmCIF source:     $INPUT_DIR"
echo "Decompressed CIF: $PDB_RENAME_DIR"
echo "Ligand SDFs:      $LIGAND_SDF_DIR"
echo "Separated SDFs:   $SEPARATED_SDF_DIR"

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: mmCIF source not found: $INPUT_DIR" >&2
        echo "Place downloaded mmCIF files in Data/mmCIF, or set MMCIF_SOURCE." >&2
    exit 1
fi

# Initialize combined error log
echo "=== Combined Error Log - $(date) ===" > "$log_file"
echo "" >> "$log_file"

# Record overall start time
SCRIPT_START_TIME=$(date +%s)
SCRIPT_START_DATE=$(date)
echo "=== Overall Script Started - $SCRIPT_START_DATE ===" >> ${TIMER_DIR}/extract_ligand_from_cif.txt

# ============================================================================
# Phase 1: Extracting ligands from PDB files
# ============================================================================
echo "=========================================="
echo "Phase 1: Extracting ligands from PDB files"
echo "=========================================="

EXTRACTION_START_TIME=0
EXTRACTION_END_TIME=0
EXTRACTION_DURATION=0
PHASE1_DONE=0
PHASE1_SKIPPED=0

# Record start time for extraction phase
echo "=== Extraction Phase Started - $(date) ===" >> ${TIMER_DIR}/extract_ligand_from_cif.txt
date >> ${TIMER_DIR}/extract_ligand_from_cif.txt
EXTRACTION_START_TIME=$(date +%s)

# Process all .gz files recursively. Process substitution keeps counters in this shell.
while IFS= read -r gz_file; do
    echo "Processing: $gz_file"

    file_name=$(basename "$gz_file")

    # Extract 4-digit PDB code
    # Handles files like: 9hhr.cif.gz, 1abc.cif.gz, pdb9kju.ent.gz
    if [[ "$file_name" =~ ^([0-9a-zA-Z]{4})\.cif\.gz$ ]]; then
        pdb_code="${BASH_REMATCH[1]}"
        new_name="${PDB_RENAME_DIR}/${pdb_code}.cif"
        input_format="cif"
    elif [[ "$file_name" =~ ^pdb([0-9a-zA-Z]{4})\.ent\.gz$ ]]; then
        pdb_code="${BASH_REMATCH[1]}"
        new_name="${PDB_RENAME_DIR}/${pdb_code}.pdb"
        input_format="pdb"
    else
        msg="Warning: Could not extract PDB code from $gz_file"
        echo "$(date): $msg" >> "$log_file"
        echo "---"
        continue
    fi

    echo "PDB code: $pdb_code (input: $input_format)"

    sdf_out="${LIGAND_SDF_DIR}/${pdb_code}.sdf"
    need_rename=false
    need_sdf=false
    [[ ! -s "$new_name" ]] && need_rename=true
    [[ ! -s "$sdf_out" ]] && need_sdf=true

    if [[ "$need_rename" == false && "$need_sdf" == false ]]; then
        PHASE1_SKIPPED=$((PHASE1_SKIPPED + 1))
        echo "Skipping $pdb_code: decompressed structure and SDF already exist"
        echo "---"
        continue
    fi

    if [[ "$need_rename" == true ]]; then
        if ! gunzip -c "$gz_file" > "$new_name"; then
            msg="Warning: Failed to extract $gz_file"
            echo "$(date): $msg" >> "$log_file"
            echo "---"
            continue
        fi
    fi

    if [[ "$need_sdf" == true ]]; then
        echo "Running unicon to convert to SDF"
        ./unicon_1.5.0/unicon -i "$new_name" --multi -o "$sdf_out"

        if [ $? -eq 0 ]; then
            echo "Successfully created $sdf_out"
            PHASE1_DONE=$((PHASE1_DONE + 1))
        else
            msg="Error running unicon for $new_name (SDF conversion)"
            echo "$(date): $msg" >> "$log_file"
        fi
    else
        PHASE1_DONE=$((PHASE1_DONE + 1))
        echo "SDF already exists for $pdb_code; refreshed decompressed structure"
    fi

    echo "---"
done < <(find "$INPUT_DIR" -name "*.gz" -type f)

# Record end time for extraction phase
EXTRACTION_END_TIME=$(date +%s)
EXTRACTION_DURATION=$((EXTRACTION_END_TIME - EXTRACTION_START_TIME))
echo "=== Extraction Phase Ended - $(date) ===" >> ${TIMER_DIR}/extract_ligand_from_cif.txt
date >> ${TIMER_DIR}/extract_ligand_from_cif.txt
echo "=== Extraction Phase Duration: ${EXTRACTION_DURATION} seconds ===" >> ${TIMER_DIR}/extract_ligand_from_cif.txt
echo "=== Extraction Phase extracted=${PHASE1_DONE} skipped=${PHASE1_SKIPPED} ===" >> ${TIMER_DIR}/extract_ligand_from_cif.txt

# ============================================================================
# Phase 2: Splitting multi-component SDF files
# ============================================================================
echo ""
echo "=========================================="
echo "Phase 2: Splitting multi-component SDF files"
echo "=========================================="

SPLITTING_START_TIME=0
SPLITTING_END_TIME=0
SPLITTING_DURATION=0
TOTAL_LIGANDS=0
TOTAL_FILES_PROCESSED=0
TOTAL_FILES_EMPTY=0
TOTAL_FILES_NO_LIGANDS=0
PHASE2_SKIPPED=0
errors=()

# Record start time for splitting phase
SPLITTING_START_TIME=$(date +%s)
echo "=== Splitting Phase Started - $(date) ===" >> ${TIMER_DIR}/extract_ligand_from_cif.txt
date >> ${TIMER_DIR}/extract_ligand_from_cif.txt

shopt -s nullglob
for INPUT_FILE in "$LIGAND_SDF_DIR"/*.sdf; do
    TOTAL_FILES_PROCESSED=$((TOTAL_FILES_PROCESSED + 1))
    LIGAND_NUM=0
    FIRST_LINE=""
    CURRENT_CONTENT=""
    INVALID_MARKERS=0

    BASE_NAME=$(basename "$INPUT_FILE" .sdf)
    echo "Processing: $BASE_NAME.sdf"

    existing_splits=("${SEPARATED_SDF_DIR}/${BASE_NAME}_"*.sdf)
    if [[ ${#existing_splits[@]} -gt 0 ]]; then
        PHASE2_SKIPPED=$((PHASE2_SKIPPED + 1))
        echo "  Skipping $BASE_NAME: ${#existing_splits[@]} separated file(s) already exist"
        echo "---"
        continue
    fi

    # Check if input file is empty
    if [ ! -s "$INPUT_FILE" ]; then
        TOTAL_FILES_EMPTY=$((TOTAL_FILES_EMPTY + 1))
        msg="Error: $INPUT_FILE is empty"
        errors+=("$msg")
        echo "$(date): $msg" >> "$log_file"
        continue
    fi

    # Process file line by line
    while IFS= read -r line || [ -n "$line" ]; do
        # Capture first non-empty line as molecule ID
        if [ -z "$FIRST_LINE" ] && [ -n "$line" ]; then
            FIRST_LINE="$line"
            IN_MOLECULE=true
        fi

        # Accumulate molecule content
        CURRENT_CONTENT="${CURRENT_CONTENT}${line}"$'\n'

        # Detect end of molecule
        if [ "$line" = '$$$$' ]; then
            if [ -n "$FIRST_LINE" ]; then
                LIGAND_NUM=$((LIGAND_NUM + 1))
                TOTAL_LIGANDS=$((TOTAL_LIGANDS + 1))

                SAFE_NAME=$(echo "$FIRST_LINE" | tr -cd '[:alnum:]_-' | cut -c1-50)
                OUTPUT_FILE="${SEPARATED_SDF_DIR}/${BASE_NAME}_${SAFE_NAME}.sdf"

                if ! printf "%s" "$CURRENT_CONTENT" > "$OUTPUT_FILE"; then
                    msg="Error: Failed to write to $OUTPUT_FILE"
                    errors+=("$msg")
                    echo "$(date): $msg" >> "$log_file"
                else
                    echo "  Created: $(basename "$OUTPUT_FILE")"
                fi

                # Reset for next molecule
                FIRST_LINE=""
                CURRENT_CONTENT=""
                IN_MOLECULE=false
            else
                INVALID_MARKERS=$((INVALID_MARKERS + 1))
                msg="Warning: Found end marker (\$\$\$\$) without valid header in $INPUT_FILE"
                errors+=("$msg")
                echo "$(date): $msg" >> "$log_file"
            fi
        fi
    done < "$INPUT_FILE"

    # Handle case where file doesn't end with $$$$
    if [ -n "$CURRENT_CONTENT" ] && [ -n "$FIRST_LINE" ]; then
        LIGAND_NUM=$((LIGAND_NUM + 1))
        TOTAL_LIGANDS=$((TOTAL_LIGANDS + 1))
        SAFE_NAME=$(echo "$FIRST_LINE" | tr -cd '[:alnum:]_-' | cut -c1-50)
        OUTPUT_FILE="${SEPARATED_SDF_DIR}/${BASE_NAME}_${SAFE_NAME}.sdf"
        if ! printf "%s" "$CURRENT_CONTENT" > "$OUTPUT_FILE"; then
            msg="Error: Failed to write final molecule to $OUTPUT_FILE"
            errors+=("$msg")
            echo "$(date): $msg" >> "$log_file"
        else
            echo "  Created: $(basename "$OUTPUT_FILE") (final molecule without $$$$ marker)"
        fi
    fi

    # Report statistics for this file
    if [ $LIGAND_NUM -eq 0 ]; then
        TOTAL_FILES_NO_LIGANDS=$((TOTAL_FILES_NO_LIGANDS + 1))
        msg="Warning: No ligands found in $INPUT_FILE"
        errors+=("$msg")
        echo "$(date): $msg" >> "$log_file"
    else
        echo "  Summary: Extracted $LIGAND_NUM ligand(s) from $BASE_NAME.sdf"
        if [ $INVALID_MARKERS -gt 0 ]; then
            echo "  Warning: Found $INVALID_MARKERS invalid marker(s) in $BASE_NAME.sdf"
        fi
    fi

    echo "---"
done
shopt -u nullglob

# Record end time for splitting phase
SPLITTING_END_TIME=$(date +%s)
SPLITTING_DURATION=$((SPLITTING_END_TIME - SPLITTING_START_TIME))
echo "=== Splitting Phase Ended - $(date) ===" >> ${TIMER_DIR}/extract_ligand_from_cif.txt
date >> ${TIMER_DIR}/extract_ligand_from_cif.txt
echo "=== Splitting Phase Duration: ${SPLITTING_DURATION} seconds ===" >> ${TIMER_DIR}/extract_ligand_from_cif.txt
echo "=== Splitting Phase Statistics ===" >> ${TIMER_DIR}/extract_ligand_from_cif.txt
echo "Total files processed: $TOTAL_FILES_PROCESSED" >> ${TIMER_DIR}/extract_ligand_from_cif.txt
echo "Files skipped (already split): $PHASE2_SKIPPED" >> ${TIMER_DIR}/extract_ligand_from_cif.txt
echo "Total ligands extracted: $TOTAL_LIGANDS" >> ${TIMER_DIR}/extract_ligand_from_cif.txt
echo "Empty files: $TOTAL_FILES_EMPTY" >> ${TIMER_DIR}/extract_ligand_from_cif.txt
echo "Files with no ligands: $TOTAL_FILES_NO_LIGANDS" >> ${TIMER_DIR}/extract_ligand_from_cif.txt
echo "Errors encountered: ${#errors[@]}" >> ${TIMER_DIR}/extract_ligand_from_cif.txt

# ============================================================================
# Final Summary
# ============================================================================
# Calculate overall script duration
SCRIPT_END_TIME=$(date +%s)
SCRIPT_END_DATE=$(date)
SCRIPT_DURATION=$((SCRIPT_END_TIME - SCRIPT_START_TIME))
echo "=== Overall Script Ended - $SCRIPT_END_DATE ===" >> ${TIMER_DIR}/extract_ligand_from_cif.txt
echo "=== Overall Script Duration: ${SCRIPT_DURATION} seconds ===" >> ${TIMER_DIR}/extract_ligand_from_cif.txt

echo ""
echo "=========================================="
echo "Processing Complete"
echo "=========================================="
echo "Phase 1: extracted ${PHASE1_DONE}, skipped ${PHASE1_SKIPPED} (already present)"
echo "Phase 2 Statistics:"
echo "  Total files processed: $TOTAL_FILES_PROCESSED"
echo "  Files skipped (already split): $PHASE2_SKIPPED"
echo "  Total ligands separated: $TOTAL_LIGANDS"
echo "  Empty files: $TOTAL_FILES_EMPTY"
echo "  Files with no ligands: $TOTAL_FILES_NO_LIGANDS"
echo "Separated SDF files saved to: $SEPARATED_SDF_DIR"

# Format duration as HH:MM:SS
format_duration() {
    local total_seconds=$1
    local hours=$((total_seconds / 3600))
    local minutes=$(((total_seconds % 3600) / 60))
    local seconds=$((total_seconds % 60))
    if [ $hours -gt 0 ]; then
        printf "%d:%02d:%02d (hours:minutes:seconds)" $hours $minutes $seconds
    elif [ $minutes -gt 0 ]; then
        printf "%d:%02d (minutes:seconds)" $minutes $seconds
    else
        printf "%d seconds" $seconds
    fi
}

echo ""
echo "=== Time Cost Summary ==="
echo "Start time: $SCRIPT_START_DATE"
echo "End time: $SCRIPT_END_DATE"
echo "Extraction phase: $(format_duration $EXTRACTION_DURATION)"
echo "Splitting phase: $(format_duration $SPLITTING_DURATION)"
echo "Total time: $(format_duration $SCRIPT_DURATION)"
echo ""

if [ ${#errors[@]} -gt 0 ]; then
    echo "=== Errors encountered (${#errors[@]} total) ==="
    echo "See $log_file for details."
else
    echo "No errors encountered!"
    echo "$(date): No errors encountered!" >> "$log_file"
fi

echo ""
echo "See $log_file for detailed error information."
