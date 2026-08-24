# Shared data-root selection. Requires PROJECT_ROOT (absolute).
#
#   Data/            default (user mmCIF in Data/mmCIF)
#   Data/test/       when MMCIF_SOURCE is Data/test/mmCIF
#   PIPELINE_DATA_DIR overrides both

pipeline_abs_path() {
    local p="$1"
    if [[ "$p" != /* ]]; then
        p="${PROJECT_ROOT}/${p}"
    fi
    realpath -m "$p"
}

select_pipeline_data_dir() {
    if [[ "${PIPELINE_TEST:-}" == "1" ]]; then
        pipeline_abs_path "${PROJECT_ROOT}/test"
        return
    fi
    if [[ -n "${PIPELINE_DATA_DIR:-}" ]]; then
        pipeline_abs_path "$PIPELINE_DATA_DIR"
        return
    fi
    local production trial trial_mmcif
    production="$(pipeline_abs_path "${PROJECT_ROOT}/Data")"
    trial="$(pipeline_abs_path "${PROJECT_ROOT}/Data/test")"
    trial_mmcif="$(pipeline_abs_path "${PROJECT_ROOT}/Data/test/mmCIF")"
    if [[ -n "${MMCIF_SOURCE:-}" ]]; then
        local mmcif
        mmcif="$(pipeline_abs_path "$MMCIF_SOURCE")"
        if [[ "$mmcif" == "$trial_mmcif" ]]; then
            printf '%s\n' "$trial"
            return
        fi
    fi
    printf '%s\n' "$production"
}
