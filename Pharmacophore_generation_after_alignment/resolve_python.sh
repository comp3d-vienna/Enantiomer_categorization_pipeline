#!/bin/bash
# Resolve a single Python 3.9 interpreter for pharmacophore steps (source this file).
#
# Use one conda env (categorize_pipeline) with CDPL, pandas, and pytz.
# CDPL may also need PYTHONPATH to CDPKit when not set by conda activate.

CDPL_PYTHONPATH_DEFAULT="${CDPL_PYTHONPATH:-/data/shared/software/CDPKit-head-RH9/Python}"
PIPELINE_CONDA_ENV="${PIPELINE_CONDA_ENV:-categorize_pipeline}"

_resolve_conda_python() {
    local env_name="$1"
    if ! command -v conda >/dev/null 2>&1; then
        return 1
    fi
    local candidate
    candidate="$(conda run -n "$env_name" which python 2>/dev/null || true)"
    if [[ -n "$candidate" && -x "$candidate" ]]; then
        echo "$candidate"
        return 0
    fi
    return 1
}

_python_has_cdpl() {
    local py="$1"
    local pp="${2:-}"
    [[ -n "$py" && -x "$py" ]] || return 1
    if [[ -n "$pp" ]]; then
        PYTHONPATH="$pp${PYTHONPATH:+:$PYTHONPATH}" "$py" -c "import CDPL.Chem; import CDPL.Pharm" 2>/dev/null
    else
        "$py" -c "import CDPL.Chem; import CDPL.Pharm" 2>/dev/null
    fi
}

_python_has_pandas() {
    local py="$1"
    [[ -n "$py" && -x "$py" ]] && "$py" -c "import pandas, pytz" 2>/dev/null
}

_python_version_ok() {
    local py="$1"
    [[ -n "$py" && -x "$py" ]] || return 1
    "$py" -c 'import sys; raise SystemExit(0 if sys.version_info[:2] == (3, 9) else 1)' 2>/dev/null
}

_ensure_cdpl_pythonpath() {
    if _python_has_cdpl "${PIPELINE_PYTHON:-}"; then
        return 0
    fi
    if [[ -d "$CDPL_PYTHONPATH_DEFAULT" ]]; then
        export PYTHONPATH="${CDPL_PYTHONPATH_DEFAULT}${PYTHONPATH:+:$PYTHONPATH}"
        _python_has_cdpl "${PIPELINE_PYTHON:-}"
    fi
}

_resolve_python_candidates() {
    local candidate
    for candidate in \
        "${PYTHON_CMD:-}" \
        "$(command -v python 2>/dev/null || true)" \
        "$(command -v python3 2>/dev/null || true)" \
        "$(_resolve_conda_python "${CONDA_DEFAULT_ENV:-}" || true)" \
        "$(_resolve_conda_python "$PIPELINE_CONDA_ENV" || true)"; do
        if [[ -n "$candidate" && -x "$candidate" ]]; then
            echo "$candidate"
        fi
    done | awk '!seen[$0]++'
}

resolve_python() {
    PIPELINE_PYTHON=""
    local candidate

    if [[ -n "${PYTHON_CMD:-}" ]]; then
        PIPELINE_PYTHON="$PYTHON_CMD"
    elif command -v python >/dev/null 2>&1 && _python_version_ok "$(command -v python)"; then
        PIPELINE_PYTHON="$(command -v python)"
    else
        while IFS= read -r candidate; do
            if _python_version_ok "$candidate"; then
                PIPELINE_PYTHON="$candidate"
                break
            fi
        done < <(_resolve_python_candidates)
    fi

    if [[ -z "$PIPELINE_PYTHON" ]]; then
        echo "Error: Python 3.9 not found for pharmacophore steps." >&2
        echo "  conda activate ${PIPELINE_CONDA_ENV}" >&2
        echo "  Or: export PYTHON_CMD=/path/to/python3.9" >&2
        return 1
    fi

    if ! _python_version_ok "$PIPELINE_PYTHON"; then
        echo "Error: pharmacophore requires Python 3.9 (found: $("$PIPELINE_PYTHON" --version 2>&1))." >&2
        return 1
    fi

    if ! _python_has_pandas "$PIPELINE_PYTHON"; then
        echo "Error: pandas/pytz not importable with: $PIPELINE_PYTHON" >&2
        echo "  conda activate ${PIPELINE_CONDA_ENV}" >&2
        echo "  python -m pip install -r requirements.txt" >&2
        return 1
    fi

    export PIPELINE_PYTHON
    _ensure_cdpl_pythonpath || {
        echo "Error: CDPL.Chem/Pharm not importable with: $PIPELINE_PYTHON" >&2
        echo "  CDPL requires PYTHONPATH to CDPKit, e.g.:" >&2
        echo "    export PYTHONPATH=${CDPL_PYTHONPATH_DEFAULT}" >&2
        echo "  Or run: conda activate ${PIPELINE_CONDA_ENV}" >&2
        return 1
    }
}
