#!/usr/bin/env python3
"""
Build a per-group feature table of maximum-relatedness RMSD and SILIRID similarity.

Reads ``feature_table.csv`` from ``whole_process.py`` and adds SILIRID similarity
of the maximum-relatedness m0–m1 pair plus the two 160-D count fingerprints.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

_CATEGORIZATION_DIR = Path(__file__).resolve().parent
if str(_CATEGORIZATION_DIR) not in sys.path:
    sys.path.insert(0, str(_CATEGORIZATION_DIR))

import paths
from Interaction_type_score_refine import (
    STRUCTURE_COUNT_RESIDUE_DEFAULT_CAP,
    build_file_signatures,
    format_silirid_fingerprint,
    min_max_similarity,
    silirid_slots,
    structure_silirid_fingerprint,
)

FEATURE_COLUMNS: tuple[str, ...] = ("maximum_relatedness_rmsd_A",)

FEATURE_TABLE_COLUMNS: tuple[str, ...] = (
    "PSG_identifier",
    "n_structures",
    "n_m0_structures",
    "n_m1_structures",
    "atom_count_unique_values",
    "atom_count_by_file",
    "maximum_relatedness_m0_file",
    "maximum_relatedness_m1_file",
    "maximum_relatedness_rmsd_A",
    "maximum_relatedness_silirid_similarity",
    "maximum_relatedness_m0_silirid",
    "maximum_relatedness_m1_silirid",
)

_PASSTHROUGH_COLUMNS: tuple[str, ...] = (
    "n_structures",
    "n_m0_structures",
    "n_m1_structures",
    "atom_count_unique_values",
    "atom_count_by_file",
)

_DROP_COLUMNS: tuple[str, ...] = (
    "category",
    "pipeline_category",
    "manual_label",
    "manual_label_join_method",
    "category_reason",
)


def psg_identifier_from_info_file(info_file: str) -> str:
    """``g000001_inchi_b55a025257_inconsistent.csv`` -> ``g000001_b55a025257``."""
    stem = Path(str(info_file).strip()).stem
    for suffix in ("_inconsistent", "_consistent"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    if "_inchi_" in stem:
        gid, hash_part = stem.split("_inchi_", 1)
        hash_part = hash_part.split("_", 1)[0]
        if gid and hash_part:
            return f"{gid}_{hash_part}"
    return stem


def _info_csv_stem_candidates(identifier: str) -> list[str]:
    """Filename stems to try when resolving a ``PSG_identifier`` to an interaction CSV."""
    raw = str(identifier).strip()
    if not raw:
        return []
    stem = Path(raw).stem if raw.endswith(".csv") else raw
    for suffix in ("_inconsistent", "_consistent"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break

    candidates: list[str] = []
    def _add(value: str) -> None:
        if value and value not in candidates:
            candidates.append(value)

    _add(stem)
    if "_inchi_" in stem:
        gid, hash_part = stem.split("_inchi_", 1)
        hash_part = hash_part.split("_", 1)[0]
        _add(f"{gid}_inchi_{hash_part}")
        _add(f"{gid}_{hash_part}")
    else:
        gid, sep, hash_part = stem.partition("_")
        if sep and gid and hash_part:
            _add(f"{gid}_inchi_{hash_part}")
    return candidates


def resolve_info_csv_path(info_dir: Path, identifier: str) -> Path | None:
    """Find the per-group interaction CSV from ``PSG_identifier`` or a leftover filename."""
    raw = str(identifier).strip()
    if not raw:
        return None
    if raw.endswith(".csv"):
        direct = info_dir / raw
        if direct.is_file():
            return direct
        canonical = _canonical_info_file(raw)
        candidate = info_dir / canonical
        if candidate.is_file():
            return candidate

    for stem in _info_csv_stem_candidates(raw):
        for suffix in (".csv", "_inconsistent.csv", "_consistent.csv"):
            candidate = info_dir / f"{stem}{suffix}"
            if candidate.is_file():
                return candidate
    return None


def normalize_group_key(group_name: str) -> list[str]:
    cleaned = group_name.strip().lstrip("{").rstrip("}").replace("'", "").strip()
    return cleaned.split("_")


def _is_invalid_category(category: object) -> bool:
    text = str(category).strip().lower()
    return text in {"invalid"}


def normalize_manual_label(value: object) -> str | None:
    """
    Canonical string form for manual labels.

    Pandas often reads label ``2`` from CSV as float ``2.0``; ``str(2.0)`` is
    ``\"2.0\"``, which must match label ``\"2\"``.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "<na>"}:
        return None
    # Float coercion: "2.0" / 2.0 -> "2"; leave "1.1", "1.2", "1.3" unchanged.
    if text.endswith(".0") and text[:-2].isdigit():
        return text[:-2]
    return text


def _group_to_info_file(group: str) -> str:
    parts = normalize_group_key(group)
    return "_".join([parts[0], "inchi"] + parts[1:]) + ".csv"


def _normalize_inchi_hash(hash_part: str) -> str:
    """Normalize inchi hash segment (e.g. float ``7813455300.0`` -> ``7813455300``)."""
    text = str(hash_part).strip()
    if text.endswith(".0") and text[:-2].replace(".", "", 1).isdigit():
        return text[:-2]
    return text


def group_key_from_info_file(info_file: str) -> str | None:
    """
    ``g000448_inchi_d463a70235.csv`` -> ``g000448_d463a70235``.

    Also accepts leftover ``_consistent`` / ``_inconsistent`` suffixes.
    """
    stem = str(info_file).strip()
    for suffix in ("_inconsistent.csv", "_consistent.csv"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    else:
        if stem.endswith(".csv"):
            stem = stem[: -len(".csv")]

    if "_inchi_" not in stem:
        return None
    gid, hash_part = stem.split("_inchi_", 1)
    if not gid or not hash_part:
        return None
    return f"{gid}_{_normalize_inchi_hash(hash_part)}"


def _canonical_info_file(info_file: str) -> str:
    """Map leftover flavored names to ``{gid}_inchi_{hash}.csv``."""
    stem = str(info_file).strip()
    for suffix in ("_inconsistent.csv", "_consistent.csv"):
        if stem.endswith(suffix):
            return stem[: -len(suffix)] + ".csv"
    return stem


def _closest_pair_silirid(
    interaction_df: pd.DataFrame,
    file_a: str,
    file_b: str,
    *,
    count_cap: int | None = STRUCTURE_COUNT_RESIDUE_DEFAULT_CAP,
) -> tuple[float, str, str]:
    """SILIRID similarity and fingerprints for one structure pair (no dedup filter)."""
    empty = (float("nan"), "", "")
    if interaction_df.empty or not file_a or not file_b:
        return empty
    try:
        sig_a = build_file_signatures(interaction_df, str(file_a).strip())
        sig_b = build_file_signatures(interaction_df, str(file_b).strip())
    except (KeyError, ValueError):
        return empty
    if not sig_a or not sig_b:
        return empty
    fp_a = structure_silirid_fingerprint(sig_a, count_cap=count_cap)
    fp_b = structure_silirid_fingerprint(sig_b, count_cap=count_cap)
    similarity = float(min_max_similarity(fp_a, fp_b))
    return similarity, format_silirid_fingerprint(fp_a), format_silirid_fingerprint(fp_b)


def manual_label_from_curation_filename(csv_path: Path) -> str:
    """
    Derive manual category from a curation CSV filename.

    Examples: ``manualcheck_11.csv`` -> ``1.1``, ``manualcheck_12.csv`` -> ``1.2``,
    ``manualcheck_13.csv`` -> ``1.3``, ``manualcheck_2.csv`` -> ``2``,
    ``manualcheck_invalid.csv`` -> ``invalid``.
    """
    stem = csv_path.stem
    prefix = "manualcheck_"
    code = stem[len(prefix) :] if stem.startswith(prefix) else stem
    if code == "invalid":
        return "invalid"
    if code == "2":
        return "2"
    if not code or not code.isdigit():
        raise ValueError(f"Cannot parse manual label from filename: {csv_path.name}")
    return ".".join(code)


def _parse_group_cell(value: object) -> str | None:
    if pd.isna(value):
        return None
    parts = normalize_group_key(str(value))
    if not parts:
        return None
    return "_".join(parts)


def load_manual_labels_from_curation_dir(
    manual_dir: Path,
    *,
    pattern: str = "*.csv",
    exclude_invalid: bool = True,
) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    """
    Build manual-label lookups from ``Data/Manual_curation/*.csv``.

    Returns ``(by_info_file, by_group_key, excluded_by_group_key)``.
    Each curated group is registered under ``{gid}_inchi_{hash}.csv``.

    When ``exclude_invalid`` is True (default), ``manualcheck_invalid.csv``
    (label ``invalid``) is omitted from the label map but recorded in
    ``excluded_by_group_key`` for diagnosis.
    """
    manual_dir = Path(manual_dir)
    if not manual_dir.is_dir():
        raise FileNotFoundError(f"Manual curation directory not found: {manual_dir}")

    csv_files = sorted(manual_dir.glob(pattern))
    if not csv_files:
        raise FileNotFoundError(f"No {pattern} files under {manual_dir}")

    by_info_file: dict[str, str] = {}
    by_group_key: dict[str, str] = {}
    excluded_by_group_key: dict[str, str] = {}
    for csv_path in csv_files:
        label = manual_label_from_curation_filename(csv_path)
        df = pd.read_csv(csv_path, sep="\t")
        if "group" not in df.columns:
            raise ValueError(f"Manual curation CSV missing 'group' column: {csv_path}")

        group_keys: set[str] = set()
        for value in df["group"].dropna().unique():
            group_key = _parse_group_cell(value)
            if group_key:
                group_keys.add(group_key)

        for group_key in group_keys:
            if exclude_invalid and _is_invalid_category(label):
                excluded_by_group_key[group_key] = label
                continue
            by_group_key[group_key] = label
            by_info_file[_group_to_info_file(group_key)] = label

    return by_info_file, by_group_key, excluded_by_group_key


def resolve_manual_label(
    info_file: str,
    *,
    by_info_file: dict[str, str] | None,
    by_group_key: dict[str, str] | None = None,
    excluded_by_group_key: dict[str, str] | None = None,
) -> tuple[str | None, str]:
    """
    Resolve manual label for a pipeline ``info_file``.

    Returns ``(label, join_method)`` where ``join_method`` explains how the label
    was found, or why it is missing.
    """
    if not by_info_file:
        return None, "manual_curation_disabled"

    info_file = str(info_file).strip()
    if info_file in by_info_file:
        return by_info_file[info_file], "exact_info_file"

    canonical = _canonical_info_file(info_file)
    if canonical != info_file and canonical in by_info_file:
        return by_info_file[canonical], f"canonical_info_file:{canonical}"

    group_key = group_key_from_info_file(info_file)
    if group_key and by_group_key and group_key in by_group_key:
        return by_group_key[group_key], f"group_key:{group_key}"

    if group_key and excluded_by_group_key and group_key in excluded_by_group_key:
        label = excluded_by_group_key[group_key]
        return None, (
            f"manual_curation_invalid_only "
            f"(manual_label={label}; group_key={group_key}; excluded as invalid)"
        )

    if group_key and by_group_key:
        gid = group_key.split("_", 1)[0]
        file_hash = _normalize_inchi_hash(group_key.split("_", 1)[1])
        for key, label in by_group_key.items():
            if key.split("_", 1)[0] != gid:
                continue
            if _normalize_inchi_hash(key.split("_", 1)[1]) == file_hash:
                return label, f"group_key_normalized:{key}"

    if group_key and by_group_key:
        return None, f"not_in_manual_curation (group_key={group_key})"

    return None, "unparseable_info_file"


def diagnose_manual_label_join(
    info_file: str,
    *,
    by_info_file: dict[str, str] | None,
    by_group_key: dict[str, str] | None = None,
    excluded_by_group_key: dict[str, str] | None = None,
) -> dict[str, object]:
    """Step-by-step diagnosis for one pipeline row."""
    label, join_method = resolve_manual_label(
        info_file,
        by_info_file=by_info_file,
        by_group_key=by_group_key,
        excluded_by_group_key=excluded_by_group_key,
    )
    group_key = group_key_from_info_file(info_file)
    canonical = _canonical_info_file(info_file)
    excluded_label = (
        excluded_by_group_key.get(group_key) if group_key and excluded_by_group_key else None
    )
    return {
        "info_file": info_file,
        "group_key": group_key,
        "canonical_info_file": canonical,
        "exact_info_file_hit": bool(by_info_file and info_file in by_info_file),
        "canonical_info_file_hit": bool(
            by_info_file and canonical != info_file and canonical in by_info_file
        ),
        "group_key_hit": bool(by_group_key and group_key and group_key in by_group_key),
        "excluded_manual_label": excluded_label,
        "manual_label": label,
        "join_method": join_method,
        "resolved": label is not None,
    }


def build_feature_table(
    groups_df: pd.DataFrame,
    info_dir: Path,
) -> pd.DataFrame:
    """
    One row per group with maximum-relatedness RMSD and SILIRID similarity.
    """
    df = groups_df.copy()
    df = df.drop(columns=list(_DROP_COLUMNS), errors="ignore")

    if df.empty:
        return pd.DataFrame(columns=list(FEATURE_TABLE_COLUMNS))

    if "PSG_identifier" not in df.columns and "info_file" in df.columns:
        df["PSG_identifier"] = df["info_file"].map(psg_identifier_from_info_file)

    required = (
        "PSG_identifier",
        "maximum_relatedness_m0_file",
        "maximum_relatedness_m1_file",
        "maximum_relatedness_rmsd_A",
    )
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(
            f"feature table missing required columns: {missing}. "
            "Run whole_process.py first."
        )

    interaction_cache: dict[str, pd.DataFrame] = {}
    rows: list[dict[str, object]] = []

    for _, row in df.sort_values("PSG_identifier").iterrows():
        psg = str(row["PSG_identifier"])
        lookup = (
            str(row["info_file"])
            if "info_file" in row.index and pd.notna(row.get("info_file"))
            else psg
        )
        file_m0 = str(row.get("maximum_relatedness_m0_file", "") or "").strip()
        file_m1 = str(row.get("maximum_relatedness_m1_file", "") or "").strip()

        similarity = float("nan")
        silirid_m0 = ""
        silirid_m1 = ""
        if file_m0 and file_m1:
            csv_path = resolve_info_csv_path(info_dir, lookup)
            if csv_path is not None:
                cache_key = csv_path.name
                if cache_key not in interaction_cache:
                    interaction_cache[cache_key] = pd.read_csv(csv_path)
                similarity, silirid_m0, silirid_m1 = _closest_pair_silirid(
                    interaction_cache[cache_key],
                    file_m0,
                    file_m1,
                )

        out_row: dict[str, object] = {
            "PSG_identifier": psg,
            "maximum_relatedness_m0_file": file_m0,
            "maximum_relatedness_m1_file": file_m1,
            "maximum_relatedness_rmsd_A": row.get("maximum_relatedness_rmsd_A"),
            "maximum_relatedness_silirid_similarity": similarity,
            "maximum_relatedness_m0_silirid": silirid_m0,
            "maximum_relatedness_m1_silirid": silirid_m1,
        }
        for col in _PASSTHROUGH_COLUMNS:
            if col in row.index and col not in out_row:
                out_row[col] = row.get(col)
        rows.append(out_row)

    table = pd.DataFrame(rows)
    if table.empty:
        return pd.DataFrame(columns=list(FEATURE_TABLE_COLUMNS))
    ordered = [c for c in FEATURE_TABLE_COLUMNS if c in table.columns]
    extra = [c for c in table.columns if c not in ordered]
    table = table[ordered + extra]
    if "PSG_identifier" in table.columns:
        table = table.sort_values("PSG_identifier")
    return table


def _write_silirid_slot_legend(out_dir: Path) -> Path:
    """Write amino-acid × feature-type order of the 160-D SILIRID count vector."""
    legend = pd.DataFrame(
        [
            {
                "slot_index": index,
                "amino_acid": aa,
                "feature_type": feature_type,
            }
            for index, (aa, feature_type) in enumerate(silirid_slots())
        ]
    )
    legend_path = out_dir / "silirid_fingerprint_slots.csv"
    legend.to_csv(legend_path, index=False)
    return legend_path


def _feature_ready_mask(table: pd.DataFrame) -> pd.Series:
    cols = [c for c in FEATURE_COLUMNS if c in table.columns]
    if not cols:
        return pd.Series(False, index=table.index)
    return table[cols].notna().all(axis=1)


def _incomplete_reason(row: pd.Series) -> str:
    def as_int(value: object) -> int:
        if value is None or pd.isna(value):
            return 0
        return int(value)

    n_m0 = as_int(row.get("n_m0_structures", 0))
    n_m1 = as_int(row.get("n_m1_structures", 0))
    if n_m0 < 1 or n_m1 < 1:
        return "missing_m0_or_m1"
    unique_vals = str(row.get("atom_count_unique_values", "") or "").strip()
    if unique_vals and "," in unique_vals:
        return "inconsistent_atom_count"
    by_file = str(row.get("atom_count_by_file", "") or "")
    if "=missing" in by_file:
        return "missing_canonical_sdf"
    return "no_computable_m0_m1_pair"


def _print_feature_completeness_summary(table: pd.DataFrame, out_dir: Path) -> None:
    """Report groups with computable closest-pair features vs incomplete rows."""
    ready = _feature_ready_mask(table)
    n_ready = int(ready.sum())
    n_total = len(table)
    n_skipped = n_total - n_ready
    print(f"  complete closest-pair features: {n_ready}/{n_total}")
    if n_skipped == 0:
        return
    print(f"  incomplete (m0–m1 metrics not computable): {n_skipped}")
    skipped = table.loc[~ready].copy()
    reasons = skipped.apply(_incomplete_reason, axis=1).value_counts()
    for reason, count in reasons.items():
        print(f"    {reason}: {int(count)}")
    incomplete_csv = out_dir / "feature_table_incomplete.csv"
    skipped.to_csv(incomplete_csv, index=False)
    print(f"  incomplete-feature rows: {incomplete_csv}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Add SILIRID similarity and fingerprints to the closest-pair feature table."
        )
    )
    parser.add_argument(
        "--in-csv",
        "--categories-csv",
        dest="in_csv",
        type=Path,
        default=Path(paths.FEATURE_TABLE_CSV),
        help="Input table from whole_process.py (default: feature_table.csv).",
    )
    parser.add_argument(
        "--info-dir",
        type=Path,
        default=Path(paths.INTERACTION_TYPE_INFO_DIR),
        help="Per-group interaction-type CSV directory.",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=None,
        help="Output CSV (default: overwrite --in-csv).",
    )

    args = parser.parse_args()

    in_csv = args.in_csv.resolve()
    if not in_csv.is_file():
        raise FileNotFoundError(f"feature table not found: {in_csv}")

    info_dir = args.info_dir.resolve()
    if not info_dir.is_dir():
        raise FileNotFoundError(f"info_dir not found: {info_dir}")

    out_csv = args.out_csv
    if out_csv is None:
        out_csv = in_csv
    out_csv = out_csv.resolve()

    groups = pd.read_csv(in_csv)

    table = build_feature_table(groups, info_dir)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_csv, index=False)
    legend_path = _write_silirid_slot_legend(out_csv.parent)
    leftover = out_csv.parent / "feature_table_missing_manual_label.csv"
    if leftover.is_file():
        leftover.unlink()

    print(f"Wrote {len(table)} row(s) to {out_csv}")
    print(
        "  SILIRID fingerprints: "
        "maximum_relatedness_m0_silirid / maximum_relatedness_m1_silirid"
    )
    print(f"  SILIRID slot order: {legend_path}")
    print(f"  RMSD feature: {', '.join(FEATURE_COLUMNS)}")
    _print_feature_completeness_summary(table, out_csv.parent)


if __name__ == "__main__":
    main()
