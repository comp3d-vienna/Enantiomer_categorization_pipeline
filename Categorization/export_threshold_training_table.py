#!/usr/bin/env python3
"""
Build a per-group training table of closest-pair RMSD and SILIRID similarity.

Exports closest-pair RMSD (from ``complex_categories.csv``), SILIRID
similarity of that same pair, and the two 160-D SILIRID count fingerprints
used for the similarity. Invalid groups are excluded. Custom RMSD scoring
uses RMSD only; similarity and fingerprints are kept in the table for reporting.
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
    SILIRID_SIMILARITY_METRIC,
    STRUCTURE_COUNT_RESIDUE_DEFAULT_CAP,
    STRUCTURE_SILIRID_LEN,
    build_file_signatures,
    format_silirid_fingerprint,
    min_max_similarity,
    silirid_slots,
    structure_silirid_fingerprint,
)

TREE_FEATURES: tuple[str, ...] = ("closest_cross_tag_rmsd_A",)


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
    ``\"2.0\"``, which must match training label ``\"2\"``.
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
    (label ``invalid``) is omitted from the training map but recorded in
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
            f"(manual_label={label}; group_key={group_key}; not used for 1.x/2 training)"
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
    pipeline_category: str | None = None,
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
        "pipeline_category": pipeline_category,
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


def build_threshold_training_table(
    categories_df: pd.DataFrame,
    info_dir: Path,
    *,
    manual_by_info_file: dict[str, str] | None = None,
    manual_by_group_key: dict[str, str] | None = None,
    excluded_manual_by_group_key: dict[str, str] | None = None,
) -> pd.DataFrame:
    """
    One row per training group with closest-pair RMSD and SILIRID similarity.
    """
    df = categories_df.copy()
    df["category"] = df["category"].astype(str).str.strip()
    df = df.loc[~df["category"].map(_is_invalid_category)].copy()

    required = (
        "info_file",
        "closest_cross_tag_file_a",
        "closest_cross_tag_file_b",
        "closest_cross_tag_rmsd_A",
    )
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(
            f"categories table missing required columns: {missing}. "
            "Run whole_process.py first."
        )

    interaction_cache: dict[str, pd.DataFrame] = {}
    rows: list[dict[str, object]] = []

    for _, row in df.sort_values("info_file").iterrows():
        info_file = str(row["info_file"])
        file_a = str(row.get("closest_cross_tag_file_a", "") or "").strip()
        file_b = str(row.get("closest_cross_tag_file_b", "") or "").strip()

        similarity = float("nan")
        silirid_a = ""
        silirid_b = ""
        if file_a and file_b:
            csv_path = info_dir / info_file
            if csv_path.is_file():
                if info_file not in interaction_cache:
                    interaction_cache[info_file] = pd.read_csv(csv_path)
                similarity, silirid_a, silirid_b = _closest_pair_silirid(
                    interaction_cache[info_file],
                    file_a,
                    file_b,
                )

        manual_label = None
        manual_join_method = ""
        row_manual = row.get("manual_label")
        if pd.notna(row_manual) and str(row_manual).strip():
            manual_label = normalize_manual_label(row_manual)
            if manual_label is not None and _is_invalid_category(manual_label):
                manual_label = None
                manual_join_method = "manual_curation_invalid"
            else:
                manual_join_method = "from_complex_categories"
        elif manual_by_info_file is not None:
            manual_label, manual_join_method = resolve_manual_label(
                info_file,
                by_info_file=manual_by_info_file,
                by_group_key=manual_by_group_key,
                excluded_by_group_key=excluded_manual_by_group_key,
            )
            manual_label = normalize_manual_label(manual_label)

        rows.append(
            {
                "info_file": info_file,
                "pipeline_category": row["category"],
                "manual_label": manual_label,
                "manual_label_join_method": manual_join_method,
                "closest_cross_tag_file_a": file_a,
                "closest_cross_tag_file_b": file_b,
                "closest_cross_tag_rmsd_A": row.get("closest_cross_tag_rmsd_A"),
                "closest_cross_tag_similarity": similarity,
                "closest_cross_tag_silirid_a": silirid_a,
                "closest_cross_tag_silirid_b": silirid_b,
                "category_reason": row.get("category_reason"),
                "n_m0_structures": row.get("n_m0_structures"),
                "n_m1_structures": row.get("n_m1_structures"),
                "similarity_metric": SILIRID_SIMILARITY_METRIC,
                "silirid_vector_length": STRUCTURE_SILIRID_LEN,
                "silirid_count_cap": STRUCTURE_COUNT_RESIDUE_DEFAULT_CAP,
            }
        )

    return pd.DataFrame(rows)


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


def _training_feature_mask(table: pd.DataFrame) -> pd.Series:
    cols = [c for c in TREE_FEATURES if c in table.columns]
    if not cols:
        return pd.Series(False, index=table.index)
    return table[cols].notna().all(axis=1)


def _print_training_readiness_summary(table: pd.DataFrame, out_dir: Path) -> None:
    """Report groups usable for custom RMSD scoring vs skipped with feature failures."""
    ready = _training_feature_mask(table)
    n_ready = int(ready.sum())
    n_total = len(table)
    n_skipped = n_total - n_ready
    print(f"  ready for custom RMSD scoring: {n_ready}/{n_total}")
    if n_skipped == 0:
        return
    print(
        "  skipped (manual labels kept; cross-tag metrics not computable): "
        f"{n_skipped}"
    )
    skipped = table.loc[~ready].copy()
    if "category_reason" in skipped.columns:
        failure = (
            skipped["category_reason"]
            .astype(str)
            .str.split(";")
            .str[-1]
            .value_counts()
            .sort_values(ascending=False)
        )
        for reason, count in failure.items():
            print(f"    {reason}: {int(count)}")
    incomplete_csv = out_dir / "threshold_training_table_incomplete_features.csv"
    skipped.to_csv(incomplete_csv, index=False)
    print(f"  incomplete-feature rows: {incomplete_csv}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Export closest-pair RMSD, SILIRID similarity, and SILIRID fingerprints "
            "(custom RMSD scoring uses RMSD only)."
        )
    )
    parser.add_argument(
        "--categories-csv",
        type=Path,
        default=Path(paths.CATEGORIES_CSV),
        help="Input complex_categories.csv from whole_process.py.",
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
        help="Output CSV (default: <categories-dir>/threshold_training_table.csv).",
    )
    parser.add_argument(
        "--manual-curation-dir",
        type=Path,
        default=Path(paths.MANUAL_CURATION_DIR),
        help="Directory of manual curation CSVs (default: Data/Manual_curation).",
    )

    args = parser.parse_args()

    categories_csv = args.categories_csv.resolve()
    if not categories_csv.is_file():
        raise FileNotFoundError(f"categories CSV not found: {categories_csv}")

    info_dir = args.info_dir.resolve()
    if not info_dir.is_dir():
        raise FileNotFoundError(f"info_dir not found: {info_dir}")

    out_csv = args.out_csv
    if out_csv is None:
        out_csv = categories_csv.parent / "threshold_training_table.csv"
    out_csv = out_csv.resolve()

    categories = pd.read_csv(categories_csv)

    manual_by_info_file: dict[str, str] | None = None
    manual_by_group_key: dict[str, str] | None = None
    excluded_manual_by_group_key: dict[str, str] | None = None
    curation_dir = args.manual_curation_dir.resolve()
    if curation_dir.is_dir():
        manual_by_info_file, manual_by_group_key, excluded_manual_by_group_key = (
            load_manual_labels_from_curation_dir(curation_dir)
        )
    else:
        print(f"Note: manual curation dir not found, skipping labels: {curation_dir}")

    table = build_threshold_training_table(
        categories,
        info_dir,
        manual_by_info_file=manual_by_info_file,
        manual_by_group_key=manual_by_group_key,
        excluded_manual_by_group_key=excluded_manual_by_group_key,
    )

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    # Keep labels as strings so category "2" is not written/read as float 2.0.
    if "manual_label" in table.columns:
        table = table.copy()
        table["manual_label"] = table["manual_label"].map(normalize_manual_label)
        table["manual_label"] = table["manual_label"].astype("string")
    table.to_csv(out_csv, index=False)
    legend_path = _write_silirid_slot_legend(out_csv.parent)

    n_manual = int(table["manual_label"].notna().sum()) if "manual_label" in table.columns else 0
    print(f"Wrote {len(table)} row(s) to {out_csv}")
    print(f"  SILIRID fingerprints: closest_cross_tag_silirid_a / closest_cross_tag_silirid_b")
    print(f"  SILIRID slot order: {legend_path}")
    print(f"  manual labels joined: {n_manual}/{len(table)}")
    print(f"  RMSD feature: {', '.join(TREE_FEATURES)}")
    _print_training_readiness_summary(table, out_csv.parent)

    if manual_by_info_file is not None and "manual_label" in table.columns:
        unlabeled = table.loc[table["manual_label"].isna()].copy()
        if not unlabeled.empty:
            diagnosis_rows = [
                diagnose_manual_label_join(
                    str(row["info_file"]),
                    by_info_file=manual_by_info_file,
                    by_group_key=manual_by_group_key,
                    excluded_by_group_key=excluded_manual_by_group_key,
                    pipeline_category=str(row.get("pipeline_category", "")),
                )
                for _, row in unlabeled.iterrows()
            ]
            diagnosis_df = pd.DataFrame(diagnosis_rows)
            diagnosis_csv = out_csv.parent / "threshold_training_table_missing_manual_label.csv"
            diagnosis_df.to_csv(diagnosis_csv, index=False)
            print(f"  missing manual_label: {len(unlabeled)} row(s) -> {diagnosis_csv}")
        else:
            print("  missing manual_label: 0")


if __name__ == "__main__":
    main()
