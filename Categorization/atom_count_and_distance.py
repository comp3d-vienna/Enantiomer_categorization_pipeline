from __future__ import annotations

import os
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem


def interaction_tsv_to_canon_sdf_name(tsv_file: str) -> str:
    """Map ``interaction_..._m0.tsv`` to the canonical ligand SDF name."""
    return "_".join(str(tsv_file).split("_")[1:]).replace(".tsv", "_canon.sdf")


def canon_sdf_path_from_interaction_tsv(tsv_file: str, ligand_structures_dir: str) -> str:
    return os.path.join(ligand_structures_dir, interaction_tsv_to_canon_sdf_name(tsv_file))


def enantiomer_tag_from_interaction_file(file_name: str) -> str | None:
    name = str(file_name)
    if name.endswith("_m0.tsv"):
        return "m0"
    if name.endswith("_m1.tsv"):
        return "m1"
    return None


def _load_atom_count(sdf_path: str) -> int | None:
    if not os.path.exists(sdf_path):
        return None
    try:
        mol = Chem.SDMolSupplier(sdf_path, removeHs=True)[0]
        if mol is None:
            return None
        return int(mol.GetNumAtoms())
    except Exception:
        return None


def enantiomer_tag_counts(files: list[str]) -> tuple[int, int]:
    n_m0 = sum(1 for f in files if enantiomer_tag_from_interaction_file(f) == "m0")
    n_m1 = sum(1 for f in files if enantiomer_tag_from_interaction_file(f) == "m1")
    return n_m0, n_m1


def has_both_enantiomer_tags(n_m0: int, n_m1: int) -> bool:
    return n_m0 >= 1 and n_m1 >= 1


def has_cross_tag_distance_pairs(dist_df: pd.DataFrame) -> bool:
    if dist_df.empty or "pair_type" not in dist_df.columns:
        return False
    return not dist_df.loc[dist_df["pair_type"] == "cross_tag"].empty


def _inconsistent_structure_files(present: list[tuple[str, int]]) -> list[str]:
    if len(present) <= 1:
        return []
    counts = [count for _, count in present]
    if len(set(counts)) == 1:
        return []
    counter = Counter(counts)
    max_freq = max(counter.values())
    if max_freq > 1:
        mode_count = counter.most_common(1)[0][0]
        return [file_name for file_name, count in present if count != mode_count]
    return [file_name for file_name, _ in present]


def _build_atom_count_entry_details(
    files: list[str],
    atom_counts: list[int | None],
) -> dict[str, object]:
    by_file_parts: list[str] = []
    missing_files: list[str] = []
    present: list[tuple[str, int]] = []

    for file_name, count in zip(files, atom_counts):
        if count is None:
            by_file_parts.append(f"{file_name}=missing")
            missing_files.append(file_name)
        else:
            by_file_parts.append(f"{file_name}={count}")
            present.append((file_name, count))

    return {
        "atom_count_by_file": "; ".join(by_file_parts),
        "atom_count_missing_files": "; ".join(missing_files),
        "atom_count_inconsistent_files": "; ".join(_inconsistent_structure_files(present)),
    }


def categorize_complex_by_atom_counts(
    info_csv_path: str | Path,
    ligand_structures_dir: str | Path,
    *,
    file_column: str = "file",
) -> dict[str, object]:
    info_csv_path = Path(info_csv_path)
    df = pd.read_csv(info_csv_path)
    files = sorted({str(x) for x in df[file_column].dropna().unique()})
    sdf_paths = [
        canon_sdf_path_from_interaction_tsv(f, str(ligand_structures_dir)) for f in files
    ]
    atom_counts = [_load_atom_count(p) for p in sdf_paths]

    n_files = len(files)
    n_missing_or_unreadable = sum(1 for c in atom_counts if c is None)
    present_counts = sorted({c for c in atom_counts if c is not None})

    is_consistent = (n_files > 0) and (n_missing_or_unreadable == 0) and (len(present_counts) == 1)
    category = 1 if is_consistent else 2
    n_m0, n_m1 = enantiomer_tag_counts(files)
    details = _build_atom_count_entry_details(files, atom_counts)

    return {
        "info_file": info_csv_path.name,
        "category": category,
        "n_structures": n_files,
        "n_m0_structures": n_m0,
        "n_m1_structures": n_m1,
        "n_missing_or_unreadable_sdf": int(n_missing_or_unreadable),
        "atom_count_unique_values": ",".join(str(c) for c in present_counts),
        "atom_count_by_file": details["atom_count_by_file"],
        "atom_count_missing_files": details["atom_count_missing_files"],
        "atom_count_inconsistent_files": details["atom_count_inconsistent_files"],
    }


def _per_atom_distances(mol1_sdf: str, mol2_sdf: str, sdf_a: str, sdf_b: str) -> pd.DataFrame:
    mol1 = Chem.SDMolSupplier(mol1_sdf, removeHs=True)[0]
    mol2 = Chem.SDMolSupplier(mol2_sdf, removeHs=True)[0]
    if mol1 is None or mol2 is None:
        raise ValueError(f"One or both SD files failed to load: {mol1_sdf} {mol2_sdf}")
    conf1 = mol1.GetConformer()
    conf2 = mol2.GetConformer()
    n_atoms = mol1.GetNumAtoms()
    n_atoms_mol2 = mol2.GetNumAtoms()
    if n_atoms != n_atoms_mol2:
        raise ValueError(
            f"Atom count mismatch: {sdf_a} ({n_atoms}) vs {sdf_b} ({n_atoms_mol2})"
        )
    rows = []
    for i in range(n_atoms):
        p1 = np.array(conf1.GetAtomPosition(i))
        p2 = np.array(conf2.GetAtomPosition(i))
        rows.append((i, mol1.GetAtomWithIdx(i).GetSymbol(), float(np.linalg.norm(p1 - p2))))
    return pd.DataFrame(rows, columns=["atom_index", "atom_type", "distance_A"])


def cross_tag_feature_failure_reason(
    cat: dict[str, object],
    errors: list[dict[str, object]],
) -> str:
    n_m0 = int(cat.get("n_m0_structures", 0) or 0)
    n_m1 = int(cat.get("n_m1_structures", 0) or 0)
    if not has_both_enantiomer_tags(n_m0, n_m1):
        return "missing_m0_or_m1"

    unique_vals = str(cat.get("atom_count_unique_values", "") or "").strip()
    if unique_vals and "," in unique_vals:
        return "inconsistent_atom_count"

    if int(cat.get("n_missing_or_unreadable_sdf", 0) or 0) > 0:
        return "missing_canonical_sdf"

    cross_errors = [e for e in errors if e.get("pair_type") == "cross_tag"]
    if cross_errors:
        messages = [str(e.get("error", "")) for e in cross_errors]
        if all("Atom count mismatch" in msg for msg in messages):
            return "cross_tag_atom_count_mismatch"
        if any("failed to load" in msg.lower() for msg in messages):
            return "sdf_load_failed"
        if all(msg == "missing_canonical_sdf" for msg in messages):
            return "missing_canonical_sdf"

    if n_m0 >= 1 and n_m1 >= 1:
        return "missing_canonical_sdf"
    return "no_computable_cross_tag_pair"


def compute_atom_distances_for_complex(
    info_csv_path: str | Path,
    ligand_structures_dir: str | Path,
    *,
    file_column: str = "file",
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    """Per-atom distances for every m0 vs m1 pair (used to pick the closest pair)."""
    info_csv_path = Path(info_csv_path)
    df = pd.read_csv(info_csv_path)
    files = sorted({str(x) for x in df[file_column].dropna().unique()})
    m0_files = [f for f in files if enantiomer_tag_from_interaction_file(f) == "m0"]
    m1_files = [f for f in files if enantiomer_tag_from_interaction_file(f) == "m1"]

    tables: list[pd.DataFrame] = []
    errors: list[dict[str, object]] = []

    for file_a_tsv, file_b_tsv in ((a, b) for a in m0_files for b in m1_files):
        sdf_a = interaction_tsv_to_canon_sdf_name(file_a_tsv)
        sdf_b = interaction_tsv_to_canon_sdf_name(file_b_tsv)
        path_a = canon_sdf_path_from_interaction_tsv(file_a_tsv, str(ligand_structures_dir))
        path_b = canon_sdf_path_from_interaction_tsv(file_b_tsv, str(ligand_structures_dir))
        if not (os.path.exists(path_a) and os.path.exists(path_b)):
            errors.append(
                {
                    "info_file": info_csv_path.name,
                    "pair_type": "cross_tag",
                    "file_a": file_a_tsv,
                    "file_b": file_b_tsv,
                    "sdf_a": sdf_a,
                    "sdf_b": sdf_b,
                    "error": "missing_canonical_sdf",
                }
            )
            continue
        try:
            df_dist = _per_atom_distances(path_a, path_b, sdf_a, sdf_b)
            df_dist.insert(0, "info_file", info_csv_path.name)
            df_dist.insert(1, "pair_type", "cross_tag")
            df_dist.insert(2, "file_a", file_a_tsv)
            df_dist.insert(3, "file_b", file_b_tsv)
            tables.append(df_dist)
        except Exception as e:
            errors.append(
                {
                    "info_file": info_csv_path.name,
                    "pair_type": "cross_tag",
                    "file_a": file_a_tsv,
                    "file_b": file_b_tsv,
                    "sdf_a": sdf_a,
                    "sdf_b": sdf_b,
                    "error": str(e),
                }
            )

    if not tables:
        return pd.DataFrame(), errors
    return pd.concat(tables, ignore_index=True), errors


def closest_cross_tag_pair_per_complex(dist_df: pd.DataFrame) -> pd.DataFrame:
    """
    Per group: m0–m1 pair with the smallest mean per-atom distance.

    RMSD of that same pair is the exported geometry feature.
    """
    if dist_df.empty:
        return pd.DataFrame()
    required = {"info_file", "pair_type", "file_a", "file_b", "distance_A"}
    missing = required - set(dist_df.columns)
    if missing:
        raise ValueError(f"closest_cross_tag_pair_per_complex missing columns: {sorted(missing)}")

    cross = dist_df.loc[dist_df["pair_type"] == "cross_tag"]
    if cross.empty:
        return pd.DataFrame()

    rows: list[dict[str, object]] = []
    for info_file, sub in cross.groupby("info_file", sort=False):
        pair_stats = (
            sub.groupby(["file_a", "file_b"], dropna=False)["distance_A"]
            .agg(pair_mean="mean", pair_rmsd=lambda d: float(np.sqrt((d ** 2).mean())))
            .reset_index()
        )
        best = pair_stats.loc[pair_stats["pair_mean"].idxmin()]
        rows.append(
            {
                "info_file": info_file,
                "closest_cross_tag_file_a": best["file_a"],
                "closest_cross_tag_file_b": best["file_b"],
                "closest_cross_tag_rmsd_A": float(best["pair_rmsd"]),
            }
        )
    return pd.DataFrame(rows)


def _finalize_and_write_categories(
    category_rows: list[dict[str, object]],
    all_distance_tables: list[pd.DataFrame],
    all_errors: list[dict[str, object]],
    out_dir: Path,
) -> None:
    categories_df = pd.DataFrame(category_rows).sort_values(["category", "info_file"])
    closest_df = (
        closest_cross_tag_pair_per_complex(pd.concat(all_distance_tables, ignore_index=True))
        if all_distance_tables
        else pd.DataFrame()
    )

    if all_errors:
        pd.DataFrame(all_errors).to_csv(out_dir / "atom_distance_errors.csv", index=False)

    categories_updated = categories_df.copy()
    categories_updated["category"] = categories_updated["category"].astype(str)
    categories_updated["closest_cross_tag_file_a"] = ""
    categories_updated["closest_cross_tag_file_b"] = ""
    categories_updated["closest_cross_tag_rmsd_A"] = float("nan")

    if not closest_df.empty:
        categories_updated = categories_updated.merge(
            closest_df, on="info_file", how="left", suffixes=("", "_best")
        )
        for col in (
            "closest_cross_tag_file_a",
            "closest_cross_tag_file_b",
            "closest_cross_tag_rmsd_A",
        ):
            best_col = f"{col}_best"
            if best_col in categories_updated.columns:
                categories_updated[col] = categories_updated[best_col].combine_first(
                    categories_updated[col]
                )
                categories_updated = categories_updated.drop(columns=[best_col])
        for col in ("closest_cross_tag_file_a", "closest_cross_tag_file_b"):
            categories_updated[col] = (
                categories_updated[col].astype(str).replace("nan", "")
            )

    categories_updated = categories_updated.sort_values(["category", "info_file"])
    categories_updated.to_csv(out_dir / "complex_categories.csv", index=False)


def process_manual_curation_categories(
    info_dir: str | Path,
    ligand_structures_dir: str | Path,
    manual_dir: str | Path,
    *,
    out_dir: str | Path,
    csv_suffix: str = ".csv",
    file_column: str = "file",
) -> None:
    """Build ``complex_categories.csv``: invalid vs training, plus closest-pair RMSD."""
    from export_threshold_training_table import (
        group_key_from_info_file,
        load_manual_labels_from_curation_dir,
    )

    info_dir = Path(info_dir)
    manual_dir = Path(manual_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    _, by_subtype, by_invalid = load_manual_labels_from_curation_dir(manual_dir)
    invalid_only = {k for k in by_invalid if k not in by_subtype}

    category_rows: list[dict[str, object]] = []
    all_distance_tables: list[pd.DataFrame] = []
    all_errors: list[dict[str, object]] = []

    for info_csv in sorted(info_dir.iterdir()):
        if not info_csv.name.endswith(csv_suffix):
            continue

        info_file = info_csv.name
        group_key = group_key_from_info_file(info_file)
        cat = categorize_complex_by_atom_counts(
            info_csv, ligand_structures_dir, file_column=file_column
        )

        if group_key in invalid_only:
            cat["category"] = "invalid"
            cat["manual_label"] = "invalid"
            cat["category_reason"] = "manual_curation_invalid"
            category_rows.append(cat)
            continue

        manual = by_subtype.get(group_key) if group_key else None
        cat["category"] = "1"
        cat["manual_label"] = manual
        cat["category_reason"] = (
            f"manual_curation_{manual}" if manual else "manual_curation_training"
        )

        dist_df, errors = compute_atom_distances_for_complex(
            info_csv, ligand_structures_dir, file_column=file_column
        )
        all_errors.extend(errors)
        if has_cross_tag_distance_pairs(dist_df):
            all_distance_tables.append(dist_df)
        else:
            failure = cross_tag_feature_failure_reason(cat, errors)
            cat["category_reason"] = f"{cat['category_reason']};{failure}"

        category_rows.append(cat)

    _finalize_and_write_categories(
        category_rows,
        all_distance_tables,
        all_errors,
        out_dir,
    )
