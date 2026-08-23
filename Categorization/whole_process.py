#!/usr/bin/env python3
"""
Binding-mode categorization from manual curation + closest-pair RMSD.

Invalid groups come only from manual curation (``manualcheck_invalid`` only groups).
All other groups are used for custom RMSD scoring (labels 1.1 / 1.2 / 1.3 / 2);
the closest m0–m1 pair (smallest mean atom distance) and its RMSD are stored.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

_CATEGORIZATION_DIR = Path(__file__).resolve().parent
if str(_CATEGORIZATION_DIR) not in sys.path:
    sys.path.insert(0, str(_CATEGORIZATION_DIR))

import paths
from atom_count_and_distance import process_manual_curation_categories
from export_threshold_training_table import load_manual_labels_from_curation_dir


@dataclass(frozen=True)
class Paths:
    info_dir: Path
    ligand_structures_dir: Path
    out_dir: Path
    manual_curation_dir: Path


def _print_manual_subtype_counts(manual_dir: Path) -> None:
    _, by_subtype, by_invalid = load_manual_labels_from_curation_dir(manual_dir)
    invalid_only = {k: v for k, v in by_invalid.items() if k not in by_subtype}

    print("\nManual curation — training label sizes:")
    if not by_subtype:
        print("  (no subtype labels found)")
    else:
        counts = pd.Series(by_subtype).value_counts().sort_index()
        for label, count in counts.items():
            print(f"  {label}: {int(count)}")
        print(f"  total: {int(counts.sum())}")

    print("\nManual curation — invalid only (excluded from custom RMSD scoring):")
    print(f"  total: {len(invalid_only)}")


def run(paths_cfg: Paths, *, report: bool = True) -> Path:
    paths_cfg.out_dir.mkdir(parents=True, exist_ok=True)

    if report:
        print(f"info_dir: {paths_cfg.info_dir.resolve()}")
        print(f"ligand_structures_dir: {paths_cfg.ligand_structures_dir.resolve()}")
        print(f"out_dir: {paths_cfg.out_dir.resolve()}")
        print(f"manual_curation_dir: {paths_cfg.manual_curation_dir.resolve()}")
        _print_manual_subtype_counts(paths_cfg.manual_curation_dir)

    process_manual_curation_categories(
        paths_cfg.info_dir,
        paths_cfg.ligand_structures_dir,
        paths_cfg.manual_curation_dir,
        out_dir=paths_cfg.out_dir,
        csv_suffix=".csv",
        file_column="file",
    )

    categories_csv = paths_cfg.out_dir / "complex_categories.csv"
    categories = pd.read_csv(categories_csv)

    if report:
        cat = categories["category"].astype(str).str.strip().str.lower()
        n_invalid = int(cat.isin(["invalid"]).sum())
        n_train = int((~cat.isin(["invalid"])).sum())
        print("\n--- Training vs invalid (manual curation) + closest-pair RMSD ---")
        print(f"Total groups: {len(categories)}")
        print(f"  Training groups: {n_train}, Invalid: {n_invalid}")
        print(
            "\nNext: export training table and score custom RMSD bins —\n"
            "  python Categorization/export_threshold_training_table.py\n"
            "  python Categorization/train_category_thresholds.py"
        )
        print(f"\nWrote {categories_csv}")

    return categories_csv


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Manual training vs invalid plus closest-pair RMSD."
        )
    )
    parser.add_argument(
        "--info-dir",
        type=Path,
        default=Path(paths.INTERACTION_TYPE_INFO_DIR),
        help="Per-group interaction-type CSV directory (pharmacophore step 3.2).",
    )
    parser.add_argument(
        "--ligand-structures-dir",
        type=Path,
        default=Path(paths.CANONICAL_LIGAND_DIR),
        help="Canonical ligand SDF directory (pharmacophore step 3.1).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(paths.CATEGORIZATION_RESULTS_DIR),
        help="Output directory for categorized results.",
    )
    parser.add_argument(
        "--manual-curation-dir",
        type=Path,
        default=Path(paths.MANUAL_CURATION_DIR),
        help="Manual curation CSV directory (default: Data/Manual_curation/).",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Do not print summary.",
    )
    args = parser.parse_args()

    manual_dir = args.manual_curation_dir.resolve()
    if not manual_dir.is_dir():
        raise FileNotFoundError(f"Manual curation directory not found: {manual_dir}")

    run(
        Paths(
            info_dir=args.info_dir,
            ligand_structures_dir=args.ligand_structures_dir,
            out_dir=args.out_dir,
            manual_curation_dir=manual_dir,
        ),
        report=not args.quiet,
    )


if __name__ == "__main__":
    main()
