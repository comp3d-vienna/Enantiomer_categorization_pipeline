#!/usr/bin/env python3
"""
Closest-pair RMSD for every paired-enantiomer group.

Each group gets the closest m0–m1 pair (smallest mean atom distance) and its RMSD.
The table is written as ``feature_table.csv``; SILIRID columns are added by
``export_feature_table.py``.
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
from atom_count_and_distance import write_closest_pair_rmsd_table


@dataclass(frozen=True)
class Paths:
    info_dir: Path
    ligand_structures_dir: Path
    out_dir: Path


def run(paths_cfg: Paths, *, report: bool = True) -> Path:
    paths_cfg.out_dir.mkdir(parents=True, exist_ok=True)

    if report:
        print(f"info_dir: {paths_cfg.info_dir.resolve()}")
        print(f"ligand_structures_dir: {paths_cfg.ligand_structures_dir.resolve()}")
        print(f"out_dir: {paths_cfg.out_dir.resolve()}")

    write_closest_pair_rmsd_table(
        paths_cfg.info_dir,
        paths_cfg.ligand_structures_dir,
        out_dir=paths_cfg.out_dir,
        csv_suffix=".csv",
        file_column="file",
    )

    feature_csv = paths_cfg.out_dir / "feature_table.csv"
    table = pd.read_csv(feature_csv)

    if report:
        print("\n--- Closest-pair RMSD ---")
        print(f"Feature groups: {len(table)}")
        if "closest_cross_tag_rmsd_A" in table.columns and not table.empty:
            n_rmsd = int(table["closest_cross_tag_rmsd_A"].notna().sum())
            print(f"  with RMSD: {n_rmsd}/{len(table)}")
        print(
            "\nNext: add SILIRID similarity —\n"
            "  python Categorization/export_feature_table.py"
        )
        print(f"\nWrote {feature_csv}")

    return feature_csv


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Closest-pair RMSD for paired enantiomer groups (writes feature_table.csv)."
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
        help="Output directory for the feature table.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Do not print summary.",
    )
    args = parser.parse_args()

    run(
        Paths(
            info_dir=args.info_dir,
            ligand_structures_dir=args.ligand_structures_dir,
            out_dir=args.out_dir,
        ),
        report=not args.quiet,
    )


if __name__ == "__main__":
    main()
