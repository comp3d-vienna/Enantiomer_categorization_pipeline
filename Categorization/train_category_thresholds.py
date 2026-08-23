#!/usr/bin/env python3
"""
Score custom RMSD bins on the full threshold training table.

Bins (not from a decision tree):
  RMSD <= 1.25        -> 1.1
  1.25 < RMSD <= 2.5  -> 1.2
  2.5  < RMSD <= 11   -> 1.3
  RMSD > 11           -> 2

Writes purity and recall for each category and for the whole dataset.
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
from export_threshold_training_table import TREE_FEATURES, normalize_manual_label

DEFAULT_LABELS: tuple[str, ...] = tuple(paths.THRESHOLD_TUNING_MANUAL_LABELS)

RMSD_FEATURE = "closest_cross_tag_rmsd_A"
FIXED_RMSD_CUTOFFS: tuple[float, float, float] = (1.25, 2.5, 11.0)
FIXED_RMSD_CATEGORY_BY_BIN: tuple[str, ...] = ("1.1", "1.2", "1.3", "2")

REPORT_NAME = "custom_rmsd_threshold_report.txt"
METRICS_CSV_NAME = "custom_rmsd_thresholds.csv"


def _default_training_csv() -> Path:
    return Path(paths.CATEGORIZATION_RESULTS_DIR) / "threshold_training_table.csv"


def _default_out_dir() -> Path:
    return Path(paths.CATEGORIZATION_RESULTS_DIR)


def load_training_data(
    csv_path: Path,
    *,
    manual_labels: tuple[str, ...],
    tree_features: tuple[str, ...] = TREE_FEATURES,
) -> pd.DataFrame:
    df = pd.read_csv(csv_path, dtype={"manual_label": "string"})
    if "manual_label" not in df.columns:
        raise ValueError(f"Missing column 'manual_label' in {csv_path}")
    df = df.loc[df["manual_label"].notna()].copy()
    df["manual_label"] = df["manual_label"].map(normalize_manual_label)
    df = df.loc[df["manual_label"].isin(manual_labels)].copy()
    if "pipeline_category" in df.columns:
        pipeline_cat = df["pipeline_category"].astype(str).str.strip()
        df = df.loc[~pipeline_cat.isin({"invalid"})].copy()
    missing = [c for c in tree_features if c not in df.columns]
    if missing:
        raise ValueError(f"Training table missing feature columns: {missing}")
    df = df.dropna(subset=list(tree_features))
    if df.empty:
        raise ValueError("No training rows after filtering labels and dropping NaN features.")
    return df


def assign_category_by_fixed_rmsd(
    rmsd: float,
    *,
    cutoffs: tuple[float, float, float] = FIXED_RMSD_CUTOFFS,
    categories: tuple[str, ...] = FIXED_RMSD_CATEGORY_BY_BIN,
) -> str:
    lo, mid, hi = cutoffs
    if rmsd <= lo:
        return categories[0]
    if rmsd <= mid:
        return categories[1]
    if rmsd <= hi:
        return categories[2]
    return categories[3]


def _rmsd_bin_label(
    category: str,
    cutoffs: tuple[float, float, float] = FIXED_RMSD_CUTOFFS,
) -> str:
    lo, mid, hi = cutoffs
    return {
        FIXED_RMSD_CATEGORY_BY_BIN[0]: f"RMSD <= {lo:g}",
        FIXED_RMSD_CATEGORY_BY_BIN[1]: f"{lo:g} < RMSD <= {mid:g}",
        FIXED_RMSD_CATEGORY_BY_BIN[2]: f"{mid:g} < RMSD <= {hi:g}",
        FIXED_RMSD_CATEGORY_BY_BIN[3]: f"RMSD > {hi:g}",
    }.get(category, category)


def evaluate_fixed_rmsd_cutoffs(
    data: pd.DataFrame,
    *,
    rmsd_col: str = RMSD_FEATURE,
    cutoffs: tuple[float, float, float] = FIXED_RMSD_CUTOFFS,
    categories: tuple[str, ...] = FIXED_RMSD_CATEGORY_BY_BIN,
) -> pd.DataFrame:
    """Purity and recall on the full labeled table (every row is assigned a bin)."""
    work = data.loc[data[rmsd_col].notna()].copy()
    work["predicted_label"] = (
        work[rmsd_col]
        .astype(float)
        .map(lambda v: assign_category_by_fixed_rmsd(v, cutoffs=cutoffs, categories=categories))
    )
    rows: list[dict[str, object]] = []
    for category in categories:
        predicted = work.loc[work["predicted_label"] == category]
        n_predicted = int(len(predicted))
        n_correct = int((predicted["manual_label"] == category).sum())
        n_true = int((work["manual_label"] == category).sum())
        rows.append(
            {
                "category": category,
                "rmsd_bin": _rmsd_bin_label(category, cutoffs),
                "n_predicted": n_predicted,
                "n_true": n_true,
                "n_correct": n_correct,
                "purity": n_correct / n_predicted if n_predicted else 0.0,
                "recall": n_correct / n_true if n_true else 0.0,
            }
        )

    n_total = int(len(work))
    n_ok = int((work["manual_label"] == work["predicted_label"]).sum())
    rows.append(
        {
            "category": "ALL",
            "rmsd_bin": "whole dataset",
            "n_predicted": n_total,
            "n_true": n_total,
            "n_correct": n_ok,
            "purity": n_ok / n_total if n_total else 0.0,
            "recall": n_ok / n_total if n_total else 0.0,
        }
    )
    return pd.DataFrame(rows)


def format_report(metrics: pd.DataFrame) -> str:
    lo, mid, hi = FIXED_RMSD_CUTOFFS
    lines = [
        "Custom RMSD thresholds",
        "=" * 60,
        (
            f"RMSD <= {lo:g} -> 1.1; {lo:g} < RMSD <= {mid:g} -> 1.2; "
            f"{mid:g} < RMSD <= {hi:g} -> 1.3; RMSD > {hi:g} -> 2"
        ),
        "",
        "Purity = n_correct / n_predicted in that bin.",
        "Recall = n_correct / n_true labels in the whole evaluated table.",
        "",
    ]
    if metrics.empty:
        lines.append("(no metrics)")
        return "\n".join(lines)

    overall = metrics.loc[metrics["category"] == "ALL"]
    per_category = metrics.loc[metrics["category"] != "ALL"]
    if not overall.empty:
        row = overall.iloc[0]
        lines.append(
            f"Whole dataset ({int(row['n_true'])} groups): "
            f"purity={float(row['purity']):.1%}, "
            f"recall={float(row['recall']):.1%} "
            f"({int(row['n_correct'])}/{int(row['n_true'])} correct)"
        )
        lines.append("")
    for _, row in per_category.iterrows():
        lines.append(
            f"Category {row['category']} ({row['rmsd_bin']}): "
            f"purity={float(row['purity']):.1%} "
            f"({int(row['n_correct'])}/{int(row['n_predicted'])}), "
            f"recall={float(row['recall']):.1%} "
            f"({int(row['n_correct'])}/{int(row['n_true'])})"
        )
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Score custom RMSD category bins on the full training table."
    )
    parser.add_argument("--training-csv", type=Path, default=_default_training_csv())
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--labels", nargs="+", default=list(DEFAULT_LABELS))
    args = parser.parse_args()

    training_csv = args.training_csv.resolve()
    if not training_csv.is_file():
        raise FileNotFoundError(f"Training CSV not found: {training_csv}")

    out_dir = (args.out_dir or _default_out_dir()).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manual_labels = tuple(str(x).strip() for x in args.labels)

    data = load_training_data(
        training_csv,
        manual_labels=manual_labels,
        tree_features=TREE_FEATURES,
    )
    raw = pd.read_csv(training_csv)
    n_skipped = len(raw) - len(data)
    print(f"Loaded {len(data)} group(s) from {training_csv}")
    if n_skipped:
        print(
            f"  ({n_skipped} curated group(s) skipped — missing RMSD; "
            "see threshold_training_table_incomplete_features.csv)"
        )

    metrics = evaluate_fixed_rmsd_cutoffs(data)
    report = format_report(metrics)
    report_path = out_dir / REPORT_NAME
    csv_path = out_dir / METRICS_CSV_NAME
    report_path.write_text(report, encoding="utf-8")
    metrics.to_csv(csv_path, index=False)

    print(report)
    print(f"Wrote {report_path}")
    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    main()
