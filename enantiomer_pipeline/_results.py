"""Build the result dict returned after each pipeline stage run."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Optional

from .config import PipelineConfig

STAGE_NAMES = {
    "preprocess": "1. Data Preprocess",
    "alignment": "2. Alignment",
    "pharmacophore": "3. Pharmacophore generation",
    "categorization": "4. Categorization",
}


def _count_files(directory: Path, pattern: str = "*") -> int:
    if not directory.is_dir():
        return 0
    return sum(1 for _ in directory.glob(pattern) if _.is_file())


def _read_text(path: Path) -> Optional[str]:
    if path.is_file():
        return path.read_text(encoding="utf-8", errors="replace").strip()
    return None


def _parse_key_value_lines(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in text.splitlines():
        if ":" in line:
            key, _, value = line.partition(":")
            result[key.strip()] = value.strip()
    return result


def _parse_grouping_summary(text: str) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for line in text.splitlines():
        m = re.match(r"Paired groups:\s*(\d+)\s*\((\d+) entries\)", line)
        if m:
            out["paired_groups"] = int(m.group(1))
            out["paired_entries"] = int(m.group(2))
        m = re.match(r"Single groups:\s*(\d+)\s*\((\d+) entries\)", line)
        if m:
            out["single_groups"] = int(m.group(1))
            out["single_entries"] = int(m.group(2))
        m = re.match(r"Skipped \(no UniProt\):\s*(\d+)", line)
        if m:
            out["skipped_no_uniprot"] = int(m.group(1))
    return out


def _line_count(path: Path) -> int:
    text = _read_text(path)
    if not text:
        return 0
    return sum(1 for line in text.splitlines() if line.strip())


def _status_from_checks(checks: list[tuple[str, bool]]) -> str:
    if all(ok for _, ok in checks):
        return "complete"
    if any(ok for _, ok in checks):
        return "partial"
    return "missing"


def _preprocess_result(config: PipelineConfig) -> dict[str, Any]:
    p = config.paths
    data = config.data_dir

    ligand_sdf = Path(p.LIGAND_SDF_DIR)
    separated = Path(p.SEPARATED_SDF_DIR)
    inchi_dir = Path(p.INCHI_DIR)
    chiral_dir = Path(p.CHIRAL_RESULTS_DIR)
    prep_dir = Path(p.PREPWIZARD_OUTPUT_DIR)

    inchi_summary = _read_text(Path(p.TIMER_DIR) / "summary_sdf_to_inchi.txt")
    chiral_summary = _read_text(chiral_dir / "summary.txt")
    db_ids = Path(p.DB_IDS_WITH_M_LAYER)

    checks = [
        ("ligand_sdfs", ligand_sdf.is_dir() and _count_files(ligand_sdf, "*.sdf") > 0),
        ("inchis", (inchi_dir / "all_inchis_dict.txt").is_file()),
        ("chiral_classification", chiral_summary is not None),
        ("filtered_enantiomers", Path(p.FILTERED_ENANTIOMER_WITH_UNIPROT_TXT).is_file()),
        ("db_ids", db_ids.is_file()),
        ("prepwizard", prep_dir.is_dir() and _count_files(prep_dir, "*_prepped.mae") > 0),
    ]
    status = _status_from_checks(checks)

    metrics: dict[str, Any] = {
        "ligand_sdf_files": _count_files(ligand_sdf, "*.sdf"),
        "separated_sdf_files": _count_files(separated, "*.sdf"),
        "inchis_converted": _line_count(inchi_dir / "all_inchis_dict.txt"),
        "db_ids_with_m_layer": _line_count(db_ids),
        "prepwizard_mae_files": _count_files(prep_dir, "*_prepped.mae"),
        "prepwizard_timeouts": _line_count(Path(p.PREPWIZARD_TIMEOUT_LOG)),
        "prepwizard_failures": _line_count(Path(p.PREPWIZARD_FAILED_LOG)),
    }
    if inchi_summary:
        metrics["inchi_conversion"] = _parse_key_value_lines(inchi_summary)
    if chiral_summary:
        metrics["chiral_classification_text"] = chiral_summary

    lines = [
        f"Data directory: {data}",
        f"Status: {status}",
        f"Ligand SDF files: {metrics['ligand_sdf_files']}",
        f"Separated SDF files: {metrics['separated_sdf_files']}",
        f"InChI entries: {metrics['inchis_converted']}",
        f"db_ids_with_m_layer lines: {metrics['db_ids_with_m_layer']}",
        f"PrepWizard MAE files: {metrics['prepwizard_mae_files']}",
    ]
    if chiral_summary:
        lines.extend(["", chiral_summary])

    return {
        "stage": "preprocess",
        "stage_name": STAGE_NAMES["preprocess"],
        "status": status,
        "data_dir": str(data),
        "metrics": metrics,
        "output_paths": {
            "ligand_sdfs": str(ligand_sdf),
            "inchis": str(inchi_dir),
            "chiral_results": str(chiral_dir),
            "prepwizard": str(prep_dir),
            "db_ids": str(db_ids),
        },
        "text": "\n".join(lines),
    }


def _alignment_result(config: PipelineConfig) -> dict[str, Any]:
    align = config.stage_dir("Alignment")
    paired = align / "paired_enantiomers"
    single = align / "single_enantiomers"
    parsed = align / "paired_enantiomers_parsed_chain"
    aligned = align / "paired_enantiomers_aligned_structures"
    ligand_extract = align / "Enantiomer_aligned_structure_ligandextract"
    grouping_summary_path = align / "grouping_summary.txt"

    grouping_text = _read_text(grouping_summary_path)
    grouping = _parse_grouping_summary(grouping_text) if grouping_text else {}

    checks = [
        ("paired_groups", paired.is_dir() and _count_files(paired, "*.csv") > 0),
        ("parsed_chains", parsed.is_dir() and _count_files(parsed, "*.pdb") > 0),
        ("aligned_structures", aligned.is_dir() and _count_files(aligned, "*.pdb") > 0),
        ("ligand_extract", ligand_extract.is_dir() and _count_files(ligand_extract, "*.sdf") > 0),
    ]
    status = _status_from_checks(checks)

    metrics = {
        "paired_group_csvs": _count_files(paired, "*.csv"),
        "single_group_csvs": _count_files(single, "*.csv"),
        "parsed_chain_pdbs": _count_files(parsed, "*.pdb"),
        "aligned_pdb_files": _count_files(aligned, "*.pdb"),
        "ligand_extract_sdfs": _count_files(ligand_extract, "*.sdf"),
        **grouping,
    }

    lines = [
        f"Data directory: {align}",
        f"Status: {status}",
        f"Paired group CSVs: {metrics['paired_group_csvs']}",
        f"Single group CSVs: {metrics['single_group_csvs']}",
        f"Parsed chain PDBs: {metrics['parsed_chain_pdbs']}",
        f"Aligned PDB files: {metrics['aligned_pdb_files']}",
        f"Ligand extract SDFs: {metrics['ligand_extract_sdfs']}",
    ]
    if grouping_text:
        lines.extend(["", grouping_text])

    return {
        "stage": "alignment",
        "stage_name": STAGE_NAMES["alignment"],
        "status": status,
        "data_dir": str(align),
        "metrics": metrics,
        "output_paths": {
            "paired_enantiomers": str(paired),
            "aligned_structures": str(aligned),
            "ligand_extract": str(ligand_extract),
            "grouping_summary": str(grouping_summary_path),
        },
        "text": "\n".join(lines),
    }


def _pharmacophore_result(config: PipelineConfig) -> dict[str, Any]:
    pharm = config.stage_dir("Pharmacophore")
    interaction = pharm / "Interaction_data"
    ph4_results = pharm / "Pharmacophore_generation_results"
    pocket_based = pharm / "paired_enantiomers_pocket_based"
    canonical = pharm / "Enantiomer_aligned_structure_ligandextract_canonical"

    checks = [
        ("interaction_data", interaction.is_dir() and _count_files(interaction, "*.tsv") > 0),
        ("pocket_based", pocket_based.is_dir() and _count_files(pocket_based, "*.csv") > 0),
        ("canonical_ligands", canonical.is_dir() and _count_files(canonical, "*.sdf") > 0),
    ]
    status = _status_from_checks(checks)

    metrics = {
        "interaction_tsv_files": _count_files(interaction, "*.tsv"),
        "ph4_result_files": _count_files(ph4_results, "*"),
        "pocket_based_csv_files": _count_files(pocket_based, "*.csv"),
        "canonical_ligand_sdfs": _count_files(canonical, "*.sdf"),
        "generation_failures": _line_count(pharm / "generation_failed_jobs.txt"),
        "generation_timeouts": _line_count(pharm / "generation_timeout.txt"),
    }

    lines = [
        f"Data directory: {pharm}",
        f"Status: {status}",
        f"Interaction TSV files: {metrics['interaction_tsv_files']}",
        f"Pocket-based CSV files: {metrics['pocket_based_csv_files']}",
        f"Canonical ligand SDFs: {metrics['canonical_ligand_sdfs']}",
    ]
    if metrics["generation_failures"]:
        lines.append(f"Generation failures: {metrics['generation_failures']}")
    if metrics["generation_timeouts"]:
        lines.append(f"Generation timeouts: {metrics['generation_timeouts']}")

    return {
        "stage": "pharmacophore",
        "stage_name": STAGE_NAMES["pharmacophore"],
        "status": status,
        "data_dir": str(pharm),
        "metrics": metrics,
        "output_paths": {
            "interaction_data": str(interaction),
            "pocket_based": str(pocket_based),
            "canonical_ligands": str(canonical),
        },
        "text": "\n".join(lines),
    }


def _categorization_result(config: PipelineConfig) -> dict[str, Any]:
    cat_dir = config.stage_dir("Categorization")
    feature_csv = cat_dir / "feature_table.csv"

    checks = [
        ("feature_table", feature_csv.is_file()),
    ]
    status = _status_from_checks(checks)

    total_groups = 0
    n_with_rmsd = 0
    if feature_csv.is_file():
        try:
            import pandas as pd

            df = pd.read_csv(feature_csv)
            total_groups = len(df)
            if "closest_cross_tag_rmsd_A" in df.columns:
                n_with_rmsd = int(df["closest_cross_tag_rmsd_A"].notna().sum())
        except Exception:
            n_with_rmsd = 0

    metrics = {
        "total_groups": total_groups,
        "groups_with_rmsd": n_with_rmsd,
    }

    lines = [
        f"Data directory: {cat_dir}",
        f"Status: {status}",
        f"Feature groups: {total_groups}",
        f"Groups with RMSD: {n_with_rmsd}",
    ]

    return {
        "stage": "categorization",
        "stage_name": STAGE_NAMES["categorization"],
        "status": status,
        "data_dir": str(cat_dir),
        "metrics": metrics,
        "output_paths": {
            "feature_table": str(feature_csv),
        },
        "text": "\n".join(lines),
    }


_COLLECTORS = {
    "preprocess": _preprocess_result,
    "alignment": _alignment_result,
    "pharmacophore": _pharmacophore_result,
    "categorization": _categorization_result,
}


def collect_stage_result(config: PipelineConfig, stage: str) -> dict[str, Any]:
    """Collect the result dict for a stage after it has been run."""
    if stage not in _COLLECTORS:
        raise ValueError(f"Unknown stage: {stage}. Choose from: {', '.join(_COLLECTORS)}")
    return _COLLECTORS[stage](config)
