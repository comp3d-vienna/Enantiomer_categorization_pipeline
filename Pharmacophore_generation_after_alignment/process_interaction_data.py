#!/usr/bin/env python3
"""
Collect pharmacophore interaction features for paired enantiomer groups.

Reads interaction .tsv files (from gen_ia_ph4s_fg.py) and paired group CSVs
from the Alignment stage. Writes one feature table per group.

Paired-only workflow: no single-enantiomer separation or end-of-pipeline merge.
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
from tqdm import tqdm

import paths

FEATURE_TYPE_DICT = {
    0: "UNKNOWN",
    1: "HYDROPHOBIC",
    2: "AROMATIC",
    3: "NEGATIVE_IONIZABLE",
    4: "POSITIVE_IONIZABLE",
    5: "H_BOND_DONOR",
    6: "H_BOND_ACCEPTOR",
    7: "HALOGEN_BOND_DONOR",
    8: "HALOGEN_BOND_ACCEPTOR",
    9: "EXCLUSION_VOLUME",
}


@dataclass
class GroupMember:
    identifier: str
    inchi: str
    m_layer: str
    inchi_hash: str


@dataclass
class PairedGroup:
    group_id: str
    uniprot_id: str
    inchi_hash: str
    members: List[GroupMember] = field(default_factory=list)


def strip_inchi_stereo_enantiomer(inchi: str) -> str:
    if isinstance(inchi, str):
        idx = inchi.find("/m")
        if idx != -1:
            return inchi[:idx]
        idx = inchi.find("/s")
        if idx != -1:
            return inchi[:idx]
    return inchi


def get_inchi_group_enantiomer(inchi: str) -> str:
    return strip_inchi_stereo_enantiomer(inchi)


def extract_identifier_from_tsv_filename(fname: str) -> Optional[str]:
    """interaction_{pdb}_{het}_{chain}_{resseq}_{m_layer}.tsv -> pdb_het_chain_resseq"""
    basename = os.path.splitext(os.path.basename(fname))[0]
    if not basename.startswith("interaction_"):
        return None
    core = basename[len("interaction_") :]
    parts = core.split("_")
    if len(parts) < 2:
        return None
    return "_".join(parts[:-1])


def get_column(df: pd.DataFrame, substring: str) -> Optional[str]:
    for col in df.columns:
        if substring in col:
            return col
    return None


def parse_atom_indices(cell) -> List[str]:
    if pd.isna(cell):
        return []
    cell_str = str(cell).strip()
    if not cell_str:
        return []
    if "," not in cell_str:
        return [cell_str]
    return [s.strip() for s in cell_str.split(",") if s.strip()]


def build_atom_to_details(df_clean, idx_col, typ_col, res_col, fg_col) -> Dict:
    atom_to_details = {}
    for _, row in df_clean.iterrows():
        atom_indices = parse_atom_indices(row[idx_col])
        ftype = FEATURE_TYPE_DICT.get(row[typ_col], str(row[typ_col]))
        fg = row[fg_col]
        pocket = row[res_col]
        for ai in atom_indices:
            atom_to_details.setdefault(ai, set()).add((ftype, pocket, fg))
    return atom_to_details


def atom_sort_key(a):
    return int(a) if str(a).isdigit() else a


def load_paired_groups(paired_dir: str) -> List[PairedGroup]:
    groups: List[PairedGroup] = []
    for csv_path in sorted(glob.glob(os.path.join(paired_dir, "*.csv"))):
        df = pd.read_csv(csv_path)
        if df.empty:
            continue
        uniprot_id = str(df["uniprot_id"].iloc[0])
        inchi_hash = str(df["inchi_hash"].iloc[0])
        group_id = ""
        if "group_id" in df.columns:
            group_id = str(df["group_id"].iloc[0]).strip()
        if not group_id:
            # Fallback: parse from filename like g000001_inchi_<hash>.csv
            base = os.path.splitext(os.path.basename(csv_path))[0]
            if "_inchi_" in base:
                group_id = base.split("_inchi_", 1)[0]
        if not group_id:
            # Last resort: unsafe (can be long), but preserves behavior
            group_id = uniprot_id
        group = PairedGroup(group_id=group_id, uniprot_id=uniprot_id, inchi_hash=inchi_hash)
        seen: Set[str] = set()
        for _, row in df.iterrows():
            ident = str(row["identifier"])
            if ident in seen:
                continue
            seen.add(ident)
            inchi = str(row["InChI"])
            m_layer = str(row.get("m_layer", ""))
            group.members.append(
                GroupMember(
                    identifier=ident,
                    inchi=inchi,
                    m_layer=m_layer,
                    inchi_hash=inchi_hash,
                )
            )
        groups.append(group)
    return groups


def index_interaction_files(interaction_dir: str) -> Dict[str, str]:
    """Map ligand identifier (without m_layer) -> interaction .tsv path."""
    index: Dict[str, str] = {}
    if not os.path.isdir(interaction_dir):
        return index
    for fname in os.listdir(interaction_dir):
        if not fname.endswith(".tsv"):
            continue
        ident = extract_identifier_from_tsv_filename(fname)
        if ident:
            index[ident] = os.path.join(interaction_dir, fname)
    return index


def load_valid_interaction_data(
    matching: List[Tuple[str, str, str]],
    skipped_log: List[str],
) -> List[Tuple[str, str, str, pd.DataFrame, Dict]]:
    """Return list of (tsv_fname, identifier, enantiomer_hash, df, atom_to_details)."""
    valid = []
    for tsv_path, identifier, enantiomer_hash in matching:
        fname = os.path.basename(tsv_path)
        try:
            df = pd.read_csv(tsv_path, sep="\t")
        except Exception as exc:
            skipped_log.append(f"Failed to read {fname}: {exc}")
            continue
        idx_col = get_column(df, "Ligand Atom Indices")
        typ_col = get_column(df, "Pharm. Feature Type")
        res_col = get_column(df, "Pocket Residues")
        fg_col = get_column(df, "Ligand Func. Group")
        if not all([idx_col, typ_col, res_col, fg_col]) or df.empty:
            skipped_log.append(f"{fname}: missing columns or empty")
            continue
        df_clean = df.dropna(subset=[idx_col])
        if df_clean.empty:
            skipped_log.append(f"{fname}: no valid atom indices")
            continue
        atom_to_details = build_atom_to_details(df_clean, idx_col, typ_col, res_col, fg_col)
        if not atom_to_details:
            skipped_log.append(f"{fname}: no parseable atom indices")
            continue
        valid.append((fname, identifier, enantiomer_hash, df_clean, atom_to_details))
    return valid


def collect_group_interactions(
    valid_data: List[Tuple[str, str, str, pd.DataFrame, Dict]],
) -> pd.DataFrame:
    """Build per-atom output rows for every structure in the group."""
    all_atom_indices = set().union(*(item[4].keys() for item in valid_data))
    output_rows = []

    for atom_idx in sorted(all_atom_indices, key=atom_sort_key):
        for fname, _identifier, enantiomer_hash, _df, atom_to_details in valid_data:
            details = atom_to_details.get(atom_idx, set())
            for ftype, pocket, functional_group in details:
                output_rows.append(
                    {
                        "atom_index": atom_idx,
                        "file": fname,
                        "feature_type": ftype,
                        "pocket_residues": pocket,
                        "functional_group": functional_group,
                        "enantiomer_short_hash": enantiomer_hash,
                    }
                )

    df_all = pd.DataFrame(output_rows)
    if not df_all.empty:
        df_all = df_all.sort_values(["atom_index", "file"])
    return df_all


def process_paired_interaction_data(
    paired_dir: Optional[str] = None,
    interaction_dir: Optional[str] = None,
    output_dir: Optional[str] = None,
    log_path: Optional[str] = None,
) -> None:
    paired_dir = paired_dir or paths.PAIRED_ENANTIOMERS_DIR
    interaction_dir = interaction_dir or paths.INTERACTION_DATA_DIR
    output_dir = output_dir or paths.PAIRED_FEATURE_OUTPUT_DIR
    log_path = log_path or paths.PROCESSING_LOG

    paths.ensure_dir(output_dir)

    groups = load_paired_groups(paired_dir)
    interaction_index = index_interaction_files(interaction_dir)
    all_tsv_files = {
        os.path.basename(p) for p in interaction_index.values()
    }

    skipped_log: List[str] = []
    files_touched: Set[str] = set()

    count_skipped_single = 0
    count_no_files = 0
    count_missing_cols = 0
    successful_groups = 0

    print(f"Paired groups to process: {len(groups)}")
    print(f"Interaction TSV files found: {len(all_tsv_files)}")

    for group in tqdm(groups, desc="Processing paired groups"):
        identifier_to_member = {m.identifier: m for m in group.members}
        matching: List[Tuple[str, str, str]] = []

        for ident, member in identifier_to_member.items():
            tsv_path = interaction_index.get(ident)
            if not tsv_path:
                continue
            enantiomer_hash = hashlib.sha1(
                get_inchi_group_enantiomer(member.inchi).encode()
            ).hexdigest()[:10]
            matching.append((tsv_path, ident, enantiomer_hash))

        if not matching:
            count_no_files += 1
            id_list = ", ".join(sorted(identifier_to_member.keys()))
            skipped_log.append(
                f"No interaction files for {group.group_id} inchi_{group.inchi_hash}: {id_list}"
            )
            continue

        valid_data = load_valid_interaction_data(matching, skipped_log)
        for item in valid_data:
            files_touched.add(item[0])

        if len(valid_data) < 2:
            count_skipped_single += 1
            skipped_log.append(
                f"Group {group.group_id}_inchi_{group.inchi_hash}: "
                f"only {len(valid_data)} valid structure(s); need >= 2 for pair comparison"
            )
            continue

        if len(valid_data) < len(matching):
            count_missing_cols += 1

        successful_groups += 1
        df_all = collect_group_interactions(valid_data)
        dest = os.path.join(output_dir, f"{group.group_id}_inchi_{group.inchi_hash}.csv")
        df_all.to_csv(dest, index=False)

    audit_lines = [
        "",
        "=" * 40,
        "FINAL AUDIT REPORT (paired enantiomers only)",
        "=" * 40,
        f"Total paired groups:                 {len(groups)}",
        f"Total interaction TSV files:         {len(all_tsv_files)}",
        f"TSV files used:                      {len(files_touched)}",
        f"Groups with output:                  {successful_groups}",
        f"Groups with no interaction files:    {count_no_files}",
        f"Groups skipped (< 2 structures):   {count_skipped_single}",
        f"Groups with some invalid TSVs:       {count_missing_cols}",
        "-" * 40,
        f"Output directory: {output_dir}",
        f"Log file: {log_path}",
    ]

    with open(log_path, "w") as f:
        if skipped_log:
            f.write("\n".join(skipped_log))
            f.write("\n\n")
        f.write("\n".join(audit_lines))

    for line in audit_lines:
        print(line)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Process pharmacophore interaction data for paired enantiomer groups."
    )
    parser.add_argument(
        "--paired-dir",
        default=paths.PAIRED_ENANTIOMERS_DIR,
        help="Directory with paired group CSVs from Alignment step 2.1",
    )
    parser.add_argument(
        "--interaction-dir",
        default=paths.INTERACTION_DATA_DIR,
        help="Directory containing interaction_*.tsv files",
    )
    parser.add_argument(
        "--output-dir",
        default=paths.PAIRED_FEATURE_OUTPUT_DIR,
        help="Output directory for per-group feature CSVs",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    process_paired_interaction_data(
        paired_dir=args.paired_dir,
        interaction_dir=args.interaction_dir,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
