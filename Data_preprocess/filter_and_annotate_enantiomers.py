#!/usr/bin/env python3
"""
Steps 1.2–1.3b of the Data Preprocess pipeline.

Replaces notebook cells in process_data.ipynb for:
  - Drug-likeness filtering of enantiomer structures
  - Syncing filtered enantiomer pair index
  - UniProt ID mapping via RCSB REST API
  - Generating db_ids_with_m_layer.txt
"""

import argparse
import os
import re

from rdkit import Chem

import paths
from utils import filter_molecules, get_uniprotID


def filter_drug_like_sdf(input_sdf, filtered_sdf, filtered_out_sdf):
    """Apply drug-likeness filters and write passed/rejected SDF files."""
    ensure = paths.ensure_parent_dir
    ensure(filtered_sdf)
    ensure(filtered_out_sdf)

    suppl = Chem.SDMolSupplier(input_sdf, removeHs=False)
    mols = [mol for mol in suppl]
    filter_molecules(mols, filtered_sdf, filtered_out_sdf)


def sync_filtered_enantiomers_txt(filtered_sdf_path, enantiomers_txt, output_txt):
    """Keep enantiomer pairs whose first entry passed the drug-likeness filter."""
    paths.ensure_parent_dir(output_txt)

    mol_name_set = set()
    for mol in Chem.SDMolSupplier(filtered_sdf_path, removeHs=False):
        if mol is None:
            continue
        if mol.HasProp("_Name"):
            mol_name_set.add(mol.GetProp("_Name").strip())

    with open(enantiomers_txt, "r") as f:
        lines = f.readlines()

    matched_groups = []
    i = 0
    while i < len(lines):
        if i + 1 >= len(lines):
            break
        line1 = lines[i].rstrip("\n")
        line2 = lines[i + 1].rstrip("\n")
        entry1_name = line1.split()[0] if line1.strip() else None
        if entry1_name and entry1_name in mol_name_set:
            matched_groups.append([line1, line2])
        i += 3

    with open(output_txt, "w") as fout:
        for idx, group in enumerate(matched_groups):
            fout.write(group[0] + "\n")
            fout.write(group[1] + "\n")
            if idx != len(matched_groups) - 1:
                fout.write("\n")

    print(f"Wrote {len(matched_groups)} filtered enantiomer pairs to {output_txt}")


def map_uniprot_ids(filtered_enantiomers_txt, output_txt):
    """Query RCSB for UniProt IDs and write annotated mapping file."""
    paths.ensure_parent_dir(output_txt)

    identifiers = set()
    with open(filtered_enantiomers_txt, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 1 or "_" not in parts[0]:
                continue
            identifiers.add(parts[0])

    print(f"Querying UniProt IDs for {len(identifiers)} identifiers...")
    results = get_uniprotID(identifiers)

    unique_lines = set()
    count_without_uniprot = 0
    with open(filtered_enantiomers_txt, "r") as in_f, open(output_txt, "w") as out_f:
        for line in in_f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 1 or "_" not in parts[0]:
                continue
            identifier = parts[0]
            inchi = parts[1]
            pdb_id = identifier.split("_")[0]
            uniprot_ids = results.get(identifier, {}).get("uniprot_ids", [])
            chain_id = results.get(identifier, {}).get("chain_id", [])
            uniprot_str = ",".join(uniprot_ids) if uniprot_ids else "NO_UNIPROT_INFO"
            chain_str = ",".join(chain_id) if chain_id else ""
            if uniprot_str == "NO_UNIPROT_INFO":
                count_without_uniprot += 1
            unique_line = f"{identifier}\t{pdb_id}\t{inchi}\t{uniprot_str}\t{chain_str}\n"
            if unique_line not in unique_lines:
                out_f.write(unique_line)
                unique_lines.add(unique_line)

    print(f"Wrote UniProt mapping to {output_txt}")
    print(f"Entries without UniProt info: {count_without_uniprot}")


def get_m_layer(inchi):
    """Extract the /m layer from an InChI string (e.g. 'm0', 'm1')."""
    if "/m" not in inchi:
        return ""
    match = re.search(r"/m([^/]+)", inchi)
    return "m" + match.group(1) if match else ""


def generate_db_ids_with_m_layer(filtered_sdf_path, output_txt):
    """Append InChI /m layer to each ligand identifier."""
    paths.ensure_parent_dir(output_txt)

    db_ids = []
    unique_pdbids = set()
    missing_m_layer = 0

    for mol in Chem.SDMolSupplier(filtered_sdf_path, removeHs=False):
        if mol is None:
            continue
        inchi_val = mol.GetProp("InChI") if mol.HasProp("InChI") else ""
        name_base = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        m_layer = get_m_layer(inchi_val)
        if not m_layer:
            missing_m_layer += 1
            print(f"No m layer for {name_base}")
        db_id_full = f"{name_base}_{m_layer}" if m_layer else name_base
        db_ids.append(db_id_full)
        if "_" in name_base:
            unique_pdbids.add(name_base.split("_")[0])
        else:
            unique_pdbids.add(name_base)

    with open(output_txt, "w") as f:
        for db_id in db_ids:
            f.write(f"{db_id}\n")

    print(f"Wrote {len(db_ids)} IDs to {output_txt}")
    print(f"Unique PDB IDs: {len(unique_pdbids)}")
    if missing_m_layer:
        print(f"Structures missing /m layer: {missing_m_layer}")


def run_all(skip_uniprot=False):
    os.makedirs(paths.DRUG_LIKENESS_DIR, exist_ok=True)

    print("=== Step 1.2: Filter drug-like enantiomer structures ===")
    filter_drug_like_sdf(
        paths.ENANTIOMER_STRUCTURE_SDF,
        paths.FILTERED_ENANTIOMER_STRUCTURE_SDF,
        paths.FILTERED_OUT_ENANTIOMER_STRUCTURE_SDF,
    )

    print("=== Step 1.2: Filter drug-like enantiomer representatives ===")
    if os.path.exists(paths.ENANTIOMER_REPRESENTATIVE_SDF):
        filter_drug_like_sdf(
            paths.ENANTIOMER_REPRESENTATIVE_SDF,
            paths.FILTERED_ENANTIOMER_REPRESENTATIVE_SDF,
            paths.FILTERED_OUT_ENANTIOMER_REPRESENTATIVE_SDF,
        )
    else:
        print(f"Skipping representatives (not found): {paths.ENANTIOMER_REPRESENTATIVE_SDF}")

    print("=== Step 1.2: Sync filtered enantiomer pair index ===")
    sync_filtered_enantiomers_txt(
        paths.FILTERED_ENANTIOMER_STRUCTURE_SDF,
        paths.ENANTIOMERS_TXT,
        paths.FILTERED_ENANTIOMERS_TXT,
    )

    if not skip_uniprot:
        print("=== Step 1.3: Map UniProt IDs ===")
        map_uniprot_ids(
            paths.FILTERED_ENANTIOMERS_TXT,
            paths.FILTERED_ENANTIOMER_WITH_UNIPROT_TXT,
        )
    else:
        print("=== Step 1.3: Skipped (--skip-uniprot) ===")

    print("=== Step 1.3b: Generate db_ids_with_m_layer.txt ===")
    generate_db_ids_with_m_layer(
        paths.FILTERED_ENANTIOMER_STRUCTURE_SDF,
        paths.DB_IDS_WITH_M_LAYER,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Filter drug-like enantiomers, map UniProt IDs, and generate db_ids_with_m_layer.txt"
    )
    parser.add_argument(
        "--skip-uniprot",
        action="store_true",
        help="Skip RCSB UniProt API queries (useful for offline testing)",
    )
    args = parser.parse_args()
    run_all(skip_uniprot=args.skip_uniprot)


if __name__ == "__main__":
    main()
