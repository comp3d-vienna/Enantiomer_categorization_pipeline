#!/usr/bin/env python3
"""
Alignment steps 2.1 and 2.2:

  2.1  Group filtered enantiomers into paired_enantiomers (m0 + m1) vs single_enantiomers,
       and convert PrepWizard MAE files to PDB for chain parsing
  2.2  Extract protein chains with BioPython and write list_with_uniprot_id.txt
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import subprocess
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

from Bio.PDB import PDBIO, PDBParser

import paths


def inchi_before_m_layer(inchi: str) -> str:
    """InChI string before the /m layer (keeps /t stereochemistry)."""
    idx = inchi.find("/m")
    return inchi[:idx] if idx != -1 else inchi


def inchi_hash(inchi_before_m: str) -> str:
    return hashlib.sha1(inchi_before_m.encode()).hexdigest()[:10]


def get_m_layer(inchi: str) -> str:
    match = re.search(r"/m(\d+)", inchi)
    return f"m{match.group(1)}" if match else ""


def load_entries(mapping_file: str) -> Tuple[List[dict], int]:
    entries = []
    skipped = 0
    with open(mapping_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                skipped += 1
                continue
            identifier, pdb_id, inchi, uniprot_id, chain_id = parts[:5]
            if not uniprot_id or uniprot_id == "NO_UNIPROT_INFO":
                skipped += 1
                continue
            entries.append(
                {
                    "identifier": identifier,
                    "pdb_id": pdb_id.lower(),
                    "inchi": inchi,
                    "uniprot_id": uniprot_id,
                    "chain_id": chain_id,
                    "m_layer": get_m_layer(inchi),
                    "inchi_before_m": inchi_before_m_layer(inchi),
                }
            )
    return entries, skipped


def group_entries(entries: List[dict]) -> Tuple[Dict, Dict]:
    """Group by (uniprot_id, inchi_hash). UniProt string kept as-is (may contain commas)."""
    groups = defaultdict(list)
    for entry in entries:
        ih = inchi_hash(entry["inchi_before_m"])
        entry["inchi_hash"] = ih
        key = (entry["uniprot_id"], ih)
        groups[key].append(entry)

    paired = {}
    single = {}
    for key, members in groups.items():
        m_layers = {m["m_layer"] for m in members if m["m_layer"]}
        if "m0" in m_layers and "m1" in m_layers:
            paired[key] = members
        else:
            single[key] = members
    return paired, single


def make_group_id_map(paired: dict, single: dict) -> Dict[Tuple[str, str], str]:
    """
    Assign stable numeric group IDs (g000001...) per (uniprot_id, inchi_hash).

    Stability rule: sort by (uniprot_id, inchi_hash) and enumerate from 1.
    """
    keys = sorted(set(paired.keys()) | set(single.keys()))
    return {key: f"g{idx:06d}" for idx, key in enumerate(keys, start=1)}


def write_group_manifest(group_id_map: Dict[Tuple[str, str], str], output_path: str) -> None:
    """Write group_id -> uniprot_id + inchi_hash for traceability."""
    paths.ensure_dir(os.path.dirname(output_path))
    with open(output_path, "w") as f:
        f.write("group_id\tuniprot_id\tinchi_hash\n")
        for (uniprot_id, ih), gid in sorted(group_id_map.items(), key=lambda x: x[1]):
            f.write(f"{gid}\t{uniprot_id}\t{ih}\n")


def write_group_csvs(groups: dict, output_dir: str, group_id_map: Dict[Tuple[str, str], str]) -> int:
    paths.ensure_dir(output_dir)
    count = 0
    fieldnames = [
        "group_id",
        "identifier",
        "pdb_id",
        "chain_id",
        "InChI",
        "m_layer",
        "uniprot_id",
        "inchi_hash",
    ]
    for (uniprot_id, ih), members in sorted(groups.items()):
        gid = group_id_map[(uniprot_id, ih)]
        fname = f"{gid}_inchi_{ih}.csv"
        out_path = os.path.join(output_dir, fname)
        with open(out_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for m in sorted(members, key=lambda x: x["identifier"]):
                writer.writerow(
                    {
                        "group_id": gid,
                        "identifier": m["identifier"],
                        "pdb_id": m["pdb_id"],
                        "chain_id": m["chain_id"],
                        "InChI": m["inchi"],
                        "m_layer": m["m_layer"],
                        "uniprot_id": m["uniprot_id"],
                        "inchi_hash": m["inchi_hash"],
                    }
                )
        count += 1
    return count


def write_grouping_summary(
    paired: dict, single: dict, skipped: int, output_path: str
) -> None:
    paths.ensure_dir(os.path.dirname(output_path))
    paired_entries = sum(len(v) for v in paired.values())
    single_entries = sum(len(v) for v in single.values())
    lines = [
        "=== Enantiomer grouping summary ===",
        f"Paired groups: {len(paired)} ({paired_entries} entries)",
        f"Single groups: {len(single)} ({single_entries} entries)",
        f"Skipped (no UniProt): {skipped}",
        "",
        "Paired = same UniProt ID (full string) + same InChI before /m, with both m0 and m1.",
    ]
    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    for line in lines:
        print(line)


def find_prep_mae(pdb_id: str) -> Optional[str]:
    mae_path = os.path.join(paths.PREPWIZARD_OUTPUT_DIR, f"{pdb_id}_prepped.mae")
    if os.path.isfile(mae_path):
        return mae_path
    return None


def find_prep_pdb(pdb_id: str) -> Optional[str]:
    import glob

    pattern = os.path.join(paths.PREP_PDB_DIR, f"{pdb_id}_prepped*.pdb")
    matches = sorted(glob.glob(pattern))
    return matches[0] if matches else None


def convert_mae_to_pdb(pdb_ids: Set[str], schrodinger: Optional[str]) -> Tuple[int, int, int]:
    """Convert PrepWizard MAE to PDB for chain parsing. Returns converted, skipped, failed."""
    if not schrodinger:
        schrodinger = os.environ.get("SCHRODINGER")
    if not schrodinger:
        print("Warning: SCHRODINGER not set; skipping MAE to PDB conversion.")
        return 0, 0, len(pdb_ids)

    structconvert = os.path.join(schrodinger, "utilities", "structconvert")
    if not os.path.isfile(structconvert):
        print(f"Warning: structconvert not found at {structconvert}")
        return 0, 0, len(pdb_ids)

    paths.ensure_dir(paths.PREP_PDB_DIR)
    converted = skipped = failed = 0

    for pdb_id in sorted(pdb_ids):
        pdb_out = os.path.join(paths.PREP_PDB_DIR, f"{pdb_id}_prepped.pdb")
        if os.path.isfile(pdb_out) and os.path.getsize(pdb_out) > 0:
            skipped += 1
            continue
        mae_in = find_prep_mae(pdb_id)
        if not mae_in:
            print(f"Warning: MAE not found for {pdb_id}")
            failed += 1
            continue
        print(f"Converting {mae_in} -> {pdb_out}")
        result = subprocess.run(
            [structconvert, mae_in, pdb_out],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0 or not os.path.isfile(pdb_out):
            print(f"Failed to convert {pdb_id}: {result.stderr.strip()}")
            failed += 1
        else:
            converted += 1

    print(f"MAE->PDB: converted={converted}, skipped={skipped}, failed={failed}")
    return converted, skipped, failed


def parse_chains(paired: Dict, group_id_map: Dict[Tuple[str, str], str]) -> Set[Tuple]:
    """
    Extract chains from prepped PDB files for all paired enantiomer entries.
    Returns set of (uniprot_id, inchi_hash, pdb_id, chain_id, group_id) written.
    """
    paths.ensure_dir(paths.PARSED_CHAIN_DIR)
    parser = PDBParser(QUIET=True)
    io = PDBIO()

    # pdb_path -> set of (uniprot_id, pdb_id, inchi_hash, chain_id)
    pdb_requests = defaultdict(set)
    for members in paired.values():
        for m in members:
            pdb_requests[m["pdb_id"]].add(
                (m["uniprot_id"], m["pdb_id"], m["inchi_hash"], m["chain_id"])
            )

    list_entries: set[tuple[str, str, str, str, str]] = set()
    missing_chain = 0
    missing_pdb = 0

    for pdb_id, requests in sorted(pdb_requests.items()):
        pdb_path = find_prep_pdb(pdb_id)
        if not pdb_path:
            print(f"Warning: prepped PDB not found for {pdb_id}")
            missing_pdb += len(requests)
            continue
        try:
            structure = parser.get_structure(pdb_id, pdb_path)
            model = next(structure.get_models())
        except Exception as exc:
            print(f"Failed to parse {pdb_path}: {exc}")
            missing_pdb += len(requests)
            continue

        for uniprot_id, pdb_id_val, ih, chain_id in requests:
            gid = group_id_map[(uniprot_id, ih)]
            out_name = f"{gid}_{ih}_{pdb_id_val}_chain_{chain_id}.pdb"
            out_path = os.path.join(paths.PARSED_CHAIN_DIR, out_name)
            if not model.has_id(chain_id):
                print(f"Chain {chain_id} not found in {pdb_path}")
                missing_chain += 1
                continue
            if not os.path.isfile(out_path):
                chain = model[chain_id]
                io.set_structure(chain)
                io.save(out_path)
            list_entries.add((uniprot_id, ih, pdb_id_val, chain_id, gid))

    with open(paths.LIST_WITH_UNIPROT_ID, "w") as f:
        for uniprot_id, ih, pdb_id, chain_id, gid in sorted(list_entries):
            f.write(f"{uniprot_id}\t{ih}\t{pdb_id}\t{chain_id}\t{gid}\n")

    print(f"Parsed chain PDBs: {len(list_entries)}")
    print(f"Missing PDB: {missing_pdb}, missing chain: {missing_chain}")
    print(f"Wrote {paths.LIST_WITH_UNIPROT_ID}")
    return list_entries


def run(skip_mae_convert: bool = False):
    if not os.path.isfile(paths.FILTERED_ENANTIOMER_WITH_UNIPROT_TXT):
        print(f"Error: input not found: {paths.FILTERED_ENANTIOMER_WITH_UNIPROT_TXT}")
        sys.exit(1)

    print("=== Step 2.1: Group enantiomers ===")
    entries, skipped = load_entries(paths.FILTERED_ENANTIOMER_WITH_UNIPROT_TXT)
    paired, single = group_entries(entries)
    group_id_map = make_group_id_map(paired, single)

    paths.ensure_dir(paths.ALIGNMENT_DIR_DATA)
    write_group_manifest(group_id_map, paths.GROUP_MANIFEST)
    print(f"Wrote group manifest: {paths.GROUP_MANIFEST}")
    n_paired = write_group_csvs(paired, paths.PAIRED_ENANTIOMERS_DIR, group_id_map)
    n_single = write_group_csvs(single, paths.SINGLE_ENANTIOMERS_DIR, group_id_map)
    write_grouping_summary(paired, single, skipped, paths.GROUPING_SUMMARY)
    print(f"Wrote {n_paired} paired CSV(s), {n_single} single CSV(s)")

    pdb_ids = {m["pdb_id"] for members in paired.values() for m in members}

    if not skip_mae_convert:
        print("Converting PrepWizard MAE to PDB")
        convert_mae_to_pdb(pdb_ids, os.environ.get("SCHRODINGER"))
    else:
        print("Skipped MAE to PDB (--skip-mae-convert)")

    print("\n=== Step 2.2: Parse chains from prepped PDB ===")
    if not paired:
        print("No paired enantiomer groups; nothing to parse.")
        return
    parse_chains(paired, group_id_map)


def main():
    parser = argparse.ArgumentParser(
        description="Alignment steps 2.1-2.2: group enantiomers, convert MAE to PDB, parse chains"
    )
    parser.add_argument(
        "--skip-mae-convert",
        action="store_true",
        help="Skip Schrödinger MAE to PDB conversion (PDB must already exist)",
    )
    args = parser.parse_args()
    run(skip_mae_convert=args.skip_mae_convert)


if __name__ == "__main__":
    main()
