import json
import re
import os
import time
from collections import defaultdict
from datetime import timedelta
from rdkit import Chem

import paths


# ============================================================================
# Configuration
# ============================================================================

inchi_path = paths.ALL_INCHIS_DICT
output_dir = paths.CHIRAL_RESULTS_DIR
sdf_path = paths.MERGED_LIGANDS_SDF
sdf_output_dir = paths.CHIRAL_SDF_DIR

# ============================================================================
# SDF Utilities
# ============================================================================

def load_mol_lookup(sdf_file):
    if not os.path.exists(sdf_file):
        print(f"SDF file not found: {sdf_file}")
        return {}
    suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
    lookup = {}
    for mol in suppl:
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else None
        if name:
            lookup[name] = mol
    print(f"Loaded {len(lookup)} molecules into lookup from {sdf_file}")
    return lookup

# ============================================================================
# Write SDF helpers
# ============================================================================

def collect_entries_from_pairs(pairs):
    """Flatten pairs into unique (key, inchi) items, preferring first occurrence."""
    entries = {}
    for key1, inchi1, key2, inchi2 in pairs:
        if key1 not in entries:
            entries[key1] = inchi1
        if key2 not in entries:
            entries[key2] = inchi2
    return list(entries.items())


def write_structures(entries, mol_lookup, out_path, key_props=None, label="structures", total_count=None):
    writer = Chem.SDWriter(out_path)
    written = 0
    missing = 0
    for key, inchi_val in entries:
        mol = mol_lookup.get(key)
        if mol is None:
            missing += 1
            continue
        if inchi_val:
            try:
                mol.SetProp("InChI", inchi_val)
            except Exception:
                pass
        if key_props and key in key_props:
            for prop, val in key_props[key].items():
                try:
                    mol.SetProp(prop, str(val))
                except Exception:
                    pass
        writer.write(mol)
        written += 1
    writer.close()
    if total_count is not None and total_count > 0:
        percentage = written / total_count * 100
        print(f"{label} written: {written} ({percentage:.2f}%) (missing: {missing}) -> {out_path}")
    else:
        print(f"{label} written: {written} (missing: {missing}) -> {out_path}")
    return out_path

# ============================================================================
# Utility Functions for InChI Parsing
# ============================================================================

def extract_t_layer(inchi):
    """Extract the /t layer content from an InChI string."""
    t_match = re.search(r'/t([^/]+)', inchi)
    return t_match.group(1) if t_match else None


def extract_m_layer(inchi):
    """Extract the /m layer content from an InChI string."""
    m_match = re.search(r'/m([^/]+)', inchi)
    return m_match.group(1) if m_match else None


def get_base_structure(inchi):
    """Extract the base structure (everything before /t layer) from an InChI string."""
    if '/t' in inchi:
        return inchi.split('/t')[0]
    else:
        return inchi


# ============================================================================
# Filtering Functions
# ============================================================================

def filter_entries(inchi_dict):
    no_stereo_count = 0
    undefined_stereo_entries = []
    remaining_entries = {}
    
    print("\n=== Step 1: Filtering entries ===")
    for key, inchi in inchi_dict.items():
        if '/t' not in inchi:
            no_stereo_count += 1
            continue
        # Check if /t layer contains "?"
        if '/t' in inchi:
            t_match = re.search(r'/t([^/]+)', inchi)
            if t_match and '?' in t_match.group(1):
                undefined_stereo_entries.append((key, inchi))
                continue
        # Keep for further processing
        remaining_entries[key] = inchi
    
    print(f"Entries without stereochemistry: {no_stereo_count}")
    print(f"Entries with undefined stereo: {len(undefined_stereo_entries)}")
    print(f"Entries remaining for grouping: {len(remaining_entries)}")
    
    return no_stereo_count, undefined_stereo_entries, remaining_entries


# ============================================================================
# Grouping Functions
# ============================================================================

def group_by_structure(remaining_entries):
    print("\n=== Step 2: Grouping by structure ===")
    structure_groups = defaultdict(list)
    
    for key, inchi in remaining_entries.items():
        base_structure = get_base_structure(inchi)
        structure_groups[base_structure].append((key, inchi))
    
    print(f"Number of structure groups: {len(structure_groups)}")
    return structure_groups


def find_enantiomer_pairs(group):
    enantiomer_pairs = []
    # Compare all pairs in the group to find all enantiomer relationships
    for i, (key1, inchi1) in enumerate(group):
        t1 = extract_t_layer(inchi1)
        m1 = extract_m_layer(inchi1)
        for j, (key2, inchi2) in enumerate(group[i+1:], start=i+1):
            t2 = extract_t_layer(inchi2)
            m2 = extract_m_layer(inchi2)
            # Check for enantiomers: same including /t, different /m
            if t1 is not None and t2 is not None and t1 == t2:
                if m1 != m2:  # Different /m layers (including one None)
                    enantiomer_pairs.append((key1, inchi1, key2, inchi2))
    
    return enantiomer_pairs

def find_diastereomer_pairs(group):
    diastereomer_pairs = []
    for i, (key1, inchi1) in enumerate(group):
        t1 = extract_t_layer(inchi1)
        # if t1 is None:
        #     continue
        for j, (key2, inchi2) in enumerate(group[i+1:], start=i+1):
            t2 = extract_t_layer(inchi2)
            # if t2 is None:
            #     continue
            # Check for diastereomers: different /t content
            if t1 != t2:
                diastereomer_pairs.append((key1, inchi1, key2, inchi2))
    return diastereomer_pairs


def classify_within_groups(structure_groups):
    print("\n=== Step 3: Classifying within groups ===")
    enantiomer_pairs = []
    diastereomer_pairs = []
    identical_entries = []
    ungrouped_count = 0
    
    for base_structure, group in structure_groups.items():
        if len(group) == 1:
            # Single entry in group - treat as identical
            key, inchi = group[0]
            identical_entries.append((key, inchi))
            continue
        # Track which entries have been paired
        paired_keys = set()
        # Group by full InChI to find identicals
        inchi_groups = defaultdict(list)
        for key, inchi in group:
            inchi_groups[inchi].append(key)
        distinct_inchis = len(inchi_groups)
        # If only one distinct InChI, mark all as identical and skip stereo classification
        if distinct_inchis == 1:
            for inchi_val, keys in inchi_groups.items():
                for key in keys:
                    identical_entries.append((key, inchi_val))
                    paired_keys.add(key)
            continue
        # If multiple distinct InChIs, skip recording identical duplicates; classify using full group
        classification_group = group
        group_enantiomer_pairs = find_enantiomer_pairs(classification_group)
        for key1, inchi1, key2, inchi2 in group_enantiomer_pairs:
            enantiomer_pairs.append((key1, inchi1, key2, inchi2))
            paired_keys.add(key1)
            paired_keys.add(key2)
        
        # Find all diastereomer pairs (can have duplicates with enantiomer pairs)
        group_diastereomer_pairs = find_diastereomer_pairs(classification_group)
        for key1, inchi1, key2, inchi2 in group_diastereomer_pairs:
            diastereomer_pairs.append((key1, inchi1, key2, inchi2))
            paired_keys.add(key1)
            paired_keys.add(key2)
        
        # All remaining entries go to ungrouped
        for key, inchi in group:
            if key not in paired_keys:
                ungrouped_count += 1
    
    print(f"Enantiomer pairs found: {len(enantiomer_pairs)}")
    print(f"Diastereomer pairs found: {len(diastereomer_pairs)}")
    print(f"Identical entries: {len(identical_entries)}")
    print(f"Ungrouped entries: {ungrouped_count}")
    
    return enantiomer_pairs, diastereomer_pairs, identical_entries


# ============================================================================
# Enantiomer Grouping Functions
# ============================================================================
def group_enantiomers(enantiomer_pairs, enantiomer_set):
    print("\n=== Grouping enantiomers ===")
    enantiomer_entry_map = {}
    unique_enantiomer_pairs = []
    
    # Build a map of all unique enantiomer pairs (using deduplicated set)
    # This ensures we only process pairs that were actually written to file
    for key1, inchi1, key2, inchi2 in enantiomer_pairs:
        pair_key = tuple(sorted([key1, key2]))
        if pair_key in enantiomer_set:
            unique_enantiomer_pairs.append((key1, key2))
            enantiomer_entry_map[key1] = (key1, inchi1)
            enantiomer_entry_map[key2] = (key2, inchi2)
    
    # Use union-find (disjoint set union) to group connected enantiomers
    # This algorithm handles transitive relationships automatically:
    # - If A-B is a pair, A and B are in the same group
    # - If B-C is a pair, B and C are in the same group
    # - Therefore, A, B, C are all in the same group
    parent = {}
    
    def find(x):
        """
        Find the root of x with path compression.
        
        Implemented iteratively with cycle detection to avoid runaway recursion
        if upstream data creates cycles (e.g., corrupted parent pointers or
        non-reflexive keys).
        """
        if x not in parent:
            parent[x] = x
        path = []
        cur = x
        seen = set()
        while True:
            if cur not in parent:
                parent[cur] = cur
            nxt = parent[cur]
            # Normal root found
            if nxt == cur:
                root = cur
                break
            seen.add(cur)
            path.append(cur)
            cur = nxt
        # Path compression
        for node in path:
            parent[node] = root
        return root
    
    def union(x, y):
        """Union two elements into the same group."""
        root_x = find(x)
        root_y = find(y)
        if root_x != root_y:
            parent[root_y] = root_x
    
    # Union all pairs - this creates groups for all connected enantiomers
    for key1, key2 in unique_enantiomer_pairs:
        union(key1, key2)
    
    # Group all keys by their root
    enantiomer_groups = defaultdict(set)
    for key in parent.keys():
        root = find(key)  # Ensure path compression is applied
        enantiomer_groups[root].add(key)
    
    print(f"Number of enantiomer groups: {len(enantiomer_groups)}")
    unique_enantiomer_structures = len({k for group in enantiomer_groups.values() for k in group})
    print(f"Unique enantiomer structures: {unique_enantiomer_structures}")
    return enantiomer_groups, enantiomer_entry_map


# ============================================================================
# File Writing Functions
# ============================================================================

def write_undefined_stereo(undefined_stereo_entries, output_dir):
    """Write undefined stereo entries to file."""
    undefined_stereo_path = os.path.join(output_dir, 'undefined_stereo.txt')
    with open(undefined_stereo_path, 'w') as f:
        for key, inchi in undefined_stereo_entries:
            f.write(f"{key}\t{inchi}\n")


def find_structure_matches_for_undefined(undefined_stereo_entries, inchi_dict):
    # Build base-structure index over all entries
    base_index = defaultdict(list)
    for key, inchi in inchi_dict.items():
        base = get_base_structure(inchi)
        base_index[base].append((key, inchi))

    matches = []
    seen = set()
    for key_u, inchi_u in undefined_stereo_entries:
        base = get_base_structure(inchi_u)
        for key_o, inchi_o in base_index.get(base, []):
            if key_o == key_u:
                continue
            pair_key = tuple(sorted([key_u, key_o]))
            if pair_key in seen:
                continue
            seen.add(pair_key)
            matches.append((key_u, inchi_u, key_o, inchi_o))
    return matches


def write_structure_matches(matches, output_dir, filename):
    """Write stereo-agnostic structure matches to a text file."""
    out_path = os.path.join(output_dir, filename)
    with open(out_path, 'w') as f:
        for key1, inchi1, key2, inchi2 in matches:
            f.write(f"{key1}\t{inchi1}\n")
            f.write(f"{key2}\t{inchi2}\n")
            f.write("\n")
    print(f"Undefined stereo structure matches (ignoring stereo) written: {len(matches)} -> {out_path}")
    return out_path


def write_enantiomer_pairs(enantiomer_pairs, output_dir):
    """Write enantiomer pairs to file (deduplicated)."""
    enantiomer_set = set()
    enantiomers_path = os.path.join(output_dir, 'enantiomers.txt')
    with open(enantiomers_path, 'w') as f:
        for key1, inchi1, key2, inchi2 in enantiomer_pairs:
            pair_key = tuple(sorted([key1, key2]))
            if pair_key not in enantiomer_set:
                enantiomer_set.add(pair_key)
                f.write(f"{key1}\t{inchi1}\n")
                f.write(f"{key2}\t{inchi2}\n")
                f.write("\n")
    print(f"Unique enantiomer pairs written: {len(enantiomer_set)}")
    return enantiomer_set

def write_enantiomer_representatives(enantiomer_groups, enantiomer_entry_map, 
                                    remaining_entries, mol_lookup, output_dir, sdf_output_dir):
    """Write enantiomer representatives to text and SDF (with group size)."""
    enantiomer_rep_path = os.path.join(output_dir, 'enantiomer_representative.txt')
    enantiomer_rep_sdf = os.path.join(sdf_output_dir, 'enantiomer_representative.sdf')
    
    rep_entries = []
    rep_props = {}
    with open(enantiomer_rep_path, 'w') as f:
        for group_id, group_keys in sorted(enantiomer_groups.items()):
            rep_key = sorted(group_keys)[0]
            rep_entry = enantiomer_entry_map.get(rep_key)
            if rep_entry is None:
                rep_inchi = remaining_entries.get(rep_key, "")
            else:
                _, rep_inchi = rep_entry
            
            group_size = len(group_keys)
            f.write(f"Group representative: {rep_key}\n")
            f.write(f"Number of entries in group: {group_size}\n")
            f.write(f"InChI: {rep_inchi}\n")
            f.write("\n")

            rep_entries.append((rep_key, rep_inchi))
            rep_props[rep_key] = {"group_size": group_size}

    write_structures(rep_entries, mol_lookup, enantiomer_rep_sdf, key_props=rep_props, label="Enantiomer representatives", total_count=len(remaining_entries))
    return enantiomer_rep_path, enantiomer_rep_sdf

def write_diastereomer_pairs(diastereomer_pairs, output_dir, filename='diastereomers.txt'):
    diastereomer_set = set()
    diastereomers_path = os.path.join(output_dir, filename)
    
    with open(diastereomers_path, 'w') as f:
        for key1, inchi1, key2, inchi2 in diastereomer_pairs:
            pair_key = tuple(sorted([key1, key2]))
            if pair_key not in diastereomer_set:
                diastereomer_set.add(pair_key)
                f.write(f"{key1}\t{inchi1}\n")
                f.write(f"{key2}\t{inchi2}\n")
                f.write("\n")
    
    return diastereomer_set

def write_identical_entries(identical_entries, output_dir):
    identical_path = os.path.join(output_dir, 'identical.txt')
    with open(identical_path, 'w') as f:
        for key, inchi in identical_entries:
            f.write(f"{key}\t{inchi}\n")
    
def write_summary(inchi_dict, no_stereo_count, undefined_stereo_entries, 
                remaining_entries, execution_time, output_dir):
    """Write summary statistics to file and print to console."""
    print("\n=== Summary ===")
    
    readable_time = str(timedelta(seconds=int(execution_time)))
    
    summary_lines = [
        "=== Chiral Compound Classification Summary ===",
        "",
        f"Total entries: {len(inchi_dict)}",
        f"Entries without stereochemistry: {no_stereo_count}, {no_stereo_count / len(inchi_dict) * 100:.2f}%",
        f"Entries with undefined stereo: {len(undefined_stereo_entries)}, {len(undefined_stereo_entries) / len(inchi_dict) * 100:.2f}%",
        f"Entries processed for grouping: {len(remaining_entries)}, {len(remaining_entries) / len(inchi_dict) * 100:.2f}%",
        "",
        f"Total script execution time: {readable_time}",
    ]
    
    for line in summary_lines:
        print(line)
    
    # Write summary to file
    summary_path = os.path.join(output_dir, 'summary.txt')
    with open(summary_path, 'w') as f:
        f.write("\n".join(summary_lines))
        f.write("\n")
    
    print(f"\nSummary written to: {summary_path}")


# ============================================================================
# Main Function
# ============================================================================

def main():
    start_time = time.time()
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(sdf_output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")
    print(f"SDF output directory: {sdf_output_dir}")

    # Load SDF molecules for structure export
    mol_lookup = load_mol_lookup(sdf_path)
    
    # Load the InChI dictionary
    print("Loading InChI dictionary...")
    with open(inchi_path, 'r') as f:
        inchi_dict = json.load(f)
    
    print(f"Total entries: {len(inchi_dict)}")
    
    # Step 1: Filter entries
    no_stereo_count, undefined_stereo_entries, remaining_entries = filter_entries(inchi_dict)
    
    # Write undefined stereo entries
    write_undefined_stereo(undefined_stereo_entries, output_dir)
    if undefined_stereo_entries:
        write_structures(undefined_stereo_entries, mol_lookup, os.path.join(sdf_output_dir, "undefined_stereo_structure.sdf"), label="Undefined stereo structures")
        # Find matches ignoring stereo for undefined entries
        undefined_matches = find_structure_matches_for_undefined(undefined_stereo_entries, inchi_dict)
        write_structure_matches(undefined_matches, output_dir, 'undefined_stereo_structure_matches.txt')
        if undefined_matches:
            undefined_match_entries = collect_entries_from_pairs(undefined_matches)
            write_structures(
                undefined_match_entries,
                mol_lookup,
                os.path.join(sdf_output_dir, "undefined_stereo_structure_matches.sdf"),
                label="Undefined stereo matches (ignoring stereo)",
                total_count=len(remaining_entries),
            )
    
    # Step 2: Group by structure
    structure_groups = group_by_structure(remaining_entries)
    
    # Step 3: Classify within groups
    enantiomer_pairs, diastereomer_pairs, identical_entries = \
        classify_within_groups(structure_groups)
    
    # Step 4: Write results
    print("\n=== Step 4: Writing results ===")
    
    # Write identical entries
    write_identical_entries(identical_entries, output_dir)
    write_structures(
        identical_entries,
        mol_lookup,
        os.path.join(sdf_output_dir, "identical_structure.sdf"),
        label="Identical structures",
        total_count=len(remaining_entries),
    )
    
    # Write enantiomer pairs
    enantiomer_set = write_enantiomer_pairs(enantiomer_pairs, output_dir)
    
    # Write enantiomer structures
    enantiomer_entries = collect_entries_from_pairs(enantiomer_pairs)
    write_structures(enantiomer_entries, mol_lookup, os.path.join(sdf_output_dir, "enantiomer_structure.sdf"), label="Enantiomer structures", total_count=len(remaining_entries))
    
    # Group enantiomers and write representatives
    enantiomer_groups, enantiomer_entry_map = group_enantiomers(enantiomer_pairs, enantiomer_set)

    write_enantiomer_representatives(
        enantiomer_groups, enantiomer_entry_map, remaining_entries, mol_lookup, output_dir, sdf_output_dir
    )
    
    # Write diastereomer pairs
    write_diastereomer_pairs(diastereomer_pairs, output_dir)
    diastereomer_entries = collect_entries_from_pairs(diastereomer_pairs)
    write_structures(diastereomer_entries, mol_lookup, os.path.join(
        sdf_output_dir, 
        "diastereomer_structure.sdf"), 
        label="Diastereomer structures", 
        total_count=len(remaining_entries)
    )
    
    # Calculate execution time
    end_time = time.time()
    execution_time = end_time - start_time
    
    write_summary(
        inchi_dict, no_stereo_count, undefined_stereo_entries, remaining_entries,
        execution_time, output_dir
    )


if __name__ == "__main__":
    main()
