from rdkit import Chem
from rdkit.Chem import SDWriter
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
import requests
from tqdm import tqdm


def filter_molecules(
    mols,
    filtered_sdf_path,
    filtered_out_sdf_path,
    mw_min=181,
    mw_max=800,
    h_acceptors_max=16,
    rot_bonds_max=20,
    small_rings_min=1,
):
    filtered_mols = []
    filtered_out_mols = []
    failed_to_load_mols = 0
    total_mols_read = 0
    for mol in mols:
        if mol is None:
            failed_to_load_mols += 1
            continue
        total_mols_read += 1
        mol_no_h = Chem.RemoveHs(mol)
        mw = float(Descriptors.MolWt(mol_no_h))
        h_acceptors = int(Descriptors.NumHAcceptors(mol))
        rot_bonds = int(rdMolDescriptors.CalcNumRotatableBonds(
            mol_no_h,
            strict=rdMolDescriptors.NumRotatableBondsOptions.StrictLinkages
        ))
        ring_count = int(Descriptors.RingCount(mol_no_h))
        mol_no_h.SetProp("MW", str(mw))
        mol_no_h.SetProp("NumHAcceptors", str(h_acceptors))
        mol_no_h.SetProp("NumRotatableBonds", str(rot_bonds))
        mol_no_h.SetProp("RingCount", str(ring_count))

        passes_criteria = (mw > mw_min and mw < mw_max and
                        h_acceptors <= h_acceptors_max and
                        rot_bonds <= rot_bonds_max and
                        ring_count >= small_rings_min)
        if passes_criteria:
            filtered_mols.append(mol_no_h)
        else:
            filtered_out_mols.append(mol_no_h)
    print(f"Number of molecules processed: {total_mols_read}")
    print(f"Number of molecules failed to load: {failed_to_load_mols}")
    print(f"Number of molecules passed all filters: {len(filtered_mols)}")
    print(f"Number of molecules filtered out: {len(filtered_out_mols)}")
    if filtered_mols:
        with SDWriter(filtered_sdf_path) as writer:
            for mol in filtered_mols:
                writer.write(mol)
        print(f"Saved filtered molecules to {filtered_sdf_path}")
    if filtered_out_mols:
        failed_counts = {
            "MW": 0,
            "NumHAcceptors": 0,
            "NumRotatableBonds": 0,
            "RingCount": 0,
        }
        for mol in filtered_out_mols:
            mw = float(mol.GetProp("MW"))
            h_acc = int(mol.GetProp("NumHAcceptors"))
            rot_b = int(mol.GetProp("NumRotatableBonds"))
            rings = int(mol.GetProp("RingCount"))
            if mw <= mw_min:
                failed_counts["MW"] += 1
            if h_acc > h_acceptors_max:
                failed_counts["NumHAcceptors"] += 1
            if rot_b > rot_bonds_max:
                failed_counts["NumRotatableBonds"] += 1
            if rings < small_rings_min:
                failed_counts["RingCount"] += 1
        print("Number of filtered out molecules failing each condition:")
        print(f"  MW <= {mw_min}: {failed_counts['MW']}")
        print(f"  NumHAcceptors > {h_acceptors_max}: {failed_counts['NumHAcceptors']}")
        print(f"  NumRotatableBonds > {rot_bonds_max}: {failed_counts['NumRotatableBonds']}")
        print(f"  RingCount < {small_rings_min}: {failed_counts['RingCount']}")

        with SDWriter(filtered_out_sdf_path) as writer:
            for mol in filtered_out_mols:
                writer.write(mol)
        print(f"Saved filtered out molecules to {filtered_out_sdf_path}")

    return filtered_mols, filtered_out_mols


def get_uniprotID(identifiers):
    results = {}
    session = requests.Session()
    base_entry_url = "https://data.rcsb.org/rest/v1/core/entry/"
    base_polymer_url = "https://data.rcsb.org/rest/v1/core/polymer_entity/"
    for identifier in tqdm(identifiers):
        pdb_id = identifier.split("_")[0]
        chain_id = identifier.split("_")[2]
        all_uniprot_ids = set()
        all_chain_ids = set()
        try:
            response = session.get(base_entry_url + pdb_id)
            response.raise_for_status()
            entry_data = response.json()
            polymer_entities = entry_data.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
            for entity_id in polymer_entities:
                try:
                    poly_resp = session.get(f"{base_polymer_url}{pdb_id}/{entity_id}")
                    poly_resp.raise_for_status()
                    entity_data = poly_resp.json()
                    chain_ids = entity_data.get("rcsb_polymer_entity_container_identifiers", {}).get("asym_ids", [])
                    for all_chains in chain_ids:
                        if chain_id in all_chains:
                            uniprot_ids = entity_data.get("rcsb_polymer_entity_container_identifiers", {}).get("uniprot_ids", [])
                            all_uniprot_ids.update(uniprot_ids)
                            all_chain_ids.update(chain_id)
                except requests.exceptions.RequestException:
                    continue
        except requests.exceptions.RequestException:
            results[identifier] = {"uniprot_ids": [], "chain_id": []}
            continue
        results[identifier] = {"pdb_id": pdb_id, "uniprot_ids": list(all_uniprot_ids), "chain_id": list(all_chain_ids)}
    session.close()
    return results
