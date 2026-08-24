# 1. Data Preprocess

Extract ligands from PDB mmCIF files, find enantiomer pairs, keep drug-like ligands, add UniProt IDs, and prepare proteins with Schrödinger.

Scripts are in this folder. Outputs go under `../Data/`. Paths: [`paths.py`](paths.py).

Place gzipped mmCIF input in **`../Data/mmCIF/`** (`*.cif.gz`). The pipeline writes decompressed files to `../Data/mmCIF_rename/`. Test data is under `../Data/test/mmCIF/` (outputs under `../Data/test/`).

## Run

From the repository root:

```bash
export SCHRODINGER=/path/to/your/schrodinger
# Optional: only if your archive is not Data/mmCIF
export MMCIF_SOURCE=/path/to/your/mmCIF          # .cif.gz; default Data/mmCIF
export MMCIF_RENAME=/path/to/your/mmCIF_rename   # decompressed .cif; default Data/mmCIF_rename
bash run_data_preprocess.sh
```

Offline: add `--skip-uniprot --skip-prepwizard`

Needs RDKit, `requests`, `tqdm`, [UniCON](https://www.zbh.uni-hamburg.de/forschung/amd/software/unicon.html) (obtain from the [NAOMI ChemBio Suite](https://software.zbh.uni-hamburg.de) and place the tree at `Data_preprocess/unicon_1.5.0/` so the binary is `unicon_1.5.0/unicon`), a licensed Schrödinger install (PrepWizard), and RCSB access for UniProt. UniCON is not included in this repository.

## Steps

| Step | Script | Output |
|------|--------|--------|
| 1.1a | `extract_split_ligands.sh` | One SDF per ligand |
| 1.1b | `convert_sdf_to_inchi.py` | `all_inchis_dict.txt` |
| 1.1c | `classify_chiral_compounds.py` | `enantiomers.txt` + `sdf/enantiomer_structure.sdf` |
| 1.2–1.3 | `filter_and_annotate_enantiomers.py` | Filtered SDF, UniProt map, `db_ids_with_m_layer.txt` |
| 1.4 | `maestro_protein_preparation.sh` | `{pdb}_prepped.mae` |

Drug-likeness (1.2): MW 181–800, HBA ≤ 16, rotatable bonds ≤ 20, ≥ 1 ring.

Identifiers look like `{pdb}_{het}_{chain}_{resseq}` (and `_m0` / `_m1` in `db_ids_with_m_layer.txt`). Enantiomers share the InChI `/t` layer and differ in `/m`.

Existing files are skipped on re-run. Delete the relevant `Data/` folders for a full redo.

## Outputs (`Data/`)

```
ligand_sdf_structure_from_mmCIF(_separated)/   # 1.1a
ligand_sdf_..._InChI/                          # 1.1b
chiral_classification_results/                 # 1.1c–1.3
db_ids_with_m_layer.txt                        # 1.3
Enantiomer_pdbstructure_maeformat/             # 1.4
Enantiomer_protein_preperation/                # 1.4 PrepWizard MAE
```

Next: [Alignment](../Alignment/README.md)
