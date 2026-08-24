# 2. Alignment

Group paired enantiomers on the same UniProt ID, align those chains, and extract ligands.

Only groups with both `m0` and `m1` are aligned. Single-enantiomer groups are listed but not processed further.

Scripts are in this folder. Outputs go under `../Data/Alignment/`. Paths: [`paths.py`](paths.py).

Needs [Data Preprocess](../Data_preprocess/README.md) through step 1.4.

## Run

From the repository root:

```bash
export SCHRODINGER=/path/to/your/schrodinger
bash run_alignment.sh
```

Skip Schrödinger: `--skip-mae-convert --skip-structalign`  
Skip ligand extract: `--skip-ligand-extract`

Needs BioPython, a licensed Schrödinger install (`structconvert`, `structalign`), and [LigandExtractor](https://www.zbh.uni-hamburg.de/forschung/amd/software/ligandextractor.html) (obtain from the [NAOMI ChemBio Suite](https://software.zbh.uni-hamburg.de) and place the tree at `Alignment/LigandExtractor_1.0.1/` so the binary is `LigandExtractor_1.0.1/LigandExtractor`). LigandExtractor is not included in this repository.

## Steps

| Step | Script | Output |
|------|--------|--------|
| 2.1–2.2 | `group_and_parse_chains.py` | Group CSVs, MAE→PDB, chain PDBs |
| 2.3 | `schrodinger.sh` | Aligned PDBs |
| 2.4 | `ligand_extractor_prep.sh` | Ligand SDFs |

Groups are `(UniProt, InChI-before-/m)` hashed to 10 hex chars, with a stable `group_id` (`g000001`, …). Filenames use `group_id`. Mapping: `group_manifest.tsv`.

Existing prep PDBs, parsed chains, and ligand SDFs are skipped. Delete `Data/Alignment/` for a full redo. Step 2.3 re-runs `structalign`.

## Outputs (`Data/Alignment/`)

```
paired_enantiomers/                          # 2.1 CSVs
single_enantiomers/                          # 2.1 audit only
group_manifest.tsv
Enantiomer_prep_pdbformat/                   # 2.1
paired_enantiomers_parsed_chain/             # 2.2
paired_enantiomers_aligned_structures/       # 2.3
Enantiomer_aligned_structure_ligandextract/  # 2.4
```

Next: [Pharmacophore generation](../Pharmacophore_generation_after_alignment/README.md)
