# 5. Maestro property

Convert aligned receptor PDBs to Maestro `.mae` format and assign user properties for downstream visualization and analysis.

**Prerequisite:** aligned structures from [Alignment](../Alignment/README.md) (step 2.3), plus the db-id list and UniProt classification file from data preprocessing.

---

## Workflow

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 5 | [`editing_meastro_property_new.py`](editing_meastro_property_new.py) | Aligned PDBs + db ids + UniProt file | MAE files with user properties |

This script does not use `run_maestro_property.sh`. That shell script only launches the older pocket-CSV editor, [`editing_meastro_property.py`](editing_meastro_property.py).

### Steps performed

1. **PDB → MAE** — convert each file in `Data/Alignment/paired_enantiomers_aligned_structures/` using Schrödinger `structconvert`. Existing outputs in `converted_mae/` are reused and not re-converted (use `--overwrite` to force reconversion).
2. **Property editing** — for each aligned structure (e.g. `rot-g000001_b55a025257_4tyb_chain_D`), assign Maestro user properties from the db-id list and the UniProt file.

### Maestro user properties

| Property | Source | Example |
|----------|--------|---------|
| `s_user_ID` | `{uniprot_field}_{pdb}_{ligand}_{chain}_{res}_mN`. The UniProt field is the whole second-to-last column, commas included. Several matches for the same PDB id and chain are stored as a list. | `['A0A7S3VNZ8,A0A7S3VSS7_9mh1_LUT_2_314_m1', ...]` |
| `s_user_group` | Group token + InChI hash from the aligned filename | `{'g000001_b55a025257'}` |

A leading `rot-` is stripped and a trailing `_copied` on a per-chain file is ignored, so `rot-g000001_b55a025257_4tyb_chain_D.pdb` and `g000001_b55a025257_4tyb_chain_D_copied.pdb` both get `{'g000001_b55a025257'}`.

`s_user_ID` is attached by PDB id and chain. The m-layer suffix (`_m0`, `_m1`, …) is removed only to find the UniProt row; the value written keeps the full db id.

Structures whose names cannot be parsed, or for which no properties can be assigned, are **not** written to `edited_mae_new/` or the combined MAE. They are listed in `unassigned_structures_new.log`.

**Reference / copied structures without group tokens**

For structures like `rot-1ks7_prepped.pdb` or `9fa6_prepped_copied.pdb`:

1. Look for matching sibling files in the aligned folder with the same PDB id. One output MAE is written. `rot-<pdb>_prepped` uses per-chain siblings that are not `_copied`. `<pdb>_prepped_copied` uses per-chain `*_copied` siblings.
2. If those siblings belong to **different groups**, `s_user_ID` and `s_user_group` are aggregated across all siblings with properties in each group. If siblings share one group (e.g. multiple chains), properties from the first matching sibling are used.
3. If no sibling exists, the structure is logged to `unassigned_structures_new.log` and skipped.

### Inputs

| Path | Source |
|------|--------|
| `Data/Alignment/paired_enantiomers_aligned_structures/` | Alignment step 2.3 |
| `Data/db_ids_with_m_layer.txt` | One db id per line, e.g. `9mh1_LUT_2_314_m1` |
| `Data/chiral_classification_results/filtered_enantiomer_with_uniprot_new.txt` | Instance id without the m-layer in column 1; UniProt is the second-to-last field |

Both text files are resolved under the pipeline `Data/` directory (`test/` when `PIPELINE_TEST=1`). Override them with `--db-ids` and `--uniprot-file`.

### Outputs

```
Data/Maestro_property/
├── converted_mae/                              # PDB -> MAE conversion
├── edited_mae_new/                             # per-structure MAE with properties
├── unassigned_structures_new.log               # structures skipped (no properties)
└── all_structures_with_properties_new.mae      # combined output (assigned only)
```

---

## Production run

Requires Schrödinger (`$SCHRODINGER`):

```bash
export SCHRODINGER=/path/to/schrodinger
$SCHRODINGER/run python3 Maestro_property/editing_meastro_property_new.py
$SCHRODINGER/run python3 Maestro_property/editing_meastro_property_new.py --skip-convert
$SCHRODINGER/run python3 Maestro_property/editing_meastro_property_new.py --overwrite
```

`--skip-convert` reuses MAE files already in `converted_mae/`. `--overwrite` re-runs `structconvert` even when those files exist.

---

## Previous pocket-CSV workflow

[`editing_meastro_property.py`](editing_meastro_property.py) still assigns `s_user_ID` and `s_user_group` from pharmacophore interaction CSVs. [`run_maestro_property.sh`](../run_maestro_property.sh) launches only that script and checks for the pocket-based CSVs and the categorization file. Its outputs are `edited_mae/`, `unassigned_structures.log`, and `all_structures_with_properties.mae`.

```bash
bash run_maestro_property.sh
bash run_maestro_property.sh --clean
bash run_maestro_property.sh --skip-convert
bash run_maestro_property.sh --overwrite
```

---

## Paths

[`paths.py`](paths.py) — respects `PIPELINE_TEST=1` for `test/` outputs.

---

## Related documentation

- [Alignment](../Alignment/README.md)
