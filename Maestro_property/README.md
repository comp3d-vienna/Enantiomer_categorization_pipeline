# 5. Maestro property

Convert aligned receptor PDBs to Maestro `.mae` format and assign user properties for downstream visualization and analysis.

**Prerequisite:** complete [Pharmacophore generation](../Pharmacophore_generation_after_alignment/README.md) (step 3).

---

## Workflow

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 5 | [`editing_meastro_property.py`](editing_meastro_property.py) | Aligned PDBs + interaction CSVs | MAE files with user properties |

### Steps performed

1. **PDB → MAE** — convert each file in `Data/Alignment/paired_enantiomers_aligned_structures/` using Schrödinger `structconvert`. Existing outputs in `converted_mae/` are reused and not re-converted (use `--overwrite` to force reconversion).
2. **Property editing** — for each aligned structure (e.g. `rot-g000001_b55a025257_4tyb_chain_D`), assign Maestro user properties from pharmacophore interaction CSVs.

### Maestro user properties

| Property | Source | Example |
|----------|--------|---------|
| `s_user_ID` | Interaction TSV instance IDs, prefixed by UniProt ID (from `group_manifest.tsv`) | `['A0A009KIJ4_4tyb_3B1_A_601_m0', ...]` |
| `s_user_group` | `group_id` + InChI hash | `{'g000001_b55a025257'}` |

Property assignment follows the reference `editing_meastro_property.py` logic: pocket CSVs are matched by `group_id` and `inchi_hash` in the filename, interaction rows are matched by PDB ID and chain in the `file` column, and `s_user_enantiomer_short_hash` is omitted.

Structures whose names cannot be parsed, or for which no properties can be assigned, are **not** written to `edited_mae/` or the combined MAE. They are listed in `unassigned_structures.log`.

**Reference / copied structures without group tokens**

For structures like `rot-1ks7_prepped.pdb` or `9fa6_prepped_copied.pdb`:

1. **Aligned folder first** — look for matching sibling files with the same PDB ID. One output MAE is written. If siblings belong to **different groups**, `s_user_ID` and `s_user_group` are aggregated across all siblings with properties in each group. If siblings share one group (e.g. multiple chains), properties from the first matching sibling are used.
2. **Pocket CSV fallback** — if no sibling exists, resolve a single identity from pharmacophore pocket-based CSVs (`file` column → PDB ID + chain).
3. If no properties can be assigned, the structure is logged to `unassigned_structures.log` and skipped.

### Inputs

| Path | Source |
|------|--------|
| `Data/Alignment/paired_enantiomers_aligned_structures/` | Alignment step 2.3 |
| `Data/Pharmacophore/paired_enantiomers_pocket_based/` | Pharmacophore step 3.2 |

### Outputs

```
Data/Maestro_property/
├── converted_mae/                         # PDB -> MAE conversion
├── edited_mae/                            # per-structure MAE with properties
├── unassigned_structures.log              # structures skipped (no properties)
└── all_structures_with_properties.mae     # combined output (assigned only)
```

---

## Production run

Requires Schrödinger (`$SCHRODINGER`):

```bash
export SCHRODINGER=/path/to/schrodinger
bash run_maestro_property.sh
bash run_maestro_property.sh --clean
bash run_maestro_property.sh --skip-convert   # skip PDB->MAE conversion, reuse converted_mae/
bash run_maestro_property.sh --overwrite   # re-run structconvert
```

Or directly:

```bash
$SCHRODINGER/run python3 Maestro_property/editing_meastro_property.py
```

---

## Paths

[`paths.py`](paths.py) — respects `PIPELINE_TEST=1` for `test/` outputs.

---

## Related documentation

- [Pharmacophore generation](../Pharmacophore_generation_after_alignment/README.md)
- [Alignment](../Alignment/README.md)
