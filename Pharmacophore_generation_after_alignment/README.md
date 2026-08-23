# 3. Pharmacophore generation

Build interaction pharmacophores from aligned paired complexes and write one feature table per group.

Needs [Alignment](../Alignment/README.md) through step 2.4.

Use one **Python 3.13** env with CDPL, pandas, and pytz (recommended: `categorize_pipeline`). Resolver: [`resolve_python.sh`](resolve_python.sh). Paths: [`paths.py`](paths.py).

## Run

```bash
conda activate categorize_pipeline
bash run_pharmacophore.sh
```

If CDPL is not on that env’s path:

```bash
export PYTHONPATH=/data/shared/software/CDPKit-head-RH9/Python
export PYTHON_CMD=$(conda run -n categorize_pipeline which python)
bash run_pharmacophore.sh
```

Test: `bash ../test/run_pharmacophore_test.sh`  
Re-run only step 3.2: `--skip-generation`

## Steps

| Step | Script | Output |
|------|--------|--------|
| 3.1 | `pharmacophore_generation_new.sh` | Canonical SDFs, interaction TSVs (and PML) |
| 3.2 | `process_interaction_data.py` | One `{group_id}_inchi_{hash}.csv` per group |

Step 3.2 needs at least two valid TSVs per group. It does not split consistent vs inconsistent features.

Existing canonical SDFs, PML, and TSVs are skipped. Delete `Data/Pharmacophore/` for a full redo.

## Outputs (`Data/Pharmacophore/` or `test/Pharmacophore/`)

```
Enantiomer_aligned_structure_ligandextract_canonical/   # 3.1
Interaction_data/                                       # 3.1 TSV
Pharmacophore_generation_results/                       # 3.1 PML
paired_enantiomers_pocket_based/                        # 3.2 CSV
```

Next: [Categorization](../Categorization/README.md)
