# 4. Categorization

Compute closest-pair **RMSD** and **SILIRID similarity** for every paired-enantiomer group.

Needs [Pharmacophore generation](../Pharmacophore_generation_after_alignment/README.md) (steps 3.1–3.2). The pipeline does not read or write manual labels.

Scripts are in this folder. Outputs go under `../Data/Categorization/`. Paths: [`paths.py`](paths.py).

## Run

From the repository root:

```bash
conda activate categorize_pipeline
bash run_categorization.sh
```

Delete `Data/Categorization/` yourself for a full redo.

Needs pandas and RDKit.

Earlier stages need [UniCON](https://www.zbh.uni-hamburg.de/forschung/amd/software/unicon.html) under `Data_preprocess/unicon_1.5.0/` and [LigandExtractor](https://www.zbh.uni-hamburg.de/forschung/amd/software/ligandextractor.html) under `Alignment/LigandExtractor_1.0.1/`. Those folders are not in the repository; download them from the [NAOMI ChemBio Suite](https://software.zbh.uni-hamburg.de) and place them there. See [repository Prerequisites](../README.md#prerequisites).

## Steps

| Step | Script | Output |
|------|--------|--------|
| 4a | `whole_process.py` | `feature_table.csv` (closest-pair RMSD + atom-count diagnostics) |
| 4b | `export_feature_table.py` | same `feature_table.csv`, with SILIRID similarity and fingerprints added |

The closest pair is the m0–m1 pose pair with the smallest mean atom distance. SILIRID similarity and the two 160-D count fingerprints for that pair (`closest_cross_tag_silirid_a` / `closest_cross_tag_silirid_b`) are written in step 4b. Slot order is `silirid_fingerprint_slots.csv` (20 amino acids × 8 feature types). Each slot count is capped at 3 (`silirid_count_cap`), matching the similarity calculation.

## Outputs (`Data/Categorization/`)

```
feature_table.csv
feature_table_incomplete.csv                   # groups missing RMSD
silirid_fingerprint_slots.csv                  # 160-D SILIRID slot order
```

This is the last pipeline stage. Closest-pair RMSD and SILIRID features are the reproducible outputs.
