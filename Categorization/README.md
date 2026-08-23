# 4. Categorization

Assign manual labels to paired groups, compute closest-pair **RMSD** and **SILIRID similarity**, and score custom RMSD bins.

Needs [Pharmacophore generation](../Pharmacophore_generation_after_alignment/README.md) (steps 3.1–3.2) and `Data/Manual_curation/` (also used in test mode).

Scripts are in this folder. Outputs go under `../Data/Categorization/` (or `../test/Categorization/` when `PIPELINE_TEST=1`). Paths: [`paths.py`](paths.py).

## Run

```bash
conda activate categorize_pipeline
bash run_categorization.sh
```

`--features-only` stops after 4a.  
Test (4a only): `bash ../test/run_categorization_test.sh`  
Delete `Data/Categorization/` yourself for a full redo.

Needs pandas and RDKit.

Earlier stages need [UniCON](https://www.zbh.uni-hamburg.de/forschung/amd/software/unicon.html) under `Data_preprocess/unicon_1.5.0/` and [LigandExtractor](https://www.zbh.uni-hamburg.de/forschung/amd/software/ligandextractor.html) under `Alignment/LigandExtractor_1.0.1/`. Those folders are not in the repository; download them from the [NAOMI ChemBio Suite](https://software.zbh.uni-hamburg.de) and place them there. See [repository Prerequisites](../README.md#prerequisites).

## Steps

| Step | Script | Output |
|------|--------|--------|
| 4a | `whole_process.py` | `complex_categories.csv` (labels + closest-pair RMSD) |
| 4b | `export_threshold_training_table.py` | `threshold_training_table.csv` (RMSD + SILIRID similarity and fingerprints) |
| 4c | `train_category_thresholds.py` | custom RMSD-bin purity and recall on the whole table |

The closest pair is the m0–m1 pose pair with the smallest mean atom distance. SILIRID similarity and the two 160-D count fingerprints for that pair (`closest_cross_tag_silirid_a` / `closest_cross_tag_silirid_b`) are written in step 4b. Slot order is `silirid_fingerprint_slots.csv` (20 amino acids × 8 feature types). Counts are capped at 3, matching the similarity calculation. Step 4c uses RMSD only.

Invalid groups are those listed only in `manualcheck_invalid.csv`. Other groups are training (`1.1` / `1.2` / `1.3` / `2`).

Custom bins: RMSD ≤ 1.25 → `1.1`; 1.25–2.5 → `1.2`; 2.5–11 → `1.3`; > 11 → `2`.

## Outputs (`Data/Categorization/` or `test/Categorization/`)

```
complex_categories.csv
threshold_training_table.csv
threshold_training_table_incomplete_features.csv   # curated but missing RMSD
silirid_fingerprint_slots.csv                      # 160-D SILIRID slot order
custom_rmsd_thresholds.csv
custom_rmsd_threshold_report.txt
```

This is the last pipeline stage. Binding-mode categories and RMSD-bin scores are the reproducible outputs.
