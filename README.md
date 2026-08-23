# enantiomer_binding_conformation

A computational pipeline for categorizing small-molecule enantiomer binding conformations.

## Badges

(Customize these badges with your own links, and check https://shields.io/ or https://badgen.net/ to see which other badges are available.)

| fair-software.eu recommendations | |
| :-- | :--  |
| (1/5) code repository              | [![github repo badge](https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/Huanni05/enantiomer_binding_conformation) |
| (2/5) license                      | [![github license badge](https://img.shields.io/github/license/Huanni05/enantiomer_binding_conformation)](https://github.com/Huanni05/enantiomer_binding_conformation) |
| (3/5) community registry           | [![RSD](https://img.shields.io/badge/rsd-enantiomer_binding_conformation-00a3e3.svg)](https://www.research-software.nl/software/enantiomer_binding_conformation) |
| (4/5) citation                     | [![DOI](https://zenodo.org/badge/DOI/<replace-with-created-DOI>.svg)](https://doi.org/<replace-with-created-DOI>)|
| (5/5) checklist                    | [![workflow cii badge](https://bestpractices.coreinfrastructure.org/projects/<replace-with-created-project-identifier>/badge)](https://bestpractices.coreinfrastructure.org/projects/<replace-with-created-project-identifier>) |
| howfairis                          | [![fair-software badge](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu) |
| **Other best practices**           | &nbsp; |
| **GitHub Actions**                 | &nbsp; |
| Build                              | [![build](https://github.com/Huanni05/enantiomer_binding_conformation/actions/workflows/build.yml/badge.svg)](https://github.com/Huanni05/enantiomer_binding_conformation/actions/workflows/build.yml) |

## Prerequisites

**Python environment** — use one **Python 3.13** conda env for all Python pipeline steps:

```bash
conda activate categorize_pipeline
python -m pip install -e .
```

`categorize_pipeline` includes CDPL (required for pharmacophore step 3.1), RDKit, pandas, and BioPython. The editable install from this clone also registers the `enantiomer-pipeline` CLI (see [Package usage](#package-usage)).

### Licensed third-party software

These tools are **not** part of this repository’s license and are **not** available from PyPI. You need your own install and a license that covers your use.

| Software | License / availability | Expected location |
| -------- | ---------------------- | ----------------- |
| [Schrödinger](https://www.schrodinger.com/) (PrepWizard, `structconvert`, `structalign`) | Commercial Schrödinger license | Set `SCHRODINGER` to your install root |
| [UniCON](https://www.zbh.uni-hamburg.de/forschung/amd/software/unicon.html) | Academic / non-commercial via the [NAOMI ChemBio Suite](https://software.zbh.uni-hamburg.de); evaluation license for non-academic users | Place the unpacked tree at `Data_preprocess/unicon_1.5.0/` so the binary is `Data_preprocess/unicon_1.5.0/unicon` |
| [LigandExtractor](https://www.zbh.uni-hamburg.de/forschung/amd/software/ligandextractor.html) | Same NAOMI ChemBio Suite terms as UniCON | Place the unpacked tree at `Alignment/LigandExtractor_1.0.1/` so the binary is `Alignment/LigandExtractor_1.0.1/LigandExtractor` |

These two tools are **not shipped** in this repository. Download them from [Universität Hamburg ZBH](https://software.zbh.uni-hamburg.de) after registration, then put the folders under `Data_preprocess/` and `Alignment/` as above. Activate each binary with your NAOMI license (`./unicon --license …` and `./LigandExtractor --license …`).

### Other dependencies

| Dependency | Stages |
| ---------- | ------ |
| RDKit, Python 3.13, `pandas`, `pytz`, `requests`, `tqdm` | Data Preprocess, Pharmacophore, Categorization |
| BioPython | Alignment |
| CDPL / CDPKit (pip, in `categorize_pipeline`) | Pharmacophore step 3.1 |
| PDB mmCIF (`Data/mmCIF/` by default) | Data Preprocess |
| Network (RCSB, UniProt APIs) | Data Preprocess |

Stage-specific setup and paths are in each stage README.

## Overview

This project develops a computational pipeline to:

1. Retrieve protein structures from the PDB that contain enantiomer pairs.
2. Filter and annotate those structures with ligand and target information.
3. Align paired enantiomer complexes on the same protein (UniProt ID).
4. Use evaluation metrics to automatically categorize enantiomer binding modes.

## Pipeline stages

| Stage | Directory | Documentation |
| ----- | --------- | ------------- |
| 1. Data Preprocess | [Data_preprocess/](Data_preprocess/) | [README.md](Data_preprocess/README.md) |
| 2. Alignment | [Alignment/](Alignment/) | [README.md](Alignment/README.md) |
| 3. Pharmacophore generation | [Pharmacophore_generation_after_alignment/](Pharmacophore_generation_after_alignment/) | [README.md](Pharmacophore_generation_after_alignment/README.md) |
| 4. Categorization | [Categorization/](Categorization/) | [README.md](Categorization/README.md) |

**Convention:** scripts live in stage directories; pipeline outputs live under [Data/](Data/). Set `PIPELINE_TEST=1` to write test outputs under [test/](test/) instead.

The [enantiomer_pipeline](enantiomer_pipeline/) package provides a unified CLI and Python API that wraps the same `run_*.sh` drivers documented below.

---

## Directory layout

```
Enantiomers_binding_conformation/
├── README.md
├── pyproject.toml
├── environment.yml               # Python 3.13 conda env
├── requirements.txt
├── enantiomer_pipeline/          # Unified Python package (CLI + API)
├── src/enantiomer_binding_conformation/
├── Data/                         # Production pipeline outputs
├── test/                         # Test sample + test outputs (PIPELINE_TEST=1)
├── Data_preprocess/              # Stage 1 scripts
├── Alignment/                    # Stage 2 scripts
├── Pharmacophore_generation_after_alignment/   # Stage 3 scripts
└── Categorization/               # Stage 4 scripts
```

---

## Installation

This pipeline is distributed from GitHub, not PyPI. Several stages need licensed third-party software that you must obtain yourself (see [Prerequisites](#prerequisites)).

Use one **Python 3.13** conda env (`categorize_pipeline`) for all Python pipeline steps.

```bash
git clone git@github.com:Huanni05/enantiomer_binding_conformation.git
cd enantiomer_binding_conformation

conda env create -f environment.yml   # skip if the env already exists
conda activate categorize_pipeline
python -m pip install -e .
```

This registers the `enantiomer-pipeline` command and the importable `enantiomer_pipeline` module from the clone. CDPKit (CDPL), RDKit, pandas, and BioPython are installed in the env.

---

## Package usage

### Command-line interface

Run one stage or the full workflow:

```bash
# Single stage (production)
enantiomer-pipeline preprocess
enantiomer-pipeline alignment
enantiomer-pipeline pharmacophore
enantiomer-pipeline categorization

# Full pipeline (stops on first failure)
enantiomer-pipeline all

# Test mode (uses test/run_*_test.sh and test/ outputs)
enantiomer-pipeline --test categorization

# Pass flags through to the underlying bash driver (note the -- separator)
enantiomer-pipeline --test preprocess -- --skip-uniprot --skip-prepwizard
enantiomer-pipeline --test alignment -- --skip-structalign
enantiomer-pipeline --test pharmacophore -- --skip-generation
```

Equivalent module invocation:

```bash
python -m enantiomer_pipeline categorization
```

Common options:

| Flag                  | Purpose                                                      |
| --------------------- | ------------------------------------------------------------ |
| `--test`              | Use `test/run_*_test.sh` drivers and write under `test/`     |
| `--project-root PATH` | Repository root (default: auto-detect from install location) |
| `--schrodinger PATH`  | Set `SCHRODINGER` for stages that need Schrödinger           |
| `--mmcif-source PATH` | Set `MMCIF_SOURCE` for preprocess (default: `Data/mmCIF`)    |
| `--python-cmd PATH`   | Set `PYTHON_CMD` for subprocess Python scripts               |
| `--json`              | Print structured result summary (metrics + output paths)     |
| `--continue-on-error` | With `all`, keep running after a failed stage                |

Each stage prints a result summary when it finishes (file counts, category breakdowns, etc.). Use `--json` for machine-readable output.

### Python API

```python
from enantiomer_pipeline import PipelineConfig, run_stage, run_all

config = PipelineConfig(
    test_mode=True,
    schrodinger="/path/to/your/schrodinger",
)

# Run one stage
result = run_stage(config, "categorization")
print(result.ok, result.returncode)
print(result.result["text"])

# Run all stages in order
results = run_all(config, stop_on_error=True)
```

`PipelineConfig.data_dir` resolves to `Data/` (production) or `test/` (when `test_mode=True`). Stage drivers are selected automatically from `run_*.sh` or `test/run_*_test.sh`.

---

## 1. Data Preprocess

Extract ligands from PDB, classify enantiomer pairs, filter drug-like compounds, map UniProt IDs, and prepare proteins with Schrödinger PrepWizard.

| Step     | Script                               | Output (high level)                                   |
| -------- | ------------------------------------ | ----------------------------------------------------- |
| 1.1a     | `extract_split_ligands.sh`           | Separated ligand SDFs                                 |
| 1.1b     | `convert_sdf_to_inchi.py`            | InChI dictionary                                      |
| 1.1c     | `classify_chiral_compounds.py`       | Enantiomer pairs + SDFs                               |
| 1.2–1.3  | `filter_and_annotate_enantiomers.py` | Filtered SDFs, UniProt map, `db_ids_with_m_layer.txt` |
| 1.4      | `maestro_protein_preparation.sh`     | Prepped MAE structures                                |

**Production run:** set `SCHRODINGER`, `MMCIF_SOURCE`, and `MMCIF_RENAME` to your local paths, then:

```bash
export SCHRODINGER=/path/to/your/schrodinger
export MMCIF_SOURCE=/path/to/your/mmCIF
export MMCIF_RENAME=/path/to/your/mmCIF_rename
bash run_data_preprocess.sh
```

If `MMCIF_*` are unset, place `.cif.gz` files in `Data/mmCIF/`.

**Test run:**

```bash
bash test/run_data_preprocess_test.sh
bash test/run_data_preprocess_test.sh --skip-uniprot --skip-prepwizard
```

Details: **[Data_preprocess/README.md](Data_preprocess/README.md)**

---

## 2. Alignment

Group paired enantiomers (same UniProt, both `m0` and `m1`), parse protein chains, align with Schrödinger, and extract ligands from aligned structures.

| Step    | Summary                                                                                           |
| ------- | ------------------------------------------------------------------------------------------------- |
| 2.1–2.2 | [group_and_parse_chains.py](Alignment/group_and_parse_chains.py) — group, MAE→PDB, parse chains |
| 2.3     | Schrödinger `structalign` — [schrodinger.sh](Alignment/schrodinger.sh) |
| 2.4     | LigandExtractor — [ligand_extractor_prep.sh](Alignment/ligand_extractor_prep.sh) |

**Production run:** set `SCHRODINGER` to your local path, then:

```bash
export SCHRODINGER=/path/to/your/schrodinger
bash run_alignment.sh
```

**Test run:**

```bash
bash test/run_alignment_test.sh
bash test/run_alignment_test.sh --skip-structalign --skip-mae-convert
```

Details: **[Alignment/README.md](Alignment/README.md)**

---

## 3. Pharmacophore generation

Generate interaction pharmacophores from aligned structures and collect interaction features for each paired enantiomer group.

| Step | Script                                                                              |
| ---- | ----------------------------------------------------------------------------------- |
| 3.1  | `pharmacophore_generation_new.sh` — interaction TSV from aligned PDBs + ligand SDFs |
| 3.2  | `process_interaction_data.py` — one feature CSV per paired group                     |

**Production run:**

```bash
conda activate categorize_pipeline
bash run_pharmacophore.sh
```

**Test run:**

```bash
bash test/run_pharmacophore_test.sh
bash test/run_pharmacophore_test.sh --skip-generation
```

Details: **[Pharmacophore_generation_after_alignment/README.md](Pharmacophore_generation_after_alignment/README.md)**

---

## 4. Categorization

Classify paired enantiomer groups into binding-mode categories using closest-pair RMSD (custom bins) and SILIRID similarity (exported with the training table).

| Step | Script |
| ---- | ------ |
| 4a | [whole_process.py](Categorization/whole_process.py) — labels + closest-pair RMSD |
| 4b | [export_threshold_training_table.py](Categorization/export_threshold_training_table.py) — RMSD + SILIRID similarity and fingerprints |
| 4c | [train_category_thresholds.py](Categorization/train_category_thresholds.py) — custom RMSD bins, purity and recall |

**Production run:**

```bash
bash run_categorization.sh
```

**Test run:**

```bash
bash test/run_categorization_test.sh
```

Details: **[Categorization/README.md](Categorization/README.md)**

---

## Credits

This package was created with [Copier](https://github.com/copier-org/copier) and the [NLeSC/python-template](https://github.com/NLeSC/python-template).
