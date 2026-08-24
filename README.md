# enantiomer_binding_conformation

A computational pipeline for categorizing small-molecule enantiomer binding conformations.

## Data

Download the data from [Zenodo](https://doi.org/<replace-with-created-DOI>) and place it under `Data/`. Put your mmCIF files (`*.cif.gz`) in `Data/mmCIF/`.

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

You need all of the following before running the pipeline. Creating the conda env is in [Installation](#installation).

| What | Detail |
| ---- | ------ |
| Python 3.13 conda env `categorize_pipeline` | From `environment.yml`. Provides RDKit, BioPython, pandas, CDPKit (CDPL), and the other Python packages. |
| [Schrödinger](https://www.schrodinger.com/) | Commercial license. PrepWizard, `structconvert`, `structalign`. Set `SCHRODINGER` to your install root. |
| [UniCON](https://www.zbh.uni-hamburg.de/forschung/amd/software/unicon.html) | Academic / non-commercial via the [NAOMI ChemBio Suite](https://software.zbh.uni-hamburg.de). Not shipped here. Place the unpacked tree at `Data_preprocess/unicon_1.5.0/` (`…/unicon` binary). Activate with `./unicon --license …`. |
| [LigandExtractor](https://www.zbh.uni-hamburg.de/forschung/amd/software/ligandextractor.html) | Same NAOMI terms as UniCON. Not shipped here. Place the unpacked tree at `Alignment/LigandExtractor_1.0.1/` (`…/LigandExtractor` binary). Activate with `./LigandExtractor --license …`. |
| PDB mmCIF input | Gzipped `*.cif.gz` in **`Data/mmCIF/`**. |
| Network | RCSB / UniProt APIs during Data Preprocess (unless you skip those steps). |

Test data is under `Data/test/mmCIF/`. Run it with:

```bash
enantiomer-pipeline --mmcif-source Data/test/mmCIF all
```

Outputs are written under `Data/test/`.

## Overview

This project develops a computational pipeline to:

1. Retrieve protein structures from the PDB that contain enantiomer pairs.
2. Filter and annotate those structures with ligand and target information.
3. Align paired enantiomer complexes on the same protein (UniProt ID).
4. Exploration of interpretable features to characterize enantiomer binding modes.

## Pipeline stages

| Stage | Directory | Documentation |
| ----- | --------- | ------------- |
| 1. Data Preprocess | [Data_preprocess/](Data_preprocess/) | [README.md](Data_preprocess/README.md) |
| 2. Alignment | [Alignment/](Alignment/) | [README.md](Alignment/README.md) |
| 3. Pharmacophore generation | [Pharmacophore_generation_after_alignment/](Pharmacophore_generation_after_alignment/) | [README.md](Pharmacophore_generation_after_alignment/README.md) |
| 4. Categorization | [Categorization/](Categorization/) | [README.md](Categorization/README.md) |

**Convention:** scripts live in stage directories; pipeline outputs live under [Data/](Data/).

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
├── Data/                         # Pipeline inputs and outputs
├── Data_preprocess/              # Stage 1 scripts
├── Alignment/                    # Stage 2 scripts
├── Pharmacophore_generation_after_alignment/   # Stage 3 scripts
└── Categorization/               # Stage 4 scripts
```

---

## Installation

```bash
git clone git@github.com:Huanni05/enantiomer_binding_conformation.git
cd enantiomer_binding_conformation

conda env create -f environment.yml
conda activate categorize_pipeline
python -m pip install -e .
```

---

## Package usage

### Command-line interface

Run one stage or the full workflow:

```bash
# Single stage
enantiomer-pipeline preprocess
enantiomer-pipeline alignment
enantiomer-pipeline pharmacophore
enantiomer-pipeline categorization

# Full pipeline (stops on first failure)
enantiomer-pipeline all

# Pass flags through to the underlying bash driver (note the -- separator)
enantiomer-pipeline preprocess -- --skip-uniprot --skip-prepwizard
enantiomer-pipeline alignment -- --skip-structalign
enantiomer-pipeline pharmacophore -- --skip-generation
```

Common options:

| Flag                  | Purpose                                                      |
| --------------------- | ------------------------------------------------------------ |
| `--project-root PATH` | Repository root (default: auto-detect from install location) |
| `--schrodinger PATH`  | Set `SCHRODINGER` for stages that need Schrödinger           |
| `--mmcif-source PATH` | Set `MMCIF_SOURCE` for preprocess (default: `Data/mmCIF`) |
| `--python-cmd PATH`   | Set `PYTHON_CMD` for subprocess Python scripts               |
| `--json`              | Print structured result summary (metrics + output paths)     |
| `--continue-on-error` | With `all`, keep running after a failed stage                |

Each stage prints a result summary when it finishes (file counts, category breakdowns, etc.). Use `--json` for machine-readable output.

### Python API

```python
from enantiomer_pipeline import PipelineConfig, run_stage, run_all

config = PipelineConfig(
    schrodinger="/path/to/your/schrodinger",
)

# Run one stage
result = run_stage(config, "categorization")
print(result.ok, result.returncode)
print(result.result["text"])

# Run all stages in order
results = run_all(config, stop_on_error=True)
```

`PipelineConfig.data_dir` is `Data/`. Stage drivers are the shared `run_*.sh` scripts.

## Credits

This package was created with [Copier](https://github.com/copier-org/copier) and the [NLeSC/python-template](https://github.com/NLeSC/python-template).
