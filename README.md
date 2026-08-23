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

### Input data (mmCIF)

Put your PDB mmCIF files in **`Data/mmCIF/`** (from the repository root). The files must be gzip-compressed **`.cif.gz`**.

```
Enantiomers_binding_conformation/
└── Data/
    └── mmCIF/                 # you create this and place *.cif.gz here
        ├── 1abc.cif.gz
        └── …
```

`Data/mmCIF_rename/` is written by the pipeline (decompressed `.cif`); you do not need to fill it. To use another archive location, set `MMCIF_SOURCE` (and optionally `MMCIF_RENAME`) or pass `--mmcif-source`. In test mode the default input is `test/mmcif/`.

### Other dependencies

| Dependency | Stages |
| ---------- | ------ |
| RDKit, Python 3.13, `pandas`, `pytz`, `requests`, `tqdm` | Data Preprocess, Pharmacophore, Categorization |
| BioPython | Alignment |
| CDPL / CDPKit (pip, in `categorize_pipeline`) | Pharmacophore step 3.1 |
| PDB mmCIF (`Data/mmCIF/*.cif.gz`) | Data Preprocess |
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

The same commands work if you call the package with Python. Put the same stage and flags after `python -m enantiomer_pipeline`:

```bash
python -m enantiomer_pipeline …
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

## Credits

This package was created with [Copier](https://github.com/copier-org/copier) and the [NLeSC/python-template](https://github.com/NLeSC/python-template).
