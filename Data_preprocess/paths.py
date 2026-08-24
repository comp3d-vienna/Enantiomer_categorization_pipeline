"""Central path configuration for pipeline data under Data/ (or Data/test/)."""

import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)


def _abs(path):
    if not os.path.isabs(path):
        path = os.path.join(PROJECT_ROOT, path)
    return os.path.normpath(os.path.abspath(path))


PRODUCTION_DATA = _abs(os.path.join(PROJECT_ROOT, "Data"))
TRIAL_DATA = _abs(os.path.join(PRODUCTION_DATA, "test"))
TRIAL_MMCIF = _abs(os.path.join(TRIAL_DATA, "mmCIF"))

IS_TEST = os.environ.get("PIPELINE_TEST", "") == "1"


def _select_data_dir():
    if IS_TEST:
        return _abs(os.path.join(PROJECT_ROOT, "test"))
    override = os.environ.get("PIPELINE_DATA_DIR")
    if override:
        return _abs(override)
    mmcif = os.environ.get("MMCIF_SOURCE")
    if mmcif and _abs(mmcif) == TRIAL_MMCIF:
        return TRIAL_DATA
    return PRODUCTION_DATA


DATA_DIR = _select_data_dir()
BASE_DIR = DATA_DIR

# Step 1.1 — ligand extraction
if IS_TEST:
    LIGAND_SDF_DIR = os.path.join(DATA_DIR, "ligand_sdf")
    SEPARATED_SDF_DIR = os.path.join(DATA_DIR, "ligand_sdf_separated")
    INCHI_DIR = os.path.join(DATA_DIR, "ligand_sdf_separated_InChI")
else:
    LIGAND_SDF_DIR = os.path.join(DATA_DIR, "ligand_sdf_structure_from_mmCIF")
    SEPARATED_SDF_DIR = os.path.join(DATA_DIR, "ligand_sdf_structure_from_mmCIF_separated")
    INCHI_DIR = os.path.join(DATA_DIR, "ligand_sdf_structure_from_mmCIF_separated_InChI")

MERGED_LIGANDS_SDF = os.path.join(INCHI_DIR, "merged_ligands.sdf")
ALL_INCHIS_DICT = os.path.join(INCHI_DIR, "all_inchis_dict.txt")
TIMER_DIR = os.path.join(BASE_DIR, "timer_mmCIF")
EXTRACT_ERROR_LOG = os.path.join(BASE_DIR, "extract_split_procedure_errors.log")

# Step 1.1c — chiral classification
CHIRAL_RESULTS_DIR = os.path.join(BASE_DIR, "chiral_classification_results")
CHIRAL_SDF_DIR = os.path.join(CHIRAL_RESULTS_DIR, "sdf")
DRUG_LIKENESS_DIR = os.path.join(CHIRAL_RESULTS_DIR, "Drug_likeness_evaluation")

# Step 1.2–1.3 — filtering and annotation
ENANTIOMER_STRUCTURE_SDF = os.path.join(CHIRAL_SDF_DIR, "enantiomer_structure.sdf")
ENANTIOMER_REPRESENTATIVE_SDF = os.path.join(CHIRAL_SDF_DIR, "enantiomer_representative.sdf")
ENANTIOMERS_TXT = os.path.join(CHIRAL_RESULTS_DIR, "enantiomers.txt")
FILTERED_ENANTIOMER_STRUCTURE_SDF = os.path.join(
    DRUG_LIKENESS_DIR, "filtered_enantiomer_structure.sdf"
)
FILTERED_OUT_ENANTIOMER_STRUCTURE_SDF = os.path.join(
    DRUG_LIKENESS_DIR, "filtered_out_enantiomer_structure.sdf"
)
FILTERED_ENANTIOMER_REPRESENTATIVE_SDF = os.path.join(
    DRUG_LIKENESS_DIR, "filtered_enantiomer_representative.sdf"
)
FILTERED_OUT_ENANTIOMER_REPRESENTATIVE_SDF = os.path.join(
    DRUG_LIKENESS_DIR, "filtered_out_enantiomer_representative.sdf"
)
FILTERED_ENANTIOMERS_TXT = os.path.join(CHIRAL_RESULTS_DIR, "filtered_enantiomers.txt")
FILTERED_ENANTIOMER_WITH_UNIPROT_TXT = os.path.join(
    CHIRAL_RESULTS_DIR, "filtered_enantiomer_with_uniprot_new.txt"
)

# Step 1.3b — structure identifiers
DB_IDS_WITH_M_LAYER = os.path.join(BASE_DIR, "db_ids_with_m_layer.txt")

# Step 1.4 — protein preparation
MAE_OUTPUT_DIR = os.path.join(BASE_DIR, "Enantiomer_pdbstructure_maeformat")
PREPWIZARD_OUTPUT_DIR = os.path.join(
    BASE_DIR, "Enantiomer_protein_preperation", "Enantiomer_prepwizard_results"
)
PREPWIZARD_TIMEOUT_LOG = os.path.join(BASE_DIR, "prepwizard_timeout_pdbids.txt")
PREPWIZARD_FAILED_LOG = os.path.join(BASE_DIR, "prepwizard_failed_pdbids.txt")

_raw_mmcif = os.environ.get("MMCIF_SOURCE")
_raw_rename = os.environ.get("MMCIF_RENAME")
MMCIF_SOURCE_DIR = _abs(_raw_mmcif) if _raw_mmcif else os.path.join(DATA_DIR, "mmCIF")
PDB_MMCIF_RENAME_DIR = (
    _abs(_raw_rename) if _raw_rename else os.path.join(DATA_DIR, "mmCIF_rename")
)


def ensure_parent_dir(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
