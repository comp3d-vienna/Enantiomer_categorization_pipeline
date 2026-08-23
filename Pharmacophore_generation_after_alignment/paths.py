"""Path configuration for the Pharmacophore generation stage."""

import importlib.util
import os

PHARMACOPHORE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(PHARMACOPHORE_DIR)

_preprocess_paths_file = os.path.join(PROJECT_ROOT, "Data_preprocess", "paths.py")
_spec = importlib.util.spec_from_file_location("preprocess_paths", _preprocess_paths_file)
preprocess_paths = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(preprocess_paths)

BASE_DIR = preprocess_paths.BASE_DIR
IS_TEST = preprocess_paths.IS_TEST

ALIGNMENT_DIR = os.path.join(BASE_DIR, "Alignment")
PHARMACOPHORE_DATA_DIR = os.path.join(BASE_DIR, "Pharmacophore")

PAIRED_ENANTIOMERS_DIR = os.path.join(ALIGNMENT_DIR, "paired_enantiomers")
ALIGNED_STRUCTURES_DIR = os.path.join(ALIGNMENT_DIR, "paired_enantiomers_aligned_structures")
LIGAND_EXTRACT_DIR = os.path.join(ALIGNMENT_DIR, "Enantiomer_aligned_structure_ligandextract")
DB_IDS_WITH_M_LAYER = preprocess_paths.DB_IDS_WITH_M_LAYER

CANONICAL_LIGAND_DIR = os.path.join(
    PHARMACOPHORE_DATA_DIR, "Enantiomer_aligned_structure_ligandextract_canonical"
)
INTERACTION_DATA_DIR = os.path.join(PHARMACOPHORE_DATA_DIR, "Interaction_data")
PH4_RESULTS_DIR = os.path.join(PHARMACOPHORE_DATA_DIR, "Pharmacophore_generation_results")

PAIRED_FEATURE_OUTPUT_DIR = os.path.join(
    PHARMACOPHORE_DATA_DIR, "paired_enantiomers_pocket_based"
)
PROCESSING_LOG = os.path.join(PHARMACOPHORE_DATA_DIR, "interaction_data_processing_log.txt")

GEN_FAILED_JOBS = os.path.join(PHARMACOPHORE_DATA_DIR, "generation_failed_jobs.txt")
GEN_RECEPTOR_NOT_FOUND = os.path.join(PHARMACOPHORE_DATA_DIR, "generation_receptor_not_found.txt")
GEN_TIMEOUT = os.path.join(PHARMACOPHORE_DATA_DIR, "generation_timeout.txt")


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
