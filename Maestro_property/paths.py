"""Path configuration for the Maestro property stage (editing_meastro_property_new)."""

import importlib.util
import os

MAESTRO_PROPERTY_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(MAESTRO_PROPERTY_DIR)

_preprocess_paths_file = os.path.join(PROJECT_ROOT, "Data_preprocess", "paths.py")
_spec = importlib.util.spec_from_file_location("preprocess_paths", _preprocess_paths_file)
preprocess_paths = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(preprocess_paths)

BASE_DIR = preprocess_paths.BASE_DIR
IS_TEST = preprocess_paths.IS_TEST

MAESTRO_PROPERTY_DATA_DIR = os.path.join(BASE_DIR, "Maestro_property")

# Alignment outputs
ALIGNMENT_DATA_DIR = os.path.join(BASE_DIR, "Alignment")
GROUP_MANIFEST = os.path.join(ALIGNMENT_DATA_DIR, "group_manifest.tsv")

# Aligned structures from the Alignment stage
ALIGNED_STRUCTURES_DIR = os.path.join(BASE_DIR, "Alignment", "paired_enantiomers_aligned_structures")

# Inputs from Data/
DATA_DIR = os.path.join(PROJECT_ROOT, "Data")
DB_IDS_WITH_M_LAYER = os.path.join(DATA_DIR, "db_ids_with_m_layer.txt")
FILTERED_ENANTIOMER_WITH_UNIPROT_TXT = os.path.join(
    DATA_DIR, "chiral_classification_results", "filtered_enantiomer_with_uniprot_new.txt"
)

# Step 5 outputs (_new suffix aligns with editing_meastro_property_new.py)
CONVERTED_MAE_DIR = os.path.join(MAESTRO_PROPERTY_DATA_DIR, "converted_mae")
EDITED_MAE_NEW_DIR = os.path.join(MAESTRO_PROPERTY_DATA_DIR, "edited_mae_new")
COMBINED_MAE_NEW_OUTPUT = os.path.join(
    MAESTRO_PROPERTY_DATA_DIR, "all_structures_with_properties_new.mae"
)
UNASSIGNED_STRUCTURES_NEW_LOG = os.path.join(
    MAESTRO_PROPERTY_DATA_DIR, "unassigned_structures_new.log"
)


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
