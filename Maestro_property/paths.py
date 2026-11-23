"""Path configuration for the Maestro property stage."""

import importlib.util
import os

MAESTRO_PROPERTY_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(MAESTRO_PROPERTY_DIR)

_preprocess_paths_file = os.path.join(PROJECT_ROOT, "Data_preprocess", "paths.py")
_spec = importlib.util.spec_from_file_location("preprocess_paths", _preprocess_paths_file)
preprocess_paths = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(preprocess_paths)

_pharmacophore_paths_file = os.path.join(
    PROJECT_ROOT, "Pharmacophore_generation_after_alignment", "paths.py"
)
_pspec = importlib.util.spec_from_file_location("pharmacophore_paths", _pharmacophore_paths_file)
pharmacophore_paths = importlib.util.module_from_spec(_pspec)
_pspec.loader.exec_module(pharmacophore_paths)

BASE_DIR = preprocess_paths.BASE_DIR
IS_TEST = preprocess_paths.IS_TEST

MAESTRO_PROPERTY_DATA_DIR = os.path.join(BASE_DIR, "Maestro_property")

# Alignment outputs needed for translating internal `group_id` -> original UniProt ID.
ALIGNMENT_DATA_DIR = os.path.join(BASE_DIR, "Alignment")
GROUP_MANIFEST = os.path.join(ALIGNMENT_DATA_DIR, "group_manifest.tsv")

# Inputs from prior stages
ALIGNED_STRUCTURES_DIR = pharmacophore_paths.ALIGNED_STRUCTURES_DIR
PAIRED_POCKET_BASED_DIR = pharmacophore_paths.PAIRED_FEATURE_OUTPUT_DIR

# Step 5 outputs
CONVERTED_MAE_DIR = os.path.join(MAESTRO_PROPERTY_DATA_DIR, "converted_mae")
EDITED_MAE_DIR = os.path.join(MAESTRO_PROPERTY_DATA_DIR, "edited_mae")
COMBINED_MAE_OUTPUT = os.path.join(
    MAESTRO_PROPERTY_DATA_DIR, "all_structures_with_properties.mae"
)
UNASSIGNED_STRUCTURES_LOG = os.path.join(
    MAESTRO_PROPERTY_DATA_DIR, "unassigned_structures.log"
)


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
