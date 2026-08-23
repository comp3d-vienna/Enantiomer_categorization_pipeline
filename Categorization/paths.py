"""Path configuration for the Categorization stage."""

import importlib.util
import os

CATEGORIZATION_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(CATEGORIZATION_DIR)

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

CATEGORIZATION_RESULTS_DIR = os.path.join(BASE_DIR, "Categorization")
CATEGORIES_CSV = os.path.join(CATEGORIZATION_RESULTS_DIR, "complex_categories.csv")

INTERACTION_TYPE_INFO_DIR = pharmacophore_paths.PAIRED_FEATURE_OUTPUT_DIR
CANONICAL_LIGAND_DIR = pharmacophore_paths.CANONICAL_LIGAND_DIR

# Manual curation spreadsheets (one CSV per category; filename encodes the label)
MANUAL_CURATION_DIR = os.path.join(preprocess_paths.DATA_DIR, "Manual_curation")

# Training labels (manualcheck_invalid excluded)
THRESHOLD_TUNING_MANUAL_LABELS = ("1.1", "1.2", "1.3", "2")
