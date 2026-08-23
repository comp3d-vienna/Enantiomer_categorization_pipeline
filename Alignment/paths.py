"""Path configuration for the Alignment stage."""

import importlib.util
import os

ALIGNMENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(ALIGNMENT_DIR)

_preprocess_paths_file = os.path.join(PROJECT_ROOT, "Data_preprocess", "paths.py")
_spec = importlib.util.spec_from_file_location("preprocess_paths", _preprocess_paths_file)
preprocess_paths = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(preprocess_paths)

BASE_DIR = preprocess_paths.BASE_DIR
IS_TEST = preprocess_paths.IS_TEST

ALIGNMENT_DIR_DATA = os.path.join(BASE_DIR, "Alignment")

FILTERED_ENANTIOMER_WITH_UNIPROT_TXT = preprocess_paths.FILTERED_ENANTIOMER_WITH_UNIPROT_TXT
PREPWIZARD_OUTPUT_DIR = preprocess_paths.PREPWIZARD_OUTPUT_DIR

PAIRED_ENANTIOMERS_DIR = os.path.join(ALIGNMENT_DIR_DATA, "paired_enantiomers")
SINGLE_ENANTIOMERS_DIR = os.path.join(ALIGNMENT_DIR_DATA, "single_enantiomers")
GROUPING_SUMMARY = os.path.join(ALIGNMENT_DIR_DATA, "grouping_summary.txt")
GROUP_MANIFEST = os.path.join(ALIGNMENT_DIR_DATA, "group_manifest.tsv")
PREP_PDB_DIR = os.path.join(ALIGNMENT_DIR_DATA, "Enantiomer_prep_pdbformat")
PARSED_CHAIN_DIR = os.path.join(ALIGNMENT_DIR_DATA, "paired_enantiomers_parsed_chain")
LIST_WITH_UNIPROT_ID = os.path.join(PARSED_CHAIN_DIR, "list_with_uniprot_id.txt")

ALIGNED_STRUCTURES_DIR = os.path.join(ALIGNMENT_DIR_DATA, "paired_enantiomers_aligned_structures")
LIGAND_EXTRACT_DIR = os.path.join(ALIGNMENT_DIR_DATA, "Enantiomer_aligned_structure_ligandextract")
LIGAND_EXTRACT_LOG = os.path.join(ALIGNMENT_DIR_DATA, "ligand_extractor_prep_errors.log")

DB_IDS_WITH_M_LAYER = preprocess_paths.DB_IDS_WITH_M_LAYER


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
