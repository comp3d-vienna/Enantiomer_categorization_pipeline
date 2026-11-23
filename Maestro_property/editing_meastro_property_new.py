#!/usr/bin/env python3
"""
Maestro property editing for aligned enantiomer structures.

1. Convert aligned PDBs to MAE format (Schrödinger structconvert).
2. Assign Maestro user properties from the m-layer db-id list and the
   UniProt classification file.

   Properties written:
   - s_user_ID: ``{uniprot_field}_{pdb}_{ligand}_{chain}_{res}_mN``
     The UniProt field is kept whole, including commas when several
     accessions are listed. Several matches for one PDB+chain become a list.
   - s_user_group: group token + InChI hash from the aligned filename
"""

from __future__ import annotations

import argparse
import ast
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

_MAESTRO_PROPERTY_DIR = Path(__file__).resolve().parent
if str(_MAESTRO_PROPERTY_DIR) not in sys.path:
    sys.path.insert(0, str(_MAESTRO_PROPERTY_DIR))

import paths

_M_LAYER_SUFFIX = re.compile(r"_m\d+$", re.IGNORECASE)


@dataclass(frozen=True)
class StructureIdentity:
    """Parsed structure identity. ``uniprot_id`` holds the short ``group_id`` token."""

    uniprot_id: str
    inchi_hash: str
    pdb_id: str
    chain_id: str

    @property
    def group_id(self) -> str:
        return f"{self.uniprot_id}_{self.inchi_hash}"

    @property
    def has_group_token(self) -> bool:
        return bool(self.uniprot_id and self.inchi_hash)

    @property
    def is_pdb_only(self) -> bool:
        return bool(self.pdb_id) and not self.has_group_token


def normalize_structure_stem(name: str) -> str:
    """Return basename without extension and without leading ``rot-``."""
    base = Path(name).name
    for ext in (".mae", ".pdb"):
        if base.endswith(ext):
            base = base[: -len(ext)]
            break
    if base.startswith("rot-"):
        base = base[4:]
    return base


def is_prepped_copied_stem(stem: str) -> bool:
    parts = stem.split("_")
    return len(parts) >= 3 and parts[1] == "prepped" and parts[2] == "copied"


def is_reference_prepped_stem(stem: str) -> bool:
    parts = stem.split("_")
    return len(parts) >= 2 and parts[1] == "prepped" and not is_prepped_copied_stem(stem)


def parse_structure_name(name: str) -> StructureIdentity | None:
    """
    Parse aligned structure names such as:

    - ``<group_id>_<inchi>_<pdb>_chain_<chain>``
    - ``<group_id>_<inchi>_<pdb>_chain_<chain>_copied``
    - ``<pdb>_prepped`` / ``rot-<pdb>_prepped`` (reference aligned structure)
    - ``<pdb>_prepped_copied`` (whole prepped protein copy)
    """
    base = normalize_structure_stem(name)
    if base.endswith("_copied") and not is_prepped_copied_stem(base):
        base = base[: -len("_copied")]

    parts = base.split("_")

    if len(parts) >= 5 and parts[3] == "chain":
        return StructureIdentity(
            uniprot_id=parts[0],
            inchi_hash=parts[1],
            pdb_id=parts[2],
            chain_id=parts[4],
        )
    if len(parts) >= 2 and parts[1] == "prepped":
        return StructureIdentity(
            uniprot_id="",
            inchi_hash="",
            pdb_id=parts[0],
            chain_id="",
        )
    return None


def is_chain_copied_stem(stem: str) -> bool:
    """True for per-chain copies such as ``g000612_..._9fa6_chain_A_copied``."""
    return stem.endswith("_copied") and not is_prepped_copied_stem(stem)


def find_sibling_stems_for_pdb_only(
    structure_name: str,
    aligned_dir: Path,
) -> list[str]:
    """
    Find aligned structure stems that share a PDB ID with a reference / copied file.

    ``rot-<pdb>_prepped`` matches per-chain aligned files without ``_copied``.
    ``<pdb>_prepped_copied`` matches per-chain ``*_copied`` files.
    """
    identity = parse_structure_name(structure_name)
    if identity is None or identity.has_group_token:
        return []

    query_stem = normalize_structure_stem(structure_name)
    pdb_key = identity.pdb_id.lower()
    query_is_prepped_copied = is_prepped_copied_stem(query_stem)
    query_is_reference_prepped = is_reference_prepped_stem(query_stem)
    if not query_is_prepped_copied and not query_is_reference_prepped:
        return []

    matches: list[str] = []
    for path in sorted(aligned_dir.glob("*.pdb")):
        sibling_stem = path.stem
        if sibling_stem == Path(structure_name).stem:
            continue
        sibling_identity = parse_structure_name(path.name)
        if sibling_identity is None or not sibling_identity.has_group_token:
            continue
        if sibling_identity.pdb_id.lower() != pdb_key:
            continue

        sibling_is_chain_copied = is_chain_copied_stem(sibling_stem)
        if query_is_prepped_copied and sibling_is_chain_copied:
            matches.append(sibling_stem)
        elif query_is_reference_prepped and not sibling_is_chain_copied:
            matches.append(sibling_stem)

    return matches


def instance_base(db_id: str) -> str:
    """Strip a trailing ``_m<digits>`` layer from a db id."""
    return _M_LAYER_SUFFIX.sub("", db_id.strip())


def pdb_chain_from_db_id(db_id: str) -> tuple[str, str] | None:
    """
    Parse ``{pdb}_{ligand}_{chain}_{res}_mN`` into ``(pdb, chain)``.

    The ligand token is a single ``_``-separated field.
    """
    parts = instance_base(db_id).split("_")
    if len(parts) < 4 or not parts[0] or not parts[2]:
        return None
    return parts[0].lower(), parts[2]


def load_uniprot_by_instance(uniprot_file: Path) -> dict[str, str]:
    """
    Map instance id (column 1, no m-layer) to the UniProt field.

    The UniProt field is the second-to-last whitespace-separated token and is
    kept whole, including commas.
    """
    if not uniprot_file.is_file():
        raise FileNotFoundError(f"UniProt classification file not found: {uniprot_file}")
    lookup: dict[str, str] = {}
    for line_no, raw in enumerate(uniprot_file.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            print(f"[WARNING] Skipping short UniProt row {line_no} in {uniprot_file.name}")
            continue
        lookup[parts[0]] = parts[-2]
    return lookup


def load_user_ids_by_pdb_chain(
    db_ids_file: Path,
    uniprot_by_instance: dict[str, str],
) -> dict[tuple[str, str], list[str]]:
    """
    Index ``{uniprot}_{full_db_id}`` values by ``(pdb_id_lower, chain)``.

    Db ids whose base instance is missing from the UniProt file are skipped.
    """
    if not db_ids_file.is_file():
        raise FileNotFoundError(f"db id file not found: {db_ids_file}")

    by_key: dict[tuple[str, str], list[str]] = {}
    n_missing_uniprot = 0
    n_malformed = 0
    for raw in db_ids_file.read_text(encoding="utf-8").splitlines():
        db_id = raw.strip()
        if not db_id or db_id.startswith("#"):
            continue
        parsed = pdb_chain_from_db_id(db_id)
        if parsed is None:
            n_malformed += 1
            print(f"[WARNING] Skipping malformed db id: {db_id}")
            continue
        base = instance_base(db_id)
        uniprot = uniprot_by_instance.get(base)
        if uniprot is None:
            n_missing_uniprot += 1
            print(f"[WARNING] No UniProt row for db id {db_id}")
            continue
        user_id = f"{uniprot}_{db_id}"
        bucket = by_key.setdefault(parsed, [])
        if user_id not in bucket:
            bucket.append(user_id)

    if n_malformed:
        print(f"[WARNING] Skipped {n_malformed} malformed db id(s)")
    if n_missing_uniprot:
        print(f"[WARNING] Skipped {n_missing_uniprot} db id(s) with no UniProt row")
    return by_key


def merge_maestro_property_dicts(
    props_list: list[dict[str, str]],
    *,
    structure_name: str = "",
) -> dict[str, str] | None:
    """Merge Maestro property dicts from multiple group identifiers into one."""
    if not props_list:
        return None
    if len(props_list) == 1:
        return props_list[0]

    instance_ids: list[str] = []
    group_keys: set[str] = set()

    for props in props_list:
        for raw_id in ast.literal_eval(props["s_user_ID"]):
            if raw_id not in instance_ids:
                instance_ids.append(raw_id)
        group_keys.update(ast.literal_eval(props["s_user_group"]))

    if not instance_ids:
        return None

    label = f" for '{structure_name}'" if structure_name else ""
    print(
        f"[INFO] Merged {len(props_list)} group identifier(s){label}: "
        f"{sorted(group_keys)}"
    )
    return {
        "s_user_ID": str(instance_ids),
        "s_user_group": str(group_keys),
    }


def build_property_values(
    identity: StructureIdentity,
    ids_by_pdb_chain: dict[tuple[str, str], list[str]],
) -> dict[str, str] | None:
    """Build Maestro properties for one group-token structure."""
    if not identity.has_group_token:
        return None
    key = (identity.pdb_id.lower(), identity.chain_id)
    instance_ids = ids_by_pdb_chain.get(key)
    if not instance_ids:
        print(
            f"[WARNING] No db ids for pdb_id={identity.pdb_id!r}, "
            f"chain_id={identity.chain_id!r} (group {identity.group_id!r})"
        )
        return None
    return {
        "s_user_ID": str(list(instance_ids)),
        "s_user_group": str({identity.group_id}),
    }


def _props_for_sibling_stem(
    stem: str,
    props_by_stem: dict[str, dict[str, str]],
    ids_by_pdb_chain: dict[tuple[str, str], list[str]],
) -> dict[str, str] | None:
    cached = props_by_stem.get(stem)
    if cached is not None:
        return cached
    sibling_identity = parse_structure_name(f"{stem}.pdb")
    if sibling_identity is None:
        return None
    return build_property_values(sibling_identity, ids_by_pdb_chain)


def properties_from_sibling_stems(
    sibling_stems: list[str],
    props_by_stem: dict[str, dict[str, str]],
    *,
    structure_name: str,
    ids_by_pdb_chain: dict[tuple[str, str], list[str]],
) -> dict[str, str] | None:
    """
    Resolve properties for pdb-only structures from aligned-folder siblings.

    When siblings belong to different groups, aggregate ``s_user_ID`` and
    ``s_user_group`` from every sibling with properties in each group (one
    output MAE). When siblings share one group (e.g. multiple chains), use
    the first sibling only.
    """
    stems_by_group: dict[str, list[str]] = {}
    for stem in sibling_stems:
        sibling_identity = parse_structure_name(f"{stem}.pdb")
        if sibling_identity is None:
            continue
        stems_by_group.setdefault(sibling_identity.group_id, []).append(stem)

    if len(stems_by_group) > 1:
        group_props: list[dict[str, str]] = []
        for group_id in sorted(stems_by_group):
            per_group_props: list[dict[str, str]] = []
            for stem in stems_by_group[group_id]:
                stem_props = _props_for_sibling_stem(
                    stem,
                    props_by_stem,
                    ids_by_pdb_chain,
                )
                if stem_props is not None:
                    per_group_props.append(stem_props)
            if not per_group_props:
                print(
                    f"[WARNING] No properties for group {group_id!r} siblings of "
                    f"'{structure_name}'"
                )
                continue
            group_merged = merge_maestro_property_dicts(
                per_group_props,
                structure_name=f"{structure_name} ({group_id})",
            )
            if group_merged is not None:
                group_props.append(group_merged)
        return merge_maestro_property_dicts(
            group_props, structure_name=structure_name
        )

    first_stem = sibling_stems[0]
    if len(sibling_stems) > 1:
        print(
            f"[INFO] Multiple aligned-folder siblings for '{structure_name}' "
            f"(same group); using '{first_stem}'"
        )
    else:
        print(
            f"[INFO] Resolved '{structure_name}' from aligned-folder sibling "
            f"'{first_stem}'"
        )
    sibling_props = _props_for_sibling_stem(
        first_stem,
        props_by_stem,
        ids_by_pdb_chain,
    )
    if sibling_props is None:
        print(
            f"[WARNING] Sibling '{first_stem}' has no properties for "
            f"'{structure_name}'"
        )
    return sibling_props


def compute_properties_for_stem(
    structure_name: str,
    ids_by_pdb_chain: dict[tuple[str, str], list[str]],
    *,
    aligned_dir: Path | None = None,
    props_by_stem: dict[str, dict[str, str]] | None = None,
) -> dict[str, str] | None:
    """
    Build Maestro user properties for one structure stem (one output MAE file).

    Group-token structures are resolved from the db-id index by PDB id and
    chain. Reference / copied structures without group tokens inherit
    properties from matching sibling files in ``aligned_dir``. Siblings from
    different groups are merged; siblings from the same group use the first
    match only.
    """
    identity = parse_structure_name(structure_name)
    if identity is None:
        return None

    if identity.has_group_token:
        return build_property_values(identity, ids_by_pdb_chain)

    if aligned_dir is not None and props_by_stem is not None:
        sibling_stems = find_sibling_stems_for_pdb_only(structure_name, aligned_dir)
        if sibling_stems:
            return properties_from_sibling_stems(
                sibling_stems,
                props_by_stem,
                structure_name=structure_name,
                ids_by_pdb_chain=ids_by_pdb_chain,
            )
    return None


def structconvert(pdb_path: Path, mae_path: Path, schrodinger: str) -> None:
    structconvert_bin = os.path.join(schrodinger, "utilities", "structconvert")
    if not os.path.isfile(structconvert_bin):
        raise FileNotFoundError(f"structconvert not found: {structconvert_bin}")
    mae_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [structconvert_bin, str(pdb_path), str(mae_path)],
        check=True,
    )


def convert_aligned_pdbs(
    aligned_dir: Path,
    converted_dir: Path,
    *,
    schrodinger: str,
    overwrite: bool = False,
) -> tuple[list[Path], int, int]:
    """
    Convert aligned PDBs to MAE.

    Existing ``.mae`` files are left unchanged unless ``overwrite`` is True.
    Returns ``(mae_paths, n_converted, n_skipped)``.
    """
    converted_dir.mkdir(parents=True, exist_ok=True)
    converted_files: list[Path] = []
    n_converted = 0
    n_skipped = 0
    for pdb_path in sorted(aligned_dir.glob("*.pdb")):
        mae_path = converted_dir / f"{pdb_path.stem}.mae"
        if mae_path.exists() and not overwrite:
            print(f"[skip convert] {mae_path.name} already exists")
            n_skipped += 1
        else:
            print(f"[convert] {pdb_path.name} -> {mae_path.name}")
            structconvert(pdb_path, mae_path, schrodinger)
            n_converted += 1
        converted_files.append(mae_path)
    return converted_files, n_converted, n_skipped


def _require_structure_module():
    try:
        from schrodinger import structure
    except ImportError as exc:
        raise SystemExit(
            "schrodinger Python module not found. Run with Schrödinger's Python, e.g.\n"
            "  $SCHRODINGER/run python3 Maestro_property/editing_meastro_property_new.py"
        ) from exc
    return structure


def _unassigned_reason(
    structure_name: str,
    identity: StructureIdentity | None,
) -> str:
    if identity is None:
        return "Could not parse structure name"
    return "No properties assigned"


def edit_mae_properties(
    converted_mae_files: list[Path],
    edited_dir: Path,
    combined_output: Path,
    unassigned_log: Path,
    *,
    aligned_dir: Path,
    ids_by_pdb_chain: dict[tuple[str, str], list[str]],
) -> tuple[int, int]:
    structure = _require_structure_module()
    edited_dir.mkdir(parents=True, exist_ok=True)
    combined_output.parent.mkdir(parents=True, exist_ok=True)
    unassigned_log.parent.mkdir(parents=True, exist_ok=True)

    props_by_stem: dict[str, dict[str, str]] = {}
    for mae_path in converted_mae_files:
        if not mae_path.is_file():
            continue
        identity = parse_structure_name(mae_path.stem)
        if identity is None or not identity.has_group_token:
            continue
        props = compute_properties_for_stem(mae_path.stem, ids_by_pdb_chain)
        if props is not None:
            props_by_stem[mae_path.stem] = props

    writer = structure.StructureWriter(str(combined_output))
    edited_count = 0
    skipped_count = 0

    with open(unassigned_log, "w", encoding="utf-8") as log_f:
        log_f.write("structure_name\treason\n")

        for mae_path in converted_mae_files:
            if not mae_path.is_file():
                continue

            structure_name = mae_path.stem
            reader = structure.StructureReader(str(mae_path))
            structures = list(reader)
            if not structures:
                reason = "No structures in MAE file"
                log_f.write(f"{structure_name}\t{reason}\n")
                print(f"[SKIP] {structure_name}: {reason}")
                skipped_count += 1
                continue

            identity = parse_structure_name(structure_name)
            if identity is not None and identity.has_group_token:
                props = props_by_stem.get(structure_name)
            else:
                props = compute_properties_for_stem(
                    structure_name,
                    ids_by_pdb_chain,
                    aligned_dir=aligned_dir,
                    props_by_stem=props_by_stem,
                )

            if props is None:
                reason = _unassigned_reason(structure_name, identity)
                log_f.write(f"{structure_name}\t{reason}\n")
                print(f"[SKIP] {structure_name}: {reason} (logged to {unassigned_log})")
                skipped_count += 1
                continue

            out_path = edited_dir / mae_path.name
            file_writer = structure.StructureWriter(str(out_path))
            for st in structures:
                for key, value in props.items():
                    st.property[key] = value
                file_writer.append(st)
                writer.append(st)
                edited_count += 1
            file_writer.close()

    writer.close()
    return edited_count, skipped_count


def run(
    *,
    aligned_dir: Path,
    converted_dir: Path,
    edited_dir: Path,
    combined_output: Path,
    unassigned_log: Path,
    db_ids_file: Path,
    uniprot_file: Path,
    schrodinger: str,
    overwrite: bool = False,
    skip_convert: bool = False,
) -> None:
    if not aligned_dir.is_dir():
        raise FileNotFoundError(f"Aligned structures directory not found: {aligned_dir}")

    uniprot_by_instance = load_uniprot_by_instance(uniprot_file)
    print(f"Loaded {len(uniprot_by_instance)} UniProt rows from {uniprot_file}")
    ids_by_pdb_chain = load_user_ids_by_pdb_chain(db_ids_file, uniprot_by_instance)
    n_ids = sum(len(values) for values in ids_by_pdb_chain.values())
    print(
        f"Indexed {n_ids} user id(s) for {len(ids_by_pdb_chain)} PDB+chain pair(s) "
        f"from {db_ids_file}"
    )

    if skip_convert:
        converted_files = sorted(converted_dir.glob("*.mae"))
        if not converted_files:
            raise FileNotFoundError(
                f"--skip-convert enabled, but no MAE files found in {converted_dir}"
            )
        print(f"[skip convert] Using {len(converted_files)} existing MAE files.")
    else:
        converted_files, n_converted, n_skipped = convert_aligned_pdbs(
            aligned_dir,
            converted_dir,
            schrodinger=schrodinger,
            overwrite=overwrite,
        )
        if n_skipped:
            print(
                f"MAE conversion: {n_converted} converted, {n_skipped} skipped "
                f"(already in {converted_dir})"
            )
        else:
            print(f"MAE conversion: {n_converted} converted in {converted_dir}")

    edited_count, skipped_count = edit_mae_properties(
        converted_files,
        edited_dir,
        combined_output,
        unassigned_log,
        aligned_dir=aligned_dir,
        ids_by_pdb_chain=ids_by_pdb_chain,
    )
    print(f"Assigned user properties to {edited_count} structure(s)")
    if skipped_count:
        print(f"Skipped {skipped_count} structure(s) without properties: {unassigned_log}")
    print(f"Edited MAE files: {edited_dir}")
    print(f"Combined output: {combined_output}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Convert aligned PDBs to MAE and assign Maestro user properties "
            "from db ids and a UniProt classification file."
        )
    )
    parser.add_argument(
        "--aligned-dir",
        type=Path,
        default=Path(paths.ALIGNED_STRUCTURES_DIR),
        help="Aligned receptor PDB directory (Alignment step 2.3).",
    )
    parser.add_argument(
        "--converted-mae-dir",
        type=Path,
        default=Path(paths.CONVERTED_MAE_DIR),
        help="Output directory for PDB -> MAE conversion.",
    )
    parser.add_argument(
        "--edited-mae-dir",
        type=Path,
        default=Path(paths.EDITED_MAE_NEW_DIR),
        help="Output directory for per-structure MAE files with properties.",
    )
    parser.add_argument(
        "--combined-output",
        type=Path,
        default=Path(paths.COMBINED_MAE_NEW_OUTPUT),
        help="Combined MAE output with all edited structures.",
    )
    parser.add_argument(
        "--unassigned-log",
        type=Path,
        default=Path(paths.UNASSIGNED_STRUCTURES_NEW_LOG),
        help="Log file for structures skipped from edited output (no properties).",
    )
    parser.add_argument(
        "--db-ids",
        type=Path,
        default=Path(paths.DB_IDS_WITH_M_LAYER),
        help="Db ids with an m-layer suffix, one per line (default: Data/db_ids_with_m_layer.txt).",
    )
    parser.add_argument(
        "--uniprot-file",
        type=Path,
        default=Path(paths.FILTERED_ENANTIOMER_WITH_UNIPROT_TXT),
        help=(
            "Classification file; UniProt is the second-to-last field "
            "(default: Data/chiral_classification_results/filtered_enantiomer_with_uniprot_new.txt)."
        ),
    )
    parser.add_argument(
        "--schrodinger",
        type=str,
        default=os.environ.get("SCHRODINGER", ""),
        help="Schrödinger installation root (default: $SCHRODINGER).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-run structconvert even when converted MAE files already exist.",
    )
    parser.add_argument(
        "--skip-convert",
        action="store_true",
        help="Skip PDB->MAE conversion; requires existing converted MAE files.",
    )
    args = parser.parse_args()

    if not args.schrodinger:
        raise SystemExit(
            "SCHRODINGER is not set. Export SCHRODINGER or pass --schrodinger."
        )

    run(
        aligned_dir=args.aligned_dir,
        converted_dir=args.converted_mae_dir,
        edited_dir=args.edited_mae_dir,
        combined_output=args.combined_output,
        unassigned_log=args.unassigned_log,
        db_ids_file=args.db_ids,
        uniprot_file=args.uniprot_file,
        schrodinger=args.schrodinger,
        overwrite=args.overwrite,
        skip_convert=args.skip_convert,
    )


if __name__ == "__main__":
    main()
