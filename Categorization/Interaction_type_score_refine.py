from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import pandas as pd


FEATURE_TYPE_INDEX = {
    "HYDROPHOBIC": 0,
    "AROMATIC": 1,
    "NEGATIVE_IONIZABLE": 2,
    "POSITIVE_IONIZABLE": 3,
    "H_BOND_DONOR": 4,
    "H_BOND_ACCEPTOR": 5,
    "HALOGEN_BOND_DONOR": 6,
    "HALOGEN_BOND_ACCEPTOR": 7,
}

N_FEATURE_TYPES = len(FEATURE_TYPE_INDEX)

AMINO_ACID_ORDER = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
]
AA_TO_POS = {aa: i for i, aa in enumerate(AMINO_ACID_ORDER)}
POCKET_RESIDUE_BITS_LEN = len(AMINO_ACID_ORDER)

SILIRID_SIMILARITY_METRIC = "structure_silirid_similarity"
STRUCTURE_SILIRID_LEN = POCKET_RESIDUE_BITS_LEN * N_FEATURE_TYPES
STRUCTURE_COUNT_RESIDUE_DEFAULT_CAP: int | None = 5
FEATURE_TYPES_IN_ORDER = [
    name for name, _ in sorted(FEATURE_TYPE_INDEX.items(), key=lambda item: item[1])
]


@dataclass(frozen=True)
class AtomSignature:
    atom_index: int
    feature_bits: str
    pocket_residue_bits: tuple[str, ...]
    pocket_residues: tuple[frozenset[str], ...]


def _empty_bits(length: int) -> str:
    return "0" * length


def _set_bit(bits: list[str], index: int) -> None:
    bits[index] = "1"


def normalize_pocket_chunk(chunk: str) -> str:
    chunk = chunk.strip()
    if not chunk:
        return ""
    parts = chunk.split("_")
    if len(parts) > 2:
        return "_".join(parts[:-1]).upper()
    return chunk.upper()


def parse_normalized_pocket_residues(pocket_residues: object) -> set[str]:
    if pocket_residues is None or (isinstance(pocket_residues, float) and pd.isna(pocket_residues)):
        return set()
    text = str(pocket_residues).strip()
    if not text:
        return set()
    residues: set[str] = set()
    for chunk in text.split(","):
        normalized = normalize_pocket_chunk(chunk)
        if normalized:
            residues.add(normalized)
    return residues


def amino_acid_codes_from_residue_set(residues: frozenset[str]) -> set[str]:
    return {
        residue.split("_", 1)[0]
        for residue in residues
        if residue.split("_", 1)[0] in AA_TO_POS
    }


def normalize_feature_type(feature_type: object) -> str | None:
    if feature_type is None or (isinstance(feature_type, float) and pd.isna(feature_type)):
        return None
    name = str(feature_type).strip().upper()
    return name if name in FEATURE_TYPE_INDEX else None


def encode_atom_rows(rows: pd.DataFrame) -> AtomSignature:
    feature_bits = list(_empty_bits(N_FEATURE_TYPES))
    pocket_residue_sets: list[set[str]] = [set() for _ in range(N_FEATURE_TYPES)]

    atom_index = int(rows["atom_index"].iloc[0])
    for _, row in rows.iterrows():
        ftype = normalize_feature_type(row["feature_type"])
        if ftype is None:
            continue
        ft_idx = FEATURE_TYPE_INDEX[ftype]
        _set_bit(feature_bits, ft_idx)
        pocket_residue_sets[ft_idx].update(
            parse_normalized_pocket_residues(row["pocket_residues"])
        )

    pocket_residues = tuple(frozenset(s) for s in pocket_residue_sets)
    pocket_residue_bits: list[str] = []
    for residues in pocket_residues:
        bits = list(_empty_bits(POCKET_RESIDUE_BITS_LEN))
        for code in amino_acid_codes_from_residue_set(residues):
            _set_bit(bits, AA_TO_POS[code])
        pocket_residue_bits.append("".join(bits))

    return AtomSignature(
        atom_index=atom_index,
        feature_bits="".join(feature_bits),
        pocket_residue_bits=tuple(pocket_residue_bits),
        pocket_residues=pocket_residues,
    )


def build_file_signatures(df: pd.DataFrame, interaction_file: str) -> dict[int, AtomSignature]:
    subset = df[df["file"] == interaction_file]
    if subset.empty:
        return {}
    return {
        int(atom_index): encode_atom_rows(atom_rows)
        for atom_index, atom_rows in subset.groupby("atom_index", sort=True)
    }


def min_max_similarity(counts_a: Sequence[int], counts_b: Sequence[int]) -> float:
    if len(counts_a) != len(counts_b):
        raise ValueError(
            f"Count vectors must have equal length ({len(counts_a)} vs {len(counts_b)})"
        )
    numerator = 0.0
    denominator = 0.0
    for a, b in zip(counts_a, counts_b):
        numerator += min(a, b)
        denominator += max(a, b)
    if denominator == 0:
        return 1.0
    return numerator / denominator


def _slot_has_feature(sig: AtomSignature, slot: int) -> bool:
    return sig.feature_bits[slot] == "1"


def _increment_count(counts: list[int], idx: int, count_cap: int | None) -> None:
    if count_cap is None or counts[idx] < count_cap:
        counts[idx] += 1


def silirid_slots() -> list[tuple[str, str]]:
    """Slot order of ``structure_silirid_fingerprint``: (amino_acid, feature_type)."""
    return [
        (aa, FEATURE_TYPES_IN_ORDER[slot])
        for aa in AMINO_ACID_ORDER
        for slot in range(N_FEATURE_TYPES)
    ]


def format_silirid_fingerprint(counts: Sequence[int]) -> str:
    return ",".join(str(int(value)) for value in counts)


def structure_silirid_fingerprint(
    sig_dict: dict[int, AtomSignature],
    *,
    count_cap: int | None = STRUCTURE_COUNT_RESIDUE_DEFAULT_CAP,
) -> tuple[int, ...]:
    counts = [0] * STRUCTURE_SILIRID_LEN
    for sig in sig_dict.values():
        for slot in range(N_FEATURE_TYPES):
            if not _slot_has_feature(sig, slot):
                continue
            for j, char in enumerate(sig.pocket_residue_bits[slot]):
                if char == "1":
                    _increment_count(counts, j * N_FEATURE_TYPES + slot, count_cap)
    return tuple(counts)


def structure_silirid_similarity(
    sig_m0: dict[int, AtomSignature],
    sig_m1: dict[int, AtomSignature],
    *,
    count_cap: int | None = STRUCTURE_COUNT_RESIDUE_DEFAULT_CAP,
) -> float:
    return min_max_similarity(
        structure_silirid_fingerprint(sig_m0, count_cap=count_cap),
        structure_silirid_fingerprint(sig_m1, count_cap=count_cap),
    )
