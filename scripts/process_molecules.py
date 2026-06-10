# pipeline/scripts/process_molecules.py
import re

import pandas as pd

from scripts.cache_manager import get_cache
from scripts.columns import CAS, CID, INCHIKEY, MOLECULE, MOLECULE_NAME, SMILES

_CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")
_CID_NUM_RE = re.compile(r"^\d+(?:\.\d+)?$")



def get_unique_cids(df: pd.DataFrame) -> list[int]:
    """Collect all non-null unique CIDs from the Molecule group."""
    if MOLECULE not in df.columns.get_level_values(0):
        raise ValueError(f"Group '{MOLECULE}' not in the DataFrame.")
    col = df[MOLECULE][CID].dropna().astype(str)
    expanded = (
        col.str.replace(r"\band\b", " ", regex=True, case=False)
        .str.split()
        .explode()
    )
    valid = expanded[expanded.str.match(r"^\d+(?:\.\d+)?$")]
    return sorted({int(float(v)) for v in valid} - {0})


def _extract_prop(props: list[dict], label: str, name: str | None = None) -> str | None:
    for prop in props:
        urn = prop.get("urn", {})
        if urn.get("label") == label:
            if name is None or urn.get("name") == name:
                return prop.get("value", {}).get("sval")
    return None


def _build_column_maps(
    unique_cids: list[int],
) -> tuple[dict[int, str], dict[int, str], dict[int, str], dict[int, str]]:
    """Read cache once per unique CID and return four cid→value dicts."""
    name_map: dict[int, str] = {}
    cas_map: dict[int, str] = {}
    inchikey_map: dict[int, str] = {}
    smiles_map: dict[int, str] = {}

    for cid in unique_cids:
        data = get_cache(str(cid), subdir="molecules")
        if not data:
            continue
        compounds = data.get("PC_Compounds", [])
        if not compounds:
            continue
        props: list[dict] = compounds[0].get("props", [])
        synonyms: list[str] = data.get("synonyms", [])

        if name := data.get("record_title"):
            name_map[cid] = name
        if key := _extract_prop(props, "InChIKey"):
            inchikey_map[cid] = key
        if smi := _extract_prop(props, "SMILES", "Isomeric") or _extract_prop(
            props, "SMILES", "Absolute"
        ):
            smiles_map[cid] = smi
        if cas := next((s for s in synonyms if _CAS_RE.match(s)), None):
            cas_map[cid] = cas

    return name_map, cas_map, inchikey_map, smiles_map


def _tokens_to_cids(tokens: list[str]) -> list[int]:
    return [
        int(float(t)) for t in tokens
        if _CID_NUM_RE.match(t) and int(float(t)) != 0
    ]


def _join_values(cids: list[int], mapping: dict[int, str]) -> str | None:
    values = [mapping[c] for c in cids if c in mapping]
    return ", ".join(values) if values else None


def enrich_molecule_columns(df: pd.DataFrame) -> None:
    """Replace Molecule Name, CAS, InChIKey and SMILES from PubChem cache.

    Multi-CID cells (e.g. '87839 and 12345') produce comma-joined values.
    """
    raw = df[MOLECULE][CID].dropna().astype(str)
    cid_lists: pd.Series = (
        raw.str.replace(r"\band\b", " ", regex=True, case=False)
        .str.split()
        .apply(_tokens_to_cids)
    )
    cid_lists = cid_lists[cid_lists.apply(len) > 0]

    unique_cids = sorted({cid for cids in cid_lists for cid in cids})
    name_map, cas_map, inchikey_map, smiles_map = _build_column_maps(unique_cids)

    for col, mapping in [
        (MOLECULE_NAME, name_map),
        (CAS, cas_map),
        (INCHIKEY, inchikey_map),
        (SMILES, smiles_map),
    ]:
        updated = cid_lists.apply(_join_values, mapping=mapping)
        mask = updated.notna()
        if mask.any():
            df.loc[updated.index[mask], (MOLECULE, col)] = updated[mask].to_numpy()
