# pipeline/scripts/process_molecules.py
import re

import pandas as pd

from scripts.cache_manager import get_cache
from scripts.columns import CAS, CID, INCHIKEY, MOLECULE, MOLECULE_NAME, SMILES

_CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")


def get_unique_cids(df: pd.DataFrame) -> list[int]:
    """Collect all non-null unique CIDs from the Molecule group."""
    if MOLECULE not in df.columns.get_level_values(0):
        raise ValueError(f"Group '{MOLECULE}' not in the DataFrame.")
    col = df[MOLECULE][CID]
    return sorted(col.dropna().astype(int).unique().tolist())


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


def enrich_molecule_columns(df: pd.DataFrame) -> None:
    """Replace Molecule Name, CAS, InChIKey and SMILES from PubChem cache."""
    cid_series = df[MOLECULE][CID].dropna().astype(int)
    unique_cids = cid_series.unique().tolist()

    name_map, cas_map, inchikey_map, smiles_map = _build_column_maps(unique_cids)

    for col, mapping in [
        (MOLECULE_NAME, name_map),
        (CAS, cas_map),
        (INCHIKEY, inchikey_map),
        (SMILES, smiles_map),
    ]:
        updated = cid_series.map(mapping)
        mask = updated.notna()
        if mask.any():
            df.loc[updated.index[mask], (MOLECULE, col)] = updated[mask].to_numpy()
