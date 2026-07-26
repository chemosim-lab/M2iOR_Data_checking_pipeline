# pipeline/scripts/process_molecules.py
import re
from pathlib import Path

import pandas as pd

from scripts.cache_manager import get_cache
from scripts.columns import CAS, CID, INCHIKEY, MIXTURE, MOLECULE, MOLECULE_NAME, SMILES

_CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")
_CID_TOKEN_RE = re.compile(r"^\d+(?:\.\d+)?$")


def parse_cid_cell(raw: str) -> list[int]:
    """Parse a CID cell into its list of CIDs.

    Accepts a single CID or several joined by 'and'/commas (e.g. '123 and 456').
    Free text such as 'mixture of 100 compounds' must be treated as not
    applicable rather than mined for a stray number, so the whole cell is
    rejected unless *every* token is numeric.
    """
    text = raw.strip()
    if not text:
        return []
    normalized = re.sub(r"\band\b", " ", text, flags=re.IGNORECASE).replace(",", " ")
    tokens = normalized.split()
    if not tokens or not all(_CID_TOKEN_RE.match(t) for t in tokens):
        return []
    return [int(float(t)) for t in tokens if int(float(t)) != 0]


def validate_cid_or_cas(df: pd.DataFrame) -> None:
    """Raise ValueError if any Molecule row has neither a usable CID nor a CAS value.

    Exception: a mixture row (Mixture == "mixture") with a Molecule Name and an
    empty CAS is kept even if its CID text isn't a parseable CID (e.g. "mixture
    of 100 compounds") -- the CID text and Name are left untouched, no error raised.
    """
    mol = df[MOLECULE]
    cid_missing = mol[CID].isna() | mol[CID].astype(str).apply(
        lambda v: len(parse_cid_cell(v)) == 0
    )
    cas_missing = mol[CAS].isna() | mol[CAS].astype(str).str.strip().isin(["", "nan"])
    both_missing = cid_missing & cas_missing

    is_mixture = mol[MIXTURE].astype(str).str.strip().str.lower() == "mixture"
    has_name = mol[MOLECULE_NAME].notna() & (
        mol[MOLECULE_NAME].astype(str).str.strip() != ""
    )
    allowed_mixture = both_missing & is_mixture & has_name

    unresolved = both_missing & ~allowed_mixture
    if unresolved.any():
        rows = unresolved[unresolved].index.tolist()
        raise ValueError(
            f"Molecule row(s) with neither CID nor CAS at index(es): {rows}"
        )


def get_unique_cids(df: pd.DataFrame) -> list[int]:
    """Collect all non-null unique CIDs from the Molecule group."""
    if MOLECULE not in df.columns.get_level_values(0):
        raise ValueError(f"Group '{MOLECULE}' not in the DataFrame.")
    col = df[MOLECULE][CID].dropna().astype(str)
    return sorted({cid for raw in col for cid in parse_cid_cell(raw)})


def _extract_prop(props: list[dict], label: str, name: str | None = None) -> str | None:
    for prop in props:
        urn = prop.get("urn", {})
        if urn.get("label") == label:
            if name is None or urn.get("name") == name:
                return prop.get("value", {}).get("sval")
    return None


def _parse_cid_cache(
    unique_cids: list[int],
) -> tuple[
    dict[int, str], dict[int, str], dict[int, str], dict[int, str], dict[int, list[str]]
]:
    """Read cache once per unique CID and return five cid→value dicts."""
    name_map: dict[int, str] = {}
    cas_map: dict[int, str] = {}
    inchikey_map: dict[int, str] = {}
    smiles_map: dict[int, str] = {}
    synonym_map: dict[int, list[str]] = {}

    for cid in unique_cids:
        data = get_cache(str(cid), subdir="molecules")
        if not data:
            continue
        compounds = data.get("PC_Compounds", [])
        if not compounds:
            continue
        props: list[dict] = compounds[0].get("props", [])
        synonyms: list[str] = data.get("synonyms", [])
        if synonyms:
            synonym_map[cid] = synonyms

        ze_name = next(
            (s for s in synonyms if re.search(r"\(Z/E\)", s, re.IGNORECASE)), None
        )
        if name := ze_name or data.get("record_title"):
            name_map[cid] = name
        if key := _extract_prop(props, "InChIKey"):
            inchikey_map[cid] = key
        if smi := _extract_prop(props, "SMILES", "Isomeric") or _extract_prop(
            props, "SMILES", "Absolute"
        ):
            smiles_map[cid] = smi
        if cas := next((s for s in synonyms if _CAS_RE.match(s)), None):
            cas_map[cid] = cas

    return name_map, cas_map, inchikey_map, smiles_map, synonym_map


def _join_values(cids: list[int], mapping: dict[int, str]) -> str | None:
    values = [mapping[c] for c in cids if c in mapping]
    return ", ".join(values) if values else None


def _extract_cid_lists(df: pd.DataFrame) -> tuple[pd.Series, list[int]]:
    raw = df[MOLECULE][CID].dropna().astype(str)
    cid_lists: pd.Series = raw.apply(parse_cid_cell)
    cid_lists = cid_lists[cid_lists.apply(len) > 0]

    unique_cids = sorted({cid for cids in cid_lists for cid in cids})
    return (cid_lists, unique_cids)


def _enrich_molecule_columns(
    df: pd.DataFrame,
    cid_lists: pd.Series,
    name_map: dict[int, str],
    cas_map: dict[int, str],
    inchikey_map: dict[int, str],
    smiles_map: dict[int, str],
) -> None:
    """Replace Molecule Name, CAS, InChIKey and SMILES from PubChem cache.

    Multi-CID cells (e.g. '87839 and 12345') produce comma-joined values.
    """

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


def _export_synonyms_csv(synonym_map: dict[int, list[str]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for cid, synonyms in synonym_map.items():
        csv_path = output_dir / f"{cid}_synonyms.csv"
        with csv_path.open("w", encoding="utf-8") as f:
            f.writelines(s + "\n" for s in synonyms)


def process_molecules(df: pd.DataFrame, output_dir: Path) -> None:
    cid_lists, unique_cids = _extract_cid_lists(df)
    name_map, cas_map, inchikey_map, smiles_map, synonym_map = _parse_cid_cache(
        unique_cids
    )
    _enrich_molecule_columns(df, cid_lists, name_map, cas_map, inchikey_map, smiles_map)
    _export_synonyms_csv(synonym_map, output_dir)
