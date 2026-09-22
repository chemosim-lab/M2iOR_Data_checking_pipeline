# pipeline/scripts/process_molecules.py
import re
from pathlib import Path

import colorlog
import pandas as pd

from scripts.cache_manager import get_cache
from scripts.columns import CAS, CID, INCHIKEY, MIXTURE, MOLECULE, MOLECULE_NAME, SMILES
from scripts.molecule_stereo import export_stereo_classification

logger = colorlog.getLogger(__name__)

_CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")
_CID_TOKEN_RE = re.compile(r"^\d+(?:\.\d+)?$")

_GREEK_TO_LATIN = {
    "α": "alpha",
    "β": "beta",
    "γ": "gamma",
    "δ": "delta",
    "ε": "epsilon",
    "ζ": "zeta",
    "η": "eta",
    "θ": "theta",
    "ι": "iota",
    "κ": "kappa",
    "λ": "lambda",
    "μ": "mu",
    "ν": "nu",
    "ξ": "xi",
    "ο": "omicron",
    "π": "pi",
    "ρ": "rho",
    "σ": "sigma",
    "ς": "sigma",
    "τ": "tau",
    "υ": "upsilon",
    "φ": "phi",
    "χ": "chi",
    "ψ": "psi",
    "ω": "omega",
}

# Chemical nomenclature (e.g. stereodescriptors like "(−)-") often uses a
# typographic minus/dash instead of a plain hyphen-minus, which PubChem's own
# names/synonyms always use.
_DASH_TO_HYPHEN = {
    "‐": "-",  # HYPHEN
    "‑": "-",  # NON-BREAKING HYPHEN
    "‒": "-",  # FIGURE DASH
    "–": "-",  # EN DASH
    "—": "-",  # EM DASH
    "−": "-",  # MINUS SIGN
}

_NAME_CHAR_REPLACEMENTS = {**_GREEK_TO_LATIN, **_DASH_TO_HYPHEN}


def _normalize_name(text: str) -> str:
    """Lowercase a molecule name, spell out Greek letters (e.g. 'β' -> 'beta')
    and normalize typographic dashes to a plain hyphen, so it compares equal to
    PubChem's Latin-spelled, plain-hyphen equivalent (e.g. '(−)-β-Elemene' vs
    '(-)-beta-Elemene')."""
    lowered = text.strip().lower()
    return "".join(_NAME_CHAR_REPLACEMENTS.get(ch, ch) for ch in lowered)


def parse_cid_cell(raw: object) -> list[int]:
    """Parse a CID cell into its list of CIDs.

    Accepts a single CID or several joined by 'and'/commas (e.g. '123 and 456').
    Free text such as 'mixture of 100 compounds' must be treated as not
    applicable rather than mined for a stray number, so the whole cell is
    rejected unless *every* token is numeric.

    `raw` may be a non-string scalar (e.g. NaN) since callers sometimes pass
    a cell through `.astype(str)` where missing values aren't guaranteed to
    become the literal string "nan" -- such values are treated as empty.
    """
    if pd.isna(raw):
        return []
    text = str(raw).strip()
    if not text:
        return []
    normalized = re.sub(r"\band\b", " ", text, flags=re.IGNORECASE).replace(",", " ")
    tokens = normalized.split()
    if not tokens or not all(_CID_TOKEN_RE.match(t) for t in tokens):
        return []
    return [int(float(t)) for t in tokens if int(float(t)) != 0]


def parse_cas_cell(raw: object) -> list[str]:
    """Parse a CAS cell into its list of CAS numbers.

    Mirrors `parse_cid_cell`: accepts a single CAS number or several joined by
    'and'/commas. Tokens that don't look like a CAS number (##-##-#) are dropped.
    """
    if pd.isna(raw):
        return []
    text = str(raw).strip()
    if not text:
        return []
    normalized = re.sub(r"\band\b", ",", text, flags=re.IGNORECASE)
    tokens = [t.strip() for t in normalized.split(",") if t.strip()]
    return [t for t in tokens if _CAS_RE.match(t)]


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


def get_unique_cas_for_cross_check(df: pd.DataFrame) -> list[str]:
    """Collect unique CAS numbers worth resolving on PubChem to cross-check
    against the row's CID.

    Only single-CID rows are considered: a mixture row's CID cell lists several
    distinct compounds, so a CAS number on that row can't be pinned to one of
    them reliably.
    """
    mol = df[MOLECULE]
    cas_values: set[str] = set()
    for idx in mol.index:
        if len(parse_cid_cell(mol.loc[idx, CID])) != 1:
            continue
        cas_values.update(parse_cas_cell(mol.loc[idx, CAS]))
    return sorted(cas_values)


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

    return name_map, cas_map, inchikey_map, smiles_map, synonym_map


def _ensure_object_dtype(df: pd.DataFrame, column: tuple[str, str]) -> None:
    """Widen an all-null (float64) column so it can hold string values.

    A column left entirely empty in the source Excel is read as float64;
    assigning strings into it then raises pandas.errors.LossySetitemError.
    """
    group, col = column
    if df[group][col].dtype != object:
        df[column] = df[column].astype(object)


def _join_values(cids: list[int], mapping: dict[int, str]) -> str | None:
    values = [mapping[c] for c in cids if c in mapping]
    return ", ".join(values) if values else None


def _aggregate_mixture(cids: list[int], mapping: dict[int, str]) -> str | None:
    """Derive the row-level Mixture label for a (possibly multi-CID) row.

    A row naming more than one CID explicitly lists distinct compounds and is
    always "mixture", regardless of each component's own stereochemistry. A
    single-CID row falls back to that compound's stereo classification
    ("sum of isomers" if its stereochemistry is undefined, "monomolecular"
    otherwise). Returns None only when a single-CID row's compound could not
    be classified (e.g. no cached SMILES).
    """
    if len(cids) > 1:
        return "mixture"
    if not cids or cids[0] not in mapping:
        return None
    return mapping[cids[0]]


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
            _ensure_object_dtype(df, (MOLECULE, col))
            df.loc[updated.index[mask], (MOLECULE, col)] = updated[mask].to_numpy()


def _enrich_mixture_column(
    df: pd.DataFrame, cid_lists: pd.Series, mixture_map: dict[int, str]
) -> None:
    """Set Mixture from the CID count and SMILES-based stereo classification.

    Rows explicitly curated as "mixture" (an actual mix of several named
    compounds, see `validate_cid_or_cas`) are left untouched.
    """
    is_explicit_mixture = (
        df[MOLECULE][MIXTURE].astype(str).str.strip().str.lower() == "mixture"
    )
    computed = cid_lists.apply(_aggregate_mixture, mapping=mixture_map)
    mask = computed.notna() & ~is_explicit_mixture.loc[computed.index]
    if mask.any():
        _ensure_object_dtype(df, (MOLECULE, MIXTURE))
        df.loc[computed.index[mask], (MOLECULE, MIXTURE)] = computed[mask].to_numpy()


def _pooled_candidates(
    cid: int, name_map: dict[int, str], synonym_map: dict[int, list[str]]
) -> set[str]:
    candidates: set[str] = set()
    if cid in name_map:
        candidates.add(_normalize_name(name_map[cid]))
    candidates.update(_normalize_name(s) for s in synonym_map.get(cid, []))
    return candidates


def reconcile_cid_cas(
    df: pd.DataFrame,
    cas_to_cid_map: dict[str, int],
    name_map: dict[int, str],
    synonym_map: dict[int, list[str]],
) -> None:
    """For single-CID rows with one or more CAS numbers, verify PubChem agrees
    that the CAS(es) and the CID refer to the same compound; correct whichever
    one the Molecule Name disagrees with.

    A CAS resolving (via `cas_to_cid_map`) to a different CID than the row's own
    is a conflict. The row's Molecule Name then decides which identifier is
    trustworthy:
      - Name matches only the CID's PubChem name/synonyms -> the CAS is wrong;
        it gets corrected downstream by `_enrich_molecule_columns` from the
        CID's own record, so nothing needs fixing here besides a warning.
      - Name matches only the CAS-resolved compound's name/synonyms -> the CID
        is wrong; it is corrected here so downstream enrichment re-derives
        Name/CAS/InChIKey/SMILES from the right compound.
      - Name matches both -> no real ambiguity, keep the CID as-is silently.
      - Name matches neither -> unresolvable, blocks the pipeline.

    Mixture rows (several CIDs in one cell) are skipped -- a CAS number can't
    be reliably pinned to one compound among several.
    """
    mol = df[MOLECULE]
    blocked: list[tuple[int, str]] = []

    for idx in mol.index:
        cids = parse_cid_cell(mol.loc[idx, CID])
        if len(cids) != 1:
            continue
        cid = cids[0]
        cas_list = parse_cas_cell(mol.loc[idx, CAS])
        if not cas_list:
            continue

        inconsistent = {
            cas: rcid
            for cas in cas_list
            if (rcid := cas_to_cid_map.get(cas)) is not None and rcid != cid
        }
        if not inconsistent:
            continue  # every resolvable CAS agrees with the table's CID

        mismatch_desc = ", ".join(
            f"{cas}→CID:{rcid}" for cas, rcid in inconsistent.items()
        )
        raw_name = mol.loc[idx, MOLECULE_NAME]
        name = _normalize_name(str(raw_name)) if pd.notna(raw_name) else ""

        cid_candidates = _pooled_candidates(cid, name_map, synonym_map)
        cas_cids = sorted(set(inconsistent.values()))
        cas_candidates: set[str] = set()
        for rcid in cas_cids:
            cas_candidates |= _pooled_candidates(rcid, name_map, synonym_map)

        match_cid = bool(name) and name in cid_candidates
        match_cas = bool(name) and name in cas_candidates

        if match_cid and match_cas:
            continue  # Name confirms both sides - no ambiguity to flag
        if match_cid:
            logger.warning(
                "  Molecule row %s: CAS(es) %s don't match the table's CID:%s - "
                "keeping CID:%s (Name %r matches it) and correcting the CAS.",
                idx,
                mismatch_desc,
                cid,
                cid,
                raw_name,
            )
        elif match_cas:
            new_cid = next(
                rcid
                for rcid in cas_cids
                if name in _pooled_candidates(rcid, name_map, synonym_map)
            )
            logger.warning(
                "  Molecule row %s: CID:%s doesn't match CAS(es) %s - Name %r matches "
                "the CAS side instead, correcting CID:%s → CID:%s.",
                idx,
                cid,
                mismatch_desc,
                raw_name,
                cid,
                new_cid,
            )
            # The CID column's dtype (string, int64...) depends on what the
            # source Excel held; widen to object first so it accepts the
            # corrected value regardless, stored as text like every other CID.
            _ensure_object_dtype(df, (MOLECULE, CID))
            df.loc[idx, (MOLECULE, CID)] = str(new_cid)
        else:
            blocked.append((idx, str(raw_name)))

    if blocked:
        rows = ", ".join(f"{idx} ({name!r})" for idx, name in blocked)
        msg = (
            "Molecule row(s) with a CID/CAS mismatch whose Name matches neither "
            f"reference's PubChem name/synonyms: {rows}"
        )
        raise ValueError(msg)


def validate_molecule_name_column(
    df: pd.DataFrame,
    cid_lists: pd.Series,
    name_map: dict[int, str],
    synonym_map: dict[int, list[str]],
) -> None:
    """Raise ValueError if a row's Molecule Name isn't among its CID(s)' PubChem
    name or synonyms.

    A multi-CID row's Name is split the same way as its CID cell (on 'and'/commas).
    When the split produces as many name parts as CIDs, each part is checked
    against its own corresponding CID's name/synonyms (positional match) --
    this catches e.g. two names swapped between the CIDs of a 2-CID mixture row,
    which a pooled check would miss since both names would still be present
    somewhere in the union. A CID with no cached name/synonyms is skipped at its
    position rather than treated as a mismatch. Otherwise (single CID, or a
    part/CID count mismatch), every part is checked against the pooled
    {record title, synonyms} of all the row's CIDs instead. Rows without a
    Molecule Name, or whose CID(s) have no cached name/synonyms to compare
    against at all, are skipped.
    """
    name_col = df[MOLECULE][MOLECULE_NAME]
    mismatches: list[tuple[int, str]] = []

    for idx, cids in cid_lists.items():
        raw_name = name_col.get(idx)
        if pd.isna(raw_name) or not str(raw_name).strip():
            continue

        per_cid_candidates = [
            _pooled_candidates(cid, name_map, synonym_map) for cid in cids
        ]
        pooled = set().union(*per_cid_candidates) if per_cid_candidates else set()
        if not pooled:
            continue  # nothing cached to compare against

        parts = [
            p.strip()
            for p in re.sub(r"\band\b", ",", str(raw_name), flags=re.IGNORECASE).split(
                ","
            )
            if p.strip()
        ]

        if len(parts) == len(cids) and len(cids) > 1:
            mismatch = any(
                candidates and _normalize_name(part) not in candidates
                for part, candidates in zip(parts, per_cid_candidates, strict=True)
            )
        else:
            mismatch = any(_normalize_name(part) not in pooled for part in parts)

        if mismatch:
            mismatches.append((idx, str(raw_name)))

    if mismatches:
        rows = ", ".join(f"{idx} ({name!r})" for idx, name in mismatches)
        msg = (
            f"Column '{MOLECULE_NAME}' has value(s) not found among the "
            f"corresponding CID's PubChem name/synonyms at row(s): {rows}"
        )
        raise ValueError(msg)


def _export_synonyms_csv(synonym_map: dict[int, list[str]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for cid, synonyms in synonym_map.items():
        csv_path = output_dir / f"{cid}_synonyms.csv"
        with csv_path.open("w", encoding="utf-8") as f:
            f.writelines(s + "\n" for s in synonyms)


def process_molecules(
    df: pd.DataFrame,
    output_dir: Path,
    images_dir: Path,
    cas_to_cid_map: dict[str, int] | None = None,
) -> None:
    cas_to_cid_map = cas_to_cid_map or {}
    _, initial_unique_cids = _extract_cid_lists(df)
    # Also cache CAS-resolved CIDs that disagree with their row's own CID, so
    # there's PubChem name/synonym data to compare the Molecule Name against.
    cache_cids = sorted(set(initial_unique_cids) | set(cas_to_cid_map.values()))
    name_map, cas_map, inchikey_map, smiles_map, synonym_map = _parse_cid_cache(
        cache_cids
    )

    reconcile_cid_cas(df, cas_to_cid_map, name_map, synonym_map)

    # Re-extract: reconcile_cid_cas may have corrected some rows' CID.
    cid_lists, unique_cids = _extract_cid_lists(df)
    validate_molecule_name_column(df, cid_lists, name_map, synonym_map)
    _enrich_molecule_columns(df, cid_lists, name_map, cas_map, inchikey_map, smiles_map)
    _export_synonyms_csv(synonym_map, output_dir)

    mixture_map = export_stereo_classification(unique_cids, smiles_map, images_dir)
    _enrich_mixture_column(df, cid_lists, mixture_map)
