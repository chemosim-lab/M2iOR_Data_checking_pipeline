# pipeline/scripts/process_molecules.py
import re
from pathlib import Path
from typing import Any

import colorlog
import pandas as pd

from scripts.cache_manager import get_cache
from scripts.columns import CAS, CID, INCHIKEY, MIXTURE, MOLECULE, MOLECULE_NAME, SMILES
from scripts.molecule_stereo import export_stereo_classification

logger = colorlog.getLogger(__name__)

_CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")
_CID_TOKEN_RE = re.compile(r"^\d+(?:\.\d+)?$")

_PUBCHEM_COMPOUND_URL = "https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
_CAS_COMMON_CHEMISTRY_DETAIL_URL = "https://commonchemistry.cas.org/detail?cas_rn={cas}"

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

# Chemical nomenclature (e.g. stereodescriptors like "(−)-") and CAS numbers
# (e.g. "2244–16-8") often use a typographic minus/dash instead of a plain
# hyphen-minus, which PubChem's own names/synonyms - and the CAS format
# itself - always use.
_DASH_TO_HYPHEN = {
    "‐": "-",  # HYPHEN
    "‑": "-",  # NON-BREAKING HYPHEN
    "‒": "-",  # FIGURE DASH
    "–": "-",  # EN DASH
    "—": "-",  # EM DASH
    "−": "-",  # MINUS SIGN
}


def normalize_dashes(text: str) -> str:
    """Replace typographic dash/minus characters with a plain hyphen-minus."""
    return "".join(_DASH_TO_HYPHEN.get(ch, ch) for ch in text)


# Locants like "4'-Ethylacetophenone" are commonly typed with a typographic
# prime or curly quote instead of a plain apostrophe.
_QUOTE_TO_APOSTROPHE = {
    "′": "'",  # PRIME
    "‘": "'",  # LEFT SINGLE QUOTATION MARK
    "’": "'",  # RIGHT SINGLE QUOTATION MARK
}

# Racemic mixtures are marked "(±)-" in the input table but PubChem/CAS spell
# it out as "+-" (e.g. "(±)-fenchone" vs "(+-)-Fenchone").
_PLUSMINUS_TO_ASCII = {
    "±": "+-",
}

_NAME_CHAR_REPLACEMENTS = {
    **_GREEK_TO_LATIN,
    **_DASH_TO_HYPHEN,
    **_QUOTE_TO_APOSTROPHE,
    **_PLUSMINUS_TO_ASCII,
}


_HTML_TAG_RE = re.compile(r"<[^>]+>")

# "(E)-"/"(Z)-" and "trans-"/"cis-" are interchangeable stereodescriptor
# conventions for the simple disubstituted alkenes these odorant names cover
# (e.g. "(E)-β-Farnesene" == PubChem's own "trans-beta-Farnesene"); collapsing
# both to the same form lets them compare equal.
_STEREO_PAREN_RE = re.compile(r"\(([ez])\)-?")


def _normalize_name(text: str) -> str:
    """Lowercase a molecule name, strip HTML markup, spell out Greek letters
    (e.g. 'β' -> 'beta'), normalize typographic dashes/quotes/± to their plain
    ASCII form, and canonicalize E/Z vs. trans/cis stereodescriptors.

    CAS Common Chemistry's synonyms wrap stereodescriptors in formatting tags
    (e.g. '(<em>R</em>)-(+)-Limonene', '<span class="text-smallcaps">D</span>
    -Limonene') that must be stripped for these to compare equal to a plain
    '(R)-(+)-Limonene' from the input table. PubChem's own names/synonyms and
    typographic dashes (e.g. '(−)-β-Elemene' vs '(-)-beta-Elemene') are
    normalized the same way."""
    lowered = _HTML_TAG_RE.sub("", text.strip().lower())
    replaced = "".join(_NAME_CHAR_REPLACEMENTS.get(ch, ch) for ch in lowered)
    replaced = _STEREO_PAREN_RE.sub(r"\1-", replaced)
    return replaced.replace("trans-", "e-").replace("cis-", "z-")


def _normalize_name_variants(text: str) -> set[str]:
    """`_normalize_name`, plus an alternate racemic-marker spelling when one
    is present. PubChem/CAS spell it inconsistently even for the same
    compound (e.g. "(+-)-Linalool" next to "(+/-)-linalool" in the same
    synonym list, while "(+/-)-Citronellal" has no "+-" form at all), so a
    name normalizing to one has to be checked against both."""
    normalized = _normalize_name(text)
    variants = {normalized}
    if "+/-" in normalized:
        variants.add(normalized.replace("+/-", "+-"))
    elif "+-" in normalized:
        variants.add(normalized.replace("+-", "+/-"))
    return variants


def parse_cid_cell(raw: object) -> list[int]:
    """Parse a CID cell into its list of CIDs.

    Accepts a single CID or several joined by 'and'/commas/'+' (e.g. '123 and 456',
    '123+456'). Free text such as 'mixture of 100 compounds' must be treated as not
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
    normalized = re.sub(r"\band\b", " ", text, flags=re.IGNORECASE).replace(
        ",", " "
    ).replace("+", " ")
    tokens = normalized.split()
    if not tokens or not all(_CID_TOKEN_RE.match(t) for t in tokens):
        return []
    return [int(float(t)) for t in tokens if int(float(t)) != 0]


def parse_cas_cell(raw: object) -> list[str]:
    """Parse a CAS cell into its list of CAS numbers.

    Mirrors `parse_cid_cell`: accepts a single CAS number or several joined by
    'and'/commas/'+'. Tokens that don't look like a CAS number (##-##-#) are
    dropped. Typographic dashes (e.g. "2244–16-8") are normalized to a plain
    hyphen first, since a CAS number is only ever typed with one but sources
    sometimes substitute an en dash or similar.
    """
    if pd.isna(raw):
        return []
    text = normalize_dashes(str(raw).strip())
    if not text:
        return []
    normalized = re.sub(r"\band\b", ",", text, flags=re.IGNORECASE).replace("+", ",")
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


def _cas_side_candidates(
    cas: str,
    rcid: int,
    cas_details_map: dict[str, dict[str, Any]],
    name_map: dict[int, str],
    synonym_map: dict[int, list[str]],
) -> set[str]:
    """Name/synonym candidates for one CAS number's identity.

    Prioritizes CAS Common Chemistry's own record for that CAS number (see
    `fetch_cas_common_chemistry_details`) - CAS itself, over PubChem's
    name/synonyms for whatever CID PubChem happens to resolve that CAS to.
    Falls back to the latter when CAS Common Chemistry has no record for it.
    """
    detail = cas_details_map.get(cas)
    if detail:
        candidates: set[str] = set()
        if detail.get("name"):
            candidates.add(_normalize_name(detail["name"]))
        candidates.update(_normalize_name(s) for s in detail.get("synonyms") or [])
        if candidates:
            return candidates
    return _pooled_candidates(rcid, name_map, synonym_map)


def reconcile_cid_cas(
    df: pd.DataFrame,
    cas_to_cid_map: dict[str, int],
    name_map: dict[int, str],
    synonym_map: dict[int, list[str]],
    cas_map: dict[int, str] | None = None,
    cas_details_map: dict[str, dict[str, Any]] | None = None,
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
      - Name matches only the CAS side's name/synonyms (priority: CAS Common
        Chemistry's own record for that CAS number, see `_cas_side_candidates`)
        -> the CID is wrong; it is corrected here so downstream enrichment
        re-derives Name/CAS/InChIKey/SMILES from the right compound.
      - Name matches both -> no real ambiguity, keep the CID as-is silently.
      - Name matches neither -> unresolvable, blocks the pipeline.

    Mixture rows (several CIDs in one cell) are skipped -- a CAS number can't
    be reliably pinned to one compound among several.
    """
    cas_map = cas_map or {}
    cas_details_map = cas_details_map or {}
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
        per_cas_candidates = {
            cas: _cas_side_candidates(cas, rcid, cas_details_map, name_map, synonym_map)
            for cas, rcid in inconsistent.items()
        }
        cas_candidates: set[str] = set().union(*per_cas_candidates.values())

        match_cid = bool(name) and name in cid_candidates
        match_cas = bool(name) and name in cas_candidates

        if match_cid and match_cas:
            continue  # Name confirms both sides - no ambiguity to flag
        if match_cid:
            corrected_cas = cas_map.get(cid)
            logger.warning(
                "  Molecule row %s: CAS(es) %s don't match the table's CID:%s - "
                "keeping CID:%s (Name %r matches it) and correcting the CAS to %s.",
                idx,
                mismatch_desc,
                cid,
                cid,
                raw_name,
                corrected_cas or "<no CAS found in PubChem's record for this CID>",
            )
        elif match_cas:
            new_cas, new_cid = next(
                (cas, inconsistent[cas])
                for cas, candidates in per_cas_candidates.items()
                if name in candidates
            )
            logger.warning(
                "  Molecule row %s: CID:%s doesn't match CAS:%s - Name %r matches "
                "the CAS side instead, correcting CID:%s → CID:%s.",
                idx,
                cid,
                new_cas,
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


def _reference_urls(cids: list[int], cas_list: list[str]) -> str:
    """Format a row's CID(s)/CAS(es) as clickable reference links, e.g.
    "CID 537747 (https://pubchem.ncbi.nlm.nih.gov/compound/537747)". Shows
    whichever of the two the row actually has - CID only, CAS only, or both."""
    refs = [
        f"CID {cid} ({_PUBCHEM_COMPOUND_URL.format(cid=cid)})" for cid in cids
    ]
    refs += [
        f"CAS {cas} ({_CAS_COMMON_CHEMISTRY_DETAIL_URL.format(cas=cas)})"
        for cas in cas_list
    ]
    return ", ".join(refs) if refs else "no CID/CAS"


def _row_name_candidates(
    cids: list[int],
    cas_list: list[str],
    cas_details_map: dict[str, dict[str, Any]],
    name_map: dict[int, str],
    synonym_map: dict[int, list[str]],
) -> tuple[set[str], list[set[str]]]:
    """Build the pooled name/synonym candidates for a row, plus each CID's own
    candidates (for the positional mixture check). CAS Common Chemistry's
    name/synonyms are folded in too, but only for single-CID rows - it's
    ambiguous which CAS maps to which CID in a mixture row."""
    per_cid_candidates = [
        _pooled_candidates(cid, name_map, synonym_map) for cid in cids
    ]
    pooled: set[str] = set()
    for candidates in per_cid_candidates:
        pooled |= candidates

    if len(cids) == 1:
        for cas in cas_list:
            pooled |= _cas_side_candidates(
                cas, cids[0], cas_details_map, name_map, synonym_map
            )

    return pooled, per_cid_candidates


def validate_molecule_name_column(
    df: pd.DataFrame,
    cid_lists: pd.Series,
    name_map: dict[int, str],
    synonym_map: dict[int, list[str]],
    cas_details_map: dict[str, dict[str, Any]] | None = None,
) -> None:
    """Raise ValueError if a row's Molecule Name isn't among its CID(s)' PubChem
    name/synonyms, or its CAS number(s)' CAS Common Chemistry name/synonyms.

    A multi-CID row's Name is split the same way as its CID cell (on 'and'/commas).
    When the split produces as many name parts as CIDs, each part is checked
    against its own corresponding CID's name/synonyms (positional match) --
    this catches e.g. two names swapped between the CIDs of a 2-CID mixture row,
    which a pooled check would miss since both names would still be present
    somewhere in the union. A CID with no cached name/synonyms is skipped at its
    position rather than treated as a mismatch. Otherwise (a part/CID count
    mismatch), every part is checked against the pooled {record title, synonyms}
    of all the row's CIDs instead. A single-CID row's Name is never split --
    chemical nomenclature routinely contains commas that aren't mixture
    separators (e.g. "1,8-cineole") -- and is checked whole against that one
    CID's name/synonyms.

    For a single-CID row with a CAS number, CAS Common Chemistry's own
    name/synonyms for that CAS (see `fetch_cas_common_chemistry_details`) are
    added to the pool alongside PubChem's: the two registries don't always
    spell a name the same way (e.g. PubChem's "trans-beta-Farnesene" vs CAS's
    "(E)-β-Farnesene"), so either source confirming the name is accepted
    rather than requiring PubChem specifically to.

    Rows without a Molecule Name, or with nothing cached to compare against
    from either source, are skipped.
    """
    cas_details_map = cas_details_map or {}
    name_col = df[MOLECULE][MOLECULE_NAME]
    cas_col = df[MOLECULE][CAS]
    mismatches: list[tuple[int, str, list[int], list[str]]] = []

    for idx, cids in cid_lists.items():
        raw_name = name_col.get(idx)
        if pd.isna(raw_name) or not str(raw_name).strip():
            continue

        cas_list = parse_cas_cell(cas_col.get(idx))
        pooled, per_cid_candidates = _row_name_candidates(
            cids, cas_list, cas_details_map, name_map, synonym_map
        )
        if not pooled:
            continue  # nothing cached to compare against

        # Only split on 'and'/commas for actual multi-CID mixture rows - a
        # single-CID name is compared whole, since chemical nomenclature
        # routinely contains commas that aren't mixture separators (e.g.
        # "1,8-cineole", "2,3-butanedione", "2,4,5-trimethylthiazole").
        if len(cids) > 1:
            parts = [
                p.strip()
                for p in re.sub(
                    r"\band\b", ",", str(raw_name), flags=re.IGNORECASE
                ).split(",")
                if p.strip()
            ]
        else:
            parts = [str(raw_name).strip()]

        if len(parts) == len(cids) and len(cids) > 1:
            mismatch = any(
                candidates and not (_normalize_name_variants(part) & candidates)
                for part, candidates in zip(parts, per_cid_candidates, strict=True)
            )
        else:
            mismatch = any(
                not (_normalize_name_variants(part) & pooled) for part in parts
            )

        if mismatch:
            mismatches.append((idx, str(raw_name), cids, cas_list))

    if mismatches:
        lines = "\n".join(
            f"  - row {idx}: {name!r} — {_reference_urls(cids, cas_list)}"
            for idx, name, cids, cas_list in mismatches
        )
        msg = (
            f"Column '{MOLECULE_NAME}' has value(s) not found among the "
            f"corresponding CID's PubChem name/synonyms:\n{lines}"
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
    cas_details_map: dict[str, dict[str, Any]] | None = None,
) -> None:
    cas_to_cid_map = cas_to_cid_map or {}
    _, initial_unique_cids = _extract_cid_lists(df)
    # Also cache CAS-resolved CIDs that disagree with their row's own CID, so
    # there's PubChem name/synonym data to compare the Molecule Name against.
    cache_cids = sorted(set(initial_unique_cids) | set(cas_to_cid_map.values()))
    name_map, cas_map, inchikey_map, smiles_map, synonym_map = _parse_cid_cache(
        cache_cids
    )

    reconcile_cid_cas(
        df,
        cas_to_cid_map,
        name_map,
        synonym_map,
        cas_map=cas_map,
        cas_details_map=cas_details_map,
    )

    # Re-extract: reconcile_cid_cas may have corrected some rows' CID.
    cid_lists, unique_cids = _extract_cid_lists(df)
    validate_molecule_name_column(
        df, cid_lists, name_map, synonym_map, cas_details_map=cas_details_map
    )
    _enrich_molecule_columns(df, cid_lists, name_map, cas_map, inchikey_map, smiles_map)
    _export_synonyms_csv(synonym_map, output_dir)

    mixture_map = export_stereo_classification(unique_cids, smiles_map, images_dir)
    _enrich_mixture_column(df, cid_lists, mixture_map)
