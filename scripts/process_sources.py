# pipeline/scripts/process_sources.py
import re

import pandas as pd

from scripts.cache_manager import get_cache
from scripts.columns import DOI, REFERENCE, SOURCE

_DOI_RE = re.compile(r"^10\.\d{4,9}/.+$", re.IGNORECASE)
_DOI_URL_PREFIX_RE = re.compile(r"^https?://doi\.org/", re.IGNORECASE)


def normalize_doi_column(df: pd.DataFrame) -> None:
    """Strip 'https://doi.org/' prefix from DOI values in place."""
    col = df[SOURCE][DOI].dropna().astype(str).str.strip()
    doi_map = {doi: _DOI_URL_PREFIX_RE.sub("", doi) for doi in col.unique()}
    normalized = col.map(doi_map)
    df.loc[normalized.index, (SOURCE, DOI)] = normalized.to_numpy()


def get_unique_dois(df: pd.DataFrame) -> list[str]:
    """Collect all non-null unique DOIs from the Source group."""
    col: pd.Series = df[SOURCE][DOI].dropna().astype(str).str.strip()
    return sorted(col.unique().tolist())


def enrich_reference_column(df: pd.DataFrame) -> None:
    """Replace Reference column values with APA-formatted references from cache."""
    col = df[SOURCE][DOI].dropna().astype(str).str.strip()
    doi_to_apa = {
        doi: data["apa"]
        for doi in col.unique()
        if (data := get_cache(doi, subdir="sources")) and data.get("apa")
    }
    if not doi_to_apa:
        return
    updated = col.map(doi_to_apa)
    mask = updated.notna()
    if mask.any():
        df.loc[updated.index[mask], (SOURCE, REFERENCE)] = updated[mask].to_numpy()


def validate_doi_column(df: pd.DataFrame) -> None:
    """Raise ValueError if any non-null DOI value does not match the DOI format."""
    col = df[SOURCE][DOI].dropna().astype(str).str.strip()
    invalid = col[~col.apply(lambda v: bool(_DOI_RE.match(v)))]
    if not invalid.empty:
        bad = invalid.unique().tolist()
        raise ValueError(
            f"Column '{DOI}' contains invalid DOI value(s): {bad}"
        )
