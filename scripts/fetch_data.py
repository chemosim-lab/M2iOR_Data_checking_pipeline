# pipeline/scripts/fetch_uniprot.py

from multiprocessing import Value
from pandas.core.frame import DataFrame
import time
from pathlib import Path
from typing import Any, cast

import pandas as pd
import requests

from scripts.cache_manager import get_missing_keys, set_cache

_UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
_UNIPROTKB_ENDPOINT_URL = "https://rest.uniprot.org/uniprotkb/{accession}"
_REQUEST_DELAY = 0.2  # seconds between requests to respect UniProt rate limits

_CACHE_FILE = Path(__file__).parent.parent / "cache" / "uniprot_sequences.json"


# def _load_disk_cache() -> dict[str, str | None]:
#     if _CACHE_FILE.exists():
#         with open(_CACHE_FILE, encoding="utf-8") as f:
#             return json.load(f)
#     return {}


# def _save_disk_cache(cache: dict[str, str | None]) -> None:
#     _CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
#     with open(_CACHE_FILE, "w", encoding="utf-8") as f:
#         json.dump(cache, f, indent=2)


def _fetch_data(accession: str) -> Any | None:
    url = _UNIPROTKB_ENDPOINT_URL.format(accession=accession)
    try:
        response = requests.get(url, timeout=10)
        json_data = response.json() if response and response.status_code == 200 else None
        return json_data if json_data else None
    except requests.RequestException as e:
        print(f"  [UniProt] Failed to fetch {accession}: {e}")
        return None


def fetch_uniprot_data(unique_ids: list[str]) -> None:
    

    to_fetch = get_missing_keys(unique_ids)

    if to_fetch:
        print(
            f"Fetching {len(to_fetch)} new accession(s) from UniProt "
            + f"({len(unique_ids) - len(to_fetch)} already cached)..."
        )
        for i, accession in enumerate(to_fetch, start=1):
            print(f"  [{i}/{len(to_fetch)}] {accession}", end=" ... ", flush=True)

            accession_data: Any | None = _fetch_data(accession)

            # Saving in the cache
            if accession_data:
                print("ok")
                set_cache(accession, accession_data)

            if i < len(to_fetch):
                time.sleep(_REQUEST_DELAY)
    else:
        print(f"All {len(unique_ids)} accession(s) found in cache, skipping API calls.")

def fetch_sequences(df: pd.DataFrame) -> pd.DataFrame:
    """
    Query UniProt for each unique accession found in Receptor["UniProt ID"] and
    Co-Receptor["UniProt ID"], then add a "sequence_ref" column right after the
    existing "Sequence" column in each group.

    Accessions already present in the on-disk cache are not re-fetched.
    The cache is persisted to disk after each run that fetches new sequences.

    Returns the augmented DataFrame (same MultiIndex structure).
    """
    groups = ["Receptor", "Co-Receptor"]
    unique_ids = get_unique_uniprot_ids(df, groups)

    cache = _load_disk_cache()
    to_fetch = [acc for acc in unique_ids if acc not in cache]

    if to_fetch:
        print(
            f"Fetching {len(to_fetch)} new accession(s) from UniProt "
            f"({len(unique_ids) - len(to_fetch)} already cached)..."
        )
        for i, accession in enumerate(to_fetch, start=1):
            print(f"  [{i}/{len(to_fetch)}] {accession}", end=" ... ", flush=True)
            seq = _fetch_sequence(accession)
            cache[accession] = seq
            print("ok" if seq else "NOT FOUND")
            if i < len(to_fetch):
                time.sleep(_REQUEST_DELAY)
        _save_disk_cache(cache)
    else:
        print(f"All {len(unique_ids)} accession(s) found in cache, skipping API calls.")

    df = df.copy()

    for group in groups:
        if group not in df.columns.get_level_values(0):
            continue

        if "sequence_ref" in df[group].columns:
            continue  # already added (e.g. called twice on same df)

        uniprot_col = df[group]["UniProt ID"].astype(str).str.strip()

        # "uniprot" when the sequence was successfully fetched, "undefined" otherwise
        database_values = uniprot_col.map(
            lambda acc: "uniprot" if cache.get(acc) is not None else "undefined"
        )
        ref_sequences = uniprot_col.map(lambda acc: cache.get(acc))

        # --- Insert "Database" right after "UniProt ID" ---
        cols = list(df.columns)
        try:
            uid_pos = cols.index((group, "UniProt ID")) + 1
        except ValueError:
            uid_pos = len(cols)
        db_col_idx = pd.MultiIndex.from_tuples([(group, "Database")])
        db_df = pd.DataFrame(
            {(group, "Database"): database_values.values},
            index=df.index,
            columns=db_col_idx,
        )
        df = pd.concat([df.iloc[:, :uid_pos], db_df, df.iloc[:, uid_pos:]], axis=1)

        # --- Insert "sequence_ref" right after "Sequence" ---
        cols = list(df.columns)
        try:
            seq_pos = cols.index((group, "Sequence")) + 1
        except ValueError:
            seq_pos = len(cols)
        ref_col_idx = pd.MultiIndex.from_tuples([(group, "sequence_ref")])
        ref_df = pd.DataFrame(
            {(group, "sequence_ref"): ref_sequences.values},
            index=df.index,
            columns=ref_col_idx,
        )
        df = pd.concat([df.iloc[:, :seq_pos], ref_df, df.iloc[:, seq_pos:]], axis=1)

    return df
