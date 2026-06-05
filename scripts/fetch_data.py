# pipeline/scripts/fetch_uniprot.py

from __future__ import annotations

import time
from multiprocessing import Value
from pathlib import Path
from typing import TYPE_CHECKING, Any

import requests

if TYPE_CHECKING:
    import pandas as pd

from scripts.cache_manager import get_missing_keys, set_cache
from scripts.fetch_blast import (
    clear_blast_cache_entries,
    count_uncached_blast_queries,
    get_failed_blast_queries,
)
from scripts.process_receptors import resolve_missing_ids_via_blast

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
        json_data = (
            response.json() if response and response.status_code == 200 else None
        )
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


def fetch_genbank_data(
    df: pd.DataFrame,
    protein_groups: list[str],
    blast_queries: list[tuple[str, str]],
) -> dict[str, str]:
    if not blast_queries:
        return {}

    new_count = count_uncached_blast_queries(blast_queries)
    failed = get_failed_blast_queries(blast_queries)

    if new_count > 0 or failed:
        print(f"\n  ! {len(blast_queries)} sequence(s) have no UniProt ID:")  # noqa: T201
        if new_count > 0:
            print(f"    {new_count} new (never queried, ~1-5 min each)")  # noqa: T201
        for _seq, sp, err in failed:
            print(f"    previously failed — {sp}: {err}")  # noqa: T201

        total = new_count + len(failed)
        answer = input(f"  Run/retry {total} BLAST lookup(s)? [y/N] ").strip().lower()
        if answer != "y":
            msg = f"Aborted by user: {total} BLAST lookup(s) required."
            raise ValueError(msg)  # noqa: TRY301

        if failed:
            clear_blast_cache_entries([(s, sp) for s, sp, _ in failed])

    return resolve_missing_ids_via_blast(df, protein_groups, blast_queries)
