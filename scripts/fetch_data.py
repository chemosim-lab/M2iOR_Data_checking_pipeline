# pipeline/scripts/fetch_uniprot.py

from __future__ import annotations

import time
from multiprocessing import Value
from pathlib import Path
from typing import TYPE_CHECKING, Any

import requests

if TYPE_CHECKING:
    import pandas as pd

from scripts.cache_manager import get_cache, get_missing_keys, set_cache
from scripts.fetch_blast import (
    clear_blast_cache_entries,
    count_uncached_blast_queries,
    get_failed_blast_queries,
)
from scripts.process_receptors import resolve_missing_ids_via_blast

_UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
_UNIPROTKB_ENDPOINT_URL = "https://rest.uniprot.org/uniprotkb/{accession}"
_NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
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
        if response.status_code == 400:
            print(
                f"  [UniProt] Not data found after fetching {accession}: {response.json()['messages']}"
            )
            return None
        else:
            return json_data if json_data else None
    except requests.RequestException as e:
        print(f"  [UniProt] Failed to fetch {accession}: {e}")
        return None


def fetch_uniprot_data(unique_ids: list[str]) -> list[str]:
    """Fetch UniProt data for each ID and cache it.
    Returns UIDs for which no data was found."""
    to_fetch = get_missing_keys(unique_ids)
    failed: list[str] = []

    if to_fetch:
        print(
            f"Fetching {len(to_fetch)} new accession(s) from UniProt "
            + f"({len(unique_ids) - len(to_fetch)} already cached)..."
        )
        for i, accession in enumerate(to_fetch, start=1):
            print(f"  [{i}/{len(to_fetch)}] {accession}", end=" ... \n", flush=False)

            accession_data: Any | None = _fetch_data(accession)

            if not accession_data:
                failed.append(accession)
                continue

            print("ok")
            set_cache(accession, accession_data)

            if i < len(to_fetch):
                time.sleep(_REQUEST_DELAY)
    else:
        print(f"All {len(unique_ids)} accession(s) found in cache, skipping API calls.")

    return failed


_HTTP_OK = 200


def _parse_gene_name_from_cds_header(header: str) -> str | None:
    """Extract gene name from a fasta_cds_aa header, e.g. [gene=Or54]."""
    for part in header.split("["):
        if part.startswith("gene="):
            return part.split("=", 1)[1].rstrip("]").strip()
    return None  # no [gene=...] tag found


def _fetch_ncbi_protein(accession: str) -> dict[str, Any] | None:
    try:
        response = requests.get(
            _NCBI_EFETCH_URL,
            params={
                "db": "nuccore",
                "id": accession,
                "rettype": "fasta_cds_aa",
                "retmode": "text",
            },
            timeout=15,
        )
        if response.status_code != _HTTP_OK:
            return None
        text = response.text.strip()
        if not text.startswith(">"):
            return None
        # Take only the first CDS entry (there may be several in the record)
        first_entry = text.split("\n>")[0]
        lines = first_entry.splitlines()
        header = lines[0][1:]  # remove leading ">"
        sequence = "".join(lines[1:])
        if not sequence:
            return None
        gene_name = _parse_gene_name_from_cds_header(header) or accession
    except requests.RequestException as e:
        print(f"  [NCBI] Failed to fetch {accession}: {e}")  # noqa: T201
        return None
    else:
        return {
            "source": "ncbi",
            "sequence": {"value": sequence},
            "genes": [{"geneName": {"value": gene_name}}],
        }


def fetch_ncbi_data(all_uids: list[str], failed_uids: list[str]) -> list[str]:
    """Try to fetch protein data from NCBI for UIDs that failed UniProt lookup.
    Returns (still_missing_uids)
    ncbi_refs is rebuilt from cache on every run for UIDs previously tagged as NCBI.
    """
    # Rebuild from cache first (covers UIDs already fetched in a previous run)
    ncbi_refs: dict[str, str] = {
        uid: data["sequence"]["value"]
        for uid in all_uids
        if (data := get_cache(uid)) and data.get("source") == "ncbi"
    }

    still_missing: list[str] = []
    for accession in failed_uids:
        if accession in ncbi_refs:
            continue  # already in cache from a previous run
        print(f"  [NCBI] {accession} ...", end=" ", flush=True)  # noqa: T201
        data = _fetch_ncbi_protein(accession)
        if data:
            print("ok")  # noqa: T201
            set_cache(accession, data)
        else:
            print("not found")  # noqa: T201
            still_missing.append(accession)
    return still_missing


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
