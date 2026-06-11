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
from scripts.columns import CAS as CAS_COL
from scripts.columns import CID as CID_COL
from scripts.columns import MOLECULE
from scripts.fetch_blast import (
    clear_blast_cache_entries,
    count_uncached_blast_queries,
    get_failed_blast_queries,
)
from scripts.process_receptors import resolve_missing_ids_via_blast

_UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
_UNIPROTKB_ENDPOINT_URL = "https://rest.uniprot.org/uniprotkb/{accession}"
_NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON"
_PUBCHEM_SYNONYMS_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
)
_PUBCHEM_VIEW_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?response_type=display"
_PUBCHEM_CAS_CID_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas}/cids/JSON"
)
_APA_URL = "https://doi.org/{doi}"
_APA_HEADERS = {"Accept": "text/x-bibliography; style=apa; locale=en-US"}
_REQUEST_DELAY = 0.2  # seconds between requests

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
    to_fetch = get_missing_keys(unique_ids, subdir="receptors")
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
            set_cache(accession, accession_data, subdir="receptors")

            if i < len(to_fetch):
                time.sleep(_REQUEST_DELAY)
    else:
        print(f"All {len(unique_ids)} accession(s) found in cache, skipping API calls.")

    return failed


_HTTP_OK = 200


def _parse_tag_from_cds_header(header: str, tag: str) -> str | None:
    """Extract a bracketed tag value from a fasta_cds_aa header, e.g. [gene=Or54]."""
    prefix = f"{tag}="
    for part in header.split("["):
        if part.startswith(prefix):
            return part.split("=", 1)[1].rstrip("]").strip()
    return None


def _parse_gene_name_from_cds_header(header: str) -> str | None:
    return _parse_tag_from_cds_header(header, "gene")


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
        protein_desc = _parse_tag_from_cds_header(header, "protein")
    except requests.RequestException as e:
        print(f"  [NCBI] Failed to fetch {accession}: {e}")  # noqa: T201
        return None
    else:
        return {
            "source": "ncbi",
            "sequence": {"value": sequence},
            "genes": [{"geneName": {"value": gene_name}}],
            "protein_description": protein_desc,
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
        if (data := get_cache(uid, subdir="receptors")) and data.get("source") == "ncbi"
    }

    still_missing: list[str] = []
    for accession in failed_uids:
        if accession in ncbi_refs:
            continue  # already in cache from a previous run
        print(f"  [NCBI] {accession} ...", end=" ", flush=True)  # noqa: T201
        data = _fetch_ncbi_protein(accession)
        if data:
            print("ok")  # noqa: T201
            set_cache(accession, data, subdir="receptors")
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

    # "no hit found" results are definitive — no point retrying them.
    no_hit = [(s, sp, err) for s, sp, err in failed if err == "no hit found"]
    retryable = [(s, sp, err) for s, sp, err in failed if err != "no hit found"]

    if new_count > 0 or failed:
        print(f"\n  ! {len(blast_queries)} sequence(s) have no UniProt ID:")  # noqa: T201
        if new_count > 0:
            print(f"    {new_count} new (never queried, ~1-5 min each)")  # noqa: T201
        for _seq, sp, err in retryable:
            print(f"    previously failed — {sp}: {err}")  # noqa: T201
        for _seq, sp, _err in no_hit:
            print(f"    no hit found on previous attempt — {sp}")  # noqa: T201

        total = new_count + len(retryable) + len(no_hit)
        answer = input(f"  Run/retry {total} BLAST lookup(s)? [y/N] ").strip().lower()
        if answer != "y":
            if new_count > 0 or retryable:
                # Cannot proceed without these — new queries or recoverable errors.
                count = new_count + len(retryable)
                msg = f"Aborted by user: {count} BLAST lookup(s) required."
                raise ValueError(msg)
            # Only "no hit found" entries: safe to skip, nothing new to find.
            return resolve_missing_ids_via_blast(df, protein_groups, blast_queries)

        if retryable:
            clear_blast_cache_entries([(s, sp) for s, sp, _ in retryable])
        if no_hit:
            clear_blast_cache_entries([(s, sp) for s, sp, _ in no_hit])

    return resolve_missing_ids_via_blast(df, protein_groups, blast_queries)


def fetch_pubchem_data(unique_cids: list[int]) -> list[int]:
    """Fetch PubChem compound data for each CID and cache it.
    Returns CIDs for which no data was found."""
    missing = set(
        get_missing_keys([str(cid) for cid in unique_cids], subdir="molecules")
    )
    to_fetch = [cid for cid in unique_cids if str(cid) in missing]
    failed: list[int] = []

    if to_fetch:
        print(
            f"Fetching {len(to_fetch)} new CID(s) from PubChem "
            f"({len(unique_cids) - len(to_fetch)} already cached)..."
        )
        for i, cid in enumerate(to_fetch, start=1):
            print(f"  [{i}/{len(to_fetch)}] CID:{cid}", end=" ... ", flush=True)
            try:
                response = requests.get(_PUBCHEM_URL.format(cid=cid), timeout=10)
                if response.status_code != _HTTP_OK:
                    print(f"not found (HTTP {response.status_code})")
                    failed.append(cid)
                    continue

                data = response.json()

                time.sleep(_REQUEST_DELAY)

                syn_response = requests.get(
                    _PUBCHEM_SYNONYMS_URL.format(cid=cid), timeout=10
                )
                if syn_response.status_code == _HTTP_OK:
                    synonyms = (
                        syn_response.json()
                        .get("InformationList", {})
                        .get("Information", [{}])[0]
                        .get("Synonym", [])
                    )
                    data["synonyms"] = synonyms

                time.sleep(_REQUEST_DELAY)

                view_response = requests.get(
                    _PUBCHEM_VIEW_URL.format(cid=cid), timeout=10
                )
                if view_response.status_code == _HTTP_OK:
                    data["record_title"] = (
                        view_response.json().get("Record", {}).get("RecordTitle")
                    )

                set_cache(str(cid), data, subdir="molecules")
                print("ok")

            except requests.RequestException as e:
                print(f"error: {e}")
                failed.append(cid)

            if i < len(to_fetch):
                time.sleep(_REQUEST_DELAY)
    else:
        print(f"All {len(unique_cids)} CID(s) found in cache, skipping API calls.")

    return failed


def fetch_cids_from_cas(df: pd.DataFrame) -> list[int]:
    """For Molecule rows with an empty CID, look up the CID via CAS on PubChem.
    Updates the CID column in df. Returns the list of newly found CIDs."""
    molecule_df = df[MOLECULE]
    no_cid = molecule_df[CID_COL].isna() | (molecule_df[CID_COL] == 0)
    missing_cid = no_cid & molecule_df[CAS_COL].notna()
    cas_series = molecule_df.loc[missing_cid, CAS_COL]

    if cas_series.empty:
        return []

    unique_cas = cas_series.unique().tolist()
    cas_to_cid: dict[str, int] = {}

    print(f"CID fallback: querying {len(unique_cas)} CAS number(s) on PubChem...")  # noqa: T201
    for i, cas in enumerate(unique_cas, start=1):
        print(f"  [{i}/{len(unique_cas)}] CAS:{cas}", end=" ... ", flush=True)  # noqa: T201
        try:
            response = requests.get(_PUBCHEM_CAS_CID_URL.format(cas=cas), timeout=10)
            if response.status_code == _HTTP_OK:
                cids = response.json().get("IdentifierList", {}).get("CID", [])
                if cids:
                    cas_to_cid[cas] = cids[0]
                    print(f"CID:{cids[0]}")  # noqa: T201
                else:
                    print("not found")  # noqa: T201
            else:
                print(f"not found (HTTP {response.status_code})")  # noqa: T201
        except requests.RequestException as e:
            print(f"error: {e}")  # noqa: T201

        if i < len(unique_cas):
            time.sleep(_REQUEST_DELAY)

    if not cas_to_cid:
        return []

    updated = cas_series.map(cas_to_cid)
    mask = updated.notna()
    df.loc[updated.index[mask], (MOLECULE, CID_COL)] = (
        updated[mask].astype(int).to_numpy()
    )

    return list(cas_to_cid.values())


def fetch_apa_references(unique_dois: list[str]) -> list[str]:
    """Fetch APA-formatted references for each DOI and cache them.
    Returns DOIs for which no reference was found."""
    to_fetch = get_missing_keys(unique_dois, subdir="sources")
    failed: list[str] = []

    if to_fetch:
        print(f"Fetching {len(to_fetch)} APA reference(s) from doi.org...")  # noqa: T201
        for i, doi in enumerate(to_fetch, start=1):
            print(f"  [{i}/{len(to_fetch)}] {doi}", end=" ... ", flush=True)  # noqa: T201
            try:
                response = requests.get(
                    _APA_URL.format(doi=doi),
                    headers=_APA_HEADERS,
                    timeout=15,
                )
                if response.status_code == _HTTP_OK:
                    set_cache(doi, {"apa": response.text.strip()}, subdir="sources")
                    print("ok")  # noqa: T201
                else:
                    print(f"not found (HTTP {response.status_code})")  # noqa: T201
                    failed.append(doi)
            except requests.RequestException as e:
                print(f"error: {e}")  # noqa: T201
                failed.append(doi)

            if i < len(to_fetch):
                time.sleep(_REQUEST_DELAY)
    else:
        print(f"All {len(unique_dois)} DOI(s) found in APA cache, skipping API calls.")  # noqa: T201

    return failed
