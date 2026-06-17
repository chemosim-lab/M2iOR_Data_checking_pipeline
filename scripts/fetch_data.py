# pipeline/scripts/fetch_uniprot.py

from __future__ import annotations

import time
from pathlib import Path
from typing import TYPE_CHECKING, Any

import colorlog
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
    get_successful_blast_queries,
)
from scripts.process_receptors import resolve_missing_accessions_via_blast

logger = colorlog.getLogger(__name__)

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


def _fetch_data(accession: str) -> Any | None:
    url = _UNIPROTKB_ENDPOINT_URL.format(accession=accession)
    try:
        response = requests.get(url, timeout=10)
        json_data = (
            response.json() if response and response.status_code == 200 else None
        )
        if response.status_code == 400:
            logger.info(
                "  [UniProt] Not data found after fetching %s: %s",
                accession,
                response.json()["messages"],
            )
            return None
        else:
            return json_data if json_data else None
    except requests.RequestException as e:
        logger.info("  [UniProt] Failed to fetch %s: %s", accession, e)
        return None


def fetch_uniprot_data(unique_accessions: list[str]) -> list[str]:
    """Fetch UniProt data for each ID and cache it.
    Returns UIDs for which no data was found."""
    to_fetch = get_missing_keys(unique_accessions, subdir="receptors")
    failed: list[str] = []

    if to_fetch:
        logger.info(
            "Fetching %d new accession(s) from UniProt (%d already cached)...",
            len(to_fetch),
            len(unique_accessions) - len(to_fetch),
        )
        for i, accession in enumerate(to_fetch, start=1):
            accession_data: Any | None = _fetch_data(accession)
            status = "ok" if accession_data else "not found"
            logger.info("  [%d/%d] %s ... %s", i, len(to_fetch), accession, status)
            if not accession_data:
                failed.append(accession)
                continue

            set_cache(accession, accession_data, subdir="receptors")

            if i < len(to_fetch):
                time.sleep(_REQUEST_DELAY)
    else:
        logger.info(
            "All %d unique accession(s) found in cache, skipping API calls.",
            len(unique_accessions),
        )

    return failed


_HTTP_OK = 200



def _fetch_ncbi_protein(accession: str) -> dict[str, Any] | None:
    try:
        response = requests.get(
            _NCBI_EFETCH_URL,
            params={
                "db": "protein",
                "id": accession,
                "rettype": "fasta",
                "retmode": "text",
            },
            timeout=15,
        )
        if response.status_code != _HTTP_OK:
            return None
        text = response.text.strip()
        if not text.startswith(">"):
            return None
        lines = text.splitlines()
        header = lines[0][1:]  # remove leading ">"
        sequence = "".join(lines[1:])
        if not sequence:
            return None
        # Protein FASTA header: ">AAT71306.1 <description> [organism]"
        parts = header.split(" ", 1)
        gene_name = accession
        protein_desc = parts[1].strip() if len(parts) > 1 else None
    except requests.RequestException as e:
        logger.info("  [NCBI] Failed to fetch %s: %s", accession, e)
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
        data = _fetch_ncbi_protein(accession)
        if data:
            logger.info("  [NCBI] %s ... ok", accession)
            set_cache(accession, data, subdir="receptors")
        else:
            logger.info("  [NCBI] %s ... not found", accession)
            still_missing.append(accession)
    return still_missing


def _print_blast_pending(
    new_count: int,
    retryable: list[tuple[str, str, str]],
    no_hit: list[tuple[str, str, str]],
    cached_ok: list[tuple[str, str]],
) -> None:
    if new_count > 0:
        logger.info("%d new (never queried, ~1-5 min each)", new_count)
    for _seq, sp, err in retryable:
        logger.info("Previously failed - (%s, %s): %s", sp, _seq[:9] + "...", err)
    for _seq, sp, _err in no_hit:
        logger.info("No hit found on previous attempt - (%s, %s)", sp, _seq[:9] + "...")
    for _seq, sp in cached_ok:
        logger.info("Cached - force re-run: (%s, %s)", sp, _seq[:9] + "...")


def _clear_blast_entries(
    retryable: list[tuple[str, str, str]],
    no_hit: list[tuple[str, str, str]],
    cached_ok: list[tuple[str, str]],
) -> None:
    if retryable:
        clear_blast_cache_entries([(s, sp) for s, sp, _ in retryable])
    if no_hit:
        clear_blast_cache_entries([(s, sp) for s, sp, _ in no_hit])
    if cached_ok:
        clear_blast_cache_entries(cached_ok)


def fetch_genbank_data(
    df: pd.DataFrame,
    protein_groups: list[str],
    blast_queries: list[tuple[str, str]],
    *,
    force_blast: bool = False,
) -> dict[str, str]:
    if not blast_queries:
        return {}

    new_count = count_uncached_blast_queries(blast_queries)
    failed: list[tuple[str, str, str]] = get_failed_blast_queries(blast_queries)
    cached_ok = get_successful_blast_queries(blast_queries) if force_blast else []

    # "no hit found" results are definitive — no point retrying them.
    no_hit = [(s, sp, err) for s, sp, err in failed if err == "no hit found"]
    retryable = [(s, sp, err) for s, sp, err in failed if err != "no hit found"]

    if new_count > 0 or failed or cached_ok:
        _print_blast_pending(new_count, retryable, no_hit, cached_ok)
        total = new_count + len(retryable) + len(no_hit) + len(cached_ok)
        answer = input(f"  Run/retry {total} BLAST lookup(s)? [y/N] ").strip().lower()
        if answer != "y":
            logger.info("Skipping BLAST lookups, continuing without new BLAST results.")
            return {}
        _clear_blast_entries(retryable, no_hit, cached_ok)

    return resolve_missing_accessions_via_blast(df, protein_groups, blast_queries)


def fetch_pubchem_data(unique_cids: list[int]) -> list[int]:
    """Fetch PubChem compound data for each CID and cache it.
    Returns CIDs for which no data was found."""
    missing = set(
        get_missing_keys([str(cid) for cid in unique_cids], subdir="molecules")
    )
    to_fetch = [cid for cid in unique_cids if str(cid) in missing]
    failed: list[int] = []

    if to_fetch:
        logger.info(
            "Fetching %d new CID(s) from PubChem (%d already cached)...",
            len(to_fetch),
            len(unique_cids) - len(to_fetch),
        )
        for i, cid in enumerate(to_fetch, start=1):
            try:
                response = requests.get(_PUBCHEM_URL.format(cid=cid), timeout=10)
                if response.status_code != _HTTP_OK:
                    logger.info(
                        "  [%d/%d] CID:%s ... not found (HTTP %d)",
                        i,
                        len(to_fetch),
                        cid,
                        response.status_code,
                    )
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
                logger.info("  [%d/%d] CID:%s ... ok", i, len(to_fetch), cid)

            except requests.RequestException as e:
                logger.info("  [%d/%d] CID:%s ... error: %s", i, len(to_fetch), cid, e)
                failed.append(cid)

            if i < len(to_fetch):
                time.sleep(_REQUEST_DELAY)
    else:
        logger.info(
            "All %d CID(s) found in cache, skipping API calls.", len(unique_cids)
        )

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

    logger.info(
        "CID fallback: querying %d CAS number(s) on PubChem...", len(unique_cas)
    )
    for i, cas in enumerate(unique_cas, start=1):
        try:
            response = requests.get(_PUBCHEM_CAS_CID_URL.format(cas=cas), timeout=10)
            if response.status_code == _HTTP_OK:
                cids = response.json().get("IdentifierList", {}).get("CID", [])
                if cids:
                    cas_to_cid[cas] = cids[0]
                    logger.info(
                        "  [%d/%d] CAS:%s ... CID:%s", i, len(unique_cas), cas, cids[0]
                    )
                else:
                    logger.info(
                        "  [%d/%d] CAS:%s ... not found", i, len(unique_cas), cas
                    )
            else:
                logger.info(
                    "  [%d/%d] CAS:%s ... not found (HTTP %d)",
                    i,
                    len(unique_cas),
                    cas,
                    response.status_code,
                )
        except requests.RequestException as e:
            logger.info("  [%d/%d] CAS:%s ... error: %s", i, len(unique_cas), cas, e)

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
        logger.info("Fetching %d APA reference(s) from doi.org...", len(to_fetch))
        for i, doi in enumerate(to_fetch, start=1):
            try:
                response = requests.get(
                    _APA_URL.format(doi=doi),
                    headers=_APA_HEADERS,
                    timeout=15,
                )
                if response.status_code == _HTTP_OK:
                    response.encoding = "utf-8"
                    set_cache(doi, {"apa": response.text.strip()}, subdir="sources")
                    logger.info("  [%d/%d] %s ... ok", i, len(to_fetch), doi)
                else:
                    logger.info(
                        "  [%d/%d] %s ... not found (HTTP %d)",
                        i,
                        len(to_fetch),
                        doi,
                        response.status_code,
                    )
                    failed.append(doi)
            except requests.RequestException as e:
                logger.info("  [%d/%d] %s ... error: %s", i, len(to_fetch), doi, e)
                failed.append(doi)

            if i < len(to_fetch):
                time.sleep(_REQUEST_DELAY)
    else:
        logger.info(
            "All %d DOI(s) found in APA cache, skipping API calls.", len(unique_dois)
        )

    return failed
