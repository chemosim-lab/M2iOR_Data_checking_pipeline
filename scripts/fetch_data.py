# pipeline/scripts/fetch_uniprot.py

from __future__ import annotations

import os
import re
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import colorlog
import requests

if TYPE_CHECKING:
    import pandas as pd

from scripts.cache_manager import get_cache, get_missing_keys, set_cache
from scripts.columns import CAS as CAS_COL
from scripts.columns import CID as CID_COL
from scripts.columns import MOLECULE, RECEPTOR_NAME, SEQUENCE, SPECIES
from scripts.fetch_blast import (
    clear_blast_cache_entries,
    get_failed_blast_queries,
    get_new_blast_queries,
    get_successful_blast_queries,
)
from scripts.process_molecules import normalize_dashes, parse_cid_cell
from scripts.process_receptors import NOT_AVAILABLE, resolve_missing_accessions_via_blast
from scripts.report import Issue, emit

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
_PUBCHEM_INCHIKEY_CID_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"
)
_CAS_COMMON_CHEMISTRY_URL = "https://commonchemistry.cas.org/api/detail"
_CAS_COMMON_CHEMISTRY_SANITY_CAS = "50-00-0"  # formaldehyde
_APA_URL = "https://doi.org/{doi}"
_APA_HEADERS = {"Accept": "text/x-bibliography; style=apa; locale=en-US"}
_REQUEST_DELAY = 0.2  # seconds between requests

_CACHE_FILE = Path(__file__).parent.parent / "cache" / "uniprot_sequences.json"

_DOTENV_PATH = Path(__file__).parent.parent / ".env"


def _load_dotenv(path: Path = _DOTENV_PATH) -> None:
    """Load KEY=VALUE lines from `path` into the environment (existing env
    vars win). Minimal by design - CAS_API_KEY is the only secret this
    project reads, so a dependency for a handful of lines isn't worth it."""
    if not path.exists():
        return
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        os.environ.setdefault(key.strip(), value.strip().strip('"').strip("'"))


def get_cas_api_key() -> str | None:
    """Return CAS_API_KEY from the environment or .env (see .env.example)."""
    _load_dotenv()
    return os.environ.get("CAS_API_KEY")


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
        if json_data and "sequence" not in json_data:
            # Entry exists but is inactive (deleted/merged/demerged) - UniProt
            # returns 200 with an `inactiveReason` instead of a sequence. Treat
            # it like "not found" so it isn't cached as a usable reference and
            # falls through to the NCBI/BLAST fallback instead.
            reason = json_data.get("inactiveReason", {}).get(
                "inactiveReasonType", "unknown"
            )
            logger.info(
                "  [UniProt] %s is inactive (%s), no sequence available",
                accession,
                reason,
            )
            return None
        return json_data
    except requests.RequestException as e:
        logger.info("  [UniProt] Failed to fetch %s: %s", accession, e)
        return None


def fetch_uniprot_data(unique_accessions: list[str]) -> list[str]:
    """Fetch UniProt data for each ID and cache it.
    Returns UIDs for which no data was found."""
    to_fetch = get_missing_keys(unique_accessions, subdir="receptors")
    failed: list[str] = []

    # The "Not available" sentinel is never a real accession - sending it to
    # UniProt/NCBI only produces a confusing "invalid format" error. Route it
    # straight to `failed` (and thus the BLAST-on-sequence fallback) instead.
    sentinel_values = [a for a in to_fetch if a.lower() == NOT_AVAILABLE.lower()]
    if sentinel_values:
        to_fetch = [a for a in to_fetch if a.lower() != NOT_AVAILABLE.lower()]
        failed.extend(sentinel_values)
        logger.info(
            '  %d accession(s) marked "%s" - skipping UniProt/NCBI lookup, '
            "falling back to BLAST on sequence.",
            len(sentinel_values),
            NOT_AVAILABLE,
        )

    already_cached = len(unique_accessions) - len(to_fetch) - len(sentinel_values)

    if to_fetch:
        logger.info(
            "Fetching %d new accession(s) from UniProt (%d already cached)...",
            len(to_fetch),
            already_cached,
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
    elif already_cached:
        logger.info(
            "All %d unique accession(s) found in cache, skipping API calls.",
            already_cached,
        )

    return failed


_HTTP_OK = 200

_NCBI_PROTEIN_ID_RE = re.compile(r'/protein_id="([^"]+)"')


def _fetch_ncbi_protein_fasta(accession: str) -> dict[str, Any] | None:
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


def _fetch_ncbi_nuccore_protein_id(accession: str) -> str | None:
    """Look up the protein_id referenced by the CDS feature of a nucleotide record.

    Some accessions we're given actually point to a nuccore (nucleotide) record
    rather than a protein one (e.g. "KM229532.1"). In that case the real protein
    sequence lives under a separate protein_id from the record's CDS feature
    (e.g. "AIO10894.1"), so fetching db=protein directly fails.
    """
    try:
        response = requests.get(
            _NCBI_EFETCH_URL,
            params={
                "db": "nuccore",
                "id": accession,
                "rettype": "gb",
                "retmode": "text",
            },
            timeout=15,
        )
        if response.status_code != _HTTP_OK:
            return None
        match = _NCBI_PROTEIN_ID_RE.search(response.text)
        return match.group(1) if match else None
    except requests.RequestException as e:
        logger.info("  [NCBI] Failed to fetch nuccore record %s: %s", accession, e)
        return None


def _fetch_ncbi_protein(accession: str) -> dict[str, Any] | None:
    data = _fetch_ncbi_protein_fasta(accession)
    if data:
        return data

    time.sleep(_REQUEST_DELAY)
    protein_id = _fetch_ncbi_nuccore_protein_id(accession)
    if not protein_id:
        return None
    logger.info(
        "  [NCBI] %s is a nucleotide accession, using its protein_id %s instead",
        accession,
        protein_id,
    )
    time.sleep(_REQUEST_DELAY)
    return _fetch_ncbi_protein_fasta(protein_id)


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
    calls_made = 0
    for accession in failed_uids:
        if accession.lower() == NOT_AVAILABLE.lower():
            # Already logged once in fetch_uniprot_data - no real ID to try here.
            still_missing.append(accession)
            continue
        if accession in ncbi_refs:
            continue  # already in cache from a previous run

        if calls_made:
            time.sleep(_REQUEST_DELAY)
        calls_made += 1

        data = _fetch_ncbi_protein(accession)
        if data:
            logger.info("  [NCBI] %s ... ok", accession)
            set_cache(accession, data, subdir="receptors")
        else:
            logger.info("  [NCBI] %s ... not found", accession)
            still_missing.append(accession)
    return still_missing


class _BlastCandidate:
    def __init__(
        self,
        category: str,
        seq: str,
        species: str,
        receptor_name: str,
        detail: str | None = None,
    ) -> None:
        self.category = category
        self.seq = seq
        self.species = species
        self.receptor_name = receptor_name
        self.detail = detail

    @property
    def key(self) -> tuple[str, str]:
        return (self.seq, self.species)


_CATEGORY_LABELS = {
    "new": "new",
    "retry": "retry",
    "no-hit": "no-hit",
    "cached": "redo",
}


def _receptor_names_by_query(
    df: pd.DataFrame, protein_groups: list[str], queries: list[tuple[str, str]]
) -> dict[tuple[str, str], str]:
    """Map each (sequence, species) query to a Receptor Name for display."""
    wanted = set(queries)
    names: dict[tuple[str, str], str] = {}
    for group in protein_groups:
        if group not in df.columns.get_level_values(0):
            continue
        sub = df[group]
        # fillna before astype(str): pandas 3.0 leaves NaN as a float instead
        # of stringifying it to "nan", which breaks downstream str ops/slicing.
        seq = sub[SEQUENCE].fillna("").astype(str).str.replace(r"\s+", "", regex=True)
        species = sub[SPECIES].fillna("").astype(str).str.strip()
        name = sub[RECEPTOR_NAME].fillna("").astype(str).str.strip()
        for s, sp, n in zip(seq, species, name, strict=False):
            key = (s, sp)
            if key in wanted and key not in names:
                names[key] = n
    return names


def _build_blast_candidates(
    df: pd.DataFrame,
    protein_groups: list[str],
    new_queries: list[tuple[str, str]],
    retryable: list[tuple[str, str, str]],
    no_hit: list[tuple[str, str, str]],
    cached_ok: list[tuple[str, str]],
) -> list[_BlastCandidate]:
    all_queries = (
        new_queries
        + [(s, sp) for s, sp, _ in retryable]
        + [(s, sp) for s, sp, _ in no_hit]
        + cached_ok
    )
    names = _receptor_names_by_query(df, protein_groups, all_queries)

    candidates = [
        _BlastCandidate("new", s, sp, names.get((s, sp), "?")) for s, sp in new_queries
    ]
    candidates += [
        _BlastCandidate("retry", s, sp, names.get((s, sp), "?"), err)
        for s, sp, err in retryable
    ]
    candidates += [
        _BlastCandidate("no-hit", s, sp, names.get((s, sp), "?"), err)
        for s, sp, err in no_hit
    ]
    candidates += [
        _BlastCandidate("cached", s, sp, names.get((s, sp), "?")) for s, sp in cached_ok
    ]
    return candidates


def _parse_selection(answer: str, count: int) -> list[int]:
    """Parse "1,3-5" style input into sorted, deduplicated 1-based indices."""
    indices: set[int] = set()
    for part in answer.split(","):
        part = part.strip()
        if not part:
            continue
        lo_s, _, hi_s = part.partition("-")
        lo, hi = int(lo_s), int(hi_s) if hi_s else int(lo_s)
        indices.update(range(lo, hi + 1))
    return sorted(i for i in indices if 1 <= i <= count)


def _select_blast_candidates(
    candidates: list[_BlastCandidate],
) -> list[_BlastCandidate]:
    """Print the numbered pending BLAST lookups and let the user pick which
    combination to (re-)run, instead of an all-or-nothing prompt."""
    for i, c in enumerate(candidates, start=1):
        seq_preview = c.seq[:12] + ("..." if len(c.seq) > 12 else "")
        detail = f" - {c.detail}" if c.detail else ""
        logger.info(
            "  [%2d] %-12s %-20s %-20s %s%s",
            i,
            _CATEGORY_LABELS[c.category],
            c.receptor_name[:20],
            c.species[:20],
            seq_preview,
            detail,
        )
    answer = (
        input(
            f'  Run which BLAST lookup(s)? (e.g. "1,3-5", "all", Enter to skip) [{len(candidates)} pending] '
        )
        .strip()
        .lower()
    )
    if not answer:
        return []
    if answer == "all":
        return candidates
    selected = _parse_selection(answer, len(candidates))
    return [candidates[i - 1] for i in selected]


# How pending BLAST lookups are handled: "ask" prompts for which to run,
# "cached" runs none (only already-cached results are applied), "run" runs
# every new or retryable one without asking.
BlastMode = Literal["ask", "cached", "run"]


def _select_blast_candidates_non_interactive(
    candidates: list[_BlastCandidate], mode: BlastMode
) -> list[_BlastCandidate]:
    if mode == "run":
        # A definitive "no hit" is never worth re-running unattended.
        return [c for c in candidates if c.category in ("new", "retry")]
    pending = [c for c in candidates if c.category in ("new", "retry")]
    if pending:
        logger.warning(
            "%d BLAST lookup(s) pending, not run (--blast cached).", len(pending)
        )
        emit(
            Issue(
                severity="warning",
                code="blast_pending",
                message=f"{len(pending)} sequence(s) without a usable accession "
                "need a BLAST lookup that wasn't run (--blast cached); their "
                "rows keep an empty accession.",
                details={
                    "queries": [
                        {
                            "receptor_name": c.receptor_name,
                            "species": c.species,
                            "category": c.category,
                            "sequence_start": c.seq[:20],
                        }
                        for c in pending
                    ]
                },
            )
        )
    return []


def fetch_genbank_data(
    df: pd.DataFrame,
    protein_groups: list[str],
    blast_queries: list[tuple[str, str]],
    *,
    force_blast: bool = False,
    mode: BlastMode = "ask",
) -> dict[str, str]:
    if not blast_queries:
        return {}

    new_queries = get_new_blast_queries(blast_queries)
    failed: list[tuple[str, str, str]] = get_failed_blast_queries(blast_queries)
    already_resolved = get_successful_blast_queries(blast_queries)
    # Only offered as re-runnable candidates under --force; otherwise they're
    # applied for free below without asking (no network call needed).
    cached_ok = already_resolved if force_blast else []

    # "no hit found" results are definitive — no point retrying them.
    no_hit = [(s, sp, err) for s, sp, err in failed if err == "no hit found"]
    retryable = [(s, sp, err) for s, sp, err in failed if err != "no hit found"]

    # Results already sitting in the cache are always applied, whether or not
    # they're also offered above as re-runnable "cached" candidates - picking
    # them there just forces a fresh BLAST query instead of reusing this.
    queries_to_resolve: set[tuple[str, str]] = set(already_resolved)

    if new_queries or failed or cached_ok:
        if already_resolved and not force_blast:
            logger.info(
                "%d already resolved from BLAST cache (applied automatically, no lookup needed)",
                len(already_resolved),
            )
        candidates = _build_blast_candidates(
            df, protein_groups, new_queries, retryable, no_hit, cached_ok
        )
        selected = (
            _select_blast_candidates(candidates)
            if mode == "ask"
            else _select_blast_candidates_non_interactive(candidates, mode)
        )
        if not selected:
            logger.info("Skipping BLAST lookups, continuing without new BLAST results.")
        else:
            to_clear = [c.key for c in selected if c.category != "new"]
            if to_clear:
                clear_blast_cache_entries(to_clear)
            queries_to_resolve.update(c.key for c in selected)
    elif already_resolved:
        logger.info(
            "All %d sequence(s) without accession resolved from BLAST cache.",
            len(already_resolved),
        )

    if not queries_to_resolve:
        return {}
    return resolve_missing_accessions_via_blast(
        df, protein_groups, list(queries_to_resolve)
    )


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
    elif unique_cids:
        logger.info(
            "All %d PubChem compound record(s) already cached, skipping API calls.",
            len(unique_cids),
        )

    return failed


_CAS_CACHE_SUBDIR = "cas_to_cid"
_CAS_COMMON_CHEMISTRY_CACHE_SUBDIR = "cas_common_chemistry"


def _cas_common_chemistry_key_is_valid(api_key: str) -> bool:
    """One cheap lookup to tell an invalid/rejected API key apart from a CAS
    number that's simply not in CAS Common Chemistry's ~500k-substance set -
    so a bad key logs once instead of a "not found" per CAS number."""
    try:
        response = requests.get(
            _CAS_COMMON_CHEMISTRY_URL,
            params={"cas_rn": _CAS_COMMON_CHEMISTRY_SANITY_CAS},
            headers={"X-API-KEY": api_key},
            timeout=10,
        )
    except requests.RequestException:
        return True  # transient network issue, not a key problem - let the batch try
    return response.status_code not in (401, 403)


def _report_cas_api_unavailable(reason: str) -> None:
    emit(
        Issue(
            severity="warning",
            code="cas_api_unavailable",
            message=f"{reason}: CAS numbers missing from the cache are checked "
            "against PubChem only, not CAS Common Chemistry.",
        )
    )


def fetch_cas_common_chemistry_details(
    unique_cas: list[str],
) -> dict[str, dict[str, Any]]:
    """Fetch (or read from cache) each CAS number's own record directly from
    CAS Common Chemistry - the registry that assigns CAS numbers - via
    https://commonchemistry.cas.org/api/detail?cas_rn={cas}.

    This is the priority source for anything CAS-related: a direct,
    unambiguous lookup by registry number, rather than indirectly resolving
    a CAS number through PubChem's fuzzy name search. Returns a
    cas -> {"name", "synonyms", "inchikey", "smiles"} mapping, omitting any
    CAS Common Chemistry has no record for (a meaningful fraction - its
    ~500k "common" substances are far fewer than PubChem's 100M+) - callers
    fall back to PubChem for those.

    Requires CAS_API_KEY (see .env.example); if it's unset or rejected,
    every CAS is skipped (logged once) so callers' PubChem fallback still runs.
    """
    if not unique_cas:
        return {}

    api_key = get_cas_api_key()
    if not api_key:
        logger.info(
            "CAS_API_KEY not set - skipping CAS Common Chemistry, falling back "
            "to PubChem for CAS number(s) (see .env.example)."
        )
        _report_cas_api_unavailable("CAS_API_KEY is not set")
        return {}

    details: dict[str, dict[str, Any]] = {
        cas: cached
        for cas in unique_cas
        if (cached := get_cache(cas, subdir=_CAS_COMMON_CHEMISTRY_CACHE_SUBDIR))
        is not None
        and not cached.get("not_found")
    }

    to_fetch = get_missing_keys(unique_cas, subdir=_CAS_COMMON_CHEMISTRY_CACHE_SUBDIR)
    if not to_fetch:
        return details

    if not _cas_common_chemistry_key_is_valid(api_key):
        logger.warning(
            "CAS Common Chemistry rejected CAS_API_KEY (401/403) - skipping it for "
            "this run, falling back to PubChem for CAS number(s)."
        )
        _report_cas_api_unavailable("CAS Common Chemistry rejected CAS_API_KEY (401/403)")
        return details

    logger.info(
        "Querying %d CAS number(s) on CAS Common Chemistry (%d already cached)...",
        len(to_fetch),
        len(unique_cas) - len(to_fetch),
    )
    headers = {"X-API-KEY": api_key}
    for i, cas in enumerate(to_fetch, start=1):
        try:
            response = requests.get(
                _CAS_COMMON_CHEMISTRY_URL,
                params={"cas_rn": cas},
                headers=headers,
                timeout=10,
            )
            if response.status_code == _HTTP_OK:
                data = response.json()
                record = {
                    "name": data.get("name"),
                    "synonyms": data.get("synonyms", []),
                    "inchikey": data.get("inchiKey"),
                    "smiles": data.get("canonicalSmile"),
                }
                set_cache(cas, record, subdir=_CAS_COMMON_CHEMISTRY_CACHE_SUBDIR)
                details[cas] = record
                logger.info("  [%d/%d] CAS:%s ... ok", i, len(to_fetch), cas)
            elif response.status_code == 404:
                set_cache(
                    cas, {"not_found": True}, subdir=_CAS_COMMON_CHEMISTRY_CACHE_SUBDIR
                )
                logger.info("  [%d/%d] CAS:%s ... not found", i, len(to_fetch), cas)
            else:
                logger.info(
                    "  [%d/%d] CAS:%s ... error (HTTP %d)",
                    i,
                    len(to_fetch),
                    cas,
                    response.status_code,
                )
        except requests.RequestException as e:
            logger.info("  [%d/%d] CAS:%s ... error: %s", i, len(to_fetch), cas, e)

        if i < len(to_fetch):
            time.sleep(_REQUEST_DELAY)

    return details


def _cid_from_inchikey(inchikey: str) -> int | None:
    try:
        response = requests.get(
            _PUBCHEM_INCHIKEY_CID_URL.format(inchikey=inchikey), timeout=10
        )
        if response.status_code == _HTTP_OK:
            cids = response.json().get("IdentifierList", {}).get("CID", [])
            return cids[0] if cids else None
    except requests.RequestException:
        pass
    return None


def _cid_via_pubchem_name_search(cas: str, i: int, total: int) -> int | None:
    """Fallback resolution: search PubChem by the CAS number as if it were a
    compound name/synonym (used when CAS Common Chemistry has no InChIKey for
    this CAS to cross-reference directly)."""
    try:
        response = requests.get(_PUBCHEM_CAS_CID_URL.format(cas=cas), timeout=10)
    except requests.RequestException as e:
        logger.info("  [%d/%d] CAS:%s ... error: %s", i, total, cas, e)
        return None
    if response.status_code != _HTTP_OK:
        logger.info(
            "  [%d/%d] CAS:%s ... not found (HTTP %d)",
            i,
            total,
            cas,
            response.status_code,
        )
        return None
    cids = response.json().get("IdentifierList", {}).get("CID", [])
    return cids[0] if cids else None


def _resolve_cas_to_cid(
    cas: str, inchikey: str | None, i: int, total: int
) -> tuple[int | None, str]:
    """Resolve one CAS number to a PubChem CID, returning (cid, source label)."""
    if inchikey:
        cid = _cid_from_inchikey(inchikey)
        if cid:
            return cid, "CAS Common Chemistry InChIKey"
    return _cid_via_pubchem_name_search(cas, i, total), "PubChem name search"


def fetch_cas_to_cid_map(unique_cas: list[str]) -> dict[str, int]:
    """Fetch (or read from cache) the PubChem CID associated with each CAS number.

    Prioritizes CAS Common Chemistry: its InChIKey for that CAS number is
    looked up on PubChem directly (an unambiguous structure match), which is
    more precise than searching PubChem by the CAS number as if it were a
    compound name/synonym. Falls back to that name search for any CAS Common
    Chemistry doesn't have a record for, or has no InChIKey for.

    Returns a cas -> cid mapping, omitting any CAS PubChem has no compound for.
    """
    if not unique_cas:
        return {}

    cas_to_cid: dict[str, int] = {
        cas: cid
        for cas in unique_cas
        if (cached := get_cache(cas, subdir=_CAS_CACHE_SUBDIR)) is not None
        and (cid := cached.get("cid")) is not None
    }

    to_fetch = get_missing_keys(unique_cas, subdir=_CAS_CACHE_SUBDIR)
    if to_fetch:
        cas_details = fetch_cas_common_chemistry_details(to_fetch)
        logger.info(
            "Querying %d CAS number(s) on PubChem (%d already cached)...",
            len(to_fetch),
            len(unique_cas) - len(to_fetch),
        )
        for i, cas in enumerate(to_fetch, start=1):
            inchikey = (cas_details.get(cas) or {}).get("inchikey")
            cid, via = _resolve_cas_to_cid(cas, inchikey, i, len(to_fetch))
            # A confirmed empty result is definitive and worth caching; a
            # non-200 or network error may be transient, so it's left
            # uncached and retried on the next run.
            set_cache(cas, {"cid": cid}, subdir=_CAS_CACHE_SUBDIR)
            if cid:
                cas_to_cid[cas] = cid
                logger.info(
                    "  [%d/%d] CAS:%s ... CID:%s (via %s)",
                    i,
                    len(to_fetch),
                    cas,
                    cid,
                    via,
                )
            else:
                logger.info("  [%d/%d] CAS:%s ... not found", i, len(to_fetch), cas)

            if i < len(to_fetch):
                time.sleep(_REQUEST_DELAY)
    elif unique_cas:
        logger.info(
            "All %d CAS→CID mapping(s) already cached, skipping API calls.",
            len(unique_cas),
        )

    return cas_to_cid


def fetch_cids_from_cas(df: pd.DataFrame) -> list[int]:
    """For Molecule rows with an empty CID, look up the CID via CAS on PubChem.
    Updates the CID column in df. Returns the list of newly found CIDs."""
    molecule_df = df[MOLECULE]
    no_cid = molecule_df[CID_COL].isna() | molecule_df[CID_COL].astype(str).apply(
        lambda v: len(parse_cid_cell(v)) == 0
    )
    missing_cid = no_cid & molecule_df[CAS_COL].notna()
    # Typographic dashes (e.g. "2244–16-8") aren't valid in a CAS number but
    # show up in source data; normalize before this is used as a lookup key
    # or sent to PubChem/CAS Common Chemistry.
    cas_series = molecule_df.loc[missing_cid, CAS_COL].astype(str).map(normalize_dashes)

    if cas_series.empty:
        return []

    unique_cas = [str(cas) for cas in cas_series.unique().tolist()]
    cas_to_cid = fetch_cas_to_cid_map(unique_cas)

    if not cas_to_cid:
        return []

    updated = cas_series.astype(str).map(cas_to_cid)
    mask = updated.notna()
    # The CID column holds text (see parse_cid_cell) - pandas 3.0's "str"
    # dtype rejects raw ints on assignment, so stringify before writing back.
    df.loc[updated.index[mask], (MOLECULE, CID_COL)] = (
        updated[mask].astype(int).astype(str).to_numpy()
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
    elif unique_dois:
        logger.info(
            "All %d DOI(s) found in APA cache, skipping API calls.", len(unique_dois)
        )

    return failed
