# scripts/fetch_blast.py

import json
import re
import time
from io import StringIO
from pathlib import Path

import requests
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML

# Required by NCBI for all Entrez/BLAST requests.
_ENTREZ_EMAIL = "andre.lanrezac@univ-cotedazur.fr"

_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
_BLAST_CACHE_FILE = Path("cache/receptors/blast_cache.json")
_POLL_INTERVAL = 10   # seconds between status polls
_BLAST_TIMEOUT = 300  # give up after 5 minutes
_ENTREZ_DELAY = 0.4   # NCBI policy: max 3 req/s without API key


def _blast_key(sequence: str, species: str) -> str:
    return f"{species.strip().lower()}|{sequence.strip()}"


def _load_blast_cache() -> dict[str, dict]:
    if _BLAST_CACHE_FILE.exists():
        with _BLAST_CACHE_FILE.open(encoding="utf-8") as f:
            return json.load(f)
    return {}


def _save_blast_cache(cache: dict[str, dict]) -> None:
    _BLAST_CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    with _BLAST_CACHE_FILE.open("w", encoding="utf-8") as f:
        json.dump(cache, f, indent=2)


def _fetch_full_sequence(accession: str) -> str | None:
    """Fetch the complete FASTA sequence for a GenBank protein accession via Entrez."""
    Entrez.email = _ENTREZ_EMAIL
    try:
        handle = Entrez.efetch(
            db="protein", id=accession, rettype="fasta", retmode="text"
        )
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq)
    except Exception as e:  # noqa: BLE001
        print(f"    [Entrez] Failed to fetch {accession}: {e}")  # noqa: T201
        return None


def _submit_blast(sequence: str, species: str) -> str | None:
    """Submit a BLASTP job to NCBI. Returns the RID on success."""
    data = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": "nr",
        "QUERY": sequence,
        "ENTREZ_QUERY": f"{species}[Organism]",
        "FORMAT_TYPE": "XML",
        "HITLIST_SIZE": "1",
        "ALIGNMENTS": "1",
    }
    try:
        response = requests.post(_BLAST_URL, data=data, timeout=30)
        response.raise_for_status()
    except requests.RequestException as e:
        print(f"    [BLAST] Submit failed: {e}")  # noqa: T201
        return None

    match = re.search(r"RID = (\w+)", response.text)
    if not match:
        print("    [BLAST] Could not find RID in response.")  # noqa: T201
        return None
    return match.group(1)


def _poll_blast(rid: str) -> str | None:
    """Poll NCBI until the job is READY. Returns the raw XML or None on failure."""
    print(f"    [BLAST] RID={rid}, polling every {_POLL_INTERVAL}s ...", flush=True)  # noqa: T201
    elapsed = 0
    while elapsed < _BLAST_TIMEOUT:
        time.sleep(_POLL_INTERVAL)
        elapsed += _POLL_INTERVAL

        try:
            resp = requests.get(
                _BLAST_URL,
                params={"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"},
                timeout=30,
            )
            resp.raise_for_status()
        except requests.RequestException as e:
            print(f"    [BLAST] Poll error: {e}")  # noqa: T201
            continue

        if "Status=WAITING" in resp.text:
            print(f"    [BLAST] Waiting... ({elapsed}s elapsed)", flush=True)  # noqa: T201
            continue
        if "Status=FAILED" in resp.text:
            print("    [BLAST] Job failed on NCBI side.")  # noqa: T201
            return None
        if "Status=UNKNOWN" in resp.text:
            print("    [BLAST] RID expired or unknown.")  # noqa: T201
            return None
        if "Status=READY" in resp.text:
            result_params = {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_TYPE": "XML",
                "HITLIST_SIZE": "1",
            }
            try:
                result = requests.get(_BLAST_URL, params=result_params, timeout=60)
                result.raise_for_status()
            except requests.RequestException as e:
                print(f"    [BLAST] Failed to retrieve results: {e}")  # noqa: T201
                return None
            else:
                return result.text

    print(f"    [BLAST] Timeout after {_BLAST_TIMEOUT}s.")  # noqa: T201
    return None


def _run_blast_query(sequence: str, species: str) -> tuple[str | None, str | None]:
    """
    Run the full BLAST pipeline.
    Returns (accession, None) on success, (None, error_reason) on failure.
    """
    rid = _submit_blast(sequence, species)
    if rid is None:
        return None, "submission failed"

    xml_text = _poll_blast(rid)
    if xml_text is None:
        return None, "timeout or NCBI error"

    try:
        blast_records = list(NCBIXML.parse(StringIO(xml_text)))
    except ValueError as e:
        return None, f"XML parse error: {e}"

    if not blast_records or not blast_records[0].alignments:
        print(f"    [BLAST] No hit found for species={species!r}")  # noqa: T201
        return None, "no hit found"

    accession = blast_records[0].alignments[0].accession
    print(f"    [BLAST] Best hit: {accession}")  # noqa: T201
    return accession, None


def count_uncached_blast_queries(queries: list[tuple[str, str]]) -> int:
    """Return how many (sequence, species) pairs have never been queried."""
    cache = _load_blast_cache()
    return sum(1 for seq, species in queries if _blast_key(seq, species) not in cache)


def get_failed_blast_queries(
    queries: list[tuple[str, str]],
) -> list[tuple[str, str, str]]:
    """
    Return entries that were previously queried but produced no result.
    Each item is (sequence, species, error_reason).
    """
    cache = _load_blast_cache()
    failed = []
    for seq, species in queries:
        entry = cache.get(_blast_key(seq, species))
        if entry is not None and not entry.get("accession"):
            error = entry.get("error", "unknown error")
            failed.append((seq, species, error))
    return failed


def clear_blast_cache_entries(queries: list[tuple[str, str]]) -> None:
    """Remove cache entries for the given queries so they can be re-run."""
    cache = _load_blast_cache()
    changed = False
    for seq, species in queries:
        key = _blast_key(seq, species)
        if key in cache:
            del cache[key]
            changed = True
    if changed:
        _save_blast_cache(cache)


def fetch_blast_reference(sequence: str, species: str) -> tuple[str, str] | None:
    """
    Resolve a (sequence, species) pair to a (genbank_accession, sequence_ref).

    Results are stored in a single file — cache/blast_cache.json — indexed by
    "{species}|{sequence}".  Negative results (no hit) are also cached so that
    repeated runs never re-query NCBI for the same pair.
    """
    cache = _load_blast_cache()
    key = _blast_key(sequence, species)
    entry = cache.get(key)

    if entry is not None:
        acc = entry.get("accession")
        seq_ref = entry.get("sequence_ref")
        return (acc, seq_ref) if (acc and seq_ref) else None  # None = negative cache

    print(f"    [BLAST] Submitting for species={species!r} ...", flush=True)  # noqa: T201
    accession, error = _run_blast_query(sequence, species)

    if accession is None:
        cache[key] = {"accession": None, "sequence_ref": None, "error": error}
        _save_blast_cache(cache)
        return None

    print(f"    [BLAST] Fetching sequence for {accession} ...", flush=True)  # noqa: T201
    time.sleep(_ENTREZ_DELAY)
    ref_sequence = _fetch_full_sequence(accession)
    if ref_sequence is None:
        return None

    cache[key] = {"accession": accession, "sequence_ref": ref_sequence}
    _save_blast_cache(cache)
    return accession, ref_sequence
