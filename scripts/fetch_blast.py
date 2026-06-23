# scripts/fetch_blast.py

import hashlib
import json
import threading
import time
from pathlib import Path
from typing import Any

from Bio import Blast, Entrez, SeqIO

from scripts.tools.get_common_name import get_taxon_id

SUBSTITUTION_MATRIX = "BLOSUM62"

# Required by NCBI for all Entrez/BLAST requests.
_ENTREZ_EMAIL = "andre.lanrezac@univ-cotedazur.fr"

_BLAST_CACHE_FILE = Path("cache/receptors/blast_cache.json")
_BLAST_XML_DIR = Path("cache/receptors/blast_xml")
_ENTREZ_DELAY = 0.4  # NCBI policy: max 3 req/s without API key
_cache_lock = threading.Lock()

Blast.email = _ENTREZ_EMAIL


def _blast_key(sequence: str, species: str) -> str:
    return f"{species.strip().lower()}|{sequence.strip()}"


def _load_blast_cache() -> dict[str, dict[str, Any]]:
    if _BLAST_CACHE_FILE.exists():
        with _BLAST_CACHE_FILE.open(encoding="utf-8") as f:
            return json.load(f)
    return {}


def _save_blast_cache(cache: dict[str, dict[str, Any]]) -> None:
    _BLAST_CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    with _BLAST_CACHE_FILE.open("w", encoding="utf-8") as f:
        json.dump(cache, f, indent=2)


def _blast_xml_path(key: str) -> Path:
    _BLAST_XML_DIR.mkdir(parents=True, exist_ok=True)
    filename = hashlib.sha256(key.encode()).hexdigest()[:16] + ".xml"
    return _BLAST_XML_DIR / filename


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


def _run_blast_query(
    sequence: str, species: str
) -> tuple[str | None, Path | None, str | None]:
    """
    Run a BLASTP search via Bio.Blast.qblast (submit + poll).
    Returns (accession, xml_path, None) on success, (None, xml_path, error) on failure.
    xml_path is None only when the query never completed.
    """
    taxon_id = get_taxon_id(species)
    msg = f"    [BLAST] Submitting for species={species!r} (txid{taxon_id}) ..."
    print(msg, flush=True)  # noqa: T201
    try:
        result_stream = Blast.qblast(
            "blastp",
            "nr",
            sequence,
            entrez_query=f"txid{taxon_id}[ORGN]",
            expect=0.05,
            hitlist_size=1,
            matrix_name=SUBSTITUTION_MATRIX,
            filter="F",
        )
    except Exception as e:  # noqa: BLE001
        print(f"    [BLAST] Query failed: {e}")  # noqa: T201
        return None, None, "submission failed"

    xml_path = _blast_xml_path(_blast_key(sequence, species))
    with xml_path.open("wb") as out_stream:
        out_stream.write(result_stream.read())
    result_stream.close()

    try:
        blast_record = Blast.read(xml_path)
    except ValueError as e:
        return None, xml_path, f"XML parse error: {e}"

    if not blast_record:
        print(f"    [BLAST] No hit found for species={species!r}")  # noqa: T201
        return None, xml_path, "no hit found"

    accession = blast_record[0].target.id
    print(f"    [BLAST] Best hit: {accession}")  # noqa: T201
    return accession, xml_path, None


def count_uncached_blast_queries(queries: list[tuple[str, str]]) -> int:
    """Return how many (sequence, species) pairs have never been queried."""
    cache = _load_blast_cache()
    return sum(1 for seq, sp in queries if _blast_key(seq, sp) not in cache)


def get_failed_blast_queries(
    queries: list[tuple[str, str]],
) -> list[tuple[str, str, str]]:
    """
    Return queries that were previously submitted but produced no result.
    Each item is (sequence, species, error_reason).
    """
    cache = _load_blast_cache()
    failed: list[tuple[str, str, str]] = []
    for seq, sp in queries:
        entry: dict[str, Any] | None = cache.get(_blast_key(seq, sp))
        if entry is not None and not entry.get("accession"):
            error = str(entry.get("error", "unknown error"))
            failed.append((seq, sp, error))
    return failed


def get_successful_blast_queries(
    queries: list[tuple[str, str]],
) -> list[tuple[str, str]]:
    """Return (sequence, species) pairs that have a successful cache entry."""
    cache = _load_blast_cache()
    return [
        (seq, sp)
        for seq, sp in queries
        if (entry := cache.get(_blast_key(seq, sp))) is not None
        and entry.get("accession")
    ]


def clear_blast_cache_entries(queries: list[tuple[str, str]]) -> None:
    """Remove cache entries for the given queries so they can be re-run."""
    cache = _load_blast_cache()
    changed = False
    for seq, sp in queries:
        key = _blast_key(seq, sp)
        if key in cache:
            del cache[key]
            changed = True
    if changed:
        _save_blast_cache(cache)


def fetch_blast_reference(sequence: str, species: str) -> tuple[str, str] | None:
    """
    Resolve a (sequence, species) pair to a (genbank_accession, sequence_ref).

    Results are stored in cache/blast_cache.json indexed by "{species}|{sequence}".
    Negative results (no hit) are also cached so repeated runs never re-query NCBI.
    """
    key = _blast_key(sequence, species)

    with _cache_lock:
        cache = _load_blast_cache()
        entry = cache.get(key)

    if entry is not None:
        acc = entry.get("accession")
        seq_ref = entry.get("sequence_ref")
        return (acc, seq_ref) if (acc and seq_ref) else None

    accession, xml_path, error = _run_blast_query(sequence, species)

    if accession is None:
        with _cache_lock:
            cache = _load_blast_cache()
            cache[key] = {
                "accession": None,
                "sequence_ref": None,
                "xml_path": str(xml_path) if xml_path else None,
                "error": error,
            }
            _save_blast_cache(cache)
        return None

    print(f"    [BLAST] Fetching sequence for {accession} ...", flush=True)  # noqa: T201
    time.sleep(_ENTREZ_DELAY)
    ref_sequence = _fetch_full_sequence(accession)
    if ref_sequence is None:
        return None

    with _cache_lock:
        cache = _load_blast_cache()
        cache[key] = {
            "accession": accession,
            "sequence_ref": ref_sequence,
            "xml_path": str(xml_path),
            "error": None,
        }
        _save_blast_cache(cache)
    return accession, ref_sequence
