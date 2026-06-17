import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request

_HTTP_TOO_MANY_REQUESTS = 429
_NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# NCBI allows 3 req/s without an API key; 0.4 s leaves a small margin.
_REQUEST_INTERVAL = 0.4


def _urlopen(url: str, retries: int = 5) -> bytes:
    """Fetch URL with exponential backoff on 429."""
    last_exc: urllib.error.HTTPError | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url) as r:  # noqa: S310
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == _HTTP_TOO_MANY_REQUESTS and attempt < retries - 1:
                time.sleep(2**attempt)
                last_exc = e
            else:
                raise
    raise last_exc  # type: ignore[misc]


def get_taxon_id(species_name: str) -> str | None:
    clean = re.sub(r"\s*\(.*?\)", "", species_name).strip()
    query = urllib.parse.quote(clean)
    url = f"{_NCBI_BASE}/esearch.fcgi?db=taxonomy&term={query}&retmode=json"
    data = json.loads(_urlopen(url))
    time.sleep(_REQUEST_INTERVAL)
    ids = data["esearchresult"]["idlist"]
    return ids[0] if ids else None


def _get_common_name_ncbi(taxon_id: str) -> str:
    url = f"{_NCBI_BASE}/esummary.fcgi?db=taxonomy&id={taxon_id}&retmode=json"
    data = json.loads(_urlopen(url))
    time.sleep(_REQUEST_INTERVAL)
    return data.get("result", {}).get(taxon_id, {}).get("commonname", "")


def _get_common_name_uniprot(taxon_id: str) -> str:
    url = f"https://rest.uniprot.org/taxonomy/{taxon_id}"
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req) as r:  # noqa: S310  # https only
            data = json.loads(r.read())
        return data.get("commonName", "")
    except urllib.error.HTTPError:
        return ""


def _get_common_name(taxon_id: str) -> str:
    return _get_common_name_ncbi(taxon_id) or _get_common_name_uniprot(taxon_id)


species_acronyms = [
    ("Anopheles gambiae", "Agam"),
    ("Apolygus lucorum", "Aluc"),
    ("Bactrocera minax", "Bmax"),
    ("Bombyx mori", "Bmor"),
    ("Dendroctonus ponderosae", "Dpon"),
    ("Drosophila melanogaster", "Dmel"),
    ("Harpegnathos saltator", "Hsal"),
    ("Helicoverpa armigera", "Harm"),
    ("Helicoverpa assulta", "Hass"),
    ("Hylobius abietis", "Habi"),
    ("Ips typographus", "Ityp"),
    ("Lampronia capitella", "Lcap"),
    ("Locusta migratoria", "Lmig"),
    ("Plutella xylostella", "Pxyl"),
    ("Spodoptera frugiperda", "Sfru"),
    ("Spodoptera littoralis", "Slit"),
]


def print_output_species_details() -> None:
    print("species name,acronym,common_name,ncbi_taxon_id")
    for species_name, acronym in species_acronyms:
        taxon_id = get_taxon_id(species_name) or ""
        common_name = _get_common_name(taxon_id) if taxon_id else ""
        print(f"{species_name},{acronym},{common_name},{taxon_id}")


if __name__ == "__main__":
    print_output_species_details()
