import json
import re
import urllib.parse
import urllib.request


def get_taxon_id(species_name: str) -> str | None:
    clean = re.sub(r"\s*\(.*?\)", "", species_name).strip()
    query = urllib.parse.quote(clean)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={query}&retmode=json"
    with urllib.request.urlopen(url) as r:
        data = json.loads(r.read())
    ids = data["esearchresult"]["idlist"]
    return ids[0] if ids else None


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
    for species, acronym in species_acronyms:
        url = f"https://rest.uniprot.org/taxonomy/search?query={urllib.parse.quote(species)}&fields=id,scientific_name,common_name&format=json&size=1"

        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req) as r:
            data = json.loads(r.read())

        hits = data.get("results", [])
        if hits:
            taxon_id = hits[0].get("taxonId", "")
            common_name = hits[0].get("commonName", "")
        else:
            taxon_id, common_name = "", ""

        print(f"{species},{acronym},{common_name},{taxon_id}")
