# M2iOR – Data Checking Pipeline

Pipeline that takes the Excel files of scientific studies (olfactory receptors, molecules, responses) as input and produces checked and enriched CSV files for the M2iOR website.

```mermaid
flowchart TD
    Root["M2iOR Pipeline"]

    Root --> Ingestion["Ingestion & validation<br/>Excel reading, structure check, cleanup"]
    Root --> Receptors["Receptors<br/>UniProt/NCBI identification, BLAST, mutations"]
    Root --> Molecules["Molecules<br/>PubChem enrichment, stereochemistry, images"]
    Root --> Responses["Experimental responses<br/>validation of measured values"]
    Root --> Sources["Bibliographic sources<br/>DOI validation, APA reference"]
    Root --> Export["Final export<br/>CSV + images for the website"]
    Root --> Cache["Cache<br/>avoids re-querying the same APIs"]
```

Full file-by-file breakdown: [docs/PIPELINE.md](docs/PIPELINE.md)
