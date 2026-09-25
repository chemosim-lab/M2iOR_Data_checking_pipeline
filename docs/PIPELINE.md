# Pipeline details

The pipeline reads one Excel file per scientific study (M2iOR double-header template), checks and enriches its data through several stages (receptors, molecules, responses, sources), then exports a final CSV together with chemical structure images for the M2iOR website. The processing state of each study is tracked so that already-successful studies are not reprocessed.

## Directory tree

```
M2iOR_Data_checking_pipeline/
├── main.py                              # Entry point: CLI and pipeline orchestration
├── registry.json                        # Persistent processing state for each study
├── excel_template.json                  # Reference schema for the expected Excel model
│
└── scripts/                              # Core business logic of the pipeline
    ├── columns.py                       # Expected column/group vocabulary for the Excel file
    ├── read_excel.py                    # Reading and structural validation of the Excel file
    ├── normalize_dataframe.py           # Generic data cleanup (whitespace, casing, renaming)
    ├── process_receptors.py             # Validation and enrichment of receptors/co-receptors
    ├── process_molecules.py             # Validation and enrichment of tested molecules
    ├── molecule_stereo.py               # Stereochemistry classification and image generation
    ├── process_responses.py             # Validation of experimental result columns
    ├── process_sources.py               # Validation and enrichment of bibliographic references
    ├── fetch_data.py                    # Calls to external APIs (UniProt, NCBI, PubChem, doi.org)
    ├── fetch_blast.py                   # Remote BLAST search to identify an unknown sequence
    ├── find_protein_mutations.py        # Sequence alignment and mutation/identity computation
    ├── registry.py                      # Tracks the processing state of each study
    ├── report.py                        # Structured validation issues and the check report
    ├── export_csv.py                    # Generates the final CSV for the website
    ├── cache_manager.py                 # Generic disk cache for external data
    ├── extract_receptors.py             # Standalone utility to extract the list of unique receptors
    ├── install_blast.sh                 # Installs the local BLAST program
    └── tools/                           # Standalone tools, outside the main flow
        ├── split_by_doi.py              # Splits a multi-study Excel file into one file per study
        ├── backfill_stereo_images.py    # Regenerates structure images without rerunning the full pipeline
        └── get_common_name.py           # Translates a species name into an NCBI taxonomy identifier
```

*(Only the files forming the functional core of the pipeline are listed here — tests, caches, the virtual environment, and ancillary development tools are omitted as they are not relevant to understanding the project.)*

## Key features

### 1. Orchestration and processing tracking
- **`main.py`** — drives the whole pipeline: command-line options to scan new files, retry failures or force reprocessing, and the per-study processing loop.
- **`scripts/registry.py`** — records for each study whether it is pending, completed, or failed, and detects files modified since their last successful run.
- **`scripts/report.py`** — validation steps raise a `ValidationError` carrying structured issues (sheet, Excel row/cell, column, value, suggestion) and report non-blocking findings through `emit()`. The pipeline runs stage by stage (structure, receptors, molecules, responses, assay, source, final checks), so one run reports every problem instead of stopping at the first. `main.py --check` prints the resulting report without exporting anything.

### 2. Excel ingestion and structural validation
- **`scripts/columns.py`** — defines the expected schema (Receptor/Co-receptor/Molecule/Response/Source groups...) and their columns.
- **`scripts/read_excel.py`** — reads the Excel template, rebuilds the double-header table, and checks that its structure matches the expected model.
- **`scripts/normalize_dataframe.py`** — cleans values (extra whitespace, casing) and renames legacy columns to their current names.

### 3. Receptors and co-receptors
- **`scripts/process_receptors.py`** — collects identifiers (UniProt/GenBank), checks the consistency of the declared species, and flags receptors whose identity with their reference sequence is too low. Receptor names are only ever taken from the Excel file, never from the UniProt/NCBI record: a differently cased `OR5` is reformatted to `Or5`, but an empty name, a name that isn't an `Or<N>`/`Orco` name, or different names sharing an accession typed in the Excel file block the export.
- **`scripts/fetch_data.py`** — fetches UniProt data, falling back to NCBI data, for each identifier.
- **`scripts/fetch_blast.py`** — species-targeted BLAST search against NCBI databases to recover an identifier from an unknown sequence.
- **`scripts/find_protein_mutations.py`** — aligns each sequence against its reference, computing mutations and percent identity.
- **`scripts/install_blast.sh`** — installs the local BLAST program required for this computation.
- **`scripts/tools/get_common_name.py`** — translates a species name into an NCBI taxonomy identifier, used to target BLAST searches.

### 4. Molecules
- **`scripts/process_molecules.py`** — checks each molecule's identifier (PubChem CID or CAS number), enriches it (name, formula, structure), and detects mixtures of compounds.
- **`scripts/molecule_stereo.py`** — stereochemistry classification via RDKit (pure molecule vs. mixture of isomers) and generation of an annotated structure image.
- **`scripts/tools/backfill_stereo_images.py`** — regenerates these images independently of a full pipeline run.
- **`scripts/fetch_data.py`** — retrieves PubChem records and resolves a CID from a CAS number.

### 5. Experimental responses
- **`scripts/process_responses.py`** — checks that result columns (response, parameter, value, unit, technique...) contain only allowed values.

### 6. Bibliographic sources
- **`scripts/process_sources.py`** — validates DOI format and replaces the reference with its full bibliographic version.
- **`scripts/fetch_data.py`** — queries doi.org to obtain this formatted reference.

### 7. Cache
- **`scripts/cache_manager.py`** — generic disk cache that avoids re-querying the same APIs on every run.
- **`scripts/fetch_blast.py`** — additionally maintains its own dedicated cache for BLAST results.

### 8. Final export
- **`scripts/export_csv.py`** — writes the final, enriched and validated table to CSV, with column names matching what the M2iOR website expects.

### 9. Ancillary tools (one-off usage, outside the main flow)
- **`scripts/tools/split_by_doi.py`** — splits a multi-study Excel file into one file per study, ahead of the pipeline.
- **`scripts/extract_receptors.py`** — extracts the list of unique receptors (standalone utility).
