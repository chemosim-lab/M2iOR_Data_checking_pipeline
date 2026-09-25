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

## Getting started

### Prerequisites

- Python >= 3.14
- [uv](https://docs.astral.sh/uv/) — the package/dependency manager this project uses

Install uv (one-time, see the [official guide](https://docs.astral.sh/uv/getting-started/installation/) for alternatives):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

### Installation

From the project root:

```bash
uv sync
```

This creates the `.venv` and installs every dependency pinned in `uv.lock`.

Receptor mutation analysis (`scripts/find_protein_mutations.py`) and the BLAST fallback (`scripts/fetch_blast.py`) require a local BLAST installation:

```bash
./scripts/install_blast.sh          # Linux/macOS
```

```powershell
.\scripts\install_blast.ps1          # Windows (PowerShell)
```

Both fetch NCBI's official prebuilt `blastp` binary and drop it into the venv (`.venv/bin` / `.venv\Scripts`) — no distro/system package, no admin rights. Re-run after any `uv sync` that recreates `.venv` from scratch.

### Running the pipeline

All commands are run through `uv run main.py`, so dependencies are resolved automatically without activating the virtual environment manually.

```bash
uv run main.py --status          # Show the processing status of every known study
uv run main.py --scan            # List Excel files, flag new/modified ones, and choose which to process
uv run main.py --pending         # Process every study still pending
uv run main.py --retry           # Retry only studies that previously failed
uv run main.py --all             # Reprocess every study, ignoring DONE/FAILED status
uv run main.py --file <path>     # Process a single Excel file
uv run main.py --force           # Reprocess every study and force BLAST re-queries
uv run main.py --force <path>    # Reprocess a single file and force BLAST re-queries
uv run main.py --check <path>    # Validate a single file without exporting anything, report every issue
```

#### Check mode and reports

`--check` runs every validation stage on one Excel file and reports all the issues it finds, without writing any CSV, image or registry entry (external API results are still cached). Its exit code is `0` when nothing blocks the export, `1` on validation errors, `2` on a pipeline crash or a usage error.

```bash
uv run main.py --check input/X.xlsx                                  # human-readable report
uv run main.py --check input/X.xlsx --report json                    # JSON report on stdout (logs on stderr)
uv run main.py --check input/X.xlsx --report json --report-file r.json
uv run main.py --file input/X.xlsx --blast cached --report json      # export, non-interactive, with its report
```

Each issue gives its severity, stage, code, sheet, Excel column/rows/cells, value and, when possible, a suggested value. Only `error` issues block the export; `warning` (needs a human look, e.g. rows the Laravel import will skip), `auto_fix` (a value the pipeline silently rewrites, e.g. a receptor name or a CID/CAS it corrects) and `info` never change what the pipeline accepts.

`--blast` controls pending BLAST lookups: `ask` prompts (default, except for `--check`), `cached` never runs one and only applies cached results (default for `--check`), `run` runs every new or retryable lookup without asking.

Note: the input directory (Excel studies) and output directory (M2iOR website data) are currently hardcoded as absolute paths at the top of [main.py](main.py) — update `M2IOR_WEB_PUBLIC_PATH` and `input_path` there if you run the pipeline on a different machine.
