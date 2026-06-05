# pipeline/main.py

import argparse
from pathlib import Path
from typing import TYPE_CHECKING, cast

from scripts.columns import (
    DATABASE,
    IDENTITY,
    SEQUENCE,
    SEQUENCE_REF,
    UNIPROT_ID,
)
from scripts.export_csv import export_to_csv
from scripts.fetch_blast import (
    clear_blast_cache_entries,
    count_uncached_blast_queries,
    get_failed_blast_queries,
)
from scripts.fetch_data import fetch_genbank_data, fetch_uniprot_data
from scripts.process_receptors import (
    add_empty_column_after,
    collect_blast_queries,
    enrich_with_reference_and_mutations,
    get_unique_uniprot_ids,
    resolve_missing_ids_via_blast,
)
from scripts.read_excel import get_raw_data_from_excel_file
from scripts.registry import ProcessingStatus, StudyFileTracker, StudyProcessingRegistry

if TYPE_CHECKING:
    from pandas import DataFrame

LARAVEL_DATA_PATH = Path(
    "/home/andre/Dev/M2iOR_migration/M2iOR/M2iOR_web_public-main/input_test"
)
REGISTRY_FILE = Path("./registry.json")


def process_raw_excel_file(excel_path: Path) -> None:
    study_id: str = excel_path.stem  # e.g. "001_Kreher_Neuron_2005"
    output_path = LARAVEL_DATA_PATH / study_id
    output_path.mkdir(parents=True, exist_ok=True)

    # Register if not already known; ignore the error if already registered
    try:
        _ = registry.register(study_id, excel_path)
    except ValueError:
        pass

    study_tracker: StudyFileTracker = registry.get(study_id)

    try:
        # 1. Read and clean the Excel file
        df: DataFrame = get_raw_data_from_excel_file(excel_path)

        protein_groups = ["Receptor", "Co-Receptor"]
        all_unique_uniprot_ids: list[str] = get_unique_uniprot_ids(
            df, groups=protein_groups, column_name="UniProt ID"
        )

        # Fetch Uniprot data from accession number (UniprotID) and store them
        # in the cache
        fetch_uniprot_data(all_unique_uniprot_ids)

        # For rows without a UniProt ID, fall back to BLAST against NCBI nr.
        # Deduplicate queries first — one API call per unique (sequence, species) pair.
        blast_queries = collect_blast_queries(df, protein_groups)
        blast_refs: dict[str, str] = fetch_genbank_data(
            df, protein_groups, blast_queries
        )

        # Add empty Sequence_ref in df for Receptor and Co-Receptor groups
        add_empty_column_after(
            df, protein_groups, after_column=SEQUENCE, new_column=SEQUENCE_REF
        )
        add_empty_column_after(
            df, protein_groups, after_column=UNIPROT_ID, new_column=DATABASE
        )
        add_empty_column_after(
            df, protein_groups, after_column=DATABASE, new_column=IDENTITY
        )

        # Add Reference sequences, identity (%) and mutations vs UniProt reference
        enrich_with_reference_and_mutations(
            df, protein_groups, all_unique_uniprot_ids, blast_refs
        )

        # 4. Export enriched dataset as CSV with two-row (group / column) header
        export_to_csv(df, output_path / f"{study_id}.csv")

        study_tracker.complete()
        print(f"✓ {study_id}")

    except ValueError as er:
        print(f"✗ {study_id} — FAILED: {er}")
        study_tracker.fail(error_message=str(er))

    registry.save(REGISTRY_FILE)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="M2iOR processing pipeline")
    _ = parser.add_argument(
        "--status", action="store_true", help="Show processing status"
    )
    _ = parser.add_argument("--file", help="Process a single Excel file")
    _ = parser.add_argument(
        "--pending", action="store_true", help="Process all pending Excel files"
    )
    _ = parser.add_argument(
        "--retry",
        action="store_true",
        help="Retry only failed files",
    )
    _ = parser.add_argument(
        "--force",
        action="store_true",
        help="Reprocess all files, including already processed ones (DONE)",
    )
    args = parser.parse_args()

    # Typed extraction of arguments (argparse returns a Namespace whose attributes are typed as Any)
    show_status = cast(bool, args.status)
    file_arg = cast(str | None, args.file)
    do_pending = cast(bool, args.pending)
    do_retry = cast(bool, args.retry)
    do_force = cast(bool, args.force)

    input_path = Path("input")

    print("Loading processing registry...")
    registry = StudyProcessingRegistry.load(REGISTRY_FILE)
    print("Done")

    if show_status:
        print(registry.summary)
    elif file_arg is not None:
        process_raw_excel_file(Path(file_arg))
    elif do_retry:
        failed_jobs: list[StudyFileTracker] = registry.all_failed()
        for job in failed_jobs:
            job.status = ProcessingStatus.PENDING
            job.error_message = None
            process_raw_excel_file(job.file_path)
    elif do_pending:
        found_files: list[StudyFileTracker] = registry.all_pending()
        print(f"Found {len(found_files)} files")
        _ = input("Press Enter to continue")
        for excel_file in found_files:
            print(f"File: {excel_file.file_path}")
            process_raw_excel_file(excel_file.file_path)
    elif do_force:
        # Reset all existing trackers to PENDING before reprocessing
        for tracker in registry.all_done() + registry.all_failed():
            tracker.status = ProcessingStatus.PENDING
            tracker.error_message = None
        found_file_paths = sorted(input_path.glob("**/*.xlsx"))
        print(f"Found {len(found_file_paths)} files — forcing full reprocess")
        _ = input("Press Enter to continue")
        for excel_file_path in found_file_paths:
            print(f"File: {excel_file_path}")
            process_raw_excel_file(excel_file_path)
    else:
        parser.print_help()
