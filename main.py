# pipeline/main.py

import argparse
from pathlib import Path
from typing import TYPE_CHECKING, cast

from scripts.columns import (
    DATABASE,
    GENE_NAME,
    IDENTITY,
    RECEPTOR_NAME,
    SEQUENCE,
    SEQUENCE_REF,
    UNIPROT_ID,
)
from scripts.export_csv import export_to_csv
from scripts.fetch_data import (
    fetch_apa_references,
    fetch_cids_from_cas,
    fetch_genbank_data,
    fetch_ncbi_data,
    fetch_pubchem_data,
    fetch_uniprot_data,
)
from scripts.process_molecules import (
    enrich_molecule_columns,
    get_unique_cids,
    validate_cid_or_cas,
)
from scripts.process_receptors import (
    add_empty_column_after,
    collect_blast_queries,
    enrich_with_reference_and_mutations,
    get_unique_uniprot_ids,
    process_receptors_name_columns,
    rename_column,
    strip_column,
)
from scripts.process_responses import (
    validate_concentration_column,
    validate_parameter_column,
    validate_responsive_column,
    validate_value_column,
)
from scripts.process_sources import (
    enrich_reference_column,
    get_unique_dois,
    normalize_doi_column,
    validate_doi_column,
)
from scripts.read_excel import get_raw_data_from_excel_file
from scripts.registry import ProcessingStatus, StudyFileTracker, StudyProcessingRegistry

if TYPE_CHECKING:
    from pandas import DataFrame

LARAVEL_DATA_PATH = Path(
    "/home/andre/Dev/M2iOR_migration/M2iOR/M2iOR_web_public-main/resources/data"
)
REGISTRY_FILE = Path("./registry.json")


def process_raw_excel_file(excel_path: Path) -> None:
    study_id: str = excel_path.stem  # e.g. "001_Kreher_Neuron_2005"
    output_path = LARAVEL_DATA_PATH
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

        # ----------------------------------------------------------------------
        # RECEPTORS AND CO-RECEPTORS -------------------------------------------
        protein_groups = ["Receptor", "Co-Receptor"]
        all_unique_uniprot_ids: list[str] = get_unique_uniprot_ids(
            df, groups=protein_groups, column_name="UniProt ID"
        )

        # Fetch Uniprot data from accession number (UniprotID) and store them
        # in the cache
        failed_uids = fetch_uniprot_data(all_unique_uniprot_ids)
        # NCBI fallback
        failed_uids = fetch_ncbi_data(all_unique_uniprot_ids, failed_uids)

        # Rename "Gene Name" columns to "Receptor Name"
        rename_column(
            df,
            protein_groups,
            old_column_name=GENE_NAME,
            new_column_name=RECEPTOR_NAME,
        )

        strip_column(df, protein_groups, UNIPROT_ID)

        # For rows without a UniProt ID, fall back to BLAST against NCBI nr.
        # Deduplicate queries first — one API call per unique (sequence, species) pair.
        blast_queries = collect_blast_queries(
            df, protein_groups, fallback_uids=failed_uids
        )
        blast_refs: dict[str, str] = fetch_genbank_data(
            df, protein_groups, blast_queries
        )

        # BLAST may have written new GenBank accessions into the UniProt ID column.
        # Re-collect all UIDs so that BLAST-discovered accessions are also named.
        all_unique_ids_after_blast: list[str] = get_unique_uniprot_ids(
            df, groups=protein_groups, column_name="UniProt ID"
        )
        blast_only_accessions = [
            uid
            for uid in all_unique_ids_after_blast
            if uid not in set(all_unique_uniprot_ids)
        ]
        if blast_only_accessions:
            fetch_ncbi_data(blast_only_accessions, blast_only_accessions)

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
            df, protein_groups, all_unique_ids_after_blast, blast_refs
        )

        # Use fetched data to check "Receptor Name" columns for Receptor and Co-Receptor
        process_receptors_name_columns(df, protein_groups, all_unique_ids_after_blast)

        # ----------------------------------------------------------------------
        # MOLECULES ------------------------------------------------------------
        validate_cid_or_cas(df)
        all_unique_cids: list[int] = get_unique_cids(df)
        failed_cids = fetch_pubchem_data(all_unique_cids)

        # Fallback: rows with empty CID but a CAS → find the CID via PubChem
        new_cids = fetch_cids_from_cas(df)
        if new_cids:
            fetch_pubchem_data(new_cids)

        enrich_molecule_columns(df)

        # ----------------------------------------------------------------------
        # RESPONSES ------------------------------------------------------------
        validate_responsive_column(df)
        validate_parameter_column(df)
        validate_value_column(df)
        validate_concentration_column(df)

        # ----------------------------------------------------------------------
        # ASSAY ----------------------------------------------------------------

        # ----------------------------------------------------------------------
        # SOURCE ---------------------------------------------------------------
        normalize_doi_column(df)
        validate_doi_column(df)
        unique_dois: list[str] = get_unique_dois(df)
        fetch_apa_references(unique_dois)
        enrich_reference_column(df)

        # Export enriched dataset as CSV with two-row (group / column) header
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
        nargs="?",
        const="__all__",
        metavar="FILE",
        help="Reprocess files, ignoring DONE/FAILED status. Without argument: reprocess all. With FILE: reprocess only that file.",
    )
    args = parser.parse_args()

    # Typed extraction of arguments (argparse returns a Namespace whose attributes are typed as Any)
    show_status = cast(bool, args.status)
    file_arg = cast(str | None, args.file)
    do_pending = cast(bool, args.pending)
    do_retry = cast(bool, args.retry)
    force_arg = cast(str | None, args.force)

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
    elif force_arg is not None:
        if force_arg == "__all__":
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
            target_path = Path(force_arg)
            study_id = target_path.stem
            if study_id in registry:
                tracker = registry.get(study_id)
                tracker.status = ProcessingStatus.PENDING
                tracker.error_message = None
            process_raw_excel_file(target_path)
    else:
        parser.print_help()
