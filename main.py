# pipeline/main.py

import argparse
import contextlib
import json
import logging
import sys
import tempfile
import warnings
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, cast

import colorlog
import pandas as pd
from pandas.core.frame import DataFrame

# The Excel template's column order (Order, Species, Receptor Name, ...) is
# deliberately not alphabetical, so the columns MultiIndex is never lexsorted -
# pandas' fast binary-search path for (group, column) indexing doesn't apply,
# and it warns on every such access even though results are still correct.
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

from scripts.columns import (
    ACCESSION,
    ASSAY,
    CANONICAL_NAME,
    CAS,
    CID,
    CO_RECEPTOR,
    DATABASE,
    DOI,
    EXPERIMENTAL_TECHNIQUE,
    GENE_NAME,
    IDENTITY,
    IDENTITY_SHORT,
    MOLECULE,
    MOLECULE_NAME,
    PARAMETER,
    RECEPTOR,
    RECEPTOR_NAME,
    REFERENCE,
    RESPONSE,
    SEQUENCE,
    SEQUENCE_REF,
    SOURCE,
    UNIPROT_ID,
)
from scripts.export_csv import export_to_csv
from scripts.fetch_data import (
    BlastMode,
    fetch_apa_references,
    fetch_cas_common_chemistry_details,
    fetch_cas_to_cid_map,
    fetch_cids_from_cas,
    fetch_genbank_data,
    fetch_ncbi_data,
    fetch_pubchem_data,
    fetch_uniprot_data,
)
from scripts.normalize_dataframe import normalize_df, rename_column
from scripts.process_molecules import (
    get_unique_cas_for_cross_check,
    get_unique_cids,
    parse_cid_cell,
    process_molecules,
    validate_cid_or_cas,
)
from scripts.process_receptors import (
    NOT_AVAILABLE,
    add_empty_column_after,
    collect_blast_queries,
    enrich_with_reference_and_mutations,
    get_unique_accessions,
    process_receptors_name_columns,
)
from scripts.process_responses import (
    validate_concentration_column,
    validate_concentration_unit_column,
    validate_experimental_technique_column,
    validate_parameter_column,
    validate_responsive_column,
    validate_stimulation_duration_unit_column,
    validate_stimulation_flux_unit_column,
    validate_unit_column,
    validate_value_column,
    validate_value_nature_column,
)
from scripts.process_sources import (
    enrich_reference_column,
    get_unique_dois,
    normalize_doi_column,
    validate_doi_column,
)
from scripts.read_excel import get_raw_data_from_excel_file
from scripts.registry import (
    ProcessingStatus,
    StudyFileTracker,
    StudyProcessingRegistry,
    compute_file_hash,
)
from scripts.report import (
    ExcelLayout,
    Issue,
    IssueCollector,
    ValidationError,
    build_report,
    emit,
    format_text_report,
)

if TYPE_CHECKING:
    from pandas import DataFrame

# The M2iOR website project is expected next to this pipeline's directory
M2IOR_WEB_PUBLIC_PATH = Path(__file__).resolve().parent.parent / "M2iOR_web_public-main"
M2IOR_DATA_INPUT_PATH = M2IOR_WEB_PUBLIC_PATH / "resources" / "data"
M2IOR_MOLECULE_IMAGES_PATH = M2IOR_WEB_PUBLIC_PATH / "public" / "images" / "molecules"
REGISTRY_FILE = Path("./registry.json")

PROTEIN_GROUPS = [RECEPTOR, CO_RECEPTOR]

PIPELINE_STAGES = [
    "structure",
    "receptors",
    "molecules",
    "responses",
    "assay",
    "source",
    "final_checks",
]

# Cells the Laravel import (M2iOR_web_public-main/app/Services/CsvImportService.php)
# requires on every row: a row missing one of them is skipped there, not
# imported. Reported as warnings, since the pipeline itself accepts them.
_IMPORT_REQUIRED_COLUMNS: list[tuple[str, str]] = [
    (RECEPTOR, ACCESSION),
    (RECEPTOR, SEQUENCE),
    (RECEPTOR, SEQUENCE_REF),
    (CO_RECEPTOR, ACCESSION),
    (CO_RECEPTOR, SEQUENCE),
    (CO_RECEPTOR, SEQUENCE_REF),
    (MOLECULE, CID),
    (ASSAY, EXPERIMENTAL_TECHNIQUE),
    (RESPONSE, PARAMETER),
]

logger = colorlog.getLogger(__name__)


@dataclass
class PipelineRun:
    df: DataFrame | None
    collector: IssueCollector
    layout: ExcelLayout | None
    exported_csv: Path | None = None

    @property
    def exit_code(self) -> int:
        """0: no blocking error; 1: validation errors; 2: pipeline crash."""
        if self.collector.crashed:
            return 2
        if self.df is None or self.collector.has_errors:
            return 1
        return 0


def _is_empty(values: pd.Series) -> pd.Series:
    # fillna before astype(str): pandas 3.0 leaves NaN as a float otherwise.
    return values.isna() | values.fillna("").astype(str).str.strip().isin(["", "nan"])


def _report_empty_rows(raw: DataFrame) -> None:
    for idx in raw.index[raw.isna().all(axis=1)]:
        emit(
            Issue(
                severity="warning",
                code="empty_row",
                message="Row is entirely empty.",
                rows=[idx],
            )
        )


def _reformat(df: DataFrame) -> DataFrame:
    df = normalize_df(df)
    rename_column(df, old_column_name=GENE_NAME, new_column_name=RECEPTOR_NAME)
    rename_column(df, old_column_name=UNIPROT_ID, new_column_name=ACCESSION)

    # Add empty Sequence_ref in df for Receptor and Co-Receptor groups
    add_empty_column_after(
        df, PROTEIN_GROUPS, after_column=SEQUENCE, new_column=SEQUENCE_REF
    )
    add_empty_column_after(
        df, PROTEIN_GROUPS, after_column=ACCESSION, new_column=DATABASE
    )
    add_empty_column_after(
        df, PROTEIN_GROUPS, after_column=DATABASE, new_column=IDENTITY
    )
    add_empty_column_after(
        df, PROTEIN_GROUPS, after_column=IDENTITY, new_column=IDENTITY_SHORT
    )
    add_empty_column_after(
        df, [MOLECULE], after_column=MOLECULE_NAME, new_column=CANONICAL_NAME
    )
    return df


def _report_unresolved_accessions(df: DataFrame, failed_accessions: list[str]) -> None:
    failed = {a for a in failed_accessions if a.lower() != NOT_AVAILABLE.lower()}
    for group in PROTEIN_GROUPS:
        accessions = df[group][ACCESSION].astype(str).str.strip()
        for idx in df.index[accessions.isin(failed)]:
            emit(
                Issue(
                    severity="warning",
                    code="accession_not_found",
                    message="Accession found in neither UniProt nor NCBI.",
                    group=group,
                    column=ACCESSION,
                    rows=[idx],
                    value=accessions[idx],
                )
            )


def _process_receptors(
    df: DataFrame, *, blast_mode: BlastMode, force_blast: bool
) -> None:
    unique_accessions: list[str] = get_unique_accessions(
        df, groups=PROTEIN_GROUPS, column_name=ACCESSION
    )

    # Fetch Uniprot data from accession number (UniprotID) and store them
    # in the cache
    failed_accessions: list[str] = fetch_uniprot_data(unique_accessions)
    # NCBI fallback
    failed_accessions = fetch_ncbi_data(unique_accessions, failed_accessions)
    _report_unresolved_accessions(df, failed_accessions)

    # For rows without a UniProt or NCBI accession but only have a sequence,
    # fall back to BLAST against NCBI nr_cluster_seq.
    # Deduplicate queries first — one API call per unique (sequence, species) pair.
    blast_queries: list[tuple[str, str]] = collect_blast_queries(
        df, PROTEIN_GROUPS, fallback_accessions=failed_accessions
    )

    # Accessions as typed in the Excel file, before BLAST fills empty ones.
    excel_accessions = {group: df[group][ACCESSION].copy() for group in PROTEIN_GROUPS}
    refseq_by_accession: dict[str, str] = fetch_genbank_data(
        df, PROTEIN_GROUPS, blast_queries, force_blast=force_blast, mode=blast_mode
    )

    # BLAST may have written new GenBank accessions into the 'Accession' column.
    # Re-collect all UIDs so that BLAST-discovered accessions are also named.
    all_unique_ids_after_blast: list[str] = get_unique_accessions(
        df, groups=PROTEIN_GROUPS, column_name=ACCESSION, verbose=False
    )
    blast_only_accessions = [
        uid for uid in all_unique_ids_after_blast if uid not in set(unique_accessions)
    ]
    if blast_only_accessions:
        fetch_ncbi_data(blast_only_accessions, blast_only_accessions)

    # Add Reference sequences, identity (%) and mutations vs UniProt reference
    enrich_with_reference_and_mutations(
        df, PROTEIN_GROUPS, all_unique_ids_after_blast, refseq_by_accession
    )

    # Use fetched data to check "Receptor Name" columns for Receptor and Co-Receptor
    process_receptors_name_columns(
        df, PROTEIN_GROUPS, all_unique_ids_after_blast, excel_accessions
    )


def _report_cid_lookups(
    df: DataFrame, cids_before: pd.Series, failed_cids: list[int]
) -> None:
    failed = set(failed_cids)
    cids_after = df[MOLECULE][CID]
    for idx in df.index:
        before = parse_cid_cell(cids_before[idx])
        if failed.intersection(before):
            emit(
                Issue(
                    severity="warning",
                    code="cid_not_found",
                    message="CID not found on PubChem.",
                    group=MOLECULE,
                    column=CID,
                    rows=[idx],
                    value=cids_before[idx],
                )
            )
        elif not before and parse_cid_cell(cids_after[idx]):
            emit(
                Issue(
                    severity="info",
                    code="cid_filled_from_cas",
                    message="Empty CID filled from the CAS number.",
                    group=MOLECULE,
                    column=CID,
                    rows=[idx],
                    suggested_value=cids_after[idx],
                    details={"cas": df[MOLECULE].loc[idx, CAS]},
                )
            )


def _process_molecules(
    df: DataFrame, collector: IssueCollector, synonyms_dir: Path, images_dir: Path
) -> None:
    collector.check(validate_cid_or_cas, df)
    all_unique_cids: list[int] = get_unique_cids(df)
    failed_cids = fetch_pubchem_data(all_unique_cids)

    # Fallback: rows with empty CID but a CAS → find the CID via PubChem
    cids_before = df[MOLECULE][CID].copy()
    new_cids = fetch_cids_from_cas(df)
    if new_cids:
        fetch_pubchem_data(new_cids)
    _report_cid_lookups(df, cids_before, failed_cids)

    # Cross-check: rows with both a (single) CID and CAS(es) → resolve each
    # CAS to a CID (prioritizing CAS Common Chemistry, see fetch_data.py)
    # so process_molecules can verify they agree with the table's CID.
    unique_cas_for_check = get_unique_cas_for_cross_check(df)
    cas_to_cid_map = fetch_cas_to_cid_map(unique_cas_for_check)
    cas_details_map = fetch_cas_common_chemistry_details(unique_cas_for_check)
    extra_cids = [cid for cid in cas_to_cid_map.values() if cid not in all_unique_cids]
    if extra_cids:
        fetch_pubchem_data(extra_cids)

    process_molecules(
        df,
        synonyms_dir,
        images_dir,
        cas_to_cid_map=cas_to_cid_map,
        cas_details_map=cas_details_map,
    )


def _process_sources(df: DataFrame) -> None:
    normalize_doi_column(df)
    validate_doi_column(df)
    unique_dois: list[str] = get_unique_dois(df)
    failed_dois = set(fetch_apa_references(unique_dois))
    dois = df[SOURCE][DOI].astype(str).str.strip()
    for idx in df.index[dois.isin(failed_dois)]:
        emit(
            Issue(
                severity="warning",
                code="doi_not_resolved",
                message="doi.org returned no reference for this DOI; the Reference "
                "column is left as is.",
                group=SOURCE,
                column=DOI,
                rows=[idx],
                value=dois[idx],
            )
        )
    enrich_reference_column(df)


def _check_row_count(df: DataFrame, input_row_count: int) -> None:
    # Safety net: no processing step should ever grow the number of rows
    # (e.g. an unintended merge/cross-product) — the CSV must have at
    # most as many data rows as the Excel input had.
    output_row_count = len(df)
    if output_row_count > input_row_count:
        msg = (
            f"Row count mismatch: CSV output would have {output_row_count} "
            f"rows, more than the {input_row_count} rows read from the "
            "Excel input."
        )
        raise ValidationError(msg, [Issue(code="row_count_grew", message=msg)])


def _report_import_blockers(df: DataFrame) -> None:
    """Flag rows the Laravel import would skip for a missing required value."""
    for group, column in _IMPORT_REQUIRED_COLUMNS:
        empty = _is_empty(df[group][column])
        for idx in df.index[empty]:
            emit(
                Issue(
                    severity="warning",
                    code="import_required_value_missing",
                    message=f"'{column}' is empty after processing: the Laravel "
                    "import skips this row.",
                    group=group,
                    column=column,
                    rows=[idx],
                )
            )
    no_reference = _is_empty(df[SOURCE][DOI]) & _is_empty(df[SOURCE][REFERENCE])
    for idx in df.index[no_reference]:
        emit(
            Issue(
                severity="warning",
                code="import_required_value_missing",
                message="Both DOI and Reference are empty: the Laravel import "
                "skips this row.",
                group=SOURCE,
                column=DOI,
                rows=[idx],
            )
        )


def run_pipeline(
    excel_path: Path,
    *,
    blast_mode: BlastMode,
    synonyms_dir: Path,
    images_dir: Path,
    force_blast: bool = False,
    catch_unexpected: bool = False,
) -> PipelineRun:
    """Read, validate and enrich one study, stage by stage.

    A stage's failure is recorded and the next independent stages still run,
    so a single run reports every problem. Nothing is exported here."""
    collector = IssueCollector(catch_unexpected=catch_unexpected)
    df: DataFrame | None = None
    layout: ExcelLayout | None = None

    with collector.stage("structure"):
        raw = get_raw_data_from_excel_file(excel_path)
        layout = ExcelLayout.from_raw(raw)
        _report_empty_rows(raw)
        df = _reformat(raw)

    if df is None:
        for stage in PIPELINE_STAGES[1:]:
            collector.skip(stage)
        return PipelineRun(df=None, collector=collector, layout=layout)

    input_row_count: int = len(df)

    with collector.stage("receptors"):
        _process_receptors(df, blast_mode=blast_mode, force_blast=force_blast)

    with collector.stage("molecules"):
        _process_molecules(df, collector, synonyms_dir, images_dir)

    with collector.stage("responses"):
        for validate in (
            validate_responsive_column,
            validate_parameter_column,
            validate_value_column,
            validate_value_nature_column,
            validate_unit_column,
            validate_concentration_column,
            validate_concentration_unit_column,
        ):
            collector.check(validate, df)

    with collector.stage("assay"):
        for validate in (
            validate_stimulation_flux_unit_column,
            validate_stimulation_duration_unit_column,
            validate_experimental_technique_column,
        ):
            collector.check(validate, df)

    with collector.stage("source"):
        _process_sources(df)

    with collector.stage("final_checks"):
        _check_row_count(df, input_row_count)
        _report_import_blockers(df)

    return PipelineRun(df=df, collector=collector, layout=layout)


def process_raw_excel_file(
    excel_path: Path,
    *,
    force_blast: bool = False,
    blast_mode: BlastMode = "ask",
    catch_unexpected: bool = False,
) -> PipelineRun:
    """Run the pipeline on one study and, if no error blocks it, export its
    CSV (and synonyms/stereo images) into the M2iOR website project."""
    study_id: str = excel_path.stem  # e.g. "001_Kreher_Neuron_2005"
    output_path = M2IOR_DATA_INPUT_PATH
    output_path.mkdir(parents=True, exist_ok=True)

    # Register if not already known; ignore the error if already registered
    try:
        _ = registry.register(study_id, excel_path)
    except ValueError:
        pass

    study_tracker: StudyFileTracker = registry.get(study_id)

    run = run_pipeline(
        excel_path,
        blast_mode=blast_mode,
        synonyms_dir=M2IOR_DATA_INPUT_PATH / "cid_synonyms",
        images_dir=M2IOR_MOLECULE_IMAGES_PATH,
        force_blast=force_blast,
        catch_unexpected=catch_unexpected,
    )

    if run.exit_code == 0 and run.df is not None:
        # Export enriched dataset as CSV with two-row (group / column) header
        run.exported_csv = output_path / f"{study_id}.csv"
        export_to_csv(run.df, run.exported_csv)

        study_tracker.complete()
        study_tracker.file_hash = compute_file_hash(excel_path)
        print(f"✓ {study_id}")
    else:
        error_message = "\n".join(run.collector.failure_messages)
        print(f"✗ {study_id} — FAILED: {error_message}")
        study_tracker.fail(error_message=error_message)

    registry.save(REGISTRY_FILE)
    return run


def check_excel_file(excel_path: Path, *, blast_mode: BlastMode) -> PipelineRun:
    """Validate one study without exporting anything: no CSV, no registry
    update, synonyms/stereo images rendered into a throwaway directory.
    External API results are still cached, as in a normal run."""
    with tempfile.TemporaryDirectory(prefix="m2ior_check_") as tmp:
        return run_pipeline(
            excel_path,
            blast_mode=blast_mode,
            synonyms_dir=Path(tmp) / "cid_synonyms",
            images_dir=Path(tmp) / "images",
            catch_unexpected=True,
        )


def _write_report(
    run: PipelineRun,
    excel_path: Path,
    *,
    mode: Literal["check", "export"],
    blast_mode: BlastMode,
    report_format: str,
    report_file: Path | None,
    stdout: Any,
) -> None:
    report = build_report(
        run.collector,
        mode=mode,
        study_id=excel_path.stem,
        excel_path=str(excel_path.resolve()),
        file_sha256=compute_file_hash(excel_path),
        layout=run.layout,
        blast_mode=blast_mode,
        exported_csv=str(run.exported_csv) if run.exported_csv else None,
    )
    if report_format == "json":
        text = json.dumps(report, indent=2, ensure_ascii=False)
    else:
        text = format_text_report(report)

    if report_file is None:
        print(text, file=stdout)
    else:
        report_file.parent.mkdir(parents=True, exist_ok=True)
        report_file.write_text(text + "\n", encoding="utf-8")
        print(f"{report['status'].upper()}: report written to {report_file}", file=sys.stderr)


def _print_internal_errors(run: PipelineRun) -> None:
    for issue in run.collector.issues:
        if issue.code == "internal_error":
            print(issue.details.get("traceback", issue.message), file=sys.stderr)


def _parse_index_selection(answer: str, count: int) -> list[int]:
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


def scan_and_select(
    registry: StudyProcessingRegistry, input_path: Path
) -> list[Path]:
    """List Excel files under input_path, newest-modified first, highlighting
    which ones are new or changed (by content hash) since their last
    successful run, then let the user pick which ones to process."""
    excel_paths = sorted(
        input_path.glob("**/*.xlsx"), key=lambda p: p.stat().st_mtime, reverse=True
    )
    if not excel_paths:
        print(f"No Excel files found under {input_path}")
        return []

    entries: list[tuple[Path, datetime, str]] = []
    for excel_path in excel_paths:
        study_id = excel_path.stem
        file_hash = compute_file_hash(excel_path)
        mtime = datetime.fromtimestamp(excel_path.stat().st_mtime)
        if study_id not in registry:
            status = "NEW"
        else:
            tracker = registry.get(study_id)
            status = "UNCHANGED" if tracker.file_hash == file_hash else "MODIFIED"
        entries.append((excel_path, mtime, status))

    to_update: list[int] = []
    for i, (excel_path, mtime, status) in enumerate(entries, start=1):
        line = f"[{i:3d}] {mtime:%Y-%m-%d %H:%M}  {excel_path.name}"
        if status == "NEW":
            logger.info("→ %s  (NEW)", line)
            to_update.append(i)
        elif status == "MODIFIED":
            logger.warning("→ %s  (MODIFIED)", line)
            to_update.append(i)
        else:
            print(f"    {line}")

    if not to_update:
        print("\nAll files are up to date.")
        return []

    default = ",".join(str(i) for i in to_update)
    print(f"\n{len(to_update)} file(s) new or modified.")
    answer = (
        input(
            f'Process which file(s)? ("1,3-5", "all", Enter for [{default}], '
            'or "none") '
        )
        .strip()
        .lower()
    )

    if answer in ("none", "n"):
        selected = []
    elif not answer:
        selected = to_update
    elif answer == "all":
        selected = list(range(1, len(entries) + 1))
    else:
        selected = _parse_index_selection(answer, len(entries))

    return [entries[i - 1][0] for i in selected]


_LOG_COLORS = {
    "DEBUG": "cyan",
    "INFO": "green",
    "WARNING": "yellow",
    "ERROR": "red",
    "CRITICAL": "bold_red",
}


class ColoredLevelFormatter(logging.Formatter):
    _default = colorlog.ColoredFormatter(
        fmt="%(log_color)s[%(levelname)s] %(message)s",
        log_colors=_LOG_COLORS,
    )
    _debug = colorlog.ColoredFormatter(
        fmt="%(log_color)s%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        log_colors=_LOG_COLORS,
    )

    def format(self, record: logging.LogRecord) -> str:
        formatter = self._debug if record.levelno == logging.DEBUG else self._default
        return formatter.format(record)


if __name__ == "__main__":
    _handler = colorlog.StreamHandler()
    _handler.setFormatter(ColoredLevelFormatter())
    logging.getLogger().addHandler(_handler)
    logging.getLogger().setLevel(logging.INFO)
    parser = argparse.ArgumentParser(description="M2iOR processing pipeline")
    _ = parser.add_argument(
        "--status", action="store_true", help="Show processing status"
    )
    _ = parser.add_argument("--file", help="Process a single Excel file")
    _ = parser.add_argument(
        "--check",
        metavar="FILE",
        help="Validate a single Excel file without exporting anything (no CSV, "
        "no registry update) and report every issue found. Exit code: 0 no "
        "error, 1 validation errors, 2 pipeline crash",
    )
    _ = parser.add_argument(
        "--report",
        choices=["text", "json"],
        help="Report format for --check (default: text) or --file (default: none)",
    )
    _ = parser.add_argument(
        "--report-file",
        metavar="PATH",
        help="Write the --report to PATH instead of stdout",
    )
    _ = parser.add_argument(
        "--blast",
        choices=["ask", "cached", "run"],
        help="Pending BLAST lookups: 'ask' prompts for which to run, 'cached' "
        "runs none (only cached results are applied), 'run' runs every new or "
        "retryable one without asking. Default: cached for --check, ask otherwise",
    )
    _ = parser.add_argument(
        "--scan",
        action="store_true",
        help="List Excel files (newest first), highlight new/modified ones "
        "by content hash, and choose which to process",
    )
    _ = parser.add_argument(
        "--pending", action="store_true", help="Process all pending Excel files"
    )
    _ = parser.add_argument(
        "--retry",
        action="store_true",
        help="Retry only failed files",
    )
    _ = parser.add_argument(
        "--all",
        action="store_true",
        help="Reprocess all files, ignoring DONE/FAILED status, like --force but without forcing BLAST re-queries.",
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
    do_scan = cast(bool, args.scan)
    do_pending = cast(bool, args.pending)
    do_retry = cast(bool, args.retry)
    do_all = cast(bool, args.all)
    force_arg = cast(str | None, args.force)
    check_arg = cast(str | None, args.check)
    report_format = cast(str | None, args.report)
    report_file_arg = cast(str | None, args.report_file)
    report_file = Path(report_file_arg) if report_file_arg else None
    blast_arg = cast(BlastMode | None, args.blast)

    input_path = Path("input/")

    single_file = check_arg or file_arg or (force_arg if force_arg != "__all__" else None)
    if (report_format or report_file) and not single_file:
        parser.error("--report/--report-file only apply to --check, --file or --force FILE")
    if single_file and not Path(single_file).is_file():
        parser.error(f"Excel file not found: {single_file}")

    if check_arg is not None:
        # stdout carries the report only: route every other print to stderr.
        real_stdout = sys.stdout
        check_blast: BlastMode = blast_arg or "cached"
        with contextlib.redirect_stdout(sys.stderr):
            check_run = check_excel_file(Path(check_arg), blast_mode=check_blast)
            _print_internal_errors(check_run)
            _write_report(
                check_run,
                Path(check_arg),
                mode="check",
                blast_mode=check_blast,
                report_format=report_format or "text",
                report_file=report_file,
                stdout=real_stdout,
            )
        sys.exit(check_run.exit_code)

    blast_mode: BlastMode = blast_arg or "ask"

    # With a report requested, stdout carries the report only.
    status_stream = sys.stderr if report_format else sys.stdout
    print("Loading processing registry...", file=status_stream)
    registry = StudyProcessingRegistry.load(REGISTRY_FILE)
    print("Done", file=status_stream)

    if show_status:
        print(registry.summary)
    elif do_scan:
        selected_paths = scan_and_select(registry, input_path)
        if selected_paths:
            print(f"\nProcessing {len(selected_paths)} file(s)...")
            for excel_file_path in selected_paths:
                print(f"File: {excel_file_path}")
                process_raw_excel_file(excel_file_path, blast_mode=blast_mode)
    elif file_arg is not None or (force_arg is not None and force_arg != "__all__"):
        target_path = Path(cast(str, file_arg or force_arg))
        force_blast = file_arg is None
        if force_blast:
            study_id = target_path.stem
            if study_id in registry:
                tracker = registry.get(study_id)
                tracker.status = ProcessingStatus.PENDING
                tracker.error_message = None
        real_stdout = sys.stdout
        # With a report requested, stdout carries the report only.
        redirect = (
            contextlib.redirect_stdout(sys.stderr)
            if report_format
            else contextlib.nullcontext()
        )
        with redirect:
            file_run = process_raw_excel_file(
                target_path,
                force_blast=force_blast,
                blast_mode=blast_mode,
                catch_unexpected=True,
            )
            _print_internal_errors(file_run)
            if report_format:
                _write_report(
                    file_run,
                    target_path,
                    mode="export",
                    blast_mode=blast_mode,
                    report_format=report_format,
                    report_file=report_file,
                    stdout=real_stdout,
                )
        sys.exit(file_run.exit_code)
    elif do_retry:
        failed_jobs: list[StudyFileTracker] = registry.all_failed()
        for job in failed_jobs:
            job.status = ProcessingStatus.PENDING
            job.error_message = None
            process_raw_excel_file(job.file_path, blast_mode=blast_mode)
    elif do_pending:
        found_files: list[StudyFileTracker] = registry.all_pending()
        print(f"Found {len(found_files)} files")
        _ = input("Press Enter to continue")
        for excel_file in found_files:
            print(f"File: {excel_file.file_path}")
            process_raw_excel_file(excel_file.file_path, blast_mode=blast_mode)
    elif do_all:
        # Reset all existing trackers to PENDING before reprocessing, without
        # forcing BLAST re-queries (unlike --force).
        for tracker in registry.all_done() + registry.all_failed():
            tracker.status = ProcessingStatus.PENDING
            tracker.error_message = None
        found_file_paths = sorted(input_path.glob("**/*.xlsx"))
        print(f"Found {len(found_file_paths)} files — reprocessing all")
        _ = input("Press Enter to continue")
        for excel_file_path in found_file_paths:
            print(f"File: {excel_file_path}")
            process_raw_excel_file(excel_file_path, blast_mode=blast_mode)
    elif force_arg is not None:
        # Reset all existing trackers to PENDING before reprocessing
        for tracker in registry.all_done() + registry.all_failed():
            tracker.status = ProcessingStatus.PENDING
            tracker.error_message = None
        found_file_paths = sorted(input_path.glob("**/*.xlsx"))
        print(f"Found {len(found_file_paths)} files — forcing full reprocess")
        _ = input("Press Enter to continue")
        for excel_file_path in found_file_paths:
            print(f"File: {excel_file_path}")
            process_raw_excel_file(
                excel_file_path, force_blast=True, blast_mode=blast_mode
            )
    else:
        parser.print_help()
