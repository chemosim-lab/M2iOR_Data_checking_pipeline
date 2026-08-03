# pipeline/scripts/process_receptors.py
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from re import Match
from typing import Any, cast

import colorlog
import pandas as pd

from scripts.cache_manager import get_cache
from scripts.columns import (
    ACCESSION,
    DATABASE,
    IDENTITY,
    IDENTITY_SHORT,
    MUTATION,
    RECEPTOR_NAME,
    SEQUENCE,
    SEQUENCE_REF,
    SPECIES,
    UNIPROT_ID,
)
from scripts.fetch_blast import fetch_blast_reference
from scripts.find_protein_mutations import align_and_annotate

logger = colorlog.getLogger(__name__)

# Some source Excel files spell out "Not available" in the Accession column
# instead of leaving it blank - treated as a recognized sentinel rather than
# a real ID: propagated as-is into Sequence_ref instead of "".
NOT_AVAILABLE = "Not available"

# Below this % identity to its reference sequence, a receptor is flagged for
# manual review - could be a mislabeled species, a wrong reference match, or
# a genuinely divergent sequence worth double-checking.
IDENTITY_WARNING_THRESHOLD = 95.0

_OR_NAME_RE = re.compile(r"^Or\d+(?:-\d+)?[a-z]?$")
_OR_NAME_RE_LOOSE = re.compile(r"^Or(\d+)(-\d+)?([a-zA-Z]?)$", re.IGNORECASE)
_ORCO_RE_LOOSE = re.compile(r"^Orco$", re.IGNORECASE)
_OR_FROM_DESC_RE = re.compile(
    r"olfactory\s+receptor\s+(?:Or)?(\d+)(-\d+)?([a-zA-Z]?)", re.IGNORECASE
)


def get_unique_accessions(
    df: pd.DataFrame, groups: list[str], column_name: str, *, verbose: bool = True
) -> list[str]:
    """Collect all non-null unique accessions (UniProt & NCBI) across the given
    column groups."""
    accessions_set: set[str] = set()
    for group in groups:
        if group not in df.columns.get_level_values(0):
            msg = f"Group {group} not in the DataFrame."
            raise ValueError(msg)
        col = df[group][column_name]
        group_accessions = set(col.dropna().astype(str).str.strip().unique())
        group_accessions.discard("")
        if verbose:
            logger.info(
                "  %s: %d unique accession(s).", group, len(group_accessions)
            )
        accessions_set.update(group_accessions)
    accessions_set.discard("")
    if verbose:
        logger.info("Found %d unique accession(s) in total.", len(accessions_set))
    return sorted(accessions_set)


def _get_normalized_receptor_name(names: list[str]) -> str | None:
    m: Match[str] | None = next(
        (m for name in names if (m := _OR_NAME_RE_LOOSE.match(name))),
        None,
    )
    if m:
        return "Or" + m.group(1) + (m.group(2) or "") + m.group(3).lower()
    if any(_ORCO_RE_LOOSE.match(name) for name in names):
        return "Orco"
    return None


def _extract_receptor_name_data(uid: str, data: dict[str, Any]) -> dict[str, Any]:
    output_data = dict[str, Any]()
    try:
        output_data["geneName"] = data["genes"][0]["geneName"]["value"]
    except KeyError as e:
        msg = f"UID:[{uid}] Unexpected UniProt cache structure: missing key {e}"
        raise ValueError(msg) from e

    try:
        synonyms_data: list[dict[str, str]] = data["synonyms"]
        output_data["synonyms"] = [syn["value"] for syn in synonyms_data]
    except KeyError:
        output_data["synonyms"] = None

    return output_data


def _extract_or_name_from_description(desc: str) -> str | None:
    """Extract a normalized Or<N>[-N][a-z] name from a protein description string."""
    m = _OR_FROM_DESC_RE.search(desc)
    if m:
        return "Or" + m.group(1) + (m.group(2) or "") + m.group(3).lower()
    return None


def _get_receptor_name(uid: str) -> str | None:

    data: dict[str, Any] | None = get_cache(uid, subdir="receptors")
    if not data:
        msg = f"UID: in _get_receptor_name [{uid}] Cache not found for this UID."
        raise ValueError(msg)
    receptor_name_data: dict[str, Any] = _extract_receptor_name_data(uid, data)

    receptor_name = receptor_name_data.get("geneName", "")
    synonyms: list[str] = receptor_name_data.get("synonyms") or []

    normalized_receptor_name: str | None = _get_normalized_receptor_name(
        [receptor_name, *synonyms]
    )
    if normalized_receptor_name:
        return normalized_receptor_name

    # NCBI fallback: extract Or<N> from protein description, or use accession
    if data.get("source") == "ncbi":
        protein_desc: str = data.get("protein_description") or ""
        if protein_desc:
            from_desc = _extract_or_name_from_description(protein_desc)
            if from_desc:
                return from_desc
        return uid  # last resort: raw accession number

    if receptor_name:
        return receptor_name

    msg = (
        f"UID:[{uid}] No Or<N>[a-z] receptor name found - "
        f"geneName={receptor_name!r}, synonyms={receptor_name_data.get('synonyms')}"
    )
    raise ValueError(msg)


def process_receptors_name_columns(
    df: pd.DataFrame, groups: list[str], all_unique_uniprot_ids: list[str]
) -> None:

    found_receptor_names_by_uid: dict[str, list[str]] = {}
    for group in groups:
        # Checking groups
        if group not in df.columns.get_level_values(0):
            msg = f"{group} missing in df"
            raise ValueError(msg)

        # Extract Uniprot ID and Receptor Name columns
        uid_recname_columns = df[group][[ACCESSION, RECEPTOR_NAME]]
        uid_recname_columns[ACCESSION] = (
            uid_recname_columns[ACCESSION].astype(str).str.strip()
        )

        # Store
        for uid, group_df in uid_recname_columns.groupby(ACCESSION):
            names: list[str] = group_df[RECEPTOR_NAME].unique()
            found_receptor_names_by_uid[str(uid)] = [
                str(name).strip() for name in names if pd.notna(name)
            ]

    receptor_name_by_uid: dict[str, str | None] = {}
    for uid in all_unique_uniprot_ids:
        found_receptor_names: list[str] = found_receptor_names_by_uid.get(uid, [""])
        normalized = _get_normalized_receptor_name(found_receptor_names)
        if normalized:
            receptor_name_by_uid[uid] = normalized
        else:
            receptor_name_by_uid[uid] = _get_receptor_name(uid)

    for group in groups:
        receptors_names = df[group][ACCESSION].map(receptor_name_by_uid)
        df.loc[:, (group, RECEPTOR_NAME)] = receptors_names.astype(
            df[group][RECEPTOR_NAME].dtype
        )

    enrich_species_column(df, groups, all_unique_uniprot_ids)


def enrich_species_column(
    df: pd.DataFrame, groups: list[str], all_unique_uniprot_ids: list[str]
) -> None:
    """Replace Species values with the organism scientific name from UniProt cache."""
    species_by_uid: dict[str, str] = {}
    for uid in all_unique_uniprot_ids:
        data = get_cache(uid, subdir="receptors")
        if not data:
            continue
        try:
            name: str = data["organism"]["scientificName"]
            if name:
                species_by_uid[uid] = name
        except (KeyError, TypeError):
            continue

    # Warn once per unique (accession, reported species) pair, not once per
    # row - the same mistyped/wrong accession is typically repeated across
    # dozens of rows.
    warned: set[tuple[str, str]] = set()

    for group in groups:
        uid_col = df[group][ACCESSION].astype(str).str.strip()
        updated = uid_col.map(species_by_uid)
        mask = updated.notna()
        if not mask.any():
            continue

        # Flag (but don't block on) a mismatch between the author-reported species
        # and the one UniProt associates with this accession - often a sign of a
        # mistyped or wrong accession number.
        reported_raw = df[group][SPECIES]
        reported = reported_raw.astype(str).str.strip()
        has_reported = reported_raw.notna() & (reported != "")
        mismatch = mask & has_reported & (reported.str.lower() != updated.str.lower())
        for idx in df.index[mismatch]:
            key = (uid_col[idx], reported[idx])
            if key in warned:
                continue
            warned.add(key)
            logger.warning(
                "  %s row %s: accession %s maps to UniProt species %r, "
                "which differs from the reported species %r - overwriting.",
                group,
                idx,
                uid_col[idx],
                updated[idx],
                reported[idx],
            )

        df.loc[updated.index[mask], (group, SPECIES)] = updated[mask].to_numpy()


def add_empty_column_after(
    df: pd.DataFrame,
    groups: list[str],
    after_column: str,
    new_column: str,
) -> None:
    for group in groups:
        if group not in df.columns.get_level_values(0):
            msg = f"{group} missing in df"
            raise ValueError(msg)

        if after_column not in df[group].columns:
            msg_0 = f'"{after_column}" column not in group "{group}"'
            raise ValueError(msg_0)

        cols = cast("list[tuple[str, str]]", list[str](df.columns))
        insert_at = cols.index((group, after_column)) + 1

        df.insert(loc=insert_at, column=(group, new_column), value=None)


def collect_blast_queries(
    df: pd.DataFrame,
    groups: list[str],
    fallback_accessions: list[str] | None = None,
) -> list[tuple[str, str]]:
    """
    Return unique (sequence, species) pairs for rows that have no UniProt or NCBI
    accession number but have both a Sequence and a Species value.
    """
    fallback_set: set[str] = set(fallback_accessions) if fallback_accessions else set()
    queries: set[tuple[str, str]] = set()
    for group in groups:
        if group not in df.columns.get_level_values(0):
            continue
        sub = df[group]
        uid_col = sub[ACCESSION].astype(str).str.strip()
        missing = sub[ACCESSION].isna() | (uid_col == "") | uid_col.isin(fallback_set)
        missing_rows = sub.loc[missing, [SEQUENCE, SPECIES]]
        seq = missing_rows[SEQUENCE].astype(str).str.replace(r"\s+", "", regex=True)
        species = missing_rows[SPECIES].astype(str).str.strip()
        has_seq = missing_rows[SEQUENCE].notna() & (seq != "")
        has_species = missing_rows[SPECIES].notna() & (species != "")

        # Sequences with no accession AND no Species can't be BLASTed either —
        # surface them so they aren't silently dropped.
        unqueryable_seqs = set(seq[has_seq & ~has_species])
        if unqueryable_seqs:
            logger.warning(
                "  %s: %d sequence(s) have no accession number in the Excel "
                "input and no Species, cannot run BLAST.",
                group,
                len(unqueryable_seqs),
            )

        valid = has_seq & has_species
        group_queries = set(zip(seq[valid], species[valid], strict=False))
        if group_queries:
            logger.warning(
                "  %s: %d sequence(s) have no accession number in the Excel input.",
                group,
                len(group_queries),
            )
        queries.update(group_queries)
    if queries:
        logger.warning(
            "%s sequence(s) have no accession number in the Excel input, in total.",
            len(queries),
        )
    return list(queries)


def resolve_missing_accessions_via_blast(
    df: pd.DataFrame,
    groups: list[str],
    queries: list[tuple[str, str]],
) -> dict[str, str]:
    """
    Run BLAST once per unique (sequence, species) pair, then fill the ACCESSION
    cells with the resulting GenBank accessions.  The DataFrame fill is vectorized
    — no per-row iteration.

    Results are persisted in cache/blast_cache.json.  Returns a transient
    {accession: sequence_ref} dict for use in enrich_with_reference_and_mutations.
    """
    # Submit up to 3 BLAST jobs concurrently — each spends most of its time polling,
    # so threads are efficient. Cache writes are protected by a lock in fetch_blast.
    accession_map: dict[tuple[str, str], str] = {}
    blast_refs: dict[str, str] = {}
    with ThreadPoolExecutor(max_workers=3) as executor:
        future_to_query = {
            executor.submit(fetch_blast_reference, seq, species): (seq, species)
            for seq, species in queries
        }
        for future in as_completed(future_to_query):
            seq, species = future_to_query[future]
            result = future.result()
            if result is not None:
                accession, seq_ref = result
                accession_map[(seq, species)] = accession
                blast_refs[accession] = seq_ref

    # Fill ACCESSION cells — vectorized per group
    for group in groups:
        if group not in df.columns.get_level_values(0):
            continue
        sub = df[group]
        missing = sub[ACCESSION].isna() | (sub[ACCESSION].astype(str).str.strip() == "")
        # Must match collect_blast_queries' normalization exactly (full
        # whitespace removal, not just trimming ends) - some source cells have
        # stray internal whitespace, and a mismatched key here silently drops
        # the fill instead of erroring.
        seq_col = sub.loc[missing, SEQUENCE].astype(str).str.replace(
            r"\s+", "", regex=True
        )
        species_col = sub.loc[missing, SPECIES].astype(str).str.strip()
        pairs = zip(seq_col, species_col, strict=False)
        new_ids = pd.Series(
            [accession_map.get((s, sp)) for s, sp in pairs],
            index=seq_col.index,
            dtype=object,
        )
        filled = new_ids.dropna()
        if not filled.empty:
            df.loc[filled.index, (group, ACCESSION)] = filled

    return blast_refs


def _resolve_seq_ref(
    uid: Any,
    blast_refs: dict[str, str],
    ncbi_uid_set: set[str],
    uniprot_cache: dict[str, str | None],
) -> tuple[str | None, str]:
    uid_str = str(uid).strip() if pd.notna(uid) else None
    if uid_str and uid_str.lower() == NOT_AVAILABLE.lower():
        return NOT_AVAILABLE, "undefined"
    if uid_str and uid_str in blast_refs:
        seq_ref, database = blast_refs[uid_str], "genbank"
    elif uid_str and uid_str in ncbi_uid_set:
        seq_ref, database = uniprot_cache.get(uid_str), "genbank"
    elif uid_str:
        seq_ref, database = uniprot_cache.get(uid_str), "uniprot"
    else:
        seq_ref, database = None, "undefined"
    if seq_ref is not None:
        seq_ref = "".join(seq_ref.split())
    return seq_ref, database


def _normalize_sequence(seq: Any, seq_ref: str | None) -> Any:
    seq_empty = pd.isna(seq) or not str(seq).strip()
    if seq_ref and seq_ref != NOT_AVAILABLE and seq_empty:
        return seq_ref
    if not seq_empty:
        return "".join(str(seq).split())
    return seq


def _compute_mutations(
    seq: Any,
    seq_ref: str | None,
) -> tuple[str | None, float | None, float | None]:
    # "Not available" is a placeholder, not a real reference sequence -
    # aligning against it would produce meaningless mutation output.
    if (
        seq_ref is None
        or seq_ref == NOT_AVAILABLE
        or pd.isna(seq)
        or not str(seq).strip()
    ):
        return None, None, None
    mut_str, pid_aln, pid_short = align_and_annotate(seq, seq_ref)
    return mut_str, pid_aln, pid_short


def enrich_with_reference_and_mutations(
    df: pd.DataFrame,
    groups: list[str],
    all_unique_uniprot_ids: list[str],
    blast_refs: dict[str, str] | None = None,
) -> None:
    """
    For each protein group, fetch reference sequences then compute sequence identity
    and mutations against each row's Sequence.

    UniProt sequences come from the per-accession disk cache.
    GenBank sequences come from blast_refs {accession: sequence_ref}, populated by
    resolve_missing_ids_via_blast from blast_cache.json.
    """
    if blast_refs is None:
        blast_refs = {}

    ncbi_uid_set: set[str] = {
        uid
        for uid in all_unique_uniprot_ids
        if (data := get_cache(uid, subdir="receptors")) and data.get("source") == "ncbi"
    }

    # Many rows across a group share the same (accession, sequence) - e.g. one
    # receptor tested against dozens of odorants - so memoize the expensive
    # pairwise alignment per unique (seq, seq_ref) pair instead of recomputing
    # it for every row.
    mutation_cache: dict[tuple[Any, str | None], tuple[str | None, float | None, float | None]] = {}

    for group in groups:
        # One disk read per unique UniProt ID (GenBank IDs come from blast_refs).
        uniprot_cache: dict[str, str | None] = {
            uid: (
                data["sequence"]["value"]
                if (data := get_cache(uid, subdir="receptors"))
                else None
            )
            for uid in all_unique_uniprot_ids
            if uid not in blast_refs
        }

        sequences: list[Any] = []
        seq_refs: list[str | None] = []
        identities: list[float | None] = []
        identities_short: list[float | None] = []
        mutations_list: list[str | None] = []
        database_list: list[str | None] = []

        for _idx, row in df[group].iterrows():
            seq_ref, database = _resolve_seq_ref(
                row[ACCESSION], blast_refs, ncbi_uid_set, uniprot_cache
            )
            seq_refs.append(seq_ref)
            database_list.append(database)

            seq = _normalize_sequence(row[SEQUENCE], seq_ref)
            sequences.append(seq)
            cache_key = (seq, seq_ref)
            is_new_pair = cache_key not in mutation_cache
            if is_new_pair:
                mutation_cache[cache_key] = _compute_mutations(seq, seq_ref)
            mut_str, pid_aln, pid_short = mutation_cache[cache_key]
            identities.append(pid_aln)
            identities_short.append(pid_short)
            mutations_list.append(mut_str)

            # Warn once per unique (sequence, reference) pair, not once per row -
            # the same receptor is typically tested against many odorants.
            if (
                is_new_pair
                and pid_aln is not None
                and pid_aln < IDENTITY_WARNING_THRESHOLD
            ):
                logger.warning(
                    "  %s: %s (%s, accession=%s) has %.2f%% identity to its "
                    "reference sequence, below the %g%% threshold.",
                    group,
                    row[RECEPTOR_NAME],
                    row[SPECIES],
                    row[ACCESSION],
                    pid_aln,
                    IDENTITY_WARNING_THRESHOLD,
                )

        df.loc[:, (group, SEQUENCE)] = sequences
        df.loc[:, (group, SEQUENCE_REF)] = seq_refs
        df.loc[:, (group, IDENTITY)] = identities
        df.loc[:, (group, IDENTITY_SHORT)] = identities_short
        df.loc[:, (group, MUTATION)] = mutations_list
        df.loc[:, (group, DATABASE)] = database_list
