# pipeline/scripts/process_receptors.py
from typing import cast

import pandas as pd

from scripts.cache_manager import get_cache
from scripts.columns import (
    DATABASE,
    IDENTITY,
    MUTATION,
    SEQUENCE,
    SEQUENCE_REF,
    SPECIES,
    UNIPROT_ID,
)
from scripts.fetch_blast import fetch_blast_reference
from scripts.find_protein_mutations import align_and_annotate


def get_unique_uniprot_ids(
    df: pd.DataFrame, groups: list[str], column_name: str
) -> list[str]:
    """Collect all non-null unique UniProt IDs across the given column groups."""
    ids: set[str] = set()
    for group in groups:
        if group not in df.columns.get_level_values(0):
            raise ValueError(f"Group {group} not in the DataFrame.")
        col = df[group][column_name]
        ids.update(col.dropna().astype(str).str.strip().unique())
    ids.discard("")
    return sorted(ids)


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


def collect_blast_queries(df: pd.DataFrame, groups: list[str]) -> list[tuple[str, str]]:
    """
    Return unique (sequence, species) pairs for rows that have no UniProt ID
    but have both a Sequence and a Species value.  Fully vectorized — no per-row
    iteration.
    """
    queries: set[tuple[str, str]] = set()
    for group in groups:
        if group not in df.columns.get_level_values(0):
            continue
        sub = df[group]
        missing = sub[UNIPROT_ID].isna() | (
            sub[UNIPROT_ID].astype(str).str.strip() == ""
        )
        candidates = sub.loc[missing, [SEQUENCE, SPECIES]].dropna()
        seq = candidates[SEQUENCE].astype(str).str.strip()
        species = candidates[SPECIES].astype(str).str.strip()
        valid = (seq != "") & (species != "")
        queries.update(zip(seq[valid], species[valid], strict=False))
    return list(queries)


def resolve_missing_ids_via_blast(
    df: pd.DataFrame,
    groups: list[str],
    queries: list[tuple[str, str]],
) -> dict[str, str]:
    """
    Run BLAST once per unique (sequence, species) pair, then fill the UNIPROT_ID
    cells with the resulting GenBank accessions.  The DataFrame fill is vectorized
    — no per-row iteration.

    Results are persisted in cache/blast_cache.json.  Returns a transient
    {accession: sequence_ref} dict for use in enrich_with_reference_and_mutations.
    """
    # One API call per unique pair (caching handled inside fetch_blast_reference)
    accession_map: dict[tuple[str, str], str] = {}
    blast_refs: dict[str, str] = {}
    for i, (seq, species) in enumerate(queries, start=1):
        print(f"  [BLAST] [{i}/{len(queries)}] species={species!r} ...", flush=True)  # noqa: T201
        result = fetch_blast_reference(seq, species)
        if result is not None:
            accession, seq_ref = result
            accession_map[(seq, species)] = accession
            blast_refs[accession] = seq_ref

    # Fill UNIPROT_ID cells — vectorized per group
    for group in groups:
        if group not in df.columns.get_level_values(0):
            continue
        sub = df[group]
        missing = sub[UNIPROT_ID].isna() | (
            sub[UNIPROT_ID].astype(str).str.strip() == ""
        )
        seq_col = sub.loc[missing, SEQUENCE].astype(str).str.strip()
        species_col = sub.loc[missing, SPECIES].astype(str).str.strip()
        pairs = zip(seq_col, species_col, strict=False)
        new_ids = pd.Series(
            [accession_map.get((s, sp)) for s, sp in pairs],
            index=seq_col.index,
            dtype=object,
        )
        filled = new_ids.dropna()
        if not filled.empty:
            df.loc[filled.index, (group, UNIPROT_ID)] = filled

    return blast_refs


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

    for group in groups:
        # One disk read per unique UniProt ID (GenBank IDs come from blast_refs).
        uniprot_cache: dict[str, str | None] = {
            uid: (data["sequence"]["value"] if (data := get_cache(uid)) else None)
            for uid in all_unique_uniprot_ids
            if uid not in blast_refs
        }

        seq_refs: list[str | None] = []
        identities: list[float | None] = []
        mutations_list: list[str | None] = []
        database_list: list[str | None] = []

        for _, row in df[group].iterrows():
            uid = row[UNIPROT_ID]
            uid_str = str(uid).strip() if pd.notna(uid) else None

            if uid_str and uid_str in blast_refs:
                seq_ref, database = blast_refs[uid_str], "genbank"
            elif uid_str:
                seq_ref, database = uniprot_cache.get(uid_str), "uniprot"
            else:
                seq_ref, database = None, "undefined"

            seq_refs.append(seq_ref)
            database_list.append(database)

            seq = row[SEQUENCE]
            if seq_ref is None or pd.isna(seq) or not str(seq).strip():
                identities.append(None)
                mutations_list.append(None)
            else:
                mut_str, identity = align_and_annotate(str(seq).strip(), seq_ref)
                identities.append(identity)
                mutations_list.append(mut_str)

        df.loc[:, (group, SEQUENCE_REF)] = seq_refs
        df.loc[:, (group, IDENTITY)] = identities
        df.loc[:, (group, MUTATION)] = mutations_list
        df.loc[:, (group, DATABASE)] = database_list
