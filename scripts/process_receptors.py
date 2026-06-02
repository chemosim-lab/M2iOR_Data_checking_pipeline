# pipeline/scripts/enrich_sequences.py
from typing import Any, Callable, cast

import pandas as pd

from scripts.cache_manager import get_cache
from scripts.columns import DATABASE, IDENTITY, MUTATION, SEQUENCE_REF, UNIPROT_ID
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


def enrich_with_reference_and_mutations(df: pd.DataFrame, groups: list[str]) -> None:
    for group in groups:
        uniprot_ids = get_unique_uniprot_ids(df, groups=[group], column_name=UNIPROT_ID)

        # One disk read per unique ID, not per row — extract sequence immediately
        tmp_cache: dict[str, str | None] = {
            uid: (data["sequence"]["value"] if (data := get_cache(uid)) else None)
            for uid in uniprot_ids
        }

        seq_refs: list[str | None] = []
        identities: list[float | None] = []
        mutations_list: list[str | None] = []
        database_list: list[str | None] = []

        for _, row in df[group].iterrows():
            uid = row[UNIPROT_ID]
            seq_ref = tmp_cache.get(str(uid).strip()) if pd.notna(uid) else None
            seq_refs.append(seq_ref)
            database_list.append("uniprot")

            seq = row["Sequence"]
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
