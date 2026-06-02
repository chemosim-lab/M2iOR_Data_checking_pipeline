# pipeline/scripts/extract_receptors.py

import pandas as pd


# Mapping: Excel column name (under "Receptor") → output CSV column name
_RECEPTOR_COLS: dict[str, str] = {
    "Order": "order",
    "Species": "species",
    "Gene Name": "gene",
    "UniProt ID": "uniprot_id",
    "Sequence": "sequence",
}


def extract_receptors(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract unique receptors from the raw multi-index DataFrame produced by
    ``get_raw_data_from_excel_file``.

    Uniqueness key: (UniProt ID, Sequence).
    The same UniProt accession can appear with different sequences when the
    study includes *variants* (e.g. mutants, chimeras), so both fields are
    required to identify a receptor unambiguously.

    Raises:
        ValueError: if any receptor row is missing a UniProt ID or a Sequence.

    Returns:
        DataFrame with columns: uniprot_id, gene, order, species, sequence
        (one row per unique receptor / variant, sorted by uniprot_id then gene).
    """
    rec: pd.DataFrame = df["Receptor"][list(_RECEPTOR_COLS.keys())].copy()

    # Normalise: treat blank / whitespace-only strings as NaN
    rec = rec.replace(r"^\s*$", pd.NA, regex=True)

    # --- Validation ---

    missing_uniprot = rec["UniProt ID"].isna()
    if missing_uniprot.any():
        n = int(missing_uniprot.sum())
        # Report the 0-based data-row indices (row 0 = 3rd row of the Excel sheet)
        bad_rows = missing_uniprot[missing_uniprot].index.tolist()
        raise ValueError(
            f"{n} receptor row(s) are missing a UniProt ID "
            f"(data row index: {bad_rows})"
        )

    missing_seq = rec["Sequence"].isna()
    if missing_seq.any():
        n = int(missing_seq.sum())
        bad_rows = missing_seq[missing_seq].index.tolist()
        raise ValueError(
            f"{n} receptor row(s) are missing a Sequence "
            f"(data row index: {bad_rows})"
        )

    # --- Deduplication ---
    # Same (UniProt ID, Sequence) → same receptor/variant
    unique_receptors: pd.DataFrame = (
        rec.drop_duplicates(subset=["UniProt ID", "Sequence"])
        .rename(columns=_RECEPTOR_COLS)
        .sort_values(["uniprot_id", "gene"])
        .reset_index(drop=True)
    )

    return unique_receptors
