# pipeline/scripts/export_csv.py

import csv
from pathlib import Path

import pandas as pd

from scripts.columns import ASSAY, CO_RECEPTOR, MOLECULE, RECEPTOR, REFERENCE, RESPONSE

# Maps MultiIndex group labels (from the Excel template) to output CSV group names.
_GROUP_RENAME: dict[str, str] = {
    RECEPTOR: "receptors",
    CO_RECEPTOR: "co_receptors",
    "Co Receptor": "co_receptors",
    MOLECULE: "compounds",
    "Compound": "compounds",
    "Compounds": "compounds",
    RESPONSE: "responses",
    "Responses": "responses",
    ASSAY: "assays",
    "Assays": "assays",
    "Assay info": "assays",
    REFERENCE: "references",
    "References": "references",
    "Source": "references",
}

# Renames individual column names for the output CSV.
_COL_RENAME: dict[str, str] = {
    "Accession": "ID",
    "sequence_ref": "Sequence_ref",
    "solvent used for dilution": "Solvent used for dilution",
}
# Case-insensitive fallback index (lower-cased key → canonical output name).
_COL_RENAME_CI: dict[str, str] = {k.lower(): v for k, v in _COL_RENAME.items()}


def export_to_csv(df: pd.DataFrame, path: Path) -> None:
    """
    Export the enriched MultiIndex DataFrame to a CSV with two header rows:
      Row 1 — group labels  (receptors, co_receptors, compounds, ...)
      Row 2 — column names  (Order, Species, ..., Order.1, Species.1, ...)
      Row 3+ — data

    Duplicate column names across groups receive a .1, .2, … suffix in row 2,
    matching standard pandas flat-column deduplication.
    Group labels with no match in _GROUP_RENAME are kept as-is;
    columns whose group is NaN (unlabelled) are written with an empty group label.
    """
    # Build the two header rows
    groups_row: list[str] = []
    seen_cols: dict[str, int] = {}
    col_row: list[str] = []

    for group, col in df.columns:
        g_str = str(group)
        groups_row.append("" if g_str == "nan" else _GROUP_RENAME.get(g_str, g_str))

        raw = str(col)
        c_str = _COL_RENAME.get(raw) or _COL_RENAME_CI.get(raw.lower(), raw)
        if c_str in seen_cols:
            seen_cols[c_str] += 1
            col_row.append(f"{c_str}.{seen_cols[c_str]}")
        else:
            seen_cols[c_str] = 0
            col_row.append(c_str)

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(groups_row)
        writer.writerow(col_row)
        for _, row in df.iterrows():
            writer.writerow(row.tolist())
