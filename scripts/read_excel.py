# pipeline/scripts/read_excel.py

from pathlib import Path

import pandas as pd

from scripts.columns import COLUMNS_BY_GROUP, GROUPS_ORDER


def get_raw_data_from_excel_file(excel_path: Path) -> pd.DataFrame:
    """
    Read the 2nd sheet of the Excel file and return the raw data as a DataFrame
    with a MultiIndex on the columns.

    Expected structure of the M2iOR template:
      - Row 0: group labels ("Receptor", "Co-Receptor", "Molecule", …)
               with NaN for merged cells within each group
      - Row 1: column names ("Order", "Species", "Gene Name", …)
      - Row 2+: data

    Returns a DataFrame with MultiIndex columns: (group, column_name).

    Raises:
        ValueError: if the file does not contain a 2nd sheet, or if groups/columns
                    do not match the expected structure.
    """
    xl = pd.ExcelFile(excel_path)

    if len(xl.sheet_names) < 2:
        raise ValueError(
            "Le fichier Excel doit contenir au moins 2 feuilles "
            + f"(feuilles trouvées : {xl.sheet_names}) : {excel_path}"
        )

    raw_data_sheet_name = xl.sheet_names[1]

    # Lire toutes la première ligne sans inférence d'en-tête
    raw_check = pd.read_excel(  # pyright: ignore[reportUnknownMemberType]
        xl, sheet_name=raw_data_sheet_name, header=None, nrows=2
    )

    # Some source files carry a "used range" far wider than the real data
    # (e.g. a fill/border applied across an entire row up to column XFD) -
    # left unbounded, openpyxl/pandas would scan millions of phantom empty
    # cells. Bound the real read to the last column actually holding a
    # group or a column label.
    last_used_col = max(
        raw_check.iloc[0].last_valid_index(),
        raw_check.iloc[1].last_valid_index(),
    )
    raw_check = raw_check.iloc[:, : last_used_col + 1]

    # Ligne 0 : labels de groupe — forward-fill pour les cellules fusionnées
    groups = raw_check.iloc[0].ffill()
    # Ligne 1 : noms de colonnes
    columns = _normalize_columns(groups, raw_check.iloc[1])

    _validate_structure(groups, columns, excel_path)

    # Lire toutes les lignes sans inférence d'en-tête
    raw = pd.read_excel(  # pyright: ignore[reportUnknownMemberType]
        xl, sheet_name=raw_data_sheet_name, header=None, usecols=range(last_used_col + 1)
    )

    # Construction du MultiIndex (groupe, colonne)
    multi_columns = pd.MultiIndex.from_arrays([groups, columns])

    # Données à partir de la ligne 2
    df = raw.iloc[2:].reset_index(drop=True)
    df.columns = multi_columns

    return df


def _normalize_columns(groups: pd.Series, columns: pd.Series) -> pd.Series:
    """Rename headers to their canonical casing (e.g. "value nature" -> "Value Nature")."""
    canonical_by_group = {
        group: {col.lower(): col for col in cols} for group, cols in COLUMNS_BY_GROUP
    }
    normalized = columns.copy()
    for idx, col in columns.items():
        if pd.isna(col):
            continue
        canonical = canonical_by_group.get(groups[idx], {}).get(str(col).lower())
        if canonical is not None:
            normalized[idx] = canonical
    return normalized


def _validate_structure(
    groups: pd.Series, columns: pd.Series, excel_path: Path
) -> None:
    actual_groups_order = list(dict.fromkeys(groups.dropna()))
    if actual_groups_order != GROUPS_ORDER:
        msg = (
            f"{excel_path}: groupes inattendus.\n"
            f"  attendu : {GROUPS_ORDER}\n"
            f"  trouvé  : {actual_groups_order}"
        )
        raise ValueError(msg)

    for group, expected_cols in COLUMNS_BY_GROUP:
        actual_cols = list(columns[groups == group])
        missing = [c for c in expected_cols if c not in actual_cols]
        if missing:
            msg = (
                f"{excel_path}: missing columns in '{group}'.\n"
                f"  missing : {missing}\n"
                f"  found   : {actual_cols}"
            )
            raise ValueError(msg)
