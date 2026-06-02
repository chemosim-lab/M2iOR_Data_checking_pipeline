# pipeline/scripts/read_excel.py

import pandas as pd
from pathlib import Path


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
        ValueError: if the file does not contain a 2nd sheet.
    """
    xl = pd.ExcelFile(excel_path)

    if len(xl.sheet_names) < 2:
        raise ValueError(
            "Le fichier Excel doit contenir au moins 2 feuilles "
            + f"(feuilles trouvées : {xl.sheet_names}) : {excel_path}"
        )

    template_data_sheet_name, raw_data_sheet_name = xl.sheet_names

    template_data_sheet_name = pd.read_excel(  # pyright: ignore[reportUnknownMemberType]
        xl, sheet_name=template_data_sheet_name, header=None
    )
    # Lire toutes les lignes sans inférence d'en-tête
    raw = pd.read_excel(xl, sheet_name=raw_data_sheet_name, header=None)  # pyright: ignore[reportUnknownMemberType]

    # Ligne 0 : labels de groupe — forward-fill pour propager les noms sur les cellules fusionnées
    groups = raw.iloc[0].ffill()

    # Ligne 1 : noms de colonnes
    columns = raw.iloc[1]

    # Construction du MultiIndex (groupe, colonne)
    multi_columns = pd.MultiIndex.from_arrays([groups, columns])

    # Données à partir de la ligne 2
    df = raw.iloc[2:].reset_index(drop=True)
    df.columns = multi_columns

    return df
