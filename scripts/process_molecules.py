# pipeline/scripts/process_molecules.py
import pandas as pd

from scripts.columns import CID, MOLECULE


def get_unique_cids(df: pd.DataFrame) -> list[int]:
    """Collect all non-null unique CIDs from the Molecule group."""
    if MOLECULE not in df.columns.get_level_values(0):
        raise ValueError(f"Group '{MOLECULE}' not in the DataFrame.")
    col = df[MOLECULE][CID]
    return sorted(col.dropna().astype(int).unique().tolist())
