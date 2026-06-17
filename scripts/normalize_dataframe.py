# pipeline/scripts/normalize_dataframe.py

from typing import Any

import pandas as pd

from scripts.columns import (
    ASSAY,
    CONCENTRATION_UNIT,
    EXPERIMENTAL_TECHNIQUE,
    MAIN_FLUX,
    MAIN_FLUX_UNIT,
    MIXTURE,
    MOLECULE,
    ODOR_DELIVERY,
    PARAMETER,
    RESPONSE,
    SOLVENT,
    STIMULATION_DURATION_UNIT,
    STIMULATION_FLUX_UNIT,
    TYPE,
    UNIT,
    VALUE_NATURE,
)


def rename_column(df: pd.DataFrame, old_column_name: str, new_column_name: str) -> None:
    df.columns = df.rename(columns={old_column_name: new_column_name}, level=1).columns


def _strip_element(x: Any) -> Any:
    return x.strip() if isinstance(x, str) else x


def _lower_element(x: Any) -> Any:
    return x.lower() if isinstance(x, str) else x


def _strip_series(s: pd.Series[Any]) -> pd.Series[Any]:
    if not pd.api.types.is_string_dtype(s):
        return s
    result = s.map(_strip_element)
    return result.astype(s.dtype)


_LOWERCASE_COLUMNS: list[tuple[str, str]] = [
    (MOLECULE, SOLVENT),
    (MOLECULE, MIXTURE),
    (RESPONSE, PARAMETER),
    (RESPONSE, UNIT),
    (RESPONSE, VALUE_NATURE),
    (RESPONSE, CONCENTRATION_UNIT),
    (ASSAY, EXPERIMENTAL_TECHNIQUE),
    (ASSAY, TYPE),
    (ASSAY, ODOR_DELIVERY),
    (ASSAY, STIMULATION_FLUX_UNIT),
    (ASSAY, MAIN_FLUX),
    (ASSAY, MAIN_FLUX_UNIT),
    (ASSAY, STIMULATION_DURATION_UNIT),
]


def normalize_df(df: pd.DataFrame) -> pd.DataFrame:
    df = df.apply(_strip_series)
    for col in _LOWERCASE_COLUMNS:
        if col in df.columns:
            df[col] = df[col].map(_lower_element)
    return df
