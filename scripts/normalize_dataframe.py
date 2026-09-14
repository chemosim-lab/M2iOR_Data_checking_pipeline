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

# Both unicode code points used as a "micro" prefix are translated to the
# ASCII "u" used by the canonical units, so any unit built on a micro-prefixed
# base ("um", "ug", "ug/ul", ...) is accepted regardless of which micro sign
# a source Excel file used (e.g. "µg/µl" or "μg/μl" both become "ug/ul").
# Greek capital Delta ("Δ") is not listed here: `_lower_element` already
# folds it to lowercase delta ("δ") before aliasing runs.
_MICRO_SIGN_TRANSLATION = str.maketrans(
    {
        "μ": "u",  # U+03BC GREEK SMALL LETTER MU
        "µ": "u",  # U+00B5 MICRO SIGN
    }
)

# Aliases let source Excel files spell a unit/value differently while still
# resolving to the canonical form expected by the `_ALLOWED_*` sets in
# scripts/process_responses.py (e.g. "%" instead of "percent", or an ASCII
# "d" used as a stand-in for "δ"). Keys must already be lowercase, stripped
# and micro-sign-translated, since aliasing runs last. Add new entries here
# rather than widening the allowed-value sets.
_VALUE_ALIASES: dict[str, str] = {
    "%": "percent",
    "%d fluorescence": "%δ fluorescence",
    "%delta fluorescence": "%δ fluorescence",
}

# Columns eligible for alias substitution, applied after lowercasing.
_ALIAS_COLUMNS: list[tuple[str, str]] = [
    (RESPONSE, UNIT),
    (RESPONSE, CONCENTRATION_UNIT),
    (ASSAY, STIMULATION_FLUX_UNIT),
    (ASSAY, STIMULATION_DURATION_UNIT),
]


def _apply_alias(x: Any) -> Any:
    if not isinstance(x, str):
        return x
    translated = x.translate(_MICRO_SIGN_TRANSLATION)
    return _VALUE_ALIASES.get(translated, translated)


def normalize_df(df: pd.DataFrame) -> pd.DataFrame:
    df = df.apply(_strip_series)
    for col in _LOWERCASE_COLUMNS:
        if col in df.columns:
            df[col] = df[col].map(_lower_element)
    for col in _ALIAS_COLUMNS:
        if col in df.columns:
            df[col] = df[col].map(_apply_alias)
    return df
