# pipeline/scripts/process_responses.py
import pandas as pd

from scripts.columns import (
    CONCENTRATION,
    PARAMETER,
    RESPONSE,
    RESPONSIVE,
    VALUE,
    VALUE_NATURE,
)

_ALLOWED_PARAMETERS = {"primary", "secondary", "dose-response"}
_ALLOWED_VALUE_NATURES = {"raw", "norm_rec", "norm_pair", "norm_other", "ec50"}


def validate_responsive_column(df: pd.DataFrame) -> None:
    """Raise ValueError if any non-null value in Responsive is not 0 or 1."""
    col = df[RESPONSE][RESPONSIVE].dropna()
    numeric: pd.Series = pd.to_numeric(col, errors="coerce")
    invalid = col[~numeric.isin([0, 1]) | numeric.isna()]
    if not invalid.empty:
        rows = invalid.index.tolist()
        raise ValueError(
            f"Column '{RESPONSIVE}' contains invalid value(s) (expected 0 or 1) "
            f"at row(s): {rows}"
        )


def _validate_float_column(df: pd.DataFrame, col_name: str) -> None:
    col = df[RESPONSE][col_name].dropna()
    numeric: pd.Series = pd.to_numeric(col, errors="coerce")
    invalid = col[numeric.isna()]
    if not invalid.empty:
        rows = invalid.index.tolist()
        raise ValueError(
            f"Column '{col_name}' contains non-float value(s) at row(s): {rows}"
        )


def validate_value_column(df: pd.DataFrame) -> None:
    """Raise ValueError if any non-null value in Value is not a float."""
    _validate_float_column(df, VALUE)


def validate_concentration_column(df: pd.DataFrame) -> None:
    """Raise ValueError if any non-null value in Concentration is not a float."""
    _validate_float_column(df, CONCENTRATION)


def validate_parameter_column(df: pd.DataFrame) -> None:
    """Raise ValueError if any non-null value in Parameter is not an allowed value."""
    col = df[RESPONSE][PARAMETER].dropna().astype(str).str.strip()
    invalid = col[~col.isin(_ALLOWED_PARAMETERS)]
    if not invalid.empty:
        bad = invalid.unique().tolist()
        raise ValueError(
            f"Column '{PARAMETER}' contains invalid value(s) {bad} "
            f"(allowed: {sorted(_ALLOWED_PARAMETERS)})"
        )


def validate_value_nature_column(df: pd.DataFrame) -> None:
    """Raise ValueError if any non-null value in Value Nature is not an allowed value."""
    col = df[RESPONSE][VALUE_NATURE].dropna().astype(str).str.strip().str.lower()
    invalid = col[~col.isin(_ALLOWED_VALUE_NATURES)]
    if not invalid.empty:
        bad = invalid.unique().tolist()
        raise ValueError(
            f"Column '{VALUE_NATURE}' contains invalid value(s) {bad} "
            f"(allowed: {sorted(_ALLOWED_VALUE_NATURES)})"
        )
