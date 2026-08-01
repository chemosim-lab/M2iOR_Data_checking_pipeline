# pipeline/scripts/process_responses.py
import pandas as pd

from scripts.columns import (
    ASSAY,
    CONCENTRATION,
    CONCENTRATION_UNIT,
    EXPERIMENTAL_TECHNIQUE,
    PARAMETER,
    RESPONSE,
    RESPONSIVE,
    STIMULATION_DURATION_UNIT,
    STIMULATION_FLUX_UNIT,
    UNIT,
    VALUE,
    VALUE_NATURE,
)

_ALLOWED_PARAMETERS = {"primary", "secondary", "dose-response"}
_ALLOWED_VALUE_NATURES = {"raw", "norm_rec", "norm_pair", "norm_other", "ec50"}
_ALLOWED_UNITS = {"na", "mol/l", "spikes/s", "v/v", "ug/ul", "um", "m", "percent"}
_ALLOWED_CONCENTRATION_UNITS = {
    "v/v",
    "m",
    "mg/ml",
    "mol/l",
    "ng/ul",
    "um",
    "ug/ul",
    "ug",
    "mg",
    "nmol",
}
_ALLOWED_STIMULATION_FLUX_UNITS = {"ml/s", "l/min", "ml/min"}
_ALLOWED_STIMULATION_DURATION_UNITS = {"ms", "s"}
_ALLOWED_EXPERIMENTAL_TECHNIQUES = {
    "two-electrode voltage clamp",
    "calcium imaging",
    "single sensillum recording",
    "fluorescence",
    "door 2.0",
    "electroantennography",
}


def _validate_allowed_values(
    df: pd.DataFrame, group: str, col_name: str, allowed: set[str]
) -> None:
    col = df[group][col_name].dropna().astype(str).str.strip().str.lower()
    invalid = col[~col.isin(allowed)]
    if not invalid.empty:
        bad = invalid.unique().tolist()
        msg = (
            f"Column '{col_name}' contains invalid value(s) {bad} "
            f"(allowed: {sorted(allowed)})"
        )
        raise ValueError(msg)


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
    """Raise ValueError if any non-null value in Value is not a float or a string."""
    col = df[RESPONSE][VALUE].dropna()
    numeric = pd.to_numeric(col, errors="coerce")
    invalid = col[numeric.isna() & ~col.astype(str).str.strip().astype(bool)]
    if not invalid.empty:
        rows = invalid.index.tolist()
        msg = f"Column '{VALUE}' contains invalid value(s) at row(s): {rows}"
        raise ValueError(msg)


def validate_concentration_column(df: pd.DataFrame) -> None:
    """Raise ValueError if any non-null value in Concentration is not a float or a string."""
    col = df[RESPONSE][CONCENTRATION].dropna()
    numeric = pd.to_numeric(col, errors="coerce")
    invalid = col[numeric.isna() & ~col.astype(str).str.strip().astype(bool)]
    if not invalid.empty:
        rows = invalid.index.tolist()
        msg = f"Column '{CONCENTRATION}' contains invalid value(s) at row(s): {rows}"
        raise ValueError(msg)


def validate_parameter_column(df: pd.DataFrame) -> None:
    _validate_allowed_values(df, RESPONSE, PARAMETER, _ALLOWED_PARAMETERS)


def validate_value_nature_column(df: pd.DataFrame) -> None:
    _validate_allowed_values(df, RESPONSE, VALUE_NATURE, _ALLOWED_VALUE_NATURES)


def validate_unit_column(df: pd.DataFrame) -> None:
    _validate_allowed_values(df, RESPONSE, UNIT, _ALLOWED_UNITS)


def validate_concentration_unit_column(df: pd.DataFrame) -> None:
    _validate_allowed_values(
        df, RESPONSE, CONCENTRATION_UNIT, _ALLOWED_CONCENTRATION_UNITS
    )


def validate_stimulation_flux_unit_column(df: pd.DataFrame) -> None:
    _validate_allowed_values(
        df, ASSAY, STIMULATION_FLUX_UNIT, _ALLOWED_STIMULATION_FLUX_UNITS
    )


def validate_stimulation_duration_unit_column(df: pd.DataFrame) -> None:
    _validate_allowed_values(
        df, ASSAY, STIMULATION_DURATION_UNIT, _ALLOWED_STIMULATION_DURATION_UNITS
    )


def validate_experimental_technique_column(df: pd.DataFrame) -> None:
    _validate_allowed_values(
        df, ASSAY, EXPERIMENTAL_TECHNIQUE, _ALLOWED_EXPERIMENTAL_TECHNIQUES
    )
